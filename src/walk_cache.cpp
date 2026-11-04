#include "walk_cache.h"

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/parallel_sort.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <type_traits>

// Exactness notes (what makes a replay, and a resume, reproduce a walk bit for bit):
//  - one step function (runWalk) for fresh walks, replays, resumes and restarts;
//  - geometric queries in double with a fixed per-triangle function and the (d2, index) rule;
//  - a walk's random numbers come only from its own pcg32 and depend only on its path so far:
//    exactly 2 draws per move, none at the end, so the state at step k is the seed advanced 2k;
//  - no dependence on threads or scheduling: per-walk RNG, and per-point sums updated in walk
//    order by one task per point after the parallel walk passes;
//  - MSVC /fp:precise (the default) and no explicit FMA, so every expression rounds the same
//    way wherever it is evaluated.

namespace {

constexpr double kPi = 3.14159265358979323846;

double secondsSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

float roundDown(double v)
{
    float f = float(v);
    if (double(f) > v) f = std::nextafter(f, -std::numeric_limits<float>::infinity());
    return f;
}

float roundUp(double v)
{
    float f = float(v);
    if (double(f) < v) f = std::nextafter(f, std::numeric_limits<float>::infinity());
    return f;
}

bool sameBits(const void* a, const void* b, size_t n) { return std::memcmp(a, b, n) == 0; }

// uniform direction on the unit sphere from exactly two draws
Vec3 sampleSphere(pcg32& rng)
{
    const double u0 = rng.nextDouble(), u1 = rng.nextDouble();
    const double z = 1.0 - 2.0*u0, r = std::sqrt(std::max(0.0, 1.0 - z*z)), phi = 2.0*kPi*u1;
    return Vec3(r*std::cos(phi), r*std::sin(phi), z);
}

void emptyBox(WalkState& st)
{
    for (int a = 0; a < 3; a++) {
        st.lo[a] = std::numeric_limits<float>::infinity();
        st.hi[a] = -std::numeric_limits<float>::infinity();
    }
}

// adds ball(c, r) to the walk's box (c = st.x while walking; a stored centre when a box is rebuilt)
void addBall(WalkState& st, const Vec3& c, double r)
{
    for (int a = 0; a < 3; a++) {
        st.lo[a] = std::min(st.lo[a], roundDown(c[a] - r));
        st.hi[a] = std::max(st.hi[a], roundUp(c[a] + r));
    }
}

// per-walk action bits of an update (WalkCache::action)
constexpr uint8_t kTested = 1, kRewalked = 2, kValue = 4;
// stored step flags (WalkPath::f): onFace + 1 in bits 0-2, "a Dirichlet face was closer" in bit 3
constexpr uint8_t kFaceMask = 7, kFaceCloser = 8;

bool countedWalk(const WalkEntry& e) { return e.status != WalkStatus::Escaped; }

WalkState makeInitialState(const QueryPoints& pts, int nWalks, size_t walk, uint32_t rngOffset)
{
    const size_t i = walk/size_t(nWalks), j = walk%size_t(nWalks);
    WalkState st;
    st.x = pts.x[i];
    st.onFace = -1;
    st.steps = 0;
    // Walk j of point i: state from a hash of (point seed, j), stream j. A hashed state (rather than
    // pcg32(seed_i, j)) keeps the walks of a point from sharing one initial state and differing only
    // in the increment, which is known to give correlated PCG streams.
    st.rng = pcg32(mixSeed(pts.seed[i], uint64_t(j)), uint64_t(j));
    if (rngOffset) st.rng.advance(int64_t(rngOffset));
    st.w = Vec3::Zero();
    emptyBox(st);
    return st;
}

enum class WalkEnd { Landed, Capped, Escaped, Stopped, Paused };

// Control 3q (2b's queue without stored steps): a walk's 2b key, accumulated by runWalk over the steps it
// evaluates, by walkKey's rule (the min of s*_k; no more s* once the key is at or below sCur). The steps are
// bitwise the stored ones, so the key is bitwise 2b's.
struct KeyScan {
    const std::vector<Vec2>* Q = nullptr;   // the future polygon of the geometry the key is for
    double sCur = 0;                        // its front
    double key = std::numeric_limits<double>::infinity();
    long long n = 0;                        // balls whose s* was computed
};

// Coupling gate (WalkParams::oddMirrorDirichlet): the walk's first (gradient) step landed at y beyond Dirichlet
// face 'face' and went on from its mirror image my (face = -1: it did not). solveMulti scores such a walk as
// Phi(y) + Phi(my) - (its terminal value).
struct OddMirror {
    int face = -1;
    Vec3 y = Vec3::Zero(), my = Vec3::Zero();
};

struct WalkControl {
    int stopStrip = -1;   // return Stopped at the first step whose ball meets this strip
    int pauseAt = -1;     // return Paused before evaluating step pauseAt
    const std::function<void(const WalkState&, const StripCut::Closest&)>* observe = nullptr;
    WalkPath* record = nullptr;   // options 2a/2b: append every evaluated step (centre, key, flags)
    KeyScan* key = nullptr;       // control 3q: fold every evaluated step into the walk's key
    OddMirror* odd = nullptr;     // coupling gate: receives the first step's odd reflection (required with oddMirror)
};

} // namespace

// =====================================================================================
// StripCut
// =====================================================================================
StripCut::StripCut(const Tet& tet, const GeometryConfig& cfg)
{
    frame_ = buildCut(tet, cfg, 1.0);   // frame (n, u, w, d, L, sMin, sMax, fullPolygon); its polygon is unused
    initFromFrame();
}

StripCut StripCut::physicalUnionReference(const Cut& frame)
{
    StripCut c;
    c.frame_ = frame;
    c.initFromFrame();
    return c;
}

void StripCut::initFromFrame()
{
    frame_.polygon.clear();
    frame_.tris.clear();
    frame_.t = 0;
    frame_.s = frame_.sMin;
    if (frame_.fullPolygon.size() < 3) return;   // the plane misses the tet: no strip can be added

    // section in (u, w) coordinates; x is exactly the p.dot(u) that buildCut's sMin/sMax come from
    std::vector<Vec2> P;
    for (const Vec3& p : frame_.fullPolygon) P.emplace_back(p.dot(frame_.u), p.dot(frame_.w));
    double a2 = 0;
    for (size_t i = 0; i < P.size(); i++) {
        const Vec2 &p = P[i], &q = P[(i + 1)%P.size()];
        a2 += p.x()*q.y() - p.y()*q.x();
    }
    if (a2 < 0) std::reverse(P.begin(), P.end());   // safeguard: planeTetSection orders CCW about n
    fullArea_ = 0.5*std::fabs(a2);

    // Lower chain: from the lower-left vertex CCW to the lower-right one; upper chain: from the
    // upper-right CCW to the upper-left, stored by increasing u. Both are u-monotone.
    const int n = int(P.size());
    int iLL = 0, iLR = 0, iUR = 0, iUL = 0;
    for (int i = 1; i < n; i++) {
        const Vec2& p = P[i];
        if (p.x() < P[iLL].x() || (p.x() == P[iLL].x() && p.y() < P[iLL].y())) iLL = i;
        if (p.x() > P[iLR].x() || (p.x() == P[iLR].x() && p.y() < P[iLR].y())) iLR = i;
        if (p.x() > P[iUR].x() || (p.x() == P[iUR].x() && p.y() > P[iUR].y())) iUR = i;
        if (p.x() < P[iUL].x() || (p.x() == P[iUL].x() && p.y() > P[iUL].y())) iUL = i;
    }
    for (int i = iLL;; i = (i + 1)%n) {
        lower_.push_back(P[i]);
        if (i == iLR) break;
    }
    for (int i = iUR;; i = (i + 1)%n) {
        upper_.push_back(P[i]);
        if (i == iUL) break;
    }
    std::reverse(upper_.begin(), upper_.end());
}

double StripCut::sOf(double t) const
{
    // the same expression as buildCut, so the fronts of both agree bit for bit
    t = std::clamp(t, 0.0, 1.0);
    return t >= 1.0 ? frame_.sMax : frame_.sMin + t*(frame_.sMax - frame_.sMin);
}

double StripCut::tOf(double s) const
{
    const double ext = frame_.sMax - frame_.sMin;
    return ext > 0 ? std::clamp((s - frame_.sMin)/ext, 0.0, 1.0) : 0.0;
}

// The two points of the front line X·u = s on the section boundary (lower and upper chain).
// One function of s for every strip, so the strips on both sides of a front share these bits.
void StripCut::frontVertices(double s, Vec2& lo, Vec2& hi) const
{
    auto onChain = [s](const std::vector<Vec2>& C) -> double {
        for (size_t k = 0; k + 1 < C.size(); k++) {
            const Vec2 &a = C[k], &b = C[k + 1];
            if (s == a.x()) return a.y();
            if (s > a.x() && s < b.x()) return a.y() + (s - a.x())/(b.x() - a.x())*(b.y() - a.y());
        }
        return s <= C.front().x() ? C.front().y() : C.back().y();
    };
    lo = Vec2(s, onChain(lower_));
    hi = Vec2(s, onChain(upper_));
}

int StripCut::addStrip(double sNew)
{
    if (lower_.empty()) return -1;
    sNew = std::min(sNew, frame_.sMax);
    const double s0 = front();
    if (!(sNew > s0)) return -1;

    // the strip polygon, CCW: lower chain left to right, then upper chain right to left
    Vec2 lo0, hi0, lo1, hi1;
    frontVertices(s0, lo0, hi0);
    frontVertices(sNew, lo1, hi1);
    std::vector<Vec2> poly;
    poly.push_back(lo0);
    for (const Vec2& v : lower_)
        if (v.x() > s0 && v.x() < sNew) poly.push_back(v);
    poly.push_back(lo1);
    poly.push_back(hi1);
    for (auto it = upper_.rbegin(); it != upper_.rend(); ++it)
        if (it->x() > s0 && it->x() < sNew) poly.push_back(*it);
    poly.push_back(hi0);
    for (size_t i = 0; i < poly.size() && poly.size() > 1;) {   // exact duplicates (a front through a vertex)
        if (poly[i] == poly[(i + 1)%poly.size()]) poly.erase(poly.begin() + i);
        else i++;
    }

    Strip S;
    S.s0 = s0;
    S.s1 = sNew;
    S.triBegin = int(tris_.size());
    S.wMin = std::numeric_limits<double>::max();
    S.wMax = -S.wMin;
    S.lo = Vec3::Constant(std::numeric_limits<double>::max());
    S.hi = -S.lo;
    for (const Vec2& v : poly) {
        S.wMin = std::min(S.wMin, v.y());
        S.wMax = std::max(S.wMax, v.y());
        Vec3 X = toWorld(v);
        S.lo = S.lo.cwiseMin(X);
        S.hi = S.hi.cwiseMax(X);
    }
    // fan triangulation with the existing rules: zero-area triangles skipped, CCW about n
    const double L = frame_.L;
    for (int i = 1; i + 1 < int(poly.size()); i++) {
        Tri T;
        T.v[0] = poly[0];
        T.v[1] = poly[i];
        T.v[2] = poly[i + 1];
        Vec2 e1 = T.v[1] - T.v[0], e2 = T.v[2] - T.v[0];
        double cr = e1.x()*e2.y() - e1.y()*e2.x();
        if (std::fabs(cr) <= 1e-12*L*L) continue;
        if (cr < 0) std::swap(T.v[1], T.v[2]);
        for (int k = 0; k < 3; k++) {
            T.e[k] = T.v[(k + 1)%3] - T.v[k];
            T.invLen2[k] = 1.0/T.e[k].squaredNorm();
        }
        T.strip = int(strips_.size());
        S.area += 0.5*std::fabs(cr);
        tris_.push_back(T);
    }
    S.triEnd = int(tris_.size());
    strips_.push_back(S);
    return int(strips_.size()) - 1;
}

void StripCut::retractTo(double sNew)
{
    while (!strips_.empty() && strips_.back().s1 > sNew) {
        tris_.resize(size_t(strips_.back().triBegin));
        strips_.pop_back();
    }
    addStrip(sNew);
}

double StripCut::area(int nStrips) const
{
    double a = 0;
    for (int m = 0; m < nStrips && m < int(strips_.size()); m++) a += strips_[m].area;
    return a;
}

Cut StripCut::frontCut(int nStrips) const
{
    Cut c = frame_;
    c.s = front(nStrips);
    c.t = tOf(c.s);
    return c;
}

// In-plane squared distance from q to triangle T, and the closest point c. Only edges that q lies
// outside of can hold the closest point of a convex polygon, so the others are skipped.
double StripCut::triDist2(const Tri& T, const Vec2& q, Vec2& c) const
{
    bool inside = true;
    double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 3; i++) {
        const Vec2 r = q - T.v[i];
        const Vec2& e = T.e[i];
        if (e.x()*r.y() - e.y()*r.x() >= 0) continue;
        inside = false;
        const double s = std::clamp(r.dot(e)*T.invLen2[i], 0.0, 1.0);
        const Vec2 p = T.v[i] + s*e;
        const double d2 = (q - p).squaredNorm();
        if (d2 < best) {
            best = d2;
            c = p;
        }
    }
    if (inside) {
        c = q;
        return 0.0;
    }
    return best;
}

void StripCut::evalStrip(int m, double h2, const Vec2& q, Closest& best) const
{
    for (int t = strips_[m].triBegin; t < strips_[m].triEnd; t++) {
        Vec2 c;
        const double d2 = h2 + triDist2(tris_[t], q, c);
        if (d2 < best.d2 || (d2 == best.d2 && t < best.tri)) {
            best.d2 = d2;
            best.tri = t;
            best.q = c;
        }
    }
}

// Lower bound on the key h^2 + e2 of every triangle of strip m: the strip lies in the box
// [s0, s1] x [wMin, wMax] of the plane. The gaps are reduced by 1e-9 L, far more than the round-off
// of any computed e2, so a strip whose bound exceeds the best key cannot hold the minimum.
double StripCut::stripLowerBound2(int m, double h2, const Vec2& q) const
{
    const Strip& S = strips_[m];
    const double slack = 1e-9*frame_.L;
    double gu = std::max(S.s0 - q.x(), q.x() - S.s1), gw = std::max(S.wMin - q.y(), q.y() - S.wMax);
    gu = std::max(0.0, gu - slack);
    gw = std::max(0.0, gw - slack);
    return h2 + (gu*gu + gw*gw);
}

StripCut::Closest StripCut::closest(const Vec3& x, int nStrips) const
{
    Closest best;
    nStrips = std::min(nStrips, int(strips_.size()));
    if (nStrips <= 0) return best;
    const double h = frame_.side(x), h2 = h*h;
    const Vec2 q(x.dot(frame_.u), x.dot(frame_.w));
    int first = 0;
    double lbFirst = std::numeric_limits<double>::infinity();
    for (int m = 0; m < nStrips; m++) {
        double lb = stripLowerBound2(m, h2, q);
        if (lb < lbFirst) {
            lbFirst = lb;
            first = m;
        }
    }
    evalStrip(first, h2, q, best);
    for (int m = 0; m < nStrips; m++)
        if (m != first && stripLowerBound2(m, h2, q) <= best.d2) evalStrip(m, h2, q, best);
    return best;
}

StripCut::Closest StripCut::closestBruteForce(const Vec3& x, int nStrips) const
{
    Closest best;
    nStrips = std::min(nStrips, int(strips_.size()));
    const double h = frame_.side(x), h2 = h*h;
    const Vec2 q(x.dot(frame_.u), x.dot(frame_.w));
    for (int m = 0; m < nStrips; m++) evalStrip(m, h2, q, best);
    return best;
}

double StripCut::stripDist2(const Vec3& x, int m) const
{
    Closest best;
    const double h = frame_.side(x);
    evalStrip(m, h*h, Vec2(x.dot(frame_.u), x.dot(frame_.w)), best);
    return best.d2;
}

bool StripCut::boxMayHitStrip(const float lo[3], const float hi[3], int m) const
{
    const Strip& S = strips_[m];
    if (S.triBegin == S.triEnd) return false;
    const double pad = 1e-9*frame_.L;   // covers every double round-off between the box, the balls and the strip
    const Vec3 bl(lo[0], lo[1], lo[2]), bh(hi[0], hi[1], hi[2]);
    for (int a = 0; a < 3; a++)   // axes x, y, z
        if (bl[a] > S.hi[a] + pad || bh[a] < S.lo[a] - pad) return false;
    const Vec3 c = 0.5*(bl + bh), e = 0.5*(bh - bl);
    auto radius = [&e](const Vec3& axis) { return std::fabs(axis.x())*e.x() + std::fabs(axis.y())*e.y() + std::fabs(axis.z())*e.z(); };
    // the cut normal: the strip lies in the plane n·X = d
    if (std::fabs(frame_.n.dot(c) - frame_.d) > radius(frame_.n) + pad) return false;
    // u: the strip's extent is exactly [s0, s1]
    double r = radius(frame_.u), cu = frame_.u.dot(c);
    if (cu - r > S.s1 + pad || cu + r < S.s0 - pad) return false;
    // w
    r = radius(frame_.w);
    double cw = frame_.w.dot(c);
    if (cw - r > S.wMax + pad || cw + r < S.wMin - pad) return false;
    return true;
}

void StripCut::mesh(int nStrips, std::vector<Vec3>& V, std::vector<std::array<int, 3>>& F,
                    std::vector<double>& triStrip) const
{
    V.clear();
    F.clear();
    triStrip.clear();
    nStrips = std::min(nStrips, int(strips_.size()));
    for (int m = 0; m < nStrips; m++)
        for (int t = strips_[m].triBegin; t < strips_[m].triEnd; t++) {
            int b = int(V.size());
            for (int k = 0; k < 3; k++) V.push_back(toWorld(tris_[t].v[k]));
            F.push_back({b, b + 1, b + 2});
            triStrip.push_back(double(m));
        }
}

// Option 2b. P ∩ {X·u >= s}, s = front(nStrips), built like a strip without a right end: the lower front
// vertex, the lower chain beyond s, the upper chain beyond s (right to left), the upper front vertex.
std::vector<Vec2> StripCut::futurePolygon(int nStrips) const
{
    std::vector<Vec2> Q;
    if (lower_.empty()) return Q;
    const double s = front(nStrips);
    Vec2 lo, hi;
    frontVertices(s, lo, hi);
    Q.push_back(lo);
    for (const Vec2& v : lower_)
        if (v.x() > s) Q.push_back(v);
    for (auto it = upper_.rbegin(); it != upper_.rend(); ++it)
        if (it->x() > s) Q.push_back(*it);
    Q.push_back(hi);
    for (size_t i = 0; i < Q.size() && Q.size() > 1;) {   // exact duplicates (a single vertex at an end)
        if (Q[i] == Q[(i + 1)%Q.size()]) Q.erase(Q.begin() + i);
        else i++;
    }
    return Q;
}

// Option 2b. The update's exact test asks whether the new strip comes closer to the step centre x than the
// step's key: fl(h2 + e2) < key2 for some triangle of the strip (h = n·x - d, h2 = h*h, e2 = in-plane squared
// distance from q = x in (u, w)). A strip passing it has a point X with |X - q|^2 < key2 - h2, up to round-off.
// The radius a of that disk is enlarged to cover the round-off many times over: the cancellation in
// key2 - h2 (a few ulps of key2; it matters when the ball barely reaches the plane, and 1e-9 key2 covers it
// 1e6 times), and the ~1e-16 L error of e2 and of the strips' interpolated front vertices (1e-9 L). With
// h2 >= key2 (1 + 1e-9) no strip can pass (fl(h2 + e2) >= h2: rounding is monotone and e2 >= 0): +inf.
// Every later strip lies in Q and has X·u <= its front s_new, so it can pass only if s_new >= the minimum of
// X·u over Q ∩ disk(q, a). That minimum is attained at the disk's leftmost point (if it lies in Q), at a
// vertex of Q inside the disk, or at an end of an edge's chord through the disk. The margin on a also keeps
// the chord ends well conditioned (a chord the uncertain true disk barely touches is >= sqrt(2 a 1e-9 L)
// long in the enlarged one), and the result is lowered by another 1e-9 L for the candidates' own round-off.
// So s* is conservative: a larger margin only adds false alarms (popped walks that the exact test rejects).
double StripCut::entryFront(const std::vector<Vec2>& Q, const Vec3& x, double key2, double bound) const
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    const double h = frame_.side(x), h2 = h*h;
    const double a2raw = key2 - h2 + 1e-9*key2;
    if (Q.empty() || !(a2raw > 0)) return inf;
    const double L = frame_.L, tol = 1e-9*L;
    const double a = std::sqrt(a2raw) + tol, a2 = a*a;
    const Vec2 q(x.dot(frame_.u), x.dot(frame_.w));
    // every candidate is >= q.x() - a (up to round-off); Q lies in X·u >= Q[0].x(), the current front
    if (q.x() - a - 2*tol >= bound || q.x() + a < Q[0].x()) return inf;
    const int n = int(Q.size());
    const Vec2 left(q.x() - a, q.y());
    bool leftInside = n >= 3;
    double best = inf;
    for (int i = 0; i < n; i++) {
        const Vec2 &A = Q[i], &B = Q[(i + 1)%n];
        if ((A - q).squaredNorm() <= a2) best = std::min(best, A.x());
        const Vec2 e = B - A, r = q - A, rl = left - A;
        const double len2 = e.squaredNorm();
        if (len2 == 0) continue;
        if (e.x()*rl.y() - e.y()*rl.x() < -tol*std::sqrt(len2)) leftInside = false;   // CCW: inside = left of every edge (lenient)
        const double cr = e.x()*r.y() - e.y()*r.x(), dperp2 = cr*cr/len2;              // q to the edge's line, squared
        if (dperp2 >= a2) continue;
        const double tf = r.dot(e)/len2, ht = std::sqrt((a2 - dperp2)/len2);         // foot of q, half chord (edge units)
        const double t0 = std::max(0.0, tf - ht), t1 = std::min(1.0, tf + ht);
        if (t0 <= t1) best = std::min({best, A.x() + t0*e.x(), A.x() + t1*e.x()});
    }
    if (leftInside) best = std::min(best, left.x());
    best -= tol;
    return best < bound ? best : inf;
}

// =====================================================================================
// Sums and estimates
// =====================================================================================
void QueryPoints::add(const std::vector<Vec3>& pts, const std::vector<char>* gradFlags, uint64_t seedBase)
{
    for (size_t i = 0; i < pts.size(); i++) {
        x.push_back(pts[i]);
        grad.push_back(gradFlags ? (*gradFlags)[i] : 0);
        seed.push_back(mixSeed(seedBase, uint64_t(i)));
    }
}

void PointSums::add(double u, const Vec3& w, double sign)
{
    n += sign > 0 ? 1 : -1;
    sH += sign*u;
    sH2 += sign*(u*u);
    for (int a = 0; a < 3; a++) {
        const double wa = w[a], wa2 = wa*wa;
        sw[a] += sign*wa;
        swu[a] += sign*(wa*u);
        sw2[a] += sign*wa2;
        sw2u[a] += sign*(wa2*u);
        sw2u2[a] += sign*(wa2*u*u);
    }
}

void PointSums::add(const PointSums& o)
{
    n += o.n;
    sH += o.sH;
    sH2 += o.sH2;
    for (int a = 0; a < 3; a++) {
        sw[a] += o.sw[a];
        swu[a] += o.swu[a];
        sw2[a] += o.sw2[a];
        sw2u[a] += o.sw2u[a];
        sw2u2[a] += o.sw2u2[a];
    }
}

Estimate EstimateD::toFloat() const
{
    Estimate e;
    e.H = float(H);
    e.Hse = float(Hse);
    e.grad = grad.cast<float>();
    e.gradSe = gradSe.cast<float>();
    e.hasGrad = hasGrad;
    e.nCounted = int(n);
    e.walkVar = float(walkVar);
    return e;
}

EstimateD estimateFromSums(const PointSums& s, bool grad)
{
    EstimateD e;
    e.n = s.n;
    e.hasGrad = grad && s.n > 0;
    if (s.n <= 0) return e;
    const double n = double(s.n);
    e.H = s.sH/n;
    e.walkVar = s.n > 1 ? std::max(0.0, (s.sH2 - n*e.H*e.H)/(n - 1)) : 0.0;
    e.Hse = std::sqrt(e.walkVar/n);
    if (e.hasGrad) {
        for (int a = 0; a < 3; a++) {
            const double g = (s.swu[a] - e.H*s.sw[a])/n;
            const double sy2 = s.sw2u2[a] - 2.0*e.H*s.sw2u[a] + e.H*e.H*s.sw2[a];   // Σ (w_a (u - H))^2
            const double var = s.n > 1 ? std::max(0.0, (sy2 - n*g*g)/(n - 1)) : 0.0;
            e.grad[a] = g;
            e.gradSe[a] = std::sqrt(var/n);
        }
    }
    return e;
}

const char* updateModeName(UpdateMode m)
{
    switch (m) {
    case UpdateMode::Incremental: return "incremental";
    case UpdateMode::FullResolve: return "full re-solve";
    case UpdateMode::BiasedRestart: return "biased restart";
    case UpdateMode::StoredBox: return "stored steps + bbox (2a)";
    case UpdateMode::StoredQueue: return "stored steps + front-key queue (2b)";
    case UpdateMode::ReplayQueue: return "seed replay + front-key queue (3q)";
    }
    return "?";
}

std::string ReuseStats::summary() const
{
    char buf[640];
    if (!full && stored) {   // options 2a/2b
        const double n = double(std::max(1LL, nWalks));
        std::snprintf(buf, sizeof buf,
                      "%lld walks: candidates %.2f%%, re-walked %.2f%% (%.1f%% of them from k = 0), value-updated %.2f%%; "
                      "%.3f s (filter %.3f, exact test %.3f, resume %.3f, new steps %.3f, keys %.3f, values %.3f, sums %.3f, "
                      "queue %.3f); %lld balls tested, %lld new steps, %lld keys",
                      nWalks, 100.0*nTested/n, 100.0*nRewalked/n,
                      nRewalked ? 100.0*double(nRewalkedAtStep0)/double(nRewalked) : 0.0, 100.0*nValueUpdated/n, tTotal,
                      tBox, tTest, tResume, tRewalk, tKeys, tValue, tSums, tQueue, ballsTested, stepsNew, keysComputed);
    } else if (full) {
        std::snprintf(buf, sizeof buf, "full solve of %lld walks: %.3f s (%lld escaped, %lld capped, %lld steps)", nWalks,
                      tTotal, escaped, capped, stepsNew);
    } else {
        const double n = double(std::max(1LL, nWalks));
        std::snprintf(buf, sizeof buf,
                      "%lld walks: tested (bbox hit) %.2f%%, re-walked %.2f%% (%.1f%% of them from k = 0), "
                      "value-updated %.2f%%; %.3f s (bbox %.3f, re-walk %.3f, values %.3f, sums %.3f); "
                      "re-walked paths: %lld escaped, %lld capped; steps replayed %lld, new %lld",
                      nWalks, 100.0*nTested/n, 100.0*nRewalked/n,
                      nRewalked ? 100.0*double(nRewalkedAtStep0)/double(nRewalked) : 0.0, 100.0*nValueUpdated/n, tTotal,
                      tBox, tRewalk, tValue, tSums, escaped, capped, stepsReplayed, stepsNew);
    }
    return buf;
}

bool sameWalk(const WalkEntry& a, const Vec3* wa, const WalkEntry& b, const Vec3* wb, bool ignoreRngOffset)
{
    if (!sameBits(a.lo, b.lo, sizeof a.lo) || !sameBits(a.hi, b.hi, sizeof a.hi)) return false;
    if (!sameBits(a.land.data(), b.land.data(), 3*sizeof(double))) return false;
    if (!sameBits(&a.value, &b.value, sizeof(double))) return false;
    if (a.steps != b.steps || a.side != b.side || a.status != b.status || a.face != b.face) return false;
    if (!ignoreRngOffset && a.rngOffset != b.rngOffset) return false;
    if ((wa == nullptr) != (wb == nullptr)) return false;
    if (wa && !sameBits(wa->data(), wb->data(), 3*sizeof(double))) return false;
    return true;
}

// =====================================================================================
// The step function
// =====================================================================================
// The domain's faces as planes n·x = c with outward unit normals (n·x <= c inside): a tet's four, in Tet::faces
// order, or the nf <= kMaxDomainPlanes faces of a convex polyhedron (the coupling gate's ring reference,
// ConvexDomain). A tet runs the same loops over its four planes, in the same order, as before nf existed.
struct FacePlanes {
    int nf = 4;
    Vec3 n[kMaxDomainPlanes];
    double c[kMaxDomainPlanes] = {0, 0, 0, 0, 0, 0, 0, 0};
    FacePlanes() = default;
    explicit FacePlanes(const Tet& tet)
    {
        for (int f = 0; f < 4; f++) {
            const auto& F = tet.faces[f];
            n[f] = (tet.v[F[1]] - tet.v[F[0]]).cross(tet.v[F[2]] - tet.v[F[0]]).normalized();
            c[f] = n[f].dot(tet.v[F[0]]);
        }
    }
    explicit FacePlanes(const ConvexDomain& dom) : nf(int(dom.n.size()))
    {
        if (nf < 4 || nf > kMaxDomainPlanes || dom.c.size() != dom.n.size())
            throw std::invalid_argument("FacePlanes: a convex domain needs 4 to kMaxDomainPlanes face planes");
        for (int f = 0; f < nf; f++) {
            n[f] = dom.n[f];
            c[f] = dom.c[f];
        }
    }
    double distInside(const Vec3& x) const
    {
        double dmin = std::numeric_limits<double>::max();
        for (int f = 0; f < nf; f++) dmin = std::min(dmin, c[f] - n[f].dot(x));
        return dmin;
    }
    bool outside(const Vec3& x, double tol) const
    {
        for (int f = 0; f < nf; f++)
            if (n[f].dot(x) - c[f] > tol) return true;
        return false;
    }
};

// Fixed inputs of the step function for one geometry version (strips [0, nStrips)).
struct WalkCtx {
    const StripCut* cut = nullptr;
    int nStrips = 0;
    FacePlanes faces;
    Cut front;              // frame with this version's front: Dirichlet data side·g(landing)
    GType g = GType::One;
    double shell = 0;       // terminate when the distance to the cut is below this
    double escTol = 0;      // a position farther than this outside a face plane is an escape
    int maxLen = 0;
    // Faces that absorb: all four with H = 0 in the experiment (WalkParams::dirichletFaces), the flagged
    // ones with faceData in the coupling gate (WalkParams::faceBC); none everywhere else (and none on a ConvexDomain).
    bool dirichletFace[kMaxDomainPlanes] = {false, false, false, false, false, false, false, false};
    bool anyDirichletFace = false;
    std::function<double(int, const Vec3&)> faceData;   // empty: H = 0 on the Dirichlet faces
    bool cutData = true;    // false: H = 0 on the cut (WalkParams::cutData)
    bool mirrorGradBall = false;   // WalkParams::mirrorGradBall (coupling gate)
    bool oddMirror = false;        // WalkParams::oddMirrorDirichlet (coupling gate, with mirrorGradBall)
    float boxLo[3], boxHi[3];   // the tet's bounding box (padded, outward): walk boxes are clipped to it

    WalkCtx(const Tet& tet, const StripCut& c, int nStrips_, GType g_, const WalkParams& wp)
        : cut(&c), nStrips(nStrips_), faces(tet), front(c.frontCut(nStrips_)), g(g_)
    {
        for (int f = 0; f < 4; f++) {
            dirichletFace[f] = wp.dirichletFaces || wp.faceBC[f] == FaceBC::Dirichlet;
            anyDirichletFace = anyDirichletFace || dirichletFace[f];
        }
        if (!wp.dirichletFaces) faceData = wp.faceData;
        cutData = wp.cutData;
        mirrorGradBall = wp.mirrorGradBall;
        oddMirror = wp.mirrorGradBall && wp.oddMirrorDirichlet;
        shell = wp.absorbShellRel*tet.L;
        escTol = 1e-12*tet.L;
        maxLen = wp.maxWalkLength;
        Vec3 lo = tet.v[0], hi = tet.v[0];
        for (const Vec3& p : tet.v) {
            lo = lo.cwiseMin(p);
            hi = hi.cwiseMax(p);
        }
        for (int a = 0; a < 3; a++) {
            boxLo[a] = roundDown(lo[a] - 1e-9*tet.L);
            boxHi[a] = roundUp(hi[a] + 1e-9*tet.L);
        }
    }
    // Coupling gate, check 4 only: the ring reference on a convex polyhedron (the octahedron), every face Neumann.
    // The walk is the tet's walk with the domain's face planes: the polyhedron is convex, so seen from inside no
    // point of its boundary is a silhouette and the star radius is the distance to the cut (starStep).
    WalkCtx(const ConvexDomain& dom, const StripCut& c, int nStrips_, GType g_, const WalkParams& wp)
        : cut(&c), nStrips(nStrips_), faces(dom), front(c.frontCut(nStrips_)), g(g_)
    {
        if (wp.dirichletFaces || wp.anyDirichletFace() || wp.faceData || wp.oddMirrorDirichlet)
            throw std::invalid_argument("WalkCtx: a ConvexDomain has Neumann faces only");
        cutData = wp.cutData;
        mirrorGradBall = wp.mirrorGradBall;
        shell = wp.absorbShellRel*dom.L;
        escTol = 1e-12*dom.L;
        maxLen = wp.maxWalkLength;
        Vec3 lo = dom.v[0], hi = dom.v[0];
        for (const Vec3& p : dom.v) {
            lo = lo.cwiseMin(p);
            hi = hi.cwiseMax(p);
        }
        for (int a = 0; a < 3; a++) {
            boxLo[a] = roundDown(lo[a] - 1e-9*dom.L);
            boxHi[a] = roundUp(hi[a] + 1e-9*dom.L);
        }
    }

    // Distance from x (inside) to the nearest Dirichlet face plane, and that face; +inf, -1 without one.
    // With all four faces Dirichlet this is bitwise FacePlanes::distInside (same terms, exact min).
    double dirichletFaceDist(const Vec3& x, int& face) const
    {
        face = -1;
        if (!anyDirichletFace) return std::numeric_limits<double>::infinity();
        double dmin = std::numeric_limits<double>::max();
        for (int f = 0; f < faces.nf; f++) {
            if (!dirichletFace[f]) continue;
            const double d = faces.c[f] - faces.n[f].dot(x);
            if (d < dmin) {
                dmin = d;
                face = f;
            }
        }
        return dmin;
    }
    // x projected onto the plane of face f
    Vec3 onFacePlane(int f, const Vec3& x) const { return x - (faces.n[f].dot(x) - faces.c[f])*faces.n[f]; }
};

namespace {

void finishWalk(const WalkCtx& c, const WalkState& st, WalkStatus status, const StripCut::Closest* q, WalkEntry* out,
                int face = -1)
{
    if (!out) return;
    for (int a = 0; a < 3; a++) {   // the strip lies in the tet, so clipping to the tet's box loses no hit
        out->lo[a] = std::max(st.lo[a], c.boxLo[a]);
        out->hi[a] = std::min(st.hi[a], c.boxHi[a]);
    }
    out->steps = st.steps;
    out->status = status;
    out->face = -1;
    if (status == WalkStatus::Landed && !q) {   // landed on Dirichlet face 'face' (experiment: H = 0 there)
        out->land = st.x;
        out->side = 0;
        out->face = int8_t(face);
        out->value = c.faceData ? c.faceData(face, c.onFacePlane(face, st.x)) : 0.0;
    } else if (status == WalkStatus::Landed) {
        out->land = c.cut->toWorld(q->q);
        // the side is decided at the walk's position, never at the projected point; a tie counts as +
        out->side = c.front.side(st.x) >= 0 ? 1 : -1;
        out->value = c.cutData ? double(out->side)*c.front.g(out->land, c.g) : 0.0;
    } else {
        out->land = st.x;
        out->side = 0;
        out->value = 0;   // capped: counted with value 0 (as the first gate); escaped: not counted
    }
}

// One step of Walk on Stars from a point x with distance R to the cut. The tet is convex, so seen
// from inside no point of its boundary is a silhouette: the star-shaped region is the whole
// ball(x, R) ∩ tet, and its radius is R itself (Zombie's computeStarRadiusForReflectingBoundary
// would return R here too, before its 0.99 shrink). On a face, the direction is flipped into the
// inward hemisphere (a flip, so the draw count stays 2). The walk goes to the first face plane
// along the direction if that is closer than R (and then lies on that face), else to the sphere.
// The same holds on any convex polyhedron (F.nf planes: the ring reference of the coupling gate).
void starStep(const FacePlanes& F, double R, Vec3 dir, WalkState& st)
{
    if (st.onFace >= 0 && F.n[st.onFace].dot(dir) > 0) dir = -dir;
    double tHit = R;
    int hit = -1;
    for (int f = 0; f < F.nf; f++) {
        if (f == st.onFace) continue;
        const double nd = F.n[f].dot(dir);
        if (nd <= 0) continue;
        const double t = std::max(0.0, (F.c[f] - F.n[f].dot(st.x))/nd);
        if (t < tHit) {
            tHit = t;
            hit = f;
        }
    }
    if (hit >= 0) {
        Vec3 y = st.x + tHit*dir;
        y = y - (F.n[hit].dot(y) - F.c[hit])*F.n[hit];   // onto the face plane (up to round-off)
        st.x = y;
        st.onFace = hit;
    } else {
        st.x = st.x + R*dir;
        st.onFace = -1;
    }
}

// Coupling gate (WalkParams::mirrorGradBall): the gradient's first ball at x may reach across the nearest
// Neumann face f into the mirror image M(T) of the tet. With zero flux through f, the even extension
// H~(y) = H(M y) of H to T u M(T) is harmonic across the open face, so the mean-value and Poisson-kernel
// identities hold on any ball B(x, R) inside T u M(T) that meets neither the cut nor its mirror image nor any
// other face plane: then B(x, R) n T and M(B(x, R) \ T) lie in T, and the disk B n plane(f) lies inside the
// face f. Sufficient: x' = M x inside every other face plane g, and
//   R <= min(dist(x, cut), dist(x', cut), dist(x, plane g), dist(x', plane g) for all g != f)
// (the Dirichlet faces are among the planes g). A landing point beyond f is reflected to M y; the weight
// stays 3 dir/R with the unreflected dir, since grad H(x) = grad H~(x). Next to a Neumann face the plain
// ball (radius ~ the face distance delta) gives a gradient variance ~ 1/delta^2, which no volume
// quadrature can integrate; the mirror ball's radius is set by the other boundaries. Returns f and sets
// R0 when the mirror ball (shrunk by kGradFirstBallShrink like the plain one) is larger, else -1.
// With c.oddMirror the face f may also be a Dirichlet face with data phi, under the same conditions: then
// H - Phi (Phi harmonic in all of space, Phi = phi on f) vanishes on the disk B n plane(f), which lies inside
// f, so its odd extension is harmonic in B (Schwarz reflection), and so is H~ = Phi + that extension:
// H~(y) = Phi(y) + Phi(M y) - H(M y) beyond f. A walk that lands beyond f goes on from M y, and solveMulti
// scores it as Phi(y) + Phi(M y) - (its terminal value) (OddMirror). Next to a Dirichlet face the plain
// ball's gradient variance grows like 1/delta (a walk from distance delta leaves the face's neighbourhood
// with probability ~ delta/l, weight 3/delta), heavy-tailed, and the variance of a split product, integrated
// over the points, is set by the few points closest to the face.
// On a convex polyhedron (ConvexDomain, Neumann faces only) the same argument holds with its nf planes.
int mirrorGradBall(const WalkCtx& c, const Vec3& x, double dCut, double& R0)
{
    int f = -1;
    double df = std::numeric_limits<double>::infinity();
    for (int g = 0; g < c.faces.nf; g++) {
        if (c.dirichletFace[g] && !c.oddMirror) continue;
        const double dg = c.faces.c[g] - c.faces.n[g].dot(x);
        if (dg < df) {
            df = dg;
            f = g;
        }
    }
    if (f < 0) return -1;
    const Vec3 xm = x + 2.0*df*c.faces.n[f];
    double R = std::min(dCut, std::sqrt(c.cut->closest(xm, c.nStrips).d2));
    for (int g = 0; g < c.faces.nf; g++) {
        if (g == f) continue;
        const double a = c.faces.c[g] - c.faces.n[g].dot(x), b = c.faces.c[g] - c.faces.n[g].dot(xm);
        if (!(b > 0)) return -1;
        R = std::min(R, std::min(a, b));
    }
    R *= kGradFirstBallShrink;
    if (!(R > R0)) return -1;
    R0 = R;
    return f;
}

// THE step function, used by fresh walks, replays, resumes and restarts alike. Continues the walk in
// st under the strips [0, c.nStrips) until it lands, hits the cap or escapes, and fills *out.
// Step k: the distance d_k to the cut decides everything (termination, landing point, radius), and
// ball(center_k, d_k) is added to the walk's box. With ctl.stopStrip >= 0 the walk is first tested
// at each step: if that strip comes closer than d_k (same key, so exactly when adding the strip
// would change d_k or the closest point) it returns Stopped with st left at that step untouched.
WalkEnd runWalk(const WalkCtx& c, bool gradPoint, const WalkControl& ctl, WalkState& st, WalkEntry* out)
{
    for (;;) {
        if (st.steps == ctl.pauseAt) return WalkEnd::Paused;
        if (st.steps >= c.maxLen) {   // no distance is evaluated at the capped position
            finishWalk(c, st, WalkStatus::Capped, nullptr, out);
            return WalkEnd::Capped;
        }
        const StripCut::Closest q = c.cut->closest(st.x, c.nStrips);
        if (ctl.observe) (*ctl.observe)(st, q);
        // Dirichlet faces (the dirichletFaces experiment, or the coupling gate's flagged faces) absorb too,
        // so they also limit the ball. When a face is closer than the cut, the step changes only if the
        // strip comes closer than that face; the 1e-12 relative margin keeps the test conservative against
        // the rounding of dF^2 vs sqrt. All faces Neumann: dF = +inf, the first gate's walk.
        int fF = -1;
        const double dF = c.dirichletFaceDist(st.x, fF);
        const double dCut = std::sqrt(q.d2);
        const bool faceCloser = dF < dCut;
        const double key2 = faceCloser ? dF*dF*(1 + 1e-12) : q.d2;
        if (ctl.stopStrip >= 0 && c.cut->stripDist2(st.x, ctl.stopStrip) < key2) return WalkEnd::Stopped;
        if (ctl.record) {   // options 2a/2b: this step's state and key, as they are now
            ctl.record->s.push_back(StoredStep{st.x, key2});
            ctl.record->f.push_back(uint8_t(uint8_t(st.onFace + 1) | (faceCloser ? kFaceCloser : 0)));
        }
        if (ctl.key && ctl.key->key > ctl.key->sCur) {   // 3q: this step's s*, exactly as walkKey computes it
            ctl.key->n++;
            ctl.key->key = std::min(ctl.key->key, c.cut->entryFront(*ctl.key->Q, st.x, key2, ctl.key->key));
        }
        const double d = faceCloser ? dF : dCut;
        addBall(st, st.x, d);
        if (d < c.shell) {
            finishWalk(c, st, WalkStatus::Landed, faceCloser ? nullptr : &q, out, fF);
            return WalkEnd::Landed;
        }
        const Vec3 dir = sampleSphere(st.rng);
        if (gradPoint && st.steps == 0) {
            // Gradient points: step 0 is the full ball of radius R0 (inside the tet and away from
            // the cut), landing at x + R0 dir. With uniform sphere sampling the Poisson kernel's
            // gradient divided by the sample density is 3 dir/R0, so w (u - c) is an estimate of
            // ∇H for any constant c. w is cached per walk: a value update changes u but not w.
            double R0 = kGradFirstBallShrink*std::min(d, c.faces.distInside(st.x));
            const int mirror = c.mirrorGradBall ? mirrorGradBall(c, st.x, dCut, R0) : -1;
            st.w = (3.0/R0)*dir;
            st.x = st.x + R0*dir;
            if (mirror >= 0) {   // a landing point beyond the mirror face is reflected back into the tet
                const double h = c.faces.n[mirror].dot(st.x) - c.faces.c[mirror];
                if (h > 0) {
                    const Vec3 y = st.x;
                    st.x = st.x - 2.0*h*c.faces.n[mirror];
                    if (c.dirichletFace[mirror]) {   // odd mirror: the value is Phi(y) + Phi(My) - u
                        if (!ctl.odd) throw std::logic_error("runWalk: the odd mirror ball needs WalkControl::odd");
                        ctl.odd->face = mirror;
                        ctl.odd->y = y;
                        ctl.odd->my = st.x;
                    }
                }
            }
            st.onFace = -1;
        } else {
            starStep(c.faces, d, dir, st);
        }
        st.steps++;
        if (st.onFace < 0 && c.faces.outside(st.x, c.escTol)) {   // round-off escape: not counted
            finishWalk(c, st, WalkStatus::Escaped, nullptr, out);
            return WalkEnd::Escaped;
        }
    }
}

// The lean step function: runWalk for solveMulti's use (no WalkControl other than the odd mirror), without the walk's
// box (only WalkCache reads it), and with allAbsorb (every face Dirichlet) without starStep's face tests: there
// R = d <= dF <= the distance dist_f to every face plane, and t = dist_f/(n_f.dir) >= dist_f >= R as long as n_f.dir <= 1
// in floating point, so starStep never hits a face and reduces to x + R dir with onFace = -1. Every other expression
// is runWalk's, so the walk is the same bit for bit. (n_f.dir can exceed 1 by an ulp only for a direction within
// ~2e-8 rad of n_f; when that face is also the nearest boundary, runWalk projects the walk onto the face where this
// one stops an ulp from it, and both end on that face next: probability ~1e-16 per step.) Used by solveMultiLean and
// the profiler (--ring-perf).
WalkEnd runWalkLean(const WalkCtx& c, bool gradPoint, OddMirror* odd, bool allAbsorb, WalkState& st, WalkEntry* out)
{
    for (;;) {
        if (st.steps >= c.maxLen) {
            finishWalk(c, st, WalkStatus::Capped, nullptr, out);
            return WalkEnd::Capped;
        }
        const StripCut::Closest q = c.cut->closest(st.x, c.nStrips);
        int fF = -1;
        const double dF = c.dirichletFaceDist(st.x, fF);
        const double dCut = std::sqrt(q.d2);
        const bool faceCloser = dF < dCut;
        const double d = faceCloser ? dF : dCut;
        if (d < c.shell) {
            finishWalk(c, st, WalkStatus::Landed, faceCloser ? nullptr : &q, out, fF);
            return WalkEnd::Landed;
        }
        const Vec3 dir = sampleSphere(st.rng);
        if (gradPoint && st.steps == 0) {
            double R0 = kGradFirstBallShrink*std::min(d, c.faces.distInside(st.x));
            const int mirror = c.mirrorGradBall ? mirrorGradBall(c, st.x, dCut, R0) : -1;
            st.w = (3.0/R0)*dir;
            st.x = st.x + R0*dir;
            if (mirror >= 0) {
                const double h = c.faces.n[mirror].dot(st.x) - c.faces.c[mirror];
                if (h > 0) {
                    const Vec3 y = st.x;
                    st.x = st.x - 2.0*h*c.faces.n[mirror];
                    if (c.dirichletFace[mirror]) {
                        odd->face = mirror;
                        odd->y = y;
                        odd->my = st.x;
                    }
                }
            }
            st.onFace = -1;
        } else if (allAbsorb) {
            st.x = st.x + d*dir;
            st.onFace = -1;
        } else {
            starStep(c.faces, d, dir, st);
        }
        st.steps++;
        if (st.onFace < 0 && c.faces.outside(st.x, c.escTol)) {
            finishWalk(c, st, WalkStatus::Escaped, nullptr, out);
            return WalkEnd::Escaped;
        }
    }
}

// 2b: a walk's key, the min over its balls of s*_k. The scan stops once the key is at or below the current
// front sCur: every later strip ends at some s_new > sCur, so such a walk is popped at the next advance
// whatever its exact key. nKeys counts the balls whose s* was computed.
double walkKey(const StripCut& cut, const std::vector<Vec2>& Q, double sCur, const WalkPath& p, long long& nKeys)
{
    double key = std::numeric_limits<double>::infinity();
    for (const StoredStep& st : p.s) {
        nKeys++;
        key = std::min(key, cut.entryFront(Q, st.c, st.key2, key));
        if (key <= sCur) break;
    }
    return key;
}

} // namespace

// =====================================================================================
// WalkCache
// =====================================================================================
void WalkCache::configure(const Tet& tet_, const QueryPoints& pts_, int nWalks_, GType g_, const WalkParams& wp_,
                          StepStore store)
{
    // the cache keeps one value per walk (side g, re-evaluated by valueUpdates), runs no walks at t = 0, and its
    // update test knows only the plain balls: none of the coupling gate's flags would be honoured
    if (wp_.usesGateFlags()) throw std::invalid_argument("WalkCache: the coupling gate's WalkParams flags need solveMulti");
    tet = tet_;
    pts = pts_;
    nWalks = std::max(1, nWalks_);
    g = g_;
    wp = wp_;
    store_ = store;
    const size_t P = pts.size(), NW = P*size_t(nWalks);
    gradSlot.assign(P, -1);
    int nGrad = 0;
    for (size_t i = 0; i < P; i++)
        if (pts.grad[i]) gradSlot[i] = nGrad++;
    walks.assign(NW, WalkEntry{});
    W.assign(size_t(nGrad)*nWalks, Vec3::Zero());
    sums.assign(P, PointSums{});
    action.assign(NW, 0);
    rewalkedFrac.assign(P, 0.0f);
    oldU_.assign(NW, 0.0);
    oldCounted_.assign(NW, 0);
    oldW_.assign(W.size(), Vec3::Zero());
    paths.clear();
    paths.shrink_to_fit();
    frontKey.clear();
    queue.clear();
    queueTmp_.clear();
    if (store == StepStore::Paths || store == StepStore::PathsAndKeys) paths.resize(NW);
    if (store == StepStore::PathsAndKeys || store == StepStore::Keys) {
        frontKey.assign(NW, -std::numeric_limits<double>::infinity());
        queue.resize(NW);
        for (size_t wi = 0; wi < NW; wi++) queue[wi] = {frontKey[wi], uint32_t(wi)};
        queueTmp_.resize(NW);
    }
    nStrips_ = 0;
}

const Vec3* WalkCache::gradWeight(size_t walk) const
{
    const int slot = gradSlot[walk/size_t(nWalks)];
    return slot < 0 ? nullptr : &W[size_t(slot)*nWalks + walk%size_t(nWalks)];
}

std::vector<EstimateD> WalkCache::estimates() const
{
    std::vector<EstimateD> e(pts.size());
    for (size_t i = 0; i < pts.size(); i++) e[i] = estimate(i);
    return e;
}

SolveStats WalkCache::solveFull(const StripCut& cut)
{
    auto t0 = std::chrono::steady_clock::now();
    SolveStats st;
    const size_t P = pts.size(), N = size_t(nWalks), NW = P*N;
    st.nPoints = int(P);
    nStrips_ = cut.nStrips();
    std::fill(action.begin(), action.end(), uint8_t(0));
    std::fill(rewalkedFrac.begin(), rewalkedFrac.end(), 0.0f);
    std::fill(sums.begin(), sums.end(), PointSums{});
    if (nStrips_ == 0) {   // t = 0: no Dirichlet boundary, no walks, H := 0 (as the first gate)
        std::fill(walks.begin(), walks.end(), WalkEntry{});
        std::fill(W.begin(), W.end(), Vec3::Zero());
        for (WalkPath& p : paths) p = WalkPath{};
        if (store_ == StepStore::PathsAndKeys || store_ == StepStore::Keys) buildKeys(cut);   // no steps: every key is +inf
        st.noCut = true;
        st.seconds = secondsSince(t0);
        return st;
    }
    // every walk is walked: "re-walked this step" is 1 everywhere after a full solve
    std::fill(action.begin(), action.end(), uint8_t(3));
    std::fill(rewalkedFrac.begin(), rewalkedFrac.end(), 1.0f);
    const WalkCtx c(tet, cut, nStrips_, g, wp);
    const bool record = store_ == StepStore::Paths || store_ == StepStore::PathsAndKeys;
    const bool keyWalks = store_ == StepStore::Keys;   // 3q: each walk accumulates its key as it walks
    const std::vector<Vec2> Q = keyWalks ? cut.futurePolygon(nStrips_) : std::vector<Vec2>{};
    tbb::parallel_for(tbb::blocked_range<size_t>(0, NW), [&](const tbb::blocked_range<size_t>& r) {
        WalkPath scratch;   // stored paths: recorded here, then copied at their exact size
        KeyScan ks;
        WalkControl ctl;
        if (record) ctl.record = &scratch;
        if (keyWalks) {
            ks.Q = &Q;
            ks.sCur = cut.front();
            ctl.key = &ks;
        }
        for (size_t wi = r.begin(); wi < r.end(); wi++) {
            const size_t i = wi/N;
            WalkState s = makeInitialState(pts, nWalks, wi, 0);
            scratch.s.clear();
            scratch.f.clear();
            ks.key = std::numeric_limits<double>::infinity();
            runWalk(c, pts.grad[i] != 0, ctl, s, &walks[wi]);
            walks[wi].rngOffset = 0;
            if (gradSlot[i] >= 0) W[size_t(gradSlot[i])*N + wi%N] = s.w;
            if (record) {
                paths[wi].s = std::vector<StoredStep>(scratch.s.begin(), scratch.s.end());
                paths[wi].f = std::vector<uint8_t>(scratch.f.begin(), scratch.f.end());
            }
            if (keyWalks) {
                frontKey[wi] = ks.key;
                queue[wi] = {ks.key, uint32_t(wi)};
            }
        }
    });
    if (store_ == StepStore::PathsAndKeys) buildKeys(cut);
    if (keyWalks) tbb::parallel_sort(queue.begin(), queue.end());
    // per-point sums, one task per point, in walk order
    std::vector<long long> nCounted(P), nCapped(P), nSteps(P);
    tbb::parallel_for(size_t(0), P, [&](size_t i) {
        PointSums& S = sums[i];
        for (size_t j = 0; j < N; j++) {
            const WalkEntry& e = walks[i*N + j];
            nSteps[i] += e.steps;
            nCapped[i] += e.status == WalkStatus::Capped;
            if (!countedWalk(e)) continue;
            nCounted[i]++;
            S.add(e.value, gradSlot[i] >= 0 ? W[size_t(gradSlot[i])*N + j] : Vec3::Zero(), 1.0);
        }
    });
    st.walksRequested = (long long)NW;
    for (size_t i = 0; i < P; i++) {
        st.walksCounted += nCounted[i];
        st.walksCapped += nCapped[i];
        st.totalSteps += nSteps[i];
    }
    st.walksEscaped = st.walksRequested - st.walksCounted;
    st.seconds = secondsSince(t0);
    return st;
}

ReuseStats WalkCache::advance(const StripCut& cut, UpdateMode mode, std::vector<ResumeRecord>* resumes)
{
    auto t0 = std::chrono::steady_clock::now();
    ReuseStats rs;
    const size_t P = pts.size(), N = size_t(nWalks), NW = P*N;
    rs.nWalks = (long long)NW;
    if (resumes) resumes->clear();
    const int m = cut.nStrips() - 1;   // the new strip
    const bool storedMode = mode == UpdateMode::StoredBox || mode == UpdateMode::StoredQueue;
    const bool replayQueue = mode == UpdateMode::ReplayQueue;
    // a cache's stored paths (and keys) must follow every update, so a cache is used with its own modes only
    if (storedMode && (store_ == StepStore::None || store_ == StepStore::Keys ||
                       (mode == UpdateMode::StoredQueue && store_ != StepStore::PathsAndKeys)))
        throw std::logic_error("WalkCache::advance: stored-step mode on a cache configured without the stored data");
    if (replayQueue && store_ != StepStore::Keys)
        throw std::logic_error("WalkCache::advance: replay-queue mode needs a cache configured with keys only");
    if (!storedMode && !replayQueue && mode != UpdateMode::FullResolve && store_ != StepStore::None)
        throw std::logic_error("WalkCache::advance: replay mode on a cache with stored steps or keys (they would go stale)");
    if (mode == UpdateMode::FullResolve || nStrips_ == 0 || m != nStrips_) {
        rs.solve = solveFull(cut);
        rs.full = true;
        rs.nTested = rs.nRewalked = rs.nWalks;
        rs.escaped = rs.solve.walksEscaped;
        rs.capped = rs.solve.walksCapped;
        rs.stepsNew = rs.solve.totalSteps;
        rs.stripArea = m >= 0 ? cut.strip(m).area : 0.0;
        rs.tTotal = secondsSince(t0);
        return rs;
    }
    if (storedMode) return advanceStored(cut, mode);
    rs.stripArea = cut.strip(m).area;
    const WalkCtx cOld(tet, cut, m, g, wp), cNew(tet, cut, m + 1, g, wp);

    // 1. bounding box of every walk against the new strip; non-intersecting walks are reused as is.
    //    3q: instead, the walks popped from the queue as in 2b; the walks after them are not touched at all
    size_t nPop = 0;
    if (replayQueue) {
        nPop = popQueue(cut.front(m + 1));
    } else {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, NW), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t wi = r.begin(); wi < r.end(); wi++)
                action[wi] = cut.boxMayHitStrip(walks[wi].lo, walks[wi].hi, m) ? kTested : 0;
        });
    }
    auto t1 = std::chrono::steady_clock::now();
    rs.tBox = std::chrono::duration<double>(t1 - t0).count();

    // 2. bbox hits: replay from the seed under the OLD geometry, testing each step's ball against the
    //    strip. At the first step k whose ball the strip enters, save the old contribution and continue
    //    the walk from that step's state (center, face, RNG, step counter k, box of balls 0..k-1) under
    //    the NEW geometry. Steps 0..k-1 are bit-identical under both geometries (the strip changes
    //    neither their distances nor their closest points), so the result is the walk a from-scratch
    //    solve with this seed computes. Biased mode restarts from the query point instead.
    //    3q: the replayed steps and then the resumed ones are folded into the walk's new key, against the
    //    future polygon after the new strip, which is what 2b's walkKey computes from the stored path.
    std::vector<size_t> hits;
    if (replayQueue) {
        hits.resize(nPop);
        for (size_t h = 0; h < nPop; h++) hits[h] = queue[h].second;
    } else {
        for (size_t wi = 0; wi < NW; wi++)
            if (action[wi]) hits.push_back(wi);
    }
    const std::vector<Vec2> Qn = replayQueue ? cut.futurePolygon(m + 1) : std::vector<Vec2>{};
    std::vector<ResumeRecord> rec(resumes ? hits.size() : 0);
    std::vector<uint8_t> recUsed(rec.size(), 0);
    std::atomic<long long> stepsReplayed{0}, replayedFalse{0}, stepsNew{0}, nEsc{0}, nCap{0}, nK0{0}, sumK{0}, nKeys{0};
    tbb::parallel_for(tbb::blocked_range<size_t>(0, hits.size()), [&](const tbb::blocked_range<size_t>& r) {
        long long locReplayed = 0, locFalse = 0, locNew = 0, locEsc = 0, locCap = 0, locK0 = 0, locSumK = 0, locKeys = 0;
        for (size_t h = r.begin(); h < r.end(); h++) {
            const size_t wi = hits[h], i = wi/N, j = wi%N;
            WalkEntry& e = walks[wi];
            const bool grad = pts.grad[i] != 0;
            WalkState s = makeInitialState(pts, nWalks, wi, e.rngOffset);
            KeyScan ks;
            ks.Q = &Qn;
            ks.sCur = cut.front(m + 1);
            WalkControl ctl;
            ctl.stopStrip = m;
            if (replayQueue) ctl.key = &ks;
            if (runWalk(cOld, grad, ctl, s, nullptr) != WalkEnd::Stopped) {   // no ball meets the strip
                locReplayed += s.steps;
                locFalse += s.steps;
                if (replayQueue) {   // a false alarm keeps its walk; its key is recomputed (the front has moved)
                    frontKey[wi] = ks.key;
                    locKeys += ks.n;
                }
                continue;
            }
            const int k = s.steps;
            locReplayed += k;
            locK0 += k == 0;
            locSumK += k;
            oldU_[wi] = e.value;
            oldCounted_[wi] = countedWalk(e);
            Vec3* w = gradSlot[i] >= 0 ? &W[size_t(gradSlot[i])*N + j] : nullptr;
            if (w) oldW_[size_t(gradSlot[i])*N + j] = *w;
            if (resumes) {
                rec[h].walk = wi;
                rec[h].k = k;
                rec[h].oldSteps = e.steps;
                rec[h].rngOffset = e.rngOffset;
                rec[h].state = s;
                recUsed[h] = 1;
            }
            uint32_t offset = e.rngOffset;
            if (mode == UpdateMode::BiasedRestart) {
                // Fig. 5 of the paper, deliberately biased: a new walk from the query point. Its random
                // numbers continue from the RNG state at step k (statistically a fresh walk); offset
                // records the draws before its step 0 so that it stays replayable from the seed.
                offset += uint32_t(2*k);
                pcg32 rng = s.rng;
                s = makeInitialState(pts, nWalks, wi, 0);
                s.rng = rng;
            }
            WalkEntry ne;
            const int stepsBefore = s.steps;
            WalkControl ctlNew;
            if (replayQueue) ctlNew.key = &ks;   // the same scan continues over the resumed steps
            WalkEnd end = runWalk(cNew, grad, ctlNew, s, &ne);
            ne.rngOffset = offset;
            locNew += s.steps - stepsBefore;
            locEsc += end == WalkEnd::Escaped;
            locCap += end == WalkEnd::Capped;
            e = ne;
            if (w) *w = s.w;
            action[wi] |= kRewalked;
            if (replayQueue) {
                frontKey[wi] = ks.key;
                locKeys += ks.n;
            }
        }
        stepsReplayed += locReplayed;
        replayedFalse += locFalse;
        nK0 += locK0;
        sumK += locSumK;
        stepsNew += locNew;
        nEsc += locEsc;
        nCap += locCap;
        nKeys += locKeys;
    });
    if (resumes)
        for (size_t h = 0; h < rec.size(); h++)
            if (recUsed[h]) resumes->push_back(rec[h]);
    rs.stepsReplayed = stepsReplayed.load();
    rs.stepsReplayedFalse = replayedFalse.load();
    rs.keysComputed = nKeys.load();
    rs.nRewalkedAtStep0 = nK0.load();
    rs.sumK = sumK.load();
    rs.stepsNew = stepsNew.load();
    rs.escaped = nEsc.load();
    rs.capped = nCap.load();
    auto t2 = std::chrono::steady_clock::now();
    rs.tRewalk = std::chrono::duration<double>(t2 - t1).count();

    // 3. value updates (g type (b) only)
    valueUpdates(cNew);
    auto t3 = std::chrono::steady_clock::now();
    rs.tValue = std::chrono::duration<double>(t3 - t2).count();

    // 4. sums
    updateSums(rs);
    nStrips_ = m + 1;
    auto t4 = std::chrono::steady_clock::now();
    rs.tSums = std::chrono::duration<double>(t4 - t3).count();

    // 5. (3q) the popped walks go back into the queue with their new keys
    if (replayQueue) pushQueue(nPop);
    auto t5 = std::chrono::steady_clock::now();
    rs.tQueue = std::chrono::duration<double>(t5 - t4).count();
    rs.tTotal = std::chrono::duration<double>(t5 - t0).count();
    rs.scratchBytes = hits.capacity()*sizeof(size_t) + rec.capacity()*sizeof(ResumeRecord) + recUsed.capacity() +
                      Qn.capacity()*sizeof(Vec2) + nPop*sizeof(queue[0]);
    return rs;
}

// 2b and 3q: the walks at the front of the queue whose key is <= sNew are this update's candidates (action =
// kTested); every other walk gets action 0 and is not touched
size_t WalkCache::popQueue(double sNew)
{
    const size_t nPop = size_t(std::upper_bound(queue.begin(), queue.end(), sNew,
                                                [](double s, const std::pair<double, uint32_t>& e) { return s < e.first; }) -
                               queue.begin());
    std::fill(action.begin(), action.end(), uint8_t(0));
    for (size_t h = 0; h < nPop; h++) action[queue[h].second] = kTested;
    return nPop;
}

// 2b and 3q: re-insert the popped walks queue[0, nPop) with their new frontKey: sort them, merge them with the
// untouched rest
void WalkCache::pushQueue(size_t nPop)
{
    std::vector<std::pair<double, uint32_t>> ins(nPop);
    for (size_t h = 0; h < nPop; h++) ins[h] = {frontKey[queue[h].second], queue[h].second};
    tbb::parallel_sort(ins.begin(), ins.end());
    std::merge(ins.begin(), ins.end(), queue.begin() + nPop, queue.end(), queueTmp_.begin());
    queue.swap(queueTmp_);
}

// 3. value updates (g type (b) only: g = 1 never changes): a reused walk keeps its landing point, and its
//    value is re-evaluated against the new front. No re-walk, no geometry needed. Shared by every
//    incremental mode (options 2a, 2b and 3), so their values are bit for bit the same.
void WalkCache::valueUpdates(const WalkCtx& cNew)
{
    if (g != GType::SqrtFront) return;
    const size_t NW = walks.size();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, NW), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t wi = r.begin(); wi < r.end(); wi++) {
            WalkEntry& e = walks[wi];
            if ((action[wi] & kRewalked) || e.status != WalkStatus::Landed) continue;
            const double v = double(e.side)*cNew.front.g(e.land, g);
            if (sameBits(&v, &e.value, sizeof v)) continue;
            oldU_[wi] = e.value;
            e.value = v;
            action[wi] |= kValue;
        }
    });
}

// 4. sums: subtract old, add new, one task per point in walk order (thread-count independent). Shared by
//    every incremental mode, so the same re-walked and value-updated walks give bitwise the same sums.
void WalkCache::updateSums(ReuseStats& rs)
{
    const size_t P = pts.size(), N = size_t(nWalks);
    std::vector<long long> nTested(P), nRe(P), nVal(P);
    tbb::parallel_for(size_t(0), P, [&](size_t i) {
        PointSums& S = sums[i];
        const int slot = gradSlot[i];
        for (size_t j = 0; j < N; j++) {
            const size_t wi = i*N + j;
            const uint8_t a = action[wi];
            if (!a) continue;
            nTested[i] += (a & kTested) != 0;
            const WalkEntry& e = walks[wi];
            const Vec3 w = slot >= 0 ? W[size_t(slot)*N + j] : Vec3::Zero();
            if (a & kRewalked) {
                nRe[i]++;
                if (oldCounted_[wi]) S.add(oldU_[wi], slot >= 0 ? oldW_[size_t(slot)*N + j] : Vec3::Zero(), -1.0);
                if (countedWalk(e)) S.add(e.value, w, 1.0);
            } else if (a & kValue) {
                nVal[i]++;
                S.add(oldU_[wi], w, -1.0);
                S.add(e.value, w, 1.0);
            }
        }
        rewalkedFrac[i] = float(double(nRe[i])/double(N));
    });
    for (size_t i = 0; i < P; i++) {
        rs.nTested += nTested[i];
        rs.nRewalked += nRe[i];
        rs.nValueUpdated += nVal[i];
    }
}

// Options 2a and 2b: every step of every walk is stored (centre, key, face, face-closer flag), so the
// update tests the stored balls directly and resumes from the stored state at the first ball the strip
// enters: no replay. The passes:
//   1. filter: 2a, the same bbox test as option 3; 2b, pop the walks whose front key is <= the new front;
//   2. exact test of each candidate's stored balls in step order, with exactly runWalk's stop criterion
//      stripDist2(c_k) < key2_k (same function, same key as the walk computed it, so the dirichletFaces
//      branch and its 1e-12 margin are included); the first ball it holds for is the invalidation step k;
//   3. resume states rebuilt from the stored step k; 4. the new steps, by runWalk, appended to the path;
//   5. (2b) new keys of every popped walk; 6./7. value updates and sums, shared with option 3;
//   8. (2b) the popped walks go back into the sorted queue.
// Why a popped walk still needs the exact test: its key is a conservative lower bound (margins on the ball
// radius and on s*), and ties are real (a ball tangent to a front edge the strip shares, whose distance
// does not drop): key <= s_new only says the strip may enter one of its balls. Only the exact test decides
// that the step changes, as it does in option 3's replay, so popped-but-unaffected walks stay bit for bit.
ReuseStats WalkCache::advanceStored(const StripCut& cut, UpdateMode mode)
{
    auto t0 = std::chrono::steady_clock::now();
    ReuseStats rs;
    rs.stored = true;
    const size_t NW = walks.size(), N = size_t(nWalks);
    rs.nWalks = (long long)NW;
    const int m = cut.nStrips() - 1;   // the new strip
    rs.stripArea = cut.strip(m).area;
    const WalkCtx cNew(tet, cut, m + 1, g, wp);
    const bool queued = mode == UpdateMode::StoredQueue;

    // 1. candidates
    std::vector<uint32_t> cand;
    size_t nPop = 0;
    if (!queued) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, NW), [&](const tbb::blocked_range<size_t>& r) {
            for (size_t wi = r.begin(); wi < r.end(); wi++)
                action[wi] = cut.boxMayHitStrip(walks[wi].lo, walks[wi].hi, m) ? kTested : 0;
        });
        for (size_t wi = 0; wi < NW; wi++)
            if (action[wi]) cand.push_back(uint32_t(wi));
    } else {
        // the queue's prefix with key <= s_new; the walks after it are not touched at all
        nPop = popQueue(cut.front(m + 1));
        cand.resize(nPop);
        for (size_t h = 0; h < nPop; h++) cand[h] = queue[h].second;
    }
    auto t1 = std::chrono::steady_clock::now();
    rs.tBox = std::chrono::duration<double>(t1 - t0).count();

    // 2. exact test of the stored balls
    std::vector<int> kHit(cand.size(), -1);
    std::atomic<long long> balls{0};
    tbb::parallel_for(tbb::blocked_range<size_t>(0, cand.size()), [&](const tbb::blocked_range<size_t>& r) {
        long long b = 0;
        for (size_t h = r.begin(); h < r.end(); h++) {
            const WalkPath& p = paths[cand[h]];
            for (size_t k = 0; k < p.s.size(); k++) {
                b++;
                if (cut.stripDist2(p.s[k].c, m) < p.s[k].key2) {
                    kHit[h] = int(k);
                    break;
                }
            }
        }
        balls += b;
    });
    rs.ballsTested = balls.load();
    std::vector<uint32_t> hits;
    std::vector<int> hitK;
    for (size_t h = 0; h < cand.size(); h++)
        if (kHit[h] >= 0) {
            hits.push_back(cand[h]);
            hitK.push_back(kHit[h]);
            rs.nRewalkedAtStep0 += kHit[h] == 0;
            rs.sumK += kHit[h];
        }
    auto t2 = std::chrono::steady_clock::now();
    rs.tTest = std::chrono::duration<double>(t2 - t1).count();

    // 3. resume states. The state runWalk has before step k is (x_k, onFace_k, k, RNG, w, box of balls
    //    0..k-1). x_k and onFace_k are stored exactly as the walk had them when it evaluated step k (steps
    //    0..k-1 are unchanged by the strip, so this is also the state the walk has under the new cut). The
    //    RNG depends only on the draw count: 2 per move, none at the end, so before step k it is the seed's
    //    generator advanced by rngOffset + 2k draws, and pcg32::advance jumps the LCG there exactly. w is set
    //    at step 0 only, so for k >= 1 it is the walk's cached W (k = 0: the initial state, w = 0). The box
    //    is rebuilt from the stored centres and radii with the same addBall: the radius is sqrt(key2) =
    //    sqrt(q.d2) (bitwise runWalk's dCut) or, with the face-closer flag, distInside(c_j) (bitwise its dF).
    //    The old contribution is saved for the sums, as in option 3.
    std::vector<WalkState> resume(hits.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, hits.size()), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t h = r.begin(); h < r.end(); h++) {
            const size_t wi = hits[h], i = wi/N, j = wi%N;
            const int k = hitK[h];
            const WalkEntry& e = walks[wi];
            const WalkPath& p = paths[wi];
            const int slot = gradSlot[i];
            oldU_[wi] = e.value;
            oldCounted_[wi] = countedWalk(e);
            if (slot >= 0) oldW_[size_t(slot)*N + j] = W[size_t(slot)*N + j];
            WalkState& s = resume[h];
            s = makeInitialState(pts, nWalks, wi, e.rngOffset);
            if (k > 0) {
                s.x = p.s[k].c;
                s.onFace = int(p.f[k] & kFaceMask) - 1;
                s.steps = k;
                s.rng.advance(int64_t(2)*k);
                if (slot >= 0) s.w = W[size_t(slot)*N + j];
                for (int q = 0; q < k; q++) {
                    const Vec3& c = p.s[q].c;
                    int fq;
                    addBall(s, c, (p.f[q] & kFaceCloser) ? cNew.dirichletFaceDist(c, fq) : std::sqrt(p.s[q].key2));
                }
            }
        }
    });
    auto t3 = std::chrono::steady_clock::now();
    rs.tResume = std::chrono::duration<double>(t3 - t2).count();

    // 4. new steps: runWalk records k, k+1, ... under the new cut, and the path keeps steps 0..k-1 and gets
    //    them appended: in place if they fit, else in an allocation of exactly the new length (no geometric
    //    growth), so a path's capacity is the longest it has been
    std::atomic<long long> stepsNew{0}, nEsc{0}, nCap{0};
    auto replaceTail = [](auto& v, size_t k, const auto& tail) {
        const size_t n = k + tail.size();
        if (n > v.capacity()) {
            std::remove_reference_t<decltype(v)> w;
            w.reserve(n);
            w.insert(w.end(), v.begin(), v.begin() + k);
            v.swap(w);
        } else {
            v.resize(k);
        }
        v.insert(v.end(), tail.begin(), tail.end());
    };
    tbb::parallel_for(tbb::blocked_range<size_t>(0, hits.size()), [&](const tbb::blocked_range<size_t>& r) {
        long long locNew = 0, locEsc = 0, locCap = 0;
        WalkPath scratch;
        WalkControl ctl;
        ctl.record = &scratch;
        for (size_t h = r.begin(); h < r.end(); h++) {
            const size_t wi = hits[h], i = wi/N;
            const int k = hitK[h];
            scratch.s.clear();
            scratch.f.clear();
            WalkState& s = resume[h];
            WalkEntry ne;
            const WalkEnd end = runWalk(cNew, pts.grad[i] != 0, ctl, s, &ne);
            WalkPath& p = paths[wi];
            replaceTail(p.s, size_t(k), scratch.s);
            replaceTail(p.f, size_t(k), scratch.f);
            ne.rngOffset = walks[wi].rngOffset;
            locNew += s.steps - k;
            locEsc += end == WalkEnd::Escaped;
            locCap += end == WalkEnd::Capped;
            walks[wi] = ne;
            if (gradSlot[i] >= 0) W[size_t(gradSlot[i])*N + wi%N] = s.w;
            action[wi] |= kRewalked;
        }
        stepsNew += locNew;
        nEsc += locEsc;
        nCap += locCap;
    });
    rs.stepsNew = stepsNew.load();
    rs.escaped = nEsc.load();
    rs.capped = nCap.load();
    auto t4 = std::chrono::steady_clock::now();
    rs.tRewalk = std::chrono::duration<double>(t4 - t3).count();

    // 5. (2b) new keys of the popped walks, re-walked or not, against the new future polygon. A false alarm
    //    keeps its path; its key is recomputed because the front has moved: the part of a ball's disk the
    //    front has passed no longer counts, so a ball that only touched the old front (an exact tie) drops out.
    if (queued) {
        const std::vector<Vec2> Q = cut.futurePolygon(m + 1);
        const double sCur = cut.front(m + 1);
        std::atomic<long long> nKeys{0};
        tbb::parallel_for(tbb::blocked_range<size_t>(0, cand.size()), [&](const tbb::blocked_range<size_t>& r) {
            long long nk = 0;
            for (size_t h = r.begin(); h < r.end(); h++)
                frontKey[cand[h]] = walkKey(cut, Q, sCur, paths[cand[h]], nk);
            nKeys += nk;
        });
        rs.keysComputed = nKeys.load();
    }
    auto t5 = std::chrono::steady_clock::now();
    rs.tKeys = std::chrono::duration<double>(t5 - t4).count();

    // 6. value updates, 7. sums (shared with option 3)
    valueUpdates(cNew);
    auto t6 = std::chrono::steady_clock::now();
    rs.tValue = std::chrono::duration<double>(t6 - t5).count();
    updateSums(rs);
    nStrips_ = m + 1;
    auto t7 = std::chrono::steady_clock::now();
    rs.tSums = std::chrono::duration<double>(t7 - t6).count();

    // 8. (2b) re-insert the popped walks: sort their new keys, merge them with the untouched rest
    if (queued) pushQueue(nPop);
    auto t8 = std::chrono::steady_clock::now();
    rs.tQueue = std::chrono::duration<double>(t8 - t7).count();
    rs.tTotal = std::chrono::duration<double>(t8 - t0).count();
    rs.scratchBytes = cand.capacity()*sizeof(uint32_t) + kHit.capacity()*sizeof(int) + hits.capacity()*sizeof(uint32_t) +
                      hitK.capacity()*sizeof(int) + resume.capacity()*sizeof(WalkState) + nPop*sizeof(queue[0]);
    return rs;
}

// 2b: every walk's key against the current future polygon, and the queue sorted by (key, walk). 3q has no
// paths and calls it only at t = 0, where no walk has a step: every key is +inf, as walkKey gives 2b then.
void WalkCache::buildKeys(const StripCut& cut)
{
    const size_t NW = walks.size();
    const std::vector<Vec2> Q = cut.futurePolygon(cut.nStrips());
    const double sCur = cut.front();
    const bool noPaths = paths.empty();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, NW), [&](const tbb::blocked_range<size_t>& r) {
        long long nk = 0;
        for (size_t wi = r.begin(); wi < r.end(); wi++) {
            frontKey[wi] = noPaths ? std::numeric_limits<double>::infinity() : walkKey(cut, Q, sCur, paths[wi], nk);
            queue[wi] = {frontKey[wi], uint32_t(wi)};
        }
    });
    tbb::parallel_sort(queue.begin(), queue.end());
}

CacheMemory WalkCache::memory() const
{
    CacheMemory mem;
    mem.entries = walks.capacity()*sizeof(WalkEntry) + action.capacity() + oldU_.capacity()*sizeof(double) +
                  oldCounted_.capacity();
    mem.gradWeights = (W.capacity() + oldW_.capacity())*sizeof(Vec3);
    mem.points = sums.capacity()*sizeof(PointSums) + rewalkedFrac.capacity()*sizeof(float) + gradSlot.capacity()*sizeof(int) +
                 pts.x.capacity()*sizeof(Vec3) + pts.grad.capacity() + pts.seed.capacity()*sizeof(uint64_t);
    mem.pathHeaders = paths.capacity()*sizeof(WalkPath);
    for (const WalkPath& p : paths) {
        mem.nSteps += (long long)p.s.size();
        mem.nStepCapacity += (long long)p.s.capacity();
        mem.stepsAllocated += p.s.capacity()*sizeof(StoredStep) + p.f.capacity();
    }
    mem.stepsPayload = size_t(mem.nSteps)*kBytesPerStoredStep;
    mem.keys = frontKey.capacity()*sizeof(double) + (queue.capacity() + queueTmp_.capacity())*sizeof(queue[0]);
    return mem;
}

// =====================================================================================
// WalkRunner (validation)
// =====================================================================================
WalkRunner::WalkRunner(const WalkCache& cache, const StripCut& cut, int nStrips)
    : cache_(cache), nStrips_(nStrips),
      ctx_(std::make_unique<WalkCtx>(cache.tet, cut, nStrips, cache.g, cache.wp))
{
}

WalkRunner::~WalkRunner() = default;

WalkState WalkRunner::initialState(size_t walk, uint32_t rngOffset) const
{
    return makeInitialState(cache_.pts, cache_.nWalks, walk, rngOffset);
}

void WalkRunner::replay(size_t walk, WalkEntry& out, Vec3& w) const
{
    const WalkEntry& e = cache_.walks[walk];
    WalkState s = initialState(walk, e.rngOffset);
    runWalk(*ctx_, cache_.pts.grad[walk/size_t(cache_.nWalks)] != 0, WalkControl{}, s, &out);
    out.rngOffset = e.rngOffset;
    w = s.w;
}

WalkState WalkRunner::runTo(size_t walk, uint32_t rngOffset, int k,
                            const std::function<void(const WalkState&, const StripCut::Closest&)>* observe) const
{
    WalkState s = initialState(walk, rngOffset);
    WalkControl ctl;
    ctl.pauseAt = k;
    ctl.observe = observe;
    runWalk(*ctx_, cache_.pts.grad[walk/size_t(cache_.nWalks)] != 0, ctl, s, nullptr);
    return s;
}

// =====================================================================================
// Sums-only solve
// =====================================================================================
std::vector<EstimateD> solveSumsOnly(const Tet& tet, const StripCut& cut, const QueryPoints& pts, int nWalks,
                                     GType g, const WalkParams& wp, SolveStats* stats)
{
    if (wp.usesGateFlags()) throw std::invalid_argument("solveSumsOnly: the coupling gate's WalkParams flags need solveMulti");
    auto t0 = std::chrono::steady_clock::now();
    const size_t P = pts.size(), N = size_t(std::max(1, nWalks));
    std::vector<EstimateD> est(P);
    SolveStats st;
    st.nPoints = int(P);
    if (cut.nStrips() == 0) {
        st.noCut = true;
        st.seconds = secondsSince(t0);
        if (stats) *stats = st;
        return est;
    }
    const WalkCtx c(tet, cut, cut.nStrips(), g, wp);
    // blocks of 64 walks of one point; block partial sums are combined in block order
    const size_t B = 64, nb = (N + B - 1)/B;
    std::vector<PointSums> part(P*nb);
    std::vector<long long> counted(P*nb), capped(P*nb), steps(P*nb);
    tbb::parallel_for(size_t(0), P*nb, [&](size_t b) {
        const size_t i = b/nb, j0 = (b%nb)*B, j1 = std::min(N, j0 + B);
        for (size_t j = j0; j < j1; j++) {
            WalkState s = makeInitialState(pts, int(N), i*N + j, 0);
            WalkEntry e;
            runWalk(c, pts.grad[i] != 0, WalkControl{}, s, &e);
            steps[b] += e.steps;
            capped[b] += e.status == WalkStatus::Capped;
            if (!countedWalk(e)) continue;
            counted[b]++;
            part[b].add(e.value, s.w, 1.0);
        }
    });
    for (size_t i = 0; i < P; i++) {
        PointSums S;
        for (size_t k = 0; k < nb; k++) {
            S.add(part[i*nb + k]);
            st.walksCounted += counted[i*nb + k];
            st.walksCapped += capped[i*nb + k];
            st.totalSteps += steps[i*nb + k];
        }
        est[i] = estimateFromSums(S, pts.grad[i] != 0);
    }
    st.walksRequested = (long long)(P*N);
    st.walksEscaped = st.walksRequested - st.walksCounted;
    st.seconds = secondsSince(t0);
    if (stats) *stats = st;
    return est;
}

// =====================================================================================
// Reuse option 1: frozen seeds, recompute everything, sums only
// =====================================================================================
std::vector<PointSums> solveSums(const Tet& tet, const StripCut& cut, const QueryPoints& pts, int nWalks, GType g,
                                 const WalkParams& wp, SolveStats* stats)
{
    if (wp.usesGateFlags()) throw std::invalid_argument("solveSums: the coupling gate's WalkParams flags need solveMulti");
    auto t0 = std::chrono::steady_clock::now();
    const size_t P = pts.size(), N = size_t(std::max(1, nWalks));
    std::vector<PointSums> sums(P);
    SolveStats st;
    st.nPoints = int(P);
    if (cut.nStrips() == 0) {   // t = 0: no walks, H := 0
        st.noCut = true;
        st.seconds = secondsSince(t0);
        if (stats) *stats = st;
        return sums;
    }
    const WalkCtx c(tet, cut, cut.nStrips(), g, wp);
    std::vector<long long> counted(P), capped(P), steps(P);
    // one task per point, its walks in walk order: the additions are exactly WalkCache::solveFull's
    tbb::parallel_for(size_t(0), P, [&](size_t i) {
        PointSums& S = sums[i];
        const bool grad = pts.grad[i] != 0;
        for (size_t j = 0; j < N; j++) {
            WalkState s = makeInitialState(pts, int(N), i*N + j, 0);
            WalkEntry e;
            runWalk(c, grad, WalkControl{}, s, &e);
            steps[i] += e.steps;
            capped[i] += e.status == WalkStatus::Capped;
            if (!countedWalk(e)) continue;
            counted[i]++;
            S.add(e.value, grad ? s.w : Vec3::Zero(), 1.0);
        }
    });
    for (size_t i = 0; i < P; i++) {
        st.walksCounted += counted[i];
        st.walksCapped += capped[i];
        st.totalSteps += steps[i];
    }
    st.walksRequested = (long long)(P*N);
    st.walksEscaped = st.walksRequested - st.walksCounted;
    st.seconds = secondsSince(t0);
    if (stats) *stats = st;
    return sums;
}

// =====================================================================================
// Several boundary data on the same walks (coupling gate)
// =====================================================================================
namespace {

// solveMulti on the walk context c (a tet's, or the ring reference's convex domain); t0: the start of the solve
MultiSums solveMultiOn(const WalkCtx& c, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                       const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats,
                       const HarmonicExtension* ext, std::chrono::steady_clock::time_point t0)
{
    if (wp.mirrorGradBall && wp.oddMirrorDirichlet && !(ext && *ext))
        throw std::invalid_argument("solveMulti: the odd mirror ball needs the data's harmonic extensions");
    MultiSums ms;
    ms.K = K;
    ms.nBatches = std::max(1, nBatches);
    ms.nPoints = pts.size();
    if (nWalks%ms.nBatches) throw std::invalid_argument("solveMulti: nWalks must be a multiple of nBatches");
    ms.walksPerBatch = nWalks/ms.nBatches;
    const size_t P = pts.size(), N = size_t(nWalks), B = size_t(ms.nBatches), NB = size_t(ms.walksPerBatch);
    ms.batch.assign(P*B*ms.stride(), 0.0);
    ms.uu.assign(P*size_t(K)*size_t(K), 0.0);
    SolveStats st;
    st.nPoints = int(P);
    if (cut.nStrips() == 0 && !c.anyDirichletFace) {   // no Dirichlet boundary at all: no walks
        st.noCut = true;
        st.seconds = secondsSince(t0);
        if (stats) *stats = st;
        return ms;
    }
    std::vector<long long> counted(P), capped(P), steps(P);
    // one task per point: its batches in order, its walks in walk order
    tbb::parallel_for(size_t(0), P, [&](size_t i) {
        const bool grad = pts.grad[i] != 0;
        std::vector<double> u(static_cast<size_t>(K)), phiY(static_cast<size_t>(K)), phiMy(static_cast<size_t>(K));
        double* uu = &ms.uu[i*size_t(K)*size_t(K)];
        OddMirror odd;
        WalkControl ctl;
        ctl.odd = &odd;
        for (size_t b = 0; b < B; b++) {
            double* s = &ms.batch[(i*B + b)*ms.stride()];
            double *su = s + 1, *sw = s + 1 + K, *swu = s + 4 + K;
            for (size_t j = b*NB; j < (b + 1)*NB; j++) {
                WalkState ws = makeInitialState(pts, int(N), i*N + j, 0);
                WalkEntry e;
                odd.face = -1;
                runWalk(c, grad, ctl, ws, &e);
                steps[i] += e.steps;
                capped[i] += e.status == WalkStatus::Capped;
                if (!countedWalk(e)) continue;
                counted[i]++;
                if (e.status == WalkStatus::Landed) {
                    WalkTerminal T;
                    T.status = e.status;
                    T.face = e.face;
                    T.side = e.side;
                    T.x = e.face >= 0 ? c.onFacePlane(e.face, e.land) : e.land;
                    values(T, u.data());
                } else {
                    std::fill(u.begin(), u.end(), 0.0);   // capped: value 0, as everywhere else
                }
                if (odd.face >= 0) {   // the first step crossed a Dirichlet face: H~(y) = Phi(y) + Phi(My) - H(My)
                    (*ext)(odd.y, phiY.data());
                    (*ext)(odd.my, phiMy.data());
                    for (int k = 0; k < K; k++) u[k] = phiY[k] + phiMy[k] - u[k];
                }
                const Vec3 w = grad ? ws.w : Vec3::Zero();
                s[0] += 1;
                for (int a = 0; a < 3; a++) sw[a] += w[a];
                for (int k = 0; k < K; k++) {
                    su[k] += u[k];
                    for (int a = 0; a < 3; a++) swu[3*k + a] += w[a]*u[k];
                    for (int l = 0; l < K; l++) uu[k*K + l] += u[k]*u[l];
                }
            }
        }
    });
    for (size_t i = 0; i < P; i++) {
        st.walksCounted += counted[i];
        st.walksCapped += capped[i];
        st.totalSteps += steps[i];
    }
    st.walksRequested = (long long)(P*N);
    st.walksEscaped = st.walksRequested - st.walksCounted;
    st.seconds = secondsSince(t0);
    if (stats) *stats = st;
    return ms;
}

} // namespace

MultiSums solveMulti(const Tet& tet, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                     const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats,
                     const HarmonicExtension* ext)
{
    auto t0 = std::chrono::steady_clock::now();
    const WalkCtx c(tet, cut, cut.nStrips(), GType::One, wp);   // values come from 'values', not from c.g
    return solveMultiOn(c, cut, pts, nWalks, nBatches, wp, K, values, stats, ext, t0);
}

MultiSums solveMulti(const ConvexDomain& dom, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                     const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats)
{
    auto t0 = std::chrono::steady_clock::now();
    const WalkCtx c(dom, cut, cut.nStrips(), GType::One, wp);   // throws if wp flags a Dirichlet face
    return solveMultiOn(c, cut, pts, nWalks, nBatches, wp, K, values, stats, nullptr, t0);
}

// =====================================================================================
// Lean multi-data solves (opt-in: the fast ring path, --ring-bench and the ring view's fast option)
// =====================================================================================
// solveMultiOn's loop with runWalkLean and per-point scratch on the stack, over the points of several jobs at once.
// The walks, their scoring and the order of every sum are solveMultiOn's, so each job's MultiSums equal solveMulti's.
std::vector<MultiSums> solveMultiLeanJobs(const std::vector<MultiJob>& jobs, std::vector<SolveStats>* stats)
{
    constexpr int kMaxK = 16;
    const auto t0 = std::chrono::steady_clock::now();
    const size_t J = jobs.size();
    std::vector<std::unique_ptr<WalkCtx>> ctx(J);
    std::vector<MultiSums> out(J);
    std::vector<char> active(J, 0), allAbsorb(J, 0);
    std::vector<size_t> first(J + 1, 0);   // global index of each active job's first point
    for (size_t j = 0; j < J; j++) {
        const MultiJob& b = jobs[j];
        if (!b.cut || !b.pts || (!b.tet == !b.dom)) throw std::invalid_argument("solveMultiLeanJobs: a job needs a cut, points and one domain");
        if (b.K < 0 || b.K > kMaxK) throw std::invalid_argument("solveMultiLeanJobs: K out of range");
        if (b.wp.mirrorGradBall && b.wp.oddMirrorDirichlet && !(b.ext && *b.ext))
            throw std::invalid_argument("solveMulti: the odd mirror ball needs the data's harmonic extensions");
        ctx[j] = b.tet ? std::make_unique<WalkCtx>(*b.tet, *b.cut, b.cut->nStrips(), GType::One, b.wp)
                       : std::make_unique<WalkCtx>(*b.dom, *b.cut, b.cut->nStrips(), GType::One, b.wp);
        const WalkCtx& c = *ctx[j];
        MultiSums& ms = out[j];
        ms.K = b.K;
        ms.nBatches = std::max(1, b.nBatches);
        ms.nPoints = b.pts->size();
        if (b.nWalks%ms.nBatches) throw std::invalid_argument("solveMulti: nWalks must be a multiple of nBatches");
        ms.walksPerBatch = b.nWalks/ms.nBatches;
        ms.batch.assign(ms.nPoints*size_t(ms.nBatches)*ms.stride(), 0.0);
        ms.uu.assign(ms.nPoints*size_t(b.K)*size_t(b.K), 0.0);
        active[j] = !(b.cut->nStrips() == 0 && !c.anyDirichletFace);
        bool all = c.anyDirichletFace;
        for (int f = 0; f < c.faces.nf; f++) all = all && c.dirichletFace[f];
        allAbsorb[j] = all;
        first[j + 1] = first[j] + (active[j] ? ms.nPoints : 0);
    }
    const size_t total = first[J];
    std::vector<long long> counted(total, 0), capped(total, 0), steps(total, 0);
    tbb::parallel_for(size_t(0), total, [&](size_t gi) {
        const size_t j = size_t(std::upper_bound(first.begin(), first.end(), gi) - first.begin()) - 1;
        const size_t i = gi - first[j];
        const MultiJob& job = jobs[j];
        const WalkCtx& c = *ctx[j];
        MultiSums& ms = out[j];
        const QueryPoints& pts = *job.pts;
        const int K = job.K;
        const size_t N = size_t(job.nWalks), B = size_t(ms.nBatches), NB = size_t(ms.walksPerBatch);
        const bool grad = pts.grad[i] != 0, absorb = allAbsorb[j] != 0;
        double u[kMaxK] = {}, phiY[kMaxK] = {}, phiMy[kMaxK] = {};   // zeroed like solveMultiOn's vectors
        double* uu = &ms.uu[i*size_t(K)*size_t(K)];
        OddMirror odd;
        for (size_t b = 0; b < B; b++) {
            double* s = &ms.batch[(i*B + b)*ms.stride()];
            double *su = s + 1, *sw = s + 1 + K, *swu = s + 4 + K;
            for (size_t w = b*NB; w < (b + 1)*NB; w++) {
                WalkState ws = makeInitialState(pts, int(N), i*N + w, 0);
                WalkEntry e;
                odd.face = -1;
                runWalkLean(c, grad, &odd, absorb, ws, &e);
                steps[gi] += e.steps;
                capped[gi] += e.status == WalkStatus::Capped;
                if (!countedWalk(e)) continue;
                counted[gi]++;
                if (e.status == WalkStatus::Landed) {
                    WalkTerminal T;
                    T.status = e.status;
                    T.face = e.face;
                    T.side = e.side;
                    T.x = e.face >= 0 ? c.onFacePlane(e.face, e.land) : e.land;
                    job.values(T, u);
                } else {
                    std::fill(u, u + K, 0.0);   // capped: value 0, as everywhere else
                }
                if (odd.face >= 0) {   // the first step crossed a Dirichlet face: H~(y) = Phi(y) + Phi(My) - H(My)
                    (*job.ext)(odd.y, phiY);
                    (*job.ext)(odd.my, phiMy);
                    for (int k = 0; k < K; k++) u[k] = phiY[k] + phiMy[k] - u[k];
                }
                const Vec3 wv = grad ? ws.w : Vec3::Zero();
                s[0] += 1;
                for (int a = 0; a < 3; a++) sw[a] += wv[a];
                for (int k = 0; k < K; k++) {
                    su[k] += u[k];
                    for (int a = 0; a < 3; a++) swu[3*k + a] += wv[a]*u[k];
                    for (int l = 0; l < K; l++) uu[k*K + l] += u[k]*u[l];
                }
            }
        }
    });
    const double sec = secondsSince(t0);
    if (stats) {
        stats->assign(J, SolveStats{});
        for (size_t j = 0; j < J; j++) {
            SolveStats& st = (*stats)[j];
            st.nPoints = int(out[j].nPoints);
            st.noCut = !active[j];
            if (active[j]) {
                for (size_t g = first[j]; g < first[j + 1]; g++) {
                    st.walksCounted += counted[g];
                    st.walksCapped += capped[g];
                    st.totalSteps += steps[g];
                }
                st.walksRequested = (long long)(out[j].nPoints*size_t(jobs[j].nWalks));
                st.walksEscaped = st.walksRequested - st.walksCounted;
            }
            st.seconds = sec;
        }
    }
    return out;
}

MultiSums solveMultiLean(const Tet& tet, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                         const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats,
                         const HarmonicExtension* ext)
{
    MultiJob jb;
    jb.tet = &tet;
    jb.cut = &cut;
    jb.pts = &pts;
    jb.nWalks = nWalks;
    jb.nBatches = nBatches;
    jb.K = K;
    jb.wp = wp;
    jb.values = values;
    jb.ext = ext;
    std::vector<SolveStats> st;
    std::vector<MultiSums> r = solveMultiLeanJobs({jb}, &st);
    if (stats) *stats = st[0];
    return std::move(r[0]);
}

MultiSums solveMultiLean(const ConvexDomain& dom, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                         const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats)
{
    MultiJob jb;
    jb.dom = &dom;
    jb.cut = &cut;
    jb.pts = &pts;
    jb.nWalks = nWalks;
    jb.nBatches = nBatches;
    jb.K = K;
    jb.wp = wp;
    jb.values = values;
    std::vector<SolveStats> st;
    std::vector<MultiSums> r = solveMultiLeanJobs({jb}, &st);
    if (stats) *stats = st[0];
    return std::move(r[0]);
}

// =====================================================================================
// Profiling only (partial_cut --ring-perf; no other mode calls this)
// =====================================================================================
#ifdef _WIN32
extern "C" __declspec(dllimport) void* __stdcall GetCurrentThread(void);
extern "C" __declspec(dllimport) unsigned long long __stdcall SetThreadAffinityMask(void* hThread, unsigned long long mask);
#endif

namespace {

struct PinThread {   // pins the calling thread to one logical processor for its lifetime
    unsigned long long old = 0;
    explicit PinThread(int cpu)
    {
#ifdef _WIN32
        if (cpu >= 0 && cpu < 64) old = SetThreadAffinityMask(GetCurrentThread(), 1ull << cpu);
#else
        (void)cpu;
#endif
    }
    ~PinThread()
    {
#ifdef _WIN32
        if (old) SetThreadAffinityMask(GetCurrentThread(), old);
#endif
    }
};

volatile double gProfileSink = 0;

} // namespace

WalkProfile profileWalks(const Tet* tet, const ConvexDomain* dom, const StripCut& cut, const QueryPoints& pts, int nWalks,
                         const WalkParams& wp, int K, const TerminalValues& values, const HarmonicExtension* ext,
                         size_t maxPoints, int cpu, bool components)
{
    using clk = std::chrono::steady_clock;
    const WalkCtx c = tet ? WalkCtx(*tet, cut, cut.nStrips(), GType::One, wp) : WalkCtx(*dom, cut, cut.nStrips(), GType::One, wp);
    if (c.oddMirror && !(ext && *ext)) throw std::invalid_argument("profileWalks: the odd mirror ball needs the extensions");
    WalkProfile pr;
    pr.cpu = cpu;
    PinThread pin(cpu);
    const size_t P = pts.size(), N = size_t(nWalks), M = std::max<size_t>(1, maxPoints);
    const size_t stride = std::max<size_t>(1, (P + M - 1)/M);
    std::vector<size_t> sel;
    for (size_t i = 0; i < P && sel.size() < M; i += stride) sel.push_back(i);
    pr.points = (long long)sel.size();
    pr.walks = pr.points*(long long)N;
    const size_t KK = size_t(K);
    std::vector<double> u(KK), phiY(KK), phiMy(KK), su(KK, 0.0), uu(KK*KK, 0.0), swu(3*KK, 0.0);
    double sw[3] = {0, 0, 0}, sn = 0, acc = 0;
    OddMirror odd;
    WalkControl ctl;
    ctl.odd = &odd;
    // solveMultiOn's scoring of one walk, verbatim
    auto score = [&](const WalkState& ws, const WalkEntry& e, bool grad) {
        if (!countedWalk(e)) return;
        if (e.status == WalkStatus::Landed) {
            WalkTerminal T;
            T.status = e.status;
            T.face = e.face;
            T.side = e.side;
            T.x = e.face >= 0 ? c.onFacePlane(e.face, e.land) : e.land;
            values(T, u.data());
        } else {
            std::fill(u.begin(), u.end(), 0.0);
        }
        if (odd.face >= 0) {
            (*ext)(odd.y, phiY.data());
            (*ext)(odd.my, phiMy.data());
            for (size_t k = 0; k < KK; k++) u[k] = phiY[k] + phiMy[k] - u[k];
        }
        const Vec3 w = grad ? ws.w : Vec3::Zero();
        sn += 1;
        for (int a = 0; a < 3; a++) sw[a] += w[a];
        for (size_t k = 0; k < KK; k++) {
            su[k] += u[k];
            for (int a = 0; a < 3; a++) swu[3*k + a] += w[a]*u[k];
            for (size_t l = 0; l < KK; l++) uu[k*KK + l] += u[k]*u[l];
        }
    };
    // warm-up (caches, branch predictors, the core's clock)
    for (size_t r = 0; r < std::min<size_t>(sel.size(), 8); r++)
        for (size_t j = 0; j < N; j++) {
            WalkState ws = makeInitialState(pts, int(N), sel[r]*N + j, 0);
            WalkEntry e;
            odd.face = -1;
            runWalk(c, pts.grad[sel[r]] != 0, ctl, ws, &e);
            score(ws, e, pts.grad[sel[r]] != 0);
        }
    auto t0 = clk::now();
    for (size_t i : sel) {
        const bool grad = pts.grad[i] != 0;
        for (size_t j = 0; j < N; j++) {
            WalkState ws = makeInitialState(pts, int(N), i*N + j, 0);
            WalkEntry e;
            odd.face = -1;
            runWalk(c, grad, ctl, ws, &e);
            score(ws, e, grad);
        }
    }
    pr.secFull = secondsSince(t0);
    t0 = clk::now();
    for (size_t i : sel) {
        const bool grad = pts.grad[i] != 0;
        for (size_t j = 0; j < N; j++) {
            WalkState ws = makeInitialState(pts, int(N), i*N + j, 0);
            WalkEntry e;
            odd.face = -1;
            runWalk(c, grad, ctl, ws, &e);
            pr.moves += e.steps;
            pr.landedCut += e.status == WalkStatus::Landed && e.face < 0;
            pr.landedFace += e.status == WalkStatus::Landed && e.face >= 0;
            pr.capped += e.status == WalkStatus::Capped;
            pr.escaped += e.status == WalkStatus::Escaped;
            pr.mirrored += odd.face >= 0;
        }
    }
    pr.secWalks = secondsSince(t0);
    pr.evals = pr.moves + pr.walks;
    // the lean step function: timed the same two ways, and compared walk by walk with runWalk
    bool allAbsorb = c.anyDirichletFace;
    for (int f = 0; f < c.faces.nf; f++) allAbsorb = allAbsorb && c.dirichletFace[f];
    pr.allAbsorb = allAbsorb;
    t0 = clk::now();
    for (size_t i : sel) {
        const bool grad = pts.grad[i] != 0;
        for (size_t j = 0; j < N; j++) {
            WalkState ws = makeInitialState(pts, int(N), i*N + j, 0);
            WalkEntry e;
            odd.face = -1;
            runWalkLean(c, grad, &odd, allAbsorb, ws, &e);
            score(ws, e, grad);
        }
    }
    pr.secFullLean = secondsSince(t0);
    long long leanMoves = 0;
    t0 = clk::now();
    for (size_t i : sel) {
        const bool grad = pts.grad[i] != 0;
        for (size_t j = 0; j < N; j++) {
            WalkState ws = makeInitialState(pts, int(N), i*N + j, 0);
            WalkEntry e;
            odd.face = -1;
            runWalkLean(c, grad, &odd, allAbsorb, ws, &e);
            leanMoves += e.steps;
        }
    }
    pr.secWalksLean = secondsSince(t0);
    acc += double(leanMoves);
    for (size_t i : sel) {
        const bool grad = pts.grad[i] != 0;
        for (size_t j = 0; j < N; j++) {
            WalkState w1 = makeInitialState(pts, int(N), i*N + j, 0), w2 = w1;
            WalkEntry e1, e2;
            odd.face = -1;
            runWalk(c, grad, ctl, w1, &e1);
            const OddMirror o1 = odd;
            odd.face = -1;
            runWalkLean(c, grad, &odd, allAbsorb, w2, &e2);
            const bool same = e1.status == e2.status && e1.steps == e2.steps && e1.side == e2.side && e1.face == e2.face &&
                              sameBits(&e1.land, &e2.land, sizeof(Vec3)) && sameBits(&w1.w, &w2.w, sizeof(Vec3)) &&
                              sameBits(&w1.x, &w2.x, sizeof(Vec3)) && o1.face == odd.face &&
                              (o1.face < 0 || (sameBits(&o1.y, &odd.y, sizeof(Vec3)) && sameBits(&o1.my, &odd.my, sizeof(Vec3))));
            pr.leanMismatch += !same;
        }
    }
    acc += sn + su[0] + sw[0] + swu[0] + uu[0];
    if (!components) {
        gProfileSink = acc;
        return pr;
    }

    // record every evaluated step state, and per walk its start, last state and terminal
    std::vector<Vec3> X;
    std::vector<int> F;
    std::vector<double> D2, R;
    X.reserve(size_t(pr.evals));
    F.reserve(size_t(pr.evals));
    D2.reserve(size_t(pr.evals));
    struct WalkRec {
        size_t i = 0, j = 0, first = 0;
        WalkEntry e;
        StripCut::Closest qLast;
        Vec3 xLast = Vec3::Zero();
    };
    std::vector<WalkRec> wr;
    wr.reserve(size_t(pr.walks));
    StripCut::Closest qPrev;
    const std::function<void(const WalkState&, const StripCut::Closest&)> obs = [&](const WalkState& st,
                                                                                     const StripCut::Closest& q) {
        X.push_back(st.x);
        F.push_back(st.onFace);
        D2.push_back(q.d2);
        qPrev = q;
    };
    WalkControl rc;
    rc.odd = &odd;
    rc.observe = &obs;
    std::vector<WalkTerminal> term;
    for (size_t i : sel) {
        const bool grad = pts.grad[i] != 0;
        for (size_t j = 0; j < N; j++) {
            WalkRec w;
            w.i = i;
            w.j = j;
            w.first = X.size();
            WalkState ws = makeInitialState(pts, int(N), i*N + j, 0);
            odd.face = -1;
            runWalk(c, grad, rc, ws, &w.e);
            w.qLast = qPrev;
            w.xLast = X.back();
            if (w.e.status == WalkStatus::Landed) {
                WalkTerminal T;
                T.status = w.e.status;
                T.face = w.e.face;
                T.side = w.e.side;
                T.x = w.e.face >= 0 ? c.onFacePlane(w.e.face, w.e.land) : w.e.land;
                term.push_back(T);
            }
            wr.push_back(w);
        }
    }
    const size_t S = X.size();
    R.resize(S);
    for (size_t s = 0; s < S; s++) {
        int f;
        const double dF = c.dirichletFaceDist(X[s], f), dC = std::sqrt(D2[s]);
        R[s] = dF < dC ? dF : dC;
    }
    auto nsPer = [](clk::time_point t, size_t n) { return n ? 1e9*secondsSince(t)/double(n) : 0.0; };
    // each part alone, over the recorded states
    t0 = clk::now();
    for (size_t s = 0; s < S; s++) {
        const StripCut::Closest q = c.cut->closest(X[s], c.nStrips);
        acc += double(q.tri) + q.q.x();
    }
    pr.nsClosest = nsPer(t0, S);
    t0 = clk::now();
    for (size_t s = 0; s < S; s++) {
        int f;
        acc += c.dirichletFaceDist(X[s], f) + double(f);
    }
    pr.nsFaceDist = nsPer(t0, S);
    t0 = clk::now();
    for (size_t s = 0; s < S; s++) acc += std::sqrt(D2[s]);
    pr.nsSqrt = nsPer(t0, S);
    {
        WalkState tmp;
        emptyBox(tmp);
        t0 = clk::now();
        for (size_t s = 0; s < S; s++) addBall(tmp, X[s], R[s]);
        pr.nsAddBall = nsPer(t0, S);
        acc += tmp.lo[0] + tmp.hi[2];
    }
    {
        pcg32 rng(0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL);
        t0 = clk::now();
        for (size_t s = 0; s < S; s++) acc += sampleSphere(rng).x();
        pr.nsSphere = nsPer(t0, S);
    }
    {
        pcg32 rng(0x1234567ULL, 0x89abcdefULL);
        std::vector<Vec3> dirs(S);
        for (size_t s = 0; s < S; s++) dirs[s] = sampleSphere(rng);
        WalkState tmp;
        t0 = clk::now();
        for (size_t s = 0; s < S; s++) {
            tmp.x = X[s];
            tmp.onFace = F[s];
            starStep(c.faces, R[s], dirs[s], tmp);
            acc += tmp.x.x() + double(tmp.onFace);
        }
        pr.nsStar = nsPer(t0, S);
    }
    t0 = clk::now();
    for (size_t s = 0; s < S; s++) acc += c.faces.outside(X[s], c.escTol) ? 1.0 : 0.0;
    pr.nsOutside = nsPer(t0, S);
    const size_t W = wr.size();
    t0 = clk::now();
    for (const WalkRec& w : wr) {
        const WalkState ws = makeInitialState(pts, int(N), w.i*N + w.j, 0);
        acc += ws.x.x();
    }
    pr.nsInit = nsPer(t0, W);
    {
        pcg32 rng(0x2468aceULL, 0x13579bdfULL);
        t0 = clk::now();
        for (const WalkRec& w : wr) {
            const Vec3& x = pts.x[w.i];
            const Vec3 dir = sampleSphere(rng);
            const double R0 = kGradFirstBallShrink*std::min(R[w.first], c.faces.distInside(x));
            const Vec3 g = (3.0/R0)*dir;
            acc += g.x() + (x + R0*dir).y();
        }
        pr.nsFirstBall = nsPer(t0, W);
    }
    if (c.mirrorGradBall) {
        t0 = clk::now();
        for (const WalkRec& w : wr) {
            double R0 = kGradFirstBallShrink*std::min(R[w.first], c.faces.distInside(pts.x[w.i]));
            acc += double(mirrorGradBall(c, pts.x[w.i], std::sqrt(D2[w.first]), R0)) + R0;
        }
        pr.nsMirror = nsPer(t0, W);
    }
    {
        WalkState tmp;
        emptyBox(tmp);
        size_t n = 0;
        t0 = clk::now();
        for (const WalkRec& w : wr) {
            if (w.e.status != WalkStatus::Landed) continue;
            WalkEntry e;
            tmp.x = w.xLast;
            tmp.steps = w.e.steps;
            finishWalk(c, tmp, WalkStatus::Landed, w.e.face >= 0 ? nullptr : &w.qLast, &e, w.e.face);
            acc += e.value;
            n++;
        }
        pr.nsFinish = nsPer(t0, n);
    }
    std::vector<double> uAll(term.size()*KK);
    t0 = clk::now();
    for (size_t t = 0; t < term.size(); t++) values(term[t], &uAll[t*KK]);
    pr.nsValues = nsPer(t0, term.size());
    {
        pcg32 rng(0x777ULL, 0x999ULL);
        std::vector<Vec3> ws(term.size());
        for (Vec3& w : ws) w = 3.0*sampleSphere(rng);
        t0 = clk::now();
        for (size_t t = 0; t < term.size(); t++) {
            const double* uk = &uAll[t*KK];
            const Vec3& w = ws[t];
            sn += 1;
            for (int a = 0; a < 3; a++) sw[a] += w[a];
            for (size_t k = 0; k < KK; k++) {
                su[k] += uk[k];
                for (int a = 0; a < 3; a++) swu[3*k + a] += w[a]*uk[k];
                for (size_t l = 0; l < KK; l++) uu[k*KK + l] += uk[k]*uk[l];
            }
        }
        pr.nsAccum = nsPer(t0, term.size());
    }
    acc += sn + su[0] + sw[0] + swu[0] + uu[0];
    gProfileSink = acc;
    return pr;
}

double timeWalksParallel(const Tet* tet, const ConvexDomain* dom, const StripCut& cut, const QueryPoints& pts, int nWalks,
                         const WalkParams& wp, bool lean, long long* moves)
{
    const WalkCtx c = tet ? WalkCtx(*tet, cut, cut.nStrips(), GType::One, wp) : WalkCtx(*dom, cut, cut.nStrips(), GType::One, wp);
    bool allAbsorb = c.anyDirichletFace;
    for (int f = 0; f < c.faces.nf; f++) allAbsorb = allAbsorb && c.dirichletFace[f];
    const size_t P = pts.size(), N = size_t(nWalks);
    std::vector<long long> steps(P, 0);
    const auto t0 = std::chrono::steady_clock::now();
    tbb::parallel_for(size_t(0), P, [&](size_t i) {
        const bool grad = pts.grad[i] != 0;
        OddMirror odd;
        WalkControl ctl;
        ctl.odd = &odd;
        for (size_t j = 0; j < N; j++) {
            WalkState ws = makeInitialState(pts, int(N), i*N + j, 0);
            WalkEntry e;
            odd.face = -1;
            if (lean) runWalkLean(c, grad, &odd, allAbsorb, ws, &e);
            else runWalk(c, grad, ctl, ws, &e);
            steps[i] += e.steps;
        }
    });
    const double sec = secondsSince(t0);
    if (moves) {
        *moves = 0;
        for (long long s : steps) *moves += s;
    }
    return sec;
}

// =====================================================================================
// General crack geometry and the mesh walks (opt-in: --ring-stress)
// =====================================================================================
CrackMesh::TriPre CrackMesh::pre(const Tri& t)
{
    TriPre r;
    r.a = t.a;
    r.b = t.b;
    r.c = t.c;
    r.ab = t.b - t.a;
    r.ac = t.c - t.a;
    r.abab = r.ab.dot(r.ab);
    r.abac = r.ab.dot(r.ac);
    r.acac = r.ac.dot(r.ac);
    return r;
}

double CrackMesh::triDist2(const TriPre& t, const Vec3& p, Vec3& q)
{
    // Ericson, Real-Time Collision Detection 5.1.5: the closest point by the Voronoi regions of the triangle, with
    // ab.bp = ab.ap - ab.ab etc. from the precomputed dot products (two dot products per call)
    const Vec3 ap = p - t.a;
    const double d1 = t.ab.dot(ap), d2 = t.ac.dot(ap);
    if (d1 <= 0 && d2 <= 0) {
        q = t.a;
        return ap.squaredNorm();
    }
    const double d3 = d1 - t.abab, d4 = d2 - t.abac;
    if (d3 >= 0 && d4 <= d3) {
        q = t.b;
        return (p - q).squaredNorm();
    }
    const double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        q = t.a + (d1/(d1 - d3))*t.ab;
        return (p - q).squaredNorm();
    }
    const double d5 = d1 - t.abac, d6 = d2 - t.acac;
    if (d6 >= 0 && d5 <= d6) {
        q = t.c;
        return (p - q).squaredNorm();
    }
    const double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        q = t.a + (d2/(d2 - d6))*t.ac;
        return (p - q).squaredNorm();
    }
    const double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        q = t.b + ((d4 - d3)/((d4 - d3) + (d5 - d6)))*(t.c - t.b);
        return (p - q).squaredNorm();
    }
    const double den = 1.0/(va + vb + vc);
    q = t.a + (vb*den)*t.ab + (vc*den)*t.ac;
    return (p - q).squaredNorm();
}

CrackMesh::CrackMesh(std::vector<Tri> tris) : tris_(std::move(tris))
{
    for (Tri& t : tris_)
        if (t.n.isZero()) {
            const Vec3 nn = (t.b - t.a).cross(t.c - t.a);
            t.n = nn/nn.norm();
        }
    if (tris_.empty()) return;
    pre_.resize(tris_.size());
    for (size_t i = 0; i < tris_.size(); i++) pre_[i] = pre(tris_[i]);
    std::vector<Vec3> cen(tris_.size());
    for (size_t i = 0; i < tris_.size(); i++) cen[i] = (tris_[i].a + tris_[i].b + tris_[i].c)/3.0;
    order_.resize(tris_.size());
    for (size_t i = 0; i < tris_.size(); i++) order_[i] = int(i);
    nodes_.reserve(2*tris_.size()/3 + 2);
    build(order_, 0, int(tris_.size()), cen);
}

int CrackMesh::build(std::vector<int>& idx, int first, int count, const std::vector<Vec3>& cen)
{
    const int id = int(nodes_.size());
    nodes_.push_back(Node{});
    Node nd;
    // the frame: the mean (area-weighted) normal of the node's triangles if they roughly agree, else the coordinate axes
    Vec3 nsum = Vec3::Zero();
    double asum = 0;
    for (int k = first; k < first + count; k++) {
        const Tri& t = tris_[size_t(idx[size_t(k)])];
        const Vec3 c = (t.b - t.a).cross(t.c - t.a);
        nsum += c.norm()*t.n;
        asum += c.norm();
    }
    if (nsum.norm() > 0.7*asum && asum > 0) {
        nd.e[2] = nsum.normalized();
        const Vec3 tt = std::fabs(nd.e[2].x()) < 0.6 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
        nd.e[0] = (tt - tt.dot(nd.e[2])*nd.e[2]).normalized();
        nd.e[1] = nd.e[2].cross(nd.e[0]);
    } else {
        nd.e[0] = Vec3(1, 0, 0);
        nd.e[1] = Vec3(0, 1, 0);
        nd.e[2] = Vec3(0, 0, 1);
    }
    for (int a = 0; a < 3; a++) {
        nd.lo[a] = std::numeric_limits<double>::infinity();
        nd.hi[a] = -std::numeric_limits<double>::infinity();
    }
    Vec3 clo = Vec3::Constant(std::numeric_limits<double>::infinity()), chi = -clo;
    for (int k = first; k < first + count; k++) {
        const Tri& t = tris_[size_t(idx[size_t(k)])];
        for (const Vec3* v : {&t.a, &t.b, &t.c})
            for (int a = 0; a < 3; a++) {
                const double c = nd.e[a].dot(*v);
                nd.lo[a] = std::min(nd.lo[a], c);
                nd.hi[a] = std::max(nd.hi[a], c);
            }
        clo = clo.cwiseMin(cen[size_t(idx[size_t(k)])]);
        chi = chi.cwiseMax(cen[size_t(idx[size_t(k)])]);
    }
    for (int a = 0; a < 3; a++) {   // padding: the box contains its triangles whatever the rounding of the projections
        nd.lo[a] -= 1e-12;
        nd.hi[a] += 1e-12;
    }
    if (count <= 4) {
        nd.first = first;
        nd.count = count;
        nodes_[size_t(id)] = nd;
        return id;
    }
    int axis = 0;
    const Vec3 ext = chi - clo;
    if (ext[1] > ext[axis]) axis = 1;
    if (ext[2] > ext[axis]) axis = 2;
    const int mid = first + count/2;
    std::nth_element(idx.begin() + first, idx.begin() + mid, idx.begin() + first + count, [&](int i, int j) {
        const double ci = cen[size_t(i)][axis], cj = cen[size_t(j)][axis];
        return ci < cj || (ci == cj && i < j);
    });
    nodes_[size_t(id)] = nd;
    const int l = build(idx, first, mid - first, cen);
    const int r = build(idx, mid, first + count - mid, cen);
    nodes_[size_t(id)].left = l;
    nodes_[size_t(id)].right = r;
    return id;
}

double CrackMesh::nodeDist2(const Node& nd, const Vec3& x)
{
    double s = 0;
    for (int a = 0; a < 3; a++) {
        const double c = nd.e[a].dot(x);
        const double v = c < nd.lo[a] ? nd.lo[a] - c : c > nd.hi[a] ? c - nd.hi[a] : 0.0;
        s += v*v;
    }
    return s;
}

CrackMesh::Hit CrackMesh::closest(const Vec3& x, int hint) const
{
    Hit best;
    if (nodes_.empty()) return best;
    if (hint >= 0 && hint < int(pre_.size())) {   // the hint's exact distance: an upper bound that prunes from the start
        best.d2 = triDist2(pre_[size_t(hint)], x, best.q);
        best.tri = hint;
    }
    struct Entry {
        int node;
        double d2;
    };
    Entry stack[128];
    int sp = 0;
    stack[sp++] = {0, nodeDist2(nodes_[0], x)};
    while (sp > 0) {
        const Entry e = stack[--sp];
        if (e.d2 > best.d2*(1.0 + 1e-12)) continue;
        const Node& nd = nodes_[size_t(e.node)];
        if (nd.left < 0) {
            for (int k = nd.first; k < nd.first + nd.count; k++) {
                const int t = order_[size_t(k)];
                if (t == best.tri) continue;
                Vec3 q;
                const double d2 = triDist2(pre_[size_t(t)], x, q);
                if (d2 < best.d2 || (d2 == best.d2 && t < best.tri)) {
                    best.d2 = d2;
                    best.tri = t;
                    best.q = q;
                }
            }
            continue;
        }
        const Node &L = nodes_[size_t(nd.left)], &R = nodes_[size_t(nd.right)];
        const double dl = nodeDist2(L, x), dr = nodeDist2(R, x);
        if (sp + 2 > 128) throw std::logic_error("CrackMesh::closest: tree too deep");
        if (dl <= dr) {   // the nearer child is popped first
            stack[sp++] = {nd.right, dr};
            stack[sp++] = {nd.left, dl};
        } else {
            stack[sp++] = {nd.left, dl};
            stack[sp++] = {nd.right, dr};
        }
    }
    return best;
}

CrackMesh::Hit CrackMesh::closestBruteForce(const Vec3& x) const
{
    Hit best;
    for (size_t t = 0; t < pre_.size(); t++) {
        Vec3 q;
        const double d2 = triDist2(pre_[t], x, q);
        if (d2 < best.d2) {   // index order: a tie keeps the lower index
            best.d2 = d2;
            best.tri = int(t);
            best.q = q;
        }
    }
    return best;
}

double CrackMesh::distance(const Vec3& x) const { return std::sqrt(closest(x).d2); }

double CrackMesh::signedDistance(const Vec3& x) const
{
    const Hit h = closest(x);
    if (h.tri < 0) return std::numeric_limits<double>::infinity();
    return side(h, x)*std::sqrt(h.d2);
}

double CrackMesh::area() const
{
    double s = 0;
    for (const Tri& t : tris_) s += 0.5*(t.b - t.a).cross(t.c - t.a).norm();
    return s;
}

void CrackMesh::clipPolygon(std::vector<Vec3>& p, std::vector<uint8_t>& f, const Vec3& n, double c, uint8_t newFlag)
{
    std::vector<Vec3> q;
    std::vector<uint8_t> g;
    const size_t m = p.size();
    for (size_t i = 0; i < m; i++) {
        const Vec3 &A = p[i], &B = p[(i + 1)%m];
        const double fa = n.dot(A) - c, fb = n.dot(B) - c;
        // the crossing point, computed from the lexicographically smaller end so that an edge shared by two polygons
        // gives the same point from both
        auto cross = [&]() -> Vec3 {
            const bool aFirst = A.x() < B.x() || (A.x() == B.x() && (A.y() < B.y() || (A.y() == B.y() && A.z() <= B.z())));
            return aFirst ? Vec3(A + (fa/(fa - fb))*(B - A)) : Vec3(B + (fb/(fb - fa))*(A - B));
        };
        if (fa <= 0) {
            q.push_back(A);
            g.push_back(f[i]);
            if (fb > 0) {
                q.push_back(cross());
                g.push_back(newFlag);
            }
        } else if (fb <= 0) {
            q.push_back(cross());
            g.push_back(f[i]);
        }
    }
    p.swap(q);
    f.swap(g);
}

void CrackMesh::fanPolygon(const std::vector<Vec3>& p, const std::vector<uint8_t>& f, int crack, const Vec3& nParent,
                           double minArea, std::vector<Tri>& out)
{
    const size_t m = p.size();
    if (m < 3) return;
    for (size_t k = 1; k + 1 < m; k++) {
        Tri t;
        t.a = p[0];
        t.b = p[k];
        t.c = p[k + 1];
        if (!(0.5*(t.b - t.a).cross(t.c - t.a).norm() > minArea)) continue;
        t.n = nParent;
        t.crack = crack;
        t.edge[0] = k == 1 ? f[0] : uint8_t(0);
        t.edge[1] = f[k];
        t.edge[2] = k + 2 == m ? f[m - 1] : uint8_t(0);
        out.push_back(t);
    }
}

CrackMesh CrackMesh::clipped(const std::vector<Vec3>& n, const std::vector<double>& c, const std::vector<uint8_t>& planeFlag,
                             double minArea) const
{
    std::vector<Tri> out;
    std::vector<Vec3> p;
    std::vector<uint8_t> f;
    for (const Tri& t : tris_) {
        p = {t.a, t.b, t.c};
        f = {t.edge[0], t.edge[1], t.edge[2]};
        for (size_t k = 0; k < n.size() && p.size() >= 3; k++) clipPolygon(p, f, n[k], c[k], planeFlag[k]);
        fanPolygon(p, f, t.crack, t.n, minArea, out);
    }
    return CrackMesh(std::move(out));
}

CrackMesh CrackMesh::mapped(const Eigen::Matrix3d& A, const Vec3& b) const
{
    const Eigen::Matrix3d AinvT = A.inverse().transpose();
    std::vector<Tri> out = tris_;
    for (Tri& t : out) {
        t.a = A*t.a + b;
        t.b = A*t.b + b;
        t.c = A*t.c + b;
        const Vec3 nn = AinvT*t.n;
        t.n = nn/nn.norm();
    }
    return CrackMesh(std::move(out));
}

namespace {

// The fixed inputs of the mesh step function (WalkCtx's, with a CrackMesh as the crack)
struct MeshCtx {
    const CrackMesh* mesh = nullptr;
    FacePlanes faces;
    bool dirichletFace[kMaxDomainPlanes] = {false, false, false, false, false, false, false, false};
    bool anyDirichletFace = false, allAbsorb = false, mirrorGradBall = false, oddMirror = false;
    double shell = 0, escTol = 0;
    int maxLen = 0;
    double dirichletFaceDist(const Vec3& x, int& face) const
    {
        face = -1;
        if (!anyDirichletFace) return std::numeric_limits<double>::infinity();
        double dmin = std::numeric_limits<double>::max();
        for (int f = 0; f < faces.nf; f++) {
            if (!dirichletFace[f]) continue;
            const double d = faces.c[f] - faces.n[f].dot(x);
            if (d < dmin) {
                dmin = d;
                face = f;
            }
        }
        return dmin;
    }
    Vec3 onFacePlane(int f, const Vec3& x) const { return x - (faces.n[f].dot(x) - faces.c[f])*faces.n[f]; }
};

MeshCtx meshCtx(const MeshJob& b)
{
    MeshCtx c;
    c.mesh = b.mesh;
    const WalkParams& wp = b.wp;
    if (b.tet) {
        c.faces = FacePlanes(*b.tet);
        for (int f = 0; f < 4; f++) {
            c.dirichletFace[f] = wp.faceBC[f] == FaceBC::Dirichlet;
            c.anyDirichletFace = c.anyDirichletFace || c.dirichletFace[f];
        }
        c.shell = wp.absorbShellRel*b.tet->L;
        c.escTol = 1e-12*b.tet->L;
    } else {
        if (wp.dirichletFaces || wp.anyDirichletFace() || wp.faceData || wp.oddMirrorDirichlet)
            throw std::invalid_argument("solveMultiMeshJobs: a ConvexDomain has Neumann faces only");
        c.faces = FacePlanes(*b.dom);
        c.shell = wp.absorbShellRel*b.dom->L;
        c.escTol = 1e-12*b.dom->L;
    }
    bool all = c.anyDirichletFace;
    for (int f = 0; f < c.faces.nf; f++) all = all && c.dirichletFace[f];
    c.allAbsorb = all;
    c.mirrorGradBall = wp.mirrorGradBall;
    c.oddMirror = wp.mirrorGradBall && wp.oddMirrorDirichlet;
    c.maxLen = wp.maxWalkLength;
    return c;
}

// mirrorGradBall with the mesh's distance at the mirror point
int mirrorGradBallMesh(const MeshCtx& c, const Vec3& x, double dCut, double& R0)
{
    int f = -1;
    double df = std::numeric_limits<double>::infinity();
    for (int g = 0; g < c.faces.nf; g++) {
        if (c.dirichletFace[g] && !c.oddMirror) continue;
        const double dg = c.faces.c[g] - c.faces.n[g].dot(x);
        if (dg < df) {
            df = dg;
            f = g;
        }
    }
    if (f < 0) return -1;
    const Vec3 xm = x + 2.0*df*c.faces.n[f];
    double R = std::min(dCut, std::sqrt(c.mesh->closest(xm).d2));
    for (int g = 0; g < c.faces.nf; g++) {
        if (g == f) continue;
        const double a = c.faces.c[g] - c.faces.n[g].dot(x), b = c.faces.c[g] - c.faces.n[g].dot(xm);
        if (!(b > 0)) return -1;
        R = std::min(R, std::min(a, b));
    }
    R *= kGradFirstBallShrink;
    if (!(R > R0)) return -1;
    R0 = R;
    return f;
}

// the lean step function with the mesh crack; fills the terminal (status, face or crack data) and returns the step count
int runWalkMesh(const MeshCtx& c, bool gradPoint, OddMirror* odd, WalkState& st, MeshTerminal& T)
{
    int hint = -1;   // the previous step's closest triangle (pruning only: the closest point does not depend on it)
    for (;;) {
        if (st.steps >= c.maxLen) {
            T.status = WalkStatus::Capped;
            T.xEnd = st.x;
            return st.steps;
        }
        const CrackMesh::Hit q = c.mesh->closest(st.x, hint);
        hint = q.tri;
        int fF = -1;
        const double dF = c.dirichletFaceDist(st.x, fF);
        const double dCut = std::sqrt(q.d2);
        const bool faceCloser = dF < dCut;
        const double d = faceCloser ? dF : dCut;
        if (d < c.shell) {
            T.status = WalkStatus::Landed;
            T.xEnd = st.x;
            if (faceCloser) {
                T.face = fF;
                T.side = 0;
                T.x = c.onFacePlane(fF, st.x);
            } else {
                T.face = -1;
                T.tri = q.tri;
                T.crack = c.mesh->tri(q.tri).crack;
                T.side = c.mesh->side(q, st.x);
                T.x = q.q;
            }
            return st.steps;
        }
        const Vec3 dir = sampleSphere(st.rng);
        if (gradPoint && st.steps == 0) {
            double R0 = kGradFirstBallShrink*std::min(d, c.faces.distInside(st.x));
            const int mirror = c.mirrorGradBall ? mirrorGradBallMesh(c, st.x, dCut, R0) : -1;
            st.w = (3.0/R0)*dir;
            st.x = st.x + R0*dir;
            if (mirror >= 0) {
                const double h = c.faces.n[mirror].dot(st.x) - c.faces.c[mirror];
                if (h > 0) {
                    const Vec3 y = st.x;
                    st.x = st.x - 2.0*h*c.faces.n[mirror];
                    if (c.dirichletFace[mirror]) {
                        odd->face = mirror;
                        odd->y = y;
                        odd->my = st.x;
                    }
                }
            }
            st.onFace = -1;
        } else if (c.allAbsorb) {
            st.x = st.x + d*dir;
            st.onFace = -1;
        } else {
            starStep(c.faces, d, dir, st);
        }
        st.steps++;
        if (st.onFace < 0 && c.faces.outside(st.x, c.escTol)) {
            T.status = WalkStatus::Escaped;
            T.xEnd = st.x;
            return st.steps;
        }
    }
}

} // namespace

std::vector<MultiSums> solveMultiMeshJobs(const std::vector<MeshJob>& jobs, std::vector<SolveStats>* stats)
{
    constexpr int kMaxK = 16;
    const auto t0 = std::chrono::steady_clock::now();
    const size_t J = jobs.size();
    std::vector<MeshCtx> ctx(J);
    std::vector<MultiSums> out(J);
    std::vector<char> active(J, 0);
    std::vector<size_t> first(J + 1, 0);
    static const CrackMesh emptyMesh;
    for (size_t j = 0; j < J; j++) {
        const MeshJob& b = jobs[j];
        if (!b.pts || (!b.tet == !b.dom)) throw std::invalid_argument("solveMultiMeshJobs: a job needs points and one domain");
        if (b.K < 0 || b.K > kMaxK) throw std::invalid_argument("solveMultiMeshJobs: K out of range");
        if (b.wp.mirrorGradBall && b.wp.oddMirrorDirichlet && !(b.ext && *b.ext))
            throw std::invalid_argument("solveMultiMeshJobs: the odd mirror ball needs the data's harmonic extensions");
        ctx[j] = meshCtx(b);
        if (!ctx[j].mesh) ctx[j].mesh = &emptyMesh;
        MultiSums& ms = out[j];
        ms.K = b.K;
        ms.nBatches = std::max(1, b.nBatches);
        ms.nPoints = b.pts->size();
        if (b.nWalks%ms.nBatches) throw std::invalid_argument("solveMultiMeshJobs: nWalks must be a multiple of nBatches");
        ms.walksPerBatch = b.nWalks/ms.nBatches;
        ms.batch.assign(ms.nPoints*size_t(ms.nBatches)*ms.stride(), 0.0);
        ms.uu.assign(ms.nPoints*size_t(b.K)*size_t(b.K), 0.0);
        active[j] = !(ctx[j].mesh->empty() && !ctx[j].anyDirichletFace);
        first[j + 1] = first[j] + (active[j] ? ms.nPoints : 0);
    }
    const size_t total = first[J];
    std::vector<long long> counted(total, 0), capped(total, 0), steps(total, 0);
    tbb::parallel_for(size_t(0), total, [&](size_t gi) {
        const size_t j = size_t(std::upper_bound(first.begin(), first.end(), gi) - first.begin()) - 1;
        const size_t i = gi - first[j];
        const MeshJob& job = jobs[j];
        const MeshCtx& c = ctx[j];
        MultiSums& ms = out[j];
        const QueryPoints& pts = *job.pts;
        const int K = job.K;
        const size_t N = size_t(job.nWalks), B = size_t(ms.nBatches), NB = size_t(ms.walksPerBatch);
        const bool grad = pts.grad[i] != 0;
        double u[kMaxK] = {}, phiY[kMaxK] = {}, phiMy[kMaxK] = {};
        double* uu = &ms.uu[i*size_t(K)*size_t(K)];
        OddMirror odd;
        for (size_t b = 0; b < B; b++) {
            double* s = &ms.batch[(i*B + b)*ms.stride()];
            double *su = s + 1, *sw = s + 1 + K, *swu = s + 4 + K;
            for (size_t w = b*NB; w < (b + 1)*NB; w++) {
                WalkState ws = makeInitialState(pts, int(N), i*N + w, 0);
                MeshTerminal T;
                odd.face = -1;
                steps[gi] += runWalkMesh(c, grad, &odd, ws, T);
                capped[gi] += T.status == WalkStatus::Capped;
                if (T.status == WalkStatus::Escaped) continue;
                counted[gi]++;
                if (T.status == WalkStatus::Landed) job.values(T, u);
                else std::fill(u, u + K, 0.0);   // capped: value 0, as everywhere else
                if (odd.face >= 0) {
                    (*job.ext)(odd.y, phiY);
                    (*job.ext)(odd.my, phiMy);
                    for (int k = 0; k < K; k++) u[k] = phiY[k] + phiMy[k] - u[k];
                }
                const Vec3 wv = grad ? ws.w : Vec3::Zero();
                s[0] += 1;
                for (int a = 0; a < 3; a++) sw[a] += wv[a];
                for (int k = 0; k < K; k++) {
                    su[k] += u[k];
                    for (int a = 0; a < 3; a++) swu[3*k + a] += wv[a]*u[k];
                    for (int l = 0; l < K; l++) uu[k*K + l] += u[k]*u[l];
                }
            }
        }
    });
    const double sec = secondsSince(t0);
    if (stats) {
        stats->assign(J, SolveStats{});
        for (size_t j = 0; j < J; j++) {
            SolveStats& st = (*stats)[j];
            st.nPoints = int(out[j].nPoints);
            st.noCut = !active[j];
            if (active[j]) {
                for (size_t g = first[j]; g < first[j + 1]; g++) {
                    st.walksCounted += counted[g];
                    st.walksCapped += capped[g];
                    st.totalSteps += steps[g];
                }
                st.walksRequested = (long long)(out[j].nPoints*size_t(jobs[j].nWalks));
                st.walksEscaped = st.walksRequested - st.walksCounted;
            }
            st.seconds = sec;
        }
    }
    return out;
}
