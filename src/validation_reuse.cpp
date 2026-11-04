// Validation of incremental walk reuse under an advancing cut. Everything printed also goes into
// the reuse section of results.txt. --compare-options (at the end) writes its own file.
#include "validation_reuse.h"

#include "validation.h"
#include "walk_cache.h"

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/parallel_for.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <memory>
#include <random>
#include <thread>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <psapi.h>   // GetProcessMemoryInfo (K32GetProcessMemoryInfo in kernel32)
#include <intrin.h>  // __cpuid
#endif

namespace {

// printf to stdout and into a string (written to results.txt at the end)
class Report {
public:
    void operator()(const char* fmt, ...)
    {
        char buf[4096];
        va_list ap;
        va_start(ap, fmt);
        std::vsnprintf(buf, sizeof buf, fmt, ap);
        va_end(ap);
        std::fputs(buf, stdout);
        std::fflush(stdout);
        text += buf;
    }
    void put(const std::string& s)   // unformatted, any length
    {
        std::fputs(s.c_str(), stdout);
        std::fflush(stdout);
        text += s;
    }
    std::string text;
};

double secondsSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

double meanOf(const std::vector<double>& v)
{
    double s = 0;
    for (double x : v) s += x;
    return v.empty() ? 0.0 : s/double(v.size());
}

double sampleStd(const std::vector<double>& v)
{
    if (v.size() < 2) return 0.0;
    double m = meanOf(v), s = 0;
    for (double x : v) s += (x - m)*(x - m);
    return std::sqrt(s/double(v.size() - 1));
}

double medianOf(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return n%2 ? v[n/2] : 0.5*(v[n/2 - 1] + v[n/2]);
}

// least-squares slope of log(y) against log(x); NaN if any y <= 0
double logLogSlope(const std::vector<double>& x, const std::vector<double>& y)
{
    std::vector<double> lx, ly;
    for (size_t i = 0; i < x.size(); i++) {
        if (y[i] <= 0) return std::nan("");
        lx.push_back(std::log(x[i]));
        ly.push_back(std::log(y[i]));
    }
    double mx = meanOf(lx), my = meanOf(ly), sxy = 0, sxx = 0;
    for (size_t i = 0; i < lx.size(); i++) {
        sxy += (lx[i] - mx)*(ly[i] - my);
        sxx += (lx[i] - mx)*(lx[i] - mx);
    }
    return sxy/sxx;
}

std::string fmt3(const Vec3& v)
{
    char buf[96];
    std::snprintf(buf, sizeof buf, "(%.4f, %.4f, %.4f)", v.x(), v.y(), v.z());
    return buf;
}

bool sameBits(const void* a, const void* b, size_t n) { return std::memcmp(a, b, n) == 0; }

bool sameState(const WalkState& a, const WalkState& b)
{
    return sameBits(a.x.data(), b.x.data(), 3*sizeof(double)) && a.onFace == b.onFace && a.steps == b.steps &&
           a.rng == b.rng && sameBits(a.w.data(), b.w.data(), 3*sizeof(double)) && sameBits(a.lo, b.lo, sizeof a.lo) &&
           sameBits(a.hi, b.hi, sizeof a.hi);
}

// distance from every point to the cut (strips [0, nStrips)); +inf without a cut
std::vector<double> cutDistances(const StripCut& cut, const std::vector<Vec3>& pts)
{
    std::vector<double> d(pts.size());
    const int ns = cut.nStrips();
    tbb::parallel_for(size_t(0), pts.size(), [&](size_t i) { d[i] = std::sqrt(cut.closest(pts[i], ns).d2); });
    return d;
}

// RMS of a difference against the RMS expected from the standard errors, and mean z^2 (z^2 uses
// its own variance se2z; points where it is 0 are left out of z^2 and counted in zeroNoise if
// their difference is not 0)
struct ErrStats {
    long long n = 0, nZ = 0, zeroNoise = 0;
    double sumD2 = 0, sumSe2 = 0, sumZ2 = 0, sumZ4 = 0;
    void add(double diff, double se2, double se2z)
    {
        n++;
        sumD2 += diff*diff;
        sumSe2 += se2;
        if (se2z > 0) {
            nZ++;
            const double z2 = diff*diff/se2z;
            sumZ2 += z2;
            sumZ4 += z2*z2;
        } else if (diff != 0) {
            zeroNoise++;
        }
    }
    void add(double diff, double se2) { add(diff, se2, se2); }
    double rms() const { return n ? std::sqrt(sumD2/double(n)) : 0.0; }
    double expected() const { return n ? std::sqrt(sumSe2/double(n)) : 0.0; }
    double z2() const { return nZ ? sumZ2/double(nZ) : 0.0; }
    // standard error of the mean z^2, from the spread of the z^2 values (no Gaussian assumption)
    double z2se() const
    {
        if (nZ < 2) return 0.0;
        const double m = z2(), var = (sumZ4/double(nZ) - m*m)*double(nZ)/double(nZ - 1);
        return std::sqrt(std::max(0.0, var)/double(nZ));
    }
};

// verdict words
const char* yesNo(bool b) { return b ? "yes" : "NO"; }

EstimateD toDouble(const Estimate& e)
{
    EstimateD d;
    d.H = e.H;
    d.Hse = e.Hse;
    d.walkVar = e.walkVar;
    d.grad = e.grad.cast<double>();
    d.gradSe = e.gradSe.cast<double>();
    d.hasGrad = e.hasGrad;
    d.n = e.nCounted;
    return d;
}

// V0: z^2 of a - b for H at the points without and with a gradient estimate, and for the grad components
struct PairZ {
    ErrStats hPlain, hGrad, grad;
};
PairZ comparePair(const std::vector<EstimateD>& a, const std::vector<EstimateD>& b, const std::vector<char>& gradPt)
{
    PairZ p;
    for (size_t k = 0; k < a.size(); k++) {
        (gradPt[k] ? p.hGrad : p.hPlain).add(a[k].H - b[k].H, a[k].Hse*a[k].Hse + b[k].Hse*b[k].Hse);
        if (a[k].hasGrad && b[k].hasGrad)
            for (int c = 0; c < 3; c++)
                p.grad.add(a[k].grad[c] - b[k].grad[c], a[k].gradSe[c]*a[k].gradSe[c] + b[k].gradSe[c]*b[k].gradSe[c]);
    }
    return p;
}

// V0 criterion: x (own - Zombie) inside [min, max] of the calibration pairs a (own - own) and c
// (Zombie - Zombie), widened by 3 combined standard errors
bool zBetween(const ErrStats& x, const ErrStats& a, const ErrStats& c)
{
    const ErrStats& lo = a.z2() <= c.z2() ? a : c;
    const ErrStats& hi = a.z2() <= c.z2() ? c : a;
    return x.z2() >= lo.z2() - 3.0*std::hypot(x.z2se(), lo.z2se()) && x.z2() <= hi.z2() + 3.0*std::hypot(x.z2se(), hi.z2se());
}

// R1: replays every cached walk from its seed under the cache's geometry and compares it bit for bit
long long replayMismatches(const WalkCache& A, const StripCut& cut)
{
    WalkRunner run(A, cut, A.nStrips());
    std::atomic<long long> bad{0};
    tbb::parallel_for(tbb::blocked_range<size_t>(0, A.walks.size()), [&](const tbb::blocked_range<size_t>& r) {
        long long b = 0;
        for (size_t wi = r.begin(); wi < r.end(); wi++) {
            WalkEntry e;
            Vec3 w;
            run.replay(wi, e, w);
            const Vec3* wc = A.gradWeight(wi);
            if (!sameWalk(e, wc ? &w : nullptr, A.walks[wi], wc, false)) b++;
        }
        bad += b;
    });
    return bad.load();
}

// Closest point of triangle abc to p in 3D (Ericson, Real-Time Collision Detection, 5.1.5). R4 uses
// it as a ball test that shares no code with StripCut's in-plane (u, w) distance.
Vec3 closestOnTriangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c)
{
    const Vec3 ab = b - a, ac = c - a, ap = p - a;
    const double d1 = ab.dot(ap), d2 = ac.dot(ap);
    if (d1 <= 0 && d2 <= 0) return a;
    const Vec3 bp = p - b;
    const double d3 = ab.dot(bp), d4 = ac.dot(bp);
    if (d3 >= 0 && d4 <= d3) return b;
    const double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) return a + (d1/(d1 - d3))*ab;
    const Vec3 cp = p - c;
    const double d5 = ab.dot(cp), d6 = ac.dot(cp);
    if (d6 >= 0 && d5 <= d6) return c;
    const double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) return a + (d2/(d2 - d6))*ac;
    const double va = d3*d6 - d5*d4;
    if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) return b + ((d4 - d3)/((d4 - d3) + (d5 - d6)))*(c - b);
    const double den = 1.0/(va + vb + vc);
    return a + (vb*den)*ab + (vc*den)*ac;
}

// R4, independent of the update's own ball test (stripDist2 against the step's d2 key): every walk
// is replayed under the old geometry (strips 0..m-1), and each step's ball(c_k, d_k) is tested
// against the 3D triangles of strip m with closestOnTriangle. Per walk: the first step whose ball
// the strip surely enters (3D distance < d_k - tol) and the first it may enter (< d_k + tol); -1: none.
struct BallHit {
    int sure = -1, maybe = -1;
};
std::vector<BallHit> independentBallHits(const WalkCache& A, const StripCut& cut, int m, double tol)
{
    std::vector<Vec3> V;
    std::vector<std::array<int, 3>> F;
    std::vector<double> triStrip;
    cut.mesh(m + 1, V, F, triStrip);
    std::vector<std::array<Vec3, 3>> T;
    for (size_t f = 0; f < F.size(); f++)
        if (int(triStrip[f]) == m) T.push_back({V[F[f][0]], V[F[f][1]], V[F[f][2]]});
    WalkRunner run(A, cut, m);
    std::vector<BallHit> hit(A.walks.size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, A.walks.size()), [&](const tbb::blocked_range<size_t>& r) {
        BallHit h;
        std::function<void(const WalkState&, const StripCut::Closest&)> obs = [&](const WalkState& st,
                                                                                  const StripCut::Closest& q) {
            if (h.sure >= 0) return;
            double D = std::numeric_limits<double>::infinity();
            for (const auto& t : T) D = std::min(D, (st.x - closestOnTriangle(st.x, t[0], t[1], t[2])).norm());
            double d = std::sqrt(q.d2);
            if (A.wp.dirichletFaces) d = std::min(d, A.tet.distToFaces(st.x));   // experiment: faces limit the ball
            if (h.maybe < 0 && D < d + tol) h.maybe = st.steps;
            if (D < d - tol) h.sure = st.steps;
        };
        for (size_t wi = r.begin(); wi < r.end(); wi++) {
            h = BallHit{};
            run.runTo(wi, A.walks[wi].rngOffset, -1, &obs);
            hit[wi] = h;
        }
    });
    return hit;
}

// log-log slope of the sample std of H over repeats against N, and its standard error from a
// bootstrap over the repeats (fixed seed; resamples with a zero std at some N are skipped).
// H[ni][r] = estimate of repeat r at Ns[ni]. NaN if a std is 0 (no noise to fit).
void ladderSlope(const std::vector<int>& Ns, const std::vector<std::vector<double>>& H, uint64_t seed, double& slope,
                 double& se)
{
    std::vector<double> nv(Ns.begin(), Ns.end()), sd;
    for (const auto& h : H) sd.push_back(sampleStd(h));
    slope = logLogSlope(nv, sd);
    se = std::nan("");
    if (std::isnan(slope)) return;
    std::mt19937_64 rng(seed);
    std::vector<double> bs;
    for (int b = 0; b < 1000; b++) {
        std::vector<double> sdB;
        for (const auto& h : H) {
            std::uniform_int_distribution<size_t> pick(0, h.size() - 1);
            std::vector<double> hb(h.size());
            for (double& x : hb) x = h[pick(rng)];
            sdB.push_back(sampleStd(hb));
        }
        double sl = logLogSlope(nv, sdB);
        if (!std::isnan(sl)) bs.push_back(sl);
    }
    if (bs.size() > 1) se = sampleStd(bs);
}

// R2: walk by walk (bitwise) and point by point
struct CacheDiff {
    long long walks = 0, mismatched = 0;
    double maxdH = 0, maxdG = 0;
};
CacheDiff compareCaches(const WalkCache& a, const WalkCache& b)
{
    CacheDiff d;
    d.walks = (long long)a.walks.size();
    std::atomic<long long> bad{0};
    tbb::parallel_for(tbb::blocked_range<size_t>(0, a.walks.size()), [&](const tbb::blocked_range<size_t>& r) {
        long long n = 0;
        for (size_t wi = r.begin(); wi < r.end(); wi++)
            if (!sameWalk(a.walks[wi], a.gradWeight(wi), b.walks[wi], b.gradWeight(wi), true)) n++;
        bad += n;
    });
    d.mismatched = bad.load();
    for (size_t i = 0; i < a.nPoints(); i++) {
        EstimateD ea = a.estimate(i), eb = b.estimate(i);
        d.maxdH = std::max(d.maxdH, std::fabs(ea.H - eb.H));
        if (ea.hasGrad && eb.hasGrad) d.maxdG = std::max(d.maxdG, (ea.grad - eb.grad).cwiseAbs().maxCoeff());
    }
    return d;
}

// R3: statistics of the invalidation step, and a direct check of every resume state: a from-scratch
// run of the walk under the NEW geometry reaches step k in exactly the replayed state, its steps
// 0..k-1 have the same distance and closest triangle under the old and the new geometry, at step k
// the new strip is strictly closer, and the RNG is exactly rngOffset + 2k draws past the seed.
struct ResumeCheck {
    long long n = 0, k0 = 0, badState = 0, badPrefix = 0, badFirst = 0, badRng = 0;
    std::vector<double> ks, kept;
};
void checkResumes(const WalkCache& A, const StripCut& cut, int m, const std::vector<ResumeRecord>& recs, ResumeCheck& rc)
{
    WalkRunner runNew(A, cut, m + 1);
    const int triBegin = cut.strip(m).triBegin;
    std::atomic<long long> badState{0}, badPrefix{0}, badFirst{0}, badRng{0};
    tbb::parallel_for(size_t(0), recs.size(), [&](size_t h) {
        const ResumeRecord& r = recs[h];
        long long prefix = 0;
        std::function<void(const WalkState&, const StripCut::Closest&)> obs =
            [&](const WalkState& st, const StripCut::Closest& qNew) {
                if (st.steps >= r.k) return;
                StripCut::Closest qOld = cut.closest(st.x, m);
                bool changed = !sameBits(&qOld.d2, &qNew.d2, sizeof(double)) || qOld.tri != qNew.tri;
                // dirichletFaces experiment: a nearer face fixes the step, so a cut change behind it is harmless
                if (changed && A.wp.dirichletFaces) {
                    const double dF = A.tet.distToFaces(st.x);
                    if (std::min(qOld.d2, qNew.d2) > dF*dF) changed = false;
                }
                if (changed) prefix++;
            };
        WalkState s = runNew.runTo(r.walk, r.rngOffset, r.k, &obs);
        badPrefix += prefix;
        if (!sameState(s, r.state)) badState++;
        StripCut::Closest qNew = cut.closest(r.state.x, m + 1), qOld = cut.closest(r.state.x, m);
        if (!(qNew.d2 < qOld.d2) || qNew.tri < triBegin) badFirst++;
        pcg32 r0 = runNew.initialState(r.walk, 0).rng;
        if (r.state.rng - r0 != int64_t(r.rngOffset) + 2*int64_t(r.k)) badRng++;
    });
    rc.badState += badState.load();
    rc.badPrefix += badPrefix.load();
    rc.badFirst += badFirst.load();
    rc.badRng += badRng.load();
    for (const ResumeRecord& r : recs) {
        rc.n++;
        rc.k0 += r.k == 0;
        rc.ks.push_back(double(r.k));
        rc.kept.push_back(r.oldSteps > 0 ? double(r.k)/double(r.oldSteps) : 0.0);
    }
}

// closest() with strip skipping against the brute force over all triangles in index order
long long closestPruningMismatches(const Tet& tet, const StripCut& cut, int nPts, uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::vector<Vec3> pts;
    while (int(pts.size()) < nPts) {   // uniform in the tet, plus points pushed onto the cut plane
        double a = U(rng), b = U(rng), c = U(rng);
        if (a + b + c > 1) continue;
        Vec3 x = tet.v[0] + a*(tet.v[1] - tet.v[0]) + b*(tet.v[2] - tet.v[0]) + c*(tet.v[3] - tet.v[0]);
        if (pts.size()%4 == 0) x -= cut.frame().side(x)*cut.frame().n;
        pts.push_back(x);
    }
    std::atomic<long long> bad{0};
    const int ns = cut.nStrips();
    tbb::parallel_for(size_t(0), pts.size(), [&](size_t i) {
        StripCut::Closest a = cut.closest(pts[i], ns), b = cut.closestBruteForce(pts[i], ns);
        if (!sameBits(&a.d2, &b.d2, sizeof(double)) || a.tri != b.tri || !sameBits(a.q.data(), b.q.data(), 2*sizeof(double)))
            bad++;
    });
    return bad.load();
}

QueryPoints makePoints(const std::vector<Vec3>& pts, const std::vector<char>* grad, uint64_t seedBase)
{
    QueryPoints q;
    q.add(pts, grad, seedBase);
    return q;
}

// pulls a probe towards the centroid until it is at least 0.01 L inside (as the first gate)
int pullInside(const Tet& tet, Vec3& x)
{
    int moved = 0;
    while (tet.distToFaces(x) < 0.01*tet.L && moved < 100) {
        x += 0.1*(tet.centroid - x);
        moved++;
    }
    return moved;
}

const char* gShortName(GType g) { return g == GType::One ? "(a) g=1" : "(b) g=min(1,sqrt(r/L))"; }

} // namespace

// =====================================================================================
bool runReuseValidation(const GeometryConfig& cfg, const Settings& s, const std::string& resultsPath, bool quick)
{
    auto tStart = std::chrono::steady_clock::now();
    Tet tet = referenceTet();
    const double L = tet.L, eps = s.epsRel*L, nearDist = 0.1*L;
    const GType gTypes[2] = {GType::One, GType::SqrtFront};
    const int N = quick ? 32 : s.nWalks;
    const int K = quick ? 6 : s.k;
    Report out;
    bool exactOk = true;
    std::vector<std::string> verdicts;   // one line per statistical check, repeated at the end

    std::time_t now = std::time(nullptr);
    std::tm tmNow{};
#ifdef _WIN32
    localtime_s(&tmNow, &now);
#else
    localtime_r(&now, &tmNow);
#endif
    char date[64];
    std::strftime(date, sizeof date, "%Y-%m-%d %H:%M:%S", &tmNow);
    StripCut frame(tet, cfg);
    const Cut& F = frame.frame();
    out("Incremental walk reuse under an advancing cut (own WoSt loop, walk_cache.h): validation%s\n",
        quick ? " -- QUICK reduced-size run (R1-R4 only, not written to results.txt)" : "");
    out("date: %s; hardware threads: %u\n", date, std::thread::hardware_concurrency());
    out("tet: v0 = %s, v1 = %s, v2 = %s, v3 = %s; L = %.6f\n", fmt3(tet.v[0]).c_str(), fmt3(tet.v[1]).c_str(),
        fmt3(tet.v[2]).c_str(), fmt3(tet.v[3]).c_str(), L);
    out("cut plane n = %s, d = %.6f; front direction u = %s; s in [%.4f, %.4f]; full cut area %.5f\n",
        fmt3(F.n).c_str(), F.d, fmt3(F.u).c_str(), F.sMin, F.sMax, frame.fullArea());
    if (frame.addStrip(frame.sOf(0.5)) < 0) {
        out("ERROR: the plane misses the tet; nothing to validate.\n");
        if (!quick) writeResultsSection(resultsPath, kResultsSectionReuse, out.text);
        return false;
    }
    out("walks: own double-precision WoSt loop; termination shell %.3g (%.0e L) on exact distances; max walk length %d;\n"
        "gradient points (interior lattice points) take a first ball of %.2f min(dist to cut, dist to faces)\n",
        s.walk.absorbShellRel*L, s.walk.absorbShellRel, s.walk.maxWalkLength, kGradFirstBallShrink);
    out("cut = list of strips; each step appends one strip [s_old, s_new]; geometry queries brute force in double\n");
    out("epsilon = %.5g (%.4g L): lattice points within epsilon of the CURRENT cut are masked out of every statistic\n"
        "(the cache holds all of them); 'near the cut' = within %.4g (0.1 L) of the current cut and not masked\n",
        eps, s.epsRel, nearDist);
    out("RNG: base seed %llu; a point set with tag T has point seeds mix(mix(base, T), i); walk j of point i uses\n"
        "pcg32(mix(seed_i, j), j), exactly 2 draws per step, none at termination\n", (unsigned long long)s.seed);

    // all lattice points (no epsilon exclusion in the cache)
    Lattice lat = buildLattice(tet, Cut{}, K, 0.0);
    int nInterior = 0;
    for (char c : lat.interior) nInterior += c;
    out("lattice k = %d: %zu points (%d interior, with grad H); N = %d walks per point unless stated\n", K, lat.pts.size(),
        nInterior, N);

    // ================================================================ V0
    if (!quick) {
        const double t0 = 0.6;
        const int NV = 1024;
        Cut cutZ = buildCut(tet, cfg, t0);
        StripCut cut0(tet, cfg);
        cut0.addStrip(cut0.sOf(t0));
        std::vector<Vec3> keepPts;
        std::vector<char> keepGrad;
        int nKeepGrad = 0;
        for (size_t i = 0; i < lat.pts.size(); i++)
            if (cutZ.distToCut(lat.pts[i]) >= eps) {   // the first gate's exclusion rule
                keepPts.push_back(lat.pts[i]);
                keepGrad.push_back(lat.interior[i]);
                nKeepGrad += lat.interior[i];
            }
        out("\n=== V0. Own WoSt loop vs Zombie (t = %.2f, N = %d, the %zu non-excluded lattice points: %zu without, %d with"
            " grad H) ===\n", t0, NV, keepPts.size(), keepPts.size() - nKeepGrad, nKeepGrad);
        out("own: strip list with one strip; Zombie: CutSolver of the first gate on buildCut(t). Here the own loop stops at\n"
            "Zombie's stop distance (shell + faceOffsetMargin) instead of its usual shell, so both stop at the same distance.\n"
            "Three pairs of independent solves (seed tags in brackets): own - Zombie, and the calibration pairs own - own and\n"
            "Zombie - Zombie. Columns: mean z^2 +- its s.e., z = difference / sqrt(s.e._1^2 + s.e._2^2).\n");
        bool v0ok = true;
        for (int gi = 0; gi < 2; gi++) {
            CutSolver zs(tet, cutZ, gTypes[gi], s.walk);
            WalkParams wpV = s.walk;
            wpV.absorbShellRel = zs.stopDistance()/L;
            SolveStats so, sz;
            auto own = [&](int tag, SolveStats* st) {
                return solveSumsOnly(tet, cut0, makePoints(keepPts, &keepGrad, mixSeed(s.seed, tag)), NV, gTypes[gi], wpV, st);
            };
            auto zombie = [&](int tag, SolveStats* st) {
                std::vector<Estimate> z;
                SolveStats ss = zs.solve(keepPts, &keepGrad, NV, mixSeed(s.seed, tag), z);
                if (st) *st = ss;
                std::vector<EstimateD> d(z.size());
                for (size_t k = 0; k < z.size(); k++) d[k] = toDouble(z[k]);
                return d;
            };
            std::vector<EstimateD> oA = own(500 + gi, &so), oB = own(520 + gi, nullptr);
            std::vector<EstimateD> zA = zombie(510 + gi, &sz), zB = zombie(530 + gi, nullptr);
            const PairZ oz = comparePair(oA, zA, keepGrad), oo = comparePair(oA, oB, keepGrad),
                        zz = comparePair(zA, zB, keepGrad);
            out("-- %s; stop distance %.3g (own's usual shell %.3g) --\n", gShortName(gTypes[gi]), zs.stopDistance(),
                s.walk.absorbShellRel*L);
            out("   pair                     H, no-grad points    H, grad points    grad H components |  RMS dH (expected)\n");
            auto row = [&](const char* name, const PairZ& p) {
                out("   %-24s %5.3f +- %5.3f      %5.3f +- %5.3f    %5.3f +- %5.3f     |  %.5f (%.5f)\n", name, p.hPlain.z2(),
                    p.hPlain.z2se(), p.hGrad.z2(), p.hGrad.z2se(), p.grad.z2(), p.grad.z2se(),
                    std::sqrt((p.hPlain.sumD2 + p.hGrad.sumD2)/double(p.hPlain.n + p.hGrad.n)),
                    std::sqrt((p.hPlain.sumSe2 + p.hGrad.sumSe2)/double(p.hPlain.n + p.hGrad.n)));
            };
            char nm[3][64];
            std::snprintf(nm[0], 64, "own - Zombie (%d, %d)", 500 + gi, 510 + gi);
            std::snprintf(nm[1], 64, "own - own (%d, %d)", 500 + gi, 520 + gi);
            std::snprintf(nm[2], 64, "Zombie - Zombie (%d, %d)", 510 + gi, 530 + gi);
            row(nm[0], oz);
            row(nm[1], oo);
            row(nm[2], zz);
            const bool ok = zBetween(oz.hPlain, oo.hPlain, zz.hPlain) && zBetween(oz.hGrad, oo.hGrad, zz.hGrad) &&
                            zBetween(oz.grad, oo.grad, zz.grad);
            v0ok = v0ok && ok;
            out("   verdict: own - Zombie lies between own - own and Zombie - Zombie (+- 3 s.e.) in all three columns: %s\n",
                yesNo(ok));
            out("   own walks:    %s\n", so.summary().c_str());
            out("   Zombie walks: %s\n", sz.summary().c_str());
        }
        out("Criterion: the z^2 of own - Zombie weighs the s.e. calibration of each solver by its variance, so it lies\n"
            "(approximately) between the own - own and Zombie - Zombie values; both solvers estimate the same H iff it does.\n"
            "Zombie's s.e. is off at the gradient points: its antithetic pairs replay the same random numbers after the\n"
            "mirrored first step (walk_on_stars.h), so the pair's two H values are positively correlated and its two\n"
            "gradient terms negatively correlated, while its s.e. treats them as independent: H z^2 above 1, grad z^2\n"
            "below 1 there. The mean walk lengths differ because Zombie shrinks every star radius by 0.99 and, with\n"
            "solveDoubleSided, keeps the tet's convex edges as silhouette candidates, so its stars shrink near edges; the\n"
            "own loop takes the full radius (a convex domain has no silhouette seen from inside).\n");
        verdicts.push_back(std::string("V0 own loop = Zombie within noise: ") + yesNo(v0ok));
    }

    // ================================================================ sweep: R1-R4 and check 1
    const double tFirst = 0.2, tLast = quick ? 0.8 : 1.0;
    const int nSteps = quick ? 6 : 16;
    out("\n=== Sweep t = %.2f -> %.2f in %d steps (lattice, N = %d): R1, R2, R3, R4%s ===\n", tFirst, tLast, nSteps, N,
        quick ? "" : " and check 1");
    // Both distances are exact up to ~1e-15 (coordinates ~1); 1e-12 L leaves a wide margin. Real ties
    // exist: a step centre on the old front's plane X·u = s_old (e.g. the lattice points with x = y when
    // the front is at t = 0.5) has its closest point on the shared front edge, at the same distance from
    // the new strip; just behind that plane the strip is only (s_old - s)^2 / 2d farther.
    const double tieTol = 1e-12*L;
    out("chain A: incremental (resume); chain B: biased restart (paper Fig. 5), started from a copy of A's initial cache;\n"
        "F: same-seed full solve on the same strip list after each step (the R2 oracle, and the timed full re-solve).\n"
        "t_update here includes collecting A's resume records for R3; check 3 and --bench-reuse time it without them.\n"
        "Both g types use the same chain seeds (tag 600), so their paths and walk sets are the same; only values differ.\n"
        "R4 is independent of the update's ball test: before each update every walk of A is replayed under the old cut and\n"
        "each step's ball(c_k, d_k) is tested against the new strip's 3D triangles (Ericson's closest point on a\n"
        "triangle). 'disagree' = walks where the update's choice (re-walked or not, and from which step k) differs from\n"
        "that test; steps whose 3D distance is within %.0e L of d_k count as ties and may go either way (a step centre\n"
        "on the old front's plane is an exact tie: its closest point lies on the front edge the strip shares).\n", tieTol/L);
    if (!quick)
        out("check 1: ref = fresh full solve with 16 N = %d walks and independent seeds (tag 1000 + 100 g + step), fresh =\n"
            "same-N solve with independent seeds (tag 2000 + 100 g + step); both keep sums only. 'expected' = RMS of the\n"
            "noise of an N-walk estimate minus the reference, sqrt(mean(v/N + v/(16 N))), v = the reference's per-walk\n"
            "variance; z^2 = (H - H_ref)^2 / (v/N + v/(16 N)) with v pooled from both estimates (16:1), so an unbiased\n"
            "estimate has mean z^2 near 1\n", 16*N);
    long long r1Bad[2] = {0, 0}, r1BadB[2] = {0, 0}, r1Walks[2] = {0, 0}, r2Bad[2] = {0, 0}, r2Walks[2] = {0, 0};
    long long r4Missed[2] = {0, 0}, r4Disagree[2] = {0, 0}, r4Ties[2] = {0, 0}, pruneBad = 0;
    double r2MaxH[2] = {0, 0}, r2MaxG[2] = {0, 0};
    ResumeCheck rc[2];
    std::vector<uint64_t> walkSetHash[2];   // per step: hash of the tested count and the re-walked (walk, k) list
    long long zeroNoise[2] = {0, 0};
    struct C1 {
        double msInc = 0, msFresh = 0, msBiased = 0, msExp = 0;
        double z2Inc = 0, z2Fresh = 0, z2Biased = 0, z2SeInc = 0, z2SeFresh = 0, z2SeBiased = 0;
    };
    std::vector<C1> c1All[2], c1Near[2];
    for (int gi = 0; gi < 2; gi++) {
        const GType g = gTypes[gi];
        StripCut cut(tet, cfg);
        cut.addStrip(cut.sOf(tFirst));
        WalkCache A;
        A.configure(tet, makePoints(lat.pts, &lat.interior, mixSeed(s.seed, 600)), N, g, s.walk);
        SolveStats s0 = A.solveFull(cut);
        WalkCache B = A;
        WalkCache Fc;
        out("\n-- %s (chain seed tag 600): initial full solve at t = %.2f: %s --\n", gShortName(g), tFirst,
            s0.summary().c_str());
        out("step    t   strip/  tested  re-walk  value-  t_update t_full  |  R1 A  R1 B |  R2 A: walks  max|dH|   max|dgrad| |"
            "  R4: missed disagree bbox FP | B vs F: walks  RMS dH\n"
            "              full                     upd      (s)     (s)   |  mism  mism |  mismatched                        |"
            "      walks   walks    rate   |  mismatched    (masked)\n");
        std::vector<std::string> c1rows;
        long long rwEsc = 0, rwCap = 0, stReplayed = 0, stNew = 0, stFull = 0;
        for (int step = 1; step <= nSteps; step++) {
            const double t = tFirst + (tLast - tFirst)*double(step)/double(nSteps);
            const int m = cut.addStrip(cut.sOf(t));
            if (m < 0) continue;
            // R4: the independent ball test on A's walks before the update (old geometry = strips 0..m-1)
            std::vector<BallHit> ind = independentBallHits(A, cut, m, tieTol);
            std::vector<ResumeRecord> recs;
            ReuseStats ra = A.advance(cut, UpdateMode::Incremental, &recs);
            B.advance(cut, UpdateMode::BiasedRestart, nullptr);
            rwEsc += ra.escaped;
            rwCap += ra.capped;
            stReplayed += ra.stepsReplayed;
            stNew += ra.stepsNew;
            std::vector<int> kOf(A.walks.size(), -1);
            uint64_t hsh = mixSeed(uint64_t(ra.nTested), uint64_t(step));
            for (const ResumeRecord& r : recs) {
                kOf[r.walk] = r.k;
                hsh = mixSeed(hsh ^ uint64_t(r.walk), uint64_t(r.k));
            }
            walkSetHash[gi].push_back(hsh);
            long long missed = 0, disagree = 0, ties = 0, bboxHits = 0, fp = 0;
            for (size_t wi = 0; wi < ind.size(); wi++) {
                const BallHit& h = ind[wi];
                const bool box = A.action[wi] & 1, rewalked = A.action[wi] & 2;
                const int k = kOf[wi];
                if (h.sure >= 0 && !box) missed++;
                if (rewalked != (k >= 0)) disagree++;
                else if (k < 0 ? h.sure >= 0 : (h.maybe < 0 || k < h.maybe || (h.sure >= 0 && k > h.sure))) disagree++;
                if (h.maybe >= 0 && h.maybe != h.sure) ties++;
                if (box) {
                    bboxHits++;
                    if (h.maybe < 0) fp++;
                }
            }
            r4Missed[gi] += missed;
            r4Disagree[gi] += disagree;
            r4Ties[gi] += ties;
            // R1 on both chains
            long long badA = replayMismatches(A, cut), badB = replayMismatches(B, cut);
            r1Bad[gi] += badA;
            r1BadB[gi] += badB;
            r1Walks[gi] += (long long)A.walks.size();
            // R2: same-seed fresh solve
            Fc.configure(tet, A.pts, N, g, s.walk);
            SolveStats sf = Fc.solveFull(cut);
            stFull += sf.totalSteps;
            CacheDiff dA = compareCaches(A, Fc), dB = compareCaches(B, Fc);
            r2Bad[gi] += dA.mismatched;
            r2Walks[gi] += dA.walks;
            r2MaxH[gi] = std::max(r2MaxH[gi], dA.maxdH);
            r2MaxG[gi] = std::max(r2MaxG[gi], dA.maxdG);
            // R3
            checkResumes(A, cut, m, recs, rc[gi]);
            // masks at the current cut
            std::vector<double> dc = cutDistances(cut, lat.pts);
            std::vector<EstimateD> eA = A.estimates(), eB = B.estimates(), eF = Fc.estimates();
            ErrStats bVsF;
            for (size_t i = 0; i < lat.pts.size(); i++)
                if (dc[i] >= eps) bVsF.add(eB[i].H - eF[i].H, 0.0);
            const double nW = double(ra.nWalks);
            out("%3d  %.3f  %.4f  %6.2f%%  %6.2f%%  %6.2f%%  %7.4f  %6.3f  | %5lld %5lld | %10lld  %.2e  %.2e | %8lld %8lld   %5.1f%% | %10lld   %.5f\n",
                step, t, cut.strip(m).area/cut.fullArea(), 100.0*ra.nTested/nW, 100.0*ra.nRewalked/nW,
                100.0*ra.nValueUpdated/nW, ra.tTotal, sf.seconds, badA, badB, dA.mismatched, dA.maxdH, dA.maxdG, missed,
                disagree, bboxHits ? 100.0*double(fp)/double(bboxHits) : 0.0, dB.mismatched, bVsF.rms());
            if (step == nSteps) pruneBad += closestPruningMismatches(tet, cut, quick ? 20000 : 200000, mixSeed(s.seed, 650 + gi));

            if (!quick) {   // check 1
                SolveStats sr, sfr;
                std::vector<EstimateD> ref = solveSumsOnly(tet, cut, makePoints(lat.pts, &lat.interior, mixSeed(s.seed, 1000 + 100*gi + step)),
                                                           16*N, g, s.walk, &sr);
                std::vector<EstimateD> fresh = solveSumsOnly(tet, cut, makePoints(lat.pts, &lat.interior, mixSeed(s.seed, 2000 + 100*gi + step)),
                                                             N, g, s.walk, &sfr);
                // Expected noise of (N-walk estimate - reference): v/N + v/(16 N), v = per-walk variance.
                // 'expected' takes v from the reference's 16 N walks. z^2 pools the two estimates' v
                // (16:1), because where almost every walk returns the same +-1 either one alone is often
                // exactly 0 (the N-walk s.e. most of the time, the reference's now and then).
                ErrStats inc[2], frs[2], bia[2];   // [all, near]
                auto pooled = [&](const EstimateD& e, const EstimateD& r) {
                    const double v = (16.0*r.walkVar + e.walkVar)/17.0;
                    return v/double(N) + v/double(16*N);
                };
                for (size_t i = 0; i < lat.pts.size(); i++) {
                    if (dc[i] < eps) continue;
                    const double se2 = ref[i].walkVar/double(N) + ref[i].Hse*ref[i].Hse;
                    for (int sub = 0; sub < 2; sub++) {
                        if (sub == 1 && dc[i] >= nearDist) continue;
                        inc[sub].add(eA[i].H - ref[i].H, se2, pooled(eA[i], ref[i]));
                        frs[sub].add(fresh[i].H - ref[i].H, se2, pooled(fresh[i], ref[i]));
                        bia[sub].add(eB[i].H - ref[i].H, se2, pooled(eB[i], ref[i]));
                    }
                }
                for (int sub = 0; sub < 2; sub++) zeroNoise[gi] += inc[sub].zeroNoise + frs[sub].zeroNoise + bia[sub].zeroNoise;
                for (int sub = 0; sub < 2; sub++) {
                    C1 c;
                    c.msInc = inc[sub].rms()*inc[sub].rms();
                    c.msFresh = frs[sub].rms()*frs[sub].rms();
                    c.msBiased = bia[sub].rms()*bia[sub].rms();
                    c.msExp = inc[sub].expected()*inc[sub].expected();
                    c.z2Inc = inc[sub].z2();
                    c.z2Fresh = frs[sub].z2();
                    c.z2Biased = bia[sub].z2();
                    c.z2SeInc = inc[sub].z2se();
                    c.z2SeFresh = frs[sub].z2se();
                    c.z2SeBiased = bia[sub].z2se();
                    (sub ? c1Near : c1All)[gi].push_back(c);
                }
                char buf[512];
                std::snprintf(buf, sizeof buf,
                              "%3d  %.3f | %4lld  %.5f  %.5f  %.5f  %.5f | %4lld  %.5f  %.5f  %.5f  %.5f | %5.2f %5.2f %5.2f | %.2f s\n",
                              step, t, inc[0].n, inc[0].rms(), frs[0].rms(), bia[0].rms(), inc[0].expected(), inc[1].n,
                              inc[1].rms(), frs[1].rms(), bia[1].rms(), inc[1].expected(), inc[1].z2(), frs[1].z2(),
                              bia[1].z2(), sr.seconds);
                c1rows.push_back(buf);
            }
        }
        out("   after the sweep: %d strips, %d triangles; R4: %lld walks with a tie step (3D distance within %.0e L of d_k)\n",
            cut.nStrips(), cut.nTris(), r4Ties[gi], tieTol/L);
        out("   chain A's re-walked paths: %lld escaped, %lld capped; the updates ran %lld replayed (old cut) + %lld new steps\n"
            "   = %.1f%% of the %lld steps of the same-seed full solves F\n", rwEsc, rwCap, stReplayed, stNew,
            stFull ? 100.0*double(stReplayed + stNew)/double(stFull) : 0.0, stFull);
        if (!quick) {
            out("\n   check 1, %s: RMS of H - H_ref over the non-masked lattice points (all) and those near the cut (near)\n",
                gShortName(g));
            out("step    t   | all: n   inc      fresh    biased   expected | near: n inc    fresh    biased   expected |"
                " near z^2: inc fresh biased | t_ref\n");
            for (const std::string& r : c1rows) out("%s", r.c_str());
            if (zeroNoise[gi])
                out("   (%lld point comparisons had a nonzero difference with zero pooled variance; left out of z^2)\n",
                    zeroNoise[gi]);
        }
    }

    // ---- exact verdicts
    out("\n=== Exact checks ===\n");
    for (int gi = 0; gi < 2; gi++) {
        out("-- %s --\n", gShortName(gTypes[gi]));
        bool r1 = r1Bad[gi] == 0 && r1BadB[gi] == 0;
        out("R1 replay exactness: %lld walks replayed per chain over %d steps; mismatches: incremental %lld, biased %lld -> %s\n",
            r1Walks[gi], nSteps, r1Bad[gi], r1BadB[gi], r1 ? "PASS" : "FAIL");
        bool r2 = r2Bad[gi] == 0;
        out("R2 resume = same-seed fresh solve: %lld walks compared, %lld mismatched; max |dH| = %.2e, max |d grad H| = %.2e -> %s\n",
            r2Walks[gi], r2Bad[gi], r2MaxH[gi], r2MaxG[gi], r2 ? "PASS" : "FAIL");
        bool r4 = r4Missed[gi] == 0 && r4Disagree[gi] == 0;
        out("R4 prefilter and ball test vs the independent 3D ball test: walks rejected by the bbox test although the strip\n"
            "   enters one of their balls: %lld; walks where the update's choice (re-walked or not, invalidation step k)\n"
            "   differs beyond ties: %lld -> %s\n", r4Missed[gi], r4Disagree[gi], r4 ? "PASS" : "FAIL");
        const ResumeCheck& c = rc[gi];
        bool r3 = c.badState == 0 && c.badPrefix == 0 && c.badFirst == 0 && c.badRng == 0;
        out("R3 resume state: %lld invalidated walks checked; resume state != from-scratch state at step k: %lld; prefix steps\n"
            "   with a different distance or closest triangle: %lld; step k not changed by the strip: %lld; RNG not at\n"
            "   rngOffset + 2k draws: %lld -> %s\n", c.n, c.badState, c.badPrefix, c.badFirst, c.badRng, r3 ? "PASS" : "FAIL");
        exactOk = exactOk && r1 && r2 && r4 && r3;
    }
    out("closest-point strip skipping vs brute force over all triangles (random points, final geometries): %lld mismatches -> %s\n",
        pruneBad, pruneBad == 0 ? "PASS" : "FAIL");
    exactOk = exactOk && pruneBad == 0;

    // ---- R3 statistics
    out("\n=== R3. Resume statistics of the invalidated walks (chain A, all sweep steps) ===\n");
    for (int gi = 0; gi < 2; gi++) {
        const ResumeCheck& c = rc[gi];
        out("%s: %lld re-walked walks; k = 0 (resume at the query point, same seed) %.2f%%; median k = %.0f; mean k = %.1f;\n"
            "   mean fraction of the walk's original steps kept (k / old step count) = %.3f\n", gShortName(gTypes[gi]), c.n,
            c.n ? 100.0*double(c.k0)/double(c.n) : 0.0, medianOf(c.ks), meanOf(c.ks), meanOf(c.kept));
    }
    const bool sameSets = walkSetHash[0] == walkSetHash[1] && rc[0].n == rc[1].n && rc[0].ks == rc[1].ks;
    out("walk sets of (a) and (b) identical at every step (same seeds; tested count and the re-walked (walk, k) list\n"
        "compared by hash): %s\n", yesNo(sameSets));

    if (quick) {
        out("\nquick run: %s, %.1f s\n", exactOk ? "all exact checks PASS" : "SOME EXACT CHECK FAILED", secondsSince(tStart));
        return exactOk;
    }

    // ---- check 1 summary
    out("\n=== 1. Unbiasedness summary (means over the %d steps; MSE = mean squared error against the 16 N reference) ===\n",
        nSteps);
    out("The steps of one chain are not independent samples: chain A is, by R2, one same-seed solve per step, so its\n"
        "errors at a point persist over the steps where its walks are reused.\n"
        "Verdict criterion on the mean z^2 (1 for an unbiased estimate with the right variance): within 1 +- 3 s.e.,\n"
        "s.e. = the mean over the steps of each step's s.e. of mean z^2 (from the spread of its z^2 values). That bounds\n"
        "the s.e. of a mean over steps however correlated the steps are. Biased chain: 'measurably worse' if above\n"
        "1 + 3 s.e. of its own z^2.\n");
    bool c1Unbiased = true, c1BiasedWorseNear = true;
    std::string c1NearVsAll;
    for (int gi = 0; gi < 2; gi++) {
        double excess[2] = {0, 0}, excessSe[2] = {0, 0};   // biased mean z^2 - 1 and its s.e., [all, near]
        for (int sub = 0; sub < 2; sub++) {
            const std::vector<C1>& v = (sub ? c1Near : c1All)[gi];
            std::vector<double> inc, fr, bi, ex, zi, zf, zb, si, sf, sb;
            for (const C1& c : v) {
                if (c.msExp <= 0) continue;   // t = 1, g = 1: H is exactly +-1 everywhere, no noise to compare
                inc.push_back(c.msInc);
                fr.push_back(c.msFresh);
                bi.push_back(c.msBiased);
                ex.push_back(c.msExp);
                zi.push_back(c.z2Inc);
                zf.push_back(c.z2Fresh);
                zb.push_back(c.z2Biased);
                si.push_back(c.z2SeInc);
                sf.push_back(c.z2SeFresh);
                sb.push_back(c.z2SeBiased);
            }
            out("%s, %s (%zu steps with noise): MSE inc %.3e, fresh %.3e, biased %.3e, expected %.3e;\n"
                "   MSE / expected: inc %.3f, fresh %.3f, biased %.3f; mean z^2: inc %.3f, fresh %.3f, biased %.3f\n",
                gShortName(gTypes[gi]), sub ? "near the cut" : "all points", inc.size(), meanOf(inc), meanOf(fr), meanOf(bi),
                meanOf(ex), meanOf(inc)/meanOf(ex), meanOf(fr)/meanOf(ex), meanOf(bi)/meanOf(ex), meanOf(zi), meanOf(zf),
                meanOf(zb));
            const bool okI = std::fabs(meanOf(zi) - 1.0) <= 3.0*meanOf(si), okF = std::fabs(meanOf(zf) - 1.0) <= 3.0*meanOf(sf);
            const bool worse = meanOf(zb) - 1.0 > 3.0*meanOf(sb);
            out("   verdict: incremental z^2 within 1 +- %.3f: %s; fresh within 1 +- %.3f: %s; biased above 1 + %.3f "
                "(measurably worse): %s\n", 3.0*meanOf(si), yesNo(okI), 3.0*meanOf(sf), yesNo(okF), 3.0*meanOf(sb), yesNo(worse));
            c1Unbiased = c1Unbiased && okI && okF;
            if (sub == 1) c1BiasedWorseNear = c1BiasedWorseNear && worse;
            excess[sub] = meanOf(zb) - 1.0;
            excessSe[sub] = meanOf(sb);
        }
        const double diff = excess[1] - excess[0], diffSe = std::hypot(excessSe[0], excessSe[1]);
        const char* rel = diff > 3.0*diffSe ? "larger" : diff < -3.0*diffSe ? "smaller" : "about the same (within 3 s.e.)";
        out("   %s: biased z^2 excess near the cut %.2f vs over all points %.2f (difference %+.2f +- %.2f): %s near the cut\n",
            gShortName(gTypes[gi]), excess[1], excess[0], diff, diffSe, rel);
        c1NearVsAll += std::string(gi ? ", " : "") + (gi ? "(b) " : "(a) ") + rel;
    }
    verdicts.push_back(std::string("1  incremental (and fresh) unbiased within noise: ") + yesNo(c1Unbiased) +
                       "; biased restart measurably worse near the cut: " + yesNo(c1BiasedWorseNear) +
                       "; its excess near the cut vs over all points: " + c1NearVsAll);

    // ================================================================ check 2
    {
        const double t2 = 0.6;
        Cut cutP = buildCut(tet, cfg, t2);
        Vec3 Fm = cutP.frontMidpoint();
        const std::vector<int> Ns = {16, 64, 256, 1024, 4096};
        struct Group {
            const char* name;
            std::vector<Vec3> probes;
            int R, tag, ctrlTag;
            bool zombieCtrl;   // control: Zombie with the first gate's seeds; else a fresh own solve
        };
        // The first gate's probes lie ahead of the front, where nearly every walk is re-walked from k = 0
        // at every step. The reuse probes lie behind the initial front (t = 0.2) close to the cut: most of
        // their walks end on the old cut and survive the 8 steps, so their ladder exercises reused (and,
        // for (b), value-updated) walks. With g = 1 their H is close to +-1 (only the few walks that pass
        // the front change sign), so they get 32 repeats to keep the sample stds away from 0.
        const Vec3 F01 = buildCut(tet, cfg, 0.1).frontMidpoint(), F02 = buildCut(tet, cfg, 0.2).frontMidpoint();
        std::vector<Group> groups = {
            {"the first gate's 4 probes", {Fm + 0.05*L*cutP.u, Fm + 0.05*L*cutP.n, Fm - 0.05*L*cutP.u + 0.03*L*cutP.n,
                                            Fm + 0.10*L*cutP.u - 0.05*L*cutP.n}, 8, 700, 100, true},
            {"4 reuse probes behind the initial front", {F01 - 0.02*L*F.n, F01 + 0.10*L*F.n, F02 - 0.05*L*F.n, F02 + 0.05*L*F.n},
             32, 760, 780, false}};
        out("\n=== 2. Convergence under reuse at t = %.2f after 8 incremental steps from t = 0.2 ===\n", t2);
        out("Per N: own seeds (tag T + 10 g + i_N), a full solve at t = 0.2, then 8 steps of 0.05. Control: for the first\n"
            "gate's probes the first gate itself (Zombie with the seeds of --validate, tags 100 + 10 g + i_N); for the reuse\n"
            "probes a fresh own solve on the same strips (tag 780 + 10 g + i_N). slope = least-squares slope of log std\n"
            "against log N, +- its s.e. from a bootstrap over the repeats (1000 resamples). Per probe: 're-walked' = mean\n"
            "fraction of its walks re-walked per step, 'kept' = fraction of its final walks never re-walked since t = 0.2.\n"
            "Verdict per group: the mean slope of its probes within 3 s.e. of -0.5 (s.e. of the mean = sqrt(sum s.e.^2)/n).\n"
            "By R2 a reused cache equals a fresh solve with its own seeds, so reuse and the fresh control have the same\n"
            "distribution and their slopes differ by noise only. Where H is near +-1 (g = 1 at the reuse probes) the noise\n"
            "comes from rare sign flips, and a bootstrap over the repeats understates the slope's s.e. there.\n");
        bool c2ok = true;
        for (size_t gr = 0; gr < groups.size(); gr++) {
            Group& G = groups[gr];
            for (Vec3& p : G.probes) pullInside(tet, p);
            const size_t P = G.probes.size();
            const int R = G.R;
            out("\n%s (T = %d, %d repeats per N):\n", G.name, G.tag, R);
            for (size_t p = 0; p < P; p++)
                out("P%zu = %s; dist to cut %.4f\n", gr*P + p + 1, fmt3(G.probes[p]).c_str(), cutP.distToCut(G.probes[p]));
            std::vector<Vec3> pts;
            for (const Vec3& p : G.probes)
                for (int r = 0; r < R; r++) pts.push_back(p);
            for (int gi = 0; gi < 2; gi++) {
                std::vector<std::vector<std::vector<double>>> H(P, std::vector<std::vector<double>>(Ns.size())), HC = H;
                std::vector<double> rew(P, 0.0), kept(P, 0.0);
                long long rewalked = 0, rewalkedK0 = 0;
                CutSolver zs(tet, cutP, gTypes[gi], s.walk);
                for (size_t ni = 0; ni < Ns.size(); ni++) {
                    const int n = Ns[ni];
                    StripCut cut(tet, cfg);
                    cut.addStrip(cut.sOf(0.2));
                    WalkCache C;
                    C.configure(tet, makePoints(pts, nullptr, mixSeed(s.seed, G.tag + 10*gi + int(ni))), n, gTypes[gi], s.walk);
                    C.solveFull(cut);
                    std::vector<uint8_t> ever(C.walks.size(), 0);
                    for (int st = 1; st <= 8; st++) {
                        cut.addStrip(cut.sOf(0.2 + 0.05*st));
                        ReuseStats rs = C.advance(cut, UpdateMode::Incremental);
                        rewalked += rs.nRewalked;
                        rewalkedK0 += rs.nRewalkedAtStep0;
                        for (size_t i = 0; i < pts.size(); i++) rew[i/R] += C.rewalkedFrac[i]/double(8*R*Ns.size());
                        for (size_t w = 0; w < ever.size(); w++) ever[w] |= C.action[w] & 2;
                    }
                    for (size_t w = 0; w < ever.size(); w++)
                        if (!ever[w]) kept[w/(size_t(n)*R)] += 1.0/double(size_t(n)*R*Ns.size());
                    std::vector<double> hc(pts.size());
                    if (G.zombieCtrl) {
                        std::vector<Estimate> ez;
                        zs.solve(pts, nullptr, n, mixSeed(s.seed, G.ctrlTag + 10*gi + int(ni)), ez);
                        for (size_t i = 0; i < pts.size(); i++) hc[i] = ez[i].H;
                    } else {
                        std::vector<EstimateD> ef = solveSumsOnly(tet, cut, makePoints(pts, nullptr, mixSeed(s.seed, G.ctrlTag + 10*gi + int(ni))),
                                                                  n, gTypes[gi], s.walk, nullptr);
                        for (size_t i = 0; i < pts.size(); i++) hc[i] = ef[i].H;
                    }
                    for (size_t p = 0; p < P; p++)
                        for (int r = 0; r < R; r++) {
                            H[p][ni].push_back(C.estimate(p*R + r).H);
                            HC[p][ni].push_back(hc[p*R + r]);
                        }
                }
                const char* ctrl = G.zombieCtrl ? "first gate" : "fresh";
                out("-- %s (re-walked walks: %.1f%% of them from k = 0) --\n", gShortName(gTypes[gi]),
                    100.0*double(rewalkedK0)/double(std::max(1LL, rewalked)));
                double sumS = 0, sumSe2 = 0, sumC = 0, sumCe2 = 0;
                int nS = 0, nC = 0;
                for (size_t p = 0; p < P; p++) {
                    out("P%zu: re-walked %.1f%% of its walks per step, kept %.1f%%\n"
                        "          N      mean        std   std*sqrt(N) | %10s: mean   std\n", gr*P + p + 1, 100.0*rew[p],
                        100.0*kept[p], ctrl);
                    for (size_t ni = 0; ni < Ns.size(); ni++) {
                        const double d = sampleStd(H[p][ni]);
                        out("     %6d   %+.5f   %.5f     %.4f    |       %+.5f  %.5f\n", Ns[ni], meanOf(H[p][ni]), d,
                            d*std::sqrt(double(Ns[ni])), meanOf(HC[p][ni]), sampleStd(HC[p][ni]));
                    }
                    double sl, se, slC, seC;
                    ladderSlope(Ns, H[p], mixSeed(s.seed, 7000 + 100*gi + 10*int(gr) + int(p)), sl, se);
                    ladderSlope(Ns, HC[p], mixSeed(s.seed, 7500 + 100*gi + 10*int(gr) + int(p)), slC, seC);
                    if (!std::isnan(sl) && !std::isnan(se)) {
                        sumS += sl;
                        sumSe2 += se*se;
                        nS++;
                    }
                    if (!std::isnan(slC) && !std::isnan(seC)) {
                        sumC += slC;
                        sumCe2 += seC*seC;
                        nC++;
                    }
                    if (std::isnan(sl))
                        out("     slope of log std vs log N: reuse n/a (std 0 at some N: no walk changed the value)\n");
                    else
                        out("     slope of log std vs log N: reuse %+.3f +- %.3f, %s %+.3f +- %.3f   (1/sqrt(N): -0.5)\n", sl, se,
                            ctrl, slC, seC);
                }
                const double mS = nS ? sumS/nS : std::nan(""), eS = nS ? std::sqrt(sumSe2)/nS : std::nan("");
                const double mC = nC ? sumC/nC : std::nan(""), eC = nC ? std::sqrt(sumCe2)/nC : std::nan("");
                const bool ok = nS > 0 && std::fabs(mS + 0.5) <= 3.0*eS;
                c2ok = c2ok && ok;
                out("   mean slope over %d probes: reuse %+.3f +- %.3f, %s %+.3f +- %.3f over %d (reuse - %s %+.3f +- %.3f);\n"
                    "   verdict: 1/sqrt(N) preserved under reuse (within 3 s.e. of -0.5): %s\n", nS, mS, eS, ctrl, mC, eC, nC,
                    ctrl, mS - mC, std::hypot(eS, eC), yesNo(ok));
            }
        }
        verdicts.push_back(std::string("2  1/sqrt(N) slope preserved under reuse (both probe groups, both g types): ") + yesNo(c2ok));
    }

    // ================================================================ check 3
    {
        const double tStart3 = 0.5;
        const int nSteps3 = 4;
        const std::vector<int> denoms = {8, 16, 32, 64};
        StripCut c0(tet, cfg);
        Slice sl = buildSlice(tet, c0.frame(), s.slice);
        out("\n=== 3. Reuse statistics vs strip size (default scene: lattice k = %d + slice, %zu points, N = %d) ===\n", K,
            lat.pts.size() + sl.verts.size(), N);
        out("full solve at t = %.2f, then %d consecutive steps of dt; means over the %d steps (so the larger dt also\n"
            "cover a later t range, where fewer points lie ahead of the front); 'first' = the first step alone, from the\n"
            "same cache at t = %.2f for every dt. t_full = same-seed full re-solve on the same geometry right after each\n"
            "update (the cache's own full-solve path). Seeds: tag 800 + g (lattice 1, slice 2 below it).\n",
            tStart3, nSteps3, nSteps3, tStart3);
        out("  dt     t range       strip/full | tested  re-walked (k=0) | first: tested re-walked | value-upd (a)  (b) |"
            " t_update (a)    (b)    | t_full  | full/update (a) (b)\n");
        std::string firstCounts, costLine;
        bool sameSets3 = true, allFaster = true, firstShrinks = true, firstSame = true;
        long long prevFirstRew = -1;
        double speedMin = std::numeric_limits<double>::infinity(), speedMax = 0;
        for (int den : denoms) {
            const double dt = 1.0/double(den);
            double area = 0, tested = 0, rew = 0, val[2] = {0, 0}, tUpd[2] = {0, 0}, tFull = 0;
            double firstTested = 0, firstRew = 0;
            long long cnt[2][5] = {};   // per g: tested, re-walked, at k = 0 (all steps), first step tested, re-walked
            for (int gi = 0; gi < 2; gi++) {
                QueryPoints q;   // the same seeds for both g types: the same paths, only the values differ
                q.add(lat.pts, &lat.interior, mixSeed(mixSeed(s.seed, 800), 1));
                q.add(sl.verts, nullptr, mixSeed(mixSeed(s.seed, 800), 2));
                StripCut cut(tet, cfg);
                cut.addStrip(cut.sOf(tStart3));
                WalkCache C, Fc;
                C.configure(tet, q, N, gTypes[gi], s.walk);
                Fc.configure(tet, q, N, gTypes[gi], s.walk);
                C.solveFull(cut);
                for (int st = 1; st <= nSteps3; st++) {
                    int m = cut.addStrip(cut.sOf(tStart3 + dt*st));
                    ReuseStats rs = C.advance(cut, UpdateMode::Incremental);
                    SolveStats sf = Fc.solveFull(cut);
                    const double nW = double(rs.nWalks);
                    cnt[gi][0] += rs.nTested;
                    cnt[gi][1] += rs.nRewalked;
                    cnt[gi][2] += rs.nRewalkedAtStep0;
                    if (st == 1) {
                        cnt[gi][3] = rs.nTested;
                        cnt[gi][4] = rs.nRewalked;
                    }
                    if (gi == 0) {
                        area += cut.strip(m).area/cut.fullArea()/nSteps3;
                        tested += rs.nTested/nW/nSteps3;
                        rew += rs.nRewalked/nW/nSteps3;
                        if (st == 1) {
                            firstTested = rs.nTested/nW;
                            firstRew = rs.nRewalked/nW;
                            char buf[160];
                            std::snprintf(buf, sizeof buf, "  1/%-3d first step: %lld tested, %lld re-walked (%lld at k = 0) of %lld walks\n",
                                          den, rs.nTested, rs.nRewalked, rs.nRewalkedAtStep0, rs.nWalks);
                            firstCounts += buf;
                        }
                    }
                    val[gi] += rs.nValueUpdated/nW/nSteps3;
                    tUpd[gi] += rs.tTotal/nSteps3;
                    tFull += sf.seconds/(2*nSteps3);
                }
            }
            for (int c = 0; c < 5; c++) sameSets3 = sameSets3 && cnt[0][c] == cnt[1][c];
            // a smaller strip should re-walk fewer walks on the same first step (from the same cache)
            if (prevFirstRew >= 0) {
                firstShrinks = firstShrinks && cnt[0][4] < prevFirstRew;
                firstSame = firstSame && cnt[0][4] == prevFirstRew;
            }
            prevFirstRew = cnt[0][4];
            {
                char buf[96];
                std::snprintf(buf, sizeof buf, "%s 1/%d: %.2f vs %.2f", costLine.empty() ? "" : ";", den,
                              0.5*(tUpd[0] + tUpd[1])/tFull, tested);
                costLine += buf;
            }
            for (int gi = 0; gi < 2; gi++) {
                speedMin = std::min(speedMin, tFull/tUpd[gi]);
                speedMax = std::max(speedMax, tFull/tUpd[gi]);
                allFaster = allFaster && tUpd[gi] < tFull;
            }
            out("  1/%-3d  %.3f-%.3f   %.4f     | %6.2f%%  %6.2f%% (%4.1f%%) |  %6.2f%%  %6.2f%%     |  %6.2f%%  %6.2f%% |"
                "  %.4f s  %.4f s | %.3f s |  %5.1fx  %5.1fx\n", den, tStart3, tStart3 + nSteps3*dt, area, 100*tested, 100*rew,
                cnt[0][1] > 0 ? 100.0*double(cnt[0][2])/double(cnt[0][1]) : 0.0, 100*firstTested, 100*firstRew, 100*val[0],
                100*val[1], tUpd[0], tUpd[1], tFull, tFull/tUpd[0], tFull/tUpd[1]);
        }
        out("(t_full is the mean over both g types; (k=0) = share of the re-walked walks invalidated at their first ball)\n"
            "tested and re-walked counts of (a) and (b) identical at every dt (same seeds, same paths), so the columns\n"
            "above stand for both: %s\n", yesNo(sameSets3));
        out("%s", firstCounts.c_str());
        out("Why the first step does not depend on dt: every step ball is empty of the old cut and touches it at its\n"
            "closest point. Inside the cut: the ball meets the plane only at that point, so no strip can enter it. On a\n"
            "lateral edge (on a tet face): the ball's disk in the plane lies beyond that edge's line, outside the section,\n"
            "so no strip can enter it. On the front edge: the disk is tangent to the front line from ahead, so every\n"
            "strip, however thin, enters it. At a front corner the disk's centre lies in the corner's normal cone, and\n"
            "the disk still contains the points of every strip next to the corner along the lateral edge, so these\n"
            "balls do not depend on dt either. So the invalidated set depends on where the front is, not on how far it\n"
            "moves, and a smaller dt only means more updates for the same front travel. The cut is the only absorbing\n"
            "boundary (the tet faces reflect), so every point ahead of the front has its closest absorbing point on\n"
            "the front edge, and its first ball is invalidated by any advance.\n");
        out("Cost: a tested walk is replayed from its seed, which costs about as much as walking it, so the update\n"
            "costs about (fraction tested) x (full solve), and a little more here because long walks have large boxes\n"
            "and are tested more often. t_update/t_full vs fraction tested (means over the steps):\n"
            "  %s\n", costLine.c_str());
        out("verdict: the first step's re-walked fraction shrinks with the strip size: %s%s\n", yesNo(firstShrinks),
            firstSame ? " (the first step's counts are identical\n   for every dt: the invalidated set does not depend on dt, as argued above)" : "");
        out("verdict: the update is faster than the same-seed full re-solve at every strip size: %s (speed-up %.1fx to %.1fx)\n",
            yesNo(allFaster), speedMin, speedMax);
        verdicts.push_back(std::string("3  first step's re-walked fraction shrinks with the strip size: ") + yesNo(firstShrinks) +
                           (firstSame ? " (identical for every dt: the invalidated set does not depend on dt here)" : ""));
        char v3[160];
        std::snprintf(v3, sizeof v3, "3  update faster than a same-N full re-solve at every strip size: %s (%.1fx to %.1fx)",
                      yesNo(allFaster), speedMin, speedMax);
        verdicts.push_back(v3);
    }

    // ================================================================ check 4
    {
        std::vector<std::pair<std::string, Vec3>> probes = {
            {"near the front at t=0.4, above the cut (F(0.4) + 0.03L n)", buildCut(tet, cfg, 0.4).frontMidpoint() + 0.03*L*F.n},
            {"near the front at t=0.6, below the cut (F(0.6) - 0.03L n)", buildCut(tet, cfg, 0.6).frontMidpoint() - 0.03*L*F.n},
            {"ahead of the front until t~0.8 (F(0.8) + 0.05L n)", buildCut(tet, cfg, 0.8).frontMidpoint() + 0.05*L*F.n},
            {"near the initial cut, far from the front (F(0.1) + 0.02L n)", buildCut(tet, cfg, 0.1).frontMidpoint() + 0.02*L*F.n},
            {"far from the cut (0.8 v3 + 0.2 centroid)", 0.8*tet.v[3] + 0.2*tet.centroid},
            {"tet centroid", tet.centroid}};
        out("\n=== 4. Temporal coherence over the sweep t = 0.2 -> 1.0 in 16 steps (6 probes, N = %d) ===\n", N);
        std::vector<Vec3> pts;
        Cut cA = buildCut(tet, cfg, 0.2), cB = buildCut(tet, cfg, 1.0);
        for (size_t p = 0; p < probes.size(); p++) {
            int moved = pullInside(tet, probes[p].second);
            pts.push_back(probes[p].second);
            out("Q%zu = %s  %s; dist to cut at t=0.2: %.4f, at t=1: %.4f%s\n", p + 1, fmt3(probes[p].second).c_str(),
                probes[p].first.c_str(), cA.distToCut(probes[p].second), cB.distToCut(probes[p].second),
                moved ? " (moved towards the centroid)" : "");
        }
        out("incremental: one cache (tag 900 + g), full solve at t = 0.2 then 16 incremental steps; fresh: a full solve\n"
            "with fresh seeds at every t (tag 3000 + 100 g + step); reference: 64 N walks, fresh seeds (tag 4000 + 100 g +\n"
            "step), to separate the true change of H from the noise. Mean |H(t_i) - H(t_(i-1))| over the 16 steps:\n");
        for (int gi = 0; gi < 2; gi++) {
            StripCut cut(tet, cfg);
            cut.addStrip(cut.sOf(0.2));
            WalkCache C;
            C.configure(tet, makePoints(pts, nullptr, mixSeed(s.seed, 900 + gi)), N, gTypes[gi], s.walk);
            C.solveFull(cut);
            std::vector<std::vector<double>> hi(pts.size()), hf(pts.size()), hr(pts.size());
            auto record = [&](int step) {
                std::vector<EstimateD> fr = solveSumsOnly(tet, cut, makePoints(pts, nullptr, mixSeed(s.seed, 3000 + 100*gi + step)), N,
                                                          gTypes[gi], s.walk, nullptr);
                std::vector<EstimateD> rf = solveSumsOnly(tet, cut, makePoints(pts, nullptr, mixSeed(s.seed, 4000 + 100*gi + step)),
                                                          64*N, gTypes[gi], s.walk, nullptr);
                for (size_t p = 0; p < pts.size(); p++) {
                    hi[p].push_back(C.estimate(p).H);
                    hf[p].push_back(fr[p].H);
                    hr[p].push_back(rf[p].H);
                }
            };
            record(0);
            for (int st = 1; st <= 16; st++) {
                cut.addStrip(cut.sOf(0.2 + 0.05*st));
                C.advance(cut, UpdateMode::Incremental);
                record(st);
            }
            out("-- %s --\n   t     ", gShortName(gTypes[gi]));
            for (size_t p = 0; p < pts.size(); p++) out(" Q%zu inc   fresh  ", p + 1);
            out("\n");
            for (int st = 0; st <= 16; st++) {
                out("  %.2f  ", 0.2 + 0.05*st);
                for (size_t p = 0; p < pts.size(); p++) out(" %+.4f %+.4f", hi[p][st], hf[p][st]);
                out("\n");
            }
            double sumI = 0, sumF = 0, sumR = 0, nzI = 0, nzF = 0;
            out("   mean |dH| (incremental / fresh / reference):");
            for (size_t p = 0; p < pts.size(); p++) {
                double di = 0, df = 0, dr = 0;
                for (int st = 1; st <= 16; st++) {
                    const double dRef = hr[p][st] - hr[p][st - 1];
                    di += std::fabs(hi[p][st] - hi[p][st - 1])/16.0;
                    df += std::fabs(hf[p][st] - hf[p][st - 1])/16.0;
                    dr += std::fabs(dRef)/16.0;
                    nzI += std::fabs(hi[p][st] - hi[p][st - 1] - dRef)/(16.0*double(pts.size()));
                    nzF += std::fabs(hf[p][st] - hf[p][st - 1] - dRef)/(16.0*double(pts.size()));
                }
                sumI += di/double(pts.size());
                sumF += df/double(pts.size());
                sumR += dr/double(pts.size());
                out("\n      Q%zu %.4f / %.4f / %.4f", p + 1, di, df, dr);
            }
            out("\n   over all 6 probes: incremental %.4f, fresh seeds %.4f (incremental / fresh = %.2f); reference %.4f\n", sumI,
                sumF, sumF > 0 ? sumI/sumF : 0.0, sumR);
            out("   noise in the step-to-step change, mean |dH - dH_ref|: incremental %.4f, fresh seeds %.4f (ratio %.2f)\n", nzI,
                nzF, nzF > 0 ? nzI/nzF : 0.0);
            const bool rawOk = sumI < sumF, noiseOk = nzI < nzF;
            out("   verdict: incremental step-to-step change smaller than with fresh seeds: %s (ratio %.2f); in its noise part:\n"
                "   %s (ratio %.2f)\n", yesNo(rawOk), sumF > 0 ? sumI/sumF : 0.0, yesNo(noiseOk), nzF > 0 ? nzI/nzF : 0.0);
            if (!rawOk)
                out("   (not met in the raw mean |dH|: the true change per step, mean |dH_ref| = %.4f, is %.1fx the incremental\n"
                    "   noise part, and |dH_true + e| averages to about |dH_true| for zero-mean noise e much smaller than\n"
                    "   dH_true, so the raw mean |dH| measures the true change in both modes, not the noise)\n", sumR,
                    nzI > 0 ? sumR/nzI : 0.0);
            char v4[200];
            std::snprintf(v4, sizeof v4, "4  %s: incremental smaller step-to-step change than fresh seeds: %s (%.2f); noise part: %s (%.2f)",
                          gShortName(gTypes[gi]), yesNo(rawOk), sumF > 0 ? sumI/sumF : 0.0, yesNo(noiseOk), nzF > 0 ? nzI/nzF : 0.0);
            verdicts.push_back(v4);
        }
    }

    double sec = secondsSince(tStart);
    out("\n=== Verdicts of the statistical checks (criteria stated with each check) ===\n");
    for (const std::string& v : verdicts) out("%s\n", v.c_str());
    out("\nexact checks R1, R2, R3, R4: %s\n", exactOk ? "all PASS" : "SOME FAILED");
    out("reuse validation wall time: %.1f s\n", sec);
    if (!writeResultsSection(resultsPath, kResultsSectionReuse, out.text)) {
        std::fprintf(stderr, "could not write %s\n", resultsPath.c_str());
        return false;
    }
    std::printf("written: %s (section \"%s\")\n", resultsPath.c_str(), kResultsSectionReuse);
    return exactOk;
}

// =====================================================================================
int runReuseBench(const GeometryConfig& cfg, const Settings& s)
{
    Tet tet = referenceTet();
    Lattice lat = buildLattice(tet, Cut{}, s.k, 0.0);
    StripCut c0(tet, cfg);
    Slice sl = buildSlice(tet, c0.frame(), s.slice);
    std::printf("bench-reuse: default scene, lattice k = %d (%zu points, all cached) + slice (%zu points), N = %d;\n"
                "full solve at t = 0.2, then 16 incremental steps of 0.05 to t = 1; t_full = same-seed full re-solve\n"
                "on the same geometry after each update\n", s.k, lat.pts.size(), sl.verts.size(), s.nWalks);
    for (GType g : {GType::One, GType::SqrtFront}) {
        QueryPoints q;
        q.add(lat.pts, &lat.interior, mixSeed(s.seed, 1));
        q.add(sl.verts, nullptr, mixSeed(s.seed, 2));
        StripCut cut(tet, cfg);
        cut.addStrip(cut.sOf(0.2));
        WalkCache C, Fc;
        C.configure(tet, q, s.nWalks, g, s.walk);
        Fc.configure(tet, q, s.nWalks, g, s.walk);
        SolveStats s0 = C.solveFull(cut);
        std::printf("\n-- %s: initial full solve: %s --\n", gTypeName(g), s0.summary().c_str());
        std::printf("step    t    strip/full  tested  re-walked  value-upd | t_box   t_rewalk t_value t_sums | t_update | t_full | full/update\n");
        double sumU = 0, sumF = 0;
        long long esc = 0, cap = 0, stRep = 0, stNew = 0, stFull = 0;
        for (int st = 1; st <= 16; st++) {
            double t = 0.2 + 0.05*st;
            int m = cut.addStrip(cut.sOf(t));
            ReuseStats rs = C.advance(cut, UpdateMode::Incremental);
            SolveStats sf = Fc.solveFull(cut);
            const double nW = double(rs.nWalks);
            sumU += rs.tTotal;
            sumF += sf.seconds;
            esc += rs.escaped;
            cap += rs.capped;
            stRep += rs.stepsReplayed;
            stNew += rs.stepsNew;
            stFull += sf.totalSteps;
            std::printf("%3d  %.3f   %.4f   %6.2f%%  %6.2f%%   %6.2f%%  | %.4f  %.4f  %.4f  %.4f | %.4f s | %.3f s | %6.1fx\n", st,
                        t, cut.strip(m).area/cut.fullArea(), 100.0*rs.nTested/nW, 100.0*rs.nRewalked/nW,
                        100.0*rs.nValueUpdated/nW, rs.tBox, rs.tRewalk, rs.tValue, rs.tSums, rs.tTotal, sf.seconds,
                        sf.seconds/rs.tTotal);
        }
        std::printf("   16 steps: incremental %.3f s in total, full re-solves %.3f s (%.1fx); last full solve: %s\n", sumU, sumF,
                    sumF/sumU, Fc.solveFull(cut).summary().c_str());
        std::printf("   re-walked paths: %lld escaped, %lld capped; the updates ran %lld replayed (old cut) + %lld new steps,\n"
                    "   %.1f%% of the %lld steps of the 16 full re-solves\n", esc, cap, stRep, stNew,
                    stFull ? 100.0*double(stRep + stNew)/double(stFull) : 0.0, stFull);
    }
    return 0;
}

// =====================================================================================
// --compare-options: the exact reuse options 1, 3, 2a and 2b side by side, and the control 3q
// =====================================================================================
namespace {

double mb(size_t bytes) { return double(bytes)/1048576.0; }

// the process's private bytes (Windows): a cross-check of the computed memory, 0 elsewhere
size_t privateBytes()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc), sizeof pmc))
        return size_t(pmc.PrivateUsage);
#endif
    return 0;
}

// The machine, for the header: the timings and their ratios depend on it and on the thread count. Windows: the
// CPU's brand string, its cores by efficiency class (on a hybrid CPU the higher class is the faster core), the
// logical processors and the L3 size; "" elsewhere.
std::string cpuInfo()
{
    std::string s;
#ifdef _WIN32
    int r[4];
    char brand[49] = {};
    __cpuid(r, int(0x80000000));
    if (unsigned(r[0]) >= 0x80000004u)
        for (int k = 0; k < 3; k++) {
            __cpuid(r, int(0x80000002 + k));
            std::memcpy(brand + 16*k, r, 16);
        }
    s = brand;
    s.erase(0, s.find_first_not_of(' '));
    DWORD len = 0;
    GetLogicalProcessorInformationEx(RelationAll, nullptr, &len);
    std::vector<char> buf(len);
    int byClass[256] = {}, nCores = 0, nLogical = 0;
    size_t l3 = 0;
    if (len && GetLogicalProcessorInformationEx(RelationAll, reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buf.data()),
                                                &len))
        for (DWORD off = 0; off < len;) {
            const auto* p = reinterpret_cast<const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buf.data() + off);
            if (p->Relationship == RelationProcessorCore) {
                nCores++;
                byClass[p->Processor.EfficiencyClass]++;
                for (WORD gi = 0; gi < p->Processor.GroupCount; gi++)
                    nLogical += std::popcount(uint64_t(p->Processor.GroupMask[gi].Mask));
            } else if (p->Relationship == RelationCache && p->Cache.Level == 3) {
                l3 += p->Cache.CacheSize;
            }
            off += p->Size;
        }
    char t[96];
    std::snprintf(t, sizeof t, "; %d cores", nCores);
    s += t;
    for (int c = 255; c >= 0; c--)
        if (byClass[c]) {
            std::snprintf(t, sizeof t, ", %d of efficiency class %d", byClass[c], c);
            s += t;
        }
    std::snprintf(t, sizeof t, "; %d logical processors; L3 %.0f MiB", nLogical, mb(l3));
    s += t;
#endif
    return s;
}

template <class F> double timeIt(F&& f)
{
    auto t0 = std::chrono::steady_clock::now();
    f();
    return secondsSince(t0);
}

double median(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n == 0 ? 0.0 : n%2 ? v[n/2] : 0.5*(v[n/2 - 1] + v[n/2]);
}

// The brute-force ball test behind the 2b soundness check: for every walk of a stored-step cache, the first
// stored ball that strip m enters under the update's own criterion stripDist2(c_k) < key2_k, applied to
// every ball of every walk with no prefilter at all; -1 if none.
std::vector<int> firstEnteredBall(const WalkCache& C, const StripCut& cut, int m)
{
    std::vector<int> E(C.paths.size(), -1);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, C.paths.size()), [&](const tbb::blocked_range<size_t>& r) {
        for (size_t wi = r.begin(); wi < r.end(); wi++) {
            const WalkPath& p = C.paths[wi];
            for (size_t k = 0; k < p.s.size(); k++)
                if (cut.stripDist2(p.s[k].c, m) < p.s[k].key2) {
                    E[wi] = int(k);
                    break;
                }
        }
    });
    return E;
}

// walks whose entry (box, landing point, value, steps, side, status, RNG offset) or gradient weight differs bitwise
long long walkMismatches(const WalkCache& a, const WalkCache& b)
{
    std::atomic<long long> bad{0};
    tbb::parallel_for(tbb::blocked_range<size_t>(0, a.walks.size()), [&](const tbb::blocked_range<size_t>& r) {
        long long n = 0;
        for (size_t wi = r.begin(); wi < r.end(); wi++)
            if (!sameWalk(a.walks[wi], a.gradWeight(wi), b.walks[wi], b.gradWeight(wi), false)) n++;
        bad += n;
    });
    return bad.load();
}

// walks whose stored path (centres, keys, flags) differs bitwise
long long pathMismatches(const WalkCache& a, const WalkCache& b)
{
    std::atomic<long long> bad{0};
    tbb::parallel_for(tbb::blocked_range<size_t>(0, a.paths.size()), [&](const tbb::blocked_range<size_t>& r) {
        long long n = 0;
        for (size_t wi = r.begin(); wi < r.end(); wi++) {
            const WalkPath &p = a.paths[wi], &q = b.paths[wi];
            if (p.s.size() != q.s.size() || p.f != q.f ||
                (!p.s.empty() && !sameBits(p.s.data(), q.s.data(), p.s.size()*sizeof(StoredStep))))
                n++;
        }
        bad += n;
    });
    return bad.load();
}

// points whose sums differ bitwise
long long sumsMismatches(const std::vector<PointSums>& a, const std::vector<PointSums>& b)
{
    long long n = 0;
    for (size_t i = 0; i < a.size(); i++) n += !sameBits(&a[i], &b[i], sizeof(PointSums));
    return n;
}

// 3q against 2b: walks whose front key differs bitwise or that one of them popped and the other did not, plus
// the queue positions that differ (both queues must be the same sequence of (key, walk))
long long keyMismatches(const WalkCache& a, const WalkCache& b)
{
    long long n = 0;
    for (size_t wi = 0; wi < a.frontKey.size(); wi++) {
        n += !sameBits(&a.frontKey[wi], &b.frontKey[wi], sizeof(double)) || ((a.action[wi] ^ b.action[wi]) & 1) != 0;
        n += !sameBits(&a.queue[wi].first, &b.queue[wi].first, sizeof(double)) || a.queue[wi].second != b.queue[wi].second;
    }
    return n;
}

// H and grad H of an option against option 1's (the same walks, summed in another order): the number of points
// whose difference exceeds a worst-case round-off bound; dH and dG are raised to the largest differences.
// The bound: a point's running sum never holds more than N + 1 terms (an update subtracts a walk's old term
// before it adds the new one), so each addition rounds by at most eps/2 (N + 1) T, T the largest |term|. After U
// updates the incremental sum has had at most N (1 + 2U) additions and option 1's N, so the two differ by at most
// A T, A = N (2 + 2U) eps/2 (N + 1). H = sH/n (|u| <= 1, T = 1): 2A/n covers the division as well. grad H =
// (Σwu - H Σw)/n (T = wMax, the largest |w| component of the point's walks over the run): the two sums, H's error
// times |Σw| <= (N + 1) wMax and the last roundings stay below (4A + (2A/n + 8 eps)(N + 1)) wMax / n. Both sides
// add the same rounded terms (the same w and u), so only the additions round differently. g = 1: every term of
// sH is an integer, so H must agree exactly. A bug in the shared value updates or sums (an old term subtracted
// with the wrong u or w) shows up at the size of one term over n, far above the bound.
long long estimateCheck(const std::vector<PointSums>& a, const std::vector<PointSums>& ref, const std::vector<char>& grad,
                        const std::vector<double>& wMax, int N, int U, bool exactH, double& dH, double& dG)
{
    const double A = double(N)*(2.0 + 2.0*double(U))*0.5*std::numeric_limits<double>::epsilon()*double(N + 1);
    long long bad = 0;
    for (size_t i = 0; i < a.size(); i++) {
        const EstimateD ea = estimateFromSums(a[i], grad[i] != 0), er = estimateFromSums(ref[i], grad[i] != 0);
        const double n = double(std::max(1LL, er.n)), bH = 2.0*A/n;
        const double eH = std::fabs(ea.H - er.H);
        dH = std::max(dH, eH);
        bool b = ea.n != er.n || ea.hasGrad != er.hasGrad || !(eH <= (exactH ? 0.0 : bH));
        if (ea.hasGrad && er.hasGrad) {
            const double eG = (ea.grad - er.grad).cwiseAbs().maxCoeff();
            dG = std::max(dG, eG);
            b = b || !(eG <= (4.0*A + (bH + 8.0*std::numeric_limits<double>::epsilon())*double(N + 1))*wMax[i]/n);
        }
        bad += b;
    }
    return bad;
}

// The options, in the order of every table. 1 keeps no per-walk state; 3, 3q, 2a and 2b are WalkCaches. 3q is
// not one of the three options asked about: it is the control that separates what 2b gains by its tighter filter
// (the front keys, which need no stored steps) from what it gains by storing every step.
enum { O1 = 0, O3 = 1, O3Q = 2, O2A = 3, O2B = 4, kOpts = 5 };
const char* const kOptName[kOpts] = {"1  frozen seeds, full recompute", "3  seed replay (WalkCache::advance)",
                                     "3q seed replay + front-key queue", "2a stored steps + bbox",
                                     "2b stored steps + front-key queue"};
const char* const kOptShort[kOpts] = {"1", "3", "3q", "2a", "2b"};

// one option's totals over a run
struct OptTotals {
    int updates = 0;
    double t = 0, tFilter = 0, tTest = 0, tResume = 0, tNew = 0, tKeys = 0, tValue = 0, tSums = 0, tQueue = 0;
    long long cand = 0, rewalked = 0, valueUpd = 0, balls = 0, stepsNew = 0, stepsReplayed = 0, replayFalse = 0,
              sumK = 0, keys = 0;
    void add(double tWall, const ReuseStats& r)
    {
        updates++;
        t += tWall;
        tFilter += r.tBox;
        tTest += r.tTest;
        tResume += r.tResume;
        tNew += r.tRewalk;
        tKeys += r.tKeys;
        tValue += r.tValue;
        tSums += r.tSums;
        tQueue += r.tQueue;
        cand += r.nTested;
        rewalked += r.nRewalked;
        valueUpd += r.nValueUpdated;
        balls += r.ballsTested;
        stepsNew += r.stepsNew;
        stepsReplayed += r.stepsReplayed;
        replayFalse += r.stepsReplayedFalse;
        sumK += r.sumK;
        keys += r.keysComputed;
    }
};

// the equivalence and soundness checks (counts must all be 0); per-option entries are indexed O3 .. O2B
struct CompareChecks {
    long long initial = 0;              // initial solves of 3q, 2a, 2b, the reference != 3's; paths, keys != the reference's, 2b's
    long long refSums = 0;              // points: option 1's sums != the same-seed solve with entries (bitwise)
    long long refRecord = 0;            // walks: recording paths changed a fresh solve's entry
    long long walks[kOpts] = {};        // walks != option 1's same-seed walks (bitwise, incl. w)
    long long paths[kOpts] = {};        // 2a, 2b: stored paths != a fresh recorded solve's (bitwise)
    long long sumsVs3[kOpts] = {};      // 3q, 2a, 2b: points whose sums != option 3's (bitwise)
    long long estimate[kOpts] = {};     // points: H or grad H off option 1's beyond the round-off bound (estimateCheck)
    long long sets[kOpts] = {};         // walks re-walked (or not) against the brute-force ball test
    long long keys3q = 0;               // 3q's keys, pops or queue != 2b's
    long long missed2b = 0;             // the strip enters a stored ball, but 2b did not pop the walk
    long long boxMissed2a = 0;          // ... but 2a's bbox test rejected the walk
    long long fullFallback = 0;         // updates that fell back to a full solve
    long long entered = 0, false2b = 0, boxFP2a = 0;   // information: entered walks, 2b false alarms, bbox false positives
    double dH[kOpts] = {}, dG[kOpts] = {};             // max |H - H_1|, |grad H - grad H_1|
    long long bad() const
    {
        long long b = initial + refSums + refRecord + keys3q + missed2b + boxMissed2a + fullFallback;
        for (int k = 0; k < kOpts; k++) b += walks[k] + paths[k] + sumsVs3[k] + estimate[k] + sets[k];
        return b;
    }
    void add(const CompareChecks& o)
    {
        initial += o.initial;
        refSums += o.refSums;
        refRecord += o.refRecord;
        keys3q += o.keys3q;
        missed2b += o.missed2b;
        boxMissed2a += o.boxMissed2a;
        fullFallback += o.fullFallback;
        entered += o.entered;
        false2b += o.false2b;
        boxFP2a += o.boxFP2a;
        for (int k = 0; k < kOpts; k++) {
            walks[k] += o.walks[k];
            paths[k] += o.paths[k];
            sumsVs3[k] += o.sumsVs3[k];
            estimate[k] += o.estimate[k];
            sets[k] += o.sets[k];
            dH[k] = std::max(dH[k], o.dH[k]);
            dG[k] = std::max(dG[k], o.dG[k]);
        }
    }
};

struct RunOut {
    OptTotals opt[kOpts];
    CompareChecks chk;
    CacheMemory mem0[kOpts], memPeak[kOpts];   // after the initial solve; after the update with the largest total + temporaries
    size_t scratchAtPeak[kOpts] = {};  // that update's temporaries
    size_t scratchMax[kOpts] = {};     // the largest temporaries of any update
    size_t mem1 = 0;                   // option 1: its per-point sums (plus the query points)
    size_t os[kOpts] = {};             // process private bytes added by configuring each cache and its initial solve
    size_t osGrowth = 0, osComputed = 0;   // largest sample after a step: private bytes over the run's start, and the
                                           // computed total of every cache alive then (the four options' and two references)
    double setup[kOpts] = {};          // each option's own same-seed initial solve at tStart (option 1: solveSums)
    double tRef = 0;                   // the per-step reference solves (validation only)
    double tRefPlain = 0, tRefRecord = 0;   // ... split: the plain same-seed full solve, the one recording paths
    long long nWalks = 0, fullSteps = 0;   // fullSteps: steps of option 1's full recomputes over the run
    std::vector<double> frac;          // per update: the share of walks re-walked (the same for every option)
    std::vector<double> stepT[kOpts];  // per update: each option's wall time
    bool first = false;                // the process's first run: its OS figures start from a heap no run has used
};

// One run: each option's own initial solve at tStart, then one strip per front in ts. Per step: the brute-force
// ball test (validation; on the reference's recorded paths, so it warms no option's data), the five timed
// updates (same clock, same threads, a random order per step from orderSeed; option 1 is a full recompute),
// then the references (a same-seed full solve with per-walk entries and one recording paths) and every
// comparison. Validation work is never inside a timed region. tables: print the per-step tables.
RunOut compareRun(Report& out, const Tet& tet, const GeometryConfig& cfg, const QueryPoints& q, int N, GType g,
                  const WalkParams& wp, double tStart, const std::vector<double>& ts, bool warmUp, uint64_t orderSeed,
                  bool tables)
{
    RunOut R;
    R.first = warmUp;
    StripCut cut(tet, cfg);
    cut.addStrip(cut.sOf(tStart));
    const size_t P = q.size();
    R.nWalks = (long long)(P*size_t(N));
    const StepStore store[kOpts] = {StepStore::None, StepStore::None, StepStore::Keys, StepStore::Paths,
                                    StepStore::PathsAndKeys};
    const UpdateMode mode[kOpts] = {UpdateMode::FullResolve, UpdateMode::Incremental, UpdateMode::ReplayQueue,
                                    UpdateMode::StoredBox, UpdateMode::StoredQueue};
    WalkCache C[kOpts], Fc, Cr;   // C[O1] stays empty: option 1 keeps no per-walk state
    if (warmUp) solveSums(tet, cut, q, N, g, wp, nullptr);   // TBB's threads and stacks, before any timing or measurement
    R.setup[O1] = timeIt([&] { solveSums(tet, cut, q, N, g, wp, nullptr); });
    const size_t mStart = privateBytes();
    size_t mPrev = mStart;
    SolveStats st0[kOpts];
    for (int k = O3; k < kOpts; k++) {
        C[k].configure(tet, q, N, g, wp, store[k]);
        R.setup[k] = timeIt([&] { st0[k] = C[k].solveFull(cut); });
        const size_t mNow = privateBytes();
        R.os[k] = mNow > mPrev ? mNow - mPrev : 0;
        mPrev = mNow;
    }
    Fc.configure(tet, q, N, g, wp, StepStore::None);
    Cr.configure(tet, q, N, g, wp, StepStore::Paths);
    Cr.solveFull(cut);   // the reference paths at tStart, for step 1's brute-force ball test
    for (int k = O3Q; k < kOpts; k++) R.chk.initial += walkMismatches(C[k], C[O3]) + sumsMismatches(C[k].sums, C[O3].sums);
    R.chk.initial += walkMismatches(Cr, C[O3]) + pathMismatches(C[O2A], Cr) + pathMismatches(C[O2B], Cr) +
                     keyMismatches(C[O3Q], C[O2B]);
    std::vector<double> wMax(P, 0.0);   // per point: the largest |w| component of its walks so far (round-off bound)
    auto raiseWMax = [&] {
        for (size_t wi = 0; wi < size_t(R.nWalks); wi++)
            if (const Vec3* w = C[O3].gradWeight(wi)) wMax[wi/size_t(N)] = std::max(wMax[wi/size_t(N)], w->cwiseAbs().maxCoeff());
    };
    raiseWMax();
    for (int k = O3; k < kOpts; k++) R.mem0[k] = R.memPeak[k] = C[k].memory();
    R.mem1 = P*sizeof(PointSums) + P*(sizeof(Vec3) + 1 + sizeof(uint64_t));
    out("initial full solve at t = %.4f: %s\n", tStart, st0[O3].summary().c_str());
    out("setup (each option's own same-seed initial solve): 1 %.3f s (sums only), 3 %.3f s, 3q %.3f s (+ keys while "
        "walking, queue), 2a %.3f s (recording every step), 2b %.3f s (+ keys and queue); initial solves bitwise "
        "identical: %s\n", R.setup[O1], R.setup[O3], R.setup[O3Q], R.setup[O2A], R.setup[O2B],
        R.chk.initial == 0 ? "yes" : "NO");

    std::string rowsT, rowsC, rowsP;
    char buf[2048];
    pcg32 ord(orderSeed, 1);
    int nUpd = 0;
    for (size_t st = 0; st < ts.size(); st++) {
        const int m = cut.addStrip(cut.sOf(ts[st]));
        if (m < 0) continue;
        nUpd++;
        // validation before the updates: the brute-force test of every stored ball, on the reference's paths
        // (bitwise those of 2a and 2b before this update: checked at the previous step and initially)
        const std::vector<int> E = firstEnteredBall(Cr, cut, m);
        long long nEntered = 0;
        for (int e : E) nEntered += e >= 0;
        const double frac = double(nEntered)/double(R.nWalks);
        // the five timed updates, in a random order per step (fixed seed): no option systematically runs first
        // or right after a given other one
        int order[kOpts] = {O1, O3, O3Q, O2A, O2B};
        for (int i = kOpts - 1; i > 0; i--) std::swap(order[i], order[ord.nextUInt(uint32_t(i + 1))]);
        std::vector<PointSums> sums1;
        SolveStats s1;
        ReuseStats r[kOpts];
        double t[kOpts] = {};
        for (int k : order) {
            if (k == O1) t[k] = timeIt([&] { sums1 = solveSums(tet, cut, q, N, g, wp, &s1); });
            else t[k] = timeIt([&] { r[k] = C[k].advance(cut, mode[k]); });
        }
        r[O1].nTested = r[O1].nRewalked = R.nWalks;   // option 1 as an "update": every walk re-walked
        r[O1].stepsNew = s1.totalSteps;
        r[O1].tRewalk = t[O1];
        R.frac.push_back(frac);
        for (int k = 0; k < kOpts; k++) {
            R.opt[k].add(t[k], r[k]);
            R.stepT[k].push_back(t[k]);
        }
        R.fullSteps += s1.totalSteps;
        // references
        const double tRefPlain = timeIt([&] { Fc.solveFull(cut); });
        const double tRefRec = timeIt([&] { Cr.solveFull(cut); });
        const double tRef = tRefPlain + tRefRec;
        R.tRef += tRef;
        R.tRefPlain += tRefPlain;
        R.tRefRecord += tRefRec;
        // checks
        raiseWMax();
        CompareChecks c;
        c.refSums = sumsMismatches(sums1, Fc.sums);
        c.refRecord = walkMismatches(Cr, Fc);
        for (int k = O3; k < kOpts; k++) {
            c.fullFallback += r[k].full;
            c.walks[k] = walkMismatches(C[k], Fc);
            if (k != O3) c.sumsVs3[k] = sumsMismatches(C[k].sums, C[O3].sums);
            c.estimate[k] = estimateCheck(C[k].sums, sums1, q.grad, wMax, N, nUpd, g == GType::One, c.dH[k], c.dG[k]);
        }
        c.paths[O2A] = pathMismatches(C[O2A], Cr);
        c.paths[O2B] = pathMismatches(C[O2B], Cr);
        c.keys3q = keyMismatches(C[O3Q], C[O2B]);
        c.entered = nEntered;
        for (size_t wi = 0; wi < E.size(); wi++) {
            const bool ent = E[wi] >= 0;
            for (int k = O3; k < kOpts; k++) c.sets[k] += ent != ((C[k].action[wi] & 2) != 0);
            c.missed2b += ent && !(C[O2B].action[wi] & 1);
            c.false2b += !ent && (C[O2B].action[wi] & 1);
            c.boxMissed2a += ent && !(C[O2A].action[wi] & 1);
            c.boxFP2a += !ent && (C[O2A].action[wi] & 1);
        }
        R.chk.add(c);
        // memory: each cache after its update plus that update's temporaries; and the process as a whole
        size_t live = Fc.memory().total() + Cr.memory().total();
        for (int k = O3; k < kOpts; k++) {
            const CacheMemory mm = C[k].memory();
            live += mm.total();
            R.scratchMax[k] = std::max(R.scratchMax[k], r[k].scratchBytes);
            if (mm.total() + r[k].scratchBytes > R.memPeak[k].total() + R.scratchAtPeak[k]) {
                R.memPeak[k] = mm;
                R.scratchAtPeak[k] = r[k].scratchBytes;
            }
        }
        const size_t pb = privateBytes();
        if (pb > mStart && pb - mStart > R.osGrowth) {
            R.osGrowth = pb - mStart;
            R.osComputed = live;
        }
        if (!tables) continue;
        const double nW = double(R.nWalks);
        std::string ordS;
        for (int k : order) ordS += std::string(ordS.empty() ? "" : " ") + kOptShort[k];
        std::snprintf(buf, sizeof buf,
                      "%4zu  %.4f  %.4f | %6.2f%%  %6.2f%%  %6.2f%%  %6.2f%% | %7.4f %7.4f %7.4f %7.4f %7.4f | %5.2fx %5.2fx "
                      "%5.2fx %5.2fx | %6.3f | %s\n",
                      st + 1, ts[st], cut.strip(m).area/cut.fullArea(), 100.0*frac, 100.0*double(r[O3].nTested)/nW,
                      100.0*double(r[O2B].nTested)/nW, 100.0*double(c.false2b)/nW, t[O1], t[O3], t[O3Q], t[O2A], t[O2B],
                      t[O1]/t[O3], t[O1]/t[O3Q], t[O1]/t[O2A], t[O1]/t[O2B], tRef, ordS.c_str());
        rowsT += buf;
        std::snprintf(buf, sizeof buf,
                      "%4zu | %lld %lld %lld %lld | %lld %lld | %lld %lld %lld | %lld %lld %lld %lld | %.1e %.1e %.1e %.1e | "
                      "%.1e %.1e %.1e %.1e | %lld %lld %lld %lld | %lld | %lld %lld | %lld %lld | %lld %lld %lld\n",
                      st + 1, c.walks[O3], c.walks[O3Q], c.walks[O2A], c.walks[O2B], c.paths[O2A], c.paths[O2B],
                      c.sumsVs3[O3Q], c.sumsVs3[O2A], c.sumsVs3[O2B], c.estimate[O3], c.estimate[O3Q], c.estimate[O2A],
                      c.estimate[O2B], c.dH[O3], c.dH[O3Q], c.dH[O2A], c.dH[O2B], c.dG[O3], c.dG[O3Q], c.dG[O2A], c.dG[O2B],
                      c.sets[O3], c.sets[O3Q], c.sets[O2A], c.sets[O2B], c.keys3q, c.missed2b, c.false2b, c.boxMissed2a,
                      c.boxFP2a, c.refSums, c.refRecord, c.fullFallback);
        rowsC += buf;
        std::snprintf(buf, sizeof buf,
                      "%4zu  opt 3  box %.4f  replay+resume+new %.4f                       values %.4f sums %.4f"
                      "               | replayed %lld (prefixes %lld, bbox false positives %lld), new %lld\n"
                      "      opt 3q pop %.4f  replay+resume+new %.4f (keys inside)         values %.4f sums %.4f"
                      "  queue %.4f | replayed %lld (prefixes %lld, false alarms %lld), new %lld, keys %lld\n"
                      "      opt 2a box %.4f  test %.4f  resume %.4f  new %.4f               values %.4f sums %.4f"
                      "               | balls %lld, new %lld\n"
                      "      opt 2b pop %.4f  test %.4f  resume %.4f  new %.4f  keys %.4f  values %.4f sums %.4f"
                      "  queue %.4f | balls %lld, new %lld, keys %lld\n"
                      "      update temporaries (MiB): 3 %.2f, 3q %.2f, 2a %.2f, 2b %.2f\n",
                      st + 1, r[O3].tBox, r[O3].tRewalk, r[O3].tValue, r[O3].tSums, r[O3].stepsReplayed, r[O3].sumK,
                      r[O3].stepsReplayedFalse, r[O3].stepsNew, r[O3Q].tBox, r[O3Q].tRewalk, r[O3Q].tValue, r[O3Q].tSums,
                      r[O3Q].tQueue, r[O3Q].stepsReplayed, r[O3Q].sumK, r[O3Q].stepsReplayedFalse, r[O3Q].stepsNew,
                      r[O3Q].keysComputed, r[O2A].tBox, r[O2A].tTest, r[O2A].tResume, r[O2A].tRewalk, r[O2A].tValue,
                      r[O2A].tSums, r[O2A].ballsTested, r[O2A].stepsNew, r[O2B].tBox, r[O2B].tTest, r[O2B].tResume,
                      r[O2B].tRewalk, r[O2B].tKeys, r[O2B].tValue, r[O2B].tSums, r[O2B].tQueue, r[O2B].ballsTested,
                      r[O2B].stepsNew, r[O2B].keysComputed, mb(r[O3].scratchBytes), mb(r[O3Q].scratchBytes),
                      mb(r[O2A].scratchBytes), mb(r[O2B].scratchBytes));
        rowsP += buf;
    }
    if (!tables) return R;
    out("per step: the share of walks re-walked; the candidates of the bbox test (options 3 and 2a) and of the queue\n"
        "(2b and 3q: the same walks, checked); 2b/3q false alarms; the update wall times; x = option 1 time / option\n"
        "time; the reference solves (validation, untimed); the order in which the updates ran\n");
    out("%4s  %6s  %6s | %7s  %7s  %7s  %7s | %7s %7s %7s %7s %7s | %6s %6s %6s %6s | %6s | %s\n", "step", "t", "strip/",
        "re-", "cand", "cand", "false", "opt 1", "opt 3", "opt 3q", "opt 2a", "opt 2b", "x 3", "x 3q", "x 2a", "x 2b", "ref",
        "order");
    out("%4s  %6s  %6s | %7s  %7s  %7s  %7s | %7s %7s %7s %7s %7s | %6s %6s %6s %6s | %6s |\n", "", "", "full", "walked",
        "bbox", "queue", "alarms", "(s)", "(s)", "(s)", "(s)", "(s)", "", "", "", "", "(s)");
    out.put(rowsT);
    out("checks per step (all counts must be 0; dH, dG = max |H - H_1|, max |grad H - grad H_1| over the points;\n"
        "'> bound' = points beyond the round-off bound of estimateCheck, with H exact for g = 1):\n"
        "step | walks != opt 1 | paths | sums != 3 | > bound    | max |dH| 3 3q 2a 2b     | max |dG| 3 3q 2a 2b     |"
        " re-walked != brute | 3q keys | 2b missed, | 2a box missed, | ref: sums1,\n"
        "     | 3 3q 2a 2b     | 2a 2b | 3q 2a 2b  | 3 3q 2a 2b |                         |                         |"
        " force: 3 3q 2a 2b  | != 2b   | false al.  | box FP         | record, fallback\n");
    out.put(rowsC);
    out("passes of each update (wall time, s):\n");
    out.put(rowsP);
    return R;
}

// "after 3.4 advances" etc.: when an option's setup + n updates (at the run's mean update time) drops below
// another's. dS = its extra setup, dU = what it saves per update.
std::string breakEven(double dS, double dU)
{
    char b[64];
    if (dS <= 0 && dU >= 0) return "from the start";
    if (dS > 0 && dU <= 0) return "never";
    std::snprintf(b, sizeof b, dS > 0 ? "after %.1f advances" : "for the first %.1f advances only", dS/dU);
    return b;
}

void printTotals(Report& out, const RunOut& R, const char* what, bool firstOfMany)
{
    out("%s: totals over %d updates (s)%s; x = option 1 time / option time\n", what, R.opt[O1].updates,
        firstOfMany ? ", first sample" : "");
    out("  option                                  total    x     | filter   test    resume  new      keys    values  "
        "sums    queue  | candidates  re-walked  balls tested  new steps    replayed  (false cand.)      sum k       keys\n");
    for (int k = 0; k < kOpts; k++) {
        const OptTotals& o = R.opt[k];
        out("  %-38s %7.3f %5.2fx  | %7.4f %7.4f %7.4f %8.4f %7.4f %7.4f %7.4f %7.4f | %10lld %10lld %12lld %11lld %11lld "
            "%14lld %10lld %10lld\n", kOptName[k], o.t, o.t > 0 ? R.opt[O1].t/o.t : 0.0, o.tFilter, o.tTest, o.tResume,
            o.tNew, o.tKeys, o.tValue, o.tSums, o.tQueue, o.cand, o.rewalked, o.balls, o.stepsNew, o.stepsReplayed,
            o.replayFalse, o.sumK, o.keys);
    }
    const OptTotals &o3 = R.opt[O3], &oq = R.opt[O3Q], &o2 = R.opt[O2A];
    const double F = double(std::max(1LL, R.fullSteps));
    out("  steps run by the updates, as a share of option 1's %lld:\n"
        "    3  %5.1f%% = new %.1f%% + replayed %.1f%% (prefixes of the re-walked walks %.1f%%, bbox false positives replayed\n"
        "                to the end %.1f%%)\n"
        "    3q %5.1f%% = new %.1f%% + replayed %.1f%% (prefixes %.1f%%, false alarms %.1f%%)\n"
        "    2a, 2b %.1f%% (new steps only; the stored prefixes they skip are %.1f%%)\n", R.fullSteps,
        100.0*double(o3.stepsNew + o3.stepsReplayed)/F, 100.0*double(o3.stepsNew)/F, 100.0*double(o3.stepsReplayed)/F,
        100.0*double(o3.sumK)/F, 100.0*double(o3.replayFalse)/F, 100.0*double(oq.stepsNew + oq.stepsReplayed)/F,
        100.0*double(oq.stepsNew)/F, 100.0*double(oq.stepsReplayed)/F, 100.0*double(oq.sumK)/F,
        100.0*double(oq.replayFalse)/F, 100.0*double(o2.stepsNew)/F, 100.0*double(o2.sumK)/F);
    const int U = std::max(1, R.opt[O1].updates);
    out("  with each option's own initial solve (setup; option 1's is its first full recompute), and when its total\n"
        "  overtakes option 3's and option 1's at this run's mean update time:\n");
    for (int k = 0; k < kOpts; k++) {
        const double tot = R.setup[k] + R.opt[k].t, tot1 = R.setup[O1] + R.opt[O1].t, tot3 = R.setup[O3] + R.opt[O3].t;
        out("    %-2s setup %.3f s + updates %.3f s = %.3f s: %5.2fx vs option 1, %5.2fx vs option 3", kOptShort[k], R.setup[k],
            R.opt[k].t, tot, tot1/tot, tot3/tot);
        if (k != O3) out("; overtakes 3 %s", breakEven(R.setup[k] - R.setup[O3], (R.opt[O3].t - R.opt[k].t)/U).c_str());
        if (k != O1) out("; overtakes 1 %s", breakEven(R.setup[k] - R.setup[O1], (R.opt[O1].t - R.opt[k].t)/U).c_str());
        out("\n");
    }
    out("  reference solves (validation, not part of any option): %.3f s = same-seed full solves %.3f s + the same recording\n"
        "  every step %.3f s (the cost of recording paths in a full solve: %+.0f%%)\n", R.tRef, R.tRefPlain, R.tRefRecord,
        R.tRefPlain > 0 ? 100.0*(R.tRefRecord/R.tRefPlain - 1.0) : 0.0);
}

// per update and option: the median wall time over the samples (the samples run the same updates)
std::vector<std::array<double, kOpts>> stepMedians(const std::vector<RunOut>& reps)
{
    std::vector<std::array<double, kOpts>> t(reps[0].frac.size());
    for (size_t u = 0; u < t.size(); u++)
        for (int k = 0; k < kOpts; k++) {
            std::vector<double> v;
            for (const RunOut& R : reps) v.push_back(R.stepT[k][u]);
            t[u][k] = median(v);
        }
    return t;
}

// the updates at which each option was slower than option ref, and the range of their re-walked shares
std::string slowerThan(const std::vector<RunOut>& reps, int ref)
{
    const std::vector<std::array<double, kOpts>> t = stepMedians(reps);
    const std::vector<double>& f = reps[0].frac;
    std::string s;
    char b[128];
    for (int k = O3; k < kOpts; k++) {
        if (k == ref) continue;
        int n = 0;
        double fMin = 1, fMax = 0;
        for (size_t u = 0; u < t.size(); u++)
            if (t[u][k] > t[u][ref]) {
                n++;
                fMin = std::min(fMin, f[u]);
                fMax = std::max(fMax, f[u]);
            }
        if (n > 1) std::snprintf(b, sizeof b, "%s%s %d of %zu (re-walked %.1f-%.1f%%)", s.empty() ? "" : ", ", kOptShort[k],
                                 n, t.size(), 100.0*fMin, 100.0*fMax);
        else if (n) std::snprintf(b, sizeof b, "%s%s 1 of %zu (re-walked %.1f%%)", s.empty() ? "" : ", ", kOptShort[k],
                                  t.size(), 100.0*fMin);
        else std::snprintf(b, sizeof b, "%s%s 0 of %zu", s.empty() ? "" : ", ", kOptShort[k], t.size());
        s += b;
    }
    return s;
}

// The totals hide a crossover when most walks are re-walked: the update times by the share of walks re-walked
// (the same for every option), and the updates at which an option lost to option 1 or to option 3. Per-update
// times are medians over the samples.
void printBands(Report& out, const std::vector<RunOut>& reps)
{
    const std::vector<std::array<double, kOpts>> t = stepMedians(reps);
    const std::vector<double>& f = reps[0].frac;
    const double lo[3] = {0.75, 0.25, 0.0}, hi[3] = {2.0, 0.75, 0.25};
    const char* name[3] = {">= 75%", "25-75%", " < 25%"};
    out("  by the share of walks re-walked in the update%s:\n", reps.size() > 1 ? " (per-update medians over the samples)" : "");
    for (int b = 0; b < 3; b++) {
        double s[kOpts] = {};
        int n = 0;
        for (size_t u = 0; u < t.size(); u++)
            if (f[u] >= lo[b] && f[u] < hi[b]) {
                n++;
                for (int k = 0; k < kOpts; k++) s[k] += t[u][k];
            }
        if (!n) continue;
        out("    %s: %2d updates: 1 %.3f s | 3 %.3f s %.2fx | 3q %.3f s %.2fx | 2a %.3f s %.2fx | 2b %.3f s %.2fx; vs 3: 3q "
            "%.2fx, 2a %.2fx, 2b %.2fx\n", name[b], n, s[O1], s[O3], s[O1]/s[O3], s[O3Q], s[O1]/s[O3Q], s[O2A],
            s[O1]/s[O2A], s[O2B], s[O1]/s[O2B], s[O3]/s[O3Q], s[O3]/s[O2A], s[O3]/s[O2B]);
    }
    out("  updates slower than option 1: %s\n", slowerThan(reps, O1).c_str());
    out("  updates slower than option 3: %s\n", slowerThan(reps, O3).c_str());
}

void printMemory(Report& out, const RunOut& R, size_t nGradWalks)
{
    const double nW = double(R.nWalks);
    out("  memory (computed from sizes and capacities, without the heap's block headers; 'temporaries' = the buffers\n"
        "  an update allocates and frees, e.g. 2a/2b's resume states; OS = process private bytes):\n");
    out("    option 1 : no per-walk state; per-point sums and query points %.2f MiB\n", mb(R.mem1));
    for (int k = O3; k < kOpts; k++) {
        const CacheMemory &a = R.mem0[k], &p = R.memPeak[k];
        const bool stored = k == O2A || k == O2B;
        const size_t fixedPerWalk = (a.entries + a.pathHeaders + a.keys)/size_t(std::max(1LL, R.nWalks));
        out("    option %-2s: %zu B per walk fixed (+%zu B per gradient-point walk, %zu of them)", kOptShort[k], fixedPerWalk,
            nGradWalks ? a.gradWeights/nGradWalks : size_t(0), nGradWalks);
        if (stored)
            out("; %zu B per stored step (%zu centre and key + 1 flags); stored steps: %lld at the start (%.1f per walk), "
                "%lld at the peak (%.1f per walk, capacity %.1f)", kBytesPerStoredStep, sizeof(StoredStep), a.nSteps,
                double(a.nSteps)/nW, p.nSteps, double(p.nSteps)/nW, double(p.nStepCapacity)/nW);
        out("\n               total %.1f MiB at the start; at the peak %.1f MiB + %.1f MiB update temporaries = %.1f MiB "
            "(largest temporaries %.1f MiB)", mb(a.total()), mb(p.total()), mb(R.scratchAtPeak[k]),
            mb(p.total() + R.scratchAtPeak[k]), mb(R.scratchMax[k]));
        if (stored) out(";\n               stored steps at the peak: %.1f MiB payload, %.1f MiB allocated", mb(p.stepsPayload),
                        mb(p.stepsAllocated));
        if (R.os[k]) out("; OS: +%.1f MiB for configuring it and its initial solve%s", mb(R.os[k]),
                         R.first ? "" : " (may read low: the heap can reuse memory an earlier run freed)");
        out("\n");
    }
    if (R.osGrowth)
        out("    process: the largest private-bytes sample after a step is %.1f MiB above the run's start, against %.1f MiB\n"
            "    computed for the four caches and the two references then (%+.0f%%: heap block headers and heap memory not\n"
            "    returned to the OS; the update temporaries are freed by then)\n", mb(R.osGrowth), mb(R.osComputed),
            R.osComputed ? 100.0*(double(R.osGrowth)/double(R.osComputed) - 1.0) : 0.0);
}

void printChecks(Report& out, const CompareChecks& c, bool ok, int samples)
{
    out("  checks over the run%s: walks != opt 1: 3 %lld, 3q %lld, 2a %lld, 2b %lld; paths: 2a %lld, 2b %lld;\n"
        "  sums != opt 3: 3q %lld, 2a %lld, 2b %lld; H/grad H beyond the round-off bound: 3 %lld, 3q %lld, 2a %lld, 2b %lld;\n"
        "  re-walked set != brute force: 3 %lld, 3q %lld, 2a %lld, 2b %lld; 3q keys/pops/queue != 2b %lld; 2b missed %lld;\n"
        "  2a box missed %lld; initial %lld, ref %lld/%lld, fallbacks %lld; max |dH| 3/3q/2a/2b %.1e/%.1e/%.1e/%.1e,\n"
        "  max |dG| %.1e/%.1e/%.1e/%.1e; entered walks %lld, 2b false alarms %lld, 2a box false positives %lld -> %s\n",
        samples > 1 ? " (all samples)" : "", c.walks[O3], c.walks[O3Q], c.walks[O2A], c.walks[O2B], c.paths[O2A],
        c.paths[O2B], c.sumsVs3[O3Q], c.sumsVs3[O2A], c.sumsVs3[O2B], c.estimate[O3], c.estimate[O3Q], c.estimate[O2A],
        c.estimate[O2B], c.sets[O3], c.sets[O3Q], c.sets[O2A], c.sets[O2B], c.keys3q, c.missed2b, c.boxMissed2a, c.initial,
        c.refSums, c.refRecord, c.fullFallback, c.dH[O3], c.dH[O3Q], c.dH[O2A], c.dH[O2B], c.dG[O3], c.dG[O3Q], c.dG[O2A],
        c.dG[O2B], c.entered, c.false2b, c.boxFP2a, ok ? "PASS" : "FAIL");
}

// "1.172 s" or, over several samples, the median and range "1.172 s [1.150, 1.201]"; a ratio ("1.17x") also gets
// " tie" when its range contains 1
std::string stat(const std::vector<RunOut>& reps, const std::function<double(const RunOut&)>& f, bool ratio)
{
    std::vector<double> v;
    for (const RunOut& R : reps) v.push_back(f(R));
    std::sort(v.begin(), v.end());
    char b[96];
    if (v.size() == 1) std::snprintf(b, sizeof b, ratio ? "%.2fx" : "%.3f s", v[0]);
    else if (ratio)
        std::snprintf(b, sizeof b, "%.2fx [%.2f, %.2f]%s", median(v), v.front(), v.back(),
                      v.front() <= 1.0 && v.back() >= 1.0 ? " tie" : "");
    else std::snprintf(b, sizeof b, "%.3f s [%.3f, %.3f]", median(v), v.front(), v.back());
    return b;
}

// a run's update times and ratios, then the same with setup, as medians [min, max] over the samples
std::string timeLines(const std::vector<RunOut>& reps)
{
    std::string L;
    for (int withSetup = 0; withSetup < 2; withSetup++) {
        auto T = [withSetup](const RunOut& R, int k) { return R.opt[k].t + (withSetup ? R.setup[k] : 0.0); };
        L += withSetup ? "\n      setup + updates: " : "updates: ";
        for (int k = 0; k < kOpts; k++) {
            L += std::string(k ? " | " : "opt ") + kOptShort[k] + " " + stat(reps, [&](const RunOut& R) { return T(R, k); }, false);
            if (k) L += " " + stat(reps, [&](const RunOut& R) { return T(R, O1)/T(R, k); }, true);
        }
        L += "; vs 3:";
        for (int k = O3Q; k < kOpts; k++)
            L += std::string(k == O3Q ? " " : ", ") + kOptShort[k] + " " +
                 stat(reps, [&](const RunOut& R) { return T(R, O3)/T(R, k); }, true);
    }
    return L;
}

std::string summaryLine(const char* what, const std::vector<RunOut>& reps, bool ok)
{
    const RunOut& R = reps[0];
    std::string L = std::string(what) + ": " + timeLines(reps) + "\n      updates slower than option 1: " +
                    slowerThan(reps, O1) + "; checks " + (ok ? "PASS" : "FAIL");
    char b[768];
    std::snprintf(b, sizeof b,
                  "\n      memory at the peak, allocated incl. the update's temporaries [OS for the initial solve%s]: opt 1 "
                  "%.2f MiB, 3 %.1f [%.1f], 3q %.1f [%.1f], 2a %.1f [%.1f], 2b %.1f [%.1f] MiB;\n      stored-step payload "
                  "at the start %.1f MiB (%lld steps x %zu B)",
                  R.first ? "" : "; may read low after the process's first run, as the heap reuses memory freed by earlier runs",
                  mb(R.mem1), mb(R.memPeak[O3].total() + R.scratchAtPeak[O3]), mb(R.os[O3]),
                  mb(R.memPeak[O3Q].total() + R.scratchAtPeak[O3Q]), mb(R.os[O3Q]),
                  mb(R.memPeak[O2A].total() + R.scratchAtPeak[O2A]), mb(R.os[O2A]),
                  mb(R.memPeak[O2B].total() + R.scratchAtPeak[O2B]), mb(R.os[O2B]), mb(R.mem0[O2A].stepsPayload),
                  R.mem0[O2A].nSteps, kBytesPerStoredStep);
    return L + b;
}

} // namespace

// =====================================================================================
bool runCompareOptions(const GeometryConfig& cfg, const Settings& s, const std::string& outPath, bool quick, int threads,
                       int repeats)
{
    auto tStart = std::chrono::steady_clock::now();
    repeats = std::max(1, repeats);
    std::unique_ptr<tbb::global_control> gc;
    if (threads > 0) gc = std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, size_t(threads));
    const int nThreads = int(tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism));
    Tet tet = referenceTet();
    Report out;
    std::time_t now = std::time(nullptr);
    std::tm tmNow{};
#ifdef _WIN32
    localtime_s(&tmNow, &now);
#else
    localtime_r(&now, &tmNow);
#endif
    char date[64];
    std::strftime(date, sizeof date, "%Y-%m-%d %H:%M:%S", &tmNow);
    const int N = quick ? 32 : s.nWalks, K = quick ? 6 : s.k;
    SliceParams sp = s.slice;
    if (quick) sp.res = 16;
    const GType gTypes[2] = {GType::One, GType::SqrtFront};
    out("Exact walk reuse under an advancing cut: options 1, 3, 2a, 2b and the control 3q compared (--compare-options)%s\n",
        quick ? " -- QUICK reduced-size run" : "");
    out("date: %s; hardware threads: %u; TBB threads used: %d%s\n", date, std::thread::hardware_concurrency(), nThreads,
        threads > 0 ? " (--threads)" : "");
    const std::string cpu = cpuInfo();
    if (!cpu.empty()) out("CPU: %s\n", cpu.c_str());
    if (repeats == 1)
        out("samples: 1 per run (--repeat R for more): the run-to-run spread is not measured, so differences of a few\n"
            "percent between options are not resolved\n");
    else
        out("samples: %d per run (--repeat), each a complete run with every check: totals are medians [min, max] over the\n"
            "samples, and a ratio whose range contains 1 is marked 'tie'; the tables show the first sample\n", repeats);
    out("tet faces: %s\n", s.walk.dirichletFaces ? "DIRICHLET, H = 0 (experiment --dirichlet-faces: a different PDE)"
                                                 : "Neumann (reflecting), the real problem");
    out("Options (all exact: the same estimator, walk for walk):\n"
        "  1  frozen seeds, full recompute: a same-seed solve keeping only the sums (solveSums), no per-walk state\n"
        "  3  seed replay (Saeed's layout, WalkCache::advance(Incremental), its production path): bbox test, replay\n"
        "     from the seed under the old cut to the first ball the strip enters, resume there\n"
        "  3q control, not one of the options asked about: option 3's replay and resume, but only for the walks popped\n"
        "     from 2b's queue; the walks accumulate their keys as they run, so nothing per step is stored. It separates\n"
        "     what 2b gains by its tighter filter from what it gains by storing every step.\n"
        "  2a every step stored (centre, key, flags): the same bbox test, then the stored balls are tested exactly\n"
        "     (stripDist2(c_k) < key2_k, runWalk's own stop test) and the walk resumes from the stored step: no replay\n"
        "  2b as 2a, but the candidates are popped from a queue sorted by each walk's entry front min_k s*_k (the\n"
        "     front position at which the known future cut surface first enters one of its balls, lowered by margins)\n"
        "  All run the one step function runWalk; value updates (g type (b)) and sums use the same code, in walk order.\n");
    out("Timing: steady_clock wall time of each option's update call, same thread count, back to back, in a random\n"
        "order per step (fixed seed; the last column of the step table). Option 1 is a full recompute. Not timed: the\n"
        "brute-force ball test before the updates (on the reference's recorded paths, so it warms no option's data) and\n"
        "the reference solves after them (a same-seed full solve with per-walk entries, which is option 1 with its walks\n"
        "kept, and one that records every path). 'updates' totals leave out each option's initial solve; 'setup +\n"
        "updates' include it (option 1's: its first full recompute). The ratios depend on the machine and on the thread\n"
        "count (--threads): compare runs at several thread counts before generalizing them.\n");
    out("Checks at every step: every option's per-walk entries (box, landing point, value, steps, side, status, w)\n"
        "bitwise equal to option 1's same-seed walks; 2a/2b's stored paths bitwise equal to a fresh recorded solve;\n"
        "sums of 3q/2a/2b bitwise equal to option 3's (same code, same walks); H and grad H of every option within a\n"
        "worst-case round-off bound of option 1's (H exactly for g = 1; the bound is derived in estimateCheck, and a\n"
        "wrong value update or sum shows up far above it); the re-walked set of every option equal to the brute-force\n"
        "ball test; 2b soundness: no walk whose stored ball the strip enters is left unpopped; 3q's keys, pops and\n"
        "queue bitwise equal to 2b's.\n");
    bool allOk = true;
    std::vector<std::string> summary;
    char what[128];
    int runId = 0;
    // one configuration: 'repeats' complete samples (every check runs in each), the first one printed in full
    auto runConfig = [&](const char* label, const QueryPoints& q, GType g, double t0, const std::vector<double>& ts,
                         size_t nGrad) {
        std::vector<RunOut> reps;
        for (int rep = 0; rep < repeats; rep++) {
            if (repeats > 1) out("[sample %d of %d]\n", rep + 1, repeats);
            reps.push_back(compareRun(out, tet, cfg, q, N, g, s.walk, t0, ts, runId == 0 && rep == 0,
                                      mixSeed(mixSeed(0x6f72646572ULL, uint64_t(runId)), uint64_t(rep)), rep == 0));
        }
        runId++;
        printTotals(out, reps[0], gShortName(g), repeats > 1);
        if (repeats > 1) out("  over the %d samples, median [min, max]: %s\n", repeats, timeLines(reps).c_str());
        printBands(out, reps);
        printMemory(out, reps[0], nGrad*size_t(N));
        CompareChecks all;
        for (const RunOut& R : reps) all.add(R.chk);
        const bool ok = all.bad() == 0;
        allOk = allOk && ok;
        printChecks(out, all, ok, repeats);
        summary.push_back(summaryLine(label, reps, ok));
    };

    // ================================================================ (i) lattice sweep
    {
        Lattice lat = buildLattice(tet, Cut{}, K, 0.0);
        size_t nGrad = 0;
        for (char c : lat.interior) nGrad += c;
        const QueryPoints q = makePoints(lat.pts, &lat.interior, mixSeed(s.seed, 600));
        std::vector<double> ts;
        for (int st = 1; st <= 16; st++) ts.push_back(0.2 + 0.8*double(st)/16.0);
        out("\n=== (i) lattice sweep t = 0.20 -> 1.00 in 16 steps: lattice k = %d, %zu points (%zu with grad H), N = %d, "
            "seed tag 600 (the validation sweep's) ===\n", K, q.size(), nGrad, N);
        for (int gi = 0; gi < 2; gi++) {
            out("\n-- %s --\n", gShortName(gTypes[gi]));
            std::snprintf(what, sizeof what, "(i) sweep %s", gTypes[gi] == GType::One ? "(a)" : "(b)");
            runConfig(what, q, gTypes[gi], 0.2, ts, nGrad);
        }
    }

    // ================================================================ (ii) check 3's schedule
    {
        Lattice lat = buildLattice(tet, Cut{}, K, 0.0);
        StripCut c0(tet, cfg);
        Slice sl = buildSlice(tet, c0.frame(), sp);
        size_t nGrad = 0;
        for (char c : lat.interior) nGrad += c;
        QueryPoints q;   // check 3's seeds
        q.add(lat.pts, &lat.interior, mixSeed(mixSeed(s.seed, 800), 1));
        q.add(sl.verts, nullptr, mixSeed(mixSeed(s.seed, 800), 2));
        out("\n=== (ii) check 3's schedule: full solve at t = 0.5, then 4 steps of dt = 1/8, 1/16, 1/32, 1/64 (one run per dt);"
            "\n    default scene: lattice k = %d + slice (res %d): %zu points (%zu with grad H), N = %d, check 3's seeds "
            "(tag 800) ===\n", K, sp.res, q.size(), nGrad, N);
        for (int gi = 0; gi < 2; gi++)
            for (int den : {8, 16, 32, 64}) {
                std::vector<double> ts;
                for (int st = 1; st <= 4; st++) ts.push_back(0.5 + double(st)/double(den));
                out("\n-- %s, dt = 1/%d --\n", gShortName(gTypes[gi]), den);
                std::snprintf(what, sizeof what, "(ii) %s dt = 1/%-2d", gTypes[gi] == GType::One ? "(a)" : "(b)", den);
                runConfig(what, q, gTypes[gi], 0.5, ts, nGrad);
            }
    }

    out("\n=== Summary (wall time summed over each run's updates; x = option 1 time / option time; 'vs 3' = option 3 time /"
        "\n    option time; %s) ===\n", repeats > 1 ? "medians [min, max] over the samples, 'tie' where a range contains 1"
                                                    : "one sample per run");
    for (const std::string& l : summary) out("%s\n", l.c_str());
    out("\nall equivalence and soundness checks: %s\n", allOk ? "PASS (0 mismatches)" : "SOME CHECK FAILED");
    out("compare-options wall time: %.1f s\n", secondsSince(tStart));
    std::ofstream f(outPath, std::ios::binary);
    f << out.text;
    f.close();
    if (!f) {
        std::fprintf(stderr, "could not write %s\n", outPath.c_str());
        return false;
    }
    std::printf("written: %s\n", outPath.c_str());
    return allOk;
}
