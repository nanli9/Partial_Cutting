// --gate coupling: can one tet's WoSt solve be glued to a neighbour through the shared face?
// See coupling_gate.h. Everything printed goes to stdout; a full (not quick) run also rewrites its
// sections of gate_report.md.
#include "coupling_gate.h"

#include "walk_cache.h"

#include "oneapi/tbb/parallel_for.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace {

constexpr double kSqrt2 = 1.41421356237309504880;
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kTwoGMax = 2.0;   // every error and noise level is divided by 2 g_max = 2, in both g types

double secondsSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// printf to stdout and into a string; a QUICK run starts every line with "QUICK "
class Report {
public:
    explicit Report(bool quick) : quick_(quick) {}
    void operator()(const char* fmt, ...)
    {
        char buf[8192];
        va_list ap;
        va_start(ap, fmt);
        std::vsnprintf(buf, sizeof buf, fmt, ap);
        va_end(ap);
        put(buf);
    }
    void put(const std::string& s)
    {
        std::string o;
        for (char ch : s) {
            if (lineStart_ && quick_) o += "QUICK ";
            o += ch;
            lineStart_ = ch == '\n';
        }
        std::fputs(o.c_str(), stdout);
        std::fflush(stdout);
        text += o;
    }
    std::string text;

private:
    bool quick_ = false, lineStart_ = true;
};

std::string fmt3(const Vec3& v, int prec = 4)
{
    char buf[128];
    std::snprintf(buf, sizeof buf, "(%.*f, %.*f, %.*f)", prec, v.x(), prec, v.y(), prec, v.z());
    return buf;
}

std::string fmtd(const char* f, double v)
{
    char buf[64];
    std::snprintf(buf, sizeof buf, f, v);
    return buf;
}

// a markdown pipe table (gate_report.md)
struct MdTable {
    std::vector<std::string> head;
    std::vector<std::vector<std::string>> rows;
    std::string str() const
    {
        std::string s = "|";
        for (const auto& h : head) s += " " + h + " |";
        s += "\n|";
        for (size_t i = 0; i < head.size(); i++) s += "---|";
        s += "\n";
        for (const auto& r : rows) {
            s += "|";
            for (const auto& c : r) s += " " + c + " |";
            s += "\n";
        }
        return s;
    }
};

// =====================================================================================
// Geometry (physical coordinates X)
// =====================================================================================
// Global vertex ids: v0, v1, v2 = the corners of F, v3 = A's apex, v4 = B's apex. The degrees of freedom
// on F are indexed by these ids, not by a tet's local index.
const std::array<Vec3, 5> kGlobal = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(0, 0, -1)};
constexpr int kFaceF = 3;   // in A and in B, F is the face opposite local vertex 3 (xi_z = 0)

// One tet of the pair in its own xi frame: X = X0 + J xi, xi on the reference tet.
struct TetFrame {
    const char* name = "";
    std::array<int, 4> global{};   // local vertex -> global id
    std::array<Vec3, 4> X;         // physical vertices, local order
    Eigen::Matrix3d J, Jinv, JinvT;
    double detJ = 0;
    Vec3 toX(const Vec3& xi) const { return X[0] + J*xi; }
    Vec3 toXi(const Vec3& x) const { return Jinv*(x - X[0]); }
    int localOf(int g) const
    {
        for (int k = 0; k < 4; k++)
            if (global[k] == g) return k;
        return -1;
    }
};

// the frame of the tet with global vertex ids global (local order) on the vertex table verts
TetFrame frameOn(const char* name, const std::array<int, 4>& global, const Vec3* verts)
{
    TetFrame T;
    T.name = name;
    T.global = global;
    for (int k = 0; k < 4; k++) T.X[k] = verts[global[k]];
    T.J.col(0) = T.X[1] - T.X[0];
    T.J.col(1) = T.X[2] - T.X[0];
    T.J.col(2) = T.X[3] - T.X[0];
    T.detJ = T.J.determinant();
    T.Jinv = T.J.inverse();
    T.JinvT = T.Jinv.transpose();
    return T;
}
TetFrame makeFrame(const char* name, const std::array<int, 4>& global) { return frameOn(name, global, kGlobal.data()); }
// A = the reference tet (xi = X); B = A mirrored across z = 0, ordered so that det J_B = +1:
// J_B = [[0,1,0],[1,0,0],[0,0,-1]], xi = (y, x, -z). Local -> global: A (0,1,2,3), B (0,2,1,4).
TetFrame frameA() { return makeFrame("A", {0, 1, 2, 3}); }
TetFrame frameB() { return makeFrame("B", {0, 2, 1, 4}); }

// The union U = A u B = {x, y >= 0, x + y + |z| <= 1} is itself a tetrahedron: conv{v1, v2, v3, v4}.
Tet unionTet() { return Tet({kGlobal[1], kGlobal[2], kGlobal[3], kGlobal[4]}); }

// The crack in physical terms: plane n.X = d (n unit, H = +g on n.X > d), front direction +y, crack region
// y < yFront. With these planes the in-plane front direction u is exactly +y, so the front is the line
// {n.X = d, y = yFront}.
struct CrackSpec {
    bool on = false;
    Vec3 n = Vec3::Zero();
    double d = 0;
    Vec3 f = Vec3(0, 1, 0);
    double yFront = 0;
    std::string name;
};
CrackSpec noCrack()
{
    CrackSpec c;
    c.name = "no crack";
    return c;
}
CrackSpec crackCheck2(double zc)   // plane z = z_c, +g above, y_f = (1 - z_c)/2
{
    CrackSpec c;
    c.on = true;
    c.n = Vec3(0, 0, 1);
    c.d = zc;
    c.yFront = (1.0 - zc)/2.0;
    c.name = "z = " + fmtd("%.2f", zc) + " (+g above), front y = " + fmtd("%.3f", c.yFront) + ", crack y < front";
    return c;
}
CrackSpec crackCheck3(double yf)   // plane x = 0.3, +g for x > 0.3
{
    CrackSpec c;
    c.on = true;
    c.n = Vec3(1, 0, 0);
    c.d = 0.3;
    c.yFront = yf;
    c.name = "x = 0.30 (+g for x > 0.3), front y = " + fmtd("%.2f", yf) + ", crack y < front";
    return c;
}

// The crack in one tet's xi frame: cfg (restVerts = the tet's physical vertices, plane and front in X)
// pulled back by buildCut on the reference tet, at the t that puts the front at y = yFront in X.
struct XiCrack {
    GeometryConfig cfg;
    double t = 0;
    Cut cut;            // in xi; empty without a crack or if the plane misses the tet
    StripCut strips;    // the same cut for the own loop: one strip [sMin, s(t)]
    Cut physFull;       // geometry only: the physical section (buildCutOnTet on the physical tet, t = 1)
    double frontErr = 0;   // max |y - yFront| over the xi front vertices mapped to X
    int nFront = 0;
};

// the vertices of a cut polygon on its front line u.x = s
std::vector<Vec3> frontVertices(const Cut& c)
{
    std::vector<Vec3> v;
    for (const Vec3& p : c.polygon)
        if (std::fabs(p.dot(c.u) - c.s) <= 1e-9*c.L) v.push_back(p);
    return v;
}

XiCrack xiCrack(const TetFrame& T, const CrackSpec& cs)
{
    XiCrack x;
    if (!cs.on) return x;
    x.cfg.restVerts = T.X;
    x.cfg.planeNormal = cs.n;
    x.cfg.planeOffset = cs.d;
    x.cfg.frontDir = cs.f;
    // t from the physical section's extent along u = +y (Cut's t is relative to that extent)
    x.physFull = buildCutOnTet(Tet(T.X), cs.n, cs.d, cs.f, 1.0);
    if (x.physFull.fullPolygon.size() < 3) return x;   // the plane misses this tet
    x.t = (cs.yFront - x.physFull.sMin)/(x.physFull.sMax - x.physFull.sMin);
    x.cut = buildCut(referenceTet(), x.cfg, x.t);
    x.strips = StripCut(referenceTet(), x.cfg);
    x.strips.addStrip(x.strips.sOf(x.t));
    for (const Vec3& p : frontVertices(x.cut)) {
        x.frontErr = std::max(x.frontErr, std::fabs(T.toX(p).y() - cs.yFront));
        x.nFront++;
    }
    return x;
}

// The crack on U in X (the reference of checks 2 and 3; geometry only in this stage). L is set to
// sqrt(2), the L of the per-tet xi frames, so that g agrees pointwise (U's own longest edge is 2).
Cut unionCut(const CrackSpec& cs, double* tOut = nullptr)
{
    const Tet U = unionTet();
    if (!cs.on) return Cut{};
    const Cut full = buildCutOnTet(U, cs.n, cs.d, cs.f, 1.0);
    const double t = (cs.yFront - full.sMin)/(full.sMax - full.sMin);
    if (tOut) *tOut = t;
    Cut c = buildCutOnTet(U, cs.n, cs.d, cs.f, t);
    c.L = kSqrt2;
    return c;
}

// outward plane of face o (opposite local vertex o) of the physical tet, unnormalized normal scaled so
// that its largest component is 1: "(a, b, c).X = e"
std::string facePlane(const std::array<Vec3, 4>& X, int o)
{
    std::array<int, 3> f{};
    int m = 0;
    for (int i = 0; i < 4; i++)
        if (i != o) f[m++] = i;
    Vec3 nrm = (X[f[1]] - X[f[0]]).cross(X[f[2]] - X[f[0]]);
    if (nrm.dot(X[o] - X[f[0]]) > 0) nrm = -nrm;
    nrm /= nrm.cwiseAbs().maxCoeff();
    char buf[96];
    std::snprintf(buf, sizeof buf, "(%g,%g,%g).X = %g", nrm.x() == 0 ? 0.0 : nrm.x(), nrm.y() == 0 ? 0.0 : nrm.y(),
                  nrm.z() == 0 ? 0.0 : nrm.z(), nrm.dot(X[f[0]]) == 0 ? 0.0 : nrm.dot(X[f[0]]));
    return buf;
}

std::string faceCorners(const TetFrame& T, int o)
{
    std::string s;
    for (int i = 0; i < 4; i++)
        if (i != o) s += (s.empty() ? "v" : " v") + std::to_string(T.global[i]);
    return s;
}

// markF: face kFaceF is F (A and B; not U, whose faces are all exterior)
std::string faceFlags(const TetFrame& T, const std::array<FaceBC, 4>& bc, const std::string& dirichletData, bool markF = true)
{
    std::string s;
    for (int o = 0; o < 4; o++) {
        s += std::string(o ? "; " : "") + "f" + std::to_string(o) + " " + facePlane(T.X, o) + " [" + faceCorners(T, o) +
             "]" + (markF && o == kFaceF ? " = F" : "") + ": " +
             (bc[o] == FaceBC::Dirichlet ? "DIRICHLET " + dirichletData : "Neumann");
    }
    return s;
}

// =====================================================================================
// Hats (physical coordinates)
// =====================================================================================
// hat of global corner g of F (g = 0, 1, 2): 1 at v_g, 0 at the other two corners, linear on F (z = 0)
double hatF(int g, const Vec3& X)
{
    return g == 0 ? 1.0 - X.x() - X.y() : g == 1 ? X.x() : X.y();
}
// barycentric coordinate of local vertex k of a physical tet (the hat of that vertex on every face)
double lambdaOf(const TetFrame& T, int k, const Vec3& X)
{
    const Vec3 xi = T.toXi(X);
    return k == 0 ? 1.0 - xi.sum() : xi[k - 1];
}

// =====================================================================================
// Quadrature point sets
// =====================================================================================
double radicalInverse(uint64_t i, uint64_t base)
{
    double f = 1.0, r = 0.0;
    while (i) {
        f /= double(base);
        r += f*double(i%base);
        i /= base;
    }
    return r;
}

// unit cube -> uniform on the simplex {a, b, c >= 0, a + b + c <= 1}: s = a+b+c has density 3s^2, then
// (a+b)/s has density 2t, then a/(a+b) is uniform (a smooth bijection, so QMC points stay well spread)
Vec3 cubeToSimplex(double u1, double u2, double u3)
{
    const double s = std::cbrt(u1), t = std::sqrt(u2);
    return Vec3(s*t*u3, s*t*(1.0 - u3), s*(1.0 - t));
}

struct Leaf {
    std::array<Vec3, 4> v;
    int depth = 0;
};
double leafVolume(const Leaf& l)
{
    return std::fabs((l.v[1] - l.v[0]).cross(l.v[2] - l.v[0]).dot(l.v[3] - l.v[0]))/6.0;
}
// Bey's red refinement: 8 children of exactly 1/8 of the volume
std::array<Leaf, 8> beyChildren(const Leaf& L)
{
    const Vec3 &x0 = L.v[0], &x1 = L.v[1], &x2 = L.v[2], &x3 = L.v[3];
    const Vec3 x01 = 0.5*(x0 + x1), x02 = 0.5*(x0 + x2), x03 = 0.5*(x0 + x3), x12 = 0.5*(x1 + x2),
               x13 = 0.5*(x1 + x3), x23 = 0.5*(x2 + x3);
    const int d = L.depth + 1;
    return {Leaf{{x0, x01, x02, x03}, d}, Leaf{{x01, x1, x12, x13}, d}, Leaf{{x02, x12, x2, x23}, d},
            Leaf{{x03, x13, x23, x3}, d}, Leaf{{x01, x02, x03, x13}, d}, Leaf{{x01, x02, x12, x13}, d},
            Leaf{{x02, x03, x13, x23}, d}, Leaf{{x02, x12, x13, x23}, d}};
}

double pointSegmentDist(const Vec3& p, const Vec3& a, const Vec3& b)
{
    const Vec3 ab = b - a;
    const double l2 = ab.squaredNorm();
    const double s = l2 > 0 ? std::clamp((p - a).dot(ab)/l2, 0.0, 1.0) : 0.0;
    return (p - (a + s*ab)).norm();
}

// One fixed, deterministic point set in one tet (built in xi on the reference tet). The tet is cut into
// leaves by Bey's refinement: 64 leaves of depth 2, and near the front (leaf centroid closer to the front
// segment than half its longest edge) depth 3 and 4. Leaf l gets m_l points, uniform in it (Halton bases 2,
// 3, 5, consecutive indices, through cubeToSimplex), so the density is q = m_l/(M V_l) on leaf l, known
// exactly, and every point of leaf l has the weight w = V_l/m_l (a stratified rule). m_l = round(c
// 8^(-(depth-2)/2)), at least mMin: the density grows by sqrt(8) per level, 8x at depth 4; c is chosen for
// about 'target' points.
// Without a front every leaf has depth 2 and the same m: the density is uniform. Then the points within eps
// of the crack are dropped: the rule integrates over T minus that eps-tube.
struct PointSet {
    std::vector<Vec3> xi, X;     // kept points, solve frame and physical
    std::vector<double> w;       // their weights (physical volume)
    std::vector<Vec3> allXi;     // every generated point (before the eps drop) and its weight
    std::vector<double> allW;
    int nLeaves = 0, leavesByDepth[5] = {0, 0, 0, 0, 0}, pointsByDepth[5] = {0, 0, 0, 0, 0};
    int nDropped = 0;
    double eps = 0;
    size_t size() const { return xi.size(); }
};

PointSet volumePoints(const TetFrame& T, const XiCrack& xc, int target, double eps, int mMin)
{
    PointSet ps;
    ps.eps = eps;
    const Tet ref = referenceTet();
    std::vector<Vec3> fr = frontVertices(xc.cut);
    const bool front = fr.size() >= 2;
    std::vector<Leaf> leaves = {Leaf{{ref.v[0], ref.v[1], ref.v[2], ref.v[3]}, 0}};
    auto refine = [&](auto pred) {
        std::vector<Leaf> out;
        for (const Leaf& l : leaves) {
            if (!pred(l)) {
                out.push_back(l);
                continue;
            }
            for (const Leaf& c : beyChildren(l)) out.push_back(c);
        }
        leaves.swap(out);
    };
    refine([](const Leaf&) { return true; });
    refine([](const Leaf&) { return true; });
    for (int level = 2; level < 4 && front; level++)
        refine([&](const Leaf& l) {
            if (l.depth != level) return false;
            const Vec3 c = 0.25*(l.v[0] + l.v[1] + l.v[2] + l.v[3]);
            double diam = 0;
            for (int i = 0; i < 4; i++)
                for (int j = i + 1; j < 4; j++) diam = std::max(diam, (l.v[i] - l.v[j]).norm());
            return pointSegmentDist(c, fr.front(), fr.back()) < 0.5*diam;
        });
    auto counts = [&](double c) {
        std::vector<int> m(leaves.size());
        for (size_t i = 0; i < leaves.size(); i++)
            m[i] = std::max(mMin, int(std::lround(c*std::pow(8.0, -0.5*(leaves[i].depth - 2)))));
        return m;
    };
    auto total = [](const std::vector<int>& m) {
        long long s = 0;
        for (int v : m) s += v;
        return s;
    };
    double lo = 1, hi = double(target);
    for (int it = 0; it < 60; it++) {   // bisection on c: the total is monotone in c
        const double mid = 0.5*(lo + hi);
        (total(counts(mid)) < target ? lo : hi) = mid;
    }
    const std::vector<int> m = counts(hi);
    uint64_t idx = 1;   // Halton index 0 is the corner (0,0,0)
    ps.nLeaves = int(leaves.size());
    for (size_t li = 0; li < leaves.size(); li++) {
        const Leaf& l = leaves[li];
        const double wl = leafVolume(l)/double(m[li])*std::fabs(T.detJ);
        ps.leavesByDepth[l.depth]++;
        ps.pointsByDepth[l.depth] += m[li];
        for (int k = 0; k < m[li]; k++, idx++) {
            const Vec3 b = cubeToSimplex(radicalInverse(idx, 2), radicalInverse(idx, 3), radicalInverse(idx, 5));
            const Vec3 p = l.v[0] + b.x()*(l.v[1] - l.v[0]) + b.y()*(l.v[2] - l.v[0]) + b.z()*(l.v[3] - l.v[0]);
            ps.allXi.push_back(p);
            ps.allW.push_back(wl);
            // J is orthogonal for A and B, so the xi distance to the crack is the physical one
            if (!xc.cut.empty() && xc.cut.distToCut(p) < eps) {
                ps.nDropped++;
                continue;
            }
            ps.xi.push_back(p);
            ps.X.push_back(T.toX(p));
            ps.w.push_back(wl);
        }
    }
    return ps;
}

// The fixed point set on F (physical, z = 0) for the trace-fit diagnostics of checks 2 and 3: 64 equal
// sub-triangles (3 midpoint subdivisions), 16 points each (Halton bases 2, 3), weight area/1024.
struct FacePointSet {
    std::vector<Vec3> X;
    std::vector<double> w;
};
FacePointSet facePoints(int perTri)
{
    std::vector<std::array<Vec3, 3>> tris = {{kGlobal[0], kGlobal[1], kGlobal[2]}};
    for (int lev = 0; lev < 3; lev++) {
        std::vector<std::array<Vec3, 3>> out;
        for (const auto& t : tris) {
            const Vec3 a = 0.5*(t[0] + t[1]), b = 0.5*(t[1] + t[2]), c = 0.5*(t[2] + t[0]);
            out.push_back({t[0], a, c});
            out.push_back({a, t[1], b});
            out.push_back({c, b, t[2]});
            out.push_back({a, b, c});
        }
        tris.swap(out);
    }
    FacePointSet fs;
    uint64_t idx = 1;
    for (const auto& t : tris) {
        const double area = 0.5*(t[1] - t[0]).cross(t[2] - t[0]).norm();
        for (int k = 0; k < perTri; k++, idx++) {
            const double r = std::sqrt(radicalInverse(idx, 2)), s = radicalInverse(idx, 3);
            fs.X.push_back((1.0 - r)*t[0] + r*(1.0 - s)*t[1] + r*s*t[2]);
            fs.w.push_back(area/double(perTri));
        }
    }
    return fs;
}

// test functions for the quadrature checks (physical X): constants and quadratics
const int kNTest = 10;
const char* kTestName[kNTest] = {"1", "x", "y", "z", "x^2", "y^2", "z^2", "xy", "yz", "xz"};
double testFn(int k, const Vec3& X)
{
    switch (k) {
    case 0: return 1.0;
    case 1: return X.x();
    case 2: return X.y();
    case 3: return X.z();
    case 4: return X.x()*X.x();
    case 5: return X.y()*X.y();
    case 6: return X.z()*X.z();
    case 7: return X.x()*X.y();
    case 8: return X.y()*X.z();
    default: return X.x()*X.z();
    }
}
// exact integral over the physical tet: the 4-point rule (exact for quadratics)
double exactTetIntegral(const TetFrame& T, int k)
{
    const double a = 0.5854101966249685, b = 0.1381966011250105;
    const double V = std::fabs(T.detJ)/6.0;
    double s = 0;
    for (int i = 0; i < 4; i++) {
        Vec3 p = Vec3::Zero();
        for (int j = 0; j < 4; j++) p += (i == j ? a : b)*T.X[j];
        s += testFn(k, p);
    }
    return V*s/4.0;
}

// Volume of the eps-tube (T within eps of the crack) and the integrals of the test functions over it, by
// QMC in a box around the cut (cut frame u, w, n; Halton with the bases given). nPts points.
void tubeIntegrals(const TetFrame& T, const Cut& cut, double eps, int nPts, const uint64_t bases[3], double out[kNTest])
{
    std::fill(out, out + kNTest, 0.0);
    if (cut.empty()) return;
    double uLo = kInf, uHi = -kInf, wLo = kInf, wHi = -kInf;
    for (const Vec3& p : cut.polygon) {
        uLo = std::min(uLo, p.dot(cut.u));
        uHi = std::max(uHi, p.dot(cut.u));
        wLo = std::min(wLo, p.dot(cut.w));
        wHi = std::max(wHi, p.dot(cut.w));
    }
    uLo -= eps, uHi += eps, wLo -= eps, wHi += eps;
    const double boxVol = (uHi - uLo)*(wHi - wLo)*2.0*eps;
    const Tet ref = referenceTet();
    const int nChunks = 64;
    std::vector<std::array<double, kNTest>> part(nChunks);
    tbb::parallel_for(0, nChunks, [&](int c) {
        std::array<double, kNTest> s{};
        for (int i = int((long long)nPts*c/nChunks); i < int((long long)nPts*(c + 1)/nChunks); i++) {
            const uint64_t idx = uint64_t(i) + 1;
            const double a = uLo + (uHi - uLo)*radicalInverse(idx, bases[0]);
            const double b = wLo + (wHi - wLo)*radicalInverse(idx, bases[1]);
            const double h = eps*(2.0*radicalInverse(idx, bases[2]) - 1.0);
            const Vec3 x = cut.d*cut.n + a*cut.u + b*cut.w + h*cut.n;
            if (ref.bary(x).minCoeff() < 0 || cut.distToCut(x) >= eps) continue;
            const Vec3 X = T.toX(x);
            for (int k = 0; k < kNTest; k++) s[k] += testFn(k, X);
        }
        part[c] = s;
    });
    for (int c = 0; c < nChunks; c++)
        for (int k = 0; k < kNTest; k++) out[k] += part[c][k];
    for (int k = 0; k < kNTest; k++) out[k] *= boxVol*std::fabs(T.detJ)/double(nPts);
}

// Prints the quadrature checks of one point set; returns the largest relative error of the kept rule. tube(bases, out):
// the eps-tube's volume and test-function integrals by QMC with the Halton bases given; hasTube: the set has one.
using TubeFn = std::function<void(const uint64_t bases[3], double out[kNTest])>;
double quadratureCheckWith(Report& out, const TetFrame& T, bool hasTube, const PointSet& ps, int nTube, MdTable* md,
                           const std::string& label, const TubeFn& tube)
{
    double exact[kNTest], tube1[kNTest], tube2[kNTest];
    const uint64_t b1[3] = {7, 11, 13}, b2[3] = {17, 19, 23};
    tube(b1, tube1);
    tube(b2, tube2);
    double maxAll = 0, maxKept = 0, maxTubeUnc = 0;
    out("  quadrature of %s: test f      exact over T    all points (rel err)     T minus tube      kept points (rel err)\n",
        label.c_str());
    for (int k = 0; k < kNTest; k++) {
        exact[k] = exactTetIntegral(T, k);
        double sAll = 0, sKept = 0;
        for (size_t i = 0; i < ps.allXi.size(); i++) sAll += ps.allW[i]*testFn(k, T.toX(ps.allXi[i]));
        for (size_t i = 0; i < ps.size(); i++) sKept += ps.w[i]*testFn(k, ps.X[i]);
        const double tube = 0.5*(tube1[k] + tube2[k]), tubeUnc = 0.5*std::fabs(tube1[k] - tube2[k]);
        const double ref = exact[k] - tube;
        const double eAll = (sAll - exact[k])/std::fabs(exact[k]), eKept = (sKept - ref)/std::fabs(ref);
        maxAll = std::max(maxAll, std::fabs(eAll));
        maxKept = std::max(maxKept, std::fabs(eKept));
        maxTubeUnc = std::max(maxTubeUnc, tubeUnc/std::fabs(ref));
        out("      %-4s  %+.8f   %+.8f (%+.1e)   %+.8f   %+.8f (%+.1e)\n", kTestName[k], exact[k], sAll, eAll, ref, sKept,
            eKept);
    }
    out("    sum of all weights = %.15f vs vol T = %.15f (exact by construction: sum of leaf volumes)\n",
        [&] { double s = 0; for (double v : ps.allW) s += v; return s; }(), std::fabs(T.detJ)/6.0);
    if (hasTube)
        out("    eps-tube: volume %.6f by QMC (two Halton sequences of %d points, half-difference %.1e relative to vol(T\\tube));\n"
            "    %d of %zu points dropped (within eps = %.5f of the crack)\n",
            0.5*(tube1[0] + tube2[0]), nTube, maxTubeUnc, ps.nDropped, ps.allXi.size(), ps.eps);
    out("    max relative quadrature error: all points %.1e, kept points %.1e\n", maxAll, maxKept);
    if (md)
        md->rows.push_back({label, std::to_string(ps.size()), std::to_string(ps.nDropped), fmtd("%.2e", maxAll),
                            !hasTube ? "(no tube)" : fmtd("%.2e", maxKept),
                            !hasTube ? "-" : fmtd("%.6f", 0.5*(tube1[0] + tube2[0]))});
    return maxKept;
}
double quadratureCheck(Report& out, const TetFrame& T, const XiCrack& xc, const PointSet& ps, int nTube, MdTable* md,
                       const std::string& label)
{
    return quadratureCheckWith(out, T, !xc.cut.empty(), ps, nTube, md, label, [&](const uint64_t b[3], double o[kNTest]) {
        tubeIntegrals(T, xc.cut, ps.eps, nTube, b, o);
    });
}

// =====================================================================================
// Fields and solves
// =====================================================================================
// One boundary-data set on a tet: +-g (or 0) on the crack, a physical function on the Dirichlet faces.
struct FieldSpec {
    std::string name;
    bool crack = false;   // H = +-g on the crack (else 0)
    GType g = GType::One;
    std::function<double(const Vec3& X)> face;   // Dirichlet faces, physical X; empty: 0
    // face, as a function on all of space, is harmonic (every face function of checks 0 and 1 and of the
    // sanity tests is linear), so it is its own harmonic extension and the gradient's first ball may cross
    // a Dirichlet face (WalkParams::oddMirrorDirichlet). false for data without a known harmonic extension.
    bool faceHarmonic = true;
};

struct OwnRun {   // one multi-data solve of the own loop (double precision)
    std::vector<FieldSpec> fields;
    MultiSums sums;
    SolveStats stats;
};

// The own loop's walk parameters in the gate: the first gate's, except
//  - the termination shell is 1e-6 L instead of 1e-4 L. A walk that ends within the shell of a Dirichlet
//    face takes the data at its projection, so H is biased by up to shell |grad H|; the quadrature points
//    come arbitrarily close to the Dirichlet faces (the crack's eps = 0.01 L exclusion keeps them away from
//    the crack only), and at 1e-4 L the bias showed in the per-point z tests within ~1e-3 of a face;
//  - the gradient's first ball may cross the nearest Neumann face (mirrorGradBall), otherwise the
//    gradient variance ~ 1/(distance to a Neumann face)^2 makes S, b and the energies unusable;
//  - with every field's face data harmonic (FieldSpec::faceHarmonic), it may also cross the nearest
//    Dirichlet face (oddMirrorDirichlet): next to a Dirichlet face the plain ball's gradient is heavy-tailed
//    (variance ~ 1/distance, from rare walks that leave the face's neighbourhood), so the noise of S and b
//    sits in the few points closest to a Dirichlet face and their jackknife s.e. varies several-fold
//    between seeds (stage 1, check 0: s.e. of S_00 from 0.0023 to 0.019 over 8 seeds, true sd 0.0065).
constexpr double kGateShellRel = 1e-6;
WalkParams gateWalkParams(const std::array<FaceBC, 4>& bc, bool mirror, bool oddMirror)
{
    WalkParams wp;
    wp.faceBC = bc;
    wp.absorbShellRel = kGateShellRel;
    wp.mirrorGradBall = mirror;
    wp.oddMirrorDirichlet = mirror && oddMirror;
    return wp;
}

// mirror: the mirror balls (both kinds, the odd one only if every field's face data is harmonic); false: the
// plain first ball everywhere. gFront: the xi cut whose front line the crack's g is measured from (check 4, where the
// physical front may lie outside the tet); null: the strips' own front (checks 0-3: the front crosses the tet or
// there is none).
OwnRun runOwn(const TetFrame& T, const XiCrack& xc, const std::array<FaceBC, 4>& bc, const std::vector<FieldSpec>& fields,
              const std::vector<Vec3>& ptsXi, int nWalks, int nBatches, uint64_t seed, bool mirror = true,
              const Cut* gFront = nullptr)
{
    OwnRun r;
    r.fields = fields;
    bool allHarmonic = true;
    for (const FieldSpec& f : fields) allHarmonic = allHarmonic && (!f.face || f.faceHarmonic);
    const WalkParams wp = gateWalkParams(bc, mirror, allHarmonic);
    QueryPoints q;
    const std::vector<char> grad(ptsXi.size(), 1);   // the gradient's full first ball at every point
    q.add(ptsXi, &grad, seed);
    const Cut front = gFront ? *gFront : xc.strips.frontCut(xc.strips.nStrips());   // g of the crack, in xi (L = sqrt(2))
    const int K = int(fields.size());
    const TerminalValues tv = [&](const WalkTerminal& t, double* u) {
        if (t.face >= 0) {
            const Vec3 X = T.toX(t.x);
            for (int k = 0; k < K; k++) u[k] = fields[k].face ? fields[k].face(X) : 0.0;
        } else {
            for (int k = 0; k < K; k++) u[k] = fields[k].crack ? double(t.side)*front.g(t.x, fields[k].g) : 0.0;
        }
    };
    // the harmonic extensions of the face data (odd mirror): the face functions themselves (0 without one)
    const HarmonicExtension ext = [&](const Vec3& xi, double* phi) {
        const Vec3 X = T.toX(xi);
        for (int k = 0; k < K; k++) phi[k] = fields[k].face ? fields[k].face(X) : 0.0;
    };
    r.sums = solveMulti(referenceTet(), xc.strips, q, nWalks, nBatches, wp, K, tv, &r.stats, &ext);
    return r;
}

// One Zombie solve (float) of one field, H (and grad H if wantGrad) at every point, once per seed. Zombie
// gets the gate's termination shell (1e-6 L, minimum star radius 1e-7 L: the first gate's ratio of 10), so
// its walks stop within stopDistance() = shell + the float face-offset margin = 6.2e-5 of the absorbing
// boundary; with the first gate's shell (1e-4 L, firstGateShell) that is 2.0e-4. The walk ends there and
// takes the data at the projected point, a bias of order stop |dH/dn| that the harmonic extension carries
// into the whole tet: at 2.0e-4 the own loop and Zombie differ measurably in check 1 (see there).
struct ZombieRun {
    std::vector<std::vector<Estimate>> est;   // [solve][point]
    SolveStats stats;
};
ZombieRun runZombie(const TetFrame& T, const XiCrack& xc, const std::array<FaceBC, 4>& bc, const FieldSpec& f,
                    const std::vector<Vec3>& ptsXi, int nWalks, bool wantGrad, const std::vector<uint64_t>& seeds,
                    bool firstGateShell = false)
{
    ZombieRun r;
    WalkParams wp;
    if (!firstGateShell) {
        wp.absorbShellRel = kGateShellRel;
        wp.reflectShellRel = 0.1*kGateShellRel;
    }
    wp.faceBC = bc;
    wp.cutData = f.crack;
    if (f.face) wp.faceData = [&T, face = f.face](int, const Vec3& xi) { return face(T.toX(xi)); };
    CutSolver solver(referenceTet(), xc.cut, f.g, wp);
    const std::vector<char> grad(ptsXi.size(), wantGrad ? 1 : 0);
    for (uint64_t s : seeds) {
        r.est.emplace_back();
        r.stats.add(solver.solve(ptsXi, &grad, nWalks, s, r.est.back()));
    }
    return r;
}

// =====================================================================================
// Estimators
// =====================================================================================
// H of the combination sum_k coef_k field_k at point i, and the variance of that mean: the per-walk
// variance of the combined value (from the K x K second moments of the same walks) over n.
void comboH(const MultiSums& s, size_t i, const std::vector<double>& coef, double& H, double& var, long long& nOut)
{
    const int K = s.K;
    double n = 0;
    std::vector<double> m(size_t(K), 0.0);
    for (int b = 0; b < s.nBatches; b++) {
        const double* d = s.batchData(i, b);
        n += d[0];
        for (int k = 0; k < K; k++) m[k] += d[1 + k];
    }
    nOut = (long long)n;
    H = 0;
    var = 0;
    if (n < 1) return;
    for (int k = 0; k < K; k++) H += coef[k]*m[k]/n;
    if (n < 2) return;
    const double* uu = s.uuData(i);
    double second = 0;
    for (int k = 0; k < K; k++)
        for (int l = 0; l < K; l++) second += coef[k]*coef[l]*uu[k*K + l];
    var = std::max(0.0, (second - n*H*H)/(n - 1))/n;
}

// grad (xi) of the combination from batch b of point i: the first-step estimator w (u - c) with the
// leave-one-out control variate c_j = mean of u over the batch's other walks, i.e.
// [sum w u - (sum u sum w - sum w u)/(n - 1)]/n. E[w] = 0 exactly (uniform directions), so this is
// unbiased for any batch size (the control variate of estimateFromSums, H of the same walks, is biased by
// a factor 1 - 1/n, which matters for 40-walk batches).
Vec3 comboBatchGrad(const MultiSums& s, size_t i, int b, const std::vector<double>& coef)
{
    const int K = s.K;
    const double* d = s.batchData(i, b);
    const double n = d[0];
    if (n < 2) return Vec3::Zero();
    const Vec3 sw(d[1 + K], d[2 + K], d[3 + K]);
    double su = 0;
    Vec3 swu = Vec3::Zero();
    for (int k = 0; k < K; k++) {
        if (coef[k] == 0) continue;
        su += coef[k]*d[1 + k];
        swu += coef[k]*Vec3(d[4 + K + 3*k], d[5 + K + 3*k], d[6 + K + 3*k]);
    }
    return (swu - (su*sw - swu)/(n - 1.0))/n;
}

std::vector<double> unitCoef(int K, int k)
{
    std::vector<double> c(size_t(K), 0.0);
    c[size_t(k)] = 1.0;
    return c;
}

// Gram matrix G_ql = integral over the kept points of grad f_q . grad f_l (physical), from B independent
// batches of per-point gradient estimates g_b. Products of two estimates from the same walks are biased by
// their covariance, so G uses only products of different batches: per point
// sum_{a != b} g_a . g_b / (B (B-1)) = (T1 . T1 - T2)/(B (B-1)), T1 = sum_b g_b, T2 = sum_b g_b . g_b.
// Standard errors: delete-one-batch jackknife over the B batches. Also the naive product of the full means
// and the two-halves split (batches 0..B/2-1 against B/2..B-1, symmetrized) as diagnostics.
struct Gram {
    int Q = 0, B = 0;
    Eigen::MatrixXd G, se, naive, halves;
    std::vector<Eigen::MatrixXd> jack;   // delete-one-batch replicates
};

Gram gramFromBatches(const std::vector<double>& w, int B, int Q, const std::function<void(size_t, int, Vec3*)>& grads)
{
    Gram R;
    R.Q = Q;
    R.B = B;
    const size_t P = w.size();
    const int nChunks = int(std::min<size_t>(64, std::max<size_t>(1, P)));
    struct Acc {
        Eigen::MatrixXd G, naive, halves;
        std::vector<Eigen::MatrixXd> jack;
    };
    std::vector<Acc> parts(static_cast<size_t>(nChunks));
    tbb::parallel_for(0, nChunks, [&](int c) {
        Acc& A = parts[size_t(c)];
        A.G = A.naive = A.halves = Eigen::MatrixXd::Zero(Q, Q);
        A.jack.assign(size_t(B), Eigen::MatrixXd::Zero(Q, Q));
        std::vector<Vec3> g(size_t(B)*Q), T1(static_cast<size_t>(Q)), h1(static_cast<size_t>(Q)), h2(static_cast<size_t>(Q));
        Eigen::MatrixXd T2(Q, Q);
        for (size_t p = P*size_t(c)/size_t(nChunks); p < P*size_t(c + 1)/size_t(nChunks); p++) {
            for (int b = 0; b < B; b++) grads(p, b, &g[size_t(b)*Q]);
            for (int q = 0; q < Q; q++) {
                T1[q] = h1[q] = h2[q] = Vec3::Zero();
                for (int b = 0; b < B; b++) {
                    T1[q] += g[size_t(b)*Q + q];
                    (2*b < B ? h1[q] : h2[q]) += g[size_t(b)*Q + q];
                }
                h1[q] /= double(B/2);
                h2[q] /= double(B - B/2);
            }
            for (int q = 0; q < Q; q++)
                for (int l = 0; l < Q; l++) {
                    double t2 = 0;
                    for (int b = 0; b < B; b++) t2 += g[size_t(b)*Q + q].dot(g[size_t(b)*Q + l]);
                    T2(q, l) = t2;
                    A.G(q, l) += w[p]*(T1[q].dot(T1[l]) - t2)/double(B*(B - 1));
                    A.naive(q, l) += w[p]*T1[q].dot(T1[l])/double(B*B);
                    A.halves(q, l) += w[p]*0.5*(h1[q].dot(h2[l]) + h2[q].dot(h1[l]));
                }
            if (B < 3) continue;
            for (int a = 0; a < B; a++)
                for (int q = 0; q < Q; q++)
                    for (int l = 0; l < Q; l++) {
                        const Vec3 &gq = g[size_t(a)*Q + q], &gl = g[size_t(a)*Q + l];
                        A.jack[a](q, l) += w[p]*((T1[q] - gq).dot(T1[l] - gl) - (T2(q, l) - gq.dot(gl)))/double((B - 1)*(B - 2));
                    }
        }
    });
    R.G = R.naive = R.halves = Eigen::MatrixXd::Zero(Q, Q);
    R.jack.assign(size_t(B), Eigen::MatrixXd::Zero(Q, Q));
    for (const Acc& A : parts) {   // in chunk order: deterministic
        R.G += A.G;
        R.naive += A.naive;
        R.halves += A.halves;
        for (int a = 0; a < B; a++) R.jack[a] += A.jack[a];
    }
    R.se = Eigen::MatrixXd::Zero(Q, Q);
    if (B >= 3) {
        Eigen::MatrixXd mean = Eigen::MatrixXd::Zero(Q, Q);
        for (const auto& J : R.jack) mean += J/double(B);
        for (const auto& J : R.jack) R.se += (J - mean).cwiseAbs2();
        R.se = (R.se*double(B - 1)/double(B)).cwiseSqrt();
    }
    return R;
}

// the Gram matrix of the combinations coefs of an own-loop run (physical gradients: J^-T grad_xi)
Gram gramOwn(const OwnRun& r, const PointSet& ps, const TetFrame& T, const std::vector<std::vector<double>>& coefs)
{
    const int Q = int(coefs.size());
    return gramFromBatches(ps.w, r.sums.nBatches, Q, [&](size_t p, int b, Vec3* g) {
        for (int q = 0; q < Q; q++) g[q] = T.JinvT*comboBatchGrad(r.sums, p, b, coefs[q]);
    });
}

// eigenvalues (ascending) of the symmetric part of a matrix, and their jackknife standard errors
void eigenWithSe(const Gram& gm, const Eigen::MatrixXd& M, const std::vector<Eigen::MatrixXd>& jack,
                 Eigen::VectorXd& ev, Eigen::VectorXd& se)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(0.5*(M + M.transpose()));
    ev = es.eigenvalues();
    se = Eigen::VectorXd::Zero(ev.size());
    const int B = gm.B;
    if (B < 3 || jack.empty()) return;
    std::vector<Eigen::VectorXd> e;
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(ev.size());
    for (const auto& J : jack) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ej(0.5*(J + J.transpose()));
        e.push_back(ej.eigenvalues());
        mean += e.back()/double(B);
    }
    for (const auto& v : e) se += (v - mean).cwiseAbs2();
    se = (se*double(B - 1)/double(B)).cwiseSqrt();
}

// physical gradients of the linear shape functions N_k (= lambda_k) of the tet's local vertices
std::array<Vec3, 4> gradN(const TetFrame& T)
{
    std::array<Vec3, 4> gN;
    gN[1] = T.JinvT.col(0);
    gN[2] = T.JinvT.col(1);
    gN[3] = T.JinvT.col(2);
    gN[0] = -(gN[1] + gN[2] + gN[3]);
    return gN;
}

// linear-tet Laplacian stiffness of the physical tet: K_ij = V grad N_i . grad N_j
Eigen::Matrix4d stiffness(const TetFrame& T)
{
    const double V = std::fabs(T.detJ)/6.0;
    const std::array<Vec3, 4> gN = gradN(T);
    Eigen::Matrix4d K;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) K(i, j) = V*gN[i].dot(gN[j]);
    return K;
}

std::string matRow(const Eigen::MatrixXd& M, int i, const char* f = "%+.5f")
{
    std::string s;
    for (int j = 0; j < M.cols(); j++) s += (j ? "  " : "") + fmtd(f, M(i, j));
    return s;
}

// z statistics of a per-point difference.
// mean z is the spec's criterion for check 1 (one estimator against itself: the direct run and the combination
// have the same per-walk distribution and the same N, so E[z] = 0 by symmetry). Between two DIFFERENT estimators
// it is biased whenever the walks' values are skewed and the two N differ: each z is normalized by its own
// sample variance, which is correlated with the sample mean (Cov(mean, s^2) = mu_3/N), so to first order
// E[z] = -(gamma/2) (1/N_a^2 - 1/N_b^2) / (1/N_a + 1/N_b)^(3/2) (gamma = the skewness): +0.008 gamma for the
// own loop (N = 2560) against Zombie (N = 1024). Measured on check 1's P_k, the same offset shows between two
// own-loop runs with N = 2560 and 1024, and none at equal N. Comparisons of different estimators therefore use
// the pooled difference instead: sum_p diff_p / sqrt(sum_p var_p) over the points, unbiased with any skewness
// (passPooled: |pooled z| < 3, and |z| < 3 at >= 99% of the points).
// Dispersion (the reference cross-check of checks 2 and 3), which sees a zero-mean discrepancy that the pooled z
// cannot. mean(diff^2 - var) is an unbiased estimate of the mean squared discrepancy between the two estimators' means
// (E diff^2 = the two variances + discrepancy^2, and var estimates the variances without bias); in the units of the
// error tables, ncRmsDiff = sqrt(max(0, mean(diff^2 - var)))/2, and ncRmsDiffBound = sqrt(max(0, mean(diff^2 - var))
// + 2 s.e.)/2 (s.e. from the spread over the points) is the largest RMS discrepancy the comparison does not exclude
// (about 2 sigma; with the estimate clamped at 0 it never falls below the resolution). mean z^2 (1 if they agree,
// s.e. ~ sqrt(2/n)) is shown too, but it is biased upward by the noise of the per-point variances (1/var is
// convex), most where nearly all of a point's walks give the same value.
struct ZStats {
    long long n = 0, below3 = 0, zeroVar = 0;
    double sumZ = 0, sumAbsZ = 0, maxAbsZ = 0, maxAbsDiffZeroVar = 0, sumDiff = 0, sumVar = 0;
    double sumZ2 = 0, sumZ4 = 0, sumEx = 0, sumEx2 = 0;   // z^2, z^4, diff^2 - var, (diff^2 - var)^2
    size_t worst = 0;
    void add(size_t i, double diff, double var)
    {
        if (!(var > 1e-24)) {   // zero sample variance, up to round-off (all walks gave the same value)
            zeroVar++;
            maxAbsDiffZeroVar = std::max(maxAbsDiffZeroVar, std::fabs(diff));
            return;
        }
        const double z = diff/std::sqrt(var), ex = diff*diff - var;
        n++;
        below3 += std::fabs(z) < 3;
        sumZ += z;
        sumAbsZ += std::fabs(z);
        sumDiff += diff;
        sumVar += var;
        sumZ2 += z*z;
        sumZ4 += z*z*z*z;
        sumEx += ex;
        sumEx2 += ex*ex;
        if (std::fabs(z) > maxAbsZ) {
            maxAbsZ = std::fabs(z);
            worst = i;
        }
    }
    double frac() const { return n ? double(below3)/double(n) : 0.0; }
    double mean() const { return n ? sumZ/double(n) : 0.0; }
    double meanAbs() const { return n ? sumAbsZ/double(n) : 0.0; }   // sqrt(2/pi) = 0.80 for independent estimates
    double meanBound() const { return n ? 3.0/std::sqrt(double(n)) : 0.0; }
    double pooled() const { return sumVar > 0 ? sumDiff/std::sqrt(sumVar) : 0.0; }
    double meanDiff() const { return n ? sumDiff/double(n) : 0.0; }        // mean difference over the points
    double meanDiffSe() const { return n ? std::sqrt(sumVar)/double(n) : 0.0; }
    bool pass() const { return n > 0 && frac() >= 0.99 && std::fabs(mean()) < meanBound(); }
    bool passPooled() const { return n > 0 && frac() >= 0.99 && std::fabs(pooled()) < 3.0; }
    double meanZ2() const { return n ? sumZ2/double(n) : 0.0; }
    double meanZ2Se() const   // from the spread of z^2 over the points
    {
        return n > 1 ? std::sqrt(std::max(0.0, sumZ4/double(n) - meanZ2()*meanZ2())/double(n)) : 0.0;
    }
    double meanEx() const { return n ? sumEx/double(n) : 0.0; }
    double meanExSe() const { return n > 1 ? std::sqrt(std::max(0.0, sumEx2/double(n) - meanEx()*meanEx())/double(n)) : 0.0; }
    double ncRmsDiff(double scale) const { return std::sqrt(std::max(0.0, meanEx()))/scale; }
    double ncRmsDiffBound(double scale) const { return std::sqrt(std::max(0.0, meanEx()) + 2.0*meanExSe())/scale; }
};

// Per-point z comparisons of two different estimators (against an exact solution, or the own loop against
// Zombie) leave out the points within kCompareBand of a Dirichlet face, for two reasons:
//  - termination bias: a walk that ends within its stop distance of a Dirichlet face takes the data at the
//    projected point, a bias of up to stop |grad H| (1.4e-6 for the own loop here, 2.0e-4 for Zombie);
//  - rare excursions: from a point at distance delta, a walk leaves the face's neighbourhood with
//    probability ~ delta/l (l ~ 0.3 here). With N delta/l of order 10 or less the per-walk values are a
//    near-constant plus a few large excursions (or none: all N walks give the same value and the sample
//    variance is 0), so the sample variance, and z, are unreliable. At delta = 0.01, N delta/l ~ 85 for the
//    own loop (N = 2560) and ~ 34 for Zombie (N = 1024).
// The band points are only counted and their largest |difference| reported. Comparisons of one estimator
// with itself (check 1's direct run against its combination of the same estimator) and every integral
// (S, b, energies) use all points.
constexpr double kCompareBand = 0.01;
double ownStopDistance() { return kGateShellRel*kSqrt2; }
double zombieStopDistance(bool firstGateShell = false)
{
    return (firstGateShell ? WalkParams{}.absorbShellRel : kGateShellRel)*kSqrt2 + faceOffsetMargin(referenceTet());
}
// distance from xi to the nearest Dirichlet face plane of the reference tet (+inf without one); J is
// orthogonal for A and B, so this is the physical distance too
double dirichletDist(const std::array<FaceBC, 4>& bc, const Vec3& xi)
{
    const Tet ref = referenceTet();
    double d = kInf;
    for (int f = 0; f < 4; f++) {
        if (bc[f] != FaceBC::Dirichlet) continue;
        const auto& F = ref.faces[f];
        const Vec3 n = (ref.v[F[1]] - ref.v[F[0]]).cross(ref.v[F[2]] - ref.v[F[0]]).normalized();
        d = std::min(d, n.dot(ref.v[F[0]] - xi));
    }
    return d;
}
// band points: their count and the largest |diff| - 3 s.e. (to compare with a bias bound)
struct BandStats {
    long long n = 0;
    double maxExcess = -kInf, maxAbsDiff = 0;
    void add(double diff, double var)
    {
        n++;
        maxAbsDiff = std::max(maxAbsDiff, std::fabs(diff));
        maxExcess = std::max(maxExcess, std::fabs(diff) - 3.0*std::sqrt(std::max(0.0, var)));
    }
};

// =====================================================================================
// Settings of the gate
// =====================================================================================
struct GateParams {
    bool quick = false;
    int nPoints = 4096;     // volume quadrature points per tet (about)
    int nFacePerTri = 16;   // F points: 64 sub-triangles x 16
    int nWalks = 2560;      // per-tet recipe runs: 10x the default N
    int nBatches = 64;      // independent batches per point (split-sample products, jackknife)
    int nZombie = 1024;     // Zombie cross-check: 4x the default N
    int nZombieBatches = 16;
    int nRef = 25600;       // the reference on U (checks 2 and 3) and on O (check 4): 100x the default N
    int nRefBatches = 64;
    int nTube = 1 << 21;    // QMC points per sequence for the eps-tube integrals
    int nRingPoints = 4096; // check 4: volume quadrature points per ring tet (about)
    double eps = 0.01*kSqrt2;
    uint64_t seed = Settings{}.seed;
    bool tipTrace = false;  // check 3 (--face-trace tip): Level B's J is the crack's own 2D front solution (levelBTrace)
};

GateParams gateParams(bool quick, bool tipTrace = false)
{
    GateParams p;
    p.quick = quick;
    p.tipTrace = tipTrace;
    if (quick) {
        p.nPoints = 512;
        p.nWalks = 256;
        p.nBatches = 16;
        p.nZombie = 128;
        p.nZombieBatches = 8;
        p.nRef = 1024;
        p.nRefBatches = 16;
        p.nTube = 1 << 17;
        p.nRingPoints = 512;
    }
    return p;
}

std::array<FaceBC, 4> onlyF()
{
    std::array<FaceBC, 4> bc{FaceBC::Neumann, FaceBC::Neumann, FaceBC::Neumann, FaceBC::Neumann};
    bc[kFaceF] = FaceBC::Dirichlet;
    return bc;
}
std::array<FaceBC, 4> allDirichlet() { return {FaceBC::Dirichlet, FaceBC::Dirichlet, FaceBC::Dirichlet, FaceBC::Dirichlet}; }
std::array<FaceBC, 4> allNeumann() { return {FaceBC::Neumann, FaceBC::Neumann, FaceBC::Neumann, FaceBC::Neumann}; }

// =====================================================================================
// Geometry setup printout
// =====================================================================================
void printCutPolygons(Report& out, const TetFrame& T, const XiCrack& xc)
{
    if (xc.cut.empty()) {
        out("    %s: no crack in this tet (the plane misses it or there is none)\n", T.name);
        return;
    }
    out("    %s: t = %.6f; xi plane n = %s, d = %.6f, front direction u = %s, front s = %.6f, L = %.6f\n", T.name, xc.t,
        fmt3(xc.cut.n).c_str(), xc.cut.d, fmt3(xc.cut.u).c_str(), xc.cut.s, xc.cut.L);
    std::string sxi, sx;
    for (const Vec3& p : xc.cut.polygon) {
        sxi += " " + fmt3(p, 3);
        sx += " " + fmt3(T.toX(p), 3);
    }
    out("       polygon (xi):%s\n       polygon (X): %s\n", sxi.c_str(), sx.c_str());
    out("       front vertices in X at y = y_f: %d vertices, max |y - y_f| = %.1e -> %s\n", xc.nFront, xc.frontErr,
        xc.nFront >= 2 && xc.frontErr < 1e-12 ? "OK" : "FAILED");
}

void printUnionPolygon(Report& out, const CrackSpec& cs)
{
    double t = 0;
    const Cut c = unionCut(cs, &t);
    std::string s;
    for (const Vec3& p : c.polygon) s += " " + fmt3(p, 3);
    double err = 0;
    int nf = 0;
    for (const Vec3& p : frontVertices(c)) {
        err = std::max(err, std::fabs(p.y() - cs.yFront));
        nf++;
    }
    out("    U (physical, the reference): t = %.6f, L overridden %.4f -> %.6f; polygon (X):%s\n", t, unionTet().L, c.L,
        s.c_str());
    out("       front vertices at y = y_f: %d, max |y - y_f| = %.1e -> %s\n", nf, err,
        nf >= 2 && err < 1e-12 ? "OK" : "FAILED");
}

// returns false if a geometric assertion failed
bool printGeometry(Report& out, const GateParams& gp)
{
    bool ok = true;
    const TetFrame A = frameA(), B = frameB();
    const Tet U = unionTet();
    out("=== Coupling gate: geometry setup (physical coordinates X) ===\n");
    out("global vertices: v0 = %s, v1 = %s, v2 = %s, v3 = %s (A's apex), v4 = %s (B's apex)\n", fmt3(kGlobal[0], 0).c_str(),
        fmt3(kGlobal[1], 0).c_str(), fmt3(kGlobal[2], 0).c_str(), fmt3(kGlobal[3], 0).c_str(), fmt3(kGlobal[4], 0).c_str());
    for (const TetFrame* T : {&A, &B}) {
        out("tet %s: local vertices 0..3 = global v%d v%d v%d v%d = %s %s %s %s; X = X0 + J xi,\n", T->name, T->global[0],
            T->global[1], T->global[2], T->global[3], fmt3(T->X[0], 0).c_str(), fmt3(T->X[1], 0).c_str(),
            fmt3(T->X[2], 0).c_str(), fmt3(T->X[3], 0).c_str());
        out("   J = [[%g,%g,%g],[%g,%g,%g],[%g,%g,%g]], det J = %+g, J^T J = I: %s; solved on the reference tet in xi\n",
            T->J(0, 0), T->J(0, 1), T->J(0, 2), T->J(1, 0), T->J(1, 1), T->J(1, 2), T->J(2, 0), T->J(2, 1), T->J(2, 2), T->detJ,
            (T->J.transpose()*T->J - Eigen::Matrix3d::Identity()).norm() == 0 ? "yes (distances, eps and g are the same in xi)"
                                                                              : "NO");
        for (int o = 0; o < 4; o++)
            out("   face f%d (opposite local vertex %d): %s, corners %s%s\n", o, o, facePlane(T->X, o).c_str(),
                faceCorners(*T, o).c_str(), o == kFaceF ? "  <- F (xi_z = 0 in this frame)" : "");
    }
    out("union U = A u B = conv{v1, v2, v3, v4} (a tetrahedron: x >= 0, y >= 0, x + y + |z| <= 1): faces");
    for (int o = 0; o < 4; o++) out(" %s", facePlane(U.v, o).c_str());
    out(";\n   longest edge L_U = %.4f (v3-v4); the crack's g uses L = sqrt(2) = %.6f in A, B and U\n", U.L, kSqrt2);
    out("shared face F = triangle z = 0 with corners v0 v1 v2; DOFs on F are indexed by global vertex (v1 and v2 swap\n"
        "   between A's local frame (local 1 = v1, 2 = v2) and B's (local 1 = v2, 2 = v1))\n");

    // hats: 1 at their corner, 0 at the others, on A's and on B's F, in X and in xi
    {
        int bad = 0, n = 0;
        for (const TetFrame* T : {&A, &B})
            for (int g = 0; g < 3; g++)
                for (int c = 0; c < 3; c++) {
                    const double want = g == c ? 1.0 : 0.0;
                    const Vec3 X = kGlobal[c], xi = T->toXi(X);
                    const double hX = hatF(g, X);                        // physical hat at the physical corner
                    const double hXi = hatF(g, T->toX(xi));              // physical hat evaluated at X = J xi
                    const double hLoc = lambdaOf(*T, T->localOf(g), X);  // the tet's own barycentric of that vertex
                    const bool onF = std::fabs(xi.z()) == 0;             // the corner lies on xi_z = 0 in this frame
                    n++;
                    if (hX != want || hXi != want || hLoc != want || !onF) bad++;
                }
        out("hat check: hat_g (g = 0, 1, 2) at the corners v0, v1, v2 of F, in A's and B's frame, evaluated in X, at\n"
            "   X = X0 + J xi, and as the tet's barycentric lambda_local(g): %d of %d (tet, hat, corner) cases wrong; every\n"
            "   corner lies on xi_z = 0 in both frames -> %s\n", bad, n, bad ? "FAILED" : "OK");
        ok = ok && bad == 0;
        out("   B's local hats in xi: lambda_1(xi) = xi_x = y = hat_2, lambda_2(xi) = xi_y = x = hat_1 (the v1/v2 swap)\n");
    }

    // crack configurations
    out("\ncrack configurations (front direction +y, crack region y < y_f; t per frame from the physical extent):\n");
    std::vector<CrackSpec> cracks = {crackCheck3(0.30), crackCheck2(0.60), crackCheck2(0.35), crackCheck3(0.15),
                                     crackCheck3(0.45)};
    const char* use[5] = {"check 1 (A alone) and check 3", "check 2", "check 2", "check 3", "check 3"};
    for (size_t ci = 0; ci < cracks.size(); ci++) {
        const CrackSpec& cs = cracks[ci];
        out("  crack %s  [%s]\n", cs.name.c_str(), use[ci]);
        const XiCrack xa = xiCrack(A, cs), xb = xiCrack(B, cs);
        printCutPolygons(out, A, xa);
        printCutPolygons(out, B, xb);
        printUnionPolygon(out, cs);
        for (const XiCrack* x : {&xa, &xb})
            if (!x->cut.empty()) ok = ok && x->nFront >= 2 && x->frontErr < 1e-12;
        // "+g" side and g at the same physical points in A's, B's and U's frames
        const Cut cu = unionCut(cs);
        std::mt19937_64 rng(mixSeed(gp.seed, 7900 + ci));
        std::uniform_real_distribution<double> Ud(0.0, 1.0);
        int nSide = 0, badSide = 0;
        double gErr = 0;
        for (int k = 0; k < 2000; k++) {
            double l[4], sl = 0;
            for (double& v : l) sl += (v = -std::log(Ud(rng)));
            const Vec3 xi(l[1]/sl, l[2]/sl, l[3]/sl);
            for (const auto& [T, x] : {std::pair<const TetFrame*, const XiCrack*>{&A, &xa}, {&B, &xb}}) {
                if (x->cut.empty()) continue;
                const Vec3 X = T->toX(xi);
                const double hX = cs.n.dot(X) - cs.d;
                if (std::fabs(hX) < 1e-9) continue;
                nSide++;
                if ((x->cut.side(xi) > 0) != (hX > 0) || (cu.side(X) > 0) != (hX > 0)) badSide++;
                for (GType g : {GType::One, GType::SqrtFront})
                    gErr = std::max(gErr, std::fabs(x->cut.g(xi, g) - cu.g(X, g)));
            }
        }
        out("    +g side at %d random points (the frame's xi cut and U's cut, each against n.X - d in X): %d disagree; |g_frame(xi) - g_U(X)| <= %.1e\n"
            "    over both g types (L = sqrt(2) everywhere) -> %s\n", nSide, badSide, gErr,
            badSide == 0 && gErr < 1e-12 ? "OK" : "FAILED");
        ok = ok && badSide == 0 && gErr < 1e-12;
    }

    // every run and its face flags
    out("\nruns and per-face flags (D = Dirichlet with the data named, N = Neumann):\n");
    out("  sanity (a) A, no crack: %s\n", faceFlags(A, onlyF(), "H = 1").c_str());
    out("  sanity (b) A, no crack: %s\n", faceFlags(A, onlyF(), "H = x + y").c_str());
    {
        std::array<FaceBC, 4> bc = allNeumann();
        bc[0] = bc[1] = FaceBC::Dirichlet;
        out("  sanity (c) A, no crack, exact solution H = x: %s\n", faceFlags(A, bc, "H = x").c_str());
    }
    out("  check 0 A, no crack (t = 0): %s\n     H0: faces 0; P_k (k = v0..v3): faces lambda_k (the hat of vertex k on every face)\n",
        faceFlags(A, allDirichlet(), "(data per field below)").c_str());
    out("  check 1 A, crack x = 0.30, y_f = 0.30: %s\n     H0: crack +-g, F 0; P_k (k = v0, v1, v2): crack 0, F hat_k; "
        "direct: crack +-g, F sum_k c_k hat_k\n", faceFlags(A, onlyF(), "(data per field)").c_str());
    out("  check 2 (z_c = 0.60 and 0.35) A: %s (Level A: H0 crack +-g, F 0; P_k crack 0, F hat_k); B: uncut, no solve (linear,\n"
        "     K only); reference U: %s\n", faceFlags(A, onlyF(), "(data per field)").c_str(),
        faceFlags(makeFrame("U", {1, 2, 3, 4}), allNeumann(), "", false).c_str());
    out("  check 3 (y_f = 0.15, 0.30 and 0.45) A: %s;\n     B (its own xi frame, the crack pulled back through J_B^-1): %s;\n"
        "     Level A, both tets: H0 crack +-g, F 0; P_k (k = v0, v1, v2) crack 0, F hat_k. Level B, both tets: H0 crack +-g, F J (the\n"
        "     same physical J in A and B); the P_k of Level A reused. Reference U: %s\n",
        faceFlags(A, onlyF(), "(data per field)").c_str(), faceFlags(B, onlyF(), "(data per field)").c_str(),
        faceFlags(makeFrame("U", {1, 2, 3, 4}), allNeumann(), "", false).c_str());
    out("geometry assertions: %s\n", ok ? "all OK" : "SOME FAILED");
    return ok;
}

// =====================================================================================
// gate_report.md
// =====================================================================================
// Every section states its findings between these markers: plain text with "- " bullets, and fixed-width tables
// in fenced blocks (fixedTable). The summary at the top of the report and the console summary at the end of a run
// are copied from them, so they repeat the sections' own numbers.
constexpr const char* kFindingsBegin = "<!-- findings -->\n";
constexpr const char* kFindingsEnd = "<!-- /findings -->\n";
std::string findingsBlock(const std::string& text)
{
    return std::string("### Findings\n\n") + kFindingsBegin + text + kFindingsEnd + "\n";
}
std::string findingsOf(const std::string& body)
{
    const size_t a = body.find(kFindingsBegin);
    if (a == std::string::npos) return "";
    const size_t b = body.find(kFindingsEnd, a);
    return b == std::string::npos ? "" : body.substr(a + std::strlen(kFindingsBegin), b - a - std::strlen(kFindingsBegin));
}
// the findings for the console: without the fences of the tables
std::string findingsPlain(const std::string& f)
{
    std::string o, line;
    std::istringstream in(f);
    while (std::getline(in, line))
        if (line != "```") o += line + "\n";
    return o;
}
// rows of cells as a fixed-width table in a fenced block (the first row is the header)
std::string fixedTable(const std::vector<std::vector<std::string>>& rows)
{
    std::vector<size_t> w;
    for (const auto& r : rows)
        for (size_t c = 0; c < r.size(); c++) {
            if (w.size() <= c) w.push_back(0);
            w[c] = std::max(w[c], r[c].size());
        }
    std::string s = "```\n";
    for (const auto& r : rows) {
        std::string line;
        for (size_t c = 0; c < r.size(); c++) line += r[c] + (c + 1 < r.size() ? std::string(w[c] - r[c].size() + 3, ' ') : "");
        s += line + "\n";
    }
    return s + "```\n";
}
std::string pct(double x, const char* f = "%.2f%%") { return fmtd(f, 100*x); }

// One marked section per part of the gate (<!-- BEGIN id --> ... <!-- END id -->); a run rewrites only the
// sections it ran and keeps the others, in a fixed order: the header and the summary, the sanity tests, the
// checks, the Polyscope views (one section per view, "view_..."), then the appendices. The summary is not
// written by a run: it is rebuilt here, every time the file is written, from the findings of the sections.
bool writeGateReport(const std::string& path, const std::vector<std::pair<std::string, std::string>>& newSecs)
{
    static const char* head[] = {"header", "summary", "sanity", "check0", "check1", "check2", "check3", "check4"};
    static const char* tail[] = {"quadrature", "geometry", "geometry2", "geometry3", "geometry4"};
    std::string old;
    {
        std::ifstream in(path, std::ios::binary);
        if (in) old.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    old.erase(std::remove(old.begin(), old.end(), '\r'), old.end());
    std::vector<std::pair<std::string, std::string>> secs;
    size_t pos = 0;
    while ((pos = old.find("<!-- BEGIN ", pos)) != std::string::npos) {
        const size_t nameEnd = old.find(" -->\n", pos);
        if (nameEnd == std::string::npos) break;
        const std::string name = old.substr(pos + 11, nameEnd - pos - 11);
        const std::string endTag = "<!-- END " + name + " -->\n";
        const size_t end = old.find(endTag, nameEnd);
        if (end == std::string::npos) {
            std::fprintf(stderr, "%s: unrecognised section markers; file left unchanged\n", path.c_str());
            return false;
        }
        secs.emplace_back(name, old.substr(nameEnd + 5, end - nameEnd - 5));
        pos = end + endTag.size();
    }
    for (const auto& ns : newSecs) {
        bool found = false;
        for (auto& s : secs)
            if (s.first == ns.first) {
                s.second = ns.second;
                found = true;
            }
        if (!found) secs.push_back(ns);
    }
    // the summary, from the findings of the sanity tests and the checks present
    std::erase_if(secs, [](const auto& s) { return s.first == "summary"; });
    std::string sum = "## Summary\n\nCopied from the \"Findings\" of the sections below when this file was written; each "
                      "section ends with the date of the run that wrote it.\n\n";
    for (const char* name : {"sanity", "check0", "check1", "check2", "check3", "check4"})
        for (const auto& s : secs)
            if (s.first == name && !findingsOf(s.second).empty()) sum += findingsOf(s.second) + "\n";
    secs.emplace_back("summary", sum);
    std::vector<std::string> order(std::begin(head), std::end(head));
    std::vector<std::string> views;
    for (const auto& s : secs)
        if (s.first.rfind("view_", 0) == 0) views.push_back(s.first);
    std::sort(views.begin(), views.end());
    order.insert(order.end(), views.begin(), views.end());
    order.insert(order.end(), std::begin(tail), std::end(tail));
    std::ofstream f(path, std::ios::binary);
    for (const std::string& name : order)
        for (const auto& s : secs)
            if (s.first == name) {
                std::string body = s.second;
                if (!body.empty() && body.back() != '\n') body += '\n';
                f << "<!-- BEGIN " << name << " -->\n" << body << "<!-- END " << name << " -->\n\n";
            }
    f.close();
    return bool(f);
}

std::string codeBlock(const std::string& s) { return "```\n" + s + (s.empty() || s.back() == '\n' ? "" : "\n") + "```\n"; }

// =====================================================================================
// Sanity tests (no crack)
// =====================================================================================
bool runSanity(Report& out, const GateParams& gp, const PointSet& ps, std::string& md)
{
    const TetFrame A = frameA();
    const XiCrack none{};
    bool ok = true;
    out("\n=== Sanity tests (tet A, no crack) ===\n");
    md += "## Sanity tests (tet A, no crack)\n\n";
    MdTable tb;
    tb.head = {"test", "faces", "solver", "result", "verdict"};
    std::string fa, fc;   // the findings of (a) and (c)

    // (a) F Dirichlet with constant 1, the other faces Neumann: every walk ends on F, so H = 1 exactly
    {
        FieldSpec one{"H", false, GType::One, [](const Vec3&) { return 1.0; }};
        const OwnRun r = runOwn(A, none, onlyF(), {one}, ps.xi, 64, 1, mixSeed(gp.seed, 7000));
        double maxDev = 0, maxVar = 0;
        for (size_t i = 0; i < ps.size(); i++) {
            double H, var;
            long long n;
            comboH(r.sums, i, {1.0}, H, var, n);
            maxDev = std::max(maxDev, std::fabs(H - 1.0));
            maxVar = std::max(maxVar, var);
        }
        const ZombieRun z = runZombie(A, none, onlyF(), one, ps.xi, 64, false, {mixSeed(gp.seed, 7001)});
        double maxDevZ = 0;
        for (const Estimate& e : z.est[0]) maxDevZ = std::max(maxDevZ, double(std::fabs(e.H - 1.0f)));
        const bool pass = maxDev < 1e-12 && maxVar == 0 && r.stats.walksCapped == 0 && maxDevZ == 0;
        ok = ok && pass;
        out("(a) F Dirichlet H = 1, other faces Neumann, %zu points, N = 64: own loop max |H - 1| = %.1e, max variance %.1e\n"
            "    (%s); Zombie (float) max |H - 1| = %.1e (%s) -> %s\n",
            ps.size(), maxDev, maxVar, r.stats.summary().c_str(), maxDevZ, z.stats.summary().c_str(), pass ? "PASS" : "FAIL");
        tb.rows.push_back({"(a) H = 1 on F", "F D, rest N", "own / Zombie (float)",
                           "max abs(H-1) = " + fmtd("%.1e", maxDev) + " / " + fmtd("%.1e", maxDevZ), pass ? "PASS" : "FAIL"});
        fa = std::string("(a) H = 1 on F, the other faces Neumann (every walk ends on F): ") + (pass ? "PASS" : "FAIL") +
             ", max abs(H - 1) = " + fmtd("%.1e", maxDev) + " (own) and " + fmtd("%.1e", maxDevZ) + " (Zombie)";
    }
    // (b) F Dirichlet with x + y, the other faces Neumann: a mixed problem whose solution is not x + y (the
    // faces x = 0 and y = 0 are not zero-flux for x + y); nothing is asserted, a few values are shown
    {
        FieldSpec xy{"H", false, GType::One, [](const Vec3& X) { return X.x() + X.y(); }};
        const OwnRun r = runOwn(A, none, onlyF(), {xy}, ps.xi, 256, 1, mixSeed(gp.seed, 7002));
        out("(b) F Dirichlet H = x + y, other faces Neumann (mixed problem, nothing asserted), N = 256 (%s):\n",
            r.stats.summary().c_str());
        for (size_t i = 0; i < ps.size(); i += std::max<size_t>(1, ps.size()/6)) {
            double H, var;
            long long n;
            comboH(r.sums, i, {1.0}, H, var, n);
            out("    X = %s: H = %.4f +- %.4f (x + y = %.4f)\n", fmt3(ps.X[i]).c_str(), H, std::sqrt(var),
                ps.X[i].x() + ps.X[i].y());
        }
        tb.rows.push_back({"(b) H = x + y on F", "F D, rest N", "own", "values printed (mixed problem)", "not asserted"});
    }
    // (c) an exact mixed problem: faces x = 0 and x + y + z = 1 Dirichlet with H = x, faces y = 0 and z = 0
    // Neumann (x has zero flux through both). H = x, grad H = (1, 0, 0), energy = vol(A) = 1/6. This tests
    // the gradient estimator and the split-sample energy next to Neumann faces.
    {
        std::array<FaceBC, 4> bc = allNeumann();
        bc[0] = bc[1] = FaceBC::Dirichlet;
        FieldSpec fx{"H", false, GType::One, [](const Vec3& X) { return X.x(); }};
        const OwnRun r = runOwn(A, none, bc, {fx}, ps.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7003));
        const double band = kCompareBand;
        ZStats zH, zG;
        BandStats bH;
        for (size_t i = 0; i < ps.size(); i++) {
            double H, var;
            long long n;
            comboH(r.sums, i, {1.0}, H, var, n);
            if (dirichletDist(bc, ps.xi[i]) < band) {
                bH.add(H - ps.X[i].x(), var);
                continue;
            }
            zH.add(i, H - ps.X[i].x(), var);
            // gradient: mean of the batch estimates, s.e. from their spread
            Vec3 m = Vec3::Zero();
            std::vector<Vec3> gb;
            if (dirichletDist(bc, ps.xi[i]) < band) continue;
            for (int b = 0; b < r.sums.nBatches; b++) {
                gb.push_back(A.JinvT*comboBatchGrad(r.sums, i, b, {1.0}));
                m += gb.back()/double(r.sums.nBatches);
            }
            Vec3 v = Vec3::Zero();
            for (const Vec3& g : gb) v += (g - m).cwiseAbs2();
            v /= double(gb.size())*double(gb.size() - 1);
            zG.add(i, m.x() - 1.0, v.x());
            zG.add(i, m.y(), v.y());
            zG.add(i, m.z(), v.z());
        }
        const Gram gm = gramOwn(r, ps, A, {{1.0}});
        const double V = 1.0/6.0;
        const double zE = (gm.G(0, 0) - V)/gm.se(0, 0);
        const bool pass = zH.passPooled() && zG.passPooled() && std::fabs(zE) < 3;
        ok = ok && pass;
        // the same problem without the mirror balls: the plain first ball next to every face
        const OwnRun rp = runOwn(A, none, bc, {fx}, ps.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7003), false);
        const Gram gp0 = gramOwn(rp, ps, A, {{1.0}});
        out("(c) exact solution H = x (faces x = 0, x+y+z = 1 Dirichlet with x; y = 0, z = 0 Neumann), N = %d in %d batches\n"
            "    (%s):\n"
            "    H vs x: |z| < 3 at %.2f%% of %lld points, pooled z = %+.2f (mean H - x = %+.1e +- %.1e), mean z = %+.4f, max |z|\n"
            "    = %.2f (%lld zero-variance points left out); %lld points within %.2f of a Dirichlet face left out of the z tests,\n"
            "    max |H - x| there %.1e\n"
            "    grad H vs (1,0,0), per component: |z| < 3 at %.2f%%, pooled z = %+.2f, mean z = %+.4f, max |z| = %.2f (s.e. from\n"
            "    the batch spread)\n"
            "    energy int |grad H|^2 = %.5f +- %.5f (split-sample, jackknife s.e.) vs 1/6 = %.5f: z = %+.2f; naive %.5f,\n"
            "    halves %.5f -> %s (pass: |z| < 3 at >= 99%% and |pooled z| < 3 for H and grad H, |z| < 3 for the energy)\n"
            "    without the mirror balls (plain first ball, same seed): energy %.5f +- %.5f, naive %.5f\n",
            gp.nWalks, gp.nBatches, r.stats.summary().c_str(), 100*zH.frac(), zH.n, zH.pooled(), zH.meanDiff(), zH.meanDiffSe(),
            zH.mean(), zH.maxAbsZ, zH.zeroVar, bH.n, band, bH.maxAbsDiff, 100*zG.frac(), zG.pooled(), zG.mean(), zG.maxAbsZ,
            gm.G(0, 0), gm.se(0, 0), V, zE, gm.naive(0, 0), gm.halves(0, 0), pass ? "PASS" : "FAIL", gp0.G(0, 0), gp0.se(0, 0),
            gp0.naive(0, 0));
        tb.rows.push_back({"(c) exact H = x", "x=0, x+y+z=1 D; y=0, z=0 N", "own",
                           "H: " + fmtd("%.2f", 100*zH.frac()) + "% abs(z)<3, pooled z " + fmtd("%+.2f", zH.pooled()) +
                               "; grad: " + fmtd("%.2f", 100*zG.frac()) + "%, pooled z " + fmtd("%+.2f", zG.pooled()) +
                               "; energy " + fmtd("%.5f", gm.G(0, 0)) + " +- " + fmtd("%.5f", gm.se(0, 0)) +
                               " vs 1/6 (without the mirror balls: " + fmtd("%.4f", gp0.G(0, 0)) + " +- " +
                               fmtd("%.4f", gp0.se(0, 0)) + ")",
                           pass ? "PASS" : "FAIL"});
        fc = std::string("(c) the exact mixed problem H = x (two faces Dirichlet, two Neumann): ") + (pass ? "PASS" : "FAIL") +
             ", H against x pooled z " + fmtd("%+.2f", zH.pooled()) + ", grad H against (1, 0, 0) pooled z " +
             fmtd("%+.2f", zG.pooled()) + ", energy " + fmtd("%.5f", gm.G(0, 0)) + " +- " + fmtd("%.5f", gm.se(0, 0)) +
             " against 1/6 (without the mirror balls of the gradient: " + fmtd("%.3f", gp0.G(0, 0)) + " +- " +
             fmtd("%.3f", gp0.se(0, 0)) + ")";
    }
    md += findingsBlock("- Sanity tests (tet A, no crack): " + std::string(ok ? "PASS" : "FAIL") + ". " + fa + "; " + fc + ".\n");
    md += tb.str() + "\n";
    return ok;
}

// =====================================================================================
// Check 0: tet A alone, no crack, all four faces Dirichlet: S = K
// =====================================================================================
bool runCheck0(Report& out, const GateParams& gp, const PointSet& ps, std::string& md)
{
    const TetFrame A = frameA();
    const XiCrack none{};
    const Eigen::Matrix4d K = stiffness(A);
    out("\n=== Check 0: machinery. Tet A alone, no crack (t = 0), all four faces Dirichlet ===\n");
    out("H0: faces 0 (so H0 = 0); P_k, k = v0..v3: faces lambda_k, so P_k = lambda_k and S = K. %zu points (uniform),\n"
        "N = %d walks per point in %d batches of %d (own loop, double), both g types (the crack is empty, so g cannot\n"
        "enter: the two runs use the same seed and must agree bit for bit)\n",
        ps.size(), gp.nWalks, gp.nBatches, gp.nWalks/gp.nBatches);
    md += "## Check 0: machinery (tet A, no crack, all four faces Dirichlet)\n\n";
    md += "S = int grad P_i . grad P_j over A, P_k = response to lambda_k on all faces; must equal K = (1/6)[[3,-1,-1,-1],"
          "[-1,1,0,0],[-1,0,1,0],[-1,0,0,1]]. Own loop (double), " + std::to_string(ps.size()) + " uniform points, N = " +
          std::to_string(gp.nWalks) + " in " + std::to_string(gp.nBatches) +
          " batches, mirror balls across the Neumann and (odd, the face data are linear) the Dirichlet faces. S from "
          "products of different batches; s.e. = delete-one-batch jackknife, checked against independent replicate runs. "
          "Pass: max abs(S - K)/s.e. < 3 over A (the spec) and over the points farther than " + fmtd("%.2f", kCompareBand) +
          " from every face (S against K vol/V there, exact since P_k = lambda_k), row sums ~ 0, S PSD, both g modes "
          "bitwise identical, and Zombie agrees. There is no reference solve here: the exact answers are K and P_k = "
          "lambda_k, and g does not enter (no crack).\n\n";
    const size_t mdIntro = md.size();   // the findings go here

    std::vector<FieldSpec> fields(5);
    OwnRun runs[2];
    const GType gts[2] = {GType::One, GType::SqrtFront};
    for (int m = 0; m < 2; m++) {
        fields[0] = FieldSpec{"H0", true, gts[m], {}};
        for (int k = 0; k < 4; k++)
            fields[1 + k] = FieldSpec{"P" + std::to_string(k), false, gts[m], [&A, k](const Vec3& X) { return lambdaOf(A, k, X); }};
        runs[m] = runOwn(A, none, allDirichlet(), fields, ps.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7100));
        out("own loop, %s: %s\n", gTypeName(gts[m]), runs[m].stats.summary().c_str());
    }
    const bool sameBits = runs[0].sums.batch == runs[1].sums.batch && runs[0].sums.uu == runs[1].sums.uu;
    out("the two g types: all per-point sums %s\n", sameBits ? "bitwise identical" : "DIFFER");
    const OwnRun& r = runs[0];

    // H0 = 0 and P_k = lambda_k at the points
    double maxH0 = 0;
    ZStats zP;
    BandStats bP;
    long long nNearFace = 0;
    const double stop = ownStopDistance(), band = kCompareBand;
    double pE2 = 0, pV = 0, pW = 0, pMax = 0;   // P_k - lambda_k over all points: quadrature-weighted e^2, variance
    for (size_t i = 0; i < ps.size(); i++) {
        double H, var;
        long long n;
        comboH(r.sums, i, unitCoef(5, 0), H, var, n);
        maxH0 = std::max(maxH0, std::fabs(H));
        const bool inBand = dirichletDist(allDirichlet(), ps.xi[i]) < band;
        for (int k = 0; k < 4; k++) {
            comboH(r.sums, i, unitCoef(5, 1 + k), H, var, n);
            const double e = H - lambdaOf(A, k, ps.X[i]);
            if (inBand) bP.add(e, var);
            else zP.add(i, e, var);
            pE2 += ps.w[i]*e*e;
            pV += ps.w[i]*var;
            pW += ps.w[i];
            pMax = std::max(pMax, std::fabs(e));
        }
        if (referenceTet().distToFaces(ps.xi[i]) < stop) nNearFace++;
    }
    const double pRms = std::sqrt(pE2/pW), pNoise = std::sqrt(pV/pW);
    out("H0: max |H0| = %.1e over the points (every walk ends on a face with value 0)\n", maxH0);
    out("P_k vs lambda_k at the points farther than %.2f from every face: |z| < 3 at %.2f%% of %lld (point, k), pooled z =\n"
        "%+.2f, mean z = %+.4f, max |z| = %.2f; %lld of them with zero variance (max |P - lambda| there %.1e); %lld (point,\n"
        "k) within %.2f of a face left out (max |P - lambda| there %.1e)\n",
        band, 100*zP.frac(), zP.n, zP.pooled(), zP.mean(), zP.maxAbsZ, zP.zeroVar, zP.maxAbsDiffZeroVar, bP.n, band,
        bP.maxAbsDiff);
    out("points within the walk termination shell (%.2e) of a face: %lld (their walks end at step 0, so their gradient\n"
        "estimate is 0)\n", stop, nNearFace);

    // S of the P_k (columns 1..4 of the run) on a weight vector: the region tests below and Zombie's point set
    std::vector<std::vector<double>> coefs;
    for (int k = 0; k < 4; k++) coefs.push_back(unitCoef(5, 1 + k));
    auto gramOn = [&](const OwnRun& run, const std::vector<double>& w) {
        return gramFromBatches(w, run.sums.nBatches, 4, [&](size_t p, int b, Vec3* g) {
            for (int q = 0; q < 4; q++) g[q] = A.JinvT*comboBatchGrad(run.sums, p, b, coefs[q]);
        });
    };
    // P_k = lambda_k, so grad P_i . grad P_j is constant: S over any subset of the points is K (its weight)/V
    const double V = std::fabs(A.detJ)/6.0;
    auto maxRatioTo = [&](const Gram& g, double vol) {
        double m = 0;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) m = std::max(m, std::fabs(g.G(i, j) - K(i, j)*vol/V)/g.se(i, j));
        return m;
    };
    const Gram gm = gramOwn(r, ps, A, coefs);
    Eigen::VectorXd ev, evSe;
    eigenWithSe(gm, gm.G, gm.jack, ev, evSe);
    double maxRatio = 0, maxDiff = 0, maxRow = 0, maxRowSe = 0, maxNaiveDiff = 0;
    for (int i = 0; i < 4; i++) {
        double row = 0;
        std::vector<double> rowJ(size_t(gm.B), 0.0);
        for (int j = 0; j < 4; j++) {
            row += gm.G(i, j);
            for (int a = 0; a < gm.B; a++) rowJ[a] += gm.jack[a](i, j);
            maxDiff = std::max(maxDiff, std::fabs(gm.G(i, j) - K(i, j)));
            maxRatio = std::max(maxRatio, std::fabs(gm.G(i, j) - K(i, j))/gm.se(i, j));
            maxNaiveDiff = std::max(maxNaiveDiff, std::fabs(gm.naive(i, j) - gm.G(i, j)));
        }
        double m = 0, v = 0;
        for (double x : rowJ) m += x/double(gm.B);
        for (double x : rowJ) v += (x - m)*(x - m);
        maxRow = std::max(maxRow, std::fabs(row));
        maxRowSe = std::max(maxRowSe, std::sqrt(v*double(gm.B - 1)/double(gm.B)));
    }
    out("\nS (own loop, split-sample)            K                                   (S - K)/s.e.\n");
    for (int i = 0; i < 4; i++) {
        std::string zr;
        for (int j = 0; j < 4; j++) zr += fmtd("  %+6.2f", (gm.G(i, j) - K(i, j))/gm.se(i, j));
        out("  %s    %s   %s\n", matRow(gm.G, i).c_str(), matRow(K, i).c_str(), zr.c_str());
    }
    out("s.e. of S (jackknife over %d batches):\n", gm.B);
    for (int i = 0; i < 4; i++) out("  %s\n", matRow(gm.se, i, "%.5f").c_str());
    out("naive S (product of the full means) - split S:\n");
    for (int i = 0; i < 4; i++) out("  %s\n", matRow(gm.naive - gm.G, i, "%+.2e").c_str());
    out("two-halves S - split S (the halves estimator, batches 0..%d vs %d..%d):\n", gm.B/2 - 1, gm.B/2, gm.B - 1);
    for (int i = 0; i < 4; i++) out("  %s\n", matRow(gm.halves - gm.G, i, "%+.2e").c_str());
    out("eigenvalues of S: %+.5f %+.5f %+.5f %+.5f (s.e. %.5f %.5f %.5f %.5f)\n",
        ev[0], ev[1], ev[2], ev[3], evSe[0], evSe[1], evSe[2], evSe[3]);
    {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> ek(K);
        out("   exact eigenvalues of K: %.5f %.5f %.5f %.5f\n", ek.eigenvalues()[0], ek.eigenvalues()[1], ek.eigenvalues()[2],
            ek.eigenvalues()[3]);
    }
    out("row sums of S: max |sum_j S_ij| = %.1e (jackknife s.e. %.1e). By construction: sum_k lambda_k = 1 on every\n"
        "face (and the extensions sum to 1), so sum_k P_k = 1 walk by walk and the row sums and the zero eigenvalue hold\n"
        "to round-off for any walks; they only confirm that the four hats sum to 1\n", maxRow, maxRowSe);
    const bool passSK = maxRatio < 3;
    const bool passRow = maxRow <= std::max(1e-12, 3*maxRowSe);
    const bool passPsd = ev[0] >= -std::max(1e-12, 3*evSe[0]);
    out("max |S - K| = %.5f, max |S - K|/s.e. = %.2f (< 3: %s); row sums ~ 0: %s; S PSD (min eigenvalue >= -3 s.e.): %s;\n"
        "max |naive - split| = %.2e\n", maxDiff, maxRatio, passSK ? "yes" : "NO", passRow ? "yes" : "NO", passPsd ? "yes" : "NO",
        maxNaiveDiff);
    // The same test by region: the points within / beyond kCompareBand of a face, S against K (region weight)/V.
    // The far region holds no point next to a face, so its s.e. does not depend on how the first ball behaves there.
    std::vector<double> wNear(ps.w), wFar(ps.w);
    double vNear = 0, vFar = 0;
    for (size_t i = 0; i < ps.size(); i++) {
        const bool nearFace = dirichletDist(allDirichlet(), ps.xi[i]) < kCompareBand;
        (nearFace ? wFar : wNear)[i] = 0;
        (nearFace ? vNear : vFar) += ps.w[i];
    }
    const Gram gn = gramOn(r, wNear), gf = gramOn(r, wFar);
    const double ratioFar = maxRatioTo(gf, vFar), ratioNear = maxRatioTo(gn, vNear);
    const bool passFar = ratioFar < 3;
    {
        std::string sn, sf;
        for (int i = 0; i < 4; i++) {
            sn += fmtd(" %.5f", gn.se(i, i));
            sf += fmtd(" %.5f", gf.se(i, i));
        }
        out("by region (S against K vol/V, exact since P_k = lambda_k): points farther than %.2f from every face (volume %.4f):\n"
            "max |S - K vol/V|/s.e. = %.2f (< 3: %s); points within %.2f of a face (volume %.4f): %.2f (reported). Jackknife\n"
            "s.e. of the diagonal: far%s; near%s\n", kCompareBand, vFar, ratioFar, passFar ? "yes" : "NO", kCompareBand, vNear,
            ratioNear, sf.c_str(), sn.c_str());
        md += "By region (S against K vol/V): farther than " + fmtd("%.2f", kCompareBand) + " from every face (volume " +
              fmtd("%.4f", vFar) + ") max abs(S - K vol/V)/s.e. = " + fmtd("%.2f", ratioFar) + "; within it (volume " +
              fmtd("%.4f", vNear) + ") " + fmtd("%.2f", ratioNear) + ". Jackknife s.e. of the diagonal: far" + sf + "; near" +
              sn + ". (Stage 1, with the plain first ball next to the Dirichlet faces, had near-face s.e. up to 0.019 on S_00.)\n\n";
    }
    // Is the jackknife s.e. right? Replicate runs with independent base seeds (the P_k only; the walks do not
    // depend on the data): the spread of S over the replicates against the jackknife s.e. of one run
    const int nRep = gp.quick ? 3 : 8;
    std::vector<Eigen::MatrixXd> repS = {gm.G}, repSe = {gm.se};
    int repPass = passSK ? 1 : 0;
    double repSec = 0;
    {
        std::vector<FieldSpec> pf(fields.begin() + 1, fields.end());
        std::vector<std::vector<double>> pc;
        for (int k = 0; k < 4; k++) pc.push_back(unitCoef(4, k));
        for (int rep = 1; rep < nRep; rep++) {
            const OwnRun rr = runOwn(A, none, allDirichlet(), pf, ps.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7300 + rep));
            repSec += rr.stats.seconds;
            const Gram g = gramOwn(rr, ps, A, pc);
            repS.push_back(g.G);
            repSe.push_back(g.se);
            double m = 0;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++) m = std::max(m, std::fabs(g.G(i, j) - K(i, j))/g.se(i, j));
            repPass += m < 3;
        }
    }
    Eigen::MatrixXd repMean = Eigen::MatrixXd::Zero(4, 4), repSd = Eigen::MatrixXd::Zero(4, 4), repRms = Eigen::MatrixXd::Zero(4, 4);
    for (int k = 0; k < nRep; k++) {
        repMean += repS[k]/double(nRep);
        repRms += repSe[k].cwiseAbs2()/double(nRep);
    }
    for (int k = 0; k < nRep; k++) repSd += (repS[k] - repMean).cwiseAbs2()/double(nRep - 1);
    repSd = repSd.cwiseSqrt();
    repRms = repRms.cwiseSqrt();
    double ratioLo = kInf, ratioHi = 0, repMeanRatio = 0;
    for (int i = 0; i < 4; i++)
        for (int j = i; j < 4; j++) {
            ratioLo = std::min(ratioLo, repSd(i, j)/repRms(i, j));
            ratioHi = std::max(ratioHi, repSd(i, j)/repRms(i, j));
            repMeanRatio = std::max(repMeanRatio, std::fabs(repMean(i, j) - K(i, j))/(repSd(i, j)/std::sqrt(double(nRep))));
        }
    out("s.e. check over %d independent replicate runs (seed tags 7100, 7301..%d; %.1f s): sd of S_ij over the runs /\n"
        "rms jackknife s.e. = %.2f to %.2f over the 10 entries (1 expected; with %d runs the sd itself scatters by ~%.0f%%);\n"
        "the full-domain test passes at %d of %d runs; mean S of the runs vs K: max |mean - K|/(sd/sqrt(%d)) = %.2f\n",
        nRep, 7300 + nRep - 1, repSec, ratioLo, ratioHi, nRep, 100.0/std::sqrt(2.0*(nRep - 1)), repPass, nRep, nRep, repMeanRatio);
    const Eigen::MatrixXd repRatio = repSd.cwiseQuotient(repRms);
    out("sd over the runs / rms jackknife s.e., per entry:\n");
    for (int i = 0; i < 4; i++) out("  %s\n", matRow(repRatio, i, "%.2f").c_str());
    md += "s.e. check over " + std::to_string(nRep) + " independent replicate runs: sd of S_ij over the runs / rms jackknife "
          "s.e. = " + fmtd("%.2f", ratioLo) + " to " + fmtd("%.2f", ratioHi) + " (1 expected, the sd of " + std::to_string(nRep) +
          " runs scatters by ~" + fmtd("%.0f", 100.0/std::sqrt(2.0*(nRep - 1))) + "%); the full-domain test passes at " +
          std::to_string(repPass) + " of " + std::to_string(nRep) + " runs; max abs(mean S - K)/(sd/sqrt(" +
          std::to_string(nRep) + ")) = " + fmtd("%.2f", repMeanRatio) + ". Per entry:\n\n" + codeBlock([&] {
        std::string s;
        for (int i = 0; i < 4; i++) s += matRow(repRatio, i, "%.2f") + "\n";
        return s;
    }()) + "\n";
    md += "Row sums and the zero eigenvalue hold by construction (sum_k lambda_k = 1 on every face, so sum_k P_k = 1 walk "
          "by walk); they only confirm that the hats sum to 1.\n\n";

    // Zombie cross-check: the same fields, gradient at every point, nZombieBatches seeds of nZombie/nZombieBatches
    // walks. The same seed for every field in a batch, so the four P_k share their walks, as in the own loop.
    out("\nZombie cross-check (Zombie works in float): N = %d per point in %d independent solves (seeds) of %d walks\n",
        gp.nZombie, gp.nZombieBatches, gp.nZombie/gp.nZombieBatches);
    std::vector<uint64_t> seeds;
    for (int b = 0; b < gp.nZombieBatches; b++) seeds.push_back(mixSeed(gp.seed, 7150 + b));
    std::vector<ZombieRun> zr;
    double zSec = 0;
    double zH0 = 0;
    for (int m = 0; m < 2; m++) {
        const ZombieRun h0 = runZombie(A, none, allDirichlet(), FieldSpec{"H0", true, gts[m], {}}, ps.xi,
                                       gp.nZombie/gp.nZombieBatches, false, {seeds[0]});
        for (const Estimate& e : h0.est[0]) zH0 = std::max(zH0, double(std::fabs(e.H)));
        zSec += h0.stats.seconds;
        out("  Zombie H0, %s: %s\n", gTypeName(gts[m]), h0.stats.summary().c_str());
    }
    for (int k = 0; k < 4; k++) {
        zr.push_back(runZombie(A, none, allDirichlet(),
                               FieldSpec{"P", false, GType::One, [&A, k](const Vec3& X) { return lambdaOf(A, k, X); }}, ps.xi,
                               gp.nZombie/gp.nZombieBatches, true, seeds));
        zSec += zr.back().stats.seconds;
        out("  Zombie P%d: %s\n", k, zr.back().stats.summary().c_str());
    }
    out("  Zombie H0 (t = 0 with Dirichlet faces runs walks now): max |H0| = %.1e\n", zH0);
    // Zombie's gradient is NaN at points within its face-offset margin of a face (its first sphere radius
    // is 0 there): such points are left out of Zombie's S (their weights set to 0) and counted
    std::vector<double> wZ(ps.w);
    long long nNaN = 0;
    double wNaN = 0, dNaN = 0;
    for (size_t p = 0; p < ps.size(); p++) {
        bool finite = true;
        for (int q = 0; q < 4; q++)
            for (int b = 0; b < gp.nZombieBatches; b++) finite = finite && zr[q].est[b][p].grad.allFinite();
        if (!finite) {
            nNaN++;
            wNaN += wZ[p];
            dNaN = std::max(dNaN, dirichletDist(allDirichlet(), ps.xi[p]));
            wZ[p] = 0;
        }
    }
    const Gram gz = gramFromBatches(wZ, gp.nZombieBatches, 4, [&](size_t p, int b, Vec3* g) {
        for (int q = 0; q < 4; q++) g[q] = wZ[p] == 0 ? Vec3::Zero() : Vec3(A.JinvT*zr[q].est[b][p].grad.cast<double>());
    });
    out("  Zombie's gradient is not finite at %lld points, at most %.1e from a face (its face-offset margin is %.1e); they\n"
        "  are left out of Zombie's S (total weight %.2e of %.4f)\n", nNaN, dNaN, faceOffsetMargin(referenceTet()), wNaN, 1.0/6.0);
    // Zombie's H against the own loop's, per point, outside kCompareBand. From separate H-only Zombie solves (one
    // per P_k, own seed): in the gradient solves above Zombie's antithetic pairs replay the same random numbers
    // after the mirrored first step, so their two H values are correlated and Zombie's s.e. of H is too small
    // (V0 of --validate-reuse)
    ZStats zHo;
    BandStats bHo;
    const double bandZ = kCompareBand;
    for (int k = 0; k < 4; k++) {
        const ZombieRun zh = runZombie(A, none, allDirichlet(),
                                       FieldSpec{"P", false, GType::One, [&A, k](const Vec3& X) { return lambdaOf(A, k, X); }},
                                       ps.xi, gp.nZombie, false, {mixSeed(gp.seed, 7170 + k)});
        zSec += zh.stats.seconds;
        for (size_t i = 0; i < ps.size(); i++) {
            const double Hz = zh.est[0][i].H, v2 = double(zh.est[0][i].Hse)*zh.est[0][i].Hse;
            double H, var;
            long long n;
            comboH(r.sums, i, unitCoef(5, 1 + k), H, var, n);
            if (dirichletDist(allDirichlet(), ps.xi[i]) < bandZ) bHo.add(H - Hz, var + v2);
            else zHo.add(i, H - Hz, var + v2);
        }
    }
    // S_Z against the own S over the same points (Zombie's NaN points left out of both) and against K (vol)/V
    const Gram gmZ = gramOn(r, wZ);
    double vZ = 0;
    for (double w : wZ) vZ += w;
    double maxZO = 0;
    const double maxZK = maxRatioTo(gz, vZ);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            maxZO = std::max(maxZO, std::fabs(gz.G(i, j) - gmZ.G(i, j))/std::hypot(gz.se(i, j), gmZ.se(i, j)));
    out("S (Zombie, split-sample over its %d seeds)    s.e. (jackknife)\n", gz.B);
    for (int i = 0; i < 4; i++) out("  %s    %s\n", matRow(gz.G, i).c_str(), matRow(gz.se, i, "%.5f").c_str());
    const bool zPass = maxZO < 3 && zHo.passPooled();
    out("Zombie vs own (the own S over the same points): max |S_Z - S_own|/combined s.e. = %.2f; max |S_Z - K vol/V|/s.e. =\n"
        "%.2f; H of P_k (H-only Zombie solves, N = %d) per point farther than %.2f from the faces: |z| < 3 at %.2f%% of %lld,\n"
        "pooled z = %+.2f (mean dH = %+.1e +- %.1e), mean z = %+.4f -> %s (agree: S within 3 combined s.e., |z| < 3 at >= 99%%\n"
        "and |pooled z| < 3); %lld (point, k) within it, max |dH| there %.1e\n",
        maxZO, maxZK, gp.nZombie, bandZ, 100*zHo.frac(), zHo.n, zHo.pooled(), zHo.meanDiff(), zHo.meanDiffSe(), zHo.mean(),
        zPass ? "agree" : "DISAGREE", bHo.n, bHo.maxAbsDiff);
    out("wall clock: own loop %.2f s per run (%.2f s both g types), Zombie %.2f s in all\n", r.stats.seconds,
        runs[0].stats.seconds + runs[1].stats.seconds, zSec);
    // the verdict: the spec's S = K test over A, the same test far from the faces, the by-construction checks,
    // and the independent cross-check (the lead's note: Zombie must agree)
    const bool pass = passSK && passFar && passRow && passPsd && sameBits && zPass;
    out("CHECK 0: %s (S = K over A: %s; far from the faces: %s; row sums / PSD: %s; g modes identical: %s; Zombie: %s)\n",
        pass ? "PASS" : "FAIL", passSK ? "pass" : "FAIL", passFar ? "pass" : "FAIL", passRow && passPsd ? "pass" : "FAIL",
        sameBits ? "yes" : "NO", zPass ? "agrees" : "DISAGREES");

    // report
    MdTable t;
    t.head = {"g mode", "solver", "max abs(S-K)", "max abs(S-K)/s.e.", "far from faces: max abs(S-K vol/V)/s.e.",
              "max abs(row sum)", "eigenvalues of S (s.e.)", "max abs(naive-split)", "wall clock (s)", "verdict"};
    for (int m = 0; m < 2; m++) {
        std::string evs;
        for (int k = 0; k < 4; k++) evs += (k ? ", " : "") + fmtd("%.4f", ev[k]) + " (" + fmtd("%.4f", evSe[k]) + ")";
        t.rows.push_back({gTypeName(gts[m]), "own loop (double), N = " + std::to_string(gp.nWalks), fmtd("%.5f", maxDiff),
                          fmtd("%.2f", maxRatio), fmtd("%.2f", ratioFar), fmtd("%.1e", maxRow) + " (by construction)", evs,
                          fmtd("%.2e", maxNaiveDiff), fmtd("%.2f", runs[m].stats.seconds), pass ? "PASS" : "FAIL"});
    }
    {
        Eigen::VectorXd evz, evzSe;
        eigenWithSe(gz, gz.G, gz.jack, evz, evzSe);
        std::string evs;
        for (int k = 0; k < 4; k++) evs += (k ? ", " : "") + fmtd("%.4f", evz[k]) + " (" + fmtd("%.4f", evzSe[k]) + ")";
        double maxDZ = 0, maxRowZ = 0;
        for (int i = 0; i < 4; i++) {
            maxRowZ = std::max(maxRowZ, std::fabs(gz.G.row(i).sum()));
            for (int j = 0; j < 4; j++) maxDZ = std::max(maxDZ, std::fabs(gz.G(i, j) - K(i, j)));
        }
        t.rows.push_back({"both (g does not enter)", "Zombie (float), N = " + std::to_string(gp.nZombie), fmtd("%.5f", maxDZ),
                          fmtd("%.2f", maxZK) + " (vs K vol/V)", "-", fmtd("%.1e", maxRowZ), evs,
                          fmtd("%.2e", (gz.naive - gz.G).cwiseAbs().maxCoeff()), fmtd("%.2f", zSec),
                          zPass ? "agrees with own" : "DISAGREES"});
    }
    md += t.str() + "\n";
    md += "S (own loop) and its s.e.:\n\n" + codeBlock([&] {
        std::string s;
        for (int i = 0; i < 4; i++) s += matRow(gm.G, i) + "    " + matRow(gm.se, i, "%.5f") + "\n";
        return s;
    }());
    md += "\nH0: max abs(H0) = " + fmtd("%.1e", maxH0) + "; P_k vs lambda_k per point: " + fmtd("%.2f", 100*zP.frac()) +
          "% with abs(z) < 3, pooled z " + fmtd("%+.2f", zP.pooled()) + " (points farther than " + fmtd("%.2f", band) +
          " from the faces); over all points RMS abs(P_k - lambda_k) " + fmtd("%.1e", pRms) + " (its noise, the RMS s.e., " +
          fmtd("%.1e", pNoise) + "), max " + fmtd("%.1e", pMax) + "; the two g modes give bitwise identical sums: " +
          (sameBits ? "yes" : "NO") + ". Zombie vs own: max abs(S_Z - S_own)/s.e. = " + fmtd("%.2f", maxZO) +
          " (same points); H of P_k: " + fmtd("%.2f", 100*zHo.frac()) + "% with abs(z) < 3, pooled z " +
          fmtd("%+.2f", zHo.pooled()) + "; Zombie works in float.\n\n";
    md += std::string("**Check 0: ") + (pass ? "PASS" : "FAIL") + "**\n";

    // findings
    std::string f = std::string("- Check 0 (machinery: tet A, no crack, all four faces Dirichlet): ") + (pass ? "PASS" : "FAIL") +
                    ". The verdict needs S = K over A within 3 s.e. (the spec), the same over the points farther than " +
                    fmtd("%.2f", kCompareBand) + " from every face, row sums ~ 0 and S PSD, identical g modes, and Zombie agreeing.\n";
    f += "- S against K: max abs(S - K) = " + fmtd("%.5f", maxDiff) + ", at most " + fmtd("%.2f", maxRatio) + " s.e. (" +
         fmtd("%.2f", ratioFar) + " far from the faces). Eigenvalues of S " + fmtd("%.4f", ev[0]) + ", " + fmtd("%.4f", ev[1]) +
         ", " + fmtd("%.4f", ev[2]) + ", " + fmtd("%.4f", ev[3]) + " (K: 0, 1/6, 1/6, 2/3), so S is PSD; the row sums (" +
         fmtd("%.1e", maxRow) + ") and the zero eigenvalue hold by construction.\n";
    f += "- The s.e. are calibrated: over " + std::to_string(nRep) + " independent runs the sd of S_ij is " + fmtd("%.2f", ratioLo) +
         " to " + fmtd("%.2f", ratioHi) + " times the jackknife s.e. (an sd from " + std::to_string(nRep) + " runs scatters by ~" +
         fmtd("%.0f", 100.0/std::sqrt(2.0*(nRep - 1))) + "%), and the S = K test passes at " + std::to_string(repPass) + " of " +
         std::to_string(nRep) + ".\n";
    f += "- H0 = 0 at every point (max " + fmtd("%.1e", maxH0) + "). P_k = lambda_k: RMS abs(P_k - lambda_k) " + fmtd("%.1e", pRms) +
         " against its noise " + fmtd("%.1e", pNoise) + "; per point abs(z) < 3 at " + fmtd("%.2f", 100*zP.frac()) + "%, pooled z " +
         fmtd("%+.2f", zP.pooled()) + ".\n";
    f += std::string("- Zombie (float, independent) ") + (zPass ? "agrees" : "DISAGREES") + ": max abs(S_Z - S_own) = " +
         fmtd("%.2f", maxZO) + " combined s.e.; H of the P_k pooled z " + fmtd("%+.2f", zHo.pooled()) + ".\n";
    f += "- Wall clock per own-loop run (" + std::to_string(ps.size()) + " points, N = " + std::to_string(gp.nWalks) + ", 5 fields): " +
         fmtd("%.2f", r.stats.seconds) + " s.\n";
    md.insert(mdIntro, findingsBlock(f));
    return pass;
}

// =====================================================================================
// Check 1: tet A alone, crack on, F Dirichlet, the other faces Neumann: superposition
// =====================================================================================
bool runCheck1(Report& out, const GateParams& gp, const PointSet& ps, const XiCrack& xa, std::string& md)
{
    const TetFrame A = frameA();
    out("\n=== Check 1: superposition. Tet A alone, crack x = 0.30, y_f = 0.30, F Dirichlet, other faces Neumann ===\n");
    out("recipe run: H0 (crack +-g, F 0) for both g types and P_k (crack 0, F hat_k), k = v0, v1, v2, all from the same\n"
        "walks; direct runs: crack +-g and F sum_k c_k hat_k, independent seeds. %zu points, N = %d in %d batches\n",
        ps.size(), gp.nWalks, gp.nBatches);
    md += "## Check 1: superposition (tet A, crack x = 0.30, y_f = 0.30, F Dirichlet, other faces Neumann)\n\n";
    md += "Recipe run: H0 (both g modes) and P_0..P_2 from the same walks. Direct runs: F data sum_k c_k hat_k, crack +-g, "
          "independent seeds. z = (H_direct - H_combo)/sqrt(var_direct + var_combo), var_combo = per-walk variance of the "
          "combined value / N. Spec criterion: abs(z) < 3 at >= 99% of points and abs(mean z) < 3/sqrt(M); a guard fails it "
          "if mean abs(z) < 0.5 (0.80 for independent runs; ~0 if the direct run reused the recipe's walks). Own loop "
          "(double), " + std::to_string(ps.size()) + " points, N = " + std::to_string(gp.nWalks) + ". There is no reference "
          "solve here; the error columns compare the direct solve with the recombined fields, divided by 2 g_max = 2.\n\n";
    const size_t mdIntro = md.size();   // the findings go here
    md += "**What this criterion can and cannot show.** H0, the P_k and the direct runs share the walk process (same "
          "geometry and face flags, the lead's note 6), and each walk's direct value is exactly sum_k coef_k u_k of its "
          "recipe values. So H_direct and H_combo are two independent samples of the same random variable and z ~ N(0,1) "
          "whatever the face code does: a wrong face flagged, or the hats evaluated in the wrong frame, would pass too "
          "(reviewer's harness). The own-loop test only confirms that the scoring is linear and the seed streams "
          "independent. The Dirichlet-face path is validated by sanity (a) and (c), check 0 (S = K and P_k = lambda_k with "
          "all faces Dirichlet) and, here with the crack, by the independent Zombie cross-check below, which is part of "
          "the verdict.\n\n";
    const GType gts[2] = {GType::One, GType::SqrtFront};
    std::vector<FieldSpec> fields = {FieldSpec{"H0 g1", true, GType::One, {}}, FieldSpec{"H0 g2", true, GType::SqrtFront, {}}};
    for (int g = 0; g < 3; g++) fields.push_back(FieldSpec{"P" + std::to_string(g), false, GType::One, [g](const Vec3& X) { return hatF(g, X); }});
    const OwnRun rec = runOwn(A, xa, onlyF(), fields, ps.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7200));
    out("recipe run: %s\n", rec.stats.summary().c_str());

    std::mt19937_64 rng(mixSeed(gp.seed, 7210));
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    bool pass = true;
    MdTable t;
    t.head = {"g mode", "c", "M", "frac abs(z)<3", "mean z", "3/sqrt(M)", "max abs(z)", "mean abs(z)",
              "H_direct - H_combo: RMS / max / 2 g_max (noise RMS)", "wall clock direct (s)", "verdict"};
    double fracLo = kInf, meanZHi = 0, mAbsLo = kInf, mAbsHi = 0, eLo = kInf, eHi = 0, nLo = kInf, nHi = 0;   // over the 6 tests
    std::vector<std::array<double, 3>> cs;
    std::vector<OwnRun> directs;
    for (int rIdx = 0; rIdx < 3; rIdx++) {
        std::array<double, 3> c{U(rng), U(rng), U(rng)};
        cs.push_back(c);
        auto faceFn = [c](const Vec3& X) { return c[0]*hatF(0, X) + c[1]*hatF(1, X) + c[2]*hatF(2, X); };
        const std::vector<FieldSpec> df = {FieldSpec{"D g1", true, GType::One, faceFn}, FieldSpec{"D g2", true, GType::SqrtFront, faceFn}};
        directs.push_back(runOwn(A, xa, onlyF(), df, ps.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7220 + rIdx)));
        const OwnRun& d = directs.back();
        out("c%d = (%+.4f, %+.4f, %+.4f): direct run %s\n", rIdx + 1, c[0], c[1], c[2], d.stats.summary().c_str());
        for (int m = 0; m < 2; m++) {
            ZStats z;
            const std::vector<double> coef = {m == 0 ? 1.0 : 0.0, m == 1 ? 1.0 : 0.0, c[0], c[1], c[2]};
            double e2 = 0, ev = 0, ew = 0, emax = 0;
            for (size_t i = 0; i < ps.size(); i++) {
                double Hc, vc, Hd, vd;
                long long nc, nd;
                comboH(rec.sums, i, coef, Hc, vc, nc);
                comboH(d.sums, i, unitCoef(2, m), Hd, vd, nd);
                z.add(i, Hd - Hc, vd + vc);
                e2 += ps.w[i]*(Hd - Hc)*(Hd - Hc);
                ev += ps.w[i]*(vd + vc);
                ew += ps.w[i];
                emax = std::max(emax, std::fabs(Hd - Hc));
            }
            const double eRms = std::sqrt(e2/ew)/kTwoGMax, eNoise = std::sqrt(ev/ew)/kTwoGMax;
            // guard: with shared walks (seeds) every difference is round-off, z ~ 0, and the criterion passes trivially
            const bool degenerate = z.meanAbs() < 0.5;
            const bool ok = z.pass() && !degenerate;
            pass = pass && ok;
            out("   %s: |z| < 3 at %.2f%% of %lld points, mean z = %+.4f (bound %.4f), max |z| = %.2f at %s, mean |z| = %.3f%s\n"
                "      -> %s\n",
                gTypeName(gts[m]), 100*z.frac(), z.n, z.mean(), z.meanBound(), z.maxAbsZ, fmt3(ps.X[z.worst]).c_str(),
                z.meanAbs(), degenerate ? " (DEGENERATE: the direct run shares the recipe's walks?)" : "", ok ? "PASS" : "FAIL");
            t.rows.push_back({gTypeName(gts[m]), "(" + fmtd("%+.3f", c[0]) + ", " + fmtd("%+.3f", c[1]) + ", " + fmtd("%+.3f", c[2]) + ")",
                              std::to_string(z.n), fmtd("%.4f", z.frac()), fmtd("%+.4f", z.mean()), fmtd("%.4f", z.meanBound()),
                              fmtd("%.2f", z.maxAbsZ), fmtd("%.3f", z.meanAbs()),
                              fmtd("%.1e", eRms) + " / " + fmtd("%.1e", emax/kTwoGMax) + " (" + fmtd("%.1e", eNoise) + ")",
                              fmtd("%.2f", d.stats.seconds), ok ? "PASS" : degenerate ? "FAIL (degenerate)" : "FAIL"});
            fracLo = std::min(fracLo, z.frac());
            meanZHi = std::max(meanZHi, std::fabs(z.mean()));
            mAbsLo = std::min(mAbsLo, z.meanAbs());
            mAbsHi = std::max(mAbsHi, z.meanAbs());
            eLo = std::min(eLo, eRms);
            eHi = std::max(eHi, eRms);
            nLo = std::min(nLo, eNoise);
            nHi = std::max(nHi, eNoise);
        }
    }
    out("(this criterion confirms only the linearity of the scoring: the face path is checked by the Zombie cross-check below)\n");
    md += t.str() + "\n";

    // S and b of A's Level-A recipe (reused by check 3 at y_f = 0.30)
    Eigen::VectorXd evS, evSSe;
    {
        std::vector<std::vector<double>> coefs;
        for (int k = 0; k < 5; k++) coefs.push_back(unitCoef(5, k));
        const Gram gm = gramOwn(rec, ps, A, coefs);
        Eigen::MatrixXd S = gm.G.block(2, 2, 3, 3);
        std::vector<Eigen::MatrixXd> jS;
        for (const auto& J : gm.jack) jS.push_back(J.block(2, 2, 3, 3));
        Eigen::VectorXd ev, evSe;
        eigenWithSe(gm, S, jS, ev, evSe);
        evS = ev;
        evSSe = evSe;
        out("\nS and b of A's recipe (F Dirichlet, crack on; check 3 reuses them), split-sample, jackknife s.e.:\n");
        for (int i = 0; i < 3; i++)
            out("  S  %s   s.e. %s   naive-split %s\n", matRow(S, i).c_str(), matRow(gm.se.block(2, 2, 3, 3), i, "%.5f").c_str(),
                matRow(gm.naive.block(2, 2, 3, 3) - S, i, "%+.1e").c_str());
        for (int m = 0; m < 2; m++)
            out("  b (%s) = %+.5f %+.5f %+.5f   s.e. %.5f %.5f %.5f;  int |grad H0|^2 = %.5f +- %.5f\n", gTypeName(gts[m]),
                gm.G(m, 2), gm.G(m, 3), gm.G(m, 4), gm.se(m, 2), gm.se(m, 3), gm.se(m, 4), gm.G(m, m), gm.se(m, m));
        out("  eigenvalues of S: %+.5f %+.5f %+.5f (s.e. %.5f %.5f %.5f)\n", ev[0], ev[1], ev[2], evSe[0], evSe[1], evSe[2]);
        md += "S and b of A's recipe (split-sample, jackknife s.e.):\n\n" + codeBlock([&] {
            std::string s;
            for (int i = 0; i < 3; i++) s += "S  " + matRow(S, i) + "    s.e. " + matRow(gm.se.block(2, 2, 3, 3), i, "%.5f") + "\n";
            for (int m = 0; m < 2; m++)
                s += std::string("b ") + (m ? "(b) " : "(a) ") + fmtd("%+.5f", gm.G(m, 2)) + " " + fmtd("%+.5f", gm.G(m, 3)) + " " +
                     fmtd("%+.5f", gm.G(m, 4)) + "    s.e. " + fmtd("%.5f", gm.se(m, 2)) + " " + fmtd("%.5f", gm.se(m, 3)) + " " +
                     fmtd("%.5f", gm.se(m, 4)) + "\n";
            s += "eigenvalues of S: " + fmtd("%.5f", ev[0]) + " " + fmtd("%.5f", ev[1]) + " " + fmtd("%.5f", ev[2]) + "  (s.e. " +
                 fmtd("%.5f", evSe[0]) + " " + fmtd("%.5f", evSe[1]) + " " + fmtd("%.5f", evSe[2]) + ")\n";
            return s;
        }());
    }

    // Zombie cross-check: every field as its own Zombie solve (H only), independent seeds. Two different estimators,
    // so the criterion is ZStats::passPooled (|z| < 3 at >= 99% of the points and |pooled z| < 3; mean z is biased
    // by the skewness of the walk values when the two N differ, see ZStats), over the points farther than
    // kCompareBand from F, and also |pooled z| < 3 over the points near the segment where the crack meets F, the
    // junction of the two Dirichlet boundaries (the path this check is about).
    std::array<Vec3, 2> junction{};
    {
        int nj = 0;
        for (const Vec3& p : xa.cut.polygon)
            if (std::fabs(p.z()) <= 1e-12 && nj < 2) junction[nj++] = A.toX(p);   // F is xi_z = 0
        if (nj < 2) throw std::logic_error("check 1: the crack does not meet F in a segment");
    }
    constexpr double kJunctionBand = 0.15;
    auto nearJunction = [&](size_t i) { return pointSegmentDist(ps.X[i], junction[0], junction[1]) < kJunctionBand; };
    const double bandZ = kCompareBand;
    long long nBand = 0, nJunction = 0;
    for (size_t i = 0; i < ps.size(); i++) {
        const bool inBand = dirichletDist(onlyF(), ps.xi[i]) < bandZ;
        nBand += inBand;
        nJunction += !inBand && nearJunction(i);
    }
    // Skew control: the recipe again with the own loop at Zombie's N and an independent seed. Own (N = nWalks)
    // against own (N = nZombie) shows the mean-z offset that unequal N and skewed walk values alone produce.
    const OwnRun recSmall = runOwn(A, xa, onlyF(), fields, ps.xi, gp.nZombie, gp.nBatches, mixSeed(gp.seed, 7400));
    out("\nZombie cross-check (Zombie works in float): each field one Zombie solve of N = %d, H only, independent seeds.\n"
        "Per point own vs Zombie over the %lld points farther than %.2f from F (%lld within it left out: Zombie's stop\n"
        "distance %.1e and the rare-excursion skew next to a Dirichlet face), and over the %lld of them within %.2f of the\n"
        "segment where the crack meets F, %s - %s. Agree: |z| < 3 at >= 99%% and |pooled z| < 3 overall, |pooled z| < 3\n"
        "near the crack-F segment. pooled z = sum(dH)/sqrt(sum var), unbiased; mean z is shown, but between N = %d and\n"
        "N = %d it is offset by the walk values' skewness: the skew control (own loop at N = %d vs N = %d, the same\n"
        "solver: %s) shows that offset without any solver difference.\n",
        gp.nZombie, (long long)ps.size() - nBand, bandZ, nBand, zombieStopDistance(), nJunction, kJunctionBand,
        fmt3(junction[0], 2).c_str(), fmt3(junction[1], 2).c_str(), gp.nWalks, gp.nZombie, gp.nWalks, gp.nZombie,
        recSmall.stats.summary().c_str());
    MdTable tz;
    tz.head = {"field", "frac abs(z)<3", "pooled z (mean dH +- s.e.)", "mean z", "near crack-F: pooled z, mean z",
               "skew control own-own: mean z, pooled z", "Zombie wall clock (s)", "verdict"};
    bool zPass = true;
    double zpLo = kInf, zpHi = -kInf, zjLo = kInf, zjHi = -kInf, zfLo = kInf;   // over the fields (findings)
    int nZ = 0;
    auto zRange = [&](const ZStats& all, const ZStats& jn) {
        zpLo = std::min(zpLo, all.pooled());
        zpHi = std::max(zpHi, all.pooled());
        zjLo = std::min(zjLo, jn.pooled());
        zjHi = std::max(zjHi, jn.pooled());
        zfLo = std::min(zfLo, all.frac());
        nZ++;
    };
    std::vector<std::vector<Estimate>> zf;
    // own (coefficients coef of run r) against estimates e at the points: overall and near the crack-F segment
    auto compare = [&](const OwnRun& r, const std::vector<double>& coef, const std::vector<Estimate>& e, ZStats& all, ZStats& jn) {
        for (size_t i = 0; i < ps.size(); i++) {
            if (dirichletDist(onlyF(), ps.xi[i]) < bandZ) continue;
            double H, var;
            long long n;
            comboH(r.sums, i, coef, H, var, n);
            const double d = H - e[i].H, v = var + double(e[i].Hse)*e[i].Hse;
            all.add(i, d, v);
            if (nearJunction(i)) jn.add(i, d, v);
        }
    };
    auto agree = [](const ZStats& all, const ZStats& jn) { return all.passPooled() && std::fabs(jn.pooled()) < 3.0; };
    for (size_t k = 0; k < fields.size(); k++) {
        const ZombieRun z = runZombie(A, xa, onlyF(), fields[k], ps.xi, gp.nZombie, false, {mixSeed(gp.seed, 7240 + k)});
        zf.push_back(z.est[0]);
        ZStats zs, zj, sk;
        compare(rec, unitCoef(5, int(k)), z.est[0], zs, zj);
        for (size_t i = 0; i < ps.size(); i++) {   // skew control: own N = nWalks vs own N = nZombie
            if (dirichletDist(onlyF(), ps.xi[i]) < bandZ) continue;
            double H1, v1, H2, v2;
            long long n;
            comboH(rec.sums, i, unitCoef(5, int(k)), H1, v1, n);
            comboH(recSmall.sums, i, unitCoef(5, int(k)), H2, v2, n);
            sk.add(i, H1 - H2, v1 + v2);
        }
        const bool ok = agree(zs, zj);
        zPass = zPass && ok;
        zRange(zs, zj);
        out("  %-6s own vs Zombie: |z| < 3 at %.2f%%, pooled z = %+.2f (mean dH %+.1e +- %.1e), mean z = %+.4f; near crack-F:\n"
            "         pooled z = %+.2f, mean z = %+.4f (bound %.3f); skew control: mean z = %+.4f, pooled z = %+.2f -> %s; Zombie %s\n",
            fields[k].name.c_str(), 100*zs.frac(), zs.pooled(), zs.meanDiff(), zs.meanDiffSe(), zs.mean(), zj.pooled(), zj.mean(),
            zj.meanBound(), sk.mean(), sk.pooled(), ok ? "agree" : "DISAGREE", z.stats.summary().c_str());
        tz.rows.push_back({fields[k].name, fmtd("%.4f", zs.frac()),
                           fmtd("%+.2f", zs.pooled()) + " (" + fmtd("%+.1e", zs.meanDiff()) + " +- " + fmtd("%.1e", zs.meanDiffSe()) + ")",
                           fmtd("%+.4f", zs.mean()), fmtd("%+.2f", zj.pooled()) + ", " + fmtd("%+.4f", zj.mean()),
                           fmtd("%+.4f", sk.mean()) + ", " + fmtd("%+.2f", sk.pooled()), fmtd("%.2f", z.stats.seconds),
                           ok ? "agree" : "DISAGREE"});
    }
    for (int rIdx = 0; rIdx < 3; rIdx++) {
        const auto c = cs[rIdx];
        auto faceFn = [c](const Vec3& X) { return c[0]*hatF(0, X) + c[1]*hatF(1, X) + c[2]*hatF(2, X); };
        for (int m = 0; m < 2; m++) {
            const ZombieRun z = runZombie(A, xa, onlyF(), FieldSpec{"D", true, gts[m], faceFn}, ps.xi, gp.nZombie, false,
                                          {mixSeed(gp.seed, 7260 + 2*rIdx + m)});
            ZStats zs, zj, zsup;
            const std::vector<double> coef = {m == 0 ? 1.0 : 0.0, m == 1 ? 1.0 : 0.0, c[0], c[1], c[2]};
            compare(directs[rIdx], unitCoef(2, m), z.est[0], zs, zj);
            for (size_t i = 0; i < ps.size(); i++) {
                // Zombie's own superposition: its direct solve against its combination of independent solves
                double comb = 0, cv = 0;
                for (int k = 0; k < 5; k++) {
                    comb += coef[k]*zf[k][i].H;
                    cv += coef[k]*coef[k]*double(zf[k][i].Hse)*zf[k][i].Hse;
                }
                zsup.add(i, z.est[0][i].H - comb, cv + double(z.est[0][i].Hse)*z.est[0][i].Hse);
            }
            const bool ok = agree(zs, zj);
            zPass = zPass && ok;
            zRange(zs, zj);
            const std::string nm = "direct c" + std::to_string(rIdx + 1) + (m ? " g2" : " g1");
            out("  %-12s own vs Zombie: |z| < 3 at %.2f%%, pooled z = %+.2f, mean z = %+.4f; near crack-F: pooled z = %+.2f\n"
                "               -> %s; Zombie's own superposition (all points): |z| < 3 at %.2f%%, pooled z = %+.2f, mean z =\n"
                "               %+.4f; Zombie %.2f s\n",
                nm.c_str(), 100*zs.frac(), zs.pooled(), zs.mean(), zj.pooled(), ok ? "agree" : "DISAGREE", 100*zsup.frac(),
                zsup.pooled(), zsup.mean(), z.stats.seconds);
            tz.rows.push_back({nm, fmtd("%.4f", zs.frac()),
                               fmtd("%+.2f", zs.pooled()) + " (" + fmtd("%+.1e", zs.meanDiff()) + " +- " + fmtd("%.1e", zs.meanDiffSe()) + ")",
                               fmtd("%+.4f", zs.mean()), fmtd("%+.2f", zj.pooled()) + ", " + fmtd("%+.4f", zj.mean()), "-",
                               fmtd("%.2f", z.stats.seconds),
                               std::string(ok ? "agree" : "DISAGREE") + "; Zombie's superposition: " + fmtd("%.4f", zsup.frac()) +
                                   ", pooled z " + fmtd("%+.2f", zsup.pooled())});
        }
    }
    // Zombie's termination bias, shown on P1: Zombie at the gate's shell (stop distance 6.2e-5) and at the first
    // gate's (2.0e-4), each from 4 independent seeds (averaged per point), compared in H over the points
    std::string termMd;
    {
        const int iP1 = 3;
        std::vector<uint64_t> sg, sf;
        for (int s = 0; s < 4; s++) {
            sg.push_back(mixSeed(gp.seed, 7280 + s));
            sf.push_back(mixSeed(gp.seed, 7290 + s));
        }
        const ZombieRun zg = runZombie(A, xa, onlyF(), fields[iP1], ps.xi, gp.nZombie, false, sg);
        const ZombieRun zf1 = runZombie(A, xa, onlyF(), fields[iP1], ps.xi, gp.nZombie, false, sf, true);
        ZStats ownG, ownF, shell;
        for (size_t i = 0; i < ps.size(); i++) {
            if (dirichletDist(onlyF(), ps.xi[i]) < bandZ) continue;
            double hg = 0, vg = 0, hf = 0, vf = 0;
            for (int s = 0; s < 4; s++) {
                hg += zg.est[s][i].H/4.0;
                vg += double(zg.est[s][i].Hse)*zg.est[s][i].Hse/16.0;
                hf += zf1.est[s][i].H/4.0;
                vf += double(zf1.est[s][i].Hse)*zf1.est[s][i].Hse/16.0;
            }
            double H, var;
            long long n;
            comboH(rec.sums, i, unitCoef(5, iP1), H, var, n);
            ownG.add(i, H - hg, var + vg);
            ownF.add(i, H - hf, var + vf);
            shell.add(i, hf - hg, vf + vg);
        }
        out("  termination bias (P1; Zombie at each shell from 4 independent seeds, N = %d each): mean H_Z(stop %.1e) -\n"
            "    H_Z(stop %.1e) = %+.1e +- %.1e (pooled z %+.2f); own - Zombie: pooled z %+.2f at the gate's shell, %+.2f at the\n"
            "    first gate's\n", gp.nZombie, zombieStopDistance(true), zombieStopDistance(), shell.meanDiff(), shell.meanDiffSe(),
            shell.pooled(), ownG.pooled(), ownF.pooled());
        termMd = "Termination bias of Zombie (P1, 4 seeds per shell): mean H_Z(stop " + fmtd("%.1e", zombieStopDistance(true)) +
              ") - H_Z(stop " + fmtd("%.1e", zombieStopDistance()) + ") = " + fmtd("%+.1e", shell.meanDiff()) + " +- " +
              fmtd("%.1e", shell.meanDiffSe()) + "; own - Zombie pooled z " + fmtd("%+.2f", ownG.pooled()) + " (gate's shell), " +
              fmtd("%+.2f", ownF.pooled()) + " (first gate's).\n\n";
    }
    out("Zombie vs own loop: %s\n", zPass ? "agree within the combined standard errors" : "DISAGREE");
    out("wall clock: own recipe run %.2f s, own direct runs %.2f s each (mean)\n", rec.stats.seconds,
        (directs[0].stats.seconds + directs[1].stats.seconds + directs[2].stats.seconds)/3.0);
    // the verdict: the spec's criterion (which only confirms linearity, see above) and the independent cross-check
    // of the face path (the lead's note: Zombie must agree)
    const bool verdict = pass && zPass;
    out("CHECK 1: %s (spec superposition criterion: %s; Zombie cross-check of the face path: %s)\n", verdict ? "PASS" : "FAIL",
        pass ? "pass" : "FAIL", zPass ? "agrees" : "DISAGREES");
    md += "\nZombie cross-check (Zombie works in float, N = " + std::to_string(gp.nZombie) + "; own loop N = " +
          std::to_string(gp.nWalks) + "). Agree: abs(z) < 3 at >= 99% and abs(pooled z) < 3 over the " +
          std::to_string(ps.size() - size_t(nBand)) + " points farther than " + fmtd("%.2f", bandZ) + " from F, and abs(pooled z) "
          "< 3 over the " + std::to_string(nJunction) + " of them within " + fmtd("%.2f", kJunctionBand) + " of the segment "
          "where the crack meets F. pooled z = sum(dH)/sqrt(sum var). Mean z is reported but not used: between unequal N it "
          "is offset by the skewness of the walk values (each z is divided by its own sample s.d., which is correlated "
          "with its mean); the skew control column compares the own loop with itself at N = " + std::to_string(gp.nWalks) +
          " and " + std::to_string(gp.nZombie) + " and shows the same kind of offset without any solver difference.\n\n" +
          tz.str() + "\n" + termMd;
    md += "Recipe run wall clock: " + fmtd("%.2f", rec.stats.seconds) + " s (" + rec.stats.summary() + ").\n\n";
    md += std::string("**Check 1: ") + (verdict ? "PASS" : "FAIL") + "** (spec superposition criterion: " +
          (pass ? "pass" : "FAIL") + "; Zombie cross-check of the face path: " + (zPass ? "agrees" : "DISAGREES") + ")\n";

    // findings
    std::string f = std::string("- Check 1 (superposition: tet A, crack x = 0.3 with y_f = 0.30, F Dirichlet, the other faces "
                                "Neumann): ") + (verdict ? "PASS" : "FAIL") + " (the spec's criterion: " + (pass ? "pass" : "FAIL") +
                    "; Zombie cross-check: " + (zPass ? "agrees" : "DISAGREES") + ").\n";
    f += "- The spec's criterion, 3 random c x 2 g modes: abs(z) < 3 at " + pct(fracLo) + " or more of " + std::to_string(ps.size()) +
         " points (needed: 99%), abs(mean z) <= " + fmtd("%.4f", meanZHi) + " (bound 3/sqrt(M) = " +
         fmtd("%.4f", 3.0/std::sqrt(double(ps.size()))) + "), mean abs(z) " + fmtd("%.3f", mAbsLo) + " to " + fmtd("%.3f", mAbsHi) +
         " (0.80 for independent runs). H_direct - H_combo has an RMS of " + fmtd("%.1e", eLo) + " to " + fmtd("%.1e", eHi) +
         " of 2 g_max, against a noise RMS of " + fmtd("%.1e", nLo) + " to " + fmtd("%.1e", nHi) + ".\n";
    f += "- What that shows: the recipe and the direct runs share one walk process, so the criterion confirms that the scoring is "
         "linear and the seed streams independent; it cannot detect a wrong face path. The face path is checked by sanity (a) "
         "and (c), check 0 and the Zombie cross-check.\n";
    f += "- Zombie (float, independent) " + std::string(zPass ? "agrees" : "DISAGREES") + " on " + std::to_string(nZ) +
         " fields (H0 x 2, P_k x 3, 6 direct solves): pooled z " + fmtd("%+.2f", zpLo) + " to " + fmtd("%+.2f", zpHi) +
         " over the points farther than " + fmtd("%.2f", bandZ) + " from F (abs(z) < 3 at " + pct(zfLo) + " or more), " +
         fmtd("%+.2f", zjLo) + " to " + fmtd("%+.2f", zjHi) + " near the segment where the crack meets F.\n";
    f += "- S of A's recipe is positive definite: eigenvalues " + fmtd("%.5f", evS[0]) + ", " + fmtd("%.5f", evS[1]) + ", " +
         fmtd("%.5f", evS[2]) + " (s.e. " + fmtd("%.5f", evSSe[0]) + ", " + fmtd("%.5f", evSSe[1]) + ", " + fmtd("%.5f", evSSe[2]) +
         ").\n";
    f += "- Wall clock per own-loop run (" + std::to_string(ps.size()) + " points, N = " + std::to_string(gp.nWalks) + "): recipe " +
         fmtd("%.2f", rec.stats.seconds) + " s (5 fields), direct " +
         fmtd("%.2f", (directs[0].stats.seconds + directs[1].stats.seconds + directs[2].stats.seconds)/3.0) + " s (2 fields).\n";
    md.insert(mdIntro, findingsBlock(f));
    return verdict;
}

// =====================================================================================
// The reference on U (checks 2 and 3)
// =====================================================================================
// The physical ground truth on the union U = A u B. U is itself a tetrahedron (unionTet), so the own loop runs on
// it directly, in X: all four faces of U Neumann (they are the exterior faces of A u B), the crack built on U by
// buildCutOnTet (the plane clipped at the front), and g with L = sqrt(2), the L of the tets' xi frames (U's own
// longest edge is 2). This is the only solve of the gate that is not in a tet's xi frame: the per-tet solves stay
// on the reference tet in xi, and buildCut keeps its guard. Walk parameters: the per-tet runs', with the same
// absolute termination shell 1e-6 sqrt(2) and the mirror ball across U's faces. Both g types come from the same
// walks.
WalkParams referenceWalkParams()
{
    WalkParams wp = gateWalkParams(allNeumann(), true, false);
    wp.absorbShellRel = kGateShellRel*kSqrt2/unionTet().L;
    return wp;
}

// U's crack as the own loop's strip list: the frame is buildCutOnTet on U at t = 1 with L = sqrt(2), and one strip
// up to the t of the front y = yFront (unionCut's t)
StripCut unionStrips(const CrackSpec& cs)
{
    Cut frame = buildCutOnTet(unionTet(), cs.n, cs.d, cs.f, 1.0);
    frame.L = kSqrt2;
    StripCut sc = StripCut::physicalUnionReference(frame);
    sc.addStrip(sc.sOf((cs.yFront - frame.sMin)/(frame.sMax - frame.sMin)));
    return sc;
}

// The reference cache (gate_cache/, ignored by git): one binary file per configuration holding the sums. Its key
// string names everything the sums depend on: the configuration, seed, N, batches, the points and their hash, the
// live values of the walk parameters, of U and of the crack frame (refKeyParams), and kRefSolverVersion; a file
// with another key is recomputed and overwritten. kRefSolverVersion is a hand-kept tag for the walk code itself
// (walk_cache.cpp's runWalk and solveMulti, the scoring of the K data): bump it when that code changes.
const char* const kRefSolverVersion =
    "own loop (walk_cache.cpp runWalk + solveMulti) v1: double; K = 2 (+-g with g types a, b)";

std::string fmtg(double v) { return fmtd("%.17g", v); }
std::string fmtg3(const Vec3& v) { return "(" + fmtg(v.x()) + ", " + fmtg(v.y()) + ", " + fmtg(v.z()) + ")"; }

// The parameters of the reference's walks as they are used, printed exactly (%.17g): U's vertices and L, the walk
// parameters (referenceWalkParams), the gradient's first-ball factor, and the crack frame of unionStrips (plane, front
// direction, extent, L of g, the strips' fronts).
std::string refKeyParams(const StripCut& sc, const WalkParams& wp)
{
    const Tet U = unionTet();
    const Cut& f = sc.frame();
    std::string k = "U";
    for (const Vec3& v : U.v) k += " " + fmtg3(v);
    k += ", L_U " + fmtg(U.L) + " | walk: absorbShellRel " + fmtg(wp.absorbShellRel) + " (absolute " +
         fmtg(wp.absorbShellRel*U.L) + "), reflectShellRel " + fmtg(wp.reflectShellRel) + ", maxWalkLength " +
         std::to_string(wp.maxWalkLength) + ", faceBC";
    for (FaceBC b : wp.faceBC) k += b == FaceBC::Dirichlet ? " D" : " N";
    k += std::string(", cutData ") + (wp.cutData ? "1" : "0") + ", dirichletFaces " + (wp.dirichletFaces ? "1" : "0") +
         ", mirrorGradBall " + (wp.mirrorGradBall ? "1" : "0") + ", oddMirrorDirichlet " + (wp.oddMirrorDirichlet ? "1" : "0") +
         ", faceData " + (wp.faceData ? "set" : "none") + ", kGradFirstBallShrink " + fmtg(kGradFirstBallShrink) +
         " | crack frame: n " + fmtg3(f.n) + ", d " + fmtg(f.d) + ", u " + fmtg3(f.u) + ", sMin " + fmtg(f.sMin) + ", sMax " +
         fmtg(f.sMax) + ", L of g " + fmtg(f.L) + ", strips " + std::to_string(sc.nStrips()) + ", front s " + fmtg(sc.front()) +
         ", g types " + gTypeName(GType::One) + " / " + gTypeName(GType::SqrtFront) + " | ";
    return k;
}

uint64_t fnv1a(const void* data, size_t n)
{
    uint64_t h = 1469598103934665603ull;
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < n; i++) {
        h ^= p[i];
        h *= 1099511628211ull;
    }
    return h;
}

struct RefRun {
    std::vector<Vec3> X;    // the points: A's, then B's, then F's (physical)
    size_t nA = 0, nB = 0, nF = 0;
    MultiSums sums;         // K = 2: +-g with g type (a), (b)
    SolveStats stats;       // of the run that made the sums (from the cache: that original run)
    bool cached = false, ok = true;
    int nPilot = 0;
    double pilotSeconds = 0, pilotExtrapolated = 0;   // the timing pilot (not run when the cache is used)
    std::string key, file;
};

bool writeRefCache(const RefRun& r)
{
    const std::string tmp = r.file + ".tmp";
    {
        std::ofstream f(tmp, std::ios::binary);
        auto put = [&](const void* p, size_t n) { f.write(static_cast<const char*>(p), std::streamsize(n)); };
        const uint64_t kl = r.key.size(), np = r.sums.nPoints, nb = r.sums.batch.size(), nu = r.sums.uu.size();
        const int32_t hdr[3] = {r.sums.K, r.sums.nBatches, r.sums.walksPerBatch};
        const int64_t st[6] = {r.stats.nPoints, r.stats.walksRequested, r.stats.walksCounted, r.stats.walksEscaped,
                               r.stats.walksCapped, r.stats.totalSteps};
        put("PCUTREF1", 8);
        put(&kl, 8);
        put(r.key.data(), r.key.size());
        put(hdr, sizeof hdr);
        put(&np, 8);
        put(st, sizeof st);
        put(&r.stats.seconds, 8);
        put(&nb, 8);
        put(r.sums.batch.data(), nb*sizeof(double));
        put(&nu, 8);
        put(r.sums.uu.data(), nu*sizeof(double));
        put("PCUTEND1", 8);
        if (!f) return false;
    }
    std::error_code ec;
    std::filesystem::remove(r.file, ec);   // a stale file (another key) is replaced
    std::filesystem::rename(tmp, r.file, ec);
    return !ec;
}

// loads the sums and the run's statistics into r if r.file exists and holds r.key
bool readRefCache(RefRun& r)
{
    std::ifstream f(r.file, std::ios::binary);
    if (!f) return false;
    auto get = [&](void* p, size_t n) { return bool(f.read(static_cast<char*>(p), std::streamsize(n))); };
    char magic[8];
    uint64_t kl = 0, np = 0, nb = 0, nu = 0;
    if (!get(magic, 8) || std::memcmp(magic, "PCUTREF1", 8) != 0 || !get(&kl, 8) || kl != r.key.size()) return false;
    std::string key(kl, '\0');
    if (!get(key.data(), kl) || key != r.key) return false;
    int32_t hdr[3];
    int64_t st[6];
    double sec = 0;
    if (!get(hdr, sizeof hdr) || !get(&np, 8) || !get(st, sizeof st) || !get(&sec, 8) || !get(&nb, 8)) return false;
    MultiSums ms;
    ms.K = hdr[0];
    ms.nBatches = hdr[1];
    ms.walksPerBatch = hdr[2];
    ms.nPoints = np;
    if (ms.K != 2 || np != r.X.size() || nb != np*size_t(ms.nBatches)*ms.stride()) return false;
    ms.batch.resize(nb);
    if (!get(ms.batch.data(), nb*sizeof(double)) || !get(&nu, 8) || nu != np*size_t(ms.K)*size_t(ms.K)) return false;
    ms.uu.resize(nu);
    if (!get(ms.uu.data(), nu*sizeof(double)) || !get(magic, 8) || std::memcmp(magic, "PCUTEND1", 8) != 0) return false;
    r.sums = std::move(ms);
    r.stats = SolveStats{};
    r.stats.nPoints = int(st[0]);
    r.stats.walksRequested = st[1];
    r.stats.walksCounted = st[2];
    r.stats.walksEscaped = st[3];
    r.stats.walksCapped = st[4];
    r.stats.totalSteps = st[5];
    r.stats.seconds = sec;
    return true;
}

// The reference of one configuration at the points of A, B and F: from the cache if its key matches; otherwise a
// timing pilot (N/100 walks per point, its own seed), the extrapolated time, and the full run, written to the cache
// (not in quick runs, which neither read nor write it). A run that the pilot puts above 3 hours is not started
// (ok = false): N is never lowered.
RefRun runReference(Report& out, const GateParams& gp, const CrackSpec& cs, const std::string& tag, uint64_t seedTag,
                    const PointSet& pa, const PointSet& pb, const FacePointSet& fps, const std::string& cacheDir)
{
    RefRun r;
    r.X = pa.X;
    r.X.insert(r.X.end(), pb.X.begin(), pb.X.end());
    r.X.insert(r.X.end(), fps.X.begin(), fps.X.end());
    r.nA = pa.size();
    r.nB = pb.size();
    r.nF = fps.X.size();
    const Tet U = unionTet();
    const StripCut sc = unionStrips(cs);
    const WalkParams wp = referenceWalkParams();
    const Cut front = sc.frontCut(sc.nStrips());
    const TerminalValues tv = [&front](const WalkTerminal& t, double* u) {   // no Dirichlet face on U: on the cut
        u[0] = double(t.side)*front.g(t.x, GType::One);
        u[1] = double(t.side)*front.g(t.x, GType::SqrtFront);
    };
    const uint64_t seedBase = mixSeed(gp.seed, seedTag);
    char buf[640];
    std::snprintf(buf, sizeof buf,
                  "coupling gate reference on U | crack %s | y_f %.17g | seed %llu | N %d | batches %d | %zu points (A %zu, "
                  "B %zu, F %zu), hash %016llx | ", cs.name.c_str(), cs.yFront, (unsigned long long)seedBase, gp.nRef,
                  gp.nRefBatches, r.X.size(), r.nA, r.nB, r.nF, (unsigned long long)fnv1a(r.X.data(), r.X.size()*sizeof(Vec3)));
    r.key = std::string(buf) + refKeyParams(sc, wp) + kRefSolverVersion;
    r.file = (std::filesystem::path(cacheDir)/("ref_" + tag + ".bin")).string();
    out("reference on U (own loop, double; N = %d per point in %d batches, both g types from the same walks; H and grad H\n"
        "at all %zu points: A %zu, B %zu, F %zu)\n", gp.nRef, gp.nRefBatches, r.X.size(), r.nA, r.nB, r.nF);
    if (!gp.quick && readRefCache(r)) {
        r.cached = true;
        out("  loaded from the cache %s (key matches: configuration, seed, N, batches, points hash, the live walk\n"
            "  parameters, U and the crack frame, solver version)\n"
            "  original run: %s\n", r.file.c_str(), r.stats.summary().c_str());
        return r;
    }
    QueryPoints q;
    const std::vector<char> grad(r.X.size(), 1);
    q.add(r.X, &grad, seedBase);
    r.nPilot = gp.nRefBatches*std::max(1, gp.nRef/(100*gp.nRefBatches));
    {
        QueryPoints qp;
        qp.add(r.X, &grad, mixSeed(gp.seed, seedTag + 100));
        SolveStats ps;
        solveMulti(U, sc, qp, r.nPilot, gp.nRefBatches, wp, 2, tv, &ps);
        r.pilotSeconds = ps.seconds;
        r.pilotExtrapolated = ps.seconds*double(gp.nRef)/double(r.nPilot);
        out("  timing pilot, N = %d: %s\n  -> extrapolated to N = %d: %.0f s\n", r.nPilot, ps.summary().c_str(), gp.nRef,
            r.pilotExtrapolated);
    }
    if (r.pilotExtrapolated > 3*3600.0) {
        r.ok = false;
        out("  NOT RUN: the extrapolated time exceeds 3 hours (N is not lowered)\n");
        return r;
    }
    r.sums = solveMulti(U, sc, q, gp.nRef, gp.nRefBatches, wp, 2, tv, &r.stats);
    out("  run: %s\n", r.stats.summary().c_str());
    if (!gp.quick) {
        std::error_code ec;
        std::filesystem::create_directories(cacheDir, ec);
        out("  %s %s\n", writeRefCache(r) ? "cached to" : "COULD NOT WRITE the cache", r.file.c_str());
    }
    return r;
}

// A field at a set of points: its value and physical gradient, and the variances of their estimates.
struct FieldAt {
    std::vector<double> H, vH;
    std::vector<Vec3> G, vG;   // vG: per component
    void resize(size_t n)
    {
        H.assign(n, 0.0);
        vH.assign(n, 0.0);
        G.assign(n, Vec3::Zero());
        vG.assign(n, Vec3::Zero());
    }
};

// The reference (g type m) at its points [off, off + n): H and its variance from comboH, grad H = the mean of the
// batch gradients (U is solved in X: physical) and its variance from their spread.
FieldAt refField(const RefRun& r, int m, size_t off, size_t n)
{
    FieldAt f;
    f.resize(n);
    const std::vector<double> coef = unitCoef(2, m);
    const int B = r.sums.nBatches;
    tbb::parallel_for(size_t(0), n, [&](size_t i) {
        long long nn;
        comboH(r.sums, off + i, coef, f.H[i], f.vH[i], nn);
        std::vector<Vec3> g(static_cast<size_t>(B));
        Vec3 mean = Vec3::Zero();
        for (int b = 0; b < B; b++) mean += (g[size_t(b)] = comboBatchGrad(r.sums, off + i, b, coef))/double(B);
        Vec3 v = Vec3::Zero();
        for (const Vec3& x : g) v += (x - mean).cwiseAbs2();
        f.G[i] = mean;
        f.vG[i] = v/(double(B)*double(B - 1));
    });
    return f;
}

// The combination sum_k coef_k f_k of an own-loop run at its points (physical gradient J^-T grad_xi; grad = the mean
// of the batch gradients), with delete-one-batch jackknife variances. coefJack[a]: the coefficients of the
// replicate without batch a (the glued c of that replicate's Gram), so that the noise of the coefficients and its
// correlation with the fields' own noise are included; empty: fixed coefficients.
// A field may sum such terms from several runs of one tet on the same points (check 3's Level B: the H0 run with
// the face function J plus the recipe run's P_k). The runs' batches are paired by index: replicate a deletes batch
// a of every run (independent runs, so batch a of all of them is one unit of the jackknife).
struct RunTerm {
    const OwnRun* run = nullptr;
    std::vector<double> coef;
    std::vector<std::vector<double>> coefJack;
};
FieldAt runField(const TetFrame& T, const std::vector<RunTerm>& terms)
{
    const size_t P = terms.front().run->sums.nPoints;
    const int B = terms.front().run->sums.nBatches;
    std::vector<std::vector<std::vector<double>>> unit(terms.size());
    for (size_t ti = 0; ti < terms.size(); ti++) {
        const MultiSums& s = terms[ti].run->sums;
        if (s.nPoints != P || s.nBatches != B) throw std::logic_error("runField: runs differ in points or batches");
        for (int k = 0; k < s.K; k++) unit[ti].push_back(unitCoef(s.K, k));
    }
    FieldAt f;
    f.resize(P);
    tbb::parallel_for(size_t(0), P, [&](size_t i) {
        std::vector<double> Ha(static_cast<size_t>(B), 0.0);
        std::vector<Vec3> Ga(static_cast<size_t>(B), Vec3::Zero());
        for (size_t ti = 0; ti < terms.size(); ti++) {
            const RunTerm& t = terms[ti];
            const MultiSums& s = t.run->sums;
            const int K = s.K;
            const std::vector<double>& coef = t.coef;
            std::vector<double> nb(static_cast<size_t>(B)), su(size_t(B)*K), SU(static_cast<size_t>(K), 0.0);
            std::vector<Vec3> gb(size_t(B)*K), GS(static_cast<size_t>(K), Vec3::Zero());
            double n = 0;
            for (int b = 0; b < B; b++) {
                const double* d = s.batchData(i, b);
                nb[b] = d[0];
                n += d[0];
                for (int k = 0; k < K; k++) {
                    su[size_t(b)*K + k] = d[1 + k];
                    SU[k] += d[1 + k];
                    gb[size_t(b)*K + k] = T.JinvT*comboBatchGrad(s, i, b, unit[ti][k]);
                    GS[k] += gb[size_t(b)*K + k];
                }
            }
            for (int k = 0; k < K; k++) {
                f.H[i] += coef[k]*SU[k]/n;
                f.G[i] += coef[k]*GS[k]/double(B);
            }
            for (int a = 0; a < B; a++) {
                const std::vector<double>& ca = t.coefJack.empty() ? coef : t.coefJack[a];
                for (int k = 0; k < K; k++) {
                    if (ca[k] == 0) continue;
                    Ha[a] += ca[k]*(SU[k] - su[size_t(a)*K + k])/(n - nb[a]);
                    Ga[a] += ca[k]*(GS[k] - gb[size_t(a)*K + k])/double(B - 1);
                }
            }
        }
        double hm = 0;
        Vec3 gm = Vec3::Zero();
        for (int a = 0; a < B; a++) {
            hm += Ha[a]/double(B);
            gm += Ga[a]/double(B);
        }
        for (int a = 0; a < B; a++) {
            f.vH[i] += (Ha[a] - hm)*(Ha[a] - hm);
            f.vG[i] += (Ga[a] - gm).cwiseAbs2();
        }
        f.vH[i] *= double(B - 1)/double(B);
        f.vG[i] *= double(B - 1)/double(B);
    });
    return f;
}

// The variance that coefficients known only up to the covariance C (estimated elsewhere, independent of this run)
// add to the field sum_k c_k f_{idx[k]} of run r: p^T C p for H, sum_kl C_kl dP_k dP_l per gradient component.
void addCoefNoise(FieldAt& f, const OwnRun& r, const TetFrame& T, const std::vector<int>& idx, const Eigen::Matrix3d& C)
{
    const MultiSums& s = r.sums;
    tbb::parallel_for(size_t(0), s.nPoints, [&](size_t i) {
        Eigen::Vector3d p;
        std::array<Vec3, 3> g;
        for (int k = 0; k < 3; k++) {
            const std::vector<double> u = unitCoef(s.K, idx[k]);
            double var;
            long long n;
            comboH(s, i, u, p[k], var, n);
            g[k] = Vec3::Zero();
            for (int b = 0; b < s.nBatches; b++) g[k] += T.JinvT*comboBatchGrad(s, i, b, u)/double(s.nBatches);
        }
        f.vH[i] += p.dot(C*p);
        for (int a = 0; a < 3; a++) {
            const Eigen::Vector3d ga(g[0][a], g[1][a], g[2][a]);
            f.vG[i][a] += ga.dot(C*ga);
        }
    });
}

// The linear interpolant on tet T of its local corner values e at points X, with the jackknife variance through
// the replicates eJack (no other noise: the field has no walks of its own).
FieldAt linearField(const TetFrame& T, const std::vector<Vec3>& X, const Eigen::Vector4d& e,
                    const std::vector<Eigen::Vector4d>& eJack)
{
    const std::array<Vec3, 4> gN = gradN(T);
    const int B = int(eJack.size());
    FieldAt f;
    f.resize(X.size());
    for (size_t i = 0; i < X.size(); i++) {
        Eigen::Vector4d lam;
        for (int k = 0; k < 4; k++) lam[k] = lambdaOf(T, k, X[i]);
        auto val = [&](const Eigen::Vector4d& v, Vec3& g) {
            g = v[0]*gN[0] + v[1]*gN[1] + v[2]*gN[2] + v[3]*gN[3];
            return lam.dot(v);
        };
        f.H[i] = val(e, f.G[i]);
        if (B < 2) continue;
        std::vector<double> Ha(static_cast<size_t>(B));
        std::vector<Vec3> Ga(static_cast<size_t>(B));
        double hm = 0;
        Vec3 gm = Vec3::Zero();
        for (int a = 0; a < B; a++) {
            Ha[a] = val(eJack[a], Ga[a]);
            hm += Ha[a]/double(B);
            gm += Ga[a]/double(B);
        }
        for (int a = 0; a < B; a++) {
            f.vH[i] += (Ha[a] - hm)*(Ha[a] - hm)*double(B - 1)/double(B);
            f.vG[i] += (Ga[a] - gm).cwiseAbs2()*double(B - 1)/double(B);
        }
    }
    return f;
}

// Errors of a field against the reference over a set of points, divided by 2 g_max = 2 in both g types (the ramp's
// local 2g vanishes at the front). RMS is quadrature-weighted, max over the points; |grad| errors are vector norms.
// With each, its noise floor: the reference's and the field's per-point s.e. (RMS over the points and max, / 2), and
// the noise-corrected RMS sqrt(max(0, mean e^2 - mean (var_ref + var_field)))/2, an estimate of the error without
// Monte Carlo noise. (Divided by kTwoGMax, defined at the top.)
struct ErrStats {
    long long n = 0;
    double w = 0, eH2 = 0, eG2 = 0, maxH = 0, maxG = 0;
    double vRefH = 0, vFldH = 0, vRefG = 0, vFldG = 0;   // quadrature-weighted sums of the variances
    double maxSeRefH = 0, maxSeFldH = 0, maxSeRefG = 0, maxSeFldG = 0;
    Vec3 worstH = Vec3::Zero(), worstG = Vec3::Zero();
    void add(double wi, const Vec3& X, double eH, const Vec3& eG, double vrH, double vfH, const Vec3& vrG, const Vec3& vfG)
    {
        n++;
        w += wi;
        eH2 += wi*eH*eH;
        eG2 += wi*eG.squaredNorm();
        vRefH += wi*vrH;
        vFldH += wi*vfH;
        vRefG += wi*vrG.sum();
        vFldG += wi*vfG.sum();
        if (std::fabs(eH) > maxH) {
            maxH = std::fabs(eH);
            worstH = X;
        }
        if (eG.norm() > maxG) {
            maxG = eG.norm();
            worstG = X;
        }
        maxSeRefH = std::max(maxSeRefH, std::sqrt(vrH));
        maxSeFldH = std::max(maxSeFldH, std::sqrt(vfH));
        maxSeRefG = std::max(maxSeRefG, std::sqrt(vrG.sum()));
        maxSeFldG = std::max(maxSeFldG, std::sqrt(vfG.sum()));
    }
    void add(double wi, const Vec3& X, const FieldAt& fld, size_t i, const FieldAt& ref, size_t j)
    {
        add(wi, X, fld.H[i] - ref.H[j], fld.G[i] - ref.G[j], ref.vH[j], fld.vH[i], ref.vG[j], fld.vG[i]);
    }
    static double rms(double s, double w) { return w > 0 ? std::sqrt(std::max(0.0, s)/w)/kTwoGMax : 0.0; }
    double rmsH() const { return rms(eH2, w); }
    double rmsG() const { return rms(eG2, w); }
    double corrH() const { return rms(eH2 - vRefH - vFldH, w); }
    double corrG() const { return rms(eG2 - vRefG - vFldG, w); }
    double seRefH() const { return rms(vRefH, w); }
    double seFldH() const { return rms(vFldH, w); }
    double seRefG() const { return rms(vRefG, w); }
    double seFldG() const { return rms(vFldG, w); }
};

// delete-one-batch jackknife s.e. of a statistic from its replicates
double jackSe(const std::vector<double>& rep)
{
    const double B = double(rep.size());
    if (rep.size() < 2) return 0.0;
    double m = 0, v = 0;
    for (double x : rep) m += x/B;
    for (double x : rep) v += (x - m)*(x - m);
    return std::sqrt(v*(B - 1.0)/B);
}

// Least-squares fit (area-weighted over the F points) of a trace on F to the three hats: c = M^-1 r with
// M_kl = sum_p w_p hat_k hat_l, r_k = sum_p w_p hat_k (H_p - J_p), and Cov(c) = M^-1 (sum_p w_p^2 var_p hat hat^T) M^-1
// (the points are independent). J: a fixed part of the trace with coefficient 1 (check 3's Level B face function;
// none: 0). Its residual: traceError.
using FaceFn = std::function<double(const Vec3& X)>;
struct TraceFit {
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
};
TraceFit fitTrace(const FacePointSet& fps, const FieldAt& refF, const FaceFn* J = nullptr)
{
    Eigen::Matrix3d M = Eigen::Matrix3d::Zero(), V = Eigen::Matrix3d::Zero();
    Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
    for (size_t p = 0; p < fps.X.size(); p++) {
        const Eigen::Vector3d h(hatF(0, fps.X[p]), hatF(1, fps.X[p]), hatF(2, fps.X[p]));
        M += fps.w[p]*h*h.transpose();
        V += fps.w[p]*fps.w[p]*refF.vH[p]*h*h.transpose();
        rhs += fps.w[p]*(refF.H[p] - (J ? (*J)(fps.X[p]) : 0.0))*h;
    }
    TraceFit tf;
    const Eigen::Matrix3d Mi = M.inverse();
    tf.c = Mi*rhs;
    tf.cov = Mi*V*Mi;
    return tf;
}
// e = H_ref - J - sum_k c_k hat_k over the F points, / 2 g_max: RMS (area-weighted), noise-corrected RMS (the
// reference's variance subtracted, and the variance of sum_k c_k hat_k: the jackknife over the replicates cJack, or
// h^T cov h), max. With near (check 3: the F points within kNearSeg of the crack's segment on F) also the same over
// the near and the far points, the near points' share of the integral of e^2 (noise-corrected: of the integral of
// e^2 - var) and their share of the area.
struct TraceErr {
    double rms = 0, nc = 0, max = 0;
    double rmsNear = 0, ncNear = 0, maxNear = 0, rmsFar = 0, ncFar = 0, maxFar = 0, shareNear = 0, shareNearNc = 0, areaNear = 0;
};
TraceErr traceError(const FacePointSet& fps, const FieldAt& refF, const Eigen::Vector3d& c,
                    const std::vector<Eigen::Vector4d>* cJack = nullptr, const Eigen::Matrix3d* cov = nullptr,
                    const FaceFn* J = nullptr, const std::vector<char>* near = nullptr)
{
    double s = 0, sv = 0, w = 0, mx = 0;
    double sN[2] = {0, 0}, svN[2] = {0, 0}, wN[2] = {0, 0}, mxN[2] = {0, 0};   // [far, near]
    for (size_t p = 0; p < fps.X.size(); p++) {
        const Eigen::Vector3d h(hatF(0, fps.X[p]), hatF(1, fps.X[p]), hatF(2, fps.X[p]));
        const double e = refF.H[p] - (J ? (*J)(fps.X[p]) : 0.0) - c.dot(h);
        double v = refF.vH[p];
        if (cJack && cJack->size() > 1) {
            std::vector<double> rep;
            for (const Eigen::Vector4d& u : *cJack) rep.push_back(u.head<3>().dot(h));
            v += jackSe(rep)*jackSe(rep);
        }
        if (cov) v += h.dot(*cov*h);
        s += fps.w[p]*e*e;
        sv += fps.w[p]*v;
        w += fps.w[p];
        mx = std::max(mx, std::fabs(e));
        if (near) {
            const int k = (*near)[p] ? 1 : 0;
            sN[k] += fps.w[p]*e*e;
            svN[k] += fps.w[p]*v;
            wN[k] += fps.w[p];
            mxN[k] = std::max(mxN[k], std::fabs(e));
        }
    }
    TraceErr te{std::sqrt(s/w)/kTwoGMax, std::sqrt(std::max(0.0, s - sv)/w)/kTwoGMax, mx/kTwoGMax};
    if (near) {
        auto rms = [](double a, double b) { return b > 0 ? std::sqrt(std::max(0.0, a)/b)/kTwoGMax : 0.0; };
        te.rmsNear = rms(sN[1], wN[1]);
        te.ncNear = rms(sN[1] - svN[1], wN[1]);
        te.maxNear = mxN[1]/kTwoGMax;
        te.rmsFar = rms(sN[0], wN[0]);
        te.ncFar = rms(sN[0] - svN[0], wN[0]);
        te.maxFar = mxN[0]/kTwoGMax;
        te.shareNear = s > 0 ? sN[1]/s : 0.0;
        const double nc1 = std::max(0.0, sN[1] - svN[1]), nc0 = std::max(0.0, sN[0] - svN[0]);
        te.shareNearNc = nc0 + nc1 > 0 ? nc1/(nc0 + nc1) : 0.0;
        te.areaNear = w > 0 ? wN[1]/w : 0.0;
    }
    return te;
}

// Zombie (float) on U at the reference's points, an independent cross-check of the reference:
//  - H at every point, one solve of N = nZombie per g type (H only, so Zombie's s.e. of H is honest);
//  - grad H at the points farther than kRefZGradBand from U's faces, nZombieBatches solves of nZombie/nZombieBatches
//    walks per g type (s.e. from the spread over the solves). Next to a face Zombie's first ball is the plain one,
//    whose gradient variance grows like 1/distance^2, hence the band.
// Agree: for H, |z| < 3 at >= 99% of the points and |pooled z| < 3 (ZStats::passPooled); for the gradient, |pooled z|
// < 3 per component (with 16 solves each per-point s.e. has only 15 degrees of freedom, so the fraction |z| < 3 is
// shown but not used).
// The verdict uses all points at once, which is blind to a discrepancy that is local or of zero mean (e.g. + in A, -
// in B). So H is also compared per region (pooled z and the dispersion of ZStats; diagnostics, not in the verdict):
// A, B, F, near F (the points of A and B with |z| < kNearF) and near the crack (the points within kRefZNearCrack of
// U's crack).
constexpr double kRefZGradBand = 0.05;
constexpr double kNearF = 0.1;           // "near F": the points of A and B with |z| < 0.1
constexpr double kRefZNearCrack = 0.1;
constexpr int kRefZRegions = 5;
const char* const kRefZRegionName[kRefZRegions] = {"A", "B", "F", "near F (A, B; abs(z) < 0.1)", "near the crack (within 0.1)"};
struct RefZombie {
    ZStats H[2], G[2][3];
    ZStats Hr[2][kRefZRegions];   // H per region
    long long nGradPts = 0, nNaN = 0;
    double seconds = 0, secondsH = 0, secondsG = 0;   // all solves; the 2 H solves; the 2 x nZombieBatches gradient solves
    SolveStats st;
    bool agreeH(int m) const { return H[m].passPooled(); }
    bool agreeG(int m) const
    {
        for (int a = 0; a < 3; a++)
            if (!(G[m][a].n > 0 && std::fabs(G[m][a].pooled()) < 3.0)) return false;
        return true;
    }
    bool agree() const { return agreeH(0) && agreeH(1) && agreeG(0) && agreeG(1); }
};
RefZombie zombieReference(const GateParams& gp, const CrackSpec& cs, const RefRun& r, uint64_t seedTag)
{
    RefZombie z;
    const Tet U = unionTet();
    const Cut cut = unionCut(cs);
    WalkParams wp;   // all four faces Neumann; the reference's absolute shell, minimum star radius 1/10 of it
    wp.absorbShellRel = kGateShellRel*kSqrt2/U.L;
    wp.reflectShellRel = 0.1*wp.absorbShellRel;
    std::vector<size_t> gi;
    std::vector<Vec3> gX;
    for (size_t i = 0; i < r.X.size(); i++)
        if (U.distToFaces(r.X[i]) > kRefZGradBand) {
            gi.push_back(i);
            gX.push_back(r.X[i]);
        }
    z.nGradPts = (long long)gi.size();
    const std::vector<char> gradFlags(gX.size(), 1);
    const GType gts[2] = {GType::One, GType::SqrtFront};
    const int S = gp.nZombieBatches, nPer = gp.nZombie/gp.nZombieBatches, Bref = r.sums.nBatches;
    std::vector<std::array<char, kRefZRegions>> inRegion(r.X.size());
    for (size_t i = 0; i < r.X.size(); i++) {
        const bool onF = i >= r.nA + r.nB;
        inRegion[i] = {char(i < r.nA), char(i >= r.nA && !onF), char(onF), char(!onF && std::fabs(r.X[i].z()) < kNearF),
                       char(cut.distToCut(r.X[i]) < kRefZNearCrack)};
    }
    for (int m = 0; m < 2; m++) {
        const CutSolver solver(U, cut, gts[m], wp);
        std::vector<Estimate> e;
        const SolveStats sH = solver.solve(r.X, nullptr, gp.nZombie, mixSeed(gp.seed, seedTag + m), e);
        z.secondsH += sH.seconds;
        z.st.add(sH);
        const std::vector<double> coef = unitCoef(2, m);
        for (size_t i = 0; i < r.X.size(); i++) {
            double H, var;
            long long n;
            comboH(r.sums, i, coef, H, var, n);
            const double diff = H - e[i].H, v = var + double(e[i].Hse)*e[i].Hse;
            z.H[m].add(i, diff, v);
            for (int k = 0; k < kRefZRegions; k++)
                if (inRegion[i][k]) z.Hr[m][k].add(i, diff, v);
        }
        std::vector<std::vector<Estimate>> eg(static_cast<size_t>(S));
        for (int s = 0; s < S; s++) {
            const SolveStats sG = solver.solve(gX, &gradFlags, nPer, mixSeed(gp.seed, seedTag + 10 + 32*m + s), eg[s]);
            z.secondsG += sG.seconds;
            z.st.add(sG);
        }
        for (size_t j = 0; j < gi.size(); j++) {
            Vec3 zm = Vec3::Zero(), zv = Vec3::Zero();
            bool finite = true;
            for (int s = 0; s < S; s++) {
                finite = finite && eg[s][j].grad.allFinite();
                zm += eg[s][j].grad.cast<double>()/double(S);
            }
            if (!finite) {
                z.nNaN++;
                continue;
            }
            for (int s = 0; s < S; s++) zv += (eg[s][j].grad.cast<double>() - zm).cwiseAbs2()/(double(S)*double(S - 1));
            Vec3 om = Vec3::Zero(), ov = Vec3::Zero();
            std::vector<Vec3> g(static_cast<size_t>(Bref));
            for (int b = 0; b < Bref; b++) om += (g[size_t(b)] = comboBatchGrad(r.sums, gi[j], b, coef))/double(Bref);
            for (const Vec3& x : g) ov += (x - om).cwiseAbs2()/(double(Bref)*double(Bref - 1));
            for (int a = 0; a < 3; a++) z.G[m][a].add(gi[j], om[a] - zm[a], ov[a] + zv[a]);
        }
    }
    z.seconds = z.st.seconds;
    return z;
}

// The cross-check's console lines and table rows (checks 2 and 3): the verdict statistics (tZ) and H per region with
// the dispersion (tZr; diagnostics). The RMS discrepancy is in the units of the error tables (/ 2 g_max).
// Check 4's verdict on H: as agreeH, and where every walk of both solvers gives the same value (K3: +-1 exactly, zero
// variance; ZStats leaves such points out of z) the values must agree exactly; with no other point that is the verdict
bool agreeHExact(const RefZombie& zr, int m)
{
    const ZStats& h = zr.H[m];
    const bool exactOk = h.zeroVar == 0 || h.maxAbsDiffZeroVar <= 1e-6;
    return exactOk && (h.n == 0 ? h.zeroVar > 0 : zr.agreeH(m));
}
// regionNames, dom, nModes, exactZeroVar: check 4's regions, its domain O, g = 1 only where the ramp is not run, and its
// verdict on H (agreeHExact)
void reportRefZombie(Report& out, const GateParams& gp, const RefZombie& zr, const std::string& cfg, MdTable& tZ, MdTable& tZr,
                     const char* const* regionNames = kRefZRegionName, const char* dom = "U", int nModes = 2,
                     bool exactZeroVar = false)
{
    const GType gts[2] = {GType::One, GType::SqrtFront};
    const char* gShort[2] = {"(a) g = 1", "(b) sqrt ramp"};
    out("Zombie cross-check of the reference (Zombie works in float; %s; per solve: H %.1f s, gradient %.2f s):\n",
        zr.st.summary().c_str(), zr.secondsH/double(nModes), zr.secondsG/double(nModes*gp.nZombieBatches));
    for (int m = 0; m < nModes; m++) {
        const ZStats& h = zr.H[m];
        const bool agree = (exactZeroVar ? agreeHExact(zr, m) : zr.agreeH(m)) && zr.agreeG(m);
        if (exactZeroVar && h.zeroVar > 0)
            out("  %s: %lld points where every walk of both solvers gives the same value (left out of z): max |dH| there %.1e\n",
                gTypeName(gts[m]), h.zeroVar, h.maxAbsDiffZeroVar);
        out("  %s: H at %lld points, N = %d: |z| < 3 at %.2f%%, pooled z = %+.2f (mean dH = %+.1e +- %.1e), mean z = %+.4f;\n"
            "    grad H at the %lld points farther than %.2f from %s's faces (%d solves of %d walks; %lld not finite):\n"
            "    pooled z x/y/z = %+.2f / %+.2f / %+.2f (|z| < 3 at %.1f%% / %.1f%% / %.1f%%) -> %s\n",
            gTypeName(gts[m]), h.n, gp.nZombie, 100*h.frac(), h.pooled(), h.meanDiff(), h.meanDiffSe(), h.mean(), zr.nGradPts,
            kRefZGradBand, dom, gp.nZombieBatches, gp.nZombie/gp.nZombieBatches, zr.nNaN, zr.G[m][0].pooled(), zr.G[m][1].pooled(),
            zr.G[m][2].pooled(), 100*zr.G[m][0].frac(), 100*zr.G[m][1].frac(), 100*zr.G[m][2].frac(),
            agree ? "agree" : "DISAGREE");
        tZ.rows.push_back({cfg, gShort[m],
                           fmtd("%.4f", h.frac()) + ", " + fmtd("%+.2f", h.pooled()) + " (" + fmtd("%+.1e", h.meanDiff()) + " +- " +
                               fmtd("%.1e", h.meanDiffSe()) + ")",
                           fmtd("%+.2f", zr.G[m][0].pooled()) + " / " + fmtd("%+.2f", zr.G[m][1].pooled()) + " / " +
                               fmtd("%+.2f", zr.G[m][2].pooled()) + " (" + fmtd("%.3f", zr.G[m][0].frac()) + " / " +
                               fmtd("%.3f", zr.G[m][1].frac()) + " / " + fmtd("%.3f", zr.G[m][2].frac()) + ")",
                           agree ? (exactZeroVar && h.zeroVar > 0 ? "agree (" + std::to_string(h.zeroVar) + " points exact: max abs(dH) " +
                                                                        fmtd("%.0e", h.maxAbsDiffZeroVar) + ")"
                                                                  : std::string("agree"))
                                 : "DISAGREE"});
        out("    H by region (diagnostic): points, pooled z, mean z^2 (s.e.; 1 if they agree), RMS discrepancy beyond the noise\n"
            "    / 2 g_max (2-s.e. bound: the largest the comparison does not exclude):\n");
        for (int k = -1; k < kRefZRegions; k++) {
            const ZStats& s = k < 0 ? h : zr.Hr[m][k];
            const char* name = k < 0 ? "all points" : regionNames[k];
            out("      %-28s %6lld   %+.2f   %.3f (%.3f)   %.1e (%.1e)\n", name, s.n, s.pooled(), s.meanZ2(), s.meanZ2Se(),
                s.ncRmsDiff(kTwoGMax), s.ncRmsDiffBound(kTwoGMax));
            tZr.rows.push_back({cfg, gShort[m], name, std::to_string(s.n), fmtd("%+.2f", s.pooled()),
                                fmtd("%.3f", s.meanZ2()) + " (" + fmtd("%.3f", s.meanZ2Se()) + ")",
                                fmtd("%.1e", s.ncRmsDiff(kTwoGMax)) + " (" + fmtd("%.1e", s.ncRmsDiffBound(kTwoGMax)) + ")"});
        }
    }
}
// the largest 2-s.e. bound on the RMS discrepancy beyond the noise over the regions of the reference check (/ 2 g_max)
double zBoundWorst(const RefZombie& zr, int m)
{
    double b = 0;
    for (int k = 0; k < kRefZRegions; k++) b = std::max(b, zr.Hr[m][k].ncRmsDiffBound(kTwoGMax));
    return b;
}
// The summary tables' "wall clock per WoSt run" entry for Zombie: per solve (the H solve at N = nZombie, one gradient
// solve of nZombie/nZombieBatches walks) and the total over its solves
std::string zombieWall(const GateParams& gp, const RefZombie& zr)
{
    const int nG = 2*gp.nZombieBatches;
    return "Zombie per solve: H (N = " + std::to_string(gp.nZombie) + ") " + fmtd("%.1f", zr.secondsH/2.0) + ", grad (" +
           std::to_string(gp.nZombie/gp.nZombieBatches) + " walks) " + fmtd("%.2f", zr.secondsG/double(nG)) + "; total of its " +
           std::to_string(2 + nG) + " solves " + fmtd("%.1f", zr.seconds);
}
// The cross-check's description in gate_report.md (checks 2 and 3)
std::string refZombieMd(const GateParams& gp, const MdTable& tZ, const MdTable& tZr)
{
    return "\n### Reference cross-check (Zombie, float)\n\nZombie on U at the reference's points, N = " + std::to_string(gp.nZombie) +
           " for H (one solve per g mode), and for grad H " + std::to_string(gp.nZombieBatches) + " solves of " +
           std::to_string(gp.nZombie/gp.nZombieBatches) + " walks at the points farther than " + fmtd("%.2f", kRefZGradBand) +
           " from U's faces. Agree: H abs(z) < 3 at >= 99% and abs(pooled z) < 3; grad abs(pooled z) < 3 per component.\n\n" +
           tZ.str() +
           "\nWhat the verdict can see: Zombie's per-point s.e. of H (N = " + std::to_string(gp.nZombie) + ") dominates the "
           "combined one, so abs(z) < 3 at >= 99% flags only per-point discrepancies of several combined s.e., and the pooled z "
           "over all points only a mean offset. A local or zero-mean discrepancy (e.g. + in A and - in B) passes both. The "
           "table below is a diagnostic, not part of the verdict: per region the pooled z, the dispersion mean z^2 (1 if the "
           "two agree and their variances are exact; s.e. from the spread of z^2 over the points; the noise of the estimated "
           "variances biases it upward, most where nearly all walks of a point give the same value), and the RMS "
           "discrepancy beyond the noise, sqrt(max(0, mean(dH^2 - var)))/2 over the points (unweighted; mean(dH^2 - var) "
           "estimates the mean squared discrepancy without bias), in the units of the error tables, with in brackets the "
           "bound sqrt(max(0, mean(dH^2 - var)) + 2 s.e.)/2: the largest RMS discrepancy between the reference and Zombie "
           "that the comparison does not exclude (about 2 sigma). Zombie's own error enters it too, so it bounds the "
           "reference's error only together with Zombie's.\n\n" + tZr.str();
}

// =====================================================================================
// Check 2: two tets, B uncut, the crack does not touch F
// =====================================================================================

// The geometry of one configuration of checks 2 and 3, printed before its first solve: the vertices of A, B and U,
// the crack polygons per tet (xi and X) and on U, the strip cuts the walks use (the front at y = y_f in X, L, and the
// same "+g" side and g in the tets' xi frames and on U at the same physical points), and the face flags of every run.
bool printRunGeometry(Report& out, const GateParams& gp, const CrackSpec& cs, const XiCrack& xa, const XiCrack& xb,
                      const StripCut& usc, const std::vector<std::string>& runs, uint64_t seedTag)
{
    const TetFrame A = frameA(), B = frameB(), Uf = makeFrame("U", {1, 2, 3, 4});
    bool ok = true;
    out("geometry of this configuration (physical coordinates X):\n");
    for (const TetFrame* T : {&A, &B, &Uf})
        out("  %s: vertices v%d v%d v%d v%d = %s %s %s %s\n", T->name, T->global[0], T->global[1], T->global[2], T->global[3],
            fmt3(T->X[0], 0).c_str(), fmt3(T->X[1], 0).c_str(), fmt3(T->X[2], 0).c_str(), fmt3(T->X[3], 0).c_str());
    out("  crack %s\n", cs.name.c_str());
    printCutPolygons(out, A, xa);
    printCutPolygons(out, B, xb);
    printUnionPolygon(out, cs);
    if (!xa.cut.empty() && !xb.cut.empty()) {
        // Check 3: the union of A's and B's crack polygons (in X) is U's. A's lies in z >= 0 and B's in z <= 0, every
        // vertex of both lies in U's polygon (so both convex polygons do), and their areas add up to U's: then A's u B's
        // is a closed subset of U's polygon with the same area, hence all of it. Also shown: the segment they share on F,
        // and random points of the plane, which lie in U's polygon exactly when they lie in A's or B's.
        const Cut cu = unionCut(cs);
        auto inPoly = [](const Cut& c, const Vec3& x) { return c.polygonContains(x) || c.distToPolygonEdges(x) < 1e-12; };
        int outside = 0, wrongSide = 0;
        std::vector<Vec3> onF[2];
        for (const auto& [T, x, side] : {std::tuple<const TetFrame*, const XiCrack*, double>{&A, &xa, 1.0}, {&B, &xb, -1.0}})
            for (const Vec3& p : x->cut.polygon) {
                const Vec3 X = T->toX(p);
                outside += !inPoly(cu, X);
                wrongSide += side*X.z() < -1e-12;
                if (std::fabs(X.z()) <= 1e-12) onF[side > 0 ? 0 : 1].push_back(X);
            }
        const double dArea = xa.cut.area() + xb.cut.area() - cu.area();
        bool segOk = onF[0].size() == 2 && onF[1].size() == 2;
        double segErr = 0;
        if (segOk)
            for (const Vec3& p : onF[0])
                segErr = std::max(segErr, std::min((p - onF[1][0]).norm(), (p - onF[1][1]).norm()));
        segOk = segOk && segErr < 1e-12;
        double uLo = kInf, uHi = -kInf, wLo = kInf, wHi = -kInf;
        for (const Vec3& p : cu.polygon) {
            uLo = std::min(uLo, p.dot(cu.u));
            uHi = std::max(uHi, p.dot(cu.u));
            wLo = std::min(wLo, p.dot(cu.w));
            wHi = std::max(wHi, p.dot(cu.w));
        }
        std::mt19937_64 rng(mixSeed(gp.seed, seedTag + 50));
        std::uniform_real_distribution<double> Ud(0.0, 1.0);
        int nTest = 0, nInU = 0, badU = 0;
        for (int k = 0; k < 20000; k++) {
            const Vec3 X = cu.d*cu.n + (uLo - 0.05 + (uHi - uLo + 0.1)*Ud(rng))*cu.u + (wLo - 0.05 + (wHi - wLo + 0.1)*Ud(rng))*cu.w;
            const Vec3 xiA = A.toXi(X), xiB = B.toXi(X);
            if (cu.distToPolygonEdges(X) < 1e-9 || xa.cut.distToPolygonEdges(xiA) < 1e-9 || xb.cut.distToPolygonEdges(xiB) < 1e-9)
                continue;   // on an edge: the side is a matter of round-off
            const bool inU = cu.polygonContains(X), inAB = xa.cut.polygonContains(xiA) || xb.cut.polygonContains(xiB);
            nTest++;
            nInU += inU;
            badU += inU != inAB;
        }
        const bool uOk = outside == 0 && wrongSide == 0 && std::fabs(dArea) < 1e-12 && segOk && badU == 0;
        ok = ok && uOk;
        out("  union of the tets' crack polygons = U's polygon: A's and B's vertices outside U's polygon: %d; on the wrong side\n"
            "    of F (A: z >= 0, B: z <= 0): %d; area A %.6f + B %.6f - U %.6f = %.1e; shared segment on F: %s - %s (A and B agree\n"
            "    to %.1e); %d random points of the plane (box around U's polygon): %d in U's polygon, %d where membership in U's\n"
            "    differs from membership in A's or B's -> %s\n", outside, wrongSide, xa.cut.area(), xb.cut.area(), cu.area(), dArea,
            onF[0].size() > 0 ? fmt3(onF[0][0], 3).c_str() : "-", onF[0].size() > 1 ? fmt3(onF[0][1], 3).c_str() : "-", segErr,
            nTest, nInU, badU, uOk ? "OK" : "FAILED");
        // The configuration is symmetric under z -> -z (plane and front direction have no z part), and X = J_B xi turns
        // that mirror into the swap xi_x <-> xi_y: B's pulled-back crack must be A's with xi_x and xi_y swapped, and then
        // B's whole xi problem (F = xi_z = 0 included) is A's up to that swap, which maps B's P_v1 onto A's P_v1.
        if (cs.n.z() == 0 && cs.f.z() == 0) {
            auto sw = [](const Vec3& v) { return Vec3(v.y(), v.x(), v.z()); };
            double e = std::max({(xb.cut.n - sw(xa.cut.n)).cwiseAbs().maxCoeff(), (xb.cut.u - sw(xa.cut.u)).cwiseAbs().maxCoeff(),
                                 std::fabs(xb.cut.d - xa.cut.d), std::fabs(xb.cut.s - xa.cut.s)});
            for (const auto& [P, Q] : {std::pair<const Cut*, const Cut*>{&xa.cut, &xb.cut}, {&xb.cut, &xa.cut}})
                for (const Vec3& p : P->polygon) {
                    double m = kInf;
                    for (const Vec3& q : Q->polygon) m = std::min(m, (sw(p) - q).norm());
                    e = std::max(e, m);
                }
            const bool swOk = e < 1e-12 && xa.cut.polygon.size() == xb.cut.polygon.size();
            ok = ok && swOk;
            out("  mirror symmetry z -> -z: B's xi crack (n, d, u, s and the polygon) is A's with xi_x and xi_y swapped, max\n"
                "    deviation %.1e -> %s (so S_B = S_A and b_B = b_A up to noise, tested after the solves)\n", e, swOk ? "OK" : "FAILED");
        }
    }
    // the strip cuts of the walks: front line at y = y_f in X, L = sqrt(2)
    const Cut fu = usc.frontCut(usc.nStrips());
    double frontErr = std::fabs(fu.frontLinePoint().y() - cs.yFront);
    std::string s = "U " + fmtd("%.1e", frontErr);
    bool lOk = fu.L == kSqrt2;
    for (const auto& [T, x] : {std::pair<const TetFrame*, const XiCrack*>{&A, &xa}, {&B, &xb}}) {
        if (x->cut.empty()) continue;
        const Cut fx = x->strips.frontCut(x->strips.nStrips());
        const double e = std::fabs(T->toX(fx.frontLinePoint()).y() - cs.yFront);
        frontErr = std::max(frontErr, e);
        s += std::string(", ") + T->name + " " + fmtd("%.1e", e);
        lOk = lOk && fx.L == kSqrt2;
    }
    // "+g" side and g of the walks' strip cuts at the same physical points (random in U; A's frame for z >= 0, B's
    // for z < 0)
    std::mt19937_64 rng(mixSeed(gp.seed, seedTag));
    std::uniform_real_distribution<double> Ud(0.0, 1.0);
    const Tet U = unionTet();
    int nCmp = 0, badSide = 0;
    double gErr = 0;
    for (int k = 0; k < 4000; k++) {
        double l[4], sl = 0;
        for (double& v : l) sl += (v = -std::log(Ud(rng)));
        const Vec3 X = (l[0]*U.v[0] + l[1]*U.v[1] + l[2]*U.v[2] + l[3]*U.v[3])/sl;
        const TetFrame& T = X.z() >= 0 ? A : B;
        const XiCrack& x = X.z() >= 0 ? xa : xb;
        if (x.cut.empty() || std::fabs(cs.n.dot(X) - cs.d) < 1e-9) continue;
        const Cut fx = x.strips.frontCut(x.strips.nStrips());
        const Vec3 xi = T.toXi(X);
        nCmp++;
        if ((fx.side(xi) > 0) != (cs.n.dot(X) - cs.d > 0) || (fu.side(X) > 0) != (cs.n.dot(X) - cs.d > 0)) badSide++;
        for (GType g : {GType::One, GType::SqrtFront}) gErr = std::max(gErr, std::fabs(fx.g(xi, g) - fu.g(X, g)));
    }
    const bool sOk = frontErr < 1e-12 && lOk && badSide == 0 && gErr < 1e-12;
    ok = ok && sOk;
    for (const XiCrack* x : {&xa, &xb})
        if (!x->cut.empty()) ok = ok && x->nFront >= 2 && x->frontErr < 1e-12;
    out("  strip cuts of the walks (the tets' in xi, U's in X): front line at y = y_f in X: max |y - y_f| = %s; L = sqrt(2):\n"
        "  %s; at %d random points of U (in the frame of the tet that holds them): %d disagree on the +g side, max |g_tet(xi) -\n"
        "  g_U(X)| = %.1e over both g types -> %s\n", s.c_str(), lOk ? "yes" : "NO", nCmp, badSide, gErr, sOk ? "OK" : "FAILED");
    out("  runs and per-face flags:\n");
    for (const std::string& r : runs) out("    %s\n", r.c_str());
    return ok;
}

// The glue of check 2, written out. Unknowns: c = (c_v0, c_v1, c_v2), A's response coefficients = the corner values
// on F (global ids), and d = c_v4, B's apex. Total energy
//   E(c, d) = (1/2) c^T S_A c + b_A^T c + const + (1/2) [c; d]^T K_B [c; d],
// K_B = B's linear-tet stiffness with its local vertices (v0, v2, v1, v4) mapped to the unknowns (v0, v1, v2, v4)
// (Kg). dE/dc = dE/dd = 0 gives the 4 x 4 system
//   [ S_A + K_cc   K_cd ] [c]   [ -b_A ]
//   [ K_dc         K_dd ] [d] = [   0  ]
// S_A is symmetrized first. S_A + K_cc is positive definite (the crack is Dirichlet in A, so S_A is) and K_dd = 1/6,
// so the system is regular. Solved with Eigen in double. Returns (c_v0, c_v1, c_v2, d); M and rhs are optional.
Eigen::Vector4d glueCheck2(const Eigen::Matrix3d& SA, const Eigen::Vector3d& bA, const Eigen::Matrix4d& Kg,
                           Eigen::Matrix4d* Mout = nullptr, Eigen::Vector4d* rhsOut = nullptr)
{
    const Eigen::Matrix3d S = 0.5*(SA + SA.transpose());
    Eigen::Matrix4d M;
    M.topLeftCorner<3, 3>() = S + Kg.topLeftCorner<3, 3>();   // S_A + K_cc
    M.topRightCorner<3, 1>() = Kg.topRightCorner<3, 1>();     // K_cd
    M.bottomLeftCorner<1, 3>() = Kg.bottomLeftCorner<1, 3>(); // K_dc
    M(3, 3) = Kg(3, 3);                                       // K_dd
    Eigen::Vector4d rhs;
    rhs << -bA, 0.0;
    if (Mout) *Mout = M;
    if (rhsOut) *rhsOut = rhs;
    return M.fullPivLu().solve(rhs);
}

// The glued (c, d) of a recipe run's Gram (fields H0 g1, H0 g2, P_v0, P_v1, P_v2) for g type m, its delete-one-batch
// replicates (the system re-solved with each replicate Gram) and their jackknife s.e.
struct Glue {
    Eigen::Vector4d u = Eigen::Vector4d::Zero(), se = Eigen::Vector4d::Zero();
    std::vector<Eigen::Vector4d> jack;
};
Glue glueRun(const Gram& gm, int m, const Eigen::Matrix4d& Kg, Eigen::Matrix4d* M = nullptr, Eigen::Vector4d* rhs = nullptr)
{
    Glue gl;
    gl.u = glueCheck2(gm.G.block(2, 2, 3, 3), Eigen::Vector3d(gm.G(m, 2), gm.G(m, 3), gm.G(m, 4)), Kg, M, rhs);
    for (const Eigen::MatrixXd& J : gm.jack)
        gl.jack.push_back(glueCheck2(J.block(2, 2, 3, 3), Eigen::Vector3d(J(m, 2), J(m, 3), J(m, 4)), Kg));
    for (int k = 0; k < 4; k++) {
        std::vector<double> rep;
        for (const Eigen::Vector4d& x : gl.jack) rep.push_back(x[k]);
        gl.se[k] = jackSe(rep);
    }
    return gl;
}

std::string fmtVec(const Eigen::VectorXd& v, const char* f = "%+.5f")
{
    std::string s = "(";
    for (int i = 0; i < v.size(); i++) s += (i ? ", " : "") + fmtd(f, v[i]);
    return s + ")";
}

// Check 2's unknowns are (c_v0, c_v1, c_v2, d = c_v4): B's stiffness on them (from B's local order v0 v2 v1 v4), and
// B's local corner values from them
int check2Unknown(int g) { return g == 4 ? 3 : g; }
Eigen::Matrix4d check2Kg(const TetFrame& B)
{
    const Eigen::Matrix4d KB = stiffness(B);
    Eigen::Matrix4d Kg = Eigen::Matrix4d::Zero();
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) Kg(check2Unknown(B.global[i]), check2Unknown(B.global[j])) = KB(i, j);
    return Kg;
}
Eigen::Vector4d check2LocalB(const TetFrame& B, const Eigen::Vector4d& u)
{
    Eigen::Vector4d e;
    for (int k = 0; k < 4; k++) e[k] = u[check2Unknown(B.global[k])];
    return e;
}

// The recipe's fields (checks 2 and 3): H0 for both g types (crack +-g, F 0) and P_v0, P_v1, P_v2 (crack 0, F hat_k)
std::vector<FieldSpec> recipeFields()
{
    std::vector<FieldSpec> fields = {FieldSpec{"H0 g1", true, GType::One, {}}, FieldSpec{"H0 g2", true, GType::SqrtFront, {}}};
    for (int g = 0; g < 3; g++)
        fields.push_back(FieldSpec{"P" + std::to_string(g), false, GType::One, [g](const Vec3& X) { return hatF(g, X); }});
    return fields;
}

bool runCheck2(Report& out, const GateParams& gp, const FacePointSet& fps, const std::string& cacheDir, std::string& md,
               std::string& geoMd)
{
    const TetFrame A = frameA(), B = frameB();
    const GType gts[2] = {GType::One, GType::SqrtFront};
    const char* gShort[2] = {"(a) g = 1", "(b) sqrt ramp"};
    out("\n=== Check 2: two tets, B uncut, the crack does not touch F ===\n");
    out("crack plane z = z_c (+g above), front {z = z_c, y = y_f}, crack region y < y_f, y_f = (1 - z_c)/2, z_c = 0.60 and 0.35.\n"
        "A: Level A recipe (F Dirichlet: H0 with crack +-g and F 0 for both g types, P_k with crack 0 and F hat_k, k = v0, v1,\n"
        "v2, all from the same walks; own loop, N = %d in %d batches). B: uncut, no WoSt: the linear interpolant of its corner\n"
        "values, energy (1/2) e^T K_B e. Glue: the 4 x 4 system for (c_v0, c_v1, c_v2, d = c_v4). Reference: one own-loop\n"
        "solve on U = A u B (N = %d), cached in %s. Errors are divided by 2 g_max = 2 in both g types; RMS quadrature-\n"
        "weighted, max over the points; energies are (1/2) int |grad H|^2 over A minus the eps-tube and B, split-sample.\n",
        gp.nWalks, gp.nBatches, gp.nRef, cacheDir.c_str());
    md += "## Check 2: two tets, B uncut, the crack does not touch F\n\n";
    md += "Crack plane z = z_c (+g above), front {z = z_c, y = y_f}, crack region y < y_f, y_f = (1 - z_c)/2. A: Level A recipe "
          "(F Dirichlet; H0 and P_v0, P_v1, P_v2 from the same walks; own loop, double, N = " + std::to_string(gp.nWalks) + " in " +
          std::to_string(gp.nBatches) + " batches). B: uncut, no WoSt, the linear interpolant of its corner values with energy "
          "(1/2) e^T K_B e. Glue: [S_A + K_cc, K_cd; K_dc, K_dd][c; d] = [-b_A; 0] (unknowns c_v0, c_v1, c_v2, d = c_v4). "
          "Reference: one own-loop solve on U = A u B in X (all four faces Neumann, the crack clipped to U, L = sqrt(2)), N = " +
          std::to_string(gp.nRef) + " in " + std::to_string(gp.nRefBatches) + " batches, cached in gate_cache/.\n\n"
          "All errors and noise levels are divided by **2 g_max = 2** in both g modes (the ramp's local 2g vanishes at the "
          "front). RMS is quadrature-weighted over the points, max is over the points; grad H errors are vector norms. "
          "\"nc\" = noise-corrected RMS, sqrt(max(0, mean e^2 - mean(var_ref + var_dec)))/2. The decomposed field's noise "
          "is the delete-one-batch jackknife over the recipe's batches, with the glued c and d recomputed in every "
          "replicate (so it includes the noise of c and d). Energies are E = (1/2) int abs(grad H)^2 over A minus the "
          "eps-tube and over B, from split-sample products (all pairs of different batches), jackknife s.e.; E_B of the "
          "decomposition is exact given c, d. Energy identity: for any H_dec that carries the crack data and is continuous "
          "across F, E_dec - E_true = (1/2) int abs(grad e)^2 + int over the eps-tube's surface of e dH_true/dn (e = H_dec - "
          "H_true; the cross term int grad H_true . grad e has no other boundary part, since H_true has zero flux on U's faces "
          "and e = 0 on the crack). Here H_dec is continuous across F for every c and d (A's face data and B's interpolant are "
          "both c . hat there), so the identity holds whatever the glue gives: (1/2) int abs(grad e)^2, shown noise-corrected "
          "next to E_dec - E_ref, checks the energy and error estimators against each other (up to the tube term), not the "
          "glue.\n\n";
    const size_t mdIntro = md.size();   // the findings go here
    struct Case {                       // the numbers behind the findings, per configuration and g mode
        std::string cfg, g;
        double rms, nc, max, noise, nearRms, nearNc, dE, seDE, fitA, gluA, fitF, gluF, meanF, dcMax, dcSig, gRms, gNc, gNoise;
        double zAll, zWorst, calLo, calHi;   // the reference check's bound (all points, worst region); replicate sd / jackknife
    };
    std::vector<Case> cases;
    int nZagree = 0, nZcases = 0;
    std::string geoTxt;   // the geometry printout of every configuration (appendix)
    MdTable tMain, tReg, tFit, tRep, tEn, tRun, tZ, tZr, tQ;
    tMain.head = {"configuration", "g mode", "H: RMS / max (A u B)", "grad H: RMS / max (A u B)",
                  "reference noise: s.e. H RMS / max; grad H RMS / max", "decomposed noise: s.e. H RMS / max; grad H RMS / max",
                  "energy E (A u B): decomposed / reference / difference", "S_A eigenvalues (s.e.)", "glued c_v0, c_v1, c_v2; d (s.e.)",
                  "wall clock per WoSt run (s)"};
    tReg.head = {"configuration", "g mode", "region", "points", "H RMS (nc)", "H max", "grad H RMS (nc)", "grad H max",
                 "ref s.e. H RMS / max", "dec s.e. H RMS / max", "ref s.e. grad RMS / max", "dec s.e. grad RMS / max"};
    tFit.head = {"configuration", "g mode", "c glued (s.e.)", "c fit (s.e.)", "fit residual on F: RMS (nc) / max",
                 "glued trace vs reference on F: RMS (nc) / max", "ref s.e. on F: RMS / max", "region", "glued: H RMS (nc) / max",
                 "fitted trace: H RMS (nc) / max", "glued: grad RMS (nc) / max", "fitted trace: grad RMS (nc) / max",
                 "fitted-trace noise: s.e. H RMS / max"};
    tRep.head = {"configuration", "g mode", "runs", "sd over the runs / rms jackknife s.e.: c_v0 / c_v1 / c_v2 / d",
                 "mean glued c over the runs +- s.e.", "c fit +- s.e.", "mean glued c - c fit +- s.e.",
                 "mean glued c's trace vs reference on F: RMS (nc) / max"};
    tEn.head = {"configuration", "g mode", "E_A dec", "E_A ref", "E_B dec (exact given c, d)", "E_B ref", "E dec", "E ref",
                "E dec - E ref", "(1/2) int abs(grad e)^2 (nc)", "E_A with the fitted trace", "naive - split: dec / ref"};
    tRun.head = {"configuration", "run", "solver", "N per point", "points", "wall clock (s)", "mean walk length",
                 "capped / escaped", "note"};
    tZ.head = {"configuration", "g mode", "H: frac abs(z)<3, pooled z (mean dH +- s.e.)",
               "grad H (points > 0.05 from U's faces): pooled z per component (frac abs(z)<3)", "verdict"};
    tZr.head = {"configuration", "g mode", "region", "points", "H pooled z", "mean z^2 (s.e.)",
                "RMS discrepancy beyond noise / 2 g_max (2-s.e. bound)"};
    tQ.head = {"point set", "kept points", "dropped (eps)", "max rel err, all points vs exact over T",
               "max rel err, kept vs T minus tube", "tube volume"};

    bool ok = true;
    const double zcs[2] = {0.60, 0.35};
    for (int ci = 0; ci < 2; ci++) {
        const CrackSpec cs = crackCheck2(zcs[ci]);
        const std::string cfg = "z_c = " + fmtd("%.2f", zcs[ci]);
        const size_t geo0 = out.text.size();
        out("\n--- check 2, %s (y_f = %.3f) ---\n", cfg.c_str(), cs.yFront);
        const XiCrack xa = xiCrack(A, cs), xb = xiCrack(B, cs);
        const StripCut usc = unionStrips(cs);
        const std::vector<std::string> runs = {
            "reference U (own loop; Zombie cross-check): " + faceFlags(makeFrame("U", {1, 2, 3, 4}), allNeumann(), "", false) +
                "; crack +-g (both g types)",
            "A recipe (H0 g1, H0 g2: crack +-g, F 0; P_v0, P_v1, P_v2: crack 0, F hat_k): " + faceFlags(A, onlyF(), "(data per field)"),
            "B: no solve (uncut): the linear interpolant of its corner values (local v0 v2 v1 v4), energy (1/2) e^T K_B e"};
        const bool geoOk = printRunGeometry(out, gp, cs, xa, xb, usc, runs, 7950 + ci);
        geoTxt += out.text.substr(geo0);
        if (!geoOk) {
            out("geometry assertions FAILED: check 2 at %s not run\n", cfg.c_str());
            ok = false;
            continue;
        }
        if (!xb.cut.empty()) throw std::logic_error("check 2: the crack must not reach B");

        // point sets
        const PointSet pa = volumePoints(A, xa, gp.nPoints, gp.eps, gp.quick ? 1 : 4);
        const PointSet pb = volumePoints(B, xb, gp.nPoints, gp.eps, gp.quick ? 1 : 4);
        for (const auto& [T, x, ps] : {std::tuple<const TetFrame*, const XiCrack*, const PointSet*>{&A, &xa, &pa}, {&B, &xb, &pb}}) {
            const std::string label = std::string(T->name) + ", " + cfg;
            out("point set %s: %d leaves (depth 2/3/4: %d/%d/%d), %zu generated, %zu kept\n", label.c_str(), ps->nLeaves,
                ps->leavesByDepth[2], ps->leavesByDepth[3], ps->leavesByDepth[4], ps->allXi.size(), ps->size());
            quadratureCheck(out, *T, *x, *ps, gp.nTube, &tQ, label);
        }

        // the reference and its Zombie cross-check
        const RefRun ref = runReference(out, gp, cs, "check2_zc" + fmtd("%.2f", zcs[ci]), 8000 + ci, pa, pb, fps, cacheDir);
        if (!ref.ok) {
            ok = false;
            continue;
        }
        const RefZombie zr = zombieReference(gp, cs, ref, 8200 + 100*ci);
        reportRefZombie(out, gp, zr, cfg, tZ, tZr);
        ok = ok && zr.agree();
        for (int m = 0; m < 2; m++) {
            nZcases++;
            nZagree += zr.agreeH(m) && zr.agreeG(m);
        }

        // A's recipe
        const std::vector<FieldSpec> fields = recipeFields();
        const OwnRun rec = runOwn(A, xa, onlyF(), fields, pa.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7500 + ci));
        out("A recipe run: %s\n", rec.stats.summary().c_str());
        std::vector<std::vector<double>> coefs;
        for (int k = 0; k < 5; k++) coefs.push_back(unitCoef(5, k));
        const Gram gm = gramOwn(rec, pa, A, coefs);
        const int Bq = gm.B;
        const Eigen::Matrix3d SA = gm.G.block(2, 2, 3, 3);
        std::vector<Eigen::MatrixXd> jS;
        for (const auto& J : gm.jack) jS.push_back(J.block(2, 2, 3, 3));
        Eigen::VectorXd ev, evSe;
        eigenWithSe(gm, SA, jS, ev, evSe);
        out("S_A (split-sample, jackknife s.e.):\n");
        for (int i = 0; i < 3; i++)
            out("  %s   s.e. %s   naive-split %s\n", matRow(SA, i).c_str(), matRow(gm.se.block(2, 2, 3, 3), i, "%.5f").c_str(),
                matRow(gm.naive.block(2, 2, 3, 3) - SA, i, "%+.1e").c_str());
        out("  eigenvalues %+.5f %+.5f %+.5f (s.e. %.5f %.5f %.5f): %s\n", ev[0], ev[1], ev[2], evSe[0], evSe[1], evSe[2],
            ev[0] > 3*evSe[0] ? "positive definite" : "NOT clearly positive definite");

        // B's stiffness on the unknowns (v0, v1, v2, v4)
        const Eigen::Matrix4d KB = stiffness(B), Kg = check2Kg(B);
        auto localB = [&B](const Eigen::Vector4d& u) { return check2LocalB(B, u); };

        // Replicate recipe runs with independent seeds (run 0 is the main run): the spread of the glued c, d over the
        // runs against one run's jackknife s.e. (the method behind the decomposed noise floor), and the mean over the
        // runs, sqrt(nRep) less noisy than one run, which shows a systematic offset of the glued c that one run's
        // noise can hide (compared with the fitted trace below)
        const int nRep = gp.quick ? 3 : 8;
        std::vector<std::array<Glue, 2>> reps(static_cast<size_t>(nRep));
        double repSec = 0;
        for (int rep = 0; rep < nRep; rep++) {
            Gram g;
            if (rep == 0) {
                g = gm;
            } else {
                const OwnRun rr = runOwn(A, xa, onlyF(), fields, pa.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7600 + 10*ci + rep));
                repSec += rr.stats.seconds;
                g = gramOwn(rr, pa, A, coefs);
            }
            for (int m = 0; m < 2; m++) reps[size_t(rep)][m] = glueRun(g, m, Kg);
        }
        out("replicate recipe runs: %d (seed tags 7500 + %d, then 7600 + %d + 1..%d; %.1f s for the extra runs)\n", nRep, ci,
            10*ci, nRep - 1, repSec);

        // the reference's energies per region (split-sample over its batches)
        std::vector<double> wA(ref.X.size(), 0.0), wB(ref.X.size(), 0.0), wAB(ref.X.size(), 0.0);
        for (size_t i = 0; i < ref.nA; i++) wA[i] = wAB[i] = pa.w[i];
        for (size_t i = 0; i < ref.nB; i++) wB[ref.nA + i] = wAB[ref.nA + i] = pb.w[i];
        const std::vector<double> unit2[2] = {unitCoef(2, 0), unitCoef(2, 1)};
        auto refGram = [&](const std::vector<double>& w) {
            return gramFromBatches(w, ref.sums.nBatches, 2, [&](size_t p, int b, Vec3* g) {
                for (int q = 0; q < 2; q++) g[q] = comboBatchGrad(ref.sums, p, b, unit2[q]);
            });
        };
        const Gram grA = refGram(wA), grB = refGram(wB), grAB = refGram(wAB);

        std::string wallRef = ref.cached ? fmtd("%.1f", ref.stats.seconds) + " (cached run)" : fmtd("%.1f", ref.stats.seconds);
        for (int m = 0; m < 2; m++) {
            out("\n%s, %s:\n", cfg.c_str(), gTypeName(gts[m]));
            const Eigen::Vector3d bA(gm.G(m, 2), gm.G(m, 3), gm.G(m, 4));
            Eigen::Matrix4d M;
            Eigen::Vector4d rhs;
            const Glue gl = glueRun(gm, m, Kg, &M, &rhs);
            const Eigen::Vector4d& u = gl.u;
            const Eigen::Vector4d& uSe = gl.se;
            const std::vector<Eigen::Vector4d>& uJ = gl.jack;
            std::vector<Eigen::Vector4d> eJ(static_cast<size_t>(Bq));
            std::vector<std::vector<double>> coefJ(static_cast<size_t>(Bq));
            for (int a = 0; a < Bq; a++) {
                eJ[a] = localB(uJ[a]);
                coefJ[a] = {m == 0 ? 1.0 : 0.0, m == 1 ? 1.0 : 0.0, uJ[a][0], uJ[a][1], uJ[a][2]};
            }
            out("  b_A = %s (s.e. %s)\n", fmtVec(bA).c_str(), fmtVec(Eigen::Vector3d(gm.se(m, 2), gm.se(m, 3), gm.se(m, 4)), "%.5f").c_str());
            out("  4 x 4 glue system [S_A + K_cc, K_cd; K_dc, K_dd] [c_v0, c_v1, c_v2, d] = [-b_A; 0]:\n");
            for (int i = 0; i < 4; i++) out("    [%s]  [%+.5f]\n", matRow(M, i).c_str(), rhs[i]);
            out("  glued: c = (%+.5f, %+.5f, %+.5f) +- (%.5f, %.5f, %.5f), d = c_v4 = %+.5f +- %.5f (jackknife over the\n"
                "  recipe's batches, the system re-solved per replicate)\n", u[0], u[1], u[2], uSe[0], uSe[1], uSe[2], u[3], uSe[3]);

            // fields at the points
            const std::vector<double> coef = {m == 0 ? 1.0 : 0.0, m == 1 ? 1.0 : 0.0, u[0], u[1], u[2]};
            const FieldAt decA = runField(A, {RunTerm{&rec, coef, coefJ}});
            const FieldAt decB = linearField(B, pb.X, localB(u), eJ);
            const FieldAt refA = refField(ref, m, 0, ref.nA), refB = refField(ref, m, ref.nA, ref.nB),
                          refF = refField(ref, m, ref.nA + ref.nB, ref.nF);
            ErrStats eA, eB, eAB, eNF, eNFA;
            for (size_t i = 0; i < ref.nA; i++) {
                eA.add(pa.w[i], pa.X[i], decA, i, refA, i);
                eAB.add(pa.w[i], pa.X[i], decA, i, refA, i);
                if (std::fabs(pa.X[i].z()) < kNearF) {
                    eNF.add(pa.w[i], pa.X[i], decA, i, refA, i);
                    eNFA.add(pa.w[i], pa.X[i], decA, i, refA, i);
                }
            }
            for (size_t i = 0; i < ref.nB; i++) {
                eB.add(pb.w[i], pb.X[i], decB, i, refB, i);
                eAB.add(pb.w[i], pb.X[i], decB, i, refB, i);
                if (std::fabs(pb.X[i].z()) < kNearF) eNF.add(pb.w[i], pb.X[i], decB, i, refB, i);
            }
            out("  errors of the glued decomposition vs the reference, / 2 g_max (nc = noise-corrected RMS):\n"
                "    region         points   H RMS (nc)          H max     grad RMS (nc)       grad max  | ref s.e. H RMS/max  dec s.e. H RMS/max  ref s.e. grad RMS/max  dec s.e. grad RMS/max\n");
            for (const auto& [name, e] : {std::pair<const char*, const ErrStats*>{"A", &eA}, {"B", &eB}, {"A u B", &eAB},
                                          {"near F (A, B)", &eNF}, {"near F (A)", &eNFA}}) {
                out("    %-14s %6lld   %.2e (%.2e)   %.2e  %.2e (%.2e)   %.2e  | %.1e / %.1e     %.1e / %.1e     %.1e / %.1e       %.1e / %.1e\n",
                    name, e->n, e->rmsH(), e->corrH(), e->maxH/kTwoGMax, e->rmsG(), e->corrG(), e->maxG/kTwoGMax, e->seRefH(),
                    e->maxSeRefH/kTwoGMax, e->seFldH(), e->maxSeFldH/kTwoGMax, e->seRefG(), e->maxSeRefG/kTwoGMax, e->seFldG(),
                    e->maxSeFldG/kTwoGMax);
                if (std::string(name) == "near F (A)") continue;
                tReg.rows.push_back({cfg, gShort[m], name, std::to_string(e->n), fmtd("%.2e", e->rmsH()) + " (" + fmtd("%.2e", e->corrH()) + ")",
                                     fmtd("%.2e", e->maxH/kTwoGMax), fmtd("%.2e", e->rmsG()) + " (" + fmtd("%.2e", e->corrG()) + ")",
                                     fmtd("%.2e", e->maxG/kTwoGMax), fmtd("%.1e", e->seRefH()) + " / " + fmtd("%.1e", e->maxSeRefH/kTwoGMax),
                                     fmtd("%.1e", e->seFldH()) + " / " + fmtd("%.1e", e->maxSeFldH/kTwoGMax),
                                     fmtd("%.1e", e->seRefG()) + " / " + fmtd("%.1e", e->maxSeRefG/kTwoGMax),
                                     fmtd("%.1e", e->seFldG()) + " / " + fmtd("%.1e", e->maxSeFldG/kTwoGMax)});
            }
            out("    largest H error at X = %s, largest grad H error at X = %s\n", fmt3(eAB.worstH, 3).c_str(), fmt3(eAB.worstG, 3).c_str());

            // energies (1/2 int |grad H|^2): A from the Gram of (H0_m, P_v0, P_v1, P_v2) with the glued c, B exact
            auto energyA = [&](const Eigen::MatrixXd& G, const Eigen::Vector3d& c) {
                const int idx[4] = {m, 2, 3, 4};
                const double v[4] = {1.0, c[0], c[1], c[2]};
                double e = 0;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++) e += v[i]*v[j]*G(idx[i], idx[j]);
                return 0.5*e;
            };
            const Eigen::Vector3d cg = u.head<3>();
            const double EdA = energyA(gm.G, cg), EdB = 0.5*localB(u).dot(KB*localB(u));
            std::vector<double> rA, rB, rT;
            for (int a = 0; a < Bq; a++) {
                rA.push_back(energyA(gm.jack[a], uJ[a].head<3>()));
                rB.push_back(0.5*eJ[a].dot(KB*eJ[a]));
                rT.push_back(rA.back() + rB.back());
            }
            const double seDA = jackSe(rA), seDB = jackSe(rB), seD = jackSe(rT);
            const double ErA = 0.5*grA.G(m, m), ErB = 0.5*grB.G(m, m), Er = 0.5*grAB.G(m, m);
            const double seRA = 0.5*grA.se(m, m), seRB = 0.5*grB.se(m, m), seR = 0.5*grAB.se(m, m);
            const double dE = EdA + EdB - Er, seDE = std::hypot(seD, seR);
            const double halfE2 = 0.5*(eAB.eG2 - eAB.vRefG - eAB.vFldG);
            const double naiveDec = energyA(gm.naive, cg) - EdA, naiveRef = 0.5*(grAB.naive(m, m) - grAB.G(m, m));

            // the trace-fit diagnostic: the reference's trace on F fitted to the three hats, A run with that trace
            // (H0 + sum_k c_fit,k P_k, by superposition, check 1), against the reference
            const TraceFit tf = fitTrace(fps, refF);
            const TraceErr fitE = traceError(fps, refF, tf.c), glE = traceError(fps, refF, cg, &uJ);
            ErrStats rf;   // the reference's own noise on F
            for (size_t p = 0; p < ref.nF; p++) rf.add(fps.w[p], fps.X[p], 0.0, Vec3::Zero(), refF.vH[p], 0.0, Vec3::Zero(), Vec3::Zero());
            FieldAt fitA = runField(A, {RunTerm{&rec, {m == 0 ? 1.0 : 0.0, m == 1 ? 1.0 : 0.0, tf.c[0], tf.c[1], tf.c[2]}, {}}});
            addCoefNoise(fitA, rec, A, {2, 3, 4}, tf.cov);
            ErrStats fA, fNFA;
            for (size_t i = 0; i < ref.nA; i++) {
                fA.add(pa.w[i], pa.X[i], fitA, i, refA, i);
                if (std::fabs(pa.X[i].z()) < kNearF) fNFA.add(pa.w[i], pa.X[i], fitA, i, refA, i);
            }
            const double EfA = energyA(gm.G, tf.c);
            Eigen::Vector3d cfSe;
            for (int k = 0; k < 3; k++) cfSe[k] = std::sqrt(tf.cov(k, k));
            out("  energy E = (1/2) int |grad H|^2 (split-sample; jackknife s.e.):\n"
                "    A minus tube: decomposed %.5f +- %.5f, reference %.5f +- %.5f\n"
                "    B:            decomposed %.5f +- %.5f (exact given c, d), reference %.5f +- %.5f\n"
                "    A u B:        decomposed %.5f +- %.5f, reference %.5f +- %.5f; difference %+.5f +- %.5f\n"
                "    (1/2) int |grad(H_dec - H_ref)|^2 (noise-corrected) = %.5f (equal to the difference for any c, d, up to the\n"
                "    eps-tube term: a consistency check of the energy and error estimators, not a test of the glue)\n"
                "    naive - split: decomposed (A) %+.1e, reference %+.1e\n",
                EdA, seDA, ErA, seRA, EdB, seDB, ErB, seRB, EdA + EdB, seD, Er, seR, dE, seDE, halfE2, naiveDec, naiveRef);
            out("  trace-fit diagnostic: least-squares fit of the reference trace on F (%zu points, area-weighted) to the 3 hats:\n"
                "    c_fit = (%+.5f, %+.5f, %+.5f) +- (%.5f, %.5f, %.5f); glued c = (%+.5f, %+.5f, %+.5f)\n"
                "    residual of the fit on F, / 2 g_max: RMS %.2e (nc %.2e), max %.2e; the glued trace vs the reference on F: RMS\n"
                "    %.2e (nc %.2e, the noise of c included), max %.2e; the reference's s.e. on F: RMS %.1e, max %.1e\n"
                "    A with the fitted trace vs the reference: A: H RMS %.2e (nc %.2e), max %.2e; grad RMS %.2e (nc %.2e), max %.2e\n"
                "                                       near F (A): H RMS %.2e (nc %.2e), max %.2e; grad RMS %.2e (nc %.2e), max %.2e\n"
                "    (its noise: s.e. H RMS %.1e, max %.1e, including c_fit's); E_A with the fitted trace %.5f (reference %.5f)\n",
                fps.X.size(), tf.c[0], tf.c[1], tf.c[2], cfSe[0], cfSe[1], cfSe[2], cg[0], cg[1], cg[2], fitE.rms, fitE.nc,
                fitE.max, glE.rms, glE.nc, glE.max,
                rf.seRefH(), rf.maxSeRefH/kTwoGMax, fA.rmsH(), fA.corrH(), fA.maxH/kTwoGMax, fA.rmsG(), fA.corrG(), fA.maxG/kTwoGMax,
                fNFA.rmsH(), fNFA.corrH(), fNFA.maxH/kTwoGMax, fNFA.rmsG(), fNFA.corrG(), fNFA.maxG/kTwoGMax, fA.seFldH(),
                fA.maxSeFldH/kTwoGMax, EfA, ErA);

            // the replicate runs: sd of the glued c, d over the runs / rms jackknife s.e., and their mean against c_fit
            Eigen::Vector4d rMean = Eigen::Vector4d::Zero(), rSd = Eigen::Vector4d::Zero(), rRms = Eigen::Vector4d::Zero();
            Eigen::Matrix3d rCov = Eigen::Matrix3d::Zero();
            for (const auto& r : reps) {
                rMean += r[m].u/double(nRep);
                rRms += r[m].se.cwiseAbs2()/double(nRep);
            }
            for (const auto& r : reps) {
                rSd += (r[m].u - rMean).cwiseAbs2()/double(nRep - 1);
                rCov += (r[m].u - rMean).head<3>()*(r[m].u - rMean).head<3>().transpose()/double(nRep - 1);
            }
            rSd = rSd.cwiseSqrt();
            rRms = rRms.cwiseSqrt();
            const Eigen::Matrix3d meanCov = rCov/double(nRep);
            const Eigen::Vector3d dc = rMean.head<3>() - tf.c;
            Eigen::Vector3d dcSe, mSe;
            for (int k = 0; k < 3; k++) {
                mSe[k] = rSd[k]/std::sqrt(double(nRep));
                dcSe[k] = std::sqrt(meanCov(k, k) + tf.cov(k, k));
            }
            const TraceErr mE = traceError(fps, refF, rMean.head<3>(), nullptr, &meanCov);
            const Eigen::Vector4d ratio = rSd.cwiseQuotient(rRms);
            out("  replicates (%d runs): sd of (c_v0, c_v1, c_v2, d) over the runs / rms jackknife s.e. = %.2f %.2f %.2f %.2f\n"
                "    (1 expected; an sd from %d runs scatters by ~%.0f%%)\n"
                "    mean glued c over the runs = %s +- %s\n"
                "    mean glued c - c_fit       = %s +- %s\n"
                "    the mean c's trace vs the reference on F, / 2 g_max: RMS %.2e (nc %.2e), max %.2e\n", nRep, ratio[0], ratio[1],
                ratio[2], ratio[3], nRep, 100.0/std::sqrt(2.0*(nRep - 1)), fmtVec(rMean.head<3>()).c_str(), fmtVec(mSe, "%.5f").c_str(),
                fmtVec(dc).c_str(), fmtVec(dcSe, "%.5f").c_str(), mE.rms, mE.nc, mE.max);
            tRep.rows.push_back({cfg, gShort[m], std::to_string(nRep),
                                 fmtd("%.2f", ratio[0]) + " / " + fmtd("%.2f", ratio[1]) + " / " + fmtd("%.2f", ratio[2]) + " / " + fmtd("%.2f", ratio[3]),
                                 fmtVec(rMean.head<3>(), "%+.4f") + " +- " + fmtVec(mSe, "%.4f"),
                                 fmtVec(tf.c, "%+.4f") + " +- " + fmtVec(cfSe, "%.4f"),
                                 fmtVec(dc, "%+.4f") + " +- " + fmtVec(dcSe, "%.4f"),
                                 fmtd("%.2e", mE.rms) + " (" + fmtd("%.2e", mE.nc) + ") / " + fmtd("%.2e", mE.max)});

            // tables
            std::string evs;
            for (int k = 0; k < 3; k++) evs += (k ? ", " : "") + fmtd("%.4f", ev[k]) + " (" + fmtd("%.4f", evSe[k]) + ")";
            tMain.rows.push_back({cfg, gShort[m], fmtd("%.2e", eAB.rmsH()) + " / " + fmtd("%.2e", eAB.maxH/kTwoGMax),
                                  fmtd("%.2e", eAB.rmsG()) + " / " + fmtd("%.2e", eAB.maxG/kTwoGMax),
                                  fmtd("%.1e", eAB.seRefH()) + " / " + fmtd("%.1e", eAB.maxSeRefH/kTwoGMax) + "; " +
                                      fmtd("%.1e", eAB.seRefG()) + " / " + fmtd("%.1e", eAB.maxSeRefG/kTwoGMax),
                                  fmtd("%.1e", eAB.seFldH()) + " / " + fmtd("%.1e", eAB.maxSeFldH/kTwoGMax) + "; " +
                                      fmtd("%.1e", eAB.seFldG()) + " / " + fmtd("%.1e", eAB.maxSeFldG/kTwoGMax),
                                  fmtd("%.5f", EdA + EdB) + " / " + fmtd("%.5f", Er) + " / " + fmtd("%+.5f", dE) + " +- " + fmtd("%.5f", seDE),
                                  evs,
                                  fmtd("%+.4f", u[0]) + ", " + fmtd("%+.4f", u[1]) + ", " + fmtd("%+.4f", u[2]) + "; " + fmtd("%+.4f", u[3]) +
                                      " (" + fmtd("%.4f", uSe[0]) + ", " + fmtd("%.4f", uSe[1]) + ", " + fmtd("%.4f", uSe[2]) + "; " +
                                      fmtd("%.4f", uSe[3]) + ")",
                                  "A recipe " + fmtd("%.1f", rec.stats.seconds) + "; reference " + wallRef + "; " + zombieWall(gp, zr)});
            for (int rgn = 0; rgn < 2; rgn++) {
                const ErrStats &g = rgn ? eNFA : eA, &f = rgn ? fNFA : fA;
                tFit.rows.push_back({cfg, gShort[m],
                                     rgn ? "" : fmtVec(cg, "%+.4f") + " (" + fmtd("%.4f", uSe[0]) + ", " + fmtd("%.4f", uSe[1]) + ", " + fmtd("%.4f", uSe[2]) + ")",
                                     rgn ? "" : fmtVec(tf.c, "%+.4f") + " (" + fmtd("%.4f", cfSe[0]) + ", " + fmtd("%.4f", cfSe[1]) + ", " + fmtd("%.4f", cfSe[2]) + ")",
                                     rgn ? "" : fmtd("%.2e", fitE.rms) + " (" + fmtd("%.2e", fitE.nc) + ") / " + fmtd("%.2e", fitE.max),
                                     rgn ? "" : fmtd("%.2e", glE.rms) + " (" + fmtd("%.2e", glE.nc) + ") / " + fmtd("%.2e", glE.max),
                                     rgn ? "" : fmtd("%.1e", rf.seRefH()) + " / " + fmtd("%.1e", rf.maxSeRefH/kTwoGMax),
                                     rgn ? "near F (A, abs(z) < 0.1)" : "A", fmtd("%.2e", g.rmsH()) + " (" + fmtd("%.2e", g.corrH()) + ") / " + fmtd("%.2e", g.maxH/kTwoGMax),
                                     fmtd("%.2e", f.rmsH()) + " (" + fmtd("%.2e", f.corrH()) + ") / " + fmtd("%.2e", f.maxH/kTwoGMax),
                                     fmtd("%.2e", g.rmsG()) + " (" + fmtd("%.2e", g.corrG()) + ") / " + fmtd("%.2e", g.maxG/kTwoGMax),
                                     fmtd("%.2e", f.rmsG()) + " (" + fmtd("%.2e", f.corrG()) + ") / " + fmtd("%.2e", f.maxG/kTwoGMax),
                                     fmtd("%.1e", f.seFldH()) + " / " + fmtd("%.1e", f.maxSeFldH/kTwoGMax)});
            }
            tEn.rows.push_back({cfg, gShort[m], fmtd("%.5f", EdA) + " +- " + fmtd("%.5f", seDA), fmtd("%.5f", ErA) + " +- " + fmtd("%.5f", seRA),
                                fmtd("%.5f", EdB) + " +- " + fmtd("%.5f", seDB), fmtd("%.5f", ErB) + " +- " + fmtd("%.5f", seRB),
                                fmtd("%.5f", EdA + EdB) + " +- " + fmtd("%.5f", seD), fmtd("%.5f", Er) + " +- " + fmtd("%.5f", seR),
                                fmtd("%+.5f", dE) + " +- " + fmtd("%.5f", seDE), fmtd("%.5f", halfE2), fmtd("%.5f", EfA),
                                fmtd("%+.1e", naiveDec) + " / " + fmtd("%+.1e", naiveRef)});
            double dcMax = 0, dcSig = 0;
            for (int k = 0; k < 3; k++) {
                dcMax = std::max(dcMax, std::fabs(dc[k]));
                dcSig = std::max(dcSig, std::fabs(dc[k])/dcSe[k]);
            }
            cases.push_back({cfg, m ? "ramp" : "g = 1", eAB.rmsH(), eAB.corrH(), eAB.maxH/kTwoGMax,
                             std::hypot(eAB.seRefH(), eAB.seFldH()), eNF.rmsH(), eNF.corrH(), dE, seDE, fA.corrH(), eA.corrH(),
                             fitE.nc, glE.nc, mE.nc, dcMax, dcSig, eAB.rmsG(), eAB.corrG(), std::hypot(eAB.seRefG(), eAB.seFldG()),
                             zr.H[m].ncRmsDiffBound(kTwoGMax), zBoundWorst(zr, m), ratio.minCoeff(), ratio.maxCoeff()});
        }
        // the runs of this configuration
        auto runRow = [&](const std::string& run, const std::string& solver, int N, const SolveStats& st, const std::string& note) {
            tRun.rows.push_back({cfg, run, solver, std::to_string(N), std::to_string(st.nPoints), fmtd("%.2f", st.seconds),
                                 fmtd("%.1f", st.meanWalkLength()), std::to_string(st.walksCapped) + " / " + std::to_string(st.walksEscaped), note});
        };
        runRow("A recipe (H0 x2, P x3)", "own loop (double)", gp.nWalks, rec.stats,
               std::to_string(nRep - 1) + " replicate runs: " + fmtd("%.1f", repSec) + " s in all");
        runRow("reference on U (both g modes)", "own loop (double)", gp.nRef, ref.stats,
               ref.cached ? "loaded from gate_cache/ (wall clock of the original run)"
                          : "timing pilot N = " + std::to_string(ref.nPilot) + ": " + fmtd("%.1f", ref.pilotSeconds) + " s, extrapolated " +
                                fmtd("%.0f", ref.pilotExtrapolated) + " s");
        runRow("Zombie cross-check of the reference (H, grad H, both g modes)", "Zombie (float)", gp.nZombie, zr.st,
               "points, wall clock and walks summed over its " + std::to_string(2*(1 + gp.nZombieBatches)) +
                   " solves (grad: " + std::to_string(gp.nZombieBatches) + " x " + std::to_string(gp.nZombie/gp.nZombieBatches) +
                   " walks)");
    }
    out("\nCHECK 2: %s (the spec sets no pass threshold for checks 2 and 3: the numbers above are the result; the verdict covers\n"
        "the geometry assertions and the Zombie cross-check of the reference)\n", ok ? "done" : "FAILED");
    md += "### Main table (A u B)\n\n" + tMain.str() + "\n### Errors by region\n\nRegions: A (minus the eps-tube), B, A u B, near F = the "
          "points of A and B with abs(z) < 0.1.\n\n" + tReg.str() + "\n### Trace-fit diagnostic\n\nThe reference's trace on F "
          "(1024 points) least-squares fitted (area-weighted) to the 3 hats; A run with that trace (H0 + sum c_fit P_k, by the "
          "superposition of check 1) and compared with the reference. Its noise includes c_fit's (from the reference's noise "
          "on F). This separates \"the 3-value face trace is too poor\" (the fitted trace's error) from \"the glue / linear B is "
          "poor\" (the glued error beyond it).\n\n" + tFit.str() + "\n### Replicate recipe runs\n\nThe A recipe repeated with "
          "independent seeds (run 1 is the main run above). The sd of the glued c, d over the runs against one run's "
          "jackknife s.e. checks the method behind the decomposed noise floor (above 1: the floor is too low; below 1: it is "
          "conservative, and the nc values above are then lower estimates). The mean over the runs is sqrt(runs) less noisy "
          "than one run: its offset from c fit is the glue's systematic error in the corner values, which one run's noise "
          "can hide.\n\n" + tRep.str() + "\n### Energies\n\nE = (1/2) int abs(grad H)^2.\n\n" +
          tEn.str() + refZombieMd(gp, tZ, tZr) + "\n### Runs\n\n" + tRun.str() + "\n### Quadrature of the point sets\n\n" + tQ.str() +
          "\n" + std::string("**Check 2: ") + (ok ? "done" : "FAILED") +
          "** (no pass threshold in the spec; the verdict covers the geometry assertions and the reference's Zombie cross-check)\n";

    // findings: generated from the cases above
    std::string f = std::string("- Check 2 (two tets, B uncut and linear, the crack z = z_c does not reach F): ") + (ok ? "done" : "FAILED") +
                    ". The spec sets no pass threshold; the geometry assertions " + (ok ? "hold" : "or the cross-check FAILED") +
                    ", and Zombie agrees with the reference in " + std::to_string(nZagree) + " of " + std::to_string(nZcases) +
                    " cases.\n";
    f += "- The glued decomposition against the reference, over A u B, divided by 2 g_max = 2. noise = the combined RMS s.e. of H, "
         "sqrt(ref^2 + dec^2), where dec includes the noise of the glued c and d; error/noise ~ 1 means the error is not resolved "
         "above the Monte Carlo noise of one run.\n\n";
    std::vector<std::vector<std::string>> rows = {{"configuration", "g", "H RMS (nc)", "H max", "noise", "error/noise",
                                                   "near F: H RMS (nc)", "E dec - E ref (s.e.)", "ref check: bound (worst region)"}};
    std::string resolved, unresolved;
    double dEsig = 0, calLo = kInf, calHi = 0;
    for (const Case& c : cases) {
        rows.push_back({c.cfg, c.g, pct(c.rms) + " (" + pct(c.nc) + ")", pct(c.max), pct(c.noise), fmtd("%.2f", c.rms/c.noise),
                        pct(c.nearRms) + " (" + pct(c.nearNc) + ")", fmtd("%+.5f", c.dE) + " (" + fmtd("%.5f", c.seDE) + ")",
                        pct(c.zAll) + " (" + pct(c.zWorst) + ")"});
        calLo = std::min(calLo, c.calLo);
        calHi = std::max(calHi, c.calHi);
        std::string& list = c.nc > 0 ? resolved : unresolved;
        list += (list.empty() ? "" : ", ") + c.cfg + " " + c.g;
        dEsig = std::max(dEsig, std::fabs(c.dE)/c.seDE);
    }
    f += fixedTable(rows) + "\n";
    double gLo = kInf, gHi = 0, gnLo = kInf, gnHi = 0, gsLo = kInf, gsHi = 0;
    for (const Case& c : cases) {
        gLo = std::min(gLo, c.gRms);
        gHi = std::max(gHi, c.gRms);
        gnLo = std::min(gnLo, c.gNc);
        gnHi = std::max(gnHi, c.gNc);
        gsLo = std::min(gsLo, c.gNoise);
        gsHi = std::max(gsHi, c.gNoise);
    }
    f += "- Resolved above the noise of one run (nc > 0): " + (resolved.empty() ? std::string("none") : resolved) +
         ". Within it (nc = 0): " + (unresolved.empty() ? std::string("none") : unresolved) + ". abs(E dec - E ref) is at most " +
         fmtd("%.2f", dEsig) + " s.e." + (dEsig < 3 ? ", so the energy comparison does not resolve the error" : "") +
         ". grad H over A u B: RMS error " + pct(gLo, "%.1f%%") + " to " + pct(gHi, "%.1f%%") + " (nc " + pct(gnLo, "%.1f%%") +
         " to " + pct(gnHi, "%.1f%%") + ") against a noise RMS of " + pct(gsLo, "%.1f%%") + " to " + pct(gsHi, "%.1f%%") + ".\n";
    f += "- What \"Zombie agrees\" can see: the last column is the largest RMS discrepancy between the reference and Zombie that "
         "their comparison does not exclude (2-s.e. bound, / 2 g_max, over all points; the worst of the regions A, B, F, near F "
         "and near the crack in brackets). Where the H error is not well above it, the cross-check does not confirm the "
         "reference at the level of the error. Noise calibration: over the replicate runs the sd of the glued c and d is " +
         fmtd("%.2f", calLo) + " to " + fmtd("%.2f", calHi) + " times their jackknife s.e. (1 expected; the sd of " +
         std::to_string(gp.quick ? 3 : 8) + " runs scatters by ~" + fmtd("%.0f", 100.0/std::sqrt(2.0*((gp.quick ? 3 : 8) - 1))) +
         "%" + (calHi < 1 ? "; below 1 the decomposed noise floor is conservative, and nc a lower estimate of the error" : "") +
         ").\n";
    f += "- Trace-fit diagnostic: A run with the best 3-hat fit of the reference trace on F, against the glued run. On F, the fit's "
         "residual against the error of the glued trace with c averaged over the replicate runs (one run's c is too noisy to "
         "show it); \"trace space\" = the fit residual's share of that error in mean square, the rest comes from the glue with "
         "the linear B. All nc, / 2 g_max.\n\n";
    rows = {{"configuration", "g", "A: H RMS, fitted / glued", "F: fit residual / glued mean c", "trace space",
             "mean glued c - c fit: max (s.e.)"}};
    double shLo = kInf, shHi = 0, sgLo = kInf, sgHi = 0;
    for (const Case& c : cases) {
        const double share = c.meanF > 0 ? c.fitF*c.fitF/(c.meanF*c.meanF) : 0.0;
        rows.push_back({c.cfg, c.g, pct(c.fitA) + " / " + pct(c.gluA), pct(c.fitF) + " / " + pct(c.meanF),
                        c.meanF > 0 ? pct(share, "%.0f%%") : "-", fmtd("%.4f", c.dcMax) + " (" + fmtd("%.1f", c.dcSig) + ")"});
        shLo = std::min(shLo, share);
        shHi = std::max(shHi, share);
        sgLo = std::min(sgLo, c.dcSig);
        sgHi = std::max(sgHi, c.dcSig);
    }
    f += fixedTable(rows) + "\n";
    f += "- The 3-hat trace space accounts for " + pct(shLo, "%.0f%%") + " to " + pct(shHi, "%.0f%%") + " of the glued trace error on "
         "F, the glue with the linear B for the rest" + (shHi < 0.5 ? ", so the glue / linear-B approximation is the larger part in "
         "every case" : shLo > 0.5 ? ", so the 3-value trace is the larger part in every case" : "") + ". In the replicate mean the "
         "glued c sits " + fmtd("%.1f", sgLo) + " to " + fmtd("%.1f", sgHi) + " s.e. from the fitted c.\n";
    md.insert(mdIntro, findingsBlock(f));
    geoMd = "## Appendix C: geometry of check 2, per configuration\n\nPrinted by the run before each configuration's first solve.\n\n" +
            codeBlock(geoTxt);
    return ok;
}

// =====================================================================================
// Check 3: two tets, both cut, the front crosses F
// =====================================================================================
// Level B's face function on F (physical X, z = 0), the same in A and in B, for the crack plane x = d:
//   J = sigma g_face bump fade,  sigma = +1 on the +g side (x > d), -1 on the other; g_face = the crack's g at the same
//   y (U's cut at (d, y, 0), so r = y_f - y); bump = max(0, 1 - |x - d|/w), w = 0.25; fade = sqrt(max(0, (y_f - y)/y_f))
//   for g = 1 (the ramp's g_face already goes to 0 at the front), 1 for the ramp; J = 0 for y >= y_f.
// J is not harmonic and has no known harmonic extension, so the runs that carry it cannot use the odd mirror ball: their
// first ball at F is the plain one (FieldSpec::faceHarmonic = false), the mirror ball across the Neumann faces only.
// With tip (--face-trace tip) J is instead the crack's own 2D solution around the front line, in the plane normal to it:
// a = y - y_f along the advance, b = x - d across the crack, theta = atan2(b, a), r = |(a, b)|, and
//   J = theta/pi (g = 1),   J = min(1, sqrt(r/L)) sin(theta/2) (ramp).
// Both are +-g on the crack (theta = +-pi), so the face data meets the crack data along their common edge on F (no
// log-divergent edge), both are 0 straight ahead of the front, and both depend on X only, so A and B see the same J.
constexpr double kBumpWidth = 0.25;
FaceFn levelBTrace(const CrackSpec& cs, GType gt, bool tip)
{
    const Cut cu = unionCut(cs);
    if (tip)
        return [L = cu.L, d = cs.d, yf = cs.yFront, gt](const Vec3& X) {
            const double a = X.y() - yf, b = X.x() - d, theta = std::atan2(b, a);
            if (gt == GType::One) return theta/std::numbers::pi;
            return std::min(1.0, std::sqrt(std::hypot(a, b)/L))*std::sin(0.5*theta);
        };
    return [cu, d = cs.d, yf = cs.yFront, gt](const Vec3& X) {
        const double x = X.x(), y = X.y();
        if (y >= yf) return 0.0;
        const double sigma = x > d ? 1.0 : x < d ? -1.0 : 0.0;
        const double bump = std::max(0.0, 1.0 - std::fabs(x - d)/kBumpWidth);
        const double fade = gt == GType::One ? std::sqrt(std::max(0.0, (yf - y)/yf)) : 1.0;
        return sigma*cu.g(Vec3(d, y, 0.0), gt)*bump*fade;
    };
}

// The F point set of a check-3 configuration: the fixed set without the points within eps of the crack's segment on F
// (a-b); the weights are unchanged, so the rule integrates over F minus that band.
FacePointSet facePointsOff(const FacePointSet& fps, const Vec3& a, const Vec3& b, double eps, int& nDropped)
{
    FacePointSet o;
    nDropped = 0;
    for (size_t p = 0; p < fps.X.size(); p++) {
        if (pointSegmentDist(fps.X[p], a, b) < eps) {
            nDropped++;
            continue;
        }
        o.X.push_back(fps.X[p]);
        o.w.push_back(fps.w[p]);
    }
    return o;
}

// The Gram of one tet's fields in check 3 over its points with weights w: the recipe run's 5 (H0 g1, H0 g2: crack +-g,
// F 0; P_v0, P_v1, P_v2: crack 0, F hat_k) at indices 0..4 and the Level-B H0 run's 2 (crack +-g, F J_g) at kH0B + 0, 1.
// The two runs are independent and their batches are paired by index, so every product is one of independent
// estimates and batch a of both runs is one unit of the jackknife.
constexpr int kH0B = 5;
Gram gramCheck3(const OwnRun& rec, const OwnRun& h0b, const std::vector<double>& w, const TetFrame& T)
{
    const int K1 = rec.sums.K, K2 = h0b.sums.K;
    if (rec.sums.nBatches != h0b.sums.nBatches || rec.sums.nPoints != h0b.sums.nPoints)
        throw std::logic_error("gramCheck3: the runs differ in batches or points");
    std::vector<std::vector<double>> u1, u2;
    for (int k = 0; k < K1; k++) u1.push_back(unitCoef(K1, k));
    for (int k = 0; k < K2; k++) u2.push_back(unitCoef(K2, k));
    return gramFromBatches(w, rec.sums.nBatches, K1 + K2, [&](size_t p, int b, Vec3* g) {
        for (int q = 0; q < K1; q++) g[q] = T.JinvT*comboBatchGrad(rec.sums, p, b, u1[q]);
        for (int q = 0; q < K2; q++) g[K1 + q] = T.JinvT*comboBatchGrad(h0b.sums, p, b, u2[q]);
    });
}

// The glue of check 3, written out. Unknowns: c = (c_v0, c_v1, c_v2), the values at the corners of F (global ids),
// shared by A and B. Each tet's P_k answers the same physical hat hat_k, evaluated at X = X0 + J xi in the tet's own
// frame (in B local corner 1 is v2 and 2 is v1: hat_1(J_B xi) = xi_y = B's lambda_2), so S_A, S_B, b_A and b_B are all
// indexed by the global corner and add directly. With H0 the tet's Level-A or Level-B particular solution
//   E(c) = E_A(c) + E_B(c) = (1/2) c^T (S_A + S_B) c + (b_A + b_B)^T c + const,
// and dE/dc = 0 gives the 3 x 3 system
//   (S_A + S_B) c = -(b_A + b_B).
// S_A and S_B are symmetrized first; both are positive definite (the crack is Dirichlet in both tets). Eigen, double.
Eigen::Vector3d glueCheck3(const Eigen::Matrix3d& SA, const Eigen::Matrix3d& SB, const Eigen::Vector3d& bA,
                           const Eigen::Vector3d& bB, Eigen::Matrix3d* Mout = nullptr, Eigen::Vector3d* rhsOut = nullptr)
{
    const Eigen::Matrix3d M = 0.5*(SA + SA.transpose()) + 0.5*(SB + SB.transpose());   // S_A + S_B
    const Eigen::Vector3d rhs = -(bA + bB);                                               // -(b_A + b_B)
    if (Mout) *Mout = M;
    if (rhsOut) *rhsOut = rhs;
    return M.fullPivLu().solve(rhs);
}

// The glued c of check 3 from the two tets' Grams (gramCheck3) with the particular solution at index h (0, 1: Level A;
// kH0B + 0, 1: Level B), its delete-one-batch replicates (the system re-solved with replicate a of both Grams) and their
// jackknife s.e.
struct Glue3 {
    Eigen::Vector3d c = Eigen::Vector3d::Zero(), se = Eigen::Vector3d::Zero();
    std::vector<Eigen::Vector3d> jack;
};
Glue3 glueRun3(const Gram& gA, const Gram& gB, int h, Eigen::Matrix3d* M = nullptr, Eigen::Vector3d* rhs = nullptr)
{
    auto S = [](const Eigen::MatrixXd& G) { return Eigen::Matrix3d(G.block(2, 2, 3, 3)); };
    auto b = [h](const Eigen::MatrixXd& G) { return Eigen::Vector3d(G(h, 2), G(h, 3), G(h, 4)); };
    Glue3 gl;
    gl.c = glueCheck3(S(gA.G), S(gB.G), b(gA.G), b(gB.G), M, rhs);
    for (size_t a = 0; a < gA.jack.size(); a++)
        gl.jack.push_back(glueCheck3(S(gA.jack[a]), S(gB.jack[a]), b(gA.jack[a]), b(gB.jack[a])));
    for (int k = 0; k < 3; k++) {
        std::vector<double> rep;
        for (const Eigen::Vector3d& x : gl.jack) rep.push_back(x[k]);
        gl.se[k] = jackSe(rep);
    }
    return gl;
}

// (1/2) int |grad(H0_h + sum_k c_k P_k)|^2 over a tet's points, from its Gram (gramCheck3 layout)
double tetEnergy(const Eigen::MatrixXd& G, int h, const Eigen::Vector3d& c)
{
    const int idx[4] = {h, 2, 3, 4};
    const double v[4] = {1.0, c[0], c[1], c[2]};
    double e = 0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) e += v[i]*v[j]*G(idx[i], idx[j]);
    return 0.5*e;
}

constexpr double kNearSeg = 0.1;   // check 3: "near the crack's segment on F" (F points and volume points)

// Level B's H0 fields (crack +-g, F = J_g for both g types; J has no harmonic extension, so the plain first ball at F)
std::vector<FieldSpec> levelBFields(const FaceFn J[2])
{
    return {FieldSpec{"H0B g1", true, GType::One, J[0], false}, FieldSpec{"H0B g2", true, GType::SqrtFront, J[1], false}};
}

// the crack's segment on F in X (A's polygon vertices on xi_z = 0, i.e. z = 0), ordered by y
std::array<Vec3, 2> crackSegmentOnF(const TetFrame& A, const XiCrack& xa)
{
    std::vector<Vec3> seg;
    for (const Vec3& p : xa.cut.polygon)
        if (std::fabs(p.z()) <= 1e-12) seg.push_back(A.toX(p));
    if (seg.size() != 2) throw std::logic_error("check 3: the crack does not meet F in a segment");
    if (seg[0].y() > seg[1].y()) std::swap(seg[0], seg[1]);
    return {seg[0], seg[1]};
}

// The decomposed field of check 3 at level L (0: H0 of the recipe run rec; 1: H0 of the Level-B run hb), g type m, in
// tet T with coefficients c (cJ: their jackknife replicates, which enter the field's variance; null: fixed coefficients)
FieldAt decField3(const TetFrame& T, const OwnRun& rec, const OwnRun& hb, int L, int m, const Eigen::Vector3d& c,
                  const std::vector<Eigen::Vector3d>* cJ)
{
    std::vector<std::vector<double>> pj;
    if (cJ)
        for (const Eigen::Vector3d& x : *cJ) pj.push_back({L == 0 && m == 0 ? 1.0 : 0.0, L == 0 && m == 1 ? 1.0 : 0.0, x[0], x[1], x[2]});
    if (L == 0) return runField(T, {RunTerm{&rec, {m == 0 ? 1.0 : 0.0, m == 1 ? 1.0 : 0.0, c[0], c[1], c[2]}, pj}});
    return runField(T, {RunTerm{&hb, unitCoef(2, m), {}}, RunTerm{&rec, {0.0, 0.0, c[0], c[1], c[2]}, pj}});
}

bool runCheck3(Report& out, const GateParams& gp, const FacePointSet& fpsAll, const std::string& cacheDir, std::string& md,
               std::string& geoMd)
{
    const TetFrame A = frameA(), B = frameB();
    const GType gts[2] = {GType::One, GType::SqrtFront};
    const char* gShort[2] = {"(a) g = 1", "(b) sqrt ramp"};
    const char* lvShort[2] = {"A (J = 0)", "B (J)"};
    out("\n=== Check 3: two tets, both cut, the front crosses F ===\n");
    out("crack plane x = 0.3 (+g for x > 0.3), front {x = 0.3, y = y_f}, crack region y < y_f, y_f = 0.15, 0.30 and 0.45: a strip\n"
        "through both tets whose front pierces F at (0.3, y_f, 0). Both tets run the recipe in their own xi frame (B's crack pulled\n"
        "back through J_B^-1): H0 (crack +-g, F 0) for both g types and P_v0, P_v1, P_v2 (crack 0, F hat_k), all from the same\n"
        "walks (own loop, N = %d in %d batches). Level A: J = 0. Level B: one more run per tet, H0 with F = J (both g types from\n"
        "the same walks); the P_k are Level A's, reused. Glue: (S_A + S_B) c = -(b_A + b_B), c = the values at v0, v1, v2.\n"
        "Reference: one own-loop solve on U (N = %d), cached in %s. Errors / 2 g_max = 2 in both g types; RMS quadrature-\n"
        "weighted, max over the points; energies (1/2) int |grad H|^2 over A and B minus their eps-tubes, split-sample.\n",
        gp.nWalks, gp.nBatches, gp.nRef, cacheDir.c_str());
    md += "## Check 3: two tets, both cut, the front crosses F\n\n";
    md += "Crack plane x = 0.3 (+g for x > 0.3), front {x = 0.3, y = y_f}, crack region y < y_f, y_f = 0.15, 0.30, 0.45: a strip "
          "through both tets whose front pierces F at (0.3, y_f, 0). Both tets run the full recipe in their own xi frame (B's "
          "crack pulled back through J_B^-1; own loop, double, N = " + std::to_string(gp.nWalks) + " in " +
          std::to_string(gp.nBatches) + " batches): H0 (crack +-g, F 0) for both g modes and P_v0, P_v1, P_v2 (crack 0, F hat_k), "
          "all from the same walks. **Level A**: J = 0, the face trace is forced continuous. **Level B**: H0 with F = J = sigma "
          "g_face bump fade (sigma = +-1 by the side of x = 0.3, g_face = the crack's g at the same y, bump = max(0, 1 - abs(x - "
          "0.3)/0.25), fade = sqrt(max(0, (y_f - y)/y_f)) for g = 1 only, J = 0 for y >= y_f); the same physical J in A and B, one "
          "more run per tet (both g modes from the same walks); the P_k are Level A's, reused. J has no harmonic extension, so "
          "its runs use the plain first ball at F (the P_k and Level A's H0 use the odd mirror ball). Glue: (S_A + S_B) c = -(b_A "
          "+ b_B) for the 3 shared corners (global ids), written out in `glueCheck3`. Reference: one own-loop solve on U = A u B "
          "in X (all four faces Neumann, the crack clipped to U, L = sqrt(2)), N = " + std::to_string(gp.nRef) + " in " +
          std::to_string(gp.nRefBatches) + " batches, cached in gate_cache/.\n\n"
          "All errors and noise levels are divided by **2 g_max = 2** in both g modes (the ramp's local 2g vanishes at the "
          "front). RMS is quadrature-weighted over the points, max is over the points; grad H errors are vector norms. \"nc\" = "
          "noise-corrected RMS, sqrt(max(0, mean e^2 - mean(var_ref + var_dec)))/2. The decomposed field's noise is the "
          "delete-one-batch jackknife over the runs' batches (batch a of all four runs of a configuration is one unit), with the "
          "glued c recomputed in every replicate. Energies E = (1/2) int abs(grad H)^2 over A and B minus their eps-tubes, from "
          "split-sample products, jackknife s.e. \"Near F\" = the points of A and B with abs(z) < 0.1; \"near seg\" = the points "
          "(and on F the F points) within 0.1 of the crack's segment on F, {x = 0.3, 0 <= y <= y_f, z = 0}.\n\n";
    const size_t mdIntro = md.size();   // the findings go here
    struct Case {                       // the numbers behind the findings, per configuration, level and g mode
        std::string yf, lv, g;
        double rms, nc, max, noise, segRms, Ed, Er, seDE, fitAB, gluAB, fitF, meanF, fitNear, slab, slabSe, edge, dqT, gRms, gNc,
            gNoise, zAll, zWorst, calLo, calHi;   // the reference check's bound (all points, worst region); replicate sd / jackknife
    };
    std::vector<Case> cases;
    int nZagree = 0, nZcases = 0;
    double evMinLo = kInf, evMinHi = 0, evMinSe = 0, quadRel = 0, quadSig = 0;   // over the configurations
    std::string geoTxt;   // the geometry printout of every configuration (appendix)
    MdTable tMain, tReg, tTr, tFit, tRep, tEn, tS, tZ, tZr, tRun, tQ;
    tMain.head = {"y_f", "level", "g mode", "H: RMS / max (A u B)", "grad H: RMS / max (A u B)",
                  "reference noise: s.e. H RMS / max; grad H RMS / max", "decomposed noise: s.e. H RMS / max; grad H RMS / max",
                  "energy E (A u B): decomposed / reference / difference", "S_A + S_B eigenvalues (s.e.)",
                  "glued c_v0, c_v1, c_v2 (s.e.)", "wall clock per WoSt run (s)"};
    tReg.head = {"y_f", "level", "g mode", "region", "points", "H RMS (nc)", "H max", "grad H RMS (nc)", "grad H max",
                 "ref s.e. H RMS / max", "dec s.e. H RMS / max", "ref s.e. grad RMS / max", "dec s.e. grad RMS / max"};
    tTr.head = {"y_f", "level", "g mode", "c glued (s.e.)", "c fit (s.e.)", "fit residual on F: RMS (nc) / max",
                "fit residual: near seg RMS (nc) / far RMS (nc); near seg share of int e^2 (nc)",
                "glued trace vs reference on F: RMS (nc) / max",
                "glued trace: near seg RMS (nc) / far RMS (nc); near seg share of int e^2 (nc)", "ref s.e. on F: RMS / max"};
    tFit.head = {"y_f", "level", "g mode", "region", "glued: H RMS (nc) / max", "fitted trace: H RMS (nc) / max",
                 "glued: grad RMS (nc) / max", "fitted trace: grad RMS (nc) / max", "fitted-trace noise: s.e. H RMS / max"};
    tRep.head = {"y_f", "level", "g mode", "runs", "sd over the runs / rms jackknife s.e.: c_v0 / c_v1 / c_v2",
                 "mean glued c over the runs +- s.e.", "c fit +- s.e.", "mean glued c - c fit +- s.e.",
                 "mean glued c's trace vs reference on F: RMS (nc) / max; near seg share (nc)",
                 "quadrature: c - c with B on A's mirrored points +- s.e.; its trace RMS on F"};
    tEn.head = {"y_f", "level", "g mode", "E_A dec", "E_A ref", "E_B dec", "E_B ref", "E dec", "E ref", "E dec - E ref",
                "(1/2) int abs(grad e)^2 (nc)", "near seg: E dec / E ref", "E with the fitted trace", "naive - split: dec / ref",
                "slab eps to 2 eps: E dec - E ref (edge model)"};
    tS.head = {"y_f", "S_A eigenvalues (s.e.)", "S_B eigenvalues (s.e.)", "S_A + S_B eigenvalues (s.e.)",
               "mirror check on common points (B on A's mirrored points vs A): max abs(dG)/s.e. S (6) / b (12) / H0 energies (4) "
               "(largest rel. dS)",
               "B's own points vs A's: the same", "quadrature, B's own points vs A's mirrored: the same"};
    tZ.head = {"y_f", "g mode", "H: frac abs(z)<3, pooled z (mean dH +- s.e.)",
               "grad H (points > 0.05 from U's faces): pooled z per component (frac abs(z)<3)", "verdict"};
    tZr.head = {"y_f", "g mode", "region", "points", "H pooled z", "mean z^2 (s.e.)",
                "RMS discrepancy beyond noise / 2 g_max (2-s.e. bound)"};
    tRun.head = {"y_f", "run", "solver", "N per point", "points", "wall clock (s)", "mean walk length", "capped / escaped", "note"};
    tQ.head = {"point set", "kept points", "dropped (eps)", "max rel err, all points vs exact over T",
               "max rel err, kept vs T minus tube", "tube volume"};
    std::string sMd;   // S_A, S_B and the b vectors per configuration (code block)

    bool ok = true;
    const double yfs[3] = {0.15, 0.30, 0.45};
    for (int ci = 0; ci < 3; ci++) {
        const CrackSpec cs = crackCheck3(yfs[ci]);
        if (cs.n != Vec3(1, 0, 0)) throw std::logic_error("check 3: levelBTrace assumes the crack plane x = d");
        const std::string cfg = fmtd("%.2f", yfs[ci]);
        const size_t geo0 = out.text.size();
        out("\n--- check 3, y_f = %s ---\n", cfg.c_str());
        const XiCrack xa = xiCrack(A, cs), xb = xiCrack(B, cs);
        const StripCut usc = unionStrips(cs);
        const Cut cu = unionCut(cs);
        const FaceFn J[2] = {levelBTrace(cs, GType::One, gp.tipTrace), levelBTrace(cs, GType::SqrtFront, gp.tipTrace)};
        const std::vector<std::string> runs = {
            "reference U (own loop; Zombie cross-check): " + faceFlags(makeFrame("U", {1, 2, 3, 4}), allNeumann(), "", false) +
                "; crack +-g (both g types)",
            "A recipe (Level A; H0 g1, H0 g2: crack +-g, F 0; P_v0, P_v1, P_v2: crack 0, F hat_k): " + faceFlags(A, onlyF(), "(data per field)"),
            "B recipe (the same fields in B's xi frame; the hats and J evaluated at X = J_B xi): " + faceFlags(B, onlyF(), "(data per field)"),
            "A Level-B H0 (H0B g1, H0B g2: crack +-g, F J_g; plain first ball at F): " + faceFlags(A, onlyF(), "H = J_g"),
            "B Level-B H0 (the same physical J_g): " + faceFlags(B, onlyF(), "H = J_g"),
            "Level B's P_k: the A and B recipes' P_v0, P_v1, P_v2, reused (not re-solved)"};
        if (!printRunGeometry(out, gp, cs, xa, xb, usc, runs, 9000 + ci)) {
            geoTxt += out.text.substr(geo0);
            out("geometry assertions FAILED: check 3 at y_f = %s not run\n", cfg.c_str());
            ok = false;
            continue;
        }
        if (xa.cut.empty() || xb.cut.empty()) throw std::logic_error("check 3: the crack must cut both tets");
        const std::array<Vec3, 2> seg = crackSegmentOnF(A, xa);
        auto nearSegment = [&](const Vec3& X) { return pointSegmentDist(X, seg[0], seg[1]) < kNearSeg; };

        // Level B's J: the same physical function through both frames, and its values at the crack's edge on F
        {
            std::mt19937_64 rng(mixSeed(gp.seed, 9060 + ci));
            std::uniform_real_distribution<double> Ud(0.0, 1.0);
            double dFrame = 0;
            for (int k = 0; k < 2000; k++) {
                double a = Ud(rng), b = Ud(rng);
                if (a + b > 1) {
                    a = 1 - a;
                    b = 1 - b;
                }
                const Vec3 X(a, b, 0.0);
                for (int m = 0; m < 2; m++)
                    for (const TetFrame* T : {&A, &B}) dFrame = std::max(dFrame, std::fabs(J[m](T->toX(T->toXi(X))) - J[m](X)));
            }
            double edge[2] = {0, 0};
            const double xp = std::nextafter(cs.d, 1.0), xm = std::nextafter(cs.d, 0.0);
            for (int k = 0; k < 100; k++) {
                const double y = cs.yFront*k/100.0;
                for (int m = 0; m < 2; m++) {
                    const double gy = cu.g(Vec3(cs.d, y, 0.0), gts[m]);
                    edge[m] = std::max({edge[m], std::fabs(J[m](Vec3(xp, y, 0.0)) - gy), std::fabs(J[m](Vec3(xm, y, 0.0)) + gy)});
                }
            }
            const bool jOk = dFrame < 1e-15;
            ok = ok && jOk;
            const double ym = 0.5*cs.yFront;
            if (gp.tipTrace)
                out("  --face-trace tip: Level B's J below is replaced by the crack's 2D front solution, J = theta/pi (g = 1) and\n"
                    "  min(1, sqrt(r/L)) sin(theta/2) (ramp), theta = atan2(x - %.2f, y - y_f), r = |(x - %.2f, y - y_f)|; the edge line\n"
                    "  below then measures this J\n", cs.d, cs.d);
            out("  Level B face function on F (physical X): J = sigma g_face bump fade; sigma = +1 for x > %.2f, -1 below; g_face = the\n"
                "  crack's g at (%.2f, y, 0) (r = y_f - y); bump = max(0, 1 - |x - %.2f|/%.2f); fade = sqrt(max(0, (y_f - y)/y_f)) for g = 1,\n"
                "  1 for the ramp; J = 0 for y >= y_f. The same physical J in A and B: max |J(X0 + J_T xi(X)) - J(X)| over 2000 random\n"
                "  points of F, both frames and both g types = %.1e -> %s. Samples at y = y_f/2: J(0.30 +- 0, .) = %+.4f / %+.4f (g = 1),\n"
                "  %+.4f / %+.4f (ramp); J(0.425, .) = %+.4f (g = 1), %+.4f (ramp). Along the crack's edge on F (x = 0.30 +- 0, 0 <= y < y_f)\n"
                "  max |J - (+-g)| = %.3f for g = 1 (J = +-sqrt((y_f - y)/y_f) there, +-1 only at y = 0) and %.1e for the ramp (J = +-g)\n",
                cs.d, cs.d, cs.d, kBumpWidth, dFrame, jOk ? "OK" : "FAILED", J[0](Vec3(xp, ym, 0.0)), J[0](Vec3(xm, ym, 0.0)),
                J[1](Vec3(xp, ym, 0.0)), J[1](Vec3(xm, ym, 0.0)), J[0](Vec3(cs.d + 0.125, ym, 0.0)), J[1](Vec3(cs.d + 0.125, ym, 0.0)),
                edge[0], edge[1]);
        }
        geoTxt += out.text.substr(geo0);

        // point sets
        const PointSet pa = volumePoints(A, xa, gp.nPoints, gp.eps, gp.quick ? 1 : 4);
        const PointSet pb = volumePoints(B, xb, gp.nPoints, gp.eps, gp.quick ? 1 : 4);
        for (const auto& [T, x, ps] : {std::tuple<const TetFrame*, const XiCrack*, const PointSet*>{&A, &xa, &pa}, {&B, &xb, &pb}}) {
            const std::string label = std::string(T->name) + ", y_f = " + cfg;
            out("point set %s: %d leaves (depth 2/3/4: %d/%d/%d), %zu generated, %zu kept\n", label.c_str(), ps->nLeaves,
                ps->leavesByDepth[2], ps->leavesByDepth[3], ps->leavesByDepth[4], ps->allXi.size(), ps->size());
            quadratureCheck(out, *T, *x, *ps, gp.nTube, &tQ, label);
        }
        int nFDrop = 0;
        const FacePointSet fps = facePointsOff(fpsAll, seg[0], seg[1], gp.eps, nFDrop);
        std::vector<char> nearF(fps.X.size());
        double swF = 0, swNear = 0;
        for (size_t p = 0; p < fps.X.size(); p++) {
            nearF[p] = nearSegment(fps.X[p]);
            swF += fps.w[p];
            swNear += nearF[p] ? fps.w[p] : 0.0;
        }
        const double band = 2.0*gp.eps*(seg[1] - seg[0]).norm() + 0.5*std::numbers::pi*gp.eps*gp.eps;
        out("F point set: %zu of %zu kept (%d within eps of the crack's segment on F, %s - %s, dropped); sum of weights %.6f vs\n"
            "  area of F minus that band %.6f (band = 2 eps y_f + pi eps^2/2); %.1f%% of the kept area lies within %.2f of the segment\n",
            fps.X.size(), fpsAll.X.size(), nFDrop, fmt3(seg[0], 2).c_str(), fmt3(seg[1], 2).c_str(), swF, 0.5 - band, 100*swNear/swF,
            kNearSeg);
        tQ.rows.push_back({"F, y_f = " + cfg, std::to_string(fps.X.size()), std::to_string(nFDrop), "-",
                           fmtd("%.2e", (swF - (0.5 - band))/(0.5 - band)) + " (sum w vs area of F minus the band)", fmtd("%.6f", band)});

        // the reference and its Zombie cross-check
        const RefRun ref = runReference(out, gp, cs, "check3_yf" + cfg, 9100 + ci, pa, pb, fps, cacheDir);
        if (!ref.ok) {
            ok = false;
            continue;
        }
        const RefZombie zr = zombieReference(gp, cs, ref, 9300 + 100*ci);
        reportRefZombie(out, gp, zr, cfg, tZ, tZr);
        ok = ok && zr.agree();
        for (int m = 0; m < 2; m++) {
            nZcases++;
            nZagree += zr.agreeH(m) && zr.agreeG(m);
        }

        // the recipes: Level A in both tets (H0 for both g types and the P_k from the same walks), and Level B's H0 in
        // both tets (one more run each, both g types from the same walks; the P_k are reused)
        const std::vector<FieldSpec> fields = recipeFields(), fieldsB = levelBFields(J);
        auto runSet = [&](uint64_t tagA, uint64_t tagB, uint64_t tagHA, uint64_t tagHB) {
            return std::array<OwnRun, 4>{runOwn(A, xa, onlyF(), fields, pa.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, tagA)),
                                         runOwn(B, xb, onlyF(), fields, pb.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, tagB)),
                                         runOwn(A, xa, onlyF(), fieldsB, pa.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, tagHA)),
                                         runOwn(B, xb, onlyF(), fieldsB, pb.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, tagHB))};
        };
        const std::array<OwnRun, 4> rs = runSet(9600 + ci, 9610 + ci, 9620 + ci, 9630 + ci);
        const OwnRun &recA = rs[0], &recB = rs[1], &hbA = rs[2], &hbB = rs[3];
        out("A recipe run: %s\nB recipe run: %s\nA Level-B H0 run: %s\nB Level-B H0 run: %s\n", recA.stats.summary().c_str(),
            recB.stats.summary().c_str(), hbA.stats.summary().c_str(), hbB.stats.summary().c_str());
        const Gram gA = gramCheck3(recA, hbA, pa.w, A), gB = gramCheck3(recB, hbB, pb.w, B);
        const int Bq = gA.B;
        auto Sof = [](const Eigen::MatrixXd& G) { return Eigen::Matrix3d(G.block(2, 2, 3, 3)); };
        const Eigen::Matrix3d SA = Sof(gA.G), SB = Sof(gB.G);
        const Eigen::Matrix3d SAB = 0.5*(SA + SA.transpose()) + 0.5*(SB + SB.transpose());
        std::vector<Eigen::MatrixXd> jSA, jSB, jSAB;
        for (int a = 0; a < Bq; a++) {
            jSA.push_back(Sof(gA.jack[a]));
            jSB.push_back(Sof(gB.jack[a]));
            const Eigen::Matrix3d ja = Sof(gA.jack[a]), jb = Sof(gB.jack[a]);
            jSAB.push_back(0.5*(ja + ja.transpose()) + 0.5*(jb + jb.transpose()));
        }
        Eigen::VectorXd evA, evASe, evB, evBSe, evAB, evABSe;
        eigenWithSe(gA, SA, jSA, evA, evASe);
        eigenWithSe(gB, SB, jSB, evB, evBSe);
        eigenWithSe(gA, SAB, jSAB, evAB, evABSe);
        auto evStr = [](const Eigen::VectorXd& e, const Eigen::VectorXd& s) {
            std::string r;
            for (int k = 0; k < 3; k++) r += (k ? ", " : "") + fmtd("%.4f", e[k]) + " (" + fmtd("%.4f", s[k]) + ")";
            return r;
        };
        std::string sBlock = "y_f = " + cfg + "\n";
        for (const auto& [name, S, gm] : {std::tuple<const char*, const Eigen::Matrix3d*, const Gram*>{"S_A", &SA, &gA}, {"S_B", &SB, &gB}}) {
            out("%s (split-sample; rows/cols v0, v1, v2; jackknife s.e.; naive - split):\n", name);
            for (int i = 0; i < 3; i++) {
                const std::string row = matRow(*S, i) + "   s.e. " + matRow(gm->se.block(2, 2, 3, 3), i, "%.5f");
                out("  %s   naive-split %s\n", row.c_str(), matRow(gm->naive.block(2, 2, 3, 3) - *S, i, "%+.1e").c_str());
                sBlock += std::string(name) + "  " + row + "\n";
            }
            for (int h : {0, 1, kH0B, kH0B + 1}) {
                const std::string bs = fmtVec(Eigen::Vector3d(gm->G(h, 2), gm->G(h, 3), gm->G(h, 4))) + " s.e. " +
                                       fmtVec(Eigen::Vector3d(gm->se(h, 2), gm->se(h, 3), gm->se(h, 4)), "%.5f");
                const char* hn = h == 0 ? "Level A, g = 1" : h == 1 ? "Level A, ramp" : h == kH0B ? "Level B, g = 1" : "Level B, ramp";
                out("  b (%s) = %s;  int |grad H0|^2 = %.5f +- %.5f\n", hn, bs.c_str(), gm->G(h, h), gm->se(h, h));
                sBlock += std::string("  b_") + (name + 2) + " (" + hn + ") = " + bs + ";  int abs(grad H0)^2 = " + fmtd("%.5f", gm->G(h, h)) +
                          " +- " + fmtd("%.5f", gm->se(h, h)) + "\n";
            }
        }
        out("eigenvalues: S_A %s; S_B %s; S_A + S_B %s\n", evStr(evA, evASe).c_str(), evStr(evB, evBSe).c_str(),
            evStr(evAB, evABSe).c_str());
        // mirror check: B's xi problem is A's with xi_x <-> xi_y swapped (the geometry printout checks the crack), and the
        // swap maps the reference tet and F onto themselves, so the walks are the same in distribution: on the same points
        // S_B = S_A, b_B = b_A and the H0 energies agree in expectation. B's own point set, though, is built in B's frame
        // (Bey's refinement is not swap-symmetric), so on its own points the comparison also contains the quadrature error,
        // which the jackknife (Monte Carlo noise only) does not. Diagnostic: B's two runs once more on A's points mirrored
        // into B (xi_B = (xi_y, xi_x, xi_z) of A's: X = A's X with z -> -z, the same weights; independent seeds). Against A
        // that is the mirror check without the quadrature; against B's own points it is the quadrature error of S and b.
        PointSet pbm;
        for (size_t i = 0; i < pa.size(); i++) {
            const Vec3 xi(pa.xi[i].y(), pa.xi[i].x(), pa.xi[i].z());
            pbm.xi.push_back(xi);
            pbm.X.push_back(B.toX(xi));
            pbm.w.push_back(pa.w[i]);
        }
        const OwnRun recBm = runOwn(B, xb, onlyF(), fields, pbm.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 9640 + ci));
        const OwnRun hbBm = runOwn(B, xb, onlyF(), fieldsB, pbm.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 9650 + ci));
        const Gram gBm = gramCheck3(recBm, hbBm, pbm.w, B);
        out("B on A's points mirrored into B: recipe run %s\n  Level-B H0 run %s\n", recBm.stats.summary().c_str(),
            hbBm.stats.summary().c_str());
        // max |G1 - G2|/s.e. over S's 6 entries, b's 12 and the 4 H0 energies; rel: the largest |dS_ij|/sqrt(S_ii S_jj)
        auto compareGrams = [](const Gram& g1, const Gram& g2, double& zS, double& zb, double& zE, double& rel) {
            zS = zb = zE = rel = 0;
            for (int i = 2; i < 5; i++)
                for (int j = i; j < 5; j++) {
                    zS = std::max(zS, std::fabs(g1.G(i, j) - g2.G(i, j))/std::hypot(g1.se(i, j), g2.se(i, j)));
                    rel = std::max(rel, std::fabs(g1.G(i, j) - g2.G(i, j))/std::sqrt(g1.G(i, i)*g1.G(j, j)));
                }
            for (int h : {0, 1, kH0B, kH0B + 1}) {
                for (int k = 2; k < 5; k++) zb = std::max(zb, std::fabs(g1.G(h, k) - g2.G(h, k))/std::hypot(g1.se(h, k), g2.se(h, k)));
                zE = std::max(zE, std::fabs(g1.G(h, h) - g2.G(h, h))/std::hypot(g1.se(h, h), g2.se(h, h)));
            }
        };
        double zS, zb, zE, rS, zSm, zbm, zEm, rSm, zSq, zbq, zEq, rSq;
        compareGrams(gB, gA, zS, zb, zE, rS);
        compareGrams(gBm, gA, zSm, zbm, zEm, rSm);
        compareGrams(gB, gBm, zSq, zbq, zEq, rSq);
        out("mirror check, max |B - A|/s.e. over S's 6 entries / the 12 entries of b (Level A and B, both g types) / the 4 H0\n"
            "energies int |grad H0|^2 (independent runs):\n"
            "  B's own points vs A's:                   %.2f / %.2f / %.2f (largest |dS_ij|/sqrt(S_ii S_jj) = %.1f%%)\n"
            "  B on A's mirrored points vs A (common):  %.2f / %.2f / %.2f (%.1f%%)  <- the mirror check proper\n"
            "  B's own points vs A's mirrored points:   %.2f / %.2f / %.2f (%.1f%%)  <- the quadrature error, beyond the noise\n",
            zS, zb, zE, 100*rS, zSm, zbm, zEm, 100*rSm, zSq, zbq, zEq, 100*rSq);
        evMinLo = std::min(evMinLo, evAB[0]);
        evMinHi = std::max(evMinHi, evAB[0]);
        evMinSe = std::max(evMinSe, evABSe[0]);
        quadRel = std::max(quadRel, rSq);
        quadSig = std::max(quadSig, zSq);
        tS.rows.push_back({cfg, evStr(evA, evASe), evStr(evB, evBSe), evStr(evAB, evABSe),
                           fmtd("%.2f", zSm) + " / " + fmtd("%.2f", zbm) + " / " + fmtd("%.2f", zEm) + " (" + fmtd("%.1f", 100*rSm) + "%)",
                           fmtd("%.2f", zS) + " / " + fmtd("%.2f", zb) + " / " + fmtd("%.2f", zE) + " (" + fmtd("%.1f", 100*rS) + "%)",
                           fmtd("%.2f", zSq) + " / " + fmtd("%.2f", zbq) + " / " + fmtd("%.2f", zEq) + " (" + fmtd("%.1f", 100*rSq) + "%)"});
        sBlock += "S_B on A's mirrored points:\n";
        for (int i = 0; i < 3; i++)
            sBlock += "S_Bm  " + matRow(Sof(gBm.G), i) + "   s.e. " + matRow(gBm.se.block(2, 2, 3, 3), i, "%.5f") + "\n";
        sMd += sBlock;

        // replicate run sets with independent seeds (set 0 is the main one): the spread of the glued c over the sets
        // against one set's jackknife s.e., and the mean over the sets against the fitted trace
        const int nRep = gp.quick ? 3 : 8;
        std::vector<std::array<std::array<Glue3, 2>, 2>> reps(static_cast<size_t>(nRep));   // [rep][level][g type]
        double repSec = 0;
        for (int rep = 0; rep < nRep; rep++) {
            if (rep == 0) {
                for (int L = 0; L < 2; L++)
                    for (int m = 0; m < 2; m++) reps[0][L][m] = glueRun3(gA, gB, L ? kH0B + m : m);
                continue;
            }
            const uint64_t t0 = 9700 + 40*ci + 4*rep;
            const std::array<OwnRun, 4> rr = runSet(t0, t0 + 1, t0 + 2, t0 + 3);
            for (const OwnRun& r : rr) repSec += r.stats.seconds;
            const Gram ga = gramCheck3(rr[0], rr[2], pa.w, A), gb = gramCheck3(rr[1], rr[3], pb.w, B);
            for (int L = 0; L < 2; L++)
                for (int m = 0; m < 2; m++) reps[size_t(rep)][L][m] = glueRun3(ga, gb, L ? kH0B + m : m);
        }
        out("replicate run sets: %d (4 runs each; seed tags 9600/9610/9620/9630 + %d, then 9700 + %d + 4 rep + 0..3; %.1f s for the\n"
            "extra runs)\n", nRep, ci, 40*ci, repSec);

        // the reference's energies (split-sample over its batches): A, B, A u B, and the points near the crack's segment
        std::vector<double> wA(ref.X.size(), 0.0), wB(ref.X.size(), 0.0), wAB(ref.X.size(), 0.0), wS(ref.X.size(), 0.0);
        std::vector<double> wSA(pa.size(), 0.0), wSB(pb.size(), 0.0);
        for (size_t i = 0; i < ref.nA; i++) {
            wA[i] = wAB[i] = pa.w[i];
            if (nearSegment(pa.X[i])) wS[i] = wSA[i] = pa.w[i];
        }
        for (size_t i = 0; i < ref.nB; i++) {
            wB[ref.nA + i] = wAB[ref.nA + i] = pb.w[i];
            if (nearSegment(pb.X[i])) wS[ref.nA + i] = wSB[i] = pb.w[i];
        }
        const std::vector<double> unit2[2] = {unitCoef(2, 0), unitCoef(2, 1)};
        auto refGram = [&](const std::vector<double>& w) {
            return gramFromBatches(w, ref.sums.nBatches, 2, [&](size_t p, int b, Vec3* g) {
                for (int q = 0; q < 2; q++) g[q] = comboBatchGrad(ref.sums, p, b, unit2[q]);
            });
        };
        const Gram grA = refGram(wA), grB = refGram(wB), grAB = refGram(wAB), grS = refGram(wS);
        const Gram gAS = gramCheck3(recA, hbA, wSA, A), gBS = gramCheck3(recB, hbB, wSB, B);   // decomposed, near the segment
        // eps-sensitivity of the energies: the same integrals over the points farther than 2 eps from the crack (the rule
        // restricted to them integrates over T minus the 2 eps-tube), so the difference is the slab between eps and 2 eps.
        // Where the face data and the crack data disagree along the edge where the crack meets F, the field near that edge
        // is the 2D solution in a right-angle wedge (crack face and F) on each side of the crack: |grad H| ~ (2/pi) jump/r.
        // Over {distance to the crack > eps} its energy is (1/pi) jump^2 ln(R/eps) per unit length and wedge, log-divergent
        // as eps -> 0, and the slab between eps and 2 eps holds (ln 2/pi) jump^2 of it. The edge is shared by both cut
        // tets, each with one wedge on either side of the crack (A above F, B below, with the same trace on F), and slabD
        // sums A and B, so the edge model counts the four wedges: 2 (ln 2/pi) int (jump_+^2 + jump_-^2) dy, jump_+- = +-g -
        // J(d +- 0, y) - sum_k c_k hat_k(d, y). The reference has no such edge (F is interior to U).
        std::vector<double> w2A(pa.size(), 0.0), w2B(pb.size(), 0.0), w2(ref.X.size(), 0.0);
        for (size_t i = 0; i < pa.size(); i++)
            if (xa.cut.distToCut(pa.xi[i]) >= 2.0*gp.eps) w2[i] = w2A[i] = pa.w[i];
        for (size_t i = 0; i < pb.size(); i++)
            if (xb.cut.distToCut(pb.xi[i]) >= 2.0*gp.eps) w2[ref.nA + i] = w2B[i] = pb.w[i];
        const Gram gA2 = gramCheck3(recA, hbA, w2A, A), gB2 = gramCheck3(recB, hbB, w2B, B), gr2 = refGram(w2);
        auto edgeModel = [&](int m, const FaceFn* Jm, const Eigen::Vector3d& c) {
            const double xp = std::nextafter(cs.d, 1.0), xm = std::nextafter(cs.d, 0.0);
            const int nY = 400;
            double s = 0;
            for (int k = 0; k < nY; k++) {
                const double y = cs.yFront*(k + 0.5)/nY;
                const Vec3 X(cs.d, y, 0.0);
                const double gy = cu.g(X, gts[m]), tau = c[0]*hatF(0, X) + c[1]*hatF(1, X) + c[2]*hatF(2, X);
                const double jp = gy - (Jm ? (*Jm)(Vec3(xp, y, 0.0)) : 0.0) - tau, jm = -gy - (Jm ? (*Jm)(Vec3(xm, y, 0.0)) : 0.0) - tau;
                s += (jp*jp + jm*jm)*cs.yFront/nY;
            }
            const int nCutTets = 2;   // A and B: two wedges each
            return nCutTets*std::log(2.0)/std::numbers::pi*s;
        };
        const auto decField = decField3;
        struct Regions {
            ErrStats A, B, AB, nearF, nearSeg;
        };
        auto regions = [&](const FieldAt& dA, const FieldAt& dB, const FieldAt& rA, const FieldAt& rB) {
            Regions e;
            for (size_t i = 0; i < pa.size(); i++) {
                e.A.add(pa.w[i], pa.X[i], dA, i, rA, i);
                e.AB.add(pa.w[i], pa.X[i], dA, i, rA, i);
                if (std::fabs(pa.X[i].z()) < kNearF) e.nearF.add(pa.w[i], pa.X[i], dA, i, rA, i);
                if (nearSegment(pa.X[i])) e.nearSeg.add(pa.w[i], pa.X[i], dA, i, rA, i);
            }
            for (size_t i = 0; i < pb.size(); i++) {
                e.B.add(pb.w[i], pb.X[i], dB, i, rB, i);
                e.AB.add(pb.w[i], pb.X[i], dB, i, rB, i);
                if (std::fabs(pb.X[i].z()) < kNearF) e.nearF.add(pb.w[i], pb.X[i], dB, i, rB, i);
                if (nearSegment(pb.X[i])) e.nearSeg.add(pb.w[i], pb.X[i], dB, i, rB, i);
            }
            return e;
        };
        auto regionList = [](const Regions& e) {
            return std::vector<std::pair<const char*, const ErrStats*>>{{"A", &e.A}, {"B", &e.B}, {"A u B", &e.AB},
                                                                        {"near F (A, B)", &e.nearF}, {"near seg (A, B)", &e.nearSeg}};
        };
        auto rmsMax = [](double rms, double nc, double mx) {
            return fmtd("%.2e", rms) + " (" + fmtd("%.2e", nc) + ") / " + fmtd("%.2e", mx);
        };

        const std::string wallRef = ref.cached ? fmtd("%.1f", ref.stats.seconds) + " (cached run)" : fmtd("%.1f", ref.stats.seconds);
        for (int m = 0; m < 2; m++) {
            const FieldAt refA = refField(ref, m, 0, ref.nA), refB = refField(ref, m, ref.nA, ref.nB),
                          refF = refField(ref, m, ref.nA + ref.nB, ref.nF);
            ErrStats rf;   // the reference's own noise on F
            for (size_t p = 0; p < ref.nF; p++) rf.add(fps.w[p], fps.X[p], 0.0, Vec3::Zero(), refF.vH[p], 0.0, Vec3::Zero(), Vec3::Zero());
            for (int L = 0; L < 2; L++) {
                const int h = L ? kH0B + m : m;
                const FaceFn* Jm = L ? &J[m] : nullptr;
                out("\ny_f = %s, Level %s, %s:\n", cfg.c_str(), L ? "B (J)" : "A (J = 0)", gTypeName(gts[m]));
                Eigen::Matrix3d M;
                Eigen::Vector3d rhs;
                const Glue3 gl = glueRun3(gA, gB, h, &M, &rhs);
                const Eigen::Vector3d& c = gl.c;
                out("  b_A = %s, b_B = %s\n", fmtVec(Eigen::Vector3d(gA.G(h, 2), gA.G(h, 3), gA.G(h, 4))).c_str(),
                    fmtVec(Eigen::Vector3d(gB.G(h, 2), gB.G(h, 3), gB.G(h, 4))).c_str());
                out("  3 x 3 glue system (S_A + S_B) [c_v0, c_v1, c_v2] = -(b_A + b_B):\n");
                for (int i = 0; i < 3; i++) out("    [%s]  [%+.5f]\n", matRow(M, i).c_str(), rhs[i]);
                out("  glued: c = (%+.5f, %+.5f, %+.5f) +- (%.5f, %.5f, %.5f) (jackknife over the batches of all four runs, the\n"
                    "  system re-solved per replicate)\n", c[0], c[1], c[2], gl.se[0], gl.se[1], gl.se[2]);

                // fields and errors per region
                const FieldAt decA = decField(A, recA, hbA, L, m, c, &gl.jack), decB = decField(B, recB, hbB, L, m, c, &gl.jack);
                const Regions e = regions(decA, decB, refA, refB);
                out("  errors of the glued decomposition vs the reference, / 2 g_max (nc = noise-corrected RMS):\n"
                    "    region           points   H RMS (nc)          H max     grad RMS (nc)       grad max  | ref s.e. H RMS/max  dec s.e. H RMS/max  ref s.e. grad RMS/max  dec s.e. grad RMS/max\n");
                for (const auto& [name, s] : regionList(e)) {
                    out("    %-16s %6lld   %.2e (%.2e)   %.2e  %.2e (%.2e)   %.2e  | %.1e / %.1e     %.1e / %.1e     %.1e / %.1e       %.1e / %.1e\n",
                        name, s->n, s->rmsH(), s->corrH(), s->maxH/kTwoGMax, s->rmsG(), s->corrG(), s->maxG/kTwoGMax, s->seRefH(),
                        s->maxSeRefH/kTwoGMax, s->seFldH(), s->maxSeFldH/kTwoGMax, s->seRefG(), s->maxSeRefG/kTwoGMax, s->seFldG(),
                        s->maxSeFldG/kTwoGMax);
                    tReg.rows.push_back({cfg, lvShort[L], gShort[m], name, std::to_string(s->n),
                                         fmtd("%.2e", s->rmsH()) + " (" + fmtd("%.2e", s->corrH()) + ")", fmtd("%.2e", s->maxH/kTwoGMax),
                                         fmtd("%.2e", s->rmsG()) + " (" + fmtd("%.2e", s->corrG()) + ")", fmtd("%.2e", s->maxG/kTwoGMax),
                                         fmtd("%.1e", s->seRefH()) + " / " + fmtd("%.1e", s->maxSeRefH/kTwoGMax),
                                         fmtd("%.1e", s->seFldH()) + " / " + fmtd("%.1e", s->maxSeFldH/kTwoGMax),
                                         fmtd("%.1e", s->seRefG()) + " / " + fmtd("%.1e", s->maxSeRefG/kTwoGMax),
                                         fmtd("%.1e", s->seFldG()) + " / " + fmtd("%.1e", s->maxSeFldG/kTwoGMax)});
                }
                out("    largest H error at X = %s, largest grad H error at X = %s\n", fmt3(e.AB.worstH, 3).c_str(),
                    fmt3(e.AB.worstG, 3).c_str());

                // energies: each tet from its Gram with the glued c; the reference from its own
                const double EdA = tetEnergy(gA.G, h, c), EdB = tetEnergy(gB.G, h, c);
                const double EdS = tetEnergy(gAS.G, h, c) + tetEnergy(gBS.G, h, c);
                std::vector<double> rA, rB, rT, rS;
                for (int a = 0; a < Bq; a++) {
                    rA.push_back(tetEnergy(gA.jack[a], h, gl.jack[a]));
                    rB.push_back(tetEnergy(gB.jack[a], h, gl.jack[a]));
                    rT.push_back(rA.back() + rB.back());
                    rS.push_back(tetEnergy(gAS.jack[a], h, gl.jack[a]) + tetEnergy(gBS.jack[a], h, gl.jack[a]));
                }
                const double seDA = jackSe(rA), seDB = jackSe(rB), seD = jackSe(rT), seDS = jackSe(rS);
                const double ErA = 0.5*grA.G(m, m), ErB = 0.5*grB.G(m, m), Er = 0.5*grAB.G(m, m), ErS = 0.5*grS.G(m, m);
                const double seRA = 0.5*grA.se(m, m), seRB = 0.5*grB.se(m, m), seR = 0.5*grAB.se(m, m), seRS = 0.5*grS.se(m, m);
                const double dE = EdA + EdB - Er, seDE = std::hypot(seD, seR);
                const double halfE2 = 0.5*(e.AB.eG2 - e.AB.vRefG - e.AB.vFldG);
                const double naiveDec = tetEnergy(gA.naive, h, c) + tetEnergy(gB.naive, h, c) - (EdA + EdB);
                const double naiveRef = 0.5*(grAB.naive(m, m) - grAB.G(m, m));
                // the slab between eps and 2 eps from the crack (decomposed, reference, their difference) and the edge model
                const double slabD = EdA + EdB - tetEnergy(gA2.G, h, c) - tetEnergy(gB2.G, h, c), slabR = Er - 0.5*gr2.G(m, m);
                std::vector<double> rSD, rSR;
                for (int a = 0; a < Bq; a++) rSD.push_back(rT[a] - tetEnergy(gA2.jack[a], h, gl.jack[a]) - tetEnergy(gB2.jack[a], h, gl.jack[a]));
                for (int a = 0; a < gr2.B; a++) rSR.push_back(0.5*(grAB.jack[a](m, m) - gr2.jack[a](m, m)));
                const double seSD = jackSe(rSD), seSR = jackSe(rSR), edge = edgeModel(m, Jm, c);

                // the trace-fit diagnostic: the reference trace on F fitted to the 3 hats (Level B: with J at coefficient 1),
                // both tets run with that trace (H0 + sum_k c_fit,k P_k, by superposition) and compared with the reference
                const TraceFit tf = fitTrace(fps, refF, Jm);
                std::vector<Eigen::Vector4d> cJ4;
                for (const Eigen::Vector3d& x : gl.jack) cJ4.push_back(Eigen::Vector4d(x[0], x[1], x[2], 0.0));
                const TraceErr fitE = traceError(fps, refF, tf.c, nullptr, nullptr, Jm, &nearF);
                const TraceErr glE = traceError(fps, refF, c, &cJ4, nullptr, Jm, &nearF);
                FieldAt fitA = decField(A, recA, hbA, L, m, tf.c, nullptr), fitB = decField(B, recB, hbB, L, m, tf.c, nullptr);
                addCoefNoise(fitA, recA, A, {2, 3, 4}, tf.cov);
                addCoefNoise(fitB, recB, B, {2, 3, 4}, tf.cov);
                const Regions f = regions(fitA, fitB, refA, refB);
                const double Efit = tetEnergy(gA.G, h, tf.c) + tetEnergy(gB.G, h, tf.c);
                Eigen::Vector3d cfSe;
                for (int k = 0; k < 3; k++) cfSe[k] = std::sqrt(tf.cov(k, k));
                out("  energy E = (1/2) int |grad H|^2 (split-sample; jackknife s.e.):\n"
                    "    A minus tube: decomposed %.5f +- %.5f, reference %.5f +- %.5f\n"
                    "    B minus tube: decomposed %.5f +- %.5f, reference %.5f +- %.5f\n"
                    "    A u B:        decomposed %.5f +- %.5f, reference %.5f +- %.5f; difference %+.5f +- %.5f\n"
                    "    near seg:     decomposed %.5f +- %.5f, reference %.5f +- %.5f (points within %.2f of the crack's segment on F)\n"
                    "    (1/2) int |grad(H_dec - H_ref)|^2 (noise-corrected) = %.5f (equal to the difference for any c up to the eps-tube\n"
                    "      term: an estimator check, not a test of the glue); with the fitted trace E = %.5f\n"
                    "    naive - split: decomposed %+.1e, reference %+.1e\n"
                    "    eps-sensitivity, the slab between eps and 2 eps from the crack: decomposed %.5f +- %.5f, reference %.5f +- %.5f,\n"
                    "      difference %+.5f +- %.5f; edge model 2 (ln 2/pi) int (jump_+^2 + jump_-^2) dy (four wedges: two per tet) = %.5f\n",
                    EdA, seDA, ErA, seRA, EdB, seDB, ErB, seRB, EdA + EdB, seD, Er, seR, dE, seDE, EdS, seDS, ErS, seRS, kNearSeg, halfE2,
                    Efit, naiveDec, naiveRef, slabD, seSD, slabR, seSR, slabD - slabR, std::hypot(seSD, seSR), edge);
                out("  trace-fit diagnostic: least-squares fit of the reference trace on F (%zu points, area-weighted) to the 3 hats%s:\n"
                    "    c_fit = (%+.5f, %+.5f, %+.5f) +- (%.5f, %.5f, %.5f); glued c = (%+.5f, %+.5f, %+.5f)\n"
                    "    residual of the fit on F, / 2 g_max: RMS %.2e (nc %.2e), max %.2e; near seg (%.0f%% of the area) RMS %.2e (nc %.2e),\n"
                    "      max %.2e; far RMS %.2e (nc %.2e); near seg share of int e^2 %.0f%% (nc %.0f%%)\n"
                    "    the glued trace vs the reference on F: RMS %.2e (nc %.2e, the noise of c included), max %.2e; near seg RMS %.2e\n"
                    "      (nc %.2e); far RMS %.2e (nc %.2e); near seg share of int e^2 %.0f%% (nc %.0f%%)\n"
                    "    the reference's s.e. on F: RMS %.1e, max %.1e\n",
                    fps.X.size(), L ? " with J at coefficient 1" : "", tf.c[0], tf.c[1], tf.c[2], cfSe[0], cfSe[1], cfSe[2], c[0], c[1],
                    c[2], fitE.rms, fitE.nc, fitE.max, 100*fitE.areaNear, fitE.rmsNear, fitE.ncNear, fitE.maxNear, fitE.rmsFar,
                    fitE.ncFar, 100*fitE.shareNear, 100*fitE.shareNearNc, glE.rms, glE.nc, glE.max, glE.rmsNear, glE.ncNear,
                    glE.rmsFar, glE.ncFar, 100*glE.shareNear, 100*glE.shareNearNc, rf.seRefH(), rf.maxSeRefH/kTwoGMax);
                out("    both tets with the fitted trace vs the reference (glued in brackets), / 2 g_max:\n");
                const auto fl = regionList(f), gl2 = regionList(e);
                for (size_t r = 0; r < fl.size(); r++) {
                    const ErrStats &fs = *fl[r].second, &gs = *gl2[r].second;
                    out("      %-16s H RMS %.2e (nc %.2e), max %.2e [%.2e (nc %.2e), %.2e]; grad RMS %.2e (nc %.2e) [%.2e (nc %.2e)]\n",
                        fl[r].first, fs.rmsH(), fs.corrH(), fs.maxH/kTwoGMax, gs.rmsH(), gs.corrH(), gs.maxH/kTwoGMax, fs.rmsG(),
                        fs.corrG(), gs.rmsG(), gs.corrG());
                    tFit.rows.push_back({cfg, lvShort[L], gShort[m], fl[r].first, rmsMax(gs.rmsH(), gs.corrH(), gs.maxH/kTwoGMax),
                                         rmsMax(fs.rmsH(), fs.corrH(), fs.maxH/kTwoGMax), rmsMax(gs.rmsG(), gs.corrG(), gs.maxG/kTwoGMax),
                                         rmsMax(fs.rmsG(), fs.corrG(), fs.maxG/kTwoGMax),
                                         fmtd("%.1e", fs.seFldH()) + " / " + fmtd("%.1e", fs.maxSeFldH/kTwoGMax)});
                }

                // the replicate sets: sd of the glued c over the sets / rms jackknife s.e., and their mean against c_fit
                Eigen::Vector3d rMean = Eigen::Vector3d::Zero(), rSd = Eigen::Vector3d::Zero(), rRms = Eigen::Vector3d::Zero();
                Eigen::Matrix3d rCov = Eigen::Matrix3d::Zero();
                for (const auto& r : reps) {
                    rMean += r[L][m].c/double(nRep);
                    rRms += r[L][m].se.cwiseAbs2()/double(nRep);
                }
                for (const auto& r : reps) {
                    rSd += (r[L][m].c - rMean).cwiseAbs2()/double(nRep - 1);
                    rCov += (r[L][m].c - rMean)*(r[L][m].c - rMean).transpose()/double(nRep - 1);
                }
                rSd = rSd.cwiseSqrt();
                rRms = rRms.cwiseSqrt();
                const Eigen::Matrix3d meanCov = rCov/double(nRep);
                const Eigen::Vector3d dc = rMean - tf.c;
                Eigen::Vector3d dcSe, mSe;
                for (int k = 0; k < 3; k++) {
                    mSe[k] = rSd[k]/std::sqrt(double(nRep));
                    dcSe[k] = std::sqrt(meanCov(k, k) + tf.cov(k, k));
                }
                const TraceErr mE = traceError(fps, refF, rMean, nullptr, &meanCov, Jm, &nearF);
                const Eigen::Vector3d ratio = rSd.cwiseQuotient(rRms);
                // the glue's sensitivity to B's quadrature: c against c glued with B's Gram on A's mirrored points (the same
                // A runs; the jackknife deletes batch a of every run, so the s.e. is that of the difference)
                const Glue3 glm = glueRun3(gA, gBm, h);
                const Eigen::Vector3d dq = c - glm.c;
                Eigen::Vector3d dqSe;
                for (int k = 0; k < 3; k++) {
                    std::vector<double> rep;
                    for (int a = 0; a < Bq; a++) rep.push_back(gl.jack[a][k] - glm.jack[a][k]);
                    dqSe[k] = jackSe(rep);
                }
                double dqT = 0, dqW = 0;
                for (size_t p = 0; p < fps.X.size(); p++) {
                    const double v = dq[0]*hatF(0, fps.X[p]) + dq[1]*hatF(1, fps.X[p]) + dq[2]*hatF(2, fps.X[p]);
                    dqT += fps.w[p]*v*v;
                    dqW += fps.w[p];
                }
                dqT = std::sqrt(dqT/dqW)/kTwoGMax;
                out("  quadrature sensitivity: c - c with B on A's mirrored points = %s +- %s; that difference's trace on F: RMS\n"
                    "    %.2e (/ 2 g_max)\n", fmtVec(dq).c_str(), fmtVec(dqSe, "%.5f").c_str(), dqT);
                out("  replicates (%d run sets): sd of (c_v0, c_v1, c_v2) over the sets / rms jackknife s.e. = %.2f %.2f %.2f (1 expected;\n"
                    "    an sd from %d sets scatters by ~%.0f%%)\n"
                    "    mean glued c over the sets = %s +- %s\n"
                    "    mean glued c - c_fit       = %s +- %s\n"
                    "    the mean c's trace vs the reference on F, / 2 g_max: RMS %.2e (nc %.2e), max %.2e; near seg share (nc) %.0f%%\n",
                    nRep, ratio[0], ratio[1], ratio[2], nRep, 100.0/std::sqrt(2.0*(nRep - 1)), fmtVec(rMean).c_str(),
                    fmtVec(mSe, "%.5f").c_str(), fmtVec(dc).c_str(), fmtVec(dcSe, "%.5f").c_str(), mE.rms, mE.nc, mE.max,
                    100*mE.shareNearNc);

                // tables
                const std::string cStr = fmtVec(c, "%+.4f") + " (" + fmtd("%.4f", gl.se[0]) + ", " + fmtd("%.4f", gl.se[1]) + ", " +
                                         fmtd("%.4f", gl.se[2]) + ")";
                tMain.rows.push_back({cfg, lvShort[L], gShort[m], fmtd("%.2e", e.AB.rmsH()) + " / " + fmtd("%.2e", e.AB.maxH/kTwoGMax),
                                      fmtd("%.2e", e.AB.rmsG()) + " / " + fmtd("%.2e", e.AB.maxG/kTwoGMax),
                                      fmtd("%.1e", e.AB.seRefH()) + " / " + fmtd("%.1e", e.AB.maxSeRefH/kTwoGMax) + "; " +
                                          fmtd("%.1e", e.AB.seRefG()) + " / " + fmtd("%.1e", e.AB.maxSeRefG/kTwoGMax),
                                      fmtd("%.1e", e.AB.seFldH()) + " / " + fmtd("%.1e", e.AB.maxSeFldH/kTwoGMax) + "; " +
                                          fmtd("%.1e", e.AB.seFldG()) + " / " + fmtd("%.1e", e.AB.maxSeFldG/kTwoGMax),
                                      fmtd("%.5f", EdA + EdB) + " / " + fmtd("%.5f", Er) + " / " + fmtd("%+.5f", dE) + " +- " + fmtd("%.5f", seDE),
                                      evStr(evAB, evABSe), cStr,
                                      "A " + fmtd("%.1f", recA.stats.seconds) + ", B " + fmtd("%.1f", recB.stats.seconds) + " (recipes); A " +
                                          fmtd("%.1f", hbA.stats.seconds) + ", B " + fmtd("%.1f", hbB.stats.seconds) + " (Level-B H0); reference " +
                                          wallRef + "; " + zombieWall(gp, zr)});
                tTr.rows.push_back({cfg, lvShort[L], gShort[m], cStr,
                                    fmtVec(tf.c, "%+.4f") + " (" + fmtd("%.4f", cfSe[0]) + ", " + fmtd("%.4f", cfSe[1]) + ", " + fmtd("%.4f", cfSe[2]) + ")",
                                    rmsMax(fitE.rms, fitE.nc, fitE.max),
                                    fmtd("%.2e", fitE.rmsNear) + " (" + fmtd("%.2e", fitE.ncNear) + ") / " + fmtd("%.2e", fitE.rmsFar) + " (" +
                                        fmtd("%.2e", fitE.ncFar) + "); " + fmtd("%.0f", 100*fitE.shareNear) + "% (" + fmtd("%.0f", 100*fitE.shareNearNc) + "%)",
                                    rmsMax(glE.rms, glE.nc, glE.max),
                                    fmtd("%.2e", glE.rmsNear) + " (" + fmtd("%.2e", glE.ncNear) + ") / " + fmtd("%.2e", glE.rmsFar) + " (" +
                                        fmtd("%.2e", glE.ncFar) + "); " + fmtd("%.0f", 100*glE.shareNear) + "% (" + fmtd("%.0f", 100*glE.shareNearNc) + "%)",
                                    fmtd("%.1e", rf.seRefH()) + " / " + fmtd("%.1e", rf.maxSeRefH/kTwoGMax)});
                tRep.rows.push_back({cfg, lvShort[L], gShort[m], std::to_string(nRep),
                                     fmtd("%.2f", ratio[0]) + " / " + fmtd("%.2f", ratio[1]) + " / " + fmtd("%.2f", ratio[2]),
                                     fmtVec(rMean, "%+.4f") + " +- " + fmtVec(mSe, "%.4f"), fmtVec(tf.c, "%+.4f") + " +- " + fmtVec(cfSe, "%.4f"),
                                     fmtVec(dc, "%+.4f") + " +- " + fmtVec(dcSe, "%.4f"),
                                     rmsMax(mE.rms, mE.nc, mE.max) + "; " + fmtd("%.0f", 100*mE.shareNearNc) + "%",
                                     fmtVec(dq, "%+.4f") + " +- " + fmtVec(dqSe, "%.4f") + "; " + fmtd("%.2e", dqT)});
                tEn.rows.push_back({cfg, lvShort[L], gShort[m], fmtd("%.5f", EdA) + " +- " + fmtd("%.5f", seDA),
                                    fmtd("%.5f", ErA) + " +- " + fmtd("%.5f", seRA), fmtd("%.5f", EdB) + " +- " + fmtd("%.5f", seDB),
                                    fmtd("%.5f", ErB) + " +- " + fmtd("%.5f", seRB), fmtd("%.5f", EdA + EdB) + " +- " + fmtd("%.5f", seD),
                                    fmtd("%.5f", Er) + " +- " + fmtd("%.5f", seR), fmtd("%+.5f", dE) + " +- " + fmtd("%.5f", seDE),
                                    fmtd("%.5f", halfE2),
                                    fmtd("%.5f", EdS) + " +- " + fmtd("%.5f", seDS) + " / " + fmtd("%.5f", ErS) + " +- " + fmtd("%.5f", seRS),
                                    fmtd("%.5f", Efit), fmtd("%+.1e", naiveDec) + " / " + fmtd("%+.1e", naiveRef),
                                    fmtd("%+.5f", slabD - slabR) + " +- " + fmtd("%.5f", std::hypot(seSD, seSR)) + " (" + fmtd("%.5f", edge) + ")"});
                cases.push_back({cfg, L ? "B (J)" : "A (J = 0)", m ? "ramp" : "g = 1", e.AB.rmsH(), e.AB.corrH(), e.AB.maxH/kTwoGMax,
                                 std::hypot(e.AB.seRefH(), e.AB.seFldH()), e.nearSeg.rmsH(), EdA + EdB, Er, seDE, f.AB.corrH(),
                                 e.AB.corrH(), fitE.nc, mE.nc, fitE.shareNearNc, slabD - slabR, std::hypot(seSD, seSR), edge, dqT,
                                 e.AB.rmsG(), e.AB.corrG(), std::hypot(e.AB.seRefG(), e.AB.seFldG()),
                                 zr.H[m].ncRmsDiffBound(kTwoGMax), zBoundWorst(zr, m), ratio.minCoeff(), ratio.maxCoeff()});
            }
        }
        // the runs of this configuration
        auto runRow = [&](const std::string& run, const std::string& solver, int N, const SolveStats& st, const std::string& note) {
            tRun.rows.push_back({cfg, run, solver, std::to_string(N), std::to_string(st.nPoints), fmtd("%.2f", st.seconds),
                                 fmtd("%.1f", st.meanWalkLength()), std::to_string(st.walksCapped) + " / " + std::to_string(st.walksEscaped), note});
        };
        const std::string repNote = std::to_string(nRep - 1) + " replicate sets of the four runs: " + fmtd("%.1f", repSec) + " s in all";
        runRow("A recipe (H0 x2, P x3)", "own loop (double)", gp.nWalks, recA.stats, repNote);
        runRow("B recipe (H0 x2, P x3), B's xi frame", "own loop (double)", gp.nWalks, recB.stats, "");
        runRow("A Level-B H0 (x2)", "own loop (double)", gp.nWalks, hbA.stats, "plain first ball at F (J not harmonic)");
        runRow("B Level-B H0 (x2)", "own loop (double)", gp.nWalks, hbB.stats, "plain first ball at F (J not harmonic)");
        runRow("B recipe on A's points mirrored into B", "own loop (double)", gp.nWalks, recBm.stats, "mirror and quadrature diagnostic");
        runRow("B Level-B H0 on A's points mirrored into B", "own loop (double)", gp.nWalks, hbBm.stats, "mirror and quadrature diagnostic");
        runRow("reference on U (both g modes)", "own loop (double)", gp.nRef, ref.stats,
               ref.cached ? "loaded from gate_cache/ (wall clock of the original run)"
                          : "timing pilot N = " + std::to_string(ref.nPilot) + ": " + fmtd("%.1f", ref.pilotSeconds) + " s, extrapolated " +
                                fmtd("%.0f", ref.pilotExtrapolated) + " s");
        runRow("Zombie cross-check of the reference (H, grad H, both g modes)", "Zombie (float)", gp.nZombie, zr.st,
               "points, wall clock and walks summed over its " + std::to_string(2*(1 + gp.nZombieBatches)) + " solves (grad: " +
                   std::to_string(gp.nZombieBatches) + " x " + std::to_string(gp.nZombie/gp.nZombieBatches) + " walks)");
    }
    out("\nCHECK 3: %s (the spec sets no pass threshold for checks 2 and 3: the numbers above are the result; the verdict covers\n"
        "the geometry assertions and the Zombie cross-check of the reference)\n", ok ? "done" : "FAILED");
    md += "### Main table (A u B)\n\n" + tMain.str() +
          "\n### Errors by region\n\nRegions: A and B (each minus its eps-tube), A u B, near F = the points with abs(z) < 0.1, near seg = "
          "the points within 0.1 of the crack's segment on F.\n\n" + tReg.str() +
          "\n### Trace on F\n\nThe reference's trace on F (the fixed 1024 points minus those within eps of the crack's segment) least-"
          "squares fitted (area-weighted) to the 3 hats, Level B with J at coefficient 1: residual e = H_ref - J - sum c hat (J = 0 "
          "at Level A). The glued trace's error is the same with the glued c. \"near seg\" = the F points within 0.1 of the segment; "
          "its share of the integral of e^2 against its share of the area (printed on the console) tells how much of the trace error "
          "sits at the crack.\n\n" + tTr.str() +
          "\n### Fitted trace against the glued c, per region\n\nBoth tets run with the fitted trace (H0 + sum c_fit P_k, by "
          "superposition) and compared with the reference; its noise includes c_fit's. This separates \"the face trace (3 hats, or "
          "3 hats + J) is too poor\" (the fitted trace's error) from \"the glue is poor\" (the glued error beyond it).\n\n" + tFit.str() +
          "\n### Replicate run sets\n\nThe four runs of a configuration repeated with independent seeds (set 1 is the main one). The sd "
          "of the glued c over the sets against one set's jackknife s.e. checks the decomposed noise floor (below 1: conservative); "
          "the mean over the sets, against c fit, shows a systematic offset of the glue that one set's noise can hide. The last "
          "column is the glue's sensitivity to B's quadrature: the main c minus c glued with B's Gram on A's mirrored points (see "
          "S matrices), and that difference's trace on F (RMS / 2 g_max).\n\n" + tRep.str() +
          "\n### Energies\n\nE = (1/2) int abs(grad H)^2 over A and B minus their eps-tubes; near seg = over the points within 0.1 of the "
          "crack's segment on F. Last column: the energy in the slab between eps and 2 eps from the crack (the same rule restricted "
          "to the points farther than 2 eps, subtracted), decomposed minus reference, and in brackets the edge model: where the face "
          "data and the crack's +-g disagree along the edge where the crack meets F, the field there is the right-angle-wedge "
          "solution abs(grad H) ~ (2/pi) jump/r on each side of the crack, whose energy outside the eps-tube grows like (1/pi) "
          "jump^2 ln(1/eps) per unit length and wedge (log-divergent as eps -> 0). The edge has four such wedges, one on each "
          "side of the crack in A and in B (the slab sums both tets), and the slab holds 2 (ln 2/pi) int (jump_+^2 + jump_-^2) dy "
          "of their energy, jump_+- = +-g - J(0.3 +- 0, y) - sum c hat(0.3, y). This is the 2D leading-order term only. The "
          "reference has no such edge. The column \"(1/2) int abs(grad e)^2\" would equal E dec - E ref for any H_dec that "
          "carries the crack data and is continuous across F, whatever c is, up to a boundary term on the eps-tube's surface "
          "(see check 2): it checks the energy and error estimators, not the glue. Here its gap to E dec - E ref is that boundary "
          "term, int over the tube's surface of e dH_true/dn (plus noise and quadrature error), which need not be small: e = "
          "H_dec - H_true vanishes on the crack but not at distance eps from it, least of all along the crack-F edge, where "
          "H_dec is singular.\n\n" +
          tEn.str() +
          "\n### S matrices\n\nEigenvalues with jackknife s.e. The mirror check compares B with A: B's xi problem is A's with xi_x and "
          "xi_y swapped (checked exactly on the crack in the geometry printout), and the swap maps the reference tet and F onto "
          "themselves, so on the same points S_B = S_A, b_B = b_A and the H0 energies agree in expectation (independent runs). B's "
          "own point set is built in B's frame by the same rule and is not the mirror image of A's, so B's own points against A's "
          "also contain the quadrature error, which the jackknife (Monte Carlo noise only) does not. Hence B's two runs are "
          "repeated on A's points mirrored into B (the same weights): against A that is the mirror check proper (common points), "
          "against B's own points it isolates the quadrature error of S and b. \"rel. dS\" = the largest abs(dS_ij)/sqrt(S_ii "
          "S_jj).\n\n" + tS.str() + "\n" + codeBlock(sMd) + refZombieMd(gp, tZ, tZr) + "\n### Runs\n\n" + tRun.str() + "\n### Quadrature of the point sets\n\n" + tQ.str() + "\n" + std::string("**Check 3: ") +
          (ok ? "done" : "FAILED") + "** (no pass threshold in the spec; the verdict covers the geometry assertions and the reference's "
          "Zombie cross-check)\n";

    // findings: generated from the cases above
    std::string f = std::string("- Check 3 (both tets cut by the crack x = 0.3, its front crossing F at (0.3, y_f, 0)): ") +
                    (ok ? "done" : "FAILED") + ". The spec sets no pass threshold; the geometry assertions " +
                    (ok ? "hold" : "or the cross-check FAILED") + ", and Zombie agrees with the reference in " +
                    std::to_string(nZagree) + " of " + std::to_string(nZcases) + " cases.\n";
    f += "- The glued decomposition against the reference, over A u B, divided by 2 g_max = 2. near seg = the points within 0.1 "
         "of the crack's segment on F; noise = the combined RMS s.e. of H, sqrt(ref^2 + dec^2), where dec includes the noise of "
         "the glued c; error/noise ~ 1 means not resolved above the Monte Carlo noise of one run.\n\n";
    std::vector<std::vector<std::string>> rows = {{"y_f", "level", "g", "H RMS (nc)", "near seg: H RMS", "H max", "noise",
                                                   "error/noise", "E dec / E ref - 1 (s.e.)", "ref check: bound (worst region)"}};
    double rLo = kInf, rHi = 0, shLo = kInf, shHi = 0, vLo = kInf, vHi = 0, calLo = kInf, calHi = 0, ebLo = kInf;
    bool allResolved = true;
    for (const Case& c : cases) {
        rows.push_back({c.yf, c.lv, c.g, pct(c.rms) + " (" + pct(c.nc) + ")", pct(c.segRms), pct(c.max, "%.1f%%"), pct(c.noise),
                        fmtd("%.1f", c.rms/c.noise), fmtd("%+.1f%%", 100*(c.Ed/c.Er - 1)) + " (" + fmtd("%.1f%%", 100*c.seDE/c.Er) + ")",
                        pct(c.zAll) + " (" + pct(c.zWorst) + ")"});
        rLo = std::min(rLo, c.rms/c.noise);
        rHi = std::max(rHi, c.rms/c.noise);
        allResolved = allResolved && c.nc > 0;
        calLo = std::min(calLo, c.calLo);
        calHi = std::max(calHi, c.calHi);
        ebLo = std::min(ebLo, c.rms/std::max(c.zAll, c.zWorst));
    }
    f += fixedTable(rows) + "\n";
    // Level B against Level A, per g mode, at each y_f (the cases come as y_f, then g, then level)
    std::string ratio[2];
    for (size_t i = 0; i + 1 < cases.size(); i += 2) {
        const int m = cases[i].g == "ramp";
        ratio[m] += (ratio[m].empty() ? "" : ", ") + fmtd("%.2f", cases[i + 1].rms/cases[i].rms);
    }
    f += "- The error/noise ratio is " + fmtd("%.1f", rLo) + " to " + fmtd("%.1f", rHi) + (allResolved ? ": every case is resolved "
         "above the noise (nc > 0), so these are decomposition errors, not Monte Carlo noise" : "") + ". Level B over Level A, H RMS "
         "over A u B at y_f = 0.15, 0.30, 0.45: g = 1 " + ratio[0] + "; ramp " + ratio[1] + ".\n";
    f += "- What \"Zombie agrees\" can see: the last column is the largest RMS discrepancy between the reference and Zombie that "
         "their comparison does not exclude (2-s.e. bound, / 2 g_max, over all points; the worst of the regions A, B, F, near F "
         "and near the crack in brackets). The H error is " + fmtd("%.1f", ebLo) + " or more times the larger of the two. Noise "
         "calibration: over the replicate sets the sd of the glued c is " + fmtd("%.2f", calLo) + " to " + fmtd("%.2f", calHi) +
         " times its jackknife s.e. (1 expected; the sd of " + std::to_string(gp.quick ? 3 : 8) + " sets scatters by ~" +
         fmtd("%.0f", 100.0/std::sqrt(2.0*((gp.quick ? 3 : 8) - 1))) + "%).\n";
    f += "- Trace-fit diagnostic: both tets run with the best fit of the reference trace on F in the same space (3 hats; Level B: "
         "3 hats + J), against the glued runs. On F, the fit's residual against the error of the glued trace with c averaged over "
         "the replicate sets; \"trace space\" = the fit residual's share of that error in mean square, the rest comes from the "
         "glue; \"near seg\" = the share of the fit residual's mean square within 0.1 of the segment. All nc, / 2 g_max.\n\n";
    rows = {{"y_f", "level", "g", "A u B: H RMS, fitted / glued", "F: fit residual / glued mean c", "trace space", "near seg"}};
    for (const Case& c : cases) {
        const double share = c.meanF > 0 ? c.fitF*c.fitF/(c.meanF*c.meanF) : 0.0;
        rows.push_back({c.yf, c.lv, c.g, pct(c.fitAB) + " / " + pct(c.gluAB), pct(c.fitF) + " / " + pct(c.meanF),
                        c.meanF > 0 ? pct(share, "%.0f%%") : "-", pct(c.fitNear, "%.0f%%")});
        shLo = std::min(shLo, share);
        shHi = std::max(shHi, share);
        if (c.gluAB > 0) {
            vLo = std::min(vLo, c.fitAB/c.gluAB);
            vHi = std::max(vHi, c.fitAB/c.gluAB);
        }
    }
    f += fixedTable(rows) + "\n";
    f += "- The trace space accounts for " + pct(shLo, "%.0f%%") + " to " + pct(shHi, "%.0f%%") + " of the glued trace error on F, "
         "and the volume error over A u B with the best trace in the same space is " + fmtd("%.2f", vLo) + " to " + fmtd("%.2f", vHi) +
         " times the glued one" + (shLo > 0.5 ? ": the limit is the face trace (3 hats, or 3 hats + J), not the glue" :
         shHi < 0.5 ? ": the glue is the larger part" : "") + ".\n";
    double gLo = kInf, gHi = 0, gsLo = kInf, gsHi = 0;
    for (const Case& c : cases) {
        gLo = std::min(gLo, c.gRms);
        gHi = std::max(gHi, c.gRms);
        gsLo = std::min(gsLo, c.gNoise);
        gsHi = std::max(gsHi, c.gNoise);
    }
    f += "- grad H over A u B: RMS error " + pct(gLo, "%.0f%%") + " to " + pct(gHi, "%.0f%%") + " of 2 g_max against a noise RMS of " +
         pct(gsLo, "%.0f%%") + " to " + pct(gsHi, "%.0f%%") + " (per case in the main table).\n";
    // the edge singularity, per level and g mode across y_f
    f += "- The edge where the crack meets F: the energy between eps and 2 eps from the crack, decomposed minus reference, at y_f = "
         "0.15, 0.30, 0.45, with the four-wedge edge model in brackets (2D leading order; log-divergent as eps -> 0 where the face "
         "data and the crack's +-g disagree):";
    for (const char* lv : {"A (J = 0)", "B (J)"})
        for (const char* g : {"g = 1", "ramp"}) {
            f += std::string(" Level ") + lv + ", " + g + ":";
            bool first = true;
            for (const Case& c : cases)
                if (c.lv == lv && c.g == g) {
                    f += std::string(first ? " " : ", ") + fmtd("%.4f", c.slab) + " (" + fmtd("%.4f", c.edge) + ")";
                    first = false;
                }
            f += std::string(lv) == "B (J)" && std::string(g) == "ramp" ? ".\n" : ";";
        }
    double dqLo = kInf, dqHi = 0;
    for (const Case& c : cases) {
        dqLo = std::min(dqLo, c.dqT);
        dqHi = std::max(dqHi, c.dqT);
    }
    f += "- S_A + S_B is positive definite: smallest eigenvalue " + fmtd("%.4f", evMinLo) + " to " + fmtd("%.4f", evMinHi) +
         " (s.e. <= " + fmtd("%.4f", evMinSe) + "). Quadrature: B's own point set and A's points mirrored into B give S differing by "
         "up to " + pct(quadRel, "%.1f%%") + " of S (up to " + fmtd("%.1f", quadSig) + " s.e.), a quadrature error the noise floors "
         "do not contain; it moves the glued trace on F by " + pct(dqLo) + " to " + pct(dqHi) + " RMS.\n";
    md.insert(mdIntro, findingsBlock(f));
    geoMd = "## Appendix D: geometry of check 3, per configuration\n\nPrinted by the run before each configuration's first solve, "
            "with Level B's face function J.\n\n" + codeBlock(geoTxt);
    return ok;
}

// =====================================================================================
// Check 4: one nodal ring
// =====================================================================================
// The ring of a mesh node c: the 8 tets (c, s_x e_x, s_y e_y, s_z e_z), one per octant s = (s_x, s_y, s_z) in {+-1}^3
// (the last two corners swapped where that makes det J_T > 0), whose union is the octahedron O = {|x| + |y| + |z| <= 1}
// (convex, 8 face planes) for any c inside O. Global vertex ids: v0 = c, v1..v6 = +e_x, -e_x, +e_y, -e_y, +e_z, -e_z.
// Local vertex 0 of every tet is c, so the face opposite it lies on a face of O and the three others are shared with a
// neighbour (12 interior faces in all).
//  - Ring 1 ("congruent"): c = 0. Every J_T is orthogonal with det +1, so xi and X agree up to a rotation and each
//    tet's xi problem is exactly the physical one: any error is decomposition error.
//  - Ring 2 ("real mesh"): c = (0.12, -0.08, 0.10). The 8 tets differ and none is congruent to the reference tet.
// O and the cracks are the same in both rings, so both rings approximate one physical problem per crack.
// Every tet runs the production recipe in its own xi frame with ALL FOUR faces Dirichlet: H0 (crack +-g, faces J) and
// P_k for its 4 corners (crack 0, faces lambda_k). O's faces are Neumann in the true problem; the glue imposes that
// weakly: the 6 outer vertex values are free unknowns beside c_v0 (the natural boundary condition of the energy).
constexpr int kRingNV = 7, kRingNT = 8;
constexpr double kNearRing = 0.1;    // check 4: "near the interior faces", "near the crack-face lines", "near the crack"
constexpr double kProbePull = 1e-3;  // the reference's vertex probes at the outer vertices: pulled this far towards 0
const char* const kRingTetName[kRingNT] = {"T+++", "T-++", "T+-+", "T--+", "T++-", "T-+-", "T+--", "T---"};
const char* const kRingVertName[kRingNV] = {"v0 = c", "v1 = +e_x", "v2 = -e_x", "v3 = +e_y", "v4 = -e_y", "v5 = +e_z",
                                            "v6 = -e_z"};
using Vec7 = Eigen::Matrix<double, kRingNV, 1>;
using Mat7 = Eigen::Matrix<double, kRingNV, kRingNV>;

struct Ring {
    int id = 1;
    std::string name;
    std::array<Vec3, kRingNV> V;
    std::vector<TetFrame> T;                      // octant o: s_x = -1 if bit 0 is set, s_y bit 1, s_z bit 2
    std::vector<std::array<int, 3>> innerFaces;   // the 12 shared faces (global ids, sorted)
};
Ring makeRing(int id)
{
    Ring R;
    R.id = id;
    R.name = id == 1 ? "ring 1 (congruent, c = 0)" : "ring 2 (real mesh, c = (0.12, -0.08, 0.10))";
    R.V = {id == 1 ? Vec3(0, 0, 0) : Vec3(0.12, -0.08, 0.10), Vec3(1, 0, 0), Vec3(-1, 0, 0), Vec3(0, 1, 0), Vec3(0, -1, 0),
           Vec3(0, 0, 1), Vec3(0, 0, -1)};
    for (int o = 0; o < kRingNT; o++) {
        std::array<int, 4> g = {0, o & 1 ? 2 : 1, o & 2 ? 4 : 3, o & 4 ? 6 : 5};
        TetFrame T = frameOn(kRingTetName[o], g, R.V.data());
        if (T.detJ < 0) {   // det J_T = s_x s_y s_z at c = 0, and its sign does not change while c stays inside O
            std::swap(g[2], g[3]);
            T = frameOn(kRingTetName[o], g, R.V.data());
        }
        R.T.push_back(T);
        for (int k = 1; k < 4; k++) {   // the faces through c: opposite local vertices 1, 2, 3
            std::array<int, 3> f{};
            int m = 0;
            for (int i = 0; i < 4; i++)
                if (i != k) f[m++] = g[i];
            std::sort(f.begin(), f.end());
            if (std::find(R.innerFaces.begin(), R.innerFaces.end(), f) == R.innerFaces.end()) R.innerFaces.push_back(f);
        }
    }
    return R;
}

// O: vertices +-e_x, +-e_y, +-e_z (global v1..v6), one outward triangle and one plane s.X <= 1 per octant s
ConvexDomain octahedron()
{
    ConvexDomain O;
    O.v = {Vec3(1, 0, 0), Vec3(-1, 0, 0), Vec3(0, 1, 0), Vec3(0, -1, 0), Vec3(0, 0, 1), Vec3(0, 0, -1)};
    for (int o = 0; o < kRingNT; o++) {
        const Vec3 s(o & 1 ? -1 : 1, o & 2 ? -1 : 1, o & 4 ? -1 : 1);
        std::array<int, 3> f = {o & 1 ? 1 : 0, o & 2 ? 3 : 2, o & 4 ? 5 : 4};
        if ((O.v[f[1]] - O.v[f[0]]).cross(O.v[f[2]] - O.v[f[0]]).dot(s) < 0) std::swap(f[1], f[2]);
        O.tris.push_back(f);
        O.n.push_back(s/std::sqrt(3.0));
        O.c.push_back(1.0/std::sqrt(3.0));
    }
    O.L = 2.0;   // diameter
    return O;
}
double distToDomainFaces(const ConvexDomain& O, const Vec3& X)
{
    double d = kInf;
    for (size_t k = 0; k < O.n.size(); k++) d = std::min(d, O.c[k] - O.n[k].dot(X));
    return d;
}

// The reference's vertex probes: c itself, and the outer vertices pulled kProbePull towards the origin (a walk cannot
// start on O's corner). The decomposed field has its vertex values exactly: they are its face data there.
std::vector<Vec3> ringProbes(const Ring& R)
{
    std::vector<Vec3> p = {R.V[0]};
    for (int k = 1; k < kRingNV; k++) p.push_back((1.0 - kProbePull)*R.V[k]);
    return p;
}

// A crack of check 4: plane n.X = d (n unit, +g on n.X > d), front direction f (u = f projected into the plane), crack
// region = the plane's section of O with u.X <= s, s = the section's extent at t. K3 is cut through (t = 1).
struct RingCrack {
    std::string name, desc;
    Vec3 n = Vec3::Zero(), f = Vec3::Zero();
    double d = 0, t = 0;
    bool through = false;
    double sJ = 0;   // the front line u.X = sJ of Level C's J: the crack's front, or (K3) 0.3 beyond O's section
    Cut cutO;        // the crack on O in X (buildCutOnDomain; L = sqrt(2), the reference tet's)
};
RingCrack ringCrack(int k)
{
    RingCrack c;
    c.name = "K" + std::to_string(k + 1);
    if (k == 1) {
        c.n = Vec3(0.3, -1.0, 0.6).normalized();
        c.d = -0.05;
        c.f = Vec3(1.0, 0.4, 0.0);
        c.desc = "n = normalize(0.3, -1, 0.6), d = -0.05, u = (1, 0.4, 0) projected into the plane; front at t = 0.30 of the "
                 "section's extent along u";
    } else {
        c.n = Vec3(1.0, 0.5, 0.2).normalized();
        c.d = 0.08;
        c.f = Vec3(0.0, 1.0, 0.0);
        c.desc = k == 0 ? "n = normalize(1, 0.5, 0.2), d = 0.08, u = (0, 1, 0) projected into the plane; front at t = 0.50 of "
                          "the section's extent along u"
                        : "K1's plane cut all the way through O (no front inside O); Level C's J has its front line 0.3 "
                          "beyond O's section";
    }
    c.t = k == 0 ? 0.5 : k == 1 ? 0.3 : 1.0;
    c.through = k == 2;
    c.cutO = buildCutOnDomain(octahedron(), c.n, c.d, c.f, kSqrt2, c.t);
    c.sJ = c.through ? c.cutO.sMax + 0.3 : c.cutO.s;
    return c;
}

// Level C's face function: the crack's own 2D front solution, one function of physical X for every tet. a = u.X - sJ
// along the advance (negative behind the front, on the crack side), b = n.X - d across the crack (positive on the +g
// side), theta = atan2(b, a), r = |(a, b)|: J = theta/pi (g = 1), J = min(1, sqrt(r/L)) sin(theta/2) (ramp), L = sqrt(2).
// Both are +-g on the crack (theta = +-pi) and harmonic away from it as long as r < L (a, b are orthonormal coordinates),
// so where a tet's crack is the whole branch cut in the tet, its H0 at Level C is J itself (tested after the runs).
FaceFn ringJ(const RingCrack& rc, GType gt)
{
    return [n = rc.n, d = rc.d, u = rc.cutO.u, s = rc.sJ, gt](const Vec3& X) {
        const double a = u.dot(X) - s, b = n.dot(X) - d, theta = std::atan2(b, a);
        if (gt == GType::One) return theta/std::numbers::pi;
        return std::min(1.0, std::sqrt(std::hypot(a, b)/kSqrt2))*std::sin(0.5*theta);
    };
}

// One ring tet's crack: the xi cut (pulled back from X by buildCut on the reference tet) at the t that puts the front at
// O's front u.X = s, and the xi frame with that front line also where it lies outside the tet (gFront: the crack data's
// g is measured from it, in xi). lines: the crack's segments on the tet's faces (faceLines) and its front segment (front),
// in X: the point set is refined around both.
struct RingTetCrack {
    XiCrack xc;          // cfg, t (clamped to [0, 1]), cut, strips, physFull, frontErr, nFront
    double tRaw = -1;    // the front's t in this tet's extent along u, unclamped (<= 0: the crack misses the tet)
    bool plane = false;  // the plane meets the tet
    Cut gFront;          // the xi frame with s_xi of O's front (the affine map keeps t, also outside [0, 1])
    Cut phys;            // the tet's crack in X (buildCutOnTet on the physical tet at t): geometry only
    std::vector<std::array<Vec3, 2>> faceLines;
    std::vector<std::array<Vec3, 2>> front;   // 0 or 1 segment
    bool cut() const { return !xc.cut.empty(); }
};
RingTetCrack ringTetCrack(const TetFrame& T, const RingCrack& rc)
{
    RingTetCrack r;
    XiCrack& x = r.xc;
    x.cfg.restVerts = T.X;
    x.cfg.planeNormal = rc.n;
    x.cfg.planeOffset = rc.d;
    x.cfg.frontDir = rc.f;
    x.physFull = buildCutOnTet(Tet(T.X), rc.n, rc.d, rc.f, 1.0);
    if (x.physFull.fullPolygon.size() < 3) return r;   // the plane misses this tet
    r.plane = true;
    r.tRaw = (rc.cutO.s - x.physFull.sMin)/(x.physFull.sMax - x.physFull.sMin);
    x.t = std::clamp(r.tRaw, 0.0, 1.0);
    x.cut = buildCut(referenceTet(), x.cfg, x.t);
    x.strips = StripCut(referenceTet(), x.cfg);
    x.strips.addStrip(x.strips.sOf(x.t));
    r.phys = buildCutOnTet(Tet(T.X), rc.n, rc.d, rc.f, x.t);
    const Cut& fr = x.strips.frame();
    r.gFront = x.strips.frontCut(x.strips.nStrips());
    r.gFront.s = fr.sMin + r.tRaw*(fr.sMax - fr.sMin);
    r.gFront.t = r.tRaw;
    if (x.cut.empty()) return r;
    const bool frontInside = r.tRaw > 0 && r.tRaw < 1;
    if (frontInside)
        for (const Vec3& p : frontVertices(x.cut)) {
            x.frontErr = std::max(x.frontErr, std::fabs(rc.cutO.u.dot(T.toX(p)) - rc.cutO.s));
            x.nFront++;
        }
    const Tet ref = referenceTet();
    const std::vector<Vec3>& P = x.cut.polygon;
    for (size_t i = 0; i < P.size(); i++) {
        const Vec3 &a = P[i], &b = P[(i + 1)%P.size()];
        const Eigen::Vector4d la = ref.bary(a), lb = ref.bary(b);
        bool onFace = false;
        for (int o = 0; o < 4; o++) onFace = onFace || (std::fabs(la[o]) <= 1e-12 && std::fabs(lb[o]) <= 1e-12);
        const bool onFront = frontInside && std::fabs(a.dot(x.cut.u) - x.cut.s) <= 1e-9*x.cut.L &&
                             std::fabs(b.dot(x.cut.u) - x.cut.s) <= 1e-9*x.cut.L;
        if (onFace) r.faceLines.push_back({T.toX(a), T.toX(b)});
        else if (onFront) r.front.push_back({T.toX(a), T.toX(b)});
    }
    return r;
}

double pointTriangleDist(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c)
{
    const Vec3 n = (b - a).cross(c - a).normalized();
    const Vec3 q = p - n.dot(p - a)*n;
    auto inner = [&](const Vec3& u, const Vec3& v) { return (v - u).cross(q - u).dot(n) >= 0; };
    if (inner(a, b) && inner(b, c) && inner(c, a)) return std::fabs(n.dot(p - a));
    return std::min({pointSegmentDist(p, a, b), pointSegmentDist(p, b, c), pointSegmentDist(p, c, a)});
}

// One ring tet's fixed point set, built in xi like volumePoints (Bey leaves, stratified Halton points, weights V_l/m_l
// times |det J|), refined to depth 3 and 4 around the crack's front segment AND the lines where the crack meets the
// tet's faces (leaf centroid closer to one of them than half the leaf's diameter, both in X): the integrands of S and b
// are ~1/r^2 at a crack-face line where the face data and the crack's disagree (check 3's quadrature error of up to 6%
// in S came from refining at the front only). Then the points within eps (physical) of O's crack are dropped.
PointSet ringPoints(const TetFrame& T, const RingTetCrack& rt, const Cut& cutO, int target, double eps, int mMin)
{
    PointSet ps;
    ps.eps = eps;
    const Tet ref = referenceTet();
    std::vector<std::array<Vec3, 2>> lines = rt.faceLines;
    lines.insert(lines.end(), rt.front.begin(), rt.front.end());
    std::vector<Leaf> leaves = {Leaf{{ref.v[0], ref.v[1], ref.v[2], ref.v[3]}, 0}};
    auto refine = [&](auto pred) {
        std::vector<Leaf> out;
        for (const Leaf& l : leaves) {
            if (!pred(l)) {
                out.push_back(l);
                continue;
            }
            for (const Leaf& c : beyChildren(l)) out.push_back(c);
        }
        leaves.swap(out);
    };
    refine([](const Leaf&) { return true; });
    refine([](const Leaf&) { return true; });
    for (int level = 2; level < 4 && !lines.empty(); level++)
        refine([&](const Leaf& l) {
            if (l.depth != level) return false;
            const Vec3 c = T.toX(0.25*(l.v[0] + l.v[1] + l.v[2] + l.v[3]));
            double diam = 0;
            for (int i = 0; i < 4; i++)
                for (int j = i + 1; j < 4; j++) diam = std::max(diam, (T.J*(l.v[i] - l.v[j])).norm());
            for (const auto& s : lines)
                if (pointSegmentDist(c, s[0], s[1]) < 0.5*diam) return true;
            return false;
        });
    auto counts = [&](double c) {
        std::vector<int> m(leaves.size());
        for (size_t i = 0; i < leaves.size(); i++)
            m[i] = std::max(mMin, int(std::lround(c*std::pow(8.0, -0.5*(leaves[i].depth - 2)))));
        return m;
    };
    auto total = [](const std::vector<int>& m) {
        long long s = 0;
        for (int v : m) s += v;
        return s;
    };
    double lo = 1, hi = double(target);
    for (int it = 0; it < 60; it++) {
        const double mid = 0.5*(lo + hi);
        (total(counts(mid)) < target ? lo : hi) = mid;
    }
    const std::vector<int> m = counts(hi);
    uint64_t idx = 1;
    ps.nLeaves = int(leaves.size());
    for (size_t li = 0; li < leaves.size(); li++) {
        const Leaf& l = leaves[li];
        const double wl = leafVolume(l)/double(m[li])*std::fabs(T.detJ);
        ps.leavesByDepth[l.depth]++;
        ps.pointsByDepth[l.depth] += m[li];
        for (int k = 0; k < m[li]; k++, idx++) {
            const Vec3 b = cubeToSimplex(radicalInverse(idx, 2), radicalInverse(idx, 3), radicalInverse(idx, 5));
            const Vec3 p = l.v[0] + b.x()*(l.v[1] - l.v[0]) + b.y()*(l.v[2] - l.v[0]) + b.z()*(l.v[3] - l.v[0]);
            ps.allXi.push_back(p);
            ps.allW.push_back(wl);
            if (cutO.distToCut(T.toX(p)) < eps) {   // physical distance to O's crack (J_T is not orthogonal in ring 2)
                ps.nDropped++;
                continue;
            }
            ps.xi.push_back(p);
            ps.X.push_back(T.toX(p));
            ps.w.push_back(wl);
        }
    }
    return ps;
}

// The eps-tube of O's crack inside one ring tet (physical X), for the quadrature check: QMC in the box spanned by the
// tet's extent along the crack frame's u and w and |n.X - d| < eps (the tube lies in it), points inside the tet and
// within eps of O's crack.
void ringTubeIntegrals(const TetFrame& T, const Cut& cutO, double eps, int nPts, const uint64_t bases[3], double out[kNTest])
{
    std::fill(out, out + kNTest, 0.0);
    double uLo = kInf, uHi = -kInf, wLo = kInf, wHi = -kInf;
    for (const Vec3& p : T.X) {
        uLo = std::min(uLo, p.dot(cutO.u));
        uHi = std::max(uHi, p.dot(cutO.u));
        wLo = std::min(wLo, p.dot(cutO.w));
        wHi = std::max(wHi, p.dot(cutO.w));
    }
    const double boxVol = (uHi - uLo)*(wHi - wLo)*2.0*eps;
    const int nChunks = 64;
    std::vector<std::array<double, kNTest>> part(nChunks);
    tbb::parallel_for(0, nChunks, [&](int c) {
        std::array<double, kNTest> s{};
        for (int i = int((long long)nPts*c/nChunks); i < int((long long)nPts*(c + 1)/nChunks); i++) {
            const uint64_t idx = uint64_t(i) + 1;
            const double a = uLo + (uHi - uLo)*radicalInverse(idx, bases[0]);
            const double b = wLo + (wHi - wLo)*radicalInverse(idx, bases[1]);
            const double h = eps*(2.0*radicalInverse(idx, bases[2]) - 1.0);
            const Vec3 X = cutO.d*cutO.n + a*cutO.u + b*cutO.w + h*cutO.n;
            const Vec3 xi = T.toXi(X);
            if (std::min({xi.x(), xi.y(), xi.z(), 1.0 - xi.sum()}) < 0 || cutO.distToCut(X) >= eps) continue;
            for (int k = 0; k < kNTest; k++) s[k] += testFn(k, X);
        }
        part[c] = s;
    });
    for (int c = 0; c < nChunks; c++)
        for (int k = 0; k < kNTest; k++) out[k] += part[c][k];
    for (int k = 0; k < kNTest; k++) out[k] *= boxVol/double(nPts);
}

// The ring reference on O at the points of the 8 tets (their kept points, tet by tet, physical X; H and grad H) and at the
// vertex probes (H only), in the same cache as checks 2 and 3. Its key names the ring, the crack (plane, front direction,
// t), seed, N, batches, the points and their hash, O's planes, the live walk parameters, the crack frame and a tag for
// the walk code on a convex domain. A pilot first (N/100 walks), as for U; a run the pilot puts above 3 hours is not
// started (N is never lowered).
const char* const kRingRefSolverVersion =
    "own loop (walk_cache.cpp runWalk + solveMulti on a ConvexDomain) v1: double; K = 2 (+-g with g types a, b)";
WalkParams ringReferenceWalkParams(const ConvexDomain& O)
{
    WalkParams wp = gateWalkParams(allNeumann(), true, false);   // the mirror ball across O's faces, no Dirichlet face
    wp.absorbShellRel = kGateShellRel*kSqrt2/O.L;               // the absolute shell of every gate solve, 1e-6 sqrt(2)
    return wp;
}
StripCut ringStrips(const RingCrack& rc)
{
    Cut frame = buildCutOnDomain(octahedron(), rc.n, rc.d, rc.f, kSqrt2, 1.0);
    StripCut sc = StripCut::physicalUnionReference(frame);
    sc.addStrip(sc.sOf(rc.t));
    return sc;
}
RefRun runRingReference(Report& out, const GateParams& gp, const ConvexDomain& O, const Ring& R, const RingCrack& rc,
                        const std::string& tag, uint64_t seedTag, const std::vector<PointSet>& ps,
                        const std::vector<Vec3>& probes, const std::string& cacheDir)
{
    RefRun r;
    for (const PointSet& p : ps) r.X.insert(r.X.end(), p.X.begin(), p.X.end());
    r.nA = r.X.size();   // check 4: nA = the tets' points, nF = the probes
    r.X.insert(r.X.end(), probes.begin(), probes.end());
    r.nF = probes.size();
    const StripCut sc = ringStrips(rc);
    const WalkParams wp = ringReferenceWalkParams(O);
    const Cut front = sc.frontCut(sc.nStrips());
    const TerminalValues tv = [&front](const WalkTerminal& t, double* u) {   // O's faces are Neumann: on the cut
        u[0] = double(t.side)*front.g(t.x, GType::One);
        u[1] = double(t.side)*front.g(t.x, GType::SqrtFront);
    };
    const uint64_t seedBase = mixSeed(gp.seed, seedTag);
    char hash[32];
    std::snprintf(hash, sizeof hash, "%016llx", (unsigned long long)fnv1a(r.X.data(), r.X.size()*sizeof(Vec3)));
    std::string key = "coupling gate check 4 reference on O | " + R.name + ", c " + fmtg3(R.V[0]) + " | crack " + rc.name +
                      ": n " + fmtg3(rc.n) + ", d " + fmtg(rc.d) + ", f " + fmtg3(rc.f) + ", t " + fmtg(rc.t) + " | seed " +
                      std::to_string(seedBase) + " | N " + std::to_string(gp.nRef) + " | batches " +
                      std::to_string(gp.nRefBatches) + " | " + std::to_string(r.X.size()) + " points (tets " +
                      std::to_string(r.nA) + ", probes " + std::to_string(r.nF) + "), hash " + hash + " | O";
    for (size_t k = 0; k < O.n.size(); k++) key += " " + fmtg3(O.n[k]) + "." + fmtg(O.c[k]);
    const Cut& f = sc.frame();
    key += ", L_O " + fmtg(O.L) + " | walk: absorbShellRel " + fmtg(wp.absorbShellRel) + " (absolute " +
           fmtg(wp.absorbShellRel*O.L) + "), maxWalkLength " + std::to_string(wp.maxWalkLength) + ", mirrorGradBall " +
           (wp.mirrorGradBall ? "1" : "0") + ", kGradFirstBallShrink " + fmtg(kGradFirstBallShrink) + " | crack frame: n " +
           fmtg3(f.n) + ", d " + fmtg(f.d) + ", u " + fmtg3(f.u) + ", sMin " + fmtg(f.sMin) + ", sMax " + fmtg(f.sMax) +
           ", L of g " + fmtg(f.L) + ", strips " + std::to_string(sc.nStrips()) + ", front s " + fmtg(sc.front()) + ", g types " +
           gTypeName(GType::One) + " / " + gTypeName(GType::SqrtFront) + " | " + kRingRefSolverVersion;
    r.key = key;
    r.file = (std::filesystem::path(cacheDir)/("ref_" + tag + ".bin")).string();
    out("reference on O (own loop, double, 8 Neumann faces; N = %d per point in %d batches, both g types from the same walks;\n"
        "H and grad H at the %zu points of the 8 tets, H at the %zu vertex probes)\n", gp.nRef, gp.nRefBatches, r.nA, r.nF);
    if (!gp.quick && readRefCache(r)) {
        r.cached = true;
        out("  loaded from the cache %s (key matches: ring, crack, seed, N, batches, points hash, O, the live walk parameters,\n"
            "  the crack frame, solver version)\n  original run: %s\n", r.file.c_str(), r.stats.summary().c_str());
        return r;
    }
    std::vector<char> grad(r.X.size(), 1);
    std::fill(grad.begin() + std::ptrdiff_t(r.nA), grad.end(), char(0));
    QueryPoints q;
    q.add(r.X, &grad, seedBase);
    r.nPilot = gp.nRefBatches*std::max(1, gp.nRef/(100*gp.nRefBatches));
    {
        QueryPoints qp;
        qp.add(r.X, &grad, mixSeed(gp.seed, seedTag + 100));
        SolveStats pst;
        solveMulti(O, sc, qp, r.nPilot, gp.nRefBatches, wp, 2, tv, &pst);
        r.pilotSeconds = pst.seconds;
        r.pilotExtrapolated = pst.seconds*double(gp.nRef)/double(r.nPilot);
        out("  timing pilot, N = %d: %s\n  -> extrapolated to N = %d: %.0f s\n", r.nPilot, pst.summary().c_str(), gp.nRef,
            r.pilotExtrapolated);
    }
    if (r.pilotExtrapolated > 3*3600.0) {
        r.ok = false;
        out("  NOT RUN: the extrapolated time exceeds 3 hours (N is not lowered)\n");
        return r;
    }
    r.sums = solveMulti(O, sc, q, gp.nRef, gp.nRefBatches, wp, 2, tv, &r.stats);
    out("  run: %s\n", r.stats.summary().c_str());
    if (!gp.quick) {
        std::error_code ec;
        std::filesystem::create_directories(cacheDir, ec);
        out("  %s %s\n", writeRefCache(r) ? "cached to" : "COULD NOT WRITE the cache", r.file.c_str());
    }
    return r;
}

// Zombie (float) on O at the ring reference's points, as zombieReference does on U: H at every point (one solve of N =
// nZombie per g type), grad H at the points farther than kRefZGradBand from O's faces (nZombieBatches solves), the same
// verdict, and H per region (inRegion; kRingZRegionName). nModes = 1: g = 1 only.
const char* const kRingZRegionName[kRefZRegions] = {"tets with the crack", "tets without it", "vertex probes",
                                                    "near interior faces (0.1)", "near the crack (0.1)"};
RefZombie zombieRingReference(const GateParams& gp, const ConvexDomain& O, const RingCrack& rc, const RefRun& r,
                              const std::vector<std::array<char, kRefZRegions>>& inRegion, int nModes, uint64_t seedTag)
{
    RefZombie z;
    WalkParams wp;   // every face Neumann; the reference's absolute shell, minimum star radius 1/10 of it
    wp.absorbShellRel = kGateShellRel*kSqrt2/O.L;
    wp.reflectShellRel = 0.1*wp.absorbShellRel;
    std::vector<size_t> gi;
    std::vector<Vec3> gX;
    for (size_t i = 0; i < r.nA; i++)
        if (distToDomainFaces(O, r.X[i]) > kRefZGradBand) {
            gi.push_back(i);
            gX.push_back(r.X[i]);
        }
    z.nGradPts = (long long)gi.size();
    const std::vector<char> gradFlags(gX.size(), 1);
    const GType gts[2] = {GType::One, GType::SqrtFront};
    const int S = gp.nZombieBatches, nPer = gp.nZombie/gp.nZombieBatches, Bref = r.sums.nBatches;
    for (int m = 0; m < nModes; m++) {
        const CutSolver solver(O, rc.cutO, gts[m], wp);
        std::vector<Estimate> e;
        const SolveStats sH = solver.solve(r.X, nullptr, gp.nZombie, mixSeed(gp.seed, seedTag + m), e);
        z.secondsH += sH.seconds;
        z.st.add(sH);
        const std::vector<double> coef = unitCoef(2, m);
        for (size_t i = 0; i < r.X.size(); i++) {
            double H, var;
            long long n;
            comboH(r.sums, i, coef, H, var, n);
            const double diff = H - e[i].H, v = var + double(e[i].Hse)*e[i].Hse;
            z.H[m].add(i, diff, v);
            for (int k = 0; k < kRefZRegions; k++)
                if (inRegion[i][k]) z.Hr[m][k].add(i, diff, v);
        }
        std::vector<std::vector<Estimate>> eg(static_cast<size_t>(S));
        for (int s = 0; s < S; s++) {
            const SolveStats sG = solver.solve(gX, &gradFlags, nPer, mixSeed(gp.seed, seedTag + 10 + 32*m + s), eg[s]);
            z.secondsG += sG.seconds;
            z.st.add(sG);
        }
        for (size_t j = 0; j < gi.size(); j++) {
            Vec3 zm = Vec3::Zero(), zv = Vec3::Zero();
            bool finite = true;
            for (int s = 0; s < S; s++) {
                finite = finite && eg[s][j].grad.allFinite();
                zm += eg[s][j].grad.cast<double>()/double(S);
            }
            if (!finite) {
                z.nNaN++;
                continue;
            }
            for (int s = 0; s < S; s++) zv += (eg[s][j].grad.cast<double>() - zm).cwiseAbs2()/(double(S)*double(S - 1));
            Vec3 om = Vec3::Zero(), ov = Vec3::Zero();
            std::vector<Vec3> g(static_cast<size_t>(Bref));
            for (int b = 0; b < Bref; b++) om += (g[size_t(b)] = comboBatchGrad(r.sums, gi[j], b, coef))/double(Bref);
            for (const Vec3& x : g) ov += (x - om).cwiseAbs2()/(double(Bref)*double(Bref - 1));
            for (int a = 0; a < 3; a++) z.G[m][a].add(gi[j], om[a] - zm[a], ov[a] + zv[a]);
        }
    }
    z.seconds = z.st.seconds;
    return z;
}

// One ring tet's two runs, all four faces Dirichlet (the production setting), on the same points:
//  - the recipe run: H0 at Level A for both g types (crack +-g, faces 0) and P_k for the tet's four corners (crack 0,
//    faces lambda_k, the hat of corner k on every face), all from the same walks; every face datum is linear, so the
//    gradient's first ball may cross the nearest Dirichlet face (odd mirror);
//  - the Level-C run: H0 for both g types with faces J_g (ringJ). J is harmonic in the tet minus its crack but not in
//    all of space (its branch cut continues into the neighbours), so the plain first ball.
constexpr int kRingP0 = 2, kRingH0C = 6;   // Gram layout: H0A g1, H0A g2, P_0..P_3, H0C g1, H0C g2
std::vector<FieldSpec> ringRecipeFields(const TetFrame& T)
{
    std::vector<FieldSpec> f = {FieldSpec{"H0A g1", true, GType::One, {}}, FieldSpec{"H0A g2", true, GType::SqrtFront, {}}};
    for (int k = 0; k < 4; k++)
        f.push_back(FieldSpec{"P" + std::to_string(k), false, GType::One, [&T, k](const Vec3& X) { return lambdaOf(T, k, X); }});
    return f;
}
std::vector<FieldSpec> ringLevelCFields(const FaceFn J[2])
{
    return {FieldSpec{"H0C g1", true, GType::One, J[0], false}, FieldSpec{"H0C g2", true, GType::SqrtFront, J[1], false}};
}

// The Gram of one ring tet's 8 fields over its points (weights w = the physical weights): phys = the physical energy,
// grad_X = J^-T grad_xi integrated over the physical tet; else the xi energy, grad_xi . grad_xi integrated over the
// reference tet (weights w/|det J|). The same gradient estimates either way; the two runs' batches are paired by index,
// so every product is one of independent estimates and batch a of both runs is one unit of the jackknife.
Gram gramRing(const OwnRun& rec, const OwnRun& hc, const std::vector<double>& w, const TetFrame& T, bool phys)
{
    const int K1 = rec.sums.K, K2 = hc.sums.K;
    if (rec.sums.nBatches != hc.sums.nBatches || rec.sums.nPoints != hc.sums.nPoints)
        throw std::logic_error("gramRing: the runs differ in batches or points");
    std::vector<std::vector<double>> u1, u2;
    for (int k = 0; k < K1; k++) u1.push_back(unitCoef(K1, k));
    for (int k = 0; k < K2; k++) u2.push_back(unitCoef(K2, k));
    std::vector<double> wx(w);
    if (!phys)
        for (double& v : wx) v /= std::fabs(T.detJ);
    const Eigen::Matrix3d M = phys ? T.JinvT : Eigen::Matrix3d::Identity();
    return gramFromBatches(wx, rec.sums.nBatches, K1 + K2, [&](size_t p, int b, Vec3* g) {
        for (int q = 0; q < K1; q++) g[q] = M*comboBatchGrad(rec.sums, p, b, u1[q]);
        for (int q = 0; q < K2; q++) g[K1 + q] = M*comboBatchGrad(hc.sums, p, b, u2[q]);
    });
}

// The glue of check 4, written out. Unknowns: c = (c_v0, ..., c_v6), the values at the ring's 7 vertices (global ids):
// the centre node v0 and the six outer vertices. O's faces are Neumann in the true problem; the glue imposes that
// weakly: the outer vertex values are free unknowns (the natural boundary condition of the energy). Tet T's P_k
// answers the hat of its local corner k (lambda_k, evaluated at X = X0 + J_T xi), so its S_T (4 x 4: S_T(i, j) =
// int grad P_i . grad P_j) and b_T (b_T(i) = int grad H0 . grad P_i) are indexed by T's local corners, local corner i
// being the global vertex T.global[i]. With c_T = (c_{T.global[0]}, ..., c_{T.global[3]}) the total energy is
//   E(c) = sum_T [ (1/2) c_T^T S_T c_T + b_T^T c_T + const ],
// and dE/dc = 0 gives the 7 x 7 system M c = r with
//   M = sum_T P_T^T S_T P_T,  r = -sum_T P_T^T b_T  (P_T: the 4 x 7 scatter matrix, P_T(i, T.global[i]) = 1),
// i.e. for every tet and every pair of its local corners i, j:
//   M(T.global[i], T.global[j]) += S_T(i, j),   r(T.global[i]) -= b_T(i).
// Each S_T is symmetrized as it is scattered. M is symmetric positive definite once the crack cuts a tet: there
// sum_k P_k is not 1 (the crack's data is 0), so 1^T S_T 1 > 0 and the constant is not in M's null space. Eigen, double.
Vec7 glueCheck4(const Ring& R, const std::vector<const Eigen::MatrixXd*>& G, int h, Mat7* Mout = nullptr, Vec7* rhsOut = nullptr)
{
    Mat7 M = Mat7::Zero();
    Vec7 r = Vec7::Zero();
    for (int t = 0; t < kRingNT; t++) {
        const Eigen::MatrixXd& g = *G[size_t(t)];
        const std::array<int, 4>& id = R.T[size_t(t)].global;
        for (int i = 0; i < 4; i++) {
            r[id[i]] -= g(h, kRingP0 + i);                                                         // -b_T(i)
            for (int j = 0; j < 4; j++) M(id[i], id[j]) += 0.5*(g(kRingP0 + i, kRingP0 + j) + g(kRingP0 + j, kRingP0 + i));   // S_T(i, j)
        }
    }
    if (Mout) *Mout = M;
    if (rhsOut) *rhsOut = r;
    return M.ldlt().solve(r);
}

// The glued vertex values of one particular solution h (0, 1: Level A; kRingH0C + 0, 1: Level C) from the 8 tets' Grams,
// its delete-one-batch replicates (the system re-solved with replicate a of every tet's Gram) and their jackknife s.e.;
// seconds: the one solve (assembly included); secondsJack: the replicates'
struct Glue4 {
    Vec7 c = Vec7::Zero(), se = Vec7::Zero(), rhs = Vec7::Zero();
    Mat7 M = Mat7::Zero();
    std::vector<Vec7> jack;
    double seconds = 0, secondsJack = 0;
};
Glue4 glueRun4(const Ring& R, const std::vector<Gram>& G, int h)
{
    Glue4 gl;
    std::vector<const Eigen::MatrixXd*> g(kRingNT);
    auto t0 = std::chrono::steady_clock::now();
    for (int t = 0; t < kRingNT; t++) g[size_t(t)] = &G[size_t(t)].G;
    gl.c = glueCheck4(R, g, h, &gl.M, &gl.rhs);
    gl.seconds = secondsSince(t0);
    t0 = std::chrono::steady_clock::now();
    for (int a = 0; a < G[0].B; a++) {
        for (int t = 0; t < kRingNT; t++) g[size_t(t)] = &G[size_t(t)].jack[size_t(a)];
        gl.jack.push_back(glueCheck4(R, g, h));
    }
    gl.secondsJack = secondsSince(t0);
    for (int k = 0; k < kRingNV; k++) {
        std::vector<double> rep;
        for (const Vec7& x : gl.jack) rep.push_back(x[k]);
        gl.se[k] = jackSe(rep);
    }
    return gl;
}
Eigen::Vector4d localValues(const TetFrame& T, const Vec7& c)
{
    return Eigen::Vector4d(c[T.global[0]], c[T.global[1]], c[T.global[2]], c[T.global[3]]);
}

// (1/2) int |grad(H0_h + sum_k c_k P_k)|^2 over a ring tet's points, from its Gram (gramRing layout)
double tetEnergy4(const Eigen::MatrixXd& G, int h, const Eigen::Vector4d& c)
{
    const int idx[5] = {h, kRingP0, kRingP0 + 1, kRingP0 + 2, kRingP0 + 3};
    const double v[5] = {1.0, c[0], c[1], c[2], c[3]};
    double e = 0;
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) e += v[i]*v[j]*G(idx[i], idx[j]);
    return 0.5*e;
}

// The decomposed field of one ring tet at level L (0: the recipe run's H0; 1: the Level-C run's H0), g type m, with its
// local corner values c (cJ: their jackknife replicates, which enter the variance; null: fixed)
FieldAt decField4(const TetFrame& T, const OwnRun& rec, const OwnRun& hc, int L, int m, const Eigen::Vector4d& c,
                  const std::vector<Eigen::Vector4d>* cJ)
{
    std::vector<std::vector<double>> pj;
    if (cJ)
        for (const Eigen::Vector4d& x : *cJ)
            pj.push_back({L == 0 && m == 0 ? 1.0 : 0.0, L == 0 && m == 1 ? 1.0 : 0.0, x[0], x[1], x[2], x[3]});
    const std::vector<double> pc = {L == 0 && m == 0 ? 1.0 : 0.0, L == 0 && m == 1 ? 1.0 : 0.0, c[0], c[1], c[2], c[3]};
    if (L == 0) return runField(T, {RunTerm{&rec, pc, pj}});
    return runField(T, {RunTerm{&hc, unitCoef(2, m), {}}, RunTerm{&rec, pc, pj}});
}

// The geometry of check 4, printed before its first solve: the rings (vertices, tets with their global ids and det J,
// J's orthogonality in ring 1, each tet's outer face on O, the volumes, the interior faces), the cracks on O (polygon,
// front, its share of the section, J's front line and r/L, the node's distance to the crack, J at the vertices), and per
// ring and crack every tet's crack (cut or not, t, the polygons in xi and X, the front at u.X = s in X, the union of the
// tets' polygons = O's, and the +g side and g of the tets' xi cuts against O's cut at random points). False if an
// assertion failed.
bool printRingGeometry(Report& out, const GateParams& gp, const ConvexDomain& O, const std::vector<Ring>& rings,
                       const std::vector<RingCrack>& cracks)
{
    bool ok = true;
    out("geometry (physical coordinates X), printed before the first solve:\n");
    out("  O = {|x| + |y| + |z| <= 1}: vertices v1..v6 = +e_x, -e_x, +e_y, -e_y, +e_z, -e_z; 8 face planes s.X <= 1, s in\n"
        "  {+-1}^3 (outward unit normals s/sqrt(3)); volume 4/3, diameter L_O = %.0f\n", O.L);
    for (const Ring& R : rings) {
        out("\n  %s: vertices", R.name.c_str());
        for (int k = 0; k < kRingNV; k++) out("%s %s %s", k ? "," : "", kRingVertName[k], fmt3(R.V[k], 2).c_str());
        out("\n");
        double vol = 0, orthMax = 0, outerErr = 0;
        bool posDet = true;
        for (const TetFrame& T : R.T) {
            vol += std::fabs(T.detJ)/6.0;
            const double orth = (T.J.transpose()*T.J - Eigen::Matrix3d::Identity()).norm();
            orthMax = std::max(orthMax, orth);
            posDet = posDet && T.detJ > 0;
            const Vec3 s = T.X[1] + T.X[2] + T.X[3];   // = (s_x, s_y, s_z): O's face s.X = 1 holds the outer face
            for (int k = 1; k < 4; k++) outerErr = std::max(outerErr, std::fabs(s.dot(T.X[k]) - 1.0));
            out("    %s: local 0..3 = v%d v%d v%d v%d, det J = %+.6f, |J^T J - I| = %.1e; face opposite local 0 on O's face "
                "(%+.0f,%+.0f,%+.0f).X = 1\n", T.name, T.global[0], T.global[1], T.global[2], T.global[3], T.detJ, orth, s.x(),
                s.y(), s.z());
        }
        const bool vOk = std::fabs(vol - 4.0/3.0) < 1e-12 && posDet && outerErr == 0 && R.innerFaces.size() == 12;
        const bool oOk = R.id != 1 || orthMax == 0;
        ok = ok && vOk && oOk;
        out("    sum of the tets' volumes %.15f (4/3 = %.15f), every det J > 0, every outer face on O, %zu interior faces -> %s\n",
            vol, 4.0/3.0, R.innerFaces.size(), vOk ? "OK" : "FAILED");
        if (R.id == 1)
            out("    ring 1: every J_T orthogonal with det +1 (max |J^T J - I| = %.1e) -> %s: xi distances, eps and g are the "
                "physical ones\n", orthMax, oOk ? "OK" : "FAILED");
        else
            out("    ring 2: max |J^T J - I| = %.2f: the tets are not congruent to the reference tet, xi distances differ from X\n",
                orthMax);
        std::string f;
        for (const auto& t : R.innerFaces) f += " v" + std::to_string(t[0]) + "v" + std::to_string(t[1]) + "v" + std::to_string(t[2]);
        out("    interior (shared) faces:%s\n", f.c_str());
    }
    for (const RingCrack& rc : cracks) {
        const Cut& c = rc.cutO;
        out("\n  crack %s: %s\n", rc.name.c_str(), rc.desc.c_str());
        out("    plane n = %s, d = %+.4f (+g on n.X > d); u = %s, w = n x u = %s\n", fmt3(c.n).c_str(), c.d, fmt3(c.u).c_str(),
            fmt3(c.w).c_str());
        double fullArea = 0;
        for (size_t i = 1; i + 1 < c.fullPolygon.size(); i++)
            fullArea += 0.5*(c.fullPolygon[i] - c.fullPolygon[0]).cross(c.fullPolygon[i + 1] - c.fullPolygon[0]).norm();
        std::string poly;
        for (const Vec3& p : c.polygon) poly += " " + fmt3(p, 3);
        out("    O's section: %zu vertices, extent along u [%+.6f, %+.6f], area %.6f; the crack u.X <= s = %+.6f (t = %.2f):\n"
            "    %zu vertices, area %.6f (%.1f%% of the section), L of g = %.6f; polygon (X):%s\n", c.fullPolygon.size(), c.sMin,
            c.sMax, fullArea, c.s, c.t, c.polygon.size(), c.area(), 100*c.area()/fullArea, c.L, poly.c_str());
        if (!rc.through) {
            const std::vector<Vec3> fv = frontVertices(c);
            out("    front line {n.X = d, u.X = s}: through %s along w; its segment in O from %s to %s (%zu polygon vertices on it)\n",
                fmt3(c.frontLinePoint()).c_str(), fv.size() > 0 ? fmt3(fv[0], 3).c_str() : "-",
                fv.size() > 1 ? fmt3(fv[1], 3).c_str() : "-", fv.size());
            ok = ok && fv.size() == 2;
        } else {
            out("    no front inside O: the whole section is cut (t = 1, s = s_max). Level C's J has its front line at u.X = %+.4f\n"
                "    (s_max + 0.3), outside O; with g = 1 and O's faces Neumann the exact solution is H = +1 / -1 on the two sides\n",
                rc.sJ);
        }
        {   // the reference walks' strip cut (ringStrips, through StripCut's chains) against this polygon: the same crack
            const StripCut sc = ringStrips(rc);
            std::mt19937_64 rng(mixSeed(gp.seed, 10900 + uint64_t(&rc - cracks.data())));
            std::uniform_real_distribution<double> Ud(-1.0, 1.0);
            double dMax = 0;
            int n = 0;
            while (n < 20000) {
                const Vec3 X(Ud(rng), Ud(rng), Ud(rng));
                if (X.cwiseAbs().sum() > 1.0) continue;
                n++;
                dMax = std::max(dMax, std::fabs(std::sqrt(sc.closest(X, sc.nStrips()).d2) - c.distToCut(X)));
            }
            const double fErr = std::fabs(sc.front() - c.s);
            const bool sOk = dMax < 1e-12 && fErr == 0 && sc.frame().L == kSqrt2 && sc.nStrips() == 1;
            ok = ok && sOk;
            out("    the reference walks' strip cut on O: distance to it minus distance to this polygon <= %.1e at %d random points\n"
                "    of O; its front s - this s = %.1e; L of g = sqrt(2): %s -> %s\n", dMax, n, fErr, sc.frame().L == kSqrt2 ? "yes" : "NO",
                sOk ? "OK" : "FAILED");
        }
        double rMax = 0;   // r of J (the distance to J's front line) is convex, so its max over O is at a vertex
        for (const Vec3& v : O.v) rMax = std::max(rMax, std::hypot(c.u.dot(v) - rc.sJ, c.n.dot(v) - c.d));
        out("    J's r (distance to its front line) over O: at most %.4f = %.3f L%s\n", rMax, rMax/kSqrt2,
            rMax < kSqrt2 ? " < L: the ramp's J is sqrt(r/L) sin(theta/2), unclipped, harmonic off the crack"
                          : ": the ramp's min(1, .) clips J somewhere in O (not harmonic there)");
        const FaceFn J[2] = {ringJ(rc, GType::One), ringJ(rc, GType::SqrtFront)};
        for (const Ring& R : rings) {
            const double hc = c.n.dot(R.V[0]) - c.d, dc = c.distToCut(R.V[0]);
            out("    %s: node c: n.c - d = %+.4f, distance to the crack %.4f (eps = %.4f)%s\n", R.name.c_str(), hc, dc, gp.eps,
                dc < gp.eps ? "  <- WITHIN eps of the crack" : dc < 2*kNearRing ? "  <- close to the crack" : "");
            std::string s;
            for (int k = 0; k < kRingNV; k++) s += fmtd(" %+.4f", J[0](R.V[k])) + "/" + fmtd("%+.4f", J[1](R.V[k]));
            out("      J at v0..v6 (g = 1 / ramp):%s\n", s.c_str());
        }
    }
    for (size_t ri = 0; ri < rings.size(); ri++)
        for (size_t ki = 0; ki < cracks.size(); ki++) {
            const Ring& R = rings[ri];
            const RingCrack& rc = cracks[ki];
            out("\n  %s, crack %s, per tet (t = the front's place in the tet's own extent along u, unclamped):\n", R.name.c_str(),
                rc.name.c_str());
            double areaSum = 0, polyErr = 0, frontErr = 0;
            int outside = 0, nCut = 0;
            bool frontOk = true;
            std::vector<RingTetCrack> rt;
            for (const TetFrame& T : R.T) rt.push_back(ringTetCrack(T, rc));
            for (size_t t = 0; t < rt.size(); t++) {
                const TetFrame& T = R.T[t];
                const RingTetCrack& r = rt[t];
                if (!r.plane) {
                    out("    %s: the plane misses it\n", T.name);
                    continue;
                }
                if (!r.cut()) {
                    out("    %s: the plane crosses it, the crack does not reach it (t = %+.4f)\n", T.name, r.tRaw);
                    continue;
                }
                nCut++;
                const Cut& x = r.xc.cut;
                out("    %s: CUT, t = %+.4f%s; xi plane n = %s, d = %+.6f, u = %s, s = %+.6f, L = %.4f; %zu crack-face segment(s)%s\n",
                    T.name, r.tRaw, r.tRaw >= 1 ? " (>= 1: the whole section)" : "", fmt3(x.n).c_str(), x.d, fmt3(x.u).c_str(), x.s,
                    x.L, r.faceLines.size(), r.front.empty() ? "" : ", the front segment");
                std::string sxi, sx, sp;
                for (const Vec3& p : x.polygon) {
                    sxi += " " + fmt3(p, 3);
                    sx += " " + fmt3(T.toX(p), 3);
                }
                out("       polygon (xi):%s\n       polygon (X): %s\n", sxi.c_str(), sx.c_str());
                if (r.tRaw > 0 && r.tRaw < 1) {
                    out("       front vertices: %d, max |u.X - s| in X = %.1e\n", r.xc.nFront, r.xc.frontErr);
                    frontOk = frontOk && r.xc.nFront == 2 && r.xc.frontErr < 1e-12;
                    frontErr = std::max(frontErr, r.xc.frontErr);
                }
                areaSum += r.phys.area();
                for (const Vec3& p : r.phys.polygon)
                    outside += !(rc.cutO.polygonContains(p) || rc.cutO.distToPolygonEdges(p) < 1e-12);
                for (const auto& [P, Q, toX] : {std::tuple<const std::vector<Vec3>*, const std::vector<Vec3>*, bool>{&x.polygon, &r.phys.polygon, true},
                                               {&r.phys.polygon, &x.polygon, false}})
                    for (const Vec3& p : *P) {
                        double m = kInf;
                        for (const Vec3& q : *Q) m = std::min(m, toX ? (T.toX(p) - q).norm() : (p - T.toX(q)).norm());
                        polyErr = std::max(polyErr, m);
                    }
            }
            // the +g side and g of the tets' xi strip cuts (gFront: the front pulled back, also outside the tet) against
            // O's cut at random points of every tet that the plane crosses
            std::mt19937_64 rng(mixSeed(gp.seed, 11000 + 10*ri + ki));
            std::uniform_real_distribution<double> Ud(0.0, 1.0);
            int nCmp = 0, badSide = 0;
            double gErr[2] = {0, 0};
            for (size_t t = 0; t < rt.size(); t++) {
                if (!rt[t].plane) continue;
                for (int k = 0; k < 1000; k++) {
                    double l[4], sl = 0;
                    for (double& v : l) sl += (v = -std::log(Ud(rng)));
                    const Vec3 xi(l[1]/sl, l[2]/sl, l[3]/sl);
                    const Vec3 X = R.T[t].toX(xi);
                    const double hX = rc.n.dot(X) - rc.d;
                    if (std::fabs(hX) < 1e-9) continue;
                    nCmp++;
                    badSide += (rt[t].gFront.side(xi) > 0) != (hX > 0) || (rc.cutO.side(X) > 0) != (hX > 0);
                    for (int m = 0; m < 2; m++) {
                        const GType gt = m ? GType::SqrtFront : GType::One;
                        gErr[m] = std::max(gErr[m], std::fabs(rt[t].gFront.g(xi, gt) - rc.cutO.g(X, gt)));
                    }
                }
            }
            const bool uOk = std::fabs(areaSum - rc.cutO.area()) < 1e-12 && outside == 0 && polyErr < 1e-12 && frontOk;
            const bool sOk = badSide == 0 && gErr[0] < 1e-12 && (R.id != 1 || gErr[1] < 1e-12);
            ok = ok && uOk && sOk && nCut > 0;
            out("    %d of 8 tets cut. Union: the tets' physical crack areas sum to %.6f, O's crack %.6f (difference %.1e); %d tet\n"
                "    polygon vertices outside O's polygon; the xi polygons mapped to X match the physical ones to %.1e; fronts at\n"
                "    u.X = s to %.1e -> %s\n", nCut, areaSum, rc.cutO.area(), areaSum - rc.cutO.area(), outside, polyErr, frontErr,
                uOk ? "OK" : "FAILED");
            out("    at %d random points of the tets the plane crosses: %d disagree on the +g side (the tet's xi cut, O's cut, n.X - d);\n"
                "    max |g_tet(xi) - g_O(X)| = %.1e for g = 1, %.1e for the ramp%s -> %s\n", nCmp, badSide, gErr[0], gErr[1],
                R.id == 1 ? "" : " (ring 2: r measured in each tet's xi is not the physical r, so the ramp's crack data would not "
                                  "match between neighbours on their shared faces; the ramp is not run in ring 2)",
                sOk ? "OK" : "FAILED");
        }
    out("\ngeometry assertions of check 4: %s\n", ok ? "all OK" : "SOME FAILED");
    return ok;
}

// The numbers behind check 4's tables and findings, per ring, crack, g mode, level and glue energy
struct Check4Case {
    int ring = 1;
    std::string crack, g, lv, en;
    bool through = false;
    double rms = 0, nc = 0, max = 0, noise = 0, intRms = 0, intNc = 0, lineRms = 0, lineNc = 0, gRms = 0, gNc = 0, gNoise = 0;
    double Ed = 0, Er = 0, seDE = 0, k3Max = 0, k3RefMax = 0, vErr = 0, vSig = 0, dcXi = -1;
    double rmsRefV = -1, lineRefV = -1, EdRefV = 0;   // the decomposition with the reference's vertex values (physical rows)
    int vWorst = 0;                                    // the vertex with the largest abs(decomposed - reference)
    double EdLine = 0, ErLine = 0, seLine = 0;         // the energy within kNearRing of the crack-face lines (physical rows)
};
struct Check4Out {
    MdTable main, vert, vdiag, k3, san, sx, run, tet, z, zr, q;
    std::vector<Check4Case> cases;
    std::string sMd;   // the glue systems of Level A and C at g = 1 (physical energy), per configuration
    int nZ = 0, nZagree = 0;
    double sxRel = 0, cxMax = 0;          // ring 1: max |S_phys - S_xi| / max |S|, max |c_phys - c_xi|
    double sanPz = 0, sanCz = 0;          // the largest |pooled z| of P_k = lambda_k and (ring 1) of H0 (Level C) = J
    double sanPfrac = 1, sanCfrac = 1;    // and the smallest fraction of |z| < 3
    double diagC = 0, diagCmax = 0;       // ring 2: the largest RMS and max of H0 (Level C) - J (xi- vs X-harmonic), / 2 g_max
    double tRecCut = 0, tRecUncut = 0, lRecCut = 0, lRecUncut = 0, tC = 0, lC = 0, glueUs = 0, glueJackMs = 0;
    int nRecCut = 0, nRecUncut = 0, nC = 0, nGlue = 0;
    double refSec = 0, refPilot = 0, zSec = 0;
    long long refPts = 0;
    bool refCached = false;
};

bool runRingConfig(Report& out, const GateParams& gp, const ConvexDomain& O, const Ring& R, const RingCrack& rc, int ki,
                   const std::string& cacheDir, Check4Out& co)
{
    const int ri = R.id - 1;
    const std::string rName = "ring " + std::to_string(R.id), cfg = rName + ", " + rc.name;
    const int nModes = R.id == 1 && !rc.through ? 2 : 1;   // the ramp: ring 1 only (r and L agree in every frame), not K3
    const char* gShort[2] = {"g = 1", "ramp"};
    const char* lvShort[2] = {"A (J = 0)", "C (J = front sol.)"};
    const char* enShort[2] = {"physical", "xi"};
    out("\n--- check 4, %s (%s; crack %s) ---\n", cfg.c_str(), R.name.c_str(), rc.desc.c_str());

    // the per-tet cracks and point sets, with their quadrature checks
    std::vector<RingTetCrack> rt;
    std::vector<PointSet> ps;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        rt.push_back(ringTetCrack(T, rc));
        // at least 1 point per leaf (checks 2 and 3: 4). A tet that the crack crosses has several crack-face lines and up to
        // ~800 depth-4 leaves around them; with 4 per leaf the minimum alone exceeds the target, every leaf gets 4 points
        // (a 64x density contrast instead of 8x) and the far field is starved (1.3% error on the quadratics in a first run)
        ps.push_back(ringPoints(T, rt.back(), rc.cutO, gp.nRingPoints, gp.eps, 1));
        const PointSet& p = ps.back();
        const std::string label = cfg + ", " + T.name;
        out("point set %s (%s): %d leaves (depth 2/3/4: %d/%d/%d), %zu generated (%d/%d/%d at depth 2/3/4), %zu kept\n",
            label.c_str(), rt.back().cut() ? "cut" : "not cut", p.nLeaves, p.leavesByDepth[2], p.leavesByDepth[3],
            p.leavesByDepth[4], p.allXi.size(), p.pointsByDepth[2], p.pointsByDepth[3], p.pointsByDepth[4], p.size());
        quadratureCheckWith(out, T, rt.back().cut() || p.nDropped > 0, p, gp.nTube, &co.q, label,
                            [&](const uint64_t b[3], double o[kNTest]) { ringTubeIntegrals(T, rc.cutO, gp.eps, gp.nTube, b, o); });
    }
    // regions (physical X): within kNearRing of an interior face, of a crack-face line (any tet's), of the crack
    std::vector<std::array<Vec3, 2>> lines;
    for (const RingTetCrack& r : rt) lines.insert(lines.end(), r.faceLines.begin(), r.faceLines.end());
    auto nearInterior = [&](const Vec3& X) {
        for (const auto& f : R.innerFaces)
            if (pointTriangleDist(X, R.V[size_t(f[0])], R.V[size_t(f[1])], R.V[size_t(f[2])]) < kNearRing) return true;
        return false;
    };
    auto nearLine = [&](const Vec3& X) {
        for (const auto& s : lines)
            if (pointSegmentDist(X, s[0], s[1]) < kNearRing) return true;
        return false;
    };
    std::vector<std::vector<char>> inInt(kRingNT), inLine(kRingNT);
    for (int t = 0; t < kRingNT; t++)
        for (const Vec3& X : ps[size_t(t)].X) {
            inInt[size_t(t)].push_back(char(nearInterior(X)));
            inLine[size_t(t)].push_back(char(nearLine(X)));
        }

    // the reference and its Zombie cross-check
    const std::vector<Vec3> probes = ringProbes(R);
    const RefRun ref = runRingReference(out, gp, O, R, rc, "check4_ring" + std::to_string(R.id) + "_" + rc.name,
                                        14000 + 10*uint64_t(ri) + uint64_t(ki), ps, probes, cacheDir);
    if (!ref.ok) return false;
    co.refSec = ref.stats.seconds;
    co.refPilot = ref.pilotExtrapolated;
    co.refCached = ref.cached;
    co.refPts = (long long)ref.X.size();
    std::vector<size_t> off(kRingNT, 0);
    for (int t = 1; t < kRingNT; t++) off[size_t(t)] = off[size_t(t - 1)] + ps[size_t(t - 1)].size();
    std::vector<std::array<char, kRefZRegions>> inRegion;
    for (int t = 0; t < kRingNT; t++)
        for (size_t p = 0; p < ps[size_t(t)].size(); p++)
            inRegion.push_back({char(rt[size_t(t)].cut()), char(!rt[size_t(t)].cut()), 0, inInt[size_t(t)][p],
                                char(rc.cutO.distToCut(ps[size_t(t)].X[p]) < kNearRing)});
    for (const Vec3& X : probes) inRegion.push_back({0, 0, 1, char(nearInterior(X)), char(rc.cutO.distToCut(X) < kNearRing)});
    const RefZombie zr = zombieRingReference(gp, O, rc, ref, inRegion, nModes, 15000 + 1000*uint64_t(ri) + 100*uint64_t(ki));
    reportRefZombie(out, gp, zr, cfg, co.z, co.zr, kRingZRegionName, "O", nModes, true);
    co.zSec = zr.seconds;
    bool zOk = true;
    for (int m = 0; m < nModes; m++) {
        const bool a = agreeHExact(zr, m) && zr.agreeG(m);
        co.nZ++;
        co.nZagree += a;
        zOk = zOk && a;
    }

    // the 16 per-tet runs (production setting: all four faces Dirichlet), each in the tet's own xi frame
    const FaceFn J[2] = {ringJ(rc, GType::One), ringJ(rc, GType::SqrtFront)};
    const std::vector<FieldSpec> fieldsC = ringLevelCFields(J);
    std::vector<OwnRun> rec, hc;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        const RingTetCrack& r = rt[size_t(t)];
        const Cut* gf = r.cut() ? &r.gFront : nullptr;
        const uint64_t tag = 12000 + 1000*uint64_t(ri) + 100*uint64_t(ki) + 10*uint64_t(t);
        rec.push_back(runOwn(T, r.xc, allDirichlet(), ringRecipeFields(T), ps[size_t(t)].xi, gp.nWalks, gp.nBatches,
                             mixSeed(gp.seed, tag), true, gf));
        hc.push_back(runOwn(T, r.xc, allDirichlet(), fieldsC, ps[size_t(t)].xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, tag + 1),
                            true, gf));
        out("%s %s (%s): recipe run %s\n          Level-C run %s\n", cfg.c_str(), T.name, r.cut() ? "cut" : "not cut",
            rec.back().stats.summary().c_str(), hc.back().stats.summary().c_str());
        const SolveStats &a = rec.back().stats, &b = hc.back().stats;
        (r.cut() ? co.tRecCut : co.tRecUncut) += a.seconds;
        (r.cut() ? co.lRecCut : co.lRecUncut) += a.meanWalkLength();
        (r.cut() ? co.nRecCut : co.nRecUncut)++;
        co.tC += b.seconds;
        co.lC += b.meanWalkLength();
        co.nC++;
        co.tet.rows.push_back({rName, rc.name, T.name, r.cut() ? "yes" : "no", r.plane ? fmtd("%+.3f", r.tRaw) : "-",
                               std::to_string(ps[size_t(t)].size()) + " (" + std::to_string(ps[size_t(t)].nDropped) + ")",
                               fmtd("%.2f", a.seconds), fmtd("%.1f", a.meanWalkLength()),
                               fmtd("%.3f", 1e6*a.seconds/double(std::max<long long>(1, a.walksRequested))),
                               fmtd("%.2f", b.seconds), fmtd("%.1f", b.meanWalkLength()),
                               std::to_string(a.walksCapped + b.walksCapped) + " / " + std::to_string(a.walksEscaped + b.walksEscaped)});
    }

    // sanity of the production recipe: where the crack does not reach a tet, H0 at Level A = 0 (every walk ends on a face
    // with value 0) and P_k = lambda_k (check 0; lambda_k is linear, so harmonic in xi as in X). In ring 1, in every tet
    // H0 at Level C = J: J is harmonic in X in the tet minus its crack and equal to the crack data there, and with J_T
    // orthogonal Delta_xi is Delta_X. In ring 2 the tet's problem is Delta_xi H = 0 with the same data, and J (harmonic in
    // X) is not harmonic in xi when J_T is not orthogonal: H0 at Level C is the xi-harmonic extension of J's face data and
    // differs from J by design; that difference is shown as a diagnostic, not tested. Per point z with the points within
    // kCompareBand of a face (in xi) left out, as in check 0, and the quadrature-weighted RMS of the difference against
    // its noise (the RMS s.e.), over all points; / 2 g_max.
    ZStats zP, zPk[4], zC[2];   // zPk: per corner (sum_k P_k = 1 walk by walk here, so z pooled over all k is 0 exactly)
    double maxH0A = 0, pE2 = 0, pV = 0, pW = 0, pMax = 0, cE2[2] = {0, 0}, cV[2] = {0, 0}, cW[2] = {0, 0}, cMax[2] = {0, 0};
    std::string uncut;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        const PointSet& p = ps[size_t(t)];
        const bool cut = rt[size_t(t)].cut();
        if (!cut) uncut += (uncut.empty() ? "" : " ") + std::string(T.name);
        for (size_t i = 0; i < p.size(); i++) {
            const bool band = dirichletDist(allDirichlet(), p.xi[i]) < kCompareBand;
            double H, var;
            long long n;
            if (!cut) {
                for (int k = 0; k < 2; k++) {
                    comboH(rec[size_t(t)].sums, i, unitCoef(6, k), H, var, n);
                    maxH0A = std::max(maxH0A, std::fabs(H));
                }
                for (int k = 0; k < 4; k++) {
                    comboH(rec[size_t(t)].sums, i, unitCoef(6, kRingP0 + k), H, var, n);
                    const double e = H - lambdaOf(T, k, p.X[i]);
                    if (!band) {
                        zP.add(i, e, var);
                        zPk[k].add(i, e, var);
                    }
                    pE2 += p.w[i]*e*e;
                    pV += p.w[i]*var;
                    pW += p.w[i];
                    pMax = std::max(pMax, std::fabs(e));
                }
            }
            for (int m = 0; m < nModes; m++) {
                comboH(hc[size_t(t)].sums, i, unitCoef(2, m), H, var, n);
                const double e = H - J[m](p.X[i]);
                if (!band) zC[m].add(i, e, var);
                cE2[m] += p.w[i]*e*e;
                cV[m] += p.w[i]*var;
                cW[m] += p.w[i];
                cMax[m] = std::max(cMax[m], std::fabs(e));
            }
        }
    }
    {
        const double pr = pW > 0 ? std::sqrt(pE2/pW) : 0.0, pn = pW > 0 ? std::sqrt(pV/pW) : 0.0;
        double pz = 0;   // the largest |pooled z| over the four corners
        for (const ZStats& z : zPk) pz = std::max(pz, std::fabs(z.pooled()));
        if (pW > 0)
            out("sanity, the tets the crack does not reach (%s): max |H0 (Level A)| = %.1e (0 exactly: every walk ends on a face\n"
                "  with 0); P_k - lambda_k: RMS %.1e against its noise %.1e, max %.1e; per point (farther than %.2f from the faces)\n"
                "  |z| < 3 at %.2f%% of %lld; pooled z per corner k = 0..3: %+.2f %+.2f %+.2f %+.2f\n", uncut.c_str(), maxH0A, pr, pn,
                pMax, kCompareBand, 100*zP.frac(), zP.n, zPk[0].pooled(), zPk[1].pooled(), zPk[2].pooled(), zPk[3].pooled());
        else
            out("sanity: the crack reaches every tet (no P_k = lambda_k test here)\n");
        std::string cs;
        for (int m = 0; m < nModes; m++) {
            const double cr = std::sqrt(cE2[m]/cW[m])/kTwoGMax, cn = std::sqrt(cV[m]/cW[m])/kTwoGMax;
            out("%s, H0 (Level C) - J over all 8 tets, %s: RMS %.1e against its noise %.1e (/ 2 g_max), max %.1e; per point |z| < 3\n"
                "  at %.2f%% of %lld, pooled z %+.2f%s\n", R.id == 1 ? "sanity" : "diagnostic (not a test)", gShort[m], cr, cn,
                cMax[m]/kTwoGMax, 100*zC[m].frac(), zC[m].n, zC[m].pooled(),
                R.id == 1 ? "" : "\n  (ring 2: H0 is xi-harmonic in each tet's frame, J is harmonic in X; they differ by design)");
            cs += std::string(m ? "; " : "") + gShort[m] + ": " + fmtd("%.1e", cr) + " (" + fmtd("%.1e", cn) + "), " +
                  fmtd("%.4f", zC[m].frac()) + ", " + fmtd("%+.2f", zC[m].pooled()) + (R.id == 1 ? "" : " (diagnostic: xi- vs X-harmonic)");
            if (R.id == 1) {
                co.sanCz = std::max(co.sanCz, std::fabs(zC[m].pooled()));
                co.sanCfrac = std::min(co.sanCfrac, zC[m].frac());
            } else {
                co.diagC = std::max(co.diagC, cr);
                co.diagCmax = std::max(co.diagCmax, cMax[m]/kTwoGMax);
            }
        }
        if (pW > 0) {
            co.sanPz = std::max(co.sanPz, pz);
            co.sanPfrac = std::min(co.sanPfrac, zP.frac());
        }
        co.san.rows.push_back({rName, rc.name, uncut.empty() ? "none" : uncut, pW > 0 ? fmtd("%.1e", maxH0A) : "-",
                               pW > 0 ? fmtd("%.1e", pr) + " (" + fmtd("%.1e", pn) + "), " + fmtd("%.1e", pMax) + "; " +
                                            fmtd("%.4f", zP.frac()) + ", " + fmtd("%.2f", pz)
                                      : "-",
                               cs});
    }

    // the Grams of every tet, physical and xi, and the reference's energy over the tets' points
    std::vector<Gram> gP, gX;
    for (int t = 0; t < kRingNT; t++) {
        gP.push_back(gramRing(rec[size_t(t)], hc[size_t(t)], ps[size_t(t)].w, R.T[size_t(t)], true));
        gX.push_back(gramRing(rec[size_t(t)], hc[size_t(t)], ps[size_t(t)].w, R.T[size_t(t)], false));
    }
    double dSx = 0, sAbs = 0;
    for (int t = 0; t < kRingNT; t++) {
        dSx = std::max(dSx, (gP[size_t(t)].G - gX[size_t(t)].G).cwiseAbs().maxCoeff());
        sAbs = std::max(sAbs, gP[size_t(t)].G.cwiseAbs().maxCoeff());
    }
    if (R.id == 1) {
        out("ring 1: physical and xi Grams of the same estimates: max |G_phys - G_xi| = %.1e (max |G| = %.3f): equal up to "
            "round-off, as J_T is orthogonal\n", dSx, sAbs);
        co.sxRel = std::max(co.sxRel, dSx/sAbs);
    }
    std::vector<double> wAll(ref.X.size(), 0.0);
    for (int t = 0; t < kRingNT; t++)
        for (size_t p = 0; p < ps[size_t(t)].size(); p++) wAll[off[size_t(t)] + p] = ps[size_t(t)].w[p];
    const std::vector<double> unit2[2] = {unitCoef(2, 0), unitCoef(2, 1)};
    const Gram gr = gramFromBatches(wAll, ref.sums.nBatches, 2, [&](size_t p, int b, Vec3* g) {
        for (int q = 0; q < 2; q++) g[q] = comboBatchGrad(ref.sums, p, b, unit2[q]);
    });
    // the same within kNearRing of the crack-face lines (where the energy of a face/crack data mismatch would sit)
    std::vector<Gram> gPL;
    std::vector<double> wLineAll(ref.X.size(), 0.0);
    for (int t = 0; t < kRingNT; t++) {
        std::vector<double> wl(ps[size_t(t)].w);
        for (size_t p = 0; p < wl.size(); p++) {
            if (!inLine[size_t(t)][p]) wl[p] = 0;
            wLineAll[off[size_t(t)] + p] = wl[p];
        }
        gPL.push_back(gramRing(rec[size_t(t)], hc[size_t(t)], wl, R.T[size_t(t)], true));
    }
    const Gram grL = gramFromBatches(wLineAll, ref.sums.nBatches, 2, [&](size_t p, int b, Vec3* g) {
        for (int q = 0; q < 2; q++) g[q] = comboBatchGrad(ref.sums, p, b, unit2[q]);
    });

    // every glue: g mode x level x energy
    double cxMax = 0, evMin = kInf;
    for (int m = 0; m < nModes; m++) {
        std::vector<FieldAt> refT;
        for (int t = 0; t < kRingNT; t++) refT.push_back(refField(ref, m, off[size_t(t)], ps[size_t(t)].size()));
        const FieldAt refP = refField(ref, m, ref.nA, ref.nF);
        for (int L = 0; L < 2; L++) {
            Vec7 cPhys = Vec7::Zero();
            for (int e = 0; e < 2; e++) {
                const int h = L ? kRingH0C + m : m;
                const Glue4 gl = glueRun4(R, e ? gX : gP, h);
                co.glueUs += 1e6*gl.seconds;
                co.glueJackMs += 1e3*gl.secondsJack;
                co.nGlue++;
                if (e == 0) cPhys = gl.c;
                const double dcXi = e ? (gl.c - cPhys).cwiseAbs().maxCoeff() : -1.0;
                if (e && R.id == 1) cxMax = std::max(cxMax, dcXi);
                Eigen::SelfAdjointEigenSolver<Mat7> es(gl.M);
                evMin = std::min(evMin, es.eigenvalues()[0]);
                out("\n%s, %s, Level %s, glue with the %s energy:\n", cfg.c_str(), gShort[m], lvShort[L], enShort[e]);
                out("  7 x 7 glue system M c = r (rows/cols v0..v6; M = sum_T scatter(S_T), r = -sum_T scatter(b_T); smallest\n"
                    "  eigenvalue of M %.4f):\n", es.eigenvalues()[0]);
                std::string sys;
                for (int i = 0; i < kRingNV; i++) {
                    const std::string row = "    [" + matRow(gl.M, i, "%+.4f") + "]  [" + fmtd("%+.5f", gl.rhs[i]) + "]";
                    out("%s\n", row.c_str());
                    sys += row + "\n";
                }
                if (m == 0 && e == 0)
                    co.sMd += cfg + ", g = 1, Level " + lvShort[L] + ", physical energy: M | r\n" + sys;
                out("  glued c = %s\n       +- %s (jackknife over the batches of all 16 runs, the system re-solved per replicate)\n"
                    "  solve %.1f us (assembly included), %d replicate solves %.2f ms\n", fmtVec(gl.c).c_str(),
                    fmtVec(gl.se, "%.5f").c_str(), 1e6*gl.seconds, int(gl.jack.size()), 1e3*gl.secondsJack);
                if (e) out("  max |c_xi - c_phys| = %.2e\n", dcXi);

                // the vertex values: the decomposed field's are its face data there, J(v) + c_v (J = 0 at Level A); the
                // reference's from its probes (c itself; the outer vertices pulled 1e-3 towards the origin)
                std::string vRow[kRingNV];
                double vErr = 0, vSig = 0;
                int vWorst = 0;
                out("  vertex values: decomposed J(v) + c_v (s.e.) / reference at the probe (s.e.):\n   ");
                for (int k = 0; k < kRingNV; k++) {
                    const double vd = (L ? J[m](R.V[size_t(k)]) : 0.0) + gl.c[k], vr = refP.H[size_t(k)];
                    const double sr = std::sqrt(refP.vH[size_t(k)]);
                    if (std::fabs(vd - vr) > vErr) vWorst = k;
                    vErr = std::max(vErr, std::fabs(vd - vr));
                    vSig = std::max(vSig, std::fabs(vd - vr)/std::hypot(gl.se[k], sr));
                    vRow[k] = fmtd("%+.4f", vd) + " (" + fmtd("%.4f", gl.se[k]) + ") / " + fmtd("%+.4f", vr) + " (" + fmtd("%.4f", sr) + ")";
                    out(" v%d %s;", k, vRow[k].c_str());
                }
                out("\n   max |decomposed - reference| = %.4f (%.1f combined s.e.)\n", vErr, vSig);

                // the fields and their errors against the reference, / 2 g_max
                ErrStats eAll, eInt, eLine;
                double k3Max = 0, k3RefMax = 0, Ed = 0, Exi = 0;
                std::vector<double> rT(size_t(gP[0].B), 0.0);
                for (int t = 0; t < kRingNT; t++) {
                    const TetFrame& T = R.T[size_t(t)];
                    const Eigen::Vector4d cT = localValues(T, gl.c);
                    std::vector<Eigen::Vector4d> cJ;
                    for (const Vec7& x : gl.jack) cJ.push_back(localValues(T, x));
                    const FieldAt dec = decField4(T, rec[size_t(t)], hc[size_t(t)], L, m, cT, &cJ);
                    const PointSet& p = ps[size_t(t)];
                    for (size_t i = 0; i < p.size(); i++) {
                        eAll.add(p.w[i], p.X[i], dec, i, refT[size_t(t)], i);
                        if (inInt[size_t(t)][i]) eInt.add(p.w[i], p.X[i], dec, i, refT[size_t(t)], i);
                        if (inLine[size_t(t)][i]) eLine.add(p.w[i], p.X[i], dec, i, refT[size_t(t)], i);
                        if (rc.through) {
                            const double sgn = rc.n.dot(p.X[i]) - rc.d > 0 ? 1.0 : -1.0;
                            k3Max = std::max(k3Max, std::fabs(dec.H[i] - sgn));
                            k3RefMax = std::max(k3RefMax, std::fabs(refT[size_t(t)].H[i] - sgn));
                        }
                    }
                    Ed += tetEnergy4(gP[size_t(t)].G, h, cT);
                    Exi += tetEnergy4(gX[size_t(t)].G, h, cT);
                    for (size_t a = 0; a < rT.size(); a++) rT[a] += tetEnergy4(gP[size_t(t)].jack[a], h, cJ[a]);
                }
                const double seD = jackSe(rT), Er = 0.5*gr.G(m, m), seR = 0.5*gr.se(m, m), seDE = std::hypot(seD, seR);
                out("  errors of the glued decomposition vs the reference, / 2 g_max (nc = noise-corrected RMS):\n"
                    "    region                    points   H RMS (nc)          H max     grad RMS (nc)       grad max  | ref s.e. H RMS/max  dec s.e. H RMS/max\n");
                for (const auto& [name, s] : {std::pair<const char*, const ErrStats*>{"the ring (8 tets)", &eAll},
                                              {"near interior faces (0.1)", &eInt}, {"near crack-face lines (0.1)", &eLine}})
                    out("    %-26s %6lld   %.2e (%.2e)   %.2e  %.2e (%.2e)   %.2e  | %.1e / %.1e     %.1e / %.1e\n", name, s->n,
                        s->rmsH(), s->corrH(), s->maxH/kTwoGMax, s->rmsG(), s->corrG(), s->maxG/kTwoGMax, s->seRefH(),
                        s->maxSeRefH/kTwoGMax, s->seFldH(), s->maxSeFldH/kTwoGMax);
                out("    largest H error at X = %s, largest grad H error at X = %s\n", fmt3(eAll.worstH, 3).c_str(),
                    fmt3(eAll.worstG, 3).c_str());
                out("  energy (1/2) int |grad_X H|^2 over the ring minus the eps-tube: decomposed %.5f +- %.5f (physical, with this\n"
                    "  glue's c), reference %.5f +- %.5f, difference %+.5f +- %.5f%s; this glue's own objective (%s energy) %.5f\n",
                    Ed, seD, Er, seR, Ed - Er, seDE, Er > 1e-12 ? (" (" + fmtd("%+.1f%%", 100*(Ed/Er - 1)) + ")").c_str() : "",
                    enShort[e], e ? Exi : Ed);
                if (rc.through)
                    out("  K3 sanity (exact solution +-1 on the two sides of the crack): max |H_dec - (+-1)| = %.4f over the ring's\n"
                        "  points (the reference: %.1e)\n", k3Max, k3RefMax);
                // diagnostic, no new walks: the same per-tet fields with the reference's own vertex values, c_v = H_ref(probe v)
                // - J(v) (J = 0 at Level A), instead of the glued ones. Its error against the reference separates the glue's
                // choice of vertex values from what the basis H0 + sum_k c_k P_k does with good ones; its energy shows what the
                // energy minimization trades. Once per g mode and level (it does not depend on the glue's energy).
                double rmsRefV = -1, lineRefV = -1, EdRefV = 0, k3RefV = 0, EdL = 0, seL = 0;
                if (e == 0) {
                    std::vector<double> rL(size_t(gP[0].B), 0.0);
                    for (int t = 0; t < kRingNT; t++) {
                        const TetFrame& T = R.T[size_t(t)];
                        EdL += tetEnergy4(gPL[size_t(t)].G, h, localValues(T, gl.c));
                        for (size_t a = 0; a < rL.size(); a++) rL[a] += tetEnergy4(gPL[size_t(t)].jack[a], h, localValues(T, gl.jack[a]));
                    }
                    seL = std::hypot(jackSe(rL), 0.5*grL.se(m, m));
                    out("  energy within %.1f of the crack-face lines: decomposed %.5f, reference %.5f, difference %+.5f +- %.5f (of the\n"
                        "    ring's %+.5f)\n", kNearRing, EdL, 0.5*grL.G(m, m), EdL - 0.5*grL.G(m, m), seL, Ed - Er);
                    Vec7 cRef;
                    for (int k = 0; k < kRingNV; k++) cRef[k] = refP.H[size_t(k)] - (L ? J[m](R.V[size_t(k)]) : 0.0);
                    ErrStats fAll, fLine;
                    for (int t = 0; t < kRingNT; t++) {
                        const TetFrame& T = R.T[size_t(t)];
                        const Eigen::Vector4d cT = localValues(T, cRef);
                        const FieldAt dv = decField4(T, rec[size_t(t)], hc[size_t(t)], L, m, cT, nullptr);
                        const PointSet& p = ps[size_t(t)];
                        for (size_t i = 0; i < p.size(); i++) {
                            fAll.add(p.w[i], p.X[i], dv, i, refT[size_t(t)], i);
                            if (inLine[size_t(t)][i]) fLine.add(p.w[i], p.X[i], dv, i, refT[size_t(t)], i);
                            if (rc.through) k3RefV = std::max(k3RefV, std::fabs(dv.H[i] - (rc.n.dot(p.X[i]) - rc.d > 0 ? 1.0 : -1.0)));
                        }
                        EdRefV += tetEnergy4(gP[size_t(t)].G, h, cT);
                    }
                    rmsRefV = fAll.rmsH();
                    lineRefV = fLine.rmsH();
                    out("  diagnostic, the reference's vertex values (c_v = H_ref(probe) - J(v)) in place of the glued ones: H RMS %.2e\n"
                        "    (nc %.2e), max %.2e; near crack-face lines RMS %.2e; energy %.5f (glued %.5f, reference %.5f)%s\n",
                        fAll.rmsH(), fAll.corrH(), fAll.maxH/kTwoGMax, fLine.rmsH(), EdRefV, Ed, Er,
                        rc.through ? ("; max |H - (+-1)| " + fmtd("%.4f", k3RefV)).c_str() : "");
                    co.vdiag.rows.push_back({rName, rc.name, gShort[m], lvShort[L],
                                             fmtd("%.2e", eAll.rmsH()) + " / " + fmtd("%.2e", eLine.rmsH()) + " / " + fmtd("%.5f", Ed),
                                             fmtd("%.2e", fAll.rmsH()) + " (" + fmtd("%.2e", fAll.corrH()) + ") / " + fmtd("%.2e", fLine.rmsH()) +
                                                 " / " + fmtd("%.5f", EdRefV),
                                             std::fabs(Er) < 1e-12 ? std::string("0 (exact)") : fmtd("%.5f", Er),
                                             fmtd("%+.5f", EdL - 0.5*grL.G(m, m)) + " +- " + fmtd("%.5f", seL) + " / " + fmtd("%+.5f", Ed - Er),
                                             rc.through ? fmtd("%.4f", k3Max) + " / " + fmtd("%.4f", k3RefV) : "-"});
                }

                // tables
                auto rmm = [](const ErrStats& s) {
                    return fmtd("%.2e", s.rmsH()) + " (" + fmtd("%.2e", s.corrH()) + ") / " + fmtd("%.2e", s.maxH/kTwoGMax);
                };
                co.main.rows.push_back({rName, rc.name, gShort[m], lvShort[L], enShort[e], rmm(eAll), rmm(eInt), rmm(eLine),
                                        fmtd("%.2e", eAll.rmsG()) + " (" + fmtd("%.2e", eAll.corrG()) + ")",
                                        fmtd("%.5f", Ed) + " / " + fmtd("%.5f", Er) + " / " + fmtd("%+.5f", Ed - Er) + " +- " + fmtd("%.5f", seDE),
                                        fmtd("%.1e", eAll.seRefH()) + " / " + fmtd("%.1e", eAll.seFldH()),
                                        fmtd("%.4f", es.eigenvalues()[0])});
                std::vector<std::string> vr = {rName, rc.name, gShort[m], lvShort[L], enShort[e]};
                for (int k = 0; k < kRingNV; k++) vr.push_back(vRow[k]);
                vr.push_back(fmtd("%.4f", vErr) + " (" + fmtd("%.1f", vSig) + ")");
                co.vert.rows.push_back(vr);
                if (rc.through)
                    co.k3.rows.push_back({rName, lvShort[L], enShort[e], fmtd("%.4f", k3Max), fmtd("%.2e", eAll.rmsH()),
                                          fmtd("%.1e", k3RefMax), fmtd("%.5f", Ed),
                                          std::fabs(Er) < 1e-12 ? "0 (" + fmtd("%.0e", Er) + ": round-off)" : fmtd("%.5f", Er)});
                Check4Case c;
                c.ring = R.id;
                c.crack = rc.name;
                c.g = gShort[m];
                c.lv = L ? "C" : "A";
                c.en = enShort[e];
                c.through = rc.through;
                c.rms = eAll.rmsH();
                c.nc = eAll.corrH();
                c.max = eAll.maxH/kTwoGMax;
                c.noise = std::hypot(eAll.seRefH(), eAll.seFldH());
                c.intRms = eInt.rmsH();
                c.intNc = eInt.corrH();
                c.lineRms = eLine.rmsH();
                c.lineNc = eLine.corrH();
                c.gRms = eAll.rmsG();
                c.gNc = eAll.corrG();
                c.gNoise = std::hypot(eAll.seRefG(), eAll.seFldG());
                c.Ed = Ed;
                c.Er = Er;
                c.seDE = seDE;
                c.k3Max = k3Max;
                c.k3RefMax = k3RefMax;
                c.vErr = vErr;
                c.vSig = vSig;
                c.dcXi = dcXi;
                c.rmsRefV = rmsRefV;
                c.vWorst = vWorst;
                c.EdLine = EdL;
                c.ErLine = 0.5*grL.G(m, m);
                c.seLine = seL;
                c.lineRefV = lineRefV;
                c.EdRefV = EdRefV;
                co.cases.push_back(c);
            }
        }
    }
    co.cxMax = std::max(co.cxMax, cxMax);
    co.sx.rows.push_back({rName, rc.name, uncut.empty() ? "none" : uncut, fmtd("%.4f", evMin),
                          R.id == 1 ? fmtd("%.1e", dSx/sAbs) + " / " + fmtd("%.1e", cxMax) : "- (not orthogonal)",
                          zOk ? "agrees" : "DISAGREES"});
    return zOk;
}

bool runCheck4(Report& out, const GateParams& gp, const std::string& cacheDir, std::string& md, std::string& geoMd)
{
    out("\n=== Check 4: one nodal ring (the 8 tets around a node, glued at its 7 vertices) ===\n");
    out("The ring of the node c: 8 tets (c, s_x e_x, s_y e_y, s_z e_z), union = the octahedron O. Ring 1: c = 0 (every tet the\n"
        "reference tet up to a rotation); ring 2: c = (0.12, -0.08, 0.10) (8 different tets, as in a real mesh). Each tet runs the\n"
        "production recipe in its own xi frame with ALL FOUR faces Dirichlet: a recipe run (H0 at Level A: crack +-g, faces 0;\n"
        "P_k, k = its 4 corners: crack 0, faces lambda_k) and a Level-C run (H0: crack +-g, faces J, the crack's own 2D front\n"
        "solution); own loop, N = %d in %d batches. Glue: the 7 x 7 system for the vertex values (O's faces Neumann, imposed\n"
        "weakly: the outer vertices are free), with S_T, b_T from the physical or the xi energy of the same estimates. Reference:\n"
        "one own-loop solve on O in X (8 Neumann faces, N = %d), cached in %s. Errors / 2 g_max = 2; energies over the ring\n"
        "minus the eps-tube.\n", gp.nWalks, gp.nBatches, gp.nRef, cacheDir.c_str());
    md += "## Check 4: one nodal ring\n\n";
    md += "The ring of a node c: the 8 tets (c, s_x e_x, s_y e_y, s_z e_z), s in {+-1}^3, ordered so that det J_T > 0; their "
          "union is the octahedron O = {abs(x) + abs(y) + abs(z) <= 1}. Global vertices v0 = c, v1..v6 = +e_x, -e_x, +e_y, -e_y, "
          "+e_z, -e_z. **Ring 1** (congruent): c = 0, every J_T orthogonal (det +1), so each tet's xi problem is the physical one "
          "and any error is decomposition error. **Ring 2** (real mesh): c = (0.12, -0.08, 0.10), the 8 tets differ and none is "
          "congruent to the reference tet. O and the cracks are the same in both rings, so both rings approximate the same "
          "physical problem.\n\n"
          "Every tet runs the production recipe in its own xi frame (the formulation: Delta_xi H = 0 on the reference tet, the "
          "crack pulled back by the tet's affine map) with **all four faces Dirichlet**, own loop, double, N = " +
          std::to_string(gp.nWalks) + " in " + std::to_string(gp.nBatches) + " batches: a recipe run (H0 at Level A: crack "
          "+-g, faces 0, both g modes; P_k for the tet's 4 corners: crack 0, faces lambda_k; all from the same walks) and a "
          "Level-C run (H0: crack +-g, faces J; both g modes from the same walks). **Level A**: J = 0. **Level C**: J = the "
          "crack's own 2D front solution, one function of physical X for all tets: theta/pi (g = 1) and min(1, sqrt(r/L)) "
          "sin(theta/2) (ramp), theta = atan2(b, a), a = u.X - s (negative behind the front), b = n.X - d, r = abs((a, b)). J is "
          "+-g on the crack and, for a planar crack with a straight front, harmonic off the crack (for the ramp while r < L), "
          "so in ring 1 the Level-C H0 is J itself in every tet (tested below; in ring 2 each tet's H0 is the xi-harmonic "
          "extension of J's data, which differs from J slightly, shown as a diagnostic): J carries the crack's whole jump, and "
          "what the glue adds (sum c_k P_k) mainly corrects J's mismatch with the Neumann outer boundary. That makes Level C much "
          "easier here than a curved crack or front would be.\n\n"
          "Glue: the 7 x 7 system M c = r for the vertex values (written out in `glueCheck4`): M(T.global[i], T.global[j]) += "
          "S_T(i, j), r(T.global[i]) -= b_T(i) over the 8 tets, each S_T symmetrized. O's faces are Neumann in the true problem "
          "and imposed weakly: the 6 outer vertex values are free unknowns (the natural boundary condition). S_T and b_T come "
          "from the same gradient estimates in two ways: **physical** (grad_X = J^-T grad_xi over the physical tet) and **xi** "
          "(grad_xi . grad_xi over the reference tet). In ring 1 they are equal up to round-off; in ring 2 they differ: the "
          "physical glue is the Ritz approximation of the physical ring problem in the per-tet basis, the xi glue minimizes the "
          "sum of the tets' xi energies. Both are compared with the same physical reference. The ramp is run in ring 1 only: in "
          "ring 2 the ramp's r measured in each tet's xi differs between neighbours, so its crack data would not match on their "
          "shared faces (the geometry printout gives the size). K3 (cut through) is run with g = 1 only.\n\n"
          "Reference: one own-loop solve on O in X (all 8 faces Neumann, mirror ball across them, the crack = the plane's section "
          "of O clipped at the front, L = sqrt(2)), N = " + std::to_string(gp.nRef) + " in " + std::to_string(gp.nRefBatches) +
          " batches, at every point of the 8 tets and at 7 vertex probes (c, and the outer vertices pulled 1e-3 towards the "
          "origin), cached in gate_cache/, cross-checked by Zombie on O. Point sets: per tet about " +
          std::to_string(gp.nRingPoints) + " points, refined around the front AND the lines where the crack meets the tet's "
          "faces, the points within eps = 0.01 sqrt(2) (physical) of the crack dropped; the same points for the decomposition "
          "and the reference.\n\n"
          "All errors and noise levels are divided by **2 g_max = 2**. RMS is quadrature-weighted, max over the points; \"nc\" "
          "= noise-corrected RMS. Regions: the ring (all points), near interior faces = within 0.1 of one of the 12 shared "
          "faces, near crack-face lines = within 0.1 of a segment where the crack meets a tet face (interior or on O). The "
          "decomposed noise is the delete-one-batch jackknife over the runs' batches (batch a of all 16 runs is one unit), with "
          "the glued c re-solved per replicate. Energies: (1/2) int abs(grad_X H)^2 over the ring minus the eps-tube, the "
          "decomposed one always physical (with the c of the glue in that row), split-sample, jackknife s.e.\n\n";
    const size_t mdIntro = md.size();   // the findings go here

    const ConvexDomain O = octahedron();
    const std::vector<Ring> rings = {makeRing(1), makeRing(2)};
    const std::vector<RingCrack> cracks = {ringCrack(0), ringCrack(1), ringCrack(2)};
    const size_t geo0 = out.text.size();
    const bool geoOk = printRingGeometry(out, gp, O, rings, cracks);
    geoMd = "## Appendix E: geometry of check 4\n\nPrinted by the run before its first solve: the rings, the cracks on O, and per ring "
            "and crack every tet's crack.\n\n" + codeBlock(out.text.substr(geo0));
    if (!geoOk) {
        out("geometry assertions FAILED: check 4 not run\n");
        md += "**Check 4: FAILED** (a geometry assertion; see Appendix E)\n";
        md.insert(mdIntro, findingsBlock("- Check 4 (one nodal ring): FAILED, a geometry assertion (Appendix E); nothing was solved.\n"));
        return false;
    }

    Check4Out co;
    co.main.head = {"ring", "crack", "g mode", "level", "glue energy", "H: RMS (nc) / max, the ring",
                    "near interior faces: RMS (nc) / max", "near crack-face lines: RMS (nc) / max", "grad H RMS (nc)",
                    "energy: decomposed / reference / difference +- s.e.", "noise: ref / dec s.e. H RMS", "M min eigenvalue"};
    co.vert.head = {"ring", "crack", "g mode", "level", "glue energy"};
    for (int k = 0; k < kRingNV; k++) co.vert.head.push_back(std::string(kRingVertName[k]) + ": dec (s.e.) / ref (s.e.)");
    co.vert.head.push_back("max abs(dec - ref) (combined s.e.)");
    co.vdiag.head = {"ring", "crack", "g mode", "level", "glued (physical): H RMS / near crack-face lines / energy",
                     "the reference's vertex values: H RMS (nc) / near crack-face lines / energy", "energy reference",
                     "energy excess (dec - ref) within 0.1 of the crack-face lines +- s.e. / over the ring",
                     "K3: max abs(H - (+-1)), glued / reference's vertex values"};
    co.k3.head = {"ring", "level", "glue energy", "max abs(H_dec - (+-1))", "H RMS / 2 g_max", "reference: max abs(H_ref - (+-1))",
                  "energy decomposed", "energy reference"};
    co.san.head = {"ring", "crack", "tets the crack does not reach", "there: max abs(H0, Level A)",
                   "there: P_k - lambda_k: RMS (noise), max; frac abs(z)<3, largest abs(pooled z) of the 4 corners",
                   "all tets: H0 (Level C) - J, / 2 g_max: RMS (noise), frac abs(z)<3, pooled z"};
    co.sx.head = {"ring", "crack", "tets not cut", "smallest eigenvalue of M (all glues)",
                  "ring 1: max abs(S_phys - S_xi)/max abs(S) / max abs(c_phys - c_xi)", "Zombie vs reference"};
    co.run.head = {"ring", "crack", "recipe run per tet: s, mean walk length (cut / not cut)", "Level-C run per tet: s, walk length",
                   "glue: one solve (us) / 64 replicate solves (ms)", "reference on O: s (points)", "Zombie cross-check: s"};
    co.tet.head = {"ring", "crack", "tet", "cut", "t", "points (dropped)", "recipe run (s)", "walk length", "us per walk (wall)",
                   "Level-C run (s)", "walk length", "capped / escaped"};
    co.z.head = {"configuration", "g mode", "H: frac abs(z)<3, pooled z (mean dH +- s.e.)",
                 "grad H (points > 0.05 from O's faces): pooled z per component (frac abs(z)<3)", "verdict"};
    co.zr.head = {"configuration", "g mode", "region", "points", "H pooled z", "mean z^2 (s.e.)",
                  "RMS discrepancy beyond noise / 2 g_max (2-s.e. bound)"};
    co.q.head = {"point set", "kept points", "dropped (eps)", "max rel err, all points vs exact over T",
                 "max rel err, kept vs T minus tube", "tube volume"};
    bool ok = true;
    for (size_t ri = 0; ri < rings.size(); ri++)
        for (size_t ki = 0; ki < cracks.size(); ki++) {
            Check4Out c1;   // this configuration's rows and numbers, merged below
            const bool cOk = runRingConfig(out, gp, O, rings[ri], cracks[ki], int(ki), cacheDir, c1);
            ok = ok && cOk;
            for (auto [dst, src] : {std::pair<MdTable*, MdTable*>{&co.main, &c1.main}, {&co.vert, &c1.vert}, {&co.vdiag, &c1.vdiag}, {&co.k3, &c1.k3},
                                    {&co.san, &c1.san}, {&co.sx, &c1.sx}, {&co.tet, &c1.tet}, {&co.z, &c1.z}, {&co.zr, &c1.zr},
                                    {&co.q, &c1.q}})
                dst->rows.insert(dst->rows.end(), src->rows.begin(), src->rows.end());
            co.cases.insert(co.cases.end(), c1.cases.begin(), c1.cases.end());
            co.sMd += c1.sMd;
            co.nZ += c1.nZ;
            co.nZagree += c1.nZagree;
            co.sxRel = std::max(co.sxRel, c1.sxRel);
            co.cxMax = std::max(co.cxMax, c1.cxMax);
            co.sanPz = std::max(co.sanPz, c1.sanPz);
            co.sanCz = std::max(co.sanCz, c1.sanCz);
            co.sanPfrac = std::min(co.sanPfrac, c1.sanPfrac);
            co.sanCfrac = std::min(co.sanCfrac, c1.sanCfrac);
            co.diagC = std::max(co.diagC, c1.diagC);
            co.diagCmax = std::max(co.diagCmax, c1.diagCmax);
            const std::string rName = "ring " + std::to_string(rings[ri].id);
            auto avg = [](double s, int n) { return n ? s/double(n) : 0.0; };
            co.run.rows.push_back({rName, cracks[ki].name,
                                   fmtd("%.2f", avg(c1.tRecCut, c1.nRecCut)) + ", " + fmtd("%.1f", avg(c1.lRecCut, c1.nRecCut)) + " / " +
                                       (c1.nRecUncut ? fmtd("%.2f", avg(c1.tRecUncut, c1.nRecUncut)) + ", " +
                                                           fmtd("%.1f", avg(c1.lRecUncut, c1.nRecUncut))
                                                     : std::string("-")),
                                   fmtd("%.2f", avg(c1.tC, c1.nC)) + ", " + fmtd("%.1f", avg(c1.lC, c1.nC)),
                                   fmtd("%.1f", avg(c1.glueUs, c1.nGlue)) + " / " + fmtd("%.2f", avg(c1.glueJackMs, c1.nGlue)),
                                   fmtd("%.0f", c1.refSec) + " (" + std::to_string(c1.refPts) + ")" +
                                       (c1.refCached ? ", loaded from gate_cache/ (wall clock of the original run)"
                                                     : ", pilot extrapolated " + fmtd("%.0f", c1.refPilot)),
                                   fmtd("%.0f", c1.zSec)});
        }
    out("\nCHECK 4: %s (no pass threshold: the numbers above are the result; the verdict covers the geometry assertions and\n"
        "the Zombie cross-check of the reference)\n", ok ? "done" : "FAILED");

    md += "### Main table (the ring)\n\n" + co.main.str() +
          "\n### Vertex values\n\nThe decomposed field's value at a vertex is its face data there, J(v) + c_v (J = 0 at Level A), "
          "with the jackknife s.e. of c_v. The reference's is its H at the probe: c itself, and the outer vertices pulled 1e-3 "
          "towards the origin (a walk cannot start on O's corner), so there it differs from the vertex value by about 1e-3 "
          "abs(grad H). Last column: the largest abs(decomposed - reference) over the 7 vertices, in combined s.e. in brackets.\n\n" +
          co.vert.str() +
          "\n### Diagnostic: the reference's vertex values in place of the glued ones\n\nThe same per-tet fields (no new walks), "
          "H = H0 + sum_k c_k P_k, with c_v = H_ref(probe v) - J(v) (J = 0 at Level A): the vertex values a perfect glue would "
          "give if the field were right at the vertices. If this is much better than the glued field, the glue picks poor "
          "vertex values; if not, the per-tet basis itself is the limit. Its energy shows what the energy minimization trades "
          "(the physical glue minimizes this same quadratic of the Grams, so its energy is the lower one by construction). The "
          "energy column splits the glued field's energy excess: the part within 0.1 of the crack-face lines against the whole "
          "ring's. H errors / 2 g_max.\n\n" + co.vdiag.str() +
          "\n### K3: the crack cut through O\n\nWith g = 1 and O's faces Neumann the exact solution is +1 on the +g side of the "
          "crack and -1 on the other (the crack separates O). The reference reproduces it exactly: every walk ends on the crack on "
          "its own side, so H_ref = +-1 with zero variance and zero energy. max abs(H_dec - (+-1)) is over the ring's sample points "
          "(not divided by 2 g_max); Level C's J has its front line 0.3 beyond O, so J is not +-1 away from the crack and the glue "
          "must correct it.\n\n" + co.k3.str() +
          "\n### Sanity of the production recipe\n\nWhere the crack does not reach a tet, the recipe must give H0 (Level A) = 0 "
          "exactly (every walk ends on a face with 0) and P_k = lambda_k (check 0; lambda_k is linear, so harmonic in xi and in "
          "X). In ring 1 the Level-C H0 must be J in every tet: J is harmonic in X off the crack and equals the crack data on it, "
          "and with J_T orthogonal Delta_xi = Delta_X. In ring 2 each tet solves Delta_xi H = 0 in its own frame, the formulation, "
          "and J is harmonic in X, not in xi, when J_T is not orthogonal: there the Level-C H0 is the xi-harmonic extension of J's "
          "data and differs from J by design. Its size is shown as a diagnostic (marked), not tested. Per point z over the points "
          "farther than 0.01 (xi) from the faces, as in check 0; RMS over all points, / 2 g_max for H0 - J.\n\n" +
          co.san.str() +
          "\n### Glue systems\n\nThe smallest eigenvalue of M over every glue of the configuration (M is symmetric positive "
          "definite once the crack cuts a tet). Ring 1: the physical and xi Grams of the same estimates, and the glued c, agree up "
          "to round-off (J_T orthogonal). Below: M and r of Level A and Level C at g = 1, physical energy.\n\n" + co.sx.str() +
          "\n" + codeBlock(co.sMd) +
          "\n### Reference cross-check (Zombie, float)\n\nZombie on O (its own octahedron mesh, 8 reflecting triangles, the crack "
          "absorbing) at the reference's points: N = " + std::to_string(gp.nZombie) + " for H at every point, and for grad H " +
          std::to_string(gp.nZombieBatches) + " solves of " + std::to_string(gp.nZombie/gp.nZombieBatches) + " walks at the "
          "points farther than " + fmtd("%.2f", kRefZGradBand) + " from O's faces; g = 1, and the ramp where it is run (ring 1, "
          "K1 and K2). Agree: H abs(z) < 3 at >= 99% and abs(pooled z) < 3; grad abs(pooled z) < 3 per component. The second "
          "table is a diagnostic per region (see check 2's description of the columns).\n\n" + co.z.str() + "\n" + co.zr.str() +
          "\n### Runs and timing\n\nPer configuration: the mean wall clock and mean walk length of the per-tet runs (all four faces "
          "Dirichlet: the production setting; the recipe run scores 6 data, the Level-C run 2, on the same walks each), the glue "
          "solve (assembly and the 7 x 7 LDLT; the 64 jackknife re-solves), the reference on O and its Zombie cross-check. "
          "16 hardware threads.\n\n" + co.run.str() + "\nPer tet:\n\n" + co.tet.str() +
          "\n### Quadrature of the point sets\n\n" + co.q.str() + "\n**Check 4: " + (ok ? "done" : "FAILED") +
          "** (no pass threshold; the verdict covers the geometry assertions and the reference's Zombie cross-check)\n";

    // findings: generated from the cases above
    auto find = [&](int ring, const std::string& crack, const std::string& g, const std::string& lv, const std::string& en) -> const Check4Case* {
        for (const Check4Case& c : co.cases)
            if (c.ring == ring && c.crack == crack && c.g == g && c.lv == lv && c.en == en) return &c;
        return nullptr;
    };
    std::string f = std::string("- Check 4 (one nodal ring: the 8 tets around a node glued at its 7 vertices, every tet with all four "
                                "faces Dirichlet, O's faces free): ") + (ok ? "done" : "FAILED") + ". No pass threshold; the "
                    "geometry assertions hold, and Zombie agrees with the reference in " + std::to_string(co.nZagree) + " of " +
                    std::to_string(co.nZ) + " cases.\n";
    for (int ring = 1; ring <= 2; ring++) {
        f += ring == 1 ? "- Ring 1 (congruent, c = 0: every tet the reference tet up to a rotation, so any error is decomposition "
                         "error). The glued decomposition against the reference on O, / 2 g_max, physical-energy glue (the xi glue "
                         "is identical here); noise = the combined RMS s.e. of H; E = (1/2) int abs(grad H)^2, decomposed / reference "
                         "- 1 (s.e.); last column: H RMS with the reference's own vertex values in place of the glued ones:\n\n"
                       : "- Ring 2 (real mesh, c = (0.12, -0.08, 0.10)), g = 1, both glues (physical / xi energy), the same columns:\n\n";
        std::vector<std::vector<std::string>> rows = {{"crack", "g", "level", "glue", "H RMS (nc)", "error/noise", "near int. faces",
                                                       "near crack-face lines", "H max", "E dec / E ref - 1 (s.e.)", "ref. vertex values"}};
        double ratioLo = kInf, ratioHi = 0, k3[2] = {0, 0}, k3v[2] = {0, 0}, vC = 0, cLo = kInf, cHi = 0, aLo = kInf, aHi = 0;
        double eA[2] = {kInf, 0}, eC[2] = {kInf, 0}, lA[2] = {kInf, 0}, lC[2] = {kInf, 0}, rv[2] = {kInf, 0}, sA[2] = {kInf, 0}, sC[2] = {kInf, 0};
        int vCk = 0;
        auto span = [](double* r, double v) {
            r[0] = std::min(r[0], v);
            r[1] = std::max(r[1], v);
        };
        auto spanStr = [](const double* r, const char* f) {
            return r[0] > r[1] ? std::string("-") : fmtd(f, r[0]) == fmtd(f, r[1]) ? fmtd(f, r[0]) : fmtd(f, r[0]) + " to " + fmtd(f, r[1]);
        };
        for (const Check4Case& c : co.cases) {
            if (c.ring != ring || (ring == 1 && c.en != "physical")) continue;
            rows.push_back({c.crack, c.g, c.lv, c.en, pct(c.rms) + " (" + pct(c.nc) + ")", fmtd("%.1f", c.rms/c.noise), pct(c.intRms),
                            pct(c.lineRms), pct(c.max, "%.1f%%"),
                            c.Er > 1e-12 ? fmtd("%+.1f%%", 100*(c.Ed/c.Er - 1)) + " (" + fmtd("%.1f%%", 100*c.seDE/c.Er) + ")"
                                         : "E dec " + fmtd("%.4f", c.Ed) + " (E ref 0)",
                            c.rmsRefV >= 0 ? pct(c.rmsRefV) : "-"});
            if (c.en != "physical") continue;
            if (c.through) {
                k3[c.lv == "C"] = std::max(k3[c.lv == "C"], c.k3Max);
                continue;
            }
            if (c.lv == "C" && c.g == "g = 1" && c.vErr > vC) {
                vC = c.vErr;
                vCk = c.vWorst;
            }
            if (c.lv == "C") {
                span(eC, 100*(c.Ed/c.Er - 1));
                span(lC, c.lineRms/c.rms);
                span(sC, 100*(c.EdLine - c.ErLine)/c.Er);
                cLo = std::min(cLo, c.rms);
                cHi = std::max(cHi, c.rms);
                if (c.rmsRefV >= 0) span(rv, c.rmsRefV/c.rms);
            } else {
                span(eA, 100*(c.Ed/c.Er - 1));
                span(lA, c.lineRms/c.rms);
                span(sA, 100*(c.EdLine - c.ErLine)/c.Er);
                aLo = std::min(aLo, c.rms);
                aHi = std::max(aHi, c.rms);
            }
            const Check4Case* a = find(ring, c.crack, c.g, "A", "physical");
            if (c.lv == "C" && a) {
                ratioLo = std::min(ratioLo, c.rms/a->rms);
                ratioHi = std::max(ratioHi, c.rms/a->rms);
            }
        }
        for (const Check4Case& c : co.cases)   // K3 with the reference's vertex values
            if (c.ring == ring && c.through && c.en == "physical") k3v[c.lv == "C"] = std::max(k3v[c.lv == "C"], c.rmsRefV);
        f += fixedTable(rows) + "\n";
        f += "  - K1 and K2 (the front inside the ring), physical glue: H RMS " + pct(aLo) + " to " + pct(aHi) + " at Level A, " + pct(cLo) +
             " to " + pct(cHi) + " at Level C (Level C / Level A " + fmtd("%.2f", ratioLo) + " to " + fmtd("%.2f", ratioHi) + "); energy " +
             spanStr(eA, "%+.0f%%") + " at Level A, " + spanStr(eC, "%+.1f%%") + " at Level C (relative to the reference's), of which " +
             spanStr(sA, "%+.0f%%") + " (Level A) and " + spanStr(sC, "%+.1f%%") + " (Level C) lie within 0.1 of the crack-face lines "
             "(same scale). Near the crack-face lines the H RMS is " +
             spanStr(lA, "%.2f") + " times the ring's at Level A and " + spanStr(lC, "%.2f") + " times at Level C.\n";
        f += "  - With the reference's own vertex values in place of the glued ones (no new walks; K1, K2, Level C) the H RMS is " +
             spanStr(rv, "%.2f") + " times the glued one. Level C, g = 1: the glued vertex values J(v) + c_v differ from the "
             "reference's probes by up to " + fmtd("%.3f", vC) + " (at " + kRingVertName[vCk] + ").\n";
        f += "  - K3 (cut through O; exact solution +-1): max abs(H_dec - (+-1)) = " + fmtd("%.3f", k3[0]) + " at Level A, " +
             fmtd("%.3f", k3[1]) + " at Level C (the reference: 0 exactly); with the reference's vertex values the H RMS is " +
             pct(k3v[0]) + " at Level A and " + pct(k3v[1]) + " at Level C.\n";
        if (ring == 1)
            f += "  - The physical and xi glues coincide, as they must with J_T orthogonal: max abs(S_phys - S_xi)/max abs(S) = " +
                 fmtd("%.1e", co.sxRel) + ", max abs(c_phys - c_xi) = " + fmtd("%.1e", co.cxMax) + ".\n";
        else {
            double dLo = kInf, dHi = 0, dcHi = 0;
            for (const Check4Case& c : co.cases)
                if (c.ring == 2 && c.en == "xi") {
                    const Check4Case* p = find(2, c.crack, c.g, c.lv, "physical");
                    if (!p) continue;
                    dLo = std::min(dLo, c.rms/p->rms);
                    dHi = std::max(dHi, c.rms/p->rms);
                    dcHi = std::max(dcHi, c.dcXi);
                }
            f += "  - The xi glue against the physical glue: H RMS ratio " + fmtd("%.2f", dLo) + " to " + fmtd("%.2f", dHi) +
                 ", the glued vertex values differ by up to " + fmtd("%.4f", dcHi) + ".\n";
        }
    }
    {
        std::vector<std::vector<std::string>> rows = {{"crack", "level", "ring 1: H RMS", "ring 2: H RMS", "ring 2 / ring 1"}};
        for (const RingCrack& rc : cracks)
            for (const char* lv : {"A", "C"}) {
                const Check4Case *a = find(1, rc.name, "g = 1", lv, "physical"), *b = find(2, rc.name, "g = 1", lv, "physical");
                if (!a || !b) continue;
                rows.push_back({rc.name, lv, pct(a->rms), pct(b->rms), fmtd("%.2f", b->rms/a->rms)});
            }
        f += "- Ring 2 against ring 1 (the same O and cracks, g = 1, physical glue):\n\n" + fixedTable(rows) + "\n";
    }
    f += "- Sanity of the production recipe (all four faces Dirichlet): where the crack misses a tet, P_k = lambda_k (largest abs(pooled "
         "z) of a corner " + fmtd("%.2f", co.sanPz) + ", abs(z) < 3 at >= " + pct(co.sanPfrac) + " of the points); in ring 1 the Level-C H0 = J in "
         "every tet (largest abs(pooled z) " + fmtd("%.2f", co.sanCz) + ", abs(z) < 3 at >= " + pct(co.sanCfrac) + "). In ring 2 the "
         "Level-C H0 is the xi-harmonic extension of J's data (each tet solves Delta_xi H = 0) while J is harmonic in X; they differ by "
         "up to " + pct(co.diagC) + " RMS (max " + pct(co.diagCmax, "%.1f%%") + "), / 2 g_max.\n";
    md.insert(mdIntro, findingsBlock(f));
    return ok;
}

} // namespace

bool runCouplingGate(int check, bool quick, const std::string& reportPath, bool tipTrace)
{
    auto t0 = std::chrono::steady_clock::now();
    const GateParams gp = gateParams(quick, tipTrace);
    Report out(quick);
    std::vector<std::pair<std::string, std::string>> secs;
    std::time_t now = std::time(nullptr);
    std::tm tmNow{};
#ifdef _WIN32
    localtime_s(&tmNow, &now);
#else
    localtime_r(&now, &tmNow);
#endif
    char date[64];
    std::strftime(date, sizeof date, "%Y-%m-%d %H:%M:%S", &tmNow);
    out("Coupling gate (--gate coupling%s)%s\n", check >= 0 ? (" --check " + std::to_string(check)).c_str() : "",
        quick ? ": QUICK development run (fewer points and walks; never use these numbers)" : "");
    out("date: %s; hardware threads: %u; base seed %llu\n", date, std::thread::hardware_concurrency(),
        (unsigned long long)gp.seed);
    out("walks: per-tet runs N = %d (10x default) in %d batches; Zombie cross-check N = %d (4x default); eps = %.5f\n\n",
        gp.nWalks, gp.nBatches, gp.nZombie, gp.eps);

    // checks 0-3 and the sanity tests share tets A and B: their geometry setup, quadrature point sets and sanity tests are
    // run (and their sections rewritten) unless the run is check 4 alone, whose ring has its own geometry and point sets
    bool ok = true;
    const std::string cacheDir = (std::filesystem::path(reportPath).parent_path()/"gate_cache").string();
    std::string md, geoMd;
    if (check != 4) {
        // geometry setup first, before any solve
        const size_t geoStart = out.text.size();
        ok = printGeometry(out, gp);
        secs.emplace_back("geometry", "## Appendix B: geometry setup\n\nPrinted by the run before the first solve: the vertices, the "
                                      "frames of A and B, U, F, the crack polygons per tet and for the union, and the per-face flags of "
                                      "every run. Appendices C and D add each configuration's own printout.\n\n" +
                                          codeBlock(out.text.substr(geoStart)));

        const TetFrame A = frameA();
        const XiCrack none{}, x1 = xiCrack(A, crackCheck3(0.30));
        // quadrature point sets
        out("\n=== Quadrature point sets ===\n");
        const PointSet ps0 = volumePoints(A, none, gp.nPoints, gp.eps, gp.quick ? 1 : 4);
        const PointSet ps1 = volumePoints(A, x1, gp.nPoints, gp.eps, gp.quick ? 1 : 4);
        const FacePointSet fps = facePoints(gp.quick ? 4 : gp.nFacePerTri);
        std::string qmd = "## Appendix A: quadrature point sets (sanity tests, checks 0 and 1)\n\nStratified rule on Bey-refined leaves "
                          "(exact leaf volumes), Halton points uniform in each leaf; weights V_leaf/m_leaf; points within eps of the "
                          "crack dropped. Checks 2 and 3 list their own point sets.\n\n";
        MdTable qt;
        qt.head = {"point set", "kept points", "dropped (eps)", "max rel err, all points vs exact over T",
                   "max rel err, kept vs T minus tube", "tube volume"};
        for (const auto& [name, ps, xc] : {std::tuple<const char*, const PointSet*, const XiCrack*>{"A, no crack (sanity, check 0)", &ps0, &none},
                                           {"A, crack x = 0.3, y_f = 0.30 (check 1)", &ps1, &x1}}) {
            out("%s: %d leaves (depth 2/3/4: %d/%d/%d), %zu points generated (%d/%d/%d at depth 2/3/4), %zu kept\n", name,
                ps->nLeaves, ps->leavesByDepth[2], ps->leavesByDepth[3], ps->leavesByDepth[4], ps->allXi.size(), ps->pointsByDepth[2],
                ps->pointsByDepth[3], ps->pointsByDepth[4], ps->size());
            quadratureCheck(out, A, *xc, *ps, gp.nTube, &qt, name);
        }
        {
            double s = 0;
            for (double v : fps.w) s += v;
            out("F point set (for the trace fits of checks 2 and 3): %zu points, sum of weights %.15f (area of F = 0.5)\n",
                fps.X.size(), s);
        }
        qmd += qt.str() + "\nF point set: " + std::to_string(fps.X.size()) + " points on 64 equal sub-triangles (area 0.5).\n";
        secs.emplace_back("quadrature", qmd);

        ok = runSanity(out, gp, ps0, md) && ok;
        secs.emplace_back("sanity", md);
        if (check < 0 || check == 0) {
            md.clear();
            ok = runCheck0(out, gp, ps0, md) && ok;
            secs.emplace_back("check0", md);
        }
        if (check < 0 || check == 1) {
            md.clear();
            ok = runCheck1(out, gp, ps1, x1, md) && ok;
            secs.emplace_back("check1", md);
        }
        if (check < 0 || check == 2) {
            md.clear();
            ok = runCheck2(out, gp, fps, cacheDir, md, geoMd) && ok;
            secs.emplace_back("check2", md);
            secs.emplace_back("geometry2", geoMd);
        }
        if (check < 0 || check == 3) {
            md.clear();
            ok = runCheck3(out, gp, fps, cacheDir, md, geoMd) && ok;
            secs.emplace_back("check3", md);
            secs.emplace_back("geometry3", geoMd);
        }
    }
    if (check < 0 || check == 4) {
        md.clear();
        ok = runCheck4(out, gp, cacheDir, md, geoMd) && ok;
        secs.emplace_back("check4", md);
        secs.emplace_back("geometry4", geoMd);
    }
    const double sec = secondsSince(t0);
    out("\ncoupling gate wall time: %.1f s; %s\n", sec, ok ? "all checks run PASS" : "SOME CHECK FAILED");

    // the console summary: the findings of the parts that ran (the summary at the top of gate_report.md repeats them)
    out("\n=== Summary (errors / 2 g_max = 2 in both g modes; nc = noise-corrected) ===\n");
    for (const auto& s : secs)
        if (!findingsOf(s.second).empty()) out.put(findingsPlain(findingsOf(s.second)));
    if (quick) return ok;
    const std::string cmd = "partial_cut --gate coupling" + (check >= 0 ? " --check " + std::to_string(check) : std::string());
    for (auto& s : secs) s.second += "\n_Written by the run of " + std::string(date) + " (`" + cmd + "`)._\n";
    secs.emplace_back("header",
        std::string("# Coupling gate report (`--gate coupling`)\n\n") +
        "Can one tet's WoSt solve be glued to a neighbour through the shared face? Tet A (the reference tet) and its mirror "
        "image B across z = 0 share the face F. Each tet is solved in its own xi frame, with F Dirichlet: a particular solution "
        "H0 (crack +-g, F = a fixed trace) and one response P_k per corner of F (crack 0, F = hat_k). The tets are glued by "
        "minimizing the total Dirichlet energy over the corner values on F, and the result is compared with a reference solve "
        "on the union U = A u B (itself a tetrahedron). Check 4 glues the 8 tets of one nodal ring (every tet with all four faces "
        "Dirichlet) at the ring's 7 vertices and compares with a reference solve on the ring, the octahedron O. The answer per "
        "check is under Summary and in each check's Findings.\n\n"
        "Last run: " + date + ", `" + cmd + "`, " + fmtd("%.0f", sec) + " s, " + (ok ? "every check run passed or was done" :
        "SOME CHECK FAILED") + ". A run rewrites only the sections of the checks it ran; each section ends with the date of "
        "its run. The Polyscope views are written by `--gate coupling --view ... --screenshot ...`.\n\n"
        "Conventions. Errors of H and grad H are divided by **2 g_max = 2** in both g modes (g = 1, and the ramp g = min(1, "
        "sqrt(r/L)) whose local 2g vanishes at the front). RMS is quadrature-weighted over the sample points, max is over the "
        "points. Noise floors are per-point standard errors (RMS and max over the points): the reference's from its batch "
        "spread, the decomposed field's from the delete-one-batch jackknife, which includes the noise of the glued corner "
        "values. \"nc\" = noise-corrected RMS, sqrt(max(0, mean e^2 - mean(var_ref + var_dec)))/2. The own loop runs in double "
        "precision; Zombie, the independent cross-check, runs in float. Checks 0 and 1 have no reference solve: they compare "
        "with exact answers (S = K, P_k = lambda_k) and with a direct solve. `--check 4` alone leaves the sections of tets A "
        "and B (sanity tests, Appendices A and B) as they are.\n\n"
        "Contents: Summary; Sanity tests; Check 0; Check 1; Check 2; Check 3; Check 4; Polyscope views; Appendix A (quadrature "
        "of the sanity tests and checks 0, 1); Appendix B (geometry setup); Appendices C, D and E (geometry per configuration "
        "of checks 2 and 3, and of check 4).\n");
    if (!writeGateReport(reportPath, secs)) {
        std::fprintf(stderr, "could not write %s\n", reportPath.c_str());
        return false;
    }
    std::printf("written: %s\n", reportPath.c_str());
    return ok;
}

bool writeGateReportSection(const std::string& path, const std::string& id, const std::string& body)
{
    return std::filesystem::exists(path) && writeGateReport(path, {{id, body}});
}

GateView gateView(int check, int config, bool tipTrace)
{
    if (check != 2 && check != 3) throw std::invalid_argument("--gate coupling --view needs --check 2 or --check 3");
    const double zcs[2] = {0.60, 0.35}, yfs[3] = {0.15, 0.30, 0.45};
    if (config < 0 || config >= (check == 2 ? 2 : 3))
        throw std::invalid_argument(check == 2 ? "--config: 0 (z_c = 0.60) or 1 (z_c = 0.35)"
                                               : "--config: 0, 1 or 2 (y_f = 0.15, 0.30, 0.45)");
    const GateParams gp = gateParams(false, tipTrace);
    Report out(false);
    const TetFrame A = frameA(), B = frameB();
    const CrackSpec cs = check == 2 ? crackCheck2(zcs[config]) : crackCheck3(yfs[config]);
    const std::string cfg = check == 2 ? "z_c = " + fmtd("%.2f", zcs[config]) : "y_f = " + fmtd("%.2f", yfs[config]);
    const XiCrack xa = xiCrack(A, cs), xb = xiCrack(B, cs);

    // the check's point sets, so that the reference's cache key (points hash) matches
    const PointSet pa = volumePoints(A, xa, gp.nPoints, gp.eps, 4), pb = volumePoints(B, xb, gp.nPoints, gp.eps, 4);
    FacePointSet fps = facePoints(gp.nFacePerTri);
    if (check == 3) {
        const std::array<Vec3, 2> seg = crackSegmentOnF(A, xa);
        int nDropped = 0;
        fps = facePointsOff(fps, seg[0], seg[1], gp.eps, nDropped);
    }
    out("coupling gate view: check %d, %s (crack %s); %zu points in A, %zu in B, %zu on F\n", check, cfg.c_str(), cs.name.c_str(),
        pa.size(), pb.size(), fps.X.size());
    const RefRun ref = check == 2 ? runReference(out, gp, cs, "check2_zc" + fmtd("%.2f", zcs[config]), 8000 + config, pa, pb, fps, "gate_cache")
                                  : runReference(out, gp, cs, "check3_yf" + fmtd("%.2f", yfs[config]), 9100 + config, pa, pb, fps, "gate_cache");
    if (!ref.ok) throw std::runtime_error("the reference on U was not run");

    GateView v;
    v.title = "coupling gate, check " + std::to_string(check) + ", " + cfg + ": crack " + cs.name;
    for (int k = 0; k < 5; k++) v.verts[k] = kGlobal[k];
    v.tets = {std::array<size_t, 4>{0, 1, 2, 3}, std::array<size_t, 4>{0, 2, 1, 4}};
    v.crack = unionCut(cs).polygon;
    v.n = cs.n;
    v.d = cs.d;
    v.X = pa.X;
    v.X.insert(v.X.end(), pb.X.begin(), pb.X.end());
    v.XF = fps.X;

    // one level's fields: the reference (once per g type), the decomposed field and the difference, with its errors
    // printed (the same numbers as gate_report.md's summary table)
    const char* gName[2] = {"g = 1", "ramp"};
    std::vector<double> wAB(pa.w);
    wAB.insert(wAB.end(), pb.w.begin(), pb.w.end());
    auto addLevel = [&](int m, const std::string& level, const FieldAt& dA, const FieldAt& dB, const std::vector<double>& trace,
                        const FieldAt& rA, const FieldAt& rB, const FieldAt& rF, bool withRef) {
        const std::string sfx = " (" + std::string(gName[m]) + ")";
        GateView::Field r{"reference H" + sfx}, dec{"decomposed H" + level + sfx}, dif{"decomposed - reference" + level + sfx};
        r.vol = rA.H;
        r.vol.insert(r.vol.end(), rB.H.begin(), rB.H.end());
        r.face = rF.H;
        dec.vol = dA.H;
        dec.vol.insert(dec.vol.end(), dB.H.begin(), dB.H.end());
        dec.face = trace;
        for (size_t i = 0; i < r.vol.size(); i++) dif.vol.push_back(dec.vol[i] - r.vol[i]);
        for (size_t p = 0; p < r.face.size(); p++) dif.face.push_back(dec.face[p] - r.face[p]);
        dif.difference = true;
        r.shot = dec.shot = dif.shot = m == 0;
        double s = 0, w = 0, mx = 0;
        for (size_t i = 0; i < dif.vol.size(); i++) {
            s += wAB[i]*dif.vol[i]*dif.vol[i];
            w += wAB[i];
            mx = std::max(mx, std::fabs(dif.vol[i]));
        }
        dif.rms = std::sqrt(s/w)/kTwoGMax;
        dif.max = mx/kTwoGMax;
        out("  decomposed - reference%s over A u B, / 2 g_max: H RMS %.2e, max %.2e\n", (level + sfx).c_str(), dif.rms, dif.max);
        if (withRef) v.fields.push_back(std::move(r));
        v.fields.push_back(std::move(dec));
        v.fields.push_back(std::move(dif));
    };
    auto traceOnF = [&](const Eigen::Vector3d& c, const FaceFn* J) {
        std::vector<double> t;
        for (const Vec3& X : fps.X) t.push_back((J ? (*J)(X) : 0.0) + c[0]*hatF(0, X) + c[1]*hatF(1, X) + c[2]*hatF(2, X));
        return t;
    };

    const std::vector<FieldSpec> fields = recipeFields();
    if (check == 2) {
        const OwnRun rec = runOwn(A, xa, onlyF(), fields, pa.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 7500 + config));
        out("A recipe run: %s\n", rec.stats.summary().c_str());
        std::vector<std::vector<double>> coefs;
        for (int k = 0; k < 5; k++) coefs.push_back(unitCoef(5, k));
        const Gram gm = gramOwn(rec, pa, A, coefs);
        const Eigen::Matrix4d Kg = check2Kg(B);
        for (int m = 0; m < 2; m++) {
            const Eigen::Vector4d u = glueRun(gm, m, Kg).u;
            out("  %s: glued c = (%+.5f, %+.5f, %+.5f), d = %+.5f\n", gName[m], u[0], u[1], u[2], u[3]);
            const FieldAt dA = runField(A, {RunTerm{&rec, {m == 0 ? 1.0 : 0.0, m == 1 ? 1.0 : 0.0, u[0], u[1], u[2]}, {}}});
            const FieldAt dB = linearField(B, pb.X, check2LocalB(B, u), {});
            addLevel(m, "", dA, dB, traceOnF(u.head<3>(), nullptr), refField(ref, m, 0, ref.nA), refField(ref, m, ref.nA, ref.nB),
                     refField(ref, m, ref.nA + ref.nB, ref.nF), true);
        }
        return v;
    }
    const FaceFn J[2] = {levelBTrace(cs, GType::One, gp.tipTrace), levelBTrace(cs, GType::SqrtFront, gp.tipTrace)};
    const std::vector<FieldSpec> fieldsB = levelBFields(J);
    const char* lvB = gp.tipTrace ? "C" : "B";   // --face-trace tip: the crack's front solution as J
    const OwnRun recA = runOwn(A, xa, onlyF(), fields, pa.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 9600 + config));
    const OwnRun recB = runOwn(B, xb, onlyF(), fields, pb.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 9610 + config));
    const OwnRun hbA = runOwn(A, xa, onlyF(), fieldsB, pa.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 9620 + config));
    const OwnRun hbB = runOwn(B, xb, onlyF(), fieldsB, pb.xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, 9630 + config));
    out("recipe and Level-B H0 runs (A, B): %.1f s\n", recA.stats.seconds + recB.stats.seconds + hbA.stats.seconds + hbB.stats.seconds);
    const Gram gA = gramCheck3(recA, hbA, pa.w, A), gB = gramCheck3(recB, hbB, pb.w, B);
    for (int m = 0; m < 2; m++) {
        const FieldAt rA = refField(ref, m, 0, ref.nA), rB = refField(ref, m, ref.nA, ref.nB), rF = refField(ref, m, ref.nA + ref.nB, ref.nF);
        for (int L = 0; L < 2; L++) {
            const Eigen::Vector3d c = glueRun3(gA, gB, L ? kH0B + m : m).c;
            out("  Level %s, %s: glued c = (%+.5f, %+.5f, %+.5f)\n", L ? lvB : "A", gName[m], c[0], c[1], c[2]);
            addLevel(m, L ? std::string(", Level ") + lvB : std::string(", Level A"), decField3(A, recA, hbA, L, m, c, nullptr), decField3(B, recB, hbB, L, m, c, nullptr),
                     traceOnF(c, L ? &J[m] : nullptr), rA, rB, rF, L == 0);
        }
    }
    return v;
}

namespace {
RingView gateRingViewFast(const RingViewParams& p, std::atomic<int>* progress);   // RingViewParams::fast (with --ring-bench, below)
}

// The interactive view of check 4: the configuration's per-tet runs and reference as in runRingConfig (the same seeds,
// so the check's own setting reproduces gate_report.md), at the sample density and walk counts asked for. Walks per
// point are rounded to whole batches: 64 batches from 1280 walks up (the check's), else 16, at least 4 walks per batch.
// p.stressCase: that --ring-stress case instead (stressRingView); p.fast: the fast method instead (gateRingViewFast);
// everything below is the default path.
RingView gateRingView(const RingViewParams& p, std::atomic<int>* progress)
{
    if (!p.stressCase.empty()) return stressRingView(p.stressCase, p.pointsPerTet, p.walksTet, p.walksR, p.walksRef, p.geometryOnly);
    if (p.fast) return gateRingViewFast(p, progress);
    if (p.ring != 1 && p.ring != 2) throw std::invalid_argument("ring view: ring 1 or 2");
    if (p.crack < 0 || p.crack > 2) throw std::invalid_argument("ring view: crack 0, 1 or 2 (K1, K2, K3)");
    const auto t0 = std::chrono::steady_clock::now();
    GateParams gp = gateParams(false);
    auto batches = [](int n) { return n >= 1280 ? 64 : 16; };
    auto whole = [](int n, int b) { return std::max(4*b, (n + b/2)/b*b); };
    const bool asCheck = p.pointsPerTet == gp.nRingPoints && p.walksTet == gp.nWalks && p.walksRef == gp.nRef;
    gp.nRingPoints = std::max(64, p.pointsPerTet);
    gp.nBatches = batches(p.walksTet);
    gp.nWalks = whole(p.walksTet, gp.nBatches);
    gp.nRefBatches = batches(p.walksRef);
    gp.nRef = whole(p.walksRef, gp.nRefBatches);
    gp.quick = !asCheck;   // runRingReference reads and writes gate_cache/ only at the check's own setting
    Report out(false);
    const ConvexDomain O = octahedron();
    const Ring R = makeRing(p.ring);
    const RingCrack rc = ringCrack(p.crack);
    const int ri = R.id - 1, ki = p.crack;
    const int nModes = R.id == 1 && !rc.through ? 2 : 1;   // the ramp: ring 1 only, not K3 (as check 4)

    RingView v;
    v.title = "check 4, " + R.name + ", crack " + rc.name + ": " + rc.desc;
    std::copy(R.V.begin(), R.V.end(), v.verts.begin());
    std::vector<RingTetCrack> rt;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        rt.push_back(ringTetCrack(T, rc));
        v.tets.push_back({size_t(T.global[0]), size_t(T.global[1]), size_t(T.global[2]), size_t(T.global[3])});
        v.tetName.push_back(T.name);
        v.tetCut.push_back(char(rt.back().cut()));
    }
    v.crack = rc.cutO.polygon;
    v.n = rc.n;
    v.d = rc.d;
    v.pointsPerTet = gp.nRingPoints;
    v.walksTet = gp.nWalks;
    v.walksRef = gp.nRef;
    if (p.geometryOnly) return v;

    std::vector<PointSet> ps;
    for (int t = 0; t < kRingNT; t++) ps.push_back(ringPoints(R.T[size_t(t)], rt[size_t(t)], rc.cutO, gp.nRingPoints, gp.eps, 1));
    const std::vector<Vec3> probes = ringProbes(R);
    const RefRun ref = runRingReference(out, gp, O, R, rc, "check4_ring" + std::to_string(R.id) + "_" + rc.name,
                                        14000 + 10*uint64_t(ri) + uint64_t(ki), ps, probes, "gate_cache");
    if (!ref.ok) throw std::runtime_error("the reference on the ring was not run (see the console)");

    // the 16 per-tet runs, all four faces Dirichlet, each in the tet's own xi frame
    const FaceFn J[2] = {ringJ(rc, GType::One), ringJ(rc, GType::SqrtFront)};
    const std::vector<FieldSpec> fieldsC = ringLevelCFields(J);
    std::vector<OwnRun> rec, hc;
    double secTets = 0, len = 0;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        const RingTetCrack& r = rt[size_t(t)];
        const Cut* gf = r.cut() ? &r.gFront : nullptr;
        const uint64_t tag = 12000 + 1000*uint64_t(ri) + 100*uint64_t(ki) + 10*uint64_t(t);
        rec.push_back(runOwn(T, r.xc, allDirichlet(), ringRecipeFields(T), ps[size_t(t)].xi, gp.nWalks, gp.nBatches,
                             mixSeed(gp.seed, tag), true, gf));
        hc.push_back(runOwn(T, r.xc, allDirichlet(), fieldsC, ps[size_t(t)].xi, gp.nWalks, gp.nBatches, mixSeed(gp.seed, tag + 1),
                            true, gf));
        secTets += rec.back().stats.seconds + hc.back().stats.seconds;
        len += 0.5*(rec.back().stats.meanWalkLength() + hc.back().stats.meanWalkLength())/kRingNT;
        out("view: %s done (%.2f s)\n", T.name, rec.back().stats.seconds + hc.back().stats.seconds);
        if (progress) ++*progress;
    }

    // the Grams (physical energy: the glue the view shows; check 4 also reports the xi one) and the points
    std::vector<Gram> G;
    for (int t = 0; t < kRingNT; t++) G.push_back(gramRing(rec[size_t(t)], hc[size_t(t)], ps[size_t(t)].w, R.T[size_t(t)], true));
    std::vector<double> w;
    std::vector<size_t> off(kRingNT, 0);
    for (int t = 0; t < kRingNT; t++) {
        off[size_t(t)] = v.X.size();
        v.X.insert(v.X.end(), ps[size_t(t)].X.begin(), ps[size_t(t)].X.end());
        w.insert(w.end(), ps[size_t(t)].w.begin(), ps[size_t(t)].w.end());
        v.tetOf.insert(v.tetOf.end(), ps[size_t(t)].size(), t);
    }
    // each vertex's basis function (P_k of its local corner in the tets around it, 0 elsewhere), and H0 (g = 1: crack
    // +-1, faces J) that the glued H adds the weighted basis functions to
    for (int k = 0; k < kRingNV; k++) {
        v.basis[size_t(k)].assign(v.X.size(), 0.0);
        v.inSupport[size_t(k)].assign(v.X.size(), 0);
    }
    GateView::Field h0{"H0 (g = 1): local solves with crack +-1 and faces = J"};
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        for (int k = 0; k < 4; k++) {
            const FieldAt f = runField(T, {RunTerm{&rec[size_t(t)], unitCoef(6, kRingP0 + k), {}}});
            for (size_t i = 0; i < f.H.size(); i++) {
                v.basis[size_t(T.global[k])][off[size_t(t)] + i] = f.H[i];
                v.inSupport[size_t(T.global[k])][off[size_t(t)] + i] = 1;
            }
        }
        const FieldAt c = runField(T, {RunTerm{&hc[size_t(t)], unitCoef(2, 0), {}}});
        h0.vol.insert(h0.vol.end(), c.H.begin(), c.H.end());
    }

    char buf[512];
    std::snprintf(buf, sizeof buf, "%zu sample points (%d per tet asked); walks per point: %d in each local solve, %d for the GT",
                  v.X.size(), gp.nRingPoints, gp.nWalks, gp.nRef);
    v.report.push_back(buf);
    std::snprintf(buf, sizeof buf, "local solves (2 per tet, 16 in all): %.1f s, mean walk length %.1f; GT (whole-ring solve): %s %.1f s",
                  secTets, len, ref.cached ? "loaded from gate_cache (original run)" : "solved,", ref.stats.seconds);
    v.report.push_back(buf);
    v.report.push_back("error of the glued H against the GT, as a fraction of the jump 2: RMS (average) / max; vertex values: glued / GT");
    const char* gName[2] = {"g = 1", "ramp g"};
    for (int m = 0; m < nModes; m++) {
        const std::string sfx = std::string(" (") + gName[m] + ")";
        GateView::Field rf{"GT" + sfx};
        for (int t = 0; t < kRingNT; t++) {
            const FieldAt f = refField(ref, m, off[size_t(t)], ps[size_t(t)].size());
            rf.vol.insert(rf.vol.end(), f.H.begin(), f.H.end());
        }
        const FieldAt rp = refField(ref, m, ref.nA, ref.nF);   // H at the vertex probes
        v.fields.push_back(rf);
        const Vec7 c = glueRun4(R, G, kRingH0C + m).c;   // the coupling: the 7 vertex weights
        GateView::Field dec{"glued H" + sfx}, dif{"error = glued H - GT" + sfx};
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = R.T[size_t(t)];
            const FieldAt f = decField4(T, rec[size_t(t)], hc[size_t(t)], 1, m, localValues(T, c), nullptr);
            dec.vol.insert(dec.vol.end(), f.H.begin(), f.H.end());
        }
        dif.difference = true;
        double s = 0, sw = 0;
        for (size_t i = 0; i < dec.vol.size(); i++) {
            dif.vol.push_back(dec.vol[i] - rf.vol[i]);
            s += w[i]*dif.vol[i]*dif.vol[i];
            sw += w[i];
            dif.max = std::max(dif.max, std::fabs(dif.vol[i]));
        }
        dif.rms = std::sqrt(s/sw)/kTwoGMax;
        dif.max /= kTwoGMax;
        std::string vt;
        for (int k = 0; k < kRingNV; k++) {   // H at vertex k is its face value there: J + c_k
            std::snprintf(buf, sizeof buf, "%s v%d %+.3f/%+.3f", k ? "," : "", k, c[k] + J[m](R.V[size_t(k)]), rp.H[size_t(k)]);
            vt += buf;
        }
        std::snprintf(buf, sizeof buf, "glued H%s: error RMS %.2f%%, max %.1f%%; vertices%s", sfx.c_str(), 100*dif.rms, 100*dif.max,
                      vt.c_str());
        v.report.push_back(buf);
        v.fields.push_back(std::move(dec));
        v.fields.push_back(std::move(dif));
    }
    v.fields.push_back(std::move(h0));
    std::snprintf(buf, sizeof buf, "total %.1f s", secondsSince(t0));
    v.report.push_back(buf);
    return v;
}

// =====================================================================================
// --ring-perf: where check 4's time goes (profiling for the performance work), the cost of the direct ring solve, and a
// first measured comparison of the current per-tet path, a control-variate prototype and the direct solve against the
// cached reference. Prints to stdout only; reads gate_cache/ and never writes it; writes no report section. Every
// existing mode is untouched (nothing here is called from them).
// =====================================================================================
namespace {

double medianOf(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n%2 ? v[n/2] : 0.5*(v[n/2 - 1] + v[n/2]);
}

// runOwn with explicit walk parameters (the shell sweep); the same scoring and harmonic extensions as runOwn
OwnRun runOwnWp(const TetFrame& T, const XiCrack& xc, const std::vector<FieldSpec>& fields, const std::vector<Vec3>& ptsXi,
                int nWalks, int nBatches, uint64_t seed, const WalkParams& wp, const Cut* gFront)
{
    OwnRun r;
    r.fields = fields;
    QueryPoints q;
    const std::vector<char> grad(ptsXi.size(), 1);
    q.add(ptsXi, &grad, seed);
    const Cut front = gFront ? *gFront : xc.strips.frontCut(xc.strips.nStrips());
    const int K = int(fields.size());
    const TerminalValues tv = [&](const WalkTerminal& t, double* u) {
        if (t.face >= 0) {
            const Vec3 X = T.toX(t.x);
            for (int k = 0; k < K; k++) u[k] = fields[k].face ? fields[k].face(X) : 0.0;
        } else {
            for (int k = 0; k < K; k++) u[k] = fields[k].crack ? double(t.side)*front.g(t.x, fields[k].g) : 0.0;
        }
    };
    const HarmonicExtension ext = [&](const Vec3& xi, double* phi) {
        const Vec3 X = T.toX(xi);
        for (int k = 0; k < K; k++) phi[k] = fields[k].face ? fields[k].face(X) : 0.0;
    };
    r.sums = solveMulti(referenceTet(), xc.strips, q, nWalks, nBatches, wp, K, tv, &r.stats, &ext);
    return r;
}

// One ring and crack at a point density, as check 4 builds them
struct PerfSetup {
    GateParams gp;
    ConvexDomain O;
    Ring R;
    RingCrack rc;
    int ri = 0, ki = 0;
    std::vector<RingTetCrack> rt;
    std::vector<PointSet> ps;
    FaceFn J[2];
    double secGeometry = 0;   // ringTetCrack + ringPoints for the 8 tets
};
void perfSetup(PerfSetup& s, int ring, int crack, int nPts)
{
    s.gp = gateParams(false);
    s.gp.nRingPoints = nPts;
    s.O = octahedron();
    s.R = makeRing(ring);
    s.rc = ringCrack(crack);
    s.ri = ring - 1;
    s.ki = crack;
    const auto t0 = std::chrono::steady_clock::now();
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = s.R.T[size_t(t)];
        s.rt.push_back(ringTetCrack(T, s.rc));
        s.ps.push_back(ringPoints(T, s.rt.back(), s.rc.cutO, nPts, s.gp.eps, 1));
    }
    s.secGeometry = secondsSince(t0);
    s.J[0] = ringJ(s.rc, GType::One);
    s.J[1] = ringJ(s.rc, GType::SqrtFront);
}
uint64_t ownTag(const PerfSetup& s, int t) { return 12000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki) + 10*uint64_t(t); }
size_t totalPoints(const PerfSetup& s)
{
    size_t n = 0;
    for (const PointSet& p : s.ps) n += p.size();
    return n;
}

// The cached ring reference of check 4 at check 4's own point sets, READ ONLY: the key is built exactly as
// runRingReference builds it; nothing is solved or written. False if the file is missing or its key differs.
bool loadRingGT(const PerfSetup& s, RefRun& r)
{
    const GateParams gp = gateParams(false);
    const Ring& R = s.R;
    const RingCrack& rc = s.rc;
    const ConvexDomain& O = s.O;
    for (const PointSet& p : s.ps) r.X.insert(r.X.end(), p.X.begin(), p.X.end());
    r.nA = r.X.size();
    const std::vector<Vec3> probes = ringProbes(R);
    r.X.insert(r.X.end(), probes.begin(), probes.end());
    r.nF = probes.size();
    const StripCut sc = ringStrips(rc);
    const WalkParams wp = ringReferenceWalkParams(O);
    const uint64_t seedBase = mixSeed(gp.seed, 14000 + 10*uint64_t(s.ri) + uint64_t(s.ki));
    char hash[32];
    std::snprintf(hash, sizeof hash, "%016llx", (unsigned long long)fnv1a(r.X.data(), r.X.size()*sizeof(Vec3)));
    std::string key = "coupling gate check 4 reference on O | " + R.name + ", c " + fmtg3(R.V[0]) + " | crack " + rc.name +
                      ": n " + fmtg3(rc.n) + ", d " + fmtg(rc.d) + ", f " + fmtg3(rc.f) + ", t " + fmtg(rc.t) + " | seed " +
                      std::to_string(seedBase) + " | N " + std::to_string(gp.nRef) + " | batches " +
                      std::to_string(gp.nRefBatches) + " | " + std::to_string(r.X.size()) + " points (tets " +
                      std::to_string(r.nA) + ", probes " + std::to_string(r.nF) + "), hash " + hash + " | O";
    for (size_t k = 0; k < O.n.size(); k++) key += " " + fmtg3(O.n[k]) + "." + fmtg(O.c[k]);
    const Cut& f = sc.frame();
    key += ", L_O " + fmtg(O.L) + " | walk: absorbShellRel " + fmtg(wp.absorbShellRel) + " (absolute " +
           fmtg(wp.absorbShellRel*O.L) + "), maxWalkLength " + std::to_string(wp.maxWalkLength) + ", mirrorGradBall " +
           (wp.mirrorGradBall ? "1" : "0") + ", kGradFirstBallShrink " + fmtg(kGradFirstBallShrink) + " | crack frame: n " +
           fmtg3(f.n) + ", d " + fmtg(f.d) + ", u " + fmtg3(f.u) + ", sMin " + fmtg(f.sMin) + ", sMax " + fmtg(f.sMax) +
           ", L of g " + fmtg(f.L) + ", strips " + std::to_string(sc.nStrips()) + ", front s " + fmtg(sc.front()) + ", g types " +
           gTypeName(GType::One) + " / " + gTypeName(GType::SqrtFront) + " | " + kRingRefSolverVersion;
    r.key = key;
    r.file = (std::filesystem::path("gate_cache")/("ref_check4_ring" + std::to_string(R.id) + "_" + rc.name + ".bin")).string();
    r.cached = readRefCache(r);
    return r.cached;
}

// The direct ring solve: the reference's solver (solveMulti on O, 8 Neumann faces, crack +-g, both g types from the same
// walks) at the tets' points (H and grad H) and the vertex probes (H), N walks per point, its own seed
RefRun runDirect(const PerfSetup& s, int N, int B, uint64_t seedTag, double shellAbs = -1)
{
    RefRun r;
    for (const PointSet& p : s.ps) r.X.insert(r.X.end(), p.X.begin(), p.X.end());
    r.nA = r.X.size();
    const std::vector<Vec3> probes = ringProbes(s.R);
    r.X.insert(r.X.end(), probes.begin(), probes.end());
    r.nF = probes.size();
    const StripCut sc = ringStrips(s.rc);
    WalkParams wp = ringReferenceWalkParams(s.O);
    if (shellAbs > 0) wp.absorbShellRel = shellAbs/s.O.L;
    const Cut front = sc.frontCut(sc.nStrips());
    const TerminalValues tv = [&front](const WalkTerminal& t, double* u) {
        u[0] = double(t.side)*front.g(t.x, GType::One);
        u[1] = double(t.side)*front.g(t.x, GType::SqrtFront);
    };
    std::vector<char> grad(r.X.size(), 1);
    std::fill(grad.begin() + std::ptrdiff_t(r.nA), grad.end(), char(0));
    QueryPoints q;
    q.add(r.X, &grad, mixSeed(s.gp.seed, seedTag));
    r.sums = solveMulti(s.O, sc, q, N, B, wp, 2, tv, &r.stats);
    return r;
}

// Level C's J built in a tet's own xi frame (the control variate of H0): the crack plane and J's front line u.X = sJ mapped
// into xi, (a, b) orthonormal xi coordinates across the front line (a along the advance, b = the signed distance to the
// plane, > 0 on the +g side), theta = atan2(b, a): J1 = theta/pi, J2 = sqrt(r/L) sin(theta/2) (no cap at 1, so J2 is
// harmonic off the half-plane). Both are xi-harmonic in the tet minus the half-plane behind the front line, which inside O
// is the crack; so H0 - J_xi has crack data 0 (g = 1) and face data J - J_xi (0 in ring 1, where J_T is a rotation).
struct XiFront {
    Vec3 n = Vec3::Zero(), u = Vec3::Zero(), p0 = Vec3::Zero();
    double dn = 0;
    double a(const Vec3& xi) const { return u.dot(xi - p0); }
    double b(const Vec3& xi) const { return n.dot(xi) - dn; }
    double r(const Vec3& xi) const { return std::hypot(a(xi), b(xi)); }
    double J1(const Vec3& xi) const { return std::atan2(b(xi), a(xi))/std::numbers::pi; }
    Vec3 grad1(const Vec3& xi) const
    {
        const double A = a(xi), B = b(xi);
        return (A*n - B*u)/(std::numbers::pi*(A*A + B*B));
    }
    double J2(const Vec3& xi) const { return std::sqrt(r(xi)/kSqrt2)*std::sin(0.5*std::atan2(b(xi), a(xi))); }
    Vec3 grad2(const Vec3& xi) const
    {
        const double A = a(xi), B = b(xi), rr = std::hypot(A, B), th = std::atan2(B, A);
        const double sn = std::sin(0.5*th), cs = std::cos(0.5*th), f = 1.0/(2.0*rr*std::sqrt(kSqrt2*rr));
        return f*((A*sn - B*cs)*u + (B*sn + A*cs)*n);
    }
};
XiFront xiFront(const TetFrame& T, const RingCrack& rc)
{
    XiFront f;
    const Vec3 nX = rc.n, uX = rc.cutO.u, wX = nX.cross(uX);
    const Vec3 jn = T.J.transpose()*nX;
    f.n = jn.normalized();
    f.dn = (rc.d - nX.dot(T.X[0]))/jn.norm();
    const Vec3 wXi = (T.Jinv*wX).normalized();
    f.u = wXi.cross(f.n);
    if (f.u.dot(T.Jinv*uX) < 0) f.u = -f.u;
    f.p0 = T.toXi(rc.d*nX + rc.sJ*uX);
    return f;
}

// One ring tet's control-variate runs: H0 = J_xi + R (R: crack 0 for g = 1, faces J - J_xi), P_k = lambda_k + Q_k (Q_k:
// crack -lambda_k, faces 0; Q = 0 where the crack misses the tet). single: R and Q on the same walks, plain first ball (R's
// face data has no known harmonic extension); split: Q alone with the odd mirror ball (extension 0), R alone with the plain
// ball, each at its own N. No walks for a field known exactly (R in ring 1, Q in an uncut tet).
struct CvRun {
    XiFront xf;
    bool cut = false, walkR = false, ramp = false;
    double maxJdiff = 0;   // max |J - J_xi| over sampled face points
    OwnRun runR, runQ;
    bool hasR = false, hasQ = false;
    int iR1 = -1, iR2 = -1, iQ0 = -1;   // in runR (single: runQ holds everything)
    const OwnRun* rRun() const { return hasR ? &runR : (hasQ && iR1 >= 0 ? &runQ : nullptr); }
    double seconds() const { return (hasR ? runR.stats.seconds : 0.0) + (hasQ ? runQ.stats.seconds : 0.0); }
    long long walks() const { return (hasR ? runR.stats.walksRequested : 0) + (hasQ ? runQ.stats.walksRequested : 0); }
    long long steps() const { return (hasR ? runR.stats.totalSteps : 0) + (hasQ ? runQ.stats.totalSteps : 0); }
};
CvRun runCv(const PerfSetup& s, int t, bool ramp, bool split, int nQ, int nR, int B, uint64_t seed)
{
    CvRun c;
    const TetFrame& T = s.R.T[size_t(t)];
    const RingTetCrack& rt = s.rt[size_t(t)];
    const PointSet& ps = s.ps[size_t(t)];
    c.xf = xiFront(T, s.rc);
    c.cut = rt.cut();
    c.ramp = ramp;
    // R needs walks iff J - J_xi is not 0 on the faces (sampled on a barycentric grid of each face)
    const Tet ref = referenceTet();
    for (int o = 0; o < 4; o++) {
        std::array<int, 3> fv{};
        int m = 0;
        for (int i = 0; i < 4; i++)
            if (i != o) fv[m++] = i;
        for (int i = 0; i <= 24; i++)
            for (int j = 0; i + j <= 24; j++) {
                const double l1 = i/24.0, l2 = j/24.0, l0 = 1.0 - l1 - l2;
                const Vec3 xi = l0*ref.v[fv[0]] + l1*ref.v[fv[1]] + l2*ref.v[fv[2]];
                if (std::fabs(c.xf.b(xi)) < 1e-9) continue;   // on the plane: atan2's sign is round-off
                const Vec3 X = T.toX(xi);
                c.maxJdiff = std::max(c.maxJdiff, std::fabs(s.J[0](X) - c.xf.J1(xi)));
                if (ramp) c.maxJdiff = std::max(c.maxJdiff, std::fabs(s.J[1](X) - c.xf.J2(xi)));
            }
    }
    c.walkR = c.maxJdiff > 1e-12;
    const Cut* gf = c.cut ? &rt.gFront : nullptr;
    const Cut front = gf ? *gf : rt.xc.strips.frontCut(rt.xc.strips.nStrips());
    const XiFront xf = c.xf;
    const FaceFn J0 = s.J[0], J1 = s.J[1];
    auto makeTv = [&T, xf, J0, J1, front](int iR1, int iR2, int iQ0) -> TerminalValues {
        return [&T, xf, J0, J1, front, iR1, iR2, iQ0](const WalkTerminal& w, double* u) {
            if (w.face >= 0) {
                const Vec3 X = T.toX(w.x);
                if (iR1 >= 0) u[iR1] = J0(X) - xf.J1(w.x);
                if (iR2 >= 0) u[iR2] = J1(X) - xf.J2(w.x);
                if (iQ0 >= 0)
                    for (int k = 0; k < 4; k++) u[iQ0 + k] = 0.0;
            } else {
                if (iR1 >= 0) u[iR1] = 0.0;   // side g - J_xi = side - side
                if (iR2 >= 0) u[iR2] = double(w.side)*(front.g(w.x, GType::SqrtFront) - std::sqrt(xf.r(w.x)/kSqrt2));
                if (iQ0 >= 0) {
                    const double l[4] = {1.0 - w.x.sum(), w.x[0], w.x[1], w.x[2]};
                    for (int k = 0; k < 4; k++) u[iQ0 + k] = -l[k];
                }
            }
        };
    };
    QueryPoints q, qR;
    const std::vector<char> grad(ps.size(), 1);
    q.add(ps.xi, &grad, seed);
    qR.add(ps.xi, &grad, mixSeed(seed, 77));
    const HarmonicExtension ext0 = [](const Vec3&, double* phi) {
        for (int k = 0; k < 4; k++) phi[k] = 0.0;
    };
    if (!split) {
        int K = 0;
        if (c.walkR) {
            c.iR1 = K++;
            if (ramp) c.iR2 = K++;
        }
        if (c.cut) {
            c.iQ0 = K;
            K += 4;
        }
        if (K > 0) {
            const WalkParams wp = gateWalkParams(allDirichlet(), true, false);   // plain first ball (Dirichlet faces only)
            c.runQ.sums = solveMulti(referenceTet(), rt.xc.strips, q, nQ, B, wp, K, makeTv(c.iR1, c.iR2, c.iQ0), &c.runQ.stats);
            c.hasQ = true;
        }
        return c;
    }
    if (c.cut) {
        const WalkParams wp = gateWalkParams(allDirichlet(), true, true);   // odd mirror ball: Q's face data 0, extension 0
        c.runQ.sums = solveMulti(referenceTet(), rt.xc.strips, q, nQ, B, wp, 4, makeTv(-1, -1, 0), &c.runQ.stats, &ext0);
        c.hasQ = true;
        c.iQ0 = 0;
    }
    if (c.walkR) {
        const int K = ramp ? 2 : 1;
        const WalkParams wp = gateWalkParams(allDirichlet(), true, false);
        c.runR.sums = solveMulti(referenceTet(), rt.xc.strips, qR, nR, B, wp, K, makeTv(0, ramp ? 1 : -1, -1), &c.runR.stats);
        c.hasR = true;
        c.iR1 = 0;
        c.iR2 = ramp ? 1 : -1;
    }
    return c;
}
// the Q run's own index of Q_0 (single: after R)
const OwnRun* cvQRun(const CvRun& c) { return c.hasQ && c.iQ0 >= 0 ? &c.runQ : nullptr; }

// Per point and batch: the batch's n, sum of u_k and gradient estimate (xi) of field k of a run
struct BatchSeries {
    std::vector<double> n, su;
    std::vector<Vec3> g;
    double N = 0, SU = 0;
    Vec3 G = Vec3::Zero();
};
void batchSeries(const MultiSums& s, size_t i, int k, const std::vector<double>& unit, BatchSeries& r)
{
    const int B = s.nBatches;
    r.n.assign(size_t(B), 0.0);
    r.su.assign(size_t(B), 0.0);
    r.g.assign(size_t(B), Vec3::Zero());
    r.N = r.SU = 0;
    r.G = Vec3::Zero();
    for (int b = 0; b < B; b++) {
        const double* d = s.batchData(i, b);
        r.n[size_t(b)] = d[0];
        r.su[size_t(b)] = d[1 + k];
        r.g[size_t(b)] = comboBatchGrad(s, i, b, unit);
        r.N += d[0];
        r.SU += d[1 + k];
        r.G += r.g[size_t(b)];
    }
}

// The Gram of one ring tet in gramRing's layout (H0 g1, H0 g2, P_0..P_3) from the control-variate runs: per batch b the
// gradient of H0 is grad J_xi (exact) + R's batch estimate, of P_k grad lambda_k (exact) + Q_k's; physical energy.
// (w: the weights per point; the overload without them uses the point set's own)
Gram cvGram(const CvRun& c, const PointSet& ps, const std::vector<double>& w, const TetFrame& T, int B)
{
    const std::array<Vec3, 4> gN = gradN(T);
    const OwnRun* rr = c.rRun();
    const OwnRun* qr = cvQRun(c);
    std::vector<std::vector<double>> uR, uQ;
    if (rr)
        for (int k = 0; k < rr->sums.K; k++) uR.push_back(unitCoef(rr->sums.K, k));
    if (qr)
        for (int k = 0; k < qr->sums.K; k++) uQ.push_back(unitCoef(qr->sums.K, k));
    return gramFromBatches(w, B, 6, [&](size_t p, int b, Vec3* g) {
        const Vec3& xi = ps.xi[p];
        Vec3 h1 = c.xf.grad1(xi), h2 = c.ramp ? c.xf.grad2(xi) : h1;
        if (rr && c.iR1 >= 0) h1 += comboBatchGrad(rr->sums, p, b, uR[size_t(c.iR1)]);
        if (rr && c.iR2 >= 0) h2 += comboBatchGrad(rr->sums, p, b, uR[size_t(c.iR2)]);
        g[0] = T.JinvT*h1;
        g[1] = T.JinvT*h2;
        for (int k = 0; k < 4; k++)
            g[2 + k] = gN[size_t(k)] + (qr ? Vec3(T.JinvT*comboBatchGrad(qr->sums, p, b, uQ[size_t(c.iQ0 + k)])) : Vec3::Zero());
    });
}
Gram cvGram(const CvRun& c, const PointSet& ps, const TetFrame& T, int B) { return cvGram(c, ps, ps.w, T, B); }

// The glued field of one ring tet from the control-variate runs: H = J_xi + R + sum_k c_k (lambda_k + Q_k), physical
// gradient, delete-one-batch jackknife variances with the replicates' coefficients cJ
FieldAt cvField(const CvRun& c, const PointSet& ps, const TetFrame& T, int m, const Eigen::Vector4d& cv,
                const std::vector<Eigen::Vector4d>& cJ, int B)
{
    const std::array<Vec3, 4> gN = gradN(T);
    const OwnRun* rr = c.rRun();
    const OwnRun* qr = cvQRun(c);
    const int iR = m == 0 ? c.iR1 : c.iR2;
    std::vector<std::vector<double>> uR, uQ;
    if (rr)
        for (int k = 0; k < rr->sums.K; k++) uR.push_back(unitCoef(rr->sums.K, k));
    if (qr)
        for (int k = 0; k < qr->sums.K; k++) uQ.push_back(unitCoef(qr->sums.K, k));
    FieldAt f;
    f.resize(ps.size());
    tbb::parallel_for(size_t(0), ps.size(), [&](size_t i) {
        const Vec3& xi = ps.xi[i];
        const double lam[4] = {1.0 - xi.sum(), xi[0], xi[1], xi[2]};
        const double h = m == 0 ? c.xf.J1(xi) : c.xf.J2(xi);
        const Vec3 gh = m == 0 ? c.xf.grad1(xi) : c.xf.grad2(xi);
        BatchSeries sR, sQ[4];
        const bool hasR = rr && iR >= 0;
        if (hasR) batchSeries(rr->sums, i, iR, uR[size_t(iR)], sR);
        if (qr)
            for (int k = 0; k < 4; k++) batchSeries(qr->sums, i, c.iQ0 + k, uQ[size_t(c.iQ0 + k)], sQ[k]);
        auto eval = [&](const Eigen::Vector4d& cc, int a, double& H, Vec3& G) {   // a < 0: all batches
            double Rv = 0, Qv[4] = {0, 0, 0, 0};
            Vec3 Rg = Vec3::Zero(), Qg[4] = {Vec3::Zero(), Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};
            if (hasR) {
                Rv = a < 0 ? sR.SU/sR.N : (sR.SU - sR.su[size_t(a)])/(sR.N - sR.n[size_t(a)]);
                Rg = a < 0 ? Vec3(sR.G/double(B)) : Vec3((sR.G - sR.g[size_t(a)])/double(B - 1));
            }
            if (qr)
                for (int k = 0; k < 4; k++) {
                    Qv[k] = a < 0 ? sQ[k].SU/sQ[k].N : (sQ[k].SU - sQ[k].su[size_t(a)])/(sQ[k].N - sQ[k].n[size_t(a)]);
                    Qg[k] = a < 0 ? Vec3(sQ[k].G/double(B)) : Vec3((sQ[k].G - sQ[k].g[size_t(a)])/double(B - 1));
                }
            H = h + Rv;
            G = T.JinvT*(gh + Rg);
            for (int k = 0; k < 4; k++) {
                H += cc[k]*(lam[k] + Qv[k]);
                G += cc[k]*(gN[size_t(k)] + T.JinvT*Qg[k]);
            }
        };
        eval(cv, -1, f.H[i], f.G[i]);
        std::vector<double> Ha(static_cast<size_t>(B));
        std::vector<Vec3> Ga(static_cast<size_t>(B));
        double hm = 0;
        Vec3 gm = Vec3::Zero();
        for (int a = 0; a < B; a++) {
            eval(cJ.empty() ? cv : cJ[size_t(a)], a, Ha[size_t(a)], Ga[size_t(a)]);
            hm += Ha[size_t(a)]/double(B);
            gm += Ga[size_t(a)]/double(B);
        }
        for (int a = 0; a < B; a++) {
            f.vH[i] += (Ha[size_t(a)] - hm)*(Ha[size_t(a)] - hm)*double(B - 1)/double(B);
            f.vG[i] += (Ga[size_t(a)] - gm).cwiseAbs2()*double(B - 1)/double(B);
        }
    });
    return f;
}

// One method's result at the ring's points against the reference
struct PerfResult {
    std::string name;
    double secWalks = 0, secPost = 0;
    long long walks = 0, steps = 0;
    ErrStats e;
    double gtGrad2 = 0;   // sum w |grad GT|^2 over the same points
    Vec7 vert = Vec7::Zero(), cse = Vec7::Zero();
    double vErr = 0;
    bool hasC = false;
};
void printPerfHeader(Report& out)
{
    out("    %-36s %8s %8s %8s %8s %8s | %7s %7s %7s %8s %8s %8s | %7s %7s\n", "method", "walks/M", "steps/M", "walk s",
        "post s", "total s", "H RMS%", "nc %", "max %", "gradRel%", "grad nc%", "vert max", "fld se%", "c se max");
}
void printPerf(Report& out, const PerfResult& r)
{
    const double g2 = std::max(1e-300, r.gtGrad2);
    out("    %-36s %8.2f %8.1f %8.2f %8.2f %8.2f | %7.3f %7.3f %7.2f %8.2f %8.2f %8.4f | %7.3f %7.4f\n", r.name.c_str(),
        1e-6*double(r.walks), 1e-6*double(r.steps), r.secWalks, r.secPost, r.secWalks + r.secPost, 100*r.e.rmsH(),
        100*r.e.corrH(), 100*r.e.maxH/kTwoGMax, 100*std::sqrt(r.e.eG2/g2),
        100*std::sqrt(std::max(0.0, r.e.eG2 - r.e.vRefG - r.e.vFldG)/g2), r.vErr, 100*r.e.seFldH(), r.hasC ? r.cse.maxCoeff() : 0.0);
}

// ---- part a: the current per-tet runs (two per tet), timing per tet at a reduced setting and at check 4's ----
void perfPartA(Report& out, int ring, int crack)
{
    out("\n=== (a) check 4's per-tet path as it is: two runs per tet (recipe run: H0 Level A x2 + P_0..P_3, odd mirror "
        "ball; Level-C run: H0 x2, plain first ball) ===\n");
    struct Setting {
        int nPts, N, B, reps;
        bool allTets;
    };
    for (const Setting st : {Setting{1024, 640, 16, 3, true}, Setting{4096, 2560, 64, 1, false}}) {
        PerfSetup s;
        perfSetup(s, ring, crack, st.nPts);
        out("\nring %d, %s, about %d points per tet (%zu kept in all; geometry + point sets %.3f s), N = %d in %d batches, "
            "%d repeat(s), 16 threads (TBB default)\n", ring, s.rc.name.c_str(), st.nPts, totalPoints(s), s.secGeometry, st.N,
            st.B, st.reps);
        out("    tet   cut  points | recipe s  len   us/walk | Level-C s  len   us/walk\n");
        double sumRec = 0, sumC = 0;
        long long walks = 0;
        std::vector<OwnRun> rec(kRingNT), hc(kRingNT);
        int nUncut = 0, nCutDone = 0;
        const std::vector<FieldSpec> fieldsC = ringLevelCFields(s.J);
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = s.R.T[size_t(t)];
            const RingTetCrack& r = s.rt[size_t(t)];
            if (!st.allTets) {   // one cut and one uncut tet
                if (r.cut() ? nCutDone++ > 0 : nUncut++ > 0) continue;
            }
            const Cut* gf = r.cut() ? &r.gFront : nullptr;
            std::vector<double> tr, tc;
            for (int rep = 0; rep < st.reps; rep++) {
                rec[size_t(t)] = runOwn(T, r.xc, allDirichlet(), ringRecipeFields(T), s.ps[size_t(t)].xi, st.N, st.B,
                                        mixSeed(s.gp.seed, ownTag(s, t)), true, gf);
                hc[size_t(t)] = runOwn(T, r.xc, allDirichlet(), fieldsC, s.ps[size_t(t)].xi, st.N, st.B,
                                       mixSeed(s.gp.seed, ownTag(s, t) + 1), true, gf);
                tr.push_back(rec[size_t(t)].stats.seconds);
                tc.push_back(hc[size_t(t)].stats.seconds);
            }
            const SolveStats &a = rec[size_t(t)].stats, &b = hc[size_t(t)].stats;
            const double mr = medianOf(tr), mc = medianOf(tc);
            sumRec += mr;
            sumC += mc;
            walks += a.walksRequested + b.walksRequested;
            out("    %-5s %-4s %6zu | %7.3f %5.1f %8.3f | %8.3f %5.1f %8.3f   (repeats rec:", T.name, r.cut() ? "yes" : "no",
                s.ps[size_t(t)].size(), mr, a.meanWalkLength(), 1e6*mr/double(a.walksRequested), mc, b.meanWalkLength(),
                1e6*mc/double(b.walksRequested));
            for (double x : tr) out(" %.3f", x);
            out("; C:");
            for (double x : tc) out(" %.3f", x);
            out(")\n");
        }
        out("    sum of the medians: recipe runs %.2f s + Level-C runs %.2f s = %.2f s (%.1f M walks)\n", sumRec, sumC, sumRec + sumC,
            1e-6*double(walks));
        if (!st.allTets) continue;
        // the rest of the path at this setting: Grams, glue, glued field (Level C, g = 1)
        auto t0 = std::chrono::steady_clock::now();
        std::vector<Gram> G;
        for (int t = 0; t < kRingNT; t++) G.push_back(gramRing(rec[size_t(t)], hc[size_t(t)], s.ps[size_t(t)].w, s.R.T[size_t(t)], true));
        const double secGram = secondsSince(t0);
        t0 = std::chrono::steady_clock::now();
        const Glue4 gl = glueRun4(s.R, G, kRingH0C);
        const double secGlue = secondsSince(t0);
        t0 = std::chrono::steady_clock::now();
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = s.R.T[size_t(t)];
            std::vector<Eigen::Vector4d> cJ;
            for (const Vec7& x : gl.jack) cJ.push_back(localValues(T, x));
            (void)decField4(T, rec[size_t(t)], hc[size_t(t)], 1, 0, localValues(T, gl.c), &cJ);
        }
        const double secField = secondsSince(t0);
        out("    after the walks: 8 Grams (gramRing, 8 x 8, jackknife) %.3f s; glue (1 + %d replicate solves) %.4f s; glued field "
            "with jackknife at all points %.3f s\n", secGram, int(gl.jack.size()), secGlue, secField);
        out("    glued c (Level C, g = 1) = %s\n", fmtVec(gl.c).c_str());
    }
}

// ---- part b: the walk stop distance (shell): walk length, time, and its bias by paired runs (same seeds) ----
void perfPartB(Report& out, int ring, int crack)
{
    out("\n=== (b) the walk stop distance: per-tet runs (all 8 tets) and the direct solve, paired runs (same seeds) against "
        "shell = 1e-6 sqrt(2) ===\n");
    PerfSetup s;
    perfSetup(s, ring, crack, 1024);
    const int N = 640, B = 16, ND = 64;
    const std::vector<FieldSpec> fieldsC = ringLevelCFields(s.J);
    struct ShellRun {
        double rel = 0;
        std::vector<OwnRun> rec, hc;
        double secRec = 0, secC = 0, lenRec = 0, lenC = 0;
        Vec7 c = Vec7::Zero();
        std::vector<FieldAt> H;
        RefRun direct;
    };
    std::vector<ShellRun> runs;
    for (double rel : {1e-6, 1e-5, 1e-4, 1e-3}) {
        ShellRun sr;
        sr.rel = rel;
        long long wr = 0, wc = 0;
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = s.R.T[size_t(t)];
            const RingTetCrack& r = s.rt[size_t(t)];
            const Cut* gf = r.cut() ? &r.gFront : nullptr;
            WalkParams wr1 = gateWalkParams(allDirichlet(), true, true), wc1 = gateWalkParams(allDirichlet(), true, false);
            wr1.absorbShellRel = wc1.absorbShellRel = rel;
            sr.rec.push_back(runOwnWp(T, r.xc, ringRecipeFields(T), s.ps[size_t(t)].xi, N, B, mixSeed(s.gp.seed, ownTag(s, t)),
                                      wr1, gf));
            sr.hc.push_back(runOwnWp(T, r.xc, fieldsC, s.ps[size_t(t)].xi, N, B, mixSeed(s.gp.seed, ownTag(s, t) + 1), wc1, gf));
            sr.secRec += sr.rec.back().stats.seconds;
            sr.secC += sr.hc.back().stats.seconds;
            sr.lenRec += double(sr.rec.back().stats.totalSteps);
            sr.lenC += double(sr.hc.back().stats.totalSteps);
            wr += sr.rec.back().stats.walksCounted;
            wc += sr.hc.back().stats.walksCounted;
        }
        sr.lenRec /= double(wr);
        sr.lenC /= double(wc);
        std::vector<Gram> G;
        for (int t = 0; t < kRingNT; t++) G.push_back(gramRing(sr.rec[size_t(t)], sr.hc[size_t(t)], s.ps[size_t(t)].w, s.R.T[size_t(t)], true));
        sr.c = glueRun4(s.R, G, kRingH0C).c;
        for (int t = 0; t < kRingNT; t++)
            sr.H.push_back(decField4(s.R.T[size_t(t)], sr.rec[size_t(t)], sr.hc[size_t(t)], 1, 0, localValues(s.R.T[size_t(t)], sr.c), nullptr));
        sr.direct = runDirect(s, ND, B, 17900 + 10*uint64_t(s.ri) + uint64_t(s.ki), rel*kSqrt2);
        runs.push_back(std::move(sr));
    }
    out("ring %d, %s, about 1024 points per tet (%zu in all), per-tet N = %d, direct N = %d, %d batches\n", ring,
        s.rc.name.c_str(), totalPoints(s), N, ND, B);
    out("    shell (abs)  | 8 tets: recipe s  len | Level-C s  len | paired vs 1e-6: max|dc|   glued H RMS%%  max%%  "
        "P_k RMS (max)        | direct: s     len   paired H RMS%% max%%\n");
    const ShellRun& r0 = runs.front();
    for (const ShellRun& sr : runs) {
        // paired differences: glued H, and P_k (recipe run, unit coefficients) against the 1e-6 run
        double e2 = 0, w = 0, emax = 0, p2 = 0, pw = 0, pmax = 0;
        for (int t = 0; t < kRingNT; t++) {
            const PointSet& p = s.ps[size_t(t)];
            for (size_t i = 0; i < p.size(); i++) {
                const double d = sr.H[size_t(t)].H[i] - r0.H[size_t(t)].H[i];
                e2 += p.w[i]*d*d;
                w += p.w[i];
                emax = std::max(emax, std::fabs(d));
                for (int k = 0; k < 4; k++) {
                    double h1, h0, v;
                    long long n;
                    comboH(sr.rec[size_t(t)].sums, i, unitCoef(6, kRingP0 + k), h1, v, n);
                    comboH(r0.rec[size_t(t)].sums, i, unitCoef(6, kRingP0 + k), h0, v, n);
                    p2 += p.w[i]*(h1 - h0)*(h1 - h0);
                    pw += p.w[i];
                    pmax = std::max(pmax, std::fabs(h1 - h0));
                }
            }
        }
        double d2 = 0, dmax = 0, dw = 0;
        {
            size_t off = 0;
            for (int t = 0; t < kRingNT; t++) {
                const PointSet& p = s.ps[size_t(t)];
                for (size_t i = 0; i < p.size(); i++) {
                    double h1, h0, v;
                    long long n;
                    comboH(sr.direct.sums, off + i, unitCoef(2, 0), h1, v, n);
                    comboH(r0.direct.sums, off + i, unitCoef(2, 0), h0, v, n);
                    d2 += p.w[i]*(h1 - h0)*(h1 - h0);
                    dw += p.w[i];
                    dmax = std::max(dmax, std::fabs(h1 - h0));
                }
                off += p.size();
            }
        }
        out("    %.1e      | %8.3f %5.1f | %8.3f %5.1f | %.2e     %.4f  %.3f  %.1e (%.1e) | %7.3f %6.1f   %.4f %.3f\n",
            sr.rel*kSqrt2, sr.secRec, sr.lenRec, sr.secC, sr.lenC, (sr.c - r0.c).cwiseAbs().maxCoeff(),
            100*std::sqrt(e2/w)/kTwoGMax, 100*emax/kTwoGMax, std::sqrt(p2/pw), pmax, sr.direct.stats.seconds,
            sr.direct.stats.meanWalkLength(), 100*std::sqrt(d2/dw)/kTwoGMax, 100*dmax/kTwoGMax);
    }
    // The paired differences above mix the shell's bias with the noise of the walks whose values the earlier stop changes.
    // Batch by batch the paired difference D_b of a field's batch means is independent across batches, so
    // mean_p (Dbar^2 - var(Dbar)) estimates the squared bias alone: noise-corrected RMS, with the signed mean and the
    // pooled z of Dbar (bias against its own noise).
    out("  noise-corrected paired bias (batch-paired differences, %d batches), / 2g: field | shell: signed mean, nc RMS (raw RMS), "
        "pooled z\n", B);
    auto pairedBias = [&](const char* name, auto&& sumsOf, int k, bool direct) {
        out("    %-30s", name);
        for (size_t ri = 1; ri < runs.size(); ri++) {
            double sw = 0, sMean = 0, sD2 = 0, sV = 0, zs = 0;
            long long nz = 0;
            size_t off = 0;
            for (int t = 0; t < kRingNT; t++) {
                const PointSet& p = s.ps[size_t(t)];
                const MultiSums& a = sumsOf(runs[ri], t);
                const MultiSums& b = sumsOf(runs[0], t);
                for (size_t i = 0; i < p.size(); i++) {
                    const size_t ii = direct ? off + i : i;
                    std::vector<double> D(static_cast<size_t>(B));
                    double Dm = 0;
                    for (int bb = 0; bb < B; bb++) {
                        const double* da = a.batchData(ii, bb);
                        const double* db = b.batchData(ii, bb);
                        D[size_t(bb)] = da[1 + k]/da[0] - db[1 + k]/db[0];
                        Dm += D[size_t(bb)]/double(B);
                    }
                    double v = 0;
                    for (double x : D) v += (x - Dm)*(x - Dm);
                    v /= double(B - 1)*double(B);
                    sw += p.w[i];
                    sMean += p.w[i]*Dm;
                    sD2 += p.w[i]*Dm*Dm;
                    sV += p.w[i]*v;
                    if (v > 0) {
                        zs += Dm/std::sqrt(v);
                        nz++;
                    }
                }
                off += p.size();
            }
            out(" | %.0e: %+.4f%%, %.4f%% (%.4f%%), z %+.1f", runs[ri].rel*kSqrt2, 100*sMean/sw/kTwoGMax,
                100*std::sqrt(std::max(0.0, sD2 - sV)/sw)/kTwoGMax, 100*std::sqrt(sD2/sw)/kTwoGMax, nz ? zs/std::sqrt(double(nz)) : 0.0);
        }
        out("\n");
    };
    pairedBias("H0 Level C g = 1 (per-tet)", [](const ShellRun& r, int t) -> const MultiSums& { return r.hc[size_t(t)].sums; }, 0, false);
    pairedBias("H0 Level C ramp (per-tet)", [](const ShellRun& r, int t) -> const MultiSums& { return r.hc[size_t(t)].sums; }, 1, false);
    pairedBias("P_0 (per-tet recipe run)", [](const ShellRun& r, int t) -> const MultiSums& { return r.rec[size_t(t)].sums; }, kRingP0, false);
    pairedBias("H g = 1 (direct on O, N = 64)", [](const ShellRun& r, int) -> const MultiSums& { return r.direct.sums; }, 0, true);
}

// ---- part c: the direct ring solve on O, cost at small N ----
void perfPartC(Report& out, int ring, int crack)
{
    out("\n=== (c) the direct ring solve (solveMulti on the octahedron O, 8 Neumann faces, crack +-g, K = 2) at the ring's "
        "points ===\n");
    for (const int nPts : {1024, 4096}) {
        PerfSetup s;
        perfSetup(s, ring, crack, nPts);
        for (const int N : {64, 256}) {
            if (nPts == 4096 && N == 256) continue;
            std::vector<double> sec;
            RefRun d;
            for (int rep = 0; rep < (nPts == 1024 ? 3 : 1); rep++) {
                d = runDirect(s, N, 16, 17000 + 10*uint64_t(s.ri) + uint64_t(s.ki));
                sec.push_back(d.stats.seconds);
            }
            const double m = medianOf(sec);
            out("    ring %d %s, %zu points (+7 probes), N = %3d: %.3f s (median of %zu:", ring, s.rc.name.c_str(), d.nA, N, m,
                sec.size());
            for (double x : sec) out(" %.3f", x);
            out("), %.2f M walks, mean walk length %.1f, %.3f us/walk wall, capped %lld, escaped %lld\n",
                1e-6*double(d.stats.walksRequested), d.stats.meanWalkLength(), 1e6*m/double(d.stats.walksRequested),
                d.stats.walksCapped, d.stats.walksEscaped);
        }
    }
}

// ---- part d: the step function's parts (single thread) and the thread scaling ----
void perfPartD(Report& out, int ring, int crack)
{
    out("\n=== (d) where a walk's time goes: profileWalks (one pinned thread; each part alone over the recorded steps) ===\n");
    PerfSetup s;
    perfSetup(s, ring, crack, 1024);
    int tCut = -1, tUncut = -1;
    for (int t = 0; t < kRingNT; t++) {
        if (s.rt[size_t(t)].cut() && tCut < 0) tCut = t;
        if (!s.rt[size_t(t)].cut() && tUncut < 0) tUncut = t;
    }
    // the per-tet runs' TerminalValues and extensions, as runOwn builds them
    struct Job {
        std::string name;
        int t = -1;   // -1: the direct solve on O
        bool recipe = true;
    };
    std::vector<Job> jobs = {{"cut tet, recipe run (K = 6, odd mirror)", tCut, true}, {"cut tet, Level-C run (K = 2, plain)", tCut, false}};
    if (tUncut >= 0) jobs.push_back({"uncut tet, recipe run (K = 6, odd mirror)", tUncut, true});
    jobs.push_back({"direct on O (K = 2, mirror across O's faces)", -1, true});
    const std::vector<FieldSpec> fieldsC = ringLevelCFields(s.J);
    const StripCut scO = ringStrips(s.rc);
    const Cut frontO = scO.frontCut(scO.nStrips());
    const TerminalValues tvO = [&frontO](const WalkTerminal& t, double* u) {
        u[0] = double(t.side)*frontO.g(t.x, GType::One);
        u[1] = double(t.side)*frontO.g(t.x, GType::SqrtFront);
    };
    std::vector<Vec3> XO;
    for (const PointSet& p : s.ps) XO.insert(XO.end(), p.X.begin(), p.X.end());
    auto profileJob = [&](const Job& j, size_t maxPts, int N, int cpu, bool comp) {
        if (j.t < 0) {
            QueryPoints q;
            const std::vector<char> g(XO.size(), 1);
            q.add(XO, &g, mixSeed(s.gp.seed, 17500));
            return profileWalks(nullptr, &s.O, scO, q, N, ringReferenceWalkParams(s.O), 2, tvO, nullptr, maxPts, cpu, comp);
        }
        const TetFrame& T = s.R.T[size_t(j.t)];
        const RingTetCrack& r = s.rt[size_t(j.t)];
        const std::vector<FieldSpec> fields = j.recipe ? ringRecipeFields(T) : fieldsC;
        bool allHarmonic = true;
        for (const FieldSpec& f : fields) allHarmonic = allHarmonic && (!f.face || f.faceHarmonic);
        const WalkParams wp = gateWalkParams(allDirichlet(), true, allHarmonic);
        QueryPoints q;
        const std::vector<char> g(s.ps[size_t(j.t)].size(), 1);
        q.add(s.ps[size_t(j.t)].xi, &g, mixSeed(s.gp.seed, ownTag(s, j.t)));
        const Cut front = r.cut() ? r.gFront : r.xc.strips.frontCut(r.xc.strips.nStrips());
        const int K = int(fields.size());
        const TerminalValues tv = [&](const WalkTerminal& t, double* u) {
            if (t.face >= 0) {
                const Vec3 X = T.toX(t.x);
                for (int k = 0; k < K; k++) u[k] = fields[k].face ? fields[k].face(X) : 0.0;
            } else {
                for (int k = 0; k < K; k++) u[k] = fields[k].crack ? double(t.side)*front.g(t.x, fields[k].g) : 0.0;
            }
        };
        const HarmonicExtension ext = [&](const Vec3& xi, double* phi) {
            const Vec3 X = T.toX(xi);
            for (int k = 0; k < K; k++) phi[k] = fields[k].face ? fields[k].face(X) : 0.0;
        };
        const Tet ref = referenceTet();
        return profileWalks(&ref, nullptr, r.xc.strips, q, N, wp, K, tv, &ext, maxPts, cpu, comp);
    };
    // per logical processor: the same short job, walks only (the CPU is hybrid: performance and efficiency cores)
    out("per logical processor (cut tet, recipe run, 16 points x 128 walks, one thread pinned): us per walk (full loop)\n   ");
    std::vector<double> perCpu;
    for (int cpu = 0; cpu < 16; cpu++) {
        double best = kInf;
        for (int rep = 0; rep < 3; rep++) {
            const WalkProfile p = profileJob(jobs[0], 16, 128, cpu, false);
            best = std::min(best, 1e6*p.secFull/double(p.walks));
        }
        perCpu.push_back(best);
        out(" cpu%d %.2f", cpu, best);
    }
    const int fast = int(std::min_element(perCpu.begin(), perCpu.end()) - perCpu.begin());
    const int slow = int(std::max_element(perCpu.begin(), perCpu.end()) - perCpu.begin());
    double sumRate = 0;
    for (double x : perCpu) sumRate += 1.0/x;
    out("\n   fastest cpu%d, slowest cpu%d; sum of the 16 single-thread rates = %.1f walks/us, i.e. an ideal 16-thread %.3f us/walk\n",
        fast, slow, sumRate, 1.0/sumRate);
    for (const Job& j : jobs) {
        for (const int cpu : {fast, slow}) {
            const WalkProfile p = profileJob(j, 128, 256, cpu, true);
            const double mv = double(p.moves)/double(p.walks), ev = double(p.evals)/double(p.walks);
            const double usFull = 1e6*p.secFull/double(p.walks), usWalk = 1e6*p.secWalks/double(p.walks);
            // modeled per walk from the parts: ev x (closest + faceDist + sqrt + addBall) + mv x (sphere + outside) + (mv - 1) x star
            // + first ball + mirror + init + finish (+ values + accum per landed walk)
            const double perEval = p.nsClosest + p.nsFaceDist + p.nsSqrt + p.nsAddBall, perMove = p.nsSphere + p.nsOutside;
            const double landedFrac = double(p.landedCut + p.landedFace)/double(p.walks);
            const double model = 1e-3*(ev*perEval + mv*perMove + (mv - 1)*p.nsStar + p.nsFirstBall + p.nsMirror + p.nsInit +
                                       landedFrac*(p.nsFinish + p.nsValues + p.nsAccum));
            out("\n  %s, cpu%d: %lld points x %lld walks; moves/walk %.2f; landed on the crack %.1f%%, on a face %.1f%%, mirrored "
                "first ball %.1f%%, capped %lld, escaped %lld\n", j.name.c_str(), cpu, p.points, p.walks/std::max(1LL, p.points), mv,
                100.0*double(p.landedCut)/double(p.walks), 100.0*double(p.landedFace)/double(p.walks),
                100.0*double(p.mirrored)/double(p.walks), p.capped, p.escaped);
            out("    full loop %.3f us/walk (%.1f ns per move); runWalk alone %.3f us/walk; terminal values + accumulation %.3f us/walk\n",
                usFull, 1e3*usFull/mv, usWalk, usFull - usWalk);
            out("    per step evaluation (x %.2f/walk): closest %.1f ns, dirichletFaceDist %.1f, sqrt %.1f, addBall %.1f\n", ev,
                p.nsClosest, p.nsFaceDist, p.nsSqrt, p.nsAddBall);
            out("    per move (x %.2f/walk): sampleSphere %.1f ns, starStep %.1f, escape test %.1f\n", mv, p.nsSphere, p.nsStar,
                p.nsOutside);
            out("    per walk: makeInitialState %.1f ns, first ball %.1f, mirrorGradBall %.1f, finishWalk %.1f, TerminalValues %.1f, "
                "accumulation %.1f\n", p.nsInit, p.nsFirstBall, p.nsMirror, p.nsFinish, p.nsValues, p.nsAccum);
            out("    sum of the parts %.3f us/walk = %.0f%% of the full loop (the rest: loop control, branches, the WalkState "
                "copies, the per-step control tests)\n", model, 100*model/usFull);
            out("    lean step function (no box%s): full loop %.3f us/walk (%.0f%% of runWalk's), walks alone %.3f us/walk; walks "
                "differing from runWalk's bit for bit: %lld of %lld\n", p.allAbsorb ? "; every face absorbs: no face tests" : "",
                1e6*p.secFullLean/double(p.walks), 100*p.secFullLean/p.secFull, 1e6*p.secWalksLean/double(p.walks), p.leanMismatch,
                p.walks);
            out("    shares of the modeled walk: closest %.0f%%, faceDist %.0f%%, sqrt %.0f%%, addBall %.0f%%, sphere %.0f%%, star %.0f%%, "
                "escape %.0f%%, per-walk parts %.0f%%\n", 100*1e-3*ev*p.nsClosest/model, 100*1e-3*ev*p.nsFaceDist/model,
                100*1e-3*ev*p.nsSqrt/model, 100*1e-3*ev*p.nsAddBall/model, 100*1e-3*mv*p.nsSphere/model,
                100*1e-3*(mv - 1)*p.nsStar/model, 100*1e-3*mv*p.nsOutside/model,
                100*1e-3*(p.nsFirstBall + p.nsMirror + p.nsInit + landedFrac*(p.nsFinish + p.nsValues + p.nsAccum))/model);
        }
    }
    // thread scaling: the same job on all threads (solveMulti) vs one thread
    out("\nthread scaling (cut tet, recipe run, all %zu points, N = 256, 16 batches): ", s.ps[size_t(tCut)].size());
    {
        const TetFrame& T = s.R.T[size_t(tCut)];
        const RingTetCrack& r = s.rt[size_t(tCut)];
        std::vector<double> sec;
        OwnRun o;
        for (int rep = 0; rep < 3; rep++) {
            o = runOwn(T, r.xc, allDirichlet(), ringRecipeFields(T), s.ps[size_t(tCut)].xi, 256, 16,
                       mixSeed(s.gp.seed, ownTag(s, tCut)), true, &r.gFront);
            sec.push_back(o.stats.seconds);
        }
        const double m = medianOf(sec);
        const WalkProfile p1 = profileJob(jobs[0], 100000, 256, fast, false);
        out("16 threads %.3f s = %.3f us/walk wall; one thread on cpu%d %.3f us/walk -> speed-up %.1fx\n", m,
            1e6*m/double(o.stats.walksRequested), fast, 1e6*p1.secFull/double(p1.walks),
            (1e6*p1.secFull/double(p1.walks))/(1e6*m/double(o.stats.walksRequested)));
    }
    // walks alone on all threads, runWalk vs the lean step function (median of 3): every tet's recipe-run walks, and the direct
    out("walks alone on 16 threads (no scoring), median of 3, runWalk vs lean step function:\n");
    {
        const WalkParams wpT = gateWalkParams(allDirichlet(), true, true);
        double sA = 0, sL = 0;
        long long walks = 0, mv = 0;
        for (int t = 0; t < kRingNT; t++) {
            const RingTetCrack& r = s.rt[size_t(t)];
            QueryPoints q;
            const std::vector<char> g(s.ps[size_t(t)].size(), 1);
            q.add(s.ps[size_t(t)].xi, &g, mixSeed(s.gp.seed, ownTag(s, t)));
            const Tet ref = referenceTet();
            std::vector<double> a, l;
            for (int rep = 0; rep < 3; rep++) {
                a.push_back(timeWalksParallel(&ref, nullptr, r.xc.strips, q, 640, wpT, false, &mv));
                l.push_back(timeWalksParallel(&ref, nullptr, r.xc.strips, q, 640, wpT, true, nullptr));
            }
            sA += medianOf(a);
            sL += medianOf(l);
            walks += (long long)q.size()*640;
        }
        out("  8 tets (all faces Dirichlet, odd mirror ball), N = 640: runWalk %.3f s (%.3f us/walk), lean %.3f s (%.3f us/walk): "
            "%.2fx\n", sA, 1e6*sA/double(walks), sL, 1e6*sL/double(walks), sA/sL);
        QueryPoints q;
        const std::vector<char> g(XO.size(), 1);
        q.add(XO, &g, mixSeed(s.gp.seed, 17500));
        std::vector<double> a, l;
        for (int rep = 0; rep < 3; rep++) {
            a.push_back(timeWalksParallel(nullptr, &s.O, scO, q, 64, ringReferenceWalkParams(s.O), false, &mv));
            l.push_back(timeWalksParallel(nullptr, &s.O, scO, q, 64, ringReferenceWalkParams(s.O), true, nullptr));
        }
        const double wO = double(q.size())*64.0;
        out("  direct on O (Neumann faces: face tests kept), N = 64: runWalk %.3f s (%.3f us/walk), lean %.3f s (%.3f us/walk): "
            "%.2fx\n", medianOf(a), 1e6*medianOf(a)/wO, medianOf(l), 1e6*medianOf(l)/wO, medianOf(a)/medianOf(l));
    }
}

// ---- part e: the current path, a control-variate prototype and the direct solve against the cached reference ----
void perfPartE(Report& out, int ring, int crack)
{
    PerfSetup s;
    perfSetup(s, ring, crack, gateParams(false).nRingPoints);
    const std::string cfg = "ring " + std::to_string(ring) + ", " + s.rc.name;
    out("\n=== (e) %s at check 4's point sets (%zu points + 7 probes): current path vs control-variate prototype vs direct, "
        "against the cached reference ===\n", cfg.c_str(), totalPoints(s));
    RefRun ref;
    auto t0 = std::chrono::steady_clock::now();
    if (!loadRingGT(s, ref)) {
        out("  the cached reference %s is missing or its key differs: skipped (nothing is solved or written)\n", ref.file.c_str());
        return;
    }
    out("  reference loaded from %s (%.1f s): N = %d per point in %d batches, original run %s\n", ref.file.c_str(),
        secondsSince(t0), ref.sums.walksPerBatch*ref.sums.nBatches, ref.sums.nBatches, ref.stats.summary().c_str());
    std::vector<size_t> off(kRingNT, 0);
    for (int t = 1; t < kRingNT; t++) off[size_t(t)] = off[size_t(t - 1)] + s.ps[size_t(t - 1)].size();
    const int m = 0;   // g = 1
    std::vector<FieldAt> refT;
    for (int t = 0; t < kRingNT; t++) refT.push_back(refField(ref, m, off[size_t(t)], s.ps[size_t(t)].size()));
    const FieldAt refP = refField(ref, m, ref.nA, ref.nF);
    double gt2 = 0, wSum = 0, gtH2 = 0;
    for (int t = 0; t < kRingNT; t++)
        for (size_t i = 0; i < s.ps[size_t(t)].size(); i++) {
            gt2 += s.ps[size_t(t)].w[i]*refT[size_t(t)].G[i].squaredNorm();
            gtH2 += s.ps[size_t(t)].w[i]*refT[size_t(t)].vH[i];
            wSum += s.ps[size_t(t)].w[i];
        }
    {
        ErrStats z;
        for (int t = 0; t < kRingNT; t++)
            for (size_t i = 0; i < s.ps[size_t(t)].size(); i++)
                z.add(s.ps[size_t(t)].w[i], s.ps[size_t(t)].X[i], 0.0, Vec3::Zero(), refT[size_t(t)].vH[i], 0.0, refT[size_t(t)].vG[i],
                      Vec3::Zero());
        out("  reference noise floor (g = 1, its batch s.e.): H RMS %.3f%% (max %.3f%%) of 2g; grad s.e. RMS / grad RMS %.2f%%; "
            "grad RMS %.4f\n", 100*z.seRefH(), 100*z.maxSeRefH/kTwoGMax, 100*std::sqrt(z.vRefG/std::max(1e-300, gt2)),
            std::sqrt(gt2/wSum));
    }
    std::vector<PerfResult> res;
    auto finish = [&](PerfResult& r, const std::vector<FieldAt>& fld, const Vec7* c, const Vec7* cse, const FieldAt* probes) {
        for (int t = 0; t < kRingNT; t++)
            for (size_t i = 0; i < s.ps[size_t(t)].size(); i++)
                r.e.add(s.ps[size_t(t)].w[i], s.ps[size_t(t)].X[i], fld[size_t(t)], i, refT[size_t(t)], i);
        r.gtGrad2 = gt2;
        for (int k = 0; k < kRingNV; k++) {
            const double v = c ? s.J[m](s.R.V[size_t(k)]) + (*c)[k] : probes->H[size_t(k)];
            r.vert[k] = v;
            r.vErr = std::max(r.vErr, std::fabs(v - refP.H[size_t(k)]));
        }
        if (cse) {
            r.cse = *cse;
            r.hasC = true;
        }
        res.push_back(r);
    };
    const int B = 16;
    // (1) the current path at N = 256: two runs per tet, physical Gram, glue, glued field with jackknife
    std::vector<OwnRun> rec(kRingNT), hc(kRingNT);
    {
        PerfResult r;
        r.name = "current, 2 runs/tet, N = 256";
        const std::vector<FieldSpec> fieldsC = ringLevelCFields(s.J);
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = s.R.T[size_t(t)];
            const RingTetCrack& rt = s.rt[size_t(t)];
            const Cut* gf = rt.cut() ? &rt.gFront : nullptr;
            rec[size_t(t)] = runOwn(T, rt.xc, allDirichlet(), ringRecipeFields(T), s.ps[size_t(t)].xi, 256, B,
                                    mixSeed(s.gp.seed, ownTag(s, t)), true, gf);
            hc[size_t(t)] = runOwn(T, rt.xc, allDirichlet(), fieldsC, s.ps[size_t(t)].xi, 256, B,
                                   mixSeed(s.gp.seed, ownTag(s, t) + 1), true, gf);
            for (const OwnRun* o : {&rec[size_t(t)], &hc[size_t(t)]}) {
                r.secWalks += o->stats.seconds;
                r.walks += o->stats.walksRequested;
                r.steps += o->stats.totalSteps;
            }
        }
        t0 = std::chrono::steady_clock::now();
        std::vector<Gram> G;
        for (int t = 0; t < kRingNT; t++) G.push_back(gramRing(rec[size_t(t)], hc[size_t(t)], s.ps[size_t(t)].w, s.R.T[size_t(t)], true));
        const Glue4 gl = glueRun4(s.R, G, kRingH0C + m);
        std::vector<FieldAt> fld;
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = s.R.T[size_t(t)];
            std::vector<Eigen::Vector4d> cJ;
            for (const Vec7& x : gl.jack) cJ.push_back(localValues(T, x));
            fld.push_back(decField4(T, rec[size_t(t)], hc[size_t(t)], 1, m, localValues(T, gl.c), &cJ));
        }
        r.secPost = secondsSince(t0);
        finish(r, fld, &gl.c, &gl.se, nullptr);
    }
    // (2) the control-variate prototype, several settings
    struct CvSetting {
        const char* name;
        bool split;
        int nQ, nR, B;
    };
    std::vector<CvRun> cvKeep;
    Vec7 cKeep = Vec7::Zero();
    for (const CvSetting cs : {CvSetting{"CV, 1 run/tet (plain), N = 256", false, 256, 256, 16},
                               CvSetting{"CV split: Q 256 (odd mirror), R 64", true, 256, 64, 16},
                               CvSetting{"CV split: Q 64 (odd mirror), R 32", true, 64, 32, 16},
                               CvSetting{"CV split: Q 32, R 16 (8 batches)", true, 32, 16, 8},
                               CvSetting{"CV split: Q 16, R 16 (8 batches)", true, 16, 16, 8}}) {
        PerfResult r;
        r.name = cs.name;
        std::vector<CvRun> cv;
        for (int t = 0; t < kRingNT; t++) {
            cv.push_back(runCv(s, t, false, cs.split, cs.nQ, cs.nR, cs.B, mixSeed(s.gp.seed, 21000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki) + 10*uint64_t(t))));
            r.secWalks += cv.back().seconds();
            r.walks += cv.back().walks();
            r.steps += cv.back().steps();
        }
        t0 = std::chrono::steady_clock::now();
        std::vector<Gram> G;
        for (int t = 0; t < kRingNT; t++) G.push_back(cvGram(cv[size_t(t)], s.ps[size_t(t)], s.R.T[size_t(t)], cs.B));
        const Glue4 gl = glueRun4(s.R, G, m);
        std::vector<FieldAt> fld;
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = s.R.T[size_t(t)];
            std::vector<Eigen::Vector4d> cJ;
            for (const Vec7& x : gl.jack) cJ.push_back(localValues(T, x));
            fld.push_back(cvField(cv[size_t(t)], s.ps[size_t(t)], T, m, localValues(T, gl.c), cJ, cs.B));
        }
        r.secPost = secondsSince(t0);
        finish(r, fld, &gl.c, &gl.se, nullptr);
        if (cvKeep.empty()) {
            out("  control-variate runs per tet (first setting): ");
            for (int t = 0; t < kRingNT; t++)
                out("%s%s %s/%s max|J-J_xi| %.1e%s", t ? "; " : "", s.R.T[size_t(t)].name, cv[size_t(t)].cut ? "cut" : "uncut",
                    cv[size_t(t)].walkR ? "R walks" : "R = 0", cv[size_t(t)].maxJdiff, cv[size_t(t)].hasQ || cv[size_t(t)].hasR ? "" : " (no walks)");
            out("\n");
            cvKeep = std::move(cv);
            cKeep = gl.c;
        }
    }
    // (3) the direct ring solve (seed independent of the reference's)
    for (const int ND : {256, 1024}) {
        PerfResult r;
        r.name = "direct on O, N = " + std::to_string(ND);
        const RefRun d = runDirect(s, ND, B, 17000 + 10*uint64_t(s.ri) + uint64_t(s.ki));
        r.secWalks = d.stats.seconds;
        r.walks = d.stats.walksRequested;
        r.steps = d.stats.totalSteps;
        t0 = std::chrono::steady_clock::now();
        std::vector<FieldAt> fld;
        for (int t = 0; t < kRingNT; t++) fld.push_back(refField(d, m, off[size_t(t)], s.ps[size_t(t)].size()));
        const FieldAt dp = refField(d, m, d.nA, d.nF);
        r.secPost = secondsSince(t0);
        finish(r, fld, nullptr, nullptr, &dp);
    }
    out("  errors vs the reference (g = 1), / 2g; gradRel = grad error RMS / GT grad RMS; vert max = max |H(v) - H_ref(probe v)| "
        "over the 7 vertices; fld se = the field's own s.e. RMS; nc = noise-corrected RMS (ref + field variances removed)\n");
    printPerfHeader(out);
    for (const PerfResult& r : res) printPerf(out, r);
    out("  vertex values v0..v6: reference");
    for (int k = 0; k < kRingNV; k++) out(" %+.4f", refP.H[size_t(k)]);
    out("\n");
    for (const PerfResult& r : res) {
        out("    %-34s", r.name.c_str());
        for (int k = 0; k < kRingNV; k++) out(" %+.4f", r.vert[k]);
        out("\n");
    }
    // the prototype targets the same quantities as the current path: H0 = J_xi + R, P_k = lambda_k + Q_k at every point,
    // within their noise (different walks), and the same glue up to noise
    {
        double zH = 0, zP = 0, dH2 = 0, vH2 = 0, dP2 = 0, vP2 = 0, jx = 0;
        long long nH = 0, nP = 0, bigH = 0, bigP = 0;
        for (int t = 0; t < kRingNT; t++) {
            const CvRun& c = cvKeep[size_t(t)];
            const PointSet& p = s.ps[size_t(t)];
            const OwnRun* rr = c.rRun();
            const OwnRun* qr = cvQRun(c);
            for (size_t i = 0; i < p.size(); i++) {
                double h, v, h2 = 0, v2 = 0;
                long long n;
                comboH(hc[size_t(t)].sums, i, unitCoef(2, 0), h, v, n);
                if (rr && c.iR1 >= 0) comboH(rr->sums, i, unitCoef(rr->sums.K, c.iR1), h2, v2, n);
                h2 += c.xf.J1(p.xi[i]);
                const double d = h2 - h, vv = v + v2;
                if (vv > 0) {
                    zH += d/std::sqrt(vv);
                    nH++;
                    bigH += std::fabs(d) > 3*std::sqrt(vv);
                }
                dH2 += p.w[i]*d*d;
                vH2 += p.w[i]*vv;
                if (ring == 1) jx = std::max(jx, std::fabs(c.xf.J1(p.xi[i]) - s.J[0](p.X[i])));
                for (int k = 0; k < 4; k++) {
                    double hp, vp, hq = 0, vq = 0;
                    comboH(rec[size_t(t)].sums, i, unitCoef(6, kRingP0 + k), hp, vp, n);
                    if (qr) comboH(qr->sums, i, unitCoef(qr->sums.K, c.iQ0 + k), hq, vq, n);
                    const double lam = k == 0 ? 1.0 - p.xi[i].sum() : p.xi[i][k - 1];
                    const double dd = lam + hq - hp, vv2 = vp + vq;
                    if (vv2 > 0) {
                        zP += dd/std::sqrt(vv2);
                        nP++;
                        bigP += std::fabs(dd) > 3*std::sqrt(vv2);
                    }
                    dP2 += p.w[i]*dd*dd;
                    vP2 += p.w[i]*vv2;
                }
            }
        }
        out("  same target check (CV N = 256 vs current N = 256, independent walks): H0: RMS diff / RMS noise = %.3f, pooled z %+.2f, "
            "|z| > 3 at %.2f%%; P_k: RMS diff / RMS noise = %.3f, pooled z %+.2f, |z| > 3 at %.2f%%; glued c max diff %.4f%s\n",
            std::sqrt(dH2/std::max(1e-300, vH2)), nH ? zH/std::sqrt(double(nH)) : 0.0, nH ? 100.0*double(bigH)/double(nH) : 0.0,
            std::sqrt(dP2/std::max(1e-300, vP2)), nP ? zP/std::sqrt(double(nP)) : 0.0, nP ? 100.0*double(bigP)/double(nP) : 0.0,
            (cKeep - res[0].vert + Vec7(s.J[0](s.R.V[0]), s.J[0](s.R.V[1]), s.J[0](s.R.V[2]), s.J[0](s.R.V[3]), s.J[0](s.R.V[4]),
                                        s.J[0](s.R.V[5]), s.J[0](s.R.V[6]))).cwiseAbs().maxCoeff(),
            ring == 1 ? (", ring 1 max |J_xi - J| at the points " + fmtd("%.1e", jx)).c_str() : "");
    }
}

// finite-difference check of J_xi's gradients (both g types) at random points of every tet
void perfCheckXiFront(Report& out, int ring, int crack)
{
    PerfSetup s;
    perfSetup(s, ring, crack, 256);
    double e1 = 0, e2 = 0, jd = 0;
    std::mt19937_64 rng(123);
    std::uniform_real_distribution<double> U(0, 1);
    for (int t = 0; t < kRingNT; t++) {
        const XiFront f = xiFront(s.R.T[size_t(t)], s.rc);
        for (int k = 0; k < 200; k++) {
            Vec3 xi(U(rng), U(rng), U(rng));
            if (xi.sum() > 1) continue;
            if (std::fabs(f.b(xi)) < 1e-3 || f.r(xi) < 1e-2) continue;
            const double h = 1e-6;
            Vec3 g1, g2;
            for (int a = 0; a < 3; a++) {
                Vec3 p = xi, q = xi;
                p[a] += h;
                q[a] -= h;
                g1[a] = (f.J1(p) - f.J1(q))/(2*h);
                g2[a] = (f.J2(p) - f.J2(q))/(2*h);
            }
            e1 = std::max(e1, (g1 - f.grad1(xi)).norm()/std::max(1e-12, f.grad1(xi).norm()));
            e2 = std::max(e2, (g2 - f.grad2(xi)).norm()/std::max(1e-12, f.grad2(xi).norm()));
            if (ring == 1) jd = std::max(jd, std::fabs(f.J1(xi) - s.J[0](s.R.T[size_t(t)].toX(xi))));
        }
    }
    out("J_xi check (ring %d, %s): finite-difference gradient max rel. error J1 %.1e, J2 %.1e%s\n", ring, s.rc.name.c_str(), e1, e2,
        ring == 1 ? (", max |J_xi - J| " + fmtd("%.1e", jd) + " (ring 1: J_T is a rotation)").c_str() : "");
}

void perfReview(Report& out, int ring, int crack, int N);   // part r / R (defined with the fast path, below)

} // namespace

bool runRingPerf(int ring, int crack, const std::string& parts)
{
    if (ring != 1 && ring != 2) throw std::invalid_argument("--ring-perf: ring 1 or 2");
    if (crack < 0 || crack > 2) throw std::invalid_argument("--ring-perf: crack (--config) 0, 1 or 2");
    Report out(false);
    const auto t0 = std::chrono::steady_clock::now();
    out("--ring-perf %s: ring %d, crack K%d; %u hardware threads\n", parts.c_str(), ring, crack + 1, std::thread::hardware_concurrency());
    for (char p : parts) {
        if (p == 'x') perfCheckXiFront(out, ring, crack);
        else if (p == 'a') perfPartA(out, ring, crack);
        else if (p == 'b') perfPartB(out, ring, crack);
        else if (p == 'c') perfPartC(out, ring, crack);
        else if (p == 'd') perfPartD(out, ring, crack);
        else if (p == 'e') perfPartE(out, ring, crack);
        else if (p == 'r') perfReview(out, ring, crack, 2560);
        else if (p == 'R') perfReview(out, ring, crack, 10240);
        else if (p == 'f') {
            for (int r = 1; r <= 2; r++)
                for (int k = 0; k < 3; k++) perfPartE(out, r, k);
        }
    }
    out("\n--ring-perf total %.1f s\n", secondsSince(t0));
    return true;
}

// =====================================================================================
// --ring-bench and the ring view's fast option
// =====================================================================================
// A, the fast per-tet path. Every tet solves in its own xi frame on the reference tet (Delta_xi H = 0, the crack pulled
// back into xi) with all four faces Dirichlet, as in check 4; the exactly known parts of its fields are taken out of the
// walks (control variates; --ring-perf parts e and f measured this design):
//   P_k = lambda_k + Q_k. Q_k: crack data -lambda_k, face data 0, odd-mirror first ball (harmonic extension 0). Q = 0
//         exactly where the crack misses the tet: no walks there.
//   H0  = J_xi + R. J_xi: Level C's J built in the tet's own xi frame (xiFront), exactly xi-harmonic off the crack and
//         +-1 on it (g = 1; the ramp's J_xi has no cap). R: crack data 0 (g = 1), face data J - J_xi, plain first ball.
//         R = 0 exactly in ring 1, where J_T is a rotation: no walks there.
// Then check 4's Gram estimator (physical energy), its 7 x 7 glue (glueRun4) and H = H0 + sum_v c_v P_v. The walks use
// solveMultiLean (the same walks bit for bit, without the walk box and, as every face absorbs, without the face tests);
// all tets' Q runs share one parallel loop and all R runs a second one.
// B, the direct ring solve: the reference's own solver (solveMulti on O, 8 Neumann faces, mirror first ball, crack +-g,
// both g types from the same walks) at N walks per point at the tets' points (H, grad H) and the vertex probes (H),
// with its own seeds, the same stop distance and solveMultiLean (no walk box; O's faces are Neumann, so the face tests
// stay).
// Opt-in: nothing here is called by check 4, the default ring view or any other existing mode. It reads gate_cache/
// (never writes it) and writes only into ring_perf/.
namespace {

// ---- A: the fast per-tet path ----
struct RingFastParams {
    int nQ = 32, nR = 16, nBatches = 8;
    double shellRel = kGateShellRel;   // walk stop distance / sqrt(2) (the reference tet's L): the gate's 1e-6 by default
    bool lean = true;                  // solveMultiLean (the same walks without box and face tests) instead of solveMulti
    bool fuse = true;                  // (lean only) all tets' Q runs in one parallel loop, all R runs in a second one
    bool ramp = false;                 // also H0 for the ramp g (R_2) and its glue (h = 1)
};
struct RingFast {
    std::vector<RingTetCrack> rt;
    std::vector<CvRun> cv;             // per tet: J_xi, which fields walk, the Q and R runs (prototype layout)
    std::vector<Gram> G;               // physical energy, layout H0 g1, H0 g2, P_0..P_3
    Glue4 gl[2];                       // h = 0 (g = 1), 1 (ramp, if asked)
    std::vector<FieldAt> F[2];         // glued H and physical grad H at every tet's points, jackknife variances
    double secSetup = 0, secQ = 0, secR = 0, secGram = 0, secGlue = 0, secField = 0, secTotal = 0;
    long long walksQ = 0, walksR = 0, stepsQ = 0, stepsR = 0, capped = 0, escaped = 0;
    int nQtets = 0, nRtets = 0;
};

// Which fields of a ring tet need walks: Q iff the crack cuts the tet; R iff J - J_xi is not 0 on the faces (sampled on a
// 25 x 25 barycentric grid per face, points within 1e-9 of the plane skipped: atan2's sign is round-off there) or, with
// the ramp, its crack data side (min(1, sqrt(r/L)) - sqrt(r/L)) is not 0 (r > L at a crack polygon vertex: r is convex).
void fastDecide(CvRun& c, const TetFrame& T, const RingTetCrack& rt, const RingCrack& rc, const FaceFn J[2], bool ramp)
{
    c.xf = xiFront(T, rc);
    c.cut = rt.cut();
    c.ramp = ramp;
    c.maxJdiff = 0;
    const Tet ref = referenceTet();
    for (int o = 0; o < 4; o++) {
        std::array<int, 3> fv{};
        int m = 0;
        for (int i = 0; i < 4; i++)
            if (i != o) fv[m++] = i;
        for (int i = 0; i <= 24; i++)
            for (int j = 0; i + j <= 24; j++) {
                const double l1 = i/24.0, l2 = j/24.0, l0 = 1.0 - l1 - l2;
                const Vec3 xi = l0*ref.v[fv[0]] + l1*ref.v[fv[1]] + l2*ref.v[fv[2]];
                if (std::fabs(c.xf.b(xi)) < 1e-9) continue;
                const Vec3 X = T.toX(xi);
                c.maxJdiff = std::max(c.maxJdiff, std::fabs(J[0](X) - c.xf.J1(xi)));
                if (ramp) c.maxJdiff = std::max(c.maxJdiff, std::fabs(J[1](X) - c.xf.J2(xi)));
            }
    }
    bool rampCrack = false;
    if (ramp && c.cut)
        for (const Vec3& v : rt.xc.cut.polygon) rampCrack = rampCrack || rt.gFront.distToFrontLine(v) > rt.gFront.L;
    c.walkR = c.maxJdiff > 1e-12 || rampCrack;
}

// The fast path on one ring: per tet the Q run (cut tets) and the R run (where needed), seeds mixSeed(seed, tagBase + 10 t)
// for Q and mixSeed(that, 77) for R (the prototype's), then the Grams, the glue and (withFields) the glued fields.
RingFast runRingFast(const Ring& R, const RingCrack& rc, const std::vector<PointSet>& ps, const FaceFn J[2],
                     const RingFastParams& p, uint64_t seed, uint64_t tagBase, bool withFields,
                     std::atomic<int>* progress = nullptr)
{
    using clk = std::chrono::steady_clock;
    const auto tStart = clk::now();
    const int B = p.nBatches;
    if (B < 3 || p.nQ%B || p.nR%B || p.nQ/B < 2 || p.nR/B < 2)
        throw std::invalid_argument("fast ring path: at least 3 batches, and nQ, nR multiples of the batch count with at least "
                                    "2 walks per batch");
    RingFast f;
    f.cv.resize(kRingNT);
    for (int t = 0; t < kRingNT; t++) {
        f.rt.push_back(ringTetCrack(R.T[size_t(t)], rc));
        fastDecide(f.cv[size_t(t)], R.T[size_t(t)], f.rt.back(), rc, J, p.ramp);
    }
    const Tet ref = referenceTet();
    WalkParams wpQ = gateWalkParams(allDirichlet(), true, true), wpR = gateWalkParams(allDirichlet(), true, false);
    wpQ.absorbShellRel = wpR.absorbShellRel = p.shellRel;
    const HarmonicExtension ext0 = [](const Vec3&, double* phi) {
        for (int k = 0; k < 4; k++) phi[k] = 0.0;
    };
    const TerminalValues tvQ = [](const WalkTerminal& w, double* u) {   // Q_k: crack -lambda_k(xi), faces 0
        if (w.face >= 0) {
            for (int k = 0; k < 4; k++) u[k] = 0.0;
            return;
        }
        const double l[4] = {1.0 - w.x.sum(), w.x[0], w.x[1], w.x[2]};
        for (int k = 0; k < 4; k++) u[k] = -l[k];
    };
    const int KR = p.ramp ? 2 : 1;
    std::vector<QueryPoints> qQ(kRingNT), qR(kRingNT);
    std::vector<TerminalValues> tvR(kRingNT);
    std::vector<MultiJob> jobsQ, jobsR;
    std::vector<int> tetQ, tetR;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        const RingTetCrack& rt = f.rt[size_t(t)];
        const CvRun& c = f.cv[size_t(t)];
        const std::vector<char> grad(ps[size_t(t)].size(), 1);
        const uint64_t sQ = mixSeed(seed, tagBase + 10*uint64_t(t));
        MultiJob jb;
        jb.tet = &ref;
        jb.cut = &rt.xc.strips;
        jb.nBatches = B;
        if (c.cut) {
            qQ[size_t(t)].add(ps[size_t(t)].xi, &grad, sQ);
            jb.pts = &qQ[size_t(t)];
            jb.nWalks = p.nQ;
            jb.K = 4;
            jb.wp = wpQ;
            jb.values = tvQ;
            jb.ext = &ext0;
            jobsQ.push_back(jb);
            tetQ.push_back(t);
        }
        if (c.walkR) {
            qR[size_t(t)].add(ps[size_t(t)].xi, &grad, mixSeed(sQ, 77));
            const Cut front = c.cut ? rt.gFront : rt.xc.strips.frontCut(rt.xc.strips.nStrips());
            tvR[size_t(t)] = [Tp = &T, xf = c.xf, J0 = J[0], J1 = J[1], front, ramp = p.ramp](const WalkTerminal& w, double* u) {
                if (w.face >= 0) {   // R's face data J - J_xi
                    const Vec3 X = Tp->toX(w.x);
                    u[0] = J0(X) - xf.J1(w.x);
                    if (ramp) u[1] = J1(X) - xf.J2(w.x);
                } else {             // its crack data: side g - J_xi = 0 (g = 1); side (min(1, sqrt(r/L)) - sqrt(r/L)) (ramp)
                    u[0] = 0.0;
                    if (ramp) u[1] = double(w.side)*(front.g(w.x, GType::SqrtFront) - std::sqrt(xf.r(w.x)/kSqrt2));
                }
            };
            jb.pts = &qR[size_t(t)];
            jb.nWalks = p.nR;
            jb.K = KR;
            jb.wp = wpR;
            jb.values = tvR[size_t(t)];
            jb.ext = nullptr;
            jobsR.push_back(jb);
            tetR.push_back(t);
        }
    }
    f.secSetup = secondsSince(tStart);
    auto run = [&](const std::vector<MultiJob>& jobs, const std::vector<int>& tets, bool isQ) {
        std::vector<MultiSums> sums(jobs.size());
        std::vector<SolveStats> st(jobs.size());
        const auto t0 = clk::now();
        if (p.lean && p.fuse) {
            if (!jobs.empty()) sums = solveMultiLeanJobs(jobs, &st);
        } else {
            for (size_t j = 0; j < jobs.size(); j++) {
                const MultiJob& b = jobs[j];
                sums[j] = p.lean ? solveMultiLean(*b.tet, *b.cut, *b.pts, b.nWalks, b.nBatches, b.wp, b.K, b.values, &st[j], b.ext)
                                 : solveMulti(*b.tet, *b.cut, *b.pts, b.nWalks, b.nBatches, b.wp, b.K, b.values, &st[j], b.ext);
            }
        }
        const double sec = secondsSince(t0);
        for (size_t j = 0; j < jobs.size(); j++) {
            CvRun& c = f.cv[size_t(tets[j])];
            OwnRun& o = isQ ? c.runQ : c.runR;
            o.sums = std::move(sums[j]);
            o.stats = st[j];
            (isQ ? f.walksQ : f.walksR) += st[j].walksRequested;
            (isQ ? f.stepsQ : f.stepsR) += st[j].totalSteps;
            f.capped += st[j].walksCapped;
            f.escaped += st[j].walksEscaped;
            if (isQ) {
                c.hasQ = true;
                c.iQ0 = 0;
            } else {
                c.hasR = true;
                c.iR1 = 0;
                c.iR2 = p.ramp ? 1 : -1;
            }
        }
        return sec;
    };
    f.secQ = run(jobsQ, tetQ, true);
    if (progress) progress->store(4);
    f.secR = run(jobsR, tetR, false);
    if (progress) progress->store(8);
    f.nQtets = int(jobsQ.size());
    f.nRtets = int(jobsR.size());
    auto t0 = clk::now();
    for (int t = 0; t < kRingNT; t++) f.G.push_back(cvGram(f.cv[size_t(t)], ps[size_t(t)], R.T[size_t(t)], B));
    f.secGram = secondsSince(t0);
    t0 = clk::now();
    const int nModes = p.ramp ? 2 : 1;
    for (int m = 0; m < nModes; m++) f.gl[m] = glueRun4(R, f.G, m);
    f.secGlue = secondsSince(t0);
    if (withFields) {
        t0 = clk::now();
        for (int m = 0; m < nModes; m++)
            for (int t = 0; t < kRingNT; t++) {
                const TetFrame& T = R.T[size_t(t)];
                std::vector<Eigen::Vector4d> cJ;
                for (const Vec7& x : f.gl[m].jack) cJ.push_back(localValues(T, x));
                f.F[m].push_back(cvField(f.cv[size_t(t)], ps[size_t(t)], T, m, localValues(T, f.gl[m].c), cJ, B));
            }
        f.secField = secondsSince(t0);
    }
    f.secTotal = secondsSince(tStart);
    return f;
}

// ---- B: the direct ring solve ----
struct RingDirect {
    RefRun run;                        // the sums at the tets' points and the probes (as the reference's)
    std::vector<FieldAt> F[2];         // per tet, per g type: H, grad H and their variances (refField)
    FieldAt P[2];                      // at the vertex probes
    double secWalks = 0, secField = 0, secTotal = 0;
};
RingDirect runRingDirect(const ConvexDomain& O, const Ring& R, const RingCrack& rc, const std::vector<PointSet>& ps, int N, int B,
                         uint64_t seed, double shellRel, bool lean, int nModes, bool withFields)
{
    using clk = std::chrono::steady_clock;
    const auto tStart = clk::now();
    RingDirect d;
    RefRun& r = d.run;
    for (const PointSet& p : ps) r.X.insert(r.X.end(), p.X.begin(), p.X.end());
    r.nA = r.X.size();
    const std::vector<Vec3> probes = ringProbes(R);
    r.X.insert(r.X.end(), probes.begin(), probes.end());
    r.nF = probes.size();
    const StripCut sc = ringStrips(rc);
    WalkParams wp = ringReferenceWalkParams(O);
    if (shellRel != kGateShellRel) wp.absorbShellRel = shellRel*kSqrt2/O.L;
    const Cut front = sc.frontCut(sc.nStrips());
    const TerminalValues tv = [&front](const WalkTerminal& t, double* u) {
        u[0] = double(t.side)*front.g(t.x, GType::One);
        u[1] = double(t.side)*front.g(t.x, GType::SqrtFront);
    };
    std::vector<char> grad(r.X.size(), 1);
    std::fill(grad.begin() + std::ptrdiff_t(r.nA), grad.end(), char(0));
    QueryPoints q;
    q.add(r.X, &grad, seed);
    r.sums = lean ? solveMultiLean(O, sc, q, N, B, wp, 2, tv, &r.stats) : solveMulti(O, sc, q, N, B, wp, 2, tv, &r.stats);
    d.secWalks = secondsSince(tStart);
    if (withFields) {
        const auto t0 = clk::now();
        for (int m = 0; m < nModes; m++) {
            size_t off = 0;
            for (const PointSet& p : ps) {
                d.F[m].push_back(refField(r, m, off, p.size()));
                off += p.size();
            }
            d.P[m] = refField(r, m, r.nA, r.nF);
        }
        d.secField = secondsSince(t0);
    }
    d.secTotal = secondsSince(tStart);
    return d;
}

// ---- check 4's current path (two runs per tet, runWalk), the baseline; at N = 2560, 64 batches and tagBase = 12000 +
// 1000 ri + 100 ki it is check 4's own run, walk for walk ----
struct RingCurrent {
    std::vector<OwnRun> rec, hc;
    std::vector<Gram> G;
    Glue4 gl[2];
    std::vector<FieldAt> F[2];
    double secSetup = 0, secWalks = 0, secGram = 0, secGlue = 0, secField = 0, secTotal = 0;
    long long walks = 0, steps = 0, capped = 0, escaped = 0;
};
RingCurrent runRingCurrent(const Ring& R, const RingCrack& rc, const std::vector<PointSet>& ps, const FaceFn J[2], int N, int B,
                           uint64_t seed, uint64_t tagBase, int nModes, bool withFields)
{
    using clk = std::chrono::steady_clock;
    const auto tStart = clk::now();
    RingCurrent c;
    std::vector<RingTetCrack> rt;
    for (const TetFrame& T : R.T) rt.push_back(ringTetCrack(T, rc));
    const std::vector<FieldSpec> fieldsC = ringLevelCFields(J);
    c.secSetup = secondsSince(tStart);
    auto t0 = clk::now();
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        const RingTetCrack& r = rt[size_t(t)];
        const Cut* gf = r.cut() ? &r.gFront : nullptr;
        const uint64_t tag = tagBase + 10*uint64_t(t);
        c.rec.push_back(runOwn(T, r.xc, allDirichlet(), ringRecipeFields(T), ps[size_t(t)].xi, N, B, mixSeed(seed, tag), true, gf));
        c.hc.push_back(runOwn(T, r.xc, allDirichlet(), fieldsC, ps[size_t(t)].xi, N, B, mixSeed(seed, tag + 1), true, gf));
        for (const OwnRun* o : {&c.rec.back(), &c.hc.back()}) {
            c.walks += o->stats.walksRequested;
            c.steps += o->stats.totalSteps;
            c.capped += o->stats.walksCapped;
            c.escaped += o->stats.walksEscaped;
        }
    }
    c.secWalks = secondsSince(t0);
    t0 = clk::now();
    for (int t = 0; t < kRingNT; t++) c.G.push_back(gramRing(c.rec[size_t(t)], c.hc[size_t(t)], ps[size_t(t)].w, R.T[size_t(t)], true));
    c.secGram = secondsSince(t0);
    t0 = clk::now();
    for (int m = 0; m < nModes; m++) c.gl[m] = glueRun4(R, c.G, kRingH0C + m);
    c.secGlue = secondsSince(t0);
    if (withFields) {
        t0 = clk::now();
        for (int m = 0; m < nModes; m++)
            for (int t = 0; t < kRingNT; t++) {
                const TetFrame& T = R.T[size_t(t)];
                std::vector<Eigen::Vector4d> cJ;
                for (const Vec7& x : c.gl[m].jack) cJ.push_back(localValues(T, x));
                c.F[m].push_back(decField4(T, c.rec[size_t(t)], c.hc[size_t(t)], 1, m, localValues(T, c.gl[m].c), &cJ));
            }
        c.secField = secondsSince(t0);
    }
    c.secTotal = secondsSince(tStart);
    return c;
}

// the glued energy (1/2) sum_T int |grad H|^2 from the tets' Grams (layout index h of H0) and the glue, with its jackknife s.e.
void glueEnergy(const Ring& R, const std::vector<Gram>& G, int h, const Glue4& gl, double& E, double& se)
{
    E = 0;
    for (int t = 0; t < kRingNT; t++) E += tetEnergy4(G[size_t(t)].G, h, localValues(R.T[size_t(t)], gl.c));
    std::vector<double> rep;
    for (size_t a = 0; a < gl.jack.size(); a++) {
        double e = 0;
        for (int t = 0; t < kRingNT; t++) e += tetEnergy4(G[size_t(t)].jack[a], h, localValues(R.T[size_t(t)], gl.jack[a]));
        rep.push_back(e);
    }
    se = rep.size() >= 3 ? jackSe(rep) : 0.0;
}

// --ring-perf r (N = 2560) / R (N = 10240), review only: the fast path (P_k = lambda_k + Q_k, H0 = J_xi + R, the skipped
// fields, lean walks) against check 4's two runs per tet (runWalk) at the same points (1024 per tet asked) and the same N
// in 64 batches, independent seeds. Per tet the Gram entries the glue uses (S_T, b_T and the H0 energy entry); the glued
// vertex values; per point (points within kCompareBand of a face left out) P_k, H0 and the glued H, value and xi / physical
// gradient. Each difference over its combined s.e. is a z: with the same targets the pooled z is N(0, 1) and mean z^2 ~ 1.
// Uncut tets: the fast S_T is exact, (sum w / V_T) K_T. Prints only.
void perfReview(Report& out, int ring, int crack, int N)
{
    PerfSetup s;
    perfSetup(s, ring, crack, 1024);
    const int B = 64;
    const int nModes = ring == 1 && !s.rc.through ? 2 : 1;
    RingFastParams fp;
    fp.nQ = fp.nR = N;
    fp.nBatches = B;
    fp.ramp = nModes == 2;
    auto t0 = std::chrono::steady_clock::now();
    const RingFast f = runRingFast(s.R, s.rc, s.ps, s.J, fp, s.gp.seed, 21000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki) + 300000, true);
    const double tf = secondsSince(t0);
    t0 = std::chrono::steady_clock::now();
    const RingCurrent cur = runRingCurrent(s.R, s.rc, s.ps, s.J, N, B, s.gp.seed, 12000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki) + 300000,
                                           nModes, true);
    const double tc = secondsSince(t0);
    out("\n=== review: fast path vs check 4's current path, ring %d %s, %zu points, N = %d in %d batches (fast %.1f s, current %.1f s)\n",
        ring, s.rc.name.c_str(), totalPoints(s), N, B, tf, tc);
    auto line = [&](const char* what, const ZStats& z) {
        out("  %-44s n %6lld: pooled z %+6.2f, mean z %+6.3f (bound %.3f), mean z^2 %.3f +- %.3f, |z| < 3 at %.2f%%, max |z| %.2f\n", what,
            z.n, z.pooled(), z.mean(), z.meanBound(), z.meanZ2(), z.meanZ2Se(), 100*z.frac(), z.maxAbsZ);
    };
    // 1. the glue's inputs
    ZStats zS, zSu, zB[2], zA[2];
    double maxKrel = 0, maxWrel = 0;
    int nUncut = 0;
    for (int t = 0; t < kRingNT; t++) {
        const Gram &gf = f.G[size_t(t)], &go = cur.G[size_t(t)];
        const bool cutT = f.cv[size_t(t)].cut;
        for (int i = 0; i < 4; i++)
            for (int j = i; j < 4; j++) {
                const int a = kRingP0 + i, b = kRingP0 + j;
                (cutT ? zS : zSu).add(size_t(t), gf.G(a, b) - go.G(a, b), gf.se(a, b)*gf.se(a, b) + go.se(a, b)*go.se(a, b));
            }
        for (int m = 0; m < nModes; m++) {
            for (int i = 0; i < 4; i++) {
                const int a = kRingP0 + i;
                zB[m].add(size_t(t), gf.G(m, a) - go.G(kRingH0C + m, a), gf.se(m, a)*gf.se(m, a) + go.se(kRingH0C + m, a)*go.se(kRingH0C + m, a));
            }
            zA[m].add(size_t(t), gf.G(m, m) - go.G(kRingH0C + m, kRingH0C + m),
                      gf.se(m, m)*gf.se(m, m) + go.se(kRingH0C + m, kRingH0C + m)*go.se(kRingH0C + m, kRingH0C + m));
        }
        if (!cutT) {   // exact: S_T = sum w gradN_i . gradN_j = (sum w / V_T) K_T
            nUncut++;
            const TetFrame& T = s.R.T[size_t(t)];
            const Eigen::Matrix4d K = stiffness(T);
            double sw = 0;
            for (double w : s.ps[size_t(t)].w) sw += w;
            const double V = std::fabs(T.detJ)/6.0;
            maxWrel = std::max(maxWrel, std::fabs(sw/V - 1.0));
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    maxKrel = std::max(maxKrel, std::fabs(gf.G(kRingP0 + i, kRingP0 + j) - sw/V*K(i, j))/K.cwiseAbs().maxCoeff());
        }
    }
    line("Gram S_T entries, cut tets", zS);
    if (zSu.n) line("Gram S_T entries, uncut tets (fast: exact)", zSu);
    for (int m = 0; m < nModes; m++) {
        line(m ? "Gram b_T entries (ramp)" : "Gram b_T entries (g = 1)", zB[m]);
        line(m ? "Gram H0.H0 entry (ramp)" : "Gram H0.H0 entry (g = 1)", zA[m]);
    }
    if (nUncut)
        out("  uncut tets (%d): fast S_T vs (sum w / V_T) K_T max rel. difference %.1e; |sum w / V_T - 1| max %.1e\n", nUncut, maxKrel, maxWrel);
    for (int m = 0; m < nModes; m++) {
        std::string zs;
        double z2 = 0;
        for (int k = 0; k < kRingNV; k++) {
            const double se = std::hypot(f.gl[m].se[k], cur.gl[m].se[k]), z = se > 0 ? (f.gl[m].c[k] - cur.gl[m].c[k])/se : 0.0;
            zs += fmtd(" %+.2f", z);
            z2 += z*z/kRingNV;
        }
        out("  glued c (%s): fast %s\n                 current %s\n    z per vertex%s, mean z^2 %.2f\n", m ? "ramp" : "g = 1",
            fmtVec(f.gl[m].c).c_str(), fmtVec(cur.gl[m].c).c_str(), zs.c_str(), z2);
    }
    // 2. per point
    ZStats zPv[4], zPg[4], zHv[2], zHg[2], zGv[2], zGg[2];
    auto gradStats = [](const MultiSums& sm, size_t i, const std::vector<double>& coef, Vec3& g, Vec3& v) {
        const int Bn = sm.nBatches;
        std::vector<Vec3> gb(static_cast<size_t>(Bn));
        g = Vec3::Zero();
        for (int b = 0; b < Bn; b++) g += (gb[size_t(b)] = comboBatchGrad(sm, i, b, coef))/double(Bn);
        v = Vec3::Zero();
        for (const Vec3& x : gb) v += (x - g).cwiseAbs2()/(double(Bn)*double(Bn - 1));
    };
    const Vec3 gLam[4] = {Vec3(-1, -1, -1), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    for (int t = 0; t < kRingNT; t++) {
        const CvRun& c = f.cv[size_t(t)];
        const PointSet& p = s.ps[size_t(t)];
        for (size_t i = 0; i < p.size(); i++) {
            const Vec3& xi = p.xi[i];
            if (dirichletDist(allDirichlet(), xi) < kCompareBand) continue;
            const double lam[4] = {1.0 - xi.sum(), xi[0], xi[1], xi[2]};
            double h, v, hq, vq;
            long long n;
            Vec3 g, gv, gq, gqv;
            for (int k = 0; k < 4; k++) {
                comboH(cur.rec[size_t(t)].sums, i, unitCoef(6, kRingP0 + k), h, v, n);
                gradStats(cur.rec[size_t(t)].sums, i, unitCoef(6, kRingP0 + k), g, gv);
                hq = vq = 0;
                gq = gqv = Vec3::Zero();
                if (c.hasQ) {
                    comboH(c.runQ.sums, i, unitCoef(4, k), hq, vq, n);
                    gradStats(c.runQ.sums, i, unitCoef(4, k), gq, gqv);
                }
                zPv[k].add(i, lam[k] + hq - h, v + vq);
                for (int a = 0; a < 3; a++) zPg[k].add(i, gLam[k][a] + gq[a] - g[a], gv[a] + gqv[a]);
            }
            for (int m = 0; m < nModes; m++) {
                comboH(cur.hc[size_t(t)].sums, i, unitCoef(2, m), h, v, n);
                gradStats(cur.hc[size_t(t)].sums, i, unitCoef(2, m), g, gv);
                const int iR = m ? c.iR2 : c.iR1;
                hq = vq = 0;
                gq = gqv = Vec3::Zero();
                if (c.hasR && iR >= 0) {
                    comboH(c.runR.sums, i, unitCoef(c.runR.sums.K, iR), hq, vq, n);
                    gradStats(c.runR.sums, i, unitCoef(c.runR.sums.K, iR), gq, gqv);
                }
                const double J = m ? c.xf.J2(xi) : c.xf.J1(xi);
                const Vec3 gJ = m ? c.xf.grad2(xi) : c.xf.grad1(xi);
                zHv[m].add(i, J + hq - h, v + vq);
                for (int a = 0; a < 3; a++) zHg[m].add(i, gJ[a] + gq[a] - g[a], gv[a] + gqv[a]);
                const FieldAt &Ff = f.F[m][size_t(t)], &Fo = cur.F[m][size_t(t)];
                zGv[m].add(i, Ff.H[i] - Fo.H[i], Ff.vH[i] + Fo.vH[i]);
                for (int a = 0; a < 3; a++) zGg[m].add(i, Ff.G[i][a] - Fo.G[i][a], Ff.vG[i][a] + Fo.vG[i][a]);
            }
        }
    }
    for (int k = 0; k < 4; k++) {
        line(("P_" + std::to_string(k) + " value: lambda + Q vs P").c_str(), zPv[k]);
        line(("P_" + std::to_string(k) + " xi gradient (3 components)").c_str(), zPg[k]);
    }
    for (int m = 0; m < nModes; m++) {
        line(m ? "H0 value (ramp): J_xi + R vs H0" : "H0 value (g = 1): J_xi + R vs H0", zHv[m]);
        line(m ? "H0 xi gradient (ramp)" : "H0 xi gradient (g = 1)", zHg[m]);
        line(m ? "glued H value (ramp)" : "glued H value (g = 1)", zGv[m]);
        line(m ? "glued H physical gradient (ramp)" : "glued H physical gradient (g = 1)", zGg[m]);
    }
}

// ---- the benchmark: one configuration ----
constexpr int kBR = 5;   // regions: all, near the front line (r < 0.1), away from it, near interior faces, near crack-face lines
const char* const kBRName[kBR] = {"all", "front_near", "front_far", "iface", "lines"};

struct BenchCfg {
    PerfSetup s;                                  // ring, crack, check 4's point sets (4096 per tet), J
    RefRun ref;
    int nModes = 1;
    std::vector<size_t> off;
    size_t nPts = 0;
    std::vector<double> w, wLine;                 // per point, all tets in order
    std::vector<std::vector<double>> wLineT;      // per tet
    std::vector<std::array<char, kBR>> reg;       // per point
    FieldAt refAll[2], refP[2];
    double gt2[2][kBR] = {}, Eref[2] = {}, ErefSe[2] = {}, ELref[2] = {}, ELrefSe[2] = {};
    ErrStats floorErr[2];                         // the reference's noise alone (error 0)
    std::string name;
};

// the energy (1/2) int |grad H|^2 over weights w of a solve in X (the reference, the direct) from its batch gradients, s.e.
void solveEnergy(const MultiSums& sums, const std::vector<double>& w, int m, double& E, double& se)
{
    const std::vector<double> unit2[2] = {unitCoef(2, 0), unitCoef(2, 1)};
    const Gram g = gramFromBatches(w, sums.nBatches, 2, [&](size_t p, int b, Vec3* gr) {
        for (int q = 0; q < 2; q++) gr[q] = comboBatchGrad(sums, p, b, unit2[q]);
    });
    E = 0.5*g.G(m, m);
    se = 0.5*g.se(m, m);
}

bool benchSetup(Report& out, BenchCfg& b, int ring, int crack)
{
    perfSetup(b.s, ring, crack, gateParams(false).nRingPoints);
    b.name = "ring " + std::to_string(ring) + " " + b.s.rc.name;
    const auto t0 = std::chrono::steady_clock::now();
    if (!loadRingGT(b.s, b.ref)) {
        out("  %s: the cached reference %s is missing or its key differs: configuration skipped (nothing is solved or written)\n",
            b.name.c_str(), b.ref.file.c_str());
        return false;
    }
    out("  reference %s loaded read-only (%.1f s): N = %d per point in %d batches; original run: %s\n", b.ref.file.c_str(),
        secondsSince(t0), b.ref.sums.walksPerBatch*b.ref.sums.nBatches, b.ref.sums.nBatches, b.ref.stats.summary().c_str());
    const PerfSetup& s = b.s;
    b.nModes = ring == 1 && !s.rc.through ? 2 : 1;
    b.off.assign(kRingNT, 0);
    for (int t = 1; t < kRingNT; t++) b.off[size_t(t)] = b.off[size_t(t - 1)] + s.ps[size_t(t - 1)].size();
    b.nPts = totalPoints(s);
    std::vector<std::array<Vec3, 2>> lines;
    for (const RingTetCrack& r : s.rt) lines.insert(lines.end(), r.faceLines.begin(), r.faceLines.end());
    const Cut& co = s.rc.cutO;
    b.wLineT.assign(kRingNT, {});
    for (int t = 0; t < kRingNT; t++) {
        const PointSet& p = s.ps[size_t(t)];
        for (size_t i = 0; i < p.size(); i++) {
            const Vec3& X = p.X[i];
            const double r = s.rc.through ? kInf : std::hypot(co.u.dot(X) - s.rc.sJ, s.rc.n.dot(X) - s.rc.d);
            bool nearI = false, nearL = false;
            for (const auto& f : s.R.innerFaces)
                nearI = nearI || pointTriangleDist(X, s.R.V[size_t(f[0])], s.R.V[size_t(f[1])], s.R.V[size_t(f[2])]) < kNearRing;
            for (const auto& l : lines) nearL = nearL || pointSegmentDist(X, l[0], l[1]) < kNearRing;
            b.reg.push_back({1, char(r < kNearRing), char(!(r < kNearRing)), char(nearI), char(nearL)});
            b.w.push_back(p.w[i]);
            b.wLine.push_back(nearL ? p.w[i] : 0.0);
            b.wLineT[size_t(t)].push_back(nearL ? p.w[i] : 0.0);
        }
    }
    for (int m = 0; m < b.nModes; m++) {
        b.refAll[m] = refField(b.ref, m, 0, b.ref.nA);
        b.refP[m] = refField(b.ref, m, b.ref.nA, b.ref.nF);
        for (size_t g = 0; g < b.nPts; g++) {
            for (int k = 0; k < kBR; k++)
                if (b.reg[g][k]) b.gt2[m][k] += b.w[g]*b.refAll[m].G[g].squaredNorm();
            b.floorErr[m].add(b.w[g], Vec3::Zero(), 0.0, Vec3::Zero(), b.refAll[m].vH[g], 0.0, b.refAll[m].vG[g], Vec3::Zero());
        }
        std::vector<double> wA(b.w), wL(b.wLine);   // the probes carry no gradient: weight 0
        wA.resize(b.ref.X.size(), 0.0);
        wL.resize(b.ref.X.size(), 0.0);
        solveEnergy(b.ref.sums, wA, m, b.Eref[m], b.ErefSe[m]);
        solveEnergy(b.ref.sums, wL, m, b.ELref[m], b.ELrefSe[m]);
        out("  %s, %s: reference noise floor (its batch s.e.): H RMS %.3f%% (max %.3f%%) of 2g; grad s.e. RMS / grad RMS %.2f%%;\n"
            "    energy over the ring minus the eps-tube %.5f +- %.5f, within %.1f of the crack-face lines %.5f +- %.5f\n",
            b.name.c_str(), m ? "ramp" : "g = 1", 100*b.floorErr[m].seRefH(), 100*b.floorErr[m].maxSeRefH/kTwoGMax,
            100*std::sqrt(b.floorErr[m].vRefG/std::max(1e-300, b.gt2[m][0])), b.Eref[m], b.ErefSe[m], kNearRing, b.ELref[m],
            b.ELrefSe[m]);
    }
    return true;
}

// One measured row: a method at one setting, one g type, one seed set
struct BenchRow {
    std::string config, method, family, label;
    int ring = 0, crack = 0, g = 0, seedRep = 0, nQ = 0, nR = 0, N = 0, B = 0;
    double shellRel = 0;
    bool lean = false, fuse = false;
    long long points = 0, walks = 0, moves = 0, capped = 0, escaped = 0;
    std::vector<double> times;   // total wall time per repeat
    double tMed = 0, tSetup = 0, tQ = 0, tR = 0, tWalks = 0, tGram = 0, tGlue = 0, tField = 0;
    ErrStats e[kBR];
    double gt2[kBR] = {};
    double E = 0, Ese = 0, Eref = 0, ErefSe = 0, EL = 0, ELse = 0, ELref = 0, ELrefSe = 0;
    bool hasE = false;
    Vec7 vert = Vec7::Zero(), vref = Vec7::Zero(), cse = Vec7::Zero();
    bool hasC = false;
    double vErr = 0;
    int nQtets = -1, nRtets = -1;
};

// errors of a method's per-tet fields (g type m) against the reference, per region; the vertex values; the energies
void fillErrors(const BenchCfg& b, int m, const std::vector<FieldAt>& F, BenchRow& row)
{
    const FieldAt& r = b.refAll[m];
    for (int t = 0; t < kRingNT; t++) {
        const FieldAt& f = F[size_t(t)];
        const PointSet& p = b.s.ps[size_t(t)];
        for (size_t i = 0; i < p.size(); i++) {
            const size_t g = b.off[size_t(t)] + i;
            for (int k = 0; k < kBR; k++)
                if (b.reg[g][k]) row.e[k].add(b.w[g], p.X[i], f.H[i] - r.H[g], f.G[i] - r.G[g], r.vH[g], f.vH[i], r.vG[g], f.vG[i]);
        }
    }
    for (int k = 0; k < kBR; k++) row.gt2[k] = b.gt2[m][k];
    for (int k = 0; k < kRingNV; k++) row.vref[k] = b.refP[m].H[size_t(k)];
    row.vErr = (row.vert - row.vref).cwiseAbs().maxCoeff();
    row.Eref = b.Eref[m];
    row.ErefSe = b.ErefSe[m];
    row.ELref = b.ELref[m];
    row.ELrefSe = b.ELrefSe[m];
    row.hasE = true;
}

double gradRel(const ErrStats& e, double gt2) { return gt2 > 1e-12 ? std::sqrt(e.eG2/gt2) : std::nan(""); }
double gradRelNc(const ErrStats& e, double gt2)
{
    return gt2 > 1e-12 ? std::sqrt(std::max(0.0, e.eG2 - e.vRefG - e.vFldG)/gt2) : std::nan("");
}

std::string csvNum(double v, const char* f = "%.6g")
{
    if (!std::isfinite(v)) return "";
    char buf[64];
    std::snprintf(buf, sizeof buf, f, v);
    return buf;
}

std::string csvHeader()
{
    std::string h = "config,ring,crack,g,method,family,label,shell_rel,lean,fuse,seed_rep,nQ,nR,N,B,points,tets_Q,tets_R,"
                    "t_total_median,t_total_min,t_total_max,n_repeats,t_setup,t_walks_Q,t_walks_R,t_walks,t_gram,t_glue,t_field,"
                    "walks,moves,mean_walk_len,capped,escaped,"
                    "H_rms,H_nc,H_max,H_se_fld,H_se_ref,grad_rel,grad_rel_nc,grad_abs,grad_abs_nc,grad_se_ref_rel";
    for (int k = 1; k < kBR; k++)
        h += std::string(",H_rms_") + kBRName[k] + ",H_nc_" + kBRName[k] + ",grad_rel_" + kBRName[k] + ",grad_rel_nc_" + kBRName[k];
    h += ",E,E_se,E_ref,E_ref_se,E_relerr,E_relerr_se,E_lines,E_lines_se,E_lines_ref,E_lines_ref_se,E_lines_relerr,"
         "E_lines_relerr_se";
    for (int k = 0; k < kRingNV; k++) h += ",v" + std::to_string(k);
    for (int k = 0; k < kRingNV; k++) h += ",vref" + std::to_string(k);
    h += ",vert_maxerr,c_se_max\n";
    return h;
}

std::string csvRow(const BenchRow& r)
{
    std::vector<std::string> c;
    auto add = [&](const std::string& s) { c.push_back(s); };
    auto num = [&](double v, const char* f = "%.6g") { c.push_back(csvNum(v, f)); };
    auto itg = [&](long long v) { c.push_back(std::to_string(v)); };
    add(r.config);
    itg(r.ring);
    add("K" + std::to_string(r.crack + 1));
    add(r.g ? "ramp" : "g1");
    add(r.method);
    add(r.family);
    add(r.label);
    num(r.shellRel, "%.3g");
    itg(r.lean);
    itg(r.fuse);
    itg(r.seedRep);
    itg(r.nQ);
    itg(r.nR);
    itg(r.N);
    itg(r.B);
    itg(r.points);
    if (r.nQtets >= 0) itg(r.nQtets); else add("");
    if (r.nRtets >= 0) itg(r.nRtets); else add("");
    const double tMin = r.times.empty() ? 0.0 : *std::min_element(r.times.begin(), r.times.end());
    const double tMax = r.times.empty() ? 0.0 : *std::max_element(r.times.begin(), r.times.end());
    num(r.tMed, "%.5f");
    num(tMin, "%.5f");
    num(tMax, "%.5f");
    itg((long long)r.times.size());
    num(r.tSetup, "%.5f");
    num(r.tQ, "%.5f");
    num(r.tR, "%.5f");
    num(r.tWalks + r.tQ + r.tR, "%.5f");   // all walks (ours: Q + R)
    num(r.tGram, "%.5f");
    num(r.tGlue, "%.6f");
    num(r.tField, "%.5f");
    itg(r.walks);
    itg(r.moves);
    num(r.walks ? double(r.moves)/double(r.walks - r.escaped) : 0.0, "%.3f");
    itg(r.capped);
    itg(r.escaped);
    const ErrStats& e = r.e[0];
    num(e.rmsH());
    num(e.corrH());
    num(e.maxH/kTwoGMax);
    num(e.seFldH());
    num(e.seRefH());
    num(gradRel(e, r.gt2[0]));
    num(gradRelNc(e, r.gt2[0]));
    num(e.rmsG());
    num(e.corrG());
    num(r.gt2[0] > 1e-12 ? std::sqrt(e.vRefG/r.gt2[0]) : std::nan(""));
    for (int k = 1; k < kBR; k++) {
        num(r.e[k].w > 0 ? r.e[k].rmsH() : std::nan(""));
        num(r.e[k].w > 0 ? r.e[k].corrH() : std::nan(""));
        num(r.e[k].w > 0 ? gradRel(r.e[k], r.gt2[k]) : std::nan(""));
        num(r.e[k].w > 0 ? gradRelNc(r.e[k], r.gt2[k]) : std::nan(""));
    }
    auto rel = [](double E, double Er) { return std::fabs(Er) > 1e-300 ? (E - Er)/Er : std::nan(""); };
    auto relSe = [](double se, double seR, double Er) { return std::fabs(Er) > 1e-300 ? std::hypot(se, seR)/std::fabs(Er) : std::nan(""); };
    for (const double v : {r.E, r.Ese, r.Eref, r.ErefSe, rel(r.E, r.Eref), relSe(r.Ese, r.ErefSe, r.Eref), r.EL, r.ELse, r.ELref,
                           r.ELrefSe, rel(r.EL, r.ELref), relSe(r.ELse, r.ELrefSe, r.ELref)})
        num(r.hasE ? v : std::nan(""));
    for (int k = 0; k < kRingNV; k++) num(r.vert[k], "%.5f");
    for (int k = 0; k < kRingNV; k++) num(r.vref[k], "%.5f");
    num(r.vErr, "%.5f");
    num(r.hasC ? r.cse.maxCoeff() : std::nan(""), "%.5f");
    std::string s;
    for (size_t i = 0; i < c.size(); i++) s += (i ? "," : "") + c[i];
    return s + "\n";
}

void printRowHeader(Report& out)
{
    out("    %-30s %-5s %7s %7s %8s %8s %8s | %6s %6s %6s %6s | %7s %7s %7s | %7s %7s | %6s %6s\n", "method", "g", "walks/M",
        "moves/M", "walks s", "post s", "total s", "H RMS%", "nc %", "max %", "fld se", "grad%", "grad nc", "|grad|%",
        "E err%", "+-", "vert", "c se");
}
void printRow(Report& out, const BenchRow& r)
{
    const ErrStats& e = r.e[0];
    const double Erel = r.hasE && std::fabs(r.Eref) > 1e-300 ? 100*(r.E - r.Eref)/r.Eref : std::nan("");
    const double Ese = r.hasE && std::fabs(r.Eref) > 1e-300 ? 100*std::hypot(r.Ese, r.ErefSe)/std::fabs(r.Eref) : std::nan("");
    const double post = r.tMed - r.tWalks - r.tQ - r.tR;
    out("    %-30s %-5s %7.2f %7.1f %8.3f %8.3f %8.3f | %6.3f %6.3f %6.2f %6.3f | %7.2f %7.2f %7.3f | %+7.2f %7.2f | %6.4f %6.4f\n",
        r.label.c_str(), r.g ? "ramp" : "g=1", 1e-6*double(r.walks), 1e-6*double(r.moves), r.tWalks + r.tQ + r.tR, post, r.tMed,
        100*e.rmsH(), 100*e.corrH(), 100*e.maxH/kTwoGMax, 100*e.seFldH(), 100*gradRel(e, r.gt2[0]), 100*gradRelNc(e, r.gt2[0]),
        100*e.rmsG(), Erel, Ese, r.vErr, r.hasC ? r.cse.maxCoeff() : 0.0);
}

// "a/b/c,d/e" -> lists of ints; "none" -> empty; "default" -> def
std::vector<std::vector<int>> parseGrid(const std::string& s, const std::string& def)
{
    const std::string src = s == "default" ? def : s;
    std::vector<std::vector<int>> out;
    if (src == "none" || src.empty()) return out;
    std::stringstream ss(src);
    std::string item;
    while (std::getline(ss, item, ',')) {
        std::vector<int> v;
        std::stringstream si(item);
        std::string tok;
        while (std::getline(si, tok, '/'))
            if (!tok.empty()) v.push_back(std::atoi(tok.c_str()));
        if (!v.empty()) out.push_back(v);
    }
    return out;
}
int defaultBatches(int n) { return n <= 32 ? 8 : n >= 1280 ? 64 : 16; }

// runs fn (which returns its total wall time) once, then again until repeats runs if the first took under 30 s
template <class Fn>
std::vector<double> timedRuns(int repeats, Fn&& fn)
{
    std::vector<double> t;
    t.push_back(fn());
    if (t[0] < 30.0)
        for (int r = 1; r < repeats; r++) t.push_back(fn());
    return t;
}
// the repeat whose total time is the median (the closest to it for an even count): its parts are the ones reported
size_t medianRun(const std::vector<double>& t)
{
    const double m = medianOf(t);
    size_t k = 0;
    for (size_t i = 1; i < t.size(); i++)
        if (std::fabs(t[i] - m) < std::fabs(t[k] - m)) k = i;
    return k;
}
struct RunParts {
    double setup = 0, q = 0, r = 0, walks = 0, gram = 0, glue = 0, field = 0;
};

// ---- the equal-time and equal-accuracy summary of one configuration and g type ----
void benchSummary(Report& out, const std::vector<BenchRow>& rows, const std::string& config, int g, std::string& csv)
{
    std::vector<std::string> fam;
    for (const BenchRow& r : rows)
        if (r.config == config && r.g == g && r.seedRep == 0 && std::find(fam.begin(), fam.end(), r.family) == fam.end())
            fam.push_back(r.family);
    if (fam.empty()) return;
    struct Pt {
        double t, rms;
        std::string label;
    };
    auto pts = [&](const std::string& f) {
        std::vector<Pt> p;
        for (const BenchRow& r : rows)
            if (r.config == config && r.g == g && r.seedRep == 0 && r.family == f) p.push_back({r.tMed, r.e[0].rmsH(), r.label});
        std::sort(p.begin(), p.end(), [](const Pt& a, const Pt& b) { return a.t < b.t; });
        return p;
    };
    const char* gn = g ? "ramp" : "g1";
    out("\n  summary %s, %s (H RMS / 2g, raw; times = medians, 16 threads)\n", config.c_str(), g ? "ramp" : "g = 1");
    out("    equal time: the best H RMS each method reaches within the budget\n      %-14s", "budget (s)");
    const double budgets[] = {0.25, 0.5, 1, 2, 5, 10, 30, 60, 120};
    for (double T : budgets) out(" %9g", T);
    out("\n");
    for (const std::string& f : fam) {
        out("      %-14s", f.c_str());
        const std::vector<Pt> p = pts(f);
        for (double T : budgets) {
            double best = kInf;
            std::string lab;
            for (const Pt& x : p)
                if (x.t <= T && x.rms < best) {
                    best = x.rms;
                    lab = x.label;
                }
            out(" %9s", std::isfinite(best) ? (fmtd("%.3f", 100*best) + "%").c_str() : "-");
            csv += config + "," + gn + "," + f + ",equal_time," + csvNum(T) + "," + (std::isfinite(best) ? csvNum(best) : "") + ",\"" +
                   lab + "\",\n";
        }
        out("\n");
    }
    out("    equal accuracy: the time each method needs (first measured setting at or below; [log-log interpolated])\n");
    const double thr[] = {0.03, 0.02, 0.015, 0.01};
    std::vector<std::vector<double>> tMeas(fam.size()), tInt(fam.size());
    for (size_t fi = 0; fi < fam.size(); fi++) {
        const std::vector<Pt> p = pts(fam[fi]);
        out("      %-14s", fam[fi].c_str());
        for (double th : thr) {
            double tm = kInf, ti = kInf;
            std::string lab;
            double bestSoFar = kInf;
            for (size_t i = 0; i < p.size(); i++) {
                if (p[i].rms <= th) {
                    tm = p[i].t;
                    lab = p[i].label;
                    ti = tm;
                    if (i > 0 && bestSoFar > th && p[i].rms < th && p[i - 1].rms > th && p[i].t > p[i - 1].t)
                        ti = std::exp(std::log(p[i - 1].t) + (std::log(th) - std::log(p[i - 1].rms))*
                                                                 (std::log(p[i].t) - std::log(p[i - 1].t))/(std::log(p[i].rms) - std::log(p[i - 1].rms)));
                    break;
                }
                bestSoFar = std::min(bestSoFar, p[i].rms);
            }
            tMeas[fi].push_back(tm);
            tInt[fi].push_back(ti);
            out("  %4.1f%%: %s", 100*th, std::isfinite(tm) ? (fmtd("%.3f s", tm) + " [" + fmtd("%.3f", ti) + "]").c_str() : "not reached");
            csv += config + "," + gn + "," + fam[fi] + ",equal_accuracy," + csvNum(th) + "," + (std::isfinite(tm) ? csvNum(tm) : "") +
                   ",\"" + lab + "\"," + (std::isfinite(ti) ? csvNum(ti) : "") + "\n";
        }
        out("\n");
    }
    // DERIVED (not measured): per node in a real mesh, a tet's local solves are shared by the 4 rings of its vertices, a
    // direct ring solve is not: ours per node = (n_r/4) (T_ours/8), direct per node = (n_r/8) T_direct (n_r tets per ring)
    const auto io = std::find(fam.begin(), fam.end(), std::string("ours")), id = std::find(fam.begin(), fam.end(), std::string("direct"));
    if (io != fam.end() && id != fam.end()) {
        const size_t a = size_t(io - fam.begin()), b = size_t(id - fam.begin());
        for (size_t k = 0; k < std::size(thr); k++) {
            if (!std::isfinite(tMeas[a][k]) || !std::isfinite(tMeas[b][k])) continue;
            const double ro = tMeas[a][k], rd = tMeas[b][k];
            out("    DERIVED per node (real mesh, n_r = 20 / 24 tets per ring, each tet in 4 rings), at %.1f%%: ours %.3f / %.3f s, direct "
                "%.3f / %.3f s; direct/ours = 4 x %.2f = %.1fx (estimate from the measured ring times)\n", 100*thr[k], 20*ro/32,
                24*ro/32, 20*rd/8, 24*rd/8, rd/ro, 4*rd/ro);
            // the first measured setting at or below the threshold overstates a method's time by up to its grid step; the
            // log-log interpolated times (each method's own) bound the ratio from the other side
            const double ri = tInt[a][k], rdi = tInt[b][k];
            if (std::isfinite(ri) && std::isfinite(rdi) && ri > 0)
                out("      same with the log-log interpolated times: ours %.3f s, direct %.3f s per ring; direct/ours = 4 x %.2f = %.1fx\n", ri,
                    rdi, rdi/ri, 4*rdi/ri);
        }
    }
}

// ---- self-test: the lean walks bit for bit, J_xi, the skip decisions, same-target z tests, determinism ----
bool sameSums(const MultiSums& a, const MultiSums& b)
{
    return a.K == b.K && a.nBatches == b.nBatches && a.walksPerBatch == b.walksPerBatch && a.nPoints == b.nPoints &&
           a.batch.size() == b.batch.size() && a.uu.size() == b.uu.size() &&
           std::memcmp(a.batch.data(), b.batch.data(), a.batch.size()*sizeof(double)) == 0 &&
           std::memcmp(a.uu.data(), b.uu.data(), a.uu.size()*sizeof(double)) == 0;
}

bool benchSelftest(Report& out)
{
    bool ok = true;
    out("\n=== self-test 1: solveMultiLean / solveMultiLeanJobs against solveMulti, bit for bit (memcmp of every batch sum and\n"
        "    the per-point second moments) ===\n");
    {
        PerfSetup s;
        perfSetup(s, 2, 0, 1024);
        const Tet ref = referenceTet();
        const int N = 64, B = 16;
        int nRuns = 0, nSame = 0;
        long long walks = 0;
        std::vector<std::vector<FieldSpec>> fieldsR(kRingNT);
        const std::vector<FieldSpec> fieldsC = ringLevelCFields(s.J);
        std::vector<QueryPoints> qs(kRingNT);
        std::vector<TerminalValues> tvs(kRingNT);
        std::vector<HarmonicExtension> exts(kRingNT);
        std::vector<Cut> fronts(kRingNT);
        std::vector<MultiJob> jobs;
        std::vector<MultiSums> plain;
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = s.R.T[size_t(t)];
            const RingTetCrack& r = s.rt[size_t(t)];
            fieldsR[size_t(t)] = ringRecipeFields(T);
            fronts[size_t(t)] = r.cut() ? r.gFront : r.xc.strips.frontCut(r.xc.strips.nStrips());
            const std::vector<char> grad(s.ps[size_t(t)].size(), 1);
            qs[size_t(t)].add(s.ps[size_t(t)].xi, &grad, mixSeed(s.gp.seed, ownTag(s, t)));
            for (int which = 0; which < 2; which++) {   // 0: the recipe run (odd mirror), 1: the Level-C run (plain first ball)
                const std::vector<FieldSpec>* fs = which ? &fieldsC : &fieldsR[size_t(t)];
                const Cut* front = &fronts[size_t(t)];
                const TetFrame* Tp = &T;
                const int K = int(fs->size());
                const TerminalValues tv = [Tp, fs, front, K](const WalkTerminal& w, double* u) {   // runOwn's scoring
                    if (w.face >= 0) {
                        const Vec3 X = Tp->toX(w.x);
                        for (int k = 0; k < K; k++) u[k] = (*fs)[size_t(k)].face ? (*fs)[size_t(k)].face(X) : 0.0;
                    } else {
                        for (int k = 0; k < K; k++) u[k] = (*fs)[size_t(k)].crack ? double(w.side)*front->g(w.x, (*fs)[size_t(k)].g) : 0.0;
                    }
                };
                const HarmonicExtension ext = [Tp, fs, K](const Vec3& xi, double* phi) {
                    const Vec3 X = Tp->toX(xi);
                    for (int k = 0; k < K; k++) phi[k] = (*fs)[size_t(k)].face ? (*fs)[size_t(k)].face(X) : 0.0;
                };
                const WalkParams wp = gateWalkParams(allDirichlet(), true, which == 0);
                SolveStats sa, sb;
                const MultiSums a = solveMulti(ref, r.xc.strips, qs[size_t(t)], N, B, wp, K, tv, &sa, &ext);
                const MultiSums b = solveMultiLean(ref, r.xc.strips, qs[size_t(t)], N, B, wp, K, tv, &sb, &ext);
                nRuns++;
                walks += sa.walksRequested;
                const bool same = sameSums(a, b) && sa.totalSteps == sb.totalSteps && sa.walksCounted == sb.walksCounted;
                nSame += same;
                if (!same) out("  MISMATCH: %s %s run\n", T.name, which ? "Level-C" : "recipe");
                if (which == 0) {   // the fused version of the recipe runs (all 8 tets in one loop), compared below
                    tvs[size_t(t)] = tv;
                    exts[size_t(t)] = ext;
                    plain.push_back(a);
                }
            }
        }
        for (int t = 0; t < kRingNT; t++) {
            MultiJob jb;
            jb.tet = &ref;
            jb.cut = &s.rt[size_t(t)].xc.strips;
            jb.pts = &qs[size_t(t)];
            jb.nWalks = N;
            jb.nBatches = B;
            jb.K = 6;
            jb.wp = gateWalkParams(allDirichlet(), true, true);
            jb.values = tvs[size_t(t)];
            jb.ext = &exts[size_t(t)];
            jobs.push_back(jb);
        }
        std::vector<SolveStats> js;
        const std::vector<MultiSums> fused = solveMultiLeanJobs(jobs, &js);
        int nFused = 0;
        for (int t = 0; t < kRingNT; t++) nFused += sameSums(fused[size_t(t)], plain[size_t(t)]);
        // the direct on O
        std::vector<Vec3> XO;
        for (const PointSet& p : s.ps) XO.insert(XO.end(), p.X.begin(), p.X.end());
        const StripCut scO = ringStrips(s.rc);
        const Cut frO = scO.frontCut(scO.nStrips());
        const TerminalValues tvO = [&frO](const WalkTerminal& t, double* u) {
            u[0] = double(t.side)*frO.g(t.x, GType::One);
            u[1] = double(t.side)*frO.g(t.x, GType::SqrtFront);
        };
        QueryPoints qO;
        const std::vector<char> gO(XO.size(), 1);
        qO.add(XO, &gO, mixSeed(s.gp.seed, 17999));
        SolveStats so1, so2;
        const MultiSums o1 = solveMulti(s.O, scO, qO, 16, 8, ringReferenceWalkParams(s.O), 2, tvO, &so1);
        const MultiSums o2 = solveMultiLean(s.O, scO, qO, 16, 8, ringReferenceWalkParams(s.O), 2, tvO, &so2);
        const bool oSame = sameSums(o1, o2) && so1.totalSteps == so2.totalSteps;
        out("  ring 2 K1, about 1024 points per tet, N = %d in %d batches: %d of %d per-tet runs identical (8 recipe runs with the\n"
            "  odd mirror ball, 8 Level-C runs with the plain one; %lld walks each way); fused (8 recipe runs in one loop): %d of 8\n"
            "  identical; direct on O (Neumann faces, N = 16, %zu points, %lld walks): %s\n", N, B, nSame, nRuns, walks, nFused,
            XO.size(), so1.walksRequested, oSame ? "identical" : "MISMATCH");
        ok = ok && nSame == nRuns && nFused == kRingNT && oSame;
    }

    out("\n=== self-test 2: J_xi (finite-difference gradients; ring 1: J_xi = J) ===\n");
    for (int ring = 1; ring <= 2; ring++)
        for (int k = 0; k < 3; k++) perfCheckXiFront(out, ring, k);

    out("\n=== self-test 3: which fields walk (ring 1: R = 0 in every tet for g = 1; Q only in cut tets) ===\n");
    for (int ring = 1; ring <= 2; ring++)
        for (int k = 0; k < 3; k++) {
            PerfSetup s;
            perfSetup(s, ring, k, 256);
            const bool rampOn = ring == 1 && k < 2;
            std::string line;
            int nR = 0, nR2 = 0, nQ = 0, nCut = 0;
            for (int t = 0; t < kRingNT; t++) {
                CvRun c, c2;
                fastDecide(c, s.R.T[size_t(t)], s.rt[size_t(t)], s.rc, s.J, false);
                nR += c.walkR;
                nQ += c.cut;
                nCut += s.rt[size_t(t)].cut();
                line += std::string(t ? " " : "") + s.R.T[size_t(t)].name + (c.cut ? ":Q" : ":-") + (c.walkR ? "R" : "-");
                if (rampOn) {
                    fastDecide(c2, s.R.T[size_t(t)], s.rt[size_t(t)], s.rc, s.J, true);
                    nR2 += c2.walkR;
                    line += c2.walkR ? "(ramp R)" : "";
                }
            }
            const bool tOk = nQ == nCut && (ring != 1 || nR == 0);
            ok = ok && tOk;
            out("  ring %d %s: Q in %d tets (cut: %d), R in %d tets%s: %s -> %s\n", ring, s.rc.name.c_str(), nQ, nCut, nR,
                rampOn ? (", ramp R in " + std::to_string(nR2)).c_str() : "", line.c_str(), tOk ? "OK" : "FAILED");
        }

    out("\n=== self-test 4: the fast path's fields against check 4's runs (same targets, independent walks, N = 256 each; points\n"
        "    within %.2f of a face left out as in check 4; pooled z per field should be N(0, 1)) ===\n", kCompareBand);
    for (const std::pair<int, int>& rk : {std::pair<int, int>{2, 0}, {2, 1}, {1, 1}}) {
        const int ring = rk.first, k = rk.second;
        PerfSetup s;
        perfSetup(s, ring, k, 1024);
        const bool rampOn = ring == 1;
        RingFastParams fp;
        fp.nQ = fp.nR = 256;
        fp.nBatches = 16;
        fp.ramp = rampOn;
        const RingFast f = runRingFast(s.R, s.rc, s.ps, s.J, fp, s.gp.seed, 21000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki) + 500000,
                                       false);
        const RingCurrent cur = runRingCurrent(s.R, s.rc, s.ps, s.J, 256, 16, s.gp.seed,
                                               12000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki), rampOn ? 2 : 1, false);
        ZStats zP[4], zH[2];
        for (int t = 0; t < kRingNT; t++) {
            const CvRun& c = f.cv[size_t(t)];
            const PointSet& p = s.ps[size_t(t)];
            for (size_t i = 0; i < p.size(); i++) {
                if (dirichletDist(allDirichlet(), p.xi[i]) < kCompareBand) continue;
                double h, v, n2;
                long long n;
                for (int q = 0; q < 4; q++) {
                    comboH(cur.rec[size_t(t)].sums, i, unitCoef(6, kRingP0 + q), h, v, n);
                    double hq = 0, vq = 0;
                    if (c.hasQ) comboH(c.runQ.sums, i, unitCoef(4, q), hq, vq, n);
                    const double lam = q == 0 ? 1.0 - p.xi[i].sum() : p.xi[i][q - 1];
                    zP[q].add(i, lam + hq - h, v + vq);
                }
                for (int m = 0; m < (rampOn ? 2 : 1); m++) {
                    comboH(cur.hc[size_t(t)].sums, i, unitCoef(2, m), h, v, n);
                    double hr = 0, vr = 0;
                    const int iR = m ? c.iR2 : c.iR1;
                    if (c.hasR && iR >= 0) comboH(c.runR.sums, i, unitCoef(c.runR.sums.K, iR), hr, vr, n);
                    n2 = (m ? c.xf.J2(p.xi[i]) : c.xf.J1(p.xi[i])) + hr;
                    zH[m].add(i, n2 - h, v + vr);
                }
            }
        }
        double zmax = 0;
        std::string zs;
        for (int q = 0; q < 4; q++) {
            zs += fmtd(" %+.2f", zP[q].pooled());
            zmax = std::max(zmax, std::fabs(zP[q].pooled()));
        }
        out("  ring %d %s: P_k = lambda_k + Q_k vs P_k: pooled z per corner%s (|z| < 3 at %.2f%% of %lld); H0 = J_xi + R vs H0 (Level C):\n"
            "    g = 1 pooled z %+.2f (|z| < 3 at %.2f%% of %lld)%s\n", ring, s.rc.name.c_str(), zs.c_str(), 100*zP[0].frac(), zP[0].n,
            zH[0].pooled(), 100*zH[0].frac(), zH[0].n,
            rampOn ? (fmtd("; ramp pooled z %+.2f", zH[1].pooled()) + fmtd(" (|z| < 3 at %.2f%%)", 100*zH[1].frac())).c_str() : "");
        zmax = std::max(zmax, std::fabs(zH[0].pooled()));
        if (rampOn) zmax = std::max(zmax, std::fabs(zH[1].pooled()));
        if (zmax > 4) out("    NOTE: a pooled |z| above 4 (largest %.2f)\n", zmax);
    }

    out("\n=== self-test 5: the fast path is deterministic and its variants agree (ring 2 K1, check 4's points, Q 32 / R 16) ===\n");
    {
        PerfSetup s;
        perfSetup(s, 2, 0, gateParams(false).nRingPoints);
        const uint64_t tag = 21000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki);
        RingFastParams fp;
        const RingFast a = runRingFast(s.R, s.rc, s.ps, s.J, fp, s.gp.seed, tag, true);
        const RingFast b = runRingFast(s.R, s.rc, s.ps, s.J, fp, s.gp.seed, tag, true);
        fp.fuse = false;
        const RingFast c = runRingFast(s.R, s.rc, s.ps, s.J, fp, s.gp.seed, tag, true);
        fp.lean = false;
        const RingFast d = runRingFast(s.R, s.rc, s.ps, s.J, fp, s.gp.seed, tag, true);
        auto same = [](const RingFast& x, const RingFast& y) {
            bool e = std::memcmp(x.gl[0].c.data(), y.gl[0].c.data(), sizeof(double)*kRingNV) == 0;
            for (int t = 0; t < kRingNT; t++) {
                e = e && x.cv[size_t(t)].hasQ == y.cv[size_t(t)].hasQ && x.cv[size_t(t)].hasR == y.cv[size_t(t)].hasR;
                if (x.cv[size_t(t)].hasQ) e = e && sameSums(x.cv[size_t(t)].runQ.sums, y.cv[size_t(t)].runQ.sums);
                if (x.cv[size_t(t)].hasR) e = e && sameSums(x.cv[size_t(t)].runR.sums, y.cv[size_t(t)].runR.sums);
                e = e && x.F[0][size_t(t)].H == y.F[0][size_t(t)].H;
            }
            return e;
        };
        const bool ab = same(a, b), ac = same(a, c), ad = same(a, d);
        out("  repeat: %s; unfused: %s; runWalk (solveMulti): %s; times %.3f / %.3f / %.3f / %.3f s; glued c = %s\n",
            ab ? "identical" : "DIFFERENT", ac ? "identical" : "DIFFERENT", ad ? "identical" : "DIFFERENT", a.secTotal, b.secTotal,
            c.secTotal, d.secTotal, fmtVec(a.gl[0].c).c_str());
        ok = ok && ab && ac && ad;
    }
    out("\nself-test: %s\n", ok ? "PASS" : "FAILED");
    return ok;
}

// fnv-1a of a double vector (bit-identity checks across builds)
uint64_t hashVec(const std::vector<double>& v) { return fnv1a(v.data(), v.size()*sizeof(double)); }

// ---- the fast ring view (RingViewParams::fast) ----
RingView gateRingViewFast(const RingViewParams& p, std::atomic<int>* progress)
{
    if (p.ring != 1 && p.ring != 2) throw std::invalid_argument("ring view: ring 1 or 2");
    if (p.crack < 0 || p.crack > 2) throw std::invalid_argument("ring view: crack 0, 1 or 2 (K1, K2, K3)");
    const auto t0 = std::chrono::steady_clock::now();
    const GateParams gp = gateParams(false);
    const int nPts = std::max(64, p.pointsPerTet);
    const int B = defaultBatches(p.walksTet);
    auto whole = [](int n, int b) { return std::max(2*b, (n + b/2)/b*b); };
    const int nQ = whole(p.walksTet, B), nR = whole(p.walksR, B);
    const int Bref = p.walksRef >= 1280 ? 64 : 16, nRef = std::max(4*Bref, (p.walksRef + Bref/2)/Bref*Bref);
    const Ring R = makeRing(p.ring);
    const RingCrack rc = ringCrack(p.crack);
    const int ri = R.id - 1, ki = p.crack;
    const int nModes = R.id == 1 && !rc.through ? 2 : 1;

    RingView v;
    v.title = "check 4, " + R.name + ", crack " + rc.name + ": " + rc.desc + " [fast method]";
    std::copy(R.V.begin(), R.V.end(), v.verts.begin());
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        v.tets.push_back({size_t(T.global[0]), size_t(T.global[1]), size_t(T.global[2]), size_t(T.global[3])});
        v.tetName.push_back(T.name);
        v.tetCut.push_back(char(ringTetCrack(T, rc).cut()));
    }
    v.crack = rc.cutO.polygon;
    v.n = rc.n;
    v.d = rc.d;
    v.pointsPerTet = nPts;
    v.walksTet = nQ;
    v.walksR = nR;
    v.walksRef = nRef;
    v.fast = true;
    if (p.geometryOnly) return v;

    PerfSetup s;
    perfSetup(s, p.ring, p.crack, nPts);
    RefRun ref;
    bool cached = false;
    if (nPts == gp.nRingPoints && nRef == gp.nRef) cached = loadRingGT(s, ref);   // read-only
    if (!cached) {   // solved here at the asked walk count, never cached
        GateParams g2 = gp;
        g2.nRingPoints = nPts;
        g2.nRef = nRef;
        g2.nRefBatches = Bref;
        g2.quick = true;
        Report out(false);
        ref = runRingReference(out, g2, s.O, s.R, s.rc, "check4_ring" + std::to_string(R.id) + "_" + rc.name,
                               14000 + 10*uint64_t(ri) + uint64_t(ki), s.ps, ringProbes(s.R), "gate_cache");
        if (!ref.ok) throw std::runtime_error("the reference on the ring was not run (see the console)");
    }
    RingFastParams fp;
    fp.nQ = nQ;
    fp.nR = nR;
    fp.nBatches = B;
    fp.ramp = nModes == 2;
    const RingFast f = runRingFast(s.R, s.rc, s.ps, s.J, fp, gp.seed, 21000 + 1000*uint64_t(ri) + 100*uint64_t(ki), true, progress);

    std::vector<double> w;
    std::vector<size_t> off(kRingNT, 0);
    for (int t = 0; t < kRingNT; t++) {
        off[size_t(t)] = v.X.size();
        v.X.insert(v.X.end(), s.ps[size_t(t)].X.begin(), s.ps[size_t(t)].X.end());
        w.insert(w.end(), s.ps[size_t(t)].w.begin(), s.ps[size_t(t)].w.end());
        v.tetOf.insert(v.tetOf.end(), s.ps[size_t(t)].size(), t);
    }
    for (int k = 0; k < kRingNV; k++) {
        v.basis[size_t(k)].assign(v.X.size(), 0.0);
        v.inSupport[size_t(k)].assign(v.X.size(), 0);
    }
    GateView::Field h0{"H0 (g = 1): local solves with crack +-1 and faces = J"};
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = s.R.T[size_t(t)];
        const CvRun& c = f.cv[size_t(t)];
        const PointSet& ps = s.ps[size_t(t)];
        for (size_t i = 0; i < ps.size(); i++) {
            const Vec3& xi = ps.xi[i];
            const double lam[4] = {1.0 - xi.sum(), xi[0], xi[1], xi[2]};
            double H, var;
            long long n;
            for (int k = 0; k < 4; k++) {   // P_k = lambda_k + Q_k
                double q = 0;
                if (c.hasQ) comboH(c.runQ.sums, i, unitCoef(4, k), q, var, n);
                v.basis[size_t(T.global[k])][off[size_t(t)] + i] = lam[k] + q;
                v.inSupport[size_t(T.global[k])][off[size_t(t)] + i] = 1;
            }
            H = c.xf.J1(xi);   // H0 = J_xi + R
            if (c.hasR) {
                double r = 0;
                comboH(c.runR.sums, i, unitCoef(c.runR.sums.K, 0), r, var, n);
                H += r;
            }
            h0.vol.push_back(H);
        }
    }
    char buf[512];
    std::snprintf(buf, sizeof buf, "%zu sample points (%d per tet asked); FAST method: walks per point Q %d (in the %d cut tets), R %d "
                  "(in %d tets), %d batches; %d for the GT", v.X.size(), nPts, nQ, f.nQtets, nR, f.nRtets, B, nRef);
    v.report.push_back(buf);
    const double len = double(f.stepsQ + f.stepsR)/double(std::max(1LL, f.walksQ + f.walksR - f.escaped));
    std::snprintf(buf, sizeof buf, "local solves (control variates P = lambda + Q, H0 = J_xi + R; lean walks): %.3f s (Q walks %.3f s, R "
                  "walks %.3f s, Grams + glue + fields %.3f s), %.2f M walks, mean walk length %.1f; GT (whole-ring solve): %s %.1f s",
                  f.secTotal, f.secQ, f.secR, f.secGram + f.secGlue + f.secField, 1e-6*double(f.walksQ + f.walksR), len,
                  cached ? "loaded from gate_cache (original run)" : "solved,", ref.stats.seconds);
    v.report.push_back(buf);
    v.report.push_back("error of the glued H against the GT, as a fraction of the jump 2: RMS (average) / max; vertex values: glued / GT");
    const char* gName[2] = {"g = 1", "ramp g"};
    for (int m = 0; m < nModes; m++) {
        const std::string sfx = std::string(" (") + gName[m] + ")";
        GateView::Field rf{"GT" + sfx};
        size_t o = 0;
        for (int t = 0; t < kRingNT; t++) {
            const FieldAt fr = refField(ref, m, o, s.ps[size_t(t)].size());
            rf.vol.insert(rf.vol.end(), fr.H.begin(), fr.H.end());
            o += s.ps[size_t(t)].size();
        }
        const FieldAt rp = refField(ref, m, ref.nA, ref.nF);
        v.fields.push_back(rf);
        const Vec7& c = f.gl[m].c;
        GateView::Field dec{"glued H" + sfx}, dif{"error = glued H - GT" + sfx};
        for (int t = 0; t < kRingNT; t++) dec.vol.insert(dec.vol.end(), f.F[m][size_t(t)].H.begin(), f.F[m][size_t(t)].H.end());
        dif.difference = true;
        double sq = 0, sw = 0;
        for (size_t i = 0; i < dec.vol.size(); i++) {
            dif.vol.push_back(dec.vol[i] - rf.vol[i]);
            sq += w[i]*dif.vol[i]*dif.vol[i];
            sw += w[i];
            dif.max = std::max(dif.max, std::fabs(dif.vol[i]));
        }
        dif.rms = std::sqrt(sq/sw)/kTwoGMax;
        dif.max /= kTwoGMax;
        std::string vt;
        for (int k = 0; k < kRingNV; k++) {
            std::snprintf(buf, sizeof buf, "%s v%d %+.3f/%+.3f", k ? "," : "", k, c[k] + s.J[m](R.V[size_t(k)]), rp.H[size_t(k)]);
            vt += buf;
        }
        std::snprintf(buf, sizeof buf, "glued H%s: error RMS %.2f%%, max %.1f%%; vertices%s", sfx.c_str(), 100*dif.rms, 100*dif.max,
                      vt.c_str());
        v.report.push_back(buf);
        v.fields.push_back(std::move(dec));
        v.fields.push_back(std::move(dif));
    }
    v.fields.push_back(std::move(h0));
    std::snprintf(buf, sizeof buf, "total %.1f s", secondsSince(t0));
    v.report.push_back(buf);
    return v;
}

} // namespace

bool runRingBench(const RingBenchParams& bp)
{
    using clk = std::chrono::steady_clock;
    const auto tStart = clk::now();
    Report out(false);
    std::time_t now = std::time(nullptr);
    std::tm tmNow{};
#ifdef _WIN32
    localtime_s(&tmNow, &now);
#else
    localtime_r(&now, &tmNow);
#endif
    char date[64], stamp[64];
    std::strftime(date, sizeof date, "%Y-%m-%d %H:%M:%S", &tmNow);
    std::strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", &tmNow);
    out("--ring-bench%s%s: %s; %u hardware threads (TBB default arena)\n", bp.selftest ? " --selftest" : "", bp.view ? " --view" : "",
        date, std::thread::hardware_concurrency());

    if (bp.view) {   // the ring view headless: its panel's lines, and hashes of its fields for bit-identity checks
        const RingView v = gateRingView(bp.viewParams);
        out("%s\n", v.title.c_str());
        for (const std::string& s : v.report) out("  %s\n", s.c_str());
        for (const GateView::Field& f : v.fields) out("  field \"%s\": %zu values, fnv1a %016llx\n", f.name.c_str(), f.vol.size(),
                                                      (unsigned long long)hashVec(f.vol));
        for (int k = 0; k < kRingNV; k++) out("  basis v%d: fnv1a %016llx\n", k, (unsigned long long)hashVec(v.basis[size_t(k)]));
        return true;
    }
    if (bp.selftest) return benchSelftest(out);

    // the grids
    const std::vector<std::vector<int>> gOurs = parseGrid(bp.ours, bp.quick ? "32/16/8,64/32/16" :
                                                          "16/16/8,32/16/8,32/64/8,64/32/16,64/256/16,128/64/16,256/128/16,1024/256/16,2560/640/64");
    const std::vector<std::vector<int>> gDirect = parseGrid(bp.direct, bp.quick ? "64,256" : "16,32,64,128,256,512,1024,2560");
    const std::vector<std::vector<int>> gCur = parseGrid(bp.current, bp.quick ? "none" : "256");
    const int repeats = std::max(1, bp.repeats);
    std::vector<std::pair<int, int>> cfgs;
    for (int r = 1; r <= 2; r++)
        for (int k = 0; k < 3; k++)
            if ((bp.ring == 0 ? !bp.quick || r == 2 : bp.ring == r) && (bp.crack < 0 ? !bp.quick || k == 0 : bp.crack == k))
                cfgs.push_back({r, k});
    const bool loose = bp.shellRel != kGateShellRel;
    std::string tag = bp.tag.empty() ? std::string("bench_") + stamp : bp.tag;
    std::error_code ec;
    std::filesystem::create_directories("ring_perf", ec);
    auto freeName = [](const std::string& base) {   // never overwrite: _2, _3, ... if the name is taken
        std::string b = base;
        for (int i = 2; std::filesystem::exists("ring_perf/" + b + "_rows.csv"); i++) b = base + "_" + std::to_string(i);
        return b;
    };
    tag = freeName(tag);
    out("output: ring_perf/%s_rows.csv, ring_perf/%s_summary.csv, ring_perf/%s_log.txt\n", tag.c_str(), tag.c_str(), tag.c_str());
    out("settings: walk stop distance %.3g sqrt(2) (both methods; the reference's: 1e-6)%s; step function: %s; repeats %d (median; runs "
        "over 30 s once); seed sets %d\n", bp.shellRel, loose ? " [LOOSE: biased, opt-in]" : "",
        bp.lean ? "lean (solveMultiLean, the same walks without box and face tests)" : "runWalk (solveMulti)", repeats, bp.seeds);
    out("  ours (nQ/nR/B):");
    for (const auto& g : gOurs) out(" %d/%d/%d", g[0], g.size() > 1 ? g[1] : g[0], g.size() > 2 ? g[2] : defaultBatches(g[0]));
    out("; direct (N/B):");
    for (const auto& g : gDirect) out(" %d/%d", g[0], g.size() > 1 ? g[1] : defaultBatches(g[0]));
    out("; current (N/B):");
    for (const auto& g : gCur) out(" %d/%d", g[0], g.size() > 1 ? g[1] : defaultBatches(g[0]));
    out("\n");

    std::vector<BenchRow> rows;
    std::string summaryCsv = "config,g,family,kind,budget_or_threshold,value,label,interpolated_time\n";
    const std::string sfxShell = loose ? fmtd(" shell %.0e", bp.shellRel) : "";
    const std::string famSfx = std::string(loose ? "-loose" : "") + (bp.lean ? "" : "-runWalk");
    for (const std::pair<int, int>& cfg : cfgs) {
        const int ring = cfg.first, crack = cfg.second;
        out("\n=== ring %d, K%d ===\n", ring, crack + 1);
        BenchCfg b;
        if (!benchSetup(out, b, ring, crack)) continue;
        const PerfSetup& s = b.s;
        out("  point sets and weights (ringTetCrack + ringPoints, 8 tets) built in %.3f s: shared by both methods, timed for neither\n",
            s.secGeometry);
        out("  %zu points in the 8 tets + 7 vertex probes; tets cut by the crack:", b.nPts);
        for (const RingTetCrack& r : s.rt) out(" %d", int(r.cut()));
        out("\n");
        printRowHeader(out);
        const size_t rowStart = rows.size();
        auto baseRow = [&](const std::string& method, int g, int rep) {
            BenchRow r;
            r.config = b.name;
            r.ring = ring;
            r.crack = crack;
            r.g = g;
            r.method = method;
            r.seedRep = rep;
            r.points = (long long)b.nPts;
            return r;
        };
        for (int rep = 0; rep < std::max(1, bp.seeds); rep++) {
            // A: ours
            for (const auto& gs : gOurs) {
                RingFastParams fp;
                fp.nQ = gs[0];
                fp.nR = gs.size() > 1 ? gs[1] : gs[0];
                fp.nBatches = gs.size() > 2 ? gs[2] : defaultBatches(fp.nQ);
                fp.shellRel = bp.shellRel;
                fp.lean = bp.lean;
                fp.fuse = bp.fuse;
                for (int m = 0; m < b.nModes; m++) {
                    fp.ramp = m == 1;
                    const uint64_t tagBase = 21000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki) + 100000*uint64_t(rep);
                    RingFast f;
                    bool determ = true;
                    Vec7 c0 = Vec7::Zero();
                    int run = 0;
                    std::vector<RunParts> parts;
                    const std::vector<double> times = timedRuns(repeats, [&] {
                        f = runRingFast(s.R, s.rc, s.ps, s.J, fp, s.gp.seed, tagBase, true);
                        if (run++ == 0) c0 = f.gl[m].c;
                        else determ = determ && std::memcmp(c0.data(), f.gl[m].c.data(), sizeof(double)*kRingNV) == 0;
                        parts.push_back({f.secSetup, f.secQ, f.secR, 0.0, f.secGram, f.secGlue, f.secField});
                        return f.secTotal;
                    });
                    const RunParts& pm = parts[medianRun(times)];
                    BenchRow r = baseRow("ours", m, rep);
                    r.family = "ours" + famSfx;
                    r.label = "ours Q" + std::to_string(fp.nQ) + "/R" + std::to_string(fp.nR) + "/B" + std::to_string(fp.nBatches) +
                              sfxShell + (rep ? " s" + std::to_string(rep) : "");
                    r.nQ = fp.nQ;
                    r.nR = fp.nR;
                    r.B = fp.nBatches;
                    r.shellRel = fp.shellRel;
                    r.lean = fp.lean;
                    r.fuse = fp.lean && fp.fuse;
                    r.times = times;
                    r.tMed = medianOf(times);
                    r.tSetup = pm.setup;
                    r.tQ = pm.q;
                    r.tR = pm.r;
                    r.tGram = pm.gram;
                    r.tGlue = pm.glue;
                    r.tField = pm.field;
                    r.walks = f.walksQ + f.walksR;
                    r.moves = f.stepsQ + f.stepsR;
                    r.capped = f.capped;
                    r.escaped = f.escaped;
                    r.nQtets = f.nQtets;
                    r.nRtets = f.nRtets;
                    for (int k = 0; k < kRingNV; k++) r.vert[k] = s.J[m](s.R.V[size_t(k)]) + f.gl[m].c[k];
                    r.cse = f.gl[m].se;
                    r.hasC = true;
                    fillErrors(b, m, f.F[m], r);
                    glueEnergy(s.R, f.G, m, f.gl[m], r.E, r.Ese);
                    std::vector<Gram> GL;   // the same Grams over the points near the crack-face lines (not timed)
                    for (int t = 0; t < kRingNT; t++)
                        GL.push_back(cvGram(f.cv[size_t(t)], s.ps[size_t(t)], b.wLineT[size_t(t)], s.R.T[size_t(t)], fp.nBatches));
                    glueEnergy(s.R, GL, m, f.gl[m], r.EL, r.ELse);
                    if (!determ) out("    WARNING: %s: the repeats gave different glued c (not deterministic)\n", r.label.c_str());
                    printRow(out, r);
                    rows.push_back(r);
                }
            }
            // B: the direct ring solve
            for (const auto& gs : gDirect) {
                const int N = gs[0], B = gs.size() > 1 ? gs[1] : defaultBatches(N);
                const uint64_t seedTag = (loose ? 17900 : 17000) + 10*uint64_t(s.ri) + uint64_t(s.ki) + 1000*uint64_t(rep);
                RingDirect d;
                std::vector<RunParts> parts;
                const std::vector<double> times = timedRuns(repeats, [&] {
                    d = runRingDirect(s.O, s.R, s.rc, s.ps, N, B, mixSeed(s.gp.seed, seedTag), bp.shellRel, bp.lean, b.nModes, true);
                    parts.push_back({0.0, 0.0, 0.0, d.secWalks, 0.0, 0.0, d.secField});
                    return d.secTotal;
                });
                const RunParts& pm = parts[medianRun(times)];
                for (int m = 0; m < b.nModes; m++) {
                    BenchRow r = baseRow("direct", m, rep);
                    r.family = "direct" + famSfx;
                    r.label = "direct N" + std::to_string(N) + "/B" + std::to_string(B) + sfxShell + (rep ? " s" + std::to_string(rep) : "");
                    r.N = N;
                    r.B = B;
                    r.shellRel = bp.shellRel;
                    r.lean = bp.lean;
                    r.times = times;
                    r.tMed = medianOf(times);
                    r.tWalks = pm.walks;
                    r.tField = pm.field;
                    r.walks = d.run.stats.walksRequested;
                    r.moves = d.run.stats.totalSteps;
                    r.capped = d.run.stats.walksCapped;
                    r.escaped = d.run.stats.walksEscaped;
                    for (int k = 0; k < kRingNV; k++) r.vert[k] = d.P[m].H[size_t(k)];
                    fillErrors(b, m, d.F[m], r);
                    std::vector<double> wA(b.w), wL(b.wLine);
                    wA.resize(d.run.X.size(), 0.0);
                    wL.resize(d.run.X.size(), 0.0);
                    solveEnergy(d.run.sums, wA, m, r.E, r.Ese);
                    solveEnergy(d.run.sums, wL, m, r.EL, r.ELse);
                    printRow(out, r);
                    rows.push_back(r);
                }
            }
            // check 4's current path (runWalk, two runs per tet; never the lean walks or the loose shell)
            for (const auto& gs : gCur) {
                const int N = gs[0], B = gs.size() > 1 ? gs[1] : defaultBatches(N);
                const uint64_t tagBase = 12000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki) + 100000*uint64_t(rep);
                RingCurrent c;
                std::vector<RunParts> parts;
                const std::vector<double> times = timedRuns(repeats, [&] {
                    c = runRingCurrent(s.R, s.rc, s.ps, s.J, N, B, s.gp.seed, tagBase, b.nModes, true);
                    parts.push_back({c.secSetup, 0.0, 0.0, c.secWalks, c.secGram, c.secGlue, c.secField});
                    return c.secTotal;
                });
                const RunParts& pm = parts[medianRun(times)];
                for (int m = 0; m < b.nModes; m++) {
                    BenchRow r = baseRow("current", m, rep);
                    r.family = "current";
                    r.label = "current N" + std::to_string(N) + "/B" + std::to_string(B) + (rep ? " s" + std::to_string(rep) : "");
                    r.N = N;
                    r.B = B;
                    r.shellRel = kGateShellRel;
                    r.times = times;
                    r.tMed = medianOf(times);
                    r.tSetup = pm.setup;
                    r.tWalks = pm.walks;
                    r.tGram = pm.gram;
                    r.tGlue = pm.glue;
                    r.tField = pm.field;
                    r.walks = c.walks;
                    r.moves = c.steps;
                    r.capped = c.capped;
                    r.escaped = c.escaped;
                    for (int k = 0; k < kRingNV; k++) r.vert[k] = s.J[m](s.R.V[size_t(k)]) + c.gl[m].c[k];
                    r.cse = c.gl[m].se;
                    r.hasC = true;
                    fillErrors(b, m, c.F[m], r);
                    glueEnergy(s.R, c.G, kRingH0C + m, c.gl[m], r.E, r.Ese);
                    std::vector<Gram> GL;
                    for (int t = 0; t < kRingNT; t++)
                        GL.push_back(gramRing(c.rec[size_t(t)], c.hc[size_t(t)], b.wLineT[size_t(t)], s.R.T[size_t(t)], true));
                    glueEnergy(s.R, GL, kRingH0C + m, c.gl[m], r.EL, r.ELse);
                    if (N == 2560 && B == 64 && rep == 0 && m == 0) {   // check 4's own run: its glue system as gate_report.md prints it
                        out("    (check 4's own setting and seeds: glue M | r, as in gate_report.md's check 4)\n");
                        for (int i = 0; i < kRingNV; i++)
                            out("      [%s]  [%s]\n", matRow(c.gl[0].M, i, "%+.4f").c_str(), fmtd("%+.5f", c.gl[0].rhs[i]).c_str());
                    }
                    printRow(out, r);
                    rows.push_back(r);
                }
            }
        }
        // the vertex values and the seed spread
        out("  vertex values v0..v6 (g = 1): reference");
        for (int k = 0; k < kRingNV; k++) out(" %+.4f", b.refP[0].H[size_t(k)]);
        out("\n");
        for (size_t i = rowStart; i < rows.size(); i++)
            if (rows[i].g == 0) {
                out("    %-30s", rows[i].label.c_str());
                for (int k = 0; k < kRingNV; k++) out(" %+.4f", rows[i].vert[k]);
                out("   (max err %.4f)\n", rows[i].vErr);
            }
        if (bp.seeds > 1) {
            out("  spread over %d seed sets (H RMS / 2g: mean, sd, min, max):\n", bp.seeds);
            std::vector<std::string> seen;
            for (size_t i = rowStart; i < rows.size(); i++) {
                const std::string key = rows[i].family + "|" + std::to_string(rows[i].nQ) + "|" + std::to_string(rows[i].nR) + "|" +
                                        std::to_string(rows[i].N) + "|" + std::to_string(rows[i].B) + "|" + std::to_string(rows[i].g);
                if (std::find(seen.begin(), seen.end(), key) != seen.end()) continue;
                seen.push_back(key);
                std::vector<double> v;
                for (size_t j = rowStart; j < rows.size(); j++)
                    if (rows[j].family == rows[i].family && rows[j].nQ == rows[i].nQ && rows[j].nR == rows[i].nR &&
                        rows[j].N == rows[i].N && rows[j].B == rows[i].B && rows[j].g == rows[i].g)
                        v.push_back(rows[j].e[0].rmsH());
                double mean = 0, sd = 0;
                for (double x : v) mean += x/double(v.size());
                for (double x : v) sd += (x - mean)*(x - mean)/std::max(1.0, double(v.size()) - 1);
                out("    %-30s %s: %.3f%% +- %.3f%% [%.3f%%, %.3f%%] (%zu)\n", rows[i].label.c_str(), rows[i].g ? "ramp" : "g=1",
                    100*mean, 100*std::sqrt(sd), 100**std::min_element(v.begin(), v.end()), 100**std::max_element(v.begin(), v.end()),
                    v.size());
            }
        }
        for (int m = 0; m < b.nModes; m++) benchSummary(out, rows, b.name, m, summaryCsv);
    }
    // the files (new names only; never gate_report.md, results*.txt or gate_cache/)
    {
        std::ofstream f("ring_perf/" + tag + "_rows.csv", std::ios::binary);
        f << csvHeader();
        for (const BenchRow& r : rows) f << csvRow(r);
    }
    {
        std::ofstream f("ring_perf/" + tag + "_summary.csv", std::ios::binary);
        f << summaryCsv;
    }
    out("\n--ring-bench total %.1f s; %zu rows written to ring_perf/%s_rows.csv\n", secondsSince(tStart), rows.size(), tag.c_str());
    {
        std::ofstream f("ring_perf/" + tag + "_log.txt", std::ios::binary);
        f << out.text;
    }
    return true;
}

// =====================================================================================
// --ring-stress: the nodal ring under other crack types (complete, curved, kinked, several cracks), our per-tet method
// against the GT (one solve on O), with ONE crack representation for both: a triangle mesh (CrackMesh, walk_cache.h)
// =====================================================================================
// Opt-in: nothing here is called by any existing mode. It reads check 4's caches in gate_cache/ read-only (A1, A2, A3),
// writes its own GT caches there under new names (ref_stress_*.bin; never a name starting with ref_check), and its
// outputs into ring_stress/. Every per-tet solve stays in the tet's own xi frame on the reference tet (Delta_xi H = 0):
// the crack mesh in X is clipped to O, each tet's part is mapped into xi by xi = J^-1 (X - X0) and clipped to the
// reference tet; the face data are physical functions evaluated at X = X0 + J xi.
#include <map>
#include <memory>
#include <set>

#include "oneapi/tbb/global_control.h"

namespace {

constexpr double kPi = std::numbers::pi;

// ---- the solid angle of an oriented triangle mesh ----
// omega(x): Van Oosterom-Strackee per triangle, signed so that it tends to +2 pi on the + side ((b - a) x (c - a)) of a
// flat triangle; harmonic off the mesh, jump 4 pi across it. grad(x): Biot-Savart over the boundary edges (interior
// edges cancel), checked against finite differences by the self-test.
struct SolidAngleMesh {
    std::vector<Vec3> V;
    std::vector<std::array<int, 3>> F;
    std::vector<std::array<int, 2>> E;   // boundary edges, oriented as in their triangle
    bool flat = false;                   // all triangles in one plane: the one-sided limits are exactly +-2 pi
    bool empty() const { return F.empty(); }
    void findBoundary()
    {
        std::map<std::pair<int, int>, std::pair<int, std::array<int, 2>>> cnt;
        for (const auto& f : F)
            for (int k = 0; k < 3; k++) {
                const int a = f[size_t(k)], b = f[size_t((k + 1)%3)];
                auto& e = cnt[{std::min(a, b), std::max(a, b)}];
                e.first++;
                e.second = {a, b};
            }
        E.clear();
        for (const auto& kv : cnt)
            if (kv.second.first == 1) E.push_back(kv.second.second);
    }
    SolidAngleMesh mapped(const Eigen::Matrix3d& A, const Vec3& b) const
    {
        SolidAngleMesh m = *this;
        for (Vec3& v : m.V) v = A*v + b;
        return m;
    }
    double omega(const Vec3& x) const
    {
        double s = 0;
        for (const auto& f : F) {
            const Vec3 r1 = V[size_t(f[0])] - x, r2 = V[size_t(f[1])] - x, r3 = V[size_t(f[2])] - x;
            const double l1 = r1.norm(), l2 = r2.norm(), l3 = r3.norm();
            const double N = r1.dot(r2.cross(r3));
            const double D = l1*l2*l3 + r1.dot(r2)*l3 + r1.dot(r3)*l2 + r2.dot(r3)*l1;
            s += std::atan2(N, D);
        }
        return -2.0*s;
    }
    Vec3 grad(const Vec3& x) const
    {
        Vec3 g = Vec3::Zero();
        for (const auto& e : E) {
            const Vec3 r1 = V[size_t(e[0])] - x, r2 = V[size_t(e[1])] - x;
            const double l1 = r1.norm(), l2 = r2.norm();
            g += r1.cross(r2)*((l1 + l2)/(l1*l2*(l1*l2 + r1.dot(r2))));
        }
        return -g;
    }
};

// ---- an indexed triangle mesh and its conversion to crack triangles ----
struct IMesh {
    std::vector<Vec3> V;
    std::vector<std::array<int, 3>> F;
    std::set<std::pair<int, int>> front;   // undirected front edges (vertex ids)
    void orient(const std::function<Vec3(const Vec3& centroid)>& plusDir)   // wind every triangle towards plusDir
    {
        for (auto& f : F) {
            const Vec3 a = V[size_t(f[0])], b = V[size_t(f[1])], c = V[size_t(f[2])];
            if ((b - a).cross(c - a).dot(plusDir((a + b + c)/3.0)) < 0) std::swap(f[1], f[2]);
        }
    }
    Vec3 centroid(const std::array<int, 3>& f) const { return (V[size_t(f[0])] + V[size_t(f[1])] + V[size_t(f[2])])/3.0; }
};
std::vector<CrackMesh::Tri> crackTris(const IMesh& m, const std::vector<int>& faces, int crack)
{
    std::vector<CrackMesh::Tri> out;
    for (int fi : faces) {
        const auto& f = m.F[size_t(fi)];
        CrackMesh::Tri t;
        t.a = m.V[size_t(f[0])];
        t.b = m.V[size_t(f[1])];
        t.c = m.V[size_t(f[2])];
        t.crack = crack;
        for (int k = 0; k < 3; k++) {
            const int a = f[size_t(k)], b = f[size_t((k + 1)%3)];
            if (m.front.count({std::min(a, b), std::max(a, b)})) t.edge[k] = CrackMesh::kEdgeFront;
        }
        out.push_back(t);
    }
    return out;
}
std::vector<int> allFaces(const IMesh& m)
{
    std::vector<int> f(m.F.size());
    for (size_t i = 0; i < f.size(); i++) f[i] = int(i);
    return f;
}
SolidAngleMesh solidOf(const IMesh& m, const std::vector<int>& faces, bool flat)
{
    SolidAngleMesh s;
    std::map<int, int> remap;
    for (int fi : faces) {
        std::array<int, 3> g{};
        for (int k = 0; k < 3; k++) {
            const int v = m.F[size_t(fi)][size_t(k)];
            auto it = remap.find(v);
            if (it == remap.end()) {
                it = remap.emplace(v, int(s.V.size())).first;
                s.V.push_back(m.V[size_t(v)]);
            }
            g[size_t(k)] = it->second;
        }
        s.F.push_back(g);
    }
    s.flat = flat;
    s.findBoundary();
    return s;
}

// zipper triangulation between two closed rings (vertex ids and their angles in [0, 2 pi), ascending from 0)
void zipRings(const std::vector<int>& A, const std::vector<double>& aa, const std::vector<int>& B, const std::vector<double>& ab,
              std::vector<std::array<int, 3>>& F)
{
    const size_t ma = A.size(), mb = B.size();
    if (ma == 1) {
        for (size_t j = 0; j < mb; j++) F.push_back({A[0], B[j], B[(j + 1)%mb]});
        return;
    }
    if (mb == 1) {
        for (size_t i = 0; i < ma; i++) F.push_back({A[i], B[0], A[(i + 1)%ma]});
        return;
    }
    size_t i = 0, j = 0;
    while (i < ma || j < mb) {
        const double ta = i + 1 < ma ? aa[i + 1] : 2*kPi + aa[0], tb = j + 1 < mb ? ab[j + 1] : 2*kPi + ab[0];
        if (i < ma && (j >= mb || ta <= tb)) {
            F.push_back({A[i%ma], B[j%mb], A[(i + 1)%ma]});
            i++;
        } else {
            F.push_back({A[i%ma], B[j%mb], B[(j + 1)%mb]});
            j++;
        }
    }
}

// a unit vector orthogonal to a (deterministic)
Vec3 orthoOf(const Vec3& a)
{
    const Vec3 t = std::fabs(a.x()) < 0.6 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    return (t - t.dot(a)*a).normalized();
}

// Flat disk (centre, unit normal n, radius rho): nRings concentric rings, the outer one (the front) with nFront vertices;
// every triangle wound towards +n; the front edges flagged
IMesh diskMesh(const Vec3& centre, const Vec3& n, double rho, int nFront, int nRings)
{
    IMesh m;
    const Vec3 e1 = orthoOf(n), e2 = n.cross(e1);
    m.V.push_back(centre);
    std::vector<int> prev = {0};
    std::vector<double> prevA = {0.0};
    for (int k = 1; k <= nRings; k++) {
        const double r = rho*double(k)/double(nRings);
        const int cnt = k == nRings ? nFront : std::max(6, int(std::lround(double(nFront)*double(k)/double(nRings))));
        std::vector<int> ring;
        std::vector<double> ang;
        for (int j = 0; j < cnt; j++) {
            const double t = 2*kPi*double(j)/double(cnt);
            ring.push_back(int(m.V.size()));
            ang.push_back(t);
            m.V.push_back(centre + r*(std::cos(t)*e1 + std::sin(t)*e2));
        }
        zipRings(prev, prevA, ring, ang, m.F);
        if (k == nRings)
            for (int j = 0; j < cnt; j++) {
                const int a = ring[size_t(j)], b = ring[size_t((j + 1)%cnt)];
                m.front.insert({std::min(a, b), std::max(a, b)});
            }
        prev = ring;
        prevA = ang;
    }
    m.orient([n](const Vec3&) { return n; });
    return m;
}
// the fan of the disk's front polygon (the same flat region as diskMesh's): the solid-angle mesh
SolidAngleMesh diskFan(const Vec3& centre, const Vec3& n, double rho, int nFront)
{
    SolidAngleMesh s;
    const Vec3 e1 = orthoOf(n), e2 = n.cross(e1);
    s.V.push_back(centre);
    for (int j = 0; j < nFront; j++) {
        const double t = 2*kPi*double(j)/double(nFront);
        s.V.push_back(centre + rho*(std::cos(t)*e1 + std::sin(t)*e2));   // the same expression as diskMesh's front ring
    }
    for (int j = 0; j < nFront; j++) {
        std::array<int, 3> f = {0, 1 + j, 1 + (j + 1)%nFront};
        const Vec3 a = s.V[0], b = s.V[size_t(f[1])], c = s.V[size_t(f[2])];
        if ((b - a).cross(c - a).dot(n) < 0) std::swap(f[1], f[2]);
        s.F.push_back(f);
    }
    s.flat = true;
    s.findBoundary();
    return s;
}

// Sphere (centre cs, radius R) as latitude-longitude rings about 'axis', with one ring exactly at cos(polar) = cosFront
// (the front ring, its edges flagged); ring spacing h where the ring comes within rFine of the origin, hCoarse elsewhere;
// wound outwards
IMesh sphereMesh(const Vec3& cs, double R, const Vec3& axis, double cosFront, double h, double hCoarse, double rFine,
                 std::vector<double>* ringCos = nullptr)
{
    const Vec3 a = axis.normalized(), e1 = orthoOf(a), e2 = a.cross(e1);
    const double thF = std::acos(std::clamp(cosFront, -1.0, 1.0));
    const Vec3 p = -cs;   // the origin relative to the centre
    const double pa = p.dot(a), pp = (p - pa*a).norm();
    auto ringDist = [&](double t) { return std::hypot(pa - R*std::cos(t), pp - R*std::sin(t)); };
    std::vector<double> th = {0.0};
    while (th.back() < kPi) {
        const double t = th.back();
        const bool fine = ringDist(t) < rFine + h || ringDist(std::min(kPi, t + hCoarse/R)) < rFine + h;
        double next = t + (fine ? h : hCoarse)/R;
        if (t < thF && next > thF - 0.3*h/R) next = thF;   // land exactly on the front ring
        if (next > kPi - 0.3*h/R) next = kPi;
        th.push_back(next);
    }
    IMesh m;
    std::vector<int> prev;
    std::vector<double> prevA;
    int frontRing = -1;
    std::vector<std::vector<int>> rings;
    for (size_t k = 0; k < th.size(); k++) {
        const double t = th[k];
        std::vector<int> ring;
        std::vector<double> ang;
        if (k == 0 || k + 1 == th.size()) {
            ring.push_back(int(m.V.size()));
            ang.push_back(0.0);
            m.V.push_back(cs + R*(k == 0 ? 1.0 : -1.0)*a);
        } else {
            const bool fine = ringDist(t) < rFine + h;
            const double cosT = t == thF ? cosFront : std::cos(t), sinT = std::sqrt(std::max(0.0, 1.0 - cosT*cosT));
            const int cnt = std::max(6, int(std::ceil(2*kPi*R*sinT/(fine ? h : hCoarse))));
            for (int j = 0; j < cnt; j++) {
                const double phi = 2*kPi*double(j)/double(cnt);
                ring.push_back(int(m.V.size()));
                ang.push_back(phi);
                m.V.push_back(cs + R*(cosT*a + sinT*(std::cos(phi)*e1 + std::sin(phi)*e2)));
            }
            if (t == thF) frontRing = int(k);
        }
        if (k > 0) zipRings(prev, prevA, ring, ang, m.F);
        prev = ring;
        prevA = ang;
        rings.push_back(ring);
        if (ringCos) ringCos->push_back(std::cos(t));
    }
    if (frontRing >= 0) {
        const std::vector<int>& r = rings[size_t(frontRing)];
        for (size_t j = 0; j < r.size(); j++) {
            const int x = r[j], y = r[(j + 1)%r.size()];
            m.front.insert({std::min(x, y), std::max(x, y)});
        }
    }
    m.orient([cs](const Vec3& c) { return Vec3(c - cs); });
    return m;
}

// A planar rectangle P0 + s ea + t eb, s in [s0, s1], t in [t0, t1]: a grid with the s-values sv (ascending) and nt
// steps in t, wound towards nPlus; the vertices with s == sFront form the front
IMesh rectMesh(const Vec3& P0, const Vec3& ea, const Vec3& eb, const std::vector<double>& sv, double t0, double t1, int nt,
               const Vec3& nPlus)
{
    IMesh m;
    const size_t ns = sv.size(), mt = size_t(nt) + 1;
    for (double s : sv)
        for (size_t j = 0; j < mt; j++) {
            const double t = t0 + (t1 - t0)*double(j)/double(nt);
            m.V.push_back(P0 + s*ea + t*eb);
        }
    for (size_t i = 0; i + 1 < ns; i++)
        for (size_t j = 0; j + 1 < mt; j++) {
            const int p = int(i*mt + j), q = int((i + 1)*mt + j), r = int((i + 1)*mt + j + 1), s = int(i*mt + j + 1);
            m.F.push_back({p, q, r});
            m.F.push_back({p, r, s});
        }
    m.orient([nPlus](const Vec3&) { return nPlus; });
    return m;
}

// A plane n.X = d bent into a sine wave along w: X = (d + A sin(2 pi t/lambda + phase)) n + s u + t w, s in sv (ascending),
// t in [-ext, ext] in steps of about h; straight along u, so every quad is planar and a column may be long; wound towards
// +n (a graph over the plane); the edges of the column s == sFront flagged as the front
IMesh waveMesh(const Vec3& n, const Vec3& u, const Vec3& w, double d, double A, double lambda, double phase,
               const std::vector<double>& sv, double sFront, double ext, double h)
{
    IMesh m;
    const int nt = int(std::ceil(2*ext/h));
    const size_t ns = sv.size(), mt = size_t(nt) + 1;
    for (double s : sv)
        for (size_t j = 0; j < mt; j++) {
            const double t = -ext + 2*ext*double(j)/double(nt);
            m.V.push_back((d + A*std::sin(2*kPi*t/lambda + phase))*n + s*u + t*w);
        }
    for (size_t i = 0; i + 1 < ns; i++)
        for (size_t j = 0; j + 1 < mt; j++) {
            const int p = int(i*mt + j), q = int((i + 1)*mt + j), r = int((i + 1)*mt + j + 1), s = int(i*mt + j + 1);
            m.F.push_back({p, q, r});
            m.F.push_back({p, r, s});
        }
    for (size_t i = 0; i < ns; i++)
        if (sv[i] == sFront)
            for (size_t j = 0; j + 1 < mt; j++) m.front.insert({int(i*mt + j), int(i*mt + j + 1)});
    m.orient([n](const Vec3&) { return n; });
    return m;
}

// ---- the cases ----
// One term of a control variate: a function exactly harmonic in every tet's xi frame off its crack, with the crack's jump
// (exactSelf: its one-sided limit on its own crack is exactly side x coef, so the walks score that instead of evaluating it)
struct CvTerm {
    enum Kind { HalfPlane, PlaneSign, MeshSign, Solid } kind = HalfPlane;
    int crack = 0;
    double coef = 1;
    bool exactSelf = true;
    Vec3 n = Vec3::Zero(), u = Vec3::Zero();   // HalfPlane: plane n.X = d, J's front line u.X = sJ; PlaneSign: n.X = d
    double d = 0, sJ = 0;
    int surf = -1;                              // MeshSign: the signed-distance surface; Solid: the solid-angle mesh
};
// One face function J (g = 1; Jr: the ramp's, empty without it) with its control variate in the fast path
struct JVariant {
    std::string name, desc;
    std::function<double(const Vec3&)> J, Jr;
    std::vector<CvTerm> cv;
};
struct StressCrack {
    std::string name;
    double g = 1;
    bool halfPlane = false;   // straight front: the ramp's r is the distance to the front line (frame)
    Cut frame;                // half-plane cracks: n, d, u, w, s (front), L = sqrt(2)
    std::function<double(const Vec3&)> dist;   // distance to the crack (Shepard weights)
};
struct StressCase {
    std::string id, plain, desc;
    int priority = 1, index = 0;
    int ringId = 2;
    Vec3 c = Vec3(0.12, -0.08, 0.10);
    std::vector<StressCrack> cracks;
    CrackMesh glob;                     // all cracks' global surfaces (crack ids, front edges), before clipping to O
    CrackMesh walk;                     // glob clipped to O: THE crack of every method
    std::vector<CrackMesh> surf;        // signed-distance surfaces (S_full) for J and the curved indicators
    std::vector<SolidAngleMesh> solid;  // solid-angle meshes (X)
    std::vector<JVariant> var;          // var[0]: the main face function
    bool ramp = false;                  // also the ramp g (straight half-plane cracks only)
    std::function<double(const Vec3&)> exact;   // complete cuts (g = 1): the exact answer, the side indicator
    int check4 = -1;                    // A1, A2, A3: check 4's point sets (4096 per tet) and its cached GT (crack K1..K3)
    bool v0Undefined = false;           // F1: the node lies on the front
    bool junction = false;              // E3: the junction line (regions)
    Vec3 jP = Vec3::Zero(), jW = Vec3::Zero();
};

// Level C's face function of a half-plane crack (frame: n, d, u; front line u.X = sJ): theta/pi and the ramp
double levelC(const Cut& f, double sJ, const Vec3& X) { return std::atan2(f.n.dot(X) - f.d, f.u.dot(X) - sJ)/kPi; }
double levelCRamp(const Cut& f, double sJ, const Vec3& X)
{
    const double a = f.u.dot(X) - sJ, b = f.n.dot(X) - f.d;
    return std::min(1.0, std::sqrt(std::hypot(a, b)/kSqrt2))*std::sin(0.5*std::atan2(b, a));
}
double halfPlaneDist(const Cut& f, const Vec3& X, bool through)
{
    const double a = f.u.dot(X) - f.s, b = f.n.dot(X) - f.d;
    return through || a <= 0 ? std::fabs(b) : std::hypot(a, b);
}

// a half-plane crack {n.X = d, u.X <= s} (through: the whole plane) as triangles: a square of the plane larger than O
// cut at the front (front edge flagged), to be clipped to O
std::vector<CrackMesh::Tri> halfPlaneTris(const Cut& f, bool through, int crack)
{
    const double R = 3.0;
    const Vec3 o = f.d*f.n;
    std::vector<Vec3> p = {o - R*f.u - R*f.w, o + R*f.u - R*f.w, o + R*f.u + R*f.w, o - R*f.u + R*f.w};
    std::vector<uint8_t> fl(4, 0);
    if (!through) CrackMesh::clipPolygon(p, fl, f.u, f.s, CrackMesh::kEdgeFront);
    std::vector<CrackMesh::Tri> out;
    CrackMesh::fanPolygon(p, fl, crack, f.n, 0.0, out);
    return out;
}

// O's 8 face planes for clipping (flag kEdgeDomain)
CrackMesh clipToO(const CrackMesh& m)
{
    const ConvexDomain O = octahedron();
    return m.clipped(O.n, O.c, std::vector<uint8_t>(O.n.size(), CrackMesh::kEdgeDomain), 1e-14);
}

// the crack frame of a plane in O: n, d, f, front at t of the section's extent (L = sqrt(2) for the ramp)
Cut stressFrame(const Vec3& n, double d, const Vec3& f, double t)
{
    return buildCutOnDomain(octahedron(), n.normalized(), d, f, kSqrt2, t);
}
Cut frameAtS(const Vec3& n, double d, const Vec3& f, double s)
{
    const Cut c1 = stressFrame(n, d, f, 1.0);
    return stressFrame(n, d, f, (s - c1.sMin)/(c1.sMax - c1.sMin));
}

const char* const kStressIds[] = {"A1", "A1s", "A2", "A3", "B1", "B2", "C1", "C2", "D1", "D2", "D3", "E1", "E2", "E3",
                                  "F1", "F2", "F3", "F4", "G1a", "G1b", "F2b", "D4", "D5"};   // new ids last: the seeds use the index
constexpr int kNStress = int(std::size(kStressIds));
const Vec3 kRing3C(0.35, -0.25, 0.20);

Ring makeRingAt(const Vec3& c, int id, const std::string& name)
{
    Ring R;
    R.id = id;
    R.name = name;
    R.V = {c, Vec3(1, 0, 0), Vec3(-1, 0, 0), Vec3(0, 1, 0), Vec3(0, -1, 0), Vec3(0, 0, 1), Vec3(0, 0, -1)};
    for (int o = 0; o < kRingNT; o++) {
        std::array<int, 4> g = {0, o & 1 ? 2 : 1, o & 2 ? 4 : 3, o & 4 ? 6 : 5};
        TetFrame T = frameOn(kRingTetName[o], g, R.V.data());
        if (T.detJ < 0) {
            std::swap(g[2], g[3]);
            T = frameOn(kRingTetName[o], g, R.V.data());
        }
        R.T.push_back(T);
        for (int k = 1; k < 4; k++) {
            std::array<int, 3> f{};
            int m = 0;
            for (int i = 0; i < 4; i++)
                if (i != k) f[m++] = g[i];
            std::sort(f.begin(), f.end());
            if (std::find(R.innerFaces.begin(), R.innerFaces.end(), f) == R.innerFaces.end()) R.innerFaces.push_back(f);
        }
    }
    return R;
}
Ring stressRing(const StressCase& sc)
{
    if (sc.ringId == 2) return makeRing(2);
    return makeRingAt(sc.c, sc.ringId, "ring 3 (distorted, c = " + fmt3(sc.c, 2) + ")");
}

// Shepard blend of several cracks' face functions: J = sum_j w_j J_j, w_j = d_j^-2 / sum_k d_k^-2 (J_j on crack j)
std::function<double(const Vec3&)> blendJ(std::vector<std::function<double(const Vec3&)>> dist,
                                          std::vector<std::function<double(const Vec3&)>> J)
{
    return [dist = std::move(dist), J = std::move(J)](const Vec3& X) {
        double sw = 0, s = 0;
        for (size_t j = 0; j < J.size(); j++) {
            const double d = dist[j](X);
            if (!(d > 1e-150)) return J[j](X);
            const double w = 1.0/(d*d);
            sw += w;
            s += w*J[j](X);
        }
        return s/sw;
    };
}

// signed in-plane distance to a closed polygon (vertices in order, in the plane of normal n): negative inside
double polygonSignedDist(const std::vector<Vec3>& P, const Vec3& n, const Vec3& X)
{
    const Vec3 e1 = orthoOf(n), e2 = n.cross(e1);
    const Vec2 x(X.dot(e1), X.dot(e2));
    double dmin = kInf;
    bool inside = false;
    const size_t m = P.size();
    for (size_t i = 0, j = m - 1; i < m; j = i++) {
        const Vec2 a(P[i].dot(e1), P[i].dot(e2)), b(P[j].dot(e1), P[j].dot(e2));
        const Vec2 ab = b - a;
        const double t = std::clamp((x - a).dot(ab)/ab.squaredNorm(), 0.0, 1.0);
        dmin = std::min(dmin, (x - (a + t*ab)).norm());
        if ((a.y() > x.y()) != (b.y() > x.y()) && x.x() < (b.x() - a.x())*(x.y() - a.y())/(b.y() - a.y()) + a.x()) inside = !inside;
    }
    return inside ? -dmin : dmin;
}

StressCase makeStressCase(const std::string& id)
{
    StressCase sc;
    sc.id = id;
    for (int k = 0; k < kNStress; k++)
        if (id == kStressIds[k]) sc.index = k;
    std::vector<CrackMesh::Tri> glob;
    const Cut K1 = stressFrame(Vec3(1.0, 0.5, 0.2), 0.08, Vec3(0, 1, 0), 0.5);   // K1: front at t = 0.5
    const Vec3 n1 = K1.n, u1 = K1.u, w1 = K1.w;
    auto addHalfPlane = [&](const std::string& name, const Cut& f, bool through, double g, double sJ) {
        StressCrack c;
        c.name = name;
        c.g = g;
        c.halfPlane = !through;
        c.frame = f;
        c.dist = [f, through](const Vec3& X) { return halfPlaneDist(f, X, through); };
        const int j = int(sc.cracks.size());
        sc.cracks.push_back(c);
        auto t = halfPlaneTris(f, through, j);
        glob.insert(glob.end(), t.begin(), t.end());
        (void)sJ;
        return j;
    };
    auto levelCVariant = [&](const Cut& f, double sJ, bool ramp, int crack, double g) {
        JVariant v;
        v.name = "levelC";
        v.desc = "Level C: J = theta/pi of the straight front" + std::string(ramp ? " (and its ramp)" : "") +
                 "; control variate J_xi (the half-plane's theta/pi in xi)";
        v.J = [f, sJ, g](const Vec3& X) { return g*levelC(f, sJ, X); };
        if (ramp) v.Jr = [f, sJ, g](const Vec3& X) { return g*levelCRamp(f, sJ, X); };
        CvTerm t;
        t.kind = CvTerm::HalfPlane;
        t.crack = crack;
        t.coef = g;
        t.n = f.n;
        t.d = f.d;
        t.u = f.u;
        t.sJ = sJ;
        v.cv.push_back(t);
        return v;
    };
    auto planeSignVariant = [&](const Cut& f, int crack) {
        JVariant v;
        v.name = "indicator";
        v.desc = "the side indicator J = sign(n.X - d) (the exact answer); control variate: the same indicator in xi (R = 0)";
        v.J = [f](const Vec3& X) { return f.n.dot(X) - f.d >= 0 ? 1.0 : -1.0; };
        CvTerm t;
        t.kind = CvTerm::PlaneSign;
        t.crack = crack;
        t.n = f.n;
        t.d = f.d;
        v.cv.push_back(t);
        return v;
    };
    auto addDisk = [&](const std::string& name, const Vec3& centre, const Vec3& n, double rho, int nFront, int nRings, double g) {
        StressCrack c;
        c.name = name;
        c.g = g;
        const int j = int(sc.cracks.size());
        const IMesh m = diskMesh(centre, n, rho, nFront, nRings);
        auto t = crackTris(m, allFaces(m), j);
        glob.insert(glob.end(), t.begin(), t.end());
        sc.solid.push_back(diskFan(centre, n, rho, nFront));
        const int s = int(sc.solid.size()) - 1;
        auto cm = std::make_shared<CrackMesh>(crackTris(m, allFaces(m), j));
        c.dist = [cm](const Vec3& X) { return cm->distance(X); };
        sc.cracks.push_back(c);
        return std::pair<int, int>(j, s);
    };
    auto solidVariant = [&](int crack, int s, double g) {
        JVariant v;
        v.name = "solid";
        v.desc = "solid angle: J = Omega(disk)/(2 pi); control variate Omega_xi of the pulled-back disk";
        const SolidAngleMesh* sm = &sc.solid[size_t(s)];
        (void)sm;
        v.J = {};   // set below (needs the final address of sc.solid)
        CvTerm t;
        t.kind = CvTerm::Solid;
        t.crack = crack;
        t.coef = g;
        t.surf = s;
        t.exactSelf = true;   // flat: the one-sided limits are exactly +-2 pi
        v.cv.push_back(t);
        return v;
    };
    const Vec3 g0 = (Vec3(0.12, -0.08, 0.10) + Vec3(1, 0, 0) + Vec3(0, 1, 0) + Vec3(0, 0, 1))/4.0;   // centroid of T+++ (ring 2)

    if (id == "A1" || id == "A1s" || id == "G1a") {
        sc.plain = id == "G1a" ? "flat crack K1 on the distorted ring 3" : "baseline K1 (one flat crack, straight front)";
        const RingCrack rc = ringCrack(0);
        const int j = addHalfPlane("K1", rc.cutO, false, 1.0, rc.sJ);
        sc.var.push_back(levelCVariant(rc.cutO, rc.sJ, true, j, 1.0));
        sc.ramp = id != "G1a";
        if (id == "A1") sc.check4 = 0;
        if (id == "G1a") {
            sc.ringId = 3;
            sc.c = kRing3C;
            sc.priority = 3;
        }
        if (id != "G1a") {   // code check: the same J with a solid-angle control variate (a big flat polygon behind the front)
            IMesh m;
            const double R = 6.0;
            const Vec3 o = rc.cutO.d*rc.cutO.n + rc.sJ*rc.cutO.u;
            m.V = {o - 2*R*rc.cutO.u - R*rc.cutO.w, o - R*rc.cutO.w, o + R*rc.cutO.w, o - 2*R*rc.cutO.u + R*rc.cutO.w};
            m.F = {{0, 1, 2}, {0, 2, 3}};
            m.orient([n = rc.cutO.n](const Vec3&) { return n; });
            sc.solid.push_back(solidOf(m, allFaces(m), true));
            JVariant v = levelCVariant(rc.cutO, rc.sJ, false, j, 1.0);
            v.name = "levelC/cvSolid";
            v.desc = "Level C with a solid-angle control variate (a 12 x 12 flat rectangle behind the front, in xi): a code "
                     "check, same answer within noise";
            v.cv[0].kind = CvTerm::Solid;
            v.cv[0].surf = int(sc.solid.size()) - 1;
            sc.var.push_back(v);
        }
        sc.desc = "half-plane n = normalize(1, 0.5, 0.2), d = 0.08, front direction (0, 1, 0), front at t = 0.5";
    } else if (id == "A2") {
        sc.plain = "baseline K2 (another flat crack, straight front)";
        const RingCrack rc = ringCrack(1);
        const int j = addHalfPlane("K2", rc.cutO, false, 1.0, rc.sJ);
        sc.var.push_back(levelCVariant(rc.cutO, rc.sJ, true, j, 1.0));
        sc.ramp = true;
        sc.check4 = 1;
        sc.desc = "half-plane n = normalize(0.3, -1, 0.6), d = -0.05, front direction (1, 0.4, 0), front at t = 0.3";
    } else if (id == "A3" || id == "B2") {
        Cut f;
        if (id == "A3") {
            sc.plain = "flat complete cut (K3)";
            const RingCrack rc = ringCrack(2);
            f = rc.cutO;
            sc.check4 = 2;
            sc.desc = "K1's plane cut all the way through O";
        } else {
            sc.plain = "flat complete cut passing 0.001 from the node";
            const Vec3 c2(0.12, -0.08, 0.10);
            f = stressFrame(n1, n1.dot(c2) - 1e-3, Vec3(0, 1, 0), 1.0);
            sc.priority = 2;
            sc.desc = "plane n1.X = n1.c - 0.001 (K1's normal), through O; the node 0.001 on its + side";
        }
        const int j = addHalfPlane(id == "A3" ? "K3" : "B2", f, true, 1.0, 0);
        const double sJ = f.sMax + 0.3;
        JVariant v = levelCVariant(f, sJ, false, j, 1.0);
        v.desc = "Level C with its front line 0.3 beyond O's section (the realistic J, check 4's K3); control variate J_xi";
        sc.var.push_back(v);
        sc.var.push_back(planeSignVariant(f, j));
        sc.exact = [f](const Vec3& X) { return f.n.dot(X) - f.d >= 0 ? 1.0 : -1.0; };
    } else if (id == "B1") {
        sc.plain = "curved complete cut (cylinder)";
        const Vec3 a = Vec3(0.2, 0.3, 1.0).normalized(), c2(0.12, -0.08, 0.10);
        Vec3 e = Vec3(-1, 0.4, 0);
        e = (e - e.dot(a)*a).normalized();
        const Vec3 q0 = c2 + 1.42*e;
        const double rho = 1.3, len = 2.5, h = 0.06;
        IMesh m;
        const Vec3 e1 = orthoOf(a), e2 = a.cross(e1);
        const int nt = int(std::ceil(2*kPi*rho/h)), nz = int(std::ceil(2*len/h));
        for (int i = 0; i <= nz; i++)
            for (int j = 0; j < nt; j++) {
                const double z = -len + 2*len*double(i)/double(nz), t = 2*kPi*double(j)/double(nt);
                m.V.push_back(q0 + z*a + rho*(std::cos(t)*e1 + std::sin(t)*e2));
            }
        for (int i = 0; i < nz; i++)
            for (int j = 0; j < nt; j++) {
                const int p = i*nt + j, q = i*nt + (j + 1)%nt, r = (i + 1)*nt + (j + 1)%nt, s = (i + 1)*nt + j;
                m.F.push_back({p, q, r});
                m.F.push_back({p, r, s});
            }
        m.orient([q0, a](const Vec3& c) { const Vec3 r = c - q0; return Vec3(r - r.dot(a)*a); });   // + = outside
        StressCrack c;
        c.name = "cylinder";
        auto t = crackTris(m, allFaces(m), 0);
        glob.insert(glob.end(), t.begin(), t.end());
        sc.surf.push_back(CrackMesh(crackTris(m, allFaces(m), 0)));
        sc.cracks.push_back(c);
        double aMax = 0;
        for (const Vec3& v : octahedron().v) aMax = std::max(aMax, a.dot(v));
        const double sPsi = aMax + 0.3;
        const CrackMesh* S = nullptr;   // bound below
        (void)S;
        JVariant th;
        th.name = "theta";
        th.desc = "realistic theta-type J = atan2(phi, psi)/pi, phi = signed distance to the tube, psi = a.X - " + fmtd("%.4f", sPsi) +
                  " (a virtual front 0.3 beyond O); control variate: the side indicator sign(phi) in xi";
        JVariant ind;
        ind.name = "indicator";
        ind.desc = "the side indicator J = sign(phi) (the exact answer); control variate the same (R = 0)";
        CvTerm cv;
        cv.kind = CvTerm::MeshSign;
        cv.surf = 0;
        th.cv.push_back(cv);
        ind.cv.push_back(cv);
        sc.var.push_back(th);
        sc.var.push_back(ind);
        sc.desc = "cylinder radius 1.3, axis direction normalize(0.2, 0.3, 1) through c + 1.42 e (the axis misses O); + side "
                  "outside; tube tessellated with h = 0.06";
        // J and exact are bound after sc is final (they capture sc.surf by pointer): see bindStressCase
        sc.jP = a;   // stash: the psi direction and offset for bindStressCase
        sc.jW = Vec3(sPsi, 0, 0);
    } else if (id == "C1" || id == "C2") {
        const bool c1 = id == "C1";
        sc.plain = c1 ? "flat crack, gently curved front" : "flat crack, strongly curved front";
        const double rho = c1 ? 1.5 : 0.45;
        const Vec3 centre = c1 ? Vec3(0.08*n1 + (K1.s - 1.5)*u1) : Vec3(0.08*n1 - 0.35*u1);
        const int nF = c1 ? 256 : 128;
        const auto js = addDisk("disk", centre, n1, rho, nF, 12, 1.0);
        sc.var.push_back(solidVariant(js.first, js.second, 1.0));
        JVariant th = solidVariant(js.first, js.second, 1.0);
        th.name = "theta";
        th.desc = "theta-type J = atan2(n.X - d, psi)/pi, psi = signed in-plane distance to the front polygon; control variate "
                  "Omega_xi of the disk";
        sc.var.push_back(th);
        sc.desc = c1 ? "K1's plane, crack = the disk of radius 1.5 tangent to K1's front (front curvature 0.67), 256 front vertices"
                     : "K1's plane, crack = the disk of radius 0.45 centred 0.35 behind O's centre (front curvature 2.2), 128 front "
                       "vertices";
        sc.jP = centre;
        sc.jW = Vec3(rho, double(nF), 0);
    } else if (id == "D1" || id == "D2" || id == "G1b") {
        const bool d2 = id == "D2";
        sc.plain = d2 ? "strongly curved crack surface (bowl, sphere radius 0.6)"
                      : id == "D1" ? "gently curved crack surface (sphere radius 2)" : "gently curved crack (sphere radius 2) on the distorted ring 3";
        Vec3 cs, axis;
        double R, sF, h;
        if (!d2) {
            R = 2.0;
            cs = (0.08 - 2.0)*n1;
            axis = u1;
            sF = K1.s;
            h = 0.06;
        } else {
            R = 0.6;
            cs = 0.6*Vec3(1, 1, 1).normalized();
            axis = Vec3(1, -1, 0).normalized();
            sF = 0.0;
            h = 0.04;
            sc.priority = 1;
        }
        const double cosF = (sF - axis.dot(cs))/R;
        const IMesh m = sphereMesh(cs, R, axis, cosF, h, 0.25, d2 ? 2.0 : 1.25);
        std::vector<int> crackF, cvF;
        for (size_t i = 0; i < m.F.size(); i++) {
            const Vec3 cen = m.centroid(m.F[i]);
            if (axis.dot(cen) - sF < 0) {
                crackF.push_back(int(i));
                if (cen.norm() < (d2 ? 1.2 + 1e-9 : 1.25)) cvF.push_back(int(i));
            }
        }
        StressCrack c;
        c.name = "sphere";
        auto t = crackTris(m, crackF, 0);
        glob.insert(glob.end(), t.begin(), t.end());
        sc.surf.push_back(CrackMesh(crackTris(m, allFaces(m), 0)));   // the full sphere: phi of J
        sc.solid.push_back(solidOf(m, cvF, false));
        auto cm = std::make_shared<CrackMesh>(crackTris(m, crackF, 0));
        c.dist = [cm](const Vec3& X) { return cm->distance(X); };
        sc.cracks.push_back(c);
        JVariant th;
        th.name = "theta";
        th.desc = "theta-type J = atan2(phi, psi)/pi, phi = signed distance to the full sphere, psi = " +
                  std::string(d2 ? "(x - y)/sqrt(2)" : "u1.X - s1") + "; control variate Omega_xi/(2 pi) of the crack-side "
                  "triangles within " + (d2 ? "1.2" : "1.25") + " of the origin (curved: not exact on the crack, evaluated at "
                  "the walk end)";
        CvTerm cv;
        cv.kind = CvTerm::Solid;
        cv.surf = 0;
        cv.exactSelf = false;
        th.cv.push_back(cv);
        sc.var.push_back(th);
        sc.desc = d2 ? "sphere radius 0.6 centred 0.6 (1,1,1)/sqrt(3); crack = sphere n O n {(x - y)/sqrt(2) < 0}; + side outside; "
                       "h = 0.04"
                     : "sphere radius 2 tangent to K1's plane at 0.08 n1 (centre -1.92 n1); crack = sphere n O n {u1.X < s1}; + "
                       "side outside; h = 0.06 within 1.25 of the origin";
        if (id == "G1b") {
            sc.ringId = 3;
            sc.c = kRing3C;
            sc.priority = 3;
        }
        sc.jP = axis;
        sc.jW = Vec3(sF, 0, 0);
    } else if (id == "D4" || id == "D5") {
        const bool d5 = id == "D5";
        sc.plain = d5 ? "strongly wavy crack (sine wave, amplitude 0.1, wavelength 0.6)"
                      : "wavy crack (sine wave, amplitude 0.1, wavelength 1)";
        // K1's plane bent into a sine wave along its front direction w1, K1's front s1: the front is a sine curve; the phase
        // puts a trough under the node (the surface 0.1 below K1's plane there, the node 0.108 above it)
        const double A = 0.1, lam = d5 ? 0.6 : 1.0, h = d5 ? 0.015 : 0.025, ext = 1.6, sF = K1.s;
        const double ph = -0.5*kPi - 2*kPi*w1.dot(Vec3(0.12, -0.08, 0.10))/lam;
        const IMesh m = waveMesh(n1, u1, w1, 0.08, A, lam, ph, {-ext, sF, ext}, sF, ext, h);
        std::vector<int> crackF, cvF;
        for (size_t i = 0; i < m.F.size(); i++) {
            const Vec3 cen = m.centroid(m.F[i]);
            if (u1.dot(cen) - sF < 0) {
                crackF.push_back(int(i));
                if (std::fabs(w1.dot(cen)) < 1.1) cvF.push_back(int(i));   // |t| > 1 is outside O (inside the unit ball)
            }
        }
        StressCrack c;
        c.name = "wave";
        auto t = crackTris(m, crackF, 0);
        glob.insert(glob.end(), t.begin(), t.end());
        sc.surf.push_back(CrackMesh(crackTris(m, allFaces(m), 0)));   // the whole sheet, past the front: phi of J
        sc.solid.push_back(solidOf(m, cvF, false));
        auto cm = std::make_shared<CrackMesh>(crackTris(m, crackF, 0));
        c.dist = [cm](const Vec3& X) { return cm->distance(X); };
        sc.cracks.push_back(c);
        JVariant th;
        th.name = "theta";
        th.desc = "theta-type J = atan2(phi, psi)/pi, phi = signed distance to the whole wavy sheet, psi = u1.X - s1; control "
                  "variate Omega_xi/(2 pi) of the crack-side strips with |w1.X| < 1.1 (curved: not exact on the crack, evaluated "
                  "at the walk end)";
        CvTerm cv;
        cv.kind = CvTerm::Solid;
        cv.surf = 0;
        cv.exactSelf = false;
        th.cv.push_back(cv);
        sc.var.push_back(th);
        sc.desc = "K1's plane bent into n1.X = 0.08 + 0.1 sin(2 pi w1.X/" + fmtd("%g", lam) + " + phase) (smallest radius of "
                  "curvature " + fmtd("%.3f", 1.0/(A*std::pow(2*kPi/lam, 2))) + ", a trough under the node), crack = {u1.X < s1} "
                  "(K1's front, now a sine curve); tessellated with h = " + fmtd("%g", h) + " along w1";
        sc.jP = u1;
        sc.jW = Vec3(sF, 0, 0);
    } else if (id == "D3") {
        sc.plain = "kinked crack (two flat pieces at 40 degrees)";
        const double sk = K1.sMin + 0.4*(K1.sMax - K1.sMin), al = 40.0*kPi/180.0, af = 0.3, ext = 4.0;
        const Vec3 Pl = 0.08*n1 + sk*u1;
        const Vec3 u2 = std::cos(al)*u1 + std::sin(al)*n1, n2 = std::cos(al)*n1 - std::sin(al)*u1;
        IMesh m1 = rectMesh(Pl, u1, w1, {-ext, -2.0, -1.0, -0.5, 0.0}, -ext, ext, 16, n1);
        IMesh m2 = rectMesh(Pl, u2, w1, {0.0, af, 0.6, 1.0, 2.0, ext}, -ext, ext, 16, n2);
        IMesh m = m1;
        const int off = int(m.V.size());
        for (const Vec3& v : m2.V) m.V.push_back(v);
        for (auto f : m2.F) m.F.push_back({f[0] + off, f[1] + off, f[2] + off});
        const double s2 = u2.dot(Pl) + af;
        std::vector<int> crackF;
        for (size_t i = 0; i < m.F.size(); i++)
            if (u2.dot(m.centroid(m.F[i])) - s2 < 0) crackF.push_back(int(i));
        // the front: facet-2 edges on u2.X = s2 (grid column s = af), between crack triangles and the rest
        for (int j = 0; j < 16; j++) {
            const int a = off + 1*17 + j, b = a + 1;
            m.front.insert({a, b});
        }
        StressCrack c;
        c.name = "kink";
        auto t = crackTris(m, crackF, 0);
        glob.insert(glob.end(), t.begin(), t.end());
        sc.surf.push_back(CrackMesh(crackTris(m, allFaces(m), 0)));
        auto cm = std::make_shared<CrackMesh>(crackTris(m, crackF, 0));
        c.dist = [cm](const Vec3& X) { return cm->distance(X); };
        sc.cracks.push_back(c);
        JVariant th;
        th.name = "theta";
        th.desc = "theta-type J = atan2(phi, psi)/pi, phi = signed distance to the kinked surface (facet 2 continued), psi = "
                  "u2.X - s_front; control variate H1 + Hf - Hl: three half-plane theta/pi terms in xi (exact, jump 2)";
        auto hp = [&](const Vec3& n, double d, const Vec3& u, double s, double coef) {
            CvTerm t2;
            t2.kind = CvTerm::HalfPlane;
            t2.n = n;
            t2.d = d;
            t2.u = u;
            t2.sJ = s;
            t2.coef = coef;
            t2.exactSelf = false;
            return t2;
        };
        th.cv.push_back(hp(n1, 0.08, u1, sk, 1.0));
        th.cv.push_back(hp(n2, n2.dot(Pl), u2, s2, 1.0));
        th.cv.push_back(hp(n2, n2.dot(Pl), u2, u2.dot(Pl), -1.0));
        sc.var.push_back(th);
        sc.desc = "facet 1 = K1's plane behind the kink line u1.X = " + fmtd("%.5f", sk) + " (t = 0.4); facet 2 turned 40 degrees "
                  "towards +n1, front 0.3 beyond the kink";
        sc.jP = u2;
        sc.jW = Vec3(s2, 0, 0);
    } else if (id == "E1") {
        sc.plain = "two parallel cracks with opposite fronts (strengths 1 and 0.5)";
        const Cut fa = stressFrame(n1, -0.14, Vec3(0, 1, 0), 0.65), fb = stressFrame(n1, 0.30, Vec3(0, -1, 0), 0.6);
        const int ja = addHalfPlane("a", fa, false, 1.0, fa.s), jb = addHalfPlane("b", fb, false, 0.5, fb.s);
        JVariant v;
        v.name = "blend";
        v.desc = "Shepard blend (weights 1/d_j^2) of the cracks' Level C functions (x g_j), and of their ramps; control variate "
                 "J_xi,a + 0.5 J_xi,b";
        v.J = blendJ({sc.cracks[0].dist, sc.cracks[1].dist},
                     {[fa](const Vec3& X) { return levelC(fa, fa.s, X); }, [fb](const Vec3& X) { return 0.5*levelC(fb, fb.s, X); }});
        v.Jr = blendJ({sc.cracks[0].dist, sc.cracks[1].dist},
                      {[fa](const Vec3& X) { return levelCRamp(fa, fa.s, X); }, [fb](const Vec3& X) { return 0.5*levelCRamp(fb, fb.s, X); }});
        v.cv = {levelCVariant(fa, fa.s, false, ja, 1.0).cv[0], levelCVariant(fb, fb.s, false, jb, 0.5).cv[0]};
        sc.var.push_back(v);
        sc.ramp = true;
        sc.desc = "crack a: n1.X = -0.14, front direction (0, 1, 0), t = 0.65, g = 1; crack b: n1.X = 0.30, front direction "
                  "(0, -1, 0), t = 0.6, g = 0.5; errors / (2 max g) = / 2";
    } else if (id == "E2" || id == "F3") {
        const bool e2 = id == "E2";
        sc.plain = e2 ? "two cracks inside the same tet" : "tiny crack inside one tet";
        sc.priority = 2;
        if (e2) {
            const auto a = addDisk("disk a", g0 + 0.05*n1, n1, 0.18, 64, 6, 1.0);
            const auto b = addDisk("disk b", g0 - 0.05*n1, n1, 0.18, 64, 6, 1.0);
            JVariant v = solidVariant(a.first, a.second, 1.0);
            v.name = "blend";
            v.desc = "Shepard blend of the two disks' Omega/(2 pi); control variate Omega_a,xi/(2 pi) + Omega_b,xi/(2 pi)";
            v.cv.push_back(solidVariant(b.first, b.second, 1.0).cv[0]);
            sc.var.push_back(v);
            sc.desc = "two parallel disks of radius 0.18, normal n1, centres g0 +- 0.05 n1 (g0 = centroid of T+++), 64 front vertices";
        } else {
            const auto a = addDisk("disk", g0, n1, 0.06, 64, 6, 1.0);
            sc.var.push_back(solidVariant(a.first, a.second, 1.0));
            JVariant th = solidVariant(a.first, a.second, 1.0);
            th.name = "theta";
            th.desc = "theta-type J = atan2(n.X - d, psi)/pi, psi = signed in-plane distance to the front polygon (fails far away); "
                      "control variate Omega_xi";
            sc.var.push_back(th);
            sc.desc = "a disk of radius 0.06, normal n1, centre g0 = centroid of T+++ (inside T+++, no face contact), 64 front "
                      "vertices";
            sc.jP = g0;
            sc.jW = Vec3(0.06, 64.0, 0);
        }
    } else if (id == "E3") {
        sc.plain = "two crossing cracks (X junction)";
        const Vec3 c2(0.12, -0.08, 0.10);
        const Cut fA = stressFrame(n1, 0.08, Vec3(0, 1, 0), 1.0);
        const Vec3 nB = Vec3(-0.4, 0.2, 1.0).normalized();
        const Cut fB = stressFrame(nB, nB.dot(c2) + 0.05, n1.cross(nB), 0.55);
        const int jA = addHalfPlane("A (through)", fA, true, 1.0, 0), jB = addHalfPlane("B", fB, false, 1.0, fB.s);
        JVariant v;
        v.name = "blend";
        v.desc = "Shepard blend of A's indicator and B's Level C; control variate sign_A,xi + J_xi,B";
        v.J = blendJ({sc.cracks[0].dist, sc.cracks[1].dist},
                     {[fA](const Vec3& X) { return fA.n.dot(X) - fA.d >= 0 ? 1.0 : -1.0; },
                      [fB](const Vec3& X) { return levelC(fB, fB.s, X); }});
        v.cv = {planeSignVariant(fA, jA).cv[0], levelCVariant(fB, fB.s, false, jB, 1.0).cv[0]};
        sc.var.push_back(v);
        sc.junction = true;
        const Vec3 lineDir = n1.cross(nB).normalized();
        // a point on both planes: least-norm solution of [n1; nB] X = [d1; dB]
        Eigen::Matrix<double, 2, 3> A;
        A.row(0) = n1.transpose();
        A.row(1) = nB.transpose();
        const Vec3 P0 = A.transpose()*(A*A.transpose()).inverse()*Eigen::Vector2d(fA.d, fB.d);
        // the junction is the half of that line on B (u_B.X <= s_B): from the point where it meets B's front, into B
        const double ul = fB.u.dot(lineDir);
        sc.jP = P0 + ((fB.s - fB.u.dot(P0))/ul)*lineDir;
        sc.jW = ul > 0 ? Vec3(-lineDir) : lineDir;
        sc.desc = "A = K1's plane through O (n1.X = 0.08); B = half-plane n_B = normalize(-0.4, 0.2, 1), d_B = n_B.c + 0.05, front "
                  "direction n1 x n_B, t = 0.55; the data disagree on the junction line (bounded solution, singular gradient)";
    } else if (id == "F1") {
        sc.plain = "front passing exactly through the node";
        sc.priority = 2;
        const Vec3 c2(0.12, -0.08, 0.10);
        const Cut f = frameAtS(n1, n1.dot(c2), Vec3(0, 1, 0), u1.dot(c2));
        const int j = addHalfPlane("F1", f, false, 1.0, f.s);
        sc.var.push_back(levelCVariant(f, f.s, false, j, 1.0));
        sc.v0Undefined = true;
        sc.desc = "K1's normal and front direction moved through c: plane n1.X = n1.c, front u1.X = u1.c (t = " + fmtd("%.4f", f.t) + ")";
    } else if (id == "F2" || id == "F2b") {
        // F2b (diagnostic): the same crack 0.03 from the face, outside the eps-tube (eps = 0.0141) that the quadrature leaves
        // out, to tell the slab's missing energy from anything else
        const double off = id == "F2" ? 1e-3 : 0.03;
        sc.plain = id == "F2" ? "crack 0.001 away from an interior face" : "crack 0.03 away from an interior face (F2's control)";
        sc.priority = id == "F2" ? 2 : 3;
        const Vec3 c2(0.12, -0.08, 0.10), ex(1, 0, 0), ey(0, 1, 0);
        Vec3 m = (ex - c2).cross(ey - c2).normalized();
        if (m.z() < 0) m = -m;
        const Vec3 fdir = c2 - 0.5*(ex + ey);
        const Vec3 cen = (c2 + ex + ey)/3.0;
        const Cut f0 = stressFrame(m, m.dot(c2) + off, fdir, 1.0);
        const Cut f = frameAtS(m, m.dot(c2) + off, fdir, f0.u.dot(cen));
        const int j = addHalfPlane(id, f, false, 1.0, f.s);
        sc.var.push_back(levelCVariant(f, f.s, false, j, 1.0));
        sc.desc = "plane m.X = m.c + " + fmtd("%g", off) + ", m = the normal of T+++'s face (c, e_x, e_y); front through that face's "
                  "centroid projected (t = " + fmtd("%.4f", f.t) + "), crack on the side of the edge e_x-e_y";
    } else if (id == "F4") {
        sc.plain = "front 0.03 from the outer boundary";
        sc.priority = 3;
        // the edge of K1's section from (u, w) ~ (0.4235, -0.5731) to (0.7707, 0.0535); front direction = its outward normal
        const Cut full = stressFrame(n1, 0.08, Vec3(0, 1, 0), 1.0);
        const std::vector<Vec3>& P = full.fullPolygon;
        size_t best = 0;
        double bd = kInf;
        for (size_t i = 0; i < P.size(); i++) {
            const Vec3 &a = P[i], &b = P[(i + 1)%P.size()];
            const double e = std::hypot(a.dot(u1) - 0.4235, a.dot(w1) + 0.5731) + std::hypot(b.dot(u1) - 0.7707, b.dot(w1) - 0.0535);
            const double e2 = std::hypot(b.dot(u1) - 0.4235, b.dot(w1) + 0.5731) + std::hypot(a.dot(u1) - 0.7707, a.dot(w1) - 0.0535);
            if (std::min(e, e2) < bd) {
                bd = std::min(e, e2);
                best = i;
            }
        }
        const Vec3 &a = P[best], &b = P[(best + 1)%P.size()];
        Vec3 nrm = (b - a).cross(n1).normalized();
        Vec3 pc = Vec3::Zero();
        for (const Vec3& p : P) pc += p/double(P.size());
        if (nrm.dot(a - pc) < 0) nrm = -nrm;
        const Cut f = frameAtS(n1, 0.08, nrm, nrm.dot(a) - 0.03);
        const int j = addHalfPlane("F4", f, false, 1.0, f.s);
        sc.var.push_back(levelCVariant(f, f.s, false, j, 1.0));
        sc.desc = "K1's plane, front direction = the outward in-plane normal " + fmt3(nrm) + " of a section edge, front 0.03 inside "
                  "that edge (t = " + fmtd("%.4f", f.t) + ")";
    } else {
        throw std::invalid_argument("--ring-stress: unknown case " + id);
    }
    sc.glob = CrackMesh(std::move(glob));
    sc.walk = clipToO(sc.glob);
    return sc;
}

// Binds the face functions that capture the case's own meshes (by address: sc must not move afterwards)
void bindStressCase(StressCase& sc)
{
    if (sc.id == "B1") {
        const CrackMesh* S = &sc.surf[0];
        const Vec3 a = sc.jP;
        const double sPsi = sc.jW.x();
        sc.var[0].J = [S, a, sPsi](const Vec3& X) { return std::atan2(S->signedDistance(X), a.dot(X) - sPsi)/kPi; };
        sc.var[1].J = [S](const Vec3& X) { return S->signedDistance(X) >= 0 ? 1.0 : -1.0; };
        sc.exact = sc.var[1].J;
        sc.jP = sc.jW = Vec3::Zero();
    } else if (sc.id == "D1" || sc.id == "D2" || sc.id == "G1b" || sc.id == "D3" || sc.id == "D4" || sc.id == "D5") {
        const CrackMesh* S = &sc.surf[0];
        const Vec3 a = sc.jP;
        const double sF = sc.jW.x();
        sc.var[0].J = [S, a, sF](const Vec3& X) { return std::atan2(S->signedDistance(X), a.dot(X) - sF)/kPi; };
        sc.jP = sc.jW = Vec3::Zero();
    } else if (sc.id == "C1" || sc.id == "C2" || sc.id == "F3") {
        const SolidAngleMesh* sm = &sc.solid[0];
        sc.var[0].J = [sm](const Vec3& X) { return sm->omega(X)/(2*kPi); };
        // theta-type: the plane's signed distance and the signed in-plane distance to the front polygon
        std::vector<Vec3> P(sm->V.begin() + 1, sm->V.end());
        const Vec3 n = (sm->V[size_t(sm->F[0][1])] - sm->V[0]).cross(sm->V[size_t(sm->F[0][2])] - sm->V[0]).normalized();
        const double d = n.dot(sm->V[0]);
        sc.var[1].J = [P, n, d](const Vec3& X) { return std::atan2(n.dot(X) - d, polygonSignedDist(P, n, X))/kPi; };
        sc.jP = sc.jW = Vec3::Zero();
    } else if (sc.id == "E2") {
        const SolidAngleMesh *a = &sc.solid[0], *b = &sc.solid[1];
        sc.var[0].J = blendJ({sc.cracks[0].dist, sc.cracks[1].dist},
                             {[a](const Vec3& X) { return a->omega(X)/(2*kPi); }, [b](const Vec3& X) { return b->omega(X)/(2*kPi); }});
    }
}

std::unique_ptr<StressCase> buildStressCase(const std::string& id)
{
    auto sc = std::make_unique<StressCase>(makeStressCase(id));
    bindStressCase(*sc);
    return sc;
}

// ---- one tet's crack in its xi frame ----
// The walk crack's triangles near the tet mapped into xi (xi = J^-1 (X - X0), det J > 0: orientation kept) and clipped to
// the reference tet; the edges on a tet face (crack-face lines) and on the front, mapped back to X for the point refinement
struct StressTet {
    CrackMesh mesh;
    bool cut = false;
    double areaX = 0;
    std::vector<std::array<Vec3, 2>> faceLines, front;
    std::vector<int> cracks;
};
StressTet stressTet(const TetFrame& T, const CrackMesh& walk)
{
    StressTet s;
    Vec3 lo = T.X[0], hi = T.X[0];
    for (const Vec3& p : T.X) {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    std::vector<CrackMesh::Tri> sel;
    for (const auto& t : walk.tris()) {
        const Vec3 tlo = t.a.cwiseMin(t.b).cwiseMin(t.c), thi = t.a.cwiseMax(t.b).cwiseMax(t.c);
        if ((tlo.array() > hi.array() + 1e-12).any() || (thi.array() < lo.array() - 1e-12).any()) continue;
        sel.push_back(t);
    }
    const CrackMesh xi = CrackMesh(std::move(sel)).mapped(T.Jinv, -T.Jinv*T.X[0]);
    const std::vector<Vec3> n = {Vec3(-1, 0, 0), Vec3(0, -1, 0), Vec3(0, 0, -1), Vec3(1, 1, 1)};
    const std::vector<double> c = {0, 0, 0, 1};
    s.mesh = xi.clipped(n, c, std::vector<uint8_t>(4, CrackMesh::kEdgeTet), 1e-14);
    s.cut = !s.mesh.empty();
    for (const auto& t : s.mesh.tris()) {
        const Vec3 A = T.toX(t.a), B = T.toX(t.b), C = T.toX(t.c);
        s.areaX += 0.5*(B - A).cross(C - A).norm();
        const Vec3* v[3] = {&A, &B, &C};
        for (int k = 0; k < 3; k++) {
            if (t.edge[k] & (CrackMesh::kEdgeDomain | CrackMesh::kEdgeTet)) s.faceLines.push_back({*v[k], *v[(k + 1)%3]});
            if (t.edge[k] & CrackMesh::kEdgeFront) s.front.push_back({*v[k], *v[(k + 1)%3]});
        }
        if (std::find(s.cracks.begin(), s.cracks.end(), t.crack) == s.cracks.end()) s.cracks.push_back(t.crack);
    }
    return s;
}

// ringPoints with the general crack: refined to depth 3-4 around every crack-face line and front segment of the tet (X),
// the points within eps (physical distance to the walk crack) dropped
PointSet stressPoints(const TetFrame& T, const StressTet& st, const CrackMesh& walk, int target, double eps)
{
    PointSet ps;
    ps.eps = eps;
    const Tet ref = referenceTet();
    std::vector<std::array<Vec3, 2>> lines = st.faceLines;
    lines.insert(lines.end(), st.front.begin(), st.front.end());
    std::vector<Leaf> leaves = {Leaf{{ref.v[0], ref.v[1], ref.v[2], ref.v[3]}, 0}};
    auto refine = [&](auto pred) {
        std::vector<Leaf> out;
        for (const Leaf& l : leaves) {
            if (!pred(l)) {
                out.push_back(l);
                continue;
            }
            for (const Leaf& c : beyChildren(l)) out.push_back(c);
        }
        leaves.swap(out);
    };
    refine([](const Leaf&) { return true; });
    refine([](const Leaf&) { return true; });
    for (int level = 2; level < 4 && !lines.empty(); level++)
        refine([&](const Leaf& l) {
            if (l.depth != level) return false;
            const Vec3 c = T.toX(0.25*(l.v[0] + l.v[1] + l.v[2] + l.v[3]));
            double diam = 0;
            for (int i = 0; i < 4; i++)
                for (int j = i + 1; j < 4; j++) diam = std::max(diam, (T.J*(l.v[i] - l.v[j])).norm());
            for (const auto& s : lines)
                if (pointSegmentDist(c, s[0], s[1]) < 0.5*diam) return true;
            return false;
        });
    auto counts = [&](double c) {
        std::vector<int> m(leaves.size());
        for (size_t i = 0; i < leaves.size(); i++) m[i] = std::max(1, int(std::lround(c*std::pow(8.0, -0.5*(leaves[i].depth - 2)))));
        return m;
    };
    auto total = [](const std::vector<int>& m) {
        long long s = 0;
        for (int v : m) s += v;
        return s;
    };
    double lo = 1, hi = double(target);
    for (int it = 0; it < 60; it++) {
        const double mid = 0.5*(lo + hi);
        (total(counts(mid)) < target ? lo : hi) = mid;
    }
    const std::vector<int> m = counts(hi);
    uint64_t idx = 1;
    ps.nLeaves = int(leaves.size());
    std::vector<Vec3> cand;
    std::vector<double> candW;
    for (size_t li = 0; li < leaves.size(); li++) {
        const Leaf& l = leaves[li];
        const double wl = leafVolume(l)/double(m[li])*std::fabs(T.detJ);
        ps.leavesByDepth[l.depth]++;
        ps.pointsByDepth[l.depth] += m[li];
        for (int k = 0; k < m[li]; k++, idx++) {
            const Vec3 b = cubeToSimplex(radicalInverse(idx, 2), radicalInverse(idx, 3), radicalInverse(idx, 5));
            cand.push_back(l.v[0] + b.x()*(l.v[1] - l.v[0]) + b.y()*(l.v[2] - l.v[0]) + b.z()*(l.v[3] - l.v[0]));
            candW.push_back(wl);
        }
    }
    std::vector<char> keep(cand.size(), 1);
    tbb::parallel_for(size_t(0), cand.size(), [&](size_t i) { keep[i] = !(walk.distance(T.toX(cand[i])) < eps); });
    for (size_t i = 0; i < cand.size(); i++) {
        ps.allXi.push_back(cand[i]);
        ps.allW.push_back(candW[i]);
        if (!keep[i]) {
            ps.nDropped++;
            continue;
        }
        ps.xi.push_back(cand[i]);
        ps.X.push_back(T.toX(cand[i]));
        ps.w.push_back(candW[i]);
    }
    return ps;
}

// Level C's theta/pi of a half-plane (plane n.X = d, front line u.X = sJ) in a tet's xi frame (xiFront with explicit data)
XiFront xiFrontOf(const TetFrame& T, const Vec3& nX, double d, const Vec3& uX, double sJ)
{
    XiFront f;
    const Vec3 wX = nX.cross(uX);
    const Vec3 jn = T.J.transpose()*nX;
    f.n = jn.normalized();
    f.dn = (d - nX.dot(T.X[0]))/jn.norm();
    const Vec3 wXi = (T.Jinv*wX).normalized();
    f.u = wXi.cross(f.n);
    if (f.u.dot(T.Jinv*uX) < 0) f.u = -f.u;
    f.p0 = T.toXi(d*nX + sJ*uX);
    return f;
}

// One variant's control variate in one tet's xi frame: the sum of its terms, each exactly harmonic in xi off its crack
struct TetCv {
    struct Term {
        CvTerm::Kind kind = CvTerm::HalfPlane;
        int crack = 0;
        double coef = 1;
        bool exactSelf = true;
        XiFront hp;
        Vec3 nX = Vec3::Zero();
        double dX = 0;
        const CrackMesh* surf = nullptr;
        SolidAngleMesh sa;   // in xi
    };
    const TetFrame* T = nullptr;
    std::vector<Term> terms;
    double term(const Term& t, const Vec3& xi, int m) const
    {
        switch (t.kind) {
        case CvTerm::HalfPlane: return t.coef*(m == 0 ? t.hp.J1(xi) : t.hp.J2(xi));
        case CvTerm::PlaneSign: return m == 0 ? t.coef*(t.nX.dot(T->toX(xi)) - t.dX >= 0 ? 1.0 : -1.0) : 0.0;
        case CvTerm::MeshSign: return m == 0 ? t.coef*(t.surf->signedDistance(T->toX(xi)) >= 0 ? 1.0 : -1.0) : 0.0;
        case CvTerm::Solid: return m == 0 ? t.coef*t.sa.omega(xi)/(2*kPi) : 0.0;
        }
        return 0.0;
    }
    double value(const Vec3& xi, int m) const
    {
        double v = 0;
        for (const Term& t : terms) v += term(t, xi, m);
        return v;
    }
    Vec3 grad(const Vec3& xi, int m) const
    {
        Vec3 g = Vec3::Zero();
        for (const Term& t : terms) {
            if (t.kind == CvTerm::HalfPlane) g += t.coef*(m == 0 ? t.hp.grad1(xi) : t.hp.grad2(xi));
            else if (t.kind == CvTerm::Solid && m == 0) g += (t.coef/(2*kPi))*t.sa.grad(xi);
        }
        return g;
    }
    // the control variate's part of the crack data: landing on crack j on side s (walk end xe, closest point q, both xi)
    double atCrack(int j, int s, const Vec3& xe, const Vec3& q, int m) const
    {
        double v = 0;
        for (const Term& t : terms) {
            if (t.crack == j && t.exactSelf) {
                if (m == 0) v += t.coef*double(s);
                else if (t.kind == CvTerm::HalfPlane) v += t.coef*double(s)*std::sqrt(t.hp.r(q)/kSqrt2);
            } else {
                v += term(t, xe, m);
            }
        }
        return v;
    }
};
TetCv tetCv(const TetFrame& T, const JVariant& v, const StressCase& sc)
{
    TetCv c;
    c.T = &T;
    for (const CvTerm& t : v.cv) {
        TetCv::Term u;
        u.kind = t.kind;
        u.crack = t.crack;
        u.coef = t.coef;
        u.exactSelf = t.exactSelf;
        if (t.kind == CvTerm::HalfPlane) u.hp = xiFrontOf(T, t.n, t.d, t.u, t.sJ);
        else if (t.kind == CvTerm::PlaneSign) {
            u.nX = t.n;
            u.dX = t.d;
        } else if (t.kind == CvTerm::MeshSign) u.surf = &sc.surf[size_t(t.surf)];
        else u.sa = sc.solid[size_t(t.surf)].mapped(T.Jinv, -T.Jinv*T.X[0]);
        c.terms.push_back(std::move(u));
    }
    return c;
}

// the crack data of the target problem: +-g_j on crack j (g = 1), times min(1, sqrt(r/L)) for the ramp (r = the physical
// distance to the crack's straight front line), at the physical point X
double stressTarget(const StressCase& sc, int j, int s, const Vec3& X, int m)
{
    const StressCrack& c = sc.cracks[size_t(j)];
    if (m == 0 || !c.halfPlane) return double(s)*c.g;
    return double(s)*c.g*c.frame.g(X, GType::SqrtFront);
}

// ---- the solve on O (the GT and the direct ring solve): the walk crack, 8 Neumann faces, mirror first ball ----
// K = 2: g = 1 and the ramp (= g = 1 for a case without the ramp, to keep the cache format). Points: the tets' points (H,
// grad H), then the 7 vertex probes (H).
RefRun stressSolveO(const StressCase& sc, const Ring& R, const std::vector<PointSet>& ps, int N, int B, uint64_t seed, int maxLen)
{
    RefRun r;
    for (const PointSet& p : ps) r.X.insert(r.X.end(), p.X.begin(), p.X.end());
    r.nA = r.X.size();
    const std::vector<Vec3> probes = ringProbes(R);
    r.X.insert(r.X.end(), probes.begin(), probes.end());
    r.nF = probes.size();
    std::vector<char> grad(r.X.size(), 1);
    std::fill(grad.begin() + std::ptrdiff_t(r.nA), grad.end(), char(0));
    QueryPoints q;
    q.add(r.X, &grad, seed);
    const ConvexDomain O = octahedron();
    MeshJob jb;
    jb.dom = &O;
    jb.mesh = &sc.walk;
    jb.pts = &q;
    jb.nWalks = N;
    jb.nBatches = B;
    jb.K = 2;
    jb.wp = ringReferenceWalkParams(O);
    jb.wp.maxWalkLength = maxLen;
    const bool ramp = sc.ramp;
    jb.values = [&sc, ramp](const MeshTerminal& w, double* u) {
        u[0] = stressTarget(sc, w.crack, w.side, w.x, 0);
        u[1] = ramp ? stressTarget(sc, w.crack, w.side, w.x, 1) : u[0];
    };
    std::vector<SolveStats> st;
    r.sums = std::move(solveMultiMeshJobs({jb}, &st)[0]);
    r.stats = st[0];
    return r;
}

uint64_t stressMeshHash(const StressCase& sc)
{
    std::vector<double> v;
    for (const auto& t : sc.walk.tris()) {
        for (const Vec3* p : {&t.a, &t.b, &t.c, &t.n}) v.insert(v.end(), {p->x(), p->y(), p->z()});
        v.push_back(double(t.crack));
        v.push_back(double(t.edge[0] + 8*t.edge[1] + 64*t.edge[2]));
    }
    for (const StressCrack& c : sc.cracks) {
        v.push_back(c.g);
        if (c.halfPlane) v.insert(v.end(), {c.frame.n.x(), c.frame.n.y(), c.frame.n.z(), c.frame.d, c.frame.u.x(), c.frame.u.y(),
                                            c.frame.u.z(), c.frame.s, c.frame.L});
    }
    v.push_back(sc.ramp ? 1.0 : 0.0);
    return fnv1a(v.data(), v.size()*sizeof(double));
}

const char* const kStressSolverVersion = "ring stress GT: own loop with a CrackMesh crack (walk_cache.cpp solveMultiMeshJobs) v1; "
                                         "K = 2 (g = 1, ramp or g = 1 again)";

// The GT of a case at its points: from the cache gate_cache/<file> if its key matches, else a timing pilot (N/100) and the
// full run, written to the cache (useCache false: neither read nor written)
RefRun stressGT(Report& out, const StressCase& sc, const Ring& R, const std::vector<PointSet>& ps, int N, int B, uint64_t seedTag,
                int maxLen, const std::string& file, bool useCache, double* pilotExtrap = nullptr, bool pilotOnly = false,
                bool skipPilot = false, bool cacheOnly = false)
{
    const GateParams gp = gateParams(false);
    const uint64_t seed = mixSeed(gp.seed, seedTag);
    RefRun r;
    for (const PointSet& p : ps) r.X.insert(r.X.end(), p.X.begin(), p.X.end());
    r.nA = r.X.size();
    const std::vector<Vec3> probes = ringProbes(R);
    r.X.insert(r.X.end(), probes.begin(), probes.end());
    r.nF = probes.size();
    const ConvexDomain O = octahedron();
    const WalkParams wp = ringReferenceWalkParams(O);
    char buf[256];
    std::snprintf(buf, sizeof buf, "%016llx | points %zu (tets %zu, probes %zu), hash %016llx", (unsigned long long)stressMeshHash(sc),
                  r.X.size(), r.nA, r.nF, (unsigned long long)fnv1a(r.X.data(), r.X.size()*sizeof(Vec3)));
    std::string key = "ring stress GT on O | case " + sc.id + " | " + R.name + ", c " + fmtg3(R.V[0]) + " | walk mesh " +
                      std::to_string(sc.walk.size()) + " triangles, hash " + buf + " | seed " + std::to_string(seed) + " | N " +
                      std::to_string(N) + " | batches " + std::to_string(B) + " | O";
    for (size_t k = 0; k < O.n.size(); k++) key += " " + fmtg3(O.n[k]) + "." + fmtg(O.c[k]);
    key += ", L_O " + fmtg(O.L) + " | walk: absorbShellRel " + fmtg(wp.absorbShellRel) + ", maxWalkLength " + std::to_string(maxLen) +
           ", mirrorGradBall " + (wp.mirrorGradBall ? "1" : "0") + ", kGradFirstBallShrink " + fmtg(kGradFirstBallShrink) + " | " +
           kStressSolverVersion;
    r.key = key;
    r.file = (std::filesystem::path("gate_cache")/file).string();
    if (file.rfind("ref_check", 0) == 0) throw std::logic_error("--ring-stress never writes a ref_check* cache");
    if (useCache && readRefCache(r)) {
        r.cached = true;
        out("  GT loaded from %s (key matches); original run: %s\n", r.file.c_str(), r.stats.summary().c_str());
        return r;
    }
    if (cacheOnly) {   // the ring view: nothing is solved here
        r.ok = false;
        return r;
    }
    r.nPilot = B*std::max(1, N/(100*B));
    if (!skipPilot) {
        const RefRun pr = stressSolveO(sc, R, ps, r.nPilot, B, mixSeed(gp.seed, seedTag + 100), maxLen);
        r.pilotSeconds = pr.stats.seconds;
        r.pilotExtrapolated = pr.stats.seconds*double(N)/double(r.nPilot);
        out("  GT timing pilot, N = %d: %s -> extrapolated to N = %d: %.0f s\n", r.nPilot, pr.stats.summary().c_str(), N,
            r.pilotExtrapolated);
    }
    if (pilotExtrap) *pilotExtrap = r.pilotExtrapolated;
    if (pilotOnly) {
        r.ok = false;
        return r;
    }
    const RefRun full = stressSolveO(sc, R, ps, N, B, seed, maxLen);
    r.sums = full.sums;
    r.stats = full.stats;
    out("  GT run: %s\n", r.stats.summary().c_str());
    if (useCache) {
        std::error_code ec;
        std::filesystem::create_directories("gate_cache", ec);
        out("  %s %s\n", writeRefCache(r) ? "GT cached to" : "COULD NOT WRITE the GT cache", r.file.c_str());
    }
    return r;
}

// ---- ours: the per-tet method on a general crack ----
struct StressOursParams {
    int nQ = 32, nR = 16, B = 8;
    bool cv = true;        // control variates (the fast path); false: plain per-tet solves (H0: crack data, faces J; P_k:
                           // crack 0, faces lambda_k) in every tet
    int maxLen = 4096;
};
struct StressOurs {
    std::vector<std::pair<int, int>> vms;          // the (variant, g mode) pairs solved
    std::vector<std::vector<TetCv>> cv;            // [variant][tet]
    std::vector<OwnRun> runP, runR;                // per tet: Q (cv) or P (plain); R (cv) or H0 (plain)
    std::vector<char> hasP, hasR;
    std::vector<std::vector<int>> comp;            // [tet][vm index]: the R component, -1: R = 0 exactly (no walks)
    std::vector<std::vector<std::vector<double>>> cvVal;   // [vm][tet][point]
    std::vector<std::vector<std::vector<Vec3>>> cvGrad;
    std::vector<std::vector<Gram>> G;              // [vm][tet], layout H0, H0, P_0..P_3 (glueCheck4's h = 0)
    std::vector<Glue4> gl;                         // [vm]
    std::vector<std::vector<FieldAt>> F;           // [vm][tet]
    std::vector<double> maxFace, maxCrack;         // [vm]: max |R's data| on the face grids / at the crack vertices
    double secSetup = 0, secP = 0, secR = 0, secGram = 0, secGlue = 0, secField = 0, secTotal = 0;
    long long walksP = 0, walksR = 0, stepsP = 0, stepsR = 0, capped = 0, escaped = 0;
    int nPtets = 0, nRtets = 0;
    int vmIndex(int v, int m) const
    {
        for (size_t k = 0; k < vms.size(); k++)
            if (vms[k].first == v && vms[k].second == m) return int(k);
        return -1;
    }
};

// The R residual of one variant and g mode in one tet: max |J - C| on a 25 x 25 barycentric grid per face (points within
// 1e-9 of the tet's crack skipped) and max |target - C| at the crack triangles' vertices and centroids (both sides)
void stressResidual(const StressTet& st, const TetFrame& T, const JVariant& v, const TetCv& cv, int m, const StressCase& sc,
                    double& maxFace, double& maxCrack)
{
    maxFace = maxCrack = 0;
    const Tet ref = referenceTet();
    for (int o = 0; o < 4; o++) {
        std::array<int, 3> fv{};
        int k = 0;
        for (int i = 0; i < 4; i++)
            if (i != o) fv[size_t(k++)] = i;
        for (int i = 0; i <= 24; i++)
            for (int j = 0; i + j <= 24; j++) {
                const double l1 = i/24.0, l2 = j/24.0, l0 = 1.0 - l1 - l2;
                const Vec3 xi = l0*ref.v[size_t(fv[0])] + l1*ref.v[size_t(fv[1])] + l2*ref.v[size_t(fv[2])];
                if (st.cut && st.mesh.distance(xi) < 1e-9) continue;
                const Vec3 X = T.toX(xi);
                maxFace = std::max(maxFace, std::fabs((m ? v.Jr(X) : v.J(X)) - cv.value(xi, m)));
            }
    }
    for (const auto& t : st.mesh.tris())
        for (const Vec3& p : {t.a, t.b, t.c, Vec3((t.a + t.b + t.c)/3.0)})
            for (int s : {-1, 1})
                maxCrack = std::max(maxCrack, std::fabs(stressTarget(sc, t.crack, s, T.toX(p), m) - cv.atCrack(t.crack, s, p, p, m)));
}

StressOurs runStressOurs(const StressCase& sc, const Ring& R, const std::vector<StressTet>& stt, const std::vector<PointSet>& ps,
                         const StressOursParams& p, uint64_t seed, uint64_t tagBase, bool withFields)
{
    using clk = std::chrono::steady_clock;
    const auto tStart = clk::now();
    const int B = p.B;
    if (B < 3 || p.nQ%B || p.nR%B || p.nQ/B < 2 || p.nR/B < 2)
        throw std::invalid_argument("--ring-stress ours: at least 3 batches, and nQ, nR multiples of the batch count with at least "
                                    "2 walks per batch");
    StressOurs o;
    const size_t nV = sc.var.size();
    for (size_t v = 0; v < nV; v++)
        for (int m = 0; m < (sc.ramp ? 2 : 1); m++)
            if (m == 0 || sc.var[v].Jr) o.vms.push_back({int(v), m});
    const size_t nVm = o.vms.size();
    o.cv.assign(nV, {});
    if (p.cv)
        for (size_t v = 0; v < nV; v++)
            for (int t = 0; t < kRingNT; t++) o.cv[v].push_back(tetCv(R.T[size_t(t)], sc.var[v], sc));
    o.comp.assign(kRingNT, std::vector<int>(nVm, -1));
    o.maxFace.assign(nVm, 0.0);
    o.maxCrack.assign(nVm, 0.0);
    std::vector<int> KR(kRingNT, 0);
    std::vector<double> rf(kRingNT*nVm, 0.0), rc(kRingNT*nVm, 0.0);   // each tet's and component's residuals (in parallel)
    if (p.cv)
        tbb::parallel_for(size_t(0), kRingNT*nVm, [&](size_t q) {
            const size_t t = q/nVm, k = q%nVm;
            stressResidual(stt[t], R.T[t], sc.var[size_t(o.vms[k].first)], o.cv[size_t(o.vms[k].first)][t], o.vms[k].second, sc, rf[q], rc[q]);
        });
    for (int t = 0; t < kRingNT; t++)
        for (size_t k = 0; k < nVm; k++) {
            const size_t q = size_t(t)*nVm + k;
            o.maxFace[k] = std::max(o.maxFace[k], rf[q]);
            o.maxCrack[k] = std::max(o.maxCrack[k], rc[q]);
            if (!p.cv || rf[q] > 1e-12 || rc[q] > 1e-12) o.comp[size_t(t)][k] = KR[size_t(t)]++;
        }
    if (*std::max_element(KR.begin(), KR.end()) > 16) throw std::invalid_argument("--ring-stress: more than 16 R components");
    const Tet ref = referenceTet();
    WalkParams wpP = gateWalkParams(allDirichlet(), true, true), wpR = gateWalkParams(allDirichlet(), true, false);
    wpP.maxWalkLength = wpR.maxWalkLength = p.maxLen;
    const HarmonicExtension ext0 = [](const Vec3&, double* phi) {
        for (int k = 0; k < 4; k++) phi[k] = 0.0;
    };
    const HarmonicExtension extL = [](const Vec3& xi, double* phi) {   // lambda_k: linear, its own harmonic extension
        phi[0] = 1.0 - xi.sum();
        phi[1] = xi[0];
        phi[2] = xi[1];
        phi[3] = xi[2];
    };
    const MeshTerminalValues tvQ = [](const MeshTerminal& w, double* u) {   // Q_k: crack -lambda_k(xi), faces 0
        if (w.face >= 0) {
            for (int k = 0; k < 4; k++) u[k] = 0.0;
            return;
        }
        const double l[4] = {1.0 - w.x.sum(), w.x[0], w.x[1], w.x[2]};
        for (int k = 0; k < 4; k++) u[k] = -l[k];
    };
    const MeshTerminalValues tvP = [](const MeshTerminal& w, double* u) {   // P_k (plain): crack 0, faces lambda_k(xi)
        if (w.face < 0) {
            for (int k = 0; k < 4; k++) u[k] = 0.0;
            return;
        }
        const double l[4] = {1.0 - w.x.sum(), w.x[0], w.x[1], w.x[2]};
        for (int k = 0; k < 4; k++) u[k] = l[k];
    };
    std::vector<QueryPoints> qP(kRingNT), qR(kRingNT);
    std::vector<MeshJob> jobsP, jobsR;
    std::vector<int> tetP, tetR;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame* T = &R.T[size_t(t)];
        const StressTet& st = stt[size_t(t)];
        const std::vector<char> grad(ps[size_t(t)].size(), 1);
        const uint64_t sP = mixSeed(seed, tagBase + 10*uint64_t(t));
        MeshJob jb;
        jb.tet = &ref;
        jb.mesh = &st.mesh;
        jb.nBatches = B;
        if (st.cut || !p.cv) {
            qP[size_t(t)].add(ps[size_t(t)].xi, &grad, sP);
            jb.pts = &qP[size_t(t)];
            jb.nWalks = p.nQ;
            jb.K = 4;
            jb.wp = wpP;
            jb.values = p.cv ? tvQ : tvP;
            jb.ext = p.cv ? &ext0 : &extL;
            jobsP.push_back(jb);
            tetP.push_back(t);
        }
        if (KR[size_t(t)] > 0) {
            qR[size_t(t)].add(ps[size_t(t)].xi, &grad, mixSeed(sP, 77));
            std::vector<std::pair<int, int>> comps;   // component -> (variant, mode)
            std::vector<const TetCv*> cvs;
            for (size_t k = 0; k < nVm; k++)
                if (o.comp[size_t(t)][k] >= 0) {
                    comps.push_back(o.vms[k]);
                    cvs.push_back(p.cv ? &o.cv[size_t(o.vms[k].first)][size_t(t)] : nullptr);
                }
            jb.pts = &qR[size_t(t)];
            jb.nWalks = p.nR;
            jb.K = int(comps.size());
            jb.wp = wpR;
            jb.ext = nullptr;
            jb.values = [T, &sc, comps, cvs](const MeshTerminal& w, double* u) {
                if (w.face >= 0) {   // faces: J - C
                    const Vec3 X = T->toX(w.x);
                    for (size_t k = 0; k < comps.size(); k++) {
                        const JVariant& v = sc.var[size_t(comps[k].first)];
                        const int m = comps[k].second;
                        u[k] = (m ? v.Jr(X) : v.J(X)) - (cvs[k] ? cvs[k]->value(w.x, m) : 0.0);
                    }
                } else {             // the crack: the target data - C (the self terms exact where known)
                    const Vec3 X = T->toX(w.x);
                    for (size_t k = 0; k < comps.size(); k++) {
                        const int m = comps[k].second;
                        u[k] = stressTarget(sc, w.crack, w.side, X, m) - (cvs[k] ? cvs[k]->atCrack(w.crack, w.side, w.xEnd, w.x, m) : 0.0);
                    }
                }
            };
            jobsR.push_back(jb);
            tetR.push_back(t);
        }
    }
    o.runP.assign(kRingNT, OwnRun{});
    o.runR.assign(kRingNT, OwnRun{});
    o.hasP.assign(kRingNT, 0);
    o.hasR.assign(kRingNT, 0);
    o.secSetup = secondsSince(tStart);
    auto run = [&](const std::vector<MeshJob>& jobs, const std::vector<int>& tets, bool isP) {
        const auto t0 = clk::now();
        std::vector<SolveStats> st;
        std::vector<MultiSums> sums;
        if (!jobs.empty()) sums = solveMultiMeshJobs(jobs, &st);
        for (size_t j = 0; j < jobs.size(); j++) {
            OwnRun& r = isP ? o.runP[size_t(tets[j])] : o.runR[size_t(tets[j])];
            r.sums = std::move(sums[j]);
            r.stats = st[j];
            (isP ? o.hasP : o.hasR)[size_t(tets[j])] = 1;
            (isP ? o.walksP : o.walksR) += st[j].walksRequested;
            (isP ? o.stepsP : o.stepsR) += st[j].totalSteps;
            o.capped += st[j].walksCapped;
            o.escaped += st[j].walksEscaped;
        }
        return secondsSince(t0);
    };
    o.secP = run(jobsP, tetP, true);
    o.secR = run(jobsR, tetR, false);
    o.nPtets = int(jobsP.size());
    o.nRtets = int(jobsR.size());
    // the control variates' values and gradients at the points (exact)
    auto t0 = clk::now();
    o.cvVal.assign(nVm, std::vector<std::vector<double>>(kRingNT));
    o.cvGrad.assign(nVm, std::vector<std::vector<Vec3>>(kRingNT));
    for (size_t k = 0; k < nVm; k++)
        for (int t = 0; t < kRingNT; t++) {
            const size_t n = ps[size_t(t)].size();
            o.cvVal[k][size_t(t)].assign(n, 0.0);
            o.cvGrad[k][size_t(t)].assign(n, Vec3::Zero());
            if (!p.cv) continue;
            const TetCv& c = o.cv[size_t(o.vms[k].first)][size_t(t)];
            const int m = o.vms[k].second;
            tbb::parallel_for(size_t(0), n, [&](size_t i) {
                o.cvVal[k][size_t(t)][i] = c.value(ps[size_t(t)].xi[i], m);
                o.cvGrad[k][size_t(t)][i] = c.grad(ps[size_t(t)].xi[i], m);
            });
        }
    // the Grams (physical energy) and the glue
    o.G.assign(nVm, {});
    const std::vector<double> u4[4] = {unitCoef(4, 0), unitCoef(4, 1), unitCoef(4, 2), unitCoef(4, 3)};
    for (size_t k = 0; k < nVm; k++)
        for (int t = 0; t < kRingNT; t++) {
            const TetFrame& T = R.T[size_t(t)];
            const std::array<Vec3, 4> gN = gradN(T);
            const OwnRun* pr = o.hasP[size_t(t)] ? &o.runP[size_t(t)] : nullptr;
            const int ci = o.comp[size_t(t)][k];
            const OwnRun* rr = ci >= 0 ? &o.runR[size_t(t)] : nullptr;
            const std::vector<double> uR = rr ? unitCoef(rr->sums.K, ci) : std::vector<double>{};
            const std::vector<Vec3>& cg = o.cvGrad[k][size_t(t)];
            const bool known = p.cv;
            o.G[k].push_back(gramFromBatches(ps[size_t(t)].w, B, 6, [&](size_t i, int b, Vec3* g) {
                Vec3 h = cg[i];
                if (rr) h += comboBatchGrad(rr->sums, i, b, uR);
                g[0] = g[1] = T.JinvT*h;
                for (int q = 0; q < 4; q++)
                    g[2 + q] = (known ? gN[size_t(q)] : Vec3::Zero()) +
                               (pr ? Vec3(T.JinvT*comboBatchGrad(pr->sums, i, b, u4[q])) : Vec3::Zero());
            }));
        }
    o.secGram = secondsSince(t0);
    t0 = clk::now();
    for (size_t k = 0; k < nVm; k++) o.gl.push_back(glueRun4(R, o.G[k], 0));
    o.secGlue = secondsSince(t0);
    if (withFields) {
        t0 = clk::now();
        o.F.assign(nVm, {});
        for (size_t k = 0; k < nVm; k++)
            for (int t = 0; t < kRingNT; t++) {
                const TetFrame& T = R.T[size_t(t)];
                const PointSet& pst = ps[size_t(t)];
                const std::array<Vec3, 4> gN = gradN(T);
                const OwnRun* pr = o.hasP[size_t(t)] ? &o.runP[size_t(t)] : nullptr;
                const int ci = o.comp[size_t(t)][k];
                const OwnRun* rr = ci >= 0 ? &o.runR[size_t(t)] : nullptr;
                const std::vector<double> uR = rr ? unitCoef(rr->sums.K, ci) : std::vector<double>{};
                const Eigen::Vector4d cc = localValues(T, o.gl[k].c);
                std::vector<Eigen::Vector4d> cJ;
                for (const Vec7& x : o.gl[k].jack) cJ.push_back(localValues(T, x));
                const bool known = p.cv;
                FieldAt f;
                f.resize(pst.size());
                tbb::parallel_for(size_t(0), pst.size(), [&](size_t i) {
                    const Vec3& xi = pst.xi[i];
                    const double lam[4] = {1.0 - xi.sum(), xi[0], xi[1], xi[2]};
                    BatchSeries sR, sQ[4];
                    if (rr) batchSeries(rr->sums, i, ci, uR, sR);
                    if (pr)
                        for (int q = 0; q < 4; q++) batchSeries(pr->sums, i, q, u4[q], sQ[q]);
                    auto eval = [&](const Eigen::Vector4d& c4, int a, double& H, Vec3& G) {
                        double Rv = 0, Qv[4] = {0, 0, 0, 0};
                        Vec3 Rg = Vec3::Zero(), Qg[4] = {Vec3::Zero(), Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};
                        if (rr) {
                            Rv = a < 0 ? sR.SU/sR.N : (sR.SU - sR.su[size_t(a)])/(sR.N - sR.n[size_t(a)]);
                            Rg = a < 0 ? Vec3(sR.G/double(B)) : Vec3((sR.G - sR.g[size_t(a)])/double(B - 1));
                        }
                        if (pr)
                            for (int q = 0; q < 4; q++) {
                                Qv[q] = a < 0 ? sQ[q].SU/sQ[q].N : (sQ[q].SU - sQ[q].su[size_t(a)])/(sQ[q].N - sQ[q].n[size_t(a)]);
                                Qg[q] = a < 0 ? Vec3(sQ[q].G/double(B)) : Vec3((sQ[q].G - sQ[q].g[size_t(a)])/double(B - 1));
                            }
                        H = o.cvVal[k][size_t(t)][i] + Rv;
                        G = T.JinvT*(o.cvGrad[k][size_t(t)][i] + Rg);
                        for (int q = 0; q < 4; q++) {
                            H += c4[q]*((known ? lam[q] : 0.0) + Qv[q]);
                            G += c4[q]*((known ? gN[size_t(q)] : Vec3::Zero()) + T.JinvT*Qg[q]);
                        }
                    };
                    eval(cc, -1, f.H[i], f.G[i]);
                    std::vector<double> Ha(static_cast<size_t>(B));
                    std::vector<Vec3> Ga(static_cast<size_t>(B));
                    double hm = 0;
                    Vec3 gm = Vec3::Zero();
                    for (int a = 0; a < B; a++) {
                        eval(cJ.empty() ? cc : cJ[size_t(a)], a, Ha[size_t(a)], Ga[size_t(a)]);
                        hm += Ha[size_t(a)]/double(B);
                        gm += Ga[size_t(a)]/double(B);
                    }
                    for (int a = 0; a < B; a++) {
                        f.vH[i] += (Ha[size_t(a)] - hm)*(Ha[size_t(a)] - hm)*double(B - 1)/double(B);
                        f.vG[i] += (Ga[size_t(a)] - gm).cwiseAbs2()*double(B - 1)/double(B);
                    }
                });
                o.F[k].push_back(std::move(f));
            }
        o.secField = secondsSince(t0);
    }
    o.secTotal = secondsSince(tStart);
    return o;
}

// H0 alone (no glue) of one variant and mode at a tet's points: C + R (cv) or the plain run, with per-point variance of H
// and per-batch gradient spread (validation of the control variates)
void stressH0(const StressOurs& o, int k, int t, const PointSet& ps, const TetFrame& T, size_t i, double& H, double& vH, Vec3& G,
              Vec3& vG)
{
    const int ci = o.comp[size_t(t)][size_t(k)];
    H = o.cvVal[size_t(k)][size_t(t)][i];
    vH = 0;
    G = o.cvGrad[size_t(k)][size_t(t)][i];
    vG = Vec3::Zero();
    if (ci < 0) {
        G = T.JinvT*G;
        return;
    }
    const MultiSums& s = o.runR[size_t(t)].sums;
    const std::vector<double> u = unitCoef(s.K, ci);
    double h, v;
    long long n;
    comboH(s, i, u, h, v, n);
    H += h;
    vH = v;
    const int Bn = s.nBatches;
    std::vector<Vec3> gb(static_cast<size_t>(Bn));
    Vec3 gm = Vec3::Zero();
    for (int b = 0; b < Bn; b++) gm += (gb[size_t(b)] = T.JinvT*comboBatchGrad(s, i, b, u))/double(Bn);
    for (const Vec3& x : gb) vG += (x - gm).cwiseAbs2()/(double(Bn)*double(Bn - 1));
    G = T.JinvT*G + gm;
    (void)ps;
}

// ---- the GT side of one case: its fields at the points and the probes, the regions, the energies ----
constexpr int kSR = 6;   // regions: all, near a front, near crack-face lines, near O's boundary, near interior faces, near the junction
const char* const kSRName[kSR] = {"all", "front", "lines", "outer", "iface", "junction"};
struct StressRef {
    std::vector<size_t> off;
    size_t nPts = 0;
    std::vector<double> w;
    std::vector<std::array<char, kSR>> reg;
    FieldAt gt[2], gtP[2];
    double gt2[2][kSR] = {}, Eref[2] = {}, ErefSe[2] = {}, gtRms[2] = {};
    ErrStats floor[2];
    int nModes = 1;
    bool exact = false, hasRun = false;
    RefRun ref;
    std::string source, check;
    double gtSeconds = 0;
};

std::vector<std::array<Vec3, 2>> walkFront(const StressCase& sc)
{
    std::vector<std::array<Vec3, 2>> f;
    for (const auto& t : sc.walk.tris()) {
        const Vec3* v[3] = {&t.a, &t.b, &t.c};
        for (int k = 0; k < 3; k++)
            if (t.edge[k] & CrackMesh::kEdgeFront) f.push_back({*v[k], *v[(k + 1)%3]});
    }
    return f;
}

void stressRefSetup(StressRef& r, const StressCase& sc, const Ring& R, const std::vector<StressTet>& stt, const std::vector<PointSet>& ps)
{
    const ConvexDomain O = octahedron();
    r.nModes = sc.ramp ? 2 : 1;
    r.off.assign(kRingNT, 0);
    for (int t = 1; t < kRingNT; t++) r.off[size_t(t)] = r.off[size_t(t - 1)] + ps[size_t(t - 1)].size();
    r.nPts = r.off.back() + ps.back().size();
    std::vector<Vec3> X;
    for (const PointSet& p : ps) {
        X.insert(X.end(), p.X.begin(), p.X.end());
        r.w.insert(r.w.end(), p.w.begin(), p.w.end());
    }
    const std::vector<std::array<Vec3, 2>> fr = walkFront(sc);
    std::vector<std::array<Vec3, 2>> lines;
    for (const StressTet& s : stt) lines.insert(lines.end(), s.faceLines.begin(), s.faceLines.end());
    r.reg.assign(r.nPts, {});
    tbb::parallel_for(size_t(0), r.nPts, [&](size_t i) {
        const Vec3& x = X[i];
        bool nf = false, nl = false, ni = false, nj = false;
        for (const auto& s : fr)
            if (pointSegmentDist(x, s[0], s[1]) < kNearRing) {
                nf = true;
                break;
            }
        for (const auto& s : lines)
            if (pointSegmentDist(x, s[0], s[1]) < kNearRing) {
                nl = true;
                break;
            }
        for (const auto& f : R.innerFaces)
            ni = ni || pointTriangleDist(x, R.V[size_t(f[0])], R.V[size_t(f[1])], R.V[size_t(f[2])]) < kNearRing;
        if (sc.junction)   // the junction half-line (E3): from B's front point jP along jW into B (4 > O's diameter)
            nj = pointSegmentDist(x, sc.jP, Vec3(sc.jP + 4.0*sc.jW)) < kNearRing;
        r.reg[i] = {1, char(nf), char(nl), char(distToDomainFaces(O, x) < kNearRing), char(ni), char(nj)};
    });
    const std::vector<Vec3> probes = ringProbes(R);
    for (int m = 0; m < r.nModes; m++) {
        if (r.hasRun) {
            r.gt[m] = refField(r.ref, m, 0, r.ref.nA);
            r.gtP[m] = refField(r.ref, m, r.ref.nA, r.ref.nF);
        } else {   // exact: the side indicator, no noise
            r.gt[m].resize(r.nPts);
            r.gtP[m].resize(probes.size());
            for (size_t i = 0; i < r.nPts; i++) r.gt[m].H[i] = sc.exact(X[i]);
            for (size_t i = 0; i < probes.size(); i++) r.gtP[m].H[i] = sc.exact(probes[i]);
        }
        double sh = 0, sw = 0;
        for (size_t g = 0; g < r.nPts; g++) {
            for (int k = 0; k < kSR; k++)
                if (r.reg[g][k]) r.gt2[m][k] += r.w[g]*r.gt[m].G[g].squaredNorm();
            r.floor[m].add(r.w[g], Vec3::Zero(), 0.0, Vec3::Zero(), r.gt[m].vH[g], 0.0, r.gt[m].vG[g], Vec3::Zero());
            sh += r.w[g]*r.gt[m].H[g]*r.gt[m].H[g];
            sw += r.w[g];
        }
        r.gtRms[m] = std::sqrt(sh/sw)/kTwoGMax;
        if (r.hasRun) {
            std::vector<double> wA(r.w);
            wA.resize(r.ref.X.size(), 0.0);
            solveEnergy(r.ref.sums, wA, m, r.Eref[m], r.ErefSe[m]);
        }
    }
}

// one measured row: a method at one setting, one face function, one g mode, one seed set
struct StressRow {
    std::string id, variant, method, label;
    int g = 0, rep = 0, nQ = 0, nR = 0, N = 0, B = 0;
    long long points = 0, walks = 0, moves = 0, capped = 0, escaped = 0;
    int tetsCut = 0, tetsP = -1, tetsR = -1;
    std::vector<double> times;
    double tMed = 0, tSetup = 0, tP = 0, tR = 0, tWalks = 0, tGram = 0, tGlue = 0, tField = 0;
    ErrStats e[kSR];
    double gt2[kSR] = {};
    double E = 0, Ese = 0, Eref = 0, ErefSe = 0;
    bool hasE = false;
    Vec7 vert = Vec7::Zero(), vref = Vec7::Zero(), cse = Vec7::Zero();
    bool hasC = false;
    double vErr = 0, gtRms = 0, maxFace = std::nan(""), maxCrack = std::nan("");
    bool exactGT = false;
};

void fillStressErrors(const StressRef& r, const StressCase& sc, int m, const std::vector<FieldAt>& F, const std::vector<PointSet>& ps,
                      StressRow& row)
{
    const FieldAt& g = r.gt[m];
    for (int t = 0; t < kRingNT; t++) {
        const FieldAt& f = F[size_t(t)];
        for (size_t i = 0; i < ps[size_t(t)].size(); i++) {
            const size_t q = r.off[size_t(t)] + i;
            for (int k = 0; k < kSR; k++)
                if (r.reg[q][k]) row.e[k].add(r.w[q], ps[size_t(t)].X[i], f.H[i] - g.H[q], f.G[i] - g.G[q], g.vH[q], f.vH[i], g.vG[q], f.vG[i]);
        }
    }
    for (int k = 0; k < kSR; k++) row.gt2[k] = r.gt2[m][k];
    for (int k = 0; k < kRingNV; k++) row.vref[k] = r.gtP[m].H[size_t(k)];
    row.vErr = 0;
    for (int k = sc.v0Undefined ? 1 : 0; k < kRingNV; k++) row.vErr = std::max(row.vErr, std::fabs(row.vert[k] - row.vref[k]));
    row.Eref = r.Eref[m];
    row.ErefSe = r.ErefSe[m];
    row.gtRms = r.gtRms[m];
    row.exactGT = !r.hasRun;
}

std::string stressCsvHeader()
{
    std::string h = "case,variant,g,method,label,seed_rep,nQ,nR,N,B,points,tets_cut,tets_P,tets_R,t_total_median,t_total_min,t_total_max,"
                    "n_repeats,t_setup,t_walks_P,t_walks_R,t_walks,t_gram,t_glue,t_field,walks,moves,mean_walk_len,capped,escaped,"
                    "H_rms,H_nc,H_max,H_se_fld,H_se_ref,H_rel_gt_rms,gt_rms,grad_rel,grad_rel_nc,grad_abs,grad_abs_nc,grad_se_ref_rel";
    for (int k = 1; k < kSR; k++) h += std::string(",H_rms_") + kSRName[k] + ",grad_rel_" + kSRName[k] + ",n_" + kSRName[k];
    h += ",E,E_se,E_ref,E_ref_se,E_relerr,E_relerr_se";
    for (int k = 0; k < kRingNV; k++) h += ",v" + std::to_string(k);
    for (int k = 0; k < kRingNV; k++) h += ",vref" + std::to_string(k);
    h += ",vert_maxerr,c_se_max,R_face_resid_max,R_crack_resid_max,exact_gt\n";
    return h;
}
std::string stressCsvRow(const StressRow& r)
{
    std::vector<std::string> c;
    auto add = [&](const std::string& s) { c.push_back(s); };
    auto num = [&](double v, const char* f = "%.6g") { c.push_back(csvNum(v, f)); };
    auto itg = [&](long long v) { c.push_back(std::to_string(v)); };
    add(r.id);
    add(r.variant);
    add(r.g ? "ramp" : "g1");
    add(r.method);
    add("\"" + r.label + "\"");
    itg(r.rep);
    itg(r.nQ);
    itg(r.nR);
    itg(r.N);
    itg(r.B);
    itg(r.points);
    itg(r.tetsCut);
    if (r.tetsP >= 0) itg(r.tetsP); else add("");
    if (r.tetsR >= 0) itg(r.tetsR); else add("");
    const double tMin = r.times.empty() ? 0.0 : *std::min_element(r.times.begin(), r.times.end());
    const double tMax = r.times.empty() ? 0.0 : *std::max_element(r.times.begin(), r.times.end());
    num(r.tMed, "%.5f");
    num(tMin, "%.5f");
    num(tMax, "%.5f");
    itg((long long)r.times.size());
    num(r.tSetup, "%.5f");
    num(r.tP, "%.5f");
    num(r.tR, "%.5f");
    num(r.tWalks + r.tP + r.tR, "%.5f");
    num(r.tGram, "%.5f");
    num(r.tGlue, "%.6f");
    num(r.tField, "%.5f");
    itg(r.walks);
    itg(r.moves);
    num(r.walks > r.escaped ? double(r.moves)/double(r.walks - r.escaped) : 0.0, "%.3f");
    itg(r.capped);
    itg(r.escaped);
    const ErrStats& e = r.e[0];
    num(e.rmsH());
    num(e.corrH());
    num(e.maxH/kTwoGMax);
    num(e.seFldH());
    num(e.seRefH());
    num(r.gtRms > 0 ? e.rmsH()/r.gtRms : std::nan(""));
    num(r.gtRms);
    num(gradRel(e, r.gt2[0]));
    num(gradRelNc(e, r.gt2[0]));
    num(e.rmsG());
    num(e.corrG());
    num(r.gt2[0] > 1e-12 ? std::sqrt(e.vRefG/r.gt2[0]) : std::nan(""));
    for (int k = 1; k < kSR; k++) {
        num(r.e[k].w > 0 ? r.e[k].rmsH() : std::nan(""));
        num(r.e[k].w > 0 ? gradRel(r.e[k], r.gt2[k]) : std::nan(""));
        itg(r.e[k].n);
    }
    auto rel = [](double E, double Er) { return std::fabs(Er) > 1e-12 ? (E - Er)/Er : std::nan(""); };
    auto relSe = [](double se, double seR, double Er) { return std::fabs(Er) > 1e-12 ? std::hypot(se, seR)/std::fabs(Er) : std::nan(""); };
    for (const double v : {r.E, r.Ese, r.Eref, r.ErefSe, rel(r.E, r.Eref), relSe(r.Ese, r.ErefSe, r.Eref)}) num(r.hasE ? v : std::nan(""));
    for (int k = 0; k < kRingNV; k++) num(r.vert[k], "%.5f");
    for (int k = 0; k < kRingNV; k++) num(r.vref[k], "%.5f");
    num(r.vErr, "%.5f");
    num(r.hasC ? r.cse.maxCoeff() : std::nan(""), "%.5f");
    num(r.maxFace, "%.3g");
    num(r.maxCrack, "%.3g");
    itg(r.exactGT);
    std::string s;
    for (size_t i = 0; i < c.size(); i++) s += (i ? "," : "") + c[i];
    return s + "\n";
}
void printStressRowHeader(Report& out)
{
    out("    %-30s %-14s %-4s %7s %8s %8s | %7s %7s %7s | %7s %7s | %7s %6s | %6s %6s\n", "method", "J", "g", "walks/M", "walks s",
        "total s", "H RMS%", "nc %", "max %", "grad%", "grad nc", "E err%", "+-", "vert", "c se");
}
void printStressRow(Report& out, const StressRow& r)
{
    const ErrStats& e = r.e[0];
    const double Erel = r.hasE && std::fabs(r.Eref) > 1e-12 ? 100*(r.E - r.Eref)/r.Eref : std::nan("");
    const double Ese = r.hasE && std::fabs(r.Eref) > 1e-12 ? 100*std::hypot(r.Ese, r.ErefSe)/std::fabs(r.Eref) : std::nan("");
    out("    %-30s %-14s %-4s %7.2f %8.3f %8.3f | %7.3f %7.3f %7.2f | %7.2f %7.2f | %+7.2f %6.2f | %6.4f %6.4f\n", r.label.c_str(),
        r.variant.c_str(), r.g ? "ramp" : "g=1", 1e-6*double(r.walks), r.tWalks + r.tP + r.tR, r.tMed, 100*e.rmsH(), 100*e.corrH(),
        100*e.maxH/kTwoGMax, 100*gradRel(e, r.gt2[0]), 100*gradRelNc(e, r.gt2[0]), Erel, Ese, r.vErr, r.hasC ? r.cse.maxCoeff() : 0.0);
}

// ---- checks of one case's geometry and face functions (printed and written to _geometry.txt) ----
// J at the crack-face lines: per line piece of every tet, 3 points, J at +-1e-7 along the face perpendicular to the line
// against side x g_j (points within 2e-3 of a front or of another crack skipped); the control variate's jump across every
// crack triangle (centroid +-1e-7 along the normal in xi) against 2 g_j (centroids within 2e-3 of a front or another crack
// skipped)
struct StressChecks {
    double maxJLine = 0, maxJrLine = 0, maxCvJump = 0;
    long long nJLine = 0, nCvJump = 0;
};
StressChecks stressChecks(const StressCase& sc, const Ring& R, const std::vector<StressTet>& stt)
{
    StressChecks ck;
    const std::vector<std::array<Vec3, 2>> fr = walkFront(sc);
    auto nearFront = [&](const Vec3& x) {
        for (const auto& s : fr)
            if (pointSegmentDist(x, s[0], s[1]) < 2e-3) return true;
        return false;
    };
    auto nearOther = [&](const Vec3& x, int j) {
        for (size_t k = 0; k < sc.cracks.size(); k++)
            if (int(k) != j && sc.cracks[k].dist && sc.cracks[k].dist(x) < 2e-3) return true;
        return false;
    };
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        const std::array<Vec3, 4> gN = gradN(T);
        std::vector<TetCv> cvs;
        for (const JVariant& v : sc.var) cvs.push_back(tetCv(T, v, sc));
        for (const auto& tr : stt[size_t(t)].mesh.tris()) {
            const Vec3 nX = (T.JinvT*tr.n).normalized();   // the triangle's normal in X (J^-T maps normals)
            const Vec3 pts[3] = {T.toX(tr.a), T.toX(tr.b), T.toX(tr.c)};
            for (int k = 0; k < 3; k++) {
                if (!(tr.edge[k] & (CrackMesh::kEdgeDomain | CrackMesh::kEdgeTet))) continue;
                const Vec3 &A = pts[k], &B = pts[(k + 1)%3];
                // the tet face holding the edge: the barycentric coordinate that vanishes at both ends
                const Vec3 xa = T.toXi(A), xb = T.toXi(B);
                const double la[4] = {1 - xa.sum(), xa[0], xa[1], xa[2]}, lb[4] = {1 - xb.sum(), xb[0], xb[1], xb[2]};
                int o = -1;
                for (int f = 0; f < 4; f++)
                    if (std::fabs(la[f]) < 1e-9 && std::fabs(lb[f]) < 1e-9) o = f;
                if (o < 0) continue;
                const Vec3 e = (B - A).normalized(), mf = gN[size_t(o)].normalized();
                Vec3 d = mf.cross(e);
                if (d.norm() < 1e-9) continue;
                d.normalize();
                for (double s : {0.25, 0.5, 0.75}) {
                    const Vec3 P = A + s*(B - A);
                    if (nearFront(P) || nearOther(P, tr.crack)) continue;
                    for (int sg : {-1, 1}) {
                        const Vec3 Q = P + sg*1e-7*d;
                        const int side = nX.dot(Q - P) >= 0 ? 1 : -1;
                        for (const JVariant& v : sc.var) {
                            ck.maxJLine = std::max(ck.maxJLine, std::fabs(v.J(Q) - stressTarget(sc, tr.crack, side, P, 0)));
                            if (v.Jr) ck.maxJrLine = std::max(ck.maxJrLine, std::fabs(v.Jr(Q) - stressTarget(sc, tr.crack, side, P, 1)));
                        }
                        ck.nJLine++;
                    }
                }
            }
            // the control variates' jump across the crack, in xi
            const Vec3 q = (tr.a + tr.b + tr.c)/3.0, Xq = T.toX(q);
            if (nearFront(Xq) || nearOther(Xq, tr.crack)) continue;
            for (const TetCv& cv : cvs) {
                const double jump = cv.value(q + 1e-7*tr.n, 0) - cv.value(q - 1e-7*tr.n, 0);
                ck.maxCvJump = std::max(ck.maxCvJump, std::fabs(jump - 2*sc.cracks[size_t(tr.crack)].g));
                ck.nCvJump++;
            }
        }
    }
    return ck;
}

// ---- validation pieces ----
// two solves of the same field at the same points: per point z of H (both g modes) and pooled z of grad H per component
struct PairZ {
    ZStats H[2], G[2][3];
};
PairZ compareSolves(const RefRun& a, const RefRun& b, int nModes)
{
    PairZ z;
    for (int m = 0; m < nModes; m++) {
        const FieldAt fa = refField(a, m, 0, a.nA), fb = refField(b, m, 0, b.nA);
        for (size_t i = 0; i < a.nA; i++) {
            z.H[m].add(i, fa.H[i] - fb.H[i], fa.vH[i] + fb.vH[i]);
            for (int c = 0; c < 3; c++) z.G[m][c].add(i, fa.G[i][c] - fb.G[i][c], fa.vG[i][c] + fb.vG[i][c]);
        }
    }
    return z;
}
std::string pairZLine(const PairZ& z, int m)
{
    const ZStats& h = z.H[m];
    char buf[512];
    std::snprintf(buf, sizeof buf, "H: pooled z %+.2f, |z| < 3 at %.2f%% (%lld points; %lld with zero variance in both, max |dH| %.1e); "
                  "grad pooled z %+.2f / %+.2f / %+.2f -> %s", h.pooled(), 100*h.frac(), h.n, h.zeroVar, h.maxAbsDiffZeroVar,
                  z.G[m][0].pooled(), z.G[m][1].pooled(), z.G[m][2].pooled(),
                  (h.n == 0 || h.passPooled()) && h.maxAbsDiffZeroVar <= 1e-12 && std::fabs(z.G[m][0].pooled()) < 3 &&
                          std::fabs(z.G[m][1].pooled()) < 3 && std::fabs(z.G[m][2].pooled()) < 3
                      ? "agree"
                      : "DISAGREE");
    return buf;
}

// Zombie (float, fcpw) on O with the walk crack absorbing (data +-g_j by the nearest triangle's crack and side), g = 1:
// H at every point (N), grad H at the points farther than kRefZGradBand from O's faces (S solves of N/S walks)
constexpr int kSZR = 5;
const char* const kSZRName[kSZR] = {"near a crack (0.1)", "near a front (0.1)", "near interior faces (0.1)", "near the junction (0.1)",
                                    "vertex probes"};
struct StressZombie {
    ZStats H, G[3], Hr[kSZR];
    long long nGradPts = 0, nNaN = 0;
    SolveStats st;
    bool agree() const
    {
        bool g = true;
        for (int a = 0; a < 3; a++) g = g && G[a].n > 0 && std::fabs(G[a].pooled()) < 3.0;
        return H.passPooled() && g;
    }
};
StressZombie zombieStress(const StressCase& sc, const Ring& R, const RefRun& r, const StressRef& sr, int N, int S, uint64_t seedTag)
{
    const GateParams gp = gateParams(false);
    StressZombie z;
    const ConvexDomain O = octahedron();
    WalkParams wp;
    wp.absorbShellRel = kGateShellRel*kSqrt2/O.L;
    wp.reflectShellRel = 0.1*wp.absorbShellRel;
    std::vector<Vec3> V;
    std::vector<std::array<int, 3>> F;
    for (const auto& t : sc.walk.tris()) {
        const int b = int(V.size());
        V.insert(V.end(), {t.a, t.b, t.c});
        F.push_back({b, b + 1, b + 2});
    }
    const CrackMesh* walk = &sc.walk;
    const StressCase* scp = &sc;
    CutSolver solver(O, V, F, [walk, scp](const Vec3& x, bool normalAligned) -> float {
        const CrackMesh::Hit h = walk->closest(x);
        const double g = scp->cracks[size_t(walk->tri(h.tri).crack)].g;
        return float(normalAligned ? g : -g);
    }, wp);
    std::vector<size_t> gi;
    std::vector<Vec3> gX;
    for (size_t i = 0; i < r.nA; i++)
        if (distToDomainFaces(O, r.X[i]) > kRefZGradBand) {
            gi.push_back(i);
            gX.push_back(r.X[i]);
        }
    z.nGradPts = (long long)gi.size();
    std::vector<Estimate> e;
    z.st.add(solver.solve(r.X, nullptr, N, mixSeed(gp.seed, seedTag), e));
    const std::vector<double> coef = unitCoef(2, 0);
    for (size_t i = 0; i < r.X.size(); i++) {
        double H, var;
        long long n;
        comboH(r.sums, i, coef, H, var, n);
        const double diff = H - e[i].H, v = var + double(e[i].Hse)*e[i].Hse;
        z.H.add(i, diff, v);
        bool inR[kSZR] = {false, false, false, false, i >= r.nA};
        if (i < r.nA) {
            inR[0] = sc.walk.distance(r.X[i]) < kNearRing;
            inR[1] = sr.reg[i][1] != 0;
            inR[2] = sr.reg[i][4] != 0;
            inR[3] = sr.reg[i][5] != 0;
        }
        for (int k = 0; k < kSZR; k++)
            if (inR[k]) z.Hr[k].add(i, diff, v);
    }
    const std::vector<char> gf(gX.size(), 1);
    std::vector<std::vector<Estimate>> eg(static_cast<size_t>(S));
    for (int s = 0; s < S; s++) z.st.add(solver.solve(gX, &gf, N/S, mixSeed(gp.seed, seedTag + 10 + s), eg[size_t(s)]));
    const int Bref = r.sums.nBatches;
    for (size_t j = 0; j < gi.size(); j++) {
        Vec3 zm = Vec3::Zero(), zv = Vec3::Zero();
        bool finite = true;
        for (int s = 0; s < S; s++) {
            finite = finite && eg[size_t(s)][j].grad.allFinite();
            zm += eg[size_t(s)][j].grad.cast<double>()/double(S);
        }
        if (!finite) {
            z.nNaN++;
            continue;
        }
        for (int s = 0; s < S; s++) zv += (eg[size_t(s)][j].grad.cast<double>() - zm).cwiseAbs2()/(double(S)*double(S - 1));
        Vec3 om = Vec3::Zero(), ov = Vec3::Zero();
        std::vector<Vec3> g(static_cast<size_t>(Bref));
        for (int b = 0; b < Bref; b++) om += (g[size_t(b)] = comboBatchGrad(r.sums, gi[j], b, coef))/double(Bref);
        for (const Vec3& x : g) ov += (x - om).cwiseAbs2()/(double(Bref)*double(Bref - 1));
        for (int a = 0; a < 3; a++) z.G[a].add(gi[j], om[a] - zm[a], ov[a] + zv[a]);
    }
    return z;
}

// ---- the driver of one case ----
struct StressSummary {
    std::string id, plain, gtLine, check, geometry;
    int priority = 1;
    bool ran = false;
    std::string failure;
};

std::string stressGtName(const std::string& id, int pts, int N, double epsRel = 0.01)
{
    return "ref_stress_" + id + "_p" + std::to_string(pts) + (N != 25600 ? "_N" + std::to_string(N) : "") +
           (epsRel != 0.01 ? "_eps" + fmtd("%g", epsRel) : "") + ".bin";
}

void runStressCase(Report& out, const RingStressParams& p, StressCase& sc, std::vector<StressRow>& rows, StressSummary& sum,
                   std::string& geoTxt)
{
    using clk = std::chrono::steady_clock;
    const auto tCase = clk::now();
    const GateParams gp = gateParams(false);
    const int k = sc.index;
    sum.id = sc.id;
    sum.plain = sc.plain;
    sum.priority = sc.priority;
    out("\n=== %s (priority %d) ===\n", sc.id.c_str(), sc.priority);
    const Ring R = stressRing(sc);
    // geometry
    auto tg = clk::now();
    std::vector<StressTet> stt;
    for (const TetFrame& T : R.T) stt.push_back(stressTet(T, sc.walk));
    double aSum = 0;
    int nCut = 0;
    std::string geo = sc.id + ": " + sc.plain + "\n  " + sc.desc + "\n  " + R.name + "; walk crack: " + std::to_string(sc.walk.size()) +
                      " triangles in O, area " + fmtd("%.6f", sc.walk.area()) + "; cracks:";
    for (const StressCrack& c : sc.cracks) geo += " " + c.name + " (g = " + fmtd("%g", c.g) + ")";
    geo += "\n  per tet (cut, triangles in xi, area in X, crack-face line pieces, front pieces):";
    for (int t = 0; t < kRingNT; t++) {
        const StressTet& s = stt[size_t(t)];
        aSum += s.areaX;
        nCut += s.cut;
        geo += std::string(" ") + kRingTetName[t] + " " + (s.cut ? "cut" : "-") + " " + std::to_string(s.mesh.size()) + " " +
               fmtd("%.4f", s.areaX) + " " + std::to_string(s.faceLines.size()) + " " + std::to_string(s.front.size()) + ";";
    }
    const double aRel = std::fabs(aSum - sc.walk.area())/std::max(1e-300, sc.walk.area());
    geo += "\n  sum of the tets' areas " + fmtd("%.12f", aSum) + " vs the walk crack's " + fmtd("%.12f", sc.walk.area()) + " (rel. diff " +
           fmtd("%.1e", aRel) + ")\n";
    const StressChecks ck = stressChecks(sc, R, stt);
    geo += "  J at the crack-face lines (+-1e-7 along the face): max |J - side g| " + fmtd("%.2e", ck.maxJLine) +
           (sc.ramp ? ", ramp " + fmtd("%.2e", ck.maxJrLine) : std::string()) + " over " + std::to_string(ck.nJLine) +
           " samples; control variates' jump across the crack (xi, +-1e-7): max |jump - 2 g| " + fmtd("%.2e", ck.maxCvJump) + " over " +
           std::to_string(ck.nCvJump) + " samples\n";
    for (const JVariant& v : sc.var) geo += "  J '" + v.name + "': " + v.desc + "\n";
    out("%s  geometry built and checked in %.2f s\n", geo.c_str(), secondsSince(tg));
    // the points and the GT
    std::vector<PointSet> ps;
    StressRef sr;
    const int target = p.quick ? 512 : p.points;
    const int Ngt = p.quick ? 1024 : p.gtWalks, Bgt = p.quick ? 16 : p.gtBatches;
    const int maxLenGT = 65536;
    auto buildPoints = [&](int tgt) {
        ps.clear();
        for (int t = 0; t < kRingNT; t++) ps.push_back(stressPoints(R.T[size_t(t)], stt[size_t(t)], sc.walk, tgt, p.epsRel*kSqrt2));
    };
    auto tGT = clk::now();
    if (sc.check4 >= 0) {
        PerfSetup s;
        perfSetup(s, 2, sc.check4, gp.nRingPoints);
        ps = s.ps;
        if (!loadRingGT(s, sr.ref)) {
            out("  check 4's cached GT %s is missing or its key differs: case skipped (nothing is solved or written)\n", sr.ref.file.c_str());
            sum.failure = "no cached GT";
            return;
        }
        sr.hasRun = true;
        sr.source = "check 4's cached GT " + sr.ref.file + " (read-only; N = 25600, 4096 points per tet)";
        sr.gtSeconds = sr.ref.stats.seconds;
        out("  GT: %s; original run: %s\n", sr.source.c_str(), sr.ref.stats.summary().c_str());
    } else if (sc.exact) {
        buildPoints(target);
        sr.source = "exact (+-1, the side indicator)";
    } else {
        int tgt = target;
        buildPoints(tgt);
        double extrap = 0;
        RefRun r = stressGT(out, sc, R, ps, Ngt, Bgt, 5000000 + 1000*uint64_t(k), maxLenGT, stressGtName(sc.id, tgt, Ngt, p.epsRel), !p.quick && p.cache, &extrap,
                            true);
        bool lowered = false;
        if (!r.cached && extrap > p.pilotLimit && tgt > 1024) {
            out("  the GT pilot puts the full GT at %.0f s > %.0f s: points per tet lowered to 1024 (the walks are not lowered)\n", extrap,
                p.pilotLimit);
            tgt = 1024;
            buildPoints(tgt);
            lowered = true;
        }
        if (!r.cached)   // (the cache at the lowered density is looked up again; the pilot is not repeated at the same points)
            r = stressGT(out, sc, R, ps, Ngt, Bgt, 5000000 + 1000*uint64_t(k), maxLenGT, stressGtName(sc.id, tgt, Ngt, p.epsRel), !p.quick && p.cache, nullptr,
                         false, !lowered);
        sr.ref = std::move(r);
        sr.hasRun = true;
        sr.source = std::string(sr.ref.cached ? "cached " : "solved ") + sr.ref.file + " (N = " + std::to_string(Ngt) + ", " +
                    std::to_string(tgt) + " points per tet asked)";
        sr.gtSeconds = sr.ref.stats.seconds;
    }
    size_t nPts = 0;
    for (const PointSet& q : ps) nPts += q.size();
    stressRefSetup(sr, sc, R, stt, ps);
    sum.gtLine = sr.source + "; GT solve " + fmtd("%.1f s", sr.gtSeconds) +
                 (sr.hasRun ? ", " + std::to_string(sr.ref.stats.walksCapped) + " capped, " + std::to_string(sr.ref.stats.walksEscaped) +
                                  " escaped walks, mean length " + fmtd("%.1f", sr.ref.stats.meanWalkLength())
                            : std::string());
    out("  %zu points in the 8 tets + 7 probes; %d tets cut; GT ready in %.1f s\n", nPts, nCut, secondsSince(tGT));
    for (int m = 0; m < sr.nModes; m++)
        out("  GT noise floor (%s): H s.e. RMS %.4f%% (max %.4f%%) of 2g; grad s.e. RMS / grad RMS %.2f%%; GT RMS %.4f; energy %.5f +- %.5f\n",
            m ? "ramp" : "g = 1", 100*sr.floor[m].seRefH(), 100*sr.floor[m].maxSeRefH/kTwoGMax,
            100*std::sqrt(sr.floor[m].vRefG/std::max(1e-300, sr.gt2[m][0])), sr.gtRms[m], sr.Eref[m], sr.ErefSe[m]);
    // validation
    std::string check;
    if (p.validate) {
        if (sc.check4 >= 0) {   // the new mesh GT code reproduces check 4's cached GT
            const int N = p.quick ? 256 : 2560, B = p.quick ? 16 : 64;
            const RefRun nr = stressSolveO(sc, R, ps, N, B, mixSeed(gp.seed, 9500000 + 1000*uint64_t(k)), 4096);
            const PairZ z = compareSolves(nr, sr.ref, sr.nModes);
            out("  flat reproduction: the mesh GT code (N = %d, new seeds, %.1f s) vs check 4's cached GT:\n", N, nr.stats.seconds);
            for (int m = 0; m < sr.nModes; m++) out("    %s %s\n", m ? "ramp " : "g = 1", pairZLine(z, m).c_str());
            check += "mesh GT (N = " + std::to_string(N) + ") vs cached: " + pairZLine(z, 0) + "; ";
        }
        if (sc.exact || sc.check4 == 2) {   // the GT solver gives exactly +-1
            const int N = p.quick ? 64 : 256, B = 16;
            const RefRun nr = stressSolveO(sc, R, ps, N, B, mixSeed(gp.seed, 9600000 + 1000*uint64_t(k)), maxLenGT);
            long long bad = 0;
            double maxd = 0;
            for (size_t i = 0; i < nr.X.size(); i++) {
                double H, var;
                long long n;
                comboH(nr.sums, i, unitCoef(2, 0), H, var, n);
                const double d = std::fabs(H - sc.exact(nr.X[i]));
                maxd = std::max(maxd, d);
                bad += d > 1e-12;
            }
            char buf[300];
            std::snprintf(buf, sizeof buf, "GT solver at N = %d: %lld of %zu points differ from the exact +-1 (max |dH| %.1e), %lld capped, "
                          "%lld escaped walks", N, bad, nr.X.size(), maxd, nr.stats.walksCapped, nr.stats.walksEscaped);
            out("  exact check: %s\n", buf);
            check += std::string(buf) + "; ";
        }
        const bool zz = p.zombie != "none" && ("," + p.zombie + ",").find("," + sc.id + ",") != std::string::npos;
        if (zz && sr.hasRun) {
            const int N = p.quick ? 128 : 1024, S = p.quick ? 8 : 16;
            const auto t0 = clk::now();
            const StressZombie z = zombieStress(sc, R, sr.ref, sr, N, S, 8000000 + 1000*uint64_t(k));
            out("  Zombie cross-check (float; %s; %.1f s): H at %lld points, N = %d: |z| < 3 at %.2f%%, pooled z %+.2f; grad at %lld points "
                "(%d x %d walks, %lld not finite): pooled z %+.2f / %+.2f / %+.2f -> %s\n", z.st.summary().c_str(), secondsSince(t0),
                z.H.n, N, 100*z.H.frac(), z.H.pooled(), z.nGradPts, S, N/S, z.nNaN, z.G[0].pooled(), z.G[1].pooled(), z.G[2].pooled(),
                z.agree() ? "agree" : "DISAGREE");
            for (int r = 0; r < kSZR; r++)
                out("    %-28s %6lld points: pooled z %+.2f, mean z^2 %.3f (%.3f), RMS discrepancy beyond noise %.1e (bound %.1e)\n",
                    kSZRName[r], z.Hr[r].n, z.Hr[r].pooled(), z.Hr[r].meanZ2(), z.Hr[r].meanZ2Se(), z.Hr[r].ncRmsDiff(kTwoGMax),
                    z.Hr[r].ncRmsDiffBound(kTwoGMax));
            check += "Zombie " + std::string(z.agree() ? "agrees" : "DISAGREES") + " (H pooled z " + fmtd("%+.2f", z.H.pooled()) + ", |z| < 3 at " +
                     fmtd("%.2f%%", 100*z.H.frac()) + "; grad pooled z " + fmtd("%+.2f", z.G[0].pooled()) + "/" + fmtd("%+.2f", z.G[1].pooled()) +
                     "/" + fmtd("%+.2f", z.G[2].pooled()) + "); ";
        }
        if (sc.id == "D1" && sr.hasRun && !p.quick) {   // the GT's own consistency: an independent solve at N = 2560
            const RefRun nr = stressSolveO(sc, R, ps, 2560, 64, mixSeed(gp.seed, 9700000 + 1000*uint64_t(k)), maxLenGT);
            const PairZ z = compareSolves(nr, sr.ref, 1);
            out("  GT re-run at N = 2560 (independent seeds, %.1f s) vs the GT: %s\n", nr.stats.seconds, pairZLine(z, 0).c_str());
            check += "GT re-run N = 2560: " + pairZLine(z, 0) + "; ";
        }
    }
    // ours
    const std::vector<std::vector<int>> gOurs = parseGrid(p.quick && p.ours == "32/16/8,32/64/8,1024/256/16" ? "8/8/4" : p.ours, "32/16/8");
    const std::vector<std::vector<int>> gDir = parseGrid(p.quick && p.direct == "16,32,64,128,256,512" ? "16,64" : p.direct, "16");
    printStressRowHeader(out);
    auto baseRow = [&](const std::string& method, const std::string& variant, int g, int rep) {
        StressRow r;
        r.id = sc.id;
        r.method = method;
        r.variant = variant;
        r.g = g;
        r.rep = rep;
        r.points = (long long)nPts;
        r.tetsCut = nCut;
        return r;
    };
    auto oursRows = [&](const StressOurs& o, const StressOursParams& op, const std::vector<double>& times, int rep, const std::string& method,
                        const std::string& lab) {
        for (size_t kk = 0; kk < o.vms.size(); kk++) {
            const int v = o.vms[kk].first, m = o.vms[kk].second;
            StressRow r = baseRow(method, sc.var[size_t(v)].name, m, rep);
            r.label = lab;
            r.nQ = op.nQ;
            r.nR = op.nR;
            r.B = op.B;
            r.times = times;
            r.tMed = medianOf(times);
            r.tSetup = o.secSetup;
            r.tP = o.secP;
            r.tR = o.secR;
            r.tGram = o.secGram;
            r.tGlue = o.secGlue;
            r.tField = o.secField;
            r.walks = o.walksP + o.walksR;
            r.moves = o.stepsP + o.stepsR;
            r.capped = o.capped;
            r.escaped = o.escaped;
            r.tetsP = o.nPtets;
            r.tetsR = o.nRtets;
            r.maxFace = o.maxFace[kk];
            r.maxCrack = o.maxCrack[kk];
            const JVariant& jv = sc.var[size_t(v)];
            for (int q = 0; q < kRingNV; q++) r.vert[q] = (m ? jv.Jr(R.V[size_t(q)]) : jv.J(R.V[size_t(q)])) + o.gl[kk].c[q];
            r.cse = o.gl[kk].se;
            r.hasC = true;
            fillStressErrors(sr, sc, m, o.F[kk], ps, r);
            glueEnergy(R, o.G[kk], 0, o.gl[kk], r.E, r.Ese);
            r.hasE = true;
            printStressRow(out, r);
            rows.push_back(r);
        }
    };
    for (size_t si = 0; si < gOurs.size(); si++) {
        StressOursParams op;
        op.nQ = gOurs[si][0];
        op.nR = gOurs[si].size() > 1 ? gOurs[si][1] : op.nQ;
        op.B = gOurs[si].size() > 2 ? gOurs[si][2] : defaultBatches(op.nQ);
        const int nSeeds = si == 0 ? std::max(1, p.seeds) : 1;
        for (int rep = 0; rep < nSeeds; rep++) {
            const uint64_t tagBase = 6000000 + 1000*uint64_t(k) + 100*uint64_t(rep) + 20000*uint64_t(si);
            StressOurs o;
            const std::vector<double> times = timedRuns(rep == 0 ? std::max(1, p.repeats) : 1, [&] {
                o = runStressOurs(sc, R, stt, ps, op, gp.seed, tagBase, true);
                return o.secTotal;
            });
            const std::string lab = "ours Q" + std::to_string(op.nQ) + "/R" + std::to_string(op.nR) + "/B" + std::to_string(op.B) +
                                    (rep ? " s" + std::to_string(rep) : "");
            oursRows(o, op, times, rep, "ours", lab);
        }
    }
    // the slow path for the control variates' check: plain per-tet solves (no control variate) at N = 2560, and H0 of the
    // control-variate path at Q32/R256, per point
    if (p.validate && (sc.id == "C2" || sc.id == "D1" || sc.id == "D3" || sc.id == "E3")) {
        StressOursParams pp;
        pp.nQ = pp.nR = p.quick ? 256 : 2560;
        pp.B = p.quick ? 16 : 64;
        pp.cv = false;
        StressOurs pl;
        const std::vector<double> tp = timedRuns(1, [&] {
            pl = runStressOurs(sc, R, stt, ps, pp, gp.seed, 9000000 + 1000*uint64_t(k), true);
            return pl.secTotal;
        });
        oursRows(pl, pp, tp, 0, "ours-plain", "plain Q" + std::to_string(pp.nQ) + "/R" + std::to_string(pp.nR) + "/B" + std::to_string(pp.B));
        StressOursParams cp;
        cp.nQ = 32;
        cp.nR = p.quick ? 64 : 256;
        cp.B = 16;
        const StressOurs cv = runStressOurs(sc, R, stt, ps, cp, gp.seed, 9100000 + 1000*uint64_t(k), false);
        ZStats zv, zg[3];
        for (int t = 0; t < kRingNT; t++)
            for (size_t i = 0; i < ps[size_t(t)].size(); i++) {
                if (dirichletDist(allDirichlet(), ps[size_t(t)].xi[i]) < kCompareBand) continue;
                double h1, v1, h2, v2;
                Vec3 g1, gv1, g2, gv2;
                stressH0(pl, 0, t, ps[size_t(t)], R.T[size_t(t)], i, h1, v1, g1, gv1);
                stressH0(cv, 0, t, ps[size_t(t)], R.T[size_t(t)], i, h2, v2, g2, gv2);
                zv.add(i, h1 - h2, v1 + v2);
                for (int a = 0; a < 3; a++) zg[a].add(i, g1[a] - g2[a], gv1[a] + gv2[a]);
            }
        char buf[400];
        std::snprintf(buf, sizeof buf, "H0 (g = 1, J '%s'): plain N = %d vs control variate Q32/R%d: value pooled z %+.2f (|z| < 3 at %.2f%%, "
                      "mean z^2 %.3f), grad pooled z %+.2f / %+.2f / %+.2f -> %s", sc.var[0].name.c_str(), pp.nR, cp.nR, zv.pooled(),
                      100*zv.frac(), zv.meanZ2(), zg[0].pooled(), zg[1].pooled(), zg[2].pooled(),
                      zv.passPooled() && std::fabs(zg[0].pooled()) < 3 && std::fabs(zg[1].pooled()) < 3 && std::fabs(zg[2].pooled()) < 3
                          ? "agree"
                          : "DISAGREE");
        out("  control-variate check: %s\n", buf);
        check += std::string(buf) + "; ";
    }
    // the direct ring solve (the GT's solver at small N)
    for (size_t di = 0; di < gDir.size(); di++) {
        const int N = gDir[di][0], B = gDir[di].size() > 1 ? gDir[di][1] : defaultBatches(N);
        RefRun d;
        const std::vector<double> times = timedRuns(1, [&] {
            const auto t0 = clk::now();
            d = stressSolveO(sc, R, ps, N, B, mixSeed(gp.seed, 7000000 + 1000*uint64_t(k) + 10*uint64_t(di)), maxLenGT);
            return secondsSince(t0);
        });
        std::vector<double> tAll = times;
        if (times[0] < 2.0)
            for (int r = 1; r < 3; r++) {
                const auto t0 = clk::now();
                d = stressSolveO(sc, R, ps, N, B, mixSeed(gp.seed, 7000000 + 1000*uint64_t(k) + 10*uint64_t(di)), maxLenGT);
                tAll.push_back(secondsSince(t0));
            }
        const auto tf = clk::now();
        std::vector<FieldAt> F[2];
        FieldAt P[2];
        for (int m = 0; m < sr.nModes; m++) {
            size_t off = 0;
            for (const PointSet& q : ps) {
                F[m].push_back(refField(d, m, off, q.size()));
                off += q.size();
            }
            P[m] = refField(d, m, d.nA, d.nF);
        }
        const double secField = secondsSince(tf);
        for (int m = 0; m < sr.nModes; m++) {
            StressRow r = baseRow("direct", "-", m, 0);
            r.label = "direct N" + std::to_string(N) + "/B" + std::to_string(B);
            r.N = N;
            r.B = B;
            r.times = tAll;
            for (double& x : r.times) x += secField;
            r.tMed = medianOf(r.times);
            r.tWalks = d.stats.seconds;
            r.tField = secField;
            r.walks = d.stats.walksRequested;
            r.moves = d.stats.totalSteps;
            r.capped = d.stats.walksCapped;
            r.escaped = d.stats.walksEscaped;
            for (int q = 0; q < kRingNV; q++) r.vert[q] = P[m].H[size_t(q)];
            fillStressErrors(sr, sc, m, F[m], ps, r);
            std::vector<double> wA(sr.w);
            wA.resize(d.X.size(), 0.0);
            solveEnergy(d.sums, wA, m, r.E, r.Ese);
            r.hasE = true;
            printStressRow(out, r);
            rows.push_back(r);
        }
    }
    sum.check = check.empty() ? "-" : check;
    sum.geometry = geo;
    sum.ran = true;
    geoTxt += geo + "  GT: " + sum.gtLine + "\n  checks: " + sum.check + "\n\n";
    out("  case %s done in %.1f s\n", sc.id.c_str(), secondsSince(tCase));
}

// ---- the per-case summary (the main comparison numbers; the report is written from these) ----
struct OursNums {
    int n = 0;
    double rms = std::nan(""), rmsSd = 0, maxp = std::nan(""), vert = std::nan(""), grad = std::nan(""), gradNc = std::nan(""),
           gradAbs = std::nan(""), E = std::nan(""), t = std::nan(""), nc = std::nan(""), relGt = std::nan("");
};
OursNums oursNums(const std::vector<StressRow>& rows, const std::string& id, const std::string& method, const std::string& variant, int g,
                  int nQ, int nR, int B)
{
    std::vector<const StressRow*> v;
    for (const StressRow& r : rows)
        if (r.id == id && r.method == method && r.variant == variant && r.g == g && r.nQ == nQ && r.nR == nR && r.B == B) v.push_back(&r);
    OursNums o;
    o.n = int(v.size());
    if (v.empty()) return o;
    auto med = [&](auto f) {
        std::vector<double> x;
        for (const StressRow* r : v) x.push_back(f(*r));
        return medianOf(x);
    };
    o.rms = med([](const StressRow& r) { return r.e[0].rmsH(); });
    double m = 0, s = 0;
    for (const StressRow* r : v) m += r->e[0].rmsH()/double(v.size());
    for (const StressRow* r : v) s += (r->e[0].rmsH() - m)*(r->e[0].rmsH() - m)/std::max(1.0, double(v.size()) - 1);
    o.rmsSd = std::sqrt(s);
    o.nc = med([](const StressRow& r) { return r.e[0].corrH(); });
    o.maxp = med([](const StressRow& r) { return r.e[0].maxH/kTwoGMax; });
    o.vert = med([](const StressRow& r) { return r.vErr; });
    o.grad = med([](const StressRow& r) { return gradRel(r.e[0], r.gt2[0]); });
    o.gradNc = med([](const StressRow& r) { return gradRelNc(r.e[0], r.gt2[0]); });
    o.gradAbs = med([](const StressRow& r) { return r.e[0].rmsG(); });
    o.E = med([](const StressRow& r) { return r.hasE && std::fabs(r.Eref) > 1e-12 ? (r.E - r.Eref)/r.Eref : std::nan(""); });
    o.t = v[0]->tMed;   // the first seed set's timing repeats
    o.relGt = med([](const StressRow& r) { return r.gtRms > 0 ? r.e[0].rmsH()/r.gtRms : std::nan(""); });
    return o;
}
// the direct solve's H RMS at time T (log-log interpolation between its measured settings; below the fastest one: that point
// scaled by sqrt(t/T), RMS ~ 1/sqrt(N); flagged)
double directAt(const std::vector<StressRow>& rows, const std::string& id, int g, double T, std::string& how)
{
    std::vector<std::pair<double, double>> p;
    for (const StressRow& r : rows)
        if (r.id == id && r.method == "direct" && r.g == g) p.push_back({r.tMed, r.e[0].rmsH()});
    std::sort(p.begin(), p.end());
    how = "-";
    if (p.empty() || !(T > 0)) return std::nan("");
    if (T <= p[0].first) {
        how = "extrapolated";
        return p[0].second*std::sqrt(p[0].first/T);
    }
    for (size_t i = 1; i < p.size(); i++)
        if (T <= p[i].first) {
            how = "interpolated";
            const double a = std::log(p[i - 1].first), b = std::log(p[i].first);
            return std::exp(std::log(p[i - 1].second) + (std::log(T) - a)*(std::log(p[i].second) - std::log(p[i - 1].second))/(b - a));
        }
    how = "beyond the largest N (its RMS)";
    return p.back().second;
}
std::string stressVerdict(const OursNums& o, bool exactIndicator)
{
    if (o.n == 0 || !std::isfinite(o.rms)) return "fails (no result)";
    if (exactIndicator) return o.maxp <= 1e-9 ? "works (exact to round-off)" : "fails (not exact: " + fmtd("%.1e", o.maxp) + ")";
    if (o.rms <= 0.02 && o.maxp <= 0.20) return "works";
    if (o.rms <= 0.05) return "degrades";
    return "fails";
}

std::string stressSummaryCsvHeader()
{
    return "case,plain_name,priority,variant,g,ran,ours_setting,n_seeds,ours_H_rms,ours_H_rms_sd,ours_H_nc,ours_H_max,ours_vert_maxerr,"
           "ours_grad_rel,ours_grad_rel_nc,ours_grad_abs,ours_E_relerr,ours_time,ours_H_rel_gt_rms,s1_H_rms,s2_H_rms,s2_H_max,"
           "s2_vert_maxerr,direct_H_rms_at_ours_time,direct_how,gt,checks,verdict,failure\n";
}

// ---- self-test ----
bool stressSelftest(Report& out)
{
    bool ok = true;
    std::mt19937_64 rng(20261103);
    std::uniform_real_distribution<double> U(-1.2, 1.2), U01(0.0, 1.0);
    // 1. the tree against brute force, bit for bit
    out("\n=== self-test 1: CrackMesh::closest (tree) against brute force, bit for bit (d2, triangle, closest point) ===\n");
    {
        std::vector<std::pair<std::string, CrackMesh>> meshes;
        for (const char* id : {"A1", "D2", "D3", "E3", "C2"}) {
            auto sc = buildStressCase(id);
            meshes.push_back({std::string(id) + " walk crack", sc->walk});
            if (std::string(id) == "D2") {
                meshes.push_back({"D2 full sphere", sc->surf[0]});
                const Ring R = stressRing(*sc);
                meshes.push_back({"D2 T+++ in xi", stressTet(R.T[0], sc->walk).mesh});
            }
        }
        {
            auto sc = buildStressCase("B1");
            meshes.push_back({"B1 tube", sc->surf[0]});
        }
        for (const auto& [name, m] : meshes) {
            long long n = 0, diff = 0;
            for (int i = 0; i < 4000; i++) {
                Vec3 x;
                if (i%2 == 0) x = Vec3(U(rng), U(rng), U(rng));
                else {   // near the mesh: a random point of a random triangle, offset along its normal
                    const auto& t = m.tri(int(rng()%m.size()));
                    double a = U01(rng), b = U01(rng);
                    if (a + b > 1) {
                        a = 1 - a;
                        b = 1 - b;
                    }
                    x = t.a + a*(t.b - t.a) + b*(t.c - t.a) + std::pow(10.0, -2 - 6*U01(rng))*(U01(rng) < 0.5 ? -1.0 : 1.0)*t.n;
                }
                const CrackMesh::Hit h1 = m.closest(x), h2 = m.closestBruteForce(x);
                n++;
                diff += !(std::memcmp(&h1.d2, &h2.d2, sizeof(double)) == 0 && h1.tri == h2.tri && std::memcmp(h1.q.data(), h2.q.data(), 3*sizeof(double)) == 0);
            }
            out("  %-22s %6zu triangles, %5d nodes: %lld points, %lld differ -> %s\n", name.c_str(), m.size(), m.nodes(), n, diff,
                diff ? "FAIL" : "ok");
            ok = ok && diff == 0;
        }
    }
    // 2. flat K1 as a mesh against StripCut
    out("\n=== self-test 2: K1 as a mesh (A1's walk crack) against check 4's StripCut on O ===\n");
    {
        auto sc = buildStressCase("A1");
        const RingCrack rc = ringCrack(0);
        const StripCut strips = ringStrips(rc);
        double maxd = 0;
        long long sideDiff = 0, n = 0;
        for (int i = 0; i < 20000; i++) {
            const Vec3 x(U(rng), U(rng), U(rng));
            const CrackMesh::Hit h = sc->walk.closest(x);
            const StripCut::Closest c = strips.closest(x, strips.nStrips());
            maxd = std::max(maxd, std::fabs(std::sqrt(h.d2) - std::sqrt(c.d2)));
            const double hgt = rc.n.dot(x) - rc.d;
            if (std::fabs(hgt) > 1e-12) sideDiff += sc->walk.side(h, x) != (hgt >= 0 ? 1 : -1);
            n++;
        }
        out("  %lld random points in O's box: max |distance difference| %.2e, side differs at %lld -> %s\n", n, maxd, sideDiff,
            maxd < 1e-14 && sideDiff == 0 ? "ok" : "FAIL");
        ok = ok && maxd < 1e-14 && sideDiff == 0;
        out("  walk crack area %.15f vs StripCut's %.15f (rel. diff %.1e)\n", sc->walk.area(), strips.area(strips.nStrips()),
            std::fabs(sc->walk.area() - strips.area(strips.nStrips()))/strips.area(strips.nStrips()));
    }
    // 3. per-tet areas, J at the crack-face lines and the control variates' jump, for every case
    out("\n=== self-test 3: every case: the tets' crack areas add up to the walk crack's; J = side g at the crack-face lines;\n"
        "    the control variates jump by exactly 2 g across the crack ===\n");
    for (int k = 0; k < kNStress; k++) {
        auto sc = buildStressCase(kStressIds[k]);
        const Ring R = stressRing(*sc);
        std::vector<StressTet> stt;
        double a = 0;
        for (const TetFrame& T : R.T) {
            stt.push_back(stressTet(T, sc->walk));
            a += stt.back().areaX;
        }
        const double aRel = std::fabs(a - sc->walk.area())/sc->walk.area();
        const StressChecks ck = stressChecks(*sc, R, stt);
        const bool good = aRel < 1e-12 && ck.maxJLine < 1e-4 && ck.maxJrLine < 1e-4 && ck.maxCvJump < 1e-4;
        out("  %-4s walk %5zu tris, area %.6f, tets' sum rel. diff %.1e | J at lines max %.1e (ramp %.1e, %lld samples) | cv jump max "
            "%.1e (%lld) -> %s\n", sc->id.c_str(), sc->walk.size(), sc->walk.area(), aRel, ck.maxJLine, ck.maxJrLine, ck.nJLine, ck.maxCvJump,
            ck.nCvJump, good ? "ok" : "FAIL");
        ok = ok && good;
    }
    // 4. the solid angle: limits, jump, Laplacian, gradient
    out("\n=== self-test 4: solid angle (C2's disk fan, D2's curved cap) ===\n");
    {
        auto c2 = buildStressCase("C2");
        auto d2 = buildStressCase("D2");
        const SolidAngleMesh& disk = c2->solid[0];
        const SolidAngleMesh& cap = d2->solid[0];
        double limErr = 0, jumpErr = 0, lapMax = 0, gradRelErr = 0;
        for (size_t f = 0; f < disk.F.size(); f += 7) {
            const Vec3 a = disk.V[size_t(disk.F[f][0])], b = disk.V[size_t(disk.F[f][1])], c = disk.V[size_t(disk.F[f][2])];
            const Vec3 q = (a + b + c)/3.0, n = (b - a).cross(c - a).normalized();
            limErr = std::max(limErr, std::max(std::fabs(disk.omega(q + 1e-9*n) - 2*kPi), std::fabs(disk.omega(q - 1e-9*n) + 2*kPi)));
        }
        for (size_t f = 0; f < cap.F.size(); f += 11) {
            const Vec3 a = cap.V[size_t(cap.F[f][0])], b = cap.V[size_t(cap.F[f][1])], c = cap.V[size_t(cap.F[f][2])];
            const Vec3 q = (a + b + c)/3.0, n = (b - a).cross(c - a).normalized();
            jumpErr = std::max(jumpErr, std::fabs(cap.omega(q + 1e-9*n) - cap.omega(q - 1e-9*n) - 4*kPi));
        }
        const CrackMesh* far[2] = {&c2->glob, &d2->surf[0]};
        for (int i = 0; i < 200; i++) {
            const Vec3 x(U(rng), U(rng), U(rng));
            for (int mi = 0; mi < 2; mi++) {
                const SolidAngleMesh* m = mi == 0 ? &disk : &cap;
                if (far[mi]->distance(x) < 0.05) continue;
                const double h = 1e-3, f0 = m->omega(x);
                double lap = 0;
                for (int a = 0; a < 3; a++) {
                    const Vec3 e = Vec3::Unit(a);
                    lap += (m->omega(x + h*e) - 2*f0 + m->omega(x - h*e))/(h*h);
                }
                const Vec3 g = m->grad(x);
                Vec3 fd;
                for (int a = 0; a < 3; a++) {
                    const Vec3 e = Vec3::Unit(a);
                    fd[a] = (m->omega(x + 1e-6*e) - m->omega(x - 1e-6*e))/2e-6;
                }
                lapMax = std::max(lapMax, std::fabs(lap)/(1.0 + g.norm()));
                gradRelErr = std::max(gradRelErr, (g - fd).norm()/(1.0 + g.norm()));
            }
        }
        const bool good = limErr < 1e-6 && jumpErr < 1e-6 && lapMax < 1e-2 && gradRelErr < 1e-5;
        out("  flat disk: one-sided limits at +-1e-9 vs +-2 pi: max err %.1e; curved cap: jump vs 4 pi max err %.1e; finite-difference "
            "Laplacian (h = 1e-3) max %.1e; Biot-Savart gradient vs central differences max rel. err %.1e -> %s\n", limErr, jumpErr, lapMax,
            gradRelErr, good ? "ok" : "FAIL");
        ok = ok && good;
    }
    // 5. determinism: the same seeds give the same sums, also with one thread
    out("\n=== self-test 5: determinism of our method on C2 (Q8/R8/B4, 256 points per tet): repeat and one thread ===\n");
    {
        auto sc = buildStressCase("C2");
        const Ring R = stressRing(*sc);
        std::vector<StressTet> stt;
        std::vector<PointSet> ps;
        for (const TetFrame& T : R.T) stt.push_back(stressTet(T, sc->walk));
        for (int t = 0; t < kRingNT; t++) ps.push_back(stressPoints(R.T[size_t(t)], stt[size_t(t)], sc->walk, 256, gateParams(false).eps));
        StressOursParams op;
        op.nQ = op.nR = 8;
        op.B = 4;
        const uint64_t seed = gateParams(false).seed;
        const StressOurs a = runStressOurs(*sc, R, stt, ps, op, seed, 123, false), b = runStressOurs(*sc, R, stt, ps, op, seed, 123, false);
        StressOurs c;
        {
            tbb::global_control gc(tbb::global_control::max_allowed_parallelism, 1);
            c = runStressOurs(*sc, R, stt, ps, op, seed, 123, false);
        }
        bool same = true;
        for (size_t k = 0; k < a.vms.size(); k++)
            same = same && std::memcmp(a.gl[k].c.data(), b.gl[k].c.data(), sizeof(double)*kRingNV) == 0 &&
                   std::memcmp(a.gl[k].c.data(), c.gl[k].c.data(), sizeof(double)*kRingNV) == 0;
        for (int t = 0; t < kRingNT; t++)
            if (a.hasR[size_t(t)]) same = same && sameSums(a.runR[size_t(t)].sums, c.runR[size_t(t)].sums);
        out("  glued c (J '%s') %s; repeat and one thread: %s\n", sc->var[0].name.c_str(), fmtVec(a.gl[0].c).c_str(), same ? "identical" : "DIFFERENT");
        ok = ok && same;
    }
    out("\nself-test: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

} // namespace

bool runRingStress(const RingStressParams& p)
{
    using clk = std::chrono::steady_clock;
    const auto tStart = clk::now();
    Report out(p.quick);
    std::time_t now = std::time(nullptr);
    std::tm tmNow{};
#ifdef _WIN32
    localtime_s(&tmNow, &now);
#else
    localtime_r(&now, &tmNow);
#endif
    char date[64], stamp[64];
    std::strftime(date, sizeof date, "%Y-%m-%d %H:%M:%S", &tmNow);
    std::strftime(stamp, sizeof stamp, "%Y%m%d_%H%M%S", &tmNow);
    out("--ring-stress%s: %s; %u hardware threads (TBB default arena)\n", p.selftest ? " --selftest" : "", date,
        std::thread::hardware_concurrency());
    if (p.selftest) return stressSelftest(out);

    // the cases
    std::vector<std::string> ids;
    const int prio[kNStress] = {1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 3, 3, 3, 3, 1, 1};
    if (p.cases == "all" || p.cases == "p1" || p.cases == "p2" || p.cases == "p3") {
        const int maxP = p.cases == "all" ? 3 : p.cases[1] - '0';
        for (int k = 0; k < kNStress; k++)
            if (p.cases == "all" ? true : prio[k] == maxP) ids.push_back(kStressIds[k]);
    } else {
        std::stringstream ss(p.cases);
        std::string item;
        while (std::getline(ss, item, ','))
            if (!item.empty()) ids.push_back(item);
    }
    for (const std::string& id : ids)
        if (std::find(std::begin(kStressIds), std::end(kStressIds), id) == std::end(kStressIds)) {
            out("unknown case %s (known: all, p1, p2, p3, or a list of", id.c_str());
            for (const char* k : kStressIds) out(" %s", k);
            out(")\n");
            return false;
        }
    std::string tag = p.tag.empty() ? std::string("stress_") + stamp : p.tag;
    std::error_code ec;
    std::filesystem::create_directories("ring_stress", ec);
    {
        std::string b = tag;
        for (int i = 2; std::filesystem::exists("ring_stress/" + b + "_rows.csv") || std::filesystem::exists("ring_stress/" + b + "_log.txt"); i++)
            b = tag + "_" + std::to_string(i);
        tag = b;
    }
    out("cases:");
    for (const std::string& id : ids) out(" %s", id.c_str());
    out("\noutput: ring_stress/%s_rows.csv, _summary.csv, _log.txt, _geometry.txt; GT caches gate_cache/ref_stress_*.bin%s\n", tag.c_str(),
        p.quick ? " (QUICK: none read or written)" : "");
    out("settings: ours (Q/R/B) %s (seed sets %d for the first, timing repeats %d: median); direct N %s; points per tet %d (A1, A2, A3: "
        "check 4's 4096); GT N = %d in %d batches, walk cap 65536; Zombie on %s; validation %s; points within %g sqrt(2) of a crack left "
        "out%s%s\n",
        p.quick ? "8/8/4" : p.ours.c_str(), p.seeds, p.repeats, p.quick ? "16,64" : p.direct.c_str(), p.quick ? 512 : p.points,
        p.quick ? 1024 : p.gtWalks, p.quick ? 16 : p.gtBatches, p.zombie.c_str(), p.validate ? "on" : "off", p.epsRel,
        p.epsRel != 0.01 ? " (NOT check 4's 0.01; A1, A2, A3 keep check 4's points)" : "", p.cache ? "" : "; GT cache off");

    std::vector<StressRow> rows;
    std::vector<StressSummary> sums;
    std::string geoTxt;
    std::vector<std::vector<int>> gOurs = parseGrid(p.quick && p.ours == "32/16/8,32/64/8,1024/256/16" ? "8/8/4" : p.ours, "32/16/8");
    auto settingOf = [&](size_t i, int& nQ, int& nR, int& B) {
        nQ = gOurs[i][0];
        nR = gOurs[i].size() > 1 ? gOurs[i][1] : nQ;
        B = gOurs[i].size() > 2 ? gOurs[i][2] : defaultBatches(nQ);
    };
    auto writeFiles = [&]() {
        {
            std::ofstream f("ring_stress/" + tag + "_rows.csv", std::ios::binary);
            f << stressCsvHeader();
            for (const StressRow& r : rows) f << stressCsvRow(r);
        }
        std::string sc = stressSummaryCsvHeader();
        for (const StressSummary& s : sums) {
            const auto cs = buildStressCase(s.id);
            for (const JVariant& v : cs->var)
                for (int g = 0; g < (cs->ramp ? 2 : 1); g++) {
                    if (g == 1 && !v.Jr) continue;
                    int nQ, nR, B;
                    settingOf(0, nQ, nR, B);
                    const OursNums o0 = oursNums(rows, s.id, "ours", v.name, g, nQ, nR, B);
                    OursNums o1, o2;
                    if (gOurs.size() > 1) {
                        settingOf(1, nQ, nR, B);
                        o1 = oursNums(rows, s.id, "ours", v.name, g, nQ, nR, B);
                    }
                    if (gOurs.size() > 2) {
                        settingOf(gOurs.size() - 1, nQ, nR, B);
                        o2 = oursNums(rows, s.id, "ours", v.name, g, nQ, nR, B);
                    }
                    std::string how;
                    const double dAt = directAt(rows, s.id, g, o0.t, how);
                    const bool exactInd = bool(cs->exact) && v.name == "indicator";
                    std::vector<std::string> c = {s.id, "\"" + s.plain + "\"", std::to_string(s.priority), v.name, g ? "ramp" : "g1",
                                                  std::to_string(int(s.ran))};
                    settingOf(0, nQ, nR, B);
                    c.push_back("Q" + std::to_string(nQ) + "/R" + std::to_string(nR) + "/B" + std::to_string(B));
                    c.push_back(std::to_string(o0.n));
                    for (double x : {o0.rms, o0.rmsSd, o0.nc, o0.maxp, o0.vert, o0.grad, o0.gradNc, o0.gradAbs, o0.E, o0.t, o0.relGt, o1.rms, o2.rms,
                                     o2.maxp, o2.vert, dAt})
                        c.push_back(csvNum(x));
                    c.push_back(how);
                    std::string gt = s.gtLine, ch = s.check;
                    std::replace(gt.begin(), gt.end(), '"', '\'');
                    std::replace(ch.begin(), ch.end(), '"', '\'');
                    c.push_back("\"" + gt + "\"");
                    c.push_back("\"" + ch + "\"");
                    c.push_back(s.ran ? stressVerdict(o0, exactInd) : "fails (" + s.failure + ")");
                    c.push_back("\"" + s.failure + "\"");
                    for (size_t i = 0; i < c.size(); i++) sc += (i ? "," : "") + c[i];
                    sc += "\n";
                }
        }
        {
            std::ofstream f("ring_stress/" + tag + "_summary.csv", std::ios::binary);
            f << sc;
        }
        {
            std::ofstream f("ring_stress/" + tag + "_geometry.txt", std::ios::binary);
            f << geoTxt;
        }
        {
            std::ofstream f("ring_stress/" + tag + "_log.txt", std::ios::binary);
            f << out.text;
        }
    };
    for (const std::string& id : ids) {
        StressSummary s;
        s.id = id;
        try {
            auto sc = buildStressCase(id);
            s.plain = sc->plain;
            s.priority = sc->priority;
            runStressCase(out, p, *sc, rows, s, geoTxt);
        } catch (const std::exception& e) {
            out("  CASE %s FAILED: %s\n", id.c_str(), e.what());
            s.failure = e.what();
        }
        sums.push_back(s);
        writeFiles();
    }
    // the console summary
    out("\n=== summary (g = 1, J = each case's main face function, ours at the first setting: median over its seed sets) ===\n");
    out("  %-4s %-52s %-14s %8s %8s %8s %8s %8s %8s %9s %s\n", "case", "plain name", "J", "H RMS%", "max %", "vertex", "grad%", "time s",
        "S_last%", "direct%", "verdict");
    for (const StressSummary& s : sums) {
        if (!s.ran) {
            out("  %-4s %-52s FAILED: %s\n", s.id.c_str(), s.plain.c_str(), s.failure.c_str());
            continue;
        }
        const auto cs = buildStressCase(s.id);
        for (const JVariant& v : cs->var) {
            int nQ, nR, B;
            settingOf(0, nQ, nR, B);
            const OursNums o0 = oursNums(rows, s.id, "ours", v.name, 0, nQ, nR, B);
            settingOf(gOurs.size() - 1, nQ, nR, B);
            const OursNums o2 = oursNums(rows, s.id, "ours", v.name, 0, nQ, nR, B);
            std::string how;
            const double dAt = directAt(rows, s.id, 0, o0.t, how);
            out("  %-4s %-52s %-14s %8.3f %8.2f %8.4f %8.2f %8.3f %8.3f %9.3f %s\n", s.id.c_str(), s.plain.substr(0, 52).c_str(), v.name.c_str(),
                100*o0.rms, 100*o0.maxp, o0.vert, 100*o0.grad, o0.t, 100*o2.rms, 100*dAt,
                stressVerdict(o0, bool(cs->exact) && v.name == "indicator").c_str());
        }
    }
    out("\n--ring-stress total %.1f s; %zu rows in ring_stress/%s_rows.csv\n", secondsSince(tStart), rows.size(), tag.c_str());
    writeFiles();
    return true;
}

int stressCachedPoints(const std::string& id)
{
    const auto sc = buildStressCase(id);   // throws on an unknown id
    if (sc->exact || sc->check4 >= 0) return -1;
    const std::string pre = "ref_stress_" + id + "_p";   // ref_stress_<id>_p<P>.bin: 25600 walks, default eps
    int best = 0;
    std::error_code ec;
    for (const auto& e : std::filesystem::directory_iterator("gate_cache", ec)) {
        const std::string f = e.path().filename().string();
        if (f.rfind(pre, 0) != 0 || f.size() <= pre.size() + 4 || f.compare(f.size() - 4, 4, ".bin") != 0) continue;
        const std::string num = f.substr(pre.size(), f.size() - 4 - pre.size());
        if (num.find_first_not_of("0123456789") == std::string::npos) best = std::max(best, std::atoi(num.c_str()));
    }
    return best;
}

RingView stressRingView(const std::string& id, int pointsPerTet, int nQ, int nR, int walksRef, bool geometryOnly)
{
    const auto t0 = std::chrono::steady_clock::now();
    if (std::find(std::begin(kStressIds), std::end(kStressIds), id) == std::end(kStressIds))
        throw std::invalid_argument("--case: unknown ring stress case " + id);
    auto scp = buildStressCase(id);
    StressCase& sc = *scp;
    const GateParams gp = gateParams(false);
    const Ring R = stressRing(sc);
    Report out(false);
    RingView v;
    v.title = "ring stress case " + sc.id + ": " + sc.plain + " (" + R.name + ") [fast method]";
    std::copy(R.V.begin(), R.V.end(), v.verts.begin());
    std::vector<StressTet> stt;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        v.tets.push_back({size_t(T.global[0]), size_t(T.global[1]), size_t(T.global[2]), size_t(T.global[3])});
        v.tetName.push_back(T.name);
        stt.push_back(stressTet(T, sc.walk));
        v.tetCut.push_back(char(stt.back().cut));
    }
    for (const auto& t : sc.walk.tris()) {
        const int b = int(v.crackV.size());
        v.crackV.insert(v.crackV.end(), {t.a, t.b, t.c});
        v.crackF.push_back({b, b + 1, b + 2});
    }
    if (!sc.cracks.empty() && sc.cracks[0].halfPlane) {   // the slice clouds: the first crack's plane (flat cracks)
        v.n = sc.cracks[0].frame.n;
        v.d = sc.cracks[0].frame.d;
    }
    if (geometryOnly) return v;
    // the points and the GT (as --ring-stress made them)
    std::vector<PointSet> ps;
    RefRun ref;
    bool haveRun = false;
    std::string gtLine;
    if (sc.check4 >= 0) {
        PerfSetup s;
        perfSetup(s, 2, sc.check4, gp.nRingPoints);
        ps = s.ps;
        haveRun = loadRingGT(s, ref);
        gtLine = haveRun ? "GT: check 4's cache " + ref.file : "";
    } else {
        for (int t = 0; t < kRingNT; t++) ps.push_back(stressPoints(R.T[size_t(t)], stt[size_t(t)], sc.walk, pointsPerTet, gp.eps));
        if (!sc.exact) {
            ref = stressGT(out, sc, R, ps, 25600, 64, 5000000 + 1000*uint64_t(sc.index), 65536, stressGtName(sc.id, pointsPerTet, 25600),
                           true, nullptr, false, false, true);
            if (ref.cached) {
                haveRun = true;
                gtLine = "GT: the --ring-stress cache " + ref.file;
            }
        }
        if (!sc.exact && !haveRun) {   // no cache at this density: solved here (not cached)
            const int B = walksRef >= 1280 ? 64 : 16, N = std::max(4*B, (walksRef + B/2)/B*B);
            ref = stressSolveO(sc, R, ps, N, B, mixSeed(gp.seed, 5000000 + 1000*uint64_t(sc.index) + 300), 65536);
            haveRun = true;
            gtLine = "GT: no cache at " + std::to_string(pointsPerTet) + " points per tet; solved here at N = " + std::to_string(N) +
                     " (not cached), " + fmtd("%.1f s", ref.stats.seconds);
        }
        if (sc.exact) gtLine = "GT: exact (+-1, the side indicator)";
    }
    StressRef sr;
    sr.ref = ref;
    sr.hasRun = haveRun;
    stressRefSetup(sr, sc, R, stt, ps);
    StressOursParams op;
    op.B = defaultBatches(nQ);
    auto whole = [&](int n) { return std::max(2*op.B, (n + op.B/2)/op.B*op.B); };
    op.nQ = whole(nQ);
    op.nR = whole(nR);
    const StressOurs o = runStressOurs(sc, R, stt, ps, op, gp.seed, 6000000 + 1000*uint64_t(sc.index), true);
    std::vector<double> w;
    for (int t = 0; t < kRingNT; t++) {
        v.X.insert(v.X.end(), ps[size_t(t)].X.begin(), ps[size_t(t)].X.end());
        w.insert(w.end(), ps[size_t(t)].w.begin(), ps[size_t(t)].w.end());
        v.tetOf.insert(v.tetOf.end(), ps[size_t(t)].size(), t);
    }
    v.pointsPerTet = sc.check4 >= 0 ? gp.nRingPoints : pointsPerTet;
    v.walksTet = op.nQ;
    v.walksR = op.nR;
    v.walksRef = haveRun ? ref.sums.walksPerBatch*ref.sums.nBatches : 0;
    v.fast = true;
    char buf[600];
    std::snprintf(buf, sizeof buf, "%zu sample points; our method: Q %d / R %d walks per point, %d batches, %.2f s (walks %.2f s); %s",
                  v.X.size(), op.nQ, op.nR, op.B, o.secTotal, o.secP + o.secR, gtLine.c_str());
    v.report.push_back(buf);
    v.report.push_back(sc.desc);
    v.report.push_back("error of the glued H against the GT, as a fraction of the jump 2: RMS / max; vertex values: glued / GT");
    GateView::Field gtF{"GT (g = 1)"};
    gtF.vol = sr.gt[0].H;
    v.fields.push_back(gtF);
    for (size_t k = 0; k < o.vms.size(); k++) {
        if (o.vms[k].second != 0) continue;
        const JVariant& jv = sc.var[size_t(o.vms[k].first)];
        const std::string sfx = " [J " + jv.name + "] (g = 1)";
        GateView::Field dec{"glued H" + sfx}, dif{"error = glued H - GT" + sfx}, h0{"H0" + sfx};
        for (int t = 0; t < kRingNT; t++) {
            dec.vol.insert(dec.vol.end(), o.F[k][size_t(t)].H.begin(), o.F[k][size_t(t)].H.end());
            for (size_t i = 0; i < ps[size_t(t)].size(); i++) {
                double H, vH;
                Vec3 G, vG;
                stressH0(o, int(k), t, ps[size_t(t)], R.T[size_t(t)], i, H, vH, G, vG);
                h0.vol.push_back(H);
            }
        }
        dif.difference = true;
        double sq = 0, sw = 0;
        for (size_t i = 0; i < dec.vol.size(); i++) {
            dif.vol.push_back(dec.vol[i] - gtF.vol[i]);
            sq += w[i]*dif.vol[i]*dif.vol[i];
            sw += w[i];
            dif.max = std::max(dif.max, std::fabs(dif.vol[i]));
        }
        dif.rms = std::sqrt(sq/sw)/kTwoGMax;
        dif.max /= kTwoGMax;
        std::string vt;
        for (int q = 0; q < kRingNV; q++) {
            std::snprintf(buf, sizeof buf, "%s v%d %+.3f/%+.3f", q ? "," : "", q, o.gl[k].c[q] + jv.J(R.V[size_t(q)]), sr.gtP[0].H[size_t(q)]);
            vt += buf;
        }
        std::snprintf(buf, sizeof buf, "glued H%s: error RMS %.2f%%, max %.1f%%; vertices%s", sfx.c_str(), 100*dif.rms, 100*dif.max, vt.c_str());
        v.report.push_back(buf);
        v.fields.push_back(std::move(dec));
        v.fields.push_back(std::move(dif));
        v.fields.push_back(std::move(h0));
    }
    // the basis functions of the main J's run: P_k = lambda_k + Q_k in every tet with corner v
    for (int q = 0; q < kRingNV; q++) {
        v.basis[size_t(q)].assign(v.X.size(), 0.0);
        v.inSupport[size_t(q)].assign(v.X.size(), 0);
    }
    size_t off = 0;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        for (size_t i = 0; i < ps[size_t(t)].size(); i++) {
            const Vec3& xi = ps[size_t(t)].xi[i];
            const double lam[4] = {1.0 - xi.sum(), xi[0], xi[1], xi[2]};
            for (int q = 0; q < 4; q++) {
                double qv = 0, var;
                long long n;
                if (o.hasP[size_t(t)]) comboH(o.runP[size_t(t)].sums, i, unitCoef(4, q), qv, var, n);
                v.basis[size_t(T.global[q])][off + i] = lam[q] + qv;
                v.inSupport[size_t(T.global[q])][off + i] = 1;
            }
        }
        off += ps[size_t(t)].size();
    }
    std::snprintf(buf, sizeof buf, "total %.1f s", secondsSince(t0));
    v.report.push_back(buf);
    return v;
}

// =====================================================================================
// --gate coupling --export-xfem --check 4: one nodal ring's glued H for an XFEM simulation (coupling_gate.h runXfemExport)
// =====================================================================================
// Our method only (the fast path), no GT: nothing is read from or written to gate_cache/, gate_report.md, ring_perf/ or
// ring_stress/, and no other mode calls this. The glue is gateRingViewFast's (its point sets, walks and seeds: the same
// c, bit for bit). The export points (the quadrature, then the self-checks' points at weight 0) are a second runRingFast
// pass with their own seeds, the same xi cuts and face data, evaluated with the glued c by cvField: H = J_xi + R + sum_k
// c_k (lambda_k + Q_k). That pass also glues its own points, printed as a diagnostic only. The render data are exact: on
// O's faces H is its face data J + sum_k c_k lambda_k, on the crack +-1.
namespace {

constexpr double kXfemFrontBand[2] = {0.15, 0.05};   // quadrature leaves within these of O's front: one more level, two
constexpr int kXfemRenderPerEdge = 24;   // the render's outer triangles (edges sqrt 2) are cut into this many per edge
constexpr int kXfemViewStride = 64;      // self-check: every 64th glue point is also an export point
constexpr double kXfemFaceOff = 1e-3;    // self-check: points this far inside O's faces (h against the face data)
constexpr double kXfemFaceCrack = 0.05;  // ... at least this far from the crack
constexpr double kXfemCornerOff = 1e-4;  // self-check: points this far from a tet's corners (h against the node value)

// motus topology/cutting/xfem_cutter.cpp positive_piece, ported: the part of a tet where a linear level set with corner
// values phi is positive, as 0, 1 or 3 sub-tets with corners in barycentric coordinates of the tet (4 positive: the tet;
// 3-1: the little tet at a lone positive corner, or the wedge under a lone negative one; 2-2: a wedge). The negative part
// is the positive part of -phi, so the two pieces are motus's +3/-1 or 2+/2- split.
int xfemPositivePiece(const Eigen::Vector4d& phi, Eigen::Vector4d (&tets)[3][4])
{
    const Eigen::Matrix4d I = Eigen::Matrix4d::Identity();
    auto crossing = [&phi, &I](int i, int j) -> Eigen::Vector4d {   // the zero of phi on edge (i, j)
        const double t = phi[i]/(phi[i] - phi[j]);
        return (1.0 - t)*I.col(i) + t*I.col(j);
    };
    auto set = [&tets](int k, const Eigen::Vector4d& a, const Eigen::Vector4d& b, const Eigen::Vector4d& c,
                       const Eigen::Vector4d& d) {
        tets[k][0] = a;
        tets[k][1] = b;
        tets[k][2] = c;
        tets[k][3] = d;
    };
    int positive = 0;
    for (int a = 0; a < 4; a++) positive += phi[a] > 0 ? 1 : 0;
    if (positive == 0) return 0;
    if (positive == 4) {
        set(0, I.col(0), I.col(1), I.col(2), I.col(3));
        return 1;
    }
    if (positive == 1 || positive == 3) {   // 3-1: L is the corner whose sign the other three do not share
        const bool lonePositive = positive == 1;
        int L = 0;
        for (int a = 0; a < 4; a++)
            if ((phi[a] > 0) == lonePositive) L = a;
        const int i = (L + 1)%4, j = (L + 2)%4, k = (L + 3)%4;
        const Eigen::Vector4d ci = crossing(L, i), cj = crossing(L, j), ck = crossing(L, k);
        if (lonePositive) {   // the corner chopped off at L
            set(0, I.col(L), ci, cj, ck);
            return 1;
        }
        set(0, I.col(i), I.col(j), I.col(k), ci);   // the wedge (ci, cj, ck) over (e_i, e_j, e_k)
        set(1, I.col(j), I.col(k), ci, cj);
        set(2, I.col(k), ci, cj, ck);
        return 3;
    }
    int p[2] = {0, 0}, n[2] = {0, 0}, np = 0, nn = 0;   // 2-2: c_ab is the crossing on edge (p_a, n_b)
    for (int a = 0; a < 4; a++) {
        if (phi[a] > 0) p[np++] = a;
        else n[nn++] = a;
    }
    const Eigen::Vector4d c00 = crossing(p[0], n[0]), c01 = crossing(p[0], n[1]), c10 = crossing(p[1], n[0]),
                          c11 = crossing(p[1], n[1]);
    set(0, I.col(p[0]), c00, c01, I.col(p[1]));   // the wedge (p0, c00, c01) to (p1, c10, c11)
    set(1, c00, c01, I.col(p[1]), c10);
    set(2, c01, I.col(p[1]), c10, c11);
    return 3;
}

// The distance from the segment [a, b] to the leaf l, 0 where the segment meets it. A point's distance to a convex set is
// convex along a segment, so a golden-section search over the segment's parameter finds its minimum (to ~1e-9 of |b - a|).
// A leaf's centroid is no proxy: the sub-tets next to the plane are slivers whose centroid can lie far from a front that
// passes through them.
double xfemLeafSegmentDist(const Leaf& l, const Vec3& a, const Vec3& b)
{
    Eigen::Matrix3d E;
    E << l.v[1] - l.v[0], l.v[2] - l.v[0], l.v[3] - l.v[0];
    const Eigen::Matrix3d Einv = E.inverse();
    auto toLeaf = [&](double t) {   // 0 inside, else the nearest face
        const Vec3 p = a + t*(b - a);
        const Vec3 lam = Einv*(p - l.v[0]);
        if (lam.minCoeff() >= 0 && lam.sum() <= 1) return 0.0;
        return std::min({pointTriangleDist(p, l.v[1], l.v[2], l.v[3]), pointTriangleDist(p, l.v[0], l.v[2], l.v[3]),
                         pointTriangleDist(p, l.v[0], l.v[1], l.v[3]), pointTriangleDist(p, l.v[0], l.v[1], l.v[2])});
    };
    constexpr double r = 0.61803398874989484820;
    double lo = 0, hi = 1, x1 = hi - r*(hi - lo), x2 = lo + r*(hi - lo), f1 = toLeaf(x1), f2 = toLeaf(x2);
    for (int it = 0; it < 44 && f1 > 0 && f2 > 0; it++) {
        if (f1 <= f2) {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - r*(hi - lo);
            f1 = toLeaf(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + r*(hi - lo);
            f2 = toLeaf(x2);
        }
    }
    return std::min({f1, f2, toLeaf(0.0), toLeaf(1.0)});
}

double leafDiameter(const Leaf& l)
{
    double d = 0;
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++) d = std::max(d, (l.v[size_t(i)] - l.v[size_t(j)]).norm());
    return d;
}

// One ring tet's quadrature rule in X: the tet split by the crack plane into the sub-tets of its + and - pieces
// (xfemPositivePiece; a tet the plane misses is one piece), every sub-tet refined level times (Bey's 8 children), then once
// more where a leaf comes within kXfemFrontBand[0] of O's front segment and once more within kXfemFrontBand[1] (the leaf's
// own distance to the segment, xfemLeafSegmentDist: every leaf the front passes through ends at depth level + 2), and the
// 4-point degree-2 Gauss rule in every leaf. The points lie strictly inside their leaves, so never on the plane, and the
// weights of the tet sum to its volume.
struct XfemQuadPoint {
    Vec3 X = Vec3::Zero();
    double w = 0;
    int side = 0;   // the piece: +1 (n.X > d) or -1
};
struct XfemTetQuad {
    std::vector<XfemQuadPoint> pts;
    int nSub[2] = {0, 0};                // sub-tets of the + and the - piece
    int nLeaves = 0, nFrontLeaves = 0;   // leaves; of them, refined for the front
    int nNear = 0, nNearShallow = 0;     // leaves within kXfemFrontBand[1] of the front; of them, not at depth level + 2
    double nearDiam = 0;                 // the largest diameter among those leaves
};
XfemTetQuad xfemTetQuad(const TetFrame& T, const RingCrack& rc, int level, const std::vector<Vec3>& front)
{
    constexpr double ga = 0.58541019662496845446, gb = 0.13819660112501051518;   // (5 + 3 sqrt 5)/20, (5 - sqrt 5)/20
    XfemTetQuad q;
    Eigen::Vector4d phi;
    for (int a = 0; a < 4; a++) phi[a] = rc.n.dot(T.X[a]) - rc.d;
    const double V = std::fabs(T.detJ)/6.0;
    auto nearFront = [&front](const Leaf& l, double band) {
        return front.size() == 2 && xfemLeafSegmentDist(l, front[0], front[1]) < band;
    };
    for (int s = 0; s < 2; s++) {
        Eigen::Vector4d sub[3][4];
        q.nSub[s] = xfemPositivePiece(s == 0 ? phi : Eigen::Vector4d(-phi), sub);
        std::vector<Leaf> leaves;
        for (int k = 0; k < q.nSub[s]; k++) {
            Leaf l;
            for (int c = 0; c < 4; c++)
                l.v[c] = sub[k][c][0]*T.X[0] + sub[k][c][1]*T.X[1] + sub[k][c][2]*T.X[2] + sub[k][c][3]*T.X[3];
            if (leafVolume(l) > 1e-15*V) leaves.push_back(l);   // a zero-volume sub-tet (a corner on the plane) carries nothing
        }
        auto refine = [&leaves](auto pred) {
            std::vector<Leaf> o;
            for (const Leaf& l : leaves) {
                if (!pred(l)) {
                    o.push_back(l);
                    continue;
                }
                for (const Leaf& ch : beyChildren(l)) o.push_back(ch);
            }
            leaves.swap(o);
        };
        for (int i = 0; i < level; i++) refine([](const Leaf&) { return true; });
        for (int e = 0; e < 2; e++)
            refine([&](const Leaf& l) { return l.depth == level + e && nearFront(l, kXfemFrontBand[e]); });
        for (const Leaf& l : leaves) {
            if (l.depth > level) q.nFrontLeaves++;
            if (nearFront(l, kXfemFrontBand[1])) {
                q.nNear++;
                q.nNearShallow += l.depth != level + 2 ? 1 : 0;
                q.nearDiam = std::max(q.nearDiam, leafDiameter(l));
            }
            const double w = 0.25*leafVolume(l);
            const Vec3 sum = l.v[0] + l.v[1] + l.v[2] + l.v[3];
            for (int c = 0; c < 4; c++) q.pts.push_back(XfemQuadPoint{(ga - gb)*l.v[c] + gb*sum, w, s == 0 ? 1 : -1});
        }
        q.nLeaves += int(leaves.size());
    }
    return q;
}

// The render data in X: each vertex has its tet, side, group (0: O's boundary, 1 / 2: the crack's + / - sheet) and
// whether it lies on the plane; each triangle its three vertices and group
struct XfemRender {
    std::vector<int> tet, side, group;
    std::vector<char> onPlane;
    std::vector<Vec3> X;
    std::vector<double> h;
    std::vector<std::array<int, 4>> tris;   // a, b, c, group
    int add(int t, int s, int g, const Vec3& x, bool plane)
    {
        tet.push_back(t);
        side.push_back(s);
        group.push_back(g);
        onPlane.push_back(char(plane));
        X.push_back(x);
        return int(X.size()) - 1;
    }
};

// O's boundary triangle of tet t (corners P, wound outwards) cut into n x n small triangles, each clipped by the crack
// plane into its parts on either side (Sutherland-Hodgman per side; a grid vertex on the plane, and every crossing of a grid
// edge, is one render vertex per side). Group 0. The vertices are shared within the face, so each side's mesh is connected.
void xfemRenderFace(XfemRender& r, int t, const std::array<Vec3, 3>& P, const RingCrack& rc, int n)
{
    std::vector<Vec3> G(size_t((n + 1)*(n + 1)), Vec3::Zero());
    std::vector<double> f(G.size(), 0.0);
    auto gid = [n](int i, int j) { return i*(n + 1) + j; };
    for (int i = 0; i <= n; i++)
        for (int j = 0; i + j <= n; j++) {
            G[size_t(gid(i, j))] = P[0] + (double(i)/n)*(P[1] - P[0]) + (double(j)/n)*(P[2] - P[0]);
            f[size_t(gid(i, j))] = rc.n.dot(G[size_t(gid(i, j))]) - rc.d;
        }
    std::map<std::array<int, 3>, int> ids;   // (grid vertex, -1, side) or (grid edge lo, hi, side) -> render vertex
    auto vertex = [&](int a, int b, int s) {
        const std::array<int, 3> key = b < 0 ? std::array<int, 3>{a, -1, s} : std::array<int, 3>{std::min(a, b), std::max(a, b), s};
        const auto it = ids.find(key);
        if (it != ids.end()) return it->second;
        const size_t lo = size_t(key[0]), hi = size_t(key[1] < 0 ? key[0] : key[1]);
        const int id = b < 0 ? r.add(t, s, 0, G[lo], f[lo] == 0) : r.add(t, s, 0, G[lo] + (f[lo]/(f[lo] - f[hi]))*(G[hi] - G[lo]), true);
        ids.emplace(key, id);
        return id;
    };
    auto clip = [&](int a, int b, int c) {
        const int v[3] = {a, b, c};
        for (int s = 1; s >= -1; s -= 2) {
            std::vector<int> poly;
            for (int e = 0; e < 3; e++) {
                const double fp = s*f[size_t(v[e])], fq = s*f[size_t(v[(e + 1)%3])];
                if (fp >= 0) poly.push_back(vertex(v[e], -1, s));
                if ((fp > 0 && fq < 0) || (fp < 0 && fq > 0)) poly.push_back(vertex(v[e], v[(e + 1)%3], s));
            }
            for (size_t k = 1; k + 1 < poly.size(); k++) {
                const Vec3 &x0 = r.X[size_t(poly[0])], &x1 = r.X[size_t(poly[k])], &x2 = r.X[size_t(poly[k + 1])];
                if ((x1 - x0).cross(x2 - x0).norm() > 1e-14) r.tris.push_back({poly[0], poly[k], poly[k + 1], 0});
            }
        }
    };
    for (int i = 0; i < n; i++)
        for (int j = 0; i + j < n; j++) {
            clip(gid(i, j), gid(i + 1, j), gid(i, j + 1));
            if (i + j < n - 1) clip(gid(i + 1, j), gid(i + 1, j + 1), gid(i, j + 1));
        }
}

// The crack inside tet t (its polygon P in X, CCW about n), fanned about its centroid, every fan triangle cut into m x m
// (m from the longest fan edge over the render spacing hR, one m for the whole polygon so that the fans' shared edges
// match), written twice: group 1 (side +1, wound with normal -n: outwards from the + piece) and group 2 (side -1, normal
// +n). Returns the fan triangles whose winding disagreed with n (0 expected).
int xfemRenderCrack(XfemRender& r, int t, const std::vector<Vec3>& P, const Vec3& nrm, double hR)
{
    const int N = int(P.size());
    if (N < 3) return 0;
    Vec3 ctr = Vec3::Zero();
    for (const Vec3& p : P) ctr += p/double(N);
    double len = 0;
    int bad = 0;
    for (int i = 0; i < N; i++) {
        len = std::max({len, (P[size_t(i)] - ctr).norm(), (P[size_t((i + 1)%N)] - P[size_t(i)]).norm()});
        if ((P[size_t(i)] - ctr).cross(P[size_t((i + 1)%N)] - ctr).dot(nrm) <= 0) bad++;
    }
    const int m = std::max(1, int(std::ceil(len/hR)));
    for (int s = 1; s >= -1; s -= 2) {
        const int g = s > 0 ? 1 : 2;
        std::map<std::array<int, 3>, int> ids;   // (spoke e, k, -1): on the segment centroid-P_e; (N + fan i, a, b): inside fan i
        auto vertex = [&](int i, int a, int b) {
            const int j = (i + 1)%N;
            std::array<int, 3> key{N + i, a, b};
            if (b == 0) key = {a == 0 ? 0 : i, a, -1};
            else if (a == 0) key = {j, b, -1};
            const auto it = ids.find(key);
            if (it != ids.end()) return it->second;
            const Vec3 x = key[2] < 0 ? Vec3(ctr + (double(key[1])/m)*(P[size_t(key[0])] - ctr))
                                      : Vec3(ctr + (double(a)/m)*(P[size_t(i)] - ctr) + (double(b)/m)*(P[size_t(j)] - ctr));
            const int id = r.add(t, s, g, x, true);
            ids.emplace(key, id);
            return id;
        };
        auto tri = [&](int a, int b, int c) {
            const Vec3 &x0 = r.X[size_t(a)], &x1 = r.X[size_t(b)], &x2 = r.X[size_t(c)];
            if ((x1 - x0).cross(x2 - x0).norm() <= 1e-14) return;
            if (s > 0) r.tris.push_back({a, c, b, g});
            else r.tris.push_back({a, b, c, g});
        };
        for (int i = 0; i < N; i++)
            for (int a = 0; a < m; a++)
                for (int b = 0; a + b < m; b++) {
                    tri(vertex(i, a, b), vertex(i, a + 1, b), vertex(i, a, b + 1));
                    if (a + b < m - 1) tri(vertex(i, a + 1, b), vertex(i, a + 1, b + 1), vertex(i, a, b + 1));
                }
    }
    return bad;
}

// printf appended to a string (the export file is built in memory, then written in one go)
void appendf(std::string& s, const char* fmt, ...)
{
    char buf[2048];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    s += buf;
}

// quantiles of v (q in [0, 1]; v is sorted here)
double xfemQuantile(std::vector<double>& v, double q)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, size_t(q*double(v.size() - 1) + 0.5))];
}

} // namespace

bool buildXfemRing(const XfemExportParams& p, XfemRingData& xd)
{
    using clk = std::chrono::steady_clock;
    const auto tStart = clk::now();
    if (p.ring != 1 && p.ring != 2) throw std::invalid_argument("--export-xfem: --ring 1 or 2");
    if (p.crack < 0 || p.crack > 2) throw std::invalid_argument("--export-xfem: --config 0, 1 or 2 (K1, K2, K3)");
    if (p.quadLevel < 0 || p.quadLevel > 4) throw std::invalid_argument("--export-xfem: --quad-level 0 to 4");
    Report out(false);
    const GateParams gp = gateParams(false);
    out("--gate coupling --export-xfem --check 4: %u hardware threads (TBB default arena)\n", std::thread::hardware_concurrency());

    // 1. the glue, as the fast ring view runs it (gateRingViewFast at its --fast default: Q 32 / R 16, 8 batches, its seeds)
    auto t0 = clk::now();
    const int nPts = std::max(64, p.gluePoints);
    const int Bg = defaultBatches(32);
    PerfSetup s;
    perfSetup(s, p.ring, p.crack, nPts);
    const Ring& R = s.R;
    const RingCrack& rc = s.rc;
    const int nModes = R.id == 1 && !rc.through ? 2 : 1;   // as the view (its ramp run also changes R's walks in ring 1)
    RingFastParams fg;
    fg.nQ = 32;
    fg.nR = 16;
    fg.nBatches = Bg;
    fg.ramp = nModes == 2;
    const RingFast glue = runRingFast(R, rc, s.ps, s.J, fg, gp.seed, 21000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki), false);
    const Vec7 c = glue.gl[0].c;
    const double secGlue = secondsSince(t0);
    out("check 4, %s, crack %s: %s\n", R.name.c_str(), rc.name.c_str(), rc.desc.c_str());
    out("glue (our fast method, as the ring view: %d points per tet asked, %zu in all; Q %d / R %d walks per point, %d batches): "
        "%.3f s (Q walks %.3f s, R walks %.3f s; setup and point sets %.3f s)\n", nPts, totalPoints(s), fg.nQ, fg.nR, Bg, secGlue,
        glue.secQ, glue.secR, s.secGeometry);

    // 2. the geometry: the tet classes, O's front, the quadrature rules and the render data
    t0 = clk::now();
    std::vector<Vec3> front;   // O's front segment: its two ends on O's boundary (none when the crack is cut through)
    if (!rc.through) front = frontVertices(rc.cutO);
    if (!rc.through && front.size() != 2) throw std::runtime_error("--export-xfem: O's crack has no front segment");
    const ConvexDomain O = octahedron();
    std::array<int, kRingNT> cls{};
    std::vector<XfemTetQuad> quad;
    for (int t = 0; t < kRingNT; t++) {
        const RingTetCrack& rt = s.rt[size_t(t)];
        cls[size_t(t)] = !rt.cut() ? 0 : rt.tRaw >= 1 ? 1 : 2;
        if (R.T[size_t(t)].global[0] != 0) throw std::logic_error("--export-xfem: a ring tet's local corner 0 is not the node");
        quad.push_back(xfemTetQuad(R.T[size_t(t)], rc, p.quadLevel, front));
    }
    XfemRender rd;
    const double hR = kSqrt2/kXfemRenderPerEdge;
    std::array<std::array<Vec3, 3>, kRingNT> outer{};   // each tet's face on O, wound outwards
    int badFans = 0;
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        std::array<Vec3, 3> P = {T.X[1], T.X[2], T.X[3]};
        if ((P[1] - P[0]).cross(P[2] - P[0]).dot(P[0] - T.X[0]) < 0) std::swap(P[1], P[2]);
        outer[size_t(t)] = P;
        xfemRenderFace(rd, t, P, rc, kXfemRenderPerEdge);
    }
    for (int t = 0; t < kRingNT; t++)
        if (cls[size_t(t)] > 0) badFans += xfemRenderCrack(rd, t, s.rt[size_t(t)].phys.polygon, rc.n, hR);
    const double secGeom = secondsSince(t0);

    // 3. the export pass: per tet the quadrature points, then (weight 0) every kXfemViewStride-th glue point, points
    // kXfemFaceOff inside O's face (away from the crack) and points kXfemCornerOff from each corner
    t0 = clk::now();
    std::vector<PointSet> pe(kRingNT);
    std::vector<size_t> nQuad(kRingNT), nView(kRingNT), nFace(kRingNT);
    std::vector<std::vector<size_t>> viewIdx(kRingNT);
    std::vector<std::vector<Vec3>> faceFoot(kRingNT);
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        PointSet& e = pe[size_t(t)];
        auto add = [&](const Vec3& X, const Vec3& xi, double w) {
            e.X.push_back(X);
            e.xi.push_back(xi);
            e.w.push_back(w);
        };
        for (const XfemQuadPoint& q : quad[size_t(t)].pts) add(q.X, T.toXi(q.X), q.w);
        nQuad[size_t(t)] = e.size();
        const PointSet& g = s.ps[size_t(t)];
        for (size_t i = 0; i < g.size(); i += kXfemViewStride) {
            viewIdx[size_t(t)].push_back(i);
            add(g.X[i], g.xi[i], 0.0);
        }
        nView[size_t(t)] = e.size();
        const std::array<Vec3, 3>& P = outer[size_t(t)];
        const Vec3 nOut = (P[1] - P[0]).cross(P[2] - P[0]).normalized();
        for (int i = 1; i < 8; i++)
            for (int j = 1; i + j < 8; j++) {
                const Vec3 foot = P[0] + (i/8.0)*(P[1] - P[0]) + (j/8.0)*(P[2] - P[0]);
                if (rc.cutO.distToCut(foot) < kXfemFaceCrack) continue;
                faceFoot[size_t(t)].push_back(foot);
                add(foot - kXfemFaceOff*nOut, T.toXi(foot - kXfemFaceOff*nOut), 0.0);
            }
        nFace[size_t(t)] = e.size();
        const Vec3 ctr = 0.25*(T.X[0] + T.X[1] + T.X[2] + T.X[3]);
        for (int k = 0; k < 4; k++) {
            const Vec3 X = T.X[size_t(k)] + kXfemCornerOff*(ctr - T.X[size_t(k)]).normalized();
            add(X, T.toXi(X), 0.0);
        }
    }
    const int Be = defaultBatches(p.walksQ);
    auto whole = [](int n, int b) { return std::max(2*b, (n + b/2)/b*b); };
    RingFastParams fe;
    fe.nQ = whole(p.walksQ, Be);
    fe.nR = whole(p.walksR, Be);
    fe.nBatches = Be;
    const RingFast ex = runRingFast(R, rc, pe, s.J, fe, gp.seed, 31000 + 1000*uint64_t(s.ri) + 100*uint64_t(s.ki), false);
    const double secExWalks = secondsSince(t0);
    t0 = clk::now();
    std::vector<FieldAt> FE, FV;   // the export pass and the glue's own field (at the glue points), both with the glued c
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        FE.push_back(cvField(ex.cv[size_t(t)], pe[size_t(t)], T, 0, localValues(T, c), {}, Be));
        FV.push_back(cvField(glue.cv[size_t(t)], s.ps[size_t(t)], T, 0, localValues(T, c), {}, Bg));
    }
    const double secFields = secondsSince(t0);
    size_t nQuadAll = 0, nCheckAll = 0;
    for (int t = 0; t < kRingNT; t++) {
        nQuadAll += nQuad[size_t(t)];
        nCheckAll += pe[size_t(t)].size() - nQuad[size_t(t)];
    }
    out("export points (the same xi cuts, face data and glued c; Q %d / R %d walks per point, %d batches): %zu quadrature + %zu "
        "check points, walks %.3f s (Q %.3f s in %d tets, R %.3f s in %d tets), fields %.3f s; geometry %.3f s\n", fe.nQ, fe.nR, Be,
        nQuadAll, nCheckAll, secExWalks, ex.secQ, ex.nQtets, ex.secR, ex.nRtets, secFields, secGeom);

    // 4. the values: h = (H/g + 1)/2 and g = grad_X h at the quadrature points, the node values, the render vertices
    auto faceH = [&](int t, const Vec3& X, int side, bool plane) {   // H's face data at X in tet t (side: on the plane)
        const TetFrame& T = R.T[size_t(t)];
        const Vec3 xi = T.toXi(X);
        const double lam[4] = {1.0 - xi.sum(), xi[0], xi[1], xi[2]};
        double H = plane ? (rc.cutO.u.dot(X) - rc.sJ < 0 ? double(side) : 0.0) : s.J[0](X);   // J = theta/pi; on the plane +-1 or 0
        for (int k = 0; k < 4; k++) H += c[T.global[size_t(k)]]*lam[k];
        return H;
    };
    std::array<double, kRingNV> Hv{}, hv{}, Jv{}, phiV{};
    for (int v = 0; v < kRingNV; v++) {
        Jv[size_t(v)] = s.J[0](R.V[size_t(v)]);
        Hv[size_t(v)] = Jv[size_t(v)] + c[v];
        hv[size_t(v)] = 0.5*(Hv[size_t(v)] + 1.0);
        phiV[size_t(v)] = rc.n.dot(R.V[size_t(v)]) - rc.d;
    }
    rd.h.resize(rd.X.size());
    for (size_t i = 0; i < rd.X.size(); i++)
        rd.h[i] = rd.group[i] == 1 ? 1.0 : rd.group[i] == 2 ? 0.0 : 0.5*(faceH(rd.tet[i], rd.X[i], rd.side[i], rd.onPlane[i]) + 1.0);

    // 5. the self-checks
    out("\nself-checks\n");
    bool ok = true;
    double wErr = 0;
    for (int t = 0; t < kRingNT; t++) {
        double sw = 0;
        for (const XfemQuadPoint& q : quad[size_t(t)].pts) sw += q.w;
        const double V = std::fabs(R.T[size_t(t)].detJ)/6.0;
        wErr = std::max(wErr, std::fabs(sw - V)/V);
    }
    ok = ok && wErr <= 1e-12;
    out("  weights: max over the tets of |sum w - V| / V = %.2e (limit 1e-12) %s\n", wErr, wErr <= 1e-12 ? "ok" : "FAILED");
    if (front.size() == 2) {   // every leaf near the front refined twice, whatever its shape
        int nNear = 0, nShallow = 0;
        double diam = 0;
        for (const XfemTetQuad& q : quad) {
            nNear += q.nNear;
            nShallow += q.nNearShallow;
            diam = std::max(diam, q.nearDiam);
        }
        ok = ok && nShallow == 0;
        out("  front refinement: %d leaves within %.2f of the front (each leaf's own distance to the segment), %d of them not at "
            "depth %d; the largest is %.3f across %s\n", nNear, kXfemFrontBand[1], nShallow, p.quadLevel + 2, diam,
            nShallow == 0 ? "ok" : "FAILED");
    }
    int badClass = 0;   // class 2 iff the tet's crack polygon has an edge on O's front
    for (int t = 0; t < kRingNT; t++) badClass += (cls[size_t(t)] == 2) != !s.rt[size_t(t)].front.empty() ? 1 : 0;
    ok = ok && badClass == 0;
    out("  classes: %d tets where class 2 and 'the tet's crack has an edge on the front' disagree %s\n", badClass,
        badClass == 0 ? "ok" : "FAILED");
    long long badSide = 0, badRender = 0;
    double minPhi = kInf, maxPlane = 0;
    for (int t = 0; t < kRingNT; t++)
        for (const XfemQuadPoint& q : quad[size_t(t)].pts) {
            const double ph = rc.n.dot(q.X) - rc.d;
            minPhi = std::min(minPhi, std::fabs(ph));
            if ((ph > 0 ? 1 : -1) != q.side) badSide++;
        }
    for (size_t i = 0; i < rd.X.size(); i++) {
        const double ph = rc.n.dot(rd.X[i]) - rc.d;
        if (rd.onPlane[i]) maxPlane = std::max(maxPlane, std::fabs(ph));
        else if ((ph > 0 ? 1 : -1) != rd.side[i]) badRender++;
    }
    ok = ok && badSide == 0 && badRender == 0 && maxPlane <= 1e-12;
    out("  sides: quadrature points with side != sign(n.X - d): %lld of %zu (closest to the plane: %.2e); render vertices off the "
        "plane with the wrong side: %lld, on it: |n.X - d| <= %.1e %s\n", badSide, nQuadAll, minPhi, badRender, maxPlane,
        badSide == 0 && badRender == 0 && maxPlane <= 1e-12 ? "ok" : "FAILED");
    long long badSheet = 0, nSheet = 0;
    for (size_t i = 0; i < rd.X.size(); i++)
        if (rd.group[i] > 0) {
            nSheet++;
            if (rd.h[i] != (rd.group[i] == 1 ? 1.0 : 0.0) || rd.side[i] != (rd.group[i] == 1 ? 1 : -1)) badSheet++;
        }
    ok = ok && badSheet == 0 && badFans == 0;
    out("  crack sheets: %lld vertices, %lld with h != 1 (+ sheet) / 0 (- sheet) or the wrong side; fan triangles wound against "
        "n: %d %s\n", nSheet, badSheet, badFans, badSheet == 0 && badFans == 0 ? "ok" : "FAILED");
    // h at the render vertices on O's faces is the face data itself; check it against the export pass kXfemFaceOff inside,
    // carried to the face to first order with its own gradient (h + off grad h . n_out)
    {
        double s2 = 0, mx = 0, z2 = 0, r2 = 0;
        long long n = 0;
        for (int t = 0; t < kRingNT; t++) {
            const std::array<Vec3, 3>& P = outer[size_t(t)];
            const Vec3 nOut = (P[1] - P[0]).cross(P[2] - P[0]).normalized();
            for (size_t k = 0; k < faceFoot[size_t(t)].size(); k++) {
                const size_t i = nView[size_t(t)] + k;
                const FieldAt& F = FE[size_t(t)];
                const double hW = 0.5*(F.H[i] + 1.0), hF = 0.5*(faceH(t, faceFoot[size_t(t)][k], 0, false) + 1.0);
                const double d = hW + kXfemFaceOff*0.5*F.G[i].dot(nOut) - hF;
                const double v = 0.25*(F.vH[i] + kXfemFaceOff*kXfemFaceOff*F.vG[i].dot(nOut.cwiseAbs2()));
                r2 += (hW - hF)*(hW - hF);
                s2 += d*d;
                mx = std::max(mx, std::fabs(d));
                z2 += v > 0 ? d*d/v : 0.0;
                n++;
            }
        }
        const double nn = double(std::max(1LL, n));
        out("  render h on O's faces = the face data (J + sum c lambda + 1)/2; the export pass %.0e inside the faces (%lld points "
            "at least %.2f from the crack), carried to the face with its gradient, against it: difference RMS %.2e, max %.2e, z "
            "RMS %.2f (without the gradient term: RMS %.2e)\n", kXfemFaceOff, n, kXfemFaceCrack, std::sqrt(s2/nn), mx,
            std::sqrt(z2/nn), std::sqrt(r2/nn));
        // where the crack meets O's faces behind the front, the face data (side + sum c lambda) is not the crack data (side)
        double gap = 0;
        long long ng = 0;
        for (size_t i = 0; i < rd.X.size(); i++)
            if (rd.group[i] == 0 && rd.onPlane[i] && rc.cutO.u.dot(rd.X[i]) - rc.sJ < 0) {
                gap = std::max(gap, std::fabs(rd.h[i] - (rd.side[i] > 0 ? 1.0 : 0.0)));
                ng++;
            }
        out("  render h on O's faces where the crack meets them (%lld vertices on the plane behind the front): max |h - h of the "
            "crack sheet| = %.4f (the face data's sum c lambda / 2: the per-tet problem's face/crack data mismatch, so the boundary "
            "and the sheets do not meet exactly in h there)\n", ng, gap);
    }
    // the node values, and the export pass kXfemCornerOff from each corner (averaged over the tets around the vertex)
    out("  node values h_v = (H(X_v) + 1)/2, H(X_v) = J(X_v) + c_v (H's face data at the vertex, every tet around it alike):\n"
        "    vertex      n.X-d    c_v       J(X_v)    H(X_v)    h_v       h %.0e from it (export pass, mean of its tets)\n",
        kXfemCornerOff);
    double jMax = 0, cornerMax = 0;
    for (int v = 0; v < kRingNV; v++) {
        double hs = 0;
        int nt = 0;
        for (int t = 0; t < kRingNT; t++) {
            const int k = R.T[size_t(t)].localOf(v);
            if (k < 0) continue;
            hs += 0.5*(FE[size_t(t)].H[nFace[size_t(t)] + size_t(k)] + 1.0);
            nt++;
        }
        hs /= double(std::max(1, nt));
        jMax = std::max(jMax, std::fabs(Jv[size_t(v)]));
        cornerMax = std::max(cornerMax, std::fabs(hs - hv[size_t(v)]));
        out("    %-10s %+.4f  %+.5f  %+.5f  %+.5f  %.5f   %.5f (%d tets)\n", kRingVertName[v], phiV[size_t(v)], c[v], Jv[size_t(v)],
            Hv[size_t(v)], hv[size_t(v)], hs, nt);
    }
    out("    J at the vertices is not 0 (max |J(X_v)| = %.3f), so h_v is not (c_v + 1)/2; the export pass next to the corners "
        "differs from h_v by at most %.4f\n", jMax, cornerMax);
    // the range of h, and the size and noise of g
    {
        std::vector<double> hq, gq, seH, seG;
        double hMin = kInf, hMax = -kInf, wOut = 0, wAll = 0;
        long long nOut = 0;
        for (int t = 0; t < kRingNT; t++)
            for (size_t i = 0; i < nQuad[size_t(t)]; i++) {
                const double h = 0.5*(FE[size_t(t)].H[i] + 1.0), w = quad[size_t(t)].pts[i].w;
                hMin = std::min(hMin, h);
                hMax = std::max(hMax, h);
                wAll += w;
                if (h < -0.05 || h > 1.05) {
                    nOut++;
                    wOut += w;
                }
                gq.push_back(0.5*FE[size_t(t)].G[i].norm());
                seH.push_back(0.5*std::sqrt(FE[size_t(t)].vH[i]));
                seG.push_back(0.5*std::sqrt(FE[size_t(t)].vG[i].sum()));
            }
        double rMin = kInf, rMax = -kInf;
        for (double h : rd.h) {
            rMin = std::min(rMin, h);
            rMax = std::max(rMax, h);
        }
        out("  h at the quadrature points: [%.4f, %.4f]; outside [-0.05, 1.05]: %lld points (%.3f%% of the volume); at the render "
            "vertices: [%.4f, %.4f]\n", hMin, hMax, nOut, 100*wOut/wAll, rMin, rMax);
        std::vector<double> g2 = gq, sh = seH, sg = seG;
        out("  |grad h| at the quadrature points: median %.3f, 90%% %.3f, 99%% %.3f, max %.3f; its s.e. (export pass): median %.4f, "
            "99%% %.4f, max %.4f; s.e. of h: median %.5f, max %.5f\n", xfemQuantile(g2, 0.5), xfemQuantile(g2, 0.9),
            xfemQuantile(g2, 0.99), xfemQuantile(g2, 1.0), xfemQuantile(sg, 0.5), xfemQuantile(sg, 0.99), xfemQuantile(sg, 1.0),
            xfemQuantile(sh, 0.5), xfemQuantile(sh, 1.0));
        if (front.size() == 2) {   // h ~ theta/(2 pi) about the front: |grad h| ~ 1/(2 pi r)
            double rMin = kInf, gAt = 0;
            std::vector<double> ratio;
            size_t m = 0;
            for (int t = 0; t < kRingNT; t++)
                for (size_t i = 0; i < nQuad[size_t(t)]; i++, m++) {
                    const double r = pointSegmentDist(quad[size_t(t)].pts[i].X, front[0], front[1]);
                    if (r < rMin) {
                        rMin = r;
                        gAt = gq[m];
                    }
                    if (r < 0.05) ratio.push_back(gq[m]*2.0*std::numbers::pi*r);
                }
            out("  near the front: the closest quadrature point is %.2e from it (|grad h| %.1f there); |grad h| 2 pi r within 0.05 of "
                "it: median %.3f (h ~ theta/(2 pi) + const: 1)\n", rMin, gAt, xfemQuantile(ratio, 0.5));
        }
        if (rc.through) {   // a crack cut through O: the exact ring solution is +-1, h the side indicator
            double e2 = 0, eMax = 0, sw = 0;
            for (int t = 0; t < kRingNT; t++)
                for (size_t i = 0; i < nQuad[size_t(t)]; i++) {
                    const XfemQuadPoint& q = quad[size_t(t)].pts[i];
                    const double e = 0.5*(FE[size_t(t)].H[i] + 1.0) - (q.side > 0 ? 1.0 : 0.0);
                    e2 += q.w*e*e;
                    sw += q.w;
                    eMax = std::max(eMax, std::fabs(e));
                }
            out("  cut through: h against the exact solution (the side indicator) at the quadrature points: RMS %.4f, max %.4f\n",
                std::sqrt(e2/sw), eMax);
        }
    }
    // the export pass against the glue's own field (what the ring view shows) at the same points, with the same c
    {
        double dh2 = 0, dhMax = 0, zh2 = 0, zhMax = 0, dg2 = 0, gg2 = 0, zg2 = 0, vh2 = 0, vg2 = 0;
        long long n = 0;
        std::vector<double> zAbs;
        for (int t = 0; t < kRingNT; t++)
            for (size_t k = 0; k < viewIdx[size_t(t)].size(); k++) {
                const size_t i = nQuad[size_t(t)] + k, j = viewIdx[size_t(t)][k];
                const FieldAt &a = FE[size_t(t)], &b = FV[size_t(t)];
                const double dh = 0.5*(a.H[i] - b.H[j]), vh = 0.25*(a.vH[i] + b.vH[j]);
                dh2 += dh*dh;
                vh2 += vh;
                dhMax = std::max(dhMax, std::fabs(dh));
                const double z = vh > 0 ? dh/std::sqrt(vh) : 0.0;
                zh2 += z*z;
                zhMax = std::max(zhMax, std::fabs(z));
                zAbs.push_back(std::fabs(z));
                const Vec3 dg = 0.5*(a.G[i] - b.G[j]), vg = 0.25*(a.vG[i] + b.vG[j]);
                dg2 += dg.squaredNorm();
                vg2 += vg.sum();
                gg2 += 0.25*b.G[j].squaredNorm();
                for (int q = 0; q < 3; q++) zg2 += vg[q] > 0 ? dg[q]*dg[q]/vg[q]/3.0 : 0.0;
                n++;
            }
        const double nn = double(std::max(1LL, n));
        long long z3 = 0;
        for (double z : zAbs) z3 += z > 3 ? 1 : 0;
        // with 4 Q walks per batch a point far from the crack can see no crack hit in any batch: its jackknife variance is then
        // about 0 and its z meaningless, so the pooled ratio sqrt(sum diff^2 / sum var) is the test (about 1: noise only)
        out("  the export pass against the ring view's glued H (the glue's own walks, Q %d / R %d) at %lld of its points (every "
            "%dth), the same c: h difference RMS %.2e, max %.2e, pooled sqrt(sum d^2 / sum var) %.2f; per point z RMS %.2f, median "
            "|z| %.2f, |z| > 3 at %lld, max |z| %.2f; grad h difference RMS / grad h RMS %.3f, pooled %.2f, z RMS %.2f (the "
            "variances: %d-batch jackknives)\n", fg.nQ, fg.nR, n, kXfemViewStride, std::sqrt(dh2/nn), dhMax,
            std::sqrt(dh2/std::max(1e-300, vh2)), std::sqrt(zh2/nn), xfemQuantile(zAbs, 0.5), z3, zhMax,
            std::sqrt(dg2/std::max(1e-300, gg2)), std::sqrt(dg2/std::max(1e-300, vg2)), std::sqrt(zg2/nn), Bg);
    }
    {
        std::string a, b;
        for (int v = 0; v < kRingNV; v++) {
            a += fmtd(" %+.4f", c[v]);
            b += fmtd(" %+.4f", ex.gl[0].c[v]);
        }
        out("  diagnostic only: the glued c (the glue's points, Q %d / R %d):%s\n                  the export points' own glue "
            "(their rule, Q %d / R %d):%s\n                  (not used: the export rule keeps no eps-tube about the crack, and "
            "|grad P_k|^2 ~ 1/r^2 at the lines where the crack meets the faces)\n", fg.nQ, fg.nR, a.c_str(), fe.nQ, fe.nR,
            b.c_str());
    }

    // 6. the tets and nodes, then the file
    out("\ntets (local corners = global vertex ids, det J > 0):\n");
    for (int t = 0; t < kRingNT; t++) {
        const TetFrame& T = R.T[size_t(t)];
        const XfemTetQuad& q = quad[size_t(t)];
        out("  %s [%d %d %d %d]: class %d (%s), %d + %d sub-tets, %d leaves (%d refined at the front), %zu quadrature points, "
            "%zu check points\n", T.name, T.global[0], T.global[1], T.global[2], T.global[3], cls[size_t(t)],
            cls[size_t(t)] == 0 ? s.rt[size_t(t)].plane ? "plane beyond the front" : "plane misses it" :
            cls[size_t(t)] == 1 ? "cut in two" : "front inside", q.nSub[0], q.nSub[1], q.nLeaves, q.nFrontLeaves,
            nQuad[size_t(t)], pe[size_t(t)].size() - nQuad[size_t(t)]);
    }
    out("node types by the per-node rule (harmonic: a class-2 tet around it; heaviside: else a class-1 tet; else none):");
    for (int v = 0; v < kRingNV; v++) {
        int best = 0;
        for (int t = 0; t < kRingNT; t++)
            if (R.T[size_t(t)].localOf(v) >= 0) best = std::max(best, cls[size_t(t)]);
        out(" v%d %s", v, best == 2 ? "harmonic" : best == 1 ? "heaviside" : "none");
    }
    double areaO = 0, areaS[3] = {0, 0, 0};
    for (const auto& tr : rd.tris)
        areaS[tr[3]] += 0.5*(rd.X[size_t(tr[1])] - rd.X[size_t(tr[0])]).cross(rd.X[size_t(tr[2])] - rd.X[size_t(tr[0])]).norm();
    for (const auto& f : O.tris) areaO += 0.5*(O.v[size_t(f[1])] - O.v[size_t(f[0])]).cross(O.v[size_t(f[2])] - O.v[size_t(f[0])]).norm();
    out("\nrender: %zu vertices, %zu triangles (O's boundary %d per edge, crack spacing %.4f); areas: boundary %.9f (O: %.9f), "
        "crack sheets %.9f and %.9f (O's crack polygon: %.9f)\n", rd.X.size(), rd.tris.size(), kXfemRenderPerEdge, hR, areaS[0],
        areaO, areaS[1], areaS[2], rc.cutO.area());
    if (!front.empty()) {
        double on = 0;
        for (const Vec3& f : front)
            for (size_t k = 0; k < O.n.size(); k++) on = std::max(on, O.n[k].dot(f) - O.c[k]);
        out("front on O: %s to %s (largest n_k.X - c_k of O's planes at its ends: %.1e, on O's boundary)\n", fmt3(front[0]).c_str(),
            fmt3(front[1]).c_str(), on);
    }

    const double secBeforeWrite = secondsSince(tStart);
    char settings[1024];
    std::snprintf(settings, sizeof settings, "our fast method (no GT); glue: %d points per tet asked (%zu), Q %d / R %d walks, %d "
                  "batches, the ring view's seeds; export points: Q %d / R %d walks, %d batches; quadrature: level %d, +1 level "
                  "for leaves within %.2f and +2 within %.3f of the front (leaf to segment), 4-point degree-2 Gauss per leaf, %zu "
                  "points; render: %d per edge; "
                  "time: glue %.2f s, export walks %.2f s, total %.2f s", nPts, totalPoints(s), fg.nQ, fg.nR, Bg, fe.nQ, fe.nR, Be,
                  p.quadLevel, kXfemFrontBand[0], kXfemFrontBand[1], nQuadAll, kXfemRenderPerEdge, secGlue, secExWalks, secBeforeWrite);
    xd = XfemRingData{};
    xd.caseName = rc.name;
    xd.ring = R.id;
    xd.g = 1;
    xd.settings = settings;
    xd.verts.assign(R.V.begin(), R.V.end());
    for (const TetFrame& T : R.T) xd.tets.push_back({T.global[0], T.global[1], T.global[2], T.global[3]});
    xd.n = rc.n;
    xd.d = rc.d;
    xd.front = front;
    xd.tetClass.assign(cls.begin(), cls.end());
    xd.nodeValue.assign(hv.begin(), hv.end());
    for (int t = 0; t < kRingNT; t++)
        for (size_t i = 0; i < nQuad[size_t(t)]; i++) {
            const XfemQuadPoint& q = quad[size_t(t)].pts[i];
            xd.quad.push_back({t, q.side, q.w, q.X, 0.5*(FE[size_t(t)].H[i] + 1.0), 0.5*FE[size_t(t)].G[i]});
        }
    for (size_t i = 0; i < rd.X.size(); i++) xd.renderV.push_back({rd.tet[i], rd.side[i], rd.X[i], rd.h[i]});
    xd.renderT = rd.tris;
    out("\n%zu quadrature points, %zu render vertices, %zu render triangles\n", xd.quad.size(), xd.renderV.size(),
        xd.renderT.size());
    out("time: glue %.3f s, geometry %.3f s, export walks %.3f s, fields %.3f s, total %.3f s\n", secGlue, secGeom, secExWalks,
        secFields, secondsSince(tStart));
    out("self-checks: %s\n", ok ? "PASS" : "FAILED");
    xd.selfChecks = ok;
    return ok;
}

bool runXfemExport(const XfemExportParams& p)
{
    if (p.out.empty()) throw std::invalid_argument("--export-xfem: --out FILE is required");
    Report out(false);
    XfemRingData xd;
    const bool ok = buildXfemRing(p, xd);
    std::string txt;
    appendf(txt, "xfem_ring_export 1\ncase %s\nring %d\ng 1\nsettings %s\n", xd.caseName.c_str(), xd.ring, xd.settings.c_str());
    appendf(txt, "vertices %zu\n", xd.verts.size());
    for (const Vec3& v : xd.verts) appendf(txt, "%.17g %.17g %.17g\n", v.x(), v.y(), v.z());
    appendf(txt, "tets %zu\n", xd.tets.size());
    for (const auto& t : xd.tets) appendf(txt, "%d %d %d %d\n", t[0], t[1], t[2], t[3]);
    appendf(txt, "plane %.17g %.17g %.17g %.17g\n", xd.n.x(), xd.n.y(), xd.n.z(), xd.d);
    appendf(txt, "front %zu\n", xd.front.size());
    for (const Vec3& v : xd.front) appendf(txt, "%.17g %.17g %.17g\n", v.x(), v.y(), v.z());
    appendf(txt, "tet_class %zu\n", xd.tetClass.size());
    for (const int c : xd.tetClass) appendf(txt, "%d\n", c);
    appendf(txt, "node_value %zu\n", xd.nodeValue.size());
    for (const double h : xd.nodeValue) appendf(txt, "%.17g\n", h);
    appendf(txt, "quad %zu\n", xd.quad.size());
    for (const XfemRingData::Quad& q : xd.quad)
        appendf(txt, "%d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", q.tet, q.side, q.w, q.X.x(), q.X.y(), q.X.z(), q.h,
                q.g.x(), q.g.y(), q.g.z());
    appendf(txt, "render_vertices %zu\n", xd.renderV.size());
    for (const XfemRingData::RenderVertex& r : xd.renderV)
        appendf(txt, "%d %d %.17g %.17g %.17g %.17g\n", r.tet, r.side, r.X.x(), r.X.y(), r.X.z(), r.h);
    appendf(txt, "render_triangles %zu\n", xd.renderT.size());
    for (const auto& tr : xd.renderT) appendf(txt, "%d %d %d %d\n", tr[0], tr[1], tr[2], tr[3]);
    appendf(txt, "end\n");
    std::error_code ec;
    const std::filesystem::path path(p.out);
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path(), ec);
    bool written = false;
    {
        std::ofstream f(path, std::ios::binary);
        f << txt;
        written = bool(f);
    }
    out("%s %s\n", written ? "wrote" : "COULD NOT WRITE", p.out.c_str());
    return ok && written;
}
