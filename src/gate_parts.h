#pragma once

// The pieces of the coupling gate (coupling_gate.cpp) that xfem_mesh.cpp reuses on a whole tet mesh: a tet's xi frame, the
// crack in that frame, the point sets, the fast path's control-variate runs with their Gram and glued field, and the XFEM
// export's quadrature and render data. The gate's file-local code lives in namespace gate (it was an unnamed namespace;
// coupling_gate.cpp opens it with a using-directive, so every name there resolves as before); the structs are defined
// here, the functions in coupling_gate.cpp. coupling_gate.h stays the gate's interface; this header is internal to
// partial_cut_core.

#include "walk_cache.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numbers>
#include <string>
#include <vector>

namespace gate {

constexpr double kSqrt2 = 1.41421356237309504880;

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
TetFrame frameOn(const char* name, const std::array<int, 4>& global, const Vec3* verts);

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
std::vector<Vec3> frontVertices(const Cut& c);

struct Leaf {
    std::array<Vec3, 4> v;
    int depth = 0;
};
double leafVolume(const Leaf& l);
// Bey's red refinement: 8 children of exactly 1/8 of the volume
std::array<Leaf, 8> beyChildren(const Leaf& L);
double pointSegmentDist(const Vec3& p, const Vec3& a, const Vec3& b);
double pointTriangleDist(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c);

// One fixed, deterministic point set in one tet (built in xi on the reference tet): coupling_gate.cpp, volumePoints.
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

// The own loop's walk parameters in the gate (coupling_gate.cpp, gateWalkParams): the termination shell 1e-6 L, the
// mirror balls of the gradient's first step.
constexpr double kGateShellRel = 1e-6;
WalkParams gateWalkParams(const std::array<FaceBC, 4>& bc, bool mirror, bool oddMirror);
std::array<FaceBC, 4> allDirichlet();

// Gram matrix G_ql = integral over the kept points of grad f_q . grad f_l (physical), from B independent batches of
// per-point gradient estimates (coupling_gate.cpp, gramFromBatches): products of different batches only, delete-one-batch
// jackknife replicates.
struct Gram {
    int Q = 0, B = 0;
    Eigen::MatrixXd G, se, naive, halves;
    std::vector<Eigen::MatrixXd> jack;   // delete-one-batch replicates
};

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
// The crack n.X = d (n unit, + side n.X > d), front direction f, front line u.X = s (u: f projected into the plane,
// normalized) in tet T: ringTetCrack's body, which calls it with O's front.
RingTetCrack tetCrackAt(const TetFrame& T, const Vec3& n, double d, const Vec3& f, double s, const Vec3& u);
// check 4's point set of one tet: Bey leaves in xi, refined around rt's crack-face lines and front, stratified Halton
// points; the points within eps (physical) of cutO dropped
PointSet ringPoints(const TetFrame& T, const RingTetCrack& rt, const Cut& cutO, int target, double eps, int mMin);

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
// Level C's theta/pi of a half-plane (plane n.X = d, front line u.X = sJ) in a tet's xi frame (xiFront with explicit data)
XiFront xiFrontOf(const TetFrame& T, const Vec3& nX, double d, const Vec3& uX, double sJ);

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
// The Gram of one tet in check 4's layout (H0 g1, H0 g2, P_0..P_3) from the control-variate runs, physical energy
Gram cvGram(const CvRun& c, const PointSet& ps, const TetFrame& T, int B);
// The glued field of one tet from the control-variate runs: H = J_xi + R + sum_k c_k (lambda_k + Q_k), physical gradient,
// delete-one-batch jackknife variances with the replicates' coefficients cJ (empty: cv fixed)
FieldAt cvField(const CvRun& c, const PointSet& ps, const TetFrame& T, int m, const Eigen::Vector4d& cv,
                const std::vector<Eigen::Vector4d>& cJ, int B);
// the fast path's batch count for n walks per point: 8 up to 32 walks, 64 from 1280, 16 between
int defaultBatches(int n);

// One tet's XFEM quadrature rule in X (coupling_gate.cpp, xfemTetQuad): the tet split by the level set with corner values
// phi into the sub-tets of its + and - pieces, refined, and refined again near the front segments; 4 Gauss points per leaf.
struct XfemQuadPoint {
    Vec3 X = Vec3::Zero();
    double w = 0;
    int side = 0;   // the piece: +1 (n.X > d) or -1
};
struct XfemTetQuad {
    std::vector<XfemQuadPoint> pts;
    int nSub[2] = {0, 0};                // sub-tets of the + and the - piece
    int nLeaves = 0, nFrontLeaves = 0;   // leaves; of them, refined for the front
    int nNear = 0, nNearShallow = 0;     // leaves within band[1] of the front; of them, not at depth level + 2
    double nearDiam = 0;                 // the largest diameter among those leaves
};
// band: the two front bands (a leaf closer than band[e] to a front segment is refined once more, e = 0, 1)
XfemTetQuad xfemTetQuad(const TetFrame& T, const Eigen::Vector4d& phi, int level,
                        const std::vector<std::array<Vec3, 2>>& front, const double band[2]);

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
// A boundary triangle of tet t (corners P, wound outwards) cut into n x n small triangles, each clipped by the plane
// nP.X = dP into its parts on either side (group 0)
void xfemRenderFace(XfemRender& r, int t, const std::array<Vec3, 3>& P, const Vec3& nP, double dP, int n);
// The crack inside tet t (its polygon P in X, CCW about nrm), fanned and cut at spacing hR, written twice: group 1 (side
// +1) and group 2 (side -1). Returns the fan triangles whose winding disagreed with nrm (0 expected).
int xfemRenderCrack(XfemRender& r, int t, const std::vector<Vec3>& P, const Vec3& nrm, double hR);

} // namespace gate
