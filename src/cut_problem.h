// Harmonic partial-cut enrichment H(ξ) on one tetrahedron, in the tet's local coordinates ξ,
// solved with Walk on Stars (Zombie). This header holds everything shared by the GUI, the
// validation checks and the command-line modes:
//   - GeometryConfig : the ONE place to change the rest tet, the cut plane and the front direction
//   - Tet, Cut       : geometry (plane ∩ tet polygon, clipped at the front, triangulated)
//   - CutSolver      : Zombie WoSt for  Δ_ξ H = 0 in the reference tet, H = ±g on the cut,
//                      ∂H/∂n = 0 on the four tet faces (the coupling gate may flag faces Dirichlet:
//                      WalkParams::faceBC)
// Every solve runs on the reference tet (referenceTet()) in ξ, never in rest or deformed space;
// a rest tet only maps the cut into ξ (pullBackCut).
//   - Lattice, Slice : query points
#pragma once

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using Vec3 = Eigen::Vector3d;   // geometry is done in double; Zombie itself works in float

// =====================================================================================
// THE ONE PLACE TO CHANGE THE EXPERIMENT GEOMETRY
// =====================================================================================
struct GeometryConfig {
    // The tet in rest space. It is NOT the solver domain: H is solved in local coordinates ξ on
    // the reference tet (0,0,0), (1,0,0), (0,1,0), (0,0,1), and restVerts only pulls the cut
    // below back into ξ through X = X0 + J ξ, J = [X1 - X0, X2 - X0, X3 - X0] (pullBackCut).
    // The vertex order is fixed: X0 is ξ = 0 and Xi is ξ = e_i. Any 4 non-coplanar points.
    // The default is the reference tet itself, so there X = ξ and the cut is used as given.
    std::array<Vec3, 4> restVerts{Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};

    // The cut plane and front direction below are in rest space (on restVerts). The numbers in
    // these comments are for the default rest tet, where rest space and ξ coincide.
    //
    // Cut plane {X : n·X = d}. planeNormal is normalized when the cut is built and
    // planeOffset is the offset along that UNIT normal. Default: x + y - z = 0.5, i.e.
    // n = (1,1,-1)/sqrt(3), d = 0.5/sqrt(3). It separates {v0, v3} from {v1, v2}, so it
    // cuts the default tet in a quad (checked at startup).
    Vec3 planeNormal{1, 1, -1};
    double planeOffset = 0.5/std::sqrt(3.0);

    // Direction in which the crack front advances (projected into the plane, normalized).
    // The cut at parameter t is the part of the plane∩tet polygon with X·u <= s(t),
    // s(t) = s_min + t (s_max - s_min) over the polygon's extent along u.
    // With (1,-1,0) the section quad has the vertices D (0,.75,.25) on edge v2-v3 at t = 0,
    // B (0,.5,0) on v0-v2 at t = 1/6, A (.5,0,0) on v0-v1 at t = 5/6 and C (.75,0,.25) on
    // v1-v3 at t = 1. The cut grows from the single vertex D towards the vertex C. The front
    // line (direction n×u) ends on the face x+y+z=1 and on the face x=0 for t < 1/6, z=0 for
    // 1/6 < t < 5/6, or y=0 for t > 5/6. At t = 1 it only touches the polygon, at C.
    Vec3 frontDir{1, -1, 0};
};

// Dirichlet data on the cut: H = +g on the side n·ξ > d, -g on the other side (the same side as n_X·X > d_X).
enum class GType : int {
    One = 0,        // (a) g = 1 (Kaufmann's simplification)
    SqrtFront = 1   // (b) g = min(1, sqrt(r/L)), r = distance to the front line, L = longest tet edge
};
const char* gTypeName(GType g);

// =====================================================================================
// Geometry
// =====================================================================================
struct Tet {
    std::array<Vec3, 4> v;
    std::array<std::array<int, 3>, 4> faces;   // outward-oriented: ((b-a)x(c-a)) points out
    Vec3 centroid;
    double L = 0;                              // longest edge
    Eigen::Matrix3d Einv;                      // x -> barycentric coordinates (l1, l2, l3)

    explicit Tet(const std::array<Vec3, 4>& verts);
    Eigen::Vector4d bary(const Vec3& x) const;          // (l0, l1, l2, l3), sum = 1
    bool contains(const Vec3& x, double tol = 0) const; // all barycentric coords >= -tol
    double distToFaces(const Vec3& x) const;            // min distance to the 4 face planes (x inside)
};

// A convex polyhedron with every face Neumann: the domain of ONE solve only, the coupling gate's ring
// reference (check 4: the octahedron, the union of a node's 8 tets, in physical X). Every per-tet solve stays
// on the reference tet in ξ. Read by CutSolver (Zombie) and the own loop's solveMulti (walk_cache.h).
constexpr int kMaxDomainPlanes = 8;
struct ConvexDomain {
    std::vector<Vec3> v;                       // vertices
    std::vector<std::array<int, 3>> tris;      // boundary triangles, outward-wound ((b-a)x(c-a) points out)
    std::vector<Vec3> n;                       // face planes n·x <= c inside, n unit and outward
    std::vector<double> c;
    double L = 0;                              // diameter: the scale of the walk shells (WalkParams *Rel)
};

// The solver domain: the reference tet in ξ. Its longest edge L = √2 is the L of g type (b), and
// epsilon, the lattice, the slice and ∇H (= ∇_ξ H) are all in ξ as well.
Tet referenceTet();
bool isReferenceTet(const Tet& tet);   // exactly the reference vertices, in order

// The cut in ξ: plane n·ξ = d (n unit) and the front direction f, not yet projected into the
// plane. pullBackCut maps the rest-space cut of cfg through X = X0 + J ξ:
//   plane:  n_ξ = Jᵀ n_X / |Jᵀ n_X|,  d_ξ = (d_X - n_X·X0) / |Jᵀ n_X|  (n_X unit; the same side is "+g");
//   front:  f_ξ = Jᵀ f_X. On the plane X·u_X is an increasing affine function of u_ξ·ξ
//           (u = f projected into the plane, normalized), so the front X·u_X = s_X(t) maps to
//           the front u_ξ·ξ = s_ξ(t) with the same t.
// For the default rest tet J = I and the result is the rest-space cut bit for bit.
struct XiCut {
    Vec3 n = Vec3::Zero();
    double d = 0;
    Vec3 f = Vec3::Zero();
};
XiCut pullBackCut(const GeometryConfig& cfg);

// Convex polygon {X : n·X = d} ∩ tet, ordered counter-clockwise about n (seen from +n).
// 'tangent' only fixes where the angular ordering starts. Empty if the plane misses the tet.
std::vector<Vec3> planeTetSection(const Tet& tet, const Vec3& n, double d, const Vec3& tangent,
                                  int* nVertsOnPlane = nullptr);

struct Cut {
    Vec3 n, u, w;              // unit plane normal, front-advance direction, front-line direction w = n×u
    double d = 0;              // plane: n·X = d
    double L = 0;              // tet longest edge (scale of g type b)
    double t = 0;              // front parameter in [0,1]
    double sMin = 0, sMax = 0; // extent of the full section along u
    double s = 0;              // front position s(t)
    std::vector<Vec3> fullPolygon;          // plane ∩ tet (t = 1)
    bool fullIsQuad = false;                // 4 edge crossings and no tet vertex on the plane
    std::vector<Vec3> polygon;              // the cut: part of fullPolygon with X·u <= s, CCW about n
    std::vector<std::array<int, 3>> tris;   // fan triangulation of polygon, each with ((b-a)x(c-a))·n > 0

    bool empty() const { return tris.empty(); }
    double side(const Vec3& x) const { return n.dot(x) - d; }   // > 0: the "+g" side
    Vec3 frontLinePoint() const { return d*n + s*u; }            // a point on the front line
    double distToFrontLine(const Vec3& x) const;
    double g(const Vec3& x, GType type) const;
    double distToCut(const Vec3& x) const;          // exact 3D distance to the cut polygon
    double distToPolygonEdges(const Vec3& x) const; // in-plane distance from x to the polygon boundary
    bool polygonContains(const Vec3& x) const;      // projection of x lies inside the cut polygon
    Vec3 frontMidpoint() const;                     // middle of the front segment (polygon centroid if none)
    double area() const;
};

// The cut at front parameter t, in ξ: cfg's cut pulled back (pullBackCut) and built on tet,
// which must be referenceTet() (throws otherwise).
Cut buildCut(const Tet& tet, const GeometryConfig& cfg, double t);
// Geometry only, on any tet: the cut of the plane n·X = d (n unit) with front direction f. buildCut uses it
// on the reference tet; the ξ check (--check-xi) also builds the rest-space cut with it, to compare
// geometry. Never solve on a cut built in rest space.
Cut buildCutOnTet(const Tet& tet, const Vec3& n, double d, const Vec3& f, double t);
// The same on a convex polyhedron (the coupling gate's ring reference, check 4: the crack in physical X on the
// octahedron), with the length scale L of g given (the ring uses √2, the reference tet's L): the section is the
// plane clipped by the domain's half-spaces, CCW about n, then cut at the front exactly as buildCutOnTet does.
Cut buildCutOnDomain(const ConvexDomain& dom, const Vec3& n, double d, const Vec3& f, double L, double t);

// =====================================================================================
// Walk on Stars solver (Zombie)
// =====================================================================================
// Boundary condition on one tet face. Everything except the coupling gate (--gate coupling) keeps all
// four faces Neumann (reflecting).
enum class FaceBC : uint8_t { Neumann = 0, Dirichlet = 1 };

struct WalkParams {
    // Walk termination shell around the cut, as a fraction of L. It is deliberately a
    // fixed, small multiple of L and much smaller than the user epsilon (the UI keeps
    // epsilon >= 10x this), so the jump check at X ± epsilon·n is not trivially exact.
    // The bias it introduces is O(shell). Walks stop within shell + faceOffsetMargin(tet)
    // of the cut (CutSolver::stopDistance()).
    double absorbShellRel = 1e-4;
    // Minimum star radius (Zombie's epsilonShellForReflectingBoundary), × L. Kept below the
    // absorbing shell: Zombie shrinks the star radius to 0.99x the distance to the cut only
    // while that distance is >= this value; below it the radius equals the distance. Before
    // the face-offset margin existed, 1e-3 put 91 of 1.72M walks on the wrong side (t = 1,
    // g = 1, N = 4096); with the margin, 0 of 6.89M at 1e-3 and 0 of 13.8M at 1e-5.
    double reflectShellRel = 1e-5;
    // Walk length cap. Capped walks still count, with terminal value 0 (biases H towards 0),
    // so they are counted and reported with every solve.
    int maxWalkLength = 4096;
    // Experiment only (--dirichlet-faces), read by the own loop of walk_cache.h; Zombie ignores it:
    // the four tet faces absorb with H = 0 instead of reflecting. A different PDE, used to test what
    // limits walk reuse; off everywhere else.
    bool dirichletFaces = false;
    // Coupling gate only (--gate coupling); the defaults are the first gate's problem, bit for bit. Read
    // by Zombie (CutSolver) and by the own loop's solveMulti (walk_cache.h) only: WalkCache, solveSums
    // and solveSumsOnly throw if any of them is set (usesGateFlags):
    //  - faceBC[f]: tet face f (Tet::faces[f], the face opposite vertex f) reflects (Neumann, the
    //    default) or absorbs (Dirichlet) with H = faceData(f, x), where x is the point of that face's
    //    plane where the walk ends, in the solver's frame (xi for the per-tet solves). No faceData: H = 0.
    //    dirichletFaces above wins over these: all four faces, H = 0.
    //  - cutData = false: H = 0 on the cut instead of +-g (the gate's response functions P_k).
    std::array<FaceBC, 4> faceBC{FaceBC::Neumann, FaceBC::Neumann, FaceBC::Neumann, FaceBC::Neumann};
    std::function<double(int face, const Vec3& x)> faceData;
    bool cutData = true;
    // Coupling gate only, own loop only (Zombie ignores both): the gradient's first ball may cross the
    // nearest face into the mirror image of the tet, so next to a face its radius is set by the other
    // boundaries instead of shrinking to the face distance. See runWalk in walk_cache.cpp.
    //  - mirrorGradBall: across a Neumann face (the even extension of H across a zero-flux plane is harmonic);
    //  - oddMirrorDirichlet (with mirrorGradBall): also across a Dirichlet face. H - Phi, Phi a harmonic
    //    function equal to the face data there, vanishes on the face, so its odd extension is harmonic;
    //    solveMulti then needs every datum's Phi (its HarmonicExtension argument).
    bool mirrorGradBall = false;
    bool oddMirrorDirichlet = false;
    bool anyDirichletFace() const
    {
        for (FaceBC b : faceBC)
            if (b == FaceBC::Dirichlet) return true;
        return false;
    }
    // any coupling-gate flag set (dirichletFaces, the experiment, is not one of them)
    bool usesGateFlags() const { return anyDirichletFace() || bool(faceData) || !cutData || mirrorGradBall || oddMirrorDirichlet; }
};

struct Estimate {
    float H = 0, Hse = 0;                          // mean and standard error (1 s.e.)
    Eigen::Vector3f grad = Eigen::Vector3f::Zero();
    Eigen::Vector3f gradSe = Eigen::Vector3f::Zero();
    bool hasGrad = false;
    int nCounted = 0;                              // walks that entered the estimate
    float walkVar = 0;                             // per-walk sample variance of H
};

struct SolveStats {
    int nPoints = 0;
    long long walksRequested = 0, walksCounted = 0;
    long long walksEscaped = 0;   // left the tet through float round-off: discarded by Zombie
    long long walksCapped = 0;    // hit maxWalkLength: counted with terminal value 0
    long long totalSteps = 0;
    double seconds = 0;
    bool noCut = false;           // t = 0 and no Dirichlet face: no Dirichlet boundary, no walks, H := 0
    void add(const SolveStats& o);
    double meanWalkLength() const { return walksCounted ? double(totalSteps)/double(walksCounted) : 0.0; }
    std::string summary() const;
};

class CutSolver {
public:
    // Builds the BVHs and the PDE for one cut. If the cut is empty (t = 0) and no face is Dirichlet
    // (WalkParams::faceBC) nothing is built and solve() returns H = 0 everywhere without running walks.
    // Dirichlet faces join the cut in the absorbing boundary; the rest of the faces reflect.
    CutSolver(const Tet& tet, const Cut& cut, GType g, const WalkParams& params);
    // The coupling gate's ring reference only (check 4): the convex polyhedron dom, every face Neumann (throws if
    // params sets any coupling-gate flag), the cut ±g as above. Escape test: exact containment in dom's planes.
    CutSolver(const ConvexDomain& dom, const Cut& cut, GType g, const WalkParams& params);
    // The ring stress test only (--ring-stress): the convex polyhedron dom, every face Neumann (as above), and the crack(s)
    // as one triangle mesh (V, F; every triangle wound towards its + side) with the Dirichlet data dirichlet(x,
    // normalAligned), x = the walk's end projected onto the mesh, normalAligned = it ended on the nearest triangle's + side.
    CutSolver(const ConvexDomain& dom, const std::vector<Vec3>& V, const std::vector<std::array<int, 3>>& F,
              std::function<float(const Vec3&, bool)> dirichlet, const WalkParams& params);
    ~CutSolver();
    CutSolver(const CutSolver&) = delete;
    CutSolver& operator=(const CutSolver&) = delete;

    // Estimates H at pts with nWalks walks each (and ∇H with Zombie's gradient estimator
    // where wantGrad[i] != 0). Point i uses its own RNG stream pcg32(mix(seed), i), so
    // results are reproducible for a given seed and independent across points and seeds.
    SolveStats solve(const std::vector<Vec3>& pts, const std::vector<char>* wantGrad, int nWalks,
                     uint64_t seed, std::vector<Estimate>& out) const;

    // walks stop at a true distance <= this from the cut: shell + faceOffsetMargin(tet)
    double stopDistance() const { return stopDistance_; }

private:
    struct Impl;                  // holds the Zombie objects (never copied or moved)
    std::unique_ptr<Impl> impl_;
    double stopDistance_ = 0;
};

// Zombie starts every step from a tet face at a point pushed off the face (by up to 2^-16,
// or 256 float ulps, per coordinate), but sizes the step's ball as 0.99x the distance to
// the cut measured at the un-pushed point. Close to the cut (within ~6e-4 L for the default
// tet) the push is larger than that 1% margin, the ball can reach across the cut and the
// walk ends on the wrong side. CutSolver therefore gives Zombie distances to the cut that
// are smaller by this margin: twice the largest push for this tet's coordinates.
double faceOffsetMargin(const Tet& tet);
double faceOffsetMargin(const std::vector<Vec3>& verts);   // the same for the convex hull of verts

uint64_t mixSeed(uint64_t a, uint64_t b);   // splitmix64-style combination, for deriving seeds

// =====================================================================================
// Query points
// =====================================================================================
// Every query point (lattice and slice) is pulled this fraction of the way towards the
// centroid before solving (a uniform 0.1% shrink), so that walks start strictly inside the
// domain: points exactly on a slanted face lie outside after float round-off and lose all
// their walks (stage 1).
constexpr double kBoundaryPull = 1e-3;

struct Lattice {
    int k = 0;
    std::vector<Vec3> pts;                   // kept points (already pulled inside)
    std::vector<std::array<int, 4>> index;   // barycentric indices (i0..i3), sum = k
    std::vector<char> interior;              // all indices >= 1 (not on a tet face)
    int nTotal = 0, nExcluded = 0;           // lattice size, points within epsilon of the cut
};
// Barycentric lattice with k subdivisions, excluding points within eps of the cut.
Lattice buildLattice(const Tet& tet, const Cut& cut, int k, double eps);

struct SliceParams {
    double offset = 0;      // signed offset of the plane from the tet centroid along its normal
    double yawDeg = 0;      // rotation of the normal about the cut normal n (towards u)
    double pitchDeg = 0;    // tilt of the normal towards n
    int res = 64;           // grid cells across the section's bounding box
};
struct Slice {
    Vec3 normal = Vec3::Zero(), origin = Vec3::Zero();
    std::vector<Vec3> verts;                 // inside the tet (pulled like the lattice)
    std::vector<std::array<int, 3>> tris;
};
// Slice plane through centroid + offset·m, m = normal from the two angles; the default
// (0, 0, 0) has normal n×u, so the cut shows up as a segment that ends at the front.
Slice buildSlice(const Tet& tet, const Cut& cut, const SliceParams& sp);

// =====================================================================================
// One full solve of the scene (what the GUI shows; --bench times exactly this)
// =====================================================================================
struct Settings {
    double t = 0.6;            // front position
    int nWalks = 256;          // walks per query point
    double epsRel = 0.01;      // user epsilon / L: lattice exclusion radius, jump-check half separation
    GType g = GType::One;
    int k = 12;                // lattice subdivisions
    SliceParams slice;
    WalkParams walk;
    uint64_t seed = 20261031;  // base seed
};

struct SceneResult {
    Cut cut;
    Lattice lattice;
    Slice slice;
    std::vector<Estimate> latticeEst, sliceEst;
    SolveStats latticeStats, sliceStats;
    double buildSeconds = 0, totalSeconds = 0;
    uint64_t seed = 0;
};
// Builds cut, lattice and slice and estimates H on both (and ∇H on the interior lattice
// points). Everything is rebuilt from scratch; nothing is reused between calls.
SceneResult solveScene(const Tet& tet, const GeometryConfig& cfg, const Settings& s, uint64_t seed);
