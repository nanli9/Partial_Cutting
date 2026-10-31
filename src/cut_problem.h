// Harmonic partial-cut enrichment H(X) on one tetrahedron, solved with Walk on Stars
// (Zombie). This header holds everything shared by the GUI, the validation checks and
// the command-line modes:
//   - GeometryConfig : the ONE place to change the tet, the cut plane and the front direction
//   - Tet, Cut       : geometry (plane ∩ tet polygon, clipped at the front, triangulated)
//   - CutSolver      : Zombie WoSt for  ΔH = 0 in the tet, H = ±g on the cut,
//                      ∂H/∂n = 0 on the four tet faces
//   - Lattice, Slice : query points
#pragma once

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using Vec3 = Eigen::Vector3d;   // geometry is done in double; Zombie itself works in float

// =====================================================================================
// THE ONE PLACE TO CHANGE THE EXPERIMENT GEOMETRY
// =====================================================================================
struct GeometryConfig {
    // Any 4 non-coplanar points (vertex order and orientation do not matter).
    std::array<Vec3, 4> tetVerts{Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};

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

// Dirichlet data on the cut: H = +g on the side n·X > d, -g on the other side.
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

Cut buildCut(const Tet& tet, const GeometryConfig& cfg, double t);

// =====================================================================================
// Walk on Stars solver (Zombie)
// =====================================================================================
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
    bool noCut = false;           // t = 0: no Dirichlet boundary, no walks, H := 0
    void add(const SolveStats& o);
    double meanWalkLength() const { return walksCounted ? double(totalSteps)/double(walksCounted) : 0.0; }
    std::string summary() const;
};

class CutSolver {
public:
    // Builds the BVHs and the PDE for one cut. If the cut is empty (t = 0) nothing is built
    // and solve() returns H = 0 everywhere without running walks.
    CutSolver(const Tet& tet, const Cut& cut, GType g, const WalkParams& params);
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
