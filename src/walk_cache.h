// Incremental Walk on Stars under a progressively growing cut, after Saeed, Lambiri, Schneider,
// Nowrouzezahrai, Kry, "Fast Cutting of Cage Based Deformation with Walk on Spheres" (SIGGRAPH 2026),
// on the single-tet gate of cut_problem.h. This header holds:
//   - StripCut  : the cut as an ordered list of strips; one advance of the front appends one strip,
//                 which is the whole CutDelta of that step (old triangles never change)
//   - the WoSt step function: our own double-precision loop, replayable walk by walk from
//                 (point seed, walk index); exact brute-force geometric queries
//   - WalkCache : per-walk cache (bounding box of the step balls + landing data) and per-point
//                 running sums; full solve, incremental update (bbox test, replay + resume, value
//                 updates), and the deliberately biased restart variant (paper Fig. 5)
//   - solveMulti: one set of walks scored against several boundary data at once, in independent batches
//                 (the coupling gate, coupling_gate.h); the step function also takes the gate's per-face
//                 Dirichlet flags and its mirror balls for the gradient (WalkParams), all off elsewhere.
//                 Only solveMulti supports them: WalkCache, solveSums and solveSumsOnly throw if one is set
//   - the reuse options compared by --compare-options: 1 = solveSums (frozen seeds, full recompute),
//                 2a/2b = WalkCache with every step stored (StepStore; exact ball test, resume without
//                 replay; 2b pops walks from a queue keyed by the front position that first enters
//                 one of their balls), 3 = WalkCache::advance(Incremental) (seed replay), and the
//                 control 3q = seed replay with 2b's queue (StepStore::Keys: keys without stored steps)
//
// Why not Zombie's WalkOnStars: all walks of a sample point draw from one shared pcg32
// (samplePt.rng), so walk j's random numbers depend on walks 0..j-1; WalkOnStars::walk() runs a
// walk to termination with no way to stop and resume it at a given step (walkStateCallback only
// observes); the gradient path draws stratified first steps for all walks up front and copies the
// RNG between antithetic pairs; the number of draws per step depends on settings; and it works in
// float with a 0.99 radius shrink and a face push. CutSolver (Zombie) stays the solver of the
// first gate; this module is separate.
#pragma once

#include "cut_problem.h"

#include <pcg32.h>   // zombie's deps/pcg32

#include <functional>
#include <limits>
#include <memory>
#include <utility>

using Vec2 = Eigen::Vector2d;

// =====================================================================================
// Cut geometry: an ordered list of strips
// =====================================================================================
// Strip m = fullPolygon ∩ {s_(m-1) <= X·u <= s_m} (s_(-1) = sMin), fan-triangulated, with its
// triangles appended to one global list (global index = append order). The cut "at version m" is
// strips 0..m-1, so the geometry before an advance is a prefix of the geometry after it and stays
// alive by construction. Everything is done in the plane's (u, w) coordinates: a point X of the
// plane is d n + (X·u) u + (X·w) w, and (u, w, n) is right-handed, so CCW in (u, w) is CCW about n
// (the existing winding rule ((b-a)x(c-a))·n > 0). The two vertices on a front line X·u = s are
// computed by one function of s, so adjacent strips share them bit for bit.
class StripCut {
public:
    struct Tri {
        Vec2 v[3];           // CCW in (u, w)
        Vec2 e[3];           // edges v[i+1] - v[i]
        double invLen2[3];   // 1/|e_i|^2
        int strip = 0;
    };
    struct Strip {
        double s0 = 0, s1 = 0;       // extent along u (exact: the front positions)
        double wMin = 0, wMax = 0;   // extent along w
        Vec3 lo = Vec3::Zero(), hi = Vec3::Zero();   // 3D bounding box of its vertices
        double area = 0;
        int triBegin = 0, triEnd = 0;
    };
    // Closest point of the cut. The key is the squared distance d2 = h^2 + e2 (h = n·x - d, e2 = in-plane
    // squared distance to the triangle). Triangles are compared lexicographically by (d2, index), so the
    // result does not depend on evaluation order, and a tie keeps the older (lower-index) triangle.
    struct Closest {
        double d2 = std::numeric_limits<double>::infinity();
        int tri = -1;
        Vec2 q = Vec2::Zero();   // closest point in (u, w)
    };

    StripCut() = default;
    StripCut(const Tet& tet, const GeometryConfig& cfg);   // no strips yet (t = 0)
    // Coupling gate only: the physical ground truth on the union U = A u B (coupling_gate.h), the one solve
    // that is not in a tet's xi frame. frame is buildCutOnTet on U at t = 1 (its L set by the caller); no
    // strips yet. The per-tet solves stay in xi and use the constructor above (buildCut keeps its guard).
    // Check 4's ring reference on the octahedron uses it too, with frame = buildCutOnDomain at t = 1.
    static StripCut physicalUnionReference(const Cut& frame);

    const Cut& frame() const { return frame_; }   // plane, u, w, L, sMin, sMax, fullPolygon (built at t = 1)
    double sOf(double t) const;
    double tOf(double s) const;
    int nStrips() const { return int(strips_.size()); }
    int nTris() const { return int(tris_.size()); }
    const Strip& strip(int m) const { return strips_[m]; }
    double front(int nStrips) const { return nStrips > 0 ? strips_[nStrips - 1].s1 : frame_.sMin; }
    double front() const { return front(nStrips()); }
    // Appends the strip [front(), sNew] and returns its index (-1 and no change if sNew <= front()).
    int addStrip(double sNew);
    // Retraction fallback: keeps the strips that end at or before sNew, then appends [front(), sNew].
    void retractTo(double sNew);
    double area(int nStrips) const;
    double fullArea() const { return fullArea_; }
    // Copy of the plane frame with the front at front(nStrips): g() of cut_problem.h on it is the
    // Dirichlet data of that geometry version.
    Cut frontCut(int nStrips) const;
    Vec3 toWorld(const Vec2& q) const { return frame_.d*frame_.n + q.x()*frame_.u + q.y()*frame_.w; }

    // Closest point over strips [0, nStrips). Strips are skipped when a safe lower bound on their
    // distance exceeds the best so far; with the (d2, index) rule that never changes the result.
    Closest closest(const Vec3& x, int nStrips) const;
    Closest closestBruteForce(const Vec3& x, int nStrips) const;   // no skipping, index order (for checks)
    double stripDist2(const Vec3& x, int m) const;                 // same key, strip m only
    // Conservative box-vs-strip test (separating axes x, y, z, n, u, w, padded): false only if the box
    // cannot meet the strip. False positives are fine.
    bool boxMayHitStrip(const float lo[3], const float hi[3], int m) const;
    // Display mesh of strips [0, nStrips) in 3D (shared front vertices are duplicated per strip).
    void mesh(int nStrips, std::vector<Vec3>& V, std::vector<std::array<int, 3>>& F, std::vector<double>& triStrip) const;

    // Reuse option 2b (the future cut surface is known): the part of the section the front has not
    // reached yet, P ∩ {X·u >= front(nStrips)}, CCW in (u, w). Every strip added later lies in it. Q[0]
    // is the lower front vertex, so Q[0].x() is the current front.
    std::vector<Vec2> futurePolygon(int nStrips) const;
    // Reuse option 2b: a lower bound s* on the front position at which a strip added inside
    // Q = futurePolygon(...) can first pass the update's exact ball test stripDist2(x) < key2: the minimum
    // of X·u over Q ∩ (the plane's disk of ball(x, sqrt(key2)), enlarged), minus a tolerance; +inf if the
    // ball cannot reach Q. Returns +inf as well when s* cannot be below 'bound' (early out for a running min).
    double entryFront(const std::vector<Vec2>& Q, const Vec3& x, double key2,
                      double bound = std::numeric_limits<double>::infinity()) const;

private:
    void initFromFrame();   // chains and area from frame_ (both constructors)
    void frontVertices(double s, Vec2& lo, Vec2& hi) const;
    double triDist2(const Tri& T, const Vec2& q, Vec2& c) const;
    void evalStrip(int m, double h2, const Vec2& q, Closest& best) const;
    double stripLowerBound2(int m, double h2, const Vec2& q) const;

    Cut frame_;
    std::vector<Vec2> lower_, upper_;   // lower and upper chains of the section, by increasing u
    double fullArea_ = 0;
    std::vector<Tri> tris_;
    std::vector<Strip> strips_;
};

// =====================================================================================
// Walks
// =====================================================================================
// Ratio of the gradient points' first ball to min(dist to cut, dist to faces), as in Zombie's
// estimateSolutionAndGradient (ZOMBIE_RADIUS_SHRINK_PERCENTAGE): the first step then lands
// strictly inside the tet, never on a face.
constexpr double kGradFirstBallShrink = 0.99;

// The complete state of a walk before its step 'steps'. Replaying a walk to step k reproduces
// this state bit for bit; resuming = continuing the walk from it under another geometry.
struct WalkState {
    Vec3 x = Vec3::Zero();
    int onFace = -1;            // tet face that x lies on, -1 inside
    int steps = 0;              // step counter k (moves made so far)
    pcg32 rng;                  // after exactly rngOffset + 2k draws (2 per move, none at the end)
    Vec3 w = Vec3::Zero();      // gradient weight 3 dir0/R0 of the first step (gradient points)
    float lo[3], hi[3];         // union of the step balls so far (rounded outward)
};

enum class WalkStatus : uint8_t { Landed = 0, Capped = 1, Escaped = 2 };

// Per-walk cache entry: the paper's bounding volume + landing data, plus the step count (only
// needed for the replay checks) and rngOffset (biased-restart mode only, 0 otherwise).
struct WalkEntry {
    float lo[3], hi[3];        // bounding box of all step balls ball(c_k, d_k), clipped to the tet's box
    Vec3 land = Vec3::Zero();  // landing point on the cut (Landed), else the last position
    double value = 0;          // side·g(land) for Landed, 0 for Capped/Escaped
    int32_t steps = 0;
    int8_t side = 0;           // sign(n·x - d) at the walk's final position (Landed)
    WalkStatus status = WalkStatus::Landed;
    int8_t face = -1;          // Landed on a Dirichlet face (WalkParams): that face; -1 = on the cut (fits the padding)
    uint32_t rngOffset = 0;    // draws consumed before this path's step 0 (biased restarts only)
};

// One seed per query point is all the per-point state besides the sums.
struct QueryPoints {
    std::vector<Vec3> x;
    std::vector<char> grad;       // estimate ∇H here (first step is a full ball, see kGradFirstBallShrink)
    std::vector<uint64_t> seed;   // walk j of point i uses pcg32(mixSeed(seed[i], j), j)
    // appends pts; point i of this group gets seed mixSeed(seedBase, i)
    void add(const std::vector<Vec3>& pts, const std::vector<char>* gradFlags, uint64_t seedBase);
    size_t size() const { return x.size(); }
};

// Running sums of one query point (double), so a walk's contribution can be subtracted without
// re-running it. u = walk value, w = first-step gradient weight (zero at non-gradient points).
struct PointSums {
    long long n = 0;              // counted walks (landed + capped)
    double sH = 0, sH2 = 0;
    double sw[3] = {0, 0, 0}, swu[3] = {0, 0, 0};                         // Σw, Σw·u
    double sw2[3] = {0, 0, 0}, sw2u[3] = {0, 0, 0}, sw2u2[3] = {0, 0, 0};  // Σw², Σw²u, Σw²u² (s.e.)
    void add(double u, const Vec3& w, double sign);   // sign +1 adds a walk, -1 removes it
    void add(const PointSums& o);
};

// Estimate in double (R2 compares them to ~1e-15).
struct EstimateD {
    double H = 0, Hse = 0, walkVar = 0;
    Vec3 grad = Vec3::Zero(), gradSe = Vec3::Zero();
    bool hasGrad = false;
    long long n = 0;
    Estimate toFloat() const;
};
// H = ΣH/n. ∇H = mean of w (u - c) with control variate c = H read from the same sums:
// (Σwu - H Σw)/n; its s.e. from Σ(w(u - H))² = Σw²u² - 2H Σw²u + H² Σw².
EstimateD estimateFromSums(const PointSums& s, bool grad);

// The reuse options of --compare-options: option 1 = solveSums (frozen seeds, full recompute, sums only);
// option 2a = StoredBox, 2b = StoredQueue (every step of every walk stored, exact ball test, resume without
// replay; need a cache configured with StepStore::Paths / PathsAndKeys); option 3 = Incremental (seed
// replay, Saeed's layout). ReplayQueue is the control 3q: option 3's replay and resume, but only for the walks
// popped from 2b's queue; the keys are accumulated by the walks themselves (StepStore::Keys, no stored steps).
enum class UpdateMode : int {
    Incremental = 0, FullResolve = 1, BiasedRestart = 2, StoredBox = 3, StoredQueue = 4, ReplayQueue = 5
};
const char* updateModeName(UpdateMode m);
enum class StepStore : int { None = 0, Paths = 1, PathsAndKeys = 2, Keys = 3 };   // what a cache keeps besides the entries

// Statistics of one update (one advance of the front).
struct ReuseStats {
    bool full = false;             // this update was a full solve (mode, empty cut, or not one strip)
    bool stored = false;           // options 2a/2b (stored steps)
    long long nWalks = 0, nTested = 0, nRewalked = 0, nValueUpdated = 0;   // nTested: bbox hits (2b, 3q: popped walks)
    long long stepsReplayed = 0;   // old-geometry steps re-run to find the invalidation step
    long long stepsReplayedFalse = 0;   // ... of them, in candidates the strip does not enter (replayed to the end)
    long long stepsNew = 0;        // steps run under the new geometry by resumed/restarted walks (full: all steps)
    long long escaped = 0, capped = 0;   // among the re-walked walks' new paths (full: all walks)
    long long nRewalkedAtStep0 = 0;      // re-walked walks invalidated at k = 0 (first ball)
    long long sumK = 0;                  // sum of the invalidation steps k of the re-walked walks
    long long ballsTested = 0;           // 2a/2b: stored balls run through the exact test
    long long keysComputed = 0;          // 2b: stored balls whose entry front s* was computed
    // wall time of each pass (s). Option 3: tBox (bbox; 3q: pop), tRewalk (replay + resume + new steps; 3q:
    // the keys too), tQueue (3q). 2a/2b: tBox (bbox test, 2b: pop), tTest (exact test of the stored balls),
    // tResume (resume states), tRewalk (new steps), tKeys and tQueue (2b: new keys, re-insertion). All:
    // tValue, tSums, tTotal.
    double tBox = 0, tRewalk = 0, tValue = 0, tSums = 0, tTotal = 0;
    double tTest = 0, tResume = 0, tKeys = 0, tQueue = 0;
    size_t scratchBytes = 0;             // the update's temporary buffers (capacities; per-task path scratch not included)
    double stripArea = 0;
    SolveStats solve;                    // full only: the full solve's statistics
    std::string summary() const;
};

// Options 2a/2b: one stored step of a walk, exactly as runWalk had it.
struct StoredStep {
    Vec3 c = Vec3::Zero();   // centre x_k
    double key2 = 0;         // the step's key: q.d2, or dF^2 (1 + 1e-12) when a Dirichlet face was closer
};
// The current path of one walk: step k = s[k] and f[k] (bits 0-2: onFace + 1; bit 3: the face was closer,
// so the ball's radius is the face distance). One entry per evaluated step, the terminal one included.
struct WalkPath {
    std::vector<StoredStep> s;
    std::vector<uint8_t> f;
};

// Memory of a cache (bytes; capacities, so what is allocated, not counting the heap's own block headers).
struct CacheMemory {
    size_t entries = 0;       // per-walk entries + update scratch (action, old value, counted flag)
    size_t gradWeights = 0;   // W and its update scratch (gradient points)
    size_t points = 0;        // per point: sums, seeds, query points, flags
    size_t pathHeaders = 0;   // 2a/2b: the two std::vector headers per walk
    size_t stepsAllocated = 0, stepsPayload = 0;   // 2a/2b: stored steps (capacity, size) x bytes per step
    size_t keys = 0;          // 2b: per-walk key + sorted queue + merge scratch
    long long nSteps = 0, nStepCapacity = 0;
    size_t total() const { return entries + gradWeights + points + pathHeaders + stepsAllocated + keys; }
};
constexpr size_t kBytesPerStoredStep = sizeof(StoredStep) + sizeof(uint8_t);

// One invalidated walk of an update (validation/R3 only).
struct ResumeRecord {
    size_t walk = 0;
    int k = 0;                 // invalidation step: steps 0..k-1 were kept
    int oldSteps = 0;          // step count of the walk before the update
    uint32_t rngOffset = 0;    // the walk's rngOffset before the update
    WalkState state;           // the replayed state at step k (resume state; restart: its RNG)
};

struct WalkCtx;   // walk_cache.cpp: the step function's fixed inputs for one geometry version

class WalkCache {
public:
    // walk j of point i is stored at i*nWalks + j. store: Paths keeps every step of every walk (options
    // 2a and 2b), PathsAndKeys also the per-walk entry front and the sorted queue of option 2b, Keys only
    // the entry fronts and the queue (control 3q).
    void configure(const Tet& tet, const QueryPoints& pts, int nWalks, GType g, const WalkParams& wp,
                   StepStore store = StepStore::None);
    // From scratch with the same seeds on all strips of cut (recording the paths, and the keys, if stored).
    // t = 0 (no strips): no walks, H := 0.
    SolveStats solveFull(const StripCut& cut);
    // The cut has exactly one more strip than the cache: incremental update, in the spec's order:
    // 1. bbox test, 2. replay with the old geometry + resume (or restart in BiasedRestart mode),
    // 3. value updates for g type (b), 4. sums. FullResolve mode, an empty old cut, or a cut that is
    // not exactly one strip ahead: full solve. resumes (optional) receives every invalidated walk.
    // StoredBox / StoredQueue: options 2a / 2b (advanceStored; resumes is not filled). ReplayQueue: 3q.
    ReuseStats advance(const StripCut& cut, UpdateMode mode, std::vector<ResumeRecord>* resumes = nullptr);

    EstimateD estimate(size_t i) const { return estimateFromSums(sums[i], pts.grad[i] != 0); }
    std::vector<EstimateD> estimates() const;
    size_t nPoints() const { return pts.size(); }
    int nStrips() const { return nStrips_; }
    StepStore store() const { return store_; }
    const Vec3* gradWeight(size_t walk) const;   // nullptr at non-gradient points
    CacheMemory memory() const;

    Tet tet = referenceTet();
    QueryPoints pts;
    int nWalks = 0;
    GType g = GType::One;
    WalkParams wp;

    std::vector<WalkEntry> walks;
    std::vector<int> gradSlot;    // gradient points: walk j's w is W[gradSlot[i]*nWalks + j]; -1 elsewhere
    std::vector<Vec3> W;          // first-step gradient weights (needed to update Σw·u by value updates)
    std::vector<PointSums> sums;
    // last update, per walk: bit 1 = bbox hit (tested; 2b: popped), 2 = re-walked, 4 = value updated (a
    // full solve sets 3 for every walk)
    std::vector<uint8_t> action;
    std::vector<float> rewalkedFrac;   // per point: fraction of its walks re-walked in the last update (1 after a full solve)
    std::vector<WalkPath> paths;       // 2a/2b: every walk's current path (empty with StepStore::None)
    std::vector<double> frontKey;      // 2b, 3q: per walk, min over its balls of the entry front s*_k
    std::vector<std::pair<double, uint32_t>> queue;   // 2b, 3q: (frontKey, walk), sorted ascending

private:
    ReuseStats advanceStored(const StripCut& cut, UpdateMode mode);
    void buildKeys(const StripCut& cut);
    size_t popQueue(double sNew);   // 2b, 3q: marks the queue's prefix with key <= sNew tested; returns its length
    void pushQueue(size_t nPop);    // 2b, 3q: puts that prefix back with the walks' new keys
    void valueUpdates(const WalkCtx& cNew);   // pass 3 of every incremental mode
    void updateSums(ReuseStats& rs);          // pass 4 of every incremental mode
    std::vector<double> oldU_;        // scratch of an update: old values
    std::vector<uint8_t> oldCounted_;
    std::vector<Vec3> oldW_;
    std::vector<std::pair<double, uint32_t>> queueTmp_;   // 2b: merge scratch
    StepStore store_ = StepStore::None;
    int nStrips_ = 0;
};

// Runs walks of a cache under one geometry version (strips [0, nStrips) of cut) with the one step
// function. Used by the validation (replays, the independent ball test, resume-state checks).
class WalkRunner {
public:
    WalkRunner(const WalkCache& cache, const StripCut& cut, int nStrips);
    ~WalkRunner();
    // walk = i*nWalks + j: the query point, step 0, pcg32(mixSeed(seed_i, j), j) advanced by rngOffset
    WalkState initialState(size_t walk, uint32_t rngOffset) const;
    // replays the walk from its seed (and its cached rngOffset) to the end
    void replay(size_t walk, WalkEntry& out, Vec3& w) const;
    // runs it from its seed (+ rngOffset) and stops before step k (as the step function leaves it;
    // k < 0: runs it to the end); observe (optional) sees every evaluated step's state and closest point
    WalkState runTo(size_t walk, uint32_t rngOffset, int k,
                    const std::function<void(const WalkState&, const StripCut::Closest&)>* observe) const;

private:
    const WalkCache& cache_;
    int nStrips_;
    std::unique_ptr<WalkCtx> ctx_;
};

// Fresh solve that keeps only the sums (no per-walk entries): references and fresh-seed comparisons.
std::vector<EstimateD> solveSumsOnly(const Tet& tet, const StripCut& cut, const QueryPoints& pts, int nWalks,
                                     GType g, const WalkParams& wp, SolveStats* stats);

// Reuse option 1 (frozen seeds, recompute everything): the same-seed solve on all strips of cut,
// keeping only the per-point sums. One task per point runs its walks in walk order and adds each as it
// ends, so the sums are bitwise those of WalkCache::solveFull (same walks, same order, same formula)
// without any per-walk state.
std::vector<PointSums> solveSums(const Tet& tet, const StripCut& cut, const QueryPoints& pts, int nWalks, GType g,
                                 const WalkParams& wp, SolveStats* stats);

// Bitwise comparison of two entries (and their gradient weights, either may be null). With
// ignoreRngOffset the RNG bookkeeping field is not compared (R2 against a fresh solve).
bool sameWalk(const WalkEntry& a, const Vec3* wa, const WalkEntry& b, const Vec3* wb, bool ignoreRngOffset);

// =====================================================================================
// Several boundary data on the same walks (coupling gate)
// =====================================================================================
// The walk paths do not depend on the Dirichlet values, only on which boundary absorbs. So one walk can be
// scored against K boundary data at once (the gate's H0 for both g types and its response functions P_k):
// a vector of K terminal values per walk.
struct WalkTerminal {
    WalkStatus status = WalkStatus::Landed;
    int face = -1;            // Landed on Dirichlet face 'face', or -1: on the cut
    int side = 0;             // Landed on the cut: +1 / -1 (decided at the walk's position, as finishWalk)
    Vec3 x = Vec3::Zero();    // Landed: the closest point of the cut, or the end point projected onto the face
};
// Writes the walk's K terminal values to u. Called for landed walks only: capped walks count with 0 in
// every component and escaped walks are not counted, as in every other solve.
using TerminalValues = std::function<void(const WalkTerminal&, double* u)>;
// WalkParams::oddMirrorDirichlet: writes Phi_k(x) for the K data, Phi_k harmonic in all of space and equal to
// datum k's face data on every Dirichlet face (the gate's face data are linear: their own extension). A walk
// whose gradient step was reflected across a Dirichlet face (y -> My) is scored Phi_k(y) + Phi_k(My) - u_k.
using HarmonicExtension = std::function<void(const Vec3& x, double* phi)>;

// Per point, the walks are split into nBatches contiguous batches of walksPerBatch walks (independent
// estimates). Per point and batch: n, Σu_k, Σw, Σw u_k (w = the gradient weight 3 dir0/R0 of a gradient
// point, 0 elsewhere). Per point: Σu_k u_l over all its walks, for the per-walk variance of any linear
// combination of the K data.
struct MultiSums {
    int K = 0, nBatches = 0, walksPerBatch = 0;
    size_t nPoints = 0;
    std::vector<double> batch;   // [point][batch][stride()]: n | Σu_k (K) | Σw (3) | Σw u_k (3 per k)
    std::vector<double> uu;      // [point][K*K]
    size_t stride() const { return size_t(4 + 4*K); }
    const double* batchData(size_t i, int b) const { return &batch[(i*size_t(nBatches) + size_t(b))*stride()]; }
    const double* uuData(size_t i) const { return &uu[i*size_t(K)*size_t(K)]; }
};
// Runs nWalks walks (a multiple of nBatches) from every point with the one step function, on all strips of
// cut and with the face conditions of wp, and scores each against values. Walk j of point i is the same
// walk as in every other solve of this module (pcg32(mixSeed(seed_i, j), j)). No strips and no Dirichlet
// face: no walks (stats->noCut), all sums 0. ext is required with wp.oddMirrorDirichlet (and mirrorGradBall).
MultiSums solveMulti(const Tet& tet, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                     const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats,
                     const HarmonicExtension* ext = nullptr);
// The same on a convex polyhedron with every face Neumann (the coupling gate's ring reference: the octahedron of
// check 4, in physical X; up to kMaxDomainPlanes faces), the one solve of this module that is not on a tet. The
// step function is runWalk with dom's face planes; a tet's walks are unchanged. Throws if wp flags a Dirichlet face.
MultiSums solveMulti(const ConvexDomain& dom, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                     const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats);

// =====================================================================================
// Lean multi-data solves (opt-in: the fast ring path of the coupling gate, --ring-bench, the ring view's fast option)
// =====================================================================================
// solveMulti with the lean step function: runWalk without the walk's bounding box (solveMulti never reads it; only
// WalkCache does) and, where every face absorbs, without starStep's face tests (a ball of radius <= the distance to the
// nearest face never reaches a face plane). The walks, their scoring and the order of every sum are solveMulti's, so the
// MultiSums are the same bit for bit (checked by --ring-bench --selftest); solveMulti itself is unchanged.
MultiSums solveMultiLean(const Tet& tet, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                         const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats,
                         const HarmonicExtension* ext = nullptr);
MultiSums solveMultiLean(const ConvexDomain& dom, const StripCut& cut, const QueryPoints& pts, int nWalks, int nBatches,
                         const WalkParams& wp, int K, const TerminalValues& values, SolveStats* stats);
// One solveMultiLean per job, all jobs' points in one parallel loop (one task per point, as solveMulti): no barrier
// between the jobs. Each job's sums are those of its own solveMultiLean call. stats[j].seconds is the wall time of the
// whole loop (shared by the jobs).
struct MultiJob {
    const Tet* tet = nullptr;            // the domain: a tet (solveMulti's first overload) or
    const ConvexDomain* dom = nullptr;   // a convex polyhedron with Neumann faces (the second); exactly one of them
    const StripCut* cut = nullptr;
    const QueryPoints* pts = nullptr;
    int nWalks = 0, nBatches = 1, K = 0; // K <= 16
    WalkParams wp;
    TerminalValues values;
    const HarmonicExtension* ext = nullptr;
};
std::vector<MultiSums> solveMultiLeanJobs(const std::vector<MultiJob>& jobs, std::vector<SolveStats>* stats);

// =====================================================================================
// Profiling only (partial_cut --ring-perf; no other mode calls this)
// =====================================================================================
// Where solveMulti's time goes, on ONE thread pinned to logical processor cpu (-1: not pinned), over every
// stride-th point of pts (at most maxPoints points), nWalks walks each, with exactly solveMulti's walks (tet if
// non-null, else dom):
//  - secFull: solveMulti's per-point loop body (walks, terminal values, odd-mirror scoring, the K-data accumulation);
//  - secWalks: the same walks, runWalk alone;
//  - ns*: each part of the step function timed alone in a tight loop over the recorded step states of those walks
//    (closest = StripCut::closest; faceDist = WalkCtx::dirichletFaceDist; addBall = the walk's box; sphere =
//    sampleSphere; star = starStep; outside = the escape test; per walk: init = makeInitialState, firstBall = the
//    gradient's first step without the mirror test, mirror = mirrorGradBall, finish = finishWalk, values = the
//    TerminalValues callback, accum = the K-data accumulation).
// componentsOnly = false skips the parts (secFull and secWalks only).
struct WalkProfile {
    long long points = 0, walks = 0, moves = 0, evals = 0;   // evals: closest() calls (moves + 1 per walk)
    long long landedCut = 0, landedFace = 0, capped = 0, escaped = 0, mirrored = 0;
    double secFull = 0, secWalks = 0;
    double nsClosest = 0, nsFaceDist = 0, nsSqrt = 0, nsAddBall = 0, nsSphere = 0, nsStar = 0, nsOutside = 0;
    double nsInit = 0, nsFirstBall = 0, nsMirror = 0, nsFinish = 0, nsValues = 0, nsAccum = 0;
    // a lean step function (profiling only): runWalk without the walk's box (solveMulti never reads it) and, where every
    // face absorbs, without starStep's face tests (the ball then never reaches a face plane): secWalksLean, secFullLean as
    // above; leanMismatch = walks whose end (status, steps, landing point, side, face, gradient weight) differs from runWalk's
    double secWalksLean = 0, secFullLean = 0;
    long long leanMismatch = 0;
    bool allAbsorb = false;
    int cpu = -1;
};
WalkProfile profileWalks(const Tet* tet, const ConvexDomain* dom, const StripCut& cut, const QueryPoints& pts, int nWalks,
                         const WalkParams& wp, int K, const TerminalValues& values, const HarmonicExtension* ext,
                         size_t maxPoints, int cpu, bool components = true);
// Profiling only: wall clock of solveMulti's walks alone (no scoring) at every point of pts on all threads (TBB, one task
// per point), with runWalk (lean = false) or the lean step function of WalkProfile; *moves: the walks' moves in all.
double timeWalksParallel(const Tet* tet, const ConvexDomain* dom, const StripCut& cut, const QueryPoints& pts, int nWalks,
                         const WalkParams& wp, bool lean, long long* moves);
