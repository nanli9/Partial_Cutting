#pragma once

#include "cut_problem.h"

#include <atomic>
#include <string>

// --gate coupling [--check 0|1|2|3|4] [--quick]: can one tet's WoSt solve be glued to a neighbour through
// the shared face? Tet A is the reference tet, tet B its mirror image across z = 0, and F (z = 0) the
// face they share. Every per-tet solve runs on referenceTet() in that tet's own xi frame (cut_problem.h):
// each face is Neumann or Dirichlet (WalkParams::faceBC), the cut is pulled back into xi, and the face
// data (hats, traces) are physical functions evaluated at X = X0 + J xi. By linearity
//   H = H0 + sum_k c_k P_k   (H0: cut +-g, Dirichlet faces = a fixed trace; P_k: cut 0, faces = hat_k),
// and the energy of a tet is (1/2) c^T S c + b^T c + const with S, b from split-sample products of the
// gradient estimates at a fixed quadrature point set (independent batches of walks, delete-one-batch
// jackknife for the standard errors). The gradient's first ball may cross the nearest face (even mirror
// across a Neumann face, odd mirror across a Dirichlet face with harmonic data). The verdicts of checks 0
// and 1 include the Zombie cross-check; check 1's own superposition test only confirms that the scoring is
// linear, since all its runs share one walk process.
//
// Implemented here: the geometry of A, B and their union U (a tetrahedron), the quadrature point sets and
// their checks, the multi-data solves of the own loop and Zombie (cross-check), the S / b / energy
// estimator, the linear-tet stiffness K, the sanity tests, checks 0 and 1, the reference and checks 2, 3.
// The reference is the physical ground truth on U: one own-loop solve in X (all four faces of U Neumann,
// the crack clipped to U, L = sqrt(2)), N = 25600, cached in gate_cache/ next to reportPath (a binary file
// per configuration, keyed by the configuration, seed, N, batches, the points' hash, the live walk
// parameters, U and the crack frame, and a tag for the walk code) and cross-checked by Zombie (the verdict,
// plus per-region dispersion diagnostics that bound a local discrepancy). Check 2 glues A's Level-A recipe to an uncut, linear B (a 4 x 4 system) and
// compares with the reference. Check 3 cuts both tets (crack x = 0.3, its front crossing F): both run the
// recipe in their own xi frame, glued by (S_A + S_B) c = -(b_A + b_B), at Level A (H0 with F = 0) and Level
// B (H0 with F = J, a fixed face function carrying the jump; the P_k reused). Checks 2 and 3 have no pass
// threshold (the spec measures accuracy).
// Check 4 glues one nodal ring: the 8 tets (c, +-e_x, +-e_y, +-e_z) around a node c, whose union is the octahedron O.
// Ring 1: c = 0 (every tet congruent to the reference tet); ring 2: c = (0.12, -0.08, 0.10) (8 different tets). Every
// tet runs the production recipe in its own xi frame with all four faces Dirichlet: H0 (crack +-g, faces J; Level A:
// J = 0, Level C: J = the crack's own 2D front solution) and P_k for its 4 corners (crack 0, faces lambda_k). Glue: the
// 7 x 7 system for the vertex values (glueCheck4; O's faces Neumann, imposed weakly: the outer vertices are free), with
// S_T and b_T from the physical or the xi energy of the same estimates. Reference: one own-loop solve on O in X (8
// Neumann faces: solveMulti on a ConvexDomain), cached in gate_cache/ and cross-checked by Zombie on O. Three cracks
// (K1, K2 with the front inside O; K3 cut through, whose exact solution is +-1). No pass threshold.
//
// Runs the sanity tests and the requested check (-1: every check; 4 alone skips the sanity tests and the setup of A and B,
// whose sections are then left as they are), prints everything, ends with a summary (the
// findings each part generates from its own numbers) and, unless quick, rewrites the sections it ran in reportPath
// (gate_report.md: header, summary rebuilt from the findings of every section present, the sanity tests and checks,
// the views, the appendices with the quadrature and the geometry printouts). quick: development only, fewer points and
// walks, every printed line marked QUICK, nothing written (not even the cache). Returns false if a check
// failed (for checks 2 and 3: a geometry assertion, or Zombie disagreeing with the reference) or the report
// could not be written. tipTrace (--face-trace tip): check 3's Level B uses the crack's own 2D front solution as
// its face function J instead of the spec's bump (levelBTrace in coupling_gate.cpp).
bool runCouplingGate(int check, bool quick, const std::string& reportPath, bool tipTrace = false);

// --gate coupling --view: one configuration of check 2 (config 0, 1: z_c = 0.60, 0.35) or check 3 (config 0, 1, 2:
// y_f = 0.15, 0.30, 0.45) at its sample points, for the Polyscope view of the union (main.cpp). The reference is
// loaded from ./gate_cache/ (solved and cached if missing); the decomposed field comes from the check's main recipe
// runs, with the check's own seeds, so it is the field behind the numbers in gate_report.md. No new walks beyond those.
struct GateView {
    std::string title;
    std::array<Vec3, 5> verts{};                  // global v0..v4
    std::array<std::array<size_t, 4>, 2> tets{};  // A and B as global ids
    std::vector<Vec3> crack;                      // U's crack polygon (X)
    Vec3 n = Vec3::Zero();                        // crack plane n.X = d, n towards the +g side
    double d = 0;
    std::vector<Vec3> X, XF;                      // the sample points of A, then B; the F points
    struct Field {
        std::string name;
        std::vector<double> vol, face;            // at X and at XF (on F the decomposed field is its trace J + sum c hat)
        bool difference = false;                  // decomposed - reference (symmetric colour range)
        bool shot = false;                        // one of the --screenshot images
        double rms = 0, max = 0;                  // a difference: its RMS (quadrature-weighted) and max over A u B, / 2 g_max
    };
    std::vector<Field> fields;
};
GateView gateView(int check, int config, bool tipTrace = false);

// --gate coupling --view --check 4: one nodal ring of check 4 (ring 1 or 2, crack K1, K2 or K3) at a chosen sample
// density, for the interactive Polyscope view (main.cpp). With the check's own setting (4096 points per tet, N = 2560
// per tet, 25600 for the reference) it is the field behind gate_report.md's check 4 and loads the reference from
// ./gate_cache/; any other setting solves the reference afresh at its own walk count and caches nothing.
struct RingViewParams {
    int ring = 2;              // 1: congruent (c = 0), 2: real mesh (c moved off-centre)
    int crack = 0;             // 0, 1, 2: K1, K2, K3
    int pointsPerTet = 4096;   // sample point density (about; refined at the front and where the crack meets the faces)
    int walksTet = 2560;       // walks per point of the per-tet runs
    int walksRef = 25600;      // walks per point of the reference solve on the whole ring
    bool geometryOnly = false; // only the ring, its crack and the title (no walks)
    // The fast method (off by default; off = check 4's path above, unchanged): every tet still solves in its own xi frame,
    // with the exactly known parts taken out of the walks (P_k = lambda_k + Q_k, H0 = J_xi + R, see --ring-bench) and the
    // lean step function; the same glue. walksTet is then Q's walks per point (cut tets), walksR R's (where J - J_xi is not
    // 0 on the faces). The reference is loaded from gate_cache/ at 4096 points and 25600 walks whatever walksTet is.
    bool fast = false;
    int walksR = 640;
};
struct RingView {
    std::string title;
    std::array<Vec3, 7> verts{};                  // v0 = the node, v1..v6 = the outer vertices
    std::vector<std::array<size_t, 4>> tets;      // the 8 tets as global ids
    std::vector<std::string> tetName;
    std::vector<char> tetCut;                     // the crack cuts the tet
    std::vector<Vec3> crack;                      // the crack on the ring (polygon, X)
    Vec3 n = Vec3::Zero();                        // crack plane n.X = d, n towards the +g side
    double d = 0;
    std::vector<Vec3> X;                          // the sample points of the 8 tets
    std::vector<int> tetOf;                       // the tet of each point
    std::vector<GateView::Field> fields;          // vol only (face unused)
    // basis[v]: vertex v's basis function at the points: in every tet with corner v its response P_k (crack 0, faces =
    // the hat of v), 0 in the other tets; inSupport[v]: the point's tet has corner v
    std::array<std::vector<double>, 7> basis;
    std::array<std::vector<char>, 7> inSupport;
    std::vector<std::string> report;              // lines for the panel: setting, timings, errors, vertex values
    int pointsPerTet = 0, walksTet = 0, walksRef = 0;   // as run (walks rounded to whole batches)
    bool fast = false;                            // made by the fast method (walksTet = Q's walks, walksR = R's)
    int walksR = 0;
};
// progress (optional): incremented after each tet's two runs (8 in all)
RingView gateRingView(const RingViewParams& p, std::atomic<int>* progress = nullptr);

// Rewrites one marked section of gate_report.md at path (the view's, id "view_..."), keeping the others, and rebuilds
// the report's summary. False if the file does not exist or cannot be written.
bool writeGateReportSection(const std::string& path, const std::string& id, const std::string& body);

// --ring-perf <parts> [--ring 1|2] [--config 0|1|2]: profiling of check 4's ring path for the performance work (prints only;
// reads gate_cache/, never writes it or gate_report.md). parts, run in order: x = J_xi checks, a = the current per-tet runs
// timed per tet, b = walk stop distance sweep (length, time, paired bias), c = the direct ring solve's cost, d = the walk
// step profile (per part, per logical processor, thread scaling), e = current path vs control-variate prototype vs direct
// solve against the cached reference (this ring and crack), f = e for all six configurations, r / R = review: the fast
// path against check 4's two runs per tet at the same points and N = 2560 / 10240 (Grams, glued c, P_k, H0, glued H z tests).
bool runRingPerf(int ring, int crack, const std::string& parts);

// --ring-bench: check 4's nodal ring solved by the fast per-tet method (A: per-tet control variates P_k = lambda_k + Q_k,
// H0 = J_xi + R in each tet's own xi frame, the lean step function, check 4's glue) and by the direct ring solve (B: the
// reference's solver on O at N walks per point, its own seeds), optionally check 4's current path, every one timed (wall
// clock, 16 threads, median of the repeats for runs under 30 s) and compared with the cached reference (read-only) at
// check 4's point sets: H and grad H errors (overall and per region), energies, vertex values. Writes
// ring_perf/<tag>_rows.csv (one row per method, setting, g type and seed), ring_perf/<tag>_summary.csv (equal time and
// equal accuracy per configuration) and ring_perf/<tag>_log.txt; never writes gate_report.md or gate_cache/.
// Lists: ours "nQ/nR[/B],...", direct and current "N[/B],..." ("none": skip; B defaults to 8 up to 32 walks, 16 below
// 1280, 64 from 1280). selftest: the lean walks against solveMulti bit for bit, J_xi, the skip decisions, same-target z
// tests and determinism. view: gateRingView(viewParams) headless (the panel's numbers on the console).
struct RingBenchParams {
    int ring = 0, crack = -1;                    // 0 / -1: all
    std::string ours = "default", direct = "default", current = "default";
    double shellRel = 1e-6;                      // walk stop distance / sqrt(2), ours and the direct (the reference's: 1e-6)
    int repeats = 3;                             // timing repeats (median), for runs under 30 s (main: 1 with --quick)
    int seeds = 1;                               // independent seed sets (the spread of the errors)
    bool lean = true, fuse = true, quick = false;
    std::string tag;                             // output name (default: bench_<date>_<time>)
    bool selftest = false, view = false;
    RingViewParams viewParams;
};
bool runRingBench(const RingBenchParams& p);
