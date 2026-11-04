// Harmonic partial-cut enrichment H(ξ) on a single tetrahedron, in its local coordinates ξ (the
// reference tet), by Walk on Stars, shown in Polyscope, with incremental walk reuse as the cut
// front advances.
//
//   partial_cut                                   GUI (own WoSt loop with the walk cache, walk_cache.h)
//   partial_cut --validate                        the first gate's four checks (Zombie), headless;
//                                                 rewrites the first-gate section of ./results.txt
//   partial_cut --bench                           one full default Zombie solve (lattice + gradient + slice), timed
//   partial_cut --validate-reuse [--quick]        the reuse validation (V0, R1-R4, checks 1-4); rewrites the
//                                                 reuse section of ./results.txt (--quick: small R1-R4 run, no file)
//   partial_cut --bench-reuse                     one default incremental sweep, with full-solve timings
//   partial_cut --compare-options [--quick] [--threads T] [--repeat R]
//                                                 the exact reuse options 1, 2a, 2b, 3 and the control 3q compared
//                                                 (equivalence, timings, memory); writes ./results_options.txt
//                                                 (--threads T: ./results_options_threadsT.txt)
//   partial_cut --check-xi                        H lives in xi: pull-back and rest-tet independence checks;
//                                                 rewrites the xi section of ./results.txt
//   partial_cut --gate coupling [--check 0|1|2|3|4] [--quick]
//                                                 the coupling gate (two tets glued through a shared face; check 4:
//                                                 the 8 tets of one nodal ring glued at its vertices);
//                                                 rewrites its sections of ./gate_report.md and caches the reference
//                                                 solves in ./gate_cache/ (--quick: writes nothing)
//   partial_cut --gate coupling --view [--check 2|3] [--config i] [--screenshot out.png [--with-ui]]
//                                                 its Polyscope view of the union (reference, decomposed, difference, a
//                                                 slice through the crack plane); the reference from ./gate_cache/;
//                                                 --screenshot also rewrites the view's section of ./gate_report.md
//   partial_cut --gate coupling --view --check 4 [--ring 1|2] [--config 0|1|2] [--fast]
//                                                 check 4's nodal ring with a panel: ring, crack, sample points per
//                                                 tet and walk counts, recomputed in the background; a "fast method"
//                                                 checkbox (off unless --fast) switches to the fast per-tet path
//   partial_cut --ring-bench [...] / --ring-bench --selftest / --ring-bench --view [--fast]
//                                                 the fast ring method vs the direct ring solve (writes ring_perf/)
//   partial_cut --screenshot out.png [--t 0.6] [--t-start 0.2 --steps 8] [--rewalked] [--with-ui]
//
// Geometry (rest tet, cut plane, front direction): GeometryConfig in cut_problem.h. The solve always
// runs on the reference tet in ξ; the rest tet only maps the cut into ξ.

#include "coupling_gate.h"
#include "cut_problem.h"
#include "validation.h"
#include "validation_reuse.h"
#include "validation_xi.h"
#include "walk_cache.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/screenshot.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include "polyscope/volume_mesh.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <future>
#include <string>
#include <system_error>
#include <vector>

namespace {

std::vector<std::array<double, 3>> toArrays(const std::vector<Vec3>& v)
{
    std::vector<std::array<double, 3>> a;
    a.reserve(v.size());
    for (const Vec3& p : v) a.push_back({p.x(), p.y(), p.z()});
    return a;
}

void printGeometry(const Tet& tet, const GeometryConfig& cfg)
{
    Cut full = buildCut(tet, cfg, 1.0);
    if (!isReferenceTet(Tet(cfg.restVerts))) {
        const auto& X = cfg.restVerts;
        std::printf("rest tet (only maps the cut into xi): (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f)\n",
                    X[0].x(), X[0].y(), X[0].z(), X[1].x(), X[1].y(), X[1].z(), X[2].x(), X[2].y(), X[2].z(),
                    X[3].x(), X[3].y(), X[3].z());
    }
    std::printf("H(xi) on the reference tet: (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f), "
                "longest edge L = %.4f\n",
                tet.v[0].x(), tet.v[0].y(), tet.v[0].z(), tet.v[1].x(), tet.v[1].y(), tet.v[1].z(),
                tet.v[2].x(), tet.v[2].y(), tet.v[2].z(), tet.v[3].x(), tet.v[3].y(), tet.v[3].z(), tet.L);
    std::printf("cut plane in xi n.xi = d: n = (%.4f, %.4f, %.4f), d = %.4f; front direction u = (%.4f, %.4f, %.4f)\n",
                full.n.x(), full.n.y(), full.n.z(), full.d, full.u.x(), full.u.y(), full.u.z());
    if (full.fullIsQuad)
        std::printf("plane-tet section is a quad (check passed)\n");
    else
        std::printf("WARNING: the plane-tet section is not a quad (%zu vertices)\n", full.fullPolygon.size());
}

// ---------------------------------------------------------------------------------------
// GUI / screenshot scene: the own WoSt loop with the walk cache. The cut is a list of strips;
// raising t appends one strip and updates the cache (incremental / full re-solve / biased restart).
// ---------------------------------------------------------------------------------------
struct App {
    GeometryConfig cfg;
    Tet tet = referenceTet();
    Settings s;
    StripCut cut{tet, cfg};
    Lattice lattice;                // ALL lattice points: epsilon only masks them (the excluded set changes with t)
    Slice slice;
    WalkCache cache, fullCache;     // fullCache: scratch for timing a same-seed full re-solve
    bool fullCacheConfigured = false;
    std::vector<EstimateD> est;     // lattice points first, then slice points
    std::vector<char> masked;       // lattice points within epsilon of the current cut

    // UI copies of the settings (ImGui wants float/int)
    float t = 0.6f, epsRel = 0.01f, sliceOffset = 0.0f, sliceYaw = 0.0f, slicePitch = 0.0f;
    int nWalks = 256, k = 12, gType = 0, mode = 0;
    bool timeFull = true;
    uint64_t seed = 0, seedCount = 0;

    SolveStats lastSolve;           // the last full solve of the cache
    ReuseStats lastUpdate;
    bool haveUpdate = false;
    double lastFullSeconds = -1;    // same-seed full re-solve after the last incremental update
    std::string lastAction;
    bool firstRegistration = true, gradCloudConfigured = false;
    std::string validationMsg;

    App()
    {
        t = float(s.t);
        epsRel = float(s.epsRel);
        nWalks = s.nWalks;
        k = s.k;
        gType = int(s.g);
        seed = s.seed;
    }

    void registerTet()
    {
        std::vector<std::array<size_t, 4>> cells = {{0, 1, 2, 3}};
        polyscope::VolumeMesh* vm =
            polyscope::registerTetMesh("tet", toArrays(std::vector<Vec3>(tet.v.begin(), tet.v.end())), cells);
        vm->setTransparency(0.25f);
        vm->setColor({0.75f, 0.75f, 0.80f});
        vm->setEdgeWidth(1.0);
    }

    void syncSettings()
    {
        s.t = t;
        s.epsRel = epsRel;
        s.nWalks = nWalks;
        s.k = k;
        s.g = GType(gType);
        s.slice.offset = sliceOffset;
        s.slice.yawDeg = sliceYaw;
        s.slice.pitchDeg = slicePitch;
    }

    void moveCutTo(double tNew)
    {
        double sNew = cut.sOf(tNew);
        if (sNew > cut.front()) cut.addStrip(sNew);
        else if (sNew < cut.front()) cut.retractTo(sNew);
    }

    // the cached points (lattice, then slice) with the current seed (lattice tag 1, slice tag 2)
    QueryPoints queryPoints() const
    {
        QueryPoints q;
        q.add(lattice.pts, &lattice.interior, mixSeed(seed, 1));
        q.add(slice.verts, nullptr, mixSeed(seed, 2));
        return q;
    }

    // GUI mode index -> update mode: 0 incremental (option 3), 1 full re-solve with a fresh seed, 2 biased
    // restart, 3 option 1 (frozen seeds: a same-seed full recompute), 4 option 2b (stored steps + queue)
    static constexpr int kModeFrozen = 3, kModeQueue = 4;
    UpdateMode updateMode() const
    {
        return mode == kModeFrozen ? UpdateMode::FullResolve : mode == kModeQueue ? UpdateMode::StoredQueue : UpdateMode(mode);
    }
    StepStore stepStore() const { return mode == kModeQueue ? StepStore::PathsAndKeys : StepStore::None; }

    // New query points and seeds, full solve on the current strip list.
    void rebuild(const char* why)
    {
        syncSettings();
        moveCutTo(s.t);
        lattice = buildLattice(tet, Cut{}, s.k, 0.0);
        slice = buildSlice(tet, cut.frame(), s.slice);
        cache.configure(tet, queryPoints(), s.nWalks, s.g, s.walk, stepStore());
        fullCacheConfigured = false;
        lastSolve = cache.solveFull(cut);
        haveUpdate = false;
        lastFullSeconds = -1;
        lastAction = std::string("full solve (") + why + ")";
        std::printf("%s: t = %.3f, N = %d, %s, k = %d, %d strips: %s\n", lastAction.c_str(), s.t, s.nWalks,
                    gTypeName(s.g), s.k, cut.nStrips(), lastSolve.summary().c_str());
        show();
    }

    // t released at a new value. Incremental, biased and option 2b modes: an update of the cache when t
    // increases, a full re-solve with the same seeds when it decreases (retraction is not incremental).
    // Option 1: a same-seed full recompute on every change. Full re-solve mode: a full solve with a fresh
    // seed on every change, as the first gate's GUI did, so the step-to-step noise that reuse avoids
    // (check 4 of the reuse validation) is visible.
    void setFront(float tNew)
    {
        syncSettings();
        const double sOld = cut.front(), sNew = cut.sOf(tNew);
        if (sNew == sOld) return;
        char buf[200];
        if (mode == int(UpdateMode::FullResolve)) {
            moveCutTo(tNew);
            seed = mixSeed(s.seed, ++seedCount);
            cache.configure(tet, queryPoints(), s.nWalks, s.g, s.walk);
            fullCacheConfigured = false;
            lastSolve = cache.solveFull(cut);
            lastUpdate = ReuseStats{};
            lastUpdate.full = true;
            lastUpdate.nWalks = lastSolve.walksRequested;
            lastUpdate.tTotal = lastSolve.seconds;
            lastUpdate.solve = lastSolve;
            haveUpdate = true;
            lastFullSeconds = lastSolve.seconds;
            std::snprintf(buf, sizeof buf, "full re-solve with a fresh seed: front s %.4f -> %.4f", sOld, sNew);
            lastAction = buf;
            std::printf("%s: %s\n", lastAction.c_str(), lastSolve.summary().c_str());
        } else if (sNew > sOld) {
            cut.addStrip(sNew);
            const UpdateMode um = updateMode();
            lastUpdate = cache.advance(cut, um);
            haveUpdate = true;
            std::snprintf(buf, sizeof buf, "%s: front s %.4f -> %.4f",
                          mode == kModeFrozen ? "frozen seeds, full recompute (option 1)"
                          : lastUpdate.full   ? "full solve (no cut before)"
                                              : updateModeName(um),
                          sOld, sNew);
            lastAction = buf;
            if (lastUpdate.full) {
                lastSolve = lastUpdate.solve;
                lastFullSeconds = lastUpdate.tTotal;
            } else if (timeFull) {   // same N, same seeds, same geometry, right after the update
                if (!fullCacheConfigured) {
                    fullCache.configure(tet, cache.pts, s.nWalks, s.g, s.walk);
                    fullCacheConfigured = true;
                }
                lastFullSeconds = fullCache.solveFull(cut).seconds;
            } else {
                lastFullSeconds = -1;
            }
            if (lastFullSeconds >= 0)
                std::printf("%s: %s; same-N full re-solve %.3f s\n", lastAction.c_str(), lastUpdate.summary().c_str(),
                            lastFullSeconds);
            else
                std::printf("%s: %s; same-N full re-solve not timed\n", lastAction.c_str(), lastUpdate.summary().c_str());
        } else {
            cut.retractTo(sNew);
            lastSolve = cache.solveFull(cut);
            haveUpdate = false;
            lastFullSeconds = lastSolve.seconds;
            std::snprintf(buf, sizeof buf, "t decreased (front s %.4f -> %.4f): full re-solve, same seeds "
                          "(retraction is not incremental)", sOld, sNew);
            lastAction = buf;
            std::printf("%s: %s\n", lastAction.c_str(), lastSolve.summary().c_str());
        }
        show();
    }

    void show()
    {
        est = cache.estimates();
        const size_t nL = lattice.pts.size();
        const double eps = s.epsRel*tet.L;
        masked.assign(nL, 0);
        for (size_t i = 0; i < nL; i++) masked[i] = std::sqrt(cut.closest(lattice.pts[i], cut.nStrips()).d2) < eps;

        // the cut: all strips, coloured by strip index (the newest has the highest)
        if (cut.nStrips() == 0) {
            polyscope::removeStructure("cut", false);
        } else {
            std::vector<Vec3> V;
            std::vector<std::array<int, 3>> Fs;
            std::vector<double> strip;
            cut.mesh(cut.nStrips(), V, Fs, strip);
            polyscope::SurfaceMesh* cm = polyscope::registerSurfaceMesh("cut", toArrays(V), Fs);
            cm->setSurfaceColor({0.95f, 0.70f, 0.15f});
            cm->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
            cm->setEdgeWidth(1.0);
            cm->addFaceScalarQuantity("strip index", strip);
        }

        // lattice: H, its standard error and the re-walked fraction at every non-masked point
        std::vector<Vec3> pts;
        std::vector<double> H, Hse, rew;
        for (size_t i = 0; i < nL; i++) {
            if (masked[i]) continue;
            pts.push_back(lattice.pts[i]);
            H.push_back(est[i].H);
            Hse.push_back(est[i].Hse);
            rew.push_back(cache.rewalkedFrac[i]);
        }
        polyscope::PointCloud* pc = polyscope::registerPointCloud("lattice", toArrays(pts));
        auto* hq = pc->addScalarQuantity("H", H);
        hq->setColorMap("coolwarm");
        hq->setMapRange({-1.0, 1.0});
        pc->addScalarQuantity("H std. error", Hse);
        auto* rq = pc->addScalarQuantity("re-walked this step", rew);
        rq->setColorMap("viridis");
        rq->setMapRange({0.0, 1.0});

        // ∇H at the interior non-masked lattice points only (its own cloud, as in the first gate)
        std::vector<Vec3> gPts;
        std::vector<std::array<double, 3>> grad;
        std::vector<double> gradSe, gnorm;
        for (size_t i = 0; i < nL; i++) {
            if (masked[i] || !est[i].hasGrad) continue;
            gPts.push_back(lattice.pts[i]);
            grad.push_back({est[i].grad.x(), est[i].grad.y(), est[i].grad.z()});
            gradSe.push_back(est[i].gradSe.norm());
            gnorm.push_back(est[i].grad.norm());
        }
        if (gPts.empty()) {
            polyscope::removeStructure("lattice grad H", false);
        } else {
            polyscope::PointCloud* gc = polyscope::registerPointCloud("lattice grad H", toArrays(gPts));
            auto* gq = gc->addVectorQuantity("grad H", grad);
            gc->addScalarQuantity("|std. error of grad H|", gradSe);
            // display scaling only (see the first gate): 90th percentile of |grad H|, never below 1/L
            std::nth_element(gnorm.begin(), gnorm.begin() + gnorm.size()*9/10, gnorm.end());
            gq->setVectorLengthRange(std::max(gnorm[gnorm.size()*9/10], 1.0/tet.L));
            if (!gradCloudConfigured) {
                gq->setEnabled(true);
                gq->setVectorLengthScale(0.6*tet.L/std::max(1, k), /*isRelative=*/false);
                gq->setVectorRadius(0.0025);
                gq->setVectorColor({0.1f, 0.1f, 0.1f});
                gc->setPointRadius(0.004);
                gc->setPointColor({0.1f, 0.1f, 0.1f});
                gc->setEnabled(false);   // hidden at startup; tick it in the structure list
                gradCloudConfigured = true;
            }
        }

        // slice: H and the re-walked fraction
        if (slice.tris.empty()) {
            polyscope::removeStructure("slice", false);
        } else {
            std::vector<double> sH, sRew;
            for (size_t v = 0; v < slice.verts.size(); v++) {
                sH.push_back(est[nL + v].H);
                sRew.push_back(cache.rewalkedFrac[nL + v]);
            }
            polyscope::SurfaceMesh* sm = polyscope::registerSurfaceMesh("slice", toArrays(slice.verts), slice.tris);
            sm->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
            auto* sq = sm->addVertexScalarQuantity("H", sH);
            sq->setColorMap("coolwarm");
            sq->setMapRange({-1.0, 1.0});
            auto* srq = sm->addVertexScalarQuantity("re-walked this step", sRew);
            srq->setColorMap("viridis");
            srq->setMapRange({0.0, 1.0});
            if (firstRegistration) sq->setEnabled(true);
        }
        if (firstRegistration) {
            hq->setEnabled(true);
            pc->setPointRadius(0.006);
            firstRegistration = false;
        }
    }

    void showRewalked()   // switch the lattice and slice colours to "re-walked this step"
    {
        if (polyscope::hasPointCloud("lattice"))
            polyscope::getPointCloud("lattice")->getQuantity("re-walked this step")->setEnabled(true);
        if (polyscope::hasSurfaceMesh("slice"))
            polyscope::getSurfaceMesh("slice")->getQuantity("re-walked this step")->setEnabled(true);
    }

    // H +- 1 s.e. at the four tet vertices (lattice points with barycentric index k, pulled 0.1% inside)
    std::string vertexHText() const
    {
        std::string out;
        char buf[64];
        for (int v = 0; v < 4; v++) {
            std::snprintf(buf, sizeof(buf), "v%d excluded  ", v);
            for (size_t i = 0; i < lattice.index.size(); i++)
                if (lattice.index[i][v] == lattice.k && !masked[i])
                    std::snprintf(buf, sizeof(buf), "v%d %+.3f+-%.3f  ", v, est[i].H, est[i].Hse);
            out += buf;
        }
        return out;
    }

    void setCamera()
    {
        Cut c = buildCut(tet, cfg, 1.0);
        Vec3 dir = (-c.w - 0.35*c.n + 0.45*c.u).normalized();
        Vec3 eye = tet.centroid + 1.8*tet.L*dir;
        polyscope::view::lookAt(glm::vec3(eye.x(), eye.y(), eye.z()),
                                glm::vec3(tet.centroid.x(), tet.centroid.y(), tet.centroid.z()));
    }

    void ui()
    {
        bool tChanged = false, rebuildNeeded = false;
        const ImGuiSliderFlags clamp = ImGuiSliderFlags_AlwaysClamp;
        ImGui::PushItemWidth(170);

        ImGui::TextUnformatted("Cut front (updates when the slider is released)");
        ImGui::SliderFloat("t (front position)", &t, 0.0f, 1.0f, "%.3f", clamp);
        tChanged |= ImGui::IsItemDeactivatedAfterEdit();
        const char* modes[] = {"incremental (resume, option 3)", "full re-solve each change (fresh seed)",
                               "biased restart (Fig. 5)", "frozen seeds, full recompute (option 1)",
                               "stored steps + front-key queue (option 2b)"};
        int newMode = mode;
        if (ImGui::Combo("update mode", &newMode, modes, 5) && newMode != mode) {
            // leaving biased mode: its cache holds restarted walks, so it is rebuilt (same seeds); entering or
            // leaving option 2b builds or drops the stored steps (a rebuild too)
            if (mode == int(UpdateMode::BiasedRestart) || mode == kModeQueue || newMode == kModeQueue) rebuildNeeded = true;
            mode = newMode;
        }
        ImGui::TextWrapped("Incremental, biased and option 2b modes: raising t appends one strip to the cut and updates "
                           "the cache; lowering t is a full re-solve with the same seeds (retraction is not incremental). "
                           "Option 1: a same-seed full recompute on every change (the estimator every exact option "
                           "reproduces). Full re-solve mode: a full solve with a fresh seed on every change.");
        ImGui::Checkbox("time a same-N full re-solve after each update", &timeFull);

        ImGui::Separator();
        ImGui::TextUnformatted("Solver (each change: new cache, full solve)");
        ImGui::SliderInt("N walks per point", &nWalks, 16, 4096, "%d", clamp | ImGuiSliderFlags_Logarithmic);
        rebuildNeeded |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SliderFloat("epsilon / L", &epsRel, float(10*s.walk.absorbShellRel), 0.1f, "%.4f",
                           clamp | ImGuiSliderFlags_Logarithmic);
        rebuildNeeded |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::TextWrapped("  epsilon = %.4g (masks lattice points; walks stop within %.2g of the cut)", epsRel*tet.L,
                    s.walk.absorbShellRel*tet.L);
        ImGui::TextUnformatted("Dirichlet data g on the cut (H = +g / -g):");
        rebuildNeeded |= ImGui::RadioButton("(a) g = 1", &gType, 0);
        ImGui::SameLine();
        rebuildNeeded |= ImGui::RadioButton("(b) g = min(1, sqrt(r/L))", &gType, 1);
        ImGui::SliderInt("k (lattice subdivisions)", &k, 1, 30, "%d", clamp);
        rebuildNeeded |= ImGui::IsItemDeactivatedAfterEdit();

        ImGui::Separator();
        ImGui::TextUnformatted("Slice plane (its points are cached too)");
        float R = 0.0f;
        for (const Vec3& v : tet.v) R = std::max(R, float((v - tet.centroid).norm()));
        ImGui::SliderFloat("offset from centroid", &sliceOffset, -R, R, "%.3f", clamp);
        rebuildNeeded |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SliderFloat("yaw (deg, about n)", &sliceYaw, -180.0f, 180.0f, "%.1f", clamp);
        rebuildNeeded |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SliderFloat("pitch (deg, towards n)", &slicePitch, -90.0f, 90.0f, "%.1f", clamp);
        rebuildNeeded |= ImGui::IsItemDeactivatedAfterEdit();

        ImGui::Separator();
        bool newSeed = ImGui::Button("Re-solve (new seed)");
        ImGui::SameLine();
        if (ImGui::Button("Run validation")) {
            std::string path = (std::filesystem::current_path()/"results.txt").string();
            auto t0 = std::chrono::steady_clock::now();
            bool ok = runValidation(cfg, s, path);
            double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            validationMsg = ok ? "first-gate validation (" + std::to_string(int(sec + 0.5)) + " s) written to " + path
                               : "validation: could not write " + path;
        }
        ImGui::SameLine();
        if (ImGui::Button("Run reuse validation")) {
            std::string path = (std::filesystem::current_path()/"results.txt").string();
            auto t0 = std::chrono::steady_clock::now();
            bool ok = runReuseValidation(cfg, s, path, false);
            double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            validationMsg = "reuse validation (" + std::to_string(int(sec + 0.5)) + " s, " +
                            (ok ? "exact checks pass" : "an exact check FAILED or the file could not be written") +
                            ") written to " + path;
        }
        ImGui::TextWrapped("The validations run on the GUI thread: the window freezes (first gate about 5 s, reuse "
                           "validation about a minute at the defaults). Both use the current N, epsilon and k.");
        if (!validationMsg.empty()) ImGui::TextWrapped("%s", validationMsg.c_str());

        ImGui::Separator();
        ImGui::TextWrapped("last: %s", lastAction.c_str());
        if (cut.nStrips() == 0) {
            ImGui::TextWrapped("t = 0: no cut. Laplace with zero Neumann data on every face has no unique "
                               "solution, so no walks are run and H is set to 0.");
        } else {
            ImGui::TextWrapped("cut: %d strips, %d triangles, area %.4f of %.4f, front at s = %.4f", cut.nStrips(), cut.nTris(),
                        cut.area(cut.nStrips()), cut.fullArea(), cut.front());
        }
        if (haveUpdate && !lastUpdate.full) {
            const double nW = double(std::max(1LL, lastUpdate.nWalks));
            ImGui::Text("last update (%lld walks, strip area %.4f):", lastUpdate.nWalks, lastUpdate.stripArea);
            ImGui::Text("  %s %.2f%%, re-walked %.2f%%, value-updated %.2f%%",
                        lastUpdate.stored ? "popped (front key <= front)" : "tested (bbox hit)", 100.0*lastUpdate.nTested/nW,
                        100.0*lastUpdate.nRewalked/nW, 100.0*lastUpdate.nValueUpdated/nW);
            if (lastUpdate.nRewalked > 0)
                ImGui::Text("  re-walked walks resumed at step k: %.1f%% at k = 0, mean k = %.2f",
                            100.0*double(lastUpdate.nRewalkedAtStep0)/double(lastUpdate.nRewalked),
                            double(lastUpdate.sumK)/double(lastUpdate.nRewalked));
            if (lastUpdate.stored) {
                ImGui::TextWrapped("  re-walked paths: %lld escaped, %lld capped; stored balls tested %lld, new steps %lld",
                                   lastUpdate.escaped, lastUpdate.capped, lastUpdate.ballsTested, lastUpdate.stepsNew);
                ImGui::TextWrapped("  update %.4f s = pop %.4f + exact test %.4f + resume %.4f + new steps %.4f + keys %.4f + "
                                   "values %.4f + sums %.4f + queue %.4f", lastUpdate.tTotal, lastUpdate.tBox, lastUpdate.tTest,
                                   lastUpdate.tResume, lastUpdate.tRewalk, lastUpdate.tKeys, lastUpdate.tValue, lastUpdate.tSums,
                                   lastUpdate.tQueue);
            } else {
                ImGui::TextWrapped("  re-walked paths: %lld escaped, %lld capped; steps replayed (old cut) %lld, new %lld",
                                   lastUpdate.escaped, lastUpdate.capped, lastUpdate.stepsReplayed, lastUpdate.stepsNew);
                ImGui::TextWrapped("  update %.4f s = bbox %.4f + re-walk %.4f + values %.4f + sums %.4f", lastUpdate.tTotal,
                                   lastUpdate.tBox, lastUpdate.tRewalk, lastUpdate.tValue, lastUpdate.tSums);
            }
            if (lastFullSeconds >= 0)
                ImGui::TextWrapped("  full re-solve, same N and seeds, same geometry, run right after: %.4f s (%.1fx)",
                            lastFullSeconds, lastFullSeconds/std::max(1e-9, lastUpdate.tTotal));
            else
                ImGui::TextUnformatted("  full re-solve time: not measured (checkbox off)");
        } else if (haveUpdate) {
            ImGui::Text("last update: full re-solve of %lld walks, %.4f s", lastUpdate.nWalks, lastUpdate.tTotal);
        }
        int nMasked = 0;
        for (char c : masked) nMasked += c;
        ImGui::TextWrapped("lattice: %zu points cached, %d within epsilon of the cut masked", lattice.pts.size(), nMasked);
        ImGui::TextWrapped("H at vertices: %s", vertexHText().c_str());
        ImGui::Text("slice: %zu points", slice.verts.size());
        if (!lastSolve.noCut)
            ImGui::TextWrapped("last full solve: %.3f s; %lld walks counted, %lld escaped, %lld capped; mean length %.1f",
                        lastSolve.seconds, lastSolve.walksCounted, lastSolve.walksEscaped, lastSolve.walksCapped,
                        lastSolve.meanWalkLength());
        ImGui::PopItemWidth();

        if (newSeed) {
            seed = mixSeed(s.seed, ++seedCount);
            rebuild("new seed");
        } else if (rebuildNeeded) {
            rebuild("settings changed");
        } else if (tChanged) {
            setFront(t);
        }
    }
};

void initPolyscope()
{
    polyscope::options::programName = "partial_cut";
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::init();
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
}

int runBench()
{
    GeometryConfig cfg;
    Tet tet = referenceTet();
    Settings s;
    printGeometry(tet, cfg);
    SceneResult r = solveScene(tet, cfg, s, s.seed);
    std::printf("bench: one full default solve: t = %.2f, N = %d, epsilon = %.4g L, %s, k = %d, slice res %d\n",
                s.t, s.nWalks, s.epsRel, gTypeName(s.g), s.k, s.slice.res);
    std::printf("  cut: %zu-gon, %zu triangles, area %.4f, front s = %.4f\n", r.cut.polygon.size(), r.cut.tris.size(),
                r.cut.area(), r.cut.s);
    int nInterior = 0;
    for (char c : r.lattice.interior) nInterior += c;
    std::printf("  lattice: %d points, %d excluded (within epsilon of the cut), %zu solved, %d with gradient\n",
                r.lattice.nTotal, r.lattice.nExcluded, r.lattice.pts.size(), nInterior);
    std::printf("  slice: %zu points, %zu triangles\n", r.slice.verts.size(), r.slice.tris.size());
    std::printf("  lattice solve: %s\n", r.latticeStats.summary().c_str());
    std::printf("  slice solve:   %s\n", r.sliceStats.summary().c_str());

    // worst per-point escape rate, and gradient noise
    double worstEsc = 0;
    Vec3 worstPt = Vec3::Zero();
    std::vector<double> gse, gnorm;
    for (size_t i = 0; i < r.latticeEst.size(); i++) {
        const Estimate& e = r.latticeEst[i];
        int req = e.hasGrad ? 2*std::max(1, s.nWalks/2) : s.nWalks;
        double esc = 1.0 - double(e.nCounted)/req;
        if (esc > worstEsc) { worstEsc = esc; worstPt = r.lattice.pts[i]; }
        if (e.hasGrad) { gse.push_back(e.gradSe.norm()); gnorm.push_back(e.grad.norm()); }
    }
    auto median = [](std::vector<double> v) {
        if (v.empty()) return 0.0;
        std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
        return v[v.size()/2];
    };
    std::printf("  worst per-point escape rate on the lattice: %.2f%% at (%.4f, %.4f, %.4f)\n", 100*worstEsc,
                worstPt.x(), worstPt.y(), worstPt.z());
    std::printf("  grad H at %zu interior points: median |grad H| = %.3f, median |s.e.| = %.3f\n", gse.size(),
                median(gnorm), median(gse));
    std::printf("  wall time: build %.3f s, total %.3f s per full solve\n", r.buildSeconds, r.totalSeconds);
    return 0;
}

// The scene with the front at t, reached from tStart by nSteps incremental updates (nSteps = 0:
// one full solve at t); showRewalked colours the lattice and the slice by "re-walked this step".
int runScreenshot(const std::string& file, double t, double tStart, int nSteps, bool showRewalked, bool withUi)
{
    App app;
    app.timeFull = false;
    printGeometry(app.tet, app.cfg);
    polyscope::options::usePrefsFile = false;   // fixed window size, nothing written to the working directory
    initPolyscope();
    polyscope::view::setWindowSize(1280, 960);
    app.registerTet();
    app.t = float(nSteps > 0 ? tStart : t);
    app.rebuild("start");
    for (int i = 1; i <= nSteps; i++) {
        app.t = float(tStart + (t - tStart)*double(i)/double(nSteps));
        app.setFront(app.t);
    }
    if (showRewalked) app.showRewalked();
    app.setCamera();
    polyscope::ScreenshotOptions opt;
    opt.transparentBackground = false;
    if (withUi) {   // also draw Polyscope's panels and the control panel
        polyscope::state::userCallback = [&app] { app.ui(); };
        opt.includeUI = true;
    }
    polyscope::screenshot(file, opt);
    std::printf("screenshot written: %s\n", file.c_str());
    return 0;
}

// --gate coupling --view: the union U = A u B of the coupling gate with the reference H, the decomposed H and their
// difference at the sample points (one point cloud; the F points carry the reference's trace and the glued trace),
// the crack, and a slice through the crack plane: the sample points within kGateSliceHalfWidth of it, one cloud per
// side (the eps-tube around the crack is empty: no sample points there). The values are the stored per-point ones of
// the check's main run (gateView): nothing is re-solved. With a file: one screenshot of the union and one of the
// slice per g = 1 field (<file stem>_<field>.png, <file stem>_slice_<field>.png), else the GUI.
constexpr double kGateSliceHalfWidth = 0.05;

std::string slug(const std::string& s)
{
    std::string o;
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c))) o += char(std::tolower(static_cast<unsigned char>(c)));
        else if (!o.empty() && o.back() != '_') o += '_';
    }
    while (!o.empty() && o.back() == '_') o.pop_back();
    return o;
}

// --gate coupling --view --check 4: one nodal ring of check 4, computed in the background (gateRingView) at the ring,
// crack, sample density and walk counts set in the panel, so the window stays live. Shown: the 8 tets, the crack, the
// 7 vertices, the sample points with every field as a quantity, and a slice through the crack plane (the points within
// kGateSliceHalfWidth of it, one cloud per side).
void showRingGeometry(const RingView& g)
{
    polyscope::removeAllStructures();
    std::vector<std::array<size_t, 4>> cells = g.tets;
    polyscope::VolumeMesh* vm = polyscope::registerTetMesh("ring (8 tets)", toArrays(std::vector<Vec3>(g.verts.begin(), g.verts.end())), cells);
    vm->setTransparency(0.15f);
    vm->setColor({0.75f, 0.75f, 0.80f});
    vm->setEdgeWidth(1.0);
    vm->addCellScalarQuantity("cut by the crack (1) or not (0)", std::vector<double>(g.tetCut.begin(), g.tetCut.end()));
    std::vector<std::array<int, 3>> tris;
    for (int k = 1; k + 1 < int(g.crack.size()); k++) tris.push_back({0, k, k + 1});
    polyscope::SurfaceMesh* cm = g.crackF.empty() ? polyscope::registerSurfaceMesh("crack", toArrays(g.crack), tris)
                                                  : polyscope::registerSurfaceMesh("crack", toArrays(g.crackV), g.crackF);
    cm->setSurfaceColor({0.95f, 0.70f, 0.15f});
    cm->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
    cm->setTransparency(0.5f);
    polyscope::PointCloud* pv = polyscope::registerPointCloud("vertices (v0 = the node)", toArrays(std::vector<Vec3>(g.verts.begin(), g.verts.end())));
    pv->setPointRadius(0.010);
    pv->setPointColor({0.1f, 0.1f, 0.1f});
}

void showRingFields(const RingView& v)
{
    std::vector<size_t> side[2];   // [0]: +g side of the crack plane, [1]: -g side, within kGateSliceHalfWidth
    for (size_t i = 0; i < v.X.size(); i++) {
        const double s = v.n.dot(v.X[i]) - v.d;
        if (s > 0 && s < kGateSliceHalfWidth) side[0].push_back(i);
        if (s < 0 && s > -kGateSliceHalfWidth) side[1].push_back(i);
    }
    polyscope::PointCloud* pc = polyscope::registerPointCloud("sample points", toArrays(v.X));
    pc->setPointRadius(0.003);
    pc->addScalarQuantity("tet", std::vector<double>(v.tetOf.begin(), v.tetOf.end()));
    polyscope::PointCloud* ps[2] = {nullptr, nullptr};
    const char* sideName[2] = {"slice, +g side of the crack plane", "slice, -g side of the crack plane"};
    for (int k = 0; k < 2; k++) {
        if (side[k].empty()) continue;
        std::vector<Vec3> X;
        for (size_t i : side[k]) X.push_back(v.X[i]);
        ps[k] = polyscope::registerPointCloud(sideName[k], toArrays(X));
        ps[k]->setPointRadius(0.005);
        ps[k]->setEnabled(false);
    }
    for (const GateView::Field& f : v.fields) {
        double r = 1.0;   // H: [-1, 1]; a difference: symmetric, its 99th percentile of |value|
        if (f.difference) {
            std::vector<double> a;
            for (double x : f.vol) a.push_back(std::fabs(x));
            std::nth_element(a.begin(), a.begin() + a.size()*99/100, a.end());
            r = std::max(1e-3, a[a.size()*99/100]);
        }
        auto add = [&](polyscope::PointCloud* c, const std::vector<double>& vals) {
            auto* q = c->addScalarQuantity(f.name, vals);
            q->setColorMap("coolwarm");
            q->setMapRange({-r, r});
            if (&f == &v.fields.front()) q->setEnabled(true);
        };
        add(pc, f.vol);
        for (int k = 0; k < 2; k++) {
            if (!ps[k]) continue;
            std::vector<double> vals;
            for (size_t i : side[k]) vals.push_back(f.vol[i]);
            add(ps[k], vals);
        }
    }
}

// One vertex's basis function (RingView::basis) as its own cloud, on the same coolwarm [-1, 1] as H (red 1, white 0),
// at the points of the tets around the vertex (supportOnly) or at every point; the vertex itself marked. k < 0: off.
const char* const kRingVertexName[7] = {"v0 (the node)", "v1 (+x)", "v2 (-x)", "v3 (+y)", "v4 (-y)", "v5 (+z)", "v6 (-z)"};
void showRingBasis(const RingView& v, int k, bool supportOnly)
{
    polyscope::removePointCloud("basis function", false);
    polyscope::removePointCloud("basis vertex", false);
    if (polyscope::hasPointCloud("sample points")) polyscope::getPointCloud("sample points")->setEnabled(k < 0);
    if (k < 0 || v.X.empty()) return;
    std::vector<Vec3> X;
    std::vector<double> vals;
    for (size_t i = 0; i < v.X.size(); i++)
        if (!supportOnly || v.inSupport[size_t(k)][i]) {
            X.push_back(v.X[i]);
            vals.push_back(v.basis[size_t(k)][i]);
        }
    polyscope::PointCloud* pc = polyscope::registerPointCloud("basis function", toArrays(X));
    pc->setPointRadius(0.003);
    auto* q = pc->addScalarQuantity(std::string("basis of ") + kRingVertexName[k], vals);
    q->setColorMap("coolwarm");
    q->setMapRange({-1.0, 1.0});
    q->setEnabled(true);
    polyscope::PointCloud* pv = polyscope::registerPointCloud("basis vertex", toArrays({v.verts[size_t(k)]}));
    pv->setPointRadius(0.018);
    pv->setPointColor({0.85f, 0.1f, 0.1f});
}

// The ring view's crack list: check 4's K1, K2, K3 (on ring 1 or 2), then --ring-stress cases (each on its own ring, by the
// fast method, the GT from its gate_cache/ref_stress_* cache)
struct RingCrackChoice {
    const char* group;   // a heading above this entry, or nullptr
    const char* label;
    int crack;           // check 4's crack (0, 1, 2), or -1: the --ring-stress case id
    const char* id;
    int ring;            // a stress case's ring (check 4's: 0, the panel's ring 1 or 2)
};
const RingCrackChoice kRingCracks[] = {
    {"flat (check 4)", "K1: flat crack", 0, "", 0},
    {nullptr, "K2: flat crack", 1, "", 0},
    {nullptr, "K3: flat cut through", 2, "", 0},
    {"curved (ring 2)", "C1: flat crack, gently curved front", -1, "C1", 2},
    {nullptr, "C2: flat crack, strongly curved front", -1, "C2", 2},
    {nullptr, "D1: curved crack (sphere, radius 2)", -1, "D1", 2},
    {nullptr, "D2: strongly curved crack (bowl, radius 0.6)", -1, "D2", 2},
    {nullptr, "D4: wavy crack (sine wave, wavelength 1)", -1, "D4", 2},
    {nullptr, "D5: strongly wavy crack (sine wave, wavelength 0.6)", -1, "D5", 2},
    {nullptr, "D3: kinked crack (40 degrees)", -1, "D3", 2},
    {nullptr, "B1: curved cut through (cylinder)", -1, "B1", 2},
    {"multiple cracks (ring 2)", "E1: two parallel cracks (g = 1 and 0.5)", -1, "E1", 2},
    {nullptr, "E2: two small cracks in one tet", -1, "E2", 2},
    {nullptr, "E3: two crossing cracks (X)", -1, "E3", 2},
    {"other placements", "B2: flat cut through, 0.001 from the node", -1, "B2", 2},
    {nullptr, "F1: front through the node", -1, "F1", 2},
    {nullptr, "F2: crack 0.001 from an interior face", -1, "F2", 2},
    {nullptr, "F2b: crack 0.03 from an interior face", -1, "F2b", 2},
    {nullptr, "F3: tiny crack inside one tet", -1, "F3", 2},
    {nullptr, "F4: front 0.03 from the outer boundary", -1, "F4", 2},
    {nullptr, "G1a: K1 on the distorted ring 3", -1, "G1a", 3},
    {nullptr, "G1b: sphere crack on the distorted ring 3", -1, "G1b", 3},
};

int runRingView(const RingViewParams& p0)
{
    initPolyscope();
    RingViewParams ui = p0;   // the panel's setting
    // the picked stress case's GT cache density (stressCachedPoints: -1 none needed, 0 no cache)
    int stressPts = ui.stressCase.empty() ? 0 : stressCachedPoints(ui.stressCase);
    auto crackLabel = [&ui]() -> std::string {
        for (const RingCrackChoice& c : kRingCracks)
            if (ui.stressCase.empty() ? c.crack == ui.crack : ui.stressCase == c.id) return c.label;
        return "stress case " + ui.stressCase;   // a case off the list (--case)
    };
    auto geometryOf = [](RingViewParams p) {
        p.geometryOnly = true;
        return gateRingView(p);
    };
    RingView v = geometryOf(ui);   // the last result (geometry only until the first run is done)
    bool have = false;
    std::string error;
    std::atomic<int> progress{0};   // tets done (8 per run)
    auto launch = [&progress](const RingViewParams& q) {
        progress = 0;
        return std::async(std::launch::async, [q, &progress] { return gateRingView(q, &progress); });
    };
    std::future<RingView> job = launch(ui);
    RingViewParams running = ui;   // the setting of the run in progress or shown
    auto started = std::chrono::steady_clock::now();
    // what the balls show: -2 the ground truth (the whole-ring solve, g = 1), -1 the glued H (Level C, physical glue,
    // g = 1), 0..6 that vertex's basis function, -3 another field picked in the list
    int show = -2;
    bool supportOnly = true;
    int fieldIdx = 0;           // the field shown on the sample points and the slices (RingView::fields)
    int gtIdx = 0, gluedIdx = 0;
    auto showField = [&](int idx) {
        if (idx < 0 || idx >= int(v.fields.size())) return;
        for (const char* c : {"sample points", "slice, +g side of the crack plane", "slice, -g side of the crack plane"})
            if (polyscope::hasPointCloud(c)) polyscope::getPointCloud(c)->getQuantity(v.fields[size_t(idx)].name)->setEnabled(true);
    };
    auto fieldLabel = [](const std::string& name) {   // what each field is, in the panel's selector
        if (name.rfind("GT", 0) == 0) return name + "   = one solve of the whole ring";
        if (name.rfind("glued H", 0) == 0) return name + "   = local solves + coupling";
        return name;
    };
    auto applyShow = [&] {
        if (show >= 0) {
            showRingBasis(v, show, supportOnly);
            return;
        }
        showRingBasis(v, -1, supportOnly);   // no basis cloud: the sample points and their field
        if (show == -2) fieldIdx = gtIdx;
        if (show == -1) fieldIdx = gluedIdx;
        showField(fieldIdx);
    };
    showRingGeometry(v);
    polyscope::view::lookAt(glm::vec3(2.2f, -3.2f, 1.8f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    std::printf("%s\ncomputing in the background (progress below)...\n", v.title.c_str());

    polyscope::state::userCallback = [&] {
        if (job.valid() && job.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
            try {
                v = job.get();
                have = true;
                error.clear();
                for (const std::string& s : v.report) std::printf("  %s\n", s.c_str());
                std::fflush(stdout);
                showRingGeometry(v);
                showRingFields(v);
                for (int i = int(v.fields.size()) - 1; i >= 0; i--) {   // the first match of each (fields are in order)
                    const std::string& n = v.fields[size_t(i)].name;
                    const bool g1 = n.size() >= 7 && n.compare(n.size() - 7, 7, "(g = 1)") == 0;
                    if (g1 && n.rfind("GT", 0) == 0) gtIdx = i;
                    if (g1 && n.rfind("glued H", 0) == 0) gluedIdx = i;
                }
                fieldIdx = std::min(fieldIdx, int(v.fields.size()) - 1);
                applyShow();
            } catch (const std::exception& e) {
                error = e.what();
            }
        }
        const bool busy = job.valid();
        const ImGuiSliderFlags logClamp = ImGuiSliderFlags_AlwaysClamp | ImGuiSliderFlags_Logarithmic;
        ImGui::PushItemWidth(180);
        ImGui::TextWrapped("%s", v.title.c_str());
        ImGui::Separator();
        ImGui::TextUnformatted("Ring and crack (applied by Recompute)");
        const bool stress = !ui.stressCase.empty();
        ImGui::PushItemWidth(-1.0f);
        if (ImGui::BeginCombo("##crack", crackLabel().c_str())) {
            for (const RingCrackChoice& c : kRingCracks) {
                if (c.group) ImGui::SeparatorText(c.group);
                const bool on = stress ? ui.stressCase == c.id : c.crack == ui.crack;
                if (ImGui::Selectable(c.label, on) && !on) {
                    ui.stressCase = c.id;
                    if (c.crack >= 0) {
                        ui.crack = c.crack;
                        if (stress) ui.pointsPerTet = 4096;   // check 4's density: its reference loads from gate_cache
                    } else {
                        stressPts = stressCachedPoints(ui.stressCase);
                        ui.pointsPerTet = stressPts > 0 ? stressPts : 2048;   // the density of the case's cached GT
                        if (!stress) {   // stress cases run the fast method only
                            ui.fast = true;
                            ui.walksTet = 32;
                            ui.walksR = 16;
                        }
                    }
                }
            }
            ImGui::EndCombo();
        }
        ImGui::PopItemWidth();
        if (!stress) {
            ImGui::RadioButton("ring 1: congruent (node at 0)", &ui.ring, 1);
            ImGui::SameLine();
            ImGui::RadioButton("ring 2: real mesh (node moved)", &ui.ring, 2);
        } else {
            int ring = 2;
            for (const RingCrackChoice& c : kRingCracks)
                if (ui.stressCase == c.id) ring = c.ring;
            ImGui::TextWrapped("on ring %d (%s); g = 1; the fast method only", ring, ring == 3 ? "distorted" : "real mesh, node moved");
        }
        ImGui::Separator();
        ImGui::TextUnformatted("Sample density and walks");
        ImGui::SliderInt("sample points per tet", &ui.pointsPerTet, 128, 16384, "%d", logClamp);
        if (!ui.fast) {
            ImGui::SliderInt("walks per point, per-tet runs", &ui.walksTet, 64, 2560, "%d", logClamp);
        } else {   // the fast method: the per-tet slider is Q's walks, and R has its own
            ImGui::SliderInt("walks per point, Q runs (fast)", &ui.walksTet, 16, 2560, "%d", logClamp);
            ImGui::SliderInt("walks per point, R runs (fast)", &ui.walksR, 16, 2560, "%d", logClamp);
        }
        ImGui::SliderInt(stress ? "walks per point, GT (when not cached)" : "walks per point, ring reference", &ui.walksRef, 256, 25600,
                         "%d", logClamp);
        if (!stress && ImGui::Button("check 4's own setting (4096 / 2560 / 25600)")) {
            ui.pointsPerTet = 4096;
            ui.walksTet = 2560;
            ui.walksRef = 25600;
        }
        if (stress && stressPts > 0 && ui.pointsPerTet != stressPts &&
            ImGui::Button(("the cached GT's density (" + std::to_string(stressPts) + " points per tet)").c_str()))
            ui.pointsPerTet = stressPts;
        ImGui::BeginDisabled(stress);
        ImGui::Checkbox("fast method (per-tet control variates + lean walks)", &ui.fast);
        ImGui::EndDisabled();
        if (ui.fast) {
            ImGui::SameLine();
            if (ImGui::Button("fast defaults (Q 32 / R 16)")) {
                ui.walksTet = 32;
                ui.walksR = 16;
            }
        }
        if (stress) {
            // the fast method's rough cost as below (curved cracks: a few times more, the solid angle at every R landing); the
            // GT: exact, the case's cache at its density, or solved here at the GT walks
            const bool cached = stressPts < 0 || stressPts == ui.pointsPerTet;
            const double pts = 8.0*ui.pointsPerTet, tTet = pts*(ui.walksTet + ui.walksR)*0.2e-6, tRef = cached ? 0.0 : pts*ui.walksRef*0.7e-6;
            const std::string gt = stressPts < 0 ? "exact or check 4's cache" : cached ? "loaded from gate_cache"
                                 : "solved here, not cached: " + std::to_string(int(tRef)) + " s" +
                                       (stressPts > 0 ? " (cached at " + std::to_string(stressPts) + " points per tet)" : "");
            ImGui::TextWrapped("estimated: about %.1f s (our runs %.1f s, longer for curved cracks; GT %s)", tTet + tRef, tTet, gt.c_str());
        } else if (!ui.fast) {
            // rough cost on 16 threads: ~0.7 us per walk on the ring, ~0.45 us per walk in a tet (two runs per tet)
            const bool asCheck = ui.pointsPerTet == 4096 && ui.walksTet == 2560 && ui.walksRef == 25600;
            const double pts = 8.0*ui.pointsPerTet, tTet = pts*ui.walksTet*2*0.45e-6, tRef = asCheck ? 0.0 : pts*ui.walksRef*0.7e-6;
            ImGui::TextWrapped("estimated: about %.0f s (per-tet runs %.0f s, reference %s)", tTet + tRef, tTet,
                               asCheck ? "loaded from gate_cache" : (std::to_string(int(tRef)) + " s").c_str());
        } else {
            // rough cost of the fast method (estimate): ~0.2 us per walk in a tet; the reference is loaded from gate_cache at
            // 4096 points and 25600 walks whatever the per-tet walks are
            const bool refCached = ui.pointsPerTet == 4096 && ui.walksRef == 25600;
            const double pts = 8.0*ui.pointsPerTet, tTet = pts*(ui.walksTet + ui.walksR)*0.2e-6, tRef = refCached ? 0.0 : pts*ui.walksRef*0.7e-6;
            ImGui::TextWrapped("estimated: about %.1f s (fast per-tet runs %.1f s, reference %s). Each tet still solves in its own xi "
                               "frame: P_k = lambda_k + Q_k and H0 = J_xi + R, with Q and R the only walks.", tTet + tRef, tTet,
                               refCached ? "loaded from gate_cache" : (std::to_string(int(tRef)) + " s").c_str());
        }
        if (busy && !running.stressCase.empty()) {   // no per-tet progress: the GT loads, then all tets run together
            ImGui::TextWrapped("solving %s... %.0f s (details in the console)", running.stressCase.c_str(),
                               std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count());
        } else if (busy) {
            const int done = progress;
            const bool asCheck = ui.fast ? ui.pointsPerTet == 4096 && ui.walksRef == 25600
                                         : ui.pointsPerTet == 4096 && ui.walksTet == 2560 && ui.walksRef == 25600;
            ImGui::TextWrapped("solving... %.0f s: per-tet runs %d / 8 done%s (details in the console)",
                               std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count(), done,
                               done == 0 && !asCheck ? " (the ring reference runs first)" : "");
            ImGui::ProgressBar(float(done)/8.0f, ImVec2(-1.0f, 0.0f));
        } else if (ImGui::Button("Recompute")) {
            v = geometryOf(ui);
            have = false;
            showRingGeometry(v);
            job = launch(ui);
            running = ui;
            started = std::chrono::steady_clock::now();
        }
        if (!error.empty()) ImGui::TextWrapped("error: %s", error.c_str());
        ImGui::Separator();
        if (have) {
            ImGui::TextUnformatted("Show (red = +1, white = 0, blue = -1), g = 1");
            const int before = show;
            const bool supportBefore = supportOnly;
            ImGui::RadioButton("GT: H from one solve of the whole ring", &show, -2);
            ImGui::RadioButton("glued H: local solves + coupling", &show, -1);
            ImGui::TextUnformatted("one vertex's basis function (used inside the coupling):");
            for (int k = 0; k < 7; k++) {
                if (k != 0 && k != 4) ImGui::SameLine();
                ImGui::RadioButton(kRingVertexName[k], &show, k);
            }
            if (ImGui::Button("previous vertex")) show = show <= 0 ? 6 : show - 1;
            ImGui::SameLine();
            if (ImGui::Button("next vertex")) show = show < 0 || show >= 6 ? 0 : show + 1;
            ImGui::SameLine();
            ImGui::Checkbox("only its tets", &supportOnly);
            if (show != before || (show >= 0 && supportOnly != supportBefore)) applyShow();
            ImGui::TextWrapped("glued H = H0 + sum over vertices of (weight c_v x basis function of v). A vertex's basis "
                               "function: in each tet around it, that tet's local solve with crack 0 and faces = the vertex's "
                               "hat; 0 in the other tets.");
            ImGui::TextUnformatted(v.crackF.empty() ? "more fields (the error, H0, the ramp g):" : "more fields (the error, H0, each face function J):");
            ImGui::PushItemWidth(-1.0f);
            if (ImGui::BeginCombo("##field", show == -3 ? fieldLabel(v.fields[size_t(fieldIdx)].name).c_str() : "(pick one)")) {
                for (int i = 0; i < int(v.fields.size()); i++)
                    if (ImGui::Selectable(fieldLabel(v.fields[size_t(i)].name).c_str(), show == -3 && i == fieldIdx)) {
                        fieldIdx = i;
                        show = -3;
                        applyShow();
                    }
                ImGui::EndCombo();
            }
            ImGui::PopItemWidth();
            ImGui::Separator();
            for (const std::string& s : v.report) ImGui::TextWrapped("%s", s.c_str());
            if (v.crackF.empty())
                ImGui::TextWrapped("The ramp g is computed in ring 1 only. The error field has its own colour range (shown in the "
                                   "Structures panel). The two slice clouds (off by default) hold the points within %.2f of the "
                                   "crack plane, one per side.", kGateSliceHalfWidth);
            else   // a stress case
                ImGui::TextWrapped("g = 1. 'glued H' uses the case's main face function J (the first in the list); the list also has "
                                   "the others. The error field has its own colour range (shown in the Structures panel).%s",
                                   v.n.isZero() ? "" : " The two slice clouds (off by default) hold the points near the first crack's plane.");
        }
        ImGui::PopItemWidth();
    };
    polyscope::show();
    if (job.valid()) {   // a run still going when the window closes finishes, and its result goes to the console
        try {
            for (const std::string& s : job.get().report) std::printf("  %s\n", s.c_str());
        } catch (const std::exception& e) {
            std::fprintf(stderr, "error: %s\n", e.what());
            return 1;
        }
    }
    if (!error.empty()) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        return 1;
    }
    return 0;
}

int runGateView(int check, int config, const std::string& shotFile, bool withUi, bool tipTrace)
{
    const GateView v = gateView(check, config, tipTrace);
    const bool shots = !shotFile.empty();
    if (shots) polyscope::options::usePrefsFile = false;   // fixed window size, nothing written to the working directory
    initPolyscope();
    if (shots) polyscope::view::setWindowSize(1280, 960);

    std::vector<std::array<size_t, 4>> cells(v.tets.begin(), v.tets.end());
    polyscope::VolumeMesh* vm = polyscope::registerTetMesh("A u B", toArrays(std::vector<Vec3>(v.verts.begin(), v.verts.end())), cells);
    vm->setTransparency(0.2f);
    vm->setColor({0.75f, 0.75f, 0.80f});
    vm->setEdgeWidth(1.0);
    std::vector<std::array<int, 3>> tris;
    for (int k = 1; k + 1 < int(v.crack.size()); k++) tris.push_back({0, k, k + 1});
    polyscope::SurfaceMesh* cm = polyscope::registerSurfaceMesh("crack (U)", toArrays(v.crack), tris);
    cm->setSurfaceColor({0.95f, 0.70f, 0.15f});
    cm->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
    cm->setTransparency(0.5f);

    // the clouds: all sample points of A and B, the F points, and the slice on either side of the crack plane
    std::vector<size_t> sideIdx[2];   // [0]: +g side, [1]: -g side
    for (size_t i = 0; i < v.X.size(); i++) {
        const double s = v.n.dot(v.X[i]) - v.d;
        if (s > 0 && s < kGateSliceHalfWidth) sideIdx[0].push_back(i);
        if (s < 0 && s > -kGateSliceHalfWidth) sideIdx[1].push_back(i);
    }
    auto pick = [](const std::vector<double>& a, const std::vector<size_t>& idx) {
        std::vector<double> o;
        for (size_t i : idx) o.push_back(a[i]);
        return o;
    };
    auto pickPts = [](const std::vector<Vec3>& a, const std::vector<size_t>& idx) {
        std::vector<Vec3> o;
        for (size_t i : idx) o.push_back(a[i]);
        return o;
    };
    polyscope::PointCloud* pcU = polyscope::registerPointCloud("A u B sample points", toArrays(v.X));
    polyscope::PointCloud* pcF = polyscope::registerPointCloud("F points (trace)", toArrays(v.XF));
    polyscope::PointCloud* pcS[2] = {nullptr, nullptr};
    const char* sideName[2] = {"slice, +g side of the crack plane", "slice, -g side of the crack plane"};
    for (int k = 0; k < 2; k++)
        if (!sideIdx[k].empty()) pcS[k] = polyscope::registerPointCloud(sideName[k], toArrays(pickPts(v.X, sideIdx[k])));
    pcU->setPointRadius(0.004);
    pcF->setPointRadius(0.004);
    for (auto* p : pcS)
        if (p) p->setPointRadius(0.007);
    std::vector<double> range;   // per field
    for (const GateView::Field& f : v.fields) {
        double r = 1.0;   // H: [-1, 1]; a difference: symmetric, its 99th percentile of |value| over the volume points
        if (f.difference) {
            std::vector<double> a;
            for (double x : f.vol) a.push_back(std::fabs(x));
            std::nth_element(a.begin(), a.begin() + a.size()*99/100, a.end());
            r = std::max(1e-3, a[a.size()*99/100]);
        }
        range.push_back(r);
        auto add = [&](polyscope::PointCloud* pc, const std::vector<double>& vals) {
            auto* q = pc->addScalarQuantity(f.name, vals);
            q->setColorMap("coolwarm");
            q->setMapRange({-r, r});
        };
        add(pcU, f.vol);
        add(pcF, f.face);
        for (int k = 0; k < 2; k++)
            if (pcS[k]) add(pcS[k], pick(f.vol, sideIdx[k]));
    }

    // cameras: the union from the +x, -y, +z side; the slice face-on from the +g side of the crack plane
    const Vec3 cU(0.25, 0.25, 0.0);
    auto lookAt = [](const Vec3& eye, const Vec3& target, const Vec3& up) {
        polyscope::view::lookAt(glm::vec3(eye.x(), eye.y(), eye.z()), glm::vec3(target.x(), target.y(), target.z()),
                                glm::vec3(up.x(), up.y(), up.z()));
    };
    const Vec3 up = std::fabs(v.n.z()) < 0.9 ? Vec3(0, 0, 1) : Vec3(0, 1, 0);
    auto unionCamera = [&] { lookAt(cU + 2.8*Vec3(0.55, -1.0, 0.4).normalized(), cU, Vec3(0, 0, 1)); };
    // the slice camera: on the normal through the slice's centroid, at a distance that frames its farthest point
    // (the default 45 degree field of view shows +-0.41 of the distance)
    Vec3 cS = v.n*v.d;
    double rS = 0.5;
    if (!sideIdx[0].empty()) {
        cS = Vec3::Zero();
        for (size_t i : sideIdx[0]) cS += v.X[i]/double(sideIdx[0].size());
        rS = 0;
        for (size_t i : sideIdx[0]) rS = std::max(rS, (v.X[i] - cS).norm());
    }
    auto sliceCamera = [&] { lookAt(cS + std::max(0.5, 2.8*rS)*v.n, cS, up); };
    auto enableOnly = [&](polyscope::PointCloud* on1, polyscope::PointCloud* on2) {
        for (auto* p : {pcU, pcF, pcS[0], pcS[1]})
            if (p) p->setEnabled(p == on1 || p == on2);
    };
    std::printf("%s\n", v.title.c_str());

    if (!shots) {
        enableOnly(pcU, pcF);
        pcU->getQuantity(v.fields.front().name)->setEnabled(true);
        pcF->getQuantity(v.fields.front().name)->setEnabled(true);
        for (auto* p : pcS)
            if (p) p->getQuantity(v.fields.front().name)->setEnabled(true);
        unionCamera();
        polyscope::state::userCallback = [&v] {
            ImGui::TextWrapped("%s", v.title.c_str());
            ImGui::TextWrapped("Values at the sample points of the check's main run (no re-solve). Pick a field in a cloud's "
                               "quantity list; the slice clouds hold the points within %.2f of the crack plane, one per side.",
                               kGateSliceHalfWidth);
        };
        polyscope::show();
        return 0;
    }
    const std::filesystem::path base(shotFile);
    const std::string stem = (base.parent_path()/base.stem()).string();
    polyscope::ScreenshotOptions opt;
    opt.transparentBackground = false;
    opt.includeUI = withUi;
    // the view's section of ./gate_report.md: the images with their colour ranges, and the errors of the fields shown
    char date[64];
    {
        const std::time_t now = std::time(nullptr);
        std::tm tmNow{};
#ifdef _WIN32
        localtime_s(&tmNow, &now);
#else
        localtime_r(&now, &tmNow);
#endif
        std::strftime(date, sizeof date, "%Y-%m-%d %H:%M:%S", &tmNow);
    }
    auto rel = [](const std::string& file) {   // the image's path from the report's folder (the working directory)
        std::error_code ec;
        const std::filesystem::path r = std::filesystem::relative(file, std::filesystem::current_path(), ec);
        return (ec || r.empty() ? std::filesystem::path(file) : r).generic_string();
    };
    char buf[512];
    std::snprintf(buf, sizeof buf, "--check %d --config %d", check, config);
    const std::string heading = v.title.rfind("coupling gate, ", 0) == 0 ? v.title.substr(15) : v.title;
    std::string md = "## Polyscope view: " + heading + "\n\nWritten by `partial_cut --gate coupling --view " + buf + " --screenshot " +
                     base.generic_string() + "` on " + date + ". Values at the sample points of the check's main runs: the "
                     "reference loaded from gate_cache/, the decomposed field from the check's recipe runs with the check's own "
                     "seeds, so these are the fields behind the tables above (nothing re-solved). Left: the union A u B (a "
                     "transparent tet mesh) with U's crack (orange), the sample points of A and B and the F points (on F the "
                     "decomposed field is its glued trace), seen from the +x, -y, +z side. Right: a slice through the crack plane, "
                     "the sample points within 0.05 of it on its +g side, seen face-on from that side (the eps-tube around the "
                     "crack holds no points). Colour map coolwarm (blue negative, red positive); H on [-1, 1], a difference on a "
                     "symmetric range at its 99th percentile of abs(value). Colour ranges are in units of H; the errors quoted "
                     "with them are divided by 2 g_max = 2, as in the tables. Without --screenshot the same command opens the "
                     "GUI with every field (both g modes; check 3 at both levels).\n\n";
    for (size_t k = 0; k < v.fields.size(); k++) {
        const GateView::Field& f = v.fields[k];
        if (!f.shot) continue;
        std::string files[2];
        for (int view = 0; view < 2; view++) {
            if (view == 1 && !pcS[0]) continue;
            for (auto* p : {pcU, pcF, pcS[0], pcS[1]})
                if (p) p->getQuantity(f.name)->setEnabled(true);
            if (view == 0) {
                enableOnly(pcU, pcF);
                unionCamera();
            } else {
                enableOnly(pcS[0], nullptr);
                sliceCamera();
            }
            const std::string file = stem + (view ? "_slice_" : "_") + slug(f.name) + ".png";
            polyscope::screenshot(file, opt);
            std::printf("screenshot written: %s\n", file.c_str());
            files[view] = rel(file);
        }
        std::snprintf(buf, sizeof buf, "colour range [%.3g, %.3g]", -range[k], range[k]);
        md += "**" + f.name + "**, " + buf;
        if (f.difference) {
            std::snprintf(buf, sizeof buf, "; over A u B, / 2 g_max: H RMS %.2e, max %.2e", f.rms, f.max);
            md += buf;
        }
        md += ".\n\n";
        for (const std::string& file : files)
            if (!file.empty()) md += "<img src=\"" + file + "\" width=\"49%\" alt=\"" + f.name + "\"> ";
        md += "\n\n";
    }
    std::snprintf(buf, sizeof buf, "view_check%d_config%d", check, config);
    if (writeGateReportSection("gate_report.md", buf, md)) std::printf("gate_report.md: section %s rewritten\n", buf);
    else std::printf("gate_report.md not found in the working directory: no report section written\n");
    return 0;
}

void printUsage()
{
    std::printf("usage:\n"
                "  partial_cut                                  Polyscope GUI (own WoSt loop + walk cache)\n"
                "  partial_cut --validate                       first gate: 4 checks (Zombie), rewrite its section of ./results.txt\n"
                "  partial_cut --bench                          time one full default Zombie solve\n"
                "  partial_cut --validate-reuse [--quick]       reuse validation (V0, R1-R4, checks 1-4), rewrite its section\n"
                "                                               of ./results.txt (--quick: small R1-R4 run, writes nothing)\n"
                "  partial_cut --bench-reuse                    one default incremental sweep with full-solve timings\n"
                "  partial_cut --compare-options [--quick] [--threads T] [--repeat R]\n"
                "                                               reuse options 1, 2a, 2b, 3 and the control 3q: equivalence\n"
                "                                               checks, timings and memory; writes ./results_options.txt\n"
                "                                               (--quick: small run, ./results_options_quick.txt); --threads\n"
                "                                               limits TBB (adds _threadsT to the file name); --repeat R:\n"
                "                                               R samples per run (medians, ranges)\n"
                "    (add --dirichlet-faces to any of the three: experiment with absorbing tet faces, H = 0 there; the\n"
                "     validation then writes ./results_dirichlet_faces.txt, the comparison\n"
                "     ./results_options_dirichlet_faces.txt)\n"
                "  partial_cut --gate coupling [--check 0|1|2|3|4] [--quick]\n"
                "                                               the coupling gate (check 4: one nodal ring); rewrite its\n"
                "                                               sections of ./gate_report.md\n"
                "                                               and cache the reference solves in ./gate_cache/\n"
                "  partial_cut --gate coupling --view [--check 2|3] [--config i] [--screenshot out.png [--with-ui]]\n"
                "                                               Polyscope view of the union: reference H, decomposed H and\n"
                "                                               their difference at the sample points, and a slice through\n"
                "                                               the crack plane (check 3, i = 0..2: y_f = 0.15, 0.30, 0.45,\n"
                "                                               default 1; check 2, i = 0, 1: z_c = 0.60, 0.35, default 0);\n"
                "                                               --screenshot writes out_<field>.png and out_slice_<field>.png\n"
                "                                               and rewrites the view's section of ./gate_report.md\n"
                "  partial_cut --gate coupling --view --check 4 [--ring 1|2] [--config 0|1|2]\n"
                "                                               check 4's nodal ring (default ring 2, K1) with a panel to\n"
                "                                               change the ring, the crack (K1, K2, K3, or a --ring-stress\n"
                "                                               case: curved, kinked, several cracks, ...), the sample points\n"
                "                                               per tet and the walk counts; recomputed in the background\n"
                "    (add --fast: start with the panel's 'fast method' checkbox on, Q 32 / R 16 walks per point)\n"
                "    (add --case ID [--points P] [--walks-tet Q] [--walks-r R] [--walks-ref N]: start on that --ring-stress\n"
                "     case, at its cached GT's density; without a cache at P the GT is solved at N walks, not cached)\n"
                "  partial_cut --ring-bench [--ring 1|2] [--config 0|1|2] [--ours LIST] [--direct LIST] [--current LIST]\n"
                "              [--shell S] [--repeat R] [--seeds S] [--no-lean] [--no-fuse] [--tag NAME] [--quick]\n"
                "                                               check 4's ring: the fast per-tet method vs the direct ring\n"
                "                                               solve (and check 4's current path), timed and compared with\n"
                "                                               the cached reference (read-only); writes ring_perf/<tag>_rows.csv,\n"
                "                                               _summary.csv and _log.txt. LIST: ours nQ/nR[/B],...; direct and\n"
                "                                               current N[/B],...; none skips. --shell: stop distance / sqrt(2)\n"
                "  partial_cut --ring-bench --selftest            lean walks vs runWalk bit for bit, J_xi, skip decisions, z tests\n"
                "  partial_cut --ring-bench --view [--fast] [--ring r] [--config k] [--points P] [--walks-tet N]\n"
                "              [--walks-r N] [--walks-ref N]   the ring view's computation headless: its panel lines\n"
                "                                               (with --case ID: a --ring-stress case's view, headless)\n"
                "  partial_cut --ring-stress [--case all|p1|p2|p3|ID[,ID...]] [--ours \"Q/R[/B],...\"] [--direct \"N[/B],...\"]\n"
                "              [--points P] [--gt-walks N] [--seeds S] [--repeat R] [--zombie IDS|none] [--no-validate]\n"
                "              [--no-cache] [--eps E] [--tag NAME] [--quick]\n"
                "                                               the nodal ring under other crack types (complete, curved, kinked,\n"
                "                                               several cracks, placement stress, ring 3): our per-tet method vs\n"
                "                                               the GT on O; writes ring_stress/<tag>_rows.csv, _summary.csv,\n"
                "                                               _log.txt, _geometry.txt and new GT caches gate_cache/ref_stress_*\n"
                "                                               (--eps E: leave out the points within E sqrt(2) of a crack,\n"
                "                                               default 0.01 as check 4; --no-cache: GTs not read or written)\n"
                "                                               (ids: A1 A1s A2 A3 B1 B2 C1 C2 D1 D2 D3 E1 E2 E3 F1 F2 F3 F4 G1a G1b F2b D4 D5)\n"
                "  partial_cut --ring-stress --selftest           crack-mesh tree, flat mesh vs StripCut, every case's J and\n"
                "                                               control-variate checks, solid angle, determinism\n"
                "  partial_cut --check-xi                     H(xi) on the reference tet: cut pull-back and rest-tet\n"
                "                                               independence checks; rewrite their section of ./results.txt\n"
                "  partial_cut --screenshot file.png [--t T] [--t-start T0 --steps K] [--rewalked] [--with-ui]\n"
                "                                               front at T (default 0.6): one full solve, or K incremental\n"
                "                                               steps from T0 (default 0.2); --rewalked shows the\n"
                "                                               're-walked this step' scalar; --with-ui adds the panels\n");
}

} // namespace

int main(int argc, char** argv)
{
    bool validate = false, bench = false, shotUi = false, validateReuse = false, quick = false, benchReuse = false;
    bool shotRewalked = false, dirichletFaces = false, compareOptions = false, checkXi = false, gateViewMode = false;
    bool tipTrace = false;
    std::string shotFile, gate, ringPerf;
    int gateCheck = -1, gateConfig = -1, gateRing = -1;
    double tShot = Settings{}.t, tStartShot = 0.2;
    int shotSteps = 0, threads = 0, repeats = 1;
    // --ring-bench (coupling_gate.h runRingBench) and the ring view's --fast option
    bool ringBench = false, repeatGiven = false, ringFast = false;
    RingBenchParams benchP;
    // --ring-stress (coupling_gate.h runRingStress); it shares --ours, --direct, --points, --seeds, --repeat, --tag, --quick and
    // --selftest with --ring-bench
    bool ringStress = false, seedsGiven = false, pointsGiven = false, walksTetGiven = false, walksRGiven = false, walksRefGiven = false;
    RingStressParams stressP;
    RingViewParams ringView;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--validate") == 0) validate = true;
        else if (std::strcmp(argv[i], "--bench") == 0) bench = true;
        else if (std::strcmp(argv[i], "--validate-reuse") == 0) validateReuse = true;
        else if (std::strcmp(argv[i], "--quick") == 0) quick = true;
        else if (std::strcmp(argv[i], "--bench-reuse") == 0) benchReuse = true;
        else if (std::strcmp(argv[i], "--compare-options") == 0) compareOptions = true;
        else if (std::strcmp(argv[i], "--check-xi") == 0) checkXi = true;
        else if (std::strcmp(argv[i], "--gate") == 0 && i + 1 < argc) gate = argv[++i];
        else if (std::strcmp(argv[i], "--check") == 0 && i + 1 < argc) gateCheck = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--view") == 0) gateViewMode = true;
        else if (std::strcmp(argv[i], "--config") == 0 && i + 1 < argc) gateConfig = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--ring") == 0 && i + 1 < argc) gateRing = std::atoi(argv[++i]);
        else if (std::strcmp(argv[i], "--face-trace") == 0 && i + 1 < argc && std::strcmp(argv[i + 1], "tip") == 0) {
            tipTrace = true;
            i++;
        }
        else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) threads = std::max(0, std::atoi(argv[++i]));
        else if (std::strcmp(argv[i], "--repeat") == 0 && i + 1 < argc) {
            repeats = std::max(1, std::atoi(argv[++i]));
            repeatGiven = true;
        }
        else if (std::strcmp(argv[i], "--screenshot") == 0 && i + 1 < argc) shotFile = argv[++i];
        else if (std::strcmp(argv[i], "--t") == 0 && i + 1 < argc) tShot = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--t-start") == 0 && i + 1 < argc) tStartShot = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--steps") == 0 && i + 1 < argc) shotSteps = std::max(0, std::atoi(argv[++i]));
        else if (std::strcmp(argv[i], "--rewalked") == 0) shotRewalked = true;
        else if (std::strcmp(argv[i], "--dirichlet-faces") == 0) dirichletFaces = true;
        else if (std::strcmp(argv[i], "--with-ui") == 0) shotUi = true;
        else if (std::strcmp(argv[i], "--ring-perf") == 0 && i + 1 < argc) ringPerf = argv[++i];
        else if (std::strcmp(argv[i], "--ring-bench") == 0) ringBench = true;
        else if (std::strcmp(argv[i], "--selftest") == 0) benchP.selftest = true;
        else if (std::strcmp(argv[i], "--ours") == 0 && i + 1 < argc) stressP.ours = benchP.ours = argv[++i];
        else if (std::strcmp(argv[i], "--direct") == 0 && i + 1 < argc) stressP.direct = benchP.direct = argv[++i];
        else if (std::strcmp(argv[i], "--current") == 0 && i + 1 < argc) benchP.current = argv[++i];
        else if (std::strcmp(argv[i], "--shell") == 0 && i + 1 < argc) benchP.shellRel = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--seeds") == 0 && i + 1 < argc) {
            benchP.seeds = std::max(1, std::atoi(argv[++i]));
            seedsGiven = true;
        }
        else if (std::strcmp(argv[i], "--tag") == 0 && i + 1 < argc) stressP.tag = benchP.tag = argv[++i];
        else if (std::strcmp(argv[i], "--ring-stress") == 0) ringStress = true;
        else if (std::strcmp(argv[i], "--case") == 0 && i + 1 < argc) stressP.cases = argv[++i];
        else if (std::strcmp(argv[i], "--gt-walks") == 0 && i + 1 < argc) stressP.gtWalks = std::max(64, std::atoi(argv[++i]));
        else if (std::strcmp(argv[i], "--zombie") == 0 && i + 1 < argc) stressP.zombie = argv[++i];
        else if (std::strcmp(argv[i], "--no-validate") == 0) stressP.validate = false;
        else if (std::strcmp(argv[i], "--no-cache") == 0) stressP.cache = false;
        else if (std::strcmp(argv[i], "--eps") == 0 && i + 1 < argc) stressP.epsRel = std::max(1e-6, std::atof(argv[++i]));
        else if (std::strcmp(argv[i], "--no-lean") == 0) benchP.lean = false;
        else if (std::strcmp(argv[i], "--no-fuse") == 0) benchP.fuse = false;
        else if (std::strcmp(argv[i], "--fast") == 0) ringFast = true;
        else if (std::strcmp(argv[i], "--points") == 0 && i + 1 < argc) {
            ringView.pointsPerTet = std::max(64, std::atoi(argv[++i]));
            pointsGiven = true;
        }
        else if (std::strcmp(argv[i], "--walks-tet") == 0 && i + 1 < argc) {
            ringView.walksTet = std::max(2, std::atoi(argv[++i]));
            walksTetGiven = true;
        }
        else if (std::strcmp(argv[i], "--walks-r") == 0 && i + 1 < argc) {
            ringView.walksR = std::max(2, std::atoi(argv[++i]));
            walksRGiven = true;
        }
        else if (std::strcmp(argv[i], "--walks-ref") == 0 && i + 1 < argc) {
            ringView.walksRef = std::max(2, std::atoi(argv[++i]));
            walksRefGiven = true;
        }
        else if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) { printUsage(); return 0; }
        else { std::fprintf(stderr, "unknown argument: %s\n", argv[i]); printUsage(); return 2; }
    }

    try {
        if (validate) {
            GeometryConfig cfg;
            Tet tet = referenceTet();
            printGeometry(tet, cfg);
            return runValidation(cfg, Settings{}, "results.txt") ? 0 : 1;
        }
        if (bench) return runBench();
        if (checkXi) return runXiCheck(Settings{}, "results.txt") ? 0 : 1;
        if (!ringPerf.empty())   // profiling only (coupling_gate.h runRingPerf)
            return runRingPerf(gateRing == 1 ? 1 : 2, gateConfig >= 0 && gateConfig <= 2 ? gateConfig : 0, ringPerf) ? 0 : 1;
        if (ringStress) {   // the ring under other crack types (coupling_gate.h runRingStress)
            stressP.quick = quick;
            stressP.selftest = benchP.selftest;
            if (repeatGiven) stressP.repeats = repeats;
            else if (quick) stressP.repeats = 1;
            if (seedsGiven) stressP.seeds = benchP.seeds;
            else if (quick) stressP.seeds = 1;
            if (pointsGiven) stressP.points = ringView.pointsPerTet;
            if (quick && stressP.zombie == "D2,E3") stressP.zombie = "none";
            return runRingStress(stressP) ? 0 : 1;
        }
        if (ringBench) {   // the fast ring method vs the direct ring solve (coupling_gate.h runRingBench)
            benchP.ring = gateRing == 1 || gateRing == 2 ? gateRing : 0;
            benchP.crack = gateConfig >= 0 && gateConfig <= 2 ? gateConfig : -1;
            benchP.quick = quick;
            if (repeatGiven) benchP.repeats = repeats;
            else if (quick) benchP.repeats = 1;   // --quick: one repeat unless --repeat says otherwise
            benchP.view = gateViewMode;
            if (gateViewMode && stressP.cases != "all") {   // a --ring-stress case's view headless: its lines and field hashes
                const int cp = stressCachedPoints(stressP.cases);   // default: the density of the case's cached GT
                const RingView v = stressRingView(stressP.cases, pointsGiven ? ringView.pointsPerTet : cp > 0 ? cp : 2048,
                                                  walksTetGiven ? ringView.walksTet : 32,
                                                  walksRGiven ? ringView.walksR : 16, walksRefGiven ? ringView.walksRef : 2560);
                std::printf("%s\n", v.title.c_str());
                for (const std::string& s : v.report) std::printf("  %s\n", s.c_str());
                for (const GateView::Field& f : v.fields) {
                    uint64_t h = 1469598103934665603ull;   // fnv1a of the values
                    const unsigned char* b = reinterpret_cast<const unsigned char*>(f.vol.data());
                    for (size_t i = 0; i < f.vol.size()*sizeof(double); i++) h = (h ^ b[i])*1099511628211ull;
                    std::printf("  field \"%s\": %zu values, fnv1a %016llx\n", f.name.c_str(), f.vol.size(), (unsigned long long)h);
                }
                return 0;
            }
            if (gateViewMode) {   // the ring view headless: ring 2, K1 unless told otherwise
                ringView.ring = benchP.ring ? benchP.ring : 2;
                ringView.crack = benchP.crack >= 0 ? benchP.crack : 0;
                ringView.fast = ringFast;
                benchP.viewParams = ringView;
            }
            return runRingBench(benchP) ? 0 : 1;
        }
        if (!gate.empty()) {
            if (gate != "coupling" || gateCheck < -1 || gateCheck > 4 || (gateViewMode && gateCheck == 4 && !shotFile.empty())) {
                std::fprintf(stderr, "usage: --gate coupling [--check 0|1|2|3|4] [--quick]\n"
                                     "       --gate coupling --view [--check 2|3] [--config i] [--screenshot out.png [--with-ui]]\n"
                                     "       --gate coupling --view --check 4 [--ring 1|2] [--config 0|1|2]\n");
                return 2;
            }
            if (gateViewMode && gateCheck == 4) {   // ring 2 (real mesh), K1 unless told otherwise; the panel changes both
                RingViewParams rp;
                if (gateRing == 1 || gateRing == 2) rp.ring = gateRing;
                if (gateConfig >= 0 && gateConfig <= 2) rp.crack = gateConfig;
                if (stressP.cases != "all") {   // --case ID: start on that --ring-stress case, at its cached GT's density
                    rp.stressCase = stressP.cases;
                    const int cp = stressCachedPoints(rp.stressCase);
                    rp.pointsPerTet = pointsGiven ? ringView.pointsPerTet : cp > 0 ? cp : 2048;
                    ringFast = true;
                }
                if (ringFast) {   // --fast: start with the fast method's checkbox on (and its default walks Q 32 / R 16)
                    rp.fast = true;
                    rp.walksTet = walksTetGiven ? ringView.walksTet : 32;
                    rp.walksR = walksRGiven ? ringView.walksR : 16;
                }
                if (walksRefGiven) rp.walksRef = ringView.walksRef;
                return runRingView(rp);
            }
            if (gateViewMode) {   // check 3 at y_f = 0.30 unless told otherwise; check 2 at z_c = 0.60
                const int c = gateCheck < 0 ? 3 : gateCheck;
                return runGateView(c, gateConfig >= 0 ? gateConfig : c == 3 ? 1 : 0, shotFile, shotUi, tipTrace);
            }
            return runCouplingGate(gateCheck, quick, "gate_report.md", tipTrace) ? 0 : 1;
        }
        // --dirichlet-faces (experiment): the own loop's tet faces absorb with H = 0; the reuse
        // validation then writes results_dirichlet_faces.txt so results.txt is never touched
        Settings reuseSettings;
        reuseSettings.walk.dirichletFaces = dirichletFaces;
        if (dirichletFaces) std::printf("EXPERIMENT: tet faces are Dirichlet (H = 0), not reflecting\n");
        if (validateReuse) {
            GeometryConfig cfg;
            Tet tet = referenceTet();
            printGeometry(tet, cfg);
            return runReuseValidation(cfg, reuseSettings, dirichletFaces ? "results_dirichlet_faces.txt" : "results.txt",
                                      quick) ? 0 : 1;
        }
        if (benchReuse) {
            GeometryConfig cfg;
            Tet tet = referenceTet();
            printGeometry(tet, cfg);
            return runReuseBench(cfg, reuseSettings);
        }
        if (compareOptions) {
            GeometryConfig cfg;
            Tet tet = referenceTet();
            printGeometry(tet, cfg);
            // --threads T adds _threadsT, so runs at several thread counts can sit side by side
            const std::string out = std::string("results_options") + (dirichletFaces ? "_dirichlet_faces" : "") +
                                    (quick ? "_quick" : "") + (threads > 0 ? "_threads" + std::to_string(threads) : "") +
                                    ".txt";
            return runCompareOptions(cfg, reuseSettings, out, quick, threads, repeats) ? 0 : 1;
        }
        if (!shotFile.empty())
            return runScreenshot(shotFile, std::clamp(tShot, 0.0, 1.0), std::clamp(tStartShot, 0.0, 1.0), shotSteps,
                                 shotRewalked, shotUi);

#ifndef NDEBUG
        // unoptimized Eigen/Zombie: ~60x slower, and the window stays blank until the first solve ends
        std::printf("WARNING: Debug build. One solve takes about a minute and the window stays blank meanwhile.\n"
                    "         Switch Visual Studio to the Release configuration (about 1 s per solve).\n");
#endif
        App app;
        printGeometry(app.tet, app.cfg);
        initPolyscope();
        app.registerTet();
        app.rebuild("start");
        app.setCamera();
        polyscope::state::userCallback = [&app] { app.ui(); };
        polyscope::show();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
