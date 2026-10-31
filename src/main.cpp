// Harmonic partial-cut enrichment H(X) on a single tetrahedron, by Walk on Stars (Zombie),
// shown in Polyscope.
//
//   partial_cut                                   GUI
//   partial_cut --validate                        the four validation checks, headless; writes ./results.txt
//   partial_cut --bench                           one full default solve (lattice + gradient + slice), timed
//   partial_cut --screenshot out.png [--t 0.6]    build the scene, solve, save a screenshot, exit
//
// Geometry (tet vertices, cut plane, front direction): GeometryConfig in cut_problem.h.

#include "cut_problem.h"
#include "validation.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/screenshot.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include "polyscope/volume_mesh.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string>
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
    std::printf("tet: (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f), longest edge L = %.4f\n",
                tet.v[0].x(), tet.v[0].y(), tet.v[0].z(), tet.v[1].x(), tet.v[1].y(), tet.v[1].z(),
                tet.v[2].x(), tet.v[2].y(), tet.v[2].z(), tet.v[3].x(), tet.v[3].y(), tet.v[3].z(), tet.L);
    std::printf("cut plane n.X = d: n = (%.4f, %.4f, %.4f), d = %.4f; front direction u = (%.4f, %.4f, %.4f)\n",
                full.n.x(), full.n.y(), full.n.z(), full.d, full.u.x(), full.u.y(), full.u.z());
    if (full.fullIsQuad)
        std::printf("plane-tet section is a quad (check passed)\n");
    else
        std::printf("WARNING: the plane-tet section is not a quad (%zu vertices)\n", full.fullPolygon.size());
}

// ---------------------------------------------------------------------------------------
// GUI / screenshot scene
// ---------------------------------------------------------------------------------------
struct App {
    GeometryConfig cfg;
    Tet tet{cfg.tetVerts};
    Settings s;

    // UI copies of the settings (ImGui wants float/int)
    float t = 0.6f, epsRel = 0.01f, sliceOffset = 0.0f, sliceYaw = 0.0f, slicePitch = 0.0f;
    int nWalks = 256, k = 12, gType = 0;

    SceneResult scene;
    uint64_t solveCount = 0;
    bool firstRegistration = true, gradCloudConfigured = false;
    std::string validationMsg;

    App()
    {
        t = float(s.t);
        epsRel = float(s.epsRel);
        nWalks = s.nWalks;
        k = s.k;
        gType = int(s.g);
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

    void solveAndShow()
    {
        s.t = t;
        s.epsRel = epsRel;
        s.nWalks = nWalks;
        s.k = k;
        s.g = GType(gType);
        s.slice.offset = sliceOffset;
        s.slice.yawDeg = sliceYaw;
        s.slice.pitchDeg = slicePitch;

        // every solve is from scratch, with a fresh seed
        scene = solveScene(tet, cfg, s, mixSeed(s.seed, solveCount++));
        std::printf("solve #%llu: t = %.3f, N = %d, epsilon = %.4g L, %s, k = %d: total %.3f s (build %.3f s)\n",
                    (unsigned long long)solveCount, s.t, s.nWalks, s.epsRel, gTypeName(s.g), s.k,
                    scene.totalSeconds, scene.buildSeconds);
        std::printf("  lattice: %s\n  slice:   %s\n", scene.latticeStats.summary().c_str(),
                    scene.sliceStats.summary().c_str());
        std::printf("  H at vertices: %s\n", vertexHText().c_str());

        const Cut& cut = scene.cut;
        if (cut.empty()) {
            polyscope::removeStructure("cut", false);
        } else {
            polyscope::SurfaceMesh* cm = polyscope::registerSurfaceMesh("cut", toArrays(cut.polygon), cut.tris);
            cm->setSurfaceColor({0.95f, 0.70f, 0.15f});
            cm->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
            cm->setEdgeWidth(1.0);
        }

        // lattice: H and its standard error at every point
        const Lattice& lat = scene.lattice;
        std::vector<double> H, Hse;
        for (const Estimate& e : scene.latticeEst) {
            H.push_back(e.H);
            Hse.push_back(e.Hse);
        }
        polyscope::PointCloud* pc = polyscope::registerPointCloud("lattice", toArrays(lat.pts));
        auto* hq = pc->addScalarQuantity("H", H);
        hq->setColorMap("coolwarm");
        hq->setMapRange({-1.0, 1.0});
        pc->addScalarQuantity("H std. error", Hse);

        // ∇H exists only at the interior lattice points (see solveScene), so it gets its own
        // point cloud holding exactly those points: a point without an estimate shows no
        // arrow, rather than a zero arrow with a zero standard error.
        std::vector<Vec3> gPts;
        std::vector<std::array<double, 3>> grad;
        std::vector<double> gradSe, gnorm;
        for (size_t i = 0; i < scene.latticeEst.size(); i++) {
            const Estimate& e = scene.latticeEst[i];
            if (!e.hasGrad) continue;
            gPts.push_back(lat.pts[i]);
            grad.push_back({e.grad.x(), e.grad.y(), e.grad.z()});
            gradSe.push_back(e.gradSe.norm());
            gnorm.push_back(e.grad.norm());
        }
        if (gPts.empty()) {   // k <= 3: no interior lattice points
            polyscope::removeStructure("lattice grad H", false);
        } else {
            polyscope::PointCloud* gc = polyscope::registerPointCloud("lattice grad H", toArrays(gPts));
            auto* gq = gc->addVectorQuantity("grad H", grad);
            gc->addScalarQuantity("|std. error of grad H|", gradSe);
            // Display scaling only: Polyscope scales arrows so the LONGEST one has a fixed
            // length, and a few noisy estimates next to the cut would shrink all others to
            // nothing. Scale to the 90th percentile of |grad H| instead (longer arrows stay
            // longer), but never below 1/L, so that a field that is ~0 everywhere (e.g. t = 1,
            // g = 1) is not blown up to full-length arrows.
            std::nth_element(gnorm.begin(), gnorm.begin() + gnorm.size()*9/10, gnorm.end());
            gq->setVectorLengthRange(std::max(gnorm[gnorm.size()*9/10], 1.0/tet.L));
            if (!gradCloudConfigured) {   // later, Polyscope's persistent options keep the user's choices
                gq->setEnabled(true);
                gq->setVectorLengthScale(0.6*tet.L/std::max(1, k), /*isRelative=*/false);   // ~0.6 lattice spacing
                gq->setVectorRadius(0.0025);
                gq->setVectorColor({0.1f, 0.1f, 0.1f});
                gc->setPointRadius(0.004);   // smaller than the lattice points, so H stays visible
                gc->setPointColor({0.1f, 0.1f, 0.1f});
                gc->setEnabled(false);       // hidden at startup; tick it in the structure list
                gradCloudConfigured = true;
            }
        }

        // slice: H
        const Slice& sl = scene.slice;
        if (sl.tris.empty()) {   // the slice plane misses the tet
            polyscope::removeStructure("slice", false);
        } else {
            std::vector<double> sH;
            for (const Estimate& e : scene.sliceEst) sH.push_back(e.H);
            polyscope::SurfaceMesh* sm = polyscope::registerSurfaceMesh("slice", toArrays(sl.verts), sl.tris);
            sm->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
            auto* sq = sm->addVertexScalarQuantity("H", sH);
            sq->setColorMap("coolwarm");
            sq->setMapRange({-1.0, 1.0});
            if (firstRegistration) sq->setEnabled(true);
        }
        if (firstRegistration) {   // later, Polyscope's persistent options keep the user's choices
            hq->setEnabled(true);
            pc->setPointRadius(0.006);
            firstRegistration = false;
        }
    }

    // H +- 1 s.e. at the four tet vertices, i.e. the lattice points with barycentric index k
    // (pulled 0.1% inside like every query point); "excluded" if within epsilon of the cut
    std::string vertexHText() const
    {
        std::string out;
        char buf[64];
        for (int v = 0; v < 4; v++) {
            std::snprintf(buf, sizeof(buf), "v%d excluded  ", v);
            for (size_t i = 0; i < scene.lattice.index.size(); i++)
                if (scene.lattice.index[i][v] == scene.lattice.k)
                    std::snprintf(buf, sizeof(buf), "v%d %+.3f+-%.3f  ", v, scene.latticeEst[i].H,
                                  scene.latticeEst[i].Hse);
            out += buf;
        }
        return out;
    }

    void setCamera()
    {
        // looks at the default slice plane (normal n x u) at an angle, from slightly ahead of
        // the front and from the -n side, so the cut shows as a strip ending at the front
        Cut c = buildCut(tet, cfg, 1.0);
        Vec3 dir = (-c.w - 0.35*c.n + 0.45*c.u).normalized();
        Vec3 eye = tet.centroid + 1.8*tet.L*dir;
        polyscope::view::lookAt(glm::vec3(eye.x(), eye.y(), eye.z()),
                                glm::vec3(tet.centroid.x(), tet.centroid.y(), tet.centroid.z()));
    }

    void ui()
    {
        bool changed = false;
        const ImGuiSliderFlags clamp = ImGuiSliderFlags_AlwaysClamp;
        ImGui::PushItemWidth(170);

        ImGui::TextUnformatted("Cut and solver (re-solves when a control is released)");
        ImGui::SliderFloat("t (front position)", &t, 0.0f, 1.0f, "%.3f", clamp);
        changed |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SliderInt("N walks per point", &nWalks, 16, 4096, "%d", clamp | ImGuiSliderFlags_Logarithmic);
        changed |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SliderFloat("epsilon / L", &epsRel, float(10*s.walk.absorbShellRel), 0.1f, "%.4f",
                           clamp | ImGuiSliderFlags_Logarithmic);
        changed |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::Text("  epsilon = %.4g (walks stop within %.2g of the cut)", epsRel*tet.L,
                    s.walk.absorbShellRel*tet.L + faceOffsetMargin(tet));
        ImGui::TextUnformatted("Dirichlet data g on the cut (H = +g / -g):");
        changed |= ImGui::RadioButton("(a) g = 1", &gType, 0);
        ImGui::SameLine();
        changed |= ImGui::RadioButton("(b) g = min(1, sqrt(r/L))", &gType, 1);
        ImGui::SliderInt("k (lattice subdivisions)", &k, 1, 30, "%d", clamp);
        changed |= ImGui::IsItemDeactivatedAfterEdit();

        ImGui::Separator();
        ImGui::TextUnformatted("Slice plane (visualization only)");
        float R = 0.0f;
        for (const Vec3& v : tet.v) R = std::max(R, float((v - tet.centroid).norm()));
        ImGui::SliderFloat("offset from centroid", &sliceOffset, -R, R, "%.3f", clamp);
        changed |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SliderFloat("yaw (deg, about n)", &sliceYaw, -180.0f, 180.0f, "%.1f", clamp);
        changed |= ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SliderFloat("pitch (deg, towards n)", &slicePitch, -90.0f, 90.0f, "%.1f", clamp);
        changed |= ImGui::IsItemDeactivatedAfterEdit();

        ImGui::Separator();
        if (ImGui::Button("Re-solve (new seed)")) changed = true;
        ImGui::SameLine();
        if (ImGui::Button("Run validation")) {
            std::string path = (std::filesystem::current_path()/"results.txt").string();
            auto t0 = std::chrono::steady_clock::now();
            bool ok = runValidation(cfg, s, path);
            double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            validationMsg = ok ? "validation (" + std::to_string(int(sec + 0.5)) + " s) written to " + path
                               : "validation: could not write " + path;
        }
        if (!validationMsg.empty()) ImGui::TextWrapped("%s", validationMsg.c_str());

        ImGui::Separator();
        const Cut& c = scene.cut;
        if (c.empty()) {
            ImGui::TextWrapped("t = 0: no cut. Laplace with zero Neumann data on every face has no unique "
                               "solution, so no walks are run and H is set to 0.");
        } else {
            ImGui::Text("cut: %zu-gon, area %.4f, front at s = %.4f in [%.4f, %.4f]", c.polygon.size(), c.area(),
                        c.s, c.sMin, c.sMax);
        }
        ImGui::Text("plane-tet section: %s", c.fullIsQuad ? "quad (ok)" : "NOT a quad");
        ImGui::Text("lattice: %zu points (%d of %d within epsilon of the cut excluded)", scene.lattice.pts.size(),
                    scene.lattice.nExcluded, scene.lattice.nTotal);
        ImGui::Text("H at vertices: %s", vertexHText().c_str());
        ImGui::Text("slice: %zu points", scene.slice.verts.size());
        ImGui::Text("last solve: %.3f s (lattice %.3f s, slice %.3f s)", scene.totalSeconds, scene.latticeStats.seconds,
                    scene.sliceStats.seconds);
        SolveStats all = scene.latticeStats;
        all.add(scene.sliceStats);
        if (!all.noCut)
            ImGui::Text("walks: %lld counted, %lld escaped, %lld hit the cap; mean length %.1f", all.walksCounted,
                        all.walksEscaped, all.walksCapped, all.meanWalkLength());
        int nGrad = 0;
        for (const Estimate& e : scene.latticeEst) nGrad += e.hasGrad;
        if (nGrad > 0)
            ImGui::TextWrapped("grad H: estimated only at the %d interior lattice points (point cloud 'lattice grad H'); "
                               "points on the tet faces have no estimate.", nGrad);
        else if (!c.empty())   // (t = 0 is explained above)
            ImGui::TextWrapped("grad H: not estimated (k = %d has no interior lattice points).", scene.lattice.k);
        ImGui::PopItemWidth();

        if (changed) solveAndShow();
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
    Tet tet(cfg.tetVerts);
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

int runScreenshot(const std::string& file, double t, bool withUi)
{
    App app;
    app.t = float(t);
    printGeometry(app.tet, app.cfg);
    polyscope::options::usePrefsFile = false;   // fixed window size, nothing written to the working directory
    initPolyscope();
    polyscope::view::setWindowSize(1280, 960);
    app.registerTet();
    app.solveAndShow();
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

void printUsage()
{
    std::printf("usage:\n"
                "  partial_cut                                  Polyscope GUI\n"
                "  partial_cut --validate                       run the 4 validation checks, write ./results.txt\n"
                "  partial_cut --bench                          time one full default solve\n"
                "  partial_cut --screenshot file.png [--t T] [--with-ui]\n"
                "                                               solve (front at T, default 0.6), save a screenshot\n"
                "                                               (--with-ui: include the control panels)\n");
}

} // namespace

int main(int argc, char** argv)
{
    bool validate = false, bench = false, shotUi = false;
    std::string shotFile;
    double tShot = Settings{}.t;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--validate") == 0) validate = true;
        else if (std::strcmp(argv[i], "--bench") == 0) bench = true;
        else if (std::strcmp(argv[i], "--screenshot") == 0 && i + 1 < argc) shotFile = argv[++i];
        else if (std::strcmp(argv[i], "--t") == 0 && i + 1 < argc) tShot = std::atof(argv[++i]);
        else if (std::strcmp(argv[i], "--with-ui") == 0) shotUi = true;
        else if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) { printUsage(); return 0; }
        else { std::fprintf(stderr, "unknown argument: %s\n", argv[i]); printUsage(); return 2; }
    }

    try {
        if (validate) {
            GeometryConfig cfg;
            Tet tet(cfg.tetVerts);
            printGeometry(tet, cfg);
            return runValidation(cfg, Settings{}, "results.txt") ? 0 : 1;
        }
        if (bench) return runBench();
        if (!shotFile.empty()) return runScreenshot(shotFile, std::clamp(tShot, 0.0, 1.0), shotUi);

#ifndef NDEBUG
        // unoptimized Eigen/Zombie: ~60x slower, and the window stays blank until the first solve ends
        std::printf("WARNING: Debug build. One solve takes about a minute and the window stays blank meanwhile.\n"
                    "         Switch Visual Studio to the Release configuration (about 1 s per solve).\n");
#endif
        App app;
        printGeometry(app.tet, app.cfg);
        initPolyscope();
        app.registerTet();
        app.solveAndShow();
        app.setCamera();
        polyscope::state::userCallback = [&app] { app.ui(); };
        polyscope::show();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
