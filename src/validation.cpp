// The four validation checks of the spec. Everything printed also goes to results.txt.
#include "validation.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <random>
#include <set>
#include <thread>

namespace {

// printf to stdout and into a string (written to results.txt at the end)
class Report {
public:
    void operator()(const char* fmt, ...)
    {
        char buf[2048];
        va_list ap;
        va_start(ap, fmt);
        std::vsnprintf(buf, sizeof buf, fmt, ap);
        va_end(ap);
        std::fputs(buf, stdout);
        text += buf;
    }
    std::string text;
};

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

} // namespace

bool runValidation(const GeometryConfig& cfg, const Settings& s, const std::string& resultsPath)
{
    auto tStart = std::chrono::steady_clock::now();
    Tet tet(cfg.tetVerts);
    const double L = tet.L, eps = s.epsRel*L;
    const double tDefault = Settings{}.t;   // checks 1 and 3 use the default t, whatever the GUI shows
    const GType gTypes[2] = {GType::One, GType::SqrtFront};
    const char* gShort[2] = {"(a) g=1", "(b) g=min(1,sqrt(r/L))"};
    Report out;
    SolveStats total;   // over every solve of this validation run

    // ---------------------------------------------------------------- header
    std::time_t now = std::time(nullptr);
    std::tm tmNow{};
#ifdef _WIN32
    localtime_s(&tmNow, &now);
#else
    localtime_r(&now, &tmNow);
#endif
    char date[64];
    std::strftime(date, sizeof date, "%Y-%m-%d %H:%M:%S", &tmNow);
    Cut cutP = buildCut(tet, cfg, tDefault);
    out("Harmonic partial-cut enrichment on one tet, Walk on Stars (Zombie): validation\n");
    out("date: %s\n", date);
    out("tet: v0 = %s, v1 = %s, v2 = %s, v3 = %s; longest edge L = %.6f\n", fmt3(tet.v[0]).c_str(),
        fmt3(tet.v[1]).c_str(), fmt3(tet.v[2]).c_str(), fmt3(tet.v[3]).c_str(), L);
    out("cut plane n.X = d: n = %s, d = %.6f; front direction u = %s; front line direction n x u = %s\n",
        fmt3(cutP.n).c_str(), cutP.d, fmt3(cutP.u).c_str(), fmt3(cutP.w).c_str());
    out("plane-tet section: %zu vertices (%s); extent along u: s in [%.4f, %.4f]\n", cutP.fullPolygon.size(),
        cutP.fullIsQuad ? "quad" : "NOT a quad", cutP.sMin, cutP.sMax);
    out("epsilon = %.5g (%.4g L): lattice exclusion radius and jump-check offset\n", eps, s.epsRel);
    const double offsetMargin = faceOffsetMargin(tet);
    out("walk termination shell = %.3g (%.0e L) on distances to the cut shrunk by the face-offset margin %.3g\n"
        "(%.2g L), so walks stop within %.3g (%.3g L) of the cut; min star radius = %.3g (%.0e L);\n"
        "max walk length = %d\n",
        s.walk.absorbShellRel*L, s.walk.absorbShellRel, offsetMargin, offsetMargin/L,
        s.walk.absorbShellRel*L + offsetMargin, s.walk.absorbShellRel + offsetMargin/L, s.walk.reflectShellRel*L,
        s.walk.reflectShellRel, s.walk.maxWalkLength);
    out("lattice k = %d, N = %d walks/point for checks 2 and 4; hardware threads: %u\n", s.k, s.nWalks,
        std::thread::hardware_concurrency());
    out("RNG: base seed %llu; solve with tag T uses seed mix(base, T); point i of a solve uses pcg32 stream i\n",
        (unsigned long long)s.seed);
    if (cutP.empty()) {
        out("ERROR: the cut is empty at t = %.2f (the plane misses the tet?); nothing to validate.\n", tDefault);
        std::ofstream(resultsPath) << out.text;
        return false;
    }

    // ---------------------------------------------------------------- 1. convergence
    {
        // Probes in the cut frame around the middle F of the front segment, at the default t.
        // They are near or ahead of the front, where walks end on both sides of the cut, so
        // the per-walk values really vary (behind the front, close to the cut, almost every
        // walk ends on the same side and the estimate is nearly deterministic).
        Vec3 F = cutP.frontMidpoint();
        struct Probe { const char* where; Vec3 x; };
        std::vector<Probe> probes = {
            {"ahead of the front, in the cut plane (F + 0.05L u)", F + 0.05*L*cutP.u},
            {"above the front edge (F + 0.05L n)", F + 0.05*L*cutP.n},
            {"above the cut, behind the front (F - 0.05L u + 0.03L n)", F - 0.05*L*cutP.u + 0.03*L*cutP.n},
            {"ahead of and below the front (F + 0.10L u - 0.05L n)", F + 0.10*L*cutP.u - 0.05*L*cutP.n}};
        out("\n=== 1. Convergence at 4 probe points (t = %.2f, front at s = %.4f) ===\n", tDefault, cutP.s);
        out("front segment midpoint F = %s\n", fmt3(F).c_str());
        for (size_t p = 0; p < probes.size(); p++) {
            int moved = 0;   // a custom geometry could put a probe outside: pull it in, and say so
            while (tet.distToFaces(probes[p].x) < 0.01*L && moved < 100) {
                probes[p].x += 0.1*(tet.centroid - probes[p].x);
                moved++;
            }
            out("P%zu = %s  %s; dist to cut %.4f, to faces %.4f%s\n", p + 1, fmt3(probes[p].x).c_str(),
                probes[p].where, cutP.distToCut(probes[p].x), tet.distToFaces(probes[p].x),
                moved ? " (moved towards the centroid to stay inside)" : "");
        }
        out("8 independent repeats per N (own RNG streams); std = sample std of the 8 repeat means;\n"
            "pred = sqrt(mean per-walk variance / N) from the same runs; distinct = distinct repeat means.\n");

        const std::vector<int> Ns = {16, 64, 256, 1024, 4096};
        const int R = 8;
        for (int gi = 0; gi < 2; gi++) {
            CutSolver solver(tet, cutP, gTypes[gi], s.walk);
            // H[p][ni][r], per-walk variance[p][ni]
            std::vector<std::vector<std::vector<double>>> H(probes.size(), std::vector<std::vector<double>>(Ns.size()));
            std::vector<std::vector<double>> walkVar(probes.size(), std::vector<double>(Ns.size(), 0.0));
            for (size_t ni = 0; ni < Ns.size(); ni++) {
                std::vector<Vec3> pts;   // all 4 probes x 8 repeats in one parallel solve
                for (const Probe& pr : probes)
                    for (int r = 0; r < R; r++) pts.push_back(pr.x);
                std::vector<Estimate> est;
                total.add(solver.solve(pts, nullptr, Ns[ni], mixSeed(s.seed, 100 + 10*gi + ni), est));
                for (size_t p = 0; p < probes.size(); p++)
                    for (int r = 0; r < R; r++) {
                        H[p][ni].push_back(est[p*R + r].H);
                        walkVar[p][ni] += est[p*R + r].walkVar/R;
                    }
            }
            out("\n-- %s (seed tags %d..%d) --\n", gShort[gi], 100 + 10*gi, 100 + 10*gi + int(Ns.size()) - 1);
            for (size_t p = 0; p < probes.size(); p++) {
                out("P%zu      N      mean        std   std*sqrt(N)   pred std   distinct\n", p + 1);
                std::vector<double> nVals, stds;
                for (size_t ni = 0; ni < Ns.size(); ni++) {
                    double m = meanOf(H[p][ni]), sd = sampleStd(H[p][ni]);
                    std::set<double> distinct(H[p][ni].begin(), H[p][ni].end());
                    out("    %6d   %+.5f   %.5f     %.4f      %.5f     %zu/%d\n", Ns[ni], m, sd,
                        sd*std::sqrt(double(Ns[ni])), std::sqrt(walkVar[p][ni]/Ns[ni]), distinct.size(), R);
                    nVals.push_back(Ns[ni]);
                    stds.push_back(sd);
                }
                double slope = logLogSlope(nVals, stds);
                if (std::isnan(slope)) out("    log-log slope of std vs N: undefined (some std = 0)\n");
                else out("    log-log slope of std vs N: %+.3f   (1/sqrt(N) gives -0.5)\n", slope);
            }
        }
    }

    // ---------------------------------------------------------------- 2. complete cut
    {
        Cut cut1 = buildCut(tet, cfg, 1.0);
        Lattice lat = buildLattice(tet, cut1, s.k, eps);
        double rMax = 0;
        for (const Vec3& v : cut1.polygon) rMax = std::max(rMax, cut1.distToFrontLine(v));   // max over a convex polygon is at a vertex
        double gbMax = std::min(1.0, std::sqrt(rMax/L));
        out("\n=== 2. Complete-cut limit (t = 1): max |H - sign((X-p).n)| over the lattice ===\n");
        out("lattice: %d points, %d within epsilon of the cut excluded, %zu used; N = %d\n", lat.nTotal,
            lat.nExcluded, lat.pts.size(), s.nWalks);
        out("data on the cut at t = 1: (a) g = 1; (b) g ranges over [0, %.4f]: the front line touches the\n"
            "section at s_max and the farthest cut point is r_max = %.4f = %.3f L from it\n", gbMax, rMax, rMax/L);
        for (int gi = 0; gi < 2; gi++) {
            CutSolver solver(tet, cut1, gTypes[gi], s.walk);
            std::vector<Estimate> est;
            SolveStats st = solver.solve(lat.pts, nullptr, s.nWalks, mixSeed(s.seed, 200 + gi), est);
            total.add(st);
            out("-- %s (seed tag %d) --\n", gShort[gi], 200 + gi);
            double wrongWalks = 0;
            for (int sgn : {+1, -1}) {
                int n = 0, over = 0;
                double maxDev = 0, sumDev = 0, maxAbsH = 0;
                Vec3 worst = Vec3::Zero();
                for (size_t i = 0; i < lat.pts.size(); i++) {
                    if ((cut1.side(lat.pts[i]) > 0 ? 1 : -1) != sgn) continue;
                    double dev = std::fabs(est[i].H - sgn);
                    n++;
                    sumDev += dev;
                    wrongWalks += 0.5*dev*est[i].nCounted;
                    if (dev > 0.1) over++;
                    maxAbsH = std::max(maxAbsH, double(std::fabs(est[i].H)));
                    if (dev > maxDev) { maxDev = dev; worst = lat.pts[i]; }
                }
                out("   %c side: %3d points, max |H %s 1| = %.5f at %s, mean = %.5f, %d points > 0.1, max |H| = %.5f\n",
                    sgn > 0 ? '+' : '-', n, sgn > 0 ? "-" : "+", maxDev, fmt3(worst).c_str(),
                    n ? sumDev/n : 0.0, over, maxAbsH);
            }
            // with g = 1 every counted walk returns exactly +1 or -1 (or 0 if capped), so
            // n|H - sign|/2 summed over the points is the number of walks that ended on the
            // other side of the cut than they started
            if (gTypes[gi] == GType::One && st.walksCapped == 0)
                out("   walks that ended on the other side of the cut: %.0f of %lld\n", wrongWalks, st.walksCounted);
            out("   walks: %s\n", st.summary().c_str());
        }
        out("note: for (b), g_b = min(1, sqrt(r/L)) < 1 wherever r < L. Here r_max = %.3f L < L, so g_b < 1 on the\n"
            "whole cut (max %.4f), and by the maximum principle |H| <= %.4f everywhere: H cannot approach\n"
            "+-1, and every point has |H - sign| >= %.4f. g_b is smallest near where the front line touches\n"
            "the polygon (its vertex at s_max).\n", rMax/L, gbMax, gbMax, 1.0 - gbMax);
    }

    // ---------------------------------------------------------------- 3. jump
    {
        const int nPairs = 20, nWalksJump = 4096;
        std::vector<double> areas;
        for (const auto& t : cutP.tris)
            areas.push_back((cutP.polygon[t[1]] - cutP.polygon[t[0]]).cross(cutP.polygon[t[2]] - cutP.polygon[t[0]]).norm());
        std::discrete_distribution<int> pickTri(areas.begin(), areas.end());
        std::uniform_real_distribution<double> U(0.0, 1.0);
        auto sampleCut = [&](std::mt19937_64& g) {   // uniform on the cut polygon
            const auto& t = cutP.tris[pickTri(g)];
            double r1 = std::sqrt(U(g)), r2 = U(g);
            return Vec3((1 - r1)*cutP.polygon[t[0]] + r1*(1 - r2)*cutP.polygon[t[1]] + r1*r2*cutP.polygon[t[2]]);
        };
        // Keep X away from the cut's edges, where H is not close to +-g: at least
        // max(0.05 L, 2 eps) from them, but at most 0.8x the largest edge distance the cut has
        // (estimated from 4000 samples with their own seed), or a large epsilon leaves no X.
        std::mt19937_64 rngDepth(mixSeed(s.seed, 301));
        double depth = 0;
        for (int i = 0; i < 4000; i++) depth = std::max(depth, cutP.distToPolygonEdges(sampleCut(rngDepth)));
        const double wantMargin = std::max(0.05*L, 2.0*eps);
        const double margin = std::min(wantMargin, 0.8*depth);
        std::mt19937_64 rng(mixSeed(s.seed, 300));
        std::vector<Vec3> X;
        for (int tries = 0; int(X.size()) < nPairs && tries < 200000; tries++) {
            Vec3 x = sampleCut(rng);
            if (cutP.distToPolygonEdges(x) < margin) continue;
            if (tet.distToFaces(x + eps*cutP.n) <= 0 || tet.distToFaces(x - eps*cutP.n) <= 0) continue;
            X.push_back(x);
        }
        out("\n=== 3. Jump across the cut (t = %.2f): %zu pairs X +- epsilon n, epsilon = %.5f, N = %d ===\n",
            tDefault, X.size(), eps, nWalksJump);
        out("X uniform on the cut, at least %.4f from its edges (sampling seed tag 300)\n", margin);
        if (margin < wantMargin)
            out("NOTE: max(0.05 L, 2 eps) = %.4f is more than this cut allows (its deepest point is %.4f from the\n"
                "edges), so the margin is capped at 0.8x that depth and all X lie near the cut's deepest point\n",
                wantMargin, depth);
        if (X.empty()) {
            out("CHECK 3 NOT RUN at this epsilon: no point X on the cut is %.4f or more from its edges with X +- epsilon n\n"
                "inside the tet, so there is no jump to report; use a smaller epsilon\n", margin);
        } else if (int(X.size()) < nPairs) {
            out("WARNING: only %zu admissible points found; the means below are over these %zu pairs\n", X.size(),
                X.size());
        }
        for (int gi = 0; gi < 2 && !X.empty(); gi++) {
            CutSolver solver(tet, cutP, gTypes[gi], s.walk);
            std::vector<Vec3> pts;
            for (const Vec3& x : X) { pts.push_back(x + eps*cutP.n); pts.push_back(x - eps*cutP.n); }
            std::vector<Estimate> est;
            SolveStats st = solver.solve(pts, nullptr, nWalksJump, mixSeed(s.seed, 310 + gi), est);
            total.add(st);
            out("-- %s (seed tag %d) --\n", gShort[gi], 310 + gi);
            out("pair  X                              r/L     H+ (s.e.)          H- (s.e.)           jump     2g       jump/2g\n");
            double sumJump = 0, sum2g = 0, sumRatio = 0;
            for (size_t i = 0; i < X.size(); i++) {
                const Estimate &ep = est[2*i], &em = est[2*i + 1];
                double jump = ep.H - em.H, twoG = 2.0*cutP.g(X[i], gTypes[gi]);
                sumJump += jump;
                sum2g += twoG;
                sumRatio += jump/twoG;
                out("%3zu   %s   %.3f   %+.4f (%.4f)   %+.4f (%.4f)   %.4f   %.4f   %.4f\n", i + 1,
                    fmt3(X[i]).c_str(), cutP.distToFrontLine(X[i])/L, ep.H, ep.Hse, em.H, em.Hse, jump, twoG, jump/twoG);
            }
            double n = double(X.size());   // > 0 here
            out("   mean jump = %.4f, mean 2g = %.4f, mean jump/2g = %.4f\n", sumJump/n, sum2g/n, sumRatio/n);
            out("   walks: %s\n", st.summary().c_str());
        }
    }

    // ---------------------------------------------------------------- 4. sweep
    {
        const double ts[] = {0.2, 0.4, 0.6, 0.8, 1.0};
        out("\n=== 4. Progressive sweep: fraction of lattice points with |H| > 0.9 (k = %d, N = %d) ===\n", s.k, s.nWalks);
        out("   t     s(t)      cut area   points used (excluded)   %-22s   %-22s   max g_b on cut\n", gShort[0],
            gShort[1]);
        long long escaped = 0, capped = 0;
        for (int ti = 0; ti < 5; ti++) {
            Cut cut = buildCut(tet, cfg, ts[ti]);
            Lattice lat = buildLattice(tet, cut, s.k, eps);
            double gbMax = 0;   // |H| <= max g on the cut (maximum principle)
            for (const Vec3& v : cut.polygon) gbMax = std::max(gbMax, cut.g(v, GType::SqrtFront));
            char cols[2][64];
            for (int gi = 0; gi < 2; gi++) {
                CutSolver solver(tet, cut, gTypes[gi], s.walk);
                std::vector<Estimate> est;
                SolveStats st = solver.solve(lat.pts, nullptr, s.nWalks, mixSeed(s.seed, 400 + 10*ti + gi), est);
                total.add(st);
                escaped += st.walksEscaped;
                capped += st.walksCapped;
                int over = 0;
                for (const Estimate& e : est) over += std::fabs(e.H) > 0.9f;
                std::snprintf(cols[gi], sizeof cols[gi], "%.4f (%d/%zu)", double(over)/double(std::max<size_t>(1, est.size())),
                              over, est.size());
            }
            out("  %.1f   %+.4f   %.5f    %4zu (%3d)               %-22s   %-22s   %.4f\n", ts[ti], cut.s, cut.area(),
                lat.pts.size(), lat.nExcluded, cols[0], cols[1], gbMax);
        }
        out("   (seed tags 400 + 10*i_t + i_g; %lld walks escaped, %lld hit the cap over these solves)\n", escaped, capped);
        out("   for (b), |H| <= max g_b on the cut (maximum principle), so |H| > 0.9 is impossible while max g_b < 0.9\n");
    }

    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - tStart).count();
    out("\nall validation solves: %s\n", total.summary().c_str());
    out("validation wall time: %.2f s\n", sec);

    std::ofstream f(resultsPath);
    f << out.text;
    f.close();
    if (!f) {
        std::fprintf(stderr, "could not write %s\n", resultsPath.c_str());
        return false;
    }
    std::printf("written: %s\n", resultsPath.c_str());
    return true;
}
