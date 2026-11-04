// --check-xi: the solve lives in local coordinates ξ on the reference tet; a rest tet only maps
// the cut into ξ. See validation_xi.h for the three checks.
#include "validation_xi.h"

#include "validation.h"
#include "walk_cache.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

const char* const kResultsSectionXi = "local coordinates xi (--check-xi)";

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
    std::string text;
};

bool sameBits(const void* a, const void* b, size_t n) { return std::memcmp(a, b, n) == 0; }
bool sameVec(const Vec3& a, const Vec3& b) { return sameBits(a.data(), b.data(), 3*sizeof(double)); }

std::string fmt3(const Vec3& v)
{
    char buf[96];
    std::snprintf(buf, sizeof buf, "(%.4f, %.4f, %.4f)", v.x(), v.y(), v.z());
    return buf;
}

Eigen::Matrix3d jacobian(const std::array<Vec3, 4>& X)
{
    Eigen::Matrix3d J;
    J.col(0) = X[1] - X[0];
    J.col(1) = X[2] - X[0];
    J.col(2) = X[3] - X[0];
    return J;
}

// every geometric field of the two cuts, bit for bit
bool sameCut(const Cut& a, const Cut& b)
{
    if (!sameVec(a.n, b.n) || !sameVec(a.u, b.u) || !sameVec(a.w, b.w)) return false;
    const std::pair<double, double> scalars[] = {{a.d, b.d}, {a.L, b.L}, {a.t, b.t}, {a.sMin, b.sMin}, {a.sMax, b.sMax}, {a.s, b.s}};
    for (const auto& [x, y] : scalars)
        if (!sameBits(&x, &y, sizeof(double))) return false;
    if (a.polygon.size() != b.polygon.size() || a.fullPolygon.size() != b.fullPolygon.size() || a.tris != b.tris)
        return false;
    for (size_t i = 0; i < a.polygon.size(); i++)
        if (!sameVec(a.polygon[i], b.polygon[i])) return false;
    for (size_t i = 0; i < a.fullPolygon.size(); i++)
        if (!sameVec(a.fullPolygon[i], b.fullPolygon[i])) return false;
    return true;
}

bool sameEstimate(const Estimate& a, const Estimate& b)
{
    return sameBits(&a.H, &b.H, sizeof(float)) && sameBits(&a.Hse, &b.Hse, sizeof(float)) &&
           sameBits(a.grad.data(), b.grad.data(), 3*sizeof(float)) &&
           sameBits(a.gradSe.data(), b.gradSe.data(), 3*sizeof(float)) && a.hasGrad == b.hasGrad &&
           a.nCounted == b.nCounted && sameBits(&a.walkVar, &b.walkVar, sizeof(float));
}

// points of the two estimate lists that differ in any bit (-1 if the lists differ in length)
long long countDiffering(const std::vector<Estimate>& a, const std::vector<Estimate>& b, double* maxDH = nullptr)
{
    if (a.size() != b.size()) return -1;
    long long bad = 0;
    double m = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if (!sameEstimate(a[i], b[i])) bad++;
        m = std::max(m, double(std::fabs(a[i].H - b[i].H)));
    }
    if (maxDH) *maxDH = m;
    return bad;
}

} // namespace

bool runXiCheck(const Settings& s, const std::string& resultsPath)
{
    Report out;
    bool ok = true;
    const Tet ref = referenceTet();
    const GeometryConfig def;

    out("H(xi) in local coordinates: the solve runs on the reference tet, a rest tet only maps the cut into xi\n");
    out("reference tet: %s %s %s %s, L = %.6f (the L of g type (b); epsilon, lattice, slice and grad H are in xi too)\n",
        fmt3(ref.v[0]).c_str(), fmt3(ref.v[1]).c_str(), fmt3(ref.v[2]).c_str(), fmt3(ref.v[3]).c_str(), ref.L);
    out("map X = X0 + J xi, J = [X1-X0, X2-X0, X3-X0]: n_xi = J^T n_X/|J^T n_X|, d_xi = (d_X - n_X.X0)/|J^T n_X|,\n"
        "front direction f_xi = J^T f_X (then projected into the xi plane); t is kept\n");

    // ---------------------------------------------------------------------------------------
    // X1: random rest tets and rest-space cuts
    // ---------------------------------------------------------------------------------------
    {
        const int nTets = 1000;
        const double ts[] = {0.05, 0.25, 0.5, 0.75, 1.0};
        std::mt19937_64 rng(mixSeed(s.seed, 4242));
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        std::normal_distribution<double> Nrm(0.0, 1.0);
        long long nCuts = 0, badCount = 0, badVert = 0, badSide = 0, nSide = 0, nInverted = 0;
        double maxErr = 0;
        for (int it = 0; it < nTets; it++) {
            GeometryConfig cfg;
            Eigen::Matrix3d J;
            do {
                for (int r = 0; r < 3; r++)
                    for (int c = 0; c < 3; c++) J(r, c) = U(rng);
            } while (std::fabs(J.determinant()) < 0.05*std::pow(J.colwise().norm().maxCoeff(), 3));
            const Vec3 X0(3*U(rng), 3*U(rng), 3*U(rng));
            cfg.restVerts = {X0, X0 + J.col(0), X0 + J.col(1), X0 + J.col(2)};
            if (J.determinant() < 0) nInverted++;
            const Tet rest(cfg.restVerts);
            // a plane through a point well inside, any orientation, an unnormalized normal
            double lam[4], sum = 0;
            for (double& l : lam) sum += (l = 0.15 + 0.85*(0.5 + 0.5*U(rng)));
            Vec3 P = Vec3::Zero();
            for (int i = 0; i < 4; i++) P += (lam[i]/sum)*cfg.restVerts[i];
            Vec3 nX(Nrm(rng), Nrm(rng), Nrm(rng));
            nX.normalize();
            cfg.planeNormal = (0.3 + 2.7*(0.5 + 0.5*U(rng)))*nX;
            cfg.planeOffset = cfg.planeNormal.normalized().dot(P);
            do {
                cfg.frontDir = Vec3(Nrm(rng), Nrm(rng), Nrm(rng));
            } while ((cfg.frontDir - cfg.frontDir.dot(nX)*nX).norm() < 0.3*cfg.frontDir.norm());

            for (double t : ts) {
                const Cut cx = buildCut(ref, cfg, t);
                const Cut cr = buildCutOnTet(rest, cfg.planeNormal.normalized(), cfg.planeOffset, cfg.frontDir, t);
                nCuts++;
                if (cx.polygon.size() != cr.polygon.size() || cx.fullIsQuad != cr.fullIsQuad) {
                    badCount++;
                    continue;
                }
                // every xi vertex, mapped to rest space, is a vertex of the rest-space cut (orientation may flip)
                for (const Vec3& q : cx.polygon) {
                    const Vec3 X = X0 + J*q;
                    double best = std::numeric_limits<double>::infinity();
                    for (const Vec3& Y : cr.polygon) best = std::min(best, (X - Y).norm());
                    maxErr = std::max(maxErr, best/rest.L);
                    if (best > 1e-10*rest.L) badVert++;
                }
                // the "+g" side: sign(n_xi.xi - d_xi) = sign(n_X.X - d_X)
                for (int k = 0; k < 4; k++) {
                    double l[4], sl = 0;
                    for (double& v : l) sl += (v = 0.5 + 0.5*U(rng));
                    const Vec3 xi(l[1]/sl, l[2]/sl, l[3]/sl);
                    const double hx = cx.side(xi), hr = cr.side(X0 + J*xi);
                    if (std::fabs(hr) < 1e-9*rest.L) continue;
                    nSide++;
                    if ((hx > 0) != (hr > 0)) badSide++;
                }
            }
        }
        const bool pass = badCount == 0 && badVert == 0 && badSide == 0;
        ok = ok && pass;
        out("\nX1 pull-back geometry: %d random rest tets (%lld with det J < 0), random planes through an inner point,\n"
            "random front directions; t = 0.05, 0.25, 0.5, 0.75, 1 (%lld cuts). The xi cut mapped by X0 + J xi against the\n"
            "cut built in rest space: vertex count or quad/triangle differs: %lld; vertices off by > 1e-10 L_rest: %lld\n"
            "(max %.2e L_rest); '+g' side differs at %lld of %lld random points -> %s\n",
            nTets, nInverted, nCuts, badCount, badVert, maxErr, badSide, nSide, pass ? "PASS" : "FAIL");
    }

    // ---------------------------------------------------------------------------------------
    // X2: the default xi cut pushed into a skewed rest tet B and pulled back
    // ---------------------------------------------------------------------------------------
    GeometryConfig cfgB;
    cfgB.restVerts = {Vec3(0.3, -0.2, 0.1), Vec3(2.1, 0.4, -0.3), Vec3(0.5, 1.7, 0.2), Vec3(-0.1, 0.6, 1.3)};
    {
        const XiCut x0 = pullBackCut(def);
        const Eigen::Matrix3d J = jacobian(cfgB.restVerts);
        const Eigen::Matrix3d JinvT = J.inverse().transpose();
        const Vec3 X0 = cfgB.restVerts[0];
        // n.xi = d with xi = J^-1 (X - X0): (J^-T n).X = d + (J^-T n).X0
        const Vec3 p = JinvT*x0.n;
        cfgB.planeNormal = p;
        cfgB.planeOffset = (x0.d + p.dot(X0))/p.norm();
        cfgB.frontDir = JinvT*x0.f;
        const Cut c0 = buildCut(ref, def, 0.6), cB = buildCut(ref, cfgB, 0.6);
        const double dn = (c0.n - cB.n).norm(), dd = std::fabs(c0.d - cB.d), du = (c0.u - cB.u).norm(),
                     dt = std::fabs(c0.s - cB.s);
        const bool pass = std::max({dn, dd, du, dt}) <= 1e-13;
        ok = ok && pass;
        out("\nX2 round trip: a skewed rest tet B = %s %s %s %s\n"
            "(det J = %.3f) carrying the default cut, pushed forward from xi. Pulled back: |dn| = %.1e, |dd| = %.1e,\n"
            "|du| = %.1e, front position at t = 0.6 |ds| = %.1e (all <= 1e-13) -> %s\n",
            fmt3(cfgB.restVerts[0]).c_str(), fmt3(cfgB.restVerts[1]).c_str(), fmt3(cfgB.restVerts[2]).c_str(),
            fmt3(cfgB.restVerts[3]).c_str(), J.determinant(), dn, dd, du, dt, pass ? "PASS" : "FAIL");
    }

    // ---------------------------------------------------------------------------------------
    // X3: rest tet C = 2 R B (R = 180 degree rotation about y) carrying the same cut moved with it.
    // Sign flips and factors of 2 are exact in floating point, so the pull-back gives bitwise B's xi cut.
    // ---------------------------------------------------------------------------------------
    {
        const Vec3 R(-1, 1, -1);
        GeometryConfig cfgC;
        for (int i = 0; i < 4; i++) cfgC.restVerts[i] = 2.0*R.cwiseProduct(cfgB.restVerts[i]);
        cfgC.planeNormal = R.cwiseProduct(cfgB.planeNormal);
        cfgC.planeOffset = 2.0*cfgB.planeOffset;
        cfgC.frontDir = R.cwiseProduct(cfgB.frontDir);

        long long cutsBad = 0;
        for (double t : {0.1, 0.3, 0.6, 0.9, 1.0})
            if (!sameCut(buildCut(ref, cfgB, t), buildCut(ref, cfgC, t))) cutsBad++;
        out("\nX3 the solve sees only the xi cut: rest tet C = 2 R B, R = rotation by 180 degrees about y, with B's cut moved\n"
            "along (C is not B: other position, orientation and size). xi cuts of B and C at t = 0.1, 0.3, 0.6, 0.9, 1\n"
            "bitwise different: %lld\n", cutsBad);
        ok = ok && cutsBad == 0;

        // Zombie (first gate): solveScene on both, all lattice and slice estimates bit for bit
        Settings sz = s;
        sz.t = 0.6;
        sz.nWalks = std::min(s.nWalks, 64);
        sz.k = 8;
        sz.slice.res = 16;
        for (GType g : {GType::One, GType::SqrtFront}) {
            sz.g = g;
            const SceneResult rB = solveScene(ref, cfgB, sz, sz.seed), rC = solveScene(ref, cfgC, sz, sz.seed);
            const SceneResult rD = solveScene(ref, def, sz, sz.seed);
            const long long dl = countDiffering(rB.latticeEst, rC.latticeEst),
                            ds = countDiffering(rB.sliceEst, rC.sliceEst);
            double mD = 0;
            const long long dDef = countDiffering(rB.latticeEst, rD.latticeEst, &mD);
            const bool pass = dl == 0 && ds == 0;
            ok = ok && pass;
            out("   Zombie, %s, t = 0.6, N = %d: B vs C, points with any bit different: lattice %lld of %zu, slice %lld of\n"
                "   %zu -> %s. (Info: B vs the default rest tet, whose xi cut agrees with B's to round-off only: %lld of %zu\n"
                "   lattice points differ, max |dH| = %.2e)\n",
                gTypeName(g), sz.nWalks, dl, rB.latticeEst.size(), ds, rB.sliceEst.size(), pass ? "PASS" : "FAIL", dDef,
                rB.latticeEst.size(), mD);
        }

        // own loop (GUI and reuse): a full solve at t = 0.2, then 4 incremental steps to t = 0.6
        const int N = std::min(s.nWalks, 64);
        const Lattice lat = buildLattice(ref, Cut{}, 8, 0.0);
        QueryPoints q;
        q.add(lat.pts, &lat.interior, mixSeed(s.seed, 4343));
        for (GType g : {GType::One, GType::SqrtFront}) {
            StripCut cutB(ref, cfgB), cutC(ref, cfgC);
            WalkCache A, B;
            A.configure(ref, q, N, g, s.walk);
            B.configure(ref, q, N, g, s.walk);
            cutB.addStrip(cutB.sOf(0.2));
            cutC.addStrip(cutC.sOf(0.2));
            A.solveFull(cutB);
            B.solveFull(cutC);
            long long bad = 0;
            auto compare = [&] {
                for (size_t i = 0; i < A.sums.size(); i++)
                    if (!sameBits(&A.sums[i], &B.sums[i], sizeof(PointSums))) bad++;
            };
            compare();
            for (int st = 1; st <= 4; st++) {
                const double t = 0.2 + 0.1*st;
                cutB.addStrip(cutB.sOf(t));
                cutC.addStrip(cutC.sOf(t));
                A.advance(cutB, UpdateMode::Incremental);
                B.advance(cutC, UpdateMode::Incremental);
                compare();
            }
            const bool pass = bad == 0;
            ok = ok && pass;
            out("   own loop, %s, N = %d, %zu lattice points, full solve at t = 0.2 + 4 incremental steps to 0.6:\n"
                "   B vs C, (point, step) pairs whose walk sums differ in any bit: %lld of %zu -> %s\n",
                gTypeName(g), N, A.sums.size(), bad, 5*A.sums.size(), pass ? "PASS" : "FAIL");
        }
    }

    out("\nxi checks X1, X2, X3: %s\n", ok ? "all PASS" : "SOME FAILED");
    if (!writeResultsSection(resultsPath, kResultsSectionXi, out.text)) {
        std::fprintf(stderr, "could not write %s\n", resultsPath.c_str());
        return false;
    }
    std::printf("written: %s (section \"%s\")\n", resultsPath.c_str(), kResultsSectionXi);
    return ok;
}
