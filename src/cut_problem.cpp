#include "cut_problem.h"

#include <zombie/zombie.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace {

double secondsSince(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// distance from p to the segment [a, b]
double pointSegmentDist(const Vec3& p, const Vec3& a, const Vec3& b)
{
    Vec3 ab = b - a;
    double len2 = ab.squaredNorm();
    double s = len2 > 0 ? std::clamp((p - a).dot(ab)/len2, 0.0, 1.0) : 0.0;
    return (p - (a + s*ab)).norm();
}

} // namespace

const char* gTypeName(GType g)
{
    return g == GType::One ? "(a) g = 1" : "(b) g = min(1, sqrt(r/L))";
}

uint64_t mixSeed(uint64_t a, uint64_t b)
{
    uint64_t z = a + 0x9E3779B97F4A7C15ull*(b + 1);
    z = (z ^ (z >> 30))*0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27))*0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

// =====================================================================================
// Tet
// =====================================================================================
Tet::Tet(const std::array<Vec3, 4>& verts) : v(verts)
{
    centroid = (v[0] + v[1] + v[2] + v[3])/4.0;
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++) L = std::max(L, (v[i] - v[j]).norm());

    Eigen::Matrix3d E;
    E.col(0) = v[1] - v[0];
    E.col(1) = v[2] - v[0];
    E.col(2) = v[3] - v[0];
    if (std::fabs(E.determinant()) < 1e-12*L*L*L) throw std::runtime_error("degenerate tetrahedron");
    Einv = E.inverse();

    // face opposite vertex o, wound so that (b-a)x(c-a) points away from o (outward)
    for (int o = 0; o < 4; o++) {
        std::array<int, 3> f{};
        int m = 0;
        for (int i = 0; i < 4; i++)
            if (i != o) f[m++] = i;
        Vec3 nrm = (v[f[1]] - v[f[0]]).cross(v[f[2]] - v[f[0]]);
        if (nrm.dot(v[o] - v[f[0]]) > 0) std::swap(f[1], f[2]);
        faces[o] = f;
    }
}

Eigen::Vector4d Tet::bary(const Vec3& x) const
{
    Vec3 b = Einv*(x - v[0]);
    return Eigen::Vector4d(1.0 - b.sum(), b[0], b[1], b[2]);
}

bool Tet::contains(const Vec3& x, double tol) const
{
    return bary(x).minCoeff() >= -tol;
}

double Tet::distToFaces(const Vec3& x) const
{
    double dmin = std::numeric_limits<double>::max();
    for (const auto& f : faces) {
        Vec3 nrm = (v[f[1]] - v[f[0]]).cross(v[f[2]] - v[f[0]]).normalized();   // outward
        dmin = std::min(dmin, (v[f[0]] - x).dot(nrm));
    }
    return dmin;
}

// =====================================================================================
// Plane sections and the cut
// =====================================================================================
std::vector<Vec3> planeTetSection(const Tet& tet, const Vec3& n, double d, const Vec3& tangent,
                                  int* nVertsOnPlane)
{
    const double tol = 1e-12*tet.L;
    double f[4];
    for (int i = 0; i < 4; i++) f[i] = n.dot(tet.v[i]) - d;

    std::vector<Vec3> pts;
    int onPlane = 0;
    for (int i = 0; i < 4; i++)
        if (std::fabs(f[i]) <= tol) { pts.push_back(tet.v[i]); onPlane++; }
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++)
            if ((f[i] > tol && f[j] < -tol) || (f[i] < -tol && f[j] > tol))
                pts.push_back(tet.v[i] + (f[i]/(f[i] - f[j]))*(tet.v[j] - tet.v[i]));
    if (nVertsOnPlane) *nVertsOnPlane = onPlane;
    if (pts.size() < 3) return {};

    // order counter-clockwise about n: angle in the basis (a, b = n×a), a×b = n
    Vec3 c = Vec3::Zero();
    for (const Vec3& p : pts) c += p;
    c /= double(pts.size());
    Vec3 a = tangent - tangent.dot(n)*n;
    if (a.norm() < 1e-9) a = n.unitOrthogonal();
    a.normalize();
    Vec3 b = n.cross(a);
    std::sort(pts.begin(), pts.end(), [&](const Vec3& p, const Vec3& q) {
        return std::atan2((p - c).dot(b), (p - c).dot(a)) < std::atan2((q - c).dot(b), (q - c).dot(a));
    });
    return pts;
}

Cut buildCut(const Tet& tet, const GeometryConfig& cfg, double t)
{
    Cut c;
    c.n = cfg.planeNormal.normalized();
    c.d = cfg.planeOffset;
    c.L = tet.L;
    c.t = std::clamp(t, 0.0, 1.0);
    Vec3 u = cfg.frontDir - cfg.frontDir.dot(c.n)*c.n;   // project the front direction into the plane
    if (u.norm() < 1e-9) throw std::runtime_error("front direction is parallel to the plane normal");
    c.u = u.normalized();
    c.w = c.n.cross(c.u);

    int onPlane = 0;
    c.fullPolygon = planeTetSection(tet, c.n, c.d, c.u, &onPlane);
    c.fullIsQuad = c.fullPolygon.size() == 4 && onPlane == 0;
    if (c.fullPolygon.empty()) return c;   // the plane misses the tet

    c.sMin = std::numeric_limits<double>::max();
    c.sMax = -c.sMin;
    for (const Vec3& p : c.fullPolygon) {
        c.sMin = std::min(c.sMin, p.dot(c.u));
        c.sMax = std::max(c.sMax, p.dot(c.u));
    }
    c.s = c.t >= 1.0 ? c.sMax : c.sMin + c.t*(c.sMax - c.sMin);
    if (c.t <= 0.0) return c;   // t = 0: no cut at all

    // clip the convex section against the half-plane X·u <= s (Sutherland-Hodgman, one plane)
    std::vector<Vec3> clipped;
    const std::vector<Vec3>& P = c.fullPolygon;
    for (size_t i = 0; i < P.size(); i++) {
        const Vec3& a = P[i];
        const Vec3& b = P[(i + 1)%P.size()];
        double fa = a.dot(c.u) - c.s, fb = b.dot(c.u) - c.s;
        if (fa <= 0) clipped.push_back(a);
        if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0)) clipped.push_back(a + (fa/(fa - fb))*(b - a));
    }
    // drop (near-)duplicate consecutive vertices produced by clipping through a vertex
    for (size_t i = 0; i < clipped.size() && clipped.size() > 1;) {
        if ((clipped[i] - clipped[(i + 1)%clipped.size()]).norm() <= 1e-9*c.L) clipped.erase(clipped.begin() + i);
        else i++;
    }
    c.polygon = clipped;

    // fan triangulation; zero-area triangles are skipped (their normal would be NaN and
    // break the side test). The polygon is CCW about n, so each triangle already has
    // ((b-a)x(c-a))·n > 0; the swap is a safeguard. This winding is what makes Zombie's
    // "normalAligned" flag mean "the walk ended on the n·X > d side".
    for (int i = 1; i + 1 < int(c.polygon.size()); i++) {
        std::array<int, 3> tri{0, i, i + 1};
        double twiceArea = (c.polygon[tri[1]] - c.polygon[tri[0]]).cross(c.polygon[tri[2]] - c.polygon[tri[0]]).dot(c.n);
        if (std::fabs(twiceArea) <= 1e-12*c.L*c.L) continue;
        if (twiceArea < 0) std::swap(tri[1], tri[2]);
        c.tris.push_back(tri);
    }
    if (c.tris.empty()) c.polygon.clear();
    return c;
}

double Cut::distToFrontLine(const Vec3& x) const
{
    Vec3 r = x - frontLinePoint();
    r -= r.dot(w)*w;
    return r.norm();
}

double Cut::g(const Vec3& x, GType type) const
{
    if (type == GType::One) return 1.0;
    return std::min(1.0, std::sqrt(distToFrontLine(x)/L));
}

bool Cut::polygonContains(const Vec3& x) const
{
    if (polygon.size() < 3) return false;
    for (size_t i = 0; i < polygon.size(); i++) {
        const Vec3& a = polygon[i];
        const Vec3& b = polygon[(i + 1)%polygon.size()];
        if ((b - a).cross(x - a).dot(n) < 0) return false;
    }
    return true;
}

double Cut::distToPolygonEdges(const Vec3& x) const
{
    Vec3 xp = x - side(x)*n;
    double dmin = std::numeric_limits<double>::max();
    for (size_t i = 0; i < polygon.size(); i++)
        dmin = std::min(dmin, pointSegmentDist(xp, polygon[i], polygon[(i + 1)%polygon.size()]));
    return dmin;
}

double Cut::distToCut(const Vec3& x) const
{
    if (empty()) return std::numeric_limits<double>::infinity();
    double h = side(x);
    if (polygonContains(x)) return std::fabs(h);
    double e = distToPolygonEdges(x);   // the closest point is on the polygon boundary
    return std::sqrt(h*h + e*e);
}

Vec3 Cut::frontMidpoint() const
{
    const std::vector<Vec3>& P = polygon.empty() ? fullPolygon : polygon;
    double wMin = std::numeric_limits<double>::max(), wMax = -wMin;
    Vec3 c = Vec3::Zero();
    for (const Vec3& p : P) {
        c += p;
        if (std::fabs(p.dot(u) - s) <= 1e-9*L) {
            wMin = std::min(wMin, p.dot(w));
            wMax = std::max(wMax, p.dot(w));
        }
    }
    if (wMin > wMax) return P.empty() ? Vec3::Zero() : Vec3(c/double(P.size()));
    return frontLinePoint() + 0.5*(wMin + wMax)*w;   // frontLinePoint()·w = 0
}

double Cut::area() const
{
    double a = 0;
    for (const auto& t : tris)
        a += 0.5*(polygon[t[1]] - polygon[t[0]]).cross(polygon[t[2]] - polygon[t[0]]).norm();
    return a;
}

// =====================================================================================
// Solver
// =====================================================================================
void SolveStats::add(const SolveStats& o)
{
    nPoints += o.nPoints;
    walksRequested += o.walksRequested;
    walksCounted += o.walksCounted;
    walksEscaped += o.walksEscaped;
    walksCapped += o.walksCapped;
    totalSteps += o.totalSteps;
    seconds += o.seconds;
    noCut = noCut || o.noCut;
}

std::string SolveStats::summary() const
{
    char buf[512];
    if (noCut) {
        std::snprintf(buf, sizeof buf, "%d points, no cut (t = 0): no walks run, H := 0", nPoints);
    } else {
        double escPct = walksRequested ? 100.0*double(walksEscaped)/double(walksRequested) : 0.0;
        std::snprintf(buf, sizeof buf,
                      "%d points, %lld walks counted of %lld, %lld escaped (%.3f%%, discarded), "
                      "%lld hit the walk-length cap, mean walk length %.1f, %.3f s",
                      nPoints, walksCounted, walksRequested, walksEscaped, escPct, walksCapped,
                      meanWalkLength(), seconds);
    }
    return buf;
}

// All Zombie objects for one cut. Zombie's query lambdas capture pointers/references to
// the handlers, to 'queries' and to the PDE, so this struct lives on the heap and is never
// copied or moved after setup; a new one is built for every cut.
struct CutSolver::Impl {
    zombie::FcpwDirichletBoundaryHandler<3> absorbing;   // the cut
    zombie::FcpwNeumannBoundaryHandler<3> reflecting;    // the four tet faces
    zombie::GeometricQueries<3> queries;
    zombie::PDE<float, 3> pde;
    zombie::WalkSettings settings{0.0f, 0.0f, 1, true};
};

double faceOffsetMargin(const Tet& tet)
{
    // offsetPointAlongDirectionImpl<3> (zombie/core/geometric_queries.h) moves coordinate i
    // by 2^-16 n_i if |x_i| < 1/32 and by int(256 n_i) float ulps otherwise. One ulp of x_i
    // is at most |x_i| 2^-23, so for |n| = 1 the push is at most max(2^-16, 2^-15 max|x_i|)
    // long, and every point of the (convex) tet has |x_i| <= the largest vertex coordinate.
    double M = 0;
    for (const Vec3& p : tet.v) M = std::max(M, p.cwiseAbs().maxCoeff());
    return 2.0*std::max(std::ldexp(1.0, -16), std::ldexp(M, -15));
}

CutSolver::CutSolver(const Tet& tet, const Cut& cut, GType g, const WalkParams& params)
{
    const double absorbShell = params.absorbShellRel*tet.L;
    const double margin = faceOffsetMargin(tet);
    stopDistance_ = absorbShell + margin;
    // t = 0: there is no Dirichlet boundary, and Laplace with pure zero-Neumann data only
    // determines H up to a constant. Nothing is built; solve() returns H = 0.
    if (cut.empty()) return;

    impl_ = std::make_unique<Impl>();
    Impl& im = *impl_;
    using zombie::Vector3;
    using zombie::Vector3i;

    std::vector<Vector3> tetV, cutV;
    std::vector<Vector3i> tetF, cutF;
    for (const Vec3& p : tet.v) tetV.push_back(p.cast<float>());
    for (const auto& f : tet.faces) tetF.emplace_back(f[0], f[1], f[2]);
    for (const Vec3& p : cut.polygon) cutV.push_back(p.cast<float>());
    for (const auto& f : cut.tris) cutF.emplace_back(f[0], f[1], f[2]);

    // bounding box, used by Zombie's default escape test (replaced below)
    Vector3 bMin = tetV[0], bMax = tetV[0];
    for (const Vector3& p : tetV) { bMin = bMin.cwiseMin(p); bMax = bMax.cwiseMax(p); }
    Vector3 pad = Vector3::Constant(1e-3f*float(tet.L));
    im.queries.domainIsWatertight = false;
    im.queries.domainMin = bMin - pad;
    im.queries.domainMax = bMax + pad;

    const float silhouettePrecision = 1e-3f;
    im.absorbing.buildAccelerationStructure(cutV, cutF, true, false, false);
    zombie::populateGeometricQueriesForDirichletBoundary<3>(im.absorbing, im.queries);
    im.reflecting.buildAccelerationStructure(
        tetV, tetF, zombie::getIgnoreCandidateSilhouetteCallback(true, silhouettePrecision), true, false, false);
    zombie::populateGeometricQueriesForNeumannBoundary<3>(im.reflecting, zombie::getBranchTraversalWeightCallback(),
                                                          im.queries);

    // Conservative distance to the cut (see faceOffsetMargin): every ball, including one
    // centred at a pushed-off face point, then stays on its own side of the cut. Zombie uses
    // the same distance for its termination test, so walks stop within shell + margin.
    // Measured at t = 1, g = 1, N = 4096: 10 of 6.89M walks ended on the wrong side without
    // it, 0 of 13.8M with it (same mean walk length).
    auto exactDist = im.queries.computeDistToAbsorbingBoundary;
    im.queries.computeDistToAbsorbingBoundary = [exactDist, m = float(margin)](const Vector3& x, bool signedDist) -> float {
        float dist = exactDist(x, signedDist);
        float shrunk = std::max(0.0f, std::fabs(dist) - m);
        return dist < 0.0f ? -shrunk : shrunk;
    };

    // Escape test: Zombie's default is "outside the bounding box", which misses walks
    // that leak out through the slanted face and lets them keep walking outside the tet
    // (stage 1: H = 0.91 instead of 1 at the apex). Use exact tet containment instead.
    // Walks that fail it are discarded by Zombie and counted as "escaped".
    Eigen::Matrix3f Einv = tet.Einv.cast<float>();
    Vector3 v0 = tetV[0];
    im.queries.outsideBoundingDomain = [Einv, v0](const Vector3& x) -> bool {
        Vector3 b = Einv*(x - v0);
        const float tol = 1e-4f;
        return b.minCoeff() < -tol || b.sum() > 1.0f + tol;
    };

    // Laplace: no source, no absorption; zero Neumann data on the faces
    im.pde.absorptionCoeff = 0.0f;
    im.pde.isSourceConstant = true;
    im.pde.areRobinConditionsPureNeumann = true;
    im.pde.areRobinCoeffsNonnegative = true;
    im.pde.source = [](const Vector3&) -> float { return 0.0f; };
    im.pde.robin = [](const Vector3&, const Vector3&, bool) -> float { return 0.0f; };
    im.pde.robinCoeff = [](const Vector3&, const Vector3&, bool) -> float { return 0.0f; };
    im.pde.hasReflectingBoundaryConditions = [](const Vector3&) -> bool { return true; };

    // Dirichlet data ±g. Zombie calls this with x = the walk's end point PROJECTED onto
    // the cut (so x lies exactly in the plane and sign((x-p)·n) would be meaningless) and
    // normalAligned = (signed distance of the unprojected end point w.r.t. the nearest cut
    // triangle's normal) > 0. Every cut triangle is wound with ((b-a)x(c-a))·n > 0, so
    // normalAligned is true exactly when the walk ended on the n·X > d side: +g there,
    // -g on the other side. g itself is evaluated at the boundary point x. Called
    // concurrently from TBB threads; the captured cut copy is only read.
    im.pde.dirichlet = [cut, g](const Vector3& x, bool normalAligned) -> float {
        float gx = float(cut.g(x.cast<double>(), g));
        return normalAligned ? gx : -gx;
    };

    // solveDoubleSided must be true: the cut is an open two-sided surface, and with
    // false Zombie always reports normalAligned = false (H = -g everywhere, stage 1).
    im.settings = zombie::WalkSettings(float(absorbShell), float(params.reflectShellRel*tet.L),
                                       params.maxWalkLength, /*solveDoubleSided=*/true);
    im.settings.silhouettePrecision = silhouettePrecision;
    // source and Neumann data are zero: skip their per-step sampling (same estimate)
    im.settings.ignoreSourceContribution = true;
    im.settings.ignoreReflectingBoundaryContribution = true;
}

CutSolver::~CutSolver() = default;

SolveStats CutSolver::solve(const std::vector<Vec3>& pts, const std::vector<char>* wantGrad, int nWalks,
                            uint64_t seed, std::vector<Estimate>& out) const
{
    auto t0 = std::chrono::steady_clock::now();
    SolveStats st;
    st.nPoints = int(pts.size());
    out.assign(pts.size(), Estimate{});
    if (!impl_) {   // t = 0
        st.noCut = true;
        st.seconds = secondsSince(t0);
        return st;
    }
    const Impl& im = *impl_;

    std::vector<zombie::SamplePoint<float, 3>> samples;
    samples.reserve(pts.size());
    const uint64_t state = mixSeed(seed, 0);
    for (size_t i = 0; i < pts.size(); i++) {
        zombie::Vector3 x = pts[i].cast<float>();
        float dA = im.queries.computeDistToAbsorbingBoundary(x, false);
        float dR = im.queries.computeDistToReflectingBoundary(x, false);
        bool grad = wantGrad && (*wantGrad)[i];
        samples.emplace_back(x, zombie::Vector3::Zero(), zombie::SampleType::InDomain,
                             grad ? zombie::EstimationQuantity::SolutionAndGradient
                                  : zombie::EstimationQuantity::Solution,
                             1.0f, dA, dR);
        // Replace Zombie's clock-based seed: point i gets its own pcg32 stream i, and the
        // state depends on the seed, so runs are reproducible and every point/solve is
        // independent. All walks of one point run on one TBB thread with this generator.
        samples.back().rng = pcg32(state, uint64_t(i));
        if (grad) st.walksRequested += 2*std::max(1, nWalks/2);   // antithetic pairs
        else st.walksRequested += dA <= im.settings.epsilonShellForAbsorbingBoundary ? 1 : nWalks;
    }

    // count walks that hit the cap; they still count, with terminal value 0
    std::atomic<long long> capped{0};
    zombie::WalkOnStars<float, 3> wost(im.queries, {},
        [&capped](zombie::WalkCompletionCode code, const zombie::WalkState<float, 3>&) -> float {
            if (code == zombie::WalkCompletionCode::ExceededMaxWalkLength) capped++;
            return 0.0f;
        });

    std::vector<zombie::SampleStatistics<float, 3>> stats(pts.size());
    std::vector<int> walks(pts.size(), nWalks);
    wost.solve(im.pde, im.settings, walks, samples, stats, /*runSingleThreaded=*/false);   // TBB over points

    for (size_t i = 0; i < pts.size(); i++) {
        const auto& s = stats[i];
        Estimate& e = out[i];
        e.nCounted = s.getSolutionEstimateCount();
        e.H = s.getEstimatedSolution();
        e.walkVar = e.nCounted > 1 ? s.getEstimatedSolutionVariance() : 0.0f;
        e.Hse = std::sqrt(e.walkVar/float(std::max(1, e.nCounted)));
        if (samples[i].estimationQuantity == zombie::EstimationQuantity::SolutionAndGradient) {
            const float* gr = s.getEstimatedGradient();
            std::vector<float> gv = s.getEstimatedGradientVariance();
            float ng = float(std::max(1, s.getGradientEstimateCount()));
            e.grad = Eigen::Vector3f(gr[0], gr[1], gr[2]);
            e.gradSe = Eigen::Vector3f(std::sqrt(gv[0]/ng), std::sqrt(gv[1]/ng), std::sqrt(gv[2]/ng));
            e.hasGrad = true;
        }
        st.walksCounted += e.nCounted;
        st.totalSteps += std::llround(double(s.getMeanWalkLength())*e.nCounted);
    }
    st.walksEscaped = st.walksRequested - st.walksCounted;
    st.walksCapped = capped.load();
    st.seconds = secondsSince(t0);
    return st;
}

// =====================================================================================
// Query points
// =====================================================================================
Lattice buildLattice(const Tet& tet, const Cut& cut, int k, double eps)
{
    Lattice lat;
    lat.k = std::max(1, k);
    const int K = lat.k;
    for (int i0 = 0; i0 <= K; i0++)
        for (int i1 = 0; i0 + i1 <= K; i1++)
            for (int i2 = 0; i0 + i1 + i2 <= K; i2++) {
                int i3 = K - i0 - i1 - i2;
                Vec3 x = (double(i0)*tet.v[0] + double(i1)*tet.v[1] + double(i2)*tet.v[2] + double(i3)*tet.v[3])/double(K);
                x += kBoundaryPull*(tet.centroid - x);   // all points: a uniform 0.1% shrink
                lat.nTotal++;
                if (cut.distToCut(x) < eps) { lat.nExcluded++; continue; }
                lat.pts.push_back(x);
                lat.index.push_back({i0, i1, i2, i3});
                lat.interior.push_back(i0 > 0 && i1 > 0 && i2 > 0 && i3 > 0);
            }
    return lat;
}

Slice buildSlice(const Tet& tet, const Cut& cut, const SliceParams& sp)
{
    const double deg = 3.14159265358979323846/180.0;
    double yaw = sp.yawDeg*deg, pitch = sp.pitchDeg*deg;
    Slice sl;
    sl.normal = (std::cos(pitch)*(std::cos(yaw)*cut.w + std::sin(yaw)*cut.u) + std::sin(pitch)*cut.n).normalized();
    sl.origin = tet.centroid + sp.offset*sl.normal;
    const Vec3& m = sl.normal;

    // in-plane axes (a, b), a×b = m; at the default orientation a = -u, b = n
    Vec3 a = cut.n.cross(m);
    if (a.norm() < 1e-6) a = cut.u.cross(m);
    a.normalize();
    Vec3 b = m.cross(a);

    std::vector<Vec3> sec = planeTetSection(tet, m, m.dot(sl.origin), a);
    if (sec.size() < 3) return sl;   // the plane misses the tet

    // section polygon in 2D (CCW in (a, b) because planeTetSection orders it CCW about m)
    std::vector<Eigen::Vector2d> P;
    Eigen::Vector2d lo(1e300, 1e300), hi(-1e300, -1e300);
    for (const Vec3& p : sec) {
        Eigen::Vector2d q((p - sl.origin).dot(a), (p - sl.origin).dot(b));
        P.push_back(q);
        lo = lo.cwiseMin(q);
        hi = hi.cwiseMax(q);
    }
    auto inside = [&](const Eigen::Vector2d& q) {
        for (size_t i = 0; i < P.size(); i++) {
            Eigen::Vector2d e = P[(i + 1)%P.size()] - P[i], r = q - P[i];
            if (e.x()*r.y() - e.y()*r.x() < 0) return false;
        }
        return true;
    };
    auto closestOnBoundary = [&](const Eigen::Vector2d& q) {
        Eigen::Vector2d best = P[0];
        double bestD = 1e300;
        for (size_t i = 0; i < P.size(); i++) {
            Eigen::Vector2d p0 = P[i], e = P[(i + 1)%P.size()] - p0;
            double s = std::clamp((q - p0).dot(e)/std::max(e.squaredNorm(), 1e-300), 0.0, 1.0);
            Eigen::Vector2d c = p0 + s*e;
            if ((q - c).norm() < bestD) { bestD = (q - c).norm(); best = c; }
        }
        return best;
    };

    // square grid over the section's bounding box; cells with at least one corner inside
    // are kept, and their outside corners are snapped onto the section boundary, so the
    // slice follows the tet faces instead of having a staircase edge
    double h = std::max(hi.x() - lo.x(), hi.y() - lo.y())/double(std::max(4, sp.res));
    int nx = std::max(1, int(std::ceil((hi.x() - lo.x())/h))), ny = std::max(1, int(std::ceil((hi.y() - lo.y())/h)));
    auto gridPt = [&](int i, int j) { return Eigen::Vector2d(lo.x() + i*h, lo.y() + j*h); };
    std::vector<char> in((nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++) in[j*(nx + 1) + i] = inside(gridPt(i, j));

    std::vector<int> id((nx + 1)*(ny + 1), -1);
    auto vertexId = [&](int i, int j) {
        int g = j*(nx + 1) + i;
        if (id[g] < 0) {
            Eigen::Vector2d q = in[g] ? gridPt(i, j) : closestOnBoundary(gridPt(i, j));
            Vec3 x = sl.origin + q.x()*a + q.y()*b;
            x += kBoundaryPull*(tet.centroid - x);
            id[g] = int(sl.verts.size());
            sl.verts.push_back(x);
        }
        return id[g];
    };
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            int g00 = j*(nx + 1) + i;
            if (!(in[g00] || in[g00 + 1] || in[g00 + nx + 1] || in[g00 + nx + 2])) continue;
            int v00 = vertexId(i, j), v10 = vertexId(i + 1, j), v11 = vertexId(i + 1, j + 1), v01 = vertexId(i, j + 1);
            sl.tris.push_back({v00, v10, v11});
            sl.tris.push_back({v00, v11, v01});
        }
    return sl;
}

// =====================================================================================
// Full scene
// =====================================================================================
SceneResult solveScene(const Tet& tet, const GeometryConfig& cfg, const Settings& s, uint64_t seed)
{
    auto t0 = std::chrono::steady_clock::now();
    SceneResult r;
    r.seed = seed;
    r.cut = buildCut(tet, cfg, s.t);
    r.lattice = buildLattice(tet, r.cut, s.k, s.epsRel*tet.L);
    r.slice = buildSlice(tet, r.cut, s.slice);
    CutSolver solver(tet, r.cut, s.g, s.walk);
    r.buildSeconds = secondsSince(t0);
    // ∇H only at interior lattice points: Zombie's gradient estimator starts with a ball of
    // radius min(dist to cut, dist to faces), and its noise grows like 1/radius, so on the
    // faces/edges/vertices (radius ~1e-3) it is useless (stage 1: s.e. ~600 at N = 256).
    r.latticeStats = solver.solve(r.lattice.pts, &r.lattice.interior, s.nWalks, mixSeed(seed, 1), r.latticeEst);
    r.sliceStats = solver.solve(r.slice.verts, nullptr, s.nWalks, mixSeed(seed, 2), r.sliceEst);
    r.totalSeconds = secondsSince(t0);
    return r;
}
