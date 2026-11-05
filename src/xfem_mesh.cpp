// The partial-cut enrichment on a whole tet mesh (xfem_mesh.h): buildXfemRing's method (coupling_gate.cpp, check 4's fast
// path) generalised from the 8-tet ring to any mesh and any set of cracks given as level sets. The per-tet pieces are the
// gate's own (gate_parts.h): tetCrackAt, ringPoints, xiFrontOf, cvGram, cvField, xfemTetQuad and the render helpers.
#include "xfem_mesh.h"

#include "gate_parts.h"
#include "walk_cache.h"

#include <Eigen/SparseCholesky>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <thread>

namespace {

using namespace gate;
using clk = std::chrono::steady_clock;

constexpr double kInf = std::numeric_limits<double>::infinity();

// Check 4's scales carried to a mesh by each tet's longest edge L (check 4's tets: sqrt 2)
constexpr double kEpsRel = 0.01;                                  // glue points within this x L of the tet's crack dropped
constexpr double kFrontBandRel[2] = {0.15/kSqrt2, 0.05/kSqrt2};   // the export quadrature's front bands (check 4: 0.15, 0.05)
constexpr double kCornerOffRel = 1e-4/kSqrt2;                     // self-check: points this far from a tet's corners
constexpr int kRenderPerEdge = 8;                                 // render triangles along the mesh's mean edge
constexpr uint64_t kSeed = 20261031;                              // the gate's seed (Settings)
constexpr int kGlueQ = 32, kGlueR = 16;                           // the glue's walks per point (check 4's fast glue)

double secondsSince(clk::time_point t0) { return std::chrono::duration<double>(clk::now() - t0).count(); }

// a unit vector orthogonal to a (deterministic)
Vec3 orthoOf(const Vec3& a)
{
    const Vec3 t = std::fabs(a.x()) < 0.6 ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
    return (t - t.dot(a)*a).normalized();
}

// A level set's corner values interpolated in tet T: the value at xi, and the physical gradient
double atXi(const Eigen::Vector4d& v, const Vec3& xi) { return (1.0 - xi.sum())*v[0] + xi[0]*v[1] + xi[1]*v[2] + xi[2]*v[3]; }
Vec3 gradOf(const TetFrame& T, const Eigen::Vector4d& v) { return T.JinvT*Vec3(v[1] - v[0], v[2] - v[0], v[3] - v[0]); }

// One crack in one tet: its linear level sets, the plane n.X = d (n unit, phi > 0 on n.X > d) and the front line u.X = s
// (u in the plane, towards growing psi) they give, the crack in the tet's xi frame (tetCrackAt), its class and, in a
// class-2 tet, the tet's piece of the front.
struct MeshTetCrack {
    Eigen::Vector4d phi = Eigen::Vector4d::Zero(), psi = Eigen::Vector4d::Zero();
    Vec3 n = Vec3::UnitX(), u = Vec3::UnitY();
    double d = 0, s = 0;
    RingTetCrack rt;
    int cls = 0;
    std::vector<Vec3> front;
    // J = atan2(phi_h, psi_h)/pi at xi (sign: phi's sign forced to it, 0 = phi's own)
    double J(const Vec3& xi, int sign = 0) const
    {
        const double f = atXi(phi, xi);
        return std::atan2(sign == 0 ? f : sign*std::fabs(f), atXi(psi, xi))/std::numbers::pi;
    }
};
MeshTetCrack meshTetCrack(const TetFrame& T, const Eigen::Vector4d& phi, const Eigen::Vector4d& psi, double L)
{
    MeshTetCrack m;
    m.phi = phi;
    m.psi = psi;
    const Vec3 gPhi = gradOf(T, phi), gPsi = gradOf(T, psi);
    if (!(gPhi.norm() > 0)) throw std::invalid_argument("phi is constant in a tet: no crack plane there");
    m.n = gPhi.normalized();
    m.d = m.n.dot(T.X[0]) - phi[0]/gPhi.norm();
    // on the plane X = d n + a u + b w: psi = psi_0 + gPsi.(X - X_0), linear in a alone when u is gPsi's in-plane part
    const Vec3 mp = gPsi - gPsi.dot(m.n)*m.n;
    if (mp.norm() > 1e-12*std::max(gPsi.norm(), 1e-300)) {
        m.u = mp.normalized();
        m.s = (gPsi.dot(T.X[0]) - psi[0] - gPsi.dot(m.n)*m.d)/mp.norm();
    } else {   // psi constant on the plane (a crack cut through: psi < 0 everywhere): the front line outside the tet
        m.u = orthoOf(m.n);
        const double onPlane = psi[0] + gPsi.dot(m.d*m.n - T.X[0]);
        double lo = kInf, hi = -kInf;
        for (const Vec3& X : T.X) {
            lo = std::min(lo, m.u.dot(X));
            hi = std::max(hi, m.u.dot(X));
        }
        m.s = onPlane < 0 ? hi + L : lo - L;
    }
    m.rt = tetCrackAt(T, m.n, m.d, m.u, m.s, m.u);
    m.cls = !m.rt.cut() ? 0 : m.rt.tRaw >= 1 ? 1 : 2;
    if (m.cls == 2) {   // the two ends of the front in the tet: the cut polygon's vertices on its front line, the farthest pair
        std::vector<Vec3> fv;
        for (const Vec3& p : frontVertices(m.rt.xc.cut)) fv.push_back(T.toX(p));
        double best = -1;
        for (size_t i = 0; i < fv.size(); i++)
            for (size_t j = i + 1; j < fv.size(); j++)
                if ((fv[i] - fv[j]).norm() > best) {
                    best = (fv[i] - fv[j]).norm();
                    m.front = {fv[i], fv[j]};
                }
    }
    return m;
}

// The fast path on a list of tets (runRingFast's per-tet part): per tet the Q run (cut tets: crack -lambda_k, faces 0, odd
// mirror ball) and the R run (where J - J_xi is not 0 on the faces: crack 0, faces J - J_xi, plain first ball) at the points
// ps, seeds mixSeed(kSeed, tagBase + 10 t) for Q and mixSeed(that, 77) for R; all Q jobs in one parallel loop, then all R.
struct MeshFast {
    std::vector<CvRun> cv;
    double secQ = 0, secR = 0;
    long long walksQ = 0, walksR = 0, capped = 0, escaped = 0;
    int nQtets = 0, nRtets = 0;
};
MeshFast runMeshFast(const std::vector<const TetFrame*>& T, const std::vector<const MeshTetCrack*>& mt,
                     const std::vector<PointSet>& ps, int nQ, int nR, int B, uint64_t tagBase)
{
    if (B < 3 || nQ%B || nR%B || nQ/B < 2 || nR/B < 2)
        throw std::invalid_argument("fast mesh path: at least 3 batches, and nQ, nR multiples of the batch count with at least "
                                    "2 walks per batch");
    const size_t N = T.size();
    MeshFast f;
    f.cv.resize(N);
    const Tet ref = referenceTet();
    for (size_t t = 0; t < N; t++) {   // fastDecide with the tet's own J: R walks iff J - J_xi is not 0 on a face
        CvRun& c = f.cv[t];
        const MeshTetCrack& m = *mt[t];
        c.xf = xiFrontOf(*T[t], m.n, m.d, m.u, m.s);
        c.cut = m.rt.cut();
        for (int o = 0; o < 4; o++) {
            std::array<int, 3> fv{};
            int k = 0;
            for (int i = 0; i < 4; i++)
                if (i != o) fv[k++] = i;
            for (int i = 0; i <= 24; i++)
                for (int j = 0; i + j <= 24; j++) {
                    const double l1 = i/24.0, l2 = j/24.0, l0 = 1.0 - l1 - l2;
                    const Vec3 xi = l0*ref.v[fv[0]] + l1*ref.v[fv[1]] + l2*ref.v[fv[2]];
                    if (std::fabs(c.xf.b(xi)) < 1e-9) continue;   // on the plane: atan2's sign is round-off
                    c.maxJdiff = std::max(c.maxJdiff, std::fabs(m.J(xi) - c.xf.J1(xi)));
                }
        }
        c.walkR = c.maxJdiff > 1e-12;
    }
    WalkParams wpQ = gateWalkParams(allDirichlet(), true, true), wpR = gateWalkParams(allDirichlet(), true, false);
    const HarmonicExtension ext0 = [](const Vec3&, double* phi) {
        for (int k = 0; k < 4; k++) phi[k] = 0.0;
    };
    const TerminalValues tvQ = [](const WalkTerminal& w, double* u) {   // Q_k: crack -lambda_k(xi), faces 0
        if (w.face >= 0) {
            for (int k = 0; k < 4; k++) u[k] = 0.0;
            return;
        }
        const double l[4] = {1.0 - w.x.sum(), w.x[0], w.x[1], w.x[2]};
        for (int k = 0; k < 4; k++) u[k] = -l[k];
    };
    std::vector<QueryPoints> qQ(N), qR(N);
    std::vector<MultiJob> jobsQ, jobsR;
    std::vector<size_t> tetQ, tetR;
    for (size_t t = 0; t < N; t++) {
        const CvRun& c = f.cv[t];
        const std::vector<char> grad(ps[t].size(), 1);
        const uint64_t sQ = mixSeed(kSeed, tagBase + 10*uint64_t(t));
        MultiJob jb;
        jb.tet = &ref;
        jb.cut = &mt[t]->rt.xc.strips;
        jb.nBatches = B;
        if (c.cut) {
            qQ[t].add(ps[t].xi, &grad, sQ);
            jb.pts = &qQ[t];
            jb.nWalks = nQ;
            jb.K = 4;
            jb.wp = wpQ;
            jb.values = tvQ;
            jb.ext = &ext0;
            jobsQ.push_back(jb);
            tetQ.push_back(t);
        }
        if (c.walkR) {
            qR[t].add(ps[t].xi, &grad, mixSeed(sQ, 77));
            jb.pts = &qR[t];
            jb.nWalks = nR;
            jb.K = 1;
            jb.wp = wpR;
            jb.values = [m = mt[t], xf = c.xf](const WalkTerminal& w, double* u) {   // faces J - J_xi, crack side g - J_xi = 0
                u[0] = w.face >= 0 ? m->J(w.x) - xf.J1(w.x) : 0.0;
            };
            jb.ext = nullptr;
            jobsR.push_back(jb);
            tetR.push_back(t);
        }
    }
    auto run = [&](const std::vector<MultiJob>& jobs, const std::vector<size_t>& tets, bool isQ) {
        std::vector<SolveStats> st(jobs.size());
        const auto t0 = clk::now();
        std::vector<MultiSums> sums;
        if (!jobs.empty()) sums = solveMultiLeanJobs(jobs, &st);
        for (size_t j = 0; j < jobs.size(); j++) {
            CvRun& c = f.cv[tets[j]];
            OwnRun& o = isQ ? c.runQ : c.runR;
            o.sums = std::move(sums[j]);
            o.stats = st[j];
            (isQ ? f.walksQ : f.walksR) += st[j].walksRequested;
            f.capped += st[j].walksCapped;
            f.escaped += st[j].walksEscaped;
            if (isQ) {
                c.hasQ = true;
                c.iQ0 = 0;
            } else {
                c.hasR = true;
                c.iR1 = 0;
            }
        }
        return secondsSince(t0);
    };
    f.secQ = run(jobsQ, tetQ, true);
    f.secR = run(jobsR, tetR, false);
    f.nQtets = int(jobsQ.size());
    f.nRtets = int(jobsR.size());
    return f;
}

// quantiles of v (q in [0, 1]; v is sorted here)
double quantile(std::vector<double>& v, double q)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, size_t(q*double(v.size() - 1) + 0.5))];
}

} // namespace

bool buildXfemMesh(const std::vector<Vec3>& verts, const std::vector<std::array<int, 4>>& tets,
                   const std::vector<XfemMeshCrack>& cracks, const XfemMeshParams& p, XfemMeshData& xd)
{
    const auto tStart = clk::now();
    const int NV = int(verts.size()), NT = int(tets.size()), NC = int(cracks.size());
    if (p.quadLevel < 0 || p.quadLevel > 4) throw std::invalid_argument("xfem mesh: quadrature level 0 to 4");
    for (const XfemMeshCrack& k : cracks)
        if (int(k.phi.size()) != NV || int(k.psi.size()) != NV)
            throw std::invalid_argument("xfem mesh: crack " + k.name + ": phi and psi need one value per vertex");
    std::printf("xfem mesh: %d vertices, %d tets, %d crack%s; %u hardware threads (TBB default arena)\n", NV, NT, NC,
                NC == 1 ? "" : "s", std::thread::hardware_concurrency());

    // 1. the tets' frames, and every crack in every tet
    std::vector<TetFrame> T(static_cast<size_t>(NT));
    std::vector<double> L(size_t(NT), 0.0);
    double edgeSum = 0;
    for (int e = 0; e < NT; e++) {
        T[size_t(e)] = frameOn("", tets[size_t(e)], verts.data());
        if (!(T[size_t(e)].detJ > 0)) throw std::invalid_argument("xfem mesh: tet " + std::to_string(e) + " is not positively oriented");
        for (int a = 0; a < 4; a++)
            for (int b = a + 1; b < 4; b++) {
                const double l = (T[size_t(e)].X[size_t(a)] - T[size_t(e)].X[size_t(b)]).norm();
                L[size_t(e)] = std::max(L[size_t(e)], l);
                edgeSum += l;
            }
    }
    const double meanEdge = edgeSum/std::max(1.0, 6.0*NT);
    auto corners = [&](const std::vector<double>& v, int e) {
        const std::array<int, 4>& t = tets[size_t(e)];
        return Eigen::Vector4d(v[size_t(t[0])], v[size_t(t[1])], v[size_t(t[2])], v[size_t(t[3])]);
    };
    std::vector<std::vector<MeshTetCrack>> mc(static_cast<size_t>(NC));
    for (int c = 0; c < NC; c++)
        for (int e = 0; e < NT; e++) {
            try {
                mc[size_t(c)].push_back(meshTetCrack(T[size_t(e)], corners(cracks[size_t(c)].phi, e), corners(cracks[size_t(c)].psi, e),
                                                     L[size_t(e)]));
            } catch (const std::invalid_argument& err) {
                throw std::invalid_argument("xfem mesh: crack " + cracks[size_t(c)].name + ", tet " + std::to_string(e) + ": " + err.what());
            }
        }

    // 2. one crack per tet, per node, and per node's support
    xd = XfemMeshData{};
    xd.tetCrack.assign(size_t(NT), -1);
    xd.tetClass.assign(size_t(NT), 0);
    xd.nodeCrack.assign(size_t(NV), -1);
    std::vector<int> nodeType(size_t(NV), 0);   // 0 none, 1 heaviside, 2 harmonic
    for (int e = 0; e < NT; e++)
        for (int c = 0; c < NC; c++) {
            if (mc[size_t(c)][size_t(e)].cls == 0) continue;
            if (xd.tetCrack[size_t(e)] >= 0)
                throw std::invalid_argument("xfem mesh: tet " + std::to_string(e) + " is cut by both " +
                                            cracks[size_t(xd.tetCrack[size_t(e)])].name + " and " + cracks[size_t(c)].name +
                                            ": one crack per tet; move the cracks apart");
            xd.tetCrack[size_t(e)] = c;
            xd.tetClass[size_t(e)] = mc[size_t(c)][size_t(e)].cls;
        }
    for (int e = 0; e < NT; e++) {
        const int c = xd.tetCrack[size_t(e)];
        if (c < 0) continue;
        for (const int v : tets[size_t(e)]) {
            if (xd.nodeCrack[size_t(v)] >= 0 && xd.nodeCrack[size_t(v)] != c)
                throw std::invalid_argument("xfem mesh: node " + std::to_string(v) + "'s support meets both " +
                                            cracks[size_t(xd.nodeCrack[size_t(v)])].name + " and " + cracks[size_t(c)].name +
                                            ": the XFEM takes one enrichment per node; move the cracks apart");
            xd.nodeCrack[size_t(v)] = c;
            nodeType[size_t(v)] = std::max(nodeType[size_t(v)], xd.tetClass[size_t(e)]);
        }
    }
    for (int e = 0; e < NT; e++) {   // a tet's enriched corners all belong to one crack: its quadrature carries one h and side
        for (const int v : tets[size_t(e)]) {
            const int c = xd.nodeCrack[size_t(v)];
            if (c < 0) continue;
            if (xd.tetCrack[size_t(e)] >= 0 && xd.tetCrack[size_t(e)] != c)
                throw std::invalid_argument("xfem mesh: tet " + std::to_string(e) + " lies in the supports of nodes enriched for " +
                                            cracks[size_t(xd.tetCrack[size_t(e)])].name + " and for " + cracks[size_t(c)].name +
                                            ": the cracks are too close for one enrichment per node; move them apart");
            xd.tetCrack[size_t(e)] = c;
        }
    }
    xd.nodePhi.resize(size_t(NV));
    for (int v = 0; v < NV; v++)
        xd.nodePhi[size_t(v)] = NC == 0 ? 1.0 : cracks[size_t(std::max(0, xd.nodeCrack[size_t(v)]))].phi[size_t(v)];
    xd.minPhi.assign(size_t(NC), kInf);
    for (int e = 0; e < NT; e++)
        if (xd.tetClass[size_t(e)] > 0) {
            const int c = xd.tetCrack[size_t(e)];
            for (const int v : tets[size_t(e)]) xd.minPhi[size_t(c)] = std::min(xd.minPhi[size_t(c)], std::fabs(cracks[size_t(c)].phi[size_t(v)]));
        }
    for (int c = 0; c < NC; c++) {
        int nCut = 0, nFront = 0, nHarm = 0, nHeav = 0, nRegion = 0;
        for (int e = 0; e < NT; e++) {
            if (xd.tetCrack[size_t(e)] != c) continue;
            nCut += xd.tetClass[size_t(e)] == 1 ? 1 : 0;
            nFront += xd.tetClass[size_t(e)] == 2 ? 1 : 0;
            bool harm = false;
            for (const int v : tets[size_t(e)]) harm = harm || nodeType[size_t(v)] == 2;
            nRegion += harm ? 1 : 0;
        }
        for (int v = 0; v < NV; v++)
            if (xd.nodeCrack[size_t(v)] == c) (nodeType[size_t(v)] == 2 ? nHarm : nHeav)++;
        std::printf("crack %s: %d tets cut in two, %d with the front inside; nodes %d harmonic, %d heaviside; region (tets with a "
                    "harmonic corner) %d tets; smallest |phi| at the cut tets' corners %.4g (%.1f%% of the mean edge %.4f)\n",
                    cracks[size_t(c)].name.c_str(), nCut, nFront, nHarm, nHeav, nRegion, xd.minPhi[size_t(c)],
                    100*xd.minPhi[size_t(c)]/meanEdge, meanEdge);
        if (nCut + nFront > 0 && !(xd.minPhi[size_t(c)] > 1e-9*meanEdge))
            throw std::invalid_argument("xfem mesh: crack " + cracks[size_t(c)].name + " passes through a node of a tet it cuts; "
                                        "move it");
    }

    // 3. per crack with a front: the region, its glue, its export pass
    bool ok = true;
    xd.nodeValue.assign(size_t(NV), 0.0);
    std::vector<std::vector<double>> cGlue(static_cast<size_t>(NC));   // per crack: the glued vertex values c_v (its region's vertices)
    std::vector<char> inRegion(size_t(NT), 0);
    double secGlue = 0, secExport = 0;
    long long walksAll = 0;
    size_t nQuadAll = 0;
    int nGluePts = 0, nRegionAll = 0;
    double wErr = 0, hMin = kInf, hMax = -kInf, cornerMax = 0;
    int nShallow = 0, badClass = 0;
    long long badSide = 0;
    std::vector<double> seH, gNorm;
    for (int e = 0; e < NT; e++) {   // class 2 iff the tet holds a piece of the front
        const int c = xd.tetCrack[size_t(e)];
        if (c < 0) continue;
        const MeshTetCrack& m = mc[size_t(c)][size_t(e)];
        badClass += (m.cls == 2) != (m.front.size() == 2) ? 1 : 0;
        if (m.cls == 2 && m.front.size() == 2) xd.front.push_back({e, m.front[0], m.front[1]});
    }
    for (int c = 0; c < NC; c++) {
        std::vector<int> region;
        for (int e = 0; e < NT; e++) {
            if (xd.tetCrack[size_t(e)] != c) continue;
            bool harm = false;
            for (const int v : tets[size_t(e)]) harm = harm || nodeType[size_t(v)] == 2;
            if (harm) region.push_back(e);
        }
        if (region.empty()) continue;
        nRegionAll += int(region.size());
        std::vector<std::array<Vec3, 2>> frontSegs;
        for (const XfemMeshData::Front& f : xd.front)
            if (xd.tetCrack[size_t(f.tet)] == c) frontSegs.push_back({f.a, f.b});
        std::vector<const TetFrame*> rT;
        std::vector<const MeshTetCrack*> rM;
        for (const int e : region) {
            inRegion[size_t(e)] = 1;
            rT.push_back(&T[size_t(e)]);
            rM.push_back(&mc[size_t(c)][size_t(e)]);
        }

        // the glue: check 4's point sets (refined around the tet's crack-face lines and front), Q 32 / R 16 walks
        auto t0 = clk::now();
        std::vector<PointSet> ps;
        for (size_t t = 0; t < region.size(); t++) {
            const int e = region[t];
            ps.push_back(ringPoints(*rT[t], rM[t]->rt, rM[t]->rt.phys, std::max(64, p.gluePoints), kEpsRel*L[size_t(e)], 1));
            nGluePts += int(ps.back().size());
        }
        const int Bg = defaultBatches(kGlueQ);
        const MeshFast glue = runMeshFast(rT, rM, ps, kGlueQ, kGlueR, Bg, 41000 + 1000*uint64_t(c));
        std::vector<int> local(size_t(NV), -1);
        std::vector<int> regionVerts;
        for (const int e : region)
            for (const int v : tets[size_t(e)])
                if (local[size_t(v)] < 0) {
                    local[size_t(v)] = int(regionVerts.size());
                    regionVerts.push_back(v);
                }
        // M c = r: M(v_i, v_j) += S_T(i, j), r(v_i) -= b_T(i), every tet's S_T symmetrized as it is scattered (glueCheck4)
        std::vector<Eigen::Triplet<double>> trip;
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(Eigen::Index(regionVerts.size()));
        for (size_t t = 0; t < region.size(); t++) {
            const Gram g = cvGram(glue.cv[t], ps[t], *rT[t], Bg);
            const std::array<int, 4>& id = rT[t]->global;
            for (int i = 0; i < 4; i++) {
                rhs[local[size_t(id[size_t(i)])]] -= g.G(0, 2 + i);
                for (int j = 0; j < 4; j++)
                    trip.emplace_back(local[size_t(id[size_t(i)])], local[size_t(id[size_t(j)])], 0.5*(g.G(2 + i, 2 + j) + g.G(2 + j, 2 + i)));
            }
        }
        Eigen::SparseMatrix<double> M(Eigen::Index(regionVerts.size()), Eigen::Index(regionVerts.size()));
        M.setFromTriplets(trip.begin(), trip.end());
        const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(M);
        const bool spd = ldlt.info() == Eigen::Success && ldlt.vectorD().minCoeff() > 0;
        ok = ok && spd;
        const Eigen::VectorXd cv = spd ? Eigen::VectorXd(ldlt.solve(rhs)) : Eigen::VectorXd::Zero(rhs.size());
        cGlue[size_t(c)].assign(size_t(NV), 0.0);
        for (size_t i = 0; i < regionVerts.size(); i++) cGlue[size_t(c)][size_t(regionVerts[i])] = cv[Eigen::Index(i)];
        const double secG = secondsSince(t0);
        secGlue += secG;
        walksAll += glue.walksQ + glue.walksR;
        size_t nPts = 0;
        for (const PointSet& s : ps) nPts += s.size();
        std::printf("crack %s: glue over %zu tets and %zu vertices (%zu points; Q %d / R %d walks per point, %d batches; Q in %d "
                    "tets, R in %d): %.2f s (Q walks %.2f s, R walks %.2f s); the %zu x %zu system %s, c in [%.4f, %.4f]\n",
                    cracks[size_t(c)].name.c_str(), region.size(), regionVerts.size(), nPts, kGlueQ, kGlueR, Bg, glue.nQtets,
                    glue.nRtets, secG, glue.secQ, glue.secR, regionVerts.size(), regionVerts.size(),
                    spd ? "positive definite" : "NOT POSITIVE DEFINITE", cv.size() ? cv.minCoeff() : 0.0,
                    cv.size() ? cv.maxCoeff() : 0.0);
        auto localValues = [&](const TetFrame& Tf) {
            const std::vector<double>& g = cGlue[size_t(c)];
            return Eigen::Vector4d(g[size_t(Tf.global[0])], g[size_t(Tf.global[1])], g[size_t(Tf.global[2])], g[size_t(Tf.global[3])]);
        };
        for (const int v : regionVerts)
            xd.nodeValue[size_t(v)] = 0.5*(std::atan2(cracks[size_t(c)].phi[size_t(v)], cracks[size_t(c)].psi[size_t(v)])/std::numbers::pi +
                                           cGlue[size_t(c)][size_t(v)] + 1.0);

        // the export pass: per region tet its quadrature, then (weight 0) a point kCornerOffRel L from each corner
        t0 = clk::now();
        std::vector<XfemTetQuad> quad;
        std::vector<PointSet> pe(region.size());
        for (size_t t = 0; t < region.size(); t++) {
            const int e = region[t];
            const TetFrame& Tf = *rT[t];
            const double band[2] = {kFrontBandRel[0]*L[size_t(e)], kFrontBandRel[1]*L[size_t(e)]};
            quad.push_back(xfemTetQuad(Tf, rM[t]->phi, p.quadLevel, frontSegs, band));
            for (const XfemQuadPoint& q : quad.back().pts) {
                pe[t].X.push_back(q.X);
                pe[t].xi.push_back(Tf.toXi(q.X));
                pe[t].w.push_back(q.w);
            }
            const Vec3 ctr = 0.25*(Tf.X[0] + Tf.X[1] + Tf.X[2] + Tf.X[3]);
            for (int k = 0; k < 4; k++) {
                const Vec3 X = Tf.X[size_t(k)] + kCornerOffRel*L[size_t(e)]*(ctr - Tf.X[size_t(k)]).normalized();
                pe[t].X.push_back(X);
                pe[t].xi.push_back(Tf.toXi(X));
                pe[t].w.push_back(0.0);
            }
        }
        const int Be = defaultBatches(p.walksQ);
        auto whole = [](int n, int b) { return std::max(2*b, (n + b/2)/b*b); };
        const MeshFast ex = runMeshFast(rT, rM, pe, whole(p.walksQ, Be), whole(p.walksR, Be), Be, 51000 + 1000*uint64_t(c));
        walksAll += ex.walksQ + ex.walksR;
        std::map<int, std::pair<double, int>> cornerH;   // vertex -> (sum of h next to it, tets)
        size_t nQuadC = 0;
        for (size_t t = 0; t < region.size(); t++) {
            const int e = region[t];
            const TetFrame& Tf = *rT[t];
            const FieldAt F = cvField(ex.cv[t], pe[t], Tf, 0, localValues(Tf), {}, Be);
            const XfemTetQuad& q = quad[t];
            double sw = 0;
            for (size_t i = 0; i < q.pts.size(); i++) {
                const XfemQuadPoint& qp = q.pts[i];
                const double h = 0.5*(F.H[i] + 1.0);
                xd.quad.push_back({e, qp.side, qp.w, qp.X, h, 0.5*F.G[i]});
                sw += qp.w;
                hMin = std::min(hMin, h);
                hMax = std::max(hMax, h);
                seH.push_back(0.5*std::sqrt(F.vH[i]));
                gNorm.push_back(0.5*F.G[i].norm());
                const double ph = atXi(rM[t]->phi, Tf.toXi(qp.X));
                if ((ph > 0 ? 1 : -1) != qp.side) badSide++;
            }
            nQuadC += q.pts.size();
            const double V = std::fabs(Tf.detJ)/6.0;
            wErr = std::max(wErr, std::fabs(sw - V)/V);
            nShallow += q.nNearShallow;
            for (int k = 0; k < 4; k++) {
                auto& acc = cornerH[Tf.global[size_t(k)]];
                acc.first += 0.5*(F.H[q.pts.size() + size_t(k)] + 1.0);
                acc.second++;
            }
        }
        for (const auto& [v, acc] : cornerH) cornerMax = std::max(cornerMax, std::fabs(acc.first/acc.second - xd.nodeValue[size_t(v)]));
        nQuadAll += nQuadC;
        const double secE = secondsSince(t0);
        secExport += secE;
        std::printf("crack %s: export pass (Q %d / R %d walks per point, %d batches): %zu quadrature points in %zu tets, %.2f s "
                    "(Q walks %.2f s, R walks %.2f s)\n", cracks[size_t(c)].name.c_str(), whole(p.walksQ, Be), whole(p.walksR, Be),
                    Be, nQuadC, region.size(), secE, ex.secQ, ex.secR);
    }

    // 4. the render data: the mesh's boundary triangles, clipped where a crack cuts their tet, and the crack sheets
    const auto tR = clk::now();
    XfemRender rd;
    std::map<std::array<int, 3>, std::pair<int, int>> faces;   // sorted face -> (count, tet * 4 + the opposite corner)
    for (int e = 0; e < NT; e++)
        for (int o = 0; o < 4; o++) {
            std::array<int, 3> f{};
            int k = 0;
            for (int a = 0; a < 4; a++)
                if (a != o) f[size_t(k++)] = tets[size_t(e)][size_t(a)];
            std::sort(f.begin(), f.end());
            auto& entry = faces[f];
            entry.first++;
            entry.second = 4*e + o;
        }
    const double hR = meanEdge/kRenderPerEdge;
    for (const auto& [f, entry] : faces) {
        if (entry.first != 1) continue;
        const int e = entry.second/4, o = entry.second%4;
        const TetFrame& Tf = T[size_t(e)];
        std::array<Vec3, 3> P{};
        int k = 0;
        for (int a = 0; a < 4; a++)
            if (a != o) P[size_t(k++)] = Tf.X[size_t(a)];
        if ((P[1] - P[0]).cross(P[2] - P[0]).dot(P[0] - Tf.X[size_t(o)]) < 0) std::swap(P[1], P[2]);
        const double longest = std::max({(P[1] - P[0]).norm(), (P[2] - P[1]).norm(), (P[0] - P[2]).norm()});
        const int c = xd.tetCrack[size_t(e)];
        const bool split = c >= 0 && xd.tetClass[size_t(e)] > 0;
        // a tet the crack does not cut is drawn whole, on side +1 (a plane no point of it reaches)
        const Vec3 nP = split ? mc[size_t(c)][size_t(e)].n : Vec3::Zero();
        const double dP = split ? mc[size_t(c)][size_t(e)].d : -1.0;
        xfemRenderFace(rd, e, P, nP, dP, std::max(1, int(std::ceil(longest/hR))));
    }
    int badFans = 0;
    for (int e = 0; e < NT; e++)
        if (xd.tetClass[size_t(e)] > 0) {
            const MeshTetCrack& m = mc[size_t(xd.tetCrack[size_t(e)])][size_t(e)];
            badFans += xfemRenderCrack(rd, e, m.rt.phys.polygon, m.n, hR);
        }
    // h: 1 / 0 on the + / - sheet; on the boundary of a region tet H's face data J + sum c lambda (J with the vertex's side
    // where the crack cuts the tet; +-1 on the plane behind the front), 0 elsewhere (no harmonic corner reads it)
    rd.h.assign(rd.X.size(), 0.0);
    for (size_t i = 0; i < rd.X.size(); i++) {
        if (rd.group[i] > 0) {
            rd.h[i] = rd.group[i] == 1 ? 1.0 : 0.0;
            continue;
        }
        const int e = rd.tet[i], c = xd.tetCrack[size_t(e)];
        if (!inRegion[size_t(e)]) continue;
        const TetFrame& Tf = T[size_t(e)];
        const MeshTetCrack& m = mc[size_t(c)][size_t(e)];
        const Vec3 xi = Tf.toXi(rd.X[i]);
        const bool split = xd.tetClass[size_t(e)] > 0;
        double H = rd.onPlane[i] && split && atXi(m.psi, xi) < 0 ? double(rd.side[i]) : m.J(xi, split ? rd.side[i] : 0);
        const double lam[4] = {1.0 - xi.sum(), xi[0], xi[1], xi[2]};
        for (int k = 0; k < 4; k++) H += cGlue[size_t(c)][size_t(Tf.global[size_t(k)])]*lam[k];
        rd.h[i] = 0.5*(H + 1.0);
    }
    for (size_t i = 0; i < rd.X.size(); i++) xd.renderV.push_back({rd.tet[i], rd.side[i], rd.X[i], rd.h[i]});
    xd.renderT = rd.tris;
    const double secRender = secondsSince(tR);

    // 5. the self-checks
    const bool wOk = wErr <= 1e-12, frontOk = nShallow == 0, classOk = badClass == 0, sideOk = badSide == 0, fanOk = badFans == 0;
    ok = ok && wOk && frontOk && classOk && sideOk && fanOk;
    std::printf("self-checks\n");
    if (nRegionAll > 0) {
        std::printf("  weights: max over the region tets of |sum w - V| / V = %.2e (limit 1e-12) %s\n", wErr, wOk ? "ok" : "FAILED");
        std::printf("  front refinement: %d leaves within the inner band of the front not at depth %d %s\n", nShallow,
                    p.quadLevel + 2, frontOk ? "ok" : "FAILED");
        std::printf("  sides: %lld of %zu quadrature points with side != the sign of their tet's phi %s\n", badSide, nQuadAll,
                    sideOk ? "ok" : "FAILED");
        std::printf("  h at the quadrature points: [%.4f, %.4f]; |grad h|: median %.2f, 99%% %.2f; s.e. of h: median %.5f, max %.5f\n",
                    hMin, hMax, quantile(gNorm, 0.5), quantile(gNorm, 0.99), quantile(seH, 0.5), quantile(seH, 1.0));
        std::printf("  h %.0e L from the corners (export pass, mean of the vertex's region tets) against h_v = (J(X_v) + c_v + 1)/2: "
                    "max difference %.4f\n", kCornerOffRel, cornerMax);
    }
    std::printf("  classes: %d tets where class 2 and 'the tet holds a piece of the front' disagree %s\n", badClass,
                classOk ? "ok" : "FAILED");
    std::printf("  crack sheets: fan triangles wound against n: %d %s\n", badFans, fanOk ? "ok" : "FAILED");
    std::printf("render: %zu vertices, %zu triangles (spacing %.4f), %.2f s\n", xd.renderV.size(), xd.renderT.size(), hR, secRender);

    char settings[1024];
    std::snprintf(settings, sizeof settings, "our fast method (no GT) on the mesh; glue: %d points per tet asked (%d), Q %d / R %d "
                  "walks; export points: Q %d / R %d walks; quadrature: level %d, front bands %.3f / %.3f of each tet's longest "
                  "edge, %zu points in %d region tets; %lld walks; time: glue %.2f s, export %.2f s, total %.2f s",
                  p.gluePoints, nGluePts, kGlueQ, kGlueR, p.walksQ, p.walksR, p.quadLevel, kFrontBandRel[0], kFrontBandRel[1],
                  nQuadAll, nRegionAll, walksAll, secGlue, secExport, secondsSince(tStart));
    xd.settings = settings;
    xd.selfChecks = ok;
    std::printf("%s\nself-checks: %s\n", settings, ok ? "PASS" : "FAILED");
    return ok;
}
