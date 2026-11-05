#pragma once

#include "cut_problem.h"

#include <array>
#include <string>
#include <vector>

// The partial-cut enrichment of --export-xfem (coupling_gate.h, buildXfemRing) on any tet mesh, for any set of cracks: the
// XFEM beam (XFEM/examples/beam) calls it in-process. The method is check 4's fast path, unchanged per tet:
//   - every crack is two level sets at the mesh's vertices, phi (the signed distance to the crack surface, > 0 on the
//     h = 1 side) and psi (the surface is cracked where psi < 0, the front is psi = 0), both linear in each tet. So in a
//     tet the crack is a flat polygon with a straight front -- the per-tet cut check 4 solves -- and a curved crack is
//     piecewise flat, continuous across faces (standard XFEM). Tet classes as in the export: 0 the crack misses the tet,
//     1 it cuts the tet in two, 2 the front is inside it;
//   - the nodes' types by the XFEM's per-node rule (harmonic: a class-2 tet in the node's support; heaviside: else a class-1
//     tet; else none). h is needed only in the region of a crack, the tets with a harmonic corner, so a crack cut through
//     the mesh everywhere needs no walks at all;
//   - in every region tet, in its own xi frame with all four faces Dirichlet: H0 = J_xi + R (J_xi: Level C's half-plane of
//     the tet's own plane and front, the control variate; R: crack 0, faces J - J_xi) and P_k = lambda_k + Q_k (Q_k: crack
//     -lambda_k, faces 0, cut tets only), at a ringPoints set of gluePoints points with Q 32 / R 16 walks. The face function
//     is J = atan2(phi_h, psi_h)/pi of the interpolated level sets: continuous across faces and +-1 on the crack;
//   - the glue: E(c) = sum_T (1/2) c_T^T S_T c_T + b_T^T c_T over the region, minimised over the vertex values c_v by one
//     sparse LDLT (check 4's 7 x 7 system grown to the region's vertices). The region's outer faces are treated as check 4
//     treats O's: their vertex values are free unknowns (the natural boundary condition);
//   - the export pass: a quadrature rule per region tet (xfemTetQuad: the tet split by its own linear phi, refined level
//     times and twice more around the crack's front polyline, 4 Gauss points per leaf) solved again with the same xi cuts,
//     face data and glued c at walksQ / walksR walks per point: h = (H + 1)/2, g = grad_X h.
// Cracks are kept apart, as the XFEM takes one enrichment per node: two cracks in one tet, a node enriched for two cracks,
// or a tet in the supports of nodes enriched for two cracks throws std::invalid_argument naming them (move the cracks
// apart). So is a node on a crack it cuts (|phi| below 1e-9 of the edge length: move the crack).
struct XfemMeshCrack {
    std::string name;
    std::vector<double> phi, psi;   // at the mesh's vertices
};
struct XfemMeshParams {
    int quadLevel = 0;              // uniform refinements of every sub-tet of the export quadrature (8 children each)
    int walksQ = 256, walksR = 256; // the export points' walks per point (Q in the cut tets, R where J - J_xi is not 0)
    int gluePoints = 4096;          // the glue's points per tet asked (check 4's ring view)
};
// What the XFEM takes, with the export file's conventions (h, side = the sign of the tet's own phi, g = grad_X h).
struct XfemMeshData {
    std::vector<int> tetCrack;      // per tet: the crack whose enrichment it carries (it cuts the tet, or a corner is
                                    // enriched for it); -1: none
    std::vector<int> tetClass;      // per tet: 0 | 1 | 2 for that crack
    std::vector<int> nodeCrack;     // per vertex: the crack it is enriched for, -1: none
    std::vector<double> nodePhi;    // per vertex: phi of that crack (of crack 0 where none)
    std::vector<double> nodeValue;  // per vertex: h at the vertex, (J(X_v) + c_v + 1)/2, in a region; else 0
    struct Quad { int tet, side; double w; Vec3 X; double h; Vec3 g; };
    std::vector<Quad> quad;         // every tet with a harmonic corner (the XFEM's quadrature path)
    struct Front { int tet; Vec3 a, b; };
    std::vector<Front> front;       // every class-2 tet's piece of its crack's front (X), where the jump is held at 0
    struct RenderVertex { int tet, side; Vec3 X; double h; };
    std::vector<RenderVertex> renderV;
    std::vector<std::array<int, 4>> renderT;   // a, b, c, group (0: the mesh's boundary, 1, 2: the cracks' + and - sheets)
    std::vector<double> minPhi;     // per crack: the smallest |phi| at the corners of the tets it cuts (class 1 or 2)
    std::string settings;
    bool selfChecks = false;        // the printed self-checks passed
};
// Prints the setting, the cracks' classes and node types, the timings and the self-checks; returns whether they passed.
bool buildXfemMesh(const std::vector<Vec3>& verts, const std::vector<std::array<int, 4>>& tets,
                   const std::vector<XfemMeshCrack>& cracks, const XfemMeshParams& p, XfemMeshData& data);
