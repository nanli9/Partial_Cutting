#pragma once

/// @file
/// A mesh's partial-cut enrichment straight from partial_cut, in this process:
/// buildXfemMesh (xfem_mesh.h), the --export-xfem solve generalised from
/// check 4's ring to a whole mesh and several cracks, handed over in the xfem
/// library's types. Kept out of the library on purpose, as live_ring.hpp is.

#include "xfem/ring_export.hpp"
#include "xfem/tet_mesh.hpp"
#include "xfem/types.hpp"

#include <functional>
#include <string>
#include <vector>

namespace xfem {

/// One crack in rest space as two level sets: phi, the signed distance to the
/// crack surface (> 0 on its + side, where h = 1), and psi, negative where the
/// surface is cracked (the front is psi = 0). Both are sampled at the mesh's
/// vertices and interpolated linearly, so in every tet the crack is flat with a
/// straight front, and a curved crack becomes piecewise flat.
struct crack_levels {
  std::string name;
  std::function<real(const vec3&)> phi, psi;
};

/// What the XFEM takes from partial_cut for a mesh (xfem_mesh.h, XfemMeshData).
struct mesh_enrichment {
  std::vector<int> tet_class;   ///< per tet, of its crack: 0 missed, 1 cut in two, 2 the front inside
  std::vector<int> tet_crack;   ///< per tet: the crack its enrichment belongs to, -1: none
  std::vector<real> node_phi;   ///< per vertex: phi of the crack it is enriched for (build_enrichment's phi)
  std::vector<real> node_value; ///< per vertex: h there
  std::vector<quadrature_point> quadrature;
  struct front_piece {
    int tet = 0;
    vec3 a = vec3::Zero(), b = vec3::Zero();
  };
  std::vector<front_piece> front;   ///< every front tet's piece of the front: the jump is held at 0 there
  std::vector<render_vertex> render_vertices;
  std::vector<render_triangle> render_triangles;   ///< group 0 the mesh's boundary, 1 / 2 the cracks' + / - sheets
  std::vector<real> min_phi;        ///< per crack: the smallest |phi| at the corners of the tets it cuts
  std::string settings;
};

/// Solves the enrichment with partial_cut's fast method at the given
/// quadrature level. Prints partial_cut's own lines (timings, self-checks);
/// throws std::runtime_error when the self-checks fail, and passes on
/// partial_cut's std::invalid_argument for cracks it refuses (two in one tet,
/// or too close for one enrichment per node).
mesh_enrichment live_mesh(const tet_mesh& mesh, const std::vector<crack_levels>& cracks, int quad_level);

}  // namespace xfem
