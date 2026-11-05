#pragma once

/// @file
/// The file partial_cut's exporter writes and this project reads: one nodal
/// ring, its flat crack, the glued harmonic field h at quadrature points in
/// every tet, and the surfaces to draw.
///
/// Text, whitespace separated, one section per keyword, counts first, 0-based
/// indices, coordinates in the ring's rest space X:
///
///     xfem_ring_export 1
///     case K1                    K1 | K2 | K3
///     ring 2
///     g 1
///     settings <free text to the end of the line>
///     vertices 7                 then x y z per line (v0 = the node)
///     tets 8                     then i j k l, det > 0
///     plane nx ny nz d           n.X = d, |n| = 1, n towards the h = 1 side
///     front F                    F = 0 or 2, then x y z: the front segment's
///                                ends on the ring boundary
///     tet_class 8                then 0 | 1 | 2
///     node_value 7               then h at the vertex
///     quad Q                     then tet side w x y z h gx gy gz
///     render_vertices R          then tet side x y z h
///     render_triangles T         then a b c group (0 outer boundary,
///                                1 crack + sheet, 2 crack - sheet)
///     end
///
/// h = (H / g + 1) / 2, g = grad_X h in rest space (J^-T grad_xi H / (2 g)),
/// side = sign(n.X - d) as +1 / -1. Quadrature points sit strictly inside the
/// leaves of a split-and-refined tessellation, never on the crack; their
/// weights sum to each tet's volume.
///
/// The reader takes the sections in any order and throws std::runtime_error
/// naming the line on anything malformed; check() then tests what the format
/// promises (orientation, weights, sides) and returns the problems as text.

#include "xfem/tet_mesh.hpp"
#include "xfem/types.hpp"

#include <string>
#include <vector>

namespace xfem {

/// One point of the exporter's quadrature: the harmonic field and its rest
/// gradient where the quadrature path evaluates F.
struct quadrature_point {
  int tet = 0;
  int side = 1;  ///< sign(n.X - d)
  real weight = 0;
  vec3 X = vec3::Zero();
  real h = 0;
  vec3 grad_h = vec3::Zero();
};

/// One vertex of the drawn surfaces, carried by the field (1) of its tet on its
/// side: vertices on the crack plane come in two copies, one per side.
struct render_vertex {
  int tet = 0;
  int side = 1;
  vec3 X = vec3::Zero();
  real h = 0;
};

struct render_triangle {
  int v[3] = {0, 0, 0};
  int group = 0;  ///< 0 outer boundary, 1 crack + sheet, 2 crack - sheet
};

struct ring_export {
  int version = 1;
  std::string case_name;
  int ring = 0;
  real g = 1;
  std::string settings;

  tet_mesh mesh;
  vec3 plane_normal = vec3::UnitZ();
  real plane_offset = 0;
  std::vector<vec3> front;  ///< empty (through cut) or the two ends of the front segment

  std::vector<int> tet_class;
  std::vector<real> node_value;
  std::vector<quadrature_point> quadrature;
  std::vector<render_vertex> render_vertices;
  std::vector<render_triangle> render_triangles;

  bool through() const { return front.empty(); }
};

ring_export read_ring_export(const std::string& path);
void write_ring_export(const std::string& path, const ring_export& data);

/// The promises of the format, tested: every problem as one line; empty means
/// the file is what the format says. `weight_tolerance` is relative to each
/// tet's volume (the format promises 1e-12).
std::vector<std::string> check(const ring_export& data, real weight_tolerance = real(1e-9));

}  // namespace xfem
