#pragma once

/// @file
/// A stand-in for partial_cut's exporter: the nodal ring of its check 4 with a
/// flat crack and a made-up h, written in the export format of ring_export.hpp.
///
/// For the self-tests, and so the demo has something to show before the real
/// files are in data/. Nothing here is partial_cut's solution:
///
///   * through: the plane cuts the whole ring (like K3) and h is the exact
///     Heaviside, 1 on the + side and 0 on the - side, gradient 0. Every node
///     is heaviside, every cut tet two-piece.
///   * partial: the crack stops at a straight front inside the ring (like K1)
///     and h is the 2D front solution theta / (2 pi) + 1/2, with
///     theta = atan2(n.X - d, u.X - s): 1 and 0 on the two faces of the crack,
///     1/2 ahead of the front, harmonic off the crack, a jump of exactly 1.
///     It is what partial_cut's Level C face function J is for g = 1, mapped
///     to h.
///
/// The ring, the planes and the quadrature rule follow the export spec: the
/// octahedron |x| + |y| + |z| <= 1 split into the 8 tets (c, +-e_x, +-e_y,
/// +-e_z) around the node c, ordered so det > 0; each tet split by the plane
/// into the +3/-1 or 2+/2- sub-tets, both pieces, refined uniformly, refined
/// further near the front, with a 4-point degree-2 Gauss rule in every leaf.

#include "xfem/ring_export.hpp"
#include "xfem/tet_mesh.hpp"
#include "xfem/types.hpp"

namespace xfem {

struct synthetic_options {
  bool through = false;                    ///< K3-like through cut (exact Heaviside) or K1-like partial
  vec3 node = vec3(0.12, -0.08, 0.10);     ///< c: partial_cut's ring 2
  vec3 normal = vec3(1.0, 0.5, 0.2);       ///< the crack plane's normal, normalised here (K1's)
  real offset = 0.08;                      ///< d in n.X = d (K1's)
  vec3 front_direction = vec3(0, 1, 0);    ///< projected into the plane: the front advance u (K1's)
  real front_t = 0.5;                      ///< the front at this fraction of the section's extent along u
  int quad_level = 1;                      ///< uniform 8-way refinements of every sub-tet
  int front_levels = 1;                    ///< extra refinements of the leaves near the front
  real front_radius = 0.15;                ///< "near the front"
  int render_resolution = 24;              ///< subdivisions per outer boundary edge
};

/// The 8 tets (c, s_x e_x, s_y e_y, s_z e_z) around c, partial_cut's makeRing:
/// v0 = c, v1..v6 = +e_x, -e_x, +e_y, -e_y, +e_z, -e_z; tet o has s_x = -1 if
/// bit 0 of o is set, s_y bit 1, s_z bit 2, its last two corners swapped where
/// that makes det > 0. Local corner 0 is always c, so the face opposite it is
/// on the octahedron's boundary.
tet_mesh nodal_ring_mesh(const vec3& node);

ring_export make_synthetic_ring(const synthetic_options& options);

}  // namespace xfem
