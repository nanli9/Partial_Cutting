#pragma once

/// @file
/// XFEM cutting, the geometry half: an oriented level set through one tet.
///
/// From motus/topology/cutting/xfem_cutter.{hpp,cpp}, the per-element part of
/// it: the sign cases, the tessellation of each side into sub-tets, the piece
/// volumes and shape integrals, the cut polygon, and the snapping. The
/// mesh-wide bookkeeping (which nodes are enriched, the numbering) lives in
/// xfem_enrichment.cpp, because here the node types are not motus's: a node
/// next to the crack front is harmonic, not Heaviside, and that is decided
/// from the tet classes the exporter wrote, not from phi alone.
///
/// One scalar per node, phi = n.X - d for the flat crack plane, interpolated by
/// the same N_i as everything else. Inside a tet phi is linear, so its zero set
/// is a plane and cutting the tet is closed form. Three cases, by how many phi
/// are positive:
///
///   4 or 0 -- uncut.
///   3-1    -- the zero set meets the three edges at the lone corner: the lone
///             corner's side is a small tet, the other side a wedge (3 tets).
///   2-2    -- the zero set is a quadrilateral and each side is a wedge.
///
/// **There is no zero case.** A node within the snap tolerance of the plane is
/// pushed onto the positive side first (snap below), so every node is plus or
/// minus and a tet is uncut, 3-1 or 2-2.
///
/// Barycentric throughout: every answer below is a fraction of the parent or a
/// point in its barycentric coordinates, so the parent's positions never
/// enter, and the same points serve the rest-space and the deformed geometry.

#include "xfem/tet_mesh.hpp"
#include "xfem/types.hpp"

#include <vector>

namespace xfem {

/// The positive piece of a tet as sub-tets, corners in barycentric coordinates
/// of the parent. Returns how many were written: 0 when nothing is positive, 1
/// for the whole tet or for the little tet at a lone positive corner, 3 for a
/// wedge -- the base under a lone negative corner, or either side of a 2-2.
/// Assumes no entry of phi is zero.
int positive_piece(const vec4& phi, vec4 (&tets)[3][4]);

/// The negative piece the same way: the positive piece of -phi. With no zero
/// in phi, the two pieces tile the tet exactly.
inline int negative_piece(const vec4& phi, vec4 (&tets)[3][4]) { return positive_piece(-phi, tets); }

/// The volume of a sub-tet given by barycentric corners, as a fraction of the
/// parent's: |det| of the 3x3 of corner differences.
real barycentric_fraction(const vec4& a, const vec4& b, const vec4& c, const vec4& d);

/// The fraction of the tet's volume on the positive side.
real positive_volume_fraction(const vec4& phi);

/// The crossing points of the zero set with the tet's edges, in barycentric
/// coordinates, wound into a polygon: 0 when uncut, 3 for 3-1, 4 for 2-2.
int cut_polygon(const vec4& phi, vec4 (&out)[4]);

/// What a cut tet hands the two-piece path: its two piece volumes and the
/// eight integrals of the corner shape functions over them. [0] / [0..3] is the
/// positive side, [1] / [4..7] the negative; corner a is the parent's own
/// corner a. motus's cut_volume and cut_shape_integral for one cut.
///
///     volume[0] + volume[1]                        = V_e
///     shape_integral[a] + shape_integral[4 + a]    = V_e / 4
///     sum_a shape_integral[4 s + a]                = volume[s]
struct cut_measure {
  real volume[2] = {0, 0};
  real shape_integral[8] = {0, 0, 0, 0, 0, 0, 0, 0};
};

/// Both sides of a tet of volume `volume` cut by phi: one tessellation, the
/// positive side measured and the negative side the remainder, so the
/// invariants above hold exactly (step 4 and 5 of motus's xfem_cutter::cut).
cut_measure measure_cut(const vec4& phi, real volume);

/// phi_i = n.X_i - d for every vertex of the mesh.
std::vector<real> sample_plane(const tet_mesh& mesh, const vec3& n, real d);

/// motus's snap: where |phi_i| < tolerance * (the shortest edge at vertex i),
/// set phi_i to +tolerance * that edge -- a definite side, never zero. Returns
/// how many nodes moved. The exporter's quadrature uses the exact plane, so
/// the default tolerance here is 1e-3 (motus: 1e-2), small enough not to fire
/// on the rings it is meant for; see xfem_enrichment.hpp.
index snap(std::vector<real>& phi, const tet_mesh& mesh, real tolerance);

}  // namespace xfem
