#pragma once

/// @file
/// The four-node linear tetrahedron.
///
/// From motus/discretization/fem/elements/tet4/tet4.hpp, cpu only. Constant
/// strain over the element, so on an uncut tet one point is exact and F is a
/// single 3x3. A cut keeps that per piece (the two-piece path); a partial cut's
/// enrichment does not, which is why the quadrature path exists at all.

#include "xfem/types.hpp"

#include <Eigen/LU>

#include <cmath>

namespace xfem {

struct tet4 {
  static constexpr int node_count = 4;

  /// The rest-shape data an element caches: enough to get from nodal positions
  /// to a deformation gradient, and to evaluate the shape functions anywhere
  /// in the tet, without touching the mesh again.
  struct rest_state {
    vec3 origin = vec3::Zero();  ///< rest position of corner 0
    mat3 inverse_reference_edges = mat3::Identity();
    real volume = 0.0;
  };

  static rest_state build_rest_state(const vec3 (&reference)[node_count]) {
    rest_state tet;

    mat3 reference_edges;
    reference_edges.col(0) = reference[1] - reference[0];
    reference_edges.col(1) = reference[2] - reference[0];
    reference_edges.col(2) = reference[3] - reference[0];

    tet.origin = reference[0];
    tet.volume = std::abs(reference_edges.determinant()) / real(6);
    tet.inverse_reference_edges = reference_edges.inverse();
    return tet;
  }

  /// F from the current nodal positions.
  static mat3 deformation_gradient(const vec3 (&positions)[node_count], const rest_state& rest) {
    mat3 deformed_edges;
    deformed_edges.col(0) = positions[1] - positions[0];
    deformed_edges.col(1) = positions[2] - positions[0];
    deformed_edges.col(2) = positions[3] - positions[0];
    return deformed_edges * rest.inverse_reference_edges;
  }

  /// The gradients of the four shape functions in rest space, one row per
  /// corner: row a is grad N_a. Constant per element, and the rows sum to zero.
  static matrix<node_count, 3> shape_gradients(const rest_state& rest) {
    // N_1, N_2, N_3 are D_m^-1 (X - X_0), and N_0 is 1 less their sum.
    matrix<node_count, 3> reference_shape_gradient;
    reference_shape_gradient << -1, -1, -1,
                                 1,  0,  0,
                                 0,  1,  0,
                                 0,  0,  1;
    return reference_shape_gradient * rest.inverse_reference_edges;
  }

  /// The four shape functions at a rest point X: its barycentric coordinates in
  /// the tet. Not in motus's tet4 (which never evaluates N away from a corner);
  /// the quadrature path and the renderer both need it.
  static vec4 shape_values(const rest_state& rest, const vec3& X) {
    const vec3 upper = rest.inverse_reference_edges * (X - rest.origin);
    return vec4(real(1) - upper.sum(), upper[0], upper[1], upper[2]);
  }
};

}  // namespace xfem
