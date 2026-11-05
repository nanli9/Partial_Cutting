#pragma once

/// @file
/// A finite-element solid cut without remeshing: XFEM, with the three element
/// paths of xfem_enrichment.hpp.
///
/// From motus/discretization/fem/xfem/xfem_system.{hpp,cpp}, simplified: one
/// cpu, no buffers or execution contexts, no assembler classes, no
/// vertex-block interface (no VBD here). What stays is motus's layout and its
/// meaning:
///
///   * The state vector is not positions. Groups [0, V) hold node positions;
///     groups [V, V + E) hold displacement *jumps*, zero at rest. Anything that
///     wants points -- a renderer, the crack opening -- folds the two together
///     through (1) first (deformed_point). Reading x as a point cloud puts every
///     enrichment at the origin and draws a star.
///   * fix_node pins u_i only, the side node i is on; the material across the
///     crack at X_i is u_i + psi a_i and stays free.
///   * The node half of the lumped mass is fem's, unchanged by cutting.
///
/// **One evaluation for all three paths.** On every path F is *linear* in the
/// element's dofs with rest-space coefficients:
///
///     F = sum_j y_j g_j^T,      j over the element's groups (4 nodes, then up
///                               to 4 enrichments), y_j their three reals
///
///     plain       one point, weight V_e:         g_a = gradN_a
///     two_piece   one "point" per piece, weight
///                 the piece volume, side s:      g_a = gradN_a,
///                                                g_{4+a} = psi_a^s gradN_a
///     quadrature  the exporter's points q:       g_a = gradN_a,
///                                                g_{4+a} = psi_a gradN_a
///                                                        + N_a gradpsi_a
///
/// The two-piece row is motus's fold, q_a = u_a + psi_a a_a, written as a
/// coefficient: F_s = sum_a q_a gradN_a^T, constant over the piece, so the
/// piece is integrated exactly by its volume and the energy is motus's sum
/// over the two pieces of V_s W(F_s). The pieces are not elements and have no
/// rest state of their own; both read the parent's gradN. The quadrature row
/// is the product rule on N_a psi_a a_a. The g_j are computed once, at
/// construction, from rest data alone, and every energy, gradient, hessian and
/// step limit afterwards is the same loop over points:
///
///     E   = sum_p w_p W(F_p)
///     g_j = sum_p w_p P(F_p) g_j                          (a 3-vector)
///     H   = sum_p w_p B_p^T C(F_p) B_p,   B_p(3c + r, 3j + r) = g_j[c]
///
/// with C projected to spd per point (per piece on the two-piece path).

#include "xfem/neo_hookean.hpp"
#include "xfem/ring_export.hpp"
#include "xfem/tet4.hpp"
#include "xfem/tet_mesh.hpp"
#include "xfem/types.hpp"
#include "xfem/xfem_enrichment.hpp"

#include <Eigen/SparseCore>

#include <array>
#include <vector>

namespace xfem {

/// One point of an element's energy: a weight and the rest coefficients g_j of
/// F = sum_j y_j g_j^T. On the two-piece path a "point" is a whole piece.
struct integration_point {
  real weight = 0;
  std::array<vec3, 8> g{};
};

struct xfem_element {
  index tet = 0;
  element_path path = element_path::plain;
  /// 4 + the number of enrichments active in this tet.
  int group_count = 4;
  /// The node groups (the tet's corners, in order), then the enrichment groups
  /// of its active corners. Only the first group_count are used.
  std::array<index, 8> groups{};
  /// For each enrichment slot j >= 4: the corner (0..3) it belongs to.
  std::array<int, 8> corner{};
  std::vector<integration_point> points;

  int width() const { return 3 * group_count; }
};

class xfem_system {
 public:
  /// The mesh, the enrichment a crack produced, the exporter's quadrature (only
  /// the quadrature-path tets read it), the material and the density. Builds
  /// the elements and the lumped mass. Throws std::invalid_argument when a
  /// quadrature-path tet has no points, or an enrichment ends with no mass.
  xfem_system(tet_mesh mesh, xfem_enrichment enrichment, const std::vector<quadrature_point>& quadrature,
              neo_hookean material, real density);

  index dofs() const { return enrichment_.dof_count(); }
  const tet_mesh& mesh() const { return mesh_; }
  const xfem_enrichment& enrichment() const { return enrichment_; }
  const neo_hookean& material() const { return material_; }
  real density() const { return density_; }
  const std::vector<xfem_element>& elements() const { return elements_; }
  const tet4::rest_state& rest(index tet) const { return rest_[static_cast<std::size_t>(tet)]; }
  index point_count() const;

  /// The lumped mass, one entry per dof: rho V_e / 4 summed per node (unchanged
  /// by cutting), and for an enrichment rho times the integral of N_a psi_a^2
  /// over its support -- motus's far-side shape integral on the two-piece path,
  /// sum_q w_q N_a psi_a^2 on the quadrature path.
  const vector_x& mass() const { return mass_; }

  /// The load of a uniform body force of g per unit mass (gravity): rho g int N_i on a node group, its lumped mass
  /// times g, and rho g int N_a psi_a on an enrichment group -- the far piece's shape integral with psi's sign there
  /// on the two-piece path, sum_q w_q N_a psi_a on the quadrature path. Not the mass times g on an enrichment: psi
  /// is signed and the mass reads psi^2. Not in motus. The integrator zeros the fixed dofs.
  vector_x body_force(const vec3& g) const;

  /// Node groups at their rest positions, enrichment groups at zero: the
  /// welded rest state.
  vector_x rest_positions() const;

  /// @name The elastic potential
  /// @{
  real potential_energy(const vector_x& x) const;
  /// Zero on the fixed dofs.
  void potential_gradient(const vector_x& x, vector_x& gradient) const;
  /// sum of the element hessians + mass_scale * M (1 / dt^2 for backward
  /// euler, 0 for the potential alone), with the fixed dofs' rows and columns
  /// cleared and 1 on their diagonal.
  void potential_hessian(const vector_x& x, bool project_spd, real mass_scale,
                         Eigen::SparseMatrix<real>& hessian) const;
  /// @}

  /// @name One element, for the self-tests
  /// @{
  real element_energy(index e, const vector_x& x) const;
  /// Local gradient and hessian, element().width() wide, ordered as the
  /// element's groups. No Dirichlet rows.
  void element_derivatives(index e, const vector_x& x, bool project_spd, vector_x& gradient,
                           matrix_x& hessian) const;
  /// The global dofs of element e's local slots.
  std::vector<index> element_dofs(index e) const;
  /// @}

  /// The largest t in (0, 1] for which x + t d keeps det F > 0 at every
  /// quadrature point and every piece; 0 when x is already inverted. motus's
  /// maximum_non_inverting_step, over points instead of pieces.
  real maximum_non_inverting_step(const vector_x& x, const vector_x& direction) const;
  real minimum_determinant(const vector_x& x) const;

  /// @name Boundary conditions
  /// @{
  void fix_dof(index dof);
  /// Hold u_i only: the side node i is on (motus's fix_node).
  void fix_node(index node);
  /// Weld node i: a_i = 0. Quiet on an unenriched node.
  void fix_node_enrichment(index node);
  void clear_fixed_dofs();
  const std::vector<index>& fixed_dofs() const { return fixed_; }
  bool is_fixed(index dof) const { return fixed_flag_[static_cast<std::size_t>(dof)] != 0; }
  /// @}

  /// @name The jump held at zero at given points (the front's constraint)
  /// Not in motus. With the front inside a tet, h jumps by its full 1 right up
  /// to the front, so |grad h| ~ 1 / (2 pi r) there and the quadrature row
  /// N_a grad h makes any jump that is nonzero on the front cost an energy
  /// that grows like log(1 / r) -- infinite in the continuum. Refining the
  /// quadrature then only drives that jump slowly to zero. These constraints
  /// put it there exactly. The jump of a tet, j(X) = sum_a N_a(X) a_a over
  /// its enriched corners, is linear in X, so zero at both ends of the
  /// front's piece in the tet makes it zero along the whole piece; and since
  /// N_a psi_a is continuous across faces, a point on a face two tets share
  /// gives the same rows from either. Each point is three scalar rows, one per
  /// component, all with the same coefficients N_a(X).
  ///
  /// Newton then works in the null space of the rows: an orthonormal basis Z
  /// (dofs x reduced), identity on every dof no row touches and on the fixed
  /// dofs, so a step is Z y with (Z^T H Z) y = -Z^T g. The rest state (a = 0)
  /// satisfies every row, and so does every state a step reaches.
  /// @{
  void zero_jump_at(index tet, const vec3& X);
  void clear_jump_constraints();
  index jump_constraint_points() const { return static_cast<index>(constraints_.size()); }
  /// The basis Z, or nullptr when there is no constraint. Rebuilt when the
  /// constraints or the fixed dofs changed since the last call.
  const Eigen::SparseMatrix<real>* constraint_basis() const;
  /// The rank of the constraint rows (all three components), so the reduced
  /// size is dofs() - constraint_rank().
  index constraint_rank() const;
  /// max over the constrained points of |j(X)|: the jump the constraints hold
  /// at zero. 0 when there are none.
  real constraint_residual(const vector_x& x) const;
  /// @}

  /// The deformed position of rest point X in tet `tet`, on `side`, where the
  /// harmonic field is h: (1), x = sum_a N_a x_a + sum_a N_a psi_a a_a over the
  /// tet's enriched corners. What the renderer draws and the crack opening
  /// measures.
  vec3 deformed_point(const vector_x& x, index tet, const vec3& X, int side, real h) const;

 private:
  mat3 deformation_gradient(const xfem_element& element, const integration_point& point,
                            const vector_x& x) const;
  void build_elements(const std::vector<quadrature_point>& quadrature);
  void build_lumped_mass(const std::vector<quadrature_point>& quadrature);

  tet_mesh mesh_;
  xfem_enrichment enrichment_;
  neo_hookean material_;
  real density_ = 1;
  std::vector<tet4::rest_state> rest_;
  std::vector<xfem_element> elements_;
  vector_x mass_;
  vector_x body_weight_;  ///< per group: rho int N_i (a node), rho int N_a psi_a (an enrichment)
  std::vector<index> fixed_;
  std::vector<char> fixed_flag_;

  // One constrained point: the enrichment groups of the tet's active corners
  // and their N_a(X).
  struct jump_constraint {
    int count = 0;
    std::array<index, 4> group{};
    std::array<real, 4> weight{};
  };
  void build_constraint_basis() const;
  std::vector<jump_constraint> constraints_;
  mutable bool basis_stale_ = true;
  mutable Eigen::SparseMatrix<real> basis_;
  mutable index rank_ = 0;
};

/// The front's constraint for a flat crack with a straight front [a, b]: the
/// segment clipped to every class-2 tet, and zero_jump_at both ends of each
/// piece. Returns the number of points. Throws std::runtime_error when the
/// pieces do not cover the segment (a front the class-2 tets do not hold). A
/// through cut has no front and nothing to call this with.
index zero_jump_on_front(xfem_system& system, const vec3& a, const vec3& b);

}  // namespace xfem
