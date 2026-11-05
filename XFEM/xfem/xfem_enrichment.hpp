#pragma once

/// @file
/// What a crack hands the simulation: which nodes are enriched and with what
/// function, which path each tet's energy takes, and the two piece volumes of
/// every tet that takes motus's two-piece path.
///
/// The displacement field is motus's (motus/discretization/fem/xfem/
/// xfem_enrichment.hpp), with one more kind of enrichment function:
///
///     u(X) = sum_i N_i(X) u_i  +  sum_{i in E} N_i(X) psi_i(X) a_i        (1)
///
///     heaviside node:  psi_i(X) = H(phi(X)) - H(phi_i)                    (2)
///     harmonic node:   psi_i(X) = h(X) - h(X_i)                           (2')
///
/// phi = n.X - d is the crack plane, H(phi) is 1 on its + side and 0 on its -
/// side. h = (H / g + 1) / 2 is partial_cut's glued harmonic field mapped to the
/// same range: 1 on the + face of the crack, 0 on the - face, so its jump across
/// the crack is 1, like the Heaviside's. Both are *shifted*: psi_i(X_i) = 0, so
/// u(X_i) = u_i and a pin on u_i still pins the material at X_i, and a_i is the
/// displacement jump at node i. a_i = 0 is the welded state.
///
/// **Three tet classes**, from the exporter (ring_export.hpp):
///
///   0 -- the crack misses the tet. The plane may still cross it, beyond the
///        front, where the material is not cracked.
///   1 -- the crack cuts the tet into two pieces.
///   2 -- the front is inside the tet: a partial cut.
///
/// **Node types are chosen per node**, so that every N_i psi_i is continuous
/// across tet faces:
///
///   harmonic   if any tet in the node's support is class 2;
///   heaviside  else, if any tet in its support is class 1;
///   none       else (not enriched).
///
/// Why per node and not per tet. In a fully cut tet the user's rule is plain
/// XFEM, +1 / -1, two pieces; in a front tet psi has to be the harmonic h. But
/// h is not a Heaviside in a fully cut tet next to the front, so switching the
/// function per *tet* would make node i's enrichment jump across the face
/// between a fully cut tet and a front tet -- a fake crack along that face. One
/// function per node cannot jump anywhere but on the crack. The rule keeps the
/// user's choice wherever it is valid: a class-1 tet whose enriched corners are
/// all heaviside runs exactly motus's two-piece path, and when the crack cuts
/// through (K3) every node is heaviside and every cut tet takes that path.
///
/// Where each psi is nonzero:
///
///   * heaviside: only in class-1 tets. In a class-0 tet the crack is absent,
///     even where the plane crosses the tet beyond the front, so psi = 0 there
///     (for a through cut this is motus's "phi has no sign change"). A class-1
///     tet's neighbour across a face the crack crosses is class 1 or 2, never
///     0, so a heaviside node's N_i psi_i stays continuous.
///   * harmonic: wherever h varies -- uncut tets included, which is why the
///     exporter writes quadrature for every tet. h is continuous across tet
///     faces (each tet's H has the same Dirichlet face data as its neighbour's),
///     so N_i (h - h_i) is too.
///
/// **Three element paths**, chosen per tet from the two lists above:
///
///   plain      no corner has a nonzero psi in this tet: a plain linear tet,
///              one F.
///   two_piece  class 1 with every corner heaviside: motus's path, the +3/-1
///              and 2+/2- cases, a constant F on each piece, energy
///              = sum over the two pieces of volume * W(F).
///   quadrature everything else: the field (1) at each of the exporter's
///              quadrature points, F_q = sum_a x_a gradN_a^T + sum_a a_a
///              (psi_a gradN_a + N_a gradpsi_a)^T, energy = sum_q w_q W(F_q).
///              F varies through the tet, hence many points. A heaviside corner
///              in such a tet takes psi from the point's side, with gradpsi = 0.
///
/// In a nodal ring every tet holds the node v0, so when the front is inside the
/// ring v0 is harmonic and every tet takes the quadrature path; when the crack
/// cuts through, the cut tets are two-piece and the rest plain.
///
/// On the front itself h is singular (|grad h| ~ 1 / (2 pi r)), so a jump that
/// does not vanish there has infinite energy in the continuum: the simulation
/// holds it at zero with linear constraints (xfem_system::zero_jump_at,
/// zero_jump_on_front), which is the limit refining the quadrature tends to.

#include "xfem/tet_mesh.hpp"
#include "xfem/types.hpp"

#include <string>
#include <vector>

namespace xfem {

enum class tet_class : int {
  missed = 0,  ///< the crack misses the tet (the plane may cross it beyond the front)
  cut = 1,     ///< the crack cuts the tet into two pieces
  front = 2,   ///< the front is inside the tet
};

enum class node_type : int {
  none = 0,       ///< not enriched
  heaviside = 1,  ///< psi = H(phi) - H(phi_i), motus's
  harmonic = 2,   ///< psi = h - h(X_i), partial_cut's h
};

enum class element_path : int {
  plain = 0,       ///< (a) one F
  two_piece = 1,   ///< (b) motus's two pieces, constant F on each
  quadrature = 2,  ///< (c) F at every quadrature point
};

/// How node types are chosen. per_node is the method; the other two exist for
/// the self-tests.
enum class enrichment_rule {
  per_node,      ///< harmonic / heaviside / none as described above
  all_harmonic,  ///< every node per_node would enrich is harmonic, so every
                 ///< tet with an enriched corner takes the quadrature path
                 ///< (with an exact Heaviside h this must reproduce two_piece)
  none,          ///< nothing enriched: plain FEM on the same mesh
};

const char* name(tet_class c);
const char* name(node_type t);
const char* name(element_path p);

/// H(phi): 1 on the + side, 0 on the - side. phi is never zero after snapping.
inline real heaviside(real phi) { return phi > real(0) ? real(1) : real(0); }

/// psi_a on one side of a cut, from the sign of phi at node a alone: (2) with
/// H(region) substituted in. Zero on the side node a is on, +-1 on the far
/// side. From motus's xfem_enrichment.hpp.
inline real xfem_psi(real phi_a, bool plus_side) {
  const bool positive = phi_a > real(0);
  if (plus_side) return positive ? real(0) : real(1);
  return positive ? real(-1) : real(0);
}

struct enrichment_options {
  enrichment_rule rule = enrichment_rule::per_node;
  /// motus's snap tolerance, in units of the shortest edge at a node (see
  /// xfem_cutter.hpp). Must be positive: it is what rules out a zero phi.
  real snap_tolerance = real(1e-3);
};

/// Flat arrays, as in motus: one entry per vertex, per enrichment, per tet,
/// and per two-piece tet.
struct xfem_enrichment {
  /// @name One entry per vertex
  /// @{
  std::vector<real> phi;        ///< n.X_i - d after snapping; never zero
  std::vector<node_type> type;  ///< heaviside / harmonic / none
  std::vector<real> node_h;     ///< h(X_i), read by harmonic nodes (the export's node_value)
  /// e_i, the enrichment this node owns, or -1. Dense: e_i counts the enriched
  /// nodes before i. Group vertex_count() + e_i holds a_i.
  std::vector<int> enrichment;
  /// @}

  /// @name One entry per enrichment
  /// @{
  std::vector<int> enrichment_node;  ///< the inverse of `enrichment`
  /// @}

  /// @name One entry per tet
  /// @{
  std::vector<tet_class> classes;
  std::vector<element_path> path;
  /// @}

  /// @name One entry per two-piece tet, motus's layout
  /// @{
  std::vector<int> cut_element;          ///< the tet
  std::vector<real> cut_volume;          ///< [2k] + side, [2k + 1] - side
  std::vector<real> cut_shape_integral;  ///< [8k + a] int_+ N_a, [8k + 4 + a] int_- N_a
  std::vector<int> cut_of;               ///< per tet: its k, or -1
  /// @}

  index snapped_nodes = 0;

  index vertex_count() const { return static_cast<index>(enrichment.size()); }
  index enrichment_count() const { return static_cast<index>(enrichment_node.size()); }
  index tet_count() const { return static_cast<index>(classes.size()); }
  index cut_count() const { return static_cast<index>(cut_element.size()); }
  index group_count() const { return vertex_count() + enrichment_count(); }
  index dof_count() const { return 3 * group_count(); }

  bool enriched(index node) const { return enrichment[static_cast<std::size_t>(node)] >= 0; }
  node_type type_of(index node) const { return type[static_cast<std::size_t>(node)]; }
  tet_class class_of(index tet) const { return classes[static_cast<std::size_t>(tet)]; }
  element_path path_of(index tet) const { return path[static_cast<std::size_t>(tet)]; }

  /// The group holding a_i, or -1 when node i is not enriched.
  index enrichment_group(index node) const {
    const int e = enrichment[static_cast<std::size_t>(node)];
    return e < 0 ? index(-1) : vertex_count() + e;
  }

  /// Whether node's psi can be nonzero in a tet of class c: a harmonic node
  /// anywhere, a heaviside node only in a class-1 tet.
  bool active(index node, tet_class c) const {
    const node_type t = type_of(node);
    return t == node_type::harmonic || (t == node_type::heaviside && c == tet_class::cut);
  }

  /// psi of corner `node` at a point of a tet of class c, on `side` (+1 / -1,
  /// the sign of n.X - d there) where the harmonic field is h. The one formula
  /// every path, the mass and the renderer read. Zero where `active` is false.
  real psi(index node, tet_class c, int side, real h) const {
    switch (type_of(node))
    {
      case node_type::harmonic: return h - node_h[static_cast<std::size_t>(node)];
      case node_type::heaviside:
        return c == tet_class::cut ? xfem_psi(phi[static_cast<std::size_t>(node)], side > 0) : real(0);
      default: return real(0);
    }
  }

  /// int_far N_a over two-piece cut k: the integral of corner a's shape
  /// function over the piece on the far side from node `node`, the only one
  /// where its psi is nonzero. Times rho, the lumped mass of a_a there.
  real far_shape_integral(index k, int a, index node) const {
    const bool node_positive = phi[static_cast<std::size_t>(node)] > real(0);
    const std::size_t base = static_cast<std::size_t>(8 * k);
    return node_positive ? cut_shape_integral[base + 4 + static_cast<std::size_t>(a)]
                         : cut_shape_integral[base + static_cast<std::size_t>(a)];
  }

  /// One line per node type and per tet path, for the stats panel.
  std::string summary() const;
};

/// The enrichment of a flat crack n.X = d through `mesh` with the exporter's
/// tet classes (0 / 1 / 2 per tet) and h at the vertices. Throws
/// std::invalid_argument when the classes contradict the plane (a class-1 or
/// class-2 tet with no sign change of phi) or the inputs are mis-sized.
xfem_enrichment build_enrichment(const tet_mesh& mesh, const vec3& n, real d,
                                 const std::vector<int>& classes, const std::vector<real>& node_h,
                                 const enrichment_options& options = {});

/// The same from phi at the vertices instead of a plane: a crack that is flat
/// in each tet but not one plane (a level set interpolated by the N_i, so a
/// curved crack becomes piecewise flat), or several cracks, each vertex taking
/// the phi of the one crack its support meets. Every per-tet rule above reads
/// phi only at the corners of a class-1 or class-2 tet, so the classes and h
/// must come from the same piecewise linear phi.
xfem_enrichment build_enrichment(const tet_mesh& mesh, std::vector<real> phi,
                                 const std::vector<int>& classes, const std::vector<real>& node_h,
                                 const enrichment_options& options = {});

}  // namespace xfem
