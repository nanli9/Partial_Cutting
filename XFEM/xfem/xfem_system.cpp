#include "xfem/xfem_system.hpp"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace xfem {

namespace {

// B for one point: d vec(F) / d(local dofs), 9 x 24 with only the first
// 3 * group_count columns used. Column-major F, as everywhere: vec(F)[3c + r]
// is F(r, c) = sum_j y_j[r] g_j[c], so its derivative by y_j[r] is g_j[c].
// tet4::deformation_gradient_jacobian is the plain path's case of this, with
// g_a = gradN_a.
matrix<9, 24> point_jacobian(const xfem_element& element, const integration_point& point) {
  matrix<9, 24> B = matrix<9, 24>::Zero();
  for (int j = 0; j < element.group_count; ++j)
  {
    for (int c = 0; c < 3; ++c)
    {
      for (int r = 0; r < 3; ++r) B(3 * c + r, 3 * j + r) = point.g[static_cast<std::size_t>(j)][c];
    }
  }
  return B;
}

// The inversion filter's arithmetic, from motus's xfem_system.cpp: polynomial
// algebra rather than anything about XFEM.

// det of the 3 x 3 whose columns are these.
inline real column_determinant(const vec3& c0, const vec3& c1, const vec3& c2) { return c0.dot(c1.cross(c2)); }

// The coefficients of det(A + t B) in t. A determinant is multilinear in its
// columns, so the t^n coefficient collects the terms that took n columns from
// B: a cubic, which is why a step limit is a root-finding problem.
void inversion_cubic(const mat3& a, const mat3& b, real (&c)[4]) {
  c[0] = column_determinant(a.col(0), a.col(1), a.col(2));
  c[1] = column_determinant(b.col(0), a.col(1), a.col(2)) + column_determinant(a.col(0), b.col(1), a.col(2)) +
         column_determinant(a.col(0), a.col(1), b.col(2));
  c[2] = column_determinant(b.col(0), b.col(1), a.col(2)) + column_determinant(b.col(0), a.col(1), b.col(2)) +
         column_determinant(a.col(0), b.col(1), b.col(2));
  c[3] = column_determinant(b.col(0), b.col(1), b.col(2));
}

inline real cubic_at(const real (&c)[4], real t) { return c[0] + t * (c[1] + t * (c[2] + t * c[3])); }

// The largest t in [0, limit] with the cubic strictly positive on [0, t),
// assuming it is positive at 0. Split at the turning points, so the cubic is
// monotone on each span, and bisect the first span whose far end is not
// positive. Errs low, which is the safe side for a filter.
real positive_until(const real (&c)[4], real limit) {
  real breaks[4];
  int count = 0;
  breaks[count++] = real(0);

  const real qa = real(3) * c[3];
  const real qb = real(2) * c[2];
  const real qc = c[1];
  const real scale = std::max({std::abs(qa), std::abs(qb), std::abs(qc), real(1)});
  const real tiny = std::numeric_limits<real>::epsilon() * scale;

  if (std::abs(qa) <= tiny)
  {
    if (std::abs(qb) > tiny)
    {
      const real t = -qc / qb;
      if (t > real(0) && t < limit) breaks[count++] = t;
    }
  }
  else
  {
    const real discriminant = qb * qb - real(4) * qa * qc;
    if (discriminant > real(0))
    {
      const real root = std::sqrt(discriminant);
      const real first = (-qb - root) / (real(2) * qa);
      const real second = (-qb + root) / (real(2) * qa);
      const real lo = std::min(first, second);
      const real hi = std::max(first, second);
      if (lo > real(0) && lo < limit) breaks[count++] = lo;
      if (hi > real(0) && hi < limit) breaks[count++] = hi;
    }
  }
  breaks[count++] = limit;
  std::sort(breaks, breaks + count);

  for (int i = 0; i + 1 < count; ++i)
  {
    real lo = breaks[i];
    real hi = breaks[i + 1];
    if (cubic_at(c, hi) > real(0)) continue;
    for (int step = 0; step < 60; ++step)
    {
      const real mid = real(0.5) * (lo + hi);
      if (cubic_at(c, mid) > real(0)) lo = mid; else hi = mid;
    }
    return lo;
  }
  return limit;
}

// Stop short of the root rather than on it: at the root the point is exactly
// degenerate and its energy is infinite.
constexpr real inversion_margin = real(0.9);

}  // namespace

xfem_system::xfem_system(tet_mesh mesh, xfem_enrichment enrichment,
                         const std::vector<quadrature_point>& quadrature, neo_hookean material, real density)
    : mesh_(std::move(mesh)), enrichment_(std::move(enrichment)), material_(material), density_(density) {
  if (enrichment_.vertex_count() != mesh_.vertex_count() || enrichment_.tet_count() != mesh_.tet_count())
  {
    throw std::invalid_argument("xfem: xfem_system: the enrichment was built for another mesh");
  }
  if (!(density_ > real(0))) throw std::invalid_argument("xfem: xfem_system: density must be positive");

  // The rest state is per tet, and a cut never touches it: both pieces of a
  // two-piece tet and every quadrature point read the parent's gradN.
  rest_.resize(static_cast<std::size_t>(mesh_.tet_count()));
  for (index e = 0; e < mesh_.tet_count(); ++e)
  {
    vec3 corners[4];
    mesh_.corners(e, corners);
    const real signed_volume = tet_volume(corners[0], corners[1], corners[2], corners[3]);
    if (!(signed_volume > real(0)))
    {
      throw std::invalid_argument("xfem: xfem_system: tet " + std::to_string(e) +
                                  " is degenerate or not positively oriented");
    }
    rest_[static_cast<std::size_t>(e)] = tet4::build_rest_state(corners);
  }

  fixed_flag_.assign(static_cast<std::size_t>(dofs()), 0);
  build_elements(quadrature);
  build_lumped_mass(quadrature);
}

index xfem_system::point_count() const {
  index count = 0;
  for (const xfem_element& element : elements_) count += static_cast<index>(element.points.size());
  return count;
}

// The g_j of every point, once. See the header: this is where the three paths
// differ, and the only place.
void xfem_system::build_elements(const std::vector<quadrature_point>& quadrature) {
  const index tets = mesh_.tet_count();
  std::vector<std::vector<const quadrature_point*>> points_of(static_cast<std::size_t>(tets));
  for (const quadrature_point& q : quadrature) points_of[static_cast<std::size_t>(q.tet)].push_back(&q);

  elements_.assign(static_cast<std::size_t>(tets), xfem_element{});
  for (index e = 0; e < tets; ++e)
  {
    xfem_element& element = elements_[static_cast<std::size_t>(e)];
    const tet4::rest_state& rest = rest_[static_cast<std::size_t>(e)];
    const matrix<4, 3> grad_n = tet4::shape_gradients(rest);
    const tet_class c = enrichment_.class_of(e);

    element.tet = e;
    element.path = enrichment_.path_of(e);
    element.group_count = 4;
    for (int a = 0; a < 4; ++a)
    {
      element.groups[static_cast<std::size_t>(a)] = mesh_.corner(e, a);
      element.corner[static_cast<std::size_t>(a)] = a;
    }

    // The enrichment slots: every corner whose psi can be nonzero here.
    if (element.path != element_path::plain)
    {
      for (int a = 0; a < 4; ++a)
      {
        const index node = mesh_.corner(e, a);
        if (!enrichment_.active(node, c)) continue;
        const std::size_t slot = static_cast<std::size_t>(element.group_count++);
        element.groups[slot] = enrichment_.enrichment_group(node);
        element.corner[slot] = a;
      }
    }

    switch (element.path)
    {
      case element_path::plain: {
        // (a) A plain linear tet: one point, the whole volume.
        integration_point p;
        p.weight = rest.volume;
        for (int a = 0; a < 4; ++a) p.g[static_cast<std::size_t>(a)] = grad_n.row(a).transpose();
        element.points.push_back(p);
        break;
      }

      case element_path::two_piece: {
        // (b) motus's two pieces. Every corner of a class-1 tet is heaviside,
        // so all four enrichment slots are there, in corner order. On the
        // piece on side s, psi_a = xfem_psi(phi_a, s) is a constant per
        // corner, and the fold q_a = x_a + psi_a a_a gives
        // F_s = sum_a x_a gradN_a^T + sum_a a_a psi_a gradN_a^T.
        const index k = enrichment_.cut_of[static_cast<std::size_t>(e)];
        for (int side = 0; side < 2; ++side)
        {
          const bool plus_side = side == 0;
          integration_point p;
          p.weight = enrichment_.cut_volume[static_cast<std::size_t>(2 * k + side)];
          for (int a = 0; a < 4; ++a)
          {
            const vec3 grad = grad_n.row(a).transpose();
            const real psi = xfem_psi(enrichment_.phi[static_cast<std::size_t>(mesh_.corner(e, a))], plus_side);
            p.g[static_cast<std::size_t>(a)] = grad;
            p.g[static_cast<std::size_t>(4 + a)] = psi * grad;
          }
          element.points.push_back(p);
        }
        break;
      }

      case element_path::quadrature: {
        // (c) The exporter's points. psi_a and its gradient at each: h - h_a
        // and grad h for a harmonic corner, the point's side and zero for a
        // heaviside one. g = psi gradN + N gradpsi is the product rule on
        // N_a psi_a a_a.
        const auto& points = points_of[static_cast<std::size_t>(e)];
        if (points.empty())
        {
          throw std::invalid_argument("xfem: xfem_system: tet " + std::to_string(e) +
                                      " takes the quadrature path but the export has no points in it");
        }
        element.points.reserve(points.size());
        for (const quadrature_point* q : points)
        {
          const vec4 N = tet4::shape_values(rest, q->X);
          integration_point p;
          p.weight = q->weight;
          for (int a = 0; a < 4; ++a) p.g[static_cast<std::size_t>(a)] = grad_n.row(a).transpose();
          for (int j = 4; j < element.group_count; ++j)
          {
            const int a = element.corner[static_cast<std::size_t>(j)];
            const index node = mesh_.corner(e, a);
            const real psi = enrichment_.psi(node, c, q->side, q->h);
            const vec3 grad_psi = enrichment_.type_of(node) == node_type::harmonic ? q->grad_h : vec3::Zero();
            p.g[static_cast<std::size_t>(j)] = psi * grad_n.row(a).transpose() + N[a] * grad_psi;
          }
          element.points.push_back(p);
        }
        break;
      }
    }
  }
}

// The lumped mass, from motus's build_lumped_mass. The node half is fem's,
// verbatim, cut or not: lumping gives corner a the integral of N_a over the
// tet, and a cut splits that integral without changing it. The enrichment half
// is rho times the integral of N_a psi_a^2 -- on the two-piece path psi^2 is 1
// on the far piece and 0 on the near one, so it is motus's far-side shape
// integral; on the quadrature path it is the quadrature of the same integrand.
// The magnitude, not the signed row sum of a consistent mass, is the point:
// what a_i carries is positive.
void xfem_system::build_lumped_mass(const std::vector<quadrature_point>& quadrature) {
  mass_ = vector_x::Zero(dofs());
  body_weight_ = vector_x::Zero(enrichment_.group_count());

  for (index e = 0; e < mesh_.tet_count(); ++e)
  {
    const real share = density_ * rest_[static_cast<std::size_t>(e)].volume / real(4);
    for (int a = 0; a < 4; ++a)
    {
      const index node = mesh_.corner(e, a);
      for (int k = 0; k < 3; ++k) mass_[3 * node + k] += share;
      body_weight_[node] += share;
    }
  }

  std::vector<std::vector<const quadrature_point*>> points_of(static_cast<std::size_t>(mesh_.tet_count()));
  for (const quadrature_point& q : quadrature) points_of[static_cast<std::size_t>(q.tet)].push_back(&q);

  for (const xfem_element& element : elements_)
  {
    const index e = element.tet;
    const tet_class c = enrichment_.class_of(e);
    for (int j = 4; j < element.group_count; ++j)
    {
      const int a = element.corner[static_cast<std::size_t>(j)];
      const index node = mesh_.corner(e, a);
      const index group = element.groups[static_cast<std::size_t>(j)];

      // integral: of N_a psi_a^2 (the mass); signed: of N_a psi_a (the body force).
      real integral = real(0);
      real signed_integral = real(0);
      if (element.path == element_path::two_piece)
      {
        integral = enrichment_.far_shape_integral(enrichment_.cut_of[static_cast<std::size_t>(e)], a, node);
        // psi on the far piece: -1 for a node on the + side, +1 for one on the - side.
        signed_integral = xfem_psi(enrichment_.phi[static_cast<std::size_t>(node)],
                                   !(enrichment_.phi[static_cast<std::size_t>(node)] > real(0))) * integral;
      }
      else
      {
        const tet4::rest_state& rest = rest_[static_cast<std::size_t>(e)];
        for (const quadrature_point* q : points_of[static_cast<std::size_t>(e)])
        {
          const real psi = enrichment_.psi(node, c, q->side, q->h);
          const real wN = q->weight * tet4::shape_values(rest, q->X)[a];
          integral += wN * psi * psi;
          signed_integral += wN * psi;
        }
      }
      for (int k = 0; k < 3; ++k) mass_[3 * group + k] += density_ * integral;
      body_weight_[group] += density_ * signed_integral;
    }
  }

  // An enrichment with no mass means its psi vanishes over its whole support:
  // a singular mass matrix in the first step, and a far better error here.
  for (index g = enrichment_.vertex_count(); g < enrichment_.group_count(); ++g)
  {
    if (mass_[3 * g] > real(0)) continue;
    throw std::invalid_argument("xfem: xfem_system: the enrichment of node " +
                                std::to_string(enrichment_.enrichment_node[static_cast<std::size_t>(
                                    g - enrichment_.vertex_count())]) +
                                " ended with zero mass; its psi vanishes over its support");
  }
}

vector_x xfem_system::body_force(const vec3& g) const {
  vector_x f = vector_x::Zero(dofs());
  for (index group = 0; group < enrichment_.group_count(); ++group) f.segment<3>(3 * group) = body_weight_[group] * g;
  return f;
}

vector_x xfem_system::rest_positions() const {
  vector_x x = vector_x::Zero(dofs());
  for (index i = 0; i < mesh_.vertex_count(); ++i) x.segment<3>(3 * i) = mesh_.vertex(i);
  return x;
}

mat3 xfem_system::deformation_gradient(const xfem_element& element, const integration_point& point,
                                       const vector_x& x) const {
  mat3 F = mat3::Zero();
  for (int j = 0; j < element.group_count; ++j)
  {
    const index group = element.groups[static_cast<std::size_t>(j)];
    F.noalias() += x.segment<3>(3 * group) * point.g[static_cast<std::size_t>(j)].transpose();
  }
  return F;
}

real xfem_system::element_energy(index e, const vector_x& x) const {
  const xfem_element& element = elements_[static_cast<std::size_t>(e)];
  real energy = real(0);
  for (const integration_point& point : element.points)
  {
    energy += point.weight * material_.energy_density(deformation_gradient(element, point, x));
  }
  return energy;
}

real xfem_system::potential_energy(const vector_x& x) const {
  real energy = real(0);
  for (index e = 0; e < static_cast<index>(elements_.size()); ++e) energy += element_energy(e, x);
  return energy;
}

std::vector<index> xfem_system::element_dofs(index e) const {
  const xfem_element& element = elements_[static_cast<std::size_t>(e)];
  std::vector<index> dofs;
  for (int j = 0; j < element.group_count; ++j)
  {
    for (int k = 0; k < 3; ++k) dofs.push_back(3 * element.groups[static_cast<std::size_t>(j)] + k);
  }
  return dofs;
}

void xfem_system::element_derivatives(index e, const vector_x& x, bool project_spd, vector_x& gradient,
                                      matrix_x& hessian) const {
  const xfem_element& element = elements_[static_cast<std::size_t>(e)];
  const int width = element.width();

  vector<24> g = vector<24>::Zero();
  matrix<24, 24> h = matrix<24, 24>::Zero();
  for (const integration_point& point : element.points)
  {
    const mat3 F = deformation_gradient(element, point, x);
    const mat3 P = material_.first_piola(F);
    for (int j = 0; j < element.group_count; ++j)
    {
      g.segment<3>(3 * j) += point.weight * (P * point.g[static_cast<std::size_t>(j)]);
    }
    const matrix<9, 24> B = point_jacobian(element, point);
    const mat9 C = material_.hessian(F, project_spd);
    h.noalias() += point.weight * (B.transpose() * (C * B));
  }
  gradient = g.head(width);
  hessian = h.topLeftCorner(width, width);
}

void xfem_system::potential_gradient(const vector_x& x, vector_x& gradient) const {
  gradient = vector_x::Zero(dofs());
  for (const xfem_element& element : elements_)
  {
    for (const integration_point& point : element.points)
    {
      const mat3 P = material_.first_piola(deformation_gradient(element, point, x));
      for (int j = 0; j < element.group_count; ++j)
      {
        const index group = element.groups[static_cast<std::size_t>(j)];
        gradient.segment<3>(3 * group) += point.weight * (P * point.g[static_cast<std::size_t>(j)]);
      }
    }
  }
  // Dirichlet rows, as in motus: a pinned dof reports no force.
  for (const index d : fixed_) gradient[d] = real(0);
}

// Element by element: the local hessian is summed over the element's points
// first and scattered once, so a quadrature tet with hundreds of points costs
// one 24 x 24 scatter. Fixed rows and columns are dropped at the scatter and
// get 1 on the diagonal (motus's apply_dirichlet_rows): the row makes the step
// zero, the column keeps the matrix symmetric, the 1 keeps it non-singular.
void xfem_system::potential_hessian(const vector_x& x, bool project_spd, real mass_scale,
                                    Eigen::SparseMatrix<real>& hessian) const {
  std::vector<Eigen::Triplet<real>> entries;
  entries.reserve(elements_.size() * 24 * 24 + static_cast<std::size_t>(dofs()));

  vector_x g;
  matrix_x h;
  for (index e = 0; e < static_cast<index>(elements_.size()); ++e)
  {
    element_derivatives(e, x, project_spd, g, h);
    const std::vector<index> local = element_dofs(e);
    for (std::size_t r = 0; r < local.size(); ++r)
    {
      if (is_fixed(local[r])) continue;
      for (std::size_t c = 0; c < local.size(); ++c)
      {
        if (is_fixed(local[c])) continue;
        entries.emplace_back(static_cast<int>(local[r]), static_cast<int>(local[c]),
                             h(static_cast<index>(r), static_cast<index>(c)));
      }
    }
  }
  for (index d = 0; d < dofs(); ++d)
  {
    const real diagonal = is_fixed(d) ? real(1) : mass_scale * mass_[d];
    entries.emplace_back(static_cast<int>(d), static_cast<int>(d), diagonal);
  }

  hessian.resize(dofs(), dofs());
  hessian.setFromTriplets(entries.begin(), entries.end());
}

// Per point rather than per tet, and that is not a detail: as motus says of
// its pieces, a thin piece or a point near the front inverts long before the
// whole tet would. Both the state and the direction go through the same g_j,
// since a dof direction moves F by sum_j d_j g_j^T -- reading d as a point
// cloud would return a step that flips a piece.
real xfem_system::maximum_non_inverting_step(const vector_x& x, const vector_x& direction) const {
  real limit = real(1);
  for (const xfem_element& element : elements_)
  {
    for (const integration_point& point : element.points)
    {
      const mat3 A = deformation_gradient(element, point, x);
      const mat3 D = deformation_gradient(element, point, direction);
      real c[4];
      inversion_cubic(A, D, c);
      if (c[0] <= real(0)) return real(0);
      const real root = positive_until(c, real(1));
      // Untouched over the whole unit step: take all of it, so Newton keeps
      // its quadratic tail.
      limit = std::min(limit, root >= real(1) ? real(1) : inversion_margin * root);
    }
  }
  return limit;
}

real xfem_system::minimum_determinant(const vector_x& x) const {
  real smallest = std::numeric_limits<real>::infinity();
  for (const xfem_element& element : elements_)
  {
    for (const integration_point& point : element.points)
    {
      smallest = std::min(smallest, deformation_gradient(element, point, x).determinant());
    }
  }
  return smallest;
}

void xfem_system::fix_dof(index dof) {
  if (dof < 0 || dof >= dofs())
  {
    throw std::invalid_argument("xfem: xfem_system: fixed dof " + std::to_string(dof) + " is outside [0, " +
                                std::to_string(dofs()) + ")");
  }
  if (fixed_flag_[static_cast<std::size_t>(dof)]) return;
  fixed_flag_[static_cast<std::size_t>(dof)] = 1;
  fixed_.push_back(dof);
  basis_stale_ = true;
}

void xfem_system::fix_node(index node) {
  for (int k = 0; k < 3; ++k) fix_dof(3 * node + k);
}

void xfem_system::fix_node_enrichment(index node) {
  const index group = enrichment_.enrichment_group(node);
  if (group < 0) return;
  for (int k = 0; k < 3; ++k) fix_dof(3 * group + k);
}

void xfem_system::clear_fixed_dofs() {
  fixed_.clear();
  std::fill(fixed_flag_.begin(), fixed_flag_.end(), 0);
  basis_stale_ = true;
}

void xfem_system::zero_jump_at(index tet, const vec3& X) {
  if (tet < 0 || tet >= mesh_.tet_count())
  {
    throw std::invalid_argument("xfem: xfem_system: zero_jump_at: tet " + std::to_string(tet) + " does not exist");
  }
  const vec4 N = tet4::shape_values(rest_[static_cast<std::size_t>(tet)], X);
  const tet_class c = enrichment_.class_of(tet);
  jump_constraint row;
  for (int a = 0; a < 4; ++a)
  {
    const index node = mesh_.corner(tet, a);
    if (!enrichment_.active(node, c)) continue;
    row.group[static_cast<std::size_t>(row.count)] = enrichment_.enrichment_group(node);
    row.weight[static_cast<std::size_t>(row.count)] = N[a];
    ++row.count;
  }
  // No enriched corner: the field has no jump here to hold.
  if (row.count == 0) return;
  constraints_.push_back(row);
  basis_stale_ = true;
}

void xfem_system::clear_jump_constraints() {
  constraints_.clear();
  basis_stale_ = true;
}

const Eigen::SparseMatrix<real>* xfem_system::constraint_basis() const {
  if (constraints_.empty()) return nullptr;
  if (basis_stale_) build_constraint_basis();
  return &basis_;
}

index xfem_system::constraint_rank() const {
  if (constraints_.empty()) return 0;
  if (basis_stale_) build_constraint_basis();
  return rank_;
}

// Per component: the rows act on the same groups with the same weights in x, y
// and z, so each component's null space is one small dense SVD over the groups
// the rows touch (a fixed dof among them is held already and left out). A
// point on a face two front tets share appears twice, so the rank is counted,
// not assumed.
void xfem_system::build_constraint_basis() const {
  std::vector<index> touched;
  for (const jump_constraint& row : constraints_)
  {
    for (int k = 0; k < row.count; ++k) touched.push_back(row.group[static_cast<std::size_t>(k)]);
  }
  std::sort(touched.begin(), touched.end());
  touched.erase(std::unique(touched.begin(), touched.end()), touched.end());

  std::vector<char> in_block(static_cast<std::size_t>(dofs()), 0);
  std::vector<std::vector<index>> free_groups(3);
  for (int r = 0; r < 3; ++r)
  {
    for (const index g : touched)
    {
      if (is_fixed(3 * g + r)) continue;
      free_groups[static_cast<std::size_t>(r)].push_back(g);
      in_block[static_cast<std::size_t>(3 * g + r)] = 1;
    }
  }

  std::vector<Eigen::Triplet<real>> entries;
  index column = 0;
  for (index d = 0; d < dofs(); ++d)
  {
    if (!in_block[static_cast<std::size_t>(d)]) entries.emplace_back(static_cast<int>(d), static_cast<int>(column++), real(1));
  }
  rank_ = 0;
  for (int r = 0; r < 3; ++r)
  {
    const std::vector<index>& groups = free_groups[static_cast<std::size_t>(r)];
    if (groups.empty()) continue;
    matrix_x C = matrix_x::Zero(static_cast<index>(constraints_.size()), static_cast<index>(groups.size()));
    for (std::size_t p = 0; p < constraints_.size(); ++p)
    {
      const jump_constraint& row = constraints_[p];
      for (int k = 0; k < row.count; ++k)
      {
        const auto at = std::lower_bound(groups.begin(), groups.end(), row.group[static_cast<std::size_t>(k)]);
        if (at == groups.end() || *at != row.group[static_cast<std::size_t>(k)]) continue;
        C(static_cast<index>(p), static_cast<index>(at - groups.begin())) += row.weight[static_cast<std::size_t>(k)];
      }
    }
    const Eigen::JacobiSVD<matrix_x, Eigen::ComputeFullV> svd(C);
    const vector_x& sigma = svd.singularValues();
    const real tolerance = real(1e-10) * std::max(sigma.size() ? sigma[0] : real(0), real(1));
    index rank = 0;
    while (rank < sigma.size() && sigma[rank] > tolerance) ++rank;
    rank_ += rank;
    const matrix_x& V = svd.matrixV();
    for (index k = rank; k < V.cols(); ++k, ++column)
    {
      for (std::size_t i = 0; i < groups.size(); ++i)
      {
        const real v = V(static_cast<index>(i), k);
        if (v != real(0)) entries.emplace_back(static_cast<int>(3 * groups[i] + r), static_cast<int>(column), v);
      }
    }
  }
  basis_.resize(dofs(), column);
  basis_.setFromTriplets(entries.begin(), entries.end());
  basis_stale_ = false;
}

real xfem_system::constraint_residual(const vector_x& x) const {
  real worst = 0;
  for (const jump_constraint& row : constraints_)
  {
    vec3 jump = vec3::Zero();
    for (int k = 0; k < row.count; ++k)
    {
      jump += row.weight[static_cast<std::size_t>(k)] * x.segment<3>(3 * row.group[static_cast<std::size_t>(k)]);
    }
    worst = std::max(worst, jump.norm());
  }
  return worst;
}

// Clip [a, b] to each front tet in barycentric coordinates: N(t) is linear in
// t along the segment, so each N_k >= 0 bounds t on one side. The front's ends
// lie on the ring's boundary, so N_0 there is 0 to rounding; a small slack
// keeps them.
index zero_jump_on_front(xfem_system& system, const vec3& a, const vec3& b) {
  const tet_mesh& mesh = system.mesh();
  const xfem_enrichment& enrichment = system.enrichment();
  constexpr real slack = real(1e-12);
  index points = 0;
  real covered = 0;
  for (index e = 0; e < mesh.tet_count(); ++e)
  {
    if (enrichment.class_of(e) != tet_class::front) continue;
    const vec4 Na = tet4::shape_values(system.rest(e), a);
    const vec4 Nb = tet4::shape_values(system.rest(e), b);
    real lo = 0, hi = 1;
    for (int k = 0; k < 4; ++k)
    {
      const real slope = Nb[k] - Na[k];
      if (std::abs(slope) < real(1e-300))
      {
        if (Na[k] < -slack) hi = -1;
        continue;
      }
      const real t = (-slack - Na[k]) / slope;
      if (slope > 0) lo = std::max(lo, t); else hi = std::min(hi, t);
    }
    if (!(hi - lo > real(1e-12))) continue;
    system.zero_jump_at(e, a + lo * (b - a));
    system.zero_jump_at(e, a + hi * (b - a));
    points += 2;
    covered += hi - lo;
  }
  if (covered < real(1) - real(1e-9))
  {
    throw std::runtime_error("xfem: zero_jump_on_front: the class-2 tets hold only " + std::to_string(covered) +
                             " of the front segment");
  }
  return points;
}

vec3 xfem_system::deformed_point(const vector_x& x, index tet, const vec3& X, int side, real h) const {
  const vec4 N = tet4::shape_values(rest_[static_cast<std::size_t>(tet)], X);
  const tet_class c = enrichment_.class_of(tet);
  vec3 p = vec3::Zero();
  for (int a = 0; a < 4; ++a)
  {
    const index node = mesh_.corner(tet, a);
    p += N[a] * x.segment<3>(3 * node);
    if (!enrichment_.active(node, c)) continue;
    const real psi = enrichment_.psi(node, c, side, h);
    p += N[a] * psi * x.segment<3>(3 * enrichment_.enrichment_group(node));
  }
  return p;
}

}  // namespace xfem
