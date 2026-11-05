// The self-tests: every one prints PASS or FAIL with its numbers, and the exit
// code is the number of failures.
//
//   xfem_tests                          all tests, on synthetic rings built in memory
//   xfem_tests --write-synthetic [DIR]  also write synthetic_partial.txt and
//                                       synthetic_through.txt (DIR: default tests/data/;
//                                       only the demo reads them, with --data tests/data)
//   xfem_tests --data DIR               look for the real exports (ring2_K*.txt) there
//                                       (default: ./data, then the project's data/)
//
// The synthetic rings (xfem/synthetic_ring.hpp) stand in for partial_cut's
// exporter: a through cut with h the exact Heaviside, and a partial cut with h
// the 2D front solution. Any real export found is checked as well: the format's
// promises, the welded and rigid states, the masses, finite differences on its
// quadrature tets, no seams, mode I and the front's constraint.

#include "xfem/backward_euler.hpp"
#include "xfem/neo_hookean.hpp"
#include "xfem/ring_export.hpp"
#include "xfem/synthetic_ring.hpp"
#include "xfem/tet4.hpp"
#include "xfem/xfem_cutter.hpp"
#include "xfem/xfem_enrichment.hpp"
#include "xfem/xfem_system.hpp"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace xfem;

namespace {

int failures = 0;
int passes = 0;

void report(bool ok, const char* name, const std::string& numbers) {
  std::printf("%s  %-58s %s\n", ok ? "PASS" : "FAIL", name, numbers.c_str());
  (ok ? passes : failures) += 1;
}

std::string fmt(const char* f, double a) {
  char buf[128];
  std::snprintf(buf, sizeof buf, f, a);
  return buf;
}
std::string fmt(const char* f, double a, double b) {
  char buf[160];
  std::snprintf(buf, sizeof buf, f, a, b);
  return buf;
}
std::string fmt(const char* f, double a, double b, double c) {
  char buf[200];
  std::snprintf(buf, sizeof buf, f, a, b, c);
  return buf;
}

const neo_hookean test_material = neo_hookean::from(1.0e4, 0.3);
constexpr real test_density = 1000;

struct built {
  ring_export data;
  xfem_enrichment enrichment;
};

xfem_system make_system(const ring_export& data, enrichment_rule rule) {
  enrichment_options options;
  options.rule = rule;
  xfem_enrichment enrichment =
      build_enrichment(data.mesh, data.plane_normal, data.plane_offset, data.tet_class, data.node_value, options);
  return xfem_system(data.mesh, std::move(enrichment), data.quadrature, test_material, test_density);
}

// A random deformed state: nodes moved by up to `node` off their rest
// positions, jumps up to `jump`, scaled down until every point keeps
// det F > 0.5.
vector_x random_state(const xfem_system& system, std::mt19937& rng, real node, real jump) {
  std::uniform_real_distribution<real> u(-1, 1);
  const vector_x rest = system.rest_positions();
  vector_x x = rest;
  const index V = system.mesh().vertex_count();
  for (index d = 0; d < system.dofs(); ++d) x[d] += (d < 3 * V ? node : jump) * u(rng);
  for (int tries = 0; tries < 40 && system.minimum_determinant(x) < 0.5; ++tries) x = rest + 0.7 * (x - rest);
  return x;
}

real max_abs(const vector_x& v) { return v.size() ? v.lpNorm<Eigen::Infinity>() : real(0); }
real max_abs(const matrix_x& m) { return m.size() ? m.cwiseAbs().maxCoeff() : real(0); }

// ---------------------------------------------------------------------------
// plain FEM, written independently of xfem_system: motus's fem_system on tet4
// ---------------------------------------------------------------------------
struct plain_fem {
  const tet_mesh& mesh;
  neo_hookean material;

  real energy(const vector_x& x) const {
    real E = 0;
    for (index e = 0; e < mesh.tet_count(); ++e)
    {
      vec3 X[4], q[4];
      mesh.corners(e, X);
      for (int a = 0; a < 4; ++a) q[a] = x.segment<3>(3 * mesh.corner(e, a));
      const tet4::rest_state rest = tet4::build_rest_state(X);
      E += rest.volume * material.energy_density(tet4::deformation_gradient(q, rest));
    }
    return E;
  }

  void derivatives(const vector_x& x, vector_x& g, matrix_x& H, bool project) const {
    const index n = 3 * mesh.vertex_count();
    g = vector_x::Zero(n);
    H = matrix_x::Zero(n, n);
    for (index e = 0; e < mesh.tet_count(); ++e)
    {
      vec3 X[4], q[4];
      mesh.corners(e, X);
      for (int a = 0; a < 4; ++a) q[a] = x.segment<3>(3 * mesh.corner(e, a));
      const tet4::rest_state rest = tet4::build_rest_state(X);
      const mat3 F = tet4::deformation_gradient(q, rest);
      const matrix<4, 3> G = tet4::shape_gradients(rest);
      // motus's tet4::deformation_gradient_jacobian
      matrix<9, 12> B = matrix<9, 12>::Zero();
      for (int node = 0; node < 4; ++node)
        for (int column = 0; column < 3; ++column)
          for (int row = 0; row < 3; ++row) B(3 * column + row, 3 * node + row) = G(node, column);
      const vector<12> ge = rest.volume * (B.transpose() * flatten(material.first_piola(F)));
      const matrix<12, 12> He = rest.volume * (B.transpose() * material.hessian(F, project) * B);
      for (int a = 0; a < 4; ++a)
      {
        const index A = mesh.corner(e, a);
        g.segment<3>(3 * A) += ge.segment<3>(3 * a);
        for (int b = 0; b < 4; ++b) H.block<3, 3>(3 * A, 3 * mesh.corner(e, b)) += He.block<3, 3>(3 * a, 3 * b);
      }
    }
  }
};

// ---------------------------------------------------------------------------
// the tests
// ---------------------------------------------------------------------------

// motus's two cutter invariants, f(phi) + f(-phi) = 1 and the shape-integral
// partition, on a 3-1, a 2-2 and an uncut tet.
void test_cutter() {
  const vec4 cases[3] = {vec4(0.3, -0.5, -0.2, -0.9), vec4(0.4, 0.7, -0.3, -0.6), vec4(0.2, 0.5, 0.1, 0.9)};
  real worst = 0;
  for (const vec4& phi : cases)
  {
    worst = std::max(worst, std::abs(positive_volume_fraction(phi) + positive_volume_fraction(-phi) - 1));
    const cut_measure m = measure_cut(phi, 2.0);
    worst = std::max(worst, std::abs(m.volume[0] + m.volume[1] - 2.0));
    for (int s = 0; s < 2; ++s)
    {
      real sum = 0;
      for (int a = 0; a < 4; ++a) sum += m.shape_integral[4 * s + a];
      worst = std::max(worst, std::abs(sum - m.volume[s]));
    }
    // the negative piece's own tessellation agrees with the remainder
    vec4 sub[3][4];
    const int count = negative_piece(phi, sub);
    real negative = 0;
    for (int k = 0; k < count; ++k) negative += barycentric_fraction(sub[k][0], sub[k][1], sub[k][2], sub[k][3]);
    worst = std::max(worst, std::abs(negative * 2.0 - m.volume[1]));
  }
  report(worst < 1e-14, "(0) cutter: f(phi) + f(-phi) = 1, partitions",
         fmt("max deviation %.2e", worst));
}

void test_export_format(const ring_export& data, const char* label) {
  const std::vector<std::string> problems = check(data, 1e-12);
  std::string numbers = std::to_string(data.quadrature.size()) + " points, " +
                        std::to_string(data.render_vertices.size()) + " render vertices";
  for (const std::string& p : problems) numbers += "; " + p;
  std::string title = std::string("(0) ") + label + ": export format promises";
  report(problems.empty(), title.c_str(), numbers);
}

void test_round_trip(const ring_export& data, const std::string& path) {
  write_ring_export(path, data);
  const ring_export back = read_ring_export(path);
  real worst = 0;
  bool same = back.quadrature.size() == data.quadrature.size() &&
              back.render_triangles.size() == data.render_triangles.size() && back.front.size() == data.front.size() &&
              back.tet_class == data.tet_class && back.case_name == data.case_name && back.settings == data.settings;
  if (same)
  {
    for (std::size_t i = 0; i < data.quadrature.size(); ++i)
    {
      worst = std::max({worst, std::abs(back.quadrature[i].weight - data.quadrature[i].weight),
                        (back.quadrature[i].X - data.quadrature[i].X).norm(),
                        (back.quadrature[i].grad_h - data.quadrature[i].grad_h).norm()});
    }
  }
  report(same && worst == 0, "(0) export file: write then read is exact", fmt("max deviation %.1e", worst));
}

// (1) No enrichment: the same numbers as plain FEM, written independently.
void test_uncut_is_plain_fem(const ring_export& data, std::mt19937& rng) {
  const xfem_system system = make_system(data, enrichment_rule::none);
  const plain_fem fem{data.mesh, test_material};
  const vector_x x = random_state(system, rng, 0.08, 0.0);

  vector_x g, gf;
  matrix_x Hf;
  system.potential_gradient(x, g);
  fem.derivatives(x, gf, Hf, true);
  Eigen::SparseMatrix<real> H;
  system.potential_hessian(x, true, 0.0, H);

  const real e = std::abs(system.potential_energy(x) - fem.energy(x)) / std::abs(fem.energy(x));
  const real eg = max_abs(vector_x(g - gf)) / max_abs(gf);
  const real eh = max_abs(matrix_x(matrix_x(H) - Hf)) / max_abs(Hf);
  report(system.dofs() == 3 * data.mesh.vertex_count() && e < 1e-13 && eg < 1e-12 && eh < 1e-12,
         "(1) uncut (no enrichment) equals plain FEM",
         fmt("energy %.1e, gradient %.1e, hessian %.1e (relative)", e, eg, eh));
}

// (2) Finite differences of each path's element energy and gradient, at random
// deformed states with nonzero jumps.
void test_finite_differences(const xfem_system& system, element_path path, std::mt19937& rng, const char* label) {
  const vector_x x = random_state(system, rng, 0.05, 0.05);
  const real h = 1e-6;
  real worst_g = 0, worst_h = 0;
  int tested = 0;
  for (index e = 0; e < static_cast<index>(system.elements().size()); ++e)
  {
    if (system.elements()[static_cast<std::size_t>(e)].path != path) continue;
    ++tested;
    vector_x g;
    matrix_x H;
    system.element_derivatives(e, x, false, g, H);
    const std::vector<index> dofs = system.element_dofs(e);
    vector_x g_fd(g.size());
    matrix_x H_fd(H.rows(), H.cols());
    for (std::size_t k = 0; k < dofs.size(); ++k)
    {
      vector_x xp = x, xm = x;
      xp[dofs[k]] += h;
      xm[dofs[k]] -= h;
      g_fd[static_cast<index>(k)] = (system.element_energy(e, xp) - system.element_energy(e, xm)) / (2 * h);
      vector_x gp, gm;
      matrix_x unused;
      system.element_derivatives(e, xp, false, gp, unused);
      system.element_derivatives(e, xm, false, gm, unused);
      H_fd.col(static_cast<index>(k)) = (gp - gm) / (2 * h);
    }
    worst_g = std::max(worst_g, max_abs(vector_x(g - g_fd)) / std::max(max_abs(g), real(1e-30)));
    worst_h = std::max(worst_h, max_abs(matrix_x(H - H_fd)) / std::max(max_abs(H), real(1e-30)));
  }
  std::string title = std::string("(2) FD, ") + xfem::name(path) + " path (" + label + ")";
  report(tested > 0 && worst_g < 1e-6 && worst_h < 1e-6, title.c_str(),
         fmt("%g tets: gradient %.1e, hessian %.1e (relative)", tested, worst_g, worst_h));
}

// (2) and the whole objective: the incremental potential's gradient and
// hessian, with inertia, a load and pins, away from the predicted point.
void test_objective_fd(const xfem_system& base, std::mt19937& rng, const char* label) {
  xfem_system system = base;
  system.fix_node(1);
  system.fix_node(4);
  const vector_x x = random_state(system, rng, 0.05, 0.05);
  const vector_x predicted = random_state(system, rng, 0.05, 0.05);
  vector_x load = vector_x::Random(system.dofs()) * 100.0;
  for (const index d : system.fixed_dofs()) load[d] = 0;
  const incremental_potential problem(system, predicted, load, 1.0 / 60);

  vector_x g;
  problem.gradient(x, g);
  Eigen::SparseMatrix<real> Hs;
  problem.hessian(x, false, Hs);
  const matrix_x H(Hs);
  const real h = 1e-6;
  vector_x g_fd = vector_x::Zero(g.size());
  matrix_x H_fd = matrix_x::Zero(H.rows(), H.cols());
  for (index k = 0; k < system.dofs(); ++k)
  {
    if (system.is_fixed(k))
    {
      H_fd(k, k) = 1;
      continue;
    }
    vector_x xp = x, xm = x;
    xp[k] += h;
    xm[k] -= h;
    g_fd[k] = (problem.value(xp) - problem.value(xm)) / (2 * h);
    vector_x gp, gm;
    problem.gradient(xp, gp);
    problem.gradient(xm, gm);
    H_fd.col(k) = (gp - gm) / (2 * h);
  }
  // A pinned dof's column is cleared in the hessian (symmetric elimination).
  for (const index d : system.fixed_dofs())
  {
    H_fd.row(d).setZero();
    H_fd(d, d) = 1;
  }
  for (const index d : system.fixed_dofs())
  {
    for (index r = 0; r < system.dofs(); ++r)
    {
      if (r != d) H_fd(r, d) = 0;
    }
  }
  const real eg = max_abs(vector_x(g - g_fd)) / max_abs(g);
  const real eh = max_abs(matrix_x(H - H_fd)) / max_abs(H);
  std::string title = std::string("(2) FD, incremental potential (") + label + ")";
  report(eg < 1e-6 && eh < 1e-6, title.c_str(), fmt("gradient %.1e, hessian %.1e (relative)", eg, eh));
}

// (3) With an exact Heaviside h, the quadrature path is the two-piece path: F
// is constant on each piece, so any rule integrates it exactly.
void test_two_piece_equals_quadrature(const ring_export& through, std::mt19937& rng) {
  const xfem_system two = make_system(through, enrichment_rule::per_node);
  const xfem_system quad = make_system(through, enrichment_rule::all_harmonic);
  index two_piece = 0, quadrature = 0;
  for (const xfem_element& e : two.elements()) two_piece += e.path == element_path::two_piece;
  for (const xfem_element& e : quad.elements()) quadrature += e.path == element_path::quadrature;

  real worst_e = 0, worst_g = 0, worst_h = 0;
  for (int trial = 0; trial < 3; ++trial)
  {
    const vector_x x = random_state(two, rng, 0.05, 0.08);
    vector_x g2, gq;
    two.potential_gradient(x, g2);
    quad.potential_gradient(x, gq);
    Eigen::SparseMatrix<real> H2, Hq;
    two.potential_hessian(x, trial % 2 == 0, 0.0, H2);
    quad.potential_hessian(x, trial % 2 == 0, 0.0, Hq);
    worst_e = std::max(worst_e, std::abs(two.potential_energy(x) - quad.potential_energy(x)) /
                                    std::abs(two.potential_energy(x)));
    worst_g = std::max(worst_g, max_abs(vector_x(g2 - gq)) / max_abs(g2));
    worst_h = std::max(worst_h, max_abs(matrix_x(matrix_x(H2) - matrix_x(Hq))) / max_abs(matrix_x(H2)));
  }
  const real worst_m = max_abs(vector_x(two.mass() - quad.mass())) / max_abs(two.mass());
  char counts[96];
  std::snprintf(counts, sizeof counts, "%d two-piece vs %d quadrature tets, %d snapped; ", int(two_piece),
                int(quadrature), int(two.enrichment().snapped_nodes));
  report(two_piece > 0 && quadrature >= two_piece && worst_e < 1e-11 && worst_g < 1e-10 && worst_h < 1e-10 &&
             worst_m < 1e-11 && two.enrichment().snapped_nodes == 0,
         "(3) exact Heaviside: two-piece == quadrature path",
         counts + fmt("energy %.1e, gradient %.1e, hessian %.1e", worst_e, worst_g, worst_h) +
             fmt(", mass %.1e", worst_m));
}

// (4) a = 0 under an affine deformation: F = A at every point and piece, so
// the energy is V W(A) and the node gradient is plain FEM's.
void test_welded_affine(const ring_export& data, const char* label) {
  const xfem_system system = make_system(data, enrichment_rule::per_node);
  const xfem_system plain = make_system(data, enrichment_rule::none);
  mat3 A;
  A << 1.08, 0.05, -0.03, 0.02, 0.93, 0.06, -0.04, 0.03, 1.05;
  const vec3 b(0.1, -0.2, 0.05);
  vector_x x = vector_x::Zero(system.dofs());
  vector_x xp = vector_x::Zero(plain.dofs());
  real volume = 0;
  for (index i = 0; i < data.mesh.vertex_count(); ++i)
  {
    x.segment<3>(3 * i) = A * data.mesh.vertex(i) + b;
    xp.segment<3>(3 * i) = x.segment<3>(3 * i);
  }
  for (index e = 0; e < data.mesh.tet_count(); ++e) volume += system.rest(e).volume;

  const real exact = volume * test_material.energy_density(A);
  const real e = std::abs(system.potential_energy(x) - exact) / exact;
  vector_x g, gp;
  system.potential_gradient(x, g);
  plain.potential_gradient(xp, gp);
  const real eg = max_abs(vector_x(g.head(gp.size()) - gp)) / max_abs(gp);
  std::string title = std::string("(4) welded a = 0, affine: plain FEM (") + label + ")";
  report(e < 1e-12 && eg < 1e-11, title.c_str(), fmt("energy vs V W(A) %.1e, node gradient %.1e", e, eg));
}

// (5) A rigid rotation with a = 0 stores nothing.
void test_rigid_rotation(const ring_export& data, const char* label) {
  const xfem_system system = make_system(data, enrichment_rule::per_node);
  const mat3 R = Eigen::AngleAxis<real>(0.7, vec3(0.3, -0.5, 0.8).normalized()).toRotationMatrix();
  vector_x x = vector_x::Zero(system.dofs());
  for (index i = 0; i < data.mesh.vertex_count(); ++i) x.segment<3>(3 * i) = R * data.mesh.vertex(i) + vec3(1, 2, 3);
  vector_x g;
  system.potential_gradient(x, g);
  // against the energy and force of a 1% stretch, for a scale
  vector_x xs = system.rest_positions() * 1.01;
  vector_x gs;
  system.potential_gradient(xs, gs);
  const real e = std::abs(system.potential_energy(x)) / system.potential_energy(xs);
  const real eg = max_abs(g) / max_abs(gs);
  std::string title = std::string("(5) rigid rotation, a = 0: zero energy (") + label + ")";
  report(e < 1e-10 && eg < 1e-10, title.c_str(), fmt("energy %.1e, force %.1e (vs a 1%% stretch)", e, eg));
}

// (6) The node half of the lumped mass is rho V, cut or not; and every
// enrichment mass is positive and at most its node's.
void test_masses(const ring_export& data, const char* label) {
  const xfem_system system = make_system(data, enrichment_rule::per_node);
  real volume = 0;
  for (index e = 0; e < data.mesh.tet_count(); ++e) volume += system.rest(e).volume;
  const index V = data.mesh.vertex_count();
  real nodes = 0, enrichments = 0;
  bool ordered = true;
  for (index i = 0; i < V; ++i) nodes += system.mass()[3 * i];
  for (index k = 0; k < system.enrichment().enrichment_count(); ++k)
  {
    const index node = system.enrichment().enrichment_node[static_cast<std::size_t>(k)];
    const real m = system.mass()[3 * (V + k)];
    enrichments += m;
    ordered = ordered && m > 0 && m <= system.mass()[3 * node] * (1 + 1e-12);
  }
  const real e = std::abs(nodes - test_density * volume) / (test_density * volume);
  std::string title = std::string("(6) lumped node masses sum to rho V (") + label + ")";
  report(e < 1e-13 && ordered, title.c_str(),
         fmt("rel %.1e; enrichments %.4g of rho V, each in (0, its node's]", e, enrichments / (test_density * volume)));
}

// (7) A free run (no pins, no load) from a random zero-mean velocity: bounded,
// the discrete energy never grows (backward euler only dissipates), and the
// node momentum is kept to the Newton tolerance.
void test_free_run(const ring_export& data, std::mt19937& rng, const char* label) {
  const xfem_system system = make_system(data, enrichment_rule::per_node);
  const index V = data.mesh.vertex_count();
  mechanical_state s;
  s.x = system.rest_positions();
  s.v = vector_x::Zero(system.dofs());
  std::uniform_real_distribution<real> u(-1, 1);
  for (index d = 0; d < system.dofs(); ++d) s.v[d] = (d < 3 * V ? 0.6 : 0.3) * u(rng);
  vec3 mean = vec3::Zero();
  real total = 0;
  for (index i = 0; i < V; ++i)
  {
    mean += system.mass()[3 * i] * s.v.segment<3>(3 * i);
    total += system.mass()[3 * i];
  }
  for (index i = 0; i < V; ++i) s.v.segment<3>(3 * i) -= mean / total;

  const auto momentum = [&] {
    vec3 p = vec3::Zero();
    for (index i = 0; i < V; ++i) p += system.mass()[3 * i] * s.v.segment<3>(3 * i);
    return p;
  };
  const auto energy = [&] { return 0.5 * s.v.dot(system.mass().cwiseProduct(s.v)) + system.potential_energy(s.x); };

  newton_options options;
  options.gradient_tolerance = 1e-9 * test_material.mu;
  backward_euler integrator(options);
  const real dt = 1.0 / 60;
  const real e0 = energy();
  const real p_scale = total * 0.6;
  real e_max = e0, drift = 0, reach = 0, det = 1;
  bool converged = true, monotone = true;
  real previous = e0;
  for (int step = 0; step < 120; ++step)
  {
    const newton_result r = integrator.step(system, s, dt, vector_x());
    converged = converged && r.converged;
    const real e = energy();
    monotone = monotone && e <= previous * (1 + 1e-9) + 1e-12;
    previous = e;
    e_max = std::max(e_max, e);
    drift = std::max(drift, momentum().norm() / p_scale);
    det = std::min(det, system.minimum_determinant(s.x));
    for (index i = 0; i < V; ++i) reach = std::max(reach, (s.x.segment<3>(3 * i) - data.mesh.vertex(i)).norm());
  }
  std::string title = std::string("(7) free run, 120 steps: bounded, dissipative (") + label + ")";
  report(converged && monotone && reach < 1.0 && drift < 1e-8 && det > 0, title.c_str(),
         fmt("energy %.3g -> %.3g, max node excursion %.3f", e0, previous, reach) +
             fmt(", momentum drift %.1e, min det F %.3f", drift, det));
}

// (8) No fake seams, read off the render data under a random state with
// nonzero jumps. Every render vertex is carried by the field of its own tet,
// so the copies of one rest point X that sit on a face two tets share (the
// outer boundary's edges, the crack polygon's edges on interior faces) must
// land at the same deformed point from either tet: that is the continuity of
// every N_i psi_i across faces the per-node rule exists for. And the two
// copies of a boundary vertex on the plane beyond the front (opposite sides,
// equal h) must land together, since the material there is not cracked.
// Relative to the largest jump. Where a boundary edge meets a crack sheet at
// the crack mouth, the boundary carries the face data and the sheet exactly
// 1 / 0, so those two separate by N psi-difference times the jump; that gap
// is printed, not tested.
void test_no_seams(const ring_export& data, std::mt19937& rng, const char* label) {
  const xfem_system system = make_system(data, enrichment_rule::per_node);
  const vector_x x = random_state(system, rng, 0.05, 0.1);
  const index V = data.mesh.vertex_count();
  const real jump = std::max(max_abs(vector_x(x.tail(x.size() - 3 * V))), real(1e-300));

  const std::size_t R = data.render_vertices.size();
  std::vector<int> group(R, -1);
  for (const render_triangle& t : data.render_triangles)
  {
    for (const int v : t.v) group[static_cast<std::size_t>(v)] = t.group;
  }
  std::vector<vec3> deformed(R);
  for (std::size_t i = 0; i < R; ++i)
  {
    const render_vertex& v = data.render_vertices[i];
    deformed[i] = system.deformed_point(x, v.tet, v.X, v.side, v.h);
  }

  // Bucket by rest position, to 1e-9.
  std::map<std::array<long long, 3>, std::vector<std::size_t>> at;
  for (std::size_t i = 0; i < R; ++i)
  {
    if (group[i] < 0) continue;
    const vec3& X = data.render_vertices[i].X;
    at[{std::llround(X[0] * 1e9), std::llround(X[1] * 1e9), std::llround(X[2] * 1e9)}].push_back(i);
  }

  real seam = 0, beyond = 0, mouth = 0;
  index seam_pairs = 0, beyond_pairs = 0, mouth_pairs = 0;
  for (const auto& [key, ids] : at)
  {
    for (std::size_t p = 0; p < ids.size(); ++p)
    {
      for (std::size_t q = p + 1; q < ids.size(); ++q)
      {
        const render_vertex& a = data.render_vertices[ids[p]];
        const render_vertex& b = data.render_vertices[ids[q]];
        const real gap = (deformed[ids[p]] - deformed[ids[q]]).norm() / jump;
        const bool boundary_pair = group[ids[p]] == 0 && group[ids[q]] == 0;
        if (a.side == b.side && group[ids[p]] == group[ids[q]] && a.tet != b.tet)
        {
          seam = std::max(seam, gap);
          ++seam_pairs;
        }
        else if (a.side != b.side && boundary_pair && std::abs(a.h - b.h) < real(0.5))
        {
          beyond = std::max(beyond, gap);
          ++beyond_pairs;
        }
        else if (a.side == b.side && (group[ids[p]] == 0) != (group[ids[q]] == 0))
        {
          mouth = std::max(mouth, gap);
          ++mouth_pairs;
        }
      }
    }
  }
  std::string title = std::string("(8) no seams along tet faces (") + label + ")";
  report(seam_pairs > 0 && seam < 1e-9 && beyond < 1e-9, title.c_str(),
         fmt("%.0f face pairs %.1e, ", double(seam_pairs), seam) +
             fmt("%.0f beyond-front pairs %.1e; mouth gap %.1e", double(beyond_pairs), beyond, mouth) +
             fmt(" (%.0f pairs, x the jump)", double(mouth_pairs)));
}

// Mode I, as the demo loads it: the outer vertices on the - side pinned (u_i
// only), the + side ones pulled along n, the demo's material, the front's jump
// held at zero or left free. With a front the pull is solved quasi-statically
// (backward euler with dt = 100 s, the load ramped over 4 steps); with no front
// the + half is a free body: 60 real steps (dt = 1/60 s). Measured on the
// crack-sheet pairs: their opening binned by the distance r to the front line,
// the smallest gap along n, and the jump on the front itself.
struct mode_one_result {
  bool converged = true;
  index iterations = 0;
  real bin_mean[5] = {0, 0, 0, 0, 0};
  index bin_count[5] = {0, 0, 0, 0, 0};
  real widest = 0, smallest_gap = 0, min_det = 0, front_jump = 0, constraint_residual = 0;
  index front_points = 0, points = 0;
};

mode_one_result solve_mode_one(const ring_export& data, real load, bool zero_front) {
  enrichment_options options;
  xfem_enrichment enrichment =
      build_enrichment(data.mesh, data.plane_normal, data.plane_offset, data.tet_class, data.node_value, options);
  xfem_system system(data.mesh, std::move(enrichment), data.quadrature, neo_hookean::from(2.0e4, 0.3), 1000);
  mode_one_result out;
  out.points = static_cast<index>(data.quadrature.size());

  const vec3 n = data.plane_normal;
  std::vector<index> plus;
  for (index v = 1; v < data.mesh.vertex_count(); ++v)  // v0 is the ring's node, inside
  {
    if (n.dot(data.mesh.vertex(v)) - data.plane_offset > 0) plus.push_back(v);
    else system.fix_node(v);
  }
  if (zero_front && !data.through()) out.front_points = zero_jump_on_front(system, data.front[0], data.front[1]);

  // The two copies of every crack-sheet vertex.
  std::vector<int> group(data.render_vertices.size(), -1);
  for (const render_triangle& t : data.render_triangles)
  {
    for (const int v : t.v) group[static_cast<std::size_t>(v)] = t.group;
  }
  std::map<std::array<long long, 4>, std::size_t> minus_copy;
  const auto key = [&](const render_vertex& v) {
    return std::array<long long, 4>{v.tet, std::llround(v.X[0] * 1e9), std::llround(v.X[1] * 1e9),
                                    std::llround(v.X[2] * 1e9)};
  };
  for (std::size_t i = 0; i < group.size(); ++i)
  {
    if (group[i] == 2) minus_copy[key(data.render_vertices[i])] = i;
  }
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  for (std::size_t i = 0; i < group.size(); ++i)
  {
    if (group[i] != 1) continue;
    const auto found = minus_copy.find(key(data.render_vertices[i]));
    if (found != minus_copy.end()) pairs.push_back({i, found->second});
  }
  const auto opening = [&](const vector_x& x, std::size_t k) {
    const render_vertex& a = data.render_vertices[pairs[k].first];
    const render_vertex& b = data.render_vertices[pairs[k].second];
    return (system.deformed_point(x, a.tet, a.X, a.side, a.h) - system.deformed_point(x, b.tet, b.X, b.side, b.h))
        .norm();
  };

  newton_options newton;
  newton.gradient_tolerance = 1e-8 * 2.0e4;
  backward_euler integrator(newton);
  mechanical_state s;
  s.x = system.rest_positions();
  s.v = vector_x::Zero(system.dofs());
  const auto force = [&](real factor) {
    vector_x f = vector_x::Zero(system.dofs());
    for (const index v : plus) f.segment<3>(3 * v) = factor * load / real(plus.size()) * n;
    return f;
  };

  if (data.through())
  {
    for (int step = 1; step <= 60; ++step)
    {
      const newton_result r = integrator.step(system, s, 1.0 / 60, force(std::min(step / 30.0, 1.0)));
      out.converged = out.converged && r.converged;
      out.iterations = std::max(out.iterations, r.iterations);
    }
  }
  else
  {
    for (int step = 1; step <= 4; ++step)
    {
      s.v.setZero();
      const newton_result r = integrator.step(system, s, 100.0, force(step / 4.0));
      out.converged = out.converged && r.converged;
      out.iterations = std::max(out.iterations, r.iterations);
    }
  }
  out.min_det = system.minimum_determinant(s.x);
  out.constraint_residual = system.constraint_residual(s.x);

  // Distance to the front line, in the plane (none for a through cut: one bin).
  const vec3 f0 = data.through() ? vec3::Zero() : data.front[0];
  const vec3 along = data.through() ? vec3::UnitX() : vec3((data.front[1] - data.front[0]).normalized());
  const real edges[] = {0.05, 0.1, 0.2, 0.4, 1e30};
  real sum[5] = {0, 0, 0, 0, 0};
  out.smallest_gap = pairs.empty() ? 0 : 1e30;
  for (std::size_t k = 0; k < pairs.size(); ++k)
  {
    const render_vertex& a = data.render_vertices[pairs[k].first];
    const render_vertex& b = data.render_vertices[pairs[k].second];
    const vec3 d = a.X - f0;
    const real r = data.through() ? real(1e29) : (d - d.dot(along) * along).norm();
    int bin = 0;
    while (r >= edges[bin]) ++bin;
    const real o = opening(s.x, k);
    sum[bin] += o;
    ++out.bin_count[bin];
    out.widest = std::max(out.widest, o);
    const vec3 apart = system.deformed_point(s.x, a.tet, a.X, a.side, a.h) - system.deformed_point(s.x, b.tet, b.X, b.side, b.h);
    out.smallest_gap = std::min(out.smallest_gap, apart.dot(n));
  }
  for (int bin = 0; bin < 5; ++bin) out.bin_mean[bin] = out.bin_count[bin] ? sum[bin] / real(out.bin_count[bin]) : 0;

  // The jump on the front, sampled along it in the front tets: + face (h = 1)
  // less - face (h = 0), sum_a N_a a_a.
  if (!data.through())
  {
    for (int i = 0; i <= 200; ++i)
    {
      const vec3 X = data.front[0] + (i / 200.0) * (data.front[1] - data.front[0]);
      for (index e = 0; e < data.mesh.tet_count(); ++e)
      {
        if (system.enrichment().class_of(e) != tet_class::front) continue;
        if (tet4::shape_values(system.rest(e), X).minCoeff() < -1e-12) continue;
        const vec3 j = system.deformed_point(s.x, e, X, 1, 1.0) - system.deformed_point(s.x, e, X, -1, 0.0);
        out.front_jump = std::max(out.front_jump, j.norm());
      }
    }
  }
  return out;
}

std::string opening_profile(const mode_one_result& m) {
  std::string profile = "mean opening by r:";
  const char* names[] = {"<0.05", "<0.1", "<0.2", "<0.4", ">=0.4"};
  for (int b = 0; b < 5; ++b)
  {
    if (m.bin_count[b] == 0) continue;
    char buf[48];
    std::snprintf(buf, sizeof buf, " %s %.4f", names[b], m.bin_mean[b]);
    profile += buf;
  }
  return profile;
}

// (9) Mode I as the demo runs it (the front's jump held at zero). With a front
// the crack must open behind it and stay shut at it: the opening has to grow
// away from the front. With no front the + half must separate with little
// strain.
void test_mode_one(const ring_export& data, real load, const char* label) {
  const mode_one_result m = solve_mode_one(data, load, true);
  std::string title = std::string("(9) mode I: ") +
                      (data.through() ? "the cut + half separates (" : "opens behind the front (") + label + ")";
  if (data.through())
  {
    report(m.converged && m.widest > 0.03 && m.min_det > 0.95, title.c_str(),
           fmt("opening %.3f after 1 s, min det F %.3f, at most %.0f Newton iterations", m.widest, m.min_det,
               double(m.iterations)));
    return;
  }
  real nearest = -1, farthest = -1;
  for (int b = 0; b < 5; ++b)
  {
    if (m.bin_count[b] == 0) continue;
    if (nearest < 0) nearest = m.bin_mean[b];
    farthest = m.bin_mean[b];
  }
  report(m.converged && m.widest > 0.01 && nearest < 0.5 * farthest && m.min_det > 0, title.c_str(),
         opening_profile(m) + fmt("; smallest gap %+.4f, min det F %.3f, at most %.0f Newton iterations",
                                  m.smallest_gap, m.min_det, double(m.iterations)));
}

// (10) The front's constraint. Its basis: orthonormal, and every column holds
// the jump at zero at every constrained point. The solve of (9): zero jump at
// the constrained points (the clipped front's ends) and so, j being linear in
// each tet, all along the front, sampled at 201 points.
void test_front_constraint(const ring_export& data, real load, const char* label) {
  xfem_system system = make_system(data, enrichment_rule::per_node);
  const index points = zero_jump_on_front(system, data.front[0], data.front[1]);
  const Eigen::SparseMatrix<real>* Z = system.constraint_basis();
  real orthonormal = 1e30, in_null_space = 1e30;
  if (Z)
  {
    const matrix_x Zd(*Z);
    orthonormal = max_abs(matrix_x(Zd.transpose() * Zd - matrix_x::Identity(Zd.cols(), Zd.cols())));
    in_null_space = 0;
    for (index k = 0; k < Zd.cols(); ++k) in_null_space = std::max(in_null_space, system.constraint_residual(Zd.col(k)));
  }
  const mode_one_result m = solve_mode_one(data, load, true);
  const real jump = std::max(m.widest, real(1e-300));
  std::string title = std::string("(10) front: jump held at 0 along it (") + label + ")";
  report(Z && points >= 4 && orthonormal < 1e-12 && in_null_space < 1e-12 && m.converged &&
             m.constraint_residual < 1e-10 * jump && m.front_jump < 1e-10 * jump,
         title.c_str(),
         fmt("%.0f points, rank %.0f; |Z^T Z - I| %.1e", double(points), double(system.constraint_rank()), orthonormal) +
             fmt(", |C Z| %.1e; after mode I: at the points %.1e, along the front %.1e (x the widest opening)",
                 in_null_space, m.constraint_residual / jump, m.front_jump / jump));
}

// (11) Refining the quadrature converges with the front's jump held at zero,
// and does not with it free: synthetic partial rings at quadrature levels 0,
// 1, 2, mode I quasi-statically, the mean opening beyond 0.4 of the front and
// the smallest gap. Free, the log(1 / r_min) energy at the front stiffens the
// jump at every level; held, the integrand is bounded and the levels agree.
void test_quadrature_convergence(real load) {
  real held[3], free_[3], held_gap[3], free_gap[3];
  index count[3];
  for (int level = 0; level < 3; ++level)
  {
    synthetic_options options;
    options.quad_level = level;
    const ring_export data = make_synthetic_ring(options);
    const mode_one_result h = solve_mode_one(data, load, true);
    const mode_one_result f = solve_mode_one(data, load, false);
    held[level] = h.bin_mean[4];
    free_[level] = f.bin_mean[4];
    held_gap[level] = h.smallest_gap;
    free_gap[level] = f.smallest_gap;
    count[level] = h.points;
  }
  const real held_change = std::abs(held[2] - held[1]) / held[2];
  const real free_change = std::abs(free_[2] - free_[1]) / free_[2];
  std::printf("      levels 0, 1, 2 (%d, %d, %d points): opening beyond 0.4 held %.4f %.4f %.4f, free %.4f %.4f %.4f;\n"
              "      smallest gap held %+.4f %+.4f %+.4f, free %+.4f %+.4f %+.4f\n",
              int(count[0]), int(count[1]), int(count[2]), held[0], held[1], held[2], free_[0], free_[1], free_[2],
              held_gap[0], held_gap[1], held_gap[2], free_gap[0], free_gap[1], free_gap[2]);
  report(held_change < 0.01 && held_change < free_change, "(11) quadrature levels 1 -> 2 converge, jump held at 0",
         fmt("opening beyond 0.4 changes by %.2f%% held, %.2f%% free", 100 * held_change, 100 * free_change));
}

}  // namespace

int main(int argc, char** argv) {
  std::string data_dir;
  std::string write_dir;
  int quad_level = -1, front_levels = -1;
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if (arg == "--data" && i + 1 < argc) data_dir = argv[++i];
    else if (arg == "--quad-level" && i + 1 < argc) quad_level = std::stoi(argv[++i]);
    else if (arg == "--front-levels" && i + 1 < argc) front_levels = std::stoi(argv[++i]);
    else if (arg == "--write-synthetic")
    {
      write_dir = (i + 1 < argc && argv[i + 1][0] != '-') ? argv[++i] : "";
      if (write_dir.empty()) write_dir = XFEM_SYNTHETIC_DIR;
    }
    else
    {
      std::fprintf(stderr, "usage: %s [--data DIR] [--write-synthetic [DIR]] [--quad-level L] [--front-levels F]\n",
                   argv[0]);
      return 2;
    }
  }
  namespace fs = std::filesystem;
  if (data_dir.empty())
  {
    data_dir = fs::exists("data") ? "data" : XFEM_DATA_DIR;
  }

  std::mt19937 rng(20261104);

  synthetic_options partial_options;
  if (quad_level >= 0) partial_options.quad_level = quad_level;
  if (front_levels >= 0) partial_options.front_levels = front_levels;
  synthetic_options through_options = partial_options;
  through_options.through = true;
  const ring_export through = make_synthetic_ring(through_options);
  const ring_export partial = make_synthetic_ring(partial_options);

  const xfem_system through_system = make_system(through, enrichment_rule::per_node);
  const xfem_system partial_system = make_system(partial, enrichment_rule::per_node);
  std::printf("synthetic through: %s\n  %s\n", through.settings.c_str(),
              through_system.enrichment().summary().c_str());
  std::printf("synthetic partial: %s\n  %s\n\n", partial.settings.c_str(),
              partial_system.enrichment().summary().c_str());

  test_cutter();
  test_export_format(through, "synthetic through");
  test_export_format(partial, "synthetic partial");
  {
    const fs::path scratch = fs::temp_directory_path() / "xfem_round_trip.txt";
    test_round_trip(partial, scratch.string());
    std::error_code ignored;
    fs::remove(scratch, ignored);
  }
  test_uncut_is_plain_fem(partial, rng);
  test_finite_differences(through_system, element_path::plain, rng, "synthetic through");
  test_finite_differences(through_system, element_path::two_piece, rng, "synthetic through");
  test_finite_differences(partial_system, element_path::quadrature, rng, "synthetic partial");
  test_objective_fd(through_system, rng, "synthetic through");
  test_objective_fd(partial_system, rng, "synthetic partial");
  test_two_piece_equals_quadrature(through, rng);
  test_welded_affine(through, "synthetic through");
  test_welded_affine(partial, "synthetic partial");
  test_rigid_rotation(through, "synthetic through");
  test_rigid_rotation(partial, "synthetic partial");
  test_masses(through, "synthetic through");
  test_masses(partial, "synthetic partial");
  test_free_run(through, rng, "synthetic through");
  test_free_run(partial, rng, "synthetic partial");
  test_no_seams(through, rng, "synthetic through");
  test_no_seams(partial, rng, "synthetic partial");
  test_mode_one(through, 400, "synthetic through");
  test_mode_one(partial, 400, "synthetic partial");
  test_front_constraint(partial, 400, "synthetic partial");
  test_quadrature_convergence(400);

  // The real exports, when the exporter has written them.
  for (const char* name : {"ring2_K1.txt", "ring2_K2.txt", "ring2_K3.txt"})
  {
    const fs::path path = fs::path(data_dir) / name;
    if (!fs::exists(path)) continue;
    std::printf("\n%s:\n", path.string().c_str());
    try
    {
      const ring_export data = read_ring_export(path.string());
      const xfem_system system = make_system(data, enrichment_rule::per_node);
      std::printf("  %s\n  %s\n", data.settings.c_str(), system.enrichment().summary().c_str());
      const std::string label = data.case_name;
      test_export_format(data, label.c_str());
      test_welded_affine(data, label.c_str());
      test_rigid_rotation(data, label.c_str());
      test_masses(data, label.c_str());
      for (const element_path p : {element_path::plain, element_path::two_piece, element_path::quadrature})
      {
        bool present = false;
        for (const xfem_element& e : system.elements()) present = present || e.path == p;
        if (present) test_finite_differences(system, p, rng, label.c_str());
      }
      test_free_run(data, rng, label.c_str());
      test_no_seams(data, rng, label.c_str());
      const real load = label == "K2" ? 3000 : 400;  // the demo's defaults
      test_mode_one(data, load, label.c_str());
      if (!data.through()) test_front_constraint(data, load, label.c_str());
    }
    catch (const std::exception& error)
    {
      report(false, name, error.what());
    }
  }

  if (!write_dir.empty())
  {
    const fs::path dir(write_dir);
    fs::create_directories(dir);
    write_ring_export((dir / "synthetic_partial.txt").string(), partial);
    write_ring_export((dir / "synthetic_through.txt").string(), through);
    std::printf("\nwrote %s and %s\n", (dir / "synthetic_partial.txt").string().c_str(),
                (dir / "synthetic_through.txt").string().c_str());
  }

  std::printf("\n%d passed, %d failed\n", passes, failures);
  return failures;
}
