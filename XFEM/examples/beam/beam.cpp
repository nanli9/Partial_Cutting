// The beam of cutting's gxfem_beam example, cracked: XFEM with partial_cut's
// enrichment on a real mesh.
//
// The beam is that example's own mesh, 1 x 0.2 x 0.2 along x with its section
// centred on the axis, 85 nodes and 203 tets (data/beam/: made by cutting's
// tetmesh from geometry::make_box with the example's TetGen switches; its
// README.md shows that every statistic the example prints matches). The
// physics is the example's too: every node at x = -0.5 is pinned to the wall
// (u_i, and a_i where it is enriched), gravity pulls the rest down -y at
// 9.81 m/s^2, neo-Hookean E = 1e6 Pa, nu = 0.4, rho = 1000 kg/m^3, dt = 1/60 s,
// backward euler, one Newton solve per frame (the xfem library's: lumped mass,
// sparse LDLT, a line search that keeps det F > 0; the example's newton takes
// the consistent mass). Gravity on an enrichment is rho g int N_a psi_a
// (xfem_system::body_force), not its mass times g. With --pull F, as in the
// example, there is no gravity: F newtons in all pull the free end along
// --pull-direction (+x unless given), shared equally over the nodes of the end
// face that lead the way (all of them along +x, the bottom edge's when pulled
// down), on each node's own dofs. The pull grows from nothing over the first
// second and is held after, and a green arrow on each of those nodes points
// the way.
//
// The beam is uncut until --cut-frame (40, as the example's). Then the case's
// cracks go in: their enrichment comes from partial_cut in this process
// (live_mesh.hpp: xfem_mesh.h's buildXfemMesh, check 4's fast method on the
// mesh, a few seconds, solved once at start-up), and the cracked system takes
// over from the nodes' positions and velocities with every new jump a_i = 0,
// so nothing moves at the cut. With a front, the jump is held at zero on it,
// tet by tet (xfem_system::zero_jump_at). No contact: a crack pushed shut
// passes through itself, so every blade comes from the top, y = +0.1, where
// the sag opens it.
//
// The cases, in rest space (beam x in [-0.5, 0.5], y and z in [-0.1, 0.1]):
//
//   notch            x = -0.245, from the top to 0.7 of the height: the
//                    example's default cut (partial); --depth d, as in the
//                    example, sets its depth
//   through          x = 0.2 all the way down: the example's --depth 1
//                    --cut-x 0.2, the free end falls away
//   curved_front     x = -0.245, 0.7 deep at z = 0 and 0.4 at the sides
//   curved_surface   a cylinder (axis z, radius 0.15) through x = -0.245 at
//                    the top, bending towards +x as it goes down, 0.6 deep
//   multi            notches at x = -0.245 (0.6 deep) and x = +0.15 (0.5 deep)
//   mixed            a notch at x = -0.245 (0.6 deep) and a cut through at
//                    x = +0.245
//
// x = -0.25 and x = +0.25 hold mesh nodes (TetGen put Steiner points there), so
// the example's -0.25 and the mixed cut's +0.25 are moved 0.005 to the side
// (4% of the mean edge, 0.123); start-up prints every crack's smallest |phi|
// at the nodes of the tets it cuts.
//
// Every frame prints a line: Newton, the energies, the tip drop (as the
// example measures it: the mean drop of the nodes at x = +0.5; with --pull,
// how far they have moved along the pull), the crack
// mouth opening per crack (the mean distance between the two sheets' copies of
// the crack's points on the top face, and its part along the crack's normal:
// negative, the faces have passed through each other), the jump on the front
// and the smallest det F. What is drawn, in the style of cutting's gxfem_beam: the body, cut
// open, every vertex carried by the field of its own tet on its own side
// (xfem_system::deformed_point), with both sheets of every crack, in the
// example's beam colour; its tets' edges, as lines carried the same way (once
// cut, each side of a crack apart, and the crack's line across a face only
// where it is cracked); a red ball on every enriched node once the cut is in;
// the pull's arrows; and the example's wall, from the camera of its clips.
//
//   xfem_beam                                    a window, the notch
//   xfem_beam --case multi                       notch | through | curved_front | curved_surface | multi | mixed
//   xfem_beam --no-window --case through --frames 60 --cut-frame 20
//   xfem_beam --no-window --depth 0.9            the notch 0.9 of the way down
//   xfem_beam --case notch --frames 120 --screenshot notch.png
//   xfem_beam --case mixed --record out/mixed    PNG frames, and an mp4 if ffmpeg is on the PATH
//   xfem_beam --csv notch.csv                    one row per frame
//   xfem_beam --quad-level 1                     the export quadrature refined once more (default 0)
//   xfem_beam --pull 100 --pull-direction 0 -1 0 no gravity: 100 N pull the free end's bottom edge down
//   xfem_beam --home-view --record out/notch     the example's own recording camera, to set beside its frames

#include "live_mesh.hpp"

#include "xfem/backward_euler.hpp"
#include "xfem/neo_hookean.hpp"
#include "xfem/tet4.hpp"
#include "xfem/tetgen_mesh.hpp"
#include "xfem/xfem_enrichment.hpp"
#include "xfem/xfem_system.hpp"

#include "imgui.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/screenshot.h"
#include "polyscope/surface_mesh.h"

#include <algorithm>
#include <array>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using namespace xfem;
namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// settings and cases
// ---------------------------------------------------------------------------

struct settings {
  double dt = 1.0 / 60;      // s
  double youngs = 1.0e6;     // Pa
  double poisson = 0.4;
  double density = 1000;     // kg / m^3, cutting's fem_options::density
  double gravity = 9.81;     // m / s^2, along -y
  double pull = 0;           // N on the free end, in place of gravity when not 0
  vec3 pull_direction = vec3::UnitX();
  double pull_ramp = 1;      // s, the pull's growth from nothing, the example's
  double tolerance = 1e-8;   // Newton: |g|_inf below this times E h^2
  int max_iterations = 50;   // the example's newton cap
  int cut_frame = 40;
  int quad_level = 0;
};

// The beam's half extents.
constexpr real half_length = 0.5, half_height = 0.1, half_width = 0.1;

// A blade straight down from the top at x = x0 to `depth` of the height, past
// both sides in z: phi = x - x0, psi = y_front - y. depth >= 1 cuts through:
// psi < 0 everywhere, so no front and no harmonic node.
crack_levels notch(const std::string& name, real x0, real depth) {
  crack_levels c;
  c.name = name;
  c.phi = [x0](const vec3& X) { return X.x() - x0; };
  const real front = half_height - depth * 2 * half_height;
  if (depth >= 1) c.psi = [](const vec3&) { return real(-1); };
  else c.psi = [front](const vec3& X) { return front - X.y(); };
  return c;
}

// The notch with a curved front: depth 0.7 at z = 0 and 0.4 at z = +-0.1,
// y_front(z) = -0.04 + 6 z^2; psi = (y_front - y) / sqrt(1 + y_front'^2), a
// distance to first order.
crack_levels curved_front(const std::string& name, real x0) {
  crack_levels c;
  c.name = name;
  c.phi = [x0](const vec3& X) { return X.x() - x0; };
  c.psi = [](const vec3& X) {
    const real z = X.z();
    const real front = real(-0.04) + 6 * z * z;
    return (front - X.y()) / std::sqrt(1 + 144 * z * z);
  };
  return c;
}

// A cylindrical crack, axis along z through (x0 + R, top): it meets the top at
// x = x0, vertical there, and bends towards +x as it goes down. phi = R - the
// distance to the axis (> 0 on the +x side); the front is the line at `depth`
// of the height, psi the distance along the cylinder's tangent there (< 0 above).
crack_levels curved_surface(const std::string& name, real x0, real radius, real depth) {
  crack_levels c;
  c.name = name;
  const real cx = x0 + radius, cy = half_height;
  c.phi = [=](const vec3& X) { return radius - std::hypot(X.x() - cx, X.y() - cy); };
  const real fy = half_height - depth * 2 * half_height;
  const real fx = cx - std::sqrt(radius * radius - (fy - cy) * (fy - cy));
  // the tangent at the front, away from the crack: the radius (fx - cx, fy - cy) / R turned by +90 degrees
  const real tx = -(fy - cy) / radius, ty = (fx - cx) / radius;
  c.psi = [=](const vec3& X) { return (X.x() - fx) * tx + (X.y() - fy) * ty; };
  return c;
}

struct crack_case {
  std::string name, description;
  std::vector<crack_levels> cracks;
};

// notch_depth is the notch case's depth, --depth (the example's flag; 0.7 by default).
std::vector<crack_case> all_cases(real notch_depth) {
  const real x_left = -0.245, x_right = 0.245;  // -0.25 and +0.25 hold nodes: see the top of the file
  char notch_text[96];
  std::snprintf(notch_text, sizeof notch_text, "x = -0.245, from the top to %g of the height (%s)", notch_depth,
                notch_depth >= 1 ? "complete" : "partial");
  return {
      {"notch", notch_text, {notch("notch", x_left, notch_depth)}},
      {"through", "x = 0.2, all the way down (complete: the end falls away)", {notch("through", 0.2, 1.0)}},
      {"curved_front", "x = -0.245, 0.7 deep at z = 0, 0.4 at the sides", {curved_front("curved_front", x_left)}},
      {"curved_surface", "a cylinder of radius 0.15 through x = -0.245 at the top, 0.6 deep",
       {curved_surface("cylinder", x_left, 0.15, 0.6)}},
      {"multi", "notches at x = -0.245 (0.6 deep) and x = +0.15 (0.5 deep)",
       {notch("left", x_left, 0.6), notch("right", 0.15, 0.5)}},
      {"mixed", "a notch at x = -0.245 (0.6 deep) and a cut through at x = +0.245",
       {notch("notch", x_left, 0.6), notch("through", x_right, 1.0)}},
  };
}

// ---------------------------------------------------------------------------
// the simulation, with no window in it
// ---------------------------------------------------------------------------

struct frame_stats {
  int frame = -1;
  double time = 0;
  bool cut = false;
  newton_result newton;
  double elastic = 0, kinetic = 0;
  double tip_moved = 0;           // along the load: down under gravity, along the pull with --pull
  std::vector<double> mouth;      // per crack: the mean distance between the sheets at the mouth
  std::vector<double> mouth_gap;  // the same along the crack's normal (< 0: the faces pass through, no contact)
  double front_jump = 0;          // the largest jump on the fronts, sum_a N_a a_a there
  double min_det = 1;
  double seconds = 0;
};

// Render vertices of one drawn group, compacted, with its triangles.
struct render_group {
  std::vector<int> ids;
  std::vector<std::array<std::size_t, 3>> faces;
};

// Render vertices joined by lines, compacted.
struct render_lines {
  std::vector<int> ids;
  std::vector<std::array<std::size_t, 2>> edges;
};

// The two sheets' copies of one point of a crack's mouth (its line on the top
// face), with its tet's shape gradients and the crack's rest normal there
// (grad phi of the tet's linear phi).
struct mouth_pair {
  int plus = 0, minus = 0;
  index tet = 0;
  matrix<4, 3> grad_n = matrix<4, 3>::Zero();
  vec3 normal = vec3::Zero();
};

class simulation {
 public:
  tet_mesh mesh;
  settings config;
  crack_case shown;
  mesh_enrichment data;            // the case's enrichment, from partial_cut
  std::unique_ptr<xfem_system> uncut, cracked;
  const xfem_system* system = nullptr;  // the one stepping: uncut until the cut, cracked after
  mechanical_state state;
  frame_stats last;
  std::vector<index> held, tip;
  std::vector<index> pulled;       // with --pull: the end's nodes that lead the way
  render_group groups[3];
  render_group body;               // the three groups as one mesh, as drawn
  render_lines edges_whole, edges_cut;  // the tets' edges before the cut and after
  std::vector<std::vector<mouth_pair>> mouth_pairs;  // per crack
  double force_scale = 1;          // E h^2
  double mean_edge = 0;
  index front_points = 0;
  int frame = 0;
  double time = 0;
  bool cut = false;
  int failed_solves = 0;

  explicit simulation(tet_mesh m) : mesh(std::move(m)) {
    for (index v = 0; v < mesh.vertex_count(); ++v)
    {
      const real x = mesh.vertex(v).x();
      if (x < -half_length + 1e-6) held.push_back(v);
      if (x > half_length - 1e-6) tip.push_back(v);
    }
    double edges = 0;
    int count = 0;
    for (const auto& t : mesh.tets)
    {
      for (int a = 0; a < 4; ++a)
      {
        for (int b = a + 1; b < 4; ++b, ++count) edges += (mesh.vertex(t[a]) - mesh.vertex(t[b])).norm();
      }
    }
    mean_edge = edges / count;
  }

  // The case's cracks: the enrichment from partial_cut, the render groups and
  // the mouth pairs. Then a reset.
  void load(const crack_case& c) {
    shown = c;
    std::printf("case %s: %s\nsolving the enrichment with partial_cut (quadrature level %d)...\n", c.name.c_str(),
                c.description.c_str(), config.quad_level);
    data = live_mesh(mesh, c.cracks, config.quad_level);
    for (std::size_t k = 0; k < c.cracks.size(); ++k)
    {
      std::printf("crack %s: smallest |phi| at the nodes of the tets it cuts %.4f = %.1f%% of the mean edge %.4f%s\n",
                  c.cracks[k].name.c_str(), data.min_phi[k], 100 * data.min_phi[k] / mean_edge, mean_edge,
                  data.min_phi[k] < 0.02 * mean_edge ? "  (under 2%: move the crack)" : "");
    }
    build_render();
    rebuild();
  }

  // The two systems from the current settings, both pinned at the wall; the
  // run back at rest, uncut, frame 0.
  void rebuild() {
    const neo_hookean material = neo_hookean::from(config.youngs, config.poisson);
    {
      enrichment_options options;
      options.rule = enrichment_rule::none;
      xfem_enrichment plain = build_enrichment(mesh, std::vector<real>(mesh.vertices.size(), real(1)),
                                               std::vector<int>(mesh.tets.size(), 0),
                                               std::vector<real>(mesh.vertices.size(), real(0)), options);
      uncut = std::make_unique<xfem_system>(mesh, std::move(plain), std::vector<quadrature_point>{}, material, config.density);
    }
    xfem_enrichment enrichment = build_enrichment(mesh, data.node_phi, data.tet_class, data.node_value);
    cracked = std::make_unique<xfem_system>(mesh, std::move(enrichment), data.quadrature, material, config.density);
    for (xfem_system* s : {uncut.get(), cracked.get()})
    {
      for (const index v : held)
      {
        s->fix_node(v);
        s->fix_node_enrichment(v);
      }
    }
    front_points = 0;
    for (const mesh_enrichment::front_piece& f : data.front)
    {
      cracked->zero_jump_at(f.tet, f.a);
      cracked->zero_jump_at(f.tet, f.b);
      front_points += 2;
    }
    force_scale = config.youngs * mean_edge * mean_edge;
    // The pulled nodes: of the end's, the ones furthest the way it is pulled,
    // so no arrow crosses the beam.
    pulled.clear();
    real leading = -DBL_MAX;
    for (const index v : tip) leading = std::max(leading, config.pull_direction.dot(mesh.vertex(v)));
    for (const index v : tip)
    {
      if (config.pull_direction.dot(mesh.vertex(v)) > leading - 1e-6) pulled.push_back(v);
    }
    reset();
  }

  vec3 load_direction() const { return config.pull != 0 ? config.pull_direction : vec3(-vec3::UnitY()); }

  // The load over the step that ends at `t`: gravity, or with --pull a share of
  // the pull on each pulled node's own dofs, all of it after the ramp.
  vector_x load(double t) const {
    if (config.pull == 0) return system->body_force(vec3(0, -config.gravity, 0));
    vector_x f = vector_x::Zero(system->dofs());
    const double share = config.pull * std::min(1.0, t / config.pull_ramp) / double(pulled.size());
    for (const index v : pulled) f.segment<3>(3 * v) += share * config.pull_direction;
    return f;
  }

  void reset() {
    system = uncut.get();
    state.x = system->rest_positions();
    state.v = vector_x::Zero(system->dofs());
    frame = 0;
    time = 0;
    cut = false;
    failed_solves = 0;
    last = frame_stats{};
    measure(last);
  }

  // The cracked system takes over: the nodes keep their positions and
  // velocities, every jump starts at 0, so the field and its energy are the
  // uncut one's.
  void cut_now() {
    if (cut) return;
    const index nodes = 3 * mesh.vertex_count();
    mechanical_state next;
    next.x = cracked->rest_positions();
    next.v = vector_x::Zero(cracked->dofs());
    next.x.head(nodes) = state.x.head(nodes);
    next.v.head(nodes) = state.v.head(nodes);
    const double before = system->potential_energy(state.x);
    system = cracked.get();
    state = next;
    cut = true;
    std::printf("frame %d: the cut goes in; elastic energy %.6e before, %.6e after\n", frame, before,
                system->potential_energy(state.x));
  }

  void step() {
    if (!cut && frame == config.cut_frame) cut_now();
    newton_options options;
    options.gradient_tolerance = config.tolerance * force_scale;
    options.max_iterations = config.max_iterations;
    backward_euler integrator(options);
    const auto start = std::chrono::steady_clock::now();
    frame_stats s;
    // Backward euler takes the load at the end of the step.
    s.newton = integrator.step(*system, state, config.dt, load(time + config.dt));
    s.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    time += config.dt;
    s.frame = frame;
    s.time = time;
    s.cut = cut;
    if (!s.newton.converged) ++failed_solves;
    measure(s);
    last = s;
    ++frame;
  }

  void measure(frame_stats& s) const {
    s.elastic = system->potential_energy(state.x);
    s.kinetic = 0.5 * state.v.dot(system->mass().cwiseProduct(state.v));
    s.min_det = system->minimum_determinant(state.x);
    s.tip_moved = 0;
    for (const index v : tip)
    {
      s.tip_moved += load_direction().dot(state.x.segment<3>(3 * v) - mesh.vertex(v)) / double(tip.size());
    }
    // The mouth: the sheets' separation, and its part along the crack's normal
    // carried by the nodes' deformation F = sum_a x_a gradN_a^T of the tet,
    // n = F^-T N normalised.
    s.mouth.assign(mouth_pairs.size(), 0.0);
    s.mouth_gap.assign(mouth_pairs.size(), 0.0);
    for (std::size_t k = 0; k < mouth_pairs.size(); ++k)
    {
      for (const mouth_pair& m : mouth_pairs[k])
      {
        const vec3 apart = deformed(m.plus) - deformed(m.minus);
        mat3 F = mat3::Zero();
        for (int a = 0; a < 4; ++a) F += state.x.segment<3>(3 * mesh.corner(m.tet, a)) * m.grad_n.row(a);
        s.mouth[k] += apart.norm();
        s.mouth_gap[k] += apart.dot((F.inverse().transpose() * m.normal).normalized());
      }
      if (mouth_pairs[k].empty()) continue;
      s.mouth[k] /= double(mouth_pairs[k].size());
      s.mouth_gap[k] /= double(mouth_pairs[k].size());
    }
    // The jump on the fronts: the + face (h = 1) less the - face (h = 0) at
    // each end of every front piece, which is sum_a N_a a_a there.
    s.front_jump = 0;
    if (!cut) return;
    for (const mesh_enrichment::front_piece& f : data.front)
    {
      for (const vec3& X : {f.a, f.b})
      {
        const vec3 j = system->deformed_point(state.x, f.tet, X, 1, 1.0) - system->deformed_point(state.x, f.tet, X, -1, 0.0);
        s.front_jump = std::max(s.front_jump, j.norm());
      }
    }
  }

  vec3 deformed(int render_id) const {
    const render_vertex& v = data.render_vertices[static_cast<std::size_t>(render_id)];
    return system->deformed_point(state.x, v.tet, v.X, v.side, v.h);
  }

  std::vector<std::array<double, 3>> positions(const std::vector<int>& ids) const {
    std::vector<std::array<double, 3>> out(ids.size());
    for (std::size_t i = 0; i < ids.size(); ++i)
    {
      const vec3 p = deformed(ids[i]);
      out[i] = {p[0], p[1], p[2]};
    }
    return out;
  }

  std::vector<vec3> pulled_positions() const {
    std::vector<vec3> out;
    for (const index v : pulled) out.push_back(state.x.segment<3>(3 * v));
    return out;
  }

  // The enriched nodes, where their own material is: the enrichment is shifted
  // (psi_a = 0 at X_a), so that is x_a, the node's own dofs.
  std::vector<std::array<double, 3>> enriched_positions() const {
    std::vector<std::array<double, 3>> out;
    for (index v = 0; v < mesh.vertex_count(); ++v)
    {
      if (cracked->enrichment().enriched(v)) out.push_back({state.x[3 * v], state.x[3 * v + 1], state.x[3 * v + 2]});
    }
    return out;
  }

  std::string describe() const {
    const xfem_enrichment& en = cracked->enrichment();
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "%s: %d dofs once cut (%d nodes, %d enrichments), %d integration points, %zu quadrature points from "
                  "partial_cut; the jump held at 0 at %d front points (rank %d); %zu nodes pinned at x = -0.5\n",
                  shown.name.c_str(), int(cracked->dofs()), int(en.vertex_count()), int(en.enrichment_count()),
                  int(cracked->point_count()), data.quadrature.size(), int(front_points), int(cracked->constraint_rank()),
                  held.size());
    std::string out = buf;
    if (config.pull != 0)
    {
      std::snprintf(buf, sizeof buf, "pull: %g N along (%g, %g, %g) on %zu nodes of the end at x = +0.5:", config.pull,
                    config.pull_direction.x(), config.pull_direction.y(), config.pull_direction.z(), pulled.size());
      out += buf;
      for (const index v : pulled) out += " " + std::to_string(v);
      std::snprintf(buf, sizeof buf, ", grown over the first %g s; no gravity\n", config.pull_ramp);
      out += buf;
    }
    return out + en.summary();
  }

 private:
  void build_render() {
    for (render_group& g : groups) g = render_group{};
    std::vector<std::array<int, 3>> local(data.render_vertices.size(), {-1, -1, -1});
    for (const render_triangle& t : data.render_triangles)
    {
      render_group& g = groups[t.group];
      std::array<std::size_t, 3> face{};
      for (int k = 0; k < 3; ++k)
      {
        int& id = local[static_cast<std::size_t>(t.v[k])][static_cast<std::size_t>(t.group)];
        if (id < 0)
        {
          id = static_cast<int>(g.ids.size());
          g.ids.push_back(t.v[k]);
        }
        face[static_cast<std::size_t>(k)] = static_cast<std::size_t>(id);
      }
      g.faces.push_back(face);
    }
    body = render_group{};
    for (const render_group& g : groups)
    {
      const std::size_t offset = body.ids.size();
      body.ids.insert(body.ids.end(), g.ids.begin(), g.ids.end());
      for (const auto& f : g.faces) body.faces.push_back({offset + f[0], offset + f[1], offset + f[2]});
    }
    edges_whole = outline(false);
    edges_cut = outline(true);

    // Pair the two copies of every crack-sheet vertex on the top face: same
    // tet, same X; per crack, from the tet's.
    mouth_pairs.assign(shown.cracks.size(), {});
    std::map<std::tuple<int, long long, long long, long long>, int> minus_copy;
    const auto key = [](const render_vertex& v) {
      const auto q = [](double x) { return static_cast<long long>(std::llround(x * 1e9)); };
      return std::make_tuple(v.tet, q(v.X[0]), q(v.X[1]), q(v.X[2]));
    };
    for (const int id : groups[2].ids) minus_copy[key(data.render_vertices[static_cast<std::size_t>(id)])] = id;
    for (const int id : groups[1].ids)
    {
      const render_vertex& v = data.render_vertices[static_cast<std::size_t>(id)];
      if (v.X.y() < half_height - 1e-9) continue;
      const auto found = minus_copy.find(key(v));
      const int crack = data.tet_crack[static_cast<std::size_t>(v.tet)];
      if (found == minus_copy.end() || crack < 0) continue;
      mouth_pair m;
      m.plus = id;
      m.minus = found->second;
      m.tet = v.tet;
      vec3 X[4];
      mesh.corners(v.tet, X);
      m.grad_n = tet4::shape_gradients(tet4::build_rest_state(X));
      for (int a = 0; a < 4; ++a) m.normal += data.node_phi[static_cast<std::size_t>(mesh.corner(v.tet, a))] * m.grad_n.row(a).transpose();
      m.normal.normalize();
      mouth_pairs[static_cast<std::size_t>(crack)].push_back(m);
    }
  }

  // The tets' edges, as cutting's beam draws them: the outline of every tet's
  // piece of the drawn surface -- each boundary face, each side of a crack
  // apart once `cut`, and the tet's piece of each crack sheet -- the edges only
  // one triangle of the piece uses. Uncut, a face's two sides are one piece, so
  // the crack's line across it is inside it, not an outline.
  render_lines outline(bool cut) const {
    using point_key = std::tuple<long long, long long, long long>;
    const auto key_of = [](const vec3& X) {
      const auto q = [](double x) { return static_cast<long long>(std::llround(x * 1e9)); };
      return point_key{q(X[0]), q(X[1]), q(X[2])};
    };
    const auto vertex = [this](int id) -> const render_vertex& {
      return data.render_vertices[static_cast<std::size_t>(id)];
    };
    // The face of its tet a boundary triangle lies on, by the corner opposite it.
    const auto face_of = [&](const render_triangle& t) {
      vec3 X[4];
      mesh.corners(vertex(t.v[0]).tet, X);
      for (int k = 0; k < 4; ++k)
      {
        const vec3 &a = X[(k + 1) % 4], &b = X[(k + 2) % 4], &c = X[(k + 3) % 4];
        const vec3 n = (b - a).cross(c - a).normalized();
        bool on = true;
        for (int i = 0; i < 3; ++i) on = on && std::abs(n.dot(vertex(t.v[i]).X - a)) < 1e-9;
        if (on) return k;
      }
      return -1;
    };

    std::map<std::array<int, 4>, int> piece_of;    // (tet, group, face, side) -> piece
    std::vector<std::array<int, 4>> pieces;
    std::vector<std::map<point_key, int>> first;   // per piece: its first render vertex at each point
    std::map<std::tuple<int, int, int>, int> uses;  // (piece, lower vertex, higher) -> triangles
    for (const render_triangle& t : data.render_triangles)
    {
      if (!cut && t.group != 0) continue;
      const render_vertex& v0 = vertex(t.v[0]);
      const std::array<int, 4> key = {v0.tet, t.group, t.group == 0 ? face_of(t) : -1, cut ? v0.side : 0};
      const auto [found, added] = piece_of.emplace(key, static_cast<int>(pieces.size()));
      if (added)
      {
        pieces.push_back(key);
        first.emplace_back();
      }
      const int piece = found->second;
      int c[3];
      for (int k = 0; k < 3; ++k)
      {
        c[k] = first[static_cast<std::size_t>(piece)].emplace(key_of(vertex(t.v[k]).X), t.v[k]).first->second;
      }
      for (int k = 0; k < 3; ++k) ++uses[{piece, std::min(c[k], c[(k + 1) % 3]), std::max(c[k], c[(k + 1) % 3])}];
    }

    // A front tet's faces are split by its crack's whole plane, past the front
    // too, where the two sides are one material: a line on the plane (its ends
    // on the other side's piece as well) is drawn only where psi < 0.
    const auto cracked = [&](int piece, int a, int b) {
      const std::array<int, 4>& k = pieces[static_cast<std::size_t>(piece)];
      if (!cut || k[1] != 0 || data.tet_class[static_cast<std::size_t>(k[0])] != 2) return true;
      const auto other = piece_of.find({k[0], 0, k[2], -k[3]});
      if (other == piece_of.end()) return true;
      const std::map<point_key, int>& there = first[static_cast<std::size_t>(other->second)];
      const vec3 &A = vertex(a).X, &B = vertex(b).X;
      if (!there.count(key_of(A)) || !there.count(key_of(B))) return true;
      return psi_at(k[0], 0.5 * (A + B)) < 0;
    };

    render_lines out;
    std::map<int, std::size_t> compact;
    const auto node = [&](int id) {
      const auto [found, added] = compact.emplace(id, out.ids.size());
      if (added) out.ids.push_back(id);
      return found->second;
    };
    for (const auto& [edge, count] : uses)
    {
      const auto [piece, a, b] = edge;
      if (count == 1 && cracked(piece, a, b)) out.edges.push_back({node(a), node(b)});
    }
    return out;
  }

  // psi of the tet's crack at X, linear in the tet as the enrichment has it.
  real psi_at(int tet, const vec3& X) const {
    vec3 C[4];
    mesh.corners(tet, C);
    mat3 D;
    D << C[1] - C[0], C[2] - C[0], C[3] - C[0];
    const vec3 xi = D.inverse() * (X - C[0]);
    const crack_levels& crack = shown.cracks[static_cast<std::size_t>(data.tet_crack[static_cast<std::size_t>(tet)])];
    return (1 - xi.sum()) * crack.psi(C[0]) + xi[0] * crack.psi(C[1]) + xi[1] * crack.psi(C[2]) + xi[2] * crack.psi(C[3]);
  }
};

void print_frame(const simulation& sim) {
  const frame_stats& s = sim.last;
  std::string mouth;
  for (std::size_t k = 0; k < s.mouth.size(); ++k)
  {
    char buf[128];
    std::snprintf(buf, sizeof buf, "  mouth %s %.5f (gap %+.5f)", sim.shown.cracks[k].name.c_str(), s.mouth[k], s.mouth_gap[k]);
    mouth += buf;
  }
  std::printf("frame %3d  t %6.3f%s  newton %2d  |g| %.2e  elastic %.4e  kinetic %.4e  tip %s %.4f%s  front jump %.1e  "
              "min det F %.3f  %.0f ms%s%s\n",
              s.frame, s.time, s.cut ? "  cut" : "     ", int(s.newton.iterations), s.newton.final_gradient_norm, s.elastic,
              s.kinetic, sim.config.pull != 0 ? "pulled" : "drop", s.tip_moved, mouth.c_str(), s.front_jump, s.min_det, 1e3 * s.seconds,
              s.newton.converged ? "" : "  NOT CONVERGED: ", s.newton.converged ? "" : s.newton.failure.c_str());
  std::fflush(stdout);
}

// One row per frame.
class csv_log {
 public:
  void open(const std::string& path, const simulation& sim) {
    if (path.empty()) return;
    const fs::path p(path);
    if (p.has_parent_path()) fs::create_directories(p.parent_path());
    out_.open(p);
    if (!out_) throw std::runtime_error("cannot write " + path);
    out_ << "frame,time,cut,newton_iterations,gradient_norm,converged,elastic,kinetic,"
         << (sim.config.pull != 0 ? "tip_pulled" : "tip_drop");
    for (const crack_levels& c : sim.shown.cracks) out_ << ",mouth_" << c.name << ",mouth_gap_" << c.name;
    out_ << ",front_jump,min_det,solve_ms\n";
  }
  void row(const simulation& sim) {
    if (!out_.is_open()) return;
    const frame_stats& s = sim.last;
    char buf[512];
    std::snprintf(buf, sizeof buf, "%d,%.6f,%d,%d,%.6e,%d,%.9e,%.9e,%.9e", s.frame, s.time, s.cut ? 1 : 0,
                  int(s.newton.iterations), s.newton.final_gradient_norm, s.newton.converged ? 1 : 0, s.elastic, s.kinetic,
                  s.tip_moved);
    out_ << buf;
    for (std::size_t k = 0; k < s.mouth.size(); ++k)
    {
      std::snprintf(buf, sizeof buf, ",%.9e,%.9e", s.mouth[k], s.mouth_gap[k]);
      out_ << buf;
    }
    std::snprintf(buf, sizeof buf, ",%.6e,%.9f,%.3f\n", s.front_jump, s.min_det, 1e3 * s.seconds);
    out_ << buf;
    out_.flush();
  }

 private:
  std::ofstream out_;
};

// ---------------------------------------------------------------------------
// the window
// ---------------------------------------------------------------------------

glm::vec3 g3(const vec3& v) { return glm::vec3(float(v[0]), float(v[1]), float(v[2])); }

// The example's wall: a plate 0.04 of the length thick whose face is the plane
// x = -0.5, reaching a quarter of the length past the beam on every side, and
// twice that below, where the beam sags.
void register_wall() {
  const double length = 2 * half_length, reach = 0.25 * length;
  const double x0 = -half_length - 0.04 * length, x1 = -half_length;
  const double y0 = -half_height - 2 * reach, y1 = half_height + reach, z = half_width + reach;
  std::vector<std::array<double, 3>> v;
  for (int k = 0; k < 8; ++k) v.push_back({k & 1 ? x1 : x0, k & 2 ? y1 : y0, k & 4 ? z : -z});
  const std::vector<std::array<std::size_t, 3>> f = {{0, 2, 3}, {0, 3, 1}, {4, 5, 7}, {4, 7, 6}, {0, 1, 5}, {0, 5, 4},
                                                     {2, 6, 7}, {2, 7, 3}, {0, 4, 6}, {0, 6, 2}, {1, 3, 7}, {1, 7, 5}};
  polyscope::SurfaceMesh* wall = polyscope::registerSurfaceMesh("wall", v, f);
  wall->setSurfaceColor({0.55f, 0.57f, 0.62f});
}

struct triangles {
  std::vector<std::array<double, 3>> points;
  std::vector<std::array<std::size_t, 3>> faces;
};

// The example's arrow: from the origin along the unit `direction`, `length`
// long, a shaft of `radius` and a cone three times as wide and a quarter of the
// length at its tip, closed at both ends.
triangles arrow(const vec3& direction, real length, real radius) {
  constexpr int sides = 16;
  const vec3 across = direction.unitOrthogonal();
  const vec3 other = direction.cross(across);
  const real pi = std::acos(real(-1));
  triangles out;
  const auto add = [&](const vec3& p) {
    out.points.push_back({p[0], p[1], p[2]});
    return out.points.size() - 1;
  };
  const auto circle = [&](real along, real r) {
    const std::size_t first = out.points.size();
    for (int k = 0; k < sides; ++k)
    {
      const real angle = 2 * pi * k / sides;
      add(along * direction + r * (std::cos(angle) * across + std::sin(angle) * other));
    }
    return first;
  };
  const real head = length / 4;
  const std::size_t base = circle(0, radius), neck = circle(length - head, radius), rim = circle(length - head, 3 * radius);
  const std::size_t tail = add(vec3::Zero()), tip = add(length * direction);
  for (std::size_t k = 0; k < sides; ++k)
  {
    const std::size_t next = (k + 1) % sides;
    out.faces.push_back({tail, base + next, base + k});  // the shaft's end
    out.faces.push_back({base + k, base + next, neck + next});  // the shaft
    out.faces.push_back({base + k, neck + next, neck + k});
    out.faces.push_back({neck + k, neck + next, rim + next});  // the cone's back
    out.faces.push_back({neck + k, rim + next, rim + k});
    out.faces.push_back({rim + k, rim + next, tip});  // the cone
  }
  return out;
}

// `shape` once at each of `at`, as one mesh.
triangles copies(const triangles& shape, const std::vector<vec3>& at) {
  triangles out;
  for (const vec3& p : at)
  {
    const std::size_t offset = out.points.size();
    for (const auto& q : shape.points) out.points.push_back({q[0] + p[0], q[1] + p[1], q[2] + p[2]});
    for (const auto& f : shape.faces) out.faces.push_back({offset + f[0], offset + f[1], offset + f[2]});
  }
  return out;
}

struct viewer {
  simulation& sim;
  std::vector<crack_case> cases;
  int case_index = 0;
  bool playing = false;
  bool edges_cut = false;  // which of the tets' edges are registered
  bool home_view = false;  // --home-view

  void register_all() {
    polyscope::removeAllStructures();
    register_wall();
    // The example's beam colour, the crack sheets too, as its cut faces are.
    polyscope::SurfaceMesh* body = polyscope::registerSurfaceMesh("beam", sim.positions(sim.body.ids), sim.body.faces);
    body->setSurfaceColor({0.85f, 0.80f, 0.72f});
    // Its red balls: 0.035 of the section across.
    polyscope::PointCloud* balls = polyscope::registerPointCloud("enriched nodes", sim.enriched_positions());
    balls->setPointColor({1.f, 0.f, 0.f});
    balls->setPointRadius(0.035 * 2 * half_height, false);
    // Its arrows, 0.6 of the section long, on the pulled nodes.
    if (sim.config.pull != 0)
    {
      const triangles drawn = arrows();
      polyscope::registerSurfaceMesh("pull", drawn.points, drawn.faces)->setSurfaceColor({0.10f, 0.60f, 0.25f});
    }
    register_edges();
    update();
  }

  triangles arrows() const {
    const real section = 2 * std::min(half_height, half_width);
    return copies(arrow(sim.config.pull_direction, 0.6 * section, 0.02 * section), sim.pulled_positions());
  }

  // The lines in Polyscope's edge colour, which the example's beam is drawn
  // with, about as wide as its edges at this camera.
  void register_edges() {
    edges_cut = sim.cut;
    const render_lines& lines = sim.cut ? sim.edges_cut : sim.edges_whole;
    polyscope::CurveNetwork* edges = polyscope::registerCurveNetwork("tet edges", sim.positions(lines.ids), lines.edges);
    edges->setColor({0.f, 0.f, 0.f});
    edges->setRadius(0.0011f, false);
  }

  void update() {
    if (sim.cut != edges_cut) register_edges();
    polyscope::getSurfaceMesh("beam")->updateVertexPositions(sim.positions(sim.body.ids));
    const render_lines& lines = sim.cut ? sim.edges_cut : sim.edges_whole;
    polyscope::getCurveNetwork("tet edges")->updateNodePositions(sim.positions(lines.ids));
    polyscope::PointCloud* balls = polyscope::getPointCloud("enriched nodes");
    balls->updatePointPositions(sim.enriched_positions());
    balls->setEnabled(sim.cut);
    if (sim.config.pull != 0) polyscope::getSurfaceMesh("pull")->updateVertexPositions(arrows().points);
  }

  // The camera of the example's clips; with --home-view, Polyscope's home view
  // of the same scene extents (the wall's and the beam's), which is what the
  // example itself records with, for frames that line up with its own.
  void camera() {
    if (home_view) polyscope::view::resetCameraToHomeView();
    else polyscope::view::lookAt(g3(vec3(1.0, 0.45, 1.05)), g3(vec3(0.1, -0.15, 0.0)), g3(vec3::UnitY()));
  }

  void load_case(int index) {
    case_index = index;
    sim.load(cases[static_cast<std::size_t>(index)]);
    std::printf("%s\n", sim.describe().c_str());
    register_all();
  }

  void ui() {
    ImGui::PushItemWidth(170);
    ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + 460);
    if (ImGui::BeginCombo("case", cases[static_cast<std::size_t>(case_index)].name.c_str()))
    {
      for (int i = 0; i < int(cases.size()); ++i)
      {
        if (ImGui::Selectable(cases[static_cast<std::size_t>(i)].name.c_str(), i == case_index) && i != case_index)
        {
          load_case(i);
        }
      }
      ImGui::EndCombo();
    }
    ImGui::TextWrapped("%s", sim.shown.description.c_str());
    if (ImGui::Button(playing ? "pause" : "play")) playing = !playing;
    ImGui::SameLine();
    if (ImGui::Button("step")) advance();
    ImGui::SameLine();
    if (ImGui::Button("reset"))
    {
      sim.reset();
      update();
    }
    ImGui::SameLine();
    if (!sim.cut && ImGui::Button("cut now")) sim.cut_now();
    ImGui::SameLine();
    if (ImGui::Button("view")) camera();
    ImGui::InputInt("cut frame", &sim.config.cut_frame);
    ImGui::PopItemWidth();

    ImGui::Separator();
    const frame_stats& s = sim.last;
    if (ImGui::TreeNode("mesh, enrichment and partial_cut"))
    {
      ImGui::TextUnformatted(sim.describe().c_str());
      ImGui::TextUnformatted(sim.data.settings.c_str());
      ImGui::TreePop();
    }
    ImGui::Text("frame %d  t = %.3f s  %s", s.frame, s.time, sim.cut ? "cut" : "uncut (the cut at frame cut frame)");
    ImGui::Text("Newton: %d iterations, |g| = %.2e (tol %.1e)%s", int(s.newton.iterations), s.newton.final_gradient_norm,
                sim.config.tolerance * sim.force_scale, s.frame < 0 || s.newton.converged ? "" : "  NOT CONVERGED");
    ImGui::Text("energy: elastic %.4e, kinetic %.4e J", s.elastic, s.kinetic);
    ImGui::Text("tip %s %.4f m", sim.config.pull != 0 ? "pulled" : "drop", s.tip_moved);
    for (std::size_t k = 0; k < s.mouth.size(); ++k)
    {
      ImGui::Text("crack mouth opening, %s: %.5f m (along its normal %+.5f%s)", sim.shown.cracks[k].name.c_str(), s.mouth[k],
                  s.mouth_gap[k], s.mouth_gap[k] < -1e-4 ? ": faces pushed through, no contact" : "");
    }
    if (!sim.data.front.empty()) ImGui::Text("jump on the fronts: %.1e (held at 0)", s.front_jump);
    ImGui::Text("min det F %.4f   step %.0f ms", s.min_det, 1e3 * s.seconds);
    ImGui::PopTextWrapPos();

    if (playing) advance();
  }

  void advance() {
    sim.step();
    print_frame(sim);
    update();
  }
};

// ---------------------------------------------------------------------------

int run(int argc, char** argv) {
  settings config;
  std::string case_name = "notch", screenshot, record, csv;
  bool window = true, home_view = false;
  int frames = 120;
  real notch_depth = 0.7;

  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    const auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument(arg + " needs a value");
      return argv[++i];
    };
    if (arg == "--no-window") window = false;
    else if (arg == "--case") case_name = next();
    else if (arg == "--frames") frames = std::stoi(next());
    else if (arg == "--cut-frame") config.cut_frame = std::stoi(next());
    else if (arg == "--depth") notch_depth = std::stod(next());
    else if (arg == "--quad-level") config.quad_level = std::stoi(next());
    else if (arg == "--screenshot") screenshot = next();
    else if (arg == "--record") record = next();
    else if (arg == "--csv") csv = next();
    else if (arg == "--home-view") home_view = true;
    else if (arg == "--dt") config.dt = std::stod(next());
    else if (arg == "--tol") config.tolerance = std::stod(next());
    else if (arg == "--pull") config.pull = std::stod(next());
    else if (arg == "--pull-direction")
    {
      for (int k = 0; k < 3; ++k) config.pull_direction[k] = std::stod(next());
      if (!(config.pull_direction.norm() > 0)) throw std::invalid_argument("--pull-direction must not be zero");
      config.pull_direction.normalize();
    }
    else
    {
      std::fprintf(stderr,
                   "usage: %s [--case notch|through|curved_front|curved_surface|multi|mixed] [--frames N] [--cut-frame N]\n"
                   "          [--depth d] [--quad-level 0..4] [--no-window] [--screenshot FILE] [--record DIR] [--csv FILE]\n"
                   "          [--dt s] [--tol t] [--pull F] [--pull-direction x y z] [--home-view]\n",
                   argv[0]);
      return 2;
    }
  }
  if (frames < 1 || config.dt <= 0 || notch_depth <= 0)
  {
    throw std::invalid_argument("--frames, --dt and --depth must be positive");
  }
  if (config.quad_level < 0 || config.quad_level > 4) throw std::invalid_argument("--quad-level: 0 to 4");

  const std::vector<crack_case> cases = all_cases(notch_depth);
  int case_index = -1;
  for (int i = 0; i < int(cases.size()); ++i)
  {
    if (cases[static_cast<std::size_t>(i)].name == case_name) case_index = i;
  }
  if (case_index < 0)
  {
    std::fprintf(stderr, "no case '%s'; the cases:", case_name.c_str());
    for (const crack_case& c : cases) std::fprintf(stderr, " %s", c.name.c_str());
    std::fprintf(stderr, "\n");
    return 2;
  }

  const fs::path mesh_prefix = fs::path(XFEM_DATA_DIR) / "beam" / "beam";
  tet_mesh mesh = read_tetgen_mesh(mesh_prefix.string());
  double volume = 0;
  for (index e = 0; e < mesh.tet_count(); ++e)
  {
    vec3 X[4];
    mesh.corners(e, X);
    volume += tet_volume(X[0], X[1], X[2], X[3]);
  }
  simulation sim(std::move(mesh));
  sim.config = config;
  std::printf("beam 1 x 0.2 x 0.2 (cutting's gxfem_beam mesh, data/beam/): %d nodes, %d tets, volume %g; %zu nodes held "
              "at x = -0.5, %zu at the tip x = +0.5\n",
              int(sim.mesh.vertex_count()), int(sim.mesh.tet_count()), volume, sim.held.size(), sim.tip.size());
  char load_text[64];
  if (config.pull != 0) std::snprintf(load_text, sizeof load_text, "no gravity: the pull");
  else std::snprintf(load_text, sizeof load_text, "gravity %.2f m/s^2 along -y", config.gravity);
  std::printf("E %.4g Pa, nu %.3f, rho %.1f kg/m^3, %s, dt %.4f s, backward euler, Newton (lumped "
              "mass) to |g|_inf < %.2e N, the cut at frame %d\n",
              config.youngs, config.poisson, config.density, load_text, config.dt,
              config.tolerance * config.youngs * sim.mean_edge * sim.mean_edge, config.cut_frame);
  sim.load(cases[static_cast<std::size_t>(case_index)]);
  std::printf("%s\n", sim.describe().c_str());
  std::printf("partial_cut: %s\n", sim.data.settings.c_str());

  csv_log log;
  log.open(csv, sim);
  const auto finish = [&]() {
    std::printf("%d of %d frames did not converge\n", sim.failed_solves, sim.frame);
    return sim.failed_solves > 0 ? 1 : 0;
  };

  if (!window && screenshot.empty() && record.empty())
  {
    for (int f = 0; f < frames; ++f)
    {
      sim.step();
      print_frame(sim);
      log.row(sim);
    }
    return finish();
  }

  polyscope::options::programName = "xfem_beam";
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::options::ssaaFactor = 2;  // the example's renderer's anti-aliasing
  if (!screenshot.empty() || !record.empty()) polyscope::options::usePrefsFile = false;
  polyscope::init();
  polyscope::view::setNavigateStyle(polyscope::NavigateStyle::Turntable);
  if (!screenshot.empty() || !record.empty()) polyscope::view::setWindowSize(1280, 720);

  viewer view{sim, cases};
  view.case_index = case_index;
  view.home_view = home_view;
  view.register_all();
  view.camera();

  if (!screenshot.empty() || !record.empty())
  {
    if (!record.empty()) fs::create_directories(record);
    polyscope::ScreenshotOptions options;
    options.transparentBackground = false;
    for (int f = 0; f < frames; ++f)
    {
      sim.step();
      print_frame(sim);
      log.row(sim);
      if (record.empty()) continue;
      view.update();
      char name[32];
      std::snprintf(name, sizeof name, "frame_%04d.png", f);
      polyscope::screenshot((fs::path(record) / name).string(), options);
    }
    if (!screenshot.empty())
    {
      view.update();
      polyscope::screenshot(screenshot, options);
      std::printf("screenshot written: %s\n", screenshot.c_str());
    }
    if (!record.empty())
    {
      // The movie beside the frames, when ffmpeg is on the PATH; the frames stay either way.
      const std::string frames_in = (fs::path(record) / "frame_%04d.png").string();
      const std::string movie = (fs::path(record) / "xfem_beam.mp4").string();
      std::printf("wrote %d frames to %s\n", frames, frames_in.c_str());
#ifdef _WIN32
      const bool have_ffmpeg = std::system("ffmpeg -version > nul 2>&1") == 0;
#else
      const bool have_ffmpeg = std::system("ffmpeg -version > /dev/null 2>&1") == 0;
#endif
      const std::string command = "ffmpeg -y -loglevel error -framerate 60 -i \"" + frames_in +
                                  "\" -pix_fmt yuv420p \"" + movie + "\"";
      if (!have_ffmpeg) std::printf("no ffmpeg on the PATH: no mp4. By hand:\n  %s\n", command.c_str());
      else if (std::system(command.c_str()) != 0) std::fprintf(stderr, "ffmpeg failed; the frames are kept. By hand:\n  %s\n", command.c_str());
      else std::printf("wrote %s\n", movie.c_str());
    }
    return finish();
  }

  polyscope::state::userCallback = [&view] { view.ui(); };
  polyscope::show();
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try
  {
    return run(argc, argv);
  }
  catch (const std::exception& error)
  {
    std::fprintf(stderr, "error: %s\n", error.what());
    return 1;
  }
}
