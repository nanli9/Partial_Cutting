// One nodal ring, cracked, pulled open: XFEM with partial_cut's enrichment.
//
// The ring is check 4's (partial_cut/README.md, "Check 4"): the octahedron
// |x| + |y| + |z| <= 1 split into the 8 tets around the node v0, with
// partial_cut's glued harmonic field h at quadrature points in every tet. h
// comes straight from partial_cut_core in this process (live_ring.hpp: the
// solve of --export-xfem, a few seconds per crack, kept once solved); with
// --data DIR it is read from export files instead. Nothing is meshed here and
// nothing is solved for h by this program: it only runs the solid.
//
// What runs, per frame: one backward-euler step of the neo-Hookean solid
// (xfem/backward_euler.hpp), Newton on the cpu with a sparse Cholesky and a
// line search that keeps det F > 0 at every quadrature point and piece. With a
// front, the jump is held at zero on it (xfem_system::zero_jump_at; off with
// --front-jump free), so Newton runs in the null space of those rows. No
// gravity. The load is external forces on the ring's outer vertices:
//
//   mode I    hold the outer vertices on the - side of the crack plane (u_i
//             only, motus's fix_node: the material across the crack at those
//             nodes stays free) and pull the + side ones along n
//   mode II   the same pins, the + side pushed in the plane, normal to the front
//   mode III  the same pins, the + side pushed along the front
//   free      no pins: + side along +n, - side along -n, equal and opposite
//
// under a ramp-and-hold, step or sine profile. With a front the sine pulses
// between no load and the full load, (1 - cos) / 2, so the crack opens and
// shuts every period without ever being pushed shut. With no front (K3, or the
// synthetic through cut) the + half is a free body once cut: under a constant
// load it drifts away, so there the sine is the default profile, a cosine
// about zero (it starts at its peak, so its mean velocity is zero too) that
// keeps the half near. It still turns slowly: the loads keep their direction
// while the half rocks, and that leaves a small mean torque. There is no
// contact: where the faces are pushed together, they pass through each other,
// and the stats print the smallest gap along n so that shows.
//
// What is drawn: the deformed outer surface and the two crack sheets, every
// vertex carried by the field of its tet on its own side,
//
//     x = sum_a N_a(X) x_a  +  sum_{a enriched} N_a(X) psi_a(X) a_a,
//
// with psi from the vertex's side and h (xfem_system::deformed_point). Drawing
// the nodes alone cannot show a crack: both sides of it share them.
//
//   xfem_nodal_ring                              a window, K1 on ring 2, solved by partial_cut here
//   xfem_nodal_ring --case K1                    K1 | K2 | K3
//   xfem_nodal_ring --ring 1 --quad-level 0      check 4's ring 1; the coarsest quadrature (live only)
//   xfem_nodal_ring --data XFEM/data             the export files there instead (partial_cut --export-xfem)
//   xfem_nodal_ring --data XFEM/tests/data       the synthetic stand-ins (xfem_tests --write-synthetic)
//   xfem_nodal_ring --no-window --case K1 --frames 120
//   xfem_nodal_ring --case K2 --frames 90 --screenshot k2.png
//   xfem_nodal_ring --mode 2 --profile sine --load 600 --dt 0.01 --E 2e4 --nu 0.3 --rho 1000
//   xfem_nodal_ring --front-jump free            the front's jump left free (default: held at zero)

#include "live_ring.hpp"

#include "xfem/backward_euler.hpp"
#include "xfem/neo_hookean.hpp"
#include "xfem/ring_export.hpp"
#include "xfem/tet4.hpp"
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
#include <cstring>
#include <filesystem>
#include <map>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using namespace xfem;
namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// settings
// ---------------------------------------------------------------------------

enum class load_mode : int { opening = 0, sliding = 1, tearing = 2, free_pair = 3 };
const char* const load_mode_names[] = {"mode I: pull + side along n", "mode II: in plane, normal to the front",
                                       "mode III: along the front", "free: equal and opposite, no pins"};

enum class load_profile : int { ramp_hold = 0, step = 1, sine = 2 };
const char* const load_profile_names[] = {"ramp and hold", "step",
                                          "sine (0 to 1 with a front, about 0 without)"};

struct settings {
  load_mode mode = load_mode::opening;
  load_profile profile = load_profile::ramp_hold;
  double load = -1;           // total force on the loaded side, N; < 0: the crack's default_load
  double ramp_time = 1.0;     // s, ramp and hold
  double period = -1;         // s, sine; < 0: the crack's default_period
  double dt = 1.0 / 60;       // s
  double youngs = 2.0e4;      // Pa
  double poisson = 0.3;
  double density = 1000;      // kg / m^3
  double tolerance = 1e-8;    // Newton: |g|_inf below this times E h^2
  int max_iterations = 40;
  // Hold the jump at zero on the front (xfem_system::zero_jump_at). Off, the
  // jump there is held only by an energy that grows like log(1 / r_min) and
  // depends on where the quadrature points fall: see the README.
  bool zero_front_jump = true;
};

// The default load per crack, chosen so that mode I opens it by about 9% (K1)
// and 13% (K2) of the ring's size at the default material, with det F above
// 0.7: K2's crack is small (its front at 30%), so the uncracked ligament
// carries most of the pull and it gets 7.5 times K1's.
double default_load(const std::string& label) { return label == "K2" ? 3000.0 : 400.0; }

// The default sine period per crack. There is no damping but backward euler's
// and no contact, so a pulse that unloads faster than the ring can follow
// swings the crack past shut and the faces overlap. Of 2 to 6 s, 3 s keeps
// K1's mean gap the closest to open (it dips to -0.0014 in 8 s; K3's free
// half stays near too), 4 s K2's (-0.0020); each period is worse for the
// other crack.
double default_period(const std::string& label) { return label == "K2" ? 4.0 : 3.0; }

// The default profile: ramp and hold with a front, the sine without one, since
// a free half under a held load flies out of the view in a few seconds.
load_profile default_profile(bool through) { return through ? load_profile::sine : load_profile::ramp_hold; }

// ---------------------------------------------------------------------------
// the cases: solved here by partial_cut, or export files in a data folder
// ---------------------------------------------------------------------------

struct case_file {
  std::string label;
  std::string path;  // an export file; empty: solved here (live_ring)
  live_ring_options live;
};

// K1, K2, K3 on one ring, solved by partial_cut_core when first shown.
std::vector<case_file> live_cases(int ring, int quad_level) {
  std::vector<case_file> out;
  for (int k = 0; k < 3; ++k) out.push_back({"K" + std::to_string(k + 1), "", {ring, k, quad_level}});
  return out;
}

std::vector<case_file> find_cases(const fs::path& dir) {
  std::vector<case_file> out;
  for (const char* k : {"K1", "K2", "K3"})
  {
    const fs::path p = dir / (std::string("ring2_") + k + ".txt");
    if (fs::exists(p)) out.push_back({k, p.string()});
  }
  if (fs::exists(dir))
  {
    std::vector<fs::path> synthetic;
    for (const auto& entry : fs::directory_iterator(dir))
    {
      const std::string name = entry.path().filename().string();
      if (name.rfind("synthetic_", 0) == 0 && entry.path().extension() == ".txt") synthetic.push_back(entry.path());
    }
    std::sort(synthetic.begin(), synthetic.end());
    for (const fs::path& p : synthetic) out.push_back({p.stem().string(), p.string()});
  }
  return out;
}

// ---------------------------------------------------------------------------
// the simulation, with no window in it
// ---------------------------------------------------------------------------

struct frame_stats {
  int frame = 0;
  double time = 0;
  double load_factor = 0;
  newton_result newton;
  double elastic = 0, kinetic = 0;
  double opening = 0;
  double normal_gap = 0, mean_gap = 0;
  double front_jump = 0;  // the largest jump on the front, sum_a N_a a_a there
  double min_det = 1;
  double seconds = 0;
};

// Render vertices of one drawn group, compacted, with its triangles.
struct render_group {
  std::vector<int> ids;  // export render vertex ids
  std::vector<std::array<std::size_t, 3>> faces;
};

class simulation {
 public:
  ring_export data;
  std::unique_ptr<xfem_system> system;
  mechanical_state state;
  settings config;
  frame_stats last;
  std::vector<double> opening_history;

  vec3 n = vec3::UnitZ(), u = vec3::UnitX(), w = vec3::UnitY();
  std::vector<index> outer, plus, minus;
  render_group groups[3];
  std::vector<std::pair<int, int>> sheet_pairs;  // (+ copy, - copy) of every crack-sheet vertex
  std::vector<std::pair<index, vec3>> front_samples;
  double ring_size = 2;
  vec3 crack_center = vec3::Zero();  // the centroid of the crack sheet's vertices, rest space
  double force_scale = 1;  // E h^2
  index front_points = 0;  // points where the jump is held at zero (0: none)
  double time = 0;
  int frame = 0;

  void load(const case_file& c) {
    std::string source = c.path;
    if (c.path.empty())
    {
      const auto key = std::make_tuple(c.live.ring, c.live.crack, c.live.quad_level);
      auto it = solved.find(key);
      if (it == solved.end())
      {
        std::printf("solving %s on ring %d with partial_cut (quadrature level %d)...\n", c.label.c_str(), c.live.ring,
                    c.live.quad_level);
        it = solved.emplace(key, live_ring(c.live)).first;
      }
      data = it->second;
      source = "partial_cut, ring " + std::to_string(c.live.ring) + ", " + c.label;
    }
    else
    {
      data = read_ring_export(c.path);
    }
    const std::vector<std::string> problems = check(data);
    for (const std::string& p : problems) std::printf("warning: %s: %s\n", source.c_str(), p.c_str());
    build_frame();
    build_render();
  }

  // The rings partial_cut has solved in this run, by (ring, crack, quadrature level).
  std::map<std::tuple<int, int, int>, ring_export> solved;

  // A new system from the current settings: material, density, pins. Starts
  // at rest, at time 0.
  void rebuild() {
    enrichment_options options;
    xfem_enrichment enrichment =
        build_enrichment(data.mesh, data.plane_normal, data.plane_offset, data.tet_class, data.node_value, options);
    system = std::make_unique<xfem_system>(data.mesh, std::move(enrichment), data.quadrature,
                                           neo_hookean::from(config.youngs, config.poisson), config.density);
    if (config.mode != load_mode::free_pair)
    {
      for (const index v : minus) system->fix_node(v);
    }
    front_points = 0;
    if (config.zero_front_jump && !data.through())
    {
      front_points = zero_jump_on_front(*system, data.front[0], data.front[1]);
    }
    double edges = 0;
    int count = 0;
    for (const auto& t : data.mesh.tets)
    {
      for (int a = 0; a < 4; ++a)
      {
        for (int b = a + 1; b < 4; ++b, ++count) edges += (data.mesh.vertex(t[a]) - data.mesh.vertex(t[b])).norm();
      }
    }
    const double h = edges / count;
    force_scale = config.youngs * h * h;
    reset();
  }

  void reset() {
    state.x = system->rest_positions();
    state.v = vector_x::Zero(system->dofs());
    time = 0;
    frame = 0;
    opening_history.clear();
    last = frame_stats{};
    measure(last);
  }

  double load_factor(double t) const {
    switch (config.profile)
    {
      case load_profile::ramp_hold: return config.ramp_time > 0 ? std::min(t / config.ramp_time, 1.0) : 1.0;
      case load_profile::step: return 1.0;
      case load_profile::sine: {
        // A free half has to see a zero-mean force, and a cracked ligament
        // never a closing one: see the top of the file.
        const double c = std::cos(2 * std::numbers::pi * t / config.period);
        return data.through() ? c : 0.5 * (1 - c);
      }
    }
    return 0;
  }

  // The load factor averaged over [t0, t1] (Simpson, 8 panels), which is the
  // step's impulse over dt. Backward euler holds the load constant over a
  // step; with the impulse exact, a free body's velocity is exact at the end
  // of every step, so the zero-mean sine really gives it a zero mean velocity.
  // Sampling the end of each step instead adds a drift of -dt F / (2 m),
  // towards closing.
  double load_factor_over(double t0, double t1) const {
    constexpr int panels = 8;
    const double h = (t1 - t0) / panels;
    double sum = load_factor(t0) + load_factor(t1);
    for (int k = 1; k < panels; ++k) sum += (k % 2 ? 4.0 : 2.0) * load_factor(t0 + k * h);
    return sum / (3.0 * panels);
  }

  // The force at time t: the load's total split evenly over the loaded side's
  // outer vertices, on their u_i (the field at X_i: psi_i(X_i) = 0).
  vector_x external_force(double t) const { return external_force_at(load_factor(t)); }

  vector_x external_force_at(double factor) const {
    vector_x f = vector_x::Zero(system->dofs());
    const double scale = config.load * factor;
    const vec3 direction = config.mode == load_mode::sliding ? u : config.mode == load_mode::tearing ? w : n;
    for (const index v : plus) f.segment<3>(3 * v) += scale / double(plus.size()) * direction;
    if (config.mode == load_mode::free_pair)
    {
      for (const index v : minus) f.segment<3>(3 * v) -= scale / double(minus.size()) * n;
    }
    return f;
  }

  void step() {
    newton_options options;
    options.gradient_tolerance = config.tolerance * force_scale;
    options.max_iterations = config.max_iterations;
    backward_euler integrator(options);

    const auto start = std::chrono::steady_clock::now();
    const double t_next = time + config.dt;
    frame_stats s;
    s.newton = integrator.step(*system, state, config.dt, external_force_at(load_factor_over(time, t_next)));
    s.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    time = t_next;
    ++frame;
    s.frame = frame;
    s.time = time;
    s.load_factor = load_factor(time);
    measure(s);
    last = s;
    opening_history.push_back(s.opening);
    if (opening_history.size() > 600) opening_history.erase(opening_history.begin());
  }

  // Energy, the largest distance between the two copies of a crack-sheet
  // vertex, the smallest and the mean separation of those copies along the
  // rest normal (negative: the faces overlap, since there is no contact), and
  // the smallest det F over every point and piece.
  void measure(frame_stats& s) const {
    s.elastic = system->potential_energy(state.x);
    s.kinetic = 0.5 * state.v.dot(system->mass().cwiseProduct(state.v));
    s.min_det = system->minimum_determinant(state.x);
    double opening = 0;
    double gap = sheet_pairs.empty() ? 0.0 : DBL_MAX;
    double sum = 0;
    for (const auto& [p, m] : sheet_pairs)
    {
      const vec3 apart = deformed(p) - deformed(m);
      opening = std::max(opening, apart.norm());
      gap = std::min(gap, apart.dot(n));
      sum += apart.dot(n);
    }
    s.opening = opening;
    s.normal_gap = gap;
    s.mean_gap = sheet_pairs.empty() ? 0.0 : sum / double(sheet_pairs.size());
    // The jump on the front: the + face (h = 1) less the - face (h = 0) at
    // each front sample, which is sum_a N_a a_a there.
    s.front_jump = 0;
    for (const auto& [tet, X] : front_samples)
    {
      const vec3 j = system->deformed_point(state.x, tet, X, 1, 1.0) - system->deformed_point(state.x, tet, X, -1, 0.0);
      s.front_jump = std::max(s.front_jump, j.norm());
    }
  }

  vec3 deformed(int render_id) const {
    const render_vertex& v = data.render_vertices[static_cast<std::size_t>(render_id)];
    return system->deformed_point(state.x, v.tet, v.X, v.side, v.h);
  }

  std::vector<std::array<double, 3>> positions(const render_group& g) const {
    std::vector<std::array<double, 3>> out(g.ids.size());
    for (std::size_t i = 0; i < g.ids.size(); ++i)
    {
      const vec3 p = deformed(g.ids[i]);
      out[i] = {p[0], p[1], p[2]};
    }
    return out;
  }

  // The front, drawn at h = 1/2: the value just ahead of a straight front.
  // Only a guide for the eye; the field is singular on the front itself.
  std::vector<std::array<double, 3>> front_positions() const {
    std::vector<std::array<double, 3>> out;
    for (const auto& [tet, X] : front_samples)
    {
      const vec3 p = system->deformed_point(state.x, tet, X, 1, 0.5);
      out.push_back({p[0], p[1], p[2]});
    }
    return out;
  }

  std::string describe() const {
    const xfem_enrichment& en = system->enrichment();
    std::string s;
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "%s: %d dofs (%d nodes, %d enrichments); %d integration points (quadrature points and "
                  "pieces) of the export's %d quadrature points\n",
                  data.case_name.c_str(), int(system->dofs()), int(en.vertex_count()), int(en.enrichment_count()),
                  int(system->point_count()), int(data.quadrature.size()));
    s += buf;
    s += en.summary() + "\n";
    if (data.through())
    {
      s += "front: none (cut through)\n";
    }
    else if (front_points > 0)
    {
      std::snprintf(buf, sizeof buf, "front: jump held at 0 at %d points (rank %d), Newton in %d of the %d dofs\n",
                    int(front_points), int(system->constraint_rank()), int(system->dofs() - system->constraint_rank()),
                    int(system->dofs()));
      s += buf;
    }
    else
    {
      s += "front: jump free (held only by the log(1 / r) energy)\n";
    }
    s += "node types:";
    for (index i = 0; i < en.vertex_count(); ++i)
    {
      std::snprintf(buf, sizeof buf, " v%d %s%s", int(i), name(en.type_of(i)),
                    std::find(minus.begin(), minus.end(), i) != minus.end() && config.mode != load_mode::free_pair
                        ? " (pinned)"
                        : std::find(plus.begin(), plus.end(), i) != plus.end() ? " (loaded)" : "");
      s += buf;
    }
    s += "\ntets (class / path / integration points):";
    for (index e = 0; e < en.tet_count(); ++e)
    {
      std::snprintf(buf, sizeof buf, " %d:%s/%s/%d", int(e), name(en.class_of(e)), name(en.path_of(e)),
                    int(system->elements()[static_cast<std::size_t>(e)].points.size()));
      s += buf;
    }
    return s;
  }

 private:
  // The crack's frame, the outer vertices on each side, and the pins.
  void build_frame() {
    n = data.plane_normal.normalized();
    crack_center = vec3::Zero();
    int sheet_vertices = 0;
    for (const render_triangle& t : data.render_triangles)
    {
      if (t.group != 1) continue;
      for (const int v : t.v) crack_center += data.render_vertices[static_cast<std::size_t>(v)].X;
      sheet_vertices += 3;
    }
    if (sheet_vertices > 0) crack_center /= sheet_vertices;
    if (data.front.size() == 2)
    {
      w = data.front[1] - data.front[0];
      w = (w - w.dot(n) * n).normalized();
      u = w.cross(n);
      // u advances the front: from the crack towards the uncracked material.
      const vec3 middle = 0.5 * (data.front[0] + data.front[1]);
      if (sheet_vertices > 0 && (middle - crack_center).dot(u) < 0) u = -u;
      w = n.cross(u);
    }
    else
    {
      // A through cut has no front; use K1's front direction (0, 1, 0) in the plane.
      const vec3 f = std::abs(n.dot(vec3::UnitY())) < 0.9 ? vec3::UnitY() : vec3::UnitX();
      u = (f - f.dot(n) * n).normalized();
      w = n.cross(u);
    }

    // Outer vertices: those on a boundary face (a face no other tet shares).
    std::map<std::array<int, 3>, int> faces;
    for (const auto& t : data.mesh.tets)
    {
      for (int skip = 0; skip < 4; ++skip)
      {
        std::array<int, 3> f{};
        int m = 0;
        for (int a = 0; a < 4; ++a)
        {
          if (a != skip) f[static_cast<std::size_t>(m++)] = t[static_cast<std::size_t>(a)];
        }
        std::sort(f.begin(), f.end());
        ++faces[f];
      }
    }
    std::vector<char> on_boundary(static_cast<std::size_t>(data.mesh.vertex_count()), 0);
    for (const auto& [f, count] : faces)
    {
      if (count != 1) continue;
      for (const int v : f) on_boundary[static_cast<std::size_t>(v)] = 1;
    }
    outer.clear();
    plus.clear();
    minus.clear();
    for (index v = 0; v < data.mesh.vertex_count(); ++v)
    {
      if (!on_boundary[static_cast<std::size_t>(v)]) continue;
      outer.push_back(v);
      (n.dot(data.mesh.vertex(v)) - data.plane_offset > 0 ? plus : minus).push_back(v);
    }

    ring_size = 0;
    for (const vec3& a : data.mesh.vertices)
    {
      for (const vec3& b : data.mesh.vertices) ring_size = std::max(ring_size, (a - b).norm());
    }

    // The front, sampled, each sample in a tet that holds it (a front tet first).
    front_samples.clear();
    if (data.front.size() == 2)
    {
      for (int i = 0; i <= 40; ++i)
      {
        const vec3 X = data.front[0] + (i / 40.0) * (data.front[1] - data.front[0]);
        index best = -1;
        double best_score = -1e30;
        for (index e = 0; e < data.mesh.tet_count(); ++e)
        {
          vec3 c[4];
          data.mesh.corners(e, c);
          const vec4 N = tet4::shape_values(tet4::build_rest_state(c), X);
          const double score = N.minCoeff() + (data.tet_class[static_cast<std::size_t>(e)] == 2 ? 1e-9 : 0);
          if (score > best_score)
          {
            best_score = score;
            best = e;
          }
        }
        front_samples.push_back({best, X});
      }
    }
  }

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

    // Pair the two copies of every crack-sheet vertex: same tet, same X.
    sheet_pairs.clear();
    std::map<std::tuple<int, long long, long long, long long>, int> minus_copy;
    const auto key = [](const render_vertex& v) {
      const auto q = [](double x) { return static_cast<long long>(std::llround(x * 1e9)); };
      return std::make_tuple(v.tet, q(v.X[0]), q(v.X[1]), q(v.X[2]));
    };
    for (const int id : groups[2].ids) minus_copy[key(data.render_vertices[static_cast<std::size_t>(id)])] = id;
    for (const int id : groups[1].ids)
    {
      const auto found = minus_copy.find(key(data.render_vertices[static_cast<std::size_t>(id)]));
      if (found != minus_copy.end()) sheet_pairs.push_back({id, found->second});
    }
  }
};

void print_frame(const simulation& sim) {
  const frame_stats& s = sim.last;
  std::printf("frame %4d  t %6.3f  load %+.3f  newton %2d  |g| %.2e  elastic %.4e  kinetic %.4e  "
              "opening %.4f (%.1f%%)  gap min %+.4f mean %+.4f  front jump %.1e  min det F %.3f  %.0f ms%s%s\n",
              s.frame, s.time, s.load_factor, int(s.newton.iterations), s.newton.final_gradient_norm, s.elastic,
              s.kinetic, s.opening, 100 * s.opening / sim.ring_size, s.normal_gap, s.mean_gap, s.front_jump, s.min_det,
              1e3 * s.seconds, s.newton.converged ? "" : "  NOT CONVERGED: ",
              s.newton.converged ? "" : s.newton.failure.c_str());
}

// ---------------------------------------------------------------------------
// the window
// ---------------------------------------------------------------------------

std::array<double, 3> a3(const vec3& v) { return {v[0], v[1], v[2]}; }
glm::vec3 g3(const vec3& v) { return glm::vec3(float(v[0]), float(v[1]), float(v[2])); }

struct viewer {
  simulation& sim;
  std::vector<case_file> cases;
  int case_index = 0;
  bool playing = false;
  bool show_points = false;
  settings pending;  // edited in the panel; material changes apply on reset

  void register_all() {
    polyscope::removeAllStructures();
    // The ring's surface shows h on polyscope's coolwarm map over [0, 1], and
    // the two crack sheets take the map's two ends: the + sheet is h = 1 (red),
    // the - sheet h = 0 (blue). Behind the front the surface's colour jumps
    // across the crack mouth; beyond it h runs smoothly through white.
    const char* names[3] = {"ring surface", "crack + sheet", "crack - sheet"};
    const glm::vec3 colours[3] = {{0.80f, 0.78f, 0.72f}, {0.706f, 0.016f, 0.150f}, {0.230f, 0.299f, 0.754f}};
    for (int k = 0; k < 3; ++k)
    {
      if (sim.groups[k].faces.empty()) continue;
      polyscope::SurfaceMesh* m = polyscope::registerSurfaceMesh(names[k], sim.positions(sim.groups[k]), sim.groups[k].faces);
      m->setSurfaceColor(colours[k]);
      m->setBackFacePolicy(polyscope::BackFacePolicy::Different);
      if (k != 0) continue;
      std::vector<double> h;
      for (const int id : sim.groups[k].ids) h.push_back(sim.data.render_vertices[static_cast<std::size_t>(id)].h);
      polyscope::SurfaceVertexScalarQuantity* field = m->addVertexScalarQuantity("h", h);
      field->setColorMap("coolwarm");
      field->setMapRange({0.0, 1.0});
      field->setEnabled(true);
      m->setTransparency(0.55f);
    }

    // The rest ring: every tet edge, once.
    std::vector<std::array<double, 3>> nodes;
    for (const vec3& v : sim.data.mesh.vertices) nodes.push_back(a3(v));
    std::vector<std::array<std::size_t, 2>> edges;
    for (const auto& t : sim.data.mesh.tets)
    {
      for (int a = 0; a < 4; ++a)
      {
        for (int b = a + 1; b < 4; ++b)
        {
          std::array<std::size_t, 2> e = {std::size_t(std::min(t[a], t[b])), std::size_t(std::max(t[a], t[b]))};
          if (std::find(edges.begin(), edges.end(), e) == edges.end()) edges.push_back(e);
        }
      }
    }
    polyscope::CurveNetwork* rest = polyscope::registerCurveNetwork("rest ring", nodes, edges);
    rest->setColor({0.55f, 0.55f, 0.55f});
    rest->setRadius(0.0015f, false);

    if (!sim.front_samples.empty())
    {
      polyscope::CurveNetwork* front = polyscope::registerCurveNetworkLine("front", sim.front_positions());
      front->setColor({0.95f, 0.75f, 0.10f});
      front->setRadius(0.008f, false);
    }

    std::vector<std::array<double, 3>> points;
    std::vector<double> h;
    for (const quadrature_point& q : sim.data.quadrature)
    {
      points.push_back(a3(q.X));
      h.push_back(q.h);
    }
    polyscope::PointCloud* cloud = polyscope::registerPointCloud("quadrature points (h)", points);
    cloud->setPointRadius(0.004, false);
    cloud->addScalarQuantity("h", h)->setEnabled(true);
    cloud->setEnabled(show_points);

    std::vector<double> role;
    for (index v = 0; v < sim.data.mesh.vertex_count(); ++v)
    {
      const bool pinned = sim.system->is_fixed(3 * v);
      const bool loaded = std::find(sim.plus.begin(), sim.plus.end(), v) != sim.plus.end() ||
                          (sim.config.mode == load_mode::free_pair &&
                           std::find(sim.minus.begin(), sim.minus.end(), v) != sim.minus.end());
      role.push_back(pinned ? 0.0 : loaded ? 1.0 : 0.5);
    }
    polyscope::PointCloud* nodes_cloud = polyscope::registerPointCloud("nodes (u_i)", nodes);
    nodes_cloud->setPointRadius(0.018, false);
    nodes_cloud->addScalarQuantity("pinned 0, free 0.5, loaded 1", role)->setEnabled(true);
    update();
  }

  void update() {
    const char* names[3] = {"ring surface", "crack + sheet", "crack - sheet"};
    for (int k = 0; k < 3; ++k)
    {
      if (sim.groups[k].faces.empty()) continue;
      polyscope::getSurfaceMesh(names[k])->updateVertexPositions(sim.positions(sim.groups[k]));
    }
    if (!sim.front_samples.empty()) polyscope::getCurveNetwork("front")->updateNodePositions(sim.front_positions());
    if (polyscope::getPointCloud("quadrature points (h)")->isEnabled())
    {
      std::vector<std::array<double, 3>> points;
      for (const quadrature_point& q : sim.data.quadrature)
      {
        points.push_back(a3(sim.system->deformed_point(sim.state.x, q.tet, q.X, q.side, q.h)));
      }
      polyscope::getPointCloud("quadrature points (h)")->updatePointPositions(points);
    }
    std::vector<std::array<double, 3>> nodes;
    std::vector<std::array<double, 3>> force;
    const vector_x f = sim.external_force(sim.time);
    // Drawn to scale: the full load on one vertex would be an arrow 0.25 long.
    const double scale = sim.config.load > 0 ? 0.25 / sim.config.load * double(std::max<std::size_t>(sim.plus.size(), 1)) : 0;
    for (index v = 0; v < sim.data.mesh.vertex_count(); ++v)
    {
      nodes.push_back(a3(sim.state.x.segment<3>(3 * v)));
      force.push_back(a3(scale * f.segment<3>(3 * v)));
    }
    polyscope::PointCloud* cloud = polyscope::getPointCloud("nodes (u_i)");
    cloud->updatePointPositions(nodes);
    polyscope::PointCloudVectorQuantity* arrows = cloud->addVectorQuantity("load", force, polyscope::VectorType::AMBIENT);
    arrows->setVectorRadius(0.008, false);
    arrows->setVectorColor({0.10f, 0.60f, 0.55f});
    arrows->setEnabled(true);
  }

  // Along the front, a little above the crack plane and a little behind the
  // front, so the crack shows as a wedge opening away from the front; aimed
  // part way from the ring's middle to the crack, so a small crack (K2) is not
  // left in a corner.
  void camera() {
    const vec3 target = 0.4 * sim.crack_center;
    const vec3 eye = target + 3.1 * (0.92 * sim.w - 0.25 * sim.u + 0.30 * sim.n).normalized();
    polyscope::view::lookAt(g3(eye), g3(target), g3(sim.n));
  }

  void load_case(int index) {
    case_index = index;
    pending.load = default_load(cases[static_cast<std::size_t>(index)].label);
    pending.period = default_period(cases[static_cast<std::size_t>(index)].label);
    sim.load(cases[static_cast<std::size_t>(index)]);
    pending.profile = default_profile(sim.data.through());
    sim.config = pending;
    sim.rebuild();
    register_all();
    camera();
  }

  void ui() {
    ImGui::PushItemWidth(170);
    ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + 440);
    if (ImGui::BeginCombo("crack", cases[static_cast<std::size_t>(case_index)].label.c_str()))
    {
      for (int i = 0; i < int(cases.size()); ++i)
      {
        if (ImGui::Selectable(cases[static_cast<std::size_t>(i)].label.c_str(), i == case_index) && i != case_index)
        {
          load_case(i);
        }
      }
      ImGui::EndCombo();
    }

    int mode = static_cast<int>(pending.mode);
    if (ImGui::Combo("load mode", &mode, load_mode_names, 4))
    {
      // The pins change with the mode, so it restarts.
      pending.mode = static_cast<load_mode>(mode);
      sim.config = pending;
      sim.rebuild();
      register_all();
    }
    int profile = static_cast<int>(pending.profile);
    if (ImGui::Combo("profile", &profile, load_profile_names, 3)) pending.profile = static_cast<load_profile>(profile);
    ImGui::InputDouble("load [N]", &pending.load, 50, 200, "%.1f");
    pending.load = std::max(0.0, pending.load);
    if (pending.profile == load_profile::ramp_hold) ImGui::InputDouble("ramp time [s]", &pending.ramp_time, 0.1, 0.5, "%.2f");
    if (pending.profile == load_profile::sine) ImGui::InputDouble("period [s]", &pending.period, 0.1, 0.5, "%.2f");
    // The load and its profile apply at once; the rest on reset.
    sim.config.load = pending.load;
    sim.config.profile = pending.profile;
    sim.config.ramp_time = std::max(0.0, pending.ramp_time);
    sim.config.period = std::max(1e-3, pending.period);

    if (ImGui::Button(playing ? "pause" : "play")) playing = !playing;
    ImGui::SameLine();
    if (ImGui::Button("step")) advance();
    ImGui::SameLine();
    if (ImGui::Button("reset"))
    {
      sim.config = pending;
      sim.rebuild();
      register_all();
      camera();
    }
    ImGui::SameLine();
    // Polyscope's own Reset View goes to its home view, not along the front.
    if (ImGui::Button("crack view")) camera();

    if (ImGui::TreeNodeEx("material and solver (apply on reset)", ImGuiTreeNodeFlags_DefaultOpen))
    {
      if (!sim.data.through()) ImGui::Checkbox("jump = 0 on the front", &pending.zero_front_jump);
      ImGui::InputDouble("dt [s]", &pending.dt, 0.001, 0.005, "%.4f");
      ImGui::InputDouble("E [Pa]", &pending.youngs, 1e3, 1e4, "%.4g");
      ImGui::InputDouble("nu", &pending.poisson, 0.01, 0.05, "%.3f");
      ImGui::InputDouble("rho [kg/m^3]", &pending.density, 50, 200, "%.1f");
      ImGui::InputDouble("Newton tol (x E h^2)", &pending.tolerance, 0, 0, "%.1e");
      pending.dt = std::max(1e-5, pending.dt);
      pending.poisson = std::clamp(pending.poisson, 0.0, 0.49);
      pending.youngs = std::max(1.0, pending.youngs);
      pending.density = std::max(1e-3, pending.density);
      pending.tolerance = std::max(1e-14, pending.tolerance);
      ImGui::TreePop();
    }
    if (ImGui::Checkbox("show quadrature points (h)", &show_points))
    {
      polyscope::getPointCloud("quadrature points (h)")->setEnabled(show_points);
      update();
    }
    ImGui::PopItemWidth();

    ImGui::Separator();
    const frame_stats& s = sim.last;
    if (ImGui::TreeNode("ring, enrichment and export"))
    {
      ImGui::TextUnformatted(sim.describe().c_str());
      ImGui::TextUnformatted(sim.data.settings.c_str());
      ImGui::TreePop();
    }
    if (sim.data.through())
    {
      ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f),
                         "No front: the crack cuts the ring through, so the + half is a free body.\n"
                         "Under a constant load it drifts away; the sine profile keeps it near.");
    }
    ImGui::Text("frame %d  t = %.3f s  load factor %+.3f", s.frame, s.time, s.load_factor);
    ImGui::Text("Newton: %d iterations, |g| = %.2e (tol %.1e)%s", int(s.newton.iterations),
                s.newton.final_gradient_norm, sim.config.tolerance * sim.force_scale,
                s.frame == 0 ? "" : s.newton.converged ? "" : "  NOT CONVERGED");
    if (!s.newton.converged && s.frame > 0) ImGui::Text("  %s", s.newton.failure.c_str());
    std::string residuals = "residuals:";
    for (const double r : s.newton.residuals)
    {
      char buf[32];
      std::snprintf(buf, sizeof buf, " %.1e", r);
      residuals += buf;
    }
    ImGui::TextUnformatted(residuals.c_str());
    ImGui::Text("energy: elastic %.4e, kinetic %.4e J", s.elastic, s.kinetic);
    ImGui::Text("max crack opening %.4f (%.1f%% of the ring's %.2f)", s.opening, 100 * s.opening / sim.ring_size,
                sim.ring_size);
    // Negative gap: faces pushed through each other (there is no contact) when
    // the mean is negative too; with the mean open, a local overlap, which
    // with the front's jump left free sits at the front (README, Known limits).
    const char* overlap = s.normal_gap >= -1e-4 ? ""
                          : s.mean_gap < 0      ? "  (faces pushed through: no contact)"
                          : sim.front_points > 0 ? "  (a local overlap, the faces open on average)"
                                                 : "  (overlap near the front: its jump is free, see README)";
    ImGui::Text("gap along n: smallest %+.4f, mean %+.4f%s", s.normal_gap, s.mean_gap, overlap);
    if (!sim.data.through())
    {
      ImGui::Text("jump on the front: %.1e%s", s.front_jump, sim.front_points > 0 ? " (held at 0)" : " (free)");
    }
    ImGui::Text("min det F %.4f   step %.0f ms", s.min_det, 1e3 * s.seconds);
    if (!sim.opening_history.empty())
    {
      std::vector<float> history(sim.opening_history.begin(), sim.opening_history.end());
      ImGui::PlotLines("opening", history.data(), int(history.size()), 0, nullptr, 0.0f, FLT_MAX, ImVec2(360, 70));
    }
    ImGui::PopTextWrapPos();

    if (playing) advance();
  }

  void advance() {
    sim.step();
    update();
  }
};

// ---------------------------------------------------------------------------

int run(int argc, char** argv) {
  settings config;
  std::string case_name, data_dir, screenshot;
  int ring = 2, quad_level = 1;
  bool window = true, with_ui = false;
  int frames = 120;
  bool frames_given = false, profile_given = false;

  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    const auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument(arg + " needs a value");
      return argv[++i];
    };
    if (arg == "--no-window") window = false;
    else if (arg == "--case") case_name = next();
    else if (arg == "--data") data_dir = next();
    else if (arg == "--ring") ring = std::stoi(next());
    else if (arg == "--quad-level") quad_level = std::stoi(next());
    else if (arg == "--frames") { frames = std::stoi(next()); frames_given = true; }
    else if (arg == "--screenshot") screenshot = next();
    else if (arg == "--with-ui") with_ui = true;
    else if (arg == "--mode") config.mode = static_cast<load_mode>(std::clamp(std::stoi(next()) - 1, 0, 3));
    else if (arg == "--free") config.mode = load_mode::free_pair;
    else if (arg == "--profile")
    {
      const std::string p = next();
      config.profile = p == "step" ? load_profile::step : p == "sine" ? load_profile::sine : load_profile::ramp_hold;
      profile_given = true;
    }
    else if (arg == "--front-jump")
    {
      const std::string p = next();
      if (p != "zero" && p != "free") throw std::invalid_argument("--front-jump: zero or free");
      config.zero_front_jump = p == "zero";
    }
    else if (arg == "--load") config.load = std::stod(next());
    else if (arg == "--ramp") config.ramp_time = std::stod(next());
    else if (arg == "--period") config.period = std::stod(next());
    else if (arg == "--dt") config.dt = std::stod(next());
    else if (arg == "--E") config.youngs = std::stod(next());
    else if (arg == "--nu") config.poisson = std::stod(next());
    else if (arg == "--rho") config.density = std::stod(next());
    else if (arg == "--tol") config.tolerance = std::stod(next());
    else
    {
      std::fprintf(stderr,
                   "usage: %s [--case K1|K2|K3] [--ring 1|2] [--quad-level 0..4]\n"
                   "          [--data DIR (export files instead; XFEM/tests/data: synthetic_partial|synthetic_through)]\n"
                   "          [--no-window] [--frames N] [--screenshot FILE [--with-ui]]\n"
                   "          [--mode 1|2|3] [--free] [--profile ramp|step|sine] [--load N] [--ramp s] [--period s]\n"
                   "          [--dt s] [--E Pa] [--nu v] [--rho kg/m^3] [--tol t] [--front-jump zero|free]\n",
                   argv[0]);
      return 2;
    }
  }

  if (ring != 1 && ring != 2) throw std::invalid_argument("--ring: 1 or 2");
  if (quad_level < 0 || quad_level > 4) throw std::invalid_argument("--quad-level: 0 to 4");
  const fs::path dir = data_dir.empty() ? fs::path(XFEM_DATA_DIR) : fs::path(data_dir);
  const std::vector<case_file> cases = data_dir.empty() ? live_cases(ring, quad_level) : find_cases(dir);
  if (cases.empty())
  {
    std::fprintf(stderr, "no export files in %s: expected ring2_K1.txt, ring2_K2.txt, ring2_K3.txt or synthetic_*.txt\n"
                         "(xfem_tests --write-synthetic writes the synthetic ones to XFEM/tests/data)\n",
                 dir.string().c_str());
    return 1;
  }
  int case_index = 0;
  if (!case_name.empty())
  {
    case_index = -1;
    for (int i = 0; i < int(cases.size()); ++i)
    {
      if (cases[static_cast<std::size_t>(i)].label == case_name) case_index = i;
    }
    if (case_index < 0)
    {
      std::fprintf(stderr, "no case '%s' in %s; found:", case_name.c_str(), dir.string().c_str());
      for (const case_file& c : cases) std::fprintf(stderr, " %s", c.label.c_str());
      std::fprintf(stderr, "\n");
      return 1;
    }
  }

  if (config.load < 0) config.load = default_load(cases[static_cast<std::size_t>(case_index)].label);
  if (config.period <= 0) config.period = default_period(cases[static_cast<std::size_t>(case_index)].label);
  simulation sim;
  sim.load(cases[static_cast<std::size_t>(case_index)]);
  if (!profile_given) config.profile = default_profile(sim.data.through());
  sim.config = config;
  sim.rebuild();
  const case_file& shown = cases[static_cast<std::size_t>(case_index)];
  const std::string source =
      shown.path.empty() ? "partial_cut, ring " + std::to_string(shown.live.ring) + ", " + shown.label : shown.path;
  std::printf("%s\n%s\n", source.c_str(), sim.describe().c_str());
  std::printf("settings: %s\n", sim.data.settings.c_str());
  std::printf("load: %s, %s (ramp %.2f s, period %.2f s), %.1f N; dt %.4f s, E %.4g Pa, nu %.3f, rho %.1f kg/m^3, "
              "Newton tol %.2e N\n",
              load_mode_names[static_cast<int>(config.mode)], load_profile_names[static_cast<int>(config.profile)],
              config.ramp_time, config.period, config.load, config.dt, config.youngs, config.poisson, config.density,
              config.tolerance * sim.force_scale);
  if (sim.data.through())
  {
    std::printf("note: no front, so the + half is a free body once cut; under a constant load it drifts\n");
  }

  if (!window)
  {
    for (int f = 0; f < frames; ++f)
    {
      sim.step();
      print_frame(sim);
    }
    return 0;
  }

  polyscope::options::programName = "xfem_nodal_ring";
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::options::transparencyMode = polyscope::TransparencyMode::Pretty;
  if (!screenshot.empty()) polyscope::options::usePrefsFile = false;
  polyscope::init();
  polyscope::view::setNavigateStyle(polyscope::NavigateStyle::Arcball);
  if (!screenshot.empty()) polyscope::view::setWindowSize(1280, 960);

  viewer view{sim, cases};
  view.case_index = case_index;
  view.pending = config;
  view.register_all();
  view.camera();

  if (!screenshot.empty())
  {
    if (!frames_given) frames = 90;
    for (int f = 0; f < frames; ++f)
    {
      sim.step();
      print_frame(sim);
    }
    view.update();
    polyscope::ScreenshotOptions options;
    options.transparentBackground = false;
    if (with_ui)
    {
      polyscope::state::userCallback = [&view] { view.ui(); };
      options.includeUI = true;
    }
    polyscope::screenshot(screenshot, options);
    std::printf("screenshot written: %s\n", screenshot.c_str());
    return 0;
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
