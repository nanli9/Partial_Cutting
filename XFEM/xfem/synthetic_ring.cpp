#include "xfem/synthetic_ring.hpp"

#include "xfem/xfem_cutter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>
#include <tuple>

namespace xfem {

namespace {

using tet_corners = std::array<vec3, 4>;

// Sutherland-Hodgman against one half-space f(X) <= 0, f affine: keep the
// vertices with f <= 0 and add the crossing on every edge whose ends are
// strictly on opposite sides.
template <class affine>
std::vector<vec3> clip(const std::vector<vec3>& polygon, affine f) {
  std::vector<vec3> out;
  for (std::size_t i = 0; i < polygon.size(); ++i)
  {
    const vec3& a = polygon[i];
    const vec3& b = polygon[(i + 1) % polygon.size()];
    const real fa = f(a);
    const real fb = f(b);
    if (fa <= real(0)) out.push_back(a);
    if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0)) out.push_back(a + (fa / (fa - fb)) * (b - a));
  }
  return out;
}

// The section of the octahedron |x| + |y| + |z| <= 1 by the plane n.X = d: a
// square of the plane larger than the octahedron, CCW in (u, w), clipped by
// its 8 half-spaces s.X <= 1 (partial_cut's buildCutOnDomain).
std::vector<vec3> octahedron_section(const vec3& n, real d, const vec3& u, const vec3& w) {
  const real R = 5;
  const vec3 o = d * n;
  std::vector<vec3> P = {o - R * u - R * w, o + R * u - R * w, o + R * u + R * w, o - R * u + R * w};
  for (int k = 0; k < 8 && !P.empty(); ++k)
  {
    const vec3 s(k & 1 ? -1 : 1, k & 2 ? -1 : 1, k & 4 ? -1 : 1);
    P = clip(P, [&](const vec3& X) { return s.dot(X) - real(1); });
  }
  return P;
}

// The crack's frame and its harmonic field. a = u.X - s along the advance
// (negative behind the front, on the crack), b = n.X - d across the crack.
struct crack_frame {
  vec3 n, u, w;
  real d = 0;
  real s = std::numeric_limits<real>::infinity();  // the front; +inf for a through cut
  bool through = false;

  real phi(const vec3& X) const { return n.dot(X) - d; }
  real distance_to_front(const vec3& X) const { return std::hypot(u.dot(X) - s, phi(X)); }

  // h and grad h at X, on `side` for a point on the plane. "On the plane" is
  // to 1e-12: a vertex the clipping put there has a phi of rounding size and
  // either sign, and atan2 would hand its + copy the - face's value.
  void field(const vec3& X, int side, real& h, vec3& grad) const {
    if (through)
    {
      h = side > 0 ? real(1) : real(0);
      grad = vec3::Zero();
      return;
    }
    const real a = u.dot(X) - s;
    real b = phi(X);
    if (std::abs(b) < real(1e-12)) b = real(0);
    real theta;
    if (b == real(0)) theta = a < 0 ? (side > 0 ? std::numbers::pi : -std::numbers::pi) : real(0);
    else theta = std::atan2(b, a);
    h = real(0.5) + theta / (real(2) * std::numbers::pi);
    const real r2 = a * a + b * b;
    grad = r2 > real(0) ? vec3((a * n - b * u) / (real(2) * std::numbers::pi * r2)) : vec3(vec3::Zero());
  }
};

// 8-way midpoint subdivision (four corner tets and the inner octahedron cut
// into four around its m02-m13 diagonal), `levels` times.
void refine(const tet_corners& t, int levels, std::vector<tet_corners>& out) {
  if (levels <= 0)
  {
    out.push_back(t);
    return;
  }
  const vec3 m01 = (t[0] + t[1]) / 2, m02 = (t[0] + t[2]) / 2, m03 = (t[0] + t[3]) / 2;
  const vec3 m12 = (t[1] + t[2]) / 2, m13 = (t[1] + t[3]) / 2, m23 = (t[2] + t[3]) / 2;
  const tet_corners children[8] = {
      {t[0], m01, m02, m03}, {m01, t[1], m12, m13}, {m02, m12, t[2], m23}, {m03, m13, m23, t[3]},
      {m02, m13, m01, m03},  {m02, m13, m03, m23},  {m02, m13, m23, m12},  {m02, m13, m12, m01}};
  for (const tet_corners& c : children) refine(c, levels - 1, out);
}

// The sub-tets of tet `corners` on one side of the plane, in X: motus's +3/-1
// and 2+/2- tessellation, or the whole tet when the plane misses it.
std::vector<tet_corners> side_pieces(const vec3 (&corners)[4], const vec4& phi, int side) {
  int positive = 0;
  for (int a = 0; a < 4; ++a) positive += phi[a] > 0 ? 1 : 0;
  std::vector<tet_corners> out;
  if (positive == 0 || positive == 4)
  {
    if ((positive == 4) == (side > 0)) out.push_back({corners[0], corners[1], corners[2], corners[3]});
    return out;
  }
  vec4 sub[3][4];
  const int count = side > 0 ? positive_piece(phi, sub) : negative_piece(phi, sub);
  for (int k = 0; k < count; ++k)
  {
    tet_corners t;
    for (int c = 0; c < 4; ++c)
    {
      t[static_cast<std::size_t>(c)] = sub[k][c][0] * corners[0] + sub[k][c][1] * corners[1] +
                                       sub[k][c][2] * corners[2] + sub[k][c][3] * corners[3];
    }
    out.push_back(t);
  }
  return out;
}

// Render vertices welded per (tet, side, position).
class render_builder {
 public:
  explicit render_builder(ring_export& out) : out_(out) {}

  int vertex(int tet, int side, const vec3& X, real h) {
    const auto q = [](real v) { return static_cast<long long>(std::llround(v * 1e9)); };
    const auto key = std::make_tuple(tet, side, q(X[0]), q(X[1]), q(X[2]));
    const auto found = index_.find(key);
    if (found != index_.end()) return found->second;
    const int id = static_cast<int>(out_.render_vertices.size());
    out_.render_vertices.push_back({tet, side, X, h});
    index_.emplace(key, id);
    return id;
  }

  void triangle(int a, int b, int c, int group) {
    if (a == b || b == c || a == c) return;
    render_triangle t;
    t.v[0] = a;
    t.v[1] = b;
    t.v[2] = c;
    t.group = group;
    out_.render_triangles.push_back(t);
  }

 private:
  ring_export& out_;
  std::map<std::tuple<int, int, long long, long long, long long>, int> index_;
};

// The k x k subdivision of triangle (A, B, C), keeping its winding.
template <class emit>
void subdivide(const vec3& A, const vec3& B, const vec3& C, int k, emit triangle) {
  const auto P = [&](int i, int j) { return vec3(A + (real(i) / k) * (B - A) + (real(j) / k) * (C - A)); };
  for (int i = 0; i < k; ++i)
  {
    for (int j = 0; i + j < k; ++j)
    {
      triangle(P(i, j), P(i + 1, j), P(i, j + 1));
      if (i + j + 1 < k) triangle(P(i + 1, j), P(i + 1, j + 1), P(i, j + 1));
    }
  }
}

}  // namespace

tet_mesh nodal_ring_mesh(const vec3& node) {
  tet_mesh mesh;
  mesh.vertices = {node, vec3(1, 0, 0), vec3(-1, 0, 0), vec3(0, 1, 0), vec3(0, -1, 0), vec3(0, 0, 1), vec3(0, 0, -1)};
  for (int o = 0; o < 8; ++o)
  {
    std::array<int, 4> g = {0, o & 1 ? 2 : 1, o & 2 ? 4 : 3, o & 4 ? 6 : 5};
    if (tet_volume(mesh.vertex(g[0]), mesh.vertex(g[1]), mesh.vertex(g[2]), mesh.vertex(g[3])) < 0)
    {
      std::swap(g[2], g[3]);
    }
    mesh.tets.push_back(g);
  }
  return mesh;
}

ring_export make_synthetic_ring(const synthetic_options& options) {
  ring_export out;
  out.case_name = options.through ? "synthetic_through" : "synthetic_partial";
  out.ring = 2;
  out.g = 1;
  out.mesh = nodal_ring_mesh(options.node);

  crack_frame crack;
  crack.n = options.normal.normalized();
  crack.d = options.offset;
  crack.through = options.through;
  const vec3 f = options.front_direction;
  crack.u = (f - f.dot(crack.n) * crack.n).normalized();
  crack.w = crack.n.cross(crack.u);
  out.plane_normal = crack.n;
  out.plane_offset = crack.d;

  // The front: at front_t of the section's extent along u, and its two ends
  // where it leaves the octahedron.
  const std::vector<vec3> section = octahedron_section(crack.n, crack.d, crack.u, crack.w);
  if (section.size() < 3) throw std::invalid_argument("xfem: make_synthetic_ring: the plane misses the ring");
  if (!options.through)
  {
    real lo = std::numeric_limits<real>::max(), hi = -lo;
    for (const vec3& p : section)
    {
      lo = std::min(lo, crack.u.dot(p));
      hi = std::max(hi, crack.u.dot(p));
    }
    crack.s = lo + options.front_t * (hi - lo);
    for (std::size_t i = 0; i < section.size(); ++i)
    {
      const vec3& a = section[i];
      const vec3& b = section[(i + 1) % section.size()];
      const real fa = crack.u.dot(a) - crack.s;
      const real fb = crack.u.dot(b) - crack.s;
      if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0)) out.front.push_back(a + (fa / (fa - fb)) * (b - a));
    }
    if (out.front.size() != 2) throw std::runtime_error("xfem: make_synthetic_ring: the front is not one segment");
  }

  const tet_mesh& mesh = out.mesh;
  const index T = mesh.tet_count();

  // The classes: from each tet's section by the plane, against the front.
  std::vector<std::vector<vec3>> crack_in(static_cast<std::size_t>(T));
  for (index e = 0; e < T; ++e)
  {
    vec3 corners[4];
    mesh.corners(e, corners);
    vec4 phi;
    for (int a = 0; a < 4; ++a) phi[a] = crack.phi(corners[a]);
    vec4 polygon[4];
    const int count = cut_polygon(phi, polygon);
    int c = 0;
    if (count > 0)
    {
      std::vector<vec3> P;
      for (int i = 0; i < count; ++i)
      {
        P.push_back(polygon[i][0] * corners[0] + polygon[i][1] * corners[1] + polygon[i][2] * corners[2] +
                    polygon[i][3] * corners[3]);
      }
      real lo = std::numeric_limits<real>::max(), hi = -lo;
      for (const vec3& p : P)
      {
        lo = std::min(lo, crack.u.dot(p));
        hi = std::max(hi, crack.u.dot(p));
      }
      if (options.through || hi <= crack.s) c = 1;
      else if (lo >= crack.s) c = 0;
      else c = 2;
      // The crack in this tet: its section behind the front.
      if (c != 0)
      {
        crack_in[static_cast<std::size_t>(e)] =
            options.through ? P : clip(P, [&](const vec3& X) { return crack.u.dot(X) - crack.s; });
      }
    }
    out.tet_class.push_back(c);
  }

  for (index v = 0; v < mesh.vertex_count(); ++v)
  {
    real h;
    vec3 grad;
    crack.field(mesh.vertex(v), crack.phi(mesh.vertex(v)) > 0 ? 1 : -1, h, grad);
    out.node_value.push_back(h);
  }

  // The quadrature: split, refine, refine near the front, Gauss in every leaf.
  constexpr real alpha = real(0.5854101966249685);
  constexpr real beta = real(0.1381966011250105);
  for (index e = 0; e < T; ++e)
  {
    vec3 corners[4];
    mesh.corners(e, corners);
    vec4 phi;
    for (int a = 0; a < 4; ++a) phi[a] = crack.phi(corners[a]);
    for (const int side : {1, -1})
    {
      for (const tet_corners& piece : side_pieces(corners, phi, side))
      {
        std::vector<tet_corners> leaves;
        refine(piece, options.quad_level, leaves);
        if (!options.through && options.front_levels > 0)
        {
          std::vector<tet_corners> refined;
          for (const tet_corners& leaf : leaves)
          {
            // Near the front if any of the leaf may be: its centroid's distance
            // less its circumradius about the centroid. The centroid alone
            // misses the slivers next to the plane that the front runs through.
            const vec3 centroid = (leaf[0] + leaf[1] + leaf[2] + leaf[3]) / 4;
            real reach = 0;
            for (const vec3& corner : leaf) reach = std::max(reach, (corner - centroid).norm());
            const bool near_front = crack.distance_to_front(centroid) - reach < options.front_radius;
            refine(leaf, near_front ? options.front_levels : 0, refined);
          }
          leaves.swap(refined);
        }
        for (const tet_corners& leaf : leaves)
        {
          const real volume = std::abs(tet_volume(leaf[0], leaf[1], leaf[2], leaf[3]));
          const vec3 sum = leaf[0] + leaf[1] + leaf[2] + leaf[3];
          for (int i = 0; i < 4; ++i)
          {
            quadrature_point q;
            q.tet = static_cast<int>(e);
            q.side = side;
            q.weight = volume / 4;
            q.X = beta * sum + (alpha - beta) * leaf[static_cast<std::size_t>(i)];
            crack.field(q.X, side, q.h, q.grad_h);
            out.quadrature.push_back(q);
          }
        }
      }
    }
  }

  // The render data. The outer boundary: the face of each tet opposite its
  // node (wound outwards), subdivided, every small triangle clipped by the
  // plane into its two sides, with a vertex on the plane written once per side.
  render_builder render(out);
  const int k = std::max(1, options.render_resolution);
  for (index e = 0; e < T; ++e)
  {
    const int tet = static_cast<int>(e);
    const vec3 A = mesh.vertex(mesh.corner(e, 1));
    const vec3 B = mesh.vertex(mesh.corner(e, 2));
    const vec3 C = mesh.vertex(mesh.corner(e, 3));
    subdivide(A, B, C, k, [&](const vec3& p, const vec3& q, const vec3& r) {
      for (const int side : {1, -1})
      {
        const std::vector<vec3> piece = clip({p, q, r}, [&](const vec3& X) { return -side * crack.phi(X); });
        if (piece.size() < 3) continue;
        std::vector<int> ids;
        for (const vec3& X : piece)
        {
          real h;
          vec3 grad;
          crack.field(X, side, h, grad);
          ids.push_back(render.vertex(tet, side, X, h));
        }
        for (std::size_t i = 1; i + 1 < ids.size(); ++i) render.triangle(ids[0], ids[i], ids[i + 1], 0);
      }
    });
  }

  // The crack in each cut or front tet, twice: the + sheet (h = 1), wound with
  // its normal to -n, out of the + material, and the - sheet (h = 0) to +n.
  const real spacing = std::sqrt(real(2)) / real(k);
  for (index e = 0; e < T; ++e)
  {
    std::vector<vec3> P = crack_in[static_cast<std::size_t>(e)];
    if (P.size() < 3) continue;
    vec3 normal = vec3::Zero();
    for (std::size_t i = 0; i < P.size(); ++i) normal += P[i].cross(P[(i + 1) % P.size()]);
    if (normal.dot(crack.n) < 0) std::reverse(P.begin(), P.end());  // now CCW about n
    const int tet = static_cast<int>(e);
    for (std::size_t i = 1; i + 1 < P.size(); ++i)
    {
      const vec3 &A = P[0], &B = P[i], &C = P[i + 1];
      const real longest = std::max({(B - A).norm(), (C - B).norm(), (A - C).norm()});
      const int m = std::max(1, static_cast<int>(std::ceil(longest / spacing)));
      subdivide(A, B, C, m, [&](const vec3& p, const vec3& q, const vec3& r) {
        // CCW about n is outward for the - side's material; reversed for +.
        render.triangle(render.vertex(tet, 1, p, 1), render.vertex(tet, 1, r, 1), render.vertex(tet, 1, q, 1), 1);
        render.triangle(render.vertex(tet, -1, p, 0), render.vertex(tet, -1, q, 0), render.vertex(tet, -1, r, 0), 2);
      });
    }
  }

  char settings[256];
  std::snprintf(settings, sizeof settings,
                "synthetic, not partial_cut's solution: %s; quad level %d (+%d within %.2f of the front), "
                "%zu points; render resolution %d",
                options.through ? "through cut, h = the exact Heaviside"
                                : "partial cut, h = theta/(2 pi) + 1/2 (the 2D front solution)",
                options.quad_level, options.through ? 0 : options.front_levels, options.front_radius,
                out.quadrature.size(), options.render_resolution);
  out.settings = settings;
  return out;
}

}  // namespace xfem
