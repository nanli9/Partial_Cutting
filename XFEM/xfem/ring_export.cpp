#include "xfem/ring_export.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace xfem {

namespace {

// Whitespace-separated tokens with their line numbers, and the rest of a line
// for `settings`, whose free text may hold anything.
class tokens {
 public:
  tokens(std::string text, std::string path) : text_(std::move(text)), path_(std::move(path)) {}

  bool at_end() {
    skip_space();
    return pos_ >= text_.size();
  }

  std::string word() {
    skip_space();
    if (pos_ >= text_.size()) fail("unexpected end of file");
    const std::size_t start = pos_;
    while (pos_ < text_.size() && !std::isspace(static_cast<unsigned char>(text_[pos_]))) ++pos_;
    return text_.substr(start, pos_ - start);
  }

  std::string rest_of_line() {
    while (pos_ < text_.size() && (text_[pos_] == ' ' || text_[pos_] == '\t')) ++pos_;
    const std::size_t start = pos_;
    while (pos_ < text_.size() && text_[pos_] != '\n' && text_[pos_] != '\r') ++pos_;
    return text_.substr(start, pos_ - start);
  }

  real number() {
    const std::string w = word();
    char* end = nullptr;
    const double v = std::strtod(w.c_str(), &end);
    if (end == w.c_str() || *end != '\0') fail("expected a number, found '" + w + "'");
    return static_cast<real>(v);
  }

  long integer() {
    const std::string w = word();
    char* end = nullptr;
    const long v = std::strtol(w.c_str(), &end, 10);
    if (end == w.c_str() || *end != '\0') fail("expected an integer, found '" + w + "'");
    return v;
  }

  index count(const std::string& section) {
    const long c = integer();
    if (c < 0) fail(section + ": a negative count");
    return static_cast<index>(c);
  }

  vec3 point() {
    const real x = number();
    const real y = number();
    const real z = number();
    return vec3(x, y, z);
  }

  [[noreturn]] void fail(const std::string& what) const {
    throw std::runtime_error("xfem: " + path_ + ", line " + std::to_string(line()) + ": " + what);
  }

  index line() const {
    return 1 + static_cast<index>(std::count(text_.begin(), text_.begin() + static_cast<std::ptrdiff_t>(
                                                                std::min(pos_, text_.size())), '\n'));
  }

 private:
  void skip_space() {
    while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) ++pos_;
  }

  std::string text_;
  std::string path_;
  std::size_t pos_ = 0;
};

}  // namespace

ring_export read_ring_export(const std::string& path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) throw std::runtime_error("xfem: cannot open " + path);
  std::stringstream buffer;
  buffer << file.rdbuf();
  tokens in(buffer.str(), path);

  ring_export out;
  if (in.word() != "xfem_ring_export") in.fail("not an xfem_ring_export file");
  out.version = static_cast<int>(in.integer());
  if (out.version != 1) in.fail("version " + std::to_string(out.version) + " (this reader knows 1)");

  bool seen_vertices = false, seen_tets = false, seen_plane = false, seen_front = false;
  bool seen_class = false, seen_value = false, seen_quad = false, seen_end = false;
  while (!in.at_end())
  {
    const std::string key = in.word();
    if (key == "end")
    {
      seen_end = true;
      break;
    }
    if (key == "case") out.case_name = in.word();
    else if (key == "ring") out.ring = static_cast<int>(in.integer());
    else if (key == "g") out.g = in.number();
    else if (key == "settings") out.settings = in.rest_of_line();
    else if (key == "vertices")
    {
      const index n = in.count(key);
      out.mesh.vertices.resize(static_cast<std::size_t>(n));
      for (vec3& v : out.mesh.vertices) v = in.point();
      seen_vertices = true;
    }
    else if (key == "tets")
    {
      const index n = in.count(key);
      out.mesh.tets.resize(static_cast<std::size_t>(n));
      for (auto& t : out.mesh.tets)
      {
        for (int& c : t) c = static_cast<int>(in.integer());
      }
      seen_tets = true;
    }
    else if (key == "plane")
    {
      out.plane_normal = in.point();
      out.plane_offset = in.number();
      seen_plane = true;
    }
    else if (key == "front")
    {
      const index n = in.count(key);
      if (n != 0 && n != 2) in.fail("front: " + std::to_string(n) + " points (0 or 2)");
      out.front.resize(static_cast<std::size_t>(n));
      for (vec3& p : out.front) p = in.point();
      seen_front = true;
    }
    else if (key == "tet_class")
    {
      const index n = in.count(key);
      out.tet_class.resize(static_cast<std::size_t>(n));
      for (int& c : out.tet_class) c = static_cast<int>(in.integer());
      seen_class = true;
    }
    else if (key == "node_value")
    {
      const index n = in.count(key);
      out.node_value.resize(static_cast<std::size_t>(n));
      for (real& h : out.node_value) h = in.number();
      seen_value = true;
    }
    else if (key == "quad")
    {
      const index n = in.count(key);
      out.quadrature.resize(static_cast<std::size_t>(n));
      for (quadrature_point& q : out.quadrature)
      {
        q.tet = static_cast<int>(in.integer());
        q.side = static_cast<int>(in.integer());
        q.weight = in.number();
        q.X = in.point();
        q.h = in.number();
        q.grad_h = in.point();
      }
      seen_quad = true;
    }
    else if (key == "render_vertices")
    {
      const index n = in.count(key);
      out.render_vertices.resize(static_cast<std::size_t>(n));
      for (render_vertex& v : out.render_vertices)
      {
        v.tet = static_cast<int>(in.integer());
        v.side = static_cast<int>(in.integer());
        v.X = in.point();
        v.h = in.number();
      }
    }
    else if (key == "render_triangles")
    {
      const index n = in.count(key);
      out.render_triangles.resize(static_cast<std::size_t>(n));
      for (render_triangle& t : out.render_triangles)
      {
        for (int& c : t.v) c = static_cast<int>(in.integer());
        t.group = static_cast<int>(in.integer());
      }
    }
    else in.fail("unknown section '" + key + "'");
  }

  if (!seen_end) in.fail("no 'end'");
  const char* missing = !seen_vertices ? "vertices"
                        : !seen_tets   ? "tets"
                        : !seen_plane  ? "plane"
                        : !seen_front  ? "front"
                        : !seen_class  ? "tet_class"
                        : !seen_value  ? "node_value"
                        : !seen_quad   ? "quad"
                                       : nullptr;
  if (missing) throw std::runtime_error("xfem: " + path + ": no '" + missing + "' section");

  // Indices are checked here rather than in check(): an out-of-range one would
  // be read out of bounds by everything downstream, not merely be wrong.
  const index V = out.mesh.vertex_count();
  const index T = out.mesh.tet_count();
  const auto bad = [&](const std::string& what) {
    throw std::runtime_error("xfem: " + path + ": " + what);
  };
  for (const auto& t : out.mesh.tets)
  {
    for (const int c : t)
    {
      if (c < 0 || c >= V) bad("a tet names vertex " + std::to_string(c));
    }
  }
  if (static_cast<index>(out.tet_class.size()) != T) bad("tet_class is not one entry per tet");
  if (static_cast<index>(out.node_value.size()) != V) bad("node_value is not one entry per vertex");
  for (const quadrature_point& q : out.quadrature)
  {
    if (q.tet < 0 || q.tet >= T) bad("a quadrature point names tet " + std::to_string(q.tet));
  }
  for (const render_vertex& v : out.render_vertices)
  {
    if (v.tet < 0 || v.tet >= T) bad("a render vertex names tet " + std::to_string(v.tet));
  }
  const index R = static_cast<index>(out.render_vertices.size());
  for (const render_triangle& t : out.render_triangles)
  {
    for (const int c : t.v)
    {
      if (c < 0 || c >= R) bad("a render triangle names vertex " + std::to_string(c));
    }
    if (t.group < 0 || t.group > 2) bad("a render triangle has group " + std::to_string(t.group));
  }
  return out;
}

void write_ring_export(const std::string& path, const ring_export& data) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("xfem: cannot write " + path);
  out << std::setprecision(17);

  const auto p3 = [&](const vec3& v) { out << v[0] << ' ' << v[1] << ' ' << v[2]; };

  out << "xfem_ring_export " << data.version << '\n';
  out << "case " << data.case_name << '\n';
  out << "ring " << data.ring << '\n';
  out << "g " << data.g << '\n';
  out << "settings " << data.settings << '\n';
  out << "vertices " << data.mesh.vertex_count() << '\n';
  for (const vec3& v : data.mesh.vertices)
  {
    p3(v);
    out << '\n';
  }
  out << "tets " << data.mesh.tet_count() << '\n';
  for (const auto& t : data.mesh.tets) out << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3] << '\n';
  out << "plane ";
  p3(data.plane_normal);
  out << ' ' << data.plane_offset << '\n';
  out << "front " << data.front.size() << '\n';
  for (const vec3& p : data.front)
  {
    p3(p);
    out << '\n';
  }
  out << "tet_class " << data.tet_class.size() << '\n';
  for (const int c : data.tet_class) out << c << '\n';
  out << "node_value " << data.node_value.size() << '\n';
  for (const real h : data.node_value) out << h << '\n';
  out << "quad " << data.quadrature.size() << '\n';
  for (const quadrature_point& q : data.quadrature)
  {
    out << q.tet << ' ' << q.side << ' ' << q.weight << ' ';
    p3(q.X);
    out << ' ' << q.h << ' ';
    p3(q.grad_h);
    out << '\n';
  }
  out << "render_vertices " << data.render_vertices.size() << '\n';
  for (const render_vertex& v : data.render_vertices)
  {
    out << v.tet << ' ' << v.side << ' ';
    p3(v.X);
    out << ' ' << v.h << '\n';
  }
  out << "render_triangles " << data.render_triangles.size() << '\n';
  for (const render_triangle& t : data.render_triangles)
  {
    out << t.v[0] << ' ' << t.v[1] << ' ' << t.v[2] << ' ' << t.group << '\n';
  }
  out << "end\n";
  if (!out) throw std::runtime_error("xfem: writing " + path + " failed");
}

std::vector<std::string> check(const ring_export& data, real weight_tolerance) {
  std::vector<std::string> problems;
  const auto say = [&](const std::string& s) { problems.push_back(s); };
  char buf[256];

  const tet_mesh& mesh = data.mesh;
  const index T = mesh.tet_count();

  if (std::abs(data.plane_normal.norm() - real(1)) > real(1e-9)) say("plane: |n| is not 1");

  std::vector<real> volume(static_cast<std::size_t>(T), real(0));
  for (index e = 0; e < T; ++e)
  {
    vec3 c[4];
    mesh.corners(e, c);
    volume[static_cast<std::size_t>(e)] = tet_volume(c[0], c[1], c[2], c[3]);
    if (!(volume[static_cast<std::size_t>(e)] > real(0)))
    {
      std::snprintf(buf, sizeof buf, "tet %d is not positively oriented (volume %.3e)", int(e),
                    volume[static_cast<std::size_t>(e)]);
      say(buf);
    }
    const int k = data.tet_class[static_cast<std::size_t>(e)];
    if (k < 0 || k > 2) say("tet " + std::to_string(e) + " has class " + std::to_string(k));
  }

  // Weights per tet, and the side of every point against the plane.
  std::vector<real> sum(static_cast<std::size_t>(T), real(0));
  std::vector<index> per_tet(static_cast<std::size_t>(T), 0);
  index wrong_side = 0, nonpositive = 0;
  for (const quadrature_point& q : data.quadrature)
  {
    sum[static_cast<std::size_t>(q.tet)] += q.weight;
    ++per_tet[static_cast<std::size_t>(q.tet)];
    if (!(q.weight > real(0))) ++nonpositive;
    const real s = data.plane_normal.dot(q.X) - data.plane_offset;
    if ((q.side != 1 && q.side != -1) || (s > real(0) ? 1 : -1) != q.side) ++wrong_side;
  }
  if (nonpositive) say(std::to_string(nonpositive) + " quadrature points with weight <= 0");
  if (wrong_side) say(std::to_string(wrong_side) + " quadrature points whose side is not sign(n.X - d)");
  for (index e = 0; e < T; ++e)
  {
    const std::size_t u = static_cast<std::size_t>(e);
    if (per_tet[u] == 0)
    {
      say("tet " + std::to_string(e) + " has no quadrature points");
      continue;
    }
    const real relative = std::abs(sum[u] - volume[u]) / std::abs(volume[u]);
    if (relative > weight_tolerance)
    {
      std::snprintf(buf, sizeof buf, "tet %d: quadrature weights sum to %.15g, volume %.15g (rel %.2e)",
                    int(e), sum[u], volume[u], relative);
      say(buf);
    }
  }

  // Render triangles: all three vertices in one tet; sides +-1.
  for (std::size_t t = 0; t < data.render_triangles.size(); ++t)
  {
    const render_triangle& tri = data.render_triangles[t];
    const int tet = data.render_vertices[static_cast<std::size_t>(tri.v[0])].tet;
    for (const int v : tri.v)
    {
      if (data.render_vertices[static_cast<std::size_t>(v)].tet != tet)
      {
        say("render triangle " + std::to_string(t) + " spans more than one tet");
        break;
      }
    }
  }
  for (const render_vertex& v : data.render_vertices)
  {
    if (v.side != 1 && v.side != -1)
    {
      say("a render vertex has side " + std::to_string(v.side));
      break;
    }
  }
  return problems;
}

}  // namespace xfem
