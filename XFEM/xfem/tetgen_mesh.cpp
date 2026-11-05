#include "xfem/tetgen_mesh.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace xfem {

namespace {

// The next line that holds data: comments cut off, blank lines skipped.
bool next_data_line(std::ifstream& in, std::string& line, int& number) {
  while (std::getline(in, line))
  {
    ++number;
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    if (line.find_first_not_of(" \t\r") != std::string::npos) return true;
  }
  return false;
}

[[noreturn]] void fail(const std::string& path, int line, const std::string& what) {
  throw std::runtime_error("xfem: read_tetgen_mesh: " + path + ", line " + std::to_string(line) + ": " + what);
}

}  // namespace

tet_mesh read_tetgen_mesh(const std::string& prefix) {
  tet_mesh mesh;
  int first = 0;

  {
    const std::string path = prefix + ".node";
    std::ifstream in(path);
    if (!in) throw std::runtime_error("xfem: read_tetgen_mesh: cannot open " + path);
    std::string line;
    int number = 0;
    if (!next_data_line(in, line, number)) fail(path, number, "no header");
    std::istringstream header(line);
    int count = 0, dimension = 0;
    if (!(header >> count >> dimension) || count < 0 || dimension != 3) fail(path, number, "expected '<nodes> 3 ...'");
    mesh.vertices.resize(static_cast<std::size_t>(count));
    for (int v = 0; v < count; ++v)
    {
      if (!next_data_line(in, line, number)) fail(path, number, "ends after " + std::to_string(v) + " nodes");
      std::istringstream values(line);
      int id = 0;
      vec3& X = mesh.vertices[static_cast<std::size_t>(v)];
      if (!(values >> id >> X[0] >> X[1] >> X[2])) fail(path, number, "expected 'id x y z'");
      if (v == 0) first = id;
      if (id != first + v) fail(path, number, "node ids are not consecutive");
    }
  }

  const std::string path = prefix + ".ele";
  std::ifstream in(path);
  if (!in) throw std::runtime_error("xfem: read_tetgen_mesh: cannot open " + path);
  std::string line;
  int number = 0;
  if (!next_data_line(in, line, number)) fail(path, number, "no header");
  std::istringstream header(line);
  int count = 0, corners = 0;
  if (!(header >> count >> corners) || count < 0 || corners != 4) fail(path, number, "expected '<tets> 4 ...'");
  mesh.tets.resize(static_cast<std::size_t>(count));
  for (int t = 0; t < count; ++t)
  {
    if (!next_data_line(in, line, number)) fail(path, number, "ends after " + std::to_string(t) + " tets");
    std::istringstream values(line);
    int id = 0;
    if (!(values >> id)) fail(path, number, "expected 'id a b c d'");
    for (int a = 0; a < 4; ++a)
    {
      int node = 0;
      if (!(values >> node)) fail(path, number, "expected 'id a b c d'");
      node -= first;
      if (node < 0 || node >= mesh.vertex_count()) fail(path, number, "a corner that is not a node");
      mesh.tets[static_cast<std::size_t>(t)][static_cast<std::size_t>(a)] = node;
    }
    vec3 X[4];
    mesh.corners(t, X);
    if (!(tet_volume(X[0], X[1], X[2], X[3]) > real(0))) fail(path, number, "a tet that is not positively oriented");
  }
  return mesh;
}

}  // namespace xfem
