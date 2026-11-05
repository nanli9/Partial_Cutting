#pragma once

/// @file
/// A tet mesh in rest space: vertices and four-corner tets, nothing else.
///
/// The host half of motus/geometry/mesh/tet_mesh.hpp, as plain std::vectors:
/// no buffers and no upload, since nothing here leaves the cpu. The mesh is
/// never edited by a cut -- that is what makes the method XFEM -- so it is
/// built once and read everywhere.

#include "xfem/types.hpp"

#include <array>
#include <vector>

namespace xfem {

struct tet_mesh {
  std::vector<vec3> vertices;
  /// Four vertex ids per tet, positively oriented: det of the three edges from
  /// corner 0 is positive. ring_export's reader and the synthetic ring both
  /// check it.
  std::vector<std::array<int, 4>> tets;

  index vertex_count() const { return static_cast<index>(vertices.size()); }
  index tet_count() const { return static_cast<index>(tets.size()); }

  int corner(index tet, int a) const {
    return tets[static_cast<std::size_t>(tet)][static_cast<std::size_t>(a)];
  }
  const vec3& vertex(index v) const { return vertices[static_cast<std::size_t>(v)]; }

  void corners(index tet, vec3 (&out)[4]) const {
    for (int a = 0; a < 4; ++a) out[a] = vertex(corner(tet, a));
  }
};

/// The signed volume of the tet (a, b, c, d): positive when (b - a, c - a,
/// d - a) is right-handed. From motus/geometry/mesh/tet_mesh.hpp.
inline real tet_volume(const vec3& a, const vec3& b, const vec3& c, const vec3& d) {
  return (b - a).dot((c - a).cross(d - a)) / real(6);
}

}  // namespace xfem
