#include "live_mesh.hpp"

#include "xfem_mesh.h"

#include <stdexcept>
#include <string>
#include <utility>

namespace xfem {

mesh_enrichment live_mesh(const tet_mesh& mesh, const std::vector<crack_levels>& cracks, int quad_level) {
  std::vector<Vec3> verts(mesh.vertices.begin(), mesh.vertices.end());
  std::vector<XfemMeshCrack> levels;
  for (const crack_levels& c : cracks)
  {
    XfemMeshCrack k;
    k.name = c.name;
    for (const vec3& X : mesh.vertices)
    {
      k.phi.push_back(c.phi(X));
      k.psi.push_back(c.psi(X));
    }
    levels.push_back(std::move(k));
  }
  XfemMeshParams p;
  p.quadLevel = quad_level;
  XfemMeshData in;
  if (!buildXfemMesh(verts, mesh.tets, levels, p, in)) throw std::runtime_error("partial_cut's self-checks failed on the mesh");

  mesh_enrichment out;
  out.tet_class = in.tetClass;
  out.tet_crack = in.tetCrack;
  out.node_phi = in.nodePhi;
  out.node_value = in.nodeValue;
  out.quadrature.reserve(in.quad.size());
  for (const XfemMeshData::Quad& q : in.quad) out.quadrature.push_back({q.tet, q.side, q.w, q.X, q.h, q.g});
  for (const XfemMeshData::Front& f : in.front) out.front.push_back({f.tet, f.a, f.b});
  out.render_vertices.reserve(in.renderV.size());
  for (const XfemMeshData::RenderVertex& r : in.renderV) out.render_vertices.push_back({r.tet, r.side, r.X, r.h});
  out.render_triangles.reserve(in.renderT.size());
  for (const auto& t : in.renderT) out.render_triangles.push_back({{t[0], t[1], t[2]}, t[3]});
  out.min_phi = in.minPhi;
  out.settings = in.settings;
  return out;
}

}  // namespace xfem
