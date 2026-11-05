#include "live_ring.hpp"

#include "coupling_gate.h"

#include <stdexcept>
#include <string>

namespace xfem {

ring_export live_ring(const live_ring_options& options) {
  XfemExportParams p;
  p.ring = options.ring;
  p.crack = options.crack;
  p.quadLevel = options.quad_level;
  XfemRingData in;
  if (!buildXfemRing(p, in)) throw std::runtime_error("partial_cut's self-checks failed for " + in.caseName);

  // Field for field what runXfemExport writes and read_ring_export reads.
  ring_export out;
  out.version = 1;
  out.case_name = in.caseName;
  out.ring = in.ring;
  out.g = in.g;
  out.settings = in.settings;
  out.mesh.vertices.assign(in.verts.begin(), in.verts.end());
  out.mesh.tets = in.tets;
  out.plane_normal = in.n;
  out.plane_offset = in.d;
  out.front.assign(in.front.begin(), in.front.end());
  out.tet_class = in.tetClass;
  out.node_value = in.nodeValue;
  out.quadrature.reserve(in.quad.size());
  for (const XfemRingData::Quad& q : in.quad) out.quadrature.push_back({q.tet, q.side, q.w, q.X, q.h, q.g});
  out.render_vertices.reserve(in.renderV.size());
  for (const XfemRingData::RenderVertex& r : in.renderV) out.render_vertices.push_back({r.tet, r.side, r.X, r.h});
  out.render_triangles.reserve(in.renderT.size());
  for (const auto& t : in.renderT) out.render_triangles.push_back({{t[0], t[1], t[2]}, t[3]});
  return out;
}

}  // namespace xfem
