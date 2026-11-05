#include "xfem/xfem_enrichment.hpp"

#include "xfem/xfem_cutter.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace xfem {

const char* name(tet_class c) {
  switch (c)
  {
    case tet_class::missed: return "missed";
    case tet_class::cut: return "cut";
    case tet_class::front: return "front";
  }
  return "?";
}

const char* name(node_type t) {
  switch (t)
  {
    case node_type::none: return "none";
    case node_type::heaviside: return "heaviside";
    case node_type::harmonic: return "harmonic";
  }
  return "?";
}

const char* name(element_path p) {
  switch (p)
  {
    case element_path::plain: return "plain";
    case element_path::two_piece: return "two-piece";
    case element_path::quadrature: return "quadrature";
  }
  return "?";
}

std::string xfem_enrichment::summary() const {
  index types[3] = {0, 0, 0};
  for (const node_type t : type) ++types[static_cast<int>(t)];
  index tets[3] = {0, 0, 0};
  for (const tet_class c : classes) ++tets[static_cast<int>(c)];
  index paths[3] = {0, 0, 0};
  for (const element_path p : path) ++paths[static_cast<int>(p)];

  std::ostringstream out;
  out << "nodes: " << types[2] << " harmonic, " << types[1] << " heaviside, " << types[0]
      << " not enriched\n"
      << "tets: " << tets[0] << " missed, " << tets[1] << " cut, " << tets[2] << " front\n"
      << "paths: " << paths[0] << " plain, " << paths[1] << " two-piece, " << paths[2]
      << " quadrature";
  return out.str();
}

// One cut, start to finish, in the order motus's xfem_cutter::cut does it, with
// the node types of xfem_enrichment.hpp in place of "every corner of a cut tet
// is enriched":
//
//   1. phi from the plane, then snap: nothing downstream sees a zero.
//   2. the tet classes, checked against phi: a class-1 or class-2 tet must be
//      crossed by the plane.
//   3. the node types, per node, from the classes of the tets around it.
//   4. the dense numbering, a prefix sum over the enriched flags.
//   5. the element paths, per tet, from the classes and the corner types.
//   6. per two-piece tet, the piece volumes and the eight shape integrals.
xfem_enrichment build_enrichment(const tet_mesh& mesh, const vec3& n, real d,
                                 const std::vector<int>& classes, const std::vector<real>& node_h,
                                 const enrichment_options& options) {
  return build_enrichment(mesh, sample_plane(mesh, n, d), classes, node_h, options);
}

xfem_enrichment build_enrichment(const tet_mesh& mesh, std::vector<real> phi,
                                 const std::vector<int>& classes, const std::vector<real>& node_h,
                                 const enrichment_options& options) {
  const index vertices = mesh.vertex_count();
  if (static_cast<index>(phi.size()) != vertices)
  {
    throw std::invalid_argument("xfem: build_enrichment: " + std::to_string(phi.size()) + " phi values for " +
                                std::to_string(vertices) + " vertices");
  }
  const index tets = mesh.tet_count();
  if (static_cast<index>(classes.size()) != tets)
  {
    throw std::invalid_argument("xfem: build_enrichment: " + std::to_string(classes.size()) +
                                " tet classes for " + std::to_string(tets) + " tets");
  }
  if (static_cast<index>(node_h.size()) != vertices)
  {
    throw std::invalid_argument("xfem: build_enrichment: " + std::to_string(node_h.size()) +
                                " node values for " + std::to_string(vertices) + " vertices");
  }
  if (!(options.snap_tolerance > real(0)))
  {
    throw std::invalid_argument("xfem: build_enrichment: snap_tolerance must be positive; it is "
                                "what rules out a node exactly on the plane");
  }

  xfem_enrichment out;

  // 1. phi, snapped.
  out.phi = std::move(phi);
  out.snapped_nodes = snap(out.phi, mesh, options.snap_tolerance);
  out.node_h = node_h;

  // 2. The classes, against the plane.
  out.classes.resize(static_cast<std::size_t>(tets));
  for (index e = 0; e < tets; ++e)
  {
    const int c = classes[static_cast<std::size_t>(e)];
    if (c < 0 || c > 2)
    {
      throw std::invalid_argument("xfem: build_enrichment: tet " + std::to_string(e) + " has class " +
                                  std::to_string(c) + "; a class is 0, 1 or 2");
    }
    out.classes[static_cast<std::size_t>(e)] = static_cast<tet_class>(c);

    int positive = 0;
    for (int a = 0; a < 4; ++a)
    {
      if (out.phi[static_cast<std::size_t>(mesh.corner(e, a))] > real(0)) ++positive;
    }
    const bool crossed = positive != 0 && positive != 4;
    if (c != 0 && !crossed)
    {
      throw std::invalid_argument("xfem: build_enrichment: tet " + std::to_string(e) + " is class " +
                                  std::to_string(c) + " but the crack plane does not cross it");
    }
  }

  // 3. The node types. per_node is the rule of the header; all_harmonic types
  // every node per_node would enrich as harmonic; none enriches nothing.
  out.type.assign(static_cast<std::size_t>(vertices), node_type::none);
  if (options.rule != enrichment_rule::none)
  {
    for (index e = 0; e < tets; ++e)
    {
      const tet_class c = out.classes[static_cast<std::size_t>(e)];
      if (c == tet_class::missed) continue;
      for (int a = 0; a < 4; ++a)
      {
        node_type& t = out.type[static_cast<std::size_t>(mesh.corner(e, a))];
        if (c == tet_class::front) t = node_type::harmonic;
        else if (t == node_type::none) t = node_type::heaviside;
      }
    }
    if (options.rule == enrichment_rule::all_harmonic)
    {
      for (node_type& t : out.type)
      {
        if (t != node_type::none) t = node_type::harmonic;
      }
    }
  }

  // 4. The dense numbering, and its inverse alongside.
  out.enrichment.assign(static_cast<std::size_t>(vertices), -1);
  for (index i = 0; i < vertices; ++i)
  {
    if (out.type[static_cast<std::size_t>(i)] == node_type::none) continue;
    out.enrichment[static_cast<std::size_t>(i)] = static_cast<int>(out.enrichment_node.size());
    out.enrichment_node.push_back(static_cast<int>(i));
  }

  // 5. The paths.
  out.path.assign(static_cast<std::size_t>(tets), element_path::plain);
  out.cut_of.assign(static_cast<std::size_t>(tets), -1);
  for (index e = 0; e < tets; ++e)
  {
    const tet_class c = out.classes[static_cast<std::size_t>(e)];
    bool all_heaviside = true;
    bool any_active = false;
    for (int a = 0; a < 4; ++a)
    {
      const index node = mesh.corner(e, a);
      if (out.type_of(node) != node_type::heaviside) all_heaviside = false;
      if (out.active(node, c)) any_active = true;
    }

    element_path& p = out.path[static_cast<std::size_t>(e)];
    if (c == tet_class::cut && all_heaviside) p = element_path::two_piece;
    else if (any_active) p = element_path::quadrature;
    else p = element_path::plain;

    if (p == element_path::two_piece)
    {
      out.cut_of[static_cast<std::size_t>(e)] = static_cast<int>(out.cut_element.size());
      out.cut_element.push_back(static_cast<int>(e));
    }
  }

  // 6. Per two-piece tet, motus's steps 4 and 5: one tessellation gives both
  // piece volumes and the eight shape integrals.
  out.cut_volume.assign(static_cast<std::size_t>(2 * out.cut_count()), real(0));
  out.cut_shape_integral.assign(static_cast<std::size_t>(8 * out.cut_count()), real(0));
  for (index k = 0; k < out.cut_count(); ++k)
  {
    const index e = out.cut_element[static_cast<std::size_t>(k)];
    vec4 phi;
    vec3 corner[4];
    for (int a = 0; a < 4; ++a)
    {
      const index node = mesh.corner(e, a);
      phi[a] = out.phi[static_cast<std::size_t>(node)];
      corner[a] = mesh.vertex(node);
    }
    const real volume = std::abs(tet_volume(corner[0], corner[1], corner[2], corner[3]));
    const cut_measure m = measure_cut(phi, volume);

    // Both pieces must have volume: a piece with none would hand its far
    // corners a massless enrichment, which motus's snapping exists to prevent.
    if (!(m.volume[0] > real(0)) || !(m.volume[1] > real(0)))
    {
      throw std::invalid_argument("xfem: build_enrichment: two-piece tet " + std::to_string(e) +
                                  " has a piece with no volume");
    }
    for (int s = 0; s < 2; ++s) out.cut_volume[static_cast<std::size_t>(2 * k + s)] = m.volume[s];
    for (int a = 0; a < 8; ++a) out.cut_shape_integral[static_cast<std::size_t>(8 * k + a)] = m.shape_integral[a];
  }

  return out;
}

}  // namespace xfem
