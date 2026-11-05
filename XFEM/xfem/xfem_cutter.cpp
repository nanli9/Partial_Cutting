#include "xfem/xfem_cutter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

// The per-element geometry of motus/topology/cutting/xfem_cutter.cpp, with its
// two invariants worth knowing before trusting a number out of it:
//
//     positive_volume_fraction  in [0, 1]  always
//     f(phi) + f(-phi)          == 1       for every phi
//
// The self-tests check both on a 3-1, a 2-2 and an uncut tet.

namespace xfem {

namespace {

// Where the zero set crosses edge (i, j): phi is linear along the edge, so it
// vanishes at t = phi_i / (phi_i - phi_j), and that point is (1 - t) e_i +
// t e_j. Only the two endpoint values enter.
vec4 edge_crossing(const vec4& phi, int i, int j) {
  const real t = phi[i] / (phi[i] - phi[j]);
  vec4 c = vec4::Zero();
  c[i] = real(1) - t;
  c[j] = t;
  return c;
}

}  // namespace

// A sub-tet with its four corners in barycentric coordinates of the parent:
// its volume, as a fraction of the parent's, is |det| of the 3x3 of corner
// differences. The differences each sum to zero, so dropping any one of the
// four coordinates leaves the same magnitude; drop the first.
real barycentric_fraction(const vec4& a, const vec4& b, const vec4& c, const vec4& d) {
  mat3 edges;
  edges.col(0) = (b - a).tail<3>();
  edges.col(1) = (c - a).tail<3>();
  edges.col(2) = (d - a).tail<3>();
  return std::abs(edges.determinant());
}

// The one tessellation: positive_volume_fraction sums its determinants,
// measure_cut integrates the shape functions over it, and the synthetic ring
// refines it into its quadrature.
int positive_piece(const vec4& phi, vec4 (&tets)[3][4]) {
  int positive = 0;
  for (int a = 0; a < 4; ++a)
  {
    if (phi[a] > real(0)) ++positive;
  }

  switch (positive)
  {
    case 0: return 0;

    case 4:
      for (int a = 0; a < 4; ++a) tets[0][a] = vec4::Unit(a);
      return 1;

    case 1:
    case 3: {
      // 3-1. The lone corner L is the one whose sign the other three do not
      // share; i, j, k are the other three.
      const bool lone_is_positive = (positive == 1);
      int L = 0;
      for (int a = 0; a < 4; ++a)
      {
        if ((phi[a] > real(0)) == lone_is_positive) L = a;
      }
      const int i = (L + 1) % 4;
      const int j = (L + 2) % 4;
      const int k = (L + 3) % 4;

      const vec4 ci = edge_crossing(phi, L, i);
      const vec4 cj = edge_crossing(phi, L, j);
      const vec4 ck = edge_crossing(phi, L, k);

      if (lone_is_positive)
      {
        // The cut chops the corner off at L: a little tet, L and its three
        // crossings.
        tets[0][0] = vec4::Unit(L);
        tets[0][1] = ci;
        tets[0][2] = cj;
        tets[0][3] = ck;
        return 1;
      }

      // What is left is a wedge, the triangle (ci, cj, ck) over the base
      // (ei, ej, ek), and a wedge is three tets.
      const vec4 ei = vec4::Unit(i);
      const vec4 ej = vec4::Unit(j);
      const vec4 ek = vec4::Unit(k);
      tets[0][0] = ei; tets[0][1] = ej; tets[0][2] = ek; tets[0][3] = ci;
      tets[1][0] = ej; tets[1][1] = ek; tets[1][2] = ci; tets[1][3] = cj;
      tets[2][0] = ek; tets[2][1] = ci; tets[2][2] = cj; tets[2][3] = ck;
      return 3;
    }

    default: break;
  }

  // 2-2, the count left. The positive pair (p0, p1) against the negative pair
  // (n0, n1); the zero set crosses the four edges joining a positive corner to
  // a negative one, and c_ab is the crossing on (p_a, n_b).
  int p[2] = {0, 0};
  int n[2] = {0, 0};
  int positives = 0;
  int negatives = 0;
  for (int a = 0; a < 4; ++a)
  {
    if (phi[a] > real(0)) p[positives++] = a;
    else n[negatives++] = a;
  }

  const vec4 c00 = edge_crossing(phi, p[0], n[0]);
  const vec4 c01 = edge_crossing(phi, p[0], n[1]);
  const vec4 c10 = edge_crossing(phi, p[1], n[0]);
  const vec4 c11 = edge_crossing(phi, p[1], n[1]);

  // The positive piece is a wedge: the triangle (p0, c00, c01) at one end,
  // (p1, c10, c11) at the other, corner over corner. Three tets, split the
  // same way as the 3-1 base. The negative piece is the same wedge with the
  // pairs swapped, and it is what f(-phi) measures.
  const vec4 ep0 = vec4::Unit(p[0]);
  const vec4 ep1 = vec4::Unit(p[1]);
  tets[0][0] = ep0; tets[0][1] = c00; tets[0][2] = c01; tets[0][3] = ep1;
  tets[1][0] = c00; tets[1][1] = c01; tets[1][2] = ep1; tets[1][3] = c10;
  tets[2][0] = c01; tets[2][1] = ep1; tets[2][2] = c10; tets[2][3] = c11;
  return 3;
}

real positive_volume_fraction(const vec4& phi) {
  vec4 tets[3][4];
  const int count = positive_piece(phi, tets);

  real fraction = real(0);
  for (int s = 0; s < count; ++s)
  {
    fraction += barycentric_fraction(tets[s][0], tets[s][1], tets[s][2], tets[s][3]);
  }
  return fraction;
}

// 3-1 gives a triangle: the three crossings on the edges at the lone corner.
// 2-2 gives a quadrilateral: the four crossings on the edges joining the pairs,
// wound so that consecutive crossings share a corner and (0, 1, 2), (0, 2, 3)
// cannot cross.
int cut_polygon(const vec4& phi, vec4 (&out)[4]) {
  int positive = 0;
  for (int a = 0; a < 4; ++a)
  {
    if (phi[a] > real(0)) ++positive;
  }
  if (positive == 0 || positive == 4) return 0;

  if (positive == 1 || positive == 3)
  {
    const bool lone_is_positive = (positive == 1);
    int L = 0;
    for (int a = 0; a < 4; ++a)
    {
      if ((phi[a] > real(0)) == lone_is_positive) L = a;
    }
    for (int i = 0; i < 3; ++i) out[i] = edge_crossing(phi, L, (L + 1 + i) % 4);
    return 3;
  }

  int p[2] = {0, 0};
  int n[2] = {0, 0};
  int positives = 0;
  int negatives = 0;
  for (int a = 0; a < 4; ++a)
  {
    if (phi[a] > real(0)) p[positives++] = a;
    else n[negatives++] = a;
  }
  out[0] = edge_crossing(phi, p[0], n[0]);
  out[1] = edge_crossing(phi, p[0], n[1]);
  out[2] = edge_crossing(phi, p[1], n[1]);
  out[3] = edge_crossing(phi, p[1], n[0]);
  return 4;
}

// Steps 4 and 5 of motus's xfem_cutter::cut, for one tet. The positive piece's
// sub-tets give its volume and its four shape integrals at once: N_a is linear,
// so over a sub-tet its integral is the sub-tet's volume times N_a at the
// centroid -- the midpoint rule, exact for a linear integrand -- and N_a at a
// barycentric point is that point's a-th coordinate. The negative side is the
// remainder of each, so the two sum back exactly.
cut_measure measure_cut(const vec4& phi, real volume) {
  vec4 sub[3][4];
  const int count = positive_piece(phi, sub);

  real plus_fraction = real(0);
  vec4 plus_integral = vec4::Zero();
  for (int s = 0; s < count; ++s)
  {
    const real fraction = barycentric_fraction(sub[s][0], sub[s][1], sub[s][2], sub[s][3]);
    const vec4 centroid = (sub[s][0] + sub[s][1] + sub[s][2] + sub[s][3]) / real(4);
    plus_fraction += fraction;
    plus_integral += fraction * centroid;
  }

  cut_measure out;
  out.volume[0] = plus_fraction * volume;
  out.volume[1] = volume - out.volume[0];
  for (int a = 0; a < 4; ++a)
  {
    out.shape_integral[a] = plus_integral[a] * volume;
    out.shape_integral[4 + a] = volume / real(4) - out.shape_integral[a];
  }
  return out;
}

// Step 1 of motus's cut, for a plane: phi_i = p.signed_distance(X_i).
std::vector<real> sample_plane(const tet_mesh& mesh, const vec3& n, real d) {
  std::vector<real> phi(static_cast<std::size_t>(mesh.vertex_count()));
  for (index i = 0; i < mesh.vertex_count(); ++i) phi[static_cast<std::size_t>(i)] = n.dot(mesh.vertex(i)) - d;
  return phi;
}

// The shortest edge at each vertex (motus's build_vertex_scale), then motus's
// snap onto the positive side.
index snap(std::vector<real>& phi, const tet_mesh& mesh, real tolerance) {
  std::vector<real> scale(static_cast<std::size_t>(mesh.vertex_count()), std::numeric_limits<real>::max());
  for (index e = 0; e < mesh.tet_count(); ++e)
  {
    for (int a = 0; a < 4; ++a)
    {
      for (int b = a + 1; b < 4; ++b)
      {
        const int i = mesh.corner(e, a);
        const int j = mesh.corner(e, b);
        const real length = (mesh.vertex(i) - mesh.vertex(j)).norm();
        scale[static_cast<std::size_t>(i)] = std::min(scale[static_cast<std::size_t>(i)], length);
        scale[static_cast<std::size_t>(j)] = std::min(scale[static_cast<std::size_t>(j)], length);
      }
    }
  }

  index snapped = 0;
  for (std::size_t i = 0; i < phi.size(); ++i)
  {
    // A vertex no tet references keeps the sentinel; give it something finite.
    const real edge = scale[i] == std::numeric_limits<real>::max() ? real(1) : scale[i];
    const real bound = tolerance * edge;
    if (std::abs(phi[i]) < bound)
    {
      phi[i] = bound;
      ++snapped;
    }
  }
  return snapped;
}

}  // namespace xfem
