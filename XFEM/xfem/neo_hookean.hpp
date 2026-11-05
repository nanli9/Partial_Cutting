#pragma once

/// @file
/// Neo-Hookean elasticity.
///
/// From motus/materials/elastic/neo_hookean.hpp: the same energy, first Piola
/// stress and 9x9 hessian, as a plain value with no adapter around it.
///
///     psi(F) = mu / 2 (tr(F^T F) - 3) - mu ln J + lambda / 2 (ln J)^2
///
/// The volumetric term goes to infinity as det F approaches zero, so a point
/// resists inversion instead of passing through it -- which is also why the
/// line search has to keep det F > 0 at every quadrature point and piece
/// (xfem_system::maximum_non_inverting_step). At det F <= 0 the energy here is
/// +infinity rather than motus's nan, so an Armijo test rejects such a trial
/// outright; the step limit means it should never be reached.

#include "xfem/types.hpp"

#include <Eigen/LU>

#include <cmath>
#include <limits>

namespace xfem {

struct neo_hookean {
  real mu = 0.0;      ///< shear modulus, Lame's second parameter
  real lambda = 0.0;  ///< Lame's first parameter

  /// Lame parameters from Young's modulus and Poisson's ratio.
  static neo_hookean from(real youngs_modulus, real poissons_ratio) {
    neo_hookean nh;
    const real E = youngs_modulus;
    const real nu = poissons_ratio;
    nh.mu = E / (2 * (1 + nu));
    nh.lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    return nh;
  }

  /// Strain energy per unit reference volume.
  real energy_density(const mat3& f) const {
    const real J = f.determinant();
    if (!(J > real(0))) return std::numeric_limits<real>::infinity();
    const real lnJ = std::log(J);
    return (mu / 2) * ((f.transpose() * f).trace() - 3) - mu * lnJ + (lambda / 2) * lnJ * lnJ;
  }

  /// The first Piola-Kirchhoff stress, P = d(psi) / d(F).
  mat3 first_piola(const mat3& f) const {
    return mu * f + (lambda * std::log(f.determinant()) - mu) * f.inverse().transpose();
  }

  /// d(P) / d(F) as a 9x9 in the column-major flattening of types.hpp,
  /// projected to positive semi-definite when asked.
  mat9 hessian(const mat3& f, bool project_spd) const {
    const real J = f.determinant();
    const real lnJ = std::log(J);

    const mat3 A = f.inverse().transpose();  // F^{-T}
    const real c = lambda * lnJ - mu;

    mat9 H;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        const int r = 3 * j + i;
        for (int k = 0; k < 3; ++k)
        {
          for (int l = 0; l < 3; ++l)
          {
            const int s = 3 * l + k;
            const real identity = (i == k && j == l) ? mu : real(0);
            const real volumetric = lambda * A(i, j) * A(k, l);
            const real inverse_term = -c * A(i, l) * A(k, j);
            H(r, s) = identity + volumetric + inverse_term;
          }
        }
      }
    }

    // Under compression this H genuinely has negative eigenvalues, and then
    // -H^-1 g points uphill and Cholesky refuses to factorise. Clamping them to
    // zero per point keeps the assembled matrix spd, because a sum of psd
    // blocks is psd -- which is why this is done here and never globally. Pass
    // false to finite-difference the real thing.
    if (project_spd) H = project_to_spd(H);
    return H;
  }
};

}  // namespace xfem
