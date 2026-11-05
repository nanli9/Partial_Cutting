#pragma once

/// @file
/// The scalar, the index, the fixed-size matrices, and the flattening
/// convention every material and every element follows.
///
/// From motus/core/types/types.hpp and motus/core/math/tensor.hpp, without the
/// device annotations: everything in this project runs on the cpu.
///
/// The convention: a 3x3 matrix flattens to a 9-vector in **column-major**
/// order, so entry (row, column) lives at index `3 * column + row`. A
/// fourth-order tensor is then a 9x9 whose (i, j) entry is
/// d(vec(a)_i) / d(vec(b)_j). Getting this wrong is silent -- the energy and the
/// gradient stay right and only Newton's convergence rate degrades -- so the
/// finite-difference self-tests are what pin it.

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cstddef>

namespace xfem {

using real = double;
using index = std::ptrdiff_t;

template <int rows, int columns>
using matrix = Eigen::Matrix<real, rows, columns>;
template <int rows>
using vector = Eigen::Matrix<real, rows, 1>;

using vec3 = vector<3>;
using vec4 = vector<4>;
using vec9 = vector<9>;
using mat3 = matrix<3, 3>;
using mat4 = matrix<4, 4>;
using mat9 = matrix<9, 9>;
using vector_x = Eigen::Matrix<real, Eigen::Dynamic, 1>;
using matrix_x = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>;

// flatten() reads a matrix's storage as it lies, which is the documented order
// only because Eigen stores column-major by default. Pin it.
static_assert((static_cast<int>(mat3::Options) & Eigen::RowMajor) == 0,
              "the flattening convention needs column-major mat3 storage");

/// Column-major flattening of a 3x3 into a 9-vector.
inline vec9 flatten(const mat3& m) { return Eigen::Map<const vec9>(m.data()); }

/// The inverse of flatten().
inline mat3 unflatten(const vec9& v) { return Eigen::Map<const mat3>(v.data()); }

/// Project a symmetric 9x9 onto the positive semi-definite cone by clamping its
/// negative eigenvalues to zero. This is what keeps Newton's assembled hessian
/// factorisable when a point is compressed. A matrix that is already positive
/// definite -- most points, most of the time -- is detected by a Cholesky that
/// costs a small fraction of the eigendecomposition, and returned as it came.
/// From motus/core/math/tensor.hpp.
inline mat9 project_to_spd(const mat9& h) {
  const Eigen::LLT<mat9> cholesky(h);
  if (cholesky.info() == Eigen::Success) return h;

  const Eigen::SelfAdjointEigenSolver<mat9> eigen(h);
  const vec9 clamped = eigen.eigenvalues().cwiseMax(real(0));
  return eigen.eigenvectors() * clamped.asDiagonal() * eigen.eigenvectors().transpose();
}

}  // namespace xfem
