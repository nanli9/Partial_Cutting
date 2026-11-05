#pragma once

/// @file
/// Backward euler as the minimisation it is, solved by Newton on the cpu.
///
/// From motus/integrators/implicit/backward_euler.{hpp,cpp},
/// motus/solvers/nonlinear/newton.{hpp,cpp} and
/// motus/solvers/nonlinear/line_search.{hpp,cpp}, folded into one small file:
/// no objective interface, no injected linear solver, no device. One step
/// minimises the incremental potential
///
///     E(x) = 1 / (2 dt^2) (x - x~)^T M (x - x~)  +  Psi(x)  -  f_ext . x,
///     x~   = x_n + dt v_n,
///
/// by Newton: solve (M / dt^2 + K) d = -g with a sparse Cholesky (LDLT), take
/// the largest step the inversion filter allows (det F > 0 at every point and
/// piece: xfem_system::maximum_non_inverting_step), backtrack from there until
/// Armijo holds, repeat. Then v = (x - x_n) / dt. No gravity and no VBD: the
/// external force is all the load there is.
///
/// With jump constraints on the system (the front's, xfem_system::zero_jump_at)
/// Newton runs in their null space: d = Z y with (Z^T (M / dt^2 + K) Z) y =
/// -Z^T g, and the convergence test reads Z^T g, the force along the allowed
/// directions. x_n, v_n and so the prediction already satisfy the constraints,
/// and so does every step.

#include "xfem/types.hpp"
#include "xfem/xfem_system.hpp"

#include <Eigen/SparseCore>

#include <string>
#include <vector>

namespace xfem {

struct mechanical_state {
  vector_x x;  ///< positions (node groups) and jumps (enrichment groups)
  vector_x v;
};

struct newton_options {
  index max_iterations = 50;
  /// Converged when |g|_inf falls below this. An absolute tolerance on a force
  /// has units, so the caller scales it to the problem (the demo uses a
  /// multiple of E h^2).
  real gradient_tolerance = real(1e-6);
  /// Converged when a step moves the solution less than this.
  real step_tolerance = real(1e-12);
  real armijo = real(1e-4);
  real shrink = real(0.5);
  index max_line_search_steps = 40;
  bool project_spd = true;
};

struct newton_result {
  index iterations = 0;
  real final_energy = 0;
  real final_gradient_norm = 0;
  bool converged = false;
  /// |g|_inf at the top of every iteration, the last one included.
  std::vector<real> residuals;
  /// Why it stopped when it did not converge; empty otherwise.
  std::string failure;
};

/// The incremental potential of one step, against a system, a predicted
/// position, an external force (zero on the fixed dofs) and dt.
class incremental_potential {
 public:
  incremental_potential(const xfem_system& system, const vector_x& predicted, const vector_x& external_force,
                        real dt);

  index dofs() const { return system_->dofs(); }
  real value(const vector_x& x) const;
  /// Zero on the fixed dofs.
  void gradient(const vector_x& x, vector_x& gradient) const;
  void hessian(const vector_x& x, bool project_spd, Eigen::SparseMatrix<real>& hessian) const;
  real step_limit(const vector_x& x, const vector_x& direction) const {
    return system_->maximum_non_inverting_step(x, direction);
  }
  /// The system's constraint basis Z (xfem_system::zero_jump_at), or nullptr:
  /// with one, Newton steps in Z's span only.
  const Eigen::SparseMatrix<real>* basis() const { return system_->constraint_basis(); }

 private:
  const xfem_system* system_;
  const vector_x* predicted_;
  const vector_x* external_force_;
  real dt_;
};

/// Minimise from x, overwriting it with the answer (the best iterate on a
/// failure; result.converged says which).
newton_result newton_solve(const incremental_potential& problem, vector_x& x, const newton_options& options);

class backward_euler {
 public:
  explicit backward_euler(newton_options options = {}) : options_(options) {}

  /// One implicit step. The solve starts from x_n (motus's
  /// initial_guess::previous, always a valid configuration); the external
  /// force is copied with every fixed dof zeroed, so a pinned dof feels no
  /// load. A solve that did not converge is still written back -- it is the
  /// best iterate found, and the result says so.
  newton_result step(const xfem_system& system, mechanical_state& state, real dt, const vector_x& external_force);

  newton_options& options() { return options_; }
  const newton_options& options() const { return options_; }

 private:
  newton_options options_;
};

}  // namespace xfem
