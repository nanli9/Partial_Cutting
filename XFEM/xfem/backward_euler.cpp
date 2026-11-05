#include "xfem/backward_euler.hpp"

#include <Eigen/SparseCholesky>

#include <cmath>
#include <stdexcept>

namespace xfem {

incremental_potential::incremental_potential(const xfem_system& system, const vector_x& predicted,
                                             const vector_x& external_force, real dt)
    : system_(&system), predicted_(&predicted), external_force_(&external_force), dt_(dt) {
  if (!(dt > 0)) throw std::invalid_argument("xfem: incremental_potential: dt must be positive");
  if (predicted.size() != system.dofs() || external_force.size() != system.dofs())
  {
    throw std::invalid_argument("xfem: incremental_potential: predicted or external_force has the wrong size");
  }
}

// motus's step 1: the inertia term pulls x towards where the last velocity was
// heading, weighted by mass; the last term is the work of a load held constant
// over the step.
real incremental_potential::value(const vector_x& x) const {
  const vector_x d = x - *predicted_;
  return real(0.5) / (dt_ * dt_) * d.dot(system_->mass().cwiseProduct(d)) + system_->potential_energy(x) -
         external_force_->dot(x);
}

// motus's step 2, with its step 5: zero the fixed dofs again, since the inertia
// put a force back on a pinned dof whose x is not its predicted.
void incremental_potential::gradient(const vector_x& x, vector_x& gradient) const {
  system_->potential_gradient(x, gradient);
  gradient += system_->mass().cwiseProduct(x - *predicted_) / (dt_ * dt_) - *external_force_;
  for (const index d : system_->fixed_dofs()) gradient[d] = real(0);
}

// M / dt^2 + K. The mass is lumped, so it only adds to the diagonal, and a
// fixed dof's row and column stay clear (motus's pinned_hessian has nothing to
// do here).
void incremental_potential::hessian(const vector_x& x, bool project_spd, Eigen::SparseMatrix<real>& hessian) const {
  system_->potential_hessian(x, project_spd, real(1) / (dt_ * dt_), hessian);
}

namespace {

// motus's filtered_line_search: start below the inversion point instead of at
// 1, then halve until Armijo holds,
//
//     f(x + t d)  <=  f(x)  +  c t (g . d).
//
// Zero when it never does, or when d is not a descent direction.
real filtered_line_search(const incremental_potential& problem, const vector_x& x, const vector_x& direction,
                          real energy, const vector_x& gradient, const newton_options& options) {
  const real limit = problem.step_limit(x, direction);
  if (!(limit > real(0))) return real(0);

  const real slope = gradient.dot(direction);
  if (!(slope < real(0))) return real(0);

  real t = std::min(limit, real(1));
  vector_x candidate(x.size());
  for (index step = 0; step < options.max_line_search_steps; ++step)
  {
    candidate = x + t * direction;
    if (problem.value(candidate) <= energy + options.armijo * t * slope) return t;
    t *= options.shrink;
  }
  return real(0);
}

}  // namespace

// motus's newton::solve: at each x, solve h d = -g, shorten the step until the
// energy actually drops, take it, repeat. The symbolic analysis is done once:
// the pattern is the system's and does not move during a solve.
newton_result newton_solve(const incremental_potential& problem, vector_x& x, const newton_options& options) {
  newton_result result;
  vector_x g(problem.dofs());
  vector_x reduced;
  Eigen::SparseMatrix<real> h, h_reduced;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>> solver;
  bool analysed = false;
  // Constrained: every vector below lives in Z's coordinates until d = Z y.
  const Eigen::SparseMatrix<real>* Z = problem.basis();

  for (index i = 0;; ++i)
  {
    const real f = problem.value(x);
    problem.gradient(x, g);
    if (Z) reduced = Z->transpose() * g;
    const vector_x& r = Z ? reduced : g;
    result.iterations = i;
    result.final_energy = f;
    result.final_gradient_norm = r.lpNorm<Eigen::Infinity>();
    result.residuals.push_back(result.final_gradient_norm);

    // The energy is the honest nan detector (motus's note): an inverted point
    // makes it infinite, and the line search never accepts one.
    if (!std::isfinite(f) || !std::isfinite(result.final_gradient_norm))
    {
      result.failure = "the energy or the gradient is not finite";
      return result;
    }
    if (result.final_gradient_norm < options.gradient_tolerance)
    {
      result.converged = true;
      return result;
    }
    if (i >= options.max_iterations)
    {
      result.failure = "out of iterations";
      return result;
    }

    problem.hessian(x, options.project_spd, h);
    if (Z) h_reduced = Z->transpose() * h * *Z;
    const Eigen::SparseMatrix<real>& system_matrix = Z ? h_reduced : h;
    if (!analysed)
    {
      solver.analyzePattern(system_matrix);
      analysed = true;
    }
    solver.factorize(system_matrix);
    if (solver.info() != Eigen::Success)
    {
      result.failure = "the hessian did not factorise";
      return result;
    }
    const vector_x y = -solver.solve(r);
    const vector_x d = Z ? vector_x(*Z * y) : y;
    if (solver.info() != Eigen::Success || !d.allFinite())
    {
      result.failure = "the linear solve failed";
      return result;
    }

    const real t = filtered_line_search(problem, x, d, f, g, options);
    if (t == real(0))
    {
      result.failure = "the line search found no decrease";
      return result;
    }
    x += t * d;

    // A step this short cannot improve x further, whatever the gradient says.
    if (t * d.norm() < options.step_tolerance)
    {
      result.iterations = i + 1;
      result.converged = true;
      return result;
    }
  }
}

// motus's backward_euler::step.
newton_result backward_euler::step(const xfem_system& system, mechanical_state& state, real dt,
                                   const vector_x& external_force) {
  const index dofs = system.dofs();
  if (state.x.size() != dofs || state.v.size() != dofs)
  {
    throw std::invalid_argument("xfem: backward_euler: the state does not have the system's dofs");
  }
  if (external_force.size() != 0 && external_force.size() != dofs)
  {
    throw std::invalid_argument("xfem: backward_euler: external_force has the wrong size");
  }

  const vector_x predicted = state.x + dt * state.v;
  const vector_x x_old = state.x;

  vector_x load = external_force.size() == 0 ? vector_x::Zero(dofs) : external_force;
  for (const index d : system.fixed_dofs()) load[d] = real(0);

  const incremental_potential problem(system, predicted, load, dt);
  newton_result result = newton_solve(problem, state.x, options_);
  state.v = (state.x - x_old) / dt;
  return result;
}

}  // namespace xfem
