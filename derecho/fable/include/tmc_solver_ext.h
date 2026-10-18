/*
 * Solver internals and utilities MICM exposes to its tests, for the C test
 * suite: the Rosenbrock error norm and alpha shift, the backward Euler
 * convergence test, the unchecked custom-parameter setter, and MICM's
 * finite-difference Jacobian verification (micm/util/jacobian_verification.hpp).
 */
#ifndef TMC_SOLVER_EXT_H
#define TMC_SOLVER_EXT_H

#include "tmc.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------- */
/* Rosenbrock and backward Euler internals                                    */
/* ------------------------------------------------------------------------- */

/* MICM's AbstractRosenbrockSolver::NormalizedError: the root mean square over
 * the state's cells and variables of errors / (atol + rtol max(|y|, |ynew|)),
 * with the state's tolerances, floored at 1e-10. y, ynew, and errors have the
 * shape and width of the state's variables; padded lanes are excluded. The
 * accumulation order (variable by variable, then cell by cell) is the solver's.
 * Returns NaN when a matrix does not match the state. */
double tmc_solver_normalized_error(const tmc_solver *s, const tmc_state *state, const tmc_dense *y, const tmc_dense *ynew, const tmc_dense *errors);

/* MICM's AbstractRosenbrockSolver::AlphaMinusJacobian: adds alpha times the
 * mass-matrix diagonal (1 for every ODE variable, 0 for an algebraic one) to
 * every diagonal element of the state's Jacobian, which holds -J, giving
 * alpha M - J in every lane; the solver's own function, as it runs before
 * factoring. */
void tmc_solver_alpha_minus_jacobian(const tmc_solver *s, tmc_state *state, double alpha);

/* MICM's AbstractBackwardEuler::IsConverged, the solver's own test: true when
 * every residual of every cell is finite and at most `small`, the variable's
 * absolute tolerance, or relative_tolerance |yn1|. residual and yn1 have one
 * row per cell and one column per variable, in the same layout
 * (absolute_tolerance has one entry per column); padded lanes are excluded. */
bool tmc_backward_euler_is_converged(const tmc_backward_euler_parameters *p, const tmc_dense *residual, const tmc_dense *yn1, const double *absolute_tolerance, double relative_tolerance);

/* MICM's State::UnsafelySetCustomRateParameters: values holds n_cells rows of
 * n_parameters values in the state's parameter order. Fails with
 * TMC_STATE_INVALID_CUSTOM_RATE_PARAM_COUNT_MULTIGRID when n_cells is not the
 * state's cell count, and TMC_STATE_INVALID_CUSTOM_RATE_PARAM_COUNT when
 * n_parameters is not its parameter count. */
int tmc_state_unsafely_set_custom_rate_parameters(tmc_state *state, const double *values, int n_cells, int n_parameters, tmc_error *err);

/* ------------------------------------------------------------------------- */
/* Jacobian verification (MICM's jacobian_verification.hpp)                  */
/* ------------------------------------------------------------------------- */

/* MICM's JacobianComparisonResult */
typedef struct tmc_jacobian_comparison
{
  bool passed;
  double max_abs_error;
  double max_rel_error;
  int worst_block;
  int worst_row;
  int worst_col;
  double worst_analytical;
  double worst_fd;
} tmc_jacobian_comparison;

/* A forcing: forcing += F(variables), for every row (block) of variables. */
typedef void (*tmc_forcing_function)(const tmc_dense *variables, tmc_dense *forcing, void *context);

/* MICM's defaults (double precision) */
#define TMC_FD_PERTURBATION 1.0e-8
#define TMC_FD_COMPARE_ATOL 1.0e-7
#define TMC_FD_COMPARE_RTOL 1.0e-7
#define TMC_FD_SPARSITY_THRESHOLD 1.0e-6

/* MICM's FiniteDifferenceJacobian: central differences with step
 * perturbation * max(1, |x|), or a forward difference where x - h < 0.
 * Returns a matrix with one row per row of variables and n_species^2 columns,
 * element (block, row * n_species + col) = d forcing_row / d x_col, in the
 * vector width of variables. The caller frees it. */
tmc_dense *tmc_finite_difference_jacobian(tmc_forcing_function forcing, void *context, const tmc_dense *variables, int n_species, double perturbation);

/* MICM's CompareJacobianToFiniteDifference: analytical holds -J (MICM's
 * convention) in its non-zero elements, which are compared with the
 * finite-difference J; an element passes when
 * |analytical - fd| <= atol + rtol max(|analytical|, |fd|). */
tmc_jacobian_comparison tmc_compare_jacobian_to_finite_difference(const tmc_sparse *analytical, const tmc_dense *fd, int n_species, double atol, double rtol);

/* MICM's CheckJacobianSparsityCompleteness: fails when a finite-difference
 * entry outside the sparsity pattern of analytical exceeds threshold in
 * magnitude; the worst entry is reported. */
tmc_jacobian_comparison tmc_check_jacobian_sparsity_completeness(const tmc_sparse *analytical, const tmc_dense *fd, int n_species, double threshold);

#ifdef __cplusplus
}
#endif

#endif /* TMC_SOLVER_EXT_H */
