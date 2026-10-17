/*
 * Constraint extensions to tmc.h: two parts of MICM's constraint API that the
 * solver otherwise uses only internally.
 */
#ifndef TMC_CONSTRAINT_EXT_H
#define TMC_CONSTRAINT_EXT_H

#include "tmc.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The Jacobian elements of the solver's constraints, as MICM's
 * ConstraintSet::NonZeroJacobianElements lists them: for each constraint, the
 * diagonal of its algebraic row, (row, row), and (row, dependency) for each
 * species it depends on, as a set ordered by row, then column. Writes up to
 * max pairs to out (which may be NULL) and returns the number of pairs.
 * Constraints of external models are not included. */
int tmc_solver_constraint_jacobian_elements(const tmc_solver *s, tmc_index_pair *out, int max);

/* MICM's AbstractRosenbrockSolver::InitializeConstraints: the Newton projection
 * of the state's algebraic variables onto the constraints, which tmc_solver_solve
 * runs before integrating, here run alone. It uses the given parameters (the
 * solver's when NULL), takes no integration step, and does not clamp. Returns
 * TMC_CONVERGED, or the failure (TMC_CONSTRAINT_INITIALIZATION_FAILED,
 * TMC_NAN_DETECTED, TMC_INF_DETECTED), after which the state's variables are
 * as they were passed in. The projection's statistics are added to *stats
 * when stats is not NULL. Rosenbrock solvers only: for another method it
 * returns TMC_CONSTRAINT_INITIALIZATION_FAILED and leaves the state alone. */
tmc_solver_state tmc_solver_initialize_constraints(const tmc_solver *s, tmc_state *state, const tmc_rosenbrock_parameters *parameters, tmc_stats *stats);

#ifdef __cplusplus
}
#endif

#endif /* TMC_CONSTRAINT_EXT_H */
