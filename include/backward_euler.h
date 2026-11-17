#ifndef BACKWARD_EULER_H
#define BACKWARD_EULER_H
#include "rosenbrock.h"

typedef struct {
    double small, h_start, time_step_reductions[5];
    size_t max_iterations;
    double relative_tolerance;
    const double *absolute_tolerance;
    const ChemSparsePlan *sparse_plan;
    ChemLuAlgorithm lu_algorithm;
} ChemBackwardEulerOptions;

ChemBackwardEulerOptions chem_backward_euler_defaults(const double *absolute_tolerance);
/* State and residual arrays have cells*variables entries in row-major order. */
int chem_backward_euler_converged(size_t cells, size_t variables, const double *residual,
    const double *state, double small, const double *absolute_tolerance, double relative_tolerance);
ChemSolveStatus chem_backward_euler_solve(size_t n, double *y, double duration,
    ChemOdeFunction rhs, ChemOdeFunction jacobian, void *context,
    const ChemBackwardEulerOptions *options, ChemSolverStats *stats);
/* Advance independent cells with shared step acceptance and convergence.
   Tolerances are per variable; cell_pattern is an optional row-major Jacobian
   mask used to construct a block-sparse plan. options->sparse_plan must be NULL. */
ChemSolveStatus chem_backward_euler_solve_grouped(size_t cells, size_t variables,
    double *y, double duration, ChemOdeFunction rhs, ChemOdeFunction jacobian,
    void *const *contexts, const unsigned char *cell_pattern,
    const ChemBackwardEulerOptions *options, ChemSolverStats *stats);
#endif
