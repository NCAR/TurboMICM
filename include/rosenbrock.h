#ifndef ROSENBROCK_H
#define ROSENBROCK_H
#include <stddef.h>

/* Autonomous ODE callbacks. Jacobian is dense and row-major. */
typedef void (*ChemOdeFunction)(void *context, const double *y, double *out);
typedef struct {
    size_t stages;
    double a[15], c[15], m[6], e[6], gamma0, error_order;
    int new_function[6];
} ChemRosenbrockMethod;
typedef struct {
    double relative_tolerance;
    const double *absolute_tolerance; /* one positive tolerance per variable */
    double h_min, h_max, h_start; /* zero chooses the reference default */
    double factor_min, factor_max, safety, rejection_factor;
    size_t max_steps;
} ChemRosenbrockOptions;
typedef struct {
    size_t steps, accepted, rejected, function_calls, jacobian_updates;
    size_t decompositions, solves;
} ChemSolverStats;
typedef enum {
    CHEM_SOLVE_OK, CHEM_SOLVE_INVALID, CHEM_SOLVE_ALLOCATION,
    CHEM_SOLVE_SINGULAR, CHEM_SOLVE_MAX_STEPS, CHEM_SOLVE_STEP_TOO_SMALL,
    CHEM_SOLVE_NONFINITE
} ChemSolveStatus;

ChemRosenbrockMethod chem_rosenbrock_three_stage(void);
ChemRosenbrockOptions chem_rosenbrock_defaults(const double *absolute_tolerance);
/* Modifies y in place and clamps negative values after solving. Each call
   resets the internal step size, as in the reference CPU solver. */
ChemSolveStatus chem_rosenbrock_solve(size_t n, double *y, double duration,
    ChemOdeFunction rhs, ChemOdeFunction jacobian, void *context,
    const ChemRosenbrockMethod *method, const ChemRosenbrockOptions *options,
    ChemSolverStats *stats);
#endif
