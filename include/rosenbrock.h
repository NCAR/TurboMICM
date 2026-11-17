#ifndef ROSENBROCK_H
#define ROSENBROCK_H
#include <stddef.h>
#include "sparse.h"
#include "linear.h"
#include "constraint.h"

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
    const ChemSparsePlan *sparse_plan; /* NULL uses dense elimination */
    ChemLuAlgorithm lu_algorithm;
    const ChemConstraintSet *constraints;
    const double *constraint_parameters; /* optional one value per constraint */
    double constraint_temperature;
    size_t constraint_init_max_iterations, constraint_init_max_backtracks;
    double constraint_init_tolerance, constraint_init_backtrack_factor, constraint_init_sufficient_decrease;
} ChemRosenbrockOptions;
typedef struct {
    size_t steps, accepted, rejected, function_calls, jacobian_updates;
    size_t decompositions, solves, constraint_init_iterations;
    double final_time;
} ChemSolverStats;
typedef enum {
    CHEM_SOLVE_OK, CHEM_SOLVE_INVALID, CHEM_SOLVE_ALLOCATION,
    CHEM_SOLVE_SINGULAR, CHEM_SOLVE_MAX_STEPS, CHEM_SOLVE_STEP_TOO_SMALL,
    CHEM_SOLVE_NONFINITE, CHEM_SOLVE_UNCONVERGED, CHEM_SOLVE_CONSTRAINT_INIT
} ChemSolveStatus;

const char *chem_solve_status_string(ChemSolveStatus status);
ChemRosenbrockMethod chem_rosenbrock_two_stage(void);
ChemRosenbrockMethod chem_rosenbrock_three_stage(void);
ChemRosenbrockMethod chem_rosenbrock_four_stage(void);
/* These parameter sets can solve ODEs and index-one DAEs. */
ChemRosenbrockMethod chem_rosenbrock_four_stage_da(void);
ChemRosenbrockMethod chem_rosenbrock_six_stage_da(void);
/* Row-major cell-by-variable arrays; tolerances are per variable. */
double chem_rosenbrock_normalized_error(size_t cells, size_t variables,
    const double *old_state, const double *new_state, const double *errors,
    const double *absolute_tolerance, double relative_tolerance);
ChemRosenbrockOptions chem_rosenbrock_defaults(const double *absolute_tolerance);
/* Project algebraic variables onto their constraint manifold; failure restores y.
   Statistics, if supplied, are accumulated rather than reset. */
ChemSolveStatus chem_rosenbrock_initialize_constraints(size_t n, double *y,
    const ChemRosenbrockOptions *options, ChemSolverStats *stats);
/* Modifies y in place and clamps negative ODE values after solving. Each call
   resets the internal step size, as in the reference CPU solver. */
ChemSolveStatus chem_rosenbrock_solve(size_t n, double *y, double duration,
    ChemOdeFunction rhs, ChemOdeFunction jacobian, void *context,
    const ChemRosenbrockMethod *method, const ChemRosenbrockOptions *options,
    ChemSolverStats *stats);
#endif
