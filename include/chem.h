#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>

/* Terms are ordered; repeated reactants have separate entries.  A reactant
 * index of -1 denotes the third body, parameterized by air density. */
typedef struct { int species; double yield; } ChemTerm;
typedef enum {
    CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_TUNNELING,
    CHEM_USER, CHEM_SURFACE, CHEM_BRANCH_ALKOXY, CHEM_BRANCH_NITRATE,
    CHEM_REVERSIBLE, CHEM_TAYLOR, CHEM_TERNARY_JPL19, CHEM_CALLBACK
} ChemRateType;
typedef struct {
    ChemRateType type;
    size_t reactant_begin, reactant_count, product_begin, product_count;
    int parameter;                 /* first custom parameter, or -1 */
    double p[22];                  /* interpreted by the rate type */
} ChemReaction;
typedef struct {
    size_t species_count, reaction_count, parameter_count;
    const char *const *species_names, *const *parameter_names;
    const double *initial_species, *initial_parameters;
    const ChemReaction *reactions;
    const ChemTerm *reactants, *products;
} ChemMechanism;
typedef struct {
    double temperature, pressure;
    double air_density; /* used when has_air_density is nonzero */
    int has_air_density;
} ChemConditions;
typedef double (*ChemRateCallback)(ChemConditions, void *context);
typedef struct { ChemRateCallback function; void *context; } ChemCallback;
ptrdiff_t chem_parameter_index(const ChemMechanism *, const char *name);

typedef enum { CHEM_VALID, CHEM_INVALID_MECHANISM, CHEM_MISSING_DIFFUSION, CHEM_MISSING_MOLECULAR_WEIGHT } ChemValidation;
ChemValidation chem_validate_mechanism(const ChemMechanism *, char *message, size_t message_size);

/* Return nonzero for invalid conditions, arguments, or non-finite results. */
int chem_rate_constants(const ChemMechanism *, ChemConditions, const double *parameters, double *rates);
/* Callback entries share the custom-parameter index space and may be replaced
 * by the caller between updates. */
int chem_rate_constants_with_callbacks(const ChemMechanism *, ChemConditions, const double *parameters,
                                       const ChemCallback *callbacks, double *rates);
void chem_forcing(const ChemMechanism *, ChemConditions, const double *rates, const double *y, double *f);
void chem_jacobian(const ChemMechanism *, ChemConditions, const double *rates, const double *y, double *j);

/* Rosenbrock parameters are data, so the step controller and formula can be
 * changed independently. Arrays use row-wise strictly lower triangular order. */
typedef enum { CHEM_LU_MOZART, CHEM_LU_DOOLITTLE } ChemLUAlgorithm;
typedef struct {
    int stages, new_function[6], max_steps;
    ChemLUAlgorithm linear_solver;
    double a[15], c[15], m[6], e[6], gamma, order;
    double factor_min, factor_max, rejection_factor, safety;
    double h_min, h_max, h_start, relative_tolerance, absolute_tolerance;
    int constraint_init_max_iterations, constraint_init_max_backtracks;
    double constraint_init_tolerance, constraint_init_backtrack_factor, constraint_init_sufficient_decrease;
    const double *absolute_tolerances; /* optional per-species override */
} ChemRosenbrockParameters;
typedef struct {
    unsigned long internal_steps, accepted, rejected, function_calls;
    unsigned long jacobian_updates, decompositions, solves, constraint_init_iterations;
    double final_time;
} ChemStatistics;
typedef enum {
    CHEM_SOLVE_SUCCESS = 0, CHEM_SOLVE_BAD_ARGUMENT = 1, CHEM_SOLVE_STEP_TOO_SMALL = 2,
    CHEM_SOLVE_SINGULAR = 3, CHEM_SOLVE_NONFINITE = 4, CHEM_SOLVE_CONSTRAINT_INITIALIZATION_FAILED = 5,
    CHEM_SOLVE_ACCEPTING_UNCONVERGED = 6
} ChemSolverStatus;
const char *chem_solver_status_string(ChemSolverStatus);
typedef struct ChemWorkspace ChemWorkspace;
typedef struct ChemConstraint ChemConstraint;

ChemRosenbrockParameters chem_rosenbrock_two_stage(void);
ChemRosenbrockParameters chem_rosenbrock_three_stage(void);
ChemRosenbrockParameters chem_rosenbrock_four_stage(void);
ChemRosenbrockParameters chem_rosenbrock_four_stage_dae(void);
ChemRosenbrockParameters chem_rosenbrock_six_stage_dae(void);
ChemWorkspace *chem_workspace_create(const ChemMechanism *);
/* Constraints and their term arrays must outlive the workspace. */
ChemWorkspace *chem_workspace_create_constrained(const ChemMechanism *, const ChemConstraint *, size_t);
void chem_workspace_free(ChemWorkspace *);
/* Initialize algebraic rows transactionally; accumulates statistics. */
int chem_initialize_constraints(ChemWorkspace *, const ChemRosenbrockParameters *,
                                ChemConditions, double *y, ChemStatistics *);
/* On failure returns nonzero; y contains the most recently accepted state. */
int chem_rosenbrock_solve(const ChemMechanism *, ChemWorkspace *,
                          const ChemRosenbrockParameters *, ChemConditions,
                          const double *rates, double dt, double *y, ChemStatistics *);

typedef struct {
    double small, h_start, reductions[5], relative_tolerance, absolute_tolerance;
    const double *absolute_tolerances;
    int max_iterations;
    ChemLUAlgorithm linear_solver;
} ChemBackwardEulerParameters;
ChemBackwardEulerParameters chem_backward_euler_defaults(void);
int chem_backward_euler_converged(size_t cells, size_t species, const double *residual,
                                  const double *y, double small, double rtol, const double *atol);
int chem_backward_euler_solve(const ChemMechanism *, ChemWorkspace *,
                              const ChemBackwardEulerParameters *, ChemConditions,
                              const double *rates, double dt, double *y, ChemStatistics *);

typedef enum { CHEM_LINEAR_CONSTRAINT, CHEM_EQUILIBRIUM_CONSTRAINT } ChemConstraintType;
struct ChemConstraint {
    const char *name;
    ChemConstraintType type;
    size_t row;
    const ChemTerm *reactants, *products;
    size_t reactant_count, product_count;
    double constant; /* linear constant, or equilibrium constant at reference temperature */
    double enthalpy, gas_constant, reference_temperature;
};
typedef enum {
    CHEM_CONSTRAINT_OK, CHEM_CONSTRAINT_BAD_ROW, CHEM_CONSTRAINT_EMPTY_REACTANTS,
    CHEM_CONSTRAINT_EMPTY_PRODUCTS, CHEM_CONSTRAINT_BAD_SPECIES,
    CHEM_CONSTRAINT_BAD_STOICHIOMETRY, CHEM_CONSTRAINT_BAD_EQUILIBRIUM
} ChemConstraintStatus;
ChemConstraintStatus chem_constraint_validate(const ChemConstraint *, size_t species_count);
/* Overwrite forcing[row] and the negative-Jacobian row; other rows are unchanged. */
ChemConstraintStatus chem_constraint_apply(const ChemConstraint *, size_t species_count,
                                          double temperature, const double *y,
                                          double *forcing, double *negative_jacobian);

#endif
