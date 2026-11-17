#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>

const char *chem_version(void);
unsigned chem_version_major(void);
unsigned chem_version_minor(void);
unsigned chem_version_patch(void);
unsigned chem_version_tweak(void);

/* Terms are ordered; repeated reactants have separate entries.  A reactant
 * index of -1 denotes air density; lower values index parameterized species. */
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
typedef struct { double constant, temperature, pressure, air_density; } ChemParameterization;
typedef struct {
    size_t species_count, reaction_count, parameter_count;
    const char *const *species_names, *const *parameter_names;
    const double *initial_species, *initial_parameters;
    const ChemReaction *reactions;
    const ChemTerm *reactants, *products;
    size_t parameterized_species_count;
    const ChemParameterization *parameterized_species;
} ChemMechanism;
typedef struct {
    double temperature, pressure;
    double air_density; /* used when has_air_density is nonzero */
    int has_air_density;
} ChemConditions;
double chem_parameterized_value(ChemParameterization, ChemConditions);
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
/* Add reaction forcing or subtract its Jacobian, leaving masked algebraic
 * rows untouched. A null mask treats every row as differential. */
int chem_process_add_forcing(const ChemMechanism *, ChemConditions, const double *rates, const double *y,
                             const unsigned char *algebraic_rows, double *forcing);
int chem_process_subtract_jacobian(const ChemMechanism *, ChemConditions, const double *rates, const double *y,
                                   const unsigned char *algebraic_rows, double *jacobian);
int chem_process_jacobian_pattern(const ChemMechanism *, unsigned char *pattern);

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
/* Generic sparse-pattern linear algebra using the same LU kernels as the
 * integrators. Pattern and matrices are row-major n-by-n arrays. Diagonal
 * entries are included automatically; zero pattern entries represent zeros. */
ChemWorkspace *chem_linear_workspace_create(size_t n, const unsigned char *pattern);
ChemWorkspace *chem_linear_workspace_create_with_initial(size_t n, const unsigned char *pattern, double initial);
int chem_linear_factor(ChemWorkspace *, const double *matrix, ChemLUAlgorithm);
/* Overwrite matrix with L below the diagonal and U on and above it. */
int chem_linear_factor_in_place(ChemWorkspace *, double *matrix, ChemLUAlgorithm);
int chem_linear_factors(const ChemWorkspace *, double *lower, double *upper);
int chem_linear_solve(const ChemWorkspace *, double *rhs);
/* Symmetric diagonal Markowitz ordering and symbolic L+U storage count. */
int chem_diagonal_markowitz_order(size_t n, const unsigned char *pattern, size_t *order);
size_t chem_linear_symbolic_nonzeros(const ChemWorkspace *);
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
    CHEM_CONSTRAINT_BAD_STOICHIOMETRY, CHEM_CONSTRAINT_BAD_EQUILIBRIUM, CHEM_CONSTRAINT_DUPLICATE_ROW
} ChemConstraintStatus;
typedef struct { const ChemConstraint *items; size_t count; } ChemConstraintSet;
ChemConstraintStatus chem_constraint_validate(const ChemConstraint *, size_t species_count);
ChemConstraintStatus chem_constraint_set_validate(ChemConstraintSet, size_t species_count);
ChemConstraintStatus chem_constraint_set_pattern(ChemConstraintSet, size_t species_count, unsigned char *pattern);
/* Cells are contiguous; temperatures, forcing, and Jacobians have one entry/block per cell. */
ChemConstraintStatus chem_constraint_set_apply(ChemConstraintSet, size_t species_count, size_t cells,
                                               const double *temperatures, const double *y,
                                               double *forcing, double *negative_jacobians);
/* Overwrite forcing[row] and the negative-Jacobian row; other rows are unchanged. */
ChemConstraintStatus chem_constraint_apply(const ChemConstraint *, size_t species_count,
                                          double temperature, const double *y,
                                          double *forcing, double *negative_jacobian);

/* Model-building metadata. Names and string values are borrowed; copies own
 * independent property arrays, while the referenced strings remain borrowed. */
typedef enum { CHEM_PROPERTY_REAL, CHEM_PROPERTY_STRING, CHEM_PROPERTY_INT, CHEM_PROPERTY_BOOL } ChemPropertyType;
typedef union { double real; const char *string; int integer; int boolean; } ChemPropertyValue;
typedef struct { const char *name; ChemPropertyType type; ChemPropertyValue value; } ChemProperty;
typedef struct {
    const char *name;
    ChemProperty *properties;
    size_t property_count, property_capacity;
    int is_parameterized;
    ChemParameterization parameterization;
} ChemSpecies;
typedef enum { CHEM_PROPERTY_OK, CHEM_PROPERTY_NOT_FOUND, CHEM_PROPERTY_INVALID_TYPE, CHEM_PROPERTY_NO_MEMORY } ChemPropertyStatus;
void chem_species_init(ChemSpecies *, const char *name);
void chem_species_free(ChemSpecies *);
ChemPropertyStatus chem_species_copy(ChemSpecies *destination, const ChemSpecies *source);
ChemPropertyStatus chem_species_set(ChemSpecies *, const char *name, ChemPropertyType, ChemPropertyValue);
ChemPropertyStatus chem_species_get(const ChemSpecies *, const char *name, ChemPropertyType, ChemPropertyValue *);
void chem_species_set_third_body(ChemSpecies *);
typedef struct { const ChemSpecies *species; int has_diffusion; double diffusion; } ChemPhaseSpecies;
void chem_phase_species_set_diffusion(ChemPhaseSpecies *, double diffusion);
typedef struct { const char *name; const ChemPhaseSpecies *species; size_t count; } ChemPhase;
typedef struct { ChemPhase gas_phase; } ChemSystem;
/* Reaction term indices refer to the owning mechanism's term arrays; phase
 * metadata and its species are borrowed. Plain assignment copies the value. */
typedef struct { ChemReaction reaction; ChemPhase phase; } ChemProcess;
size_t chem_phase_state_size(ChemPhase);
/* Returns zero on success; output receives phase.name or the unqualified name. */
int chem_phase_state_name(ChemPhase, size_t state_index, int qualified, char *output, size_t capacity);
size_t chem_system_state_size(ChemSystem);
int chem_system_state_name(ChemSystem, size_t state_index, const size_t *order, char *output, size_t capacity);

#endif
