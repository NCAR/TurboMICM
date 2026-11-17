#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>

/* Version of the reference API represented by this port. */
const char *chem_version(void);
unsigned chem_version_major(void);
unsigned chem_version_minor(void);
unsigned chem_version_patch(void);
unsigned chem_version_tweak(void);

/* All mechanism data are immutable and may be shared between solver states.
 * A caller may initialize these structs in C or link a generated definition. */
typedef enum {
    CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_TERNARY_JPL19,
    CHEM_BRANCHED, CHEM_TUNNELING, CHEM_TAYLOR, CHEM_REVERSIBLE,
    CHEM_USER, CHEM_SURFACE, CHEM_CALLBACK
} ChemRateKind;

typedef double (*ChemRateFunction)(double temperature, double pressure, double air_density, void *context);
typedef struct { double constant, temperature, pressure, air_density; } ChemParameterizedFunction;

typedef struct {
    const char *name;
    double initial, absolute_tolerance, molecular_weight, diffusion;
} ChemSpecies;

typedef struct { int species; double yield; } ChemProduct;
typedef struct { int species; double coefficient; } ChemTerm;

typedef enum { CHEM_LINEAR_CONSTRAINT, CHEM_EQUILIBRIUM_CONSTRAINT } ChemConstraintKind;
typedef struct {
    const char *name;
    ChemConstraintKind kind;
    int algebraic_species;
    const ChemTerm *reactants, *products;
    size_t reactant_count, product_count;
    double constant;
    double equilibrium_constant, enthalpy, gas_constant, reference_temperature;
} ChemConstraint;

typedef struct {
    ChemRateKind kind;
    int parameter;                 /* index, or -1; surface uses index and index+1 */
    double a, b, c, d, e;
    double k0_a, k0_b, k0_c, k0_d;
    double kinf_a, kinf_b, kinf_c, kinf_d;
    double kint_a, kint_b, kint_c, kint_d;
    double fc, n;
    double branch_x, branch_y, branch_a0, branch_n;
    int nitrate_branch;
    double diffusion, molecular_weight, probability;
    size_t coefficient_count;
    double coefficients[16];
    const char *label;
    ChemRateFunction function;
    void *context;
    /* Parameterized reactants multiply the rate by their concentration. */
    size_t parameterization_count;
    const ChemParameterizedFunction *parameterizations;
} ChemRate;

typedef struct {
    size_t reactant_start, reactant_count, product_start, product_count;
    unsigned air_density_power;
    ChemRate rate;
} ChemReaction;

typedef struct {
    size_t species_count, reaction_count, parameter_count;
    const ChemSpecies *species;
    const ChemReaction *reactions;
    const int *reactants;
    const ChemProduct *products;
    const char *const *parameter_names;
    const double *parameter_defaults;
    double initial_temperature, initial_pressure;
} ChemMechanism;

/* Stable rate-type index into the mechanism's reaction array. Entries for
 * kind k occupy offsets[k] through offsets[k+1] in reaction_indices. */
typedef struct { size_t offsets[CHEM_CALLBACK+2]; size_t *reaction_indices; } ChemRateIndex;
int chem_rate_index_init(ChemRateIndex *index, const ChemMechanism *mechanism);
void chem_rate_index_destroy(ChemRateIndex *index);
int chem_branched_constants(const ChemRate *rate, double *k0, double *z);
int chem_surface_speed_factor(const ChemRate *rate, double *factor);

typedef struct {
    unsigned long internal_steps, accepted, rejected, function_calls;
    unsigned long jacobian_updates, decompositions, solves, constraint_init_iterations;
    double final_time;
} ChemStats;

/* The arrays a and c store the strictly lower triangular matrices in row order.
 * The maximum of six stages matches the supported Rosenbrock parameter sets. */
typedef struct ChemRosenbrock ChemRosenbrock;
typedef double (*ChemStepController)(double error, double step, const ChemRosenbrock *parameters, void *context);

struct ChemRosenbrock {
    size_t stages;
    double a[15], c[15], m[6], e[6], gamma[6];
    int new_function[6];
    double order, relative_tolerance, h_min, h_max, h_start;
    double factor_min, factor_max, safety, rejection_factor;
    unsigned max_steps;
    unsigned constraint_max_iterations, constraint_max_backtracks;
    double constraint_tolerance, constraint_backtrack_factor, constraint_sufficient_decrease;
    ChemStepController step_controller;
    void *step_controller_context;
};

typedef struct {
    double small, h_start, relative_tolerance;
    unsigned max_iterations;
    double time_step_reductions[5];
} ChemBackwardEuler;

typedef struct ChemState ChemState;
typedef struct { const char *name; const double *values; size_t count; } ChemNamedValues;

/* Matrix entries use sparse row order; values[entry*cells+cell] holds one
 * value per cell. lookup[row*species_count+column] is -1 for structural zero.
 * factor must add alpha only to differential diagonals before factorizing. */
typedef struct {
    size_t species_count, cells, nonzeros;
    const size_t *row_start;
    const int *columns, *lookup;
    const unsigned char *algebraic;
    double *values;
} ChemLinearSystem;
typedef struct {
    int (*factor)(ChemLinearSystem *system, double alpha, void *context);
    int (*solve)(const ChemLinearSystem *system, double *rhs, void *context);
} ChemLinearSolver;

/* Extension callbacks operate on the same species-major arrays as the solver.
 * The mechanism must include all species and parameters used by the models.
 * Descriptors, dependencies, and callback contexts are borrowed and must outlive
 * all states created from them. Finalize may cache sparse entry indices (not
 * cell-dependent value offsets); use a separate context for different plans.
 * A constraint callback receives a pointer to the selected cell and reads
 * values[species * stride]; its gradient, when non-NULL, is a zero-filled
 * species_count-element array. Declare every species with a possible nonzero
 * gradient in dependencies. */
typedef struct { size_t row, column; } ChemJacobianEntry;
typedef struct {
    const char *name;
    int algebraic_species;
    const int *dependencies;
    size_t dependency_count;
    double (*evaluate)(const ChemState *state, size_t cell, const double *values,
                       size_t stride, double *gradient, void *context);
    void *context;
} ChemExternalConstraint;
typedef struct {
    const ChemJacobianEntry *jacobian_pattern;
    size_t jacobian_count;
    const ChemExternalConstraint *constraints;
    size_t constraint_count;
    int (*finalize)(const ChemMechanism *mechanism, const ChemLinearSystem *pattern, void *context);
    int (*update)(ChemState *state, void *context);
    void (*add_forcing)(const ChemState *state, const double *values, double *forcing, void *context);
    void (*subtract_jacobian)(const ChemState *state, const double *values,
                              ChemLinearSystem *negative_jacobian, void *context);
    void *context;
    /* Optional additions for chem_mechanism_compose. */
    const ChemSpecies *species;
    size_t species_count;
    const char *const *parameter_names;
    const double *parameter_defaults;
    size_t parameter_count;
} ChemExternalModel;
typedef struct {
    ChemMechanism mechanism;
    ChemSpecies *species_storage;
    const char **parameter_name_storage;
    double *parameter_default_storage;
} ChemMechanismComposition;
/* Copy base metadata and append model-owned declarations into an empty result.
 * String data, reactions and contexts are borrowed; destroy after every
 * associated state. The result is unchanged on failure. */
int chem_mechanism_compose(ChemMechanismComposition *result, const ChemMechanism *base,
                           const ChemExternalModel *models, size_t model_count);
void chem_mechanism_composition_destroy(ChemMechanismComposition *composition);

ChemRosenbrock chem_rosenbrock_three_stage(void);
ChemRosenbrock chem_rosenbrock_two_stage(void);
ChemRosenbrock chem_rosenbrock_four_stage(void);
ChemRosenbrock chem_rosenbrock_four_stage_dae(void);
ChemRosenbrock chem_rosenbrock_six_stage_dae(void);
ChemBackwardEuler chem_backward_euler_defaults(void);
ChemRate chem_rate_defaults(ChemRateKind kind);
int chem_rate_calculate(const ChemRate *rate, double temperature, double pressure,
                        double air_density, const double *parameters, size_t parameter_count,
                        double *result);
ChemState *chem_state_create(const ChemMechanism *mechanism, size_t cells);
ChemState *chem_state_create_constrained(const ChemMechanism *mechanism, const ChemConstraint *constraints,
                                         size_t constraint_count, size_t cells);
ChemState *chem_state_create_extended(const ChemMechanism *mechanism, const ChemConstraint *constraints,
                                     size_t constraint_count, const ChemExternalModel *models,
                                     size_t model_count, size_t cells);
int chem_validate_constraints(const ChemMechanism *mechanism, const ChemConstraint *constraints, size_t count);
/* pattern is a zeroed-or-overwritten species_count-square row-major byte array. */
int chem_constraint_pattern(const ChemMechanism *mechanism, const ChemConstraint *constraints,
                            size_t count, unsigned char *pattern);
int chem_constraint_evaluate(const ChemConstraint *constraint, const double *values, size_t species_count,
                             double temperature, double *residual, double *gradient);
int chem_constraints_apply(const ChemState *state, size_t cell, const double *values,
                           double *forcing, double *negative_jacobian);
int chem_initialize_constraints(ChemState *state, const ChemRosenbrock *parameters, ChemStats *stats);
/* Share the immutable symbolic linear-solver plan with an existing state. */
ChemState *chem_state_create_like(const ChemState *state, size_t cells);
/* Deep-copy mutable state; immutable mechanism data and callback contexts remain borrowed. */
ChemState *chem_state_clone(const ChemState *state);
void chem_state_free(ChemState *state);
size_t chem_cell_count(const ChemState *state);
/* Each returned array has chem_cell_count() entries. */
double *chem_species_values(ChemState *state, size_t species);
const double *chem_species_values_const(const ChemState *state, size_t species);
double *chem_parameter_values(ChemState *state, size_t parameter);
const double *chem_parameter_values_const(const ChemState *state, size_t parameter);
double *chem_temperatures(ChemState *state);
double *chem_pressures(ChemState *state);
double *chem_absolute_tolerances(ChemState *state);
/* Rate constants from the most recent update, one entry per cell. */
const double *chem_rate_constants(const ChemState *state, size_t reaction);
double *chem_constraint_constants(ChemState *state, size_t constraint);
int chem_find_constraint(const ChemState *state, const char *name);
int chem_is_algebraic(const ChemState *state, size_t species);
/* By default air density is calculated from temperature and pressure at rate
 * update. An explicit density can be installed independently for each cell. */
int chem_set_air_density(ChemState *state, size_t cell, double density);
int chem_use_ideal_air_density(ChemState *state, size_t cell);
int chem_find_species(const ChemMechanism *mechanism, const char *name);
int chem_find_parameter(const ChemMechanism *mechanism, const char *name);
/* Named inputs must contain exactly one value per cell. Bulk setters validate
 * every entry before changing any values. Row-major parameters are indexed
 * as values[cell * parameter_count + parameter]. */
int chem_set_concentrations(ChemState *state, const ChemNamedValues *items, size_t count);
int chem_set_parameters(ChemState *state, const ChemNamedValues *items, size_t count);
int chem_set_parameters_row_major(ChemState *state, const double *values, size_t cells, size_t parameter_count);
int chem_override_rate_function(ChemState *state, const char *label, ChemRateFunction function, void *context);
/* NULL restores the default in-place Mozart solver; callback context is borrowed. */
int chem_set_linear_solver(ChemState *state, const ChemLinearSolver *solver, void *context);
ChemLinearSolver chem_doolittle_linear_solver(void);
ChemLinearSolver chem_mozart_linear_solver(void);
/* Copy and factor a matrix without modifying the input. Both views must have
 * the same sparse pattern and distinct value arrays. */
int chem_linear_factor_copy(const ChemLinearSolver *solver, const ChemLinearSystem *input,
                            ChemLinearSystem *factors, double alpha, void *context);
int chem_linear_add_alpha(ChemLinearSystem *system, double alpha);
int chem_update_rates(ChemState *state);
/* Negative Jacobian, in dense row-major form for one cell, and forcing, for
 * diagnostics and tests. The solver itself uses sparse storage. */
int chem_evaluate(const ChemState *state, size_t cell, double *forcing, double *negative_jacobian);
int chem_solve(ChemState *state, const ChemRosenbrock *parameters, double dt, ChemStats *stats);
/* Arrays are species-major across cells, as returned by chem_species_values. */
int chem_rosenbrock_normalized_error(size_t species_count, size_t cells, const double *old_values,
                                     const double *new_values, const double *errors,
                                     const double *absolute_tolerances, double relative_tolerance, double *result);
int chem_backward_euler_converged(size_t species_count, size_t cells, const double *residual,
                                  const double *values, const double *absolute_tolerances,
                                  double relative_tolerance, double small);
int chem_solve_backward_euler(ChemState *state, const ChemBackwardEuler *parameters, double dt, ChemStats *stats);
const char *chem_error_string(int code);

enum { CHEM_OK = 0, CHEM_INVALID = 1, CHEM_MEMORY = 2, CHEM_SINGULAR = 3,
       CHEM_TOO_MANY_STEPS = 4, CHEM_STEP_TOO_SMALL = 5, CHEM_NONFINITE = 6,
       CHEM_MISSING_DIFFUSION = 7, CHEM_MISSING_MOLECULAR_WEIGHT = 8,
       CHEM_UNCONVERGED = 9, CHEM_CONSTRAINT_INITIALIZATION_FAILED = 10 };

#endif
