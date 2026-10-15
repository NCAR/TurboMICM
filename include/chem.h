#ifndef CHEM_H
#define CHEM_H
#include <stddef.h>

/* All rates use concentrations in mol m^-3, temperature K, pressure Pa. */
typedef enum {
    CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_TERNARY_JPL19,
    CHEM_TUNNELING, CHEM_BRANCHED_ALKOXY, CHEM_BRANCHED_NITRATE,
    CHEM_TAYLOR, CHEM_REVERSIBLE, CHEM_USER, CHEM_SURFACE, CHEM_CALLBACK
} ChemRateType;
typedef double (*ChemRateFunction)(double temperature, double pressure,
                                   double air_density, const double *custom, void *context);
typedef struct {
    ChemRateType type;
    /* Arrhenius/Taylor: A B C D E, then polynomial coefficients.
       Troe/ternary: k0 A B C, kinf A B C, Fc N.
       JPL19: k0 A B C D, kinf A B C D, kint A B C D, Fc N.
       Tunnel: A B C. Branched: X Y a0 n. Reversible: A C kr.
       User: scaling. Surface: diffusion, molecular weight, probability. */
    double p[24];
    int parameter, coefficient_count, third_body_count;
    ChemRateFunction function;
    void *context;
    const char *name;
} ChemRate;
typedef struct { int species; double yield; } ChemProduct;
typedef struct {
    size_t reactant_count, product_count;
    const int *reactants;
    const ChemProduct *products;
    ChemRate rate;
} ChemReaction;
typedef struct ChemMechanism ChemMechanism;
typedef struct ChemState ChemState;
typedef enum { CHEM_LINEAR_CONSTRAINT, CHEM_EQUILIBRIUM_CONSTRAINT, CHEM_CALLBACK_CONSTRAINT } ChemConstraintType;
/* A constraint callback returns G and optionally fills its dense gradient.
   y and custom are contiguous values for one cell; gradient has species_count entries. */
typedef double (*ChemConstraintFunction)(double temperature, const double *y,
    const double *custom, double *gradient, void *context);
typedef struct {
    ChemConstraintType type;
    int algebraic_species;
    size_t reactant_count, product_count;
    const ChemProduct *reactants, *products; /* linear terms use reactants */
    double constant; /* linear */
    double equilibrium_constant, enthalpy, gas_constant, reference_temperature;
    int parameter; /* equilibrium parameter column, -1 to calculate directly */
    ChemConstraintFunction function;
    void *context;
    const char *name;
} ChemConstraint;
typedef struct { size_t row, column; } ChemMatrixIndex;
typedef void (*ChemModelFunction)(ChemState *, const double *y, double *output, void *context);
typedef struct {
    size_t nonzero_count;
    const ChemMatrixIndex *nonzeros;
    ChemModelFunction forcing, negative_jacobian;
    void (*update)(ChemState *, void *context);
    void *context;
} ChemExternalModel;
typedef struct {
    size_t species_count, reaction_count, parameter_count;
    const char *const *species_names;
    const char *const *parameter_names;
    const double *initial_concentrations, *initial_parameters, *absolute_tolerances;
    double temperature, pressure;
    const ChemReaction *reactions;
    size_t constraint_count;
    const ChemConstraint *constraints;
    size_t external_model_count;
    const ChemExternalModel *external_models;
} ChemMechanismData;
typedef enum { CHEM_OK, CHEM_INVALID, CHEM_MEMORY, CHEM_MAX_STEPS,
               CHEM_SMALL_STEP, CHEM_NONFINITE, CHEM_SINGULAR, CHEM_UNCONVERGED, CHEM_INITIALIZATION_FAILED } ChemStatus;
typedef struct {
    size_t internal_steps, accepted, rejected, function_calls,
           jacobian_updates, decompositions, solves;
    double final_time;
    size_t constraint_init_iterations;
} ChemStats;
typedef double (*ChemStepController)(double error, double order, double safety,
                                     double factor_min, double factor_max, void *context);
typedef struct {
    ChemStatus (*factor)(const ChemMechanism *, size_t cells, double *matrix);
    void (*solve)(const ChemMechanism *, size_t cells, const double *matrix, double *rhs);
} ChemLinearSolver;
typedef struct {
    double small, h_start;
    size_t max_iterations, reduction_count;
    double reductions[5];
} ChemBackwardEuler;
typedef struct {
    int stages, new_function[6];
    size_t max_steps;
    double a[15], c[15], m[6], e[6], alpha[6], gamma[6], order;
    double factor_min, factor_max, rejection_factor, safety;
    double h_min, h_max, h_start;
    ChemStepController controller;
    void *controller_context;
    size_t constraint_init_max_iterations, constraint_init_max_backtracks;
    double constraint_init_tolerance, constraint_init_backtrack_factor, constraint_init_sufficient_decrease;
} ChemRosenbrock;

/* Data and all pointer members must outlive the mechanism. Construction validates
   indices and precomputes the sparse fill pattern. Separate states are independent. */
ChemStatus chem_mechanism_create(const ChemMechanismData *, ChemMechanism **);
void chem_mechanism_free(ChemMechanism *);
/* Callback rates can be replaced by name before updating state parameters. */
ChemRate *chem_rate_by_name(ChemMechanism *, const char *name);
const ChemMechanismData *chem_mechanism_data(const ChemMechanism *);
int chem_is_algebraic(const ChemMechanism *, size_t species);
int chem_constraint_nonzero(const ChemMechanism *, size_t row, size_t column);
size_t chem_jacobian_size(const ChemMechanism *);
ChemStatus chem_state_create(const ChemMechanism *, size_t cells, ChemState **);
void chem_state_free(ChemState *);
/* Species-major and parameter-major arrays, with cells contiguous. */
size_t chem_cell_count(const ChemState *);
double *chem_concentrations(ChemState *);
double *chem_parameters(ChemState *);
double *chem_temperatures(ChemState *);
double *chem_pressures(ChemState *);
double *chem_air_densities(ChemState *);
ChemStatus chem_calculate_air_density(ChemState *);
double *chem_absolute_tolerances(ChemState *);
void chem_relative_tolerance(ChemState *, double);
ChemStatus chem_update_rates(ChemState *);
const double *chem_rate_constants(const ChemState *);
void chem_forcing(const ChemState *, const double *y, double *forcing);
void chem_alpha_minus_jacobian(const ChemMechanism *, size_t cells, double *matrix, double alpha);
void chem_negative_jacobian(const ChemState *, const double *y, double *jacobian);
/* Diagnostic row/column to storage index, SIZE_MAX if structural zero. */
size_t chem_jacobian_index(const ChemMechanism *, size_t row, size_t column);
void chem_constraint_forcing(const ChemState *, const double *y, double *forcing);
void chem_constraint_negative_jacobian(const ChemState *, const double *y, double *jacobian);
ChemConstraint chem_constraint_default(ChemConstraintType);
ChemStatus chem_initialize_constraints(ChemState *, const ChemRosenbrock *, ChemStats *);
ChemRosenbrock chem_rosenbrock_parameters(int stages); /* 2, 3, 4, -4 (DA), 6 (DA) */
ChemLinearSolver chem_linear_mozart(void);
ChemLinearSolver chem_linear_doolittle(void);
ChemStatus chem_set_linear_solver(ChemState *, ChemLinearSolver);
ChemBackwardEuler chem_backward_euler_parameters(void);
int chem_backward_euler_converged(const ChemState *, const ChemBackwardEuler *, const double *residual, const double *y);
ChemStatus chem_solve_backward_euler(ChemState *, const ChemBackwardEuler *, double dt, ChemStats *);
double chem_normalized_error(const ChemState *, const double *old, const double *candidate, const double *error);
ChemStatus chem_solve(ChemState *, const ChemRosenbrock *, double dt, ChemStats *);
const char *chem_status_string(ChemStatus);
ChemRate chem_rate_default(ChemRateType);
double chem_rate_evaluate(const ChemRate *, double t, double p, double air, const double *custom);
#endif
