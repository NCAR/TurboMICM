#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>
#include <stdint.h>

#define CHEM_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CHEM_PI 3.14159265358979323846
#define CHEM_NONE SIZE_MAX

typedef enum {
    CHEM_OK, CHEM_INVALID, CHEM_MEMORY, CHEM_NOT_FOUND, CHEM_SIZE,
    CHEM_SINGULAR, CHEM_MAX_STEPS, CHEM_SMALL_STEP, CHEM_NAN, CHEM_INF,
    CHEM_CONSTRAINT_INIT, CHEM_UNCONVERGED
} chem_status;
const char *chem_status_string(chem_status status);

typedef struct { double temperature, pressure, air_density; } chem_conditions;
void chem_ideal_air_density(chem_conditions *conditions);

typedef enum {
    CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_TERNARY_JPL19,
    CHEM_TUNNELING, CHEM_BRANCHED, CHEM_TAYLOR, CHEM_REVERSIBLE,
    CHEM_USER, CHEM_SURFACE, CHEM_LAMBDA
} chem_rate_type;
typedef double (*chem_rate_function)(const chem_conditions *, const double *, void *);
typedef struct { double A, B, C, D, E; } chem_arrhenius;
typedef struct { double k0_A, k0_B, k0_C, kinf_A, kinf_B, kinf_C, Fc, N; } chem_falloff;
typedef struct {
    double k0_A, k0_B, k0_C, k0_D, kinf_A, kinf_B, kinf_C, kinf_D;
    double kint_A, kint_B, kint_C, kint_D, Fc, N;
} chem_jpl19;
typedef struct { double A, B, C; } chem_tunneling;
typedef struct { int nitrate; double X, Y, a0; size_t n; double k0, z; } chem_branched;
typedef struct { chem_arrhenius arrhenius; size_t n; double coefficients[16]; } chem_taylor;
typedef struct { double A, C, k_r; } chem_reversible;
typedef struct { double scaling; } chem_user;
typedef struct { double diffusion, molecular_weight, probability, speed_factor; } chem_surface;
typedef struct { chem_rate_function function; void *context; } chem_lambda;
typedef struct {
    chem_rate_type type;
    const char *label;
    union {
        chem_arrhenius arrhenius;
        chem_falloff falloff;
        chem_jpl19 jpl19;
        chem_tunneling tunneling;
        chem_branched branched;
        chem_taylor taylor;
        chem_reversible reversible;
        chem_user user;
        chem_surface surface;
        chem_lambda lambda;
    } p;
} chem_rate;
chem_rate chem_rate_default(chem_rate_type type);
chem_status chem_rate_prepare(chem_rate *rate);
double chem_rate_calculate(const chem_rate *rate, const chem_conditions *conditions, const double *custom);

typedef enum { CHEM_SPECIES_VARIABLE, CHEM_SPECIES_THIRD_BODY, CHEM_SPECIES_PARAMETER } chem_species_type;
typedef struct {
    const char *name;
    chem_species_type type;
    double absolute_tolerance, initial;
    double molecular_weight, diffusion;
    chem_rate_function parameter;
    void *context;
} chem_species;
typedef struct { size_t species; double yield; } chem_product;
typedef struct {
    size_t n_reactants, n_products;
    const size_t *reactants;
    const chem_product *products;
    chem_rate rate;
} chem_reaction;
typedef struct { const char *name; double value; } chem_named_value;

struct chem_model;
/* Kernels add to the destination. Columns are contiguous arrays of cells. */
typedef void (*chem_kernel)(const struct chem_model *, size_t, const double *, const double *, double *);
typedef struct {
    size_t n_species, n_reactions;
    const chem_species *species;
    const chem_reaction *reactions;
    size_t n_defaults;
    const chem_named_value *defaults;
    chem_conditions conditions;
    chem_kernel forcing, jacobian;
    int preserve_order;
} chem_mechanism;

typedef enum { CHEM_LU_MOZART, CHEM_LU_DOOLITTLE } chem_lu_algorithm;
typedef struct { size_t target, left, right; unsigned kind; } chem_lu_op;
typedef struct {
    size_t n, nnz, n_ops;
    size_t *map, *row_start, *column, *diagonal;
    chem_lu_op *ops;
    chem_lu_algorithm algorithm;
} chem_lu_plan;
chem_status chem_lu_create(chem_lu_plan *plan, size_t n, const unsigned char *pattern, chem_lu_algorithm algorithm);
void chem_lu_destroy(chem_lu_plan *plan);
chem_status chem_markowitz(size_t n, const unsigned char *pattern, size_t *permutation);
void chem_lu_factor(const chem_lu_plan *plan, size_t cells, double *matrix, double *scratch);
void chem_lu_solve(const chem_lu_plan *plan, size_t cells, const double *matrix, double *rhs);

typedef struct {
    size_t reaction, independent, n_reactants, n_products;
    size_t *reactants, *reactant_flat, *product_flat;
    double *yields;
    size_t diagonal;
} chem_jac_job;
typedef struct chem_model {
    const chem_mechanism *mechanism; /* Borrowed immutable description; must outlive model. */
    size_t n_species, n_reactions, n_custom, n_jobs;
    size_t *spec_to_var, *var_to_spec, *custom_offset;
    char **custom_names;
    chem_rate *rates;
    chem_jac_job *jobs;
    chem_lu_plan lu;
    chem_kernel forcing, jacobian;
    unsigned char *algebraic;
} chem_model;
chem_status chem_model_create(chem_model *model, const chem_mechanism *mechanism, chem_lu_algorithm algorithm);
void chem_model_destroy(chem_model *model);
size_t chem_species_index(const chem_model *model, const char *name);
size_t chem_custom_index(const chem_model *model, const char *name);
void chem_forcing(const chem_model *model, size_t cells, const double *rates, const double *y, double *forcing);
void chem_jacobian(const chem_model *model, size_t cells, const double *rates, const double *y, double *negative_jacobian);

typedef struct {
    const chem_model *model;
    size_t cells;
    double *y, *rates, *custom, *atol;
    chem_conditions *conditions;
    double rtol;
    double *work, *jacobian, *scratch;
} chem_state;
chem_status chem_state_create(chem_state *state, const chem_model *model, size_t cells);
void chem_state_destroy(chem_state *state);
chem_status chem_state_copy(chem_state *destination, const chem_state *source);
chem_status chem_state_set_species(chem_state *state, const char *name, const double *values, size_t count);
chem_status chem_state_set_custom(chem_state *state, const char *name, const double *values, size_t count);
void chem_calculate_rates(chem_state *state);

typedef struct {
    uint64_t internal_steps, accepted, rejected, function_calls, jacobian_updates, decompositions, solves;
    double final_time;
} chem_stats;
void chem_stats_add(chem_stats *total, const chem_stats *step);
typedef struct { chem_status status; chem_stats stats; } chem_result;
typedef enum { CHEM_ROS2, CHEM_ROS3, CHEM_ROS4, CHEM_RODA3, CHEM_RODA4 } chem_ros_method;
typedef struct chem_ros_parameters {
    size_t stages, max_steps;
    double round_off, factor_min, factor_max, rejection_factor, safety_factor;
    double h_min, h_max, h_start, order;
    size_t init_max_iterations, init_max_backtracks;
    double init_tolerance, init_backtrack_factor, init_sufficient_decrease;
    int new_function[6];
    double a[15], c[15], m[6], e[6], alpha[6], gamma[6];
} chem_ros_parameters;
chem_ros_parameters chem_ros_parameters_default(chem_ros_method method);
/* Return an unbounded proposed next h; acceptance/rejection policy stays in the solver. */
typedef double (*chem_step_controller)(double h, double error, const chem_ros_parameters *, void *);
double chem_ros_controller(double h, double error, const chem_ros_parameters *parameters, void *context);
double chem_normalized_error(const chem_state *state, const double *y, const double *ynew, const double *error);
chem_result chem_ros_solve(chem_state *state, double dt, const chem_ros_parameters *parameters,
                           chem_step_controller controller, void *context);

typedef struct {
    double small, h_start;
    size_t max_steps;
    double reductions[5];
} chem_be_parameters;
chem_be_parameters chem_be_parameters_default(void);
int chem_be_converged(const chem_state *state, const chem_be_parameters *parameters,
                      const double *residual, const double *y);
chem_result chem_be_solve(chem_state *state, double dt, const chem_be_parameters *parameters);
chem_status chem_model_set_lambda(chem_model *model, const char *name, chem_rate_function function, void *context);

typedef struct { chem_status code; char message[256]; } chem_error;
chem_status chem_model_create_checked(chem_model *model, const chem_mechanism *mechanism,
                                      chem_lu_algorithm algorithm, chem_error *error);

#endif
