#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>
#include <stdint.h>

/* All concentrations are mol m^-3, pressure Pa, temperature K, time s. */
#define CK_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CK_AVOGADRO 6.02214076e23
#define CK_PI 3.14159265358979323846

typedef enum
{
  CK_OK,
  CK_INVALID,
  CK_MEMORY,
  CK_UNKNOWN,
  CK_MAX_STEPS,
  CK_STEP_SMALL,
  CK_NAN,
  CK_INF,
  CK_CONSTRAINT_FAILED,
  CK_SINGULAR,
  CK_UNCONVERGED
} ck_status;
const char *ck_version(void);
int ck_version_major(void);
int ck_version_minor(void);
int ck_version_patch(void);
int ck_version_tweak(void);
const char *ck_status_string(ck_status status);

typedef struct
{
  double temperature, pressure, air_density;
} ck_conditions;
void ck_ideal_air_density(ck_conditions *conditions);
typedef double (*ck_rate_function)(const ck_conditions *, const double *, void *);
typedef enum
{
  CK_ARRHENIUS,
  CK_TROE,
  CK_TERNARY,
  CK_TERNARY_JPL19,
  CK_TUNNELING,
  CK_BRANCH_ALKOXY,
  CK_BRANCH_NITRATE,
  CK_TAYLOR,
  CK_REVERSIBLE,
  CK_USER,
  CK_SURFACE,
  CK_LAMBDA
} ck_rate_type;
typedef struct
{
  ck_rate_type type;
  double A, B, C, D, E;
  double k0_A, k0_B, k0_C, k0_D, kinf_A, kinf_B, kinf_C, kinf_D;
  double kint_A, kint_B, kint_C, kint_D, Fc, N;
  double X, Y, a0, n, k0, z;
  double probability, diffusion, molecular_weight, scaling, reverse;
  double coefficients[16];
  size_t ncoefficients;
  const char *label;
  ck_rate_function function;
  void *context;
} ck_rate;
ck_rate ck_rate_defaults(ck_rate_type type);
ck_status ck_rate_prepare(ck_rate *rate);
double ck_rate_calculate(const ck_rate *rate, const ck_conditions *conditions, const double *custom);

typedef enum
{
  CK_PROPERTY_REAL,
  CK_PROPERTY_STRING,
  CK_PROPERTY_INT,
  CK_PROPERTY_BOOL
} ck_property_type;
typedef struct
{
  const char *name;
  ck_property_type type;
  union
  {
    double real;
    const char *string;
    int integer;
    int boolean;
  } value;
} ck_property;
typedef struct
{
  const char *name;
  double molecular_weight, diffusion, absolute_tolerance, initial;
  int parameterized; /* 1 = ideal air density; 2 = callback */
  ck_rate_function function;
  void *context;
  size_t nproperties;
  const ck_property *properties;
} ck_species;
ck_status ck_species_property(const ck_species *species, const char *name, ck_property_type type,
                              void *value);
double ck_species_parameter(const ck_species *species, const ck_conditions *conditions);
typedef struct
{
  const char *name;
  size_t nspecies;
  const ck_species *species;
} ck_phase;
size_t ck_phase_state_size(const ck_phase *phase);
ck_status ck_phase_unique_name(const ck_phase *phase, size_t variable, char *buffer, size_t capacity);
typedef struct
{
  size_t species;
  double yield;
} ck_product;
typedef struct
{
  size_t nreactants;
  const size_t *reactants;
  size_t nproducts;
  const ck_product *products;
  ck_rate rate;
} ck_reaction;
typedef enum
{
  CK_LINEAR_CONSTRAINT,
  CK_EQUILIBRIUM_CONSTRAINT
} ck_constraint_type;
typedef struct
{
  ck_constraint_type type;
  const char *name;
  size_t algebraic_species;
  size_t nreactants, nproducts;
  const ck_product *reactants, *products; /* yield is coefficient/exponent */
  double constant, equilibrium_ref, enthalpy, gas_constant, reference_temperature;
} ck_constraint;
ck_status ck_constraint_validate(const ck_constraint *constraint);
typedef struct
{
  const char *name;
  double value;
} ck_named_value;
typedef struct ck_model ck_model;
typedef struct ck_state ck_state;
typedef struct
{
  const char *row, *column;
  int constraint;
} ck_named_jacobian_entry;
typedef struct
{
  const char *name;
  size_t nspecies;
  const ck_species *species;
  size_t ncustom;
  const char *const *custom_names;
  size_t nalgebraic;
  const char *const *algebraic_names;
  size_t njacobian;
  const ck_named_jacobian_entry *jacobian;
  void *context;
  void *(*clone_context)(const void *);
  void (*destroy_context)(void *);
  ck_status (*finalize)(const ck_model *, void *);
  ck_status (*update)(ck_state *, void *);
  void (*forcing)(const ck_state *, const double *, double *, void *);
  void (*jacobian_terms)(const ck_state *, const double *, double *, void *);
  ck_status (*update_constraints)(ck_state *, void *);
  ck_status (*initialize_constraints)(ck_state *, void *);
  void (*constraint_residual)(const ck_state *, const double *, double *, void *);
  void (*constraint_jacobian)(const ck_state *, const double *, double *, void *);
} ck_external_model;

typedef struct
{
  const char *name;
  size_t nspecies, nreactions;
  const ck_species *species;
  const ck_reaction *reactions;
  double temperature, pressure;
  size_t ncustom_defaults;
  const ck_named_value *custom_defaults;
  size_t nconstraints;
  const ck_constraint *constraints;
  size_t nexternal;
  const ck_external_model *external;
} ck_mechanism;

/* The header/species descriptors and all bound numerical structures are copied.
   Strings, property tables, original reaction arrays and initial-default arrays
   are borrowed immutable data and must outlive the model and its states.
   External contexts are borrowed unless clone_context/destroy_context are set.
   A model may be shared by threads; each thread needs an independent state. */
ck_status ck_model_create(const ck_mechanism *description, ck_model **model);
ck_status ck_model_create_with_order(const ck_mechanism *description, int reorder, ck_model **model);
ck_status ck_markowitz_order(size_t n, const unsigned char *pattern, size_t *permutation);
void ck_model_free(ck_model *model);
size_t ck_model_nspecies(const ck_model *model);
size_t ck_model_nreactions(const ck_model *model);
size_t ck_model_ncustom(const ck_model *model);
const char *ck_species_name(const ck_model *model, size_t variable);
const char *ck_custom_name(const ck_model *model, size_t parameter);
ck_status ck_species_index(const ck_model *model, const char *name, size_t *index);
ck_status ck_custom_index(const ck_model *model, const char *name, size_t *index);

/* A state is a group of cells sharing adaptive step size and error norm.
   Arrays are species-major: y[variable * cells + cell]. */
struct ck_state
{
  const ck_model *model;
  size_t cells;
  double *y, *rates, *custom;
  ck_conditions *conditions;
  double relative_tolerance, *absolute_tolerance;
  void *workspace;
};
ck_status ck_state_create(const ck_model *model, size_t cells, ck_state **state);
ck_status ck_state_clone(const ck_state *source, ck_state **destination);
void ck_state_free(ck_state *state);
ck_status ck_state_set(ck_state *state, const char *name, const double *values, size_t count);
ck_status ck_state_set_custom(ck_state *state, const char *name, const double *values, size_t count);
ck_status ck_state_tolerances(ck_state *state, double relative, const double *absolute, size_t count);
typedef struct
{
  const char *name;
  const double *values;
  size_t count;
} ck_named_values;
ck_status ck_state_set_many(ck_state *state, const ck_named_values *values, size_t count, int custom);
ck_status ck_state_set_custom_rows(ck_state *state, const double *values, size_t rows, size_t columns);
ck_status ck_state_move(ck_state **source, ck_state **destination);
ck_status ck_model_set_lambda(ck_model *model, const char *name, ck_rate_function function, void *context);
ck_status ck_model_validate_species_used(const ck_mechanism *description);
ck_status ck_update_rates(ck_state *state);
/* Accumulate chemical reactions only, preserving masked algebraic rows. */
void ck_add_forcing(const ck_state *state, const double *y, double *forcing);
void ck_subtract_jacobian(const ck_state *state, const double *y, double *negative_jacobian);
/* Clear then accumulate chemical and external processes. */
void ck_forcing(const ck_state *state, const double *y, double *forcing);
void ck_jacobian(const ck_state *state, const double *y, double *negative_jacobian);
void ck_constraint_forcing(const ck_state *state, const double *y, double *forcing);
void ck_constraint_jacobian(const ck_state *state, const double *y, double *negative_jacobian);
size_t ck_model_nconstraints(const ck_model *model);
int ck_is_algebraic(const ck_model *model, size_t variable);
size_t ck_jacobian_size(const ck_model *model);
size_t ck_jacobian_index(const ck_model *model, size_t row, size_t column);

/* Sparse symbolic fill is independent of numeric values. Both algorithms use
   unit lower diagonal and retain U's diagonal. No pivoting is introduced. */
typedef enum
{
  CK_LU_MOZART,
  CK_LU_DOOLITTLE
} ck_lu_algorithm;
ck_status ck_linear_model_create(size_t order, const unsigned char *pattern, ck_model **model);
ck_status ck_factor_separate(const ck_model *model, size_t cells, const double *matrix, double *lower,
                             double *upper, ck_lu_algorithm algorithm);
void ck_linear_solve_separate(const ck_model *model, size_t cells, const double *lower, const double *upper,
                              double *rhs);
void ck_alpha_minus_jacobian(const ck_model *model, size_t cells, double *matrix, double alpha);
void ck_matrix_vector(const ck_model *model, size_t cells, const double *matrix, const double *x,
                      double *result);
void ck_factor(const ck_model *model, size_t cells, double *matrix, ck_lu_algorithm algorithm);
void ck_linear_solve(const ck_model *model, size_t cells, const double *lu, double *rhs);

typedef struct
{
  void (*factor)(const ck_model *, size_t, double *, ck_lu_algorithm, void *);
  void (*solve)(const ck_model *, size_t, const double *, double *, void *);
  void *context;
} ck_linear_operations;

typedef struct ck_rosenbrock_parameters
{
  size_t stages, max_steps;
  double round_off, factor_min, factor_max, rejection_decrease, safety;
  double h_min, h_max, h_start, order;
  const ck_linear_operations *linear;
  double (*step_factor)(const struct ck_rosenbrock_parameters *, double, void *);
  void *controller_context;
  int new_function[6];
  double a[15], c[15], m[6], e[6], alpha[6], gamma[6];
  ck_lu_algorithm lu;
  size_t constraint_init_max_iterations, constraint_init_max_backtracks;
  double constraint_init_tolerance, constraint_init_backtrack_factor, constraint_init_sufficient_decrease;
} ck_rosenbrock_parameters;
/* methods: 2, 3, 4, 43 (four-stage DAE), 6 (six-stage DAE). */
ck_status ck_rosenbrock_parameters_init(int method, ck_rosenbrock_parameters *parameters);
typedef struct
{
  uint64_t internal_steps, accepted, rejected, function_calls, jacobian_updates, decompositions, solves;
  uint64_t constraint_init_iterations;
  double final_time;
} ck_stats;
void ck_stats_add(ck_stats *destination, const ck_stats *source);
double ck_normalized_error(const ck_state *state, const double *y, const double *next, const double *errors);
double ck_next_step_factor(const ck_rosenbrock_parameters *parameters, double error);
typedef struct
{
  double small, h_start;
  const ck_linear_operations *linear;
  size_t max_iterations, nreductions;
  double reductions[16];
  ck_lu_algorithm lu;
} ck_backward_euler_parameters;
ck_backward_euler_parameters ck_backward_euler_defaults(void);
int ck_backward_euler_converged(const ck_state *state, const double *residual, const double *y, double small);
ck_status ck_backward_euler_solve(ck_state *state, double dt, const ck_backward_euler_parameters *parameters,
                                  ck_stats *stats);
ck_status ck_initialize_constraints(ck_state *state, const ck_rosenbrock_parameters *parameters,
                                    ck_stats *stats);
void ck_post_solve_clamp(ck_state *state);
ck_status ck_solve_raw(ck_state *state, double dt, const ck_rosenbrock_parameters *parameters,
                       ck_stats *stats);
ck_status ck_solve(ck_state *state, double dt, const ck_rosenbrock_parameters *parameters, ck_stats *stats);

#endif
