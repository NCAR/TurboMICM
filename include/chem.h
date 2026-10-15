#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>
#include <stdint.h>

/* All concentrations are mol m^-3, pressure Pa, temperature K, time s. */
#define CK_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CK_AVOGADRO 6.02214076e23
#define CK_PI 3.14159265358979323846

typedef enum {
  CK_OK, CK_INVALID, CK_MEMORY, CK_UNKNOWN, CK_MAX_STEPS, CK_STEP_SMALL,
  CK_NAN, CK_INF, CK_CONSTRAINT_FAILED, CK_SINGULAR
} ck_status;
const char *ck_status_string(ck_status status);

typedef struct { double temperature, pressure, air_density; } ck_conditions;
void ck_ideal_air_density(ck_conditions *conditions);
typedef double (*ck_rate_function)(const ck_conditions *, const double *, void *);
typedef enum {
  CK_ARRHENIUS, CK_TROE, CK_TERNARY, CK_TERNARY_JPL19, CK_TUNNELING,
  CK_BRANCH_ALKOXY, CK_BRANCH_NITRATE, CK_TAYLOR, CK_REVERSIBLE,
  CK_USER, CK_SURFACE, CK_LAMBDA
} ck_rate_type;
typedef struct {
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

typedef struct {
  const char *name;
  double molecular_weight, diffusion, absolute_tolerance, initial;
  int parameterized; /* 1 = ideal air density; 2 = callback */
  ck_rate_function function;
  void *context;
} ck_species;
typedef struct { size_t species; double yield; } ck_product;
typedef struct {
  size_t nreactants;
  const size_t *reactants;
  size_t nproducts;
  const ck_product *products;
  ck_rate rate;
} ck_reaction;
typedef struct { const char *name; double value; } ck_named_value;
typedef struct {
  const char *name;
  size_t nspecies, nreactions;
  const ck_species *species;
  const ck_reaction *reactions;
  double temperature, pressure;
  size_t ncustom_defaults;
  const ck_named_value *custom_defaults;
} ck_mechanism;

typedef struct ck_model ck_model;
typedef struct ck_state ck_state;
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
struct ck_state {
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
ck_status ck_update_rates(ck_state *state);
void ck_forcing(const ck_state *state, const double *y, double *forcing);
void ck_jacobian(const ck_state *state, const double *y, double *negative_jacobian);
size_t ck_jacobian_size(const ck_model *model);
size_t ck_jacobian_index(const ck_model *model, size_t row, size_t column);

/* Sparse symbolic fill is independent of numeric values. Both algorithms use
   unit lower diagonal and retain U's diagonal. No pivoting is introduced. */
typedef enum { CK_LU_MOZART, CK_LU_DOOLITTLE } ck_lu_algorithm;
void ck_factor(const ck_model *model, size_t cells, double *matrix, ck_lu_algorithm algorithm);
void ck_linear_solve(const ck_model *model, size_t cells, const double *lu, double *rhs);

typedef struct {
  size_t stages, max_steps;
  double round_off, factor_min, factor_max, rejection_decrease, safety;
  double h_min, h_max, h_start, order;
  int new_function[6];
  double a[15], c[15], m[6], e[6], alpha[6], gamma[6];
  ck_lu_algorithm lu;
} ck_rosenbrock_parameters;
/* methods: 2, 3, 4, 43 (four-stage DAE), 6 (six-stage DAE). */
ck_status ck_rosenbrock_parameters_init(int method, ck_rosenbrock_parameters *parameters);
typedef struct {
  uint64_t internal_steps, accepted, rejected, function_calls, jacobian_updates, decompositions, solves;
  double final_time;
} ck_stats;
void ck_stats_add(ck_stats *destination, const ck_stats *source);
double ck_normalized_error(const ck_state *state, const double *y, const double *next, const double *errors);
double ck_next_step_factor(const ck_rosenbrock_parameters *parameters, double error);
ck_status ck_solve(ck_state *state, double dt, const ck_rosenbrock_parameters *parameters, ck_stats *stats);

#endif
