#ifndef CHEM_H
#define CHEM_H

#include "chem_constraints.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  CHEM_RATE_ARRHENIUS = 0,
  CHEM_RATE_TROE,
  CHEM_RATE_TERNARY,
  CHEM_RATE_TERNARY_JPL19,
  CHEM_RATE_BRANCHED,
  CHEM_RATE_TUNNELING,
  CHEM_RATE_TAYLOR,
  CHEM_RATE_REVERSIBLE,
  CHEM_RATE_USER,
  CHEM_RATE_SURFACE,
  CHEM_RATE_CALLBACK
} chem_rate_kind;

enum { CHEM_RATE_PARAMETER_COUNT = 24, CHEM_SPECIES_THIRD_BODY = 1u };
enum { CHEM_ROSENBROCK_MAX_STAGES = 6, CHEM_ROSENBROCK_MAX_COEFFICIENTS = 15 };

typedef enum {
  CHEM_LU_MOZART_IN_PLACE = 0,
  CHEM_LU_DOOLITTLE = 1
} chem_lu_method;

/* Row-major dense LU routines. The in-place form stores the lower multipliers
   below the diagonal and U on and above it; the separate form leaves matrix
   unchanged and writes unit-diagonal L and upper-triangular U. These methods do
   not pivot, require disjoint matrix/output arrays, and return -1 for invalid
   input or a zero/non-finite pivot. The separate form returns -2 on allocation
   failure. */
int chem_lu_factor_dense_in_place(double *combined_lu, size_t order, chem_lu_method method);
int chem_lu_factor_dense(
    const double *matrix,
    size_t order,
    double *lower,
    double *upper,
    chem_lu_method method);
int chem_lu_solve_dense_in_place(const double *combined_lu, size_t order, double *right_hand_side);

/* Coordinate sparse LU utilities. Each matrix entry has one unique (row,
   column) pair. Factorization stores lower multipliers and U in `values`; the
   diagonal is stored as its reciprocal to support direct back substitution.
   The supplied sparsity pattern must contain every fill entry produced by
   elimination. These routines do not pivot. */
int chem_lu_factor_sparse_in_place(
    size_t order, size_t nonzero_count, const size_t *rows, const size_t *columns, double *values);
int chem_lu_solve_sparse_in_place(
    size_t order, size_t nonzero_count, const size_t *rows, const size_t *columns,
    const double *factor_values, double *right_hand_side);
int chem_sparse_matrix_vector(
    size_t row_count, size_t column_count, size_t nonzero_count,
    const size_t *rows, const size_t *columns, const double *values,
    const double *vector, double *result);

typedef struct {
  const char *name;
  double molecular_weight;
  double diffusion_coefficient;
  double initial_value;
  double absolute_tolerance;
  unsigned flags;
} chem_species;

/* Owning species metadata used when assembling mechanisms through C code.
   Initialize owning structs with {0}; names and properties are copied by all
   init/copy functions. */
typedef enum {
  CHEM_PROPERTY_STRING = 0,
  CHEM_PROPERTY_DOUBLE,
  CHEM_PROPERTY_INTEGER,
  CHEM_PROPERTY_BOOLEAN
} chem_property_type;

typedef struct {
  double constant;
  double temperature;
  double pressure;
  double air_density;
  int enabled;
} chem_species_parameterization;

typedef struct {
  char *name;
  chem_property_type type;
  union {
    char *string_value;
    double double_value;
    int integer_value;
    int boolean_value;
  } value;
} chem_species_property;

typedef struct {
  char *name;
  size_t property_count;
  chem_species_property *properties;
  chem_species_parameterization parameterization;
} chem_species_info;

typedef struct {
  chem_species_info species;
  double diffusion_coefficient;
  int has_diffusion_coefficient;
  double density;
  int has_density;
} chem_phase_species;

typedef struct {
  char *name;
  size_t species_count;
  chem_phase_species *species;
} chem_phase;

typedef struct {
  char **items;
  size_t count;
} chem_string_list;

typedef struct {
  chem_phase gas_phase;
} chem_system;

/* C equivalents for MICM's species, phase, and system value objects. Return
   0 on success, -1 for invalid input or missing properties, -2 for a property
   type mismatch, and -3 for allocation failure. */
int chem_species_info_init(chem_species_info *species, const char *name);
void chem_species_info_destroy(chem_species_info *species);
int chem_species_info_copy(chem_species_info *destination, const chem_species_info *source);
int chem_species_info_set_string(chem_species_info *species, const char *name, const char *value);
int chem_species_info_set_double(chem_species_info *species, const char *name, double value);
int chem_species_info_set_integer(chem_species_info *species, const char *name, int value);
int chem_species_info_set_boolean(chem_species_info *species, const char *name, int value);
int chem_species_info_get_property(
    const chem_species_info *species,
    const char *name,
    chem_property_type type,
    void *value);
int chem_species_info_set_parameterization(
    chem_species_info *species,
    chem_species_parameterization parameterization);
void chem_species_info_set_third_body(chem_species_info *species);
int chem_species_info_is_parameterized(const chem_species_info *species);
double chem_species_info_parameterized_value(
    const chem_species_info *species,
    double temperature,
    double pressure,
    double air_density);

int chem_phase_species_init(
    chem_phase_species *phase_species,
    const chem_species_info *species,
    const double *diffusion_coefficient,
    const double *density);
int chem_phase_species_set_diffusion_coefficient(chem_phase_species *phase_species, double value);
int chem_phase_species_set_density(chem_phase_species *phase_species, double value);
void chem_phase_species_destroy(chem_phase_species *phase_species);
int chem_phase_init(chem_phase *phase, const char *name, const chem_phase_species *species, size_t species_count);
void chem_phase_destroy(chem_phase *phase);
size_t chem_phase_state_size(const chem_phase *phase);
int chem_phase_unique_names(const chem_phase *phase, chem_string_list *names);
int chem_system_init(chem_system *system, const chem_phase *gas_phase);
void chem_system_destroy(chem_system *system);
size_t chem_system_state_size(const chem_system *system);
typedef size_t (*chem_name_reorder)(const char *const *names, size_t count, size_t output_index, void *context);
int chem_system_unique_names(
    const chem_system *system,
    chem_name_reorder reorder,
    void *context,
    chem_string_list *names);
void chem_string_list_destroy(chem_string_list *list);

typedef struct {
  unsigned species;
  double coefficient;
} chem_term;

/* Callback rates receive cell conditions and caller-owned context. The context
   must outlive the mechanism and be thread-safe for parallel integrations. */
typedef double (*chem_rate_callback)(double temperature, double pressure, double air_density, void *context);

typedef struct {
  chem_rate_kind kind;
  const char *label;
  double p[CHEM_RATE_PARAMETER_COUNT];
  const double *coefficients;
  size_t coefficient_count;
  unsigned custom_parameter;
  unsigned surface_species;
  unsigned air_density_power;
  unsigned branch;
  chem_rate_callback callback;
  void *context;
} chem_rate;

typedef struct {
  unsigned reactant_offset;
  unsigned reactant_count;
  unsigned product_offset;
  unsigned product_count;
  chem_rate rate;
} chem_reaction;

typedef struct {
  const char *name;
  double initial_value;
} chem_parameter;

typedef struct {
  const char *name;
  size_t species_count;
  size_t reaction_count;
  size_t reactant_term_count;
  size_t product_term_count;
  size_t parameter_count;
  const chem_species *species;
  const chem_reaction *reactions;
  const chem_term *reactant_terms;
  const chem_term *product_terms;
  const chem_parameter *parameters;
  double temperature;
  double pressure;
  double relative_tolerance;
} chem_mechanism;

typedef struct {
  unsigned long internal_steps;
  unsigned long accepted;
  unsigned long rejected;
  unsigned long function_calls;
  unsigned long jacobian_updates;
  unsigned long decompositions;
  unsigned long solves;
} chem_statistics;

/* Mutable, cell-major state storage. The state owns its name copies and arrays;
   initialize an empty instance with {0} before calling init/copy/destroy. */
typedef struct {
  size_t cell_count;
  size_t species_count;
  size_t parameter_count;
  size_t rate_constant_count;
  char **species_names;
  char **parameter_names;
  double *concentrations;
  double *custom_parameters;
  double *rate_constants;
  double relative_tolerance;
  double *absolute_tolerances;
} chem_state;

typedef struct {
  const char *name;
  const double *values;
  size_t value_count;
} chem_named_values;

/* Return 0 on success, -1 for invalid dimensions/names, and -2 on allocation
   failure. Scalar concentration/parameter setters require a one-cell state;
   vector setters take one value per cell. */
int chem_state_init(
    chem_state *state,
    const char *const *species_names,
    size_t species_count,
    size_t rate_constant_count,
    const char *const *parameter_names,
    size_t parameter_count,
    size_t cell_count,
    double relative_tolerance,
    const double *absolute_tolerances);
void chem_state_destroy(chem_state *state);
int chem_state_copy(chem_state *destination, const chem_state *source);
void chem_state_move(chem_state *destination, chem_state *source);
int chem_state_set_concentration(chem_state *state, const char *species, double value);
int chem_state_set_concentrations(chem_state *state, const char *species, const double *values, size_t value_count);
int chem_state_set_all_concentrations(chem_state *state, const chem_named_values *values, size_t value_count);
int chem_state_get_concentration(const chem_state *state, const char *species, size_t cell, double *value);
int chem_state_set_custom_parameter(chem_state *state, const char *parameter, double value);
int chem_state_set_custom_parameters(chem_state *state, const char *parameter, const double *values, size_t value_count);
int chem_state_set_all_custom_parameters(chem_state *state, const chem_named_values *values, size_t value_count);
int chem_state_set_custom_parameters_raw(chem_state *state, const double *values, size_t rows, size_t columns);

typedef struct {
  size_t stage_count;
  unsigned char new_function_evaluation[CHEM_ROSENBROCK_MAX_STAGES];
  /* Strictly lower-triangular row storage: offset i*(i-1)/2 + j for stage i, j<i. */
  double a[CHEM_ROSENBROCK_MAX_COEFFICIENTS];
  double c[CHEM_ROSENBROCK_MAX_COEFFICIENTS];
  double m[CHEM_ROSENBROCK_MAX_STAGES];
  double e[CHEM_ROSENBROCK_MAX_STAGES];
  double gamma0;
  chem_lu_method lu_method;
  double estimator_order;
  double round_off;
  double factor_min;
  double factor_max;
  double rejection_factor_decrease;
  double safety_factor;
  /* Step limits are fractions of the requested interval. */
  double h_start_fraction;
  double h_min_fraction;
  double h_max_fraction;
  size_t max_internal_steps;
  size_t constraint_init_max_updates;
  double constraint_init_tolerance;
  size_t constraint_init_max_backtracks;
  double constraint_init_backtrack_factor;
  double constraint_init_sufficient_decrease;
} chem_rosenbrock_parameters;

const chem_mechanism *chem_generated_mechanism(void);

/* Default MICM three-stage Rosenbrock coefficients and error controls. */
chem_rosenbrock_parameters chem_rosenbrock_default_parameters(void);
chem_rosenbrock_parameters chem_rosenbrock_two_stage_parameters(void);
chem_rosenbrock_parameters chem_rosenbrock_four_stage_parameters(void);
chem_rosenbrock_parameters chem_rosenbrock_four_stage_da_parameters(void);
chem_rosenbrock_parameters chem_rosenbrock_six_stage_da_parameters(void);

/* Evaluate one reaction's rate coefficient for one cell. Parameters are in
   descriptor order; pass NULL when the mechanism has no custom parameters.
   Returns 0 on success, -1 for invalid input, and -2 when a surface-rate
   species lacks its required molecular weight or diffusion coefficient. */
int chem_evaluate_rate(
    const chem_mechanism *mechanism,
    size_t reaction_index,
    const double *parameters,
    double temperature,
    double pressure,
    double *rate_constant);

/* Apply already evaluated reaction rates to one or more cell-major states.
   y and forcing use [cell * species_count + species]; rate_constants use
   [cell * reaction_count + reaction]. Input and output arrays must not overlap.
   Return 0 on success, -1 for invalid descriptors or dimensions, and -2 on
   allocation failure. */
int chem_evaluate_forcing(
    const chem_mechanism *mechanism,
    size_t cell_count,
    const double *y,
    const double *rate_constants,
    double *forcing);

/* Add reaction forcing to an existing cell-major buffer. A non-NULL mask has
   one byte per species; zero entries suppress updates to those species rows,
   while masked concentrations still participate in reaction rates. Inputs
   and output must not overlap. Return 0 on success, -1 for invalid descriptors
   or dimensions, and -2 on allocation failure. */
int chem_accumulate_forcing(
    const chem_mechanism *mechanism,
    size_t cell_count,
    const double *y,
    const double *rate_constants,
    const unsigned char *differential_species,
    double *forcing);

/* Add -df/dy to a dense cell-major Jacobian in [cell][row][column] order.
   A row mask follows chem_accumulate_forcing semantics. Inputs and output must
   not overlap. Return 0 on success and -1 for invalid descriptors or
   dimensions. */
int chem_accumulate_negative_jacobian(
    const chem_mechanism *mechanism,
    size_t cell_count,
    const double *y,
    const double *rate_constants,
    const unsigned char *differential_species,
    double *negative_jacobian);

/* Integrate each row of y in place for one chemistry interval.  State layout
   is cell-major: y[cell * species_count + species].  Parameters use the same
   cell-major convention.  A separate call owns a separate adaptive history. */
int chem_integrate_group(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    chem_statistics *statistics);

int chem_integrate_group_steps(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    size_t interval_count,
    chem_statistics *statistics);

int chem_integrate_group_with_parameters(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    const chem_rosenbrock_parameters *solver_parameters,
    chem_statistics *statistics);

int chem_integrate_group_steps_with_parameters(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    size_t interval_count,
    const chem_rosenbrock_parameters *solver_parameters,
    chem_statistics *statistics);

/* Integrate a constrained DAE system with algebraic rows supplied by a
   chem_constraint_set. constraint_parameters uses [cell][state_parameter]
   (for example, precomputed K_eq values); the mechanism's rate parameters
   retain their usual cell-major layout. */
int chem_integrate_group_with_constraints(
    const chem_mechanism *mechanism,
    const chem_constraint_set *constraints,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *constraint_parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    const chem_rosenbrock_parameters *solver_parameters,
    chem_statistics *statistics);

#ifdef __cplusplus
}
#endif

#endif
