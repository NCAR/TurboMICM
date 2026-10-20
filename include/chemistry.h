#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>

#define CHEM_VERSION_MAJOR 0
#define CHEM_VERSION_MINOR 1
#define CHEM_VERSION_PATCH 0
#define CHEM_VERSION_TWEAK 0

typedef enum {
  CHEM_ARRHENIUS,
  CHEM_PHOTOLYSIS,
  CHEM_USER_DEFINED,
  CHEM_TROE,
  CHEM_SURFACE,
  CHEM_TERNARY_ACTIVATION,
  CHEM_TUNNELING,
  CHEM_BRANCHED_ALKOXY,
  CHEM_BRANCHED_NITRATE,
  CHEM_TERNARY_ACTIVATION_JPL19,
  CHEM_TAYLOR_SERIES,
  CHEM_REVERSIBLE,
  CHEM_LAMBDA
} chem_rate_kind;

typedef double (*chem_rate_function)(double temperature, double pressure,
                                     double air_density,
                                     const double *rate_parameters);

const char *chem_version(void);
int chem_version_major(void);
int chem_version_minor(void);
int chem_version_patch(void);
int chem_version_tweak(void);

typedef struct {
  size_t species;
  double coefficient;
} chem_term;

typedef struct {
  chem_rate_kind kind;
  const char *label;
  const chem_term *reactants;
  size_t reactant_count;
  const chem_term *products;
  size_t product_count;
  size_t surface_species;
  double p[24];
  chem_rate_function rate_function;
} chem_reaction;

typedef enum {
  CHEM_PROPERTY_DOUBLE,
  CHEM_PROPERTY_INTEGER,
  CHEM_PROPERTY_BOOLEAN,
  CHEM_PROPERTY_STRING
} chem_property_kind;

typedef union {
  double real_value;
  int integer_value;
  int boolean_value;
  const char *string_value;
} chem_property_value;

typedef struct {
  const char *name;
  chem_property_kind kind;
  chem_property_value value;
} chem_property;

typedef struct {
  const char *name;
  int third_body;
  double molecular_weight;
  double diffusion;
  double initial_concentration;
  int density_parameterized;
  double density_factor;
  const chem_property *properties;
  size_t property_count;
  double temperature_parameterization_factor;
} chem_species;

typedef enum {
  CHEM_PROPERTY_STATUS_OK,
  CHEM_PROPERTY_STATUS_NOT_FOUND,
  CHEM_PROPERTY_STATUS_TYPE_MISMATCH,
  CHEM_PROPERTY_STATUS_INVALID
} chem_property_status;

typedef struct {
  const chem_species *species;
  int has_diffusion_coefficient;
  double diffusion_coefficient;
} chem_phase_species;

typedef struct {
  const char *name;
  const chem_phase_species *species;
  size_t species_count;
} chem_phase;

typedef struct {
  const chem_phase *gas_phase;
} chem_system;

int chem_species_is_parameterized(const chem_species *species);
double chem_species_parameterize(const chem_species *species,
                                 double air_density);
double chem_species_parameterize_conditions(const chem_species *species,
                                           double temperature,
                                           double air_density);
int chem_species_set_third_body(chem_species *species);
chem_property_status chem_species_get_double_property(
    const chem_species *species, const char *name, double *value);
chem_property_status chem_species_get_integer_property(
    const chem_species *species, const char *name, int *value);
chem_property_status chem_species_get_boolean_property(
    const chem_species *species, const char *name, int *value);
chem_property_status chem_species_get_string_property(
    const chem_species *species, const char *name, const char **value);
int chem_phase_species_init(chem_phase_species *phase_species,
                            const chem_species *species,
                            int has_diffusion_coefficient,
                            double diffusion_coefficient);
int chem_phase_species_set_diffusion_coefficient(
    chem_phase_species *phase_species, double diffusion_coefficient);
size_t chem_phase_state_size(const chem_phase *phase);
int chem_phase_unique_name(const chem_phase *phase, size_t species_index,
                           char *buffer, size_t buffer_size);
int chem_system_init(chem_system *system, const chem_phase *gas_phase);
size_t chem_system_state_size(const chem_system *system);
int chem_system_unique_name(const chem_system *system, size_t species_index,
                            char *buffer, size_t buffer_size);

/* Non-owning, cell-major state view. The caller owns the names and arrays. */
typedef struct {
  const char *const *species_names;
  size_t species_count;
  const char *const *custom_parameter_labels;
  size_t custom_parameter_count;
  size_t cell_count;
  double *concentrations;
  double *custom_parameters;
} chem_state_view;

typedef struct {
  const char *label;
  const double *values;
  size_t value_count;
} chem_state_series;

int chem_state_view_init(chem_state_view *state,
                         const char *const *species_names,
                         size_t species_count,
                         const char *const *custom_parameter_labels,
                         size_t custom_parameter_count, size_t cell_count,
                         double *concentrations, double *custom_parameters);
int chem_state_set_concentration(chem_state_view *state, const char *species,
                                 const double *values, size_t value_count);
int chem_state_get_concentration(const chem_state_view *state,
                                 const char *species, size_t cell,
                                 double *value);
int chem_state_set_concentrations(chem_state_view *state,
                                  const chem_state_series *series,
                                  size_t series_count);
int chem_state_set_custom_parameter(chem_state_view *state, const char *label,
                                    const double *values,
                                    size_t value_count);
int chem_state_set_custom_parameters(chem_state_view *state,
                                     const chem_state_series *series,
                                     size_t series_count);
int chem_state_set_custom_parameter_matrix(chem_state_view *state,
                                           const double *values,
                                           size_t row_count,
                                           size_t column_count);

typedef struct {
  const chem_species *species;
  size_t species_count;
  const chem_reaction *reactions;
  size_t reaction_count;
  const char *const *custom_labels;
  const double *default_custom_parameters;
  size_t custom_parameter_count;
  double initial_temperature;
  double initial_pressure;
} chem_mechanism;

typedef enum {
  CHEM_CONSTRAINT_LINEAR,
  CHEM_CONSTRAINT_EQUILIBRIUM
} chem_constraint_kind;

typedef struct {
  size_t species;
  double coefficient;
} chem_constraint_term;

typedef struct {
  chem_constraint_kind kind;
  const char *name;
  size_t algebraic_species;
  const chem_constraint_term *terms;
  size_t term_count;
  double constant;
  const chem_constraint_term *reactants;
  size_t reactant_count;
  const chem_constraint_term *products;
  size_t product_count;
  double equilibrium_constant_ref;
  double enthalpy_change;
  double reference_temperature;
} chem_constraint;

typedef struct {
  const chem_constraint *constraints;
  size_t constraint_count;
  size_t species_count;
} chem_constraint_set;

int chem_constraint_set_validate(const chem_constraint_set *set);
int chem_constraint_set_algebraic_mask(const chem_constraint_set *set,
                                       unsigned char *mask);
int chem_constraint_set_jacobian_pattern(const chem_constraint_set *set,
                                         unsigned char *pattern);
int chem_constraint_set_add_forcing(const chem_constraint_set *set,
                                    const double *states,
                                    const double *temperatures,
                                    size_t cell_count, double *forcing);
int chem_constraint_set_add_negative_jacobian(
    const chem_constraint_set *set, const double *states,
    const double *temperatures, size_t cell_count,
    double *negative_jacobians);
int chem_forcing_add_with_algebraic_mask(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, double air_density,
    const unsigned char *algebraic_mask, double *forcing);
int chem_negative_jacobian_add_with_algebraic_mask(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, const unsigned char *algebraic_mask,
    double *negative_jacobian);

/* A non-owning rectangular view over row-major or strided matrix storage. */
typedef struct {
  double *data;
  size_t rows;
  size_t columns;
  size_t row_stride;
  size_t column_stride;
} chem_matrix_view;

typedef struct {
  const double *data;
  size_t rows;
  size_t columns;
  size_t row_stride;
  size_t column_stride;
} chem_const_matrix_view;

typedef void (*chem_matrix_element_function)(size_t row, size_t column,
                                              double *value, void *context);
typedef void (*chem_matrix_row_function)(
    size_t row, double *destination, size_t destination_columns,
    size_t destination_stride, const double *first, size_t first_columns,
    size_t first_stride, const double *second, size_t second_columns,
    size_t second_stride, void *context);
typedef int (*chem_matrix_predicate)(double value, void *context);
typedef double (*chem_matrix_combine)(double accumulator, double value,
                                      void *context);

int chem_matrix_view_init(chem_matrix_view *view, double *data,
                          size_t rows, size_t columns,
                          size_t row_stride, size_t column_stride);
int chem_const_matrix_view_init(chem_const_matrix_view *view,
                                const double *data, size_t rows,
                                size_t columns, size_t row_stride,
                                size_t column_stride);
int chem_const_matrix_view_get(const chem_const_matrix_view *view,
                               size_t row, size_t column, double *value);
int chem_const_matrix_view_copy_row_to_vector(
    const chem_const_matrix_view *view, size_t row, double *vector,
    size_t vector_size);

typedef struct {
  double *data;
  size_t rows;
  size_t columns;
  size_t group_width;
} chem_vector_matrix_view;

typedef struct {
  const double *data;
  size_t rows;
  size_t columns;
  size_t group_width;
} chem_const_vector_matrix_view;

int chem_vector_matrix_view_storage_size(size_t rows, size_t columns,
                                         size_t group_width,
                                         size_t *storage_size);
int chem_vector_matrix_view_init(chem_vector_matrix_view *view, double *data,
                                 size_t rows, size_t columns,
                                 size_t group_width);
int chem_const_vector_matrix_view_init(chem_const_vector_matrix_view *view,
                                       const double *data, size_t rows,
                                       size_t columns, size_t group_width);
size_t chem_vector_matrix_view_group_size(const chem_vector_matrix_view *view);
size_t chem_vector_matrix_view_number_of_groups(
    const chem_vector_matrix_view *view);
size_t chem_vector_matrix_view_row_stride(const chem_vector_matrix_view *view);
size_t chem_vector_matrix_view_column_stride(const chem_vector_matrix_view *view);
int chem_vector_matrix_view_fill(chem_vector_matrix_view *view, double value);
int chem_vector_matrix_view_copy(chem_vector_matrix_view *destination,
                                 const chem_vector_matrix_view *source);
int chem_vector_matrix_view_get(const chem_vector_matrix_view *view,
                                size_t row, size_t column, double *value);
int chem_vector_matrix_view_set(chem_vector_matrix_view *view,
                                size_t row, size_t column, double value);
int chem_const_vector_matrix_view_get(
    const chem_const_vector_matrix_view *view, size_t row, size_t column,
    double *value);
int chem_vector_matrix_view_copy_to_vector(const chem_vector_matrix_view *view,
                                           double *vector,
                                           size_t vector_size);
int chem_vector_matrix_view_copy_from_vector(chem_vector_matrix_view *view,
                                             const double *vector,
                                             size_t vector_size);
int chem_const_vector_matrix_view_copy_to_vector(
    const chem_const_vector_matrix_view *view, double *vector,
    size_t vector_size);
int chem_vector_matrix_view_copy_row_to_vector(
    const chem_vector_matrix_view *view, size_t row, double *vector,
    size_t vector_size);
int chem_vector_matrix_view_copy_row_from_vector(
    chem_vector_matrix_view *view, size_t row, const double *vector,
    size_t vector_size);
int chem_vector_matrix_view_copy_column_to_vector(
    const chem_vector_matrix_view *view, size_t column, double *vector,
    size_t vector_size);
int chem_vector_matrix_view_copy_column_from_vector(
    chem_vector_matrix_view *view, size_t column, const double *vector,
    size_t vector_size);
int chem_vector_matrix_view_axpy(chem_vector_matrix_view *destination,
                                 const chem_vector_matrix_view *source,
                                 double alpha);
int chem_vector_matrix_view_for_each(chem_vector_matrix_view *view,
                                    chem_matrix_element_function function,
                                    void *context);
int chem_vector_matrix_view_for_each_row(
    chem_vector_matrix_view *destination,
    const chem_vector_matrix_view *first,
    const chem_vector_matrix_view *second,
    chem_matrix_row_function function, void *context);
int chem_vector_matrix_view_min(const chem_vector_matrix_view *view,
                                double *value);
int chem_vector_matrix_view_max(const chem_vector_matrix_view *view,
                                double *value);
int chem_matrix_view_get(const chem_matrix_view *view, size_t row,
                         size_t column, double *value);
int chem_matrix_view_set(chem_matrix_view *view, size_t row,
                         size_t column, double value);
int chem_matrix_view_fill(chem_matrix_view *view, double value);
int chem_matrix_view_for_each(chem_matrix_view *view,
                              chem_matrix_element_function function,
                              void *context);
int chem_matrix_view_for_each_row(chem_matrix_view *destination,
                                  const chem_matrix_view *first,
                                  const chem_matrix_view *second,
                                  chem_matrix_row_function function,
                                  void *context);
int chem_matrix_view_for_each_row_const_sources(
    chem_matrix_view *destination, const chem_const_matrix_view *first,
    const chem_const_matrix_view *second, chem_matrix_row_function function,
    void *context);
int chem_matrix_view_copy_to_vector(const chem_matrix_view *view,
                                    double *vector, size_t vector_size);
int chem_matrix_view_copy_from_vector(chem_matrix_view *view,
                                      const double *vector,
                                      size_t vector_size);
int chem_matrix_view_copy_row_to_vector(const chem_matrix_view *view,
                                        size_t row, double *vector,
                                        size_t vector_size);
int chem_matrix_view_copy_row_from_vector(chem_matrix_view *view,
                                          size_t row, const double *vector,
                                          size_t vector_size);
int chem_matrix_view_copy_column_to_vector(const chem_matrix_view *view,
                                           size_t column, double *vector,
                                           size_t vector_size);
int chem_matrix_view_copy_column_from_vector(chem_matrix_view *view,
                                             size_t column,
                                             const double *vector,
                                             size_t vector_size);
int chem_matrix_view_axpy(chem_matrix_view *destination,
                          const chem_matrix_view *source, double alpha);
int chem_matrix_view_min(const chem_matrix_view *view, double *value);
int chem_matrix_view_max(const chem_matrix_view *view, double *value);
int chem_matrix_view_sum(const chem_matrix_view *view, double *value);
int chem_matrix_view_reduce_or(const chem_matrix_view *view,
                               chem_matrix_predicate predicate, void *context,
                               int initial, int *value);
int chem_matrix_view_reduce_and(const chem_matrix_view *view,
                                chem_matrix_predicate predicate, void *context,
                                int initial, int *value);
int chem_matrix_view_reduce_strict(const chem_const_matrix_view *view,
                                   chem_matrix_combine combine, void *context,
                                   double initial,
                                   double *value);

/* Non-owning block CSR matrix; values are block-major over the shared pattern. */
typedef struct {
  size_t dimension;
  size_t block_count;
  size_t nonzero_count;
  const size_t *row_offsets;
  const size_t *column_indices;
  double *values;
  int compressed_by_column;
} chem_csr_matrix;

int chem_csr_matrix_init(chem_csr_matrix *matrix, size_t dimension,
                         size_t block_count, size_t nonzero_count,
                         const size_t *row_offsets,
                         const size_t *column_indices, double *values);
int chem_csr_matrix_init_csc(chem_csr_matrix *matrix, size_t dimension,
                             size_t block_count, size_t nonzero_count,
                             const size_t *column_offsets,
                             const size_t *row_indices, double *values);
int chem_csr_matrix_get(const chem_csr_matrix *matrix, size_t block,
                        size_t row, size_t column, double *value);
int chem_csr_matrix_set(chem_csr_matrix *matrix, size_t block, size_t row,
                        size_t column, double value);
int chem_csr_matrix_vector_index(const chem_csr_matrix *matrix, size_t block,
                                 size_t row, size_t column, size_t *index);
int chem_csr_matrix_is_structural_zero(const chem_csr_matrix *matrix,
                                       size_t row, size_t column,
                                       int *is_zero);
int chem_csr_matrix_fill(chem_csr_matrix *matrix, double value);
int chem_csr_matrix_add_to_diagonal(chem_csr_matrix *matrix, double value);
typedef struct {
  double *data;
  size_t count;
  size_t stride;
} chem_csr_value_view;
typedef struct {
  const double *data;
  size_t count;
  size_t stride;
} chem_csr_const_value_view;
typedef void (*chem_csr_block_function)(size_t block, double *destination,
                                       const chem_csr_const_value_view *sources,
                                       size_t source_count, void *context);
int chem_csr_matrix_get_value_view(const chem_csr_matrix *matrix, size_t row,
                                   size_t column,
                                   chem_csr_value_view *view);
int chem_csr_matrix_get_const_value_view(const chem_csr_matrix *matrix,
                                        size_t row, size_t column,
                                        chem_csr_const_value_view *view);
int chem_csr_matrix_apply_blocks(chem_csr_value_view *destination,
                                 const chem_csr_const_value_view *sources,
                                 size_t source_count,
                                 chem_csr_block_function function,
                                 void *context);
int chem_csr_value_view_fill(chem_csr_value_view *view, double value);
int chem_csr_value_view_copy(chem_csr_value_view *destination,
                             const chem_csr_const_value_view *source);
int chem_csr_value_view_copy_to_vector(
    const chem_csr_const_value_view *view, double *vector,
    size_t vector_size);
int chem_csr_value_view_copy_from_vector(chem_csr_value_view *view,
                                         const double *vector,
                                         size_t vector_size);

double chem_arrhenius(double a, double b, double c, double d, double e,
                      double temperature, double air_density);
double chem_reaction_rate(const chem_reaction *reaction, double temperature,
                          double air_density, const double *rate_parameters);
int chem_forcing(const chem_mechanism *mechanism, const double *concentration,
                 double temperature, double air_density,
                 const double *rate_parameters, double *forcing);
int chem_jacobian(const chem_mechanism *mechanism, const double *concentration,
                  double temperature, double air_density,
                  const double *rate_parameters, double *jacobian);
int chem_rate_constants(const chem_mechanism *mechanism, double temperature,
                        double pressure, double air_density,
                        const double *rate_parameters,
                        double *rate_constants);
int chem_forcing_with_constants(const chem_mechanism *mechanism,
                                const double *concentration,
                                const double *rate_constants, double air_density,
                                double *forcing);
int chem_jacobian_with_constants(const chem_mechanism *mechanism,
                                 const double *concentration,
                                 const double *rate_constants,
                                 double air_density, double *jacobian);
int chem_negative_jacobian_with_constants(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, double *negative_jacobian);

typedef struct {
  size_t internal_steps;
  size_t accepted;
  size_t rejected;
  size_t function_calls;
  size_t jacobian_updates;
  size_t decompositions;
  size_t solves;
  size_t constraint_init_iterations;
} chem_solver_stats;

/* Lower-triangular Rosenbrock coefficients use row-wise packed indexing:
   A(i,j) and C(i,j) are at i*(i-1)/2+j for 0 <= j < i < stages. */
typedef struct {
  size_t stages;
  double a[15];
  double c[15];
  double m[6];
  double e[6];
  double gamma[6];
  unsigned char new_function_evaluation[6];
  double estimator_order;
} chem_rosenbrock_parameters;

typedef enum {
  CHEM_ROSENBROCK_THREE_STAGE,
  CHEM_ROSENBROCK_TWO_STAGE,
  CHEM_ROSENBROCK_FOUR_STAGE,
  CHEM_ROSENBROCK_FOUR_STAGE_DAE,
  CHEM_ROSENBROCK_SIX_STAGE_DAE
} chem_rosenbrock_parameter_set;

chem_rosenbrock_parameters chem_rosenbrock_builtin_parameters(
    chem_rosenbrock_parameter_set parameter_set);

typedef struct {
  double relative_tolerance;
  double absolute_tolerance;
  size_t max_steps;
} chem_solver_options;

typedef struct {
  void *implementation;
} chem_solver_workspace;

chem_solver_options chem_solver_default_options(void);
int chem_solver_workspace_init(chem_solver_workspace *workspace,
                              const chem_mechanism *mechanism);
int chem_solver_workspace_init_group(chem_solver_workspace *workspace,
                                    const chem_mechanism *mechanism,
                                    size_t cell_capacity);
void chem_solver_workspace_destroy(chem_solver_workspace *workspace);
int chem_rosenbrock3_solve_with_workspace(
    const chem_mechanism *mechanism, chem_solver_workspace *workspace,
    double *state, double temperature, double pressure, double air_density,
    const double *rate_parameters, double duration,
    const chem_solver_options *options, chem_solver_stats *statistics);
int chem_rosenbrock_solve_with_parameters(
    const chem_mechanism *mechanism, chem_solver_workspace *workspace,
    double *state, double temperature, double pressure, double air_density,
    const double *rate_parameters, double duration,
    const chem_rosenbrock_parameters *parameters,
    const chem_solver_options *options, chem_solver_stats *statistics);
int chem_rosenbrock_solve_constrained_with_parameters(
    const chem_mechanism *mechanism, chem_solver_workspace *workspace,
    double *state, double temperature, double pressure, double air_density,
    const double *rate_parameters, double duration,
    const chem_rosenbrock_parameters *parameters,
    const chem_constraint_set *constraints,
    const chem_solver_options *options, chem_solver_stats *statistics);
int chem_rosenbrock_solve_constrained_group_with_parameters(
    const chem_mechanism *mechanism, chem_solver_workspace *workspace,
    double *states, const double *temperatures, const double *pressures,
    const double *air_densities, const double *rate_parameters,
    size_t cell_count, double duration,
    const chem_rosenbrock_parameters *parameters,
    const chem_constraint_set *constraints,
    const chem_solver_options *options, chem_solver_stats *statistics);
int chem_rosenbrock3_solve_group_with_workspace(
    const chem_mechanism *mechanism, chem_solver_workspace *workspace,
    double *states, const double *temperatures, const double *pressures,
    const double *air_densities, const double *rate_parameters,
    size_t cell_count, double duration, const chem_solver_options *options,
    chem_solver_stats *statistics);
int chem_rosenbrock3_solve(const chem_mechanism *mechanism, double *state,
                           double temperature, double air_density,
                           const double *rate_parameters, double duration,
                           const chem_solver_options *options,
                           chem_solver_stats *statistics);

#endif
