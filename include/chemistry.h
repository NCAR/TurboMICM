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

typedef struct {
  const char *name;
  int third_body;
  double molecular_weight;
  double diffusion;
  double initial_concentration;
  int density_parameterized;
  double density_factor;
} chem_species;

int chem_species_is_parameterized(const chem_species *species);
double chem_species_parameterize(const chem_species *species,
                                 double air_density);

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
int chem_matrix_view_get(const chem_matrix_view *view, size_t row,
                         size_t column, double *value);
int chem_matrix_view_set(chem_matrix_view *view, size_t row,
                         size_t column, double value);
int chem_matrix_view_fill(chem_matrix_view *view, double value);
int chem_matrix_view_for_each(chem_matrix_view *view,
                              chem_matrix_element_function function,
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
int chem_matrix_view_axpy(chem_matrix_view *destination,
                          const chem_matrix_view *source, double alpha);
int chem_matrix_view_min(const chem_matrix_view *view, double *value);
int chem_matrix_view_max(const chem_matrix_view *view, double *value);
int chem_matrix_view_sum(const chem_matrix_view *view, double *value);

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
