#ifndef CHEM_KINETICS_H
#define CHEM_KINETICS_H

#include <stddef.h>
#include <stdint.h>
#include "constraint.h"

typedef enum {
  CHEM_ARRHENIUS = 0,
  CHEM_TROE,
  CHEM_TERNARY,
  CHEM_TERNARY_JPL19,
  CHEM_BRANCHED,
  CHEM_TUNNELING,
  CHEM_TAYLOR,
  CHEM_REVERSIBLE,
  CHEM_USER_DEFINED,
  CHEM_SURFACE,
  CHEM_CALLBACK,
  CHEM_RATE_KIND_COUNT
} chem_rate_kind;

typedef double (*chem_rate_callback)(double temperature, double pressure,
                                     double air_density,
                                     const double *custom_parameters,
                                     void *context);

typedef enum {
  CHEM_SPECIES_UNPARAMETERIZED = 0,
  CHEM_SPECIES_CONSTANT,
  CHEM_SPECIES_TEMPERATURE,
  CHEM_SPECIES_AIR_DENSITY
} chem_species_parameter_kind;

typedef enum {
  CHEM_PROPERTY_DOUBLE = 1,
  CHEM_PROPERTY_STRING,
  CHEM_PROPERTY_INTEGER,
  CHEM_PROPERTY_BOOLEAN
} chem_species_property_kind;

typedef union {
  double real_value;
  const char *string_value;
  int integer_value;
  int boolean_value;
} chem_species_property_value;

typedef struct {
  char *name;
  chem_species_property_kind kind;
  chem_species_property_value value;
} chem_species_property;

typedef struct {
  const char *name;
  double molecular_weight;
  double diffusion_coefficient;
  double initial_value;
  double absolute_tolerance;
  chem_species_property *properties;
  size_t property_count;
  size_t property_capacity;
  chem_species_parameter_kind parameter_kind;
  double parameter_c0;
  double parameter_c_temperature;
  double parameter_c_air_density;
  int owns_properties;
} chem_species;

typedef struct {
  int kind;
  unsigned reactant_offset;
  unsigned reactant_count;
  unsigned product_offset;
  unsigned product_count;
  unsigned body_order;
  unsigned coefficient_offset;
  unsigned coefficient_count;
  int custom0;
  int custom1;
  int surface_species;
  double p[16];
  chem_rate_callback callback;
  void *callback_context;
} chem_reaction;

typedef struct {
  const char *name;
  double default_value;
} chem_custom_parameter;

struct chem_specialized_kernels;

typedef struct {
  const char *name;
  size_t species_count;
  size_t reaction_count;
  size_t reactant_count;
  size_t product_count;
  size_t custom_count;
  double default_temperature;
  double default_pressure;
  const chem_species *species;
  const chem_reaction *reactions;
  const uint32_t *reactants;
  const uint32_t *products;
  const double *yields;
  const double *coefficients;
  const chem_custom_parameter *custom_parameters;
  const struct chem_specialized_kernels *kernels;
} chem_mechanism;

typedef struct {
  const chem_mechanism *mechanism;
  size_t cells;
  double *concentrations;
  double *temperature;
  double *pressure;
  double *custom_parameters;
} chem_state;

/* Generated mechanism used by the TS1 benchmark target. */
extern const chem_mechanism chem_model;

typedef struct {
  uint64_t internal_steps;
  uint64_t accepted;
  uint64_t rejected;
  uint64_t function_calls;
  uint64_t jacobian_updates;
  uint64_t decompositions;
  uint64_t solves;
  uint64_t constraint_init_iterations;
} chem_solver_stats;

typedef enum {
  CHEM_SOLVER_CONVERGED = 0,
  CHEM_SOLVER_RUNNING,
  CHEM_SOLVER_FAILED,
  CHEM_SOLVER_CONSTRAINT_INITIALIZATION_FAILED
} chem_solver_state;

const char *chem_solver_state_string(chem_solver_state state);

typedef struct {
  double small;
  size_t max_iterations;
} chem_backward_euler_parameters;

int chem_backward_euler_is_converged(
    const chem_backward_euler_parameters *parameters,
    const double *residuals, const double *new_values,
    const double *absolute_tolerances, double relative_tolerance,
    size_t cell_count, size_t species_count);
int chem_backward_euler_step(
    const chem_mechanism *mechanism, size_t cell_count,
    const double *rate_constants, double *concentrations, double time_step,
    const chem_backward_euler_parameters *parameters,
    double relative_tolerance, const double *absolute_tolerances,
    chem_solver_stats *stats);

typedef struct chem_sparse_pattern {
  size_t n;
  size_t nnz;
  size_t *row_offsets;
  uint32_t *columns;
  size_t *column_offsets;
  uint32_t *rows;
  size_t *diagonal;
  size_t *position;
  double *mass_diagonal;
} chem_sparse_pattern;

/* Optional build-time generated kernels for a fixed mechanism. */
typedef struct chem_specialized_kernels {
  size_t species_count;
  size_t pattern_nonzeros;
  void (*rates)(double temperature, double pressure,
                const double *custom_parameters, double *rate_constants);
  void (*custom_rates)(double temperature, double pressure,
                       const double *custom_parameters,
                       double *rate_constants);
  int (*evaluate)(const double *concentrations, const double *rate_constants,
                  double *forcing, double *negative_jacobian,
                  const chem_sparse_pattern *pattern);
  int (*factor)(const chem_sparse_pattern *pattern, double alpha,
                const double *negative_jacobian, double *lu, void *context);
  void (*solve)(const chem_sparse_pattern *pattern, const double *lu,
                double *right_hand_side, void *context);
} chem_specialized_kernels;

typedef double (*chem_step_controller)(double error, double order,
                                       double safety_factor, void *context);
typedef int (*chem_lu_factor_callback)(const chem_sparse_pattern *pattern,
                                      double alpha,
                                      const double *negative_jacobian,
                                      double *lu, void *context);
typedef void (*chem_lu_solve_callback)(const chem_sparse_pattern *pattern,
                                       const double *lu,
                                       double *right_hand_side,
                                       void *context);

typedef struct {
  size_t stages;
  unsigned char new_function_evaluation[6];
  double a[15];
  double c[15];
  double m[6];
  double e[6];
  double gamma[6];
  double estimator_order;
  double relative_tolerance;
  double default_absolute_tolerance;
  double factor_min;
  double factor_max;
  double rejection_factor_decrease;
  double safety_factor;
  double h_min_fraction;
  double h_start_fraction;
  double round_off;
  size_t max_number_of_steps;
  size_t constraint_init_max_iterations;
  double constraint_init_tolerance;
  size_t constraint_init_max_backtracks;
  double constraint_init_backtrack_factor;
  double constraint_init_sufficient_decrease;
  chem_step_controller controller;
  chem_lu_factor_callback factor;
  chem_lu_solve_callback solve;
  void *callback_context;
} chem_rosenbrock_parameters;

void chem_rosenbrock_default(chem_rosenbrock_parameters *parameters);
void chem_rosenbrock_two_stage(chem_rosenbrock_parameters *parameters);
void chem_rosenbrock_four_stage(chem_rosenbrock_parameters *parameters);
void chem_rosenbrock_four_stage_dae(chem_rosenbrock_parameters *parameters);
void chem_rosenbrock_six_stage_dae(chem_rosenbrock_parameters *parameters);

typedef struct chem_solver chem_solver;

/* Construction caches environment-dependent rate constants and sparse fill. */
chem_solver *chem_solver_create(const chem_mechanism *mechanism,
                                size_t cell_count,
                                const double *temperature,
                                const double *pressure,
                                const double *custom_parameters,
                                const chem_rosenbrock_parameters *parameters,
                                size_t cell_group_size);
void chem_solver_destroy(chem_solver *solver);
/* Integrate one chemistry interval over the persistent cell-major state. */
int chem_solver_step(chem_solver *solver, double *concentrations,
                     double time_step, chem_solver_stats *stats);
/* Refresh CHEM_USER_DEFINED rates; custom_parameters is cell-major [cell][custom]. */
int chem_solver_update_custom_rates(chem_solver *solver,
                                   const double *temperature,
                                   const double *pressure,
                                   const double *custom_parameters);
int chem_solver_set_constraints(chem_solver *solver,
                                const chem_constraint_set *constraints,
                                const double *equilibrium_constants,
                                const double *constraint_temperatures);
/* Constraint constants use cell-major [cell][constraint] order; temperatures are per cell. */
double chem_solver_absolute_tolerance(const chem_solver *solver, size_t species);
double chem_solver_relative_tolerance(const chem_solver *solver);
double chem_normalized_error(const double *old_values,
                             const double *new_values,
                             const double *estimated_errors,
                             const chem_species *species, size_t cell_count,
                             size_t species_count, double relative_tolerance,
                             double default_absolute_tolerance);

double chem_rate_constant(const chem_mechanism *mechanism,
                          const chem_reaction *reaction, double temperature,
                          double pressure, const double *custom_parameters);
double chem_effective_rate_constant(const chem_mechanism *mechanism,
                                    const chem_reaction *reaction,
                                    double temperature, double pressure,
                                    const double *custom_parameters);
int chem_rate_store_layout(const chem_mechanism *mechanism,
                           size_t counts[CHEM_RATE_KIND_COUNT],
                           size_t offsets[CHEM_RATE_KIND_COUNT + 1],
                           size_t *ordered_reaction_indices);
int chem_mechanism_forcing(const chem_mechanism *mechanism,
                           const double *concentrations,
                           const double *rate_constants, double *forcing);
int chem_mechanism_jacobian_dense(const chem_mechanism *mechanism,
                                  const double *concentrations,
                                  const double *rate_constants,
                                  double *jacobian);

int chem_species_init(chem_species *species, const char *name);
int chem_species_copy(chem_species *destination, const chem_species *source);
void chem_species_destroy(chem_species *species);
int chem_species_set_property(chem_species *species, const char *name,
                              chem_species_property_kind kind,
                              chem_species_property_value value);
int chem_species_get_property(const chem_species *species, const char *name,
                              chem_species_property_kind kind, void *value);
int chem_species_set_parameterization(chem_species *species,
                                      chem_species_parameter_kind kind,
                                      double c0, double c_temperature,
                                      double c_air_density);
int chem_species_set_third_body(chem_species *species);
int chem_species_is_parameterized(const chem_species *species);
double chem_species_parameterize(const chem_species *species,
                                 double temperature, double air_density);

int chem_state_init(chem_state *state, const chem_mechanism *mechanism,
                    size_t cell_count);
int chem_state_copy(chem_state *destination, const chem_state *source);
void chem_state_move(chem_state *destination, chem_state *source);
void chem_state_destroy(chem_state *state);
size_t chem_state_species_index(const chem_state *state, const char *name);
size_t chem_state_custom_parameter_index(const chem_state *state,
                                         const char *name);
int chem_state_set_concentration(chem_state *state, const char *species,
                                 size_t cell, double value);
int chem_state_get_concentration(const chem_state *state, const char *species,
                                 size_t cell, double *value);
int chem_state_set_species_values(chem_state *state, const char *species,
                                 const double *values, size_t value_count);
int chem_state_set_cell_values(chem_state *state, size_t cell,
                               const double *values, size_t value_count);
int chem_state_set_custom_parameter(chem_state *state, const char *name,
                                    size_t cell, double value);
int chem_state_set_custom_parameter_values(chem_state *state, const char *name,
                                           const double *values,
                                           size_t value_count);
int chem_state_set_custom_parameters(chem_state *state, const double *values,
                                     size_t cell_count,
                                     size_t parameter_count);

#endif
