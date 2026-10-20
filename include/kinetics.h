#ifndef CHEM_KINETICS_H
#define CHEM_KINETICS_H

#include <stddef.h>
#include <stdint.h>

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
  CHEM_CALLBACK
} chem_rate_kind;

typedef double (*chem_rate_callback)(double temperature, double pressure,
                                     double air_density,
                                     const double *custom_parameters,
                                     void *context);

typedef struct {
  const char *name;
  double molecular_weight;
  double diffusion_coefficient;
  double initial_value;
  double absolute_tolerance;
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
} chem_solver_stats;

typedef struct chem_sparse_pattern {
  size_t n;
  size_t nnz;
  size_t *row_offsets;
  uint32_t *columns;
  size_t *column_offsets;
  uint32_t *rows;
  size_t *diagonal;
  size_t *position;
} chem_sparse_pattern;

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
  chem_step_controller controller;
  chem_lu_factor_callback factor;
  chem_lu_solve_callback solve;
  void *callback_context;
} chem_rosenbrock_parameters;

void chem_rosenbrock_default(chem_rosenbrock_parameters *parameters);

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

double chem_rate_constant(const chem_mechanism *mechanism,
                          const chem_reaction *reaction, double temperature,
                          double pressure, const double *custom_parameters);

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
