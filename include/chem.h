/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>
#include <stdio.h>

#define CH_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CH_AVOGADRO 6.02214076e23

typedef enum {
  CH_OK = 0, CH_INVALID = -1, CH_MEMORY = -2, CH_NOT_FOUND = -3,
  CH_SIZE = -4, CH_FROZEN = -5, CH_SINGULAR = -6,
  CH_MAX_STEPS = -7, CH_SMALL_STEP = -8, CH_NAN = -9, CH_INF = -10,
  CH_CONSTRAINT_INIT = -11
} ch_status;
const char *ch_status_string(int status);

typedef struct { double temperature, pressure, air_density; } ch_conditions;
void ch_ideal_air_density(ch_conditions *conditions);

typedef enum {
  CH_ARRHENIUS, CH_TROE, CH_TERNARY, CH_TERNARY_JPL19, CH_TUNNELING,
  CH_BRANCHED, CH_SURFACE, CH_USER, CH_TAYLOR, CH_REVERSIBLE, CH_LAMBDA
} ch_rate_type;
typedef double (*ch_rate_function)(const ch_conditions *, const double *, void *);

/* Parameter order is documented in DESIGN.md. Initialize with ch_rate_default. */
typedef struct {
  int type;
  double p[32];
  int custom_index, coefficient_count;
  ch_rate_function function;
  void *context;
} ch_rate;
ch_rate ch_rate_default(int type);
void ch_rate_prepare(ch_rate *rate);
double ch_rate_calculate(const ch_rate *, const ch_conditions *, const double *custom);

typedef struct ch_model ch_model;
typedef struct ch_state ch_state;

ch_model *ch_model_create(void);
void ch_model_free(ch_model *);
/* Species IDs are stable builder IDs; state indices are assigned at finalize. */
int ch_model_add_species(ch_model *, const char *name, double atol, int third_body);
int ch_model_set_species_property(ch_model *, int species, const char *key, double value);
int ch_model_add_parameter(ch_model *, const char *name, double initial_value);
int ch_model_add_reaction(ch_model *, size_t nr, const int *reactants,
                          size_t np, const int *products, const double *yields,
                          const ch_rate *rate);
int ch_model_finalize(ch_model *, int reorder, int ignore_unused);
size_t ch_model_species_count(const ch_model *);
size_t ch_model_reaction_count(const ch_model *);
size_t ch_model_parameter_count(const ch_model *);
const char *ch_model_species_name(const ch_model *, size_t index);
const char *ch_model_parameter_name(const ch_model *, size_t index);
int ch_model_species_index(const ch_model *, const char *name);
int ch_model_parameter_index(const ch_model *, const char *name);
ch_rate *ch_model_rate(ch_model *, size_t reaction);

/* A state owns one independent adaptive solver group, in species-major order. */
ch_state *ch_state_create(const ch_model *, size_t cells);
ch_state *ch_state_clone(const ch_state *);
void ch_state_free(ch_state *);
size_t ch_state_cells(const ch_state *);
double *ch_state_values(ch_state *);
double *ch_state_parameters(ch_state *);
double *ch_state_rates(ch_state *);
ch_conditions *ch_state_conditions(ch_state *);
int ch_state_set(ch_state *, const char *name, const double *values, size_t count);
int ch_state_set_parameter(ch_state *, const char *name, const double *, size_t count);
int ch_state_tolerances(ch_state *, double rtol, const double *atol, size_t count);
int ch_update_rates(ch_state *);
void ch_force(const ch_state *, const double *y, double *forcing);
/* Jacobian values are stored by symbolic index, then cell. */
void ch_jacobian(const ch_state *, const double *y, double *negative_jacobian);
size_t ch_model_jacobian_size(const ch_model *);
int ch_model_jacobian_index(const ch_model *, size_t row, size_t column);

typedef struct {
  size_t internal_steps, accepted, rejected, function_calls, jacobian_updates;
  size_t decompositions, solves, constraint_init_iterations;
  double final_time;
} ch_stats;
typedef struct ch_rosenbrock_parameters ch_rosenbrock_parameters;
typedef double (*ch_step_controller)(double error, const ch_rosenbrock_parameters *, void *);
struct ch_rosenbrock_parameters {
  int stages, max_steps, new_function[6];
  double round_off, factor_min, factor_max, rejection_factor, safety_factor;
  double h_min, h_max, h_start, estimator_order;
  double a[15], c[15], m[6], e[6], alpha[6], gamma[6];
  int init_max_iterations, init_max_backtracks;
  double init_tolerance, init_backtrack_factor, init_sufficient_decrease;
  ch_step_controller controller;
  void *controller_context;
};
/* Parameter sets: 2, 3, 4, 43 (four-stage DAE), 6 (six-stage DAE). */
ch_rosenbrock_parameters ch_rosenbrock_default(int method);
double ch_default_controller(double, const ch_rosenbrock_parameters *, void *);
double ch_normalized_error(const ch_state *, const double *, const double *, const double *);
int ch_solve(ch_state *, double time_step, const ch_rosenbrock_parameters *, ch_stats *);
void ch_stats_add(ch_stats *, const ch_stats *);

/* Symbolic sparse LU can also be used independently of a chemistry model. */
typedef struct ch_sparse ch_sparse;
ch_sparse *ch_sparse_create(size_t n, const unsigned char *pattern);
void ch_sparse_free(ch_sparse *);
size_t ch_sparse_size(const ch_sparse *);
int ch_sparse_index(const ch_sparse *, size_t row, size_t col);
int ch_markowitz(size_t n, const unsigned char *pattern, int *permutation);
int ch_lu_factor(const ch_sparse *, size_t cells, double *a, int algorithm);
void ch_lu_solve(const ch_sparse *, size_t cells, const double *a, double *x);

#endif
