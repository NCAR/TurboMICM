/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>
#include <stdio.h>

#define CH_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CH_AVOGADRO 6.02214076e23
#define CH_VERSION_MAJOR 1
#define CH_VERSION_MINOR 0
#define CH_VERSION_PATCH 0
#define CH_VERSION_TWEAK 0
const char *ch_version(void);
int ch_version_component(int index);

typedef enum {
  CH_OK = 0, CH_INVALID = -1, CH_MEMORY = -2, CH_NOT_FOUND = -3,
  CH_SIZE = -4, CH_FROZEN = -5, CH_SINGULAR = -6,
  CH_MAX_STEPS = -7, CH_SMALL_STEP = -8, CH_NAN = -9, CH_INF = -10,
  CH_CONSTRAINT_INIT = -11, CH_UNCONVERGED = -12
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
typedef struct ch_species_spec ch_species_spec;
typedef struct ch_phase ch_phase;
enum { CH_PROPERTY_REAL, CH_PROPERTY_STRING, CH_PROPERTY_INT, CH_PROPERTY_BOOL };
/* Metadata is owned and copied. Strings returned by getters are borrowed. */
ch_species_spec *ch_species_create(const char *name);
ch_species_spec *ch_species_clone(const ch_species_spec *);
void ch_species_free(ch_species_spec *);
const char *ch_species_name(const ch_species_spec *);
size_t ch_species_property_count(const ch_species_spec *, int type);
int ch_species_set_property(ch_species_spec *, const char *key, int type, const void *value);
int ch_species_get_property(const ch_species_spec *, const char *key, int type, void *value);
int ch_species_parameterize(ch_species_spec *, int enabled, double constant, double temperature,
                            double pressure, double air_density);
int ch_species_is_parameterized(const ch_species_spec *);
double ch_species_parameter_value(const ch_species_spec *, const ch_conditions *);
ch_phase *ch_phase_create(const char *name);
void ch_phase_free(ch_phase *);
const char *ch_phase_name(const ch_phase *);
int ch_phase_add(ch_phase *, const ch_species_spec *);
size_t ch_phase_species_count(const ch_phase *);
size_t ch_phase_state_size(const ch_phase *);
const ch_species_spec *ch_phase_species(const ch_phase *, size_t);
const char *ch_phase_state_name(const ch_phase *, size_t index, int qualified);

ch_model *ch_model_create(void);
ch_model *ch_model_clone_builder(const ch_model *);
void ch_model_free(ch_model *);
/* Species IDs are stable builder IDs; state indices are assigned at finalize. */
int ch_model_add_species(ch_model *, const char *name, double atol, int third_body);
int ch_model_add_species_spec(ch_model *, const ch_species_spec *);
int ch_model_add_phase(ch_model *, const ch_phase *, int gas_phase);
const ch_species_spec *ch_model_species_info(const ch_model *, size_t builder_id);
size_t ch_model_builder_species_count(const ch_model *);
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
typedef struct {
  size_t reactants, products, parameterized_reactants;
  const int *reactant_ids, *product_ids, *parameterized_ids;
  const double *yields;
  const ch_rate *rate;
} ch_reaction_info;
int ch_model_reaction_info(const ch_model *, size_t reaction, ch_reaction_info *);
int ch_model_rate_group(const ch_model *, int type, size_t *offset, size_t *count);
const ch_rate *ch_model_rate_in_group(const ch_model *, int type, size_t index);
int ch_model_species_used(const ch_model *, size_t builder_id);
int ch_model_kinetic_pattern(const ch_model *, unsigned char *pattern, size_t size);

int ch_model_add_linear(ch_model *, const char *name, int algebraic_species,
                        size_t n, const int *species, const double *coefficients, double constant);
int ch_model_add_equilibrium(ch_model *, const char *name, int algebraic_species,
                             size_t nr, const int *reactants, const double *reactant_stoichiometry,
                             size_t np, const int *products, const double *product_stoichiometry,
                             double equilibrium_constant, double enthalpy);
size_t ch_model_constraint_count(const ch_model *);
int ch_model_is_algebraic(const ch_model *, size_t variable);
void ch_constraint_force(const ch_state *, const double *y, double *forcing);
void ch_constraint_jacobian(const ch_state *, const double *y, double *negative_jacobian);

/* External callbacks add to the supplied arrays; no global mutable data needed.
 * Pattern and algebraic IDs refer to builder species IDs and are copied.
 * bind runs once after symbolic setup. Context remains owned by the caller. */
typedef void (*ch_external_eval)(const ch_state *, const double *, double *, void *);
typedef void (*ch_external_update)(ch_state *, void *);
typedef struct {
  size_t elements, algebraic_count;
  const int *rows, *columns, *algebraic_species;
  void (*bind)(const ch_model *, void *);
  ch_external_update update, initialize_constraints;
  ch_external_eval force, jacobian, constraint_force, constraint_jacobian;
  void *context;
} ch_external;
int ch_model_add_external(ch_model *, const ch_external *);
const ch_model *ch_state_model(const ch_state *);

/* A state owns one independent adaptive solver group, in species-major order. */
ch_state *ch_state_create(const ch_model *, size_t cells);
ch_state *ch_state_clone(const ch_state *);
int ch_state_copy(ch_state *destination, const ch_state *source);
void ch_state_move(ch_state **destination, ch_state **source);
void ch_state_free(ch_state *);
size_t ch_state_cells(const ch_state *);
double *ch_state_values(ch_state *);
double *ch_state_parameters(ch_state *);
double *ch_state_rates(ch_state *);
ch_conditions *ch_state_conditions(ch_state *);
int ch_state_set(ch_state *, const char *name, const double *values, size_t count);
int ch_state_set_parameter(ch_state *, const char *name, const double *, size_t count);
int ch_state_set_exact(ch_state *, const char *name, const double *, size_t count);
int ch_state_set_parameter_exact(ch_state *, const char *name, const double *, size_t count);
int ch_state_set_parameter_matrix(ch_state *, const double *row_major_values, size_t rows, size_t columns);
double ch_state_relative_tolerance(const ch_state *);
const double *ch_state_absolute_tolerances(const ch_state *);
int ch_state_tolerances(ch_state *, double rtol, const double *atol, size_t count);
int ch_update_rates(ch_state *);
void ch_force(const ch_state *, const double *y, double *forcing);
int ch_add_force(const ch_state *, const double *y, double *forcing);
/* Jacobian values are stored by symbolic index, then cell. */
void ch_jacobian(const ch_state *, const double *y, double *negative_jacobian);
int ch_add_negative_jacobian(const ch_state *, const double *y, double *negative_jacobian);
void ch_alpha_minus_jacobian(const ch_state *, double alpha, double *negative_jacobian);
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
/* Independent states may share an immutable model. Each state must occur once.
 * threads <= 0 selects OpenMP's default. Status and statistics have count entries. */
int ch_solve_many(ch_state *const *states, size_t count, double time_step,
                  const ch_rosenbrock_parameters *, int threads, int *status, ch_stats *statistics);
int ch_initialize_constraints(ch_state *, const ch_rosenbrock_parameters *, ch_stats *);
void ch_stats_add(ch_stats *, const ch_stats *);

/* Experimental ODE integrators. Legacy ch_solve and its ABI are unchanged.
 * A workspace borrows one state, which must outlive it. It keeps a step-size
 * prediction across calls, but invalidates matrices at every call (rates may
 * have changed). Algebraic constraints are handled by the legacy solvers only. */
typedef struct ch_integrator ch_integrator;
typedef struct {
  int method;                    /* 2, 3, 4, 43, 6, 431: Rosenbrock; 102: SDIRK2 */
  int controller;                /* 0: integral, 1: PI, 2: H211b */
  int norm;                      /* 0: RMS, 1: max cell RMS, 2: infinity */
  int jacobian_steps, warm_start, max_steps, newton_iterations;
  double h_start, h_max, safety, factor_min, factor_max;
  double pi_beta, lu_band, positivity, newton_tolerance, filter_b, filter_k;
} ch_integrator_options;
typedef struct {
  size_t jacobian_reuses, factor_reuses, positivity_rejections;
  size_t nonfinite_rejections, clipped_values, nonlinear_iterations, nonlinear_failures;
  double clipped_amount, maximum_clip;
} ch_integrator_stats;
ch_integrator_options ch_integrator_default(int method);
ch_integrator *ch_integrator_create(ch_state *, const ch_integrator_options *);
void ch_integrator_free(ch_integrator *);
void ch_integrator_reset(ch_integrator *);
int ch_integrator_step(ch_integrator *, double dt, ch_stats *, ch_integrator_stats *);
void ch_integrator_stats_add(ch_integrator_stats *, const ch_integrator_stats *);

typedef struct {
  double small, h_start;
  int max_iterations, reduction_count;
  double reductions[16];
} ch_backward_euler_parameters;
ch_backward_euler_parameters ch_backward_euler_default(void);
int ch_backward_euler_converged(const ch_state *, const double *residual, const double *y,
                                const ch_backward_euler_parameters *);
int ch_solve_backward_euler(ch_state *, double time_step, const ch_backward_euler_parameters *, ch_stats *);

/* Symbolic sparse LU can also be used independently of a chemistry model. */
typedef struct ch_sparse ch_sparse;
typedef int (*ch_factor_function)(const ch_sparse *, size_t, double *, double *scratch, void *);
typedef void (*ch_linear_function)(const ch_sparse *, size_t, const double *, double *, void *);
int ch_state_set_lu(ch_state *, int algorithm);
int ch_state_set_linear_solver(ch_state *, ch_factor_function, ch_linear_function, void *);
ch_sparse *ch_sparse_create(size_t n, const unsigned char *pattern);
void ch_sparse_free(ch_sparse *);
size_t ch_sparse_size(const ch_sparse *);
size_t ch_sparse_dimension(const ch_sparse *);
/* Borrowed CSR metadata, including symbolic fill; valid until sparse_free. */
const int *ch_sparse_row_offsets(const ch_sparse *);
const int *ch_sparse_column_indices(const ch_sparse *);
const int *ch_sparse_diagonal_indices(const ch_sparse *);
int ch_sparse_index(const ch_sparse *, size_t row, size_t col);
int ch_markowitz(size_t n, const unsigned char *pattern, int *permutation);
int ch_lu_factor(const ch_sparse *, size_t cells, double *a, int algorithm);
int ch_lu_factor_copy(const ch_sparse *, size_t cells, const double *a, double *lu, int algorithm);
void ch_lu_solve(const ch_sparse *, size_t cells, const double *a, double *x);
/* Matrix and vectors use index-major, then cell layout; output must not alias input. */
void ch_sparse_multiply(const ch_sparse *, size_t cells, const double *a, const double *x, double *output);

typedef void (*ch_forcing_function)(size_t species, size_t cells, const double *y, double *forcing, void *context);
typedef struct {
  int passed;
  double max_absolute_error, max_relative_error;
  size_t worst_cell, worst_row, worst_column;
  double worst_analytical, worst_finite_difference;
} ch_jacobian_comparison;
/* Dense derivative layout: (row * species + column) * cells + cell. */
int ch_finite_difference(size_t species, size_t cells, const double *y, double perturbation,
                         ch_forcing_function, void *context, double *jacobian);
ch_jacobian_comparison ch_compare_jacobian(size_t species, size_t cells, const unsigned char *pattern,
                                          const double *negative_analytical, const double *finite_difference,
                                          double atol, double rtol);
ch_jacobian_comparison ch_check_sparsity(size_t species, size_t cells, const unsigned char *pattern,
                                        const double *finite_difference, double threshold);

#endif
