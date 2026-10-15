#ifndef CHEM_H
#define CHEM_H

#include <stddef.h>
#include <stdint.h>

#define CHEM_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CHEM_PI 3.14159265358979323846

typedef enum {
  CHEM_OK, CHEM_INVALID, CHEM_MEMORY, CHEM_UNKNOWN, CHEM_SIZE,
  CHEM_MAX_STEPS, CHEM_SMALL_STEP, CHEM_NAN, CHEM_INF, CHEM_INIT_FAILED, CHEM_UNCONVERGED
} ChemStatus;
const char *chem_status_string(ChemStatus status);

typedef struct { double temperature, pressure, air_density; } ChemConditions;
void chem_ideal_air_density(ChemConditions *conditions);
typedef double (*ChemRateFunction)(const ChemConditions *, const double *, void *);

/* The enum order is the stable reaction order used by the reference solver. */
typedef enum {
  CHEM_ARRHENIUS, CHEM_TROE, CHEM_TERNARY, CHEM_JPL19, CHEM_BRANCHED,
  CHEM_TUNNELING, CHEM_TAYLOR, CHEM_REVERSIBLE, CHEM_USER, CHEM_SURFACE, CHEM_LAMBDA
} ChemRateType;
typedef struct {
  ChemRateType type;
  double A, B, C, D, E;
  double k0_A, k0_B, k0_C, k0_D;
  double kinf_A, kinf_B, kinf_C, kinf_D;
  double kint_A, kint_B, kint_C, kint_D;
  double Fc, N;
  double X, Y, a0, n, k0, z;
  int branch;
  double diffusion, molecular_weight, probability;
  double coefficients[10];
  size_t ncoefficients;
  char *label;
  size_t custom_index;
  ChemRateFunction function;
  void *context;
  double mean_free_speed_factor;
} ChemRate;
ChemRate chem_rate(ChemRateType type);
double chem_rate_value(const ChemRate *rate, const ChemConditions *conditions, const double *custom);

typedef struct {
  char *name;
  double initial, absolute_tolerance, molecular_weight, diffusion;
  int parameterized; /* 0: variable, 1: ideal air density, 2: callback */
  ChemRateFunction parameterization;
  void *context;
  size_t variable;
  struct ChemProperty *properties;
  size_t nproperties;
} ChemSpecies;

typedef struct {
  size_t nreactants, nproducts;
  size_t *reactants, *products;
  double *yields;
  ChemRate rate;
} ChemReaction;

typedef enum { CHEM_LINEAR, CHEM_EQUILIBRIUM, CHEM_CUSTOM_CONSTRAINT } ChemConstraintType;
typedef double (*ChemConstraintFunction)(const double *, const double *, void *);
typedef void (*ChemConstraintJacobian)(const double *, const double *, double *, void *);
typedef struct {
  ChemConstraintType type;
  size_t algebraic, nterms, nreactants;
  size_t *species;
  double *coefficients;
  double constant;
  ChemConstraintFunction residual;
  ChemConstraintJacobian jacobian;
  void *context;
  char *name;
  double delta_enthalpy, reference_temperature;
  size_t custom_index;
} ChemConstraint;

typedef struct { size_t first, count; } ChemRange;
typedef struct { size_t target, left, right; } ChemTriple;
typedef struct { size_t reaction, column, omitted, nreactants, nproducts; } ChemJacJob;
typedef struct ChemModel ChemModel;
typedef struct ChemState ChemState;
typedef void (*ChemKernel)(const ChemModel *, size_t, size_t, const double *, const double *, double *);
typedef void (*ChemExternalKernel)(const ChemModel *, const ChemState *, const double *, double *, void *);
struct ChemModel {
  size_t nspecies, nvariables, nreactions, ncustom, nconstraints;
  ChemSpecies *species;
  ChemReaction *reactions;
  ChemConstraint *constraints;
  char **custom_names;
  double *custom_defaults;
  double default_temperature, default_pressure;
  size_t nnz;
  size_t *variable_species, *diagonal;
  int *index;
  unsigned char *algebraic;
  ChemRange *lower, *upper, *scale, *updates;
  size_t *lower_ids, *upper_ids, *scale_ids;
  size_t *lower_columns, *upper_columns;
  ChemTriple *update_ops;
  size_t njobs;
  ChemJacJob *jobs;
  ChemKernel force_kernel, jacobian_kernel;
  ChemExternalKernel external_force, external_jacobian;
  void *external_context;
  int finalized;
  unsigned char *declared_sparsity;
  int ignore_unused;
  void (*setup)(const ChemModel *model, void *context);
  void *setup_context;
};
ChemModel *chem_model_create(void);
void chem_model_free(ChemModel *model);
ChemStatus chem_model_species(ChemModel *model, const char *name, double initial, int parameterized, size_t *id);
ChemStatus chem_model_reaction(ChemModel *model, const size_t *reactants, size_t nreactants,
                             const size_t *products, const double *yields, size_t nproducts, ChemRate rate);
ChemStatus chem_model_constraint(ChemModel *model, const ChemConstraint *constraint);
ChemStatus chem_constraint_sparsity(const ChemModel *model, unsigned char *row_major_pattern, size_t order);
ChemStatus chem_model_finalize(ChemModel *model);
ChemStatus chem_model_sparsity(ChemModel *model, const unsigned char *row_major_pattern, size_t order);
ChemStatus chem_markowitz(size_t order, const unsigned char *pattern, size_t *permutation);
ChemStatus chem_model_reorder(ChemModel *model);
ChemRate *chem_model_lambda(ChemModel *model, const char *label);
size_t chem_rate_offset(const ChemModel *model, ChemRateType type);
size_t chem_species_id(const ChemModel *model, const char *name);
size_t chem_variable_id(const ChemModel *model, const char *name);
size_t chem_custom_id(const ChemModel *model, const char *name);

struct ChemState {
  const ChemModel *model;
  size_t cells, stride;
  double *variables, *custom, *rates;
  ChemConditions *conditions;
  double *absolute_tolerance, relative_tolerance;
  double *jacobian, *work, *new_y, *forcing, *error, *stages[6];
};
ChemState *chem_state_create(const ChemModel *model, size_t cells);
ChemState *chem_state_copy(const ChemState *state);
void chem_state_free(ChemState *state);
ChemStatus chem_state_set(ChemState *state, const char *species, const double *values, size_t count);
ChemStatus chem_state_set_custom(ChemState *state, const char *label, const double *values, size_t count);
ChemStatus chem_state_custom_array(ChemState *state, const double *row_major_values, size_t cells, size_t parameters);
ChemStatus chem_state_tolerances(ChemState *state, double relative, const double *absolute, size_t count);
void chem_calculate_rates(ChemState *state);
void chem_forcing(const ChemState *state, const double *variables, double *forcing);
void chem_jacobian(const ChemState *state, const double *variables, double *negative_jacobian);
void chem_constraints(const ChemState *state, const double *variables, double *forcing, double *negative_jacobian);
void chem_generic_forcing(const ChemModel *, size_t, size_t, const double *, const double *, double *);
void chem_generic_jacobian(const ChemModel *, size_t, size_t, const double *, const double *, double *);

/* Sparse unpivoted LU. Values are structure-of-arrays, nnz * stride. */
void chem_lu_mozart(const ChemModel *model, size_t cells, size_t stride, double *matrix);
void chem_lu_doolittle(const ChemModel *model, size_t cells, size_t stride, double *matrix);
void chem_linear_solve(const ChemModel *model, size_t cells, size_t stride, const double *matrix, double *rhs);
ChemStatus chem_lu_separate(const ChemModel *, size_t cells, size_t stride, const double *matrix,
                            double *lower, double *upper, int doolittle);
void chem_linear_solve_separate(const ChemModel *, size_t cells, size_t stride, const double *lower,
                                const double *upper, double *rhs);

typedef struct {
  int stages, new_function[6];
  double a[15], c[15], m[6], e[6], alpha[6], gamma[6], order;
  size_t max_steps;
  double roundoff, factor_min, factor_max, rejection_factor, safety;
  double h_min, h_max, h_start;
  size_t init_iterations, init_backtracks;
  double init_tolerance, init_backtrack_factor, init_sufficient_decrease;
} ChemRosenbrock;
typedef struct {
  uint64_t internal_steps, accepted, rejected, function_calls, jacobian_updates, decompositions, solves;
  uint64_t constraint_init_iterations;
  double final_time;
} ChemStats;
typedef double (*ChemController)(double error, const ChemRosenbrock *parameters, void *context);
typedef void (*ChemFactor)(const ChemModel *, size_t, size_t, double *);
typedef void (*ChemSolve)(const ChemModel *, size_t, size_t, const double *, double *);
typedef struct {
  ChemRosenbrock parameters;
  ChemController controller;
  void *controller_context;
  ChemFactor factor;
  ChemSolve solve;
} ChemSolver;
/* Parameter sets: 2, 3, 4, 5 (four-stage DAE), 6 (six-stage DAE). */
ChemRosenbrock chem_rosenbrock_parameters(int method);
ChemSolver chem_solver(int method);
double chem_default_controller(double error, const ChemRosenbrock *parameters, void *context);
double chem_normalized_error(const ChemState *state, const double *old_y, const double *new_y, const double *errors);
void chem_alpha_minus_jacobian(ChemState *state, double alpha);
ChemStatus chem_initialize_constraints(ChemSolver *solver, ChemState *state, ChemStats *stats);
ChemStatus chem_solve(ChemSolver *solver, double time_step, ChemState *state, ChemStats *stats);

typedef struct {
  double small, h_start;
  size_t max_steps, nreductions;
  double reductions[5];
  ChemFactor factor;
  ChemSolve solve;
} ChemEuler;
ChemEuler chem_euler_parameters(void);
int chem_euler_converged(const ChemState *state, double small, const double *residual, const double *variables);
ChemStatus chem_backward_euler(const ChemEuler *parameters, double time_step, ChemState *state, ChemStats *stats);

#endif
