/* The library's pieces for the Python prototypes, through ctypes: a TS1
 * solver with the generated kernels, states, rate constants, the forcing,
 * the shifted Jacobian factored in the state, the triangular solves, the
 * dense Jacobian, and the Rosenbrock parameter sets, with concentrations in
 * row-major (cell, variable) arrays. A development tool, built by
 * `make proto` into build/libtmc_proto; not part of the library or the
 * benchmarks. */
#include "internal.h"
#include "ts1.h"

#include <math.h>

typedef struct
{
  tmc_solver *solver;
  tmc_system system;
} proto_solver;

typedef struct
{
  const proto_solver *owner;
  tmc_state *state;
  tmc_dense *y, *f; /* scratch, vector ordered */
  tmc_sparse *lu[TMC_MAX_EXTRAPOLATION]; /* stored factorizations (slots), made on demand */
  double *recip[TMC_MAX_EXTRAPOLATION];
} proto_state;

proto_solver *proto_open(int generic)
{
  proto_solver *ps = tmc_calloc(1, sizeof *ps);
  ps->system = ts1_system();
  tmc_reaction *reactions;
  int n = ts1_reactions(&reactions);
  tmc_builder *b = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(b, &ps->system);
  tmc_builder_set_reactions(b, reactions, n);
  tmc_builder_set_ignore_unused_species(b, true);
  if (!generic)
    tmc_builder_set_kernels(b, ts1_kernels());
  tmc_error err;
  ps->solver = tmc_builder_build(b, &err);
  tmc_builder_free(b);
  tmc_reactions_free(reactions, n);
  if (!ps->solver)
  {
    tmc_system_free(&ps->system);
    free(ps);
    return NULL;
  }
  return ps;
}

void proto_close(proto_solver *ps)
{
  if (!ps)
    return;
  tmc_solver_free(ps->solver);
  tmc_system_free(&ps->system);
  free(ps);
}

int proto_n_vars(const proto_solver *ps) { return ps->solver->n_vars; }
int proto_n_custom(const proto_solver *ps) { return ps->solver->n_custom; }
int proto_n_reactions(const proto_solver *ps) { return ps->solver->n_reactions; }
const char *proto_variable_name(const proto_solver *ps, int i) { return ps->solver->variable_names[i]; }
const char *proto_custom_label(const proto_solver *ps, int i) { return ps->solver->custom_labels[i]; }
int proto_uses_kernels(const proto_solver *ps) { return ps->solver->kernels != NULL; }

proto_state *proto_state_create(const proto_solver *ps, int n_cells)
{
  proto_state *s = tmc_calloc(1, sizeof *s);
  s->owner = ps;
  s->state = tmc_solver_create_state(ps->solver, n_cells);
  ts1_initial_state(s->state);
  s->y = tmc_dense_create(n_cells, ps->solver->n_vars, s->state->L, 0.0);
  s->f = tmc_dense_create(n_cells, ps->solver->n_vars, s->state->L, 0.0);
  return s;
}

void proto_state_free(proto_state *s)
{
  if (!s)
    return;
  tmc_state_free(s->state);
  tmc_dense_free(s->y);
  tmc_dense_free(s->f);
  for (int i = 0; i < TMC_MAX_EXTRAPOLATION; ++i)
  {
    tmc_sparse_free(s->lu[i]);
    tmc_free_doubles(s->recip[i]);
  }
  free(s);
}

void proto_set_conditions(proto_state *s, int cell, double temperature, double pressure)
{
  tmc_conditions *c = &s->state->conditions[cell];
  c->temperature = temperature;
  c->pressure = pressure;
  tmc_conditions_calculate_ideal_air_density(c);
}

void proto_set_variable(proto_state *s, int cell, int v, double value) { tmc_dense_set(s->state->variables, cell, v, value); }
double proto_get_variable(const proto_state *s, int cell, int v) { return tmc_dense_get(s->state->variables, cell, v); }
void proto_set_custom(proto_state *s, int cell, int k, double value) { tmc_dense_set(s->state->custom_rate_parameters, cell, k, value); }
double proto_get_custom(const proto_state *s, int cell, int k) { return tmc_dense_get(s->state->custom_rate_parameters, cell, k); }
double proto_rate_constant(const proto_state *s, int cell, int r) { return tmc_dense_get(s->state->rate_constants, cell, r); }
void proto_set_tolerances(proto_state *s, double rtol, double atol)
{
  s->state->relative_tolerance = rtol;
  for (int v = 0; v < s->state->n_vars; ++v)
    s->state->absolute_tolerance[v] = atol;
}
double proto_absolute_tolerance(const proto_state *s, int v) { return s->state->absolute_tolerance[v]; }
void proto_update_rate_constants(proto_state *s) { tmc_solver_update_state_parameters(s->owner->solver, s->state); }
void proto_update_custom_rate_constants(proto_state *s) { tmc_solver_update_custom_rate_constants(s->owner->solver, s->state); }
double proto_rate_change(const proto_state *s) { return s->state->rc_change; }
void proto_reset_rate_change(proto_state *s) { s->state->rc_change = 0.0; }

/* row-major (cell, variable) <-> the state's vector ordering */
static void load(tmc_dense *d, const double *rows)
{
  for (int c = 0; c < d->n_rows; ++c)
    for (int v = 0; v < d->n_cols; ++v)
      tmc_dense_set(d, c, v, rows[(size_t)c * d->n_cols + v]);
}

static void store(const tmc_dense *d, double *rows)
{
  for (int c = 0; c < d->n_rows; ++c)
    for (int v = 0; v < d->n_cols; ++v)
      rows[(size_t)c * d->n_cols + v] = tmc_dense_get(d, c, v);
}

void proto_get_variables(const proto_state *s, double *rows) { store(s->state->variables, rows); }
void proto_set_variables(proto_state *s, const double *rows) { load(s->state->variables, rows); }

/* f = F(y) */
void proto_forcing(proto_state *s, const double *y, double *f)
{
  load(s->y, y);
  tmc_rosenbrock_forcing(s->owner->solver, s->state, s->y, s->f);
  store(s->f, f);
}

/* alpha I - J(y), factored into the state */
void proto_factor(proto_state *s, const double *y, double alpha)
{
  tmc_stats stats;
  memset(&stats, 0, sizeof stats);
  load(s->y, y);
  tmc_rosenbrock_factor(s->owner->solver, s->state, s->y, alpha, &stats);
}

/* x = (alpha I - J)^-1 x with the state's factorization */
void proto_solve(proto_state *s, double *x)
{
  tmc_stats stats;
  memset(&stats, 0, sizeof stats);
  load(s->f, x);
  tmc_rosenbrock_linear_solve(s->owner->solver, s->state, s->f, &stats);
  store(s->f, x);
}

/* alpha I - J(y) factored into slot `slot` (0 to 7), kept until overwritten:
 * the state's own matrix is factored as proto_factor does, then copied. */
void proto_factor_slot(proto_state *s, int slot, const double *y, double alpha)
{
  const tmc_solver *sv = s->owner->solver;
  const size_t n_recip = (size_t)s->state->n_blocks * (size_t)s->state->L * (size_t)sv->n_vars;
  if (!s->lu[slot])
  {
    s->lu[slot] = tmc_lu_matrix(sv->lu, s->state->n_cells, 0.0);
    s->recip[slot] = tmc_alloc_doubles(n_recip);
  }
  proto_factor(s, y, alpha);
  memcpy(s->lu[slot]->data, s->state->jacobian->data, s->state->jacobian->size * sizeof(double));
  memcpy(s->recip[slot], s->state->reciprocals, n_recip * sizeof(double));
}

/* x = (alpha I - J)^-1 x with slot `slot`'s factorization */
void proto_solve_slot(proto_state *s, int slot, double *x)
{
  const tmc_solver *sv = s->owner->solver;
  tmc_stats stats;
  memset(&stats, 0, sizeof stats);
  load(s->f, x);
  if (sv->kernels)
    sv->kernels->lu_solve(s->f->data, s->lu[slot]->data, s->recip[slot], s->state->n_blocks);
  else
    tmc_lu_solve_in_place(sv->lu, s->lu[slot], s->recip[slot], s->f);
  store(s->f, x);
}

/* The dense Jacobian dF/dy of one cell, n_vars x n_vars row-major, from the
 * generic sparse assembly (the kernels' J agrees with it to roundoff). */
void proto_dense_jacobian(proto_state *s, int cell, const double *y, double *J)
{
  const tmc_solver *sv = s->owner->solver;
  const int n = sv->n_vars;
  load(s->y, y);
  tmc_sparse *neg = tmc_lu_matrix(sv->lu, s->state->n_cells, 0.0);
  tmc_eval_jacobian(sv, s->state, s->y, neg); /* -J in the LU pattern, zero in the fill */
  for (int i = 0; i < n * n; ++i)
    J[i] = 0.0;
  for (int e = 0; e < neg->nnz; ++e)
    J[neg->row_of[e] * n + neg->col_of[e]] = -tmc_sparse_get(neg, cell, neg->row_of[e], neg->col_of[e]);
  tmc_sparse_free(neg);
}

/* A parameter set by name: stages, a[15], c[15], m[6], e[6], gamma, the
 * new-function-evaluation flags, and the estimator order; 0 if unknown. */
int proto_tableau(const char *name, int *stages, double *a, double *c, double *m, double *e, double *gamma, int *new_f, double *order)
{
  tmc_rosenbrock_parameters p;
  if (!strcmp(name, "ros2"))
    p = tmc_rosenbrock_two_stage();
  else if (!strcmp(name, "ros3") || !strcmp(name, "micm3"))
    p = tmc_rosenbrock_three_stage();
  else if (!strcmp(name, "ros4"))
    p = tmc_rosenbrock_four_stage();
  else if (!strcmp(name, "rodas3"))
    p = tmc_rosenbrock_four_stage_dae();
  else if (!strcmp(name, "rodas4"))
    p = tmc_rosenbrock_six_stage_dae();
  else if (!strcmp(name, "ros34pw2"))
    p = tmc_rosenbrock_ros34pw2();
  else
    return 0;
  *stages = p.stages;
  for (int i = 0; i < 15; ++i)
  {
    a[i] = p.a[i];
    c[i] = p.c[i];
  }
  for (int i = 0; i < TMC_MAX_STAGES; ++i)
  {
    m[i] = p.m[i];
    e[i] = p.e[i];
    new_f[i] = p.new_function_evaluation[i];
  }
  *gamma = p.gamma[0];
  *order = p.estimator_of_local_order;
  return 1;
}

/* MICM's loop constants, so the prototype can mirror them. */
void proto_loop_constants(double *factor_min, double *factor_max, double *rejection_factor_decrease, double *safety_factor, int *max_steps, double *round_off)
{
  tmc_rosenbrock_parameters p = tmc_rosenbrock_three_stage();
  *factor_min = p.factor_min;
  *factor_max = p.factor_max;
  *rejection_factor_decrease = p.rejection_factor_decrease;
  *safety_factor = p.safety_factor;
  *max_steps = p.max_number_of_steps;
  *round_off = p.round_off;
}
