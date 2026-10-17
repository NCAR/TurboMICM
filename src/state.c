/* Solver states: a group of grid cells with its concentrations, conditions,
 * rate constants, parameters, tolerances, and the solver's work arrays. */
#include "internal.h"

tmc_state *tmc_solver_create_state(const tmc_solver *s, int n_cells)
{
  tmc_state *st = tmc_calloc(1, sizeof *st);
  st->solver = s;
  st->n_cells = n_cells;
  st->L = s->width;
  st->n_blocks = (n_cells + st->L - 1) / st->L;
  st->n_vars = s->n_vars;
  st->n_reactions = s->n_reactions;
  st->n_custom = s->n_custom;
  int L = st->L;
  st->variables = tmc_dense_create(n_cells, s->n_vars, L, 0.0);
  st->rate_constants = tmc_dense_create(n_cells, s->n_reactions, L, 0.0);
  st->custom_rate_parameters = tmc_dense_create(n_cells, s->n_custom, L, 0.0);
  st->conditions = tmc_calloc((size_t)st->n_blocks * (size_t)L, sizeof *st->conditions);
  st->jacobian = tmc_lu_matrix(s->lu, n_cells, 0.0);
  st->reciprocals = tmc_calloc((size_t)st->n_blocks * (size_t)L * (size_t)s->n_vars, sizeof(double));
  st->relative_tolerance = 1.0e-6;
  st->absolute_tolerance = tmc_malloc((size_t)s->n_vars * sizeof(double));
  memcpy(st->absolute_tolerance, s->absolute_tolerance, (size_t)s->n_vars * sizeof(double));
  if (s->method == TMC_METHOD_ROSENBROCK)
  {
    st->ynew = tmc_dense_create(n_cells, s->n_vars, L, 0.0);
    st->initial_forcing = tmc_dense_create(n_cells, s->n_vars, L, 0.0);
    st->yerror = tmc_dense_create(n_cells, s->n_vars, L, 0.0);
    for (int i = 0; i < s->ros.stages; ++i)
      st->K[i] = tmc_dense_create(n_cells, s->n_vars, L, 0.0);
  }
  else
  {
    st->yn = tmc_dense_create(n_cells, s->n_vars, L, 0.0);
    st->yn1 = tmc_dense_create(n_cells, s->n_vars, L, 0.0);
    st->forcing = tmc_dense_create(n_cells, s->n_vars, L, 0.0);
  }
  return st;
}

void tmc_state_free(tmc_state *st)
{
  if (!st)
    return;
  tmc_dense_free(st->variables);
  tmc_dense_free(st->rate_constants);
  tmc_dense_free(st->custom_rate_parameters);
  free(st->conditions);
  tmc_sparse_free(st->jacobian);
  free(st->reciprocals);
  free(st->absolute_tolerance);
  tmc_dense_free(st->ynew);
  tmc_dense_free(st->initial_forcing);
  tmc_dense_free(st->yerror);
  for (int i = 0; i < TMC_MAX_STAGES; ++i)
    tmc_dense_free(st->K[i]);
  tmc_dense_free(st->yn);
  tmc_dense_free(st->yn1);
  tmc_dense_free(st->forcing);
  free(st);
}

tmc_state *tmc_state_copy(const tmc_state *src)
{
  tmc_state *st = tmc_solver_create_state(src->solver, src->n_cells);
  tmc_dense_copy_into(st->variables, src->variables);
  tmc_dense_copy_into(st->rate_constants, src->rate_constants);
  tmc_dense_copy_into(st->custom_rate_parameters, src->custom_rate_parameters);
  memcpy(st->conditions, src->conditions, (size_t)st->n_blocks * (size_t)st->L * sizeof *st->conditions);
  memcpy(st->jacobian->data, src->jacobian->data, st->jacobian->size * sizeof(double));
  st->relative_tolerance = src->relative_tolerance;
  memcpy(st->absolute_tolerance, src->absolute_tolerance, (size_t)st->n_vars * sizeof(double));
  return st;
}

int tmc_state_num_cells(const tmc_state *st)
{
  return st->n_cells;
}
int tmc_state_num_variables(const tmc_state *st)
{
  return st->n_vars;
}
int tmc_state_num_blocks(const tmc_state *st)
{
  return st->n_blocks;
}
int tmc_state_width(const tmc_state *st)
{
  return st->L;
}
tmc_dense *tmc_state_variables(tmc_state *st)
{
  return st->variables;
}
tmc_dense *tmc_state_rate_constants(tmc_state *st)
{
  return st->rate_constants;
}
tmc_dense *tmc_state_custom_rate_parameters(tmc_state *st)
{
  return st->custom_rate_parameters;
}
tmc_sparse *tmc_state_jacobian(tmc_state *st)
{
  return st->jacobian;
}
tmc_conditions *tmc_state_conditions(tmc_state *st)
{
  return st->conditions;
}

double tmc_state_get(const tmc_state *st, int cell, int variable)
{
  return tmc_dense_get(st->variables, cell, variable);
}

void tmc_state_set(tmc_state *st, int cell, int variable, double value)
{
  tmc_dense_set(st->variables, cell, variable, value);
}

int tmc_state_variable_index(const tmc_state *st, const char *name)
{
  for (int i = 0; i < st->n_vars; ++i)
    if (strcmp(st->solver->variable_names[i], name) == 0)
      return i;
  return -1;
}

const char *tmc_state_variable_name(const tmc_state *st, int index)
{
  return st->solver->variable_names[index];
}

int tmc_state_custom_parameter_index(const tmc_state *st, const char *label)
{
  for (int i = 0; i < st->n_custom; ++i)
    if (strcmp(st->solver->custom_labels[i], label) == 0)
      return i;
  return -1;
}

int tmc_state_set_concentration(tmc_state *st, const char *name, double value, tmc_error *err)
{
  int v = tmc_state_variable_index(st, name);
  if (v < 0)
    return tmc_error_set(err, TMC_ERROR_STATE, TMC_STATE_UNKNOWN_SPECIES, "Unknown species '%s'", name);
  for (int c = 0; c < st->n_cells; ++c)
    tmc_dense_set(st->variables, c, v, value);
  return TMC_OK;
}

int tmc_state_set_concentrations(tmc_state *st, const char *name, const double *values, int n, tmc_error *err)
{
  int v = tmc_state_variable_index(st, name);
  if (v < 0)
    return tmc_error_set(err, TMC_ERROR_STATE, TMC_STATE_UNKNOWN_SPECIES, "Unknown species '%s'", name);
  if (n != st->n_cells)
    return tmc_error_set(err, TMC_ERROR_STATE, TMC_STATE_INVALID_CONCENTRATION_COUNT_MULTIGRID, "Expected %d concentrations for '%s', got %d", st->n_cells, name, n);
  for (int c = 0; c < n; ++c)
    tmc_dense_set(st->variables, c, v, values[c]);
  return TMC_OK;
}

int tmc_state_set_custom_rate_parameter(tmc_state *st, const char *label, double value, tmc_error *err)
{
  int p = tmc_state_custom_parameter_index(st, label);
  if (p < 0)
    return tmc_error_set(err, TMC_ERROR_STATE, TMC_STATE_UNKNOWN_RATE_CONSTANT_PARAMETER, "Unknown rate constant parameter '%s'", label);
  for (int c = 0; c < st->n_cells; ++c)
    tmc_dense_set(st->custom_rate_parameters, c, p, value);
  return TMC_OK;
}

int tmc_state_set_custom_rate_parameters(tmc_state *st, const char *label, const double *values, int n, tmc_error *err)
{
  int p = tmc_state_custom_parameter_index(st, label);
  if (p < 0)
    return tmc_error_set(err, TMC_ERROR_STATE, TMC_STATE_UNKNOWN_RATE_CONSTANT_PARAMETER, "Unknown rate constant parameter '%s'", label);
  if (n != st->n_cells)
    return tmc_error_set(err, TMC_ERROR_STATE, TMC_STATE_INVALID_CUSTOM_RATE_PARAM_COUNT_MULTIGRID, "Expected %d values for '%s', got %d", st->n_cells, label, n);
  for (int c = 0; c < n; ++c)
    tmc_dense_set(st->custom_rate_parameters, c, p, values[c]);
  return TMC_OK;
}

void tmc_state_set_relative_tolerance(tmc_state *st, double rtol)
{
  st->relative_tolerance = rtol;
}

double tmc_state_relative_tolerance(const tmc_state *st)
{
  return st->relative_tolerance;
}

int tmc_state_set_absolute_tolerances(tmc_state *st, const double *atol, int n, tmc_error *err)
{
  if (n != st->n_vars)
    return tmc_error_set(err, TMC_ERROR_SOLVER, TMC_SOLVER_INVALID_TOLERANCE_SIZE, "Expected %d absolute tolerances, got %d", st->n_vars, n);
  memcpy(st->absolute_tolerance, atol, (size_t)n * sizeof(double));
  return TMC_OK;
}

const double *tmc_state_absolute_tolerances(const tmc_state *st)
{
  return st->absolute_tolerance;
}

void tmc_state_clamp(tmc_state *st)
{
  tmc_dense_max(st->variables, 0.0);
}
