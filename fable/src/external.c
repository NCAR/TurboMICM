/* External models: dispatch of their callbacks (MICM's RatesBundle and
 * ConstraintBundle) and the name maps they receive. */
#include "internal.h"

int tmc_name_map_lookup(const tmc_name_map *map, const char *name)
{
  for (int i = 0; i < map->n; ++i)
    if (strcmp(map->names[i], name) == 0)
      return i;
  return -1;
}

int tmc_name_map_size(const tmc_name_map *map)
{
  return map->n;
}

const char *tmc_name_map_name(const tmc_name_map *map, int index)
{
  return map->names[index];
}

static int compare_names(const void *a, const void *b)
{
  return strcmp(a, b);
}

int tmc_external_names(const tmc_external_model *m, int (*cb)(void *, char (*)[TMC_NAME_MAX], int), char (**out)[TMC_NAME_MAX])
{
  *out = NULL;
  if (!cb)
    return 0;
  int n = cb(m->context, NULL, 0);
  if (n <= 0)
    return 0;
  char (*names)[TMC_NAME_MAX] = tmc_calloc((size_t)n, sizeof *names);
  int got = cb(m->context, names, n);
  if (got < n)
    n = got;
  /* std::set order and uniqueness */
  qsort(names, (size_t)n, sizeof *names, compare_names);
  int k = 0;
  for (int i = 0; i < n; ++i)
    if (k == 0 || strcmp(names[k - 1], names[i]) != 0)
      memcpy(names[k++], names[i], TMC_NAME_MAX);
  *out = names;
  return k;
}

static bool has_processes(const tmc_external_model *m)
{
  return m->add_forcing_terms != NULL;
}

static bool has_constraints(const tmc_external_model *m)
{
  return m->constraint_algebraic_variable_names != NULL;
}

void tmc_external_update_parameters(const tmc_solver *s, const tmc_conditions *conditions, int n_cells, tmc_dense *params)
{
  for (int i = 0; i < s->n_external; ++i)
  {
    const tmc_external_model *m = &s->externals[i];
    if (has_processes(m) && m->update_state_parameters)
      m->update_state_parameters(m->context, conditions, n_cells, params);
  }
}

void tmc_external_add_forcing(const tmc_solver *s, const tmc_dense *params, const tmc_dense *y, tmc_dense *forcing)
{
  for (int i = 0; i < s->n_external; ++i)
  {
    const tmc_external_model *m = &s->externals[i];
    if (has_processes(m))
      m->add_forcing_terms(m->context, params, y, forcing);
  }
}

void tmc_external_subtract_jacobian(const tmc_solver *s, const tmc_dense *params, const tmc_dense *y, tmc_sparse *jacobian)
{
  for (int i = 0; i < s->n_external; ++i)
  {
    const tmc_external_model *m = &s->externals[i];
    if (has_processes(m) && m->subtract_jacobian_terms)
      m->subtract_jacobian_terms(m->context, params, y, jacobian);
  }
}

void tmc_external_update_constraint_parameters(const tmc_solver *s, const tmc_conditions *conditions, int n_cells, tmc_dense *params)
{
  for (int i = 0; i < s->n_external; ++i)
  {
    const tmc_external_model *m = &s->externals[i];
    if (has_constraints(m) && s->external_constraint_active[i] && m->update_constraint_state_parameters)
      m->update_constraint_state_parameters(m->context, conditions, n_cells, params);
  }
}

void tmc_external_add_constraint_residual(const tmc_solver *s, const tmc_dense *params, const tmc_dense *y, tmc_dense *forcing)
{
  for (int i = 0; i < s->n_external; ++i)
  {
    const tmc_external_model *m = &s->externals[i];
    if (has_constraints(m) && s->external_constraint_active[i] && m->add_constraint_residual)
      m->add_constraint_residual(m->context, params, y, forcing);
  }
}

void tmc_external_subtract_constraint_jacobian(const tmc_solver *s, const tmc_dense *params, const tmc_dense *y, tmc_sparse *jacobian)
{
  for (int i = 0; i < s->n_external; ++i)
  {
    const tmc_external_model *m = &s->externals[i];
    if (has_constraints(m) && s->external_constraint_active[i] && m->subtract_constraint_jacobian)
      m->subtract_constraint_jacobian(m->context, params, y, jacobian);
  }
}

void tmc_external_initialize_constraint_parameters(const tmc_solver *s, const tmc_dense *y, tmc_dense *params)
{
  for (int i = 0; i < s->n_external; ++i)
  {
    const tmc_external_model *m = &s->externals[i];
    if (has_constraints(m) && s->external_constraint_active[i] && m->initialize_constraint_parameters)
      m->initialize_constraint_parameters(m->context, y, params);
  }
}
