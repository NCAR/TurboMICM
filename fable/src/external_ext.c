/* Additions to the external-model API (include/tmc_external_ext.h): name
 * maps built by hand, and the constraint half of the forcing and Jacobian
 * (MICM's ConstraintBundle). */
#include "internal.h"

#include "tmc_external_ext.h"

tmc_name_map *tmc_name_map_create(const char *const *names, int n)
{
  tmc_name_map *map = tmc_calloc(1, sizeof *map);
  map->n = n > 0 ? n : 0;
  map->names = tmc_calloc((size_t)(map->n > 0 ? map->n : 1), sizeof *map->names);
  for (int i = 0; i < map->n; ++i)
    tmc_strlcpy(map->names[i], names[i], TMC_NAME_MAX);
  return map;
}

void tmc_name_map_free(tmc_name_map *map)
{
  if (!map)
    return;
  free(map->names);
  free(map);
}

/* The same sequence as the tail of tmc_eval_forcing (process.c). */
void tmc_solver_add_constraint_terms(const tmc_solver *s, const tmc_state *state, const tmc_dense *variables, tmc_dense *forcing)
{
  if (s->n_constraints)
    tmc_constraints_add_forcing(s, variables, state->custom_rate_parameters, forcing);
  if (s->n_external)
    tmc_external_add_constraint_residual(s, state->custom_rate_parameters, variables, forcing);
}

/* The same sequence as the tail of tmc_eval_jacobian (process.c). */
void tmc_solver_subtract_constraint_jacobian_terms(const tmc_solver *s, const tmc_state *state, const tmc_dense *variables, tmc_sparse *jacobian)
{
  if (s->n_constraints)
    tmc_constraints_subtract_jacobian(s, variables, state->custom_rate_parameters, jacobian);
  if (s->n_external)
    tmc_external_subtract_constraint_jacobian(s, state->custom_rate_parameters, variables, jacobian);
}
