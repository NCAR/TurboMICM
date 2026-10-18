/* Constraint extensions (include/tmc_constraint_ext.h): the constraints'
 * Jacobian elements, and the constraint initialization run on its own. */
#include "tmc_constraint_ext.h"

#include "internal.h"

#include <string.h>

int tmc_solver_constraint_jacobian_elements(const tmc_solver *s, tmc_index_pair *out, int max)
{
  const int n = s->n_vars;
  unsigned char *present = tmc_calloc((size_t)n * (size_t)n, 1);
  for (int k = 0; k < s->n_constraints; ++k)
  {
    const tmc_constraint_plan *plan = &s->constraint_plans[k];
    present[(size_t)plan->row * n + plan->row] = 1;
    for (int i = 0; i < plan->n_dependencies; ++i)
      present[(size_t)plan->row * n + plan->dependencies[i]] = 1;
  }
  int count = 0;
  for (int row = 0; row < n; ++row)
    for (int col = 0; col < n; ++col)
      if (present[(size_t)row * n + col])
      {
        if (out && count < max)
          out[count] = (tmc_index_pair){ row, col };
        ++count;
      }
  free(present);
  return count;
}

tmc_solver_state tmc_solver_initialize_constraints(const tmc_solver *s, tmc_state *state, const tmc_rosenbrock_parameters *parameters, tmc_stats *stats)
{
  if (s->method != TMC_METHOD_ROSENBROCK)
    return TMC_CONSTRAINT_INITIALIZATION_FAILED;
  tmc_stats projection;
  memset(&projection, 0, sizeof projection);
  tmc_solver_state result = tmc_rosenbrock_initialize_constraints(s, state, parameters ? parameters : &s->ros, &projection);
  if (stats)
  {
    stats->function_calls += projection.function_calls;
    stats->jacobian_updates += projection.jacobian_updates;
    stats->number_of_steps += projection.number_of_steps;
    stats->accepted += projection.accepted;
    stats->rejected += projection.rejected;
    stats->decompositions += projection.decompositions;
    stats->solves += projection.solves;
    stats->constraint_init_iterations += projection.constraint_init_iterations;
  }
  return result;
}
