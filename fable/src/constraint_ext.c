/* Constraint extensions (include/tmc_constraint_ext.h): the constraints'
 * Jacobian elements, and the constraint initialization run on its own. */
#include "tmc_constraint_ext.h"

#include "internal.h"

#include <float.h>

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
  tmc_result r;
  if (s->n_algebraic == 0)
  {
    /* MICM's loop finds a zero residual and returns after one iteration */
    memset(&r, 0, sizeof r);
    r.state = TMC_CONVERGED;
    r.stats.constraint_init_iterations = 1;
  }
  else
  {
    /* The Rosenbrock solve projects first and then steps while
     * present_time - time_step + round_off <= 0, which a zero time step and a
     * positive round-off make false at once. The projection reads neither.
     * A shallow copy of the solver carries the parameters, so the caller's
     * solver is not changed. */
    tmc_solver projection = *s;
    if (parameters)
      projection.ros = *parameters;
    if (!(projection.ros.round_off > 0.0))
      projection.ros.round_off = DBL_EPSILON;
    r = tmc_rosenbrock_solve(&projection, 0.0, state);
  }
  if (stats)
  {
    stats->function_calls += r.stats.function_calls;
    stats->jacobian_updates += r.stats.jacobian_updates;
    stats->number_of_steps += r.stats.number_of_steps;
    stats->accepted += r.stats.accepted;
    stats->rejected += r.stats.rejected;
    stats->decompositions += r.stats.decompositions;
    stats->solves += r.stats.solves;
    stats->constraint_init_iterations += r.stats.constraint_init_iterations;
  }
  return r.state;
}
