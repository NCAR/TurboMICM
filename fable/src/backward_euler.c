/* Backward Euler integrator (MICM's AbstractBackwardEuler::Solve): Newton
 * iterations on each step, halving the step on non-convergence and doubling
 * it after two successes. */
#include "internal.h"

#include <math.h>

/* MICM's convergence test: every residual is finite and either below `small`,
 * below the variable's absolute tolerance, or below rtol |y|. Padded lanes
 * are excluded. */
static bool is_converged(const tmc_state *st, const tmc_backward_euler_parameters *p, const tmc_dense *residual, const tmc_dense *yn1)
{
  const int L = st->L, n_vars = st->n_vars, n_cells = st->n_cells, n_blocks = st->n_blocks;
  const double rtol = st->relative_tolerance;
  for (int v = 0; v < n_vars; ++v)
  {
    const double atol = st->absolute_tolerance[v];
    for (int block = 0; block < n_blocks; ++block)
    {
      const size_t base = ((size_t)block * n_vars + v) * L;
      const int lanes = block == n_blocks - 1 ? n_cells - block * L : L;
      for (int lane = 0; lane < lanes; ++lane)
      {
        double r = residual->data[base + lane];
        if (!isfinite(r))
          return false;
        double ar = fabs(r);
        if (ar > p->small && ar > atol && ar > rtol * fabs(yn1->data[base + lane]))
          return false;
      }
    }
  }
  return true;
}

static bool any_nonfinite(const tmc_dense *m, bool *nan_found, bool *inf_found)
{
  *nan_found = *inf_found = false;
  for (size_t i = 0; i < m->size; ++i)
  {
    if (isnan(m->data[i]))
      *nan_found = true;
    else if (isinf(m->data[i]))
      *inf_found = true;
  }
  return *nan_found || *inf_found;
}

tmc_result tmc_backward_euler_solve(tmc_solver *s, double time_step, tmc_state *st)
{
  const tmc_backward_euler_parameters *p = &s->be;
  tmc_result result;
  memset(&result, 0, sizeof result);
  result.state = TMC_NOT_YET_CALLED;
  tmc_dense *Yn = st->yn;
  tmc_dense *Yn1 = st->variables;
  tmc_dense *forcing = st->forcing;
  const int max_iter = p->max_number_of_steps;
  const int n_reductions = p->time_step_reductions_count;

  double H = p->h_start == 0.0 ? time_step : p->h_start;
  double present_time = 0.0;
  int n_successful_integrations = 0;
  int n_convergence_failures = 0;

  tmc_dense_copy_into(Yn, Yn1);
  while (present_time < time_step)
  {
    result.state = TMC_RUNNING;
    bool converged = false;
    int iterations = 0;
    do
    {
      result.stats.number_of_steps++;
      tmc_eval_forcing(s, st, Yn1, forcing);
      result.stats.function_calls++;
      tmc_eval_jacobian(s, st, Yn1, st->jacobian);
      result.stats.jacobian_updates++;
      tmc_sparse_add_to_diagonal(st->jacobian, 1.0 / H);
      tmc_lu_decompose_in_place(s->lu, st->jacobian, st->reciprocals);
      result.stats.decompositions++;
      {
        double *restrict f = forcing->data;
        const double *restrict y1 = Yn1->data, *restrict y0 = Yn->data;
        for (size_t i = 0; i < forcing->size; ++i)
          f[i] -= (y1[i] - y0[i]) / H;
      }
      tmc_lu_solve_in_place(s->lu, st->jacobian, NULL, forcing);
      result.stats.solves++;
      {
        double *restrict y1 = Yn1->data;
        const double *restrict f = forcing->data;
        for (size_t i = 0; i < Yn1->size; ++i)
        {
          double v = y1[i] + f[i];
          y1[i] = 0.0 > v ? 0.0 : v;
        }
      }
      if (iterations++ == 0)
        continue;
      converged = is_converged(st, p, forcing, Yn1);
    } while (!converged && iterations < max_iter);

    if (!converged)
    {
      result.stats.rejected++;
      n_successful_integrations = 0;
      if (n_convergence_failures >= n_reductions)
      {
        present_time += H;
        result.state = TMC_ACCEPTING_UNCONVERGED_INTEGRATION;
        bool nan_found, inf_found;
        any_nonfinite(Yn1, &nan_found, &inf_found);
        if (nan_found)
        {
          result.state = TMC_NAN_DETECTED;
          break;
        }
        if (inf_found)
          result.state = TMC_INF_DETECTED;
        break;
      }
      tmc_dense_copy_into(Yn1, Yn);
      H *= p->time_step_reductions[n_convergence_failures++];
    }
    else
    {
      result.stats.accepted++;
      result.state = TMC_CONVERGED;
      present_time += H;
      tmc_dense_copy_into(Yn, Yn1);
      n_successful_integrations++;
      if (n_successful_integrations >= 2)
      {
        n_successful_integrations = 0;
        H *= 2.0;
      }
    }
    H = fmin(H, time_step - present_time);
  }
  result.stats.final_time = present_time;
  return result;
}
