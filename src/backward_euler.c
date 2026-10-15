/* Newton iteration, positivity, failure reductions and growth follow MICM. */
#include "internal.h"
ck_backward_euler_parameters ck_backward_euler_defaults(void)
{
  ck_backward_euler_parameters p = {.small=1e-40, .max_iterations=11, .nreductions=5,
    .reductions={0.5,0.5,0.5,0.5,0.1}, .lu=CK_LU_MOZART};
  return p;
}
int ck_backward_euler_converged(const ck_state *s, const double *residual, const double *y, double small)
{
  for (size_t j = 0; j < s->model->n; ++j) for (size_t c = 0; c < s->cells; ++c) {
    size_t i = j * s->cells + c;
    if (!isfinite(residual[i])) return 0;
    if (fabs(residual[i]) > small && fabs(residual[i]) > s->absolute_tolerance[j] &&
        fabs(residual[i]) > s->relative_tolerance * fabs(y[i])) return 0;
  }
  return 1;
}
ck_status ck_backward_euler_solve(ck_state *s, double dt, const ck_backward_euler_parameters *p, ck_stats *stats)
{
  if (!s || !p || !stats || !(dt >= 0) || !isfinite(dt) || !p->max_iterations || p->nreductions > 16) return CK_INVALID;
  memset(stats, 0, sizeof *stats);
  const ck_model *m = s->model;
  ck_workspace *w = s->workspace;
  size_t size = m->n * s->cells;
  double h = p->h_start == 0 ? dt : p->h_start, time = 0;
  size_t successes = 0, failures = 0;
  memcpy(w->next, s->y, size * sizeof(double));
  ck_status status = CK_OK;
  while (time < dt) {
    int converged = 0;
    size_t iteration = 0;
    do {
      ++stats->internal_steps;
      ck_forcing(s, s->y, w->initial); ++stats->function_calls;
      ck_jacobian(s, s->y, w->jacobian); ++stats->jacobian_updates;
      for (size_t j = 0; j < m->n; ++j)
        for (size_t c = 0; c < s->cells; ++c) w->jacobian[m->diagonal[j] * s->cells + c] += 1 / h;
      ck_factor(m, s->cells, w->jacobian, p->lu); ++stats->decompositions;
      #pragma omp simd
      for (size_t i = 0; i < size; ++i) w->initial[i] -= (s->y[i] - w->next[i]) / h;
      ck_linear_solve(m, s->cells, w->jacobian, w->initial); ++stats->solves;
      #pragma omp simd
      for (size_t i = 0; i < size; ++i) {
        double next = s->y[i] + w->initial[i]; s->y[i] = 0 > next ? 0 : next;
      }
      if (iteration++ == 0) continue;
      converged = ck_backward_euler_converged(s, w->initial, s->y, p->small);
    } while (!converged && iteration < p->max_iterations);
    if (!converged) {
      ++stats->rejected; successes = 0;
      if (failures >= p->nreductions) {
        time += h; status = CK_UNCONVERGED;
        int inf = 0;
        for (size_t i = 0; i < size; ++i) {
          if (isnan(s->y[i])) { status = CK_NAN; break; }
          if (isinf(s->y[i])) inf = 1;
        }
        if (status != CK_NAN && inf) status = CK_INF;
        break;
      }
      memcpy(s->y, w->next, size * sizeof(double)); h *= p->reductions[failures++];
    } else {
      ++stats->accepted; time += h;
      memcpy(w->next, s->y, size * sizeof(double));
      if (++successes >= 2) { successes = 0; h *= 2; }
    }
    h = fmin(h, dt - time);
  }
  stats->final_time = time; return status;
}
