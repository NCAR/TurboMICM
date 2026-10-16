/* Backward Euler iteration and step reductions ported from MICM.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

ch_backward_euler_parameters ch_backward_euler_default(void) {
  return (ch_backward_euler_parameters){1e-40,0.0,11,5,{0.5,0.5,0.5,0.5,0.1}};
}
int ch_backward_euler_converged(const ch_state *s, const double *residual, const double *y,
                                const ch_backward_euler_parameters *p) {
  for (size_t i = 0; i < s->model->nv; ++i) for (size_t c = 0; c < s->cells; ++c) {
    size_t k = i*s->cells+c;
    double r = fabs(residual[k]);
    if (!isfinite(r) || (r > p->small && r > s->atol[i] && r > s->rtol*fabs(y[k]))) return 0;
  }
  return 1;
}
int ch_solve_backward_euler(ch_state *s, double dt, const ch_backward_euler_parameters *p, ch_stats *stats) {
  if (!s || !p || !stats || !isfinite(dt) || dt < 0 || p->max_iterations < 1 ||
      p->reduction_count < 0 || p->reduction_count > 16) return CH_INVALID;
  if (s->model->nalgebraic) return CH_INVALID; /* MICM's BE is an ODE method. */
  *stats = (ch_stats){0};
  const ch_model *m = s->model;
  size_t n = m->nv*s->cells;
  double h = p->h_start == 0.0 ? dt : p->h_start, time = 0.0;
  int successes = 0, failures = 0, status = CH_OK;
  double *old = s->new_y, *f = s->forcing;
  memcpy(old,s->y,n*sizeof(double));
  while (time < dt) {
    int converged = 0;
    for (int iteration = 0; iteration < p->max_iterations; ++iteration) {
      ++stats->internal_steps;
      ch_force(s,s->y,f); ++stats->function_calls;
      ch_jacobian(s,s->y,s->jacobian); ++stats->jacobian_updates;
      for (size_t i = 0; i < m->nv; ++i)
        for (size_t c = 0; c < s->cells; ++c) s->jacobian[(size_t)m->sparse->diagonal[i]*s->cells+c] += 1.0/h;
      status = ch_lu_factor(m->sparse,s->cells,s->jacobian,0);
      if (status) goto done;
      ++stats->decompositions;
      #pragma omp simd
      for (size_t i = 0; i < n; ++i) f[i] -= (s->y[i]-old[i])/h;
      ch_lu_solve(m->sparse,s->cells,s->jacobian,f); ++stats->solves;
      #pragma omp simd
      for (size_t i = 0; i < n; ++i) { double next = s->y[i]+f[i]; s->y[i] = 0.0 > next ? 0.0 : next; }
      if (!iteration) continue;
      converged = ch_backward_euler_converged(s,f,s->y,p);
      if (converged) break;
    }
    if (!converged) {
      ++stats->rejected; successes = 0;
      if (failures >= p->reduction_count) {
        time += h; status = CH_MAX_STEPS;
        for (size_t i = 0; i < n; ++i) {
          if (isnan(s->y[i])) { status = CH_NAN; break; }
          if (isinf(s->y[i])) status = CH_INF;
        }
        break;
      }
      memcpy(s->y,old,n*sizeof(double)); h *= p->reductions[failures++];
    } else {
      ++stats->accepted; time += h; memcpy(old,s->y,n*sizeof(double));
      if (++successes >= 2) { successes = 0; h *= 2.0; }
    }
    h = fmin(h,dt-time);
    if (time < dt && (h <= 0 || time+h == time)) { status = CH_SMALL_STEP; break; }
  }
done:
  stats->final_time = time; ch_post_solve(s); return status;
}
