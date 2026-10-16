/* Rosenbrock integration and error controller ported from MICM.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

static void axpy(size_t n, double scale, const double *restrict x, double *restrict y) {
  #pragma omp simd
  for (size_t i = 0; i < n; ++i) y[i] += scale*x[i];
}
double ch_normalized_error(const ch_state *s, const double *y, const double *new_y, const double *errors) {
  double sum = 0.0;
  for (size_t i = 0; i < s->model->nv; ++i) {
    double atol = s->atol[i];
    for (size_t c = 0; c < s->cells; ++c) {
      size_t k = i*s->cells+c;
      double ymax = fabs(y[k]) > fabs(new_y[k]) ? fabs(y[k]) : fabs(new_y[k]);
      double e = errors[k]/(atol+s->rtol*ymax);
      sum += e*e;
    }
  }
  double error = sqrt(sum/(s->model->nv*s->cells));
  return error < 1e-10 ? 1e-10 : error;
}
double ch_default_controller(double error, const ch_rosenbrock_parameters *p, void *context) {
  (void)context;
  return fmin(p->factor_max,fmax(p->factor_min,p->safety_factor/pow(error,1.0/p->estimator_order)));
}
void ch_stats_add(ch_stats *a, const ch_stats *b) {
  a->internal_steps += b->internal_steps; a->accepted += b->accepted; a->rejected += b->rejected;
  a->function_calls += b->function_calls; a->jacobian_updates += b->jacobian_updates;
  a->decompositions += b->decompositions; a->solves += b->solves;
  a->constraint_init_iterations += b->constraint_init_iterations; a->final_time += b->final_time;
}
int ch_solve(ch_state *s, double dt, const ch_rosenbrock_parameters *p, ch_stats *stats) {
  if (!s || !p || !stats || !isfinite(dt) || dt < 0 || p->stages < 1 || p->stages > 6 ||
      p->gamma[0] == 0 || p->estimator_order <= 0 || !p->controller) return CH_INVALID;
  *stats = (ch_stats){0};
  const ch_model *m = s->model;
  size_t n = m->nv*s->cells;
  const double h_min = p->h_min == 0 ? 1e-15*dt : p->h_min;
  const double h_max = p->h_max == 0 ? dt : fmin(dt,p->h_max);
  const double h_start = p->h_start == 0 ? 1e-6*dt : fmin(h_max,p->h_start);
  double h = fmin(fmax(h_min,fabs(h_start)),fabs(h_max)), time = 0.0;
  int reject_last = 0, reject_more = 0, status = CH_OK;
  if (m->nalgebraic) {
    for (size_t i = 0; i < m->nexternal; ++i) {
      const ch_external *e = &m->externals[i];
      if (e->algebraic_count && e->initialize_constraints) e->initialize_constraints(s,e->context);
    }
    status = ch_initialize_constraints(s,p,stats);
    if (status) goto done;
  }
  while (time-dt+p->round_off <= 0) {
    if (stats->internal_steps > (size_t)p->max_steps) { status = CH_MAX_STEPS; break; }
    if (time+0.1*h == time || h <= fmin(p->round_off,h_min)) { status = CH_SMALL_STEP; break; }
    h = fmin(h,fabs(dt-time));
    ch_force(s,s->y,s->forcing); ++stats->function_calls;
    if (m->nalgebraic) ch_constraint_force(s,s->y,s->forcing);
    ch_jacobian(s,s->y,s->jacobian); ++stats->jacobian_updates;
    if (m->nalgebraic) ch_constraint_jacobian(s,s->y,s->jacobian);
    int accepted = 0;
    while (!accepted) {
      double alpha = 1.0/(h*p->gamma[0]);
      for (size_t i = 0; i < m->nv; ++i) {
        if (m->algebraic[i]) continue;
        double *restrict diagonal = s->jacobian+(size_t)m->sparse->diagonal[i]*s->cells;
        #pragma omp simd
        for (size_t c = 0; c < s->cells; ++c) diagonal[c] += alpha;
      }
      status = ch_lu_factor(m->sparse,s->cells,s->jacobian,0);
      if (status) goto done;
      ++stats->decompositions;
      for (int stage = 0; stage < p->stages; ++stage) {
        double *k = s->stages+(size_t)stage*n;
        int offset = stage*(stage-1)/2;
        if (!stage) memcpy(k,s->forcing,n*sizeof(double));
        else if (p->new_function[stage]) {
          memcpy(s->new_y,s->y,n*sizeof(double));
          for (int j = 0; j < stage; ++j) axpy(n,p->a[offset+j],s->stages+(size_t)j*n,s->new_y);
          ch_force(s,s->new_y,k); ++stats->function_calls;
          if (m->nalgebraic) ch_constraint_force(s,s->new_y,k);
        }
        if (stage+1 < p->stages && !p->new_function[stage+1]) memcpy(k+n,k,n*sizeof(double));
        for (int j = 0; j < stage; ++j) {
          double scale = p->c[offset+j]/h;
          const double *previous = s->stages+(size_t)j*n;
          if (!m->nalgebraic) axpy(n,scale,previous,k);
          else for (size_t i = 0; i < m->nv; ++i) if (!m->algebraic[i])
            axpy(s->cells,scale,previous+i*s->cells,k+i*s->cells);
        }
        ch_lu_solve(m->sparse,s->cells,s->jacobian,k); ++stats->solves;
      }
      memcpy(s->new_y,s->y,n*sizeof(double)); memset(s->error,0,n*sizeof(double));
      for (int stage = 0; stage < p->stages; ++stage) axpy(n,p->m[stage],s->stages+(size_t)stage*n,s->new_y);
      for (int stage = 0; stage < p->stages; ++stage) axpy(n,p->e[stage],s->stages+(size_t)stage*n,s->error);
      if (m->nalgebraic) for (size_t i = 0; i < m->nv; ++i) if (m->algebraic[i])
        for (size_t c = 0; c < s->cells; ++c) s->error[i*s->cells+c] = s->new_y[i*s->cells+c]-s->y[i*s->cells+c];
      double error = ch_normalized_error(s,s->y,s->new_y,s->error);
      double hnew = h*p->controller(error,p,p->controller_context);
      ++stats->internal_steps;
      if (isnan(error)) { status = CH_NAN; goto done; }
      if (isinf(error)) { status = CH_INF; goto done; }
      if (error < 1 || h < h_min) {
        ++stats->accepted; time += h;
        double *swap = s->y; s->y = s->new_y; s->new_y = swap;
        hnew = fmax(h_min,fmin(hnew,h_max));
        if (reject_last) hnew = fmin(hnew,h);
        reject_last = reject_more = 0; h = hnew; accepted = 1;
      } else {
        if (reject_more) hnew = h*p->rejection_factor;
        reject_more = reject_last; reject_last = 1; h = hnew;
        if (stats->accepted >= 1) ++stats->rejected;
        ch_jacobian(s,s->y,s->jacobian); ++stats->jacobian_updates;
        if (m->nalgebraic) ch_constraint_jacobian(s,s->y,s->jacobian);
        /* A failed trial cannot make progress once its step underflows. */
        if (h == 0 || !isfinite(h)) { status = CH_SMALL_STEP; goto done; }
      }
    }
  }
done:
  stats->final_time = time;
  ch_post_solve(s);
  return status;
}
