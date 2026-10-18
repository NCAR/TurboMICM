/* Adaptive ODE integration. Equations and sources: DESIGN.md, REFERENCES.md.
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <float.h>

struct ch_integrator {
  ch_state *state;
  ch_integrator_options options;
  ch_rosenbrock_parameters table;
  double *raw_jacobian;
  double next_h;
};

ch_integrator_options ch_integrator_default(int method) {
  ch_integrator_options o = {method,0,1,1,1,100000,8,
                            0.01,0.0,0.9,0.2,6.0,0.4,1.0,0.1,0.03,4.0,0.0};
  return o;
}

/* RODAS3.1's published mathematical table (KPP documentation, 2026).
 * No implementation code was imported. Other tables are the retained Tier 1
 * coefficients, avoiding the shortened numbers in the web documentation. */
static ch_rosenbrock_parameters table(int method) {
  ch_rosenbrock_parameters p = ch_rosenbrock_default(method == 431 ? 43 : method);
  if (method == 431) {
    p.gamma[0] = 0.515; p.gamma[1] = 1.628546001287715;
    p.a[1] = p.a[3] = p.m[0] = 0.646601929740551;
    p.a[2] = p.a[4] = p.m[1] = 0.409567801987914;
    p.c[0] = 4.198495621784201;
    p.c[1] = 3.711590161613010; p.c[2] = -1.787771994729384;
    p.c[3] = 4.458898153216104; p.c[4] = -2.024095448516552;
    p.c[5] = -2.626700600119396;
  }
  return p;
}

ch_integrator *ch_integrator_create(ch_state *s, const ch_integrator_options *o) {
  if (!s || !o || s->model->nalgebraic ||
      (o->method != 102 && !table(o->method).stages) ||
      o->controller < 0 || o->controller > 2 || o->norm < 0 || o->norm > 2 ||
      o->jacobian_steps < 1 || o->warm_start < 0 || o->warm_start > 1 ||
      o->max_steps < 1 || o->newton_iterations < 1 ||
      !isfinite(o->h_start) || o->h_start <= 0 || !isfinite(o->h_max) || o->h_max < 0 ||
      !isfinite(o->safety) || o->safety <= 0 || o->safety >= 1 ||
      !isfinite(o->factor_min) || o->factor_min <= 0 || o->factor_min >= 1 ||
      !isfinite(o->factor_max) || o->factor_max < 1 ||
      !isfinite(o->pi_beta) || o->pi_beta < 0 || o->pi_beta > 1 ||
      !isfinite(o->lu_band) || o->lu_band < 1 ||
      !isfinite(o->positivity) || o->positivity < 0 ||
      !isfinite(o->newton_tolerance) || o->newton_tolerance <= 0 || o->newton_tolerance >= 1 ||
      !isfinite(o->filter_b) || o->filter_b < 1 || !isfinite(o->filter_k) || o->filter_k < 0)
    return NULL;
  ch_integrator *w = calloc(1,sizeof(*w));
  if (!w) return NULL;
  w->raw_jacobian = ch_allocate(s->model->sparse->nnz*s->cells,sizeof(double));
  if (!w->raw_jacobian) { free(w); return NULL; }
  w->state = s; w->options = *o; w->table = table(o->method);
  return w;
}
void ch_integrator_free(ch_integrator *w) {
  if (w) { free(w->raw_jacobian); free(w); }
}
void ch_integrator_reset(ch_integrator *w) { if (w) w->next_h = 0; }
void ch_integrator_stats_add(ch_integrator_stats *a, const ch_integrator_stats *b) {
  a->jacobian_reuses += b->jacobian_reuses; a->factor_reuses += b->factor_reuses;
  a->positivity_rejections += b->positivity_rejections;
  a->nonfinite_rejections += b->nonfinite_rejections; a->clipped_values += b->clipped_values;
  a->nonlinear_iterations += b->nonlinear_iterations; a->nonlinear_failures += b->nonlinear_failures;
  a->clipped_amount += b->clipped_amount;
  a->maximum_clip = fmax(a->maximum_clip,b->maximum_clip);
}

static void add(size_t n, double a, const double *restrict x, double *restrict y) {
  #pragma omp simd
  for (size_t i = 0; i < n; ++i) y[i] += a*x[i];
}

/* Nonfinite errors must not be hidden by fmax or comparisons with NaN. */
static double norm(const ch_integrator *w, const double *v, const double *candidate) {
  const ch_state *s = w->state;
  size_t nv = s->model->nv, nc = s->cells;
  double sum = 0, maximum = 0;
  int bad = 0;
  double *restrict cell_sum = s->scratch;
  memset(cell_sum,0,nc*sizeof(double));
  for (size_t i = 0; i < nv; ++i) {
    const double *restrict old = s->y+i*nc;
    const double *restrict next = candidate+i*nc;
    const double *restrict errors = v+i*nc;
    const double atol = s->atol[i], rtol = s->rtol;
    if (w->options.norm == 2) {
      #pragma omp simd reduction(max:maximum) reduction(|:bad)
      for (size_t c = 0; c < nc; ++c) {
        double e = errors[c] == 0 ? 0 : fabs(errors[c])/(atol+rtol*fmax(fabs(old[c]),fabs(next[c])));
        bad |= !isfinite(e) || !isfinite(next[c]);
        maximum = fmax(maximum,e);
      }
    } else {
      /* Each cell still accumulates species in the original order. SIMD is
       * across independent cells, so this is not a reassociated reduction. */
      #pragma omp simd reduction(|:bad)
      for (size_t c = 0; c < nc; ++c) {
        double e = errors[c] == 0 ? 0 : fabs(errors[c])/(atol+rtol*fmax(fabs(old[c]),fabs(next[c])));
        bad |= !isfinite(e) || !isfinite(next[c]);
        cell_sum[c] += e*e;
      }
    }
  }
  if (bad) return INFINITY;
  for (size_t c = 0; c < nc; ++c) {
    sum += cell_sum[c];
    if (w->options.norm == 1) maximum = fmax(maximum,sqrt(cell_sum[c]/nv));
  }
  return w->options.norm == 0 ? sqrt(sum/(nv*nc)) : maximum;
}

static void rosenbrock_trial(ch_integrator *w, double h, ch_stats *stats, int evaluate_initial_force) {
  ch_state *s = w->state;
  const ch_rosenbrock_parameters *p = &w->table;
  size_t n = s->model->nv*s->cells;
  if (evaluate_initial_force) { ch_force(s,s->y,s->forcing); ++stats->function_calls; }
  for (int stage = 0; stage < p->stages; ++stage) {
    double *k = s->stages+(size_t)stage*n;
    int offset = stage*(stage-1)/2;
    if (!stage) memcpy(k,s->forcing,n*sizeof(double));
    else if (p->new_function[stage]) {
      memcpy(s->new_y,s->y,n*sizeof(double));
      for (int j = 0; j < stage; ++j) add(n,p->a[offset+j],s->stages+(size_t)j*n,s->new_y);
      ch_force(s,s->new_y,k); ++stats->function_calls;
    }
    if (stage+1 < p->stages && !p->new_function[stage+1]) memcpy(k+n,k,n*sizeof(double));
    for (int j = 0; j < stage; ++j) add(n,p->c[offset+j]/h,s->stages+(size_t)j*n,k);
    ch_linear_state(s,k); ++stats->solves;
  }
  memcpy(s->new_y,s->y,n*sizeof(double)); memset(s->error,0,n*sizeof(double));
  for (int stage = 0; stage < p->stages; ++stage) add(n,p->m[stage],s->stages+(size_t)stage*n,s->new_y);
  for (int stage = 0; stage < p->stages; ++stage) add(n,p->e[stage],s->stages+(size_t)stage*n,s->error);
}

static int sdirk_trial(ch_integrator *w, double h, ch_stats *stats, ch_integrator_stats *extra) {
  ch_state *s = w->state;
  size_t n = s->model->nv*s->cells;
  const double gamma = 1.0-1.0/sqrt(2.0), alpha = 1.0/(gamma*h);
  double *first_increment = s->stages, *base = s->stages+n;
  memcpy(s->new_y,s->y,n*sizeof(double)); memcpy(base,s->y,n*sizeof(double));
  for (int stage = 0; stage < 2; ++stage) {
    int converged = 0;
    for (int iteration = 0; iteration < w->options.newton_iterations; ++iteration) {
      ch_force(s,s->new_y,s->forcing); ++stats->function_calls;
      for (size_t k = 0; k < n; ++k) s->forcing[k] -= alpha*(s->new_y[k]-base[k]);
      ch_linear_state(s,s->forcing); ++stats->solves; ++extra->nonlinear_iterations;
      add(n,1.0,s->forcing,s->new_y);
      double error = norm(w,s->forcing,s->new_y);
      if (!isfinite(error)) break;
      if (error <= w->options.newton_tolerance) { converged = 1; break; }
    }
    if (!converged) { ++extra->nonlinear_failures; return 0; }
    if (!stage) {
      for (size_t k = 0; k < n; ++k) {
        first_increment[k] = s->new_y[k]-s->y[k];
        base[k] = s->y[k]+((1-gamma)/gamma)*first_increment[k];
      }
    }
  }
  for (size_t k = 0; k < n; ++k)
    s->error[k] = first_increment[k]/gamma-(s->new_y[k]-s->y[k]);
  return 1;
}

int ch_integrator_step(ch_integrator *w, double dt, ch_stats *stats, ch_integrator_stats *extra) {
  if (!w || !stats || !extra || !isfinite(dt) || dt < 0) return CH_INVALID;
  *stats = (ch_stats){0}; *extra = (ch_integrator_stats){0};
  ch_state *s = w->state;
  const ch_integrator_options *o = &w->options;
  size_t n = s->model->nv*s->cells, nj = s->model->sparse->nnz*s->cells;
  for (size_t k = 0; k < n; ++k) if (!isfinite(s->y[k]) || s->y[k] < 0) return CH_INVALID;
  if (dt == 0) return CH_OK;
  double hmax = o->h_max > 0 ? fmin(dt,o->h_max) : dt;
  double h = fmin(hmax,o->warm_start && w->next_h > 0 ? w->next_h : o->h_start);
  double time = 0, factor_h = 0, previous_error = 1, previous_h = h;
  int age = o->jacobian_steps, have_factor = 0, have_force = 0, rejected = 0, status = CH_OK;
  double gamma = o->method == 102 ? 1.0-1.0/sqrt(2.0) : w->table.gamma[0];
  double order = o->method == 102 ? 2.0 : w->table.estimator_order;
  while (time < dt) {
    if (stats->internal_steps >= (size_t)o->max_steps) { status = CH_MAX_STEPS; break; }
    double proposed = h;
    h = fmin(h,dt-time);
    if (h <= 0 || time+h == time) { status = CH_SMALL_STEP; break; }
    int fresh = age >= o->jacobian_steps;
    if (fresh) {
      ch_jacobian(s,s->y,w->raw_jacobian); ++stats->jacobian_updates;
      age = 0; have_factor = 0;
    } else ++extra->jacobian_reuses;
    if (!have_factor || factor_h != h) {
      memcpy(s->jacobian,w->raw_jacobian,nj*sizeof(double));
      ch_alpha_minus_jacobian(s,1.0/(gamma*h),s->jacobian);
      status = ch_factor_state(s); ++stats->decompositions;
      if (status && status != CH_SINGULAR) break;
      have_factor = status == CH_OK; factor_h = h;
    } else ++extra->factor_reuses;
    ++stats->internal_steps;
    int trial_ok = have_factor;
    if (trial_ok) {
      if (o->method == 102) trial_ok = sdirk_trial(w,h,stats,extra);
      else { rosenbrock_trial(w,h,stats,!have_force); have_force = 1; }
    }
    int negative = 0;
    size_t clips = 0; double clipped_amount = 0, maximum_clip = 0;
    if (trial_ok) for (size_t i = 0; i < s->model->nv; ++i)
      for (size_t c = 0; c < s->cells; ++c) {
        size_t k = i*s->cells+c;
        if (s->new_y[k] < -o->positivity*s->atol[i]) negative = 1;
        else if (s->new_y[k] < 0) {
          double correction = -s->new_y[k];
          /* The actual candidate includes the projection. Charge its change to
           * the embedded error before deciding whether to accept this trial. */
          s->error[k] += correction; s->new_y[k] = 0;
          ++clips; clipped_amount += correction; maximum_clip = fmax(maximum_clip,correction);
        }
      }
    double error = trial_ok ? norm(w,s->error,s->new_y) : INFINITY;
    if (!isfinite(error)) ++extra->nonfinite_rejections;
    if (negative) ++extra->positivity_rejections;
    double safe_error = fmax(error,1e-12);
    double factor = o->safety*pow(safe_error,-1.0/order);
    if (o->controller == 1 && stats->accepted && !rejected && error <= 1)
      factor = o->safety*pow(safe_error,-(1.0-o->pi_beta)/order)*pow(previous_error,o->pi_beta/order);
    if (o->controller == 2 && stats->accepted && !rejected && error <= 1) {
      double k = o->filter_k > 0 ? o->filter_k : order;
      factor = pow(o->safety,2.0/(o->filter_b*k))*pow(safe_error*previous_error,-1.0/(o->filter_b*k))*
               pow(h/previous_h,-1.0/o->filter_b);
    }
    factor = fmin(o->factor_max,fmax(o->factor_min,factor));
    if (isfinite(error) && error <= 1 && !negative) {
      extra->clipped_values += clips; extra->clipped_amount += clipped_amount;
      extra->maximum_clip = fmax(extra->maximum_clip,maximum_clip);
      double *swap = s->y; s->y = s->new_y; s->new_y = swap;
      have_force = 0;
      ++stats->accepted; ++age;
      time = h >= dt-time ? dt : time+h;
      if (rejected) factor = fmin(1.0,factor);
      double next = fmin(hmax,h*factor);
      /* Hold h only inside the controller's allowed growth window. This
       * permits exact reuse of a factorization of a lagged Jacobian. */
      if (o->lu_band > 1 && age < o->jacobian_steps &&
          factor >= 1.0/o->lu_band && factor <= o->lu_band) next = h;
      w->next_h = time == dt && h < proposed && error < 0.5 ? fmax(next,proposed) : next;
      previous_h = h; h = next; rejected = 0; previous_error = fmax(error,1e-12); status = CH_OK;
    } else {
      ++stats->rejected; rejected = 1;
      /* A stale matrix gets one retry with the exact current Jacobian. */
      if (fresh) h *= negative || !isfinite(error) ? o->factor_min : fmin(0.9,factor);
      age = o->jacobian_steps; have_factor = 0; status = CH_OK;
    }
  }
  stats->final_time = time;
  if (status) w->next_h = 0;
  ch_post_solve(s);
  return status;
}
