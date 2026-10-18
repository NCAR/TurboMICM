/* Adaptive ODE experiments, independently implemented from the formulas cited
   in REFERENCES.md. Negative-Jacobian storage follows the existing C library. */
#include "internal.h"
#include "chem_integrator.h"
#include <float.h>

struct ck_integrator
{
  ck_state *state;
  ck_integrator_options o;
  ck_rosenbrock_parameters p;
  double *jacobian, *table[3];
  double proposed;
};

ck_integrator_options ck_integrator_defaults(void)
{
  return (ck_integrator_options){.method = 4, .controller = 0, .norm = 1,
      .jacobian_steps = 1, .max_steps = 100000, .start = 0.01, .safety = 0.9,
      .growth = 6, .shrink = 0.1, .deadband = 0, .warm_start = 1, .positivity = 1, .reuse_rhs = 1};
}
void ck_integrator_free(ck_integrator *v)
{
  if (!v) return;
  free(v->jacobian);
  for (size_t j = 0; j < 3; ++j) free(v->table[j]);
  free(v);
}
ck_status ck_integrator_create(ck_state *s, const ck_integrator_options *o, ck_integrator **out)
{
  if (!out) return CK_INVALID;
  *out = NULL;
  if (!s || !o || s->model->nconstraints || !o->jacobian_steps || !o->max_steps ||
      o->controller < 0 || o->controller > 1 || o->norm < 0 || o->norm > 2 ||
      (o->warm_start != 0 && o->warm_start != 1) ||
      (o->positivity != 0 && o->positivity != 1) ||
      (o->reuse_rhs != 0 && o->reuse_rhs != 1) ||
      !(s->relative_tolerance >= 0) || !isfinite(s->relative_tolerance) ||
      !(o->start > 0) || !isfinite(o->start) || !(o->safety > 0 && o->safety <= 1) ||
      !(o->growth >= 1) || !isfinite(o->growth) || !(o->shrink > 0 && o->shrink < 1) ||
      !(o->deadband >= 0 && o->deadband < 1)) return CK_INVALID;
  for (size_t j = 0; j < s->model->n; ++j)
    if (!(s->absolute_tolerance[j] > 0) || !isfinite(s->absolute_tolerance[j])) return CK_INVALID;
  ck_integrator *v = ck_allocate(1, sizeof *v);
  if (!v) return CK_MEMORY;
  v->state = s;
  v->o = *o;
  if (o->method != 12 && o->method != 13 && ck_rosenbrock_parameters_init(o->method, &v->p))
  { ck_integrator_free(v); return CK_INVALID; }
  v->jacobian = ck_allocate(s->cells * s->model->nnz, sizeof(double));
  for (size_t j = 0; j < 3; ++j)
    v->table[j] = ck_allocate(s->cells * s->model->n, sizeof(double));
  if (!v->jacobian || !v->table[0] || !v->table[1] || !v->table[2])
  { ck_integrator_free(v); return CK_MEMORY; }
  *out = v;
  return CK_OK;
}
static void add(size_t n, double a, const double *restrict x, double *restrict y)
{
#pragma omp simd
  for (size_t i = 0; i < n; ++i) y[i] += a * x[i];
}
static double norm(const ck_integrator *v, const double *next, const double *error)
{
  const ck_state *s = v->state;
  size_t n = s->model->n, cells = s->cells;
  double total = 0, worst = 0, negative = 0;
  for (size_t c = 0; c < cells; ++c)
  {
    double sum = 0;
    for (size_t j = 0; j < n; ++j)
    {
      size_t i = j * cells + c;
      if (!isfinite(next[i]) || !isfinite(error[i])) return INFINITY;
      double scale = s->absolute_tolerance[j] + s->relative_tolerance * fmax(fabs(s->y[i]), fabs(next[i]));
      double e = error[i] / scale;
      sum += e * e;
      if (v->o.norm == 2) worst = fmax(worst, fabs(e));
      if (v->o.positivity) negative = fmax(negative, -next[i] / scale);
    }
    total += sum;
    if (v->o.norm == 1) worst = fmax(worst, sqrt(sum / n));
  }
  if (v->o.norm == 0) worst = sqrt(total / (n * cells));
  return fmax(negative, fmax(worst, 1e-12));
}
static void factor(ck_integrator *v, double alpha, ck_stats *st)
{
  ck_state *s = v->state;
  ck_workspace *w = s->workspace;
  memcpy(w->jacobian, v->jacobian, s->cells * s->model->nnz * sizeof(double));
  ck_alpha_minus_jacobian(s->model, s->cells, w->jacobian, alpha);
  ck_factor(s->model, s->cells, w->jacobian, CK_LU_MOZART);
  ++st->decompositions;
}
static void rosenbrock(ck_integrator *v, double h, int refactor, ck_stats *st)
{
  ck_state *s = v->state;
  const ck_model *m = s->model;
  ck_workspace *w = s->workspace;
  const ck_rosenbrock_parameters *p = &v->p;
  size_t n = m->n * s->cells;
  if (refactor) factor(v, 1 / (h * p->gamma[0]), st);
  for (size_t k = 0; k < p->stages; ++k)
  {
    size_t off = k * (k - 1) / 2;
    if (k == 0) memcpy(w->stage[k], w->initial, n * sizeof(double));
    else if (p->new_function[k])
    {
      memcpy(w->next, s->y, n * sizeof(double));
      for (size_t j = 0; j < k; ++j) add(n, p->a[off+j], w->stage[j], w->next);
      ck_forcing(s, w->next, w->stage[k]);
      ++st->function_calls;
    }
    if (k+1 < p->stages && !p->new_function[k+1])
      memcpy(w->stage[k+1], w->stage[k], n * sizeof(double));
    for (size_t j = 0; j < k; ++j) add(n, p->c[off+j] / h, w->stage[j], w->stage[k]);
    ck_linear_solve(m, s->cells, w->jacobian, w->stage[k]);
    ++st->solves;
  }
  memcpy(w->next, s->y, n * sizeof(double));
  memset(w->error, 0, n * sizeof(double));
  for (size_t k = 0; k < p->stages; ++k)
  {
    add(n, p->m[k], w->stage[k], w->next);
    add(n, p->e[k], w->stage[k], w->error);
  }
}
static void extrapolate(ck_integrator *v, double h, ck_stats *st)
{
  ck_state *s = v->state;
  ck_workspace *w = s->workspace;
  size_t size = s->cells * s->model->n;
  int order = v->o.method - 10;
  /* Richardson tableau at n=1,2,3 substeps, frozen J over the macro step.
     Each sequence reuses its LU exactly. The final two columns estimate LTE. */
  for (int k = 0; k < order; ++k)
  {
    int steps = k+1;
    double sub = h / steps;
    factor(v, 1 / sub, st);
    memcpy(w->next, s->y, size * sizeof(double));
    for (int j = 0; j < steps; ++j)
    {
      if (j == 0) memcpy(w->stage[0], w->initial, size * sizeof(double));
      else { ck_forcing(s, w->next, w->stage[0]); ++st->function_calls; }
      ck_linear_solve(s->model, s->cells, w->jacobian, w->stage[0]);
      ++st->solves;
      add(size, 1, w->stage[0], w->next);
    }
    for (int j = 0; j < k; ++j)
    {
      double ratio = (double)steps / (k-j) - 1;
      for (size_t i = 0; i < size; ++i)
      {
        double old = w->next[i];
        w->next[i] += (old - v->table[j][i]) / ratio;
        v->table[j][i] = old;
      }
    }
    if (k == order-1)
      for (size_t i = 0; i < size; ++i) w->error[i] = w->next[i] - v->table[k-1][i];
    memcpy(v->table[k], w->next, size * sizeof(double));
  }
}
ck_status ck_integrator_advance(ck_integrator *v, double dt, ck_stats *st)
{
  if (!v || !st || !(dt >= 0) || !isfinite(dt)) return CK_INVALID;
  memset(st, 0, sizeof *st);
  ck_state *s = v->state;
  ck_workspace *w = s->workspace;
  ck_status status = ck_external_initialize(s);
  if (status || dt == 0) return status;
  double h = fmin(dt, v->o.warm_start && v->proposed > 0 ? v->proposed : v->o.start);
  double time = 0, previous = 1, factored_h = 0;
  size_t age = v->o.jacobian_steps;
  int rejected = 0, fresh = 0, rhs_valid = 0;
  while (time < dt)
  {
    if (st->internal_steps >= v->o.max_steps) return CK_MAX_STEPS;
    h = fmin(h, dt-time);
    if (!(h > 0) || time + h == time) return CK_STEP_SMALL;
    if (!rhs_valid || !v->o.reuse_rhs)
    {
      ck_forcing(s, s->y, w->initial);
      ++st->function_calls;
      rhs_valid = 1;
    }
    if (age >= v->o.jacobian_steps)
    {
      ck_jacobian(s, s->y, v->jacobian);
      ++st->jacobian_updates;
      age = 0; fresh = 1; factored_h = 0;
    }
    if ((v->o.method == 12 || v->o.method == 13)) extrapolate(v, h, st);
    else rosenbrock(v, h, h != factored_h, st);
    factored_h = h;
    double e = norm(v, w->next, w->error);
    ++st->internal_steps;
    double order = (v->o.method == 12 || v->o.method == 13) ? v->o.method - 10 : v->p.order;
    double factor_h = v->o.safety * pow(e, -1 / order);
    if (v->o.controller && !rejected && st->accepted)
      factor_h = v->o.safety * pow(e, -0.7 / order) * pow(previous, 0.4 / order);
    factor_h = fmin(v->o.growth, fmax(v->o.shrink, factor_h));
    if (e <= 1)
    {
      ++st->accepted;
      time = h >= dt-time ? dt : time+h;
      double *swap = s->y; s->y = w->next; w->next = swap;
      rhs_valid = 0;
      if (v->o.positivity) ck_post_solve_clamp(s);
      ++age; fresh = 0;
      if (rejected) factor_h = fmin(1, factor_h);
      previous = fmax(e, 1e-12); rejected = 0;
      /* Holding h also permits exact LU reuse when J has not expired. */
      if (fabs(factor_h-1) <= v->o.deadband) factor_h = 1;
      v->proposed = h * factor_h;
      h = v->proposed;
    }
    else
    {
      ++st->rejected;
      /* Retry a stale Jacobian with a current one before shrinking h. */
      if (!fresh) age = v->o.jacobian_steps;
      else h *= fmin(0.5, factor_h);
      rejected = 1;
    }
    st->final_time = time;
  }
  ck_post_solve_clamp(s);
  return CK_OK;
}
