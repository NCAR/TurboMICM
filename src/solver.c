/* Rosenbrock formulas and step control translated from MICM, Apache-2.0.
   See NOTICE for attribution. The solver operates only through chemistry and
   linear algebra interfaces, keeping parameter tables/controller replaceable. */
#include "internal.h"

void ck_stats_add(ck_stats *d, const ck_stats *s)
{
  d->internal_steps += s->internal_steps;
  d->accepted += s->accepted;
  d->rejected += s->rejected;
  d->function_calls += s->function_calls;
  d->jacobian_updates += s->jacobian_updates;
  d->decompositions += s->decompositions;
  d->solves += s->solves;
  d->constraint_init_iterations += s->constraint_init_iterations;
}
static void axpy(size_t count, double factor, const double *restrict x, double *restrict y)
{
#pragma omp simd
  for (size_t i = 0; i < count; ++i)
    y[i] += factor * x[i];
}
double ck_normalized_error(const ck_state *s, const double *y, const double *next, const double *errors)
{
  double sum = 0;
  for (size_t j = 0; j < s->model->n; ++j)
  {
    /* Deliberately ordered reduction: no fast-math or SIMD reassociation. */
    for (size_t c = 0; c < s->cells; ++c)
    {
      size_t i = j * s->cells + c;
      double ymax = fabs(y[i]) > fabs(next[i]) ? fabs(y[i]) : fabs(next[i]);
      double scaled = errors[i] / (s->absolute_tolerance[j] + s->relative_tolerance * ymax);
      sum += scaled * scaled;
    }
  }
  double error = sqrt(sum / (s->cells * s->model->n));
  return error > 1e-10 || isnan(error) ? error : 1e-10;
}
double ck_next_step_factor(const ck_rosenbrock_parameters *p, double error)
{
  return fmin(p->factor_max, fmax(p->factor_min, p->safety / pow(error, 1 / p->order)));
}
ck_status ck_solve_raw(ck_state *s, double dt, const ck_rosenbrock_parameters *p, ck_stats *stats)
{
  if (!s || !p || !stats || !(dt >= 0) || !isfinite(dt) || p->stages < 1 || p->stages > 6 ||
      !(p->gamma[0] > 0))
    return CK_INVALID;
  memset(stats, 0, sizeof *stats);
  ck_status init = ck_initialize_constraints(s, p, stats);
  if (init)
    return init;
  if (dt == 0)
    return CK_OK;
  ck_workspace *w = s->workspace;
  const ck_model *m = s->model;
  size_t size = m->n * s->cells;
  const double h_min = p->h_min == 0 ? 1e-15 * dt : p->h_min;
  const double h_max = p->h_max == 0 ? dt : fmin(dt, p->h_max);
  const double h_start = p->h_start == 0 ? 1e-6 * dt : fmin(h_max, p->h_start);
  double h = fmin(fmax(h_min, fabs(h_start)), fabs(h_max)), time = 0;
  int reject_last = 0, reject_more = 0;
  ck_status status = CK_OK;
  while ((time - dt + p->round_off) <= 0)
  {
    if (stats->internal_steps > p->max_steps)
    {
      status = CK_MAX_STEPS;
      break;
    }
    if (time + 0.1 * h == time || h <= fmin(p->round_off, h_min))
    {
      status = CK_STEP_SMALL;
      break;
    }
    h = fmin(h, fabs(dt - time));
    ck_forcing(s, s->y, w->initial);
    ck_constraint_forcing(s, s->y, w->initial);
    ++stats->function_calls;
    ck_jacobian(s, s->y, w->jacobian);
    ck_constraint_jacobian(s, s->y, w->jacobian);
    ++stats->jacobian_updates;
    int accepted = 0;
    while (!accepted)
    {
      double alpha = 1.0 / (h * p->gamma[0]);
      ck_alpha_minus_jacobian(m, s->cells, w->jacobian, alpha);
      ck_dispatch_factor(p->linear, m, s->cells, w->jacobian, p->lu);
      ++stats->decompositions;
      for (size_t stage = 0; stage < p->stages; ++stage)
      {
        size_t offset = stage * (stage - 1) / 2;
        if (stage == 0)
          memcpy(w->stage[stage], w->initial, size * sizeof(double));
        else if (p->new_function[stage])
        {
          memcpy(w->next, s->y, size * sizeof(double));
          for (size_t j = 0; j < stage; ++j)
            axpy(size, p->a[offset + j], w->stage[j], w->next);
          ck_forcing(s, w->next, w->stage[stage]);
          ck_constraint_forcing(s, w->next, w->stage[stage]);
          ++stats->function_calls;
        }
        if (stage + 1 < p->stages && !p->new_function[stage + 1])
          memcpy(w->stage[stage + 1], w->stage[stage], size * sizeof(double));
        for (size_t j = 0; j < stage; ++j)
        {
          if (!m->nconstraints)
            axpy(size, p->c[offset + j] / h, w->stage[j], w->stage[stage]);
          else
            for (size_t v = 0; v < m->n; ++v)
              if (!m->algebraic[v])
                axpy(s->cells, p->c[offset + j] / h, w->stage[j] + v * s->cells,
                     w->stage[stage] + v * s->cells);
        }
        ck_dispatch_solve(p->linear, m, s->cells, w->jacobian, w->stage[stage]);
        ++stats->solves;
      }
      memcpy(w->next, s->y, size * sizeof(double));
      for (size_t stage = 0; stage < p->stages; ++stage)
        axpy(size, p->m[stage], w->stage[stage], w->next);
      memset(w->error, 0, size * sizeof(double));
      for (size_t stage = 0; stage < p->stages; ++stage)
        axpy(size, p->e[stage], w->stage[stage], w->error);
      if (m->nconstraints)
        for (size_t v = 0; v < m->n; ++v)
          if (m->algebraic[v])
            for (size_t c = 0; c < s->cells; ++c)
              w->error[v * s->cells + c] = w->next[v * s->cells + c] - s->y[v * s->cells + c];
      double error = ck_normalized_error(s, s->y, w->next, w->error);
      double next_h = h * (p->step_factor ? p->step_factor(p, error, p->controller_context)
                                          : ck_next_step_factor(p, error));
      ++stats->internal_steps;
      if (isnan(error))
      {
        status = CK_NAN;
        goto done;
      }
      if (isinf(error))
      {
        status = CK_INF;
        goto done;
      }
      if (error < 1 || h < h_min)
      {
        ++stats->accepted;
        time += h;
        double *swap = s->y;
        s->y = w->next;
        w->next = swap;
        next_h = fmax(h_min, fmin(next_h, h_max));
        if (reject_last)
          next_h = fmin(next_h, h);
        reject_last = reject_more = 0;
        h = next_h;
        accepted = 1;
      }
      else
      {
        if (reject_more)
          next_h = h * p->rejection_decrease;
        reject_more = reject_last;
        reject_last = 1;
        h = next_h;
        if (stats->accepted >= 1)
          ++stats->rejected;
        ck_jacobian(s, s->y, w->jacobian);
        ck_constraint_jacobian(s, s->y, w->jacobian);
        ++stats->jacobian_updates;
      }
    }
  }
done:
  stats->final_time = time;
  return status;
}

void ck_post_solve_clamp(ck_state *s)
{
  for (size_t j = 0; j < s->model->n; ++j)
    if (!ck_is_algebraic(s->model, j))
    {
      double *restrict y = s->y + j * s->cells;
#pragma omp simd
      for (size_t c = 0; c < s->cells; ++c)
        y[c] = y[c] > 0 ? y[c] : 0;
    }
}
ck_status ck_solve(ck_state *s, double dt, const ck_rosenbrock_parameters *p, ck_stats *stats)
{
  if (!s || !p || !stats)
    return CK_INVALID;
  memset(stats, 0, sizeof *stats);
  ck_status status = ck_external_initialize(s);
  if (status)
    return status;
  status = ck_solve_raw(s, dt, p, stats);
  if (s)
    ck_post_solve_clamp(s);
  return status;
}
