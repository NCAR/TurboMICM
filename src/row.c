/* A Rosenbrock-Wanner loop for the scenario solver. It runs MICM's stage
 * formulas (tmc_rosenbrock_stages) with the choices MICM's loop fixes:
 * the first step carried over from the state's previous solve, a choice of
 * error norms and step-size controllers, clamping after each accepted step,
 * steps rounded to tile the solve interval, and reuse of a factorization
 * across steps and solves, which W-methods allow without losing their order. */
#include "internal.h"

#include <math.h>

static const double H_MIN_FRACTION = 1.0e-15;  /* of the solve interval, as MICM */
static const double H_START_FRACTION = 1.0e-6; /* MICM's first step */
static const double SAME_STEP = 1.0e-12;       /* relative difference of steps counted as equal */

double tmc_row_quantize(double h, double interval)
{
  if (h >= interval)
    return interval;
  double n = ceil(interval / h - 1.0e-9);
  if (n < 1.0)
    n = 1.0;
  return interval / n;
}

double tmc_row_initial_step(const tmc_rosenbrock_parameters *p, const tmc_state *st, double time_step, double h_min, double h_max)
{
  double H;
  if (p->warm_start && st->h_next > 0.0)
    H = st->h_next;
  else if (p->h_start != 0.0)
    H = p->h_start;
  else
    H = (p->h_start_factor > 0.0 ? p->h_start_factor : H_START_FRACTION) * time_step;
  H = fmin(fmax(h_min, fabs(H)), h_max);
  if (p->quantize)
    H = tmc_row_quantize(H, time_step);
  return H;
}

/* The factor Hnew / H: MICM's elementary rule, Gustafsson's predictive
 * controller (the smaller of the elementary and predicted factors, as RODAS
 * and RADAU5 take it), or Soderlind's H211b filter; all within MICM's bounds. */
double tmc_row_step_factor(const tmc_rosenbrock_parameters *p, double error, double H, const tmc_state *st)
{
  const double order = p->estimator_of_local_order;
  double fac = p->safety_factor / pow(error, 1.0 / order);
  const bool history = st->err_prev > 0.0 && st->h_prev > 0.0;
  if (p->controller_type == TMC_CONTROLLER_PREDICTIVE && history)
  {
    const double predicted = p->safety_factor * (H / st->h_prev) * pow(st->err_prev / (error * error), 1.0 / order);
    if (predicted < fac)
      fac = predicted;
  }
  else if (p->controller_type == TMC_CONTROLLER_H211B && history)
  {
    /* h_new / h = (1 / err)^(1 / (b k)) (1 / err_prev)^(1 / (b k)) (h / h_prev)^(-1 / b) */
    const double b = p->h211_b > 0.0 ? p->h211_b : 4.0, k = p->h211_k > 0.0 ? p->h211_k : order;
    fac = p->safety_factor * pow(error * st->err_prev, -1.0 / (b * k)) * pow(H / st->h_prev, -1.0 / b);
  }
  if (fac < p->factor_min)
    fac = p->factor_min;
  if (fac > p->factor_max)
    fac = p->factor_max;
  return fac;
}

/* MICM's norm (the root mean square over the group), the largest of the
 * cells' root mean squares, or the largest scaled error; floored at 1e-10. */
double tmc_row_error_norm(const tmc_state *st, const tmc_dense *y, const tmc_dense *ynew, const tmc_dense *errors, int norm)
{
  if (norm == TMC_NORM_RMS)
    return tmc_rosenbrock_normalized_error(st, y, ynew, errors);
  const int L = st->L, n_vars = st->n_vars, n_cells = st->n_cells, n_blocks = st->n_blocks;
  const double rtol = st->relative_tolerance;
  double worst = 0.0;
  for (int block = 0; block < n_blocks; ++block)
  {
    const int lanes = block == n_blocks - 1 ? n_cells - block * L : L;
    for (int lane = 0; lane < lanes; ++lane)
    {
      double acc = 0.0, max_e = 0.0;
      for (int v = 0; v < n_vars; ++v)
      {
        const size_t i = ((size_t)block * n_vars + v) * L + lane;
        const double a = fabs(y->data[i]), b = fabs(ynew->data[i]);
        const double e = fabs(errors->data[i]) / (st->absolute_tolerance[v] + rtol * (a > b ? a : b));
        acc += e * e;
        if (e > max_e || isnan(e))
          max_e = e;
      }
      const double cell = norm == TMC_NORM_MAX ? max_e : sqrt(acc / (double)n_vars);
      if (cell > worst || isnan(cell))
        worst = cell;
    }
  }
  return worst > 1.0e-10 ? worst : 1.0e-10;
}

tmc_result tmc_row_solve(const tmc_solver *s, double time_step, tmc_state *st)
{
  const tmc_rosenbrock_parameters *p = &s->ros;
  tmc_result result;
  memset(&result, 0, sizeof result);
  result.state = TMC_RUNNING;
  tmc_stats *stats = &result.stats;
  tmc_dense *Y = st->variables, *Ynew = st->ynew, *Yerror = st->yerror, *F0 = st->initial_forcing;

  const double h_min = p->h_min == 0.0 ? H_MIN_FRACTION * time_step : p->h_min;
  const double h_max = p->h_max == 0.0 ? time_step : fmin(time_step, p->h_max);
  double H = tmc_row_initial_step(p, st, time_step, h_min, h_max);
  double proposed = H; /* the controller's untruncated proposal, carried to the next solve */
  double t = 0.0;
  bool reject_last = false, reject_more = false;

  while (t - time_step + p->round_off <= 0.0 && result.state == TMC_RUNNING)
  {
    if (stats->number_of_steps > p->max_number_of_steps)
    {
      result.state = TMC_CONVERGENCE_EXCEEDED_MAX_STEPS;
      break;
    }
    if (t + 0.1 * H == t || H <= fmin(p->round_off, h_min))
    {
      result.state = TMC_STEP_SIZE_TOO_SMALL;
      break;
    }
    H = fmin(H, fabs(time_step - t));

    /* The step's forcing, and its matrix unless the one in the state may serve again. */
    const bool reuse = p->max_matrix_age > 0 && st->lu_h > 0.0 && fabs(st->lu_h - H) <= SAME_STEP * H &&
                       st->lu_age < p->max_matrix_age && st->rc_change <= p->matrix_rate_change;
    if (reuse)
    {
      tmc_rosenbrock_forcing(s, st, Y, F0);
      stats->matrix_reuses += 1;
    }
    else if (s->kernels && s->kernels->forcing_jacobian)
    {
      s->kernels->forcing_jacobian(Y->data, st->rate_constants->data, 1.0 / (H * p->gamma[0]), F0->data, st->jacobian->data, st->n_blocks);
      s->kernels->lu_decompose(st->jacobian->data, st->reciprocals, st->n_blocks);
      stats->jacobian_updates += 1;
      stats->decompositions += 1;
      st->lu_h = H;
      st->lu_age = 0;
      st->rc_change = 0.0;
    }
    else
    {
      tmc_rosenbrock_forcing(s, st, Y, F0);
      tmc_rosenbrock_factor(s, st, Y, 1.0 / (H * p->gamma[0]), stats);
      st->lu_h = H;
      st->lu_age = 0;
      st->rc_change = 0.0;
    }
    stats->function_calls += 1;
    bool factored = true;

    bool accepted = false;
    while (!accepted)
    {
      if (!factored)
      {
        /* after a rejection: the matrix for the new H, at the same Y */
        tmc_rosenbrock_factor(s, st, Y, 1.0 / (H * p->gamma[0]), stats);
        st->lu_h = H;
        st->lu_age = 0;
        st->rc_change = 0.0;
      }
      factored = false;
      tmc_rosenbrock_stages(s, st, p, H, stats);
      tmc_rosenbrock_combine_solution_and_error(Ynew, Yerror, Y, p->stages, p->m, p->e, st->K);
      if (p->filter_error)
      {
        /* the estimate seen through the step's own matrix, as RADAU5 filters its
         * first-step estimate: damps the stiff components a transient excites */
        const double alpha = 1.0 / (H * p->gamma[0]);
        tmc_rosenbrock_linear_solve(s, st, Yerror, stats);
        double *restrict e = Yerror->data;
        const size_t len = Yerror->size;
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
          e[i] *= alpha;
      }
      const double error = tmc_row_error_norm(st, Y, Ynew, Yerror, p->norm);
      stats->number_of_steps += 1;
      if (isnan(error))
      {
        result.state = TMC_NAN_DETECTED;
        break;
      }
      if (isinf(error))
      {
        result.state = TMC_INF_DETECTED;
        break;
      }
      const double fac = tmc_row_step_factor(p, error, H, st);
      if (error < 1.0 || H < h_min)
      {
        stats->accepted += 1;
        t += H;
        tmc_dense_swap(Y, Ynew);
        if (p->clamp_internal)
          tmc_dense_max(Y, 0.0);
        st->lu_age += 1;
        st->err_prev = error;
        st->h_prev = H;
        double Hnew = H * fac;
        if (reject_last)
          Hnew = fmin(Hnew, H);
        if (p->freeze_band > 1.0 && Hnew < H * p->freeze_band && Hnew > H / p->freeze_band)
          Hnew = H; /* keep the step, and with it the factorization */
        Hnew = fmax(h_min, fmin(Hnew, h_max));
        if (p->quantize)
          Hnew = tmc_row_quantize(Hnew, time_step);
        reject_last = reject_more = false;
        proposed = Hnew;
        H = Hnew;
        accepted = true;
      }
      else
      {
        double Hnew = H * fac;
        if (reject_more)
          Hnew = H * p->rejection_factor_decrease;
        reject_more = reject_last;
        reject_last = true;
        if (p->quantize)
          Hnew = tmc_row_quantize(Hnew, time_step);
        H = Hnew;
        stats->rejected += 1;
        st->lu_h = 0.0; /* the in-place factorization is for the old H */
      }
    }
  }

  if (result.state == TMC_RUNNING)
    result.state = TMC_CONVERGED;
  st->h_next = proposed;
  stats->final_time = t;
  return result;
}
