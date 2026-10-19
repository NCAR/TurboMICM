/* Extrapolation of the linearly implicit Euler method (Deuflhard; Hairer and
 * Wanner, Solving ODEs II, IV.9). One attempt with step H integrates from Y
 * with n_j substeps of H / n_j for each entry of a step-number sequence,
 * the Jacobian evaluated once at Y and the matrix n_j / H I - J factored once
 * per entry, and combines the results by Aitken-Neville extrapolation in
 * powers of H. The error estimate is the difference of the last two entries
 * of the tableau's last row. The step-size rules are the ROW loop's. */
#include "internal.h"

#include <math.h>

static const double H_MIN_FRACTION = 1.0e-15;

static int sequence_entry(int kind, int j)
{
  switch (kind)
  {
    case TMC_SEQUENCE_SHIFTED: return j + 2;      /* 2, 3, 4, 5, ... */
    case TMC_SEQUENCE_EVEN: return 2 * (j + 1);   /* 2, 4, 6, 8, ... */
    default: return j + 1;                        /* 1, 2, 3, 4, ... (harmonic) */
  }
}

/* -J of Y, unshifted, into the state's copy. */
static void assemble_jacobian(const tmc_solver *s, tmc_state *st, const tmc_dense *y, tmc_stats *stats)
{
  if (s->kernels)
    s->kernels->jacobian(y->data, st->rate_constants->data, 0.0, st->jacobian_copy->data, st->n_blocks);
  else
    tmc_eval_jacobian(s, st, y, st->jacobian_copy);
  stats->jacobian_updates += 1;
}

/* alpha I - J from the copy, factored in place. */
static void factor_shifted(const tmc_solver *s, tmc_state *st, double alpha, tmc_stats *stats)
{
  memcpy(st->jacobian->data, st->jacobian_copy->data, st->jacobian->size * sizeof(double));
  tmc_rosenbrock_alpha_minus_jacobian(s, st, alpha);
  if (s->kernels)
    s->kernels->lu_decompose(st->jacobian->data, st->reciprocals, st->n_blocks);
  else
    tmc_lu_decompose_in_place(s->lu, st->jacobian, st->reciprocals);
  stats->decompositions += 1;
}

/* y += x */
static void add_into(tmc_dense *y, const tmc_dense *x)
{
  double *restrict a = y->data;
  const double *restrict b = x->data;
  const size_t len = y->size;
#pragma omp simd
  for (size_t i = 0; i < len; ++i)
    a[i] += b[i];
}

tmc_result tmc_extrapolation_solve(const tmc_solver *s, double time_step, tmc_state *st)
{
  const tmc_rosenbrock_parameters *p = &s->ros;
  tmc_result result;
  memset(&result, 0, sizeof result);
  result.state = TMC_RUNNING;
  tmc_stats *stats = &result.stats;
  const int k = p->extrapolation_order < 2 ? 2 : (p->extrapolation_order > TMC_MAX_EXTRAPOLATION ? TMC_MAX_EXTRAPOLATION : p->extrapolation_order);
  tmc_dense *Y = st->variables, *Ynew = st->ynew, *Yerror = st->yerror, *F0 = st->initial_forcing, *rhs = st->K[0];
  tmc_dense **prev = st->extrap_prev, **cur = st->extrap_cur;

  const double h_min = p->h_min == 0.0 ? H_MIN_FRACTION * time_step : p->h_min;
  const double h_max = p->h_max == 0.0 ? time_step : fmin(time_step, p->h_max);
  double H = tmc_row_initial_step(p, st, time_step, h_min, h_max);
  double proposed = H;
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

    tmc_rosenbrock_forcing(s, st, Y, F0);
    stats->function_calls += 1;
    assemble_jacobian(s, st, Y, stats);

    bool accepted = false;
    while (!accepted)
    {
      for (int j = 0; j < k; ++j)
      {
        const int n_j = sequence_entry(p->extrapolation_sequence, j);
        factor_shifted(s, st, (double)n_j / H, stats);
        /* n_j linearly implicit Euler substeps: y += (alpha I - J)^-1 F(y), alpha = n_j / H */
        tmc_dense_copy_into(Ynew, Y);
        for (int m = 0; m < n_j; ++m)
        {
          if (m == 0)
            tmc_dense_copy_into(rhs, F0);
          else
          {
            tmc_rosenbrock_forcing(s, st, Ynew, rhs);
            stats->function_calls += 1;
          }
          tmc_rosenbrock_linear_solve(s, st, rhs, stats);
          add_into(Ynew, rhs);
        }
        /* the tableau's row j: T_j0 = Ynew, T_jl = T_j,l-1 + (T_j,l-1 - T_j-1,l-1) / (n_j / n_j-l - 1) */
        tmc_dense_copy_into(cur[0], Ynew);
        for (int l = 1; l <= j; ++l)
        {
          const double ratio = (double)n_j / (double)sequence_entry(p->extrapolation_sequence, j - l);
          const double w = 1.0 / (ratio - 1.0);
          double *restrict out = cur[l]->data;
          const double *restrict a = cur[l - 1]->data, *restrict b = prev[l - 1]->data;
          const size_t len = cur[l]->size;
#pragma omp simd
          for (size_t i = 0; i < len; ++i)
            out[i] = a[i] + w * (a[i] - b[i]);
        }
        if (j < k - 1)
          for (int l = 0; l <= j; ++l)
            tmc_dense_swap(cur[l], prev[l]);
      }
      /* the result and its error estimate: the last row's last two entries */
      tmc_dense_copy_into(Ynew, cur[k - 1]);
      {
        double *restrict e = Yerror->data;
        const double *restrict a = cur[k - 1]->data, *restrict b = cur[k - 2]->data;
        const size_t len = Yerror->size;
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
          e[i] = a[i] - b[i];
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
        st->err_prev = error;
        st->h_prev = H;
        double Hnew = H * fac;
        if (reject_last)
          Hnew = fmin(Hnew, H);
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
      }
    }
  }

  if (result.state == TMC_RUNNING)
    result.state = TMC_CONVERGED;
  st->h_next = proposed;
  st->lu_h = 0.0; /* the factorization in the state is a substep's, not reusable */
  stats->final_time = t;
  return result;
}
