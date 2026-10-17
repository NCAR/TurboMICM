/* The Rosenbrock integrator (MICM's AbstractRosenbrockSolver::Solve), for
 * ODE systems. Each state is one group of cells with one step size. */
#include "internal.h"

#include <math.h>

static const double DEFAULT_H_MIN = 1.0e-15;  /* relative to the time step */
static const double DEFAULT_H_START = 1.0e-6; /* relative to the time step */

/* The forcing F(y) into `forcing` (overwritten). */
static void forcing(const tmc_solver *s, tmc_state *st, const tmc_dense *y, tmc_dense *f)
{
  if (s->kernels)
    s->kernels->forcing(y->data, st->rate_constants->data, f->data, st->n_blocks);
  else
  {
    tmc_dense_fill(f, 0.0);
    tmc_generic_forcing(s, st->rate_constants, y, f);
  }
}

/* alpha I - J(y) into the state's LU matrix, then its factorization. */
static void factor(const tmc_solver *s, tmc_state *st, const tmc_dense *y, double alpha, tmc_stats *stats)
{
  if (s->kernels)
  {
    s->kernels->jacobian(y->data, st->rate_constants->data, alpha, st->jacobian->data, st->n_blocks);
    s->kernels->lu_decompose(st->jacobian->data, st->reciprocals, st->n_blocks);
  }
  else
  {
    tmc_sparse_fill(st->jacobian, 0.0);
    tmc_generic_jacobian(s, st->rate_constants, y, st->jacobian);
    /* AlphaMinusJacobian: alpha on the diagonal of every block */
    const int L = st->L, nnz = st->jacobian->nnz;
    for (int g = 0; g < st->jacobian->n_groups; ++g)
    {
      double *J = st->jacobian->data + (size_t)g * (size_t)nnz * L;
      for (int i = 0; i < st->n_vars; ++i)
      {
        double *d = J + (size_t)s->lu->mz_diag[i] * L;
        for (int c = 0; c < L; ++c)
          d[c] += alpha;
      }
    }
    tmc_lu_decompose_in_place(s->lu, st->jacobian, st->reciprocals);
  }
  stats->jacobian_updates += 1;
  stats->decompositions += 1;
}

static void solve(const tmc_solver *s, tmc_state *st, tmc_dense *x, tmc_stats *stats)
{
  if (s->kernels)
    s->kernels->lu_solve(x->data, st->jacobian->data, st->reciprocals, st->n_blocks);
  else
    tmc_lu_solve_in_place(s->lu, st->jacobian, NULL, x);
  stats->solves += 1;
}

/* out = base + sum_j coef[j] k[j], element by element with one fused
 * multiply-add per term in order, as MICM's sequence of Axpy calls computes
 * it. base may be NULL (zero) or out itself. */
static void combine(tmc_dense *out, const tmc_dense *base, int n, const double *coef, tmc_dense *const *k)
{
  const size_t len = out->size;
  double *restrict o = out->data;
  const double *restrict b = base ? base->data : NULL;
  const double *restrict k0 = n > 0 ? k[0]->data : NULL;
  const double *restrict k1 = n > 1 ? k[1]->data : NULL;
  const double *restrict k2 = n > 2 ? k[2]->data : NULL;
  switch (n)
  {
    case 1:
      if (b)
      {
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
          o[i] = b[i] + coef[0] * k0[i];
      }
      else
      {
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
          o[i] = coef[0] * k0[i];
      }
      return;
    case 2:
      if (b)
      {
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
        {
          double t = b[i] + coef[0] * k0[i];
          o[i] = t + coef[1] * k1[i];
        }
      }
      else
      {
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
        {
          double t = coef[0] * k0[i];
          o[i] = t + coef[1] * k1[i];
        }
      }
      return;
    case 3:
      if (b)
      {
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
        {
          double t = b[i] + coef[0] * k0[i];
          t = t + coef[1] * k1[i];
          o[i] = t + coef[2] * k2[i];
        }
      }
      else
      {
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
        {
          double t = coef[0] * k0[i];
          t = t + coef[1] * k1[i];
          o[i] = t + coef[2] * k2[i];
        }
      }
      return;
    default:
      if (b && b != o)
        memcpy(o, b, len * sizeof(double));
      else if (!b)
        memset(o, 0, len * sizeof(double));
      for (int j = 0; j < n; ++j)
      {
        const double *restrict kj = k[j]->data;
        const double a = coef[j];
#pragma omp simd
        for (size_t i = 0; i < len; ++i)
          o[i] += a * kj[i];
      }
      return;
  }
}

/* MICM's normalized error: the root mean square over all cells and variables
 * of error / (atol + rtol max(|y|, |ynew|)), accumulated variable by variable
 * and cell by cell, padded lanes excluded, floored at 1e-10. */
static double normalized_error(const tmc_state *st, const tmc_dense *y, const tmc_dense *ynew, const tmc_dense *errors)
{
  const int L = st->L, n_vars = st->n_vars, n_cells = st->n_cells;
  const double rtol = st->relative_tolerance;
  const int n_blocks = st->n_blocks;
  double acc = 0.0;
  for (int v = 0; v < n_vars; ++v)
  {
    const double atol = st->absolute_tolerance[v];
    for (int block = 0; block < n_blocks; ++block)
    {
      const size_t base = ((size_t)block * n_vars + v) * L;
      const int lanes = block == n_blocks - 1 ? n_cells - block * L : L;
      for (int lane = 0; lane < lanes; ++lane)
      {
        double a = fabs(y->data[base + lane]), b = fabs(ynew->data[base + lane]);
        double ymax = a > b ? a : b;
        double e = errors->data[base + lane] / (atol + rtol * ymax);
        acc += e * e;
      }
    }
  }
  size_t N = (size_t)n_cells * (size_t)n_vars;
  if (N == 0)
    N = 1;
  double error = sqrt(acc / (double)N);
  return error > 1.0e-10 ? error : 1.0e-10;
}

tmc_result tmc_rosenbrock_solve(tmc_solver *s, double time_step, tmc_state *st)
{
  const tmc_rosenbrock_parameters *p = &s->ros;
  tmc_result result;
  memset(&result, 0, sizeof result);
  result.state = TMC_RUNNING;
  tmc_stats *stats = &result.stats;

  tmc_dense *Y = st->variables;
  tmc_dense *Ynew = st->ynew;
  tmc_dense *initial_forcing = st->initial_forcing;
  tmc_dense *Yerror = st->yerror;
  tmc_dense **K = st->K;

  const double h_min = p->h_min == 0.0 ? DEFAULT_H_MIN * time_step : p->h_min;
  const double h_max = p->h_max == 0.0 ? time_step : fmin(time_step, p->h_max);
  const double h_start = p->h_start == 0.0 ? DEFAULT_H_START * time_step : fmin(h_max, p->h_start);
  double H = fmin(fmax(h_min, fabs(h_start)), fabs(h_max));

  double present_time = 0.0;
  bool reject_last_h = false;
  bool reject_more_h = false;

  while ((present_time - time_step + p->round_off) <= 0 && result.state == TMC_RUNNING)
  {
    if (stats->number_of_steps > p->max_number_of_steps)
    {
      result.state = TMC_CONVERGENCE_EXCEEDED_MAX_STEPS;
      break;
    }
    if (((present_time + 0.1 * H) == present_time) || (H <= fmin(p->round_off, h_min)))
    {
      result.state = TMC_STEP_SIZE_TOO_SMALL;
      break;
    }
    H = fmin(H, fabs(time_step - present_time));

    forcing(s, st, Y, initial_forcing);
    stats->function_calls += 1;

    bool accepted = false;
    while (!accepted)
    {
      const double alpha = 1.0 / (H * p->gamma[0]);
      factor(s, st, Y, alpha, stats);

      /* The stages. Stage s solves (alpha I - J) K_s = F(Y + sum_j a_sj K_j) + sum_j (c_sj / H) K_j.
       * A stage that reuses the previous stage's forcing gets it before the
       * c terms are added, as MICM's K[stage + 1].Copy(K[stage]) does. */
      for (int stage = 0; stage < p->stages; ++stage)
      {
        const int stage_combinations = stage * (stage - 1) / 2;
        const bool reused_by_next = stage + 1 < p->stages && !p->new_function_evaluation[stage + 1];
        tmc_dense *raw = reused_by_next ? K[stage + 1] : K[stage]; /* where the raw forcing goes */
        if (stage == 0)
        {
          if (raw != initial_forcing)
            tmc_dense_copy_into(raw, initial_forcing);
        }
        else if (p->new_function_evaluation[stage])
        {
          combine(Ynew, Y, stage, &p->a[stage_combinations], K);
          forcing(s, st, Ynew, raw);
          stats->function_calls += 1;
        }
        else if (raw != K[stage])
        {
          tmc_dense_copy_into(raw, K[stage]);
        }
        if (stage > 0)
        {
          double c_over_h[TMC_MAX_STAGES];
          for (int j = 0; j < stage; ++j)
            c_over_h[j] = p->c[stage_combinations + j] / H;
          combine(K[stage], raw, stage, c_over_h, K);
        }
        else if (raw != K[0])
        {
          tmc_dense_copy_into(K[0], raw);
        }
        solve(s, st, K[stage], stats);
      }

      combine(Ynew, Y, p->stages, p->m, K);
      combine(Yerror, NULL, p->stages, p->e, K);

      const double error = normalized_error(st, Y, Ynew, Yerror);
      const double fac = s->controller(error, p, s->controller_context);
      double Hnew = H * fac;
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
      if ((error < 1) || (H < h_min))
      {
        stats->accepted += 1;
        present_time = present_time + H;
        tmc_dense_swap(Y, Ynew);
        Hnew = fmax(h_min, fmin(Hnew, h_max));
        if (reject_last_h)
          Hnew = fmin(Hnew, H);
        reject_last_h = false;
        reject_more_h = false;
        H = Hnew;
        accepted = true;
      }
      else
      {
        if (reject_more_h)
          Hnew = H * p->rejection_factor_decrease;
        reject_more_h = reject_last_h;
        reject_last_h = true;
        H = Hnew;
        if (stats->accepted >= 1)
          stats->rejected += 1;
        /* the in-place factorization overwrote the Jacobian; the next attempt
         * rebuilds it (counted when it does) */
      }
    }
  }

  if (result.state == TMC_RUNNING)
    result.state = TMC_CONVERGED;
  stats->final_time = present_time;
  return result;
}
