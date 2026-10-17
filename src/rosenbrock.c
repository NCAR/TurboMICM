/* The Rosenbrock integrator (MICM's AbstractRosenbrockSolver::Solve), for
 * ODE systems. Each state is one group of cells with one step size. */
#include "internal.h"

#include <math.h>

static const double DEFAULT_H_MIN = 1.0e-15;  /* relative to the time step */
static const double DEFAULT_H_START = 1.0e-6; /* relative to the time step */

/* The forcing F(y) into `forcing` (overwritten). */
static void forcing(const tmc_solver *s, tmc_state *st, const tmc_dense *y, tmc_dense *f)
{
  tmc_eval_forcing(s, st, y, f);
}

/* AlphaMinusJacobian: alpha M on the diagonal of every block, M the mass
 * matrix (0 for algebraic rows). */
static void alpha_minus_jacobian(const tmc_solver *s, tmc_state *st, double alpha)
{
  const int L = st->L, nnz = st->jacobian->nnz;
  for (int g = 0; g < st->jacobian->n_groups; ++g)
  {
    double *J = st->jacobian->data + (size_t)g * (size_t)nnz * L;
    for (int i = 0; i < st->n_vars; ++i)
    {
      const double scaled_alpha = alpha * s->mass_diagonal[i];
      double *d = J + (size_t)s->lu->mz_diag[i] * L;
      for (int c = 0; c < L; ++c)
        d[c] += scaled_alpha;
    }
  }
}

/* alpha M - J(y) into the state's LU matrix, then its factorization. */
static void factor(const tmc_solver *s, tmc_state *st, const tmc_dense *y, double alpha, tmc_stats *stats)
{
  if (s->kernels)
  {
    s->kernels->jacobian(y->data, st->rate_constants->data, alpha, st->jacobian->data, st->n_blocks);
    s->kernels->lu_decompose(st->jacobian->data, st->reciprocals, st->n_blocks);
  }
  else
  {
    tmc_eval_jacobian(s, st, y, st->jacobian);
    alpha_minus_jacobian(s, st, alpha);
    tmc_lu_decompose_in_place(s->lu, st->jacobian, st->reciprocals);
  }
  stats->jacobian_updates += 1;
  stats->decompositions += 1;
}

/* K_stage += (c / H) K_j on the differential rows only (MICM's mass_coupling). */
static void add_coupling_masked(const tmc_solver *s, tmc_dense *k_stage, double c_over_h, const tmc_dense *k_j)
{
  const int L = k_stage->L, n_vars = s->n_vars;
  for (int g = 0; g < k_stage->n_groups; ++g)
    for (int v = 0; v < n_vars; ++v)
    {
      if (s->mass_diagonal[v] == 0.0)
        continue;
      double *ks = k_stage->data + ((size_t)g * n_vars + v) * L;
      const double *kj = k_j->data + ((size_t)g * n_vars + v) * L;
      for (int c = 0; c < L; ++c)
        ks[c] += c_over_h * kj[c];
    }
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

/* Ynew = Y + sum m_j K_j and Yerror = sum e_j K_j in one pass over the stage
 * vectors, each as MICM's sequence of fused multiply-adds. */
static void combine_solution_and_error(tmc_dense *ynew, tmc_dense *yerror, const tmc_dense *y, int stages, const double *m, const double *e, tmc_dense *const *k)
{
  if (stages != 3)
  {
    combine(ynew, y, stages, m, k);
    combine(yerror, NULL, stages, e, k);
    return;
  }
  const size_t len = ynew->size;
  double *restrict out = ynew->data, *restrict err = yerror->data;
  const double *restrict base = y->data;
  const double *restrict k0 = k[0]->data, *restrict k1 = k[1]->data, *restrict k2 = k[2]->data;
  const double m0 = m[0], m1 = m[1], m2 = m[2], e0 = e[0], e1 = e[1], e2 = e[2];
#pragma omp simd
  for (size_t i = 0; i < len; ++i)
  {
    double a = k0[i], b = k1[i], c = k2[i];
    double t = base[i] + m0 * a;
    t = t + m1 * b;
    out[i] = t + m2 * c;
    double u = e0 * a;
    u = u + e1 * b;
    err[i] = u + e2 * c;
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

/* ---- Newton initialization of the algebraic variables (MICM's InitializeConstraints) ---- */

typedef struct algebraic_check
{
  bool nan, inf;
  double max_residual;
} algebraic_check;

/* Over the algebraic rows of the active lanes: non-finite flags and the largest magnitude. */
static algebraic_check check_algebraic(const tmc_solver *s, const tmc_state *st, const tmc_dense *m)
{
  algebraic_check r = { false, false, 0.0 };
  const int L = st->L, n_vars = s->n_vars;
  for (int a = 0; a < s->n_algebraic; ++a)
  {
    int v = s->algebraic_ids[a];
    for (int block = 0; block < st->n_blocks; ++block)
    {
      const int lanes = block == st->n_blocks - 1 ? st->n_cells - block * L : L;
      const double *x = m->data + ((size_t)block * n_vars + v) * L;
      for (int lane = 0; lane < lanes; ++lane)
      {
        double val = x[lane];
        if (isnan(val))
          r.nan = true;
        else if (isinf(val))
          r.inf = true;
        else
        {
          double av = fabs(val);
          if (av > r.max_residual)
            r.max_residual = av;
        }
      }
    }
  }
  return r;
}

/* max over algebraic rows and active lanes of |delta| / (atol + rtol max(|z|, |z + delta|)) */
static double weighted_correction(const tmc_solver *s, const tmc_state *st, const tmc_dense *y, const tmc_dense *delta)
{
  double q_max = 0.0;
  const int L = st->L, n_vars = s->n_vars;
  const double rtol = st->relative_tolerance;
  for (int a = 0; a < s->n_algebraic; ++a)
  {
    int v = s->algebraic_ids[a];
    const double atol = st->absolute_tolerance[v];
    for (int block = 0; block < st->n_blocks; ++block)
    {
      const int lanes = block == st->n_blocks - 1 ? st->n_cells - block * L : L;
      size_t base = ((size_t)block * n_vars + v) * L;
      for (int lane = 0; lane < lanes; ++lane)
      {
        double z = y->data[base + lane], correction = delta->data[base + lane];
        double corrected = z + correction;
        double z_max = fabs(z) > fabs(corrected) ? fabs(z) : fabs(corrected);
        double scale = atol + rtol * z_max;
        double q = scale > 0.0 ? fabs(correction) / scale : (correction == 0.0 ? 0.0 : INFINITY);
        if (q > q_max)
          q_max = q;
      }
    }
  }
  return q_max;
}

/* y[alg] += delta[alg] where delta is finite */
static void apply_update(const tmc_solver *s, const tmc_state *st, tmc_dense *y, const tmc_dense *delta)
{
  const int L = st->L, n_vars = s->n_vars;
  for (int a = 0; a < s->n_algebraic; ++a)
  {
    int v = s->algebraic_ids[a];
    for (int g = 0; g < y->n_groups; ++g)
    {
      size_t base = ((size_t)g * n_vars + v) * L;
      for (int lane = 0; lane < L; ++lane)
      {
        double d = delta->data[base + lane];
        if (!isnan(d) && !isinf(d))
          y->data[base + lane] += d;
      }
    }
  }
}

/* The constraint residuals G(y) into delta (overwritten) */
static void constraint_residuals(const tmc_solver *s, tmc_state *st, const tmc_dense *y, tmc_dense *delta)
{
  tmc_dense_fill(delta, 0.0);
  if (s->n_constraints)
    tmc_constraints_add_forcing(s, y, st->custom_rate_parameters, delta);
  if (s->n_external)
    tmc_external_add_constraint_residual(s, st->custom_rate_parameters, y, delta);
}

static void constraint_jacobian_factor(const tmc_solver *s, tmc_state *st, const tmc_dense *y, tmc_stats *stats)
{
  tmc_sparse_fill(st->jacobian, 0.0);
  if (s->n_constraints)
    tmc_constraints_subtract_jacobian(s, y, st->custom_rate_parameters, st->jacobian);
  if (s->n_external)
    tmc_external_subtract_constraint_jacobian(s, st->custom_rate_parameters, y, st->jacobian);
  stats->jacobian_updates += 1;
  alpha_minus_jacobian(s, st, 1.0);
  tmc_lu_decompose_in_place(s->lu, st->jacobian, st->reciprocals);
  stats->decompositions += 1;
}

static tmc_solver_state initialize_constraints(const tmc_solver *s, tmc_state *st, const tmc_rosenbrock_parameters *p, tmc_stats *stats)
{
  tmc_dense *Y = st->variables;
  tmc_dense *delta = st->initial_forcing;   /* residual and Newton correction */
  tmc_dense *original = st->K[0];           /* rollback copy */
  tmc_dense *candidate = st->ynew;
  tmc_dense *trial_delta = st->yerror;
  bool modified = false;
#define RESTORE_AND_RETURN(status)                   \
  do                                                 \
  {                                                  \
    if (modified)                                    \
      tmc_dense_copy_into(Y, original);              \
    return (status);                                 \
  } while (0)

  for (int update = 0; update <= p->constraint_init_max_iterations; ++update)
  {
    constraint_residuals(s, st, Y, delta);
    algebraic_check chk = check_algebraic(s, st, delta);
    if (chk.nan)
      RESTORE_AND_RETURN(TMC_NAN_DETECTED);
    if (chk.inf)
      RESTORE_AND_RETURN(TMC_INF_DETECTED);
    stats->constraint_init_iterations += 1;
    if (chk.max_residual == 0.0)
      return TMC_CONVERGED;

    constraint_jacobian_factor(s, st, Y, stats);
    solve(s, st, delta, stats);
    chk = check_algebraic(s, st, delta);
    if (chk.nan)
      RESTORE_AND_RETURN(TMC_NAN_DETECTED);
    if (chk.inf)
      RESTORE_AND_RETURN(TMC_INF_DETECTED);

    double max_correction = weighted_correction(s, st, Y, delta);
    if (max_correction <= p->constraint_init_tolerance)
      return TMC_CONVERGED;
    if (update == p->constraint_init_max_iterations)
      break;

    if (!modified)
    {
      tmc_dense_copy_into(original, Y);
      modified = true;
    }
    if (p->constraint_init_max_backtracks > 0)
    {
      const double current_norm = max_correction;
      double step = 1.0;
      bool accepted = false;
      for (int backtrack = 0;; ++backtrack, step *= p->constraint_init_backtrack_factor)
      {
        if (backtrack > p->constraint_init_max_backtracks)
          RESTORE_AND_RETURN(TMC_CONSTRAINT_INITIALIZATION_FAILED);
        tmc_dense_fill(trial_delta, 0.0);
        tmc_dense_axpy(trial_delta, step, delta);
        tmc_dense_copy_into(candidate, Y);
        apply_update(s, st, candidate, trial_delta);
        chk = check_algebraic(s, st, candidate);
        if (chk.nan || chk.inf)
          continue;
        constraint_residuals(s, st, candidate, trial_delta);
        chk = check_algebraic(s, st, trial_delta);
        if (chk.nan || chk.inf)
          continue;
        if (chk.max_residual == 0.0)
        {
          accepted = true;
          break;
        }
        solve(s, st, trial_delta, stats);
        chk = check_algebraic(s, st, trial_delta);
        if (chk.nan || chk.inf)
          continue;
        max_correction = weighted_correction(s, st, candidate, trial_delta);
        if (max_correction <= p->constraint_init_tolerance || max_correction < (1 - p->constraint_init_sufficient_decrease * step) * current_norm)
        {
          accepted = true;
          break;
        }
      }
      if (accepted)
      {
        tmc_dense_copy_into(Y, candidate);
        continue; /* recheck convergence at the accepted state with a fresh Jacobian */
      }
    }
    apply_update(s, st, Y, delta);
    chk = check_algebraic(s, st, delta);
    if (chk.nan)
      RESTORE_AND_RETURN(TMC_NAN_DETECTED);
    if (chk.inf)
      RESTORE_AND_RETURN(TMC_INF_DETECTED);
  }
  RESTORE_AND_RETURN(TMC_CONSTRAINT_INITIALIZATION_FAILED);
#undef RESTORE_AND_RETURN
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

  if (s->n_algebraic)
  {
    tmc_solver_state init = initialize_constraints(s, st, p, stats);
    if (init != TMC_CONVERGED)
    {
      result.state = init;
      return result;
    }
  }

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

    bool factored = false;
    if (s->kernels && s->kernels->forcing_jacobian)
    {
      /* the forcing and the shifted Jacobian of Y in one pass, for the first attempt's H */
      const double alpha = 1.0 / (H * p->gamma[0]);
      s->kernels->forcing_jacobian(Y->data, st->rate_constants->data, alpha, initial_forcing->data, st->jacobian->data, st->n_blocks);
      s->kernels->lu_decompose(st->jacobian->data, st->reciprocals, st->n_blocks);
      stats->jacobian_updates += 1;
      stats->decompositions += 1;
      factored = true;
    }
    else
      forcing(s, st, Y, initial_forcing);
    stats->function_calls += 1;

    bool accepted = false;
    while (!accepted)
    {
      const double alpha = 1.0 / (H * p->gamma[0]);
      if (!factored)
        factor(s, st, Y, alpha, stats);
      factored = false;

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
        if (stage > 0 && s->n_algebraic == 0)
        {
          double c_over_h[TMC_MAX_STAGES];
          for (int j = 0; j < stage; ++j)
            c_over_h[j] = p->c[stage_combinations + j] / H;
          combine(K[stage], raw, stage, c_over_h, K);
        }
        else if (stage > 0)
        {
          if (raw != K[stage])
            tmc_dense_copy_into(K[stage], raw);
          for (int j = 0; j < stage; ++j)
            add_coupling_masked(s, K[stage], p->c[stage_combinations + j] / H, K[j]);
        }
        else if (raw != K[0])
        {
          tmc_dense_copy_into(K[0], raw);
        }
        solve(s, st, K[stage], stats);
      }

      combine_solution_and_error(Ynew, Yerror, Y, p->stages, p->m, p->e, K);
      /* For DAE systems the embedded error is near zero on algebraic rows;
       * MICM uses the step change there instead. */
      if (s->n_algebraic)
        tmc_set_algebraic_errors(s, Yerror, Y, Ynew);

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
