/* Solver internals and utilities MICM exposes to its tests (tmc_solver_ext.h). */
#include "internal.h"
#include "tmc_solver_ext.h"

#include <math.h>

/* ---- Rosenbrock and backward Euler internals ------------------------------ */

static bool shaped_like_variables(const tmc_state *st, const tmc_dense *m)
{
  return m->n_rows == st->n_cells && m->n_cols == st->n_vars && m->L == st->L;
}

double tmc_solver_normalized_error(const tmc_solver *s, const tmc_state *st, const tmc_dense *y, const tmc_dense *ynew, const tmc_dense *errors)
{
  (void)s;
  if (!shaped_like_variables(st, y) || !shaped_like_variables(st, ynew) || !shaped_like_variables(st, errors))
    return NAN;
  /* The accumulation of src/rosenbrock.c: variable by variable, block by
   * block, lane by lane, padded lanes excluded. */
  const int L = st->L, n_vars = st->n_vars, n_cells = st->n_cells, n_blocks = st->n_blocks;
  const double rtol = st->relative_tolerance;
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

void tmc_solver_alpha_minus_jacobian(const tmc_solver *s, tmc_state *st, double alpha)
{
  /* As the generic path of src/rosenbrock.c shifts the diagonal before factoring. */
  tmc_sparse *J = st->jacobian;
  const int L = J->L, nnz = J->nnz;
  for (int g = 0; g < J->n_groups; ++g)
  {
    double *block = J->data + (size_t)g * (size_t)nnz * L;
    for (int i = 0; i < s->n_vars; ++i)
    {
      double *d = block + (size_t)s->lu->mz_diag[i] * L;
      for (int c = 0; c < L; ++c)
        d[c] += alpha;
    }
  }
}

bool tmc_backward_euler_is_converged(const tmc_backward_euler_parameters *p, const tmc_dense *residual, const tmc_dense *yn1, const double *absolute_tolerance, double relative_tolerance)
{
  /* The test of src/backward_euler.c with the tolerances passed in. */
  const int L = residual->L, n_vars = residual->n_cols, n_cells = residual->n_rows, n_blocks = residual->n_groups;
  for (int v = 0; v < n_vars; ++v)
  {
    const double atol = absolute_tolerance[v];
    for (int block = 0; block < n_blocks; ++block)
    {
      const size_t base = ((size_t)block * n_vars + v) * L;
      const int lanes = block == n_blocks - 1 ? n_cells - block * L : L;
      for (int lane = 0; lane < lanes; ++lane)
      {
        double r = residual->data[base + lane];
        if (!isfinite(r))
          return false;
        double ar = fabs(r);
        if (ar > p->small && ar > atol && ar > relative_tolerance * fabs(yn1->data[tmc_dense_index(yn1, block * L + lane, v)]))
          return false;
      }
    }
  }
  return true;
}

int tmc_state_unsafely_set_custom_rate_parameters(tmc_state *st, const double *values, int n_cells, int n_parameters, tmc_error *err)
{
  if (n_cells != st->n_cells)
    return tmc_error_set(err, TMC_ERROR_STATE, TMC_STATE_INVALID_CUSTOM_RATE_PARAM_COUNT_MULTIGRID, "Custom rate parameter vector size does not match the number of grid cells");
  if (n_parameters != st->n_custom)
    return tmc_error_set(err, TMC_ERROR_STATE, TMC_STATE_INVALID_CUSTOM_RATE_PARAM_COUNT, "Custom rate parameter count does not match the number of rate constants");
  for (int c = 0; c < n_cells; ++c)
    tmc_dense_set_row(st->custom_rate_parameters, c, values + (size_t)c * (size_t)n_parameters);
  return TMC_OK;
}

/* ---- Jacobian verification ------------------------------------------------- */

tmc_dense *tmc_finite_difference_jacobian(tmc_forcing_function forcing, void *context, const tmc_dense *variables, int n_species, double perturbation)
{
  const int n_blocks = variables->n_rows, n_cols = variables->n_cols, L = variables->L;
  tmc_dense *result = tmc_dense_create(n_blocks, n_species * n_species, L, 0.0);
  tmc_dense *plus = tmc_dense_create(n_blocks, n_cols, L, 0.0);
  tmc_dense *minus = tmc_dense_create(n_blocks, n_cols, L, 0.0);
  tmc_dense *f_plus = tmc_dense_create(n_blocks, n_cols, L, 0.0);
  tmc_dense *f_minus = tmc_dense_create(n_blocks, n_cols, L, 0.0);
  tmc_dense *f_base = NULL; /* F(x), computed when a forward difference needs it */
  for (int col = 0; col < n_species; ++col)
  {
    tmc_dense_copy_into(plus, variables);
    tmc_dense_copy_into(minus, variables);
    for (int block = 0; block < n_blocks; ++block)
    {
      double x = tmc_dense_get(variables, block, col);
      double h = perturbation * fmax(1.0, fabs(x));
      tmc_dense_set(plus, block, col, x + h);
      tmc_dense_set(minus, block, col, x - h);
    }
    tmc_dense_fill(f_plus, 0.0);
    tmc_dense_fill(f_minus, 0.0);
    forcing(plus, f_plus, context);
    forcing(minus, f_minus, context);
    for (int block = 0; block < n_blocks; ++block)
    {
      double x = tmc_dense_get(variables, block, col);
      double h = perturbation * fmax(1.0, fabs(x));
      if (x - h >= 0.0)
      {
        double two_h = 2.0 * h;
        for (int row = 0; row < n_species; ++row)
          tmc_dense_set(result, block, row * n_species + col, (tmc_dense_get(f_plus, block, row) - tmc_dense_get(f_minus, block, row)) / two_h);
      }
      else
      {
        if (!f_base)
        {
          f_base = tmc_dense_create(n_blocks, n_cols, L, 0.0);
          forcing(variables, f_base, context);
        }
        for (int row = 0; row < n_species; ++row)
          tmc_dense_set(result, block, row * n_species + col, (tmc_dense_get(f_plus, block, row) - tmc_dense_get(f_base, block, row)) / h);
      }
    }
  }
  tmc_dense_free(plus);
  tmc_dense_free(minus);
  tmc_dense_free(f_plus);
  tmc_dense_free(f_minus);
  tmc_dense_free(f_base);
  return result;
}

static tmc_jacobian_comparison comparison_init(void)
{
  tmc_jacobian_comparison r;
  memset(&r, 0, sizeof r);
  r.passed = true;
  return r;
}

tmc_jacobian_comparison tmc_compare_jacobian_to_finite_difference(const tmc_sparse *analytical, const tmc_dense *fd, int n_species, double atol, double rtol)
{
  tmc_jacobian_comparison result = comparison_init();
  for (int block = 0; block < analytical->n_blocks; ++block)
    for (int row = 0; row < n_species; ++row)
      for (int col = 0; col < n_species; ++col)
      {
        if (tmc_sparse_is_zero(analytical, row, col))
          continue;
        /* analytical holds -(df/dx) */
        double analytical_value = -tmc_sparse_get(analytical, block, row, col);
        double fd_value = tmc_dense_get(fd, block, row * n_species + col);
        double abs_error = fabs(analytical_value - fd_value);
        double scale = fmax(fabs(analytical_value), fabs(fd_value));
        double rel_error = scale > 0.0 ? abs_error / scale : 0.0;
        if (abs_error > result.max_abs_error)
        {
          result.max_abs_error = abs_error;
          result.worst_block = block;
          result.worst_row = row;
          result.worst_col = col;
          result.worst_analytical = analytical_value;
          result.worst_fd = fd_value;
        }
        if (rel_error > result.max_rel_error)
          result.max_rel_error = rel_error;
        if (abs_error > atol + rtol * scale)
          result.passed = false;
      }
  return result;
}

tmc_jacobian_comparison tmc_check_jacobian_sparsity_completeness(const tmc_sparse *analytical, const tmc_dense *fd, int n_species, double threshold)
{
  tmc_jacobian_comparison result = comparison_init();
  for (int block = 0; block < analytical->n_blocks; ++block)
    for (int row = 0; row < n_species; ++row)
      for (int col = 0; col < n_species; ++col)
      {
        if (!tmc_sparse_is_zero(analytical, row, col))
          continue;
        double fd_value = fabs(tmc_dense_get(fd, block, row * n_species + col));
        if (fd_value > threshold)
        {
          result.passed = false;
          if (fd_value > result.max_abs_error)
          {
            result.max_abs_error = fd_value;
            result.worst_block = block;
            result.worst_row = row;
            result.worst_col = col;
            result.worst_analytical = 0.0;
            result.worst_fd = fd_value;
          }
        }
      }
  return result;
}
