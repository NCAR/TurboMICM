#include "internal.h"

/* Arithmetic within each cell preserves reaction/derivative order. Each vector
   loop has contiguous independent cells, with no indirect scatter in SIMD loops. */
void ck_add_forcing(const ck_state *s, const double *restrict y, double *restrict f)
{
  const ck_model *m = s->model;
  const size_t cells = s->cells;
  double *restrict rate = ((ck_workspace *)s->workspace)->scratch;
  for (size_t r = 0; r < m->nr; ++r)
  {
    memcpy(rate, s->rates + r * cells, cells * sizeof(double));
    for (size_t j = 0; j < m->nreactants[r]; ++j)
    {
      const double *restrict v = y + m->reactants[r][j] * cells;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        rate[c] *= v[c];
    }
    for (size_t j = 0; j < m->nreactants[r]; ++j)
    {
      if (m->algebraic[m->reactants[r][j]])
        continue;
      double *restrict v = f + m->reactants[r][j] * cells;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        v[c] -= rate[c];
    }
    for (size_t j = 0; j < m->nproducts[r]; ++j)
    {
      if (m->algebraic[m->products[r][j].species])
        continue;
      double *restrict v = f + m->products[r][j].species * cells;
      double yield = m->products[r][j].yield;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        v[c] += yield * rate[c];
    }
  }
}
void ck_forcing(const ck_state *s, const double *restrict y, double *restrict f)
{
  const ck_model *m = s->model;
  const size_t cells = s->cells;
  memset(f, 0, m->n * cells * sizeof(double));
  ck_add_forcing(s, y, f);
  for (size_t k = 0; k < m->nexternal; ++k)
    if (m->external[k].forcing)
      m->external[k].forcing(s, y, f, m->external[k].context);
  if (m->nexternal && m->nconstraints)
    for (size_t j = 0; j < m->n; ++j)
      if (m->algebraic[j])
        memset(f + j * cells, 0, cells * sizeof(double));
}
void ck_subtract_jacobian(const ck_state *s, const double *restrict y, double *restrict a)
{
  const ck_model *m = s->model;
  const size_t cells = s->cells;
  double *restrict rate = ((ck_workspace *)s->workspace)->scratch;
  for (size_t job = 0; job < m->njobs; ++job)
  {
    const ck_derivative *d = &m->jobs[job];
    memcpy(rate, s->rates + d->reaction * cells, cells * sizeof(double));
    for (size_t j = 0; j < d->nreactants; ++j)
    {
      const double *restrict v = y + d->reactants[j] * cells;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        rate[c] *= v[c];
    }
    for (size_t j = 0; j <= d->nreactants; ++j)
    {
      size_t row = j < d->nreactants ? d->reactants[j] : d->independent;
      if (m->algebraic[row])
        continue;
      double *restrict v = a + d->targets[j] * cells;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        v[c] += rate[c];
    }
    for (size_t j = 0; j < d->nproducts; ++j)
    {
      if (m->algebraic[d->products[j].species])
        continue;
      double *restrict v = a + d->targets[d->nreactants + 1 + j] * cells;
      double yield = d->products[j].yield;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        v[c] -= yield * rate[c];
    }
  }
}
void ck_jacobian(const ck_state *s, const double *restrict y, double *restrict a)
{
  const ck_model *m = s->model;
  const size_t cells = s->cells;
  memset(a, 0, m->nnz * cells * sizeof(double));
  ck_subtract_jacobian(s, y, a);
  for (size_t k = 0; k < m->nexternal; ++k)
    if (m->external[k].jacobian_terms)
      m->external[k].jacobian_terms(s, y, a, m->external[k].context);
  if (m->nexternal && m->nconstraints)
    for (size_t j = 0; j < m->n; ++j)
      if (m->algebraic[j])
        for (size_t p = m->row_start[j]; p < m->row_start[j + 1]; ++p)
          memset(a + p * cells, 0, cells * sizeof(double));
}
void ck_factor(const ck_model *m, size_t cells, double *restrict a, ck_lu_algorithm algorithm)
{
  if (algorithm == CK_LU_MOZART)
  {
    for (size_t i = 0; i < m->n; ++i)
    {
      const double *restrict diagonal = a + m->diagonal[i] * cells;
      /* Divide all below-diagonal values by the same reciprocal, as Mozart. */
      for (size_t j = m->divide_start[i]; j < m->divide_start[i + 1]; ++j)
      {
        double *restrict v = a + m->divides[j] * cells;
#pragma omp simd
        for (size_t c = 0; c < cells; ++c)
          v[c] *= 1.0 / diagonal[c];
      }
      for (size_t j = m->update_start[i]; j < m->update_start[i + 1]; ++j)
      {
        ck_operation op = m->ops[j];
        double *restrict target = a + op.target * cells;
        const double *restrict left = a + op.left * cells, *restrict right = a + op.right * cells;
#pragma omp simd
        for (size_t c = 0; c < cells; ++c)
          target[c] -= left[c] * right[c];
      }
    }
  }
  else
  {
    /* Doolittle updates each upper row, then each lower column. */
    const size_t n = m->n;
    for (size_t i = 0; i < n; ++i)
    {
      for (size_t p = m->row_start[i]; p < m->row_start[i + 1]; ++p)
      {
        size_t k = m->columns[p];
        if (k < i)
          continue;
        double *restrict v = a + p * cells;
        for (size_t j = 0; j < i; ++j)
        {
          size_t l = m->map[i * n + j], u = m->map[j * n + k];
          if (l == SIZE_MAX || u == SIZE_MAX)
            continue;
          const double *restrict left = a + l * cells, *restrict right = a + u * cells;
#pragma omp simd
          for (size_t c = 0; c < cells; ++c)
            v[c] -= left[c] * right[c];
        }
      }
      for (size_t k = i + 1; k < n; ++k)
      {
        size_t off = m->map[k * n + i];
        if (off == SIZE_MAX)
          continue;
        double *restrict v = a + off * cells;
        for (size_t j = 0; j < i; ++j)
        {
          size_t l = m->map[k * n + j], u = m->map[j * n + i];
          if (l == SIZE_MAX || u == SIZE_MAX)
            continue;
          const double *restrict left = a + l * cells, *restrict right = a + u * cells;
#pragma omp simd
          for (size_t c = 0; c < cells; ++c)
            v[c] -= left[c] * right[c];
        }
        const double *restrict diag = a + m->diagonal[i] * cells;
#pragma omp simd
        for (size_t c = 0; c < cells; ++c)
          v[c] /= diag[c];
      }
    }
  }
}
void ck_linear_solve(const ck_model *m, size_t cells, const double *restrict a, double *restrict x)
{
  for (size_t i = 0; i < m->n; ++i)
  {
    double *restrict xi = x + i * cells;
    for (size_t p = m->row_start[i]; p < m->diagonal[i]; ++p)
    {
      const double *restrict v = a + p * cells, *restrict xj = x + m->columns[p] * cells;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        xi[c] -= v[c] * xj[c];
    }
  }
  for (size_t i = m->n; i-- > 0;)
  {
    double *restrict xi = x + i * cells;
    for (size_t p = m->diagonal[i] + 1; p < m->row_start[i + 1]; ++p)
    {
      const double *restrict v = a + p * cells, *restrict xj = x + m->columns[p] * cells;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        xi[c] -= v[c] * xj[c];
    }
    const double *restrict diag = a + m->diagonal[i] * cells;
#pragma omp simd
    for (size_t c = 0; c < cells; ++c)
      xi[c] /= diag[c];
  }
}

ck_status ck_factor_separate(const ck_model *m, size_t cells, const double *a, double *l, double *u,
                             ck_lu_algorithm algorithm)
{
  if (!m || !cells || !a || !l || !u)
    return CK_INVALID;
  double *packed = ck_allocate(m->nnz * cells, sizeof(double));
  if (!packed)
    return CK_MEMORY;
  memcpy(packed, a, m->nnz * cells * sizeof(double));
  ck_factor(m, cells, packed, algorithm);
  memset(l, 0, m->nnz * cells * sizeof(double));
  memset(u, 0, m->nnz * cells * sizeof(double));
  for (size_t i = 0; i < m->n; ++i)
  {
    for (size_t p = m->row_start[i]; p < m->row_start[i + 1]; ++p)
    {
      if (m->columns[p] < i)
        memcpy(l + p * cells, packed + p * cells, cells * sizeof(double));
      else
        memcpy(u + p * cells, packed + p * cells, cells * sizeof(double));
    }
    for (size_t c = 0; c < cells; ++c)
      l[m->diagonal[i] * cells + c] = 1;
  }
  free(packed);
  return CK_OK;
}
void ck_linear_solve_separate(const ck_model *m, size_t cells, const double *restrict l,
                              const double *restrict u, double *restrict x)
{
  for (size_t i = 0; i < m->n; ++i)
  {
    double *restrict xi = x + i * cells;
    for (size_t p = m->row_start[i]; p < m->diagonal[i]; ++p)
    {
      const double *restrict v = l + p * cells, *restrict xj = x + m->columns[p] * cells;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        xi[c] -= v[c] * xj[c];
    }
    const double *restrict diagonal = l + m->diagonal[i] * cells;
#pragma omp simd
    for (size_t c = 0; c < cells; ++c)
      xi[c] /= diagonal[c];
  }
  for (size_t i = m->n; i-- > 0;)
  {
    double *restrict xi = x + i * cells;
    for (size_t p = m->diagonal[i] + 1; p < m->row_start[i + 1]; ++p)
    {
      const double *restrict v = u + p * cells, *restrict xj = x + m->columns[p] * cells;
#pragma omp simd
      for (size_t c = 0; c < cells; ++c)
        xi[c] -= v[c] * xj[c];
    }
    const double *restrict diagonal = u + m->diagonal[i] * cells;
#pragma omp simd
    for (size_t c = 0; c < cells; ++c)
      xi[c] /= diagonal[c];
  }
}

void ck_alpha_minus_jacobian(const ck_model *m, size_t cells, double *restrict a, double alpha)
{
  for (size_t i = 0; i < m->n; ++i)
  {
    double scaled_alpha = alpha * !ck_is_algebraic(m, i);
    double *restrict diagonal = a + m->diagonal[i] * cells;
#pragma omp simd
    for (size_t c = 0; c < cells; ++c)
      diagonal[c] += scaled_alpha;
  }
}

void ck_dispatch_factor(const ck_linear_operations *ops, const ck_model *m, size_t cells, double *matrix,
                        ck_lu_algorithm algorithm)
{
  if (ops && ops->factor)
    ops->factor(m, cells, matrix, algorithm, ops->context);
  else
    ck_factor(m, cells, matrix, algorithm);
}
void ck_dispatch_solve(const ck_linear_operations *ops, const ck_model *m, size_t cells, const double *lu,
                       double *rhs)
{
  if (ops && ops->solve)
    ops->solve(m, cells, lu, rhs, ops->context);
  else
    ck_linear_solve(m, cells, lu, rhs);
}

void ck_matrix_vector(const ck_model *m, size_t cells, const double *matrix, const double *x, double *result)
{
  memset(result, 0, m->n * cells * sizeof(double));
  for (size_t row = 0; row < m->n; ++row)
    for (size_t p = m->row_start[row]; p < m->row_start[row + 1]; ++p)
      for (size_t c = 0; c < cells; ++c)
        result[row * cells + c] += matrix[p * cells + c] * x[m->columns[p] * cells + c];
}
