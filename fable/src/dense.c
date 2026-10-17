/* Dense matrices in vector ordering (MICM's Matrix is the L = 1 case of its
 * VectorMatrix<L>). */
#include "internal.h"

tmc_dense *tmc_dense_create(int n_rows, int n_cols, int L, double initial_value)
{
  if (L < 1)
    L = 1;
  tmc_dense *m = tmc_calloc(1, sizeof *m);
  m->n_rows = n_rows;
  m->n_cols = n_cols;
  m->L = L;
  m->n_groups = (n_rows + L - 1) / L;
  m->size = (size_t)m->n_groups * (size_t)L * (size_t)n_cols;
  m->data = tmc_alloc_doubles(m->size);
  for (size_t i = 0; i < m->size; ++i)
    m->data[i] = initial_value;
  return m;
}

tmc_dense *tmc_dense_copy(const tmc_dense *src)
{
  tmc_dense *m = tmc_dense_create(src->n_rows, src->n_cols, src->L, 0.0);
  memcpy(m->data, src->data, m->size * sizeof(double));
  return m;
}

void tmc_dense_free(tmc_dense *m)
{
  if (!m)
    return;
  tmc_free_doubles(m->data);
  free(m);
}

size_t tmc_dense_index(const tmc_dense *m, int row, int col)
{
  return ((size_t)(row / m->L) * (size_t)m->n_cols + (size_t)col) * (size_t)m->L + (size_t)(row % m->L);
}

double tmc_dense_get(const tmc_dense *m, int row, int col)
{
  return m->data[tmc_dense_index(m, row, col)];
}

void tmc_dense_set(tmc_dense *m, int row, int col, double value)
{
  m->data[tmc_dense_index(m, row, col)] = value;
}

void tmc_dense_fill(tmc_dense *m, double value)
{
  for (size_t i = 0; i < m->size; ++i)
    m->data[i] = value;
}

void tmc_dense_copy_into(tmc_dense *dst, const tmc_dense *src)
{
  memcpy(dst->data, src->data, dst->size * sizeof(double));
}

void tmc_dense_swap(tmc_dense *a, tmc_dense *b)
{
  double *t = a->data;
  a->data = b->data;
  b->data = t;
}

void tmc_dense_axpy(tmc_dense *y, double alpha, const tmc_dense *x)
{
  double *restrict yd = y->data;
  const double *restrict xd = x->data;
  size_t n = y->size;
#pragma omp simd
  for (size_t i = 0; i < n; ++i)
    yd[i] += alpha * xd[i];
}

/* y = std::max(y, x) and y = std::min(y, x) as MICM writes them: the
 * comparison is on x, so a NaN or a -0.0 in y stays (MICM's post-solve clamp
 * leaves a failed solve's NaN visible). */
void tmc_dense_max(tmc_dense *m, double x)
{
  for (size_t i = 0; i < m->size; ++i)
    m->data[i] = m->data[i] < x ? x : m->data[i];
}

void tmc_dense_min(tmc_dense *m, double x)
{
  for (size_t i = 0; i < m->size; ++i)
    m->data[i] = x < m->data[i] ? x : m->data[i];
}

size_t tmc_dense_row_stride(const tmc_dense *m)
{
  return m->L == 1 ? (size_t)m->n_cols : 1;
}

size_t tmc_dense_column_stride(const tmc_dense *m)
{
  return (size_t)m->L;
}

void tmc_dense_get_row(const tmc_dense *m, int row, double *out)
{
  for (int j = 0; j < m->n_cols; ++j)
    out[j] = m->data[tmc_dense_index(m, row, j)];
}

void tmc_dense_set_row(tmc_dense *m, int row, const double *values)
{
  for (int j = 0; j < m->n_cols; ++j)
    m->data[tmc_dense_index(m, row, j)] = values[j];
}
