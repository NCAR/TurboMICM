/* Block-diagonal sparse matrices with one sparsity pattern per block, in CSR
 * or CSC element numbering and vector ordering of width L. */
#include "internal.h"

struct tmc_sparse_builder
{
  int n;
  int n_blocks;
  double initial_value;
  unsigned char *present; /* n x n */
};

tmc_sparse_builder *tmc_sparse_builder_create(int n)
{
  tmc_sparse_builder *b = tmc_calloc(1, sizeof *b);
  b->n = n;
  b->n_blocks = 1;
  b->initial_value = 0.0;
  b->present = tmc_calloc((size_t)n * (size_t)n, 1);
  return b;
}

void tmc_sparse_builder_free(tmc_sparse_builder *b)
{
  if (!b)
    return;
  free(b->present);
  free(b);
}

int tmc_sparse_builder_with_element(tmc_sparse_builder *b, int row, int col, tmc_error *err)
{
  if (row < 0 || col < 0 || row >= b->n || col >= b->n)
    return tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Element (%d, %d) is out of range for a %d x %d block", row, col, b->n, b->n);
  b->present[(size_t)row * (size_t)b->n + (size_t)col] = 1;
  return TMC_OK;
}

void tmc_sparse_builder_set_number_of_blocks(tmc_sparse_builder *b, int n_blocks)
{
  b->n_blocks = n_blocks;
}

void tmc_sparse_builder_set_initial_value(tmc_sparse_builder *b, double value)
{
  b->initial_value = value;
}

int tmc_sparse_builder_number_of_elements(const tmc_sparse_builder *b)
{
  int count = 0;
  for (size_t i = 0; i < (size_t)b->n * (size_t)b->n; ++i)
    count += b->present[i];
  return count;
}

static tmc_sparse *from_pattern(int n, const unsigned char *present, tmc_sparse_ordering ordering, int L, int n_blocks, double initial_value, bool indexing_only)
{
  if (L < 1)
    L = 1;
  tmc_sparse *m = tmc_calloc(1, sizeof *m);
  m->n = n;
  m->n_blocks = n_blocks;
  m->L = L;
  m->n_groups = (n_blocks + L - 1) / L;
  m->ordering = ordering;
  int nnz = 0;
  for (size_t i = 0; i < (size_t)n * (size_t)n; ++i)
    nnz += present[i] ? 1 : 0;
  m->nnz = nnz;
  m->row_of = tmc_malloc((size_t)nnz * sizeof(int));
  m->col_of = tmc_malloc((size_t)nnz * sizeof(int));
  m->lookup = tmc_malloc((size_t)n * (size_t)n * sizeof(int));
  for (size_t i = 0; i < (size_t)n * (size_t)n; ++i)
    m->lookup[i] = -1;
  /* number the elements in the chosen order */
  int e = 0;
  if (ordering == TMC_SPARSE_CSR)
  {
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        if (present[(size_t)i * n + j])
        {
          m->row_of[e] = i;
          m->col_of[e] = j;
          m->lookup[(size_t)i * n + j] = e;
          ++e;
        }
  }
  else
  {
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        if (present[(size_t)i * n + j])
        {
          m->row_of[e] = i;
          m->col_of[e] = j;
          m->lookup[(size_t)i * n + j] = e;
          ++e;
        }
  }
  /* row and column indexes over the element ids */
  m->row_start = tmc_calloc((size_t)n + 1, sizeof(int));
  m->col_start = tmc_calloc((size_t)n + 1, sizeof(int));
  m->row_cols = tmc_malloc((size_t)nnz * sizeof(int));
  m->row_elem = tmc_malloc((size_t)nnz * sizeof(int));
  m->col_rows = tmc_malloc((size_t)nnz * sizeof(int));
  m->col_elem = tmc_malloc((size_t)nnz * sizeof(int));
  m->diag = tmc_malloc((size_t)n * sizeof(int));
  int k = 0;
  for (int i = 0; i < n; ++i)
  {
    m->row_start[i] = k;
    for (int j = 0; j < n; ++j)
    {
      int id = m->lookup[(size_t)i * n + j];
      if (id >= 0)
      {
        m->row_cols[k] = j;
        m->row_elem[k] = id;
        ++k;
      }
    }
    m->diag[i] = m->lookup[(size_t)i * n + i];
  }
  m->row_start[n] = k;
  k = 0;
  for (int j = 0; j < n; ++j)
  {
    m->col_start[j] = k;
    for (int i = 0; i < n; ++i)
    {
      int id = m->lookup[(size_t)i * n + j];
      if (id >= 0)
      {
        m->col_rows[k] = i;
        m->col_elem[k] = id;
        ++k;
      }
    }
  }
  m->col_start[n] = k;
  m->size = (size_t)m->n_groups * (size_t)L * (size_t)nnz;
  if (!indexing_only)
  {
    m->data = tmc_malloc(m->size * sizeof(double));
    for (size_t i = 0; i < m->size; ++i)
      m->data[i] = initial_value;
  }
  return m;
}

tmc_sparse *tmc_sparse_build(const tmc_sparse_builder *b, tmc_sparse_ordering ordering, int L, bool indexing_only)
{
  return from_pattern(b->n, b->present, ordering, L, b->n_blocks, b->initial_value, indexing_only);
}

tmc_sparse *tmc_sparse_create_like(const tmc_sparse *pattern, int n_blocks, double initial_value)
{
  unsigned char *present = tmc_calloc((size_t)pattern->n * (size_t)pattern->n, 1);
  for (int e = 0; e < pattern->nnz; ++e)
    present[(size_t)pattern->row_of[e] * pattern->n + pattern->col_of[e]] = 1;
  tmc_sparse *m = from_pattern(pattern->n, present, pattern->ordering, pattern->L, n_blocks, initial_value, false);
  free(present);
  return m;
}

tmc_sparse *tmc_sparse_copy(const tmc_sparse *src)
{
  tmc_sparse *m = tmc_sparse_create_like(src, src->n_blocks, 0.0);
  if (src->data)
    memcpy(m->data, src->data, m->size * sizeof(double));
  else
  {
    free(m->data);
    m->data = NULL;
  }
  return m;
}

void tmc_sparse_free(tmc_sparse *m)
{
  if (!m)
    return;
  free(m->row_of);
  free(m->col_of);
  free(m->row_start);
  free(m->row_cols);
  free(m->row_elem);
  free(m->col_start);
  free(m->col_rows);
  free(m->col_elem);
  free(m->diag);
  free(m->lookup);
  free(m->data);
  free(m);
}

bool tmc_sparse_is_zero(const tmc_sparse *m, int row, int col)
{
  if (row < 0 || col < 0 || row >= m->n || col >= m->n)
    return true;
  return m->lookup[(size_t)row * m->n + col] < 0;
}

int tmc_sparse_element(const tmc_sparse *m, int row, int col)
{
  if (row < 0 || col < 0 || row >= m->n || col >= m->n)
    return -1;
  return m->lookup[(size_t)row * m->n + col];
}

size_t tmc_sparse_data_index(const tmc_sparse *m, int block, int element)
{
  return ((size_t)(block / m->L) * (size_t)m->nnz + (size_t)element) * (size_t)m->L + (size_t)(block % m->L);
}

size_t tmc_sparse_index(const tmc_sparse *m, int block, int row, int col, tmc_error *err)
{
  if (row < 0 || col < 0 || row >= m->n || col >= m->n)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Element (%d, %d) is out of range", row, col);
    return (size_t)-1;
  }
  int e = m->lookup[(size_t)row * m->n + col];
  if (e < 0)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ZERO_ELEMENT_ACCESS, "Element (%d, %d) is zero", row, col);
    return (size_t)-1;
  }
  if (block < 0 || block >= m->n_blocks)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Block %d is out of range", block);
    return (size_t)-1;
  }
  return tmc_sparse_data_index(m, block, e);
}

double tmc_sparse_get(const tmc_sparse *m, int block, int row, int col)
{
  int e = m->lookup[(size_t)row * m->n + col];
  return e < 0 ? 0.0 : m->data[tmc_sparse_data_index(m, block, e)];
}

void tmc_sparse_set(tmc_sparse *m, int block, int row, int col, double value)
{
  int e = m->lookup[(size_t)row * m->n + col];
  if (e >= 0)
    m->data[tmc_sparse_data_index(m, block, e)] = value;
}

void tmc_sparse_fill(tmc_sparse *m, double value)
{
  for (size_t i = 0; i < m->size; ++i)
    m->data[i] = value;
}

void tmc_sparse_add_to_diagonal(tmc_sparse *m, double value)
{
  for (int b = 0; b < m->n_blocks; ++b)
    for (int i = 0; i < m->n; ++i)
      if (m->diag[i] >= 0)
        m->data[tmc_sparse_data_index(m, b, m->diag[i])] += value;
}

size_t tmc_sparse_flat_block_size(const tmc_sparse *m)
{
  return (size_t)m->nnz;
}
