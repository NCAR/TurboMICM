/* MICM's sparse container interface beyond what the solver uses: VectorIndex
 * with MICM's checks, a checked IsZero, DiagonalIndices, the group geometry of
 * the vector orderings, and printing. See include/tmc_sparse_ext.h. */
#include "internal.h"
#include "tmc_sparse_ext.h"

#include <stdio.h>

/* ------------------------------------------------------------------------- */
/* Indices                                                                   */
/* ------------------------------------------------------------------------- */

static bool in_block(const tmc_sparse *m, int row, int column)
{
  return row >= 0 && column >= 0 && row < m->n && column < m->n;
}

size_t tmc_sparse_vector_index(const tmc_sparse *m, int block, int row, int column, tmc_error *err)
{
  if (!in_block(m, row, column) || block < 0 || block >= m->n_blocks)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ELEMENT_OUT_OF_RANGE,
                  "Element (%d, %d) of block %d is out of range for %d blocks of %d x %d", row, column, block, m->n_blocks, m->n, m->n);
    return (size_t)-1;
  }
  int e = tmc_sparse_element(m, row, column);
  if (e < 0)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ZERO_ELEMENT_ACCESS, "Element (%d, %d) is always zero", row, column);
    return (size_t)-1;
  }
  return tmc_sparse_data_index(m, block, e);
}

size_t tmc_sparse_vector_index_one_block(const tmc_sparse *m, int row, int column, tmc_error *err)
{
  if (m->n_blocks != 1)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_MISSING_BLOCK_INDEX,
                  "Matrix has %d blocks; use the (block, row, column) form of the vector index instead", m->n_blocks);
    return (size_t)-1;
  }
  return tmc_sparse_vector_index(m, 0, row, column, err);
}

bool tmc_sparse_is_zero_checked(const tmc_sparse *m, int row, int column, tmc_error *err)
{
  if (!in_block(m, row, column))
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Element (%d, %d) is out of range for a %d x %d block", row, column, m->n, m->n);
    return true;
  }
  return tmc_sparse_is_zero(m, row, column);
}

int tmc_sparse_diagonal_indices(const tmc_sparse *m, int block, size_t *indices, tmc_error *err)
{
  if (block < 0 || block >= m->n_blocks)
  {
    tmc_error_set(err, TMC_ERROR_MATRIX, TMC_MATRIX_ELEMENT_OUT_OF_RANGE, "Block %d is out of range for %d blocks", block, m->n_blocks);
    return -1;
  }
  int count = 0;
  for (int i = 0; i < m->n; ++i)
    if (m->diag[i] >= 0)
    {
      if (indices)
        indices[count] = tmc_sparse_data_index(m, block, m->diag[i]);
      ++count;
    }
  return count;
}

/* ------------------------------------------------------------------------- */
/* Group geometry                                                            */
/* ------------------------------------------------------------------------- */

int tmc_sparse_group_width(const tmc_sparse *m)
{
  return m->L;
}

int tmc_sparse_group_size(const tmc_sparse *m)
{
  return m->L * m->nnz;
}

int tmc_sparse_number_of_groups(const tmc_sparse *m, int n_blocks)
{
  if (n_blocks <= 0)
    return 0;
  return (n_blocks + m->L - 1) / m->L;
}

/* ------------------------------------------------------------------------- */
/* Printing                                                                  */
/* ------------------------------------------------------------------------- */

/* Where the text goes: a buffer of `size` characters (NUL included) that
 * keeps the start of a text too long for it, or a stream. `length` counts the
 * whole text either way. */
typedef struct sink
{
  char *buffer;
  size_t size;
  FILE *stream;
  size_t length;
  bool failed;
} sink;

static void put(sink *s, const char *text)
{
  for (const char *c = text; *c; ++c, ++s->length)
  {
    if (s->stream)
    {
      if (!s->failed && fputc(*c, s->stream) == EOF)
        s->failed = true;
    }
    else if (s->length + 1 < s->size)
      s->buffer[s->length] = *c;
  }
}

static void put_value(sink *s, double value)
{
  char text[64];
  snprintf(text, sizeof text, "%g", value);
  put(s, text);
}

static void put_int(sink *s, int value)
{
  char text[32];
  snprintf(text, sizeof text, "%d", value);
  put(s, text);
}

static void put_block_header(sink *s, int block)
{
  put(s, "Block ");
  put_int(s, block);
  put(s, "\n");
}

/* MICM's operator<<: every element of every block, zeros included. */
static void write_matrix(const tmc_sparse *m, sink *s)
{
  for (int b = 0; b < m->n_blocks; ++b)
  {
    put_block_header(s, b);
    for (int i = 0; i < m->n; ++i)
      for (int j = 0; j < m->n; ++j)
      {
        if (tmc_sparse_is_zero(m, i, j))
          put(s, "0");
        else
          put_value(s, tmc_sparse_get(m, b, i, j));
        put(s, j + 1 < m->n ? "," : "\n");
      }
  }
}

/* MICM's PrintNonZeroElements: the non-zero elements, rows then columns. */
static void write_non_zero(const tmc_sparse *m, sink *s)
{
  for (int b = 0; b < m->n_blocks; ++b)
  {
    put_block_header(s, b);
    for (int i = 0; i < m->n; ++i)
      for (int j = 0; j < m->n; ++j)
        if (!tmc_sparse_is_zero(m, i, j))
        {
          put_int(s, i);
          put(s, ", ");
          put_int(s, j);
          put(s, ", ");
          put_value(s, tmc_sparse_get(m, b, i, j));
          put(s, "\n");
        }
  }
}

static size_t format_with(void (*write)(const tmc_sparse *, sink *), const tmc_sparse *m, char *buffer, size_t size)
{
  sink s = { buffer, size, NULL, 0, false };
  write(m, &s);
  if (size > 0)
    buffer[s.length < size ? s.length : size - 1] = '\0';
  return s.length;
}

static int print_with(void (*write)(const tmc_sparse *, sink *), const tmc_sparse *m, FILE *stream)
{
  sink s = { NULL, 0, stream, 0, false };
  write(m, &s);
  if (s.failed)
    return -1;
  return (int)s.length;
}

size_t tmc_sparse_format(const tmc_sparse *m, char *buffer, size_t size)
{
  return format_with(write_matrix, m, buffer, size);
}

int tmc_sparse_print(const tmc_sparse *m, FILE *stream)
{
  return print_with(write_matrix, m, stream);
}

size_t tmc_sparse_format_non_zero(const tmc_sparse *m, char *buffer, size_t size)
{
  return format_with(write_non_zero, m, buffer, size);
}

int tmc_sparse_print_non_zero(const tmc_sparse *m, FILE *stream)
{
  return print_with(write_non_zero, m, stream);
}
