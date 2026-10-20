#include "matrix.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int matrix_valid(const double *data, size_t rows, size_t columns,
                        size_t row_stride, size_t column_stride)
{
  return (rows == 0 || columns == 0 || data != NULL) &&
         (rows <= 1 || row_stride != 0) &&
         (columns <= 1 || column_stride != 0);
}

static int same_shape(size_t rows_a, size_t columns_a, size_t rows_b,
                      size_t columns_b)
{
  return rows_a == rows_b && columns_a == columns_b;
}

int chem_matrix_view_init(chem_matrix_view *view, double *data, size_t rows,
                          size_t columns, size_t row_stride,
                          size_t column_stride)
{
  if (!view || !matrix_valid(data, rows, columns, row_stride, column_stride))
    return -1;
  *view = (chem_matrix_view){data, rows, columns, row_stride, column_stride};
  return 0;
}

int chem_const_matrix_view_init(chem_const_matrix_view *view, const double *data,
                                size_t rows, size_t columns,
                                size_t row_stride, size_t column_stride)
{
  if (!view || !matrix_valid(data, rows, columns, row_stride, column_stride))
    return -1;
  *view = (chem_const_matrix_view){data, rows, columns, row_stride, column_stride};
  return 0;
}

double *chem_matrix_at(chem_matrix_view *view, size_t row, size_t column)
{
  if (!view || !view->data || row >= view->rows || column >= view->columns)
    return NULL;
  return view->data + row * view->row_stride + column * view->column_stride;
}

const double *chem_const_matrix_at(const chem_const_matrix_view *view,
                                   size_t row, size_t column)
{
  if (!view || !view->data || row >= view->rows || column >= view->columns)
    return NULL;
  return view->data + row * view->row_stride + column * view->column_stride;
}

int chem_matrix_fill(chem_matrix_view view, double value)
{
  size_t row, column;
  if (!matrix_valid(view.data, view.rows, view.columns, view.row_stride,
                    view.column_stride))
    return -1;
  for (row = 0; row < view.rows; ++row)
    for (column = 0; column < view.columns; ++column)
      view.data[row * view.row_stride + column * view.column_stride] = value;
  return 0;
}

int chem_matrix_copy(chem_matrix_view destination,
                     chem_const_matrix_view source)
{
  size_t row, column;
  if (!matrix_valid(destination.data, destination.rows, destination.columns,
                    destination.row_stride, destination.column_stride) ||
      !matrix_valid(source.data, source.rows, source.columns, source.row_stride,
                    source.column_stride) ||
      !same_shape(destination.rows, destination.columns, source.rows,
                  source.columns))
    return -1;
  for (row = 0; row < source.rows; ++row)
    for (column = 0; column < source.columns; ++column)
      destination.data[row * destination.row_stride + column * destination.column_stride] =
          source.data[row * source.row_stride + column * source.column_stride];
  return 0;
}

int chem_matrix_axpy(chem_matrix_view destination, double alpha,
                     chem_const_matrix_view source)
{
  size_t row, column;
  if (!matrix_valid(destination.data, destination.rows, destination.columns,
                    destination.row_stride, destination.column_stride) ||
      !matrix_valid(source.data, source.rows, source.columns, source.row_stride,
                    source.column_stride) ||
      !same_shape(destination.rows, destination.columns, source.rows,
                  source.columns))
    return -1;
  for (row = 0; row < source.rows; ++row)
    for (column = 0; column < source.columns; ++column)
      destination.data[row * destination.row_stride + column * destination.column_stride] +=
          alpha * source.data[row * source.row_stride + column * source.column_stride];
  return 0;
}

int chem_matrix_transform(chem_matrix_view destination,
                          chem_const_matrix_view source,
                          chem_unary_function function, void *context)
{
  size_t row, column;
  if (!function || !matrix_valid(destination.data, destination.rows,
                                 destination.columns, destination.row_stride,
                                 destination.column_stride) ||
      !matrix_valid(source.data, source.rows, source.columns, source.row_stride,
                    source.column_stride) ||
      !same_shape(destination.rows, destination.columns, source.rows,
                  source.columns))
    return -1;
  for (row = 0; row < source.rows; ++row)
    for (column = 0; column < source.columns; ++column)
      destination.data[row * destination.row_stride + column * destination.column_stride] =
          function(source.data[row * source.row_stride + column * source.column_stride], context);
  return 0;
}

int chem_matrix_transform2(chem_matrix_view destination,
                           chem_const_matrix_view left,
                           chem_const_matrix_view right,
                           chem_binary_function function, void *context)
{
  size_t row, column;
  if (!function || !matrix_valid(destination.data, destination.rows,
                                 destination.columns, destination.row_stride,
                                 destination.column_stride) ||
      !matrix_valid(left.data, left.rows, left.columns, left.row_stride,
                    left.column_stride) ||
      !matrix_valid(right.data, right.rows, right.columns, right.row_stride,
                    right.column_stride) ||
      !same_shape(destination.rows, destination.columns, left.rows, left.columns) ||
      !same_shape(destination.rows, destination.columns, right.rows, right.columns))
    return -1;
  for (row = 0; row < left.rows; ++row)
    for (column = 0; column < left.columns; ++column)
      destination.data[row * destination.row_stride + column * destination.column_stride] =
          function(left.data[row * left.row_stride + column * left.column_stride],
                   right.data[row * right.row_stride + column * right.column_stride], context);
  return 0;
}

static int matrix_reduce(chem_const_matrix_view view, double *value, int kind)
{
  size_t row, column, count = 0;
  double result = 0.0;
  if (!value || !matrix_valid(view.data, view.rows, view.columns, view.row_stride,
                              view.column_stride) || !view.rows || !view.columns)
    return -1;
  for (row = 0; row < view.rows; ++row) {
    for (column = 0; column < view.columns; ++column) {
      double x = view.data[row * view.row_stride + column * view.column_stride];
      if (kind == 0) {
        result += x;
      } else if (!count || (kind == 1 && x < result) || (kind == 2 && x > result)) {
        result = x;
      }
      ++count;
    }
  }
  *value = result;
  return 0;
}

int chem_matrix_sum(chem_const_matrix_view view, double *value)
{
  return matrix_reduce(view, value, 0);
}

int chem_matrix_min(chem_const_matrix_view view, double *value)
{
  return matrix_reduce(view, value, 1);
}

int chem_matrix_max(chem_const_matrix_view view, double *value)
{
  return matrix_reduce(view, value, 2);
}

int chem_vector_view_init(chem_vector_view *view, double *data, size_t length,
                          size_t stride)
{
  if (!view || (length && (!data || !stride)))
    return -1;
  *view = (chem_vector_view){data, length, stride};
  return 0;
}

int chem_const_vector_view_init(chem_const_vector_view *view,
                                const double *data, size_t length,
                                size_t stride)
{
  if (!view || (length && (!data || !stride)))
    return -1;
  *view = (chem_const_vector_view){data, length, stride};
  return 0;
}

int chem_vector_copy(chem_vector_view destination,
                     chem_const_vector_view source)
{
  size_t i;
  if (destination.length != source.length ||
      (destination.length && (!destination.data || !source.data ||
                              !destination.stride || !source.stride)))
    return -1;
  for (i = 0; i < source.length; ++i)
    destination.data[i * destination.stride] = source.data[i * source.stride];
  return 0;
}

int chem_vector_fill(chem_vector_view destination, double value)
{
  size_t i;
  if (destination.length && (!destination.data || !destination.stride))
    return -1;
  for (i = 0; i < destination.length; ++i)
    destination.data[i * destination.stride] = value;
  return 0;
}

int chem_vector_axpy(chem_vector_view destination, double alpha,
                     chem_const_vector_view source)
{
  size_t i;
  if (destination.length != source.length ||
      (destination.length && (!destination.data || !source.data ||
                              !destination.stride || !source.stride)))
    return -1;
  for (i = 0; i < source.length; ++i)
    destination.data[i * destination.stride] += alpha * source.data[i * source.stride];
  return 0;
}

typedef struct {
  uint32_t major;
  uint32_t minor;
} chem_sparse_pair;

static int compare_sparse_pair(const void *left, const void *right)
{
  const chem_sparse_pair *a = (const chem_sparse_pair *)left;
  const chem_sparse_pair *b = (const chem_sparse_pair *)right;
  if (a->major != b->major)
    return a->major < b->major ? -1 : 1;
  if (a->minor != b->minor)
    return a->minor < b->minor ? -1 : 1;
  return 0;
}

static size_t sparse_value_index(const chem_sparse_matrix *matrix,
                                 size_t block, size_t entry)
{
  return ((block / matrix->lanes) * matrix->nonzeros + entry) * matrix->lanes +
         block % matrix->lanes;
}

void chem_sparse_matrix_destroy(chem_sparse_matrix *matrix)
{
  if (!matrix)
    return;
  free(matrix->offsets);
  free(matrix->indices);
  free(matrix->values);
  *matrix = (chem_sparse_matrix){0};
}

int chem_sparse_matrix_init(chem_sparse_matrix *matrix, size_t rows,
                            size_t columns, const uint32_t *entry_rows,
                            const uint32_t *entry_columns, size_t nonzeros,
                            size_t blocks, size_t lanes,
                            chem_sparse_format format)
{
  const size_t major_count = format == CHEM_SPARSE_CSR ? rows : columns;
  chem_sparse_pair *pairs = NULL;
  size_t value_count, i;
  if (!matrix || (format != CHEM_SPARSE_CSR && format != CHEM_SPARSE_CSC) ||
      !lanes || (nonzeros && (!entry_rows || !entry_columns)) ||
      major_count == SIZE_MAX || blocks > SIZE_MAX - (lanes - 1) ||
      nonzeros > SIZE_MAX / sizeof(*matrix->indices) ||
      nonzeros > SIZE_MAX / sizeof(chem_sparse_pair))
    return -1;
  *matrix = (chem_sparse_matrix){0};
  matrix->rows = rows;
  matrix->columns = columns;
  matrix->nonzeros = nonzeros;
  matrix->blocks = blocks;
  matrix->lanes = lanes;
  matrix->groups = blocks ? (blocks + lanes - 1) / lanes : 0;
  matrix->format = format;
  matrix->offsets = (size_t *)calloc(major_count + 1, sizeof(*matrix->offsets));
  matrix->indices = (uint32_t *)malloc((nonzeros ? nonzeros : 1) * sizeof(*matrix->indices));
  pairs = (chem_sparse_pair *)malloc((nonzeros ? nonzeros : 1) * sizeof(*pairs));
  if (!matrix->offsets || !matrix->indices || !pairs)
    goto allocation_failure;
  for (i = 0; i < nonzeros; ++i) {
    const uint32_t row = entry_rows[i];
    const uint32_t column = entry_columns[i];
    if ((size_t)row >= rows || (size_t)column >= columns)
      goto invalid_pattern;
    pairs[i].major = format == CHEM_SPARSE_CSR ? row : column;
    pairs[i].minor = format == CHEM_SPARSE_CSR ? column : row;
  }
  qsort(pairs, nonzeros, sizeof(*pairs), compare_sparse_pair);
  for (i = 0; i < nonzeros; ++i) {
    if (i && pairs[i].major == pairs[i - 1].major &&
        pairs[i].minor == pairs[i - 1].minor)
      goto invalid_pattern;
    ++matrix->offsets[(size_t)pairs[i].major + 1];
    matrix->indices[i] = pairs[i].minor;
  }
  for (i = 0; i < major_count; ++i)
    matrix->offsets[i + 1] += matrix->offsets[i];
  if (matrix->groups && nonzeros > SIZE_MAX / matrix->groups)
    goto allocation_failure;
  value_count = matrix->groups * nonzeros;
  if (value_count && lanes > SIZE_MAX / value_count)
    goto allocation_failure;
  value_count *= lanes;
  matrix->values = (double *)calloc(value_count ? value_count : 1,
                                    sizeof(*matrix->values));
  free(pairs);
  if (!matrix->values)
    goto allocation_failure;
  return 0;

invalid_pattern:
  free(pairs);
  chem_sparse_matrix_destroy(matrix);
  return -1;
allocation_failure:
  free(pairs);
  chem_sparse_matrix_destroy(matrix);
  return -1;
}

static size_t sparse_find_entry(const chem_sparse_matrix *matrix,
                                size_t row, size_t column)
{
  size_t major;
  uint32_t minor;
  size_t low, high;
  if (!matrix || row >= matrix->rows || column >= matrix->columns)
    return SIZE_MAX;
  major = matrix->format == CHEM_SPARSE_CSR ? row : column;
  minor = (uint32_t)(matrix->format == CHEM_SPARSE_CSR ? column : row);
  low = matrix->offsets[major];
  high = matrix->offsets[major + 1];
  while (low < high) {
    const size_t middle = low + (high - low) / 2;
    if (matrix->indices[middle] < minor)
      low = middle + 1;
    else
      high = middle;
  }
  return low < matrix->offsets[major + 1] && matrix->indices[low] == minor
             ? low
             : SIZE_MAX;
}

int chem_sparse_matrix_contains(const chem_sparse_matrix *matrix, size_t row,
                                size_t column)
{
  if (!matrix || row >= matrix->rows || column >= matrix->columns)
    return -1;
  return sparse_find_entry(matrix, row, column) == SIZE_MAX ? 0 : 1;
}

int chem_sparse_matrix_get(const chem_sparse_matrix *matrix, size_t block,
                           size_t row, size_t column, double *value)
{
  size_t entry;
  if (!matrix || !value || block >= matrix->blocks || row >= matrix->rows ||
      column >= matrix->columns)
    return -1;
  entry = sparse_find_entry(matrix, row, column);
  *value = entry == SIZE_MAX ? 0.0 : matrix->values[sparse_value_index(matrix, block, entry)];
  return 0;
}

int chem_sparse_matrix_set(chem_sparse_matrix *matrix, size_t block,
                           size_t row, size_t column, double value)
{
  size_t entry;
  if (!matrix || block >= matrix->blocks || row >= matrix->rows ||
      column >= matrix->columns)
    return -1;
  entry = sparse_find_entry(matrix, row, column);
  if (entry == SIZE_MAX)
    return -1;
  matrix->values[sparse_value_index(matrix, block, entry)] = value;
  return 0;
}

int chem_sparse_matrix_fill(chem_sparse_matrix *matrix, double value)
{
  size_t i, count;
  if (!matrix || matrix->groups > SIZE_MAX / (matrix->nonzeros ? matrix->nonzeros : 1))
    return -1;
  count = matrix->groups * matrix->nonzeros;
  if (count && matrix->lanes > SIZE_MAX / count)
    return -1;
  count *= matrix->lanes;
  for (i = 0; i < count; ++i)
    matrix->values[i] = value;
  return 0;
}

int chem_sparse_matrix_add_diagonal(chem_sparse_matrix *matrix, double value)
{
  size_t i, count;
  if (!matrix)
    return -1;
  count = matrix->rows < matrix->columns ? matrix->rows : matrix->columns;
  for (i = 0; i < count; ++i)
    if (sparse_find_entry(matrix, i, i) == SIZE_MAX)
      return -1;
  for (i = 0; i < count; ++i) {
    size_t block;
    for (block = 0; block < matrix->blocks; ++block) {
      const size_t entry = sparse_find_entry(matrix, i, i);
      matrix->values[sparse_value_index(matrix, block, entry)] += value;
    }
  }
  return 0;
}

static int sparse_same_pattern(const chem_sparse_matrix *a,
                               const chem_sparse_matrix *b)
{
  size_t major, entry;
  size_t a_major;
  if (!a || !b || a->rows != b->rows || a->columns != b->columns ||
      a->nonzeros != b->nonzeros)
    return 0;
  a_major = a->format == CHEM_SPARSE_CSR ? a->rows : a->columns;
  for (major = 0; major < a_major; ++major) {
    for (entry = a->offsets[major]; entry < a->offsets[major + 1]; ++entry) {
      const size_t row = a->format == CHEM_SPARSE_CSR ? major : a->indices[entry];
      const size_t column = a->format == CHEM_SPARSE_CSR ? a->indices[entry] : major;
      if (chem_sparse_matrix_contains(b, row, column) != 1)
        return 0;
    }
  }
  return 1;
}

int chem_sparse_matrix_copy(chem_sparse_matrix *destination,
                            const chem_sparse_matrix *source)
{
  size_t block, major, entry;
  const size_t major_count = source && source->format == CHEM_SPARSE_CSR
                                 ? source->rows
                                 : source ? source->columns : 0;
  if (!destination || !source || destination->blocks != source->blocks ||
      !sparse_same_pattern(destination, source))
    return -1;
  for (block = 0; block < source->blocks; ++block) {
    for (major = 0; major < major_count; ++major) {
      for (entry = source->offsets[major]; entry < source->offsets[major + 1]; ++entry) {
        const size_t row = source->format == CHEM_SPARSE_CSR ? major : source->indices[entry];
        const size_t column = source->format == CHEM_SPARSE_CSR ? source->indices[entry] : major;
        const double value = source->values[sparse_value_index(source, block, entry)];
        if (chem_sparse_matrix_set(destination, block, row, column, value) != 0)
          return -1;
      }
    }
  }
  return 0;
}

int chem_sparse_matrix_axpy(chem_sparse_matrix *destination, double alpha,
                            const chem_sparse_matrix *source)
{
  size_t block, major, entry;
  const size_t major_count = source && source->format == CHEM_SPARSE_CSR
                                 ? source->rows
                                 : source ? source->columns : 0;
  if (!destination || !source || destination->blocks != source->blocks ||
      !sparse_same_pattern(destination, source))
    return -1;
  for (block = 0; block < source->blocks; ++block) {
    for (major = 0; major < major_count; ++major) {
      for (entry = source->offsets[major]; entry < source->offsets[major + 1]; ++entry) {
        const size_t row = source->format == CHEM_SPARSE_CSR ? major : source->indices[entry];
        const size_t column = source->format == CHEM_SPARSE_CSR ? source->indices[entry] : major;
        const double source_value = source->values[sparse_value_index(source, block, entry)];
        const size_t destination_entry = sparse_find_entry(destination, row, column);
        destination->values[sparse_value_index(destination, block, destination_entry)] +=
            alpha * source_value;
      }
    }
  }
  return 0;
}

int chem_sparse_matrix_transform(chem_sparse_matrix *matrix,
                                 chem_unary_function function, void *context)
{
  size_t i, count;
  if (!matrix || !function || matrix->groups > SIZE_MAX / (matrix->nonzeros ? matrix->nonzeros : 1))
    return -1;
  count = matrix->groups * matrix->nonzeros;
  if (count && matrix->lanes > SIZE_MAX / count)
    return -1;
  count *= matrix->lanes;
  for (i = 0; i < count; ++i)
    matrix->values[i] = function(matrix->values[i], context);
  return 0;
}

int chem_sparse_matrix_validate_blocks(const chem_sparse_matrix *matrix,
                                       size_t block_count)
{
  return matrix && matrix->blocks == block_count ? 0 : -1;
}

int chem_lu_decompose_dense(const double *matrix, size_t order,
                            double *lower, double *upper)
{
  size_t i, j, k;
  if ((!matrix || !lower || !upper) && order)
    return -1;
  if (order && (order > SIZE_MAX / order ||
                order * order > SIZE_MAX / sizeof(*lower)))
    return -1;
  if (order) {
    memset(lower, 0, order * order * sizeof(*lower));
    memset(upper, 0, order * order * sizeof(*upper));
  }
  for (i = 0; i < order; ++i) {
    for (j = 0; j < order; ++j) {
      if (i <= j) {
        double value = matrix[i * order + j];
        for (k = 0; k < i; ++k)
          value -= lower[i * order + k] * upper[k * order + j];
        upper[i * order + j] = value;
      } else {
        double value = matrix[i * order + j];
        for (k = 0; k < j; ++k)
          value -= lower[i * order + k] * upper[k * order + j];
        if (upper[j * order + j] == 0.0 || !isfinite(upper[j * order + j]))
          return -1;
        lower[i * order + j] = value / upper[j * order + j];
      }
    }
    lower[i * order + i] = 1.0;
  }
  return 0;
}

int chem_lu_factor_dense(double *matrix, size_t order)
{
  size_t k, i, j;
  if (!matrix && order)
    return -1;
  if (order && order > SIZE_MAX / order)
    return -1;
  for (k = 0; k < order; ++k) {
    const double pivot = matrix[k * order + k];
    if (pivot == 0.0 || !isfinite(pivot))
      return -1;
    for (i = k + 1; i < order; ++i) {
      matrix[i * order + k] /= pivot;
      for (j = k + 1; j < order; ++j)
        matrix[i * order + j] -= matrix[i * order + k] * matrix[k * order + j];
    }
  }
  return 0;
}

int chem_lu_solve_dense(const double *factor, size_t order,
                        double *right_hand_side)
{
  size_t i, j;
  if ((!factor || !right_hand_side) && order)
    return -1;
  for (i = 0; i < order; ++i) {
    for (j = 0; j < i; ++j)
      right_hand_side[i] -= factor[i * order + j] * right_hand_side[j];
  }
  for (i = order; i-- > 0;) {
    for (j = i + 1; j < order; ++j)
      right_hand_side[i] -= factor[i * order + j] * right_hand_side[j];
    if (factor[i * order + i] == 0.0 || !isfinite(factor[i * order + i]))
      return -1;
    right_hand_side[i] /= factor[i * order + i];
  }
  return 0;
}
