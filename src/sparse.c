#include "chemistry.h"

#include <math.h>

static int csr_valid(const chem_csr_matrix *matrix) {
  if (!matrix || !matrix->row_offsets ||
      (matrix->nonzero_count && !matrix->column_indices) ||
      (matrix->block_count && matrix->nonzero_count && !matrix->values) ||
      matrix->dimension == (size_t)-1 ||
      (matrix->block_count && matrix->nonzero_count >
          (size_t)-1 / matrix->block_count) ||
      matrix->row_offsets[0] != 0 ||
      matrix->row_offsets[matrix->dimension] != matrix->nonzero_count)
    return 0;
  for (size_t row = 0; row < matrix->dimension; ++row) {
    const size_t begin = matrix->row_offsets[row];
    const size_t end = matrix->row_offsets[row + 1];
    if (begin > end || end > matrix->nonzero_count) return 0;
    for (size_t k = begin; k < end; ++k)
      if (matrix->column_indices[k] >= matrix->dimension ||
          (k > begin && matrix->column_indices[k - 1] >=
                            matrix->column_indices[k])) return 0;
  }
  return 1;
}

int chem_csr_matrix_init(chem_csr_matrix *matrix, size_t dimension,
                         size_t block_count, size_t nonzero_count,
                         const size_t *row_offsets,
                         const size_t *column_indices, double *values) {
  if (!matrix) return 1;
  chem_csr_matrix candidate = {dimension, block_count, nonzero_count,
                               row_offsets, column_indices, values};
  if (!csr_valid(&candidate)) return 1;
  *matrix = candidate;
  return 0;
}

int chem_csr_matrix_vector_index(const chem_csr_matrix *matrix, size_t block,
                                 size_t row, size_t column, size_t *index) {
  if (!csr_valid(matrix) || !index || block >= matrix->block_count ||
      row >= matrix->dimension || column >= matrix->dimension) return 1;
  size_t lower = matrix->row_offsets[row];
  size_t upper = matrix->row_offsets[row + 1];
  while (lower < upper) {
    const size_t middle = lower + (upper - lower) / 2;
    if (matrix->column_indices[middle] < column) lower = middle + 1;
    else upper = middle;
  }
  if (lower >= matrix->row_offsets[row + 1] ||
      matrix->column_indices[lower] != column) return 1;
  *index = block * matrix->nonzero_count + lower;
  return 0;
}

static int csr_entry(const chem_csr_matrix *matrix, size_t row, size_t column,
                     size_t *entry) {
  if (!csr_valid(matrix) || !entry || row >= matrix->dimension ||
      column >= matrix->dimension) return 1;
  size_t lower = matrix->row_offsets[row];
  size_t upper = matrix->row_offsets[row + 1];
  while (lower < upper) {
    const size_t middle = lower + (upper - lower) / 2;
    if (matrix->column_indices[middle] < column) lower = middle + 1;
    else upper = middle;
  }
  if (lower >= matrix->row_offsets[row + 1] ||
      matrix->column_indices[lower] != column) return 1;
  *entry = lower;
  return 0;
}

int chem_csr_matrix_get_value_view(const chem_csr_matrix *matrix, size_t row,
                                   size_t column, chem_csr_value_view *view) {
  size_t entry = 0;
  if (!view || csr_entry(matrix,row,column,&entry)) return 1;
  view->data = matrix->values ? matrix->values + entry : NULL;
  view->count = matrix->block_count;
  view->stride = matrix->nonzero_count;
  return 0;
}

int chem_csr_matrix_get_const_value_view(const chem_csr_matrix *matrix,
                                        size_t row, size_t column,
                                        chem_csr_const_value_view *view) {
  size_t entry = 0;
  if (!view || csr_entry(matrix,row,column,&entry)) return 1;
  view->data = matrix->values ? matrix->values + entry : NULL;
  view->count = matrix->block_count;
  view->stride = matrix->nonzero_count;
  return 0;
}

int chem_csr_matrix_apply_blocks(chem_csr_value_view *destination,
                                 const chem_csr_const_value_view *sources,
                                 size_t source_count,
                                 chem_csr_block_function function,
                                 void *context) {
  if (!destination || !function || (source_count && !sources) ||
      (destination->count && !destination->data)) return 1;
  for (size_t i = 0; i < source_count; ++i)
    if (sources[i].count != destination->count ||
        (sources[i].count && !sources[i].data)) return 1;
  for (size_t block = 0; block < destination->count; ++block)
    function(block,destination->data + block*destination->stride,sources,
             source_count,context);
  return 0;
}

int chem_csr_matrix_get(const chem_csr_matrix *matrix, size_t block,
                        size_t row, size_t column, double *value) {
  size_t index = 0;
  if (!value || chem_csr_matrix_vector_index(matrix,block,row,column,&index))
    return 1;
  *value = matrix->values[index];
  return 0;
}

int chem_csr_matrix_set(chem_csr_matrix *matrix, size_t block, size_t row,
                        size_t column, double value) {
  size_t index = 0;
  if (chem_csr_matrix_vector_index(matrix,block,row,column,&index)) return 1;
  matrix->values[index] = value;
  return 0;
}

int chem_csr_matrix_is_structural_zero(const chem_csr_matrix *matrix,
                                       size_t row, size_t column,
                                       int *is_zero) {
  if (!csr_valid(matrix) || !is_zero || row >= matrix->dimension ||
      column >= matrix->dimension) return 1;
  size_t lower = matrix->row_offsets[row];
  size_t upper = matrix->row_offsets[row + 1];
  while (lower < upper) {
    const size_t middle = lower + (upper - lower) / 2;
    if (matrix->column_indices[middle] < column) lower = middle + 1;
    else upper = middle;
  }
  *is_zero = lower >= matrix->row_offsets[row + 1] ||
             matrix->column_indices[lower] != column;
  return 0;
}

int chem_csr_matrix_fill(chem_csr_matrix *matrix, double value) {
  if (!csr_valid(matrix) || !isfinite(value)) return 1;
  const size_t count = matrix->block_count * matrix->nonzero_count;
  for (size_t k = 0; k < count; ++k) matrix->values[k] = value;
  return 0;
}

int chem_csr_matrix_add_to_diagonal(chem_csr_matrix *matrix, double value) {
  if (!csr_valid(matrix) || !isfinite(value)) return 1;
  for (size_t row = 0; row < matrix->dimension; ++row) {
    size_t lower = matrix->row_offsets[row];
    size_t upper = matrix->row_offsets[row + 1];
    while (lower < upper) {
      const size_t middle = lower + (upper - lower) / 2;
      if (matrix->column_indices[middle] < row) lower = middle + 1;
      else upper = middle;
    }
    if (lower < matrix->row_offsets[row + 1] &&
        matrix->column_indices[lower] == row)
      for (size_t block = 0; block < matrix->block_count; ++block)
        matrix->values[block*matrix->nonzero_count + lower] += value;
  }
  return 0;
}
