#include "chem_sparse_matrix.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int sparse_valid(const chem_sparse_matrix *matrix) {
  return matrix && matrix->dimension && matrix->lanes &&
      (matrix->nonzero_count == 0u || (matrix->rows && matrix->columns)) &&
      (matrix->storage_length == 0u || matrix->values) &&
      (matrix->order == CHEM_SPARSE_ROW_COMPRESSED || matrix->order == CHEM_SPARSE_COLUMN_COMPRESSED) &&
      (matrix->layout == CHEM_SPARSE_STANDARD_BLOCKS || matrix->layout == CHEM_SPARSE_GROUPED_BLOCKS);
}

static int coordinate_before(chem_sparse_order order, size_t row_a, size_t column_a, size_t row_b, size_t column_b) {
  return order == CHEM_SPARSE_ROW_COMPRESSED
      ? row_a < row_b || (row_a == row_b && column_a < column_b)
      : column_a < column_b || (column_a == column_b && row_a < row_b);
}

static size_t find_entry(const chem_sparse_matrix *matrix, size_t row, size_t column) {
  size_t first = 0u, last = matrix->nonzero_count;
  while (first < last) {
    const size_t middle = first + (last - first) / 2u;
    if (coordinate_before(matrix->order, matrix->rows[middle], matrix->columns[middle], row, column))
      first = middle + 1u;
    else
      last = middle;
  }
  if (first < matrix->nonzero_count && matrix->rows[first] == row && matrix->columns[first] == column)
    return first;
  return SIZE_MAX;
}

static size_t value_index(const chem_sparse_matrix *matrix, size_t block, size_t entry) {
  if (matrix->layout == CHEM_SPARSE_STANDARD_BLOCKS) return block * matrix->nonzero_count + entry;
  return ((block / matrix->lanes * matrix->nonzero_count + entry) * matrix->lanes) + block % matrix->lanes;
}

void chem_sparse_matrix_destroy(chem_sparse_matrix *matrix) {
  if (!matrix) return;
  free(matrix->rows);
  free(matrix->columns);
  free(matrix->values);
  memset(matrix, 0, sizeof(*matrix));
}

int chem_sparse_matrix_init(
    chem_sparse_matrix *matrix, size_t dimension,
    const size_t *rows, const size_t *columns, size_t nonzero_count,
    size_t block_count, size_t lanes, chem_sparse_order order,
    chem_sparse_layout layout, double initial_value) {
  if (!matrix || !dimension || !lanes || (nonzero_count && (!rows || !columns)) ||
      (order != CHEM_SPARSE_ROW_COMPRESSED && order != CHEM_SPARSE_COLUMN_COMPRESSED) ||
      (layout != CHEM_SPARSE_STANDARD_BLOCKS && layout != CHEM_SPARSE_GROUPED_BLOCKS) ||
      (layout == CHEM_SPARSE_STANDARD_BLOCKS && lanes != 1u) ||
      nonzero_count > SIZE_MAX / sizeof(size_t)) return -1;

  chem_sparse_matrix result = { 0 };
  result.dimension = dimension;
  result.nonzero_count = nonzero_count;
  result.block_count = block_count;
  result.lanes = lanes;
  result.order = order;
  result.layout = layout;
  result.rows = nonzero_count ? malloc(nonzero_count * sizeof(*result.rows)) : NULL;
  result.columns = nonzero_count ? malloc(nonzero_count * sizeof(*result.columns)) : NULL;
  if (nonzero_count && (!result.rows || !result.columns)) {
    chem_sparse_matrix_destroy(&result);
    return -2;
  }
  for (size_t entry = 0; entry < nonzero_count; ++entry) {
    if (rows[entry] >= dimension || columns[entry] >= dimension) {
      chem_sparse_matrix_destroy(&result);
      return -1;
    }
    result.rows[entry] = rows[entry];
    result.columns[entry] = columns[entry];
  }
  for (size_t i = 1u; i < nonzero_count; ++i) {
    const size_t row = result.rows[i], column = result.columns[i];
    size_t j = i;
    while (j > 0u && coordinate_before(order, row, column, result.rows[j - 1u], result.columns[j - 1u])) {
      result.rows[j] = result.rows[j - 1u];
      result.columns[j] = result.columns[j - 1u];
      --j;
    }
    result.rows[j] = row;
    result.columns[j] = column;
  }
  for (size_t entry = 1u; entry < nonzero_count; ++entry) {
    if (result.rows[entry - 1u] == result.rows[entry] && result.columns[entry - 1u] == result.columns[entry]) {
      chem_sparse_matrix_destroy(&result);
      return -1;
    }
  }

  if (layout == CHEM_SPARSE_STANDARD_BLOCKS) {
    if (nonzero_count && block_count > SIZE_MAX / nonzero_count) {
      chem_sparse_matrix_destroy(&result);
      return -1;
    }
    result.storage_length = block_count * nonzero_count;
  } else {
    const size_t groups = block_count / lanes + (block_count % lanes != 0u);
    if (nonzero_count && groups > SIZE_MAX / nonzero_count) {
      chem_sparse_matrix_destroy(&result);
      return -1;
    }
    result.storage_length = groups * nonzero_count;
    if (result.storage_length && lanes > SIZE_MAX / result.storage_length) {
      chem_sparse_matrix_destroy(&result);
      return -1;
    }
    result.storage_length *= lanes;
  }
  if (result.storage_length > SIZE_MAX / sizeof(*result.values)) {
    chem_sparse_matrix_destroy(&result);
    return -1;
  }
  result.values = result.storage_length ? malloc(result.storage_length * sizeof(*result.values)) : NULL;
  if (result.storage_length && !result.values) {
    chem_sparse_matrix_destroy(&result);
    return -2;
  }
  for (size_t i = 0; i < result.storage_length; ++i) result.values[i] = initial_value;
  *matrix = result;
  return 0;
}

int chem_sparse_matrix_vector_index(
    const chem_sparse_matrix *matrix, size_t row, size_t column, size_t *index) {
  if (!sparse_valid(matrix) || !index || row >= matrix->dimension || column >= matrix->dimension) return -1;
  const size_t position = find_entry(matrix, row, column);
  if (position == SIZE_MAX) return -1;
  *index = position;
  return 0;
}

int chem_sparse_matrix_is_zero(
    const chem_sparse_matrix *matrix, size_t row, size_t column, int *is_zero) {
  if (!sparse_valid(matrix) || !is_zero || row >= matrix->dimension || column >= matrix->dimension) return -1;
  *is_zero = find_entry(matrix, row, column) == SIZE_MAX;
  return 0;
}

int chem_sparse_matrix_get(
    const chem_sparse_matrix *matrix, size_t block, size_t row, size_t column, double *value) {
  if (!sparse_valid(matrix) || !value || block >= matrix->block_count ||
      row >= matrix->dimension || column >= matrix->dimension) return -1;
  const size_t entry = find_entry(matrix, row, column);
  *value = entry == SIZE_MAX ? 0.0 : matrix->values[value_index(matrix, block, entry)];
  return 0;
}

int chem_sparse_matrix_set(
    chem_sparse_matrix *matrix, size_t block, size_t row, size_t column, double value) {
  if (!sparse_valid(matrix) || block >= matrix->block_count ||
      row >= matrix->dimension || column >= matrix->dimension) return -1;
  const size_t entry = find_entry(matrix, row, column);
  if (entry == SIZE_MAX) return -1;
  matrix->values[value_index(matrix, block, entry)] = value;
  return 0;
}

static int block_group_view_parameters(
    const chem_sparse_matrix *matrix, size_t group, size_t entry,
    size_t *offset, size_t *length, size_t *stride) {
  if (!sparse_valid(matrix) || entry >= matrix->nonzero_count) return -1;
  const size_t groups = matrix->block_count / matrix->lanes + (matrix->block_count % matrix->lanes != 0u);
  if (group >= groups) return -1;
  const size_t first_block = group * matrix->lanes;
  *length = matrix->block_count - first_block;
  if (*length > matrix->lanes) *length = matrix->lanes;
  if (matrix->layout == CHEM_SPARSE_STANDARD_BLOCKS) {
    *offset = first_block * matrix->nonzero_count + entry;
    *stride = matrix->nonzero_count;
  } else {
    *offset = (group * matrix->nonzero_count + entry) * matrix->lanes;
    *stride = 1u;
  }
  return 0;
}

int chem_sparse_matrix_block_group_view(
    chem_sparse_matrix *matrix, size_t group, size_t entry, chem_group_view *view) {
  if (!view) return -1;
  size_t offset, length, stride;
  if (block_group_view_parameters(matrix, group, entry, &offset, &length, &stride) != 0) return -1;
  return chem_group_view_init(matrix->values, matrix->storage_length, offset, length, stride, view);
}

int chem_sparse_matrix_const_block_group_view(
    const chem_sparse_matrix *matrix, size_t group, size_t entry, chem_const_group_view *view) {
  if (!view) return -1;
  size_t offset, length, stride;
  if (block_group_view_parameters(matrix, group, entry, &offset, &length, &stride) != 0) return -1;
  return chem_const_group_view_init(matrix->values, matrix->storage_length, offset, length, stride, view);
}

int chem_sparse_matrix_fill(chem_sparse_matrix *matrix, double value) {
  if (!sparse_valid(matrix)) return -1;
  for (size_t i = 0; i < matrix->storage_length; ++i) matrix->values[i] = value;
  return 0;
}

int chem_sparse_matrix_copy(chem_sparse_matrix *destination, const chem_sparse_matrix *source) {
  if (!sparse_valid(destination) || !sparse_valid(source) ||
      destination->dimension != source->dimension || destination->nonzero_count != source->nonzero_count ||
      destination->block_count != source->block_count || destination->lanes != source->lanes ||
      destination->layout != source->layout || destination->order != source->order) return -1;
  for (size_t entry = 0; entry < source->nonzero_count; ++entry)
    if (destination->rows[entry] != source->rows[entry] || destination->columns[entry] != source->columns[entry]) return -1;
  if (destination->storage_length != source->storage_length) return -1;
  if (destination->storage_length)
    memmove(destination->values, source->values, destination->storage_length * sizeof(*source->values));
  return 0;
}

int chem_sparse_matrix_add_to_diagonal(chem_sparse_matrix *matrix, double value) {
  if (!sparse_valid(matrix)) return -1;
  for (size_t entry = 0; entry < matrix->nonzero_count; ++entry) {
    if (matrix->rows[entry] != matrix->columns[entry]) continue;
    for (size_t block = 0; block < matrix->block_count; ++block)
      matrix->values[value_index(matrix, block, entry)] += value;
  }
  return 0;
}

int chem_sparse_matrix_print(const chem_sparse_matrix *matrix, FILE *stream) {
  if (!sparse_valid(matrix) || !stream) return -1;
  for (size_t block = 0; block < matrix->block_count; ++block) {
    if (fprintf(stream, "Block %zu\n", block) < 0) return -1;
    for (size_t row = 0; row < matrix->dimension; ++row) {
      for (size_t column = 0; column < matrix->dimension; ++column) {
        if (column && fputc(',', stream) == EOF) return -1;
        double value;
        if (chem_sparse_matrix_get(matrix, block, row, column, &value) != 0 || fprintf(stream, "%g", value) < 0) return -1;
      }
      if (fputc('\n', stream) == EOF) return -1;
    }
  }
  return 0;
}

int chem_sparse_matrix_print_nonzero(const chem_sparse_matrix *matrix, FILE *stream) {
  if (!sparse_valid(matrix) || !stream) return -1;
  for (size_t block = 0; block < matrix->block_count; ++block) {
    if (fprintf(stream, "Block %zu\n", block) < 0) return -1;
    for (size_t entry = 0; entry < matrix->nonzero_count; ++entry)
      if (fprintf(stream, "%zu, %zu, %g\n", matrix->rows[entry], matrix->columns[entry],
              matrix->values[value_index(matrix, block, entry)]) < 0) return -1;
  }
  return 0;
}
