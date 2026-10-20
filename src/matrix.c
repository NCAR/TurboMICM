#include "chem_matrix.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int matrix_valid(const chem_matrix *matrix) {
  return matrix && matrix->lanes != 0u &&
      (matrix->layout == CHEM_MATRIX_ROW_MAJOR || matrix->layout == CHEM_MATRIX_GROUPED) &&
      (matrix->storage_length == 0u || matrix->data != NULL);
}

static size_t matrix_index(const chem_matrix *matrix, size_t row, size_t column) {
  if (matrix->layout == CHEM_MATRIX_ROW_MAJOR) return row * matrix->columns + column;
  return (row / matrix->lanes * matrix->columns + column) * matrix->lanes + row % matrix->lanes;
}

static int matrix_same_shape(const chem_matrix *a, const chem_matrix *b) {
  return matrix_valid(a) && matrix_valid(b) && a->rows == b->rows && a->columns == b->columns;
}

static int logical_element_count(size_t rows, size_t columns, size_t *count) {
  if (columns != 0u && rows > SIZE_MAX / columns) return -1;
  *count = rows * columns;
  return 0;
}

static int matrix_allocate_temporary(const chem_matrix *matrix, double **temporary, size_t *count) {
  if (logical_element_count(matrix->rows, matrix->columns, count) != 0 || *count > SIZE_MAX / sizeof(double)) return -1;
  *temporary = *count ? malloc(*count * sizeof(**temporary)) : NULL;
  return *count && !*temporary ? -2 : 0;
}

int chem_matrix_init(
    chem_matrix *matrix, size_t rows, size_t columns, size_t lanes,
    chem_matrix_layout layout, double initial_value) {
  if (!matrix || !lanes ||
      (layout != CHEM_MATRIX_ROW_MAJOR && layout != CHEM_MATRIX_GROUPED) ||
      (layout == CHEM_MATRIX_ROW_MAJOR && lanes != 1u)) return -1;
  size_t storage_length = 0u;
  if (layout == CHEM_MATRIX_ROW_MAJOR) {
    if (columns && rows > SIZE_MAX / columns) return -1;
    storage_length = rows * columns;
  } else {
    const size_t groups = rows / lanes + (rows % lanes != 0u);
    if (columns && groups > SIZE_MAX / columns) return -1;
    storage_length = groups * columns;
    if (storage_length && lanes > SIZE_MAX / storage_length) return -1;
    storage_length *= lanes;
  }
  if (storage_length > SIZE_MAX / sizeof(double)) return -1;
  double *data = storage_length ? malloc(storage_length * sizeof(*data)) : NULL;
  if (storage_length && !data) return -2;
  for (size_t i = 0; i < storage_length; ++i) data[i] = initial_value;
  *matrix = (chem_matrix){ rows, columns, lanes, storage_length, layout, data };
  return 0;
}

void chem_matrix_destroy(chem_matrix *matrix) {
  if (!matrix) return;
  free(matrix->data);
  memset(matrix, 0, sizeof(*matrix));
}

int chem_matrix_get(const chem_matrix *matrix, size_t row, size_t column, double *value) {
  if (!matrix_valid(matrix) || !value || row >= matrix->rows || column >= matrix->columns) return -1;
  *value = matrix->data[matrix_index(matrix, row, column)];
  return 0;
}

int chem_matrix_set(chem_matrix *matrix, size_t row, size_t column, double value) {
  if (!matrix_valid(matrix) || row >= matrix->rows || column >= matrix->columns) return -1;
  matrix->data[matrix_index(matrix, row, column)] = value;
  return 0;
}

int chem_matrix_row_view(chem_matrix *matrix, size_t row, chem_group_view *view) {
  if (!matrix_valid(matrix) || !view || row >= matrix->rows) return -1;
  if (matrix->columns == 0u)
    return chem_group_view_init(matrix->data, matrix->storage_length, 0u, 0u, 1u, view);
  const size_t offset = matrix_index(matrix, row, 0u);
  const size_t stride = matrix->layout == CHEM_MATRIX_ROW_MAJOR ? 1u : matrix->lanes;
  return chem_group_view_init(matrix->data, matrix->storage_length, offset, matrix->columns, stride, view);
}

int chem_matrix_const_row_view(const chem_matrix *matrix, size_t row, chem_const_group_view *view) {
  if (!matrix_valid(matrix) || !view || row >= matrix->rows) return -1;
  if (matrix->columns == 0u)
    return chem_const_group_view_init(matrix->data, matrix->storage_length, 0u, 0u, 1u, view);
  const size_t offset = matrix_index(matrix, row, 0u);
  const size_t stride = matrix->layout == CHEM_MATRIX_ROW_MAJOR ? 1u : matrix->lanes;
  return chem_const_group_view_init(matrix->data, matrix->storage_length, offset, matrix->columns, stride, view);
}

int chem_matrix_column_group_view(chem_matrix *matrix, size_t group, size_t column, chem_group_view *view) {
  if (!matrix_valid(matrix) || !view || column >= matrix->columns) return -1;
  const size_t groups = matrix->rows / matrix->lanes + (matrix->rows % matrix->lanes != 0u);
  if (group >= groups) return -1;
  const size_t first_row = group * matrix->lanes;
  size_t count = matrix->rows - first_row;
  if (count > matrix->lanes) count = matrix->lanes;
  const size_t offset = matrix_index(matrix, first_row, column);
  const size_t stride = matrix->layout == CHEM_MATRIX_ROW_MAJOR ? matrix->columns : 1u;
  return chem_group_view_init(matrix->data, matrix->storage_length, offset, count, stride, view);
}

int chem_matrix_const_column_group_view(
    const chem_matrix *matrix, size_t group, size_t column, chem_const_group_view *view) {
  if (!matrix_valid(matrix) || !view || column >= matrix->columns) return -1;
  const size_t groups = matrix->rows / matrix->lanes + (matrix->rows % matrix->lanes != 0u);
  if (group >= groups) return -1;
  const size_t first_row = group * matrix->lanes;
  size_t count = matrix->rows - first_row;
  if (count > matrix->lanes) count = matrix->lanes;
  const size_t offset = matrix_index(matrix, first_row, column);
  const size_t stride = matrix->layout == CHEM_MATRIX_ROW_MAJOR ? matrix->columns : 1u;
  return chem_const_group_view_init(matrix->data, matrix->storage_length, offset, count, stride, view);
}

int chem_matrix_fill(chem_matrix *matrix, double value) {
  if (!matrix_valid(matrix)) return -1;
  for (size_t i = 0; i < matrix->storage_length; ++i) matrix->data[i] = value;
  return 0;
}

int chem_matrix_copy(chem_matrix *destination, const chem_matrix *source) {
  if (!matrix_same_shape(destination, source)) return -1;
  for (size_t row = 0; row < source->rows; ++row)
    for (size_t column = 0; column < source->columns; ++column)
      destination->data[matrix_index(destination, row, column)] = source->data[matrix_index(source, row, column)];
  return 0;
}

int chem_matrix_assign_row(chem_matrix *matrix, size_t row, const double *values, size_t value_count) {
  if (!matrix_valid(matrix) || row >= matrix->rows || value_count < matrix->columns ||
      (matrix->columns && !values)) return -1;
  for (size_t column = 0; column < matrix->columns; ++column)
    matrix->data[matrix_index(matrix, row, column)] = values[column];
  return 0;
}

int chem_matrix_read_row(const chem_matrix *matrix, size_t row, double *values, size_t value_count) {
  if (!matrix_valid(matrix) || row >= matrix->rows || value_count < matrix->columns ||
      (matrix->columns && !values)) return -1;
  for (size_t column = 0; column < matrix->columns; ++column)
    values[column] = matrix->data[matrix_index(matrix, row, column)];
  return 0;
}

int chem_matrix_axpy(chem_matrix *destination, double alpha, const chem_matrix *source) {
  if (!matrix_same_shape(destination, source)) return -1;
  for (size_t row = 0; row < source->rows; ++row)
    for (size_t column = 0; column < source->columns; ++column) {
      const size_t dest_index = matrix_index(destination, row, column);
      const double value = source->data[matrix_index(source, row, column)];
      destination->data[dest_index] += alpha * value;
    }
  return 0;
}

int chem_matrix_max_scalar(chem_matrix *matrix, double lower_bound) {
  if (!matrix_valid(matrix)) return -1;
  for (size_t i = 0; i < matrix->storage_length; ++i)
    if (matrix->data[i] < lower_bound) matrix->data[i] = lower_bound;
  return 0;
}

int chem_matrix_min_scalar(chem_matrix *matrix, double upper_bound) {
  if (!matrix_valid(matrix)) return -1;
  for (size_t i = 0; i < matrix->storage_length; ++i)
    if (matrix->data[i] > upper_bound) matrix->data[i] = upper_bound;
  return 0;
}

int chem_matrix_min_max(const chem_matrix *matrix, double *minimum, double *maximum) {
  if (!matrix_valid(matrix) || !minimum || !maximum || matrix->rows == 0u || matrix->columns == 0u) return -1;
  *minimum = *maximum = matrix->data[matrix_index(matrix, 0u, 0u)];
  for (size_t row = 0; row < matrix->rows; ++row)
    for (size_t column = 0; column < matrix->columns; ++column) {
      const double value = matrix->data[matrix_index(matrix, row, column)];
      if (value < *minimum) *minimum = value;
      if (value > *maximum) *maximum = value;
    }
  return 0;
}

int chem_matrix_reduce(
    const chem_matrix *matrix, double *accumulator, chem_matrix_reducer reducer, void *context) {
  if (!matrix_valid(matrix) || !accumulator || !reducer) return -1;
  for (size_t row = 0; row < matrix->rows; ++row)
    for (size_t column = 0; column < matrix->columns; ++column)
      reducer(matrix->data[matrix_index(matrix, row, column)], accumulator, context);
  return 0;
}

int chem_matrix_transform_unary(
    chem_matrix *destination, const chem_matrix *source,
    chem_matrix_unary_function function, void *context) {
  if (!matrix_same_shape(destination, source) || !function) return -1;
  double *temporary = NULL;
  size_t count = 0u;
  const int allocation_status = matrix_allocate_temporary(source, &temporary, &count);
  if (allocation_status != 0) return allocation_status;
  for (size_t row = 0; row < source->rows; ++row)
    for (size_t column = 0; column < source->columns; ++column) {
      const size_t logical = row * source->columns + column;
      temporary[logical] = function(source->data[matrix_index(source, row, column)], context);
    }
  for (size_t row = 0; row < destination->rows; ++row)
    for (size_t column = 0; column < destination->columns; ++column)
      destination->data[matrix_index(destination, row, column)] = temporary[row * destination->columns + column];
  free(temporary);
  return 0;
}

int chem_matrix_transform_binary(
    chem_matrix *destination, const chem_matrix *left, const chem_matrix *right,
    chem_matrix_binary_function function, void *context) {
  if (!matrix_same_shape(destination, left) || !matrix_same_shape(destination, right) || !function) return -1;
  double *temporary = NULL;
  size_t count = 0u;
  const int allocation_status = matrix_allocate_temporary(left, &temporary, &count);
  if (allocation_status != 0) return allocation_status;
  for (size_t row = 0; row < left->rows; ++row)
    for (size_t column = 0; column < left->columns; ++column) {
      const size_t logical = row * left->columns + column;
      temporary[logical] = function(
          left->data[matrix_index(left, row, column)], right->data[matrix_index(right, row, column)], context);
    }
  for (size_t row = 0; row < destination->rows; ++row)
    for (size_t column = 0; column < destination->columns; ++column)
      destination->data[matrix_index(destination, row, column)] = temporary[row * destination->columns + column];
  free(temporary);
  return 0;
}

int chem_matrix_transform_columns(
    chem_matrix *destination, size_t destination_column,
    const chem_matrix *left, size_t left_column,
    const chem_matrix *right, size_t right_column,
    chem_matrix_binary_function function, void *context) {
  if (!matrix_valid(destination) || !matrix_valid(left) || !matrix_valid(right) || !function ||
      destination->rows != left->rows || destination->rows != right->rows ||
      destination_column >= destination->columns || left_column >= left->columns || right_column >= right->columns ||
      destination->rows > SIZE_MAX / sizeof(double)) return -1;
  double *temporary = destination->rows ? malloc(destination->rows * sizeof(*temporary)) : NULL;
  if (destination->rows && !temporary) return -2;
  for (size_t row = 0; row < destination->rows; ++row)
    temporary[row] = function(
        left->data[matrix_index(left, row, left_column)], right->data[matrix_index(right, row, right_column)], context);
  for (size_t row = 0; row < destination->rows; ++row)
    destination->data[matrix_index(destination, row, destination_column)] = temporary[row];
  free(temporary);
  return 0;
}

int chem_matrix_transform_column_vector(
    chem_matrix *destination, size_t destination_column,
    const chem_matrix *matrix, size_t matrix_column,
    const double *vector, size_t vector_length,
    chem_matrix_binary_function function, void *context) {
  if (!matrix_valid(destination) || !matrix_valid(matrix) || !function ||
      destination->rows != matrix->rows || destination_column >= destination->columns ||
      matrix_column >= matrix->columns || destination->rows != vector_length ||
      (vector_length && !vector) || destination->rows > SIZE_MAX / sizeof(double)) return -1;
  double *temporary = destination->rows ? malloc(destination->rows * sizeof(*temporary)) : NULL;
  if (destination->rows && !temporary) return -2;
  for (size_t row = 0; row < destination->rows; ++row)
    temporary[row] = function(matrix->data[matrix_index(matrix, row, matrix_column)], vector[row], context);
  for (size_t row = 0; row < destination->rows; ++row)
    destination->data[matrix_index(destination, row, destination_column)] = temporary[row];
  free(temporary);
  return 0;
}

int chem_matrix_print(const chem_matrix *matrix, FILE *stream) {
  if (!matrix_valid(matrix) || !stream) return -1;
  for (size_t row = 0; row < matrix->rows; ++row) {
    for (size_t column = 0; column < matrix->columns; ++column) {
      if (column && fputc(',', stream) == EOF) return -1;
      if (fprintf(stream, "%g", matrix->data[matrix_index(matrix, row, column)]) < 0) return -1;
    }
    if (fputc('\n', stream) == EOF) return -1;
  }
  return 0;
}
