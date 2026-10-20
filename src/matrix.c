#include "chemistry.h"

#include <math.h>

static int valid_view(const chem_matrix_view *view) {
  if (view == NULL) return 0;
  if (view->rows == 0 || view->columns == 0) return 1;
  if (view->data == NULL || view->row_stride == 0 || view->column_stride == 0)
    return 0;
  if (view->rows - 1 > (size_t)-1 / view->row_stride ||
      view->columns - 1 > (size_t)-1 / view->column_stride) return 0;
  const size_t row_offset = (view->rows - 1) * view->row_stride;
  const size_t column_offset = (view->columns - 1) * view->column_stride;
  return row_offset <= (size_t)-1 - column_offset;
}

static int valid_const_view(const chem_const_matrix_view *view) {
  if (view == NULL) return 0;
  if (view->rows == 0 || view->columns == 0) return 1;
  if (view->data == NULL || view->row_stride == 0 || view->column_stride == 0)
    return 0;
  if (view->rows - 1 > (size_t)-1 / view->row_stride ||
      view->columns - 1 > (size_t)-1 / view->column_stride) return 0;
  const size_t row_offset = (view->rows - 1) * view->row_stride;
  const size_t column_offset = (view->columns - 1) * view->column_stride;
  return row_offset <= (size_t)-1 - column_offset;
}

int chem_matrix_view_init(chem_matrix_view *view, double *data,
                          size_t rows, size_t columns,
                          size_t row_stride, size_t column_stride) {
  if (view == NULL || (rows && columns && data == NULL) ||
      (rows && columns && (!row_stride || !column_stride)) ||
      (rows && columns && rows - 1 > (size_t)-1 / row_stride) ||
      (rows && columns && columns - 1 > (size_t)-1 / column_stride) ||
      (rows && columns && (rows - 1) * row_stride >
          (size_t)-1 - (columns - 1) * column_stride)) return 1;
  *view = (chem_matrix_view){data, rows, columns, row_stride, column_stride};
  return 0;
}

int chem_const_matrix_view_init(chem_const_matrix_view *view,
                                const double *data, size_t rows,
                                size_t columns, size_t row_stride,
                                size_t column_stride) {
  if (view == NULL || (rows && columns && data == NULL) ||
      (rows && columns && (!row_stride || !column_stride)) ||
      (rows && columns && rows - 1 > (size_t)-1 / row_stride) ||
      (rows && columns && columns - 1 > (size_t)-1 / column_stride) ||
      (rows && columns && (rows - 1) * row_stride >
          (size_t)-1 - (columns - 1) * column_stride)) return 1;
  *view = (chem_const_matrix_view){data, rows, columns, row_stride, column_stride};
  return 0;
}

int chem_const_matrix_view_get(const chem_const_matrix_view *view,
                               size_t row, size_t column, double *value) {
  if (!valid_const_view(view) || value == NULL || row >= view->rows ||
      column >= view->columns) return 1;
  *value = view->data[row * view->row_stride + column * view->column_stride];
  return 0;
}

int chem_const_matrix_view_copy_row_to_vector(
    const chem_const_matrix_view *view, size_t row, double *vector,
    size_t vector_size) {
  if (!valid_const_view(view) || row >= view->rows ||
      vector_size != view->columns || (vector_size && vector == NULL)) return 1;
  for (size_t j = 0; j < view->columns; ++j)
    vector[j] = view->data[row*view->row_stride + j*view->column_stride];
  return 0;
}

int chem_matrix_view_get(const chem_matrix_view *view, size_t row,
                         size_t column, double *value) {
  if (!valid_view(view) || value == NULL || row >= view->rows ||
      column >= view->columns) return 1;
  *value = view->data[row * view->row_stride + column * view->column_stride];
  return 0;
}

int chem_matrix_view_set(chem_matrix_view *view, size_t row,
                         size_t column, double value) {
  if (!valid_view(view) || row >= view->rows || column >= view->columns) return 1;
  view->data[row * view->row_stride + column * view->column_stride] = value;
  return 0;
}

int chem_matrix_view_fill(chem_matrix_view *view, double value) {
  if (!valid_view(view)) return 1;
  for (size_t i = 0; i < view->rows; ++i)
    for (size_t j = 0; j < view->columns; ++j)
      view->data[i * view->row_stride + j * view->column_stride] = value;
  return 0;
}

int chem_matrix_view_for_each(chem_matrix_view *view,
                              chem_matrix_element_function function,
                              void *context) {
  if (!valid_view(view) || function == NULL) return 1;
  for (size_t i = 0; i < view->rows; ++i)
    for (size_t j = 0; j < view->columns; ++j)
      function(i, j, &view->data[i * view->row_stride + j * view->column_stride],
               context);
  return 0;
}

static int vector_dimensions(const chem_matrix_view *view, size_t vector_size) {
  return valid_view(view) && view->rows <= (size_t)-1 / (view->columns ? view->columns : 1) &&
      view->rows * view->columns == vector_size;
}

int chem_matrix_view_copy_to_vector(const chem_matrix_view *view,
                                    double *vector, size_t vector_size) {
  if (!vector_dimensions(view, vector_size) || (vector_size && vector == NULL)) return 1;
  size_t k = 0;
  for (size_t i = 0; i < view->rows; ++i)
    for (size_t j = 0; j < view->columns; ++j)
      vector[k++] = view->data[i * view->row_stride + j * view->column_stride];
  return 0;
}

int chem_matrix_view_copy_from_vector(chem_matrix_view *view,
                                      const double *vector,
                                      size_t vector_size) {
  if (!vector_dimensions(view, vector_size) || (vector_size && vector == NULL)) return 1;
  size_t k = 0;
  for (size_t i = 0; i < view->rows; ++i)
    for (size_t j = 0; j < view->columns; ++j)
      view->data[i * view->row_stride + j * view->column_stride] = vector[k++];
  return 0;
}

int chem_matrix_view_copy_row_to_vector(const chem_matrix_view *view,
                                        size_t row, double *vector,
                                        size_t vector_size) {
  if (!valid_view(view) || row >= view->rows || vector_size != view->columns ||
      (vector_size && vector == NULL)) return 1;
  for (size_t j = 0; j < view->columns; ++j)
    vector[j] = view->data[row*view->row_stride + j*view->column_stride];
  return 0;
}

int chem_matrix_view_copy_row_from_vector(chem_matrix_view *view,
                                          size_t row, const double *vector,
                                          size_t vector_size) {
  if (!valid_view(view) || row >= view->rows || vector_size != view->columns ||
      (vector_size && vector == NULL)) return 1;
  for (size_t j = 0; j < view->columns; ++j)
    view->data[row*view->row_stride + j*view->column_stride] = vector[j];
  return 0;
}

int chem_matrix_view_copy_column_to_vector(const chem_matrix_view *view,
                                           size_t column, double *vector,
                                           size_t vector_size) {
  if (!valid_view(view) || column >= view->columns || vector_size != view->rows ||
      (vector_size && vector == NULL)) return 1;
  for (size_t i = 0; i < view->rows; ++i)
    vector[i] = view->data[i*view->row_stride + column*view->column_stride];
  return 0;
}

int chem_matrix_view_copy_column_from_vector(chem_matrix_view *view,
                                             size_t column,
                                             const double *vector,
                                             size_t vector_size) {
  if (!valid_view(view) || column >= view->columns || vector_size != view->rows ||
      (vector_size && vector == NULL)) return 1;
  for (size_t i = 0; i < view->rows; ++i)
    view->data[i*view->row_stride + column*view->column_stride] = vector[i];
  return 0;
}

int chem_matrix_view_axpy(chem_matrix_view *destination,
                          const chem_matrix_view *source, double alpha) {
  if (!valid_view(destination) || !valid_view(source) ||
      destination->rows != source->rows ||
      destination->columns != source->columns || !isfinite(alpha)) return 1;
  for (size_t i = 0; i < destination->rows; ++i)
    for (size_t j = 0; j < destination->columns; ++j)
      destination->data[i * destination->row_stride + j * destination->column_stride] +=
          alpha * source->data[i * source->row_stride + j * source->column_stride];
  return 0;
}

static int matrix_reduce(const chem_matrix_view *view, double *value, int mode) {
  if (!valid_view(view) || value == NULL || view->rows == 0 || view->columns == 0) return 1;
  double result = view->data[0];
  for (size_t i = 0; i < view->rows; ++i)
    for (size_t j = 0; j < view->columns; ++j) {
      const double x = view->data[i * view->row_stride + j * view->column_stride];
      if (mode == 0) result += (i || j) ? x : 0.0;
      else if (mode == 1 && x < result) result = x;
      else if (mode == 2 && x > result) result = x;
    }
  *value = result;
  return 0;
}

int chem_matrix_view_sum(const chem_matrix_view *view, double *value) {
  return matrix_reduce(view, value, 0);
}

int chem_matrix_view_min(const chem_matrix_view *view, double *value) {
  return matrix_reduce(view, value, 1);
}

int chem_matrix_view_max(const chem_matrix_view *view, double *value) {
  return matrix_reduce(view, value, 2);
}
