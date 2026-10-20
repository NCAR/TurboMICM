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

int chem_vector_matrix_view_storage_size(size_t rows, size_t columns,
                                         size_t group_width,
                                         size_t *storage_size) {
  if (storage_size == NULL || group_width == 0 ||
      columns > (size_t)-1/group_width) return 1;
  const size_t group_size = columns*group_width;
  const size_t groups = rows/group_width + (rows%group_width != 0);
  if (group_size && groups > (size_t)-1/group_size) return 1;
  *storage_size = groups*group_size;
  return 0;
}

int chem_vector_matrix_view_init(chem_vector_matrix_view *view, double *data,
                                 size_t rows, size_t columns,
                                 size_t group_width) {
  size_t storage_size = 0;
  if (view == NULL || group_width == 0 ||
      chem_vector_matrix_view_storage_size(rows,columns,group_width,
                                           &storage_size) ||
      (storage_size && data == NULL)) return 1;
  *view = (chem_vector_matrix_view){data,rows,columns,group_width};
  return 0;
}

int chem_const_vector_matrix_view_init(chem_const_vector_matrix_view *view,
                                       const double *data, size_t rows,
                                       size_t columns, size_t group_width) {
  size_t storage_size = 0;
  if (view == NULL || group_width == 0 ||
      chem_vector_matrix_view_storage_size(rows,columns,group_width,
                                           &storage_size) ||
      (storage_size && data == NULL)) return 1;
  *view = (chem_const_vector_matrix_view){data,rows,columns,group_width};
  return 0;
}

size_t chem_vector_matrix_view_group_size(const chem_vector_matrix_view *view) {
  if (view == NULL || view->group_width == 0 ||
      view->columns > (size_t)-1/view->group_width) return 0;
  return view->columns*view->group_width;
}

size_t chem_vector_matrix_view_number_of_groups(
    const chem_vector_matrix_view *view) {
  if (view == NULL || view->group_width == 0) return 0;
  return view->rows/view->group_width + (view->rows%view->group_width != 0);
}

size_t chem_vector_matrix_view_row_stride(const chem_vector_matrix_view *view) {
  return view && view->group_width ? 1 : 0;
}

size_t chem_vector_matrix_view_column_stride(const chem_vector_matrix_view *view) {
  return view ? view->group_width : 0;
}

static size_t vector_matrix_offset(size_t row, size_t column,
                                   size_t columns, size_t group_width) {
  return (row/group_width)*(columns*group_width) + column*group_width +
      row%group_width;
}

int chem_vector_matrix_view_get(const chem_vector_matrix_view *view,
                                size_t row, size_t column, double *value) {
  if (view == NULL || view->data == NULL || value == NULL ||
      row >= view->rows || column >= view->columns || view->group_width == 0)
    return 1;
  *value = view->data[vector_matrix_offset(row,column,view->columns,
                                           view->group_width)];
  return 0;
}

int chem_vector_matrix_view_set(chem_vector_matrix_view *view,
                                size_t row, size_t column, double value) {
  if (view == NULL || view->data == NULL || row >= view->rows ||
      column >= view->columns || view->group_width == 0) return 1;
  view->data[vector_matrix_offset(row,column,view->columns,view->group_width)] = value;
  return 0;
}

int chem_vector_matrix_view_fill(chem_vector_matrix_view *view, double value) {
  size_t storage_size = 0;
  if (view == NULL || view->group_width == 0 ||
      chem_vector_matrix_view_storage_size(view->rows,view->columns,
                                           view->group_width,&storage_size) ||
      (storage_size && view->data == NULL)) return 1;
  for (size_t row = 0; row < view->rows; ++row)
    for (size_t column = 0; column < view->columns; ++column)
      view->data[vector_matrix_offset(row,column,view->columns,
                                      view->group_width)] = value;
  return 0;
}

int chem_const_vector_matrix_view_get(
    const chem_const_vector_matrix_view *view, size_t row, size_t column,
    double *value) {
  if (view == NULL || view->data == NULL || value == NULL ||
      row >= view->rows || column >= view->columns || view->group_width == 0)
    return 1;
  *value = view->data[vector_matrix_offset(row,column,view->columns,
                                           view->group_width)];
  return 0;
}

static int vector_matrix_view_valid(const chem_vector_matrix_view *view) {
  size_t size = 0;
  return view && !chem_vector_matrix_view_storage_size(
      view->rows,view->columns,view->group_width,&size) &&
      (!size || view->data != NULL);
}

static int const_vector_matrix_view_valid(
    const chem_const_vector_matrix_view *view) {
  size_t size = 0;
  return view && !chem_vector_matrix_view_storage_size(
      view->rows,view->columns,view->group_width,&size) &&
      (!size || view->data != NULL);
}

static int vector_matrix_vector_size(size_t rows, size_t columns,
                                     size_t vector_size) {
  return columns == 0 ? rows == 0 && vector_size == 0
      : rows <= (size_t)-1/columns && rows*columns == vector_size;
}

int chem_vector_matrix_view_copy_to_vector(const chem_vector_matrix_view *view,
                                           double *vector,
                                           size_t vector_size) {
  if (!vector_matrix_view_valid(view) ||
      !vector_matrix_vector_size(view->rows,view->columns,vector_size) ||
      (vector_size && vector == NULL)) return 1;
  size_t k = 0;
  for (size_t row = 0; row < view->rows; ++row)
    for (size_t column = 0; column < view->columns; ++column)
      vector[k++] = view->data[vector_matrix_offset(
          row,column,view->columns,view->group_width)];
  return 0;
}

int chem_vector_matrix_view_copy_from_vector(chem_vector_matrix_view *view,
                                             const double *vector,
                                             size_t vector_size) {
  if (!vector_matrix_view_valid(view) ||
      !vector_matrix_vector_size(view->rows,view->columns,vector_size) ||
      (vector_size && vector == NULL)) return 1;
  size_t k = 0;
  for (size_t row = 0; row < view->rows; ++row)
    for (size_t column = 0; column < view->columns; ++column)
      view->data[vector_matrix_offset(row,column,view->columns,
                                      view->group_width)] = vector[k++];
  return 0;
}

int chem_const_vector_matrix_view_copy_to_vector(
    const chem_const_vector_matrix_view *view, double *vector,
    size_t vector_size) {
  if (!const_vector_matrix_view_valid(view) ||
      !vector_matrix_vector_size(view->rows,view->columns,vector_size) ||
      (vector_size && vector == NULL)) return 1;
  size_t k = 0;
  for (size_t row = 0; row < view->rows; ++row)
    for (size_t column = 0; column < view->columns; ++column)
      vector[k++] = view->data[vector_matrix_offset(row,column,view->columns,
                                                    view->group_width)];
  return 0;
}

int chem_vector_matrix_view_axpy(chem_vector_matrix_view *destination,
                                 const chem_vector_matrix_view *source,
                                 double alpha) {
  if (!vector_matrix_view_valid(destination) ||
      !vector_matrix_view_valid(source) || destination->rows != source->rows ||
      destination->columns != source->columns || !isfinite(alpha)) return 1;
  for (size_t row = 0; row < destination->rows; ++row)
    for (size_t column = 0; column < destination->columns; ++column)
      destination->data[vector_matrix_offset(row,column,destination->columns,
                                             destination->group_width)] +=
          alpha*source->data[vector_matrix_offset(row,column,source->columns,
                                                  source->group_width)];
  return 0;
}

int chem_vector_matrix_view_for_each(chem_vector_matrix_view *view,
                                    chem_matrix_element_function function,
                                    void *context) {
  if (!vector_matrix_view_valid(view) || function == NULL) return 1;
  for (size_t row = 0; row < view->rows; ++row)
    for (size_t column = 0; column < view->columns; ++column)
      function(row,column,&view->data[vector_matrix_offset(
          row,column,view->columns,view->group_width)],context);
  return 0;
}

static int vector_matrix_reduce(const chem_vector_matrix_view *view,
                                double *value, int maximum) {
  if (!vector_matrix_view_valid(view) || value == NULL || view->rows == 0 ||
      view->columns == 0) return 1;
  double result = view->data[vector_matrix_offset(0,0,view->columns,
                                                   view->group_width)];
  for (size_t row = 0; row < view->rows; ++row)
    for (size_t column = 0; column < view->columns; ++column) {
      const double item = view->data[vector_matrix_offset(
          row,column,view->columns,view->group_width)];
      if (maximum ? item > result : item < result) result = item;
    }
  *value = result;
  return 0;
}

int chem_vector_matrix_view_min(const chem_vector_matrix_view *view,
                                double *value) {
  return vector_matrix_reduce(view,value,0);
}

int chem_vector_matrix_view_max(const chem_vector_matrix_view *view,
                                double *value) {
  return vector_matrix_reduce(view,value,1);
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
