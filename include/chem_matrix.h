#ifndef CHEM_MATRIX_H
#define CHEM_MATRIX_H

#include "chem_views.h"

#include <stddef.h>
#include <stdio.h>

typedef enum {
  CHEM_MATRIX_ROW_MAJOR = 0,
  CHEM_MATRIX_GROUPED = 1
} chem_matrix_layout;

typedef struct {
  size_t rows;
  size_t columns;
  size_t lanes;
  size_t storage_length;
  chem_matrix_layout layout;
  double *data;
} chem_matrix;

typedef double (*chem_matrix_unary_function)(double value, void *context);
typedef double (*chem_matrix_binary_function)(double left, double right, void *context);
typedef void (*chem_matrix_reducer)(double value, double *accumulator, void *context);

/* Matrices own their storage. Initialize an empty instance with {0}; init
   returns -1 for invalid dimensions/layout and -2 on allocation failure.
   Use lanes=1 for row-major storage; grouped storage keeps each lane group
   contiguous by column, padding only its final group. Empty row or column
   dimensions are valid. Destroy releases storage and clears the struct. */
int chem_matrix_init(
    chem_matrix *matrix, size_t rows, size_t columns, size_t lanes,
    chem_matrix_layout layout, double initial_value);
void chem_matrix_destroy(chem_matrix *matrix);
int chem_matrix_get(const chem_matrix *matrix, size_t row, size_t column, double *value);
int chem_matrix_set(chem_matrix *matrix, size_t row, size_t column, double value);
int chem_matrix_row_view(chem_matrix *matrix, size_t row, chem_group_view *view);
int chem_matrix_const_row_view(const chem_matrix *matrix, size_t row, chem_const_group_view *view);
int chem_matrix_column_group_view(chem_matrix *matrix, size_t group, size_t column, chem_group_view *view);
int chem_matrix_const_column_group_view(
    const chem_matrix *matrix, size_t group, size_t column, chem_const_group_view *view);

int chem_matrix_fill(chem_matrix *matrix, double value);
int chem_matrix_copy(chem_matrix *destination, const chem_matrix *source);
int chem_matrix_assign_row(chem_matrix *matrix, size_t row, const double *values, size_t value_count);
int chem_matrix_read_row(const chem_matrix *matrix, size_t row, double *values, size_t value_count);
int chem_matrix_axpy(chem_matrix *destination, double alpha, const chem_matrix *source);
int chem_matrix_max_scalar(chem_matrix *matrix, double lower_bound);
int chem_matrix_min_scalar(chem_matrix *matrix, double upper_bound);
int chem_matrix_min_max(const chem_matrix *matrix, double *minimum, double *maximum);
int chem_matrix_reduce(
    const chem_matrix *matrix, double *accumulator, chem_matrix_reducer reducer, void *context);
int chem_matrix_transform_unary(
    chem_matrix *destination, const chem_matrix *source,
    chem_matrix_unary_function function, void *context);
int chem_matrix_transform_binary(
    chem_matrix *destination, const chem_matrix *left, const chem_matrix *right,
    chem_matrix_binary_function function, void *context);
int chem_matrix_transform_columns(
    chem_matrix *destination, size_t destination_column,
    const chem_matrix *left, size_t left_column,
    const chem_matrix *right, size_t right_column,
    chem_matrix_binary_function function, void *context);
int chem_matrix_transform_column_vector(
    chem_matrix *destination, size_t destination_column,
    const chem_matrix *matrix, size_t matrix_column,
    const double *vector, size_t vector_length,
    chem_matrix_binary_function function, void *context);
int chem_matrix_print(const chem_matrix *matrix, FILE *stream);

#endif
