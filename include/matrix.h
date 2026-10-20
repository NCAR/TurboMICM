#ifndef CHEM_MATRIX_H
#define CHEM_MATRIX_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
  double *data;
  size_t rows;
  size_t columns;
  size_t row_stride;
  size_t column_stride;
} chem_matrix_view;

typedef struct {
  const double *data;
  size_t rows;
  size_t columns;
  size_t row_stride;
  size_t column_stride;
} chem_const_matrix_view;

typedef struct {
  double *data;
  size_t length;
  size_t stride;
} chem_vector_view;

typedef struct {
  const double *data;
  size_t length;
  size_t stride;
} chem_const_vector_view;

typedef double (*chem_unary_function)(double value, void *context);
typedef double (*chem_binary_function)(double left, double right, void *context);

typedef enum {
  CHEM_SPARSE_CSR = 0,
  CHEM_SPARSE_CSC = 1
} chem_sparse_format;

typedef struct {
  size_t rows;
  size_t columns;
  size_t nonzeros;
  size_t blocks;
  size_t lanes;
  size_t groups;
  chem_sparse_format format;
  size_t *offsets;
  uint32_t *indices;
  double *values;
} chem_sparse_matrix;

int chem_matrix_view_init(chem_matrix_view *view, double *data, size_t rows,
                          size_t columns, size_t row_stride,
                          size_t column_stride);
int chem_const_matrix_view_init(chem_const_matrix_view *view, const double *data,
                                size_t rows, size_t columns,
                                size_t row_stride, size_t column_stride);
double *chem_matrix_at(chem_matrix_view *view, size_t row, size_t column);
const double *chem_const_matrix_at(const chem_const_matrix_view *view,
                                   size_t row, size_t column);
int chem_matrix_fill(chem_matrix_view view, double value);
int chem_matrix_copy(chem_matrix_view destination,
                     chem_const_matrix_view source);
int chem_matrix_axpy(chem_matrix_view destination, double alpha,
                     chem_const_matrix_view source);
int chem_matrix_transform(chem_matrix_view destination,
                          chem_const_matrix_view source,
                          chem_unary_function function, void *context);
int chem_matrix_transform2(chem_matrix_view destination,
                           chem_const_matrix_view left,
                           chem_const_matrix_view right,
                           chem_binary_function function, void *context);
int chem_matrix_sum(chem_const_matrix_view view, double *value);
int chem_matrix_min(chem_const_matrix_view view, double *value);
int chem_matrix_max(chem_const_matrix_view view, double *value);

int chem_vector_view_init(chem_vector_view *view, double *data, size_t length,
                          size_t stride);
int chem_const_vector_view_init(chem_const_vector_view *view,
                                const double *data, size_t length,
                                size_t stride);
int chem_vector_copy(chem_vector_view destination,
                     chem_const_vector_view source);
int chem_vector_fill(chem_vector_view destination, double value);
int chem_vector_axpy(chem_vector_view destination, double alpha,
                     chem_const_vector_view source);

int chem_sparse_matrix_init(chem_sparse_matrix *matrix, size_t rows,
                            size_t columns, const uint32_t *entry_rows,
                            const uint32_t *entry_columns, size_t nonzeros,
                            size_t blocks, size_t lanes,
                            chem_sparse_format format);
void chem_sparse_matrix_destroy(chem_sparse_matrix *matrix);
int chem_sparse_matrix_contains(const chem_sparse_matrix *matrix, size_t row,
                                size_t column);
int chem_sparse_matrix_get(const chem_sparse_matrix *matrix, size_t block,
                           size_t row, size_t column, double *value);
int chem_sparse_matrix_set(chem_sparse_matrix *matrix, size_t block,
                           size_t row, size_t column, double value);
int chem_sparse_matrix_fill(chem_sparse_matrix *matrix, double value);
int chem_sparse_matrix_add_diagonal(chem_sparse_matrix *matrix, double value);
int chem_sparse_matrix_copy(chem_sparse_matrix *destination,
                            const chem_sparse_matrix *source);
int chem_sparse_matrix_axpy(chem_sparse_matrix *destination, double alpha,
                            const chem_sparse_matrix *source);
int chem_sparse_matrix_transform(chem_sparse_matrix *matrix,
                                 chem_unary_function function, void *context);
int chem_sparse_matrix_validate_blocks(const chem_sparse_matrix *matrix,
                                       size_t block_count);

int chem_lu_decompose_dense(const double *matrix, size_t order,
                            double *lower, double *upper);
int chem_lu_factor_dense(double *matrix, size_t order);
int chem_lu_solve_dense(const double *factor, size_t order,
                        double *right_hand_side);

#endif
