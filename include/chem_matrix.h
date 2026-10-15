#ifndef CHEM_MATRIX_H
#define CHEM_MATRIX_H
#include "chem.h"
#include <stdio.h>

/* width=0: ordinary row-major. width>0: groups of width rows, then columns,
   then contiguous lanes. All allocation sizes and indices are checked. */
typedef struct
{
  size_t rows, columns, width, padded_rows, size, element_size;
  void *data;
} ck_matrix;
ck_status ck_matrix_create(size_t rows, size_t columns, size_t width, size_t element_size,
                           const void *initial, ck_matrix *matrix);
void ck_matrix_free(ck_matrix *matrix);
ck_status ck_matrix_clone(const ck_matrix *source, ck_matrix *destination);
size_t ck_matrix_index(const ck_matrix *matrix, size_t row, size_t column);
void *ck_matrix_at(ck_matrix *matrix, size_t row, size_t column);
const void *ck_matrix_const_at(const ck_matrix *matrix, size_t row, size_t column);
ck_status ck_matrix_set_row(ck_matrix *matrix, size_t row, const void *values, size_t count);
ck_status ck_matrix_get_row(const ck_matrix *matrix, size_t row, void *values, size_t count);
ck_status ck_matrix_from_rows(size_t rows, const void *const *values, const size_t *lengths, size_t width,
                              size_t element_size, ck_matrix *matrix);
void ck_matrix_fill(ck_matrix *matrix, const void *value);
ck_status ck_matrix_axpy(ck_matrix *y, double alpha, const ck_matrix *x);
void ck_matrix_bound(ck_matrix *matrix, double value, int lower_bound);
ck_status ck_matrix_print(const ck_matrix *matrix, FILE *output);

typedef enum
{
  CK_SPARSE_ROW,
  CK_SPARSE_COLUMN
} ck_sparse_order;
typedef struct
{
  size_t order, nnz;
  ck_sparse_order ordering;
  size_t *map, *rows, *columns;
  ck_matrix values;
} ck_sparse;
ck_status ck_sparse_create(size_t order, size_t blocks, size_t width, ck_sparse_order ordering,
                           const unsigned char *pattern, size_t element_size, const void *initial,
                           ck_sparse *matrix);
ck_status ck_sparse_from_entries(size_t order, size_t blocks, size_t width, ck_sparse_order ordering,
                                 const size_t *rows, const size_t *columns, size_t count, size_t element_size,
                                 const void *initial, ck_sparse *matrix);
void ck_sparse_free(ck_sparse *matrix);
size_t ck_sparse_index(const ck_sparse *matrix, size_t block, size_t row, size_t column);
void *ck_sparse_at(ck_sparse *matrix, size_t block, size_t row, size_t column);
const void *ck_sparse_const_at(const ck_sparse *matrix, size_t block, size_t row, size_t column);
ck_status ck_sparse_diagonal(ck_sparse *matrix, double alpha);
ck_status ck_sparse_print(const ck_sparse *matrix, FILE *output, int nonzero_only);

/* A view represents one matrix column/block element, or a caller-owned vector.
   strict operations visit exactly rows; non-strict operations may visit padding.
   writable=0 enforces const access at runtime, replacing C++ signature checks. */
typedef struct
{
  void *data;
  size_t rows, capacity, element_size, width, group_stride, lane_stride, offset;
  int writable;
} ck_view;
ck_status ck_matrix_column(const ck_matrix *matrix, size_t column, int writable, ck_view *view);
ck_status ck_sparse_block_view(const ck_sparse *matrix, size_t flat_index, int writable, ck_view *view);
ck_status ck_vector_view(void *data, size_t rows, size_t element_size, int writable, ck_view *view);
const void *ck_view_const_at(const ck_view *view, size_t row);
void *ck_view_at(ck_view *view, size_t row);
/* Precompute one contiguous group, including padding in its capacity. */
ck_status ck_view_group(ck_view source, size_t group, ck_view *view);
ck_status ck_view_fill(ck_view destination, const void *value, int strict);
ck_status ck_view_copy(ck_view destination, ck_view source, int strict);
/* Mutable slots must be obtained through ck_view_at; the callback returns a
   status so an error can abort iteration without exceptions or longjmp. */
typedef ck_status (*ck_view_callback)(ck_view *views, size_t count, size_t row, void *context);
ck_status ck_view_apply(ck_view *views, size_t count, int strict, ck_view_callback function, void *context);
typedef ck_status (*ck_reduce_callback)(const void *value, void *accumulator, void *context);
ck_status ck_view_reduce(ck_view view, int strict, ck_reduce_callback function, void *accumulator,
                         void *context);

typedef enum
{
  CK_ARGUMENT_MATRIX,
  CK_ARGUMENT_SPARSE,
  CK_ARGUMENT_VECTOR
} ck_argument_kind;
typedef struct
{
  ck_argument_kind kind;
  const ck_matrix *matrix;
  const ck_sparse *sparse;
  ck_view vector;
  int writable;
} ck_argument;
ck_argument ck_matrix_argument(const ck_matrix *matrix, int writable);
ck_argument ck_sparse_argument(const ck_sparse *matrix, int writable);
ck_argument ck_vector_argument(ck_view vector);
typedef ck_status (*ck_array_callback)(const ck_argument *arguments, size_t count, void *context);
typedef struct
{
  size_t count;
  ck_argument *prototypes;
  ck_array_callback callback;
  void *context;
} ck_array_function;
ck_status ck_array_function_create(const ck_argument *prototypes, size_t count, ck_array_callback callback,
                                   void *context, ck_array_function *function);
ck_status ck_array_function_invoke(const ck_array_function *function, const ck_argument *arguments,
                                   size_t count);
void ck_array_function_free(ck_array_function *function);
#endif
