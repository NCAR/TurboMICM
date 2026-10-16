/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_ARRAY_H
#define CHEM_ARRAY_H

#include "chem.h"

/* Dense rows and sparse blocks share an array of independent records.
 * group=0 gives ordinary row/block-major storage; group>0 interleaves that
 * many records. Sparse order is 0 (rows) or 1 (columns). */
enum { CH_ARRAY_REAL, CH_ARRAY_INT };
enum { CH_ARRAY_DENSE, CH_ARRAY_SPARSE, CH_ARRAY_VECTOR };
enum { CH_VIEW_COLUMN, CH_VIEW_BLOCK, CH_VIEW_VECTOR, CH_VIEW_TEMPORARY };
typedef struct ch_array ch_array;
int ch_pattern_set(size_t dimension, unsigned char *pattern, size_t row, size_t column);
ch_array *ch_array_dense(size_t rows, size_t columns, size_t group, int type, double initial);
ch_array *ch_array_sparse(size_t blocks, size_t dimension, const unsigned char *pattern,
                          int order, size_t group, int type, double initial);
ch_array *ch_array_vector(size_t size, int type, double initial);
ch_array *ch_array_clone(const ch_array *);
void ch_array_free(ch_array *);
size_t ch_array_records(const ch_array *);
size_t ch_array_elements(const ch_array *);
size_t ch_array_storage_size(const ch_array *);
size_t ch_array_group_size(const ch_array *);
size_t ch_array_groups(const ch_array *);
size_t ch_array_row_stride(const ch_array *);
size_t ch_array_column_stride(const ch_array *);
void *ch_array_data(ch_array *);
int ch_array_element(const ch_array *, size_t row, size_t column);
int ch_array_index(const ch_array *, size_t record, size_t element, size_t *index);
int ch_array_get(const ch_array *, size_t record, size_t element, double *value);
int ch_array_set(ch_array *, size_t record, size_t element, double value);
int ch_array_set_row(ch_array *, size_t row, const double *, size_t count);
int ch_array_fill(ch_array *, double);
int ch_array_copy(ch_array *, const ch_array *);
int ch_array_axpy(ch_array *, double alpha, const ch_array *);
int ch_array_max(ch_array *, double);
int ch_array_min(ch_array *, double);
int ch_array_add_diagonal(ch_array *, double);
int ch_array_print(const ch_array *, FILE *);
/* Return bytes required (excluding NUL). Output is always NUL terminated when
 * capacity>0; a NULL destination queries the length. */
size_t ch_array_format(const ch_array *, char *destination, size_t capacity);
size_t ch_array_format_dense(const ch_array *, char *destination, size_t capacity);

typedef struct {
  void *data;
  size_t size, capacity, stride;
  int type, read_only, category;
} ch_view;
typedef struct {
  ch_array *array;
  size_t start, size, capacity;
  int read_only;
} ch_group;
int ch_array_group(ch_array *, size_t group, int read_only, ch_group *);
int ch_group_view(const ch_group *, size_t element, ch_view *);
int ch_group_sparse_view(const ch_group *, size_t row, size_t column, ch_view *);
ch_view ch_temporary_view(double *data, size_t size, size_t capacity);
double ch_view_get(const ch_view *, size_t);
int ch_view_set(ch_view *, size_t, double);
int ch_view_fill(ch_view *, double, int strict);
int ch_view_copy(ch_view *, const ch_view *, int strict);
/* values are gathered before invoking the callback, then writable views are
 * scattered in argument order. Callbacks can capture arbitrary caller context. */
typedef void (*ch_element_function)(double *values, size_t count, void *context);
int ch_views_map(size_t count, ch_view *, ch_element_function, void *context, int strict);
typedef void (*ch_reduce_function)(double value, double *accumulator, void *context);
int ch_view_reduce(const ch_view *, ch_reduce_function, double *accumulator, void *context, int strict);
int ch_array_map(size_t count, ch_array *const *, const unsigned char *read_only,
                  ch_element_function, void *context);

/* Reusable dispatch captures dimensions and storage policies, not array data.
 * Record counts may change between calls; all arguments must agree at invocation.
 * A callback receives a bounded group from each matrix/vector. */
typedef struct ch_array_function ch_array_function;
typedef int (*ch_group_function)(size_t count, ch_group *groups, void *context);
int ch_array_function_create(size_t count, ch_array *const *examples,
                             const unsigned char *read_only, ch_group_function,
                             void *context, ch_array_function **result);
int ch_array_function_call(const ch_array_function *, size_t count, ch_array *const *);
void ch_array_function_free(ch_array_function *);

#endif
