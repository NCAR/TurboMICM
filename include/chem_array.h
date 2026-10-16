#ifndef CHEM_ARRAY_H
#define CHEM_ARRAY_H
#include "chem.h"
#include <stdio.h>

/* A zero vector length selects row-major storage. A positive length selects
 * groups of column-contiguous lanes, with initialized padding in the last group.
 * Integer arrays retain native integer storage; numeric helpers convert to double. */
typedef enum { CHEM_DOUBLE, CHEM_INT, CHEM_INDEX } chem_element_type;
typedef enum { CHEM_DENSE_ARRAY, CHEM_VECTOR_ARRAY, CHEM_SPARSE_ARRAY } chem_array_kind;
typedef struct {
    size_t rows, columns, vector_length, size;
    chem_element_type type;
    chem_array_kind kind;
    void *data;
} chem_array;
size_t chem_element_size(chem_element_type type);
chem_status chem_array_create(chem_array *array, size_t rows, size_t columns,
                              size_t vector_length, chem_element_type type, double initial);
chem_status chem_vector_create(chem_array *array, size_t count, size_t vector_length,
                               chem_element_type type, double initial);
void chem_array_destroy(chem_array *array);
chem_status chem_array_clone(chem_array *destination, const chem_array *source);
chem_status chem_array_import(chem_array *array, size_t rows, const size_t *lengths,
                              const double *const *values, size_t vector_length);
size_t chem_array_index(const chem_array *array, size_t row, size_t column);
size_t chem_array_row_stride(const chem_array *array);
size_t chem_array_column_stride(const chem_array *array);
size_t chem_array_groups(const chem_array *array);
double chem_array_get(const chem_array *array, size_t row, size_t column);
chem_status chem_array_set(chem_array *array, size_t row, size_t column, double value);
chem_status chem_array_assign_row(chem_array *array, size_t row, const double *values, size_t count);
chem_status chem_array_export_row(const chem_array *array, size_t row, double *values, size_t count);
void chem_array_fill(chem_array *array, double value);
void chem_array_max(chem_array *array, double value);
void chem_array_min(chem_array *array, double value);
chem_status chem_array_axpy(chem_array *y, double alpha, const chem_array *x);
chem_status chem_array_print(FILE *file, const chem_array *array);

typedef enum {
    CHEM_COLUMN_VIEW, CHEM_VECTOR_VIEW, CHEM_TEMPORARY_VIEW, CHEM_ELEMENT_VIEW,
    CHEM_RAW_COLUMN_VIEW, CHEM_RAW_ELEMENT_VIEW
} chem_view_kind;
typedef struct {
    const void *data;
    void *writable;
    size_t stride, count, active;
    chem_element_type type;
    chem_view_kind kind;
    chem_status status;
} chem_span;
typedef struct {
    const chem_array *array;
    size_t index, start, width, active;
    int writable;
} chem_array_group;
chem_status chem_array_get_group(const chem_array *array, size_t index, int writable, chem_array_group *group);
chem_span chem_group_column(chem_array_group group, size_t column);
chem_span chem_group_const_column(chem_array_group group, size_t column);
chem_span chem_group_vector(chem_array_group group);
chem_span chem_span_temporary(double *storage, size_t width, size_t active);
chem_span chem_span_vector(const void *data, void *writable, size_t count, chem_element_type type);
chem_span chem_span_const(chem_span span);
double chem_span_get(chem_span span, size_t index);
chem_status chem_span_set(chem_span span, size_t index, double value);
chem_status chem_span_fill(chem_span destination, double value);
chem_status chem_span_copy(chem_span destination, chem_span source);
typedef double (*chem_map_function)(const double *values, size_t count, void *context);
chem_status chem_span_map(chem_span destination, size_t n_sources, const chem_span *sources,
                          chem_map_function function, void *context);
typedef double (*chem_reduce_function)(double value, double accumulator, void *context);
chem_status chem_span_reduce(chem_span source, int strict, double *accumulator,
                             chem_reduce_function function, void *context);

typedef struct { const chem_array *array; int writable; } chem_array_argument;
chem_array_argument chem_array_mutable(chem_array *array);
chem_array_argument chem_array_const(const chem_array *array);
typedef chem_status (*chem_group_function)(size_t count, const chem_array_group *groups, void *context);
typedef struct {
    size_t count;
    size_t *columns, *vector_lengths;
    chem_element_type *types;
    chem_array_kind *kinds;
    int *writable;
    chem_group_function function;
    void *context;
} chem_array_function;
/* Creation records types, layouts and column counts, never row counts.
 * Invocation checks all arguments before writing any of them. */
chem_status chem_array_function_create(chem_array_function *function, size_t count,
                                       const chem_array_argument *arguments,
                                       chem_group_function callback, void *context);
chem_status chem_array_function_apply(const chem_array_function *function, size_t count,
                                      const chem_array_argument *arguments);
void chem_array_function_destroy(chem_array_function *function);

/* Raw views keep a matrix and logical slot; binding resolves the group once. */
typedef struct { const chem_array *array; size_t slot; int writable; chem_view_kind kind; } chem_raw_view;
chem_raw_view chem_array_column(const chem_array *array, size_t column, int writable);
chem_span chem_bind_view(chem_raw_view view, chem_array_group group);
#endif
