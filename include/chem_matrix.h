/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_MATRIX_H
#define CHEM_MATRIX_H
#include "chem.h"
#include <stdio.h>

typedef enum { CHEM_MATRIX_REAL, CHEM_MATRIX_INT } ChemMatrixType;
/* group == 0 selects row-major storage; group > 0 selects groups of
 * contiguous column vectors. Padding belongs to the allocation, not the shape.
 * Initialize owning matrices with {0}; free before reusing them as outputs. */
typedef struct {
    size_t rows, columns, group, size;
    ChemMatrixType type;
    void *data;
} ChemMatrix;
ChemStatus chem_matrix_create(size_t rows, size_t columns, size_t group,
                              ChemMatrixType type, double initial, ChemMatrix *out);
void chem_matrix_free(ChemMatrix *);
ChemStatus chem_matrix_copy(const ChemMatrix *, ChemMatrix *out);
ChemStatus chem_matrix_from_rows(size_t rows, const size_t *lengths,
    const void *const *values, size_t group, ChemMatrixType, ChemMatrix *out);
size_t chem_matrix_index(const ChemMatrix *, size_t row, size_t column);
size_t chem_matrix_row_stride(const ChemMatrix *);
size_t chem_matrix_column_stride(const ChemMatrix *);
size_t chem_matrix_group_count(const ChemMatrix *);
size_t chem_matrix_group_size(const ChemMatrix *);
void *chem_matrix_at(ChemMatrix *, size_t row, size_t column);
const void *chem_matrix_const_at(const ChemMatrix *, size_t row, size_t column);
double chem_matrix_get(const ChemMatrix *, size_t row, size_t column);
ChemStatus chem_matrix_set(ChemMatrix *, size_t row, size_t column, double);
ChemStatus chem_matrix_copy_row(const ChemMatrix *, size_t row, void *out, size_t count);
ChemStatus chem_matrix_assign_row(ChemMatrix *, size_t row, const void *, size_t count);
ChemStatus chem_matrix_fill(ChemMatrix *, double);
ChemStatus chem_matrix_max(ChemMatrix *, double);
ChemStatus chem_matrix_min(ChemMatrix *, double);
ChemStatus chem_matrix_axpy(ChemMatrix *y, double alpha, const ChemMatrix *x);
/* Each callback returns the new destination element. Input matrices must have
 * the same logical shape; their layouts may differ. */
typedef double (*ChemMatrixElementFunction)(double destination,
    const double *inputs, size_t count, void *context);
ChemStatus chem_matrix_for_each(ChemMatrix *, size_t input_count,
    const ChemMatrix *const *inputs, ChemMatrixElementFunction, void *context);
ChemStatus chem_matrix_print(const ChemMatrix *, FILE *);

/* Borrowed arrays and views never own their storage. Const constructors retain
 * only a const pointer; all write operations reject such views. */
typedef enum { CHEM_ARRAY_MATRIX, CHEM_ARRAY_VECTOR, CHEM_ARRAY_SPARSE } ChemArrayKind;
typedef struct {
    const void *data;
    void *writable;
    size_t rows, columns, group;
    ChemMatrixType type;
    ChemArrayKind kind;
    int read_only;
} ChemArray;
typedef enum {
    CHEM_VIEW_VECTOR, CHEM_VIEW_DENSE, CHEM_VIEW_SPARSE,
    CHEM_VIEW_GROUPED_DENSE, CHEM_VIEW_GROUPED_SPARSE
} ChemViewKind;
typedef struct {
    const void *data;
    void *writable;
    size_t length, padded_length, group, block_stride, stride;
    ChemMatrixType type;
    ChemViewKind kind;
} ChemView;
ChemArray chem_array_matrix(ChemMatrix *);
ChemArray chem_array_const_matrix(const ChemMatrix *);
ChemArray chem_array_vector(void *, size_t length, ChemMatrixType);
ChemArray chem_array_const_vector(const void *, size_t length, ChemMatrixType);
ChemStatus chem_array_column(ChemArray, size_t column, ChemView *);
ChemView chem_view_const(ChemView);
/* Select one group and precompute its base pointer. Logical length excludes
 * tail padding. Grouped fill/copy and non-strict row loops include padding. */
ChemStatus chem_view_group(ChemView raw, size_t group_index, ChemView *out);
double chem_view_get(ChemView, size_t row);
ChemStatus chem_view_set(ChemView, size_t row, double);
ChemStatus chem_view_fill(ChemView, double);
ChemStatus chem_view_copy(ChemView destination, ChemView source);
/* The row callback sees a snapshot of all arguments. Its updated values are
 * written to mutable views; const arguments are never written. Logical rows
 * only are visited for raw views, so external vectors need no padding. When
 * all arguments are grouped views, the default loop includes their padding. */
typedef void (*ChemRowFunction)(double *values, size_t count, void *context);
ChemStatus chem_views_for_each(size_t count, const ChemView *, ChemRowFunction, void *context);
ChemStatus chem_views_for_each_strict(size_t count, const ChemView *, ChemRowFunction, void *context);
typedef double (*ChemReduceFunction)(double value, double accumulator, void *context);
ChemStatus chem_view_reduce(ChemView, int strict, ChemReduceFunction, void *context, double *accumulator);
/* Function creation fixes argument kinds, scalar types and mutability. Dense
 * functions fix columns; functions with sparse arguments allow different
 * structures but require matching matrix group widths. Row counts can change
 * and must agree only at invocation time. Empty calls skip the callback. Context is
 * borrowed, providing the C equivalent of captured lambda state. */
typedef ChemStatus (*ChemArrayCallback)(const ChemArray *, size_t count, void *context);
typedef struct ChemArrayFunction ChemArrayFunction;
ChemStatus chem_array_function_create(size_t count, const ChemArray *prototypes,
    ChemArrayCallback, void *context, ChemArrayFunction **out);
ChemStatus chem_array_function_call(const ChemArrayFunction *, size_t count, const ChemArray *);
void chem_array_function_free(ChemArrayFunction *);
#endif
