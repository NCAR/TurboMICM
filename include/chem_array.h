#ifndef CHEM_ARRAY_H
#define CHEM_ARRAY_H
#include "chem.h"

/* Rows are cells (blocks for sparse arrays); columns are dense columns or sparse elements.
 * A group interleaves adjacent cells. A group of one gives ordinary row-major storage. */
typedef struct {
  size_t rows, columns, group, allocated;
  double *values;
  size_t order;
  int *index;
  int column_order;
} ChemArray;
typedef struct {
  ChemArray *array;
  double *vector;
  size_t column, start, count;
  int readonly;
} ChemView;
typedef void (*ChemElementFunction)(double *values, size_t count, void *context);
typedef void (*ChemReduceFunction)(double value, double *accumulator, void *context);
typedef struct { ChemArray *array; double *vector; size_t count; int readonly; } ChemArgument;
typedef void (*ChemArrayFunction)(ChemArgument *arguments, size_t count, void *context);
typedef struct ChemFunction ChemFunction;

ChemArray *chem_array_dense(size_t rows, size_t columns, size_t group, double initial);
ChemStatus chem_array_from_rows(const double *const *rows, const size_t *lengths, size_t count, size_t group, ChemArray **result);
ChemArray *chem_array_sparse(size_t blocks, size_t order, size_t group, int column_order,
                            const size_t *row_ids, const size_t *column_ids, size_t elements, double initial);
void chem_array_free(ChemArray *);
size_t chem_array_index(const ChemArray *, size_t row, size_t column);
size_t chem_sparse_index(const ChemArray *, size_t block, size_t row, size_t column);
ChemStatus chem_array_get(const ChemArray *, size_t row, size_t column, double *value);
ChemStatus chem_array_set(ChemArray *, size_t row, size_t column, double value);
ChemStatus chem_array_row(ChemArray *, size_t row, double *values, size_t count, int write);
void chem_array_fill(ChemArray *, double value);
void chem_array_bound(ChemArray *, double value, int maximum);
ChemStatus chem_array_copy(ChemArray *destination, const ChemArray *source);
ChemStatus chem_array_axpy(ChemArray *destination, double alpha, const ChemArray *source);
void chem_array_add_diagonal(ChemArray *, double value);
size_t chem_array_format(const ChemArray *, int nonzero_only, char *buffer, size_t capacity);
ChemView chem_array_view(ChemArray *, size_t column, int readonly);
ChemView chem_vector_view(double *, size_t count, int readonly);
ChemStatus chem_view_fill(ChemView *, double value);
ChemStatus chem_view_copy(ChemView *destination, const ChemView *source);
ChemStatus chem_view_for_each(ChemView *views, size_t count, ChemElementFunction, void *context);
/* Reduction operators: 0=sum, 1=max, 2=logical OR, 3=logical AND. */
ChemStatus chem_view_reduce(const ChemView *, int operation, int strict, double *accumulator);
ChemStatus chem_view_reduce_function(const ChemView *, int strict, double *accumulator, ChemReduceFunction, void *context);
ChemStatus chem_function_create(const ChemArgument *, size_t count, ChemArrayFunction, void *, ChemFunction **);
ChemStatus chem_function_invoke(const ChemFunction *, ChemArgument *, size_t count);
void chem_function_free(ChemFunction *);
#endif
