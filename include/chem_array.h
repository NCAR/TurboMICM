#ifndef CHEM_ARRAY_H
#define CHEM_ARRAY_H
#include "chem.h"

/* Runtime element types replace C++ matrix templates. Group length zero is
   row-major; positive lengths interleave that many rows, padding the last group. */
typedef enum { CHEM_REAL, CHEM_INTEGER, CHEM_BOOLEAN } ChemElementType;
typedef struct ChemArray ChemArray;
typedef struct {
    void *data;
    size_t length, capacity, block_length, block_stride, element_stride;
    ChemElementType type;
    int read_only;
} ChemView;
typedef void (*ChemElementFunction)(size_t count, void *const *elements, void *context);
typedef void (*ChemReductionFunction)(size_t count, const void *const *elements, void *accumulator, void *context);
typedef struct { ChemArray *array; int read_only; } ChemArrayArgument;
typedef ChemStatus (*ChemArrayFunctionBody)(size_t count, const ChemArrayArgument *, void *context);
typedef struct ChemArrayFunction ChemArrayFunction;

ChemStatus chem_array_create(size_t rows, size_t columns, size_t group_length,
    ChemElementType, const void *initial_value, ChemArray **);
void chem_array_free(ChemArray *);
ChemStatus chem_array_clone(const ChemArray *, ChemArray **);
/* Move transfers ownership and clears the source pointer; destination must be NULL. */
ChemStatus chem_array_move(ChemArray **destination, ChemArray **source);
size_t chem_array_rows(const ChemArray *);
size_t chem_array_columns(const ChemArray *);
size_t chem_array_group_length(const ChemArray *);
size_t chem_array_groups(const ChemArray *);
size_t chem_array_storage_size(const ChemArray *);
size_t chem_array_index(const ChemArray *, size_t row, size_t column);
void *chem_array_data(ChemArray *);
const void *chem_array_const_data(const ChemArray *);
ChemStatus chem_array_get(const ChemArray *, size_t row, size_t column, void *value);
ChemStatus chem_array_set(ChemArray *, size_t row, size_t column, const void *value);
ChemStatus chem_array_fill(ChemArray *, const void *value);
ChemStatus chem_array_copy(ChemArray *, const ChemArray *);
ChemStatus chem_array_assign_row(ChemArray *, size_t row, const void *values, size_t count);
ChemStatus chem_array_read_row(const ChemArray *, size_t row, void *values, size_t count);
ChemStatus chem_array_from_rows(size_t rows, const size_t *row_lengths,
    const void *const *row_data, size_t group_length, ChemElementType, ChemArray **);
ChemStatus chem_array_axpy(ChemArray *, double alpha, const ChemArray *);
ChemStatus chem_array_limit(ChemArray *, double value, int lower_bound);
/* Buffer formatting follows snprintf: result is the required length, excluding NUL. */
size_t chem_array_format(const ChemArray *, char *buffer, size_t capacity);
ChemStatus chem_array_row(const ChemArray *, size_t row, int read_only, ChemView *);
ChemStatus chem_array_column(const ChemArray *, size_t column, int read_only, ChemView *);
ChemStatus chem_array_group_column(const ChemArray *, size_t group, size_t column, int read_only, ChemView *);
ChemView chem_vector_view(void *data, size_t length, ChemElementType, int read_only);
void *chem_view_element(const ChemView *, size_t index);
ChemStatus chem_view_fill(ChemView, const void *value);
ChemStatus chem_view_copy(ChemView destination, ChemView source);
/* Apply visits logical rows only, retaining aliasing between arguments. The
   writable mask declares which arguments the callback may modify. */
ChemStatus chem_views_apply(size_t count, const ChemView *, size_t writable_mask, ChemElementFunction, void *context);
/* strict excludes padding; non-strict includes capacity common to all views. */
ChemStatus chem_views_reduce(size_t count, const ChemView *, int strict, ChemReductionFunction, void *accumulator, void *context);
ChemStatus chem_array_function_create(size_t count, const ChemArrayArgument *, ChemArrayFunctionBody, void *context, ChemArrayFunction **);
ChemStatus chem_array_function_call(const ChemArrayFunction *, size_t count, const ChemArrayArgument *);
void chem_array_function_free(ChemArrayFunction *);
#endif
