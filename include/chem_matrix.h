#ifndef CHEM_MATRIX_H
#define CHEM_MATRIX_H

#include <stddef.h>

/* Matrix elements may be any trivially copyable C type. With group_size zero,
 * storage is row-major; otherwise groups store each column in adjacent lanes. */
typedef struct {
    size_t rows, columns, group_size, element_size, slots;
    unsigned char *data;
} ChemMatrix;

int chem_matrix_init(ChemMatrix *matrix, size_t rows, size_t columns, size_t group_size,
                     size_t element_size, const void *initial_value);
int chem_matrix_from_rows(ChemMatrix *matrix, const void *const *rows, const size_t *lengths,
                          size_t row_count, size_t group_size, size_t element_size);
void chem_matrix_destroy(ChemMatrix *matrix);
void *chem_matrix_at(ChemMatrix *matrix, size_t row, size_t column);
const void *chem_matrix_at_const(const ChemMatrix *matrix, size_t row, size_t column);
size_t chem_matrix_row_stride(const ChemMatrix *matrix);
size_t chem_matrix_column_stride(const ChemMatrix *matrix);
int chem_matrix_get_row(const ChemMatrix *matrix, size_t row, void *values, size_t count);
int chem_matrix_set_row(ChemMatrix *matrix, size_t row, const void *values, size_t count);
int chem_matrix_get_column(const ChemMatrix *matrix, size_t column, void *values, size_t count);
int chem_matrix_set_column(ChemMatrix *matrix, size_t column, const void *values, size_t count);
int chem_matrix_fill(ChemMatrix *matrix, const void *value);
int chem_matrix_fill_column(ChemMatrix *matrix, size_t column, const void *value);
int chem_matrix_copy(ChemMatrix *destination, const ChemMatrix *source);
int chem_matrix_copy_column(ChemMatrix *destination, size_t dst_column,
                            const ChemMatrix *source, size_t src_column);
int chem_matrix_axpy(ChemMatrix *destination, double alpha, const ChemMatrix *source);
int chem_matrix_bound(ChemMatrix *matrix, double value, int minimum);
/* Format double-valued matrices as comma-separated rows, including newlines. */
int chem_matrix_format(const ChemMatrix *matrix, char *buffer, size_t capacity);

typedef void (*ChemMatrixReducer)(const void *value, void *accumulator, void *context);
/* A non-strict reduction includes padded lanes, while strict visits real rows. */
int chem_matrix_reduce_column(const ChemMatrix *matrix, size_t column, int strict,
                              ChemMatrixReducer reducer, void *accumulator, void *context);

typedef int (*ChemMatrixFunctionBody)(ChemMatrix *const *matrices, void *const *vectors, void *context);
typedef struct {
    size_t matrix_count, vector_count, *columns;
    ChemMatrixFunctionBody body;
    void *context;
} ChemMatrixFunction;
/* Creation records column counts; row and vector lengths are checked at call
 * time, so a function may be reused with a different number of rows. */
int chem_matrix_function_init(ChemMatrixFunction *function, ChemMatrix *const *matrices,
                              size_t matrix_count, size_t vector_count,
                              ChemMatrixFunctionBody body, void *context);
void chem_matrix_function_destroy(ChemMatrixFunction *function);
int chem_matrix_function_call(const ChemMatrixFunction *function, ChemMatrix *const *matrices,
                              void *const *vectors, const size_t *vector_lengths);

#endif
