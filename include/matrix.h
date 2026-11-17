#ifndef CHEM_MATRIX_H
#define CHEM_MATRIX_H
#include <stddef.h>
#include <stdio.h>

typedef struct {
    size_t rows, columns, width, groups, storage_size;
    int grouped;
    double *data;
} ChemMatrix;
typedef void (*ChemElementFunction)(double *out, const double *a, const double *b, void *context);

int chem_matrix_init(ChemMatrix *, size_t rows, size_t columns, size_t width, int grouped, double initial);
void chem_matrix_free(ChemMatrix *);
int chem_matrix_from_rows(ChemMatrix *, size_t rows, const double *const *values, const size_t *lengths, size_t width, int grouped);
double *chem_matrix_at(ChemMatrix *, size_t row, size_t column);
const double *chem_matrix_const_at(const ChemMatrix *, size_t row, size_t column);
size_t chem_matrix_row_stride(const ChemMatrix *);
size_t chem_matrix_column_stride(const ChemMatrix *);
int chem_matrix_get_row(const ChemMatrix *, size_t row, double *values, size_t count);
int chem_matrix_set_row(ChemMatrix *, size_t row, const double *values, size_t count);
int chem_matrix_copy(ChemMatrix *, const ChemMatrix *);
void chem_matrix_fill(ChemMatrix *, double);
void chem_matrix_max(ChemMatrix *, double);
void chem_matrix_min(ChemMatrix *, double);
int chem_matrix_axpy(ChemMatrix *, double, const ChemMatrix *);
int chem_matrix_foreach(ChemMatrix *, const ChemMatrix *, const ChemMatrix *, ChemElementFunction, void *);
int chem_matrix_print(const ChemMatrix *, FILE *);
typedef void (*ChemMatrixReducer)(double value, void *context);
int chem_matrix_reduce_column(const ChemMatrix *, size_t column, int strict, ChemMatrixReducer, void *context);
int chem_matrix_fill_column(ChemMatrix *, size_t column, double value);
int chem_matrix_copy_column(ChemMatrix *destination, size_t to_column, const ChemMatrix *source, size_t from_column);
int chem_matrix_column_to_vector(const ChemMatrix *, size_t column, double *out, size_t size);
int chem_matrix_vector_to_column(ChemMatrix *, size_t column, const double *values, size_t size);

typedef struct { void *data; size_t size; } ChemVector;
typedef int (*ChemRowFunction)(size_t row, ChemMatrix *const *matrices, ChemVector *const *vectors, void *context);
typedef struct {
    size_t matrix_count, vector_count;
    size_t *columns;
    ChemRowFunction function;
    void *context;
} ChemMatrixFunction;
int chem_matrix_function_init(ChemMatrixFunction *, ChemMatrix *const *prototypes, size_t matrix_count,
                              size_t vector_count, ChemRowFunction, void *context);
void chem_matrix_function_free(ChemMatrixFunction *);
int chem_matrix_function_call(const ChemMatrixFunction *, ChemMatrix *const *matrices, size_t matrix_count,
                              ChemVector *const *vectors, size_t vector_count);

typedef struct { size_t rows, columns, width, groups; int *data; } ChemIntMatrix;
typedef void (*ChemIntRowFunction)(int *output, int input, void *context);
int chem_int_matrix_init(ChemIntMatrix *, size_t rows, size_t columns, size_t width, int initial);
void chem_int_matrix_free(ChemIntMatrix *);
int *chem_int_matrix_at(ChemIntMatrix *, size_t row, size_t column);
int chem_int_matrix_apply_vector(ChemIntMatrix *, size_t column, const int *, size_t size, ChemIntRowFunction, void *);
#endif
