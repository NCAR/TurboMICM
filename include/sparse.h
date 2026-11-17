#ifndef CHEM_SPARSE_H
#define CHEM_SPARSE_H
#include <stddef.h>
#include <stdio.h>
#include "matrix.h"

typedef struct { size_t order, count; unsigned char *pattern; } ChemSparseBuilder;
typedef struct {
    size_t order, blocks, width, groups, nonzeros, storage_size;
    int column_major;
    size_t *indices;
    double *data;
} ChemSparseMatrix;
int chem_sparse_builder_init(ChemSparseBuilder *, size_t order);
void chem_sparse_builder_free(ChemSparseBuilder *);
int chem_sparse_builder_add(ChemSparseBuilder *, size_t row, size_t column);
int chem_sparse_init(ChemSparseMatrix *, const ChemSparseBuilder *, size_t blocks, size_t width, int column_major, double initial);
void chem_sparse_free(ChemSparseMatrix *);
int chem_sparse_index(const ChemSparseMatrix *, size_t block, size_t row, size_t column, size_t *index);
int chem_sparse_index_single(const ChemSparseMatrix *, size_t row, size_t column, size_t *index);
int chem_sparse_is_zero(const ChemSparseMatrix *, size_t row, size_t column, int *zero);
int chem_sparse_diagonal_indices(const ChemSparseMatrix *, size_t block, size_t *indices, size_t capacity, size_t *count);
double *chem_sparse_at(ChemSparseMatrix *, size_t block, size_t row, size_t column);
const double *chem_sparse_const_at(const ChemSparseMatrix *, size_t block, size_t row, size_t column);
void chem_sparse_fill(ChemSparseMatrix *, double);
void chem_sparse_add_diagonal(ChemSparseMatrix *, double);
int chem_sparse_print(const ChemSparseMatrix *, FILE *, int nonzero_only);
double *chem_sparse_at_vector_index(ChemSparseMatrix *, size_t block, size_t block_zero_index);
int chem_sparse_fill_block_view(ChemSparseMatrix *, size_t row, size_t column, double value);
int chem_sparse_block_to_vector(const ChemSparseMatrix *, size_t row, size_t column, double *values, size_t size);
int chem_sparse_vector_to_block(ChemSparseMatrix *, size_t row, size_t column, const double *values, size_t size);
int chem_sparse_copy_block_view(ChemSparseMatrix *out, size_t out_row, size_t out_column,
                                const ChemSparseMatrix *in, size_t in_row, size_t in_column);

typedef int (*ChemBlockFunction)(size_t block, ChemSparseMatrix *const *sparse,
                                 ChemMatrix *const *dense, ChemVector *const *vectors, void *context);
typedef struct {
    size_t width, sparse_count, dense_count, vector_count;
    ChemBlockFunction function;
    void *context;
} ChemSparseFunction;
int chem_sparse_function_init(ChemSparseFunction *, ChemSparseMatrix *const *sparse, size_t sparse_count,
                              ChemMatrix *const *dense, size_t dense_count, size_t vector_count,
                              ChemBlockFunction, void *context);
int chem_sparse_function_call(const ChemSparseFunction *, ChemSparseMatrix *const *sparse, size_t sparse_count,
                              ChemMatrix *const *dense, size_t dense_count, ChemVector *const *vectors, size_t vector_count);
#endif
