/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_SPARSE_MATRIX_H
#define CHEM_SPARSE_MATRIX_H
#include "chem_matrix.h"

typedef enum { CHEM_SPARSE_COLUMN, CHEM_SPARSE_ROW } ChemSparseOrdering;
typedef struct { size_t row,column; } ChemSparseElement;
typedef struct ChemSparseBuilder ChemSparseBuilder;
typedef struct {
    size_t order;
    ChemSparseOrdering ordering;
    ChemSparseElement *elements; /* sorted, unique structural entries */
    ChemMatrix values;          /* logical rows are blocks, columns are entries */
} ChemSparseMatrix;
ChemStatus chem_sparse_builder_create(size_t order,ChemSparseBuilder **);
void chem_sparse_builder_free(ChemSparseBuilder *);
ChemStatus chem_sparse_builder_add(ChemSparseBuilder *,size_t row,size_t column);
ChemStatus chem_sparse_builder_blocks(ChemSparseBuilder *,size_t blocks);
void chem_sparse_builder_initial(ChemSparseBuilder *,double value);
size_t chem_sparse_builder_elements(const ChemSparseBuilder *);
/* group=0 is standard block-major; group>0 groups contiguous block lanes.
 * No diagonal or fill is added. Initialize outputs with {0} and free before reuse. */
ChemStatus chem_sparse_matrix_build(const ChemSparseBuilder *,size_t group,ChemSparseOrdering,ChemMatrixType,ChemSparseMatrix *);
void chem_sparse_matrix_free(ChemSparseMatrix *);
ChemStatus chem_sparse_matrix_copy(const ChemSparseMatrix *,ChemSparseMatrix *);
ChemStatus chem_sparse_matrix_index(const ChemSparseMatrix *,size_t block,size_t row,size_t column,size_t *out);
ChemStatus chem_sparse_matrix_single_index(const ChemSparseMatrix *,size_t row,size_t column,size_t *out);
ChemStatus chem_sparse_matrix_is_zero(const ChemSparseMatrix *,size_t row,size_t column,int *out);
double chem_sparse_matrix_get(const ChemSparseMatrix *,size_t block,size_t row,size_t column);
ChemStatus chem_sparse_matrix_set(ChemSparseMatrix *,size_t block,size_t row,size_t column,double);
ChemStatus chem_sparse_matrix_diagonal_indices(const ChemSparseMatrix *,size_t block,size_t *out,size_t capacity,size_t *count);
ChemStatus chem_sparse_matrix_add_diagonal(ChemSparseMatrix *,double value);
ChemStatus chem_sparse_matrix_print(const ChemSparseMatrix *,FILE *,int nonzero_only);
ChemArray chem_array_sparse_matrix(ChemSparseMatrix *);
ChemArray chem_array_const_sparse_matrix(const ChemSparseMatrix *);
/* vector_index is the raw offset in block zero, as returned by index(). */
ChemStatus chem_sparse_block_view(ChemArray,size_t vector_index,ChemView *);
#endif
