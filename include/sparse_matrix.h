#ifndef CHEM_SPARSE_MATRIX_H
#define CHEM_SPARSE_MATRIX_H
#include "matrix.h"

typedef struct { size_t row, column; } ChemSparseEntry;
typedef enum { CHEM_SPARSE_CSR, CHEM_SPARSE_CSC } ChemSparseOrdering;
/* All blocks have the same sparsity; values use ChemMatrix group layout. */
typedef struct {
    size_t dimension, nonzeros;
    ChemSparseOrdering ordering;
    size_t *offsets, *indices;
    ChemMatrix values;
} ChemSparseMatrix;
/* Duplicate entries are coalesced. Initialize with {0}; invalid coordinates
   return an error and leave the existing matrix unchanged. */
ChemMatrixStatus chem_sparse_matrix_init(ChemSparseMatrix *m,size_t dimension,size_t blocks,size_t width,
        ChemSparseOrdering ordering,const ChemSparseEntry *entries,size_t count,double initial_value);
void chem_sparse_matrix_destroy(ChemSparseMatrix *m);
ChemMatrixStatus chem_sparse_matrix_clone(ChemSparseMatrix *dst,const ChemSparseMatrix *src);
/* Returns SIZE_MAX when the entry is absent or the coordinates are invalid. */
size_t chem_sparse_matrix_entry(const ChemSparseMatrix *m,size_t row,size_t column);
size_t chem_sparse_matrix_index(const ChemSparseMatrix *m,size_t block,size_t row,size_t column);
double *chem_sparse_matrix_at(ChemSparseMatrix *m,size_t block,size_t row,size_t column);
const double *chem_sparse_matrix_const_at(const ChemSparseMatrix *m,size_t block,size_t row,size_t column);
ChemMatrixStatus chem_sparse_matrix_diagonal_indices(const ChemSparseMatrix *m,size_t block,
        size_t *indices,size_t capacity,size_t *count);
/* Borrow a view across all blocks of one structural entry. The index form
   accepts a physical storage index for block zero. */
ChemMatrixStatus chem_sparse_matrix_block(ChemSparseMatrix *m,size_t row,size_t column,ChemColumnView *view);
ChemMatrixStatus chem_sparse_matrix_block_by_index(ChemSparseMatrix *m,size_t block_zero_index,ChemColumnView *view);
ChemMatrixStatus chem_sparse_matrix_group_block(ChemSparseMatrix *m,size_t group,size_t row,size_t column,ChemGroupSpan *span);
ChemMatrixStatus chem_sparse_matrix_const_group_block(const ChemSparseMatrix *m,size_t group,size_t row,size_t column,ChemConstGroupSpan *span);
void chem_sparse_matrix_fill(ChemSparseMatrix *m,double value);
void chem_sparse_matrix_add_diagonal(ChemSparseMatrix *m,double value);
ChemMatrixStatus chem_sparse_matrix_print(const ChemSparseMatrix *m,FILE *stream,int nonzeros_only);
#endif
