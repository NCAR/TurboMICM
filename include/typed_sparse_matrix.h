#ifndef CHEM_TYPED_SPARSE_MATRIX_H
#define CHEM_TYPED_SPARSE_MATRIX_H
#include "sparse_matrix.h"
#include "typed_matrix.h"

/* Sparse blocks with values of any trivially copyable C type. Pointers stored
   as values remain caller-owned. The shared sparsity and values are owned. */
typedef struct { ChemSparseMatrix structure; ChemTypedMatrix values; } ChemTypedSparseMatrix;
ChemMatrixStatus chem_typed_sparse_init(ChemTypedSparseMatrix *matrix,size_t dimension,size_t blocks,size_t width,
    ChemSparseOrdering ordering,const ChemSparseEntry *entries,size_t count,size_t element_size,const void *initial_value);
ChemMatrixStatus chem_typed_sparse_clone(ChemTypedSparseMatrix *destination,const ChemTypedSparseMatrix *source);
void chem_typed_sparse_destroy(ChemTypedSparseMatrix *matrix);
void *chem_typed_sparse_at(ChemTypedSparseMatrix *matrix,size_t block,size_t row,size_t column);
const void *chem_typed_sparse_const_at(const ChemTypedSparseMatrix *matrix,size_t block,size_t row,size_t column);
#endif
