#ifndef CHEM_SPARSE_H
#define CHEM_SPARSE_H

#include "chem_matrix.h"
#include <stddef.h>

typedef struct { size_t row, column; } ChemSparseCoordinate;
typedef enum { CHEM_SPARSE_ROW, CHEM_SPARSE_COLUMN } ChemSparseOrder;
typedef enum { CHEM_NUMBER_INT, CHEM_NUMBER_DOUBLE } ChemNumberType;
typedef struct {
    size_t order, nonzeros;
    ChemSparseOrder ordering;
    size_t *lookup;
    ChemSparseCoordinate *coordinates;
    /* Each row is a block, each column a structural nonzero. */
    ChemMatrix storage;
} ChemSparseMatrix;

int chem_sparse_init(ChemSparseMatrix *matrix, size_t order, size_t blocks,
                     const ChemSparseCoordinate *coordinates, size_t count,
                     ChemSparseOrder ordering, size_t group_size, size_t element_size,
                     const void *initial_value);
void chem_sparse_destroy(ChemSparseMatrix *matrix);
size_t chem_sparse_entry(const ChemSparseMatrix *matrix, size_t row, size_t column);
size_t chem_sparse_slot(const ChemSparseMatrix *matrix, size_t block, size_t row, size_t column);
void *chem_sparse_at(ChemSparseMatrix *matrix, size_t block, size_t row, size_t column);
const void *chem_sparse_at_const(const ChemSparseMatrix *matrix, size_t block, size_t row, size_t column);
/* A block-zero slot addresses a given nonzero across blocks. */
void *chem_sparse_block_slot(ChemSparseMatrix *matrix, size_t block, size_t block_zero_slot);
int chem_sparse_group_block(ChemSparseMatrix *matrix, size_t group, size_t row, size_t column, ChemGroupView *view);
int chem_sparse_const_group_block(const ChemSparseMatrix *matrix, size_t group, size_t row, size_t column,
                                  ChemConstGroupView *view);
int chem_sparse_fill(ChemSparseMatrix *matrix, const void *value);
int chem_sparse_add_diagonal(ChemSparseMatrix *matrix, const void *value, ChemNumberType type);
int chem_sparse_format(const ChemSparseMatrix *matrix, char *buffer, size_t capacity,
                       ChemNumberType type, int nonzeros_only);

#endif
