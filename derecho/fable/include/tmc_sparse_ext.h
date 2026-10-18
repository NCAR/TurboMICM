/*
 * tmc: the rest of MICM's sparse container interface.
 *
 * MICM's SparseMatrix<T, SparseMatrixStandardOrderingCompressedSparseRow> is
 * tmc_sparse with TMC_SPARSE_CSR and L = 1, its ...CompressedSparseColumn
 * variant TMC_SPARSE_CSC, and the SparseMatrixVectorOrdering...<L> variants
 * the same orderings with width L (tmc.h). The solver needs only the
 * operations in tmc.h; this header adds the parts of MICM's container
 * interface that it does not use: VectorIndex with MICM's checks (including
 * the two-argument form of a one-block matrix), a checked IsZero,
 * DiagonalIndices, the group geometry of the vector orderings, and printing.
 *
 * MICM's Function facility over sparse matrices (GroupView, BlockView,
 * BlockVariable, ForEachBlock, Fill, Copy) is tmc_array_function_* in
 * tmc_dense_ext.h, whose tmc_array_group_element takes an element id as
 * tmc_sparse_element returns it: MICM's GetBlockView(VectorIndex(0, row,
 * column)) is tmc_array_group_element(g, arg, tmc_sparse_element(m, row,
 * column)).
 *
 * Data order. Element ids number the non-zero elements of a block row by row
 * (CSR) or column by column (CSC), each row or column ascending; element e of
 * block b is stored at ((b / L) * nnz + e) * L + b % L (tmc_sparse_data_index).
 * At L = 1 that is MICM's standard orderings, block * nnz + e; at L > 1 it is
 * MICM's vector orderings, e * L + b % L + (b / L) * L * nnz, in which the
 * block count is padded up to a multiple of L.
 */
#ifndef TMC_SPARSE_EXT_H
#define TMC_SPARSE_EXT_H

#include "tmc.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* MICM's VectorIndex(block, row, column): the data index of a non-zero
 * element, with MICM's checks in MICM's order: TMC_MATRIX_ELEMENT_OUT_OF_RANGE
 * when row, column, or block is out of range, then
 * TMC_MATRIX_ZERO_ELEMENT_ACCESS when (row, column) is a structural zero.
 * Returns (size_t)-1 and sets err (category TMC_ERROR_MATRIX) on failure. */
size_t tmc_sparse_vector_index(const tmc_sparse *m, int block, int row, int column, tmc_error *err);

/* MICM's VectorIndex(row, column): the index in block 0 of a one-block matrix.
 * A matrix with any other block count fails with
 * TMC_MATRIX_MISSING_BLOCK_INDEX before the element is looked at. */
size_t tmc_sparse_vector_index_one_block(const tmc_sparse *m, int row, int column, tmc_error *err);

/* MICM's IsZero with its range check: true when (row, column) is a structural
 * zero. An element out of range is reported as zero, with
 * TMC_MATRIX_ELEMENT_OUT_OF_RANGE set in err. */
bool tmc_sparse_is_zero_checked(const tmc_sparse *m, int row, int column, tmc_error *err);

/* MICM's DiagonalIndices(block): the data indices of the non-zero diagonal
 * elements of `block`, by row, written to indices (which needs room for m->n
 * values) unless it is NULL. Returns the count, or -1 with
 * TMC_MATRIX_ELEMENT_OUT_OF_RANGE when the block does not exist. */
int tmc_sparse_diagonal_indices(const tmc_sparse *m, int block, size_t *indices, tmc_error *err);

/* Group geometry, as MICM's vector orderings report it: GroupVectorSize (the
 * width L), GroupSize (the values per group of blocks, L * nnz), and
 * NumberOfGroups(n_blocks) (the groups of L blocks that n_blocks blocks
 * occupy, a partial last group included). */
int tmc_sparse_group_width(const tmc_sparse *m);
int tmc_sparse_group_size(const tmc_sparse *m);
int tmc_sparse_number_of_groups(const tmc_sparse *m, int n_blocks);

/* MICM's stream output: for each block, a line "Block <b>", then each row of
 * the block with its values formatted as %g (the C++ stream default) and its
 * structural zeros as 0, separated by commas, each row ending in a newline.
 * tmc_sparse_format writes at most size - 1 characters and a terminating NUL
 * (nothing when size is 0) and returns the length of the whole text, as
 * snprintf does. tmc_sparse_print writes to a stream and returns the number
 * of characters written, or a negative value on an output error. */
size_t tmc_sparse_format(const tmc_sparse *m, char *buffer, size_t size);
int tmc_sparse_print(const tmc_sparse *m, FILE *stream);

/* MICM's PrintNonZeroElements: for each block, a line "Block <b>", then one
 * line "<row>, <column>, <value>" per non-zero element, rows then columns
 * ascending whatever the element numbering. The same conventions as
 * tmc_sparse_format and tmc_sparse_print. */
size_t tmc_sparse_format_non_zero(const tmc_sparse *m, char *buffer, size_t size);
int tmc_sparse_print_non_zero(const tmc_sparse *m, FILE *stream);

#ifdef __cplusplus
}
#endif

#endif
