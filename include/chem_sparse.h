#ifndef CHEM_SPARSE_H
#define CHEM_SPARSE_H
#include "chem_array.h"

typedef enum { CHEM_CSR, CHEM_CSC } ChemCompression;
typedef struct ChemSparse ChemSparse;
/* Duplicate coordinates are coalesced. All blocks have the same pattern. */
ChemStatus chem_sparse_create(size_t order, size_t blocks, size_t group_length,
    ChemCompression, ChemElementType, size_t count, const ChemMatrixIndex *coordinates,
    const void *initial_value, ChemSparse **);
void chem_sparse_free(ChemSparse *);
ChemStatus chem_sparse_clone(const ChemSparse *, ChemSparse **);
size_t chem_sparse_order(const ChemSparse *);
size_t chem_sparse_nonzeros(const ChemSparse *);
ChemCompression chem_sparse_compression(const ChemSparse *);
/* Values are a blocks x nonzeros array. Its columns are sparse element views.
   This borrowed array can participate in the same functions as dense arrays. */
ChemArray *chem_sparse_values(ChemSparse *);
const size_t *chem_sparse_outer_starts(const ChemSparse *);
const size_t *chem_sparse_inner_indices(const ChemSparse *);
size_t chem_sparse_flat_index(const ChemSparse *, size_t row, size_t column);
size_t chem_sparse_index(const ChemSparse *, size_t block, size_t row, size_t column);
size_t chem_sparse_single_index(const ChemSparse *, size_t row, size_t column);
ChemStatus chem_sparse_get(const ChemSparse *, size_t block, size_t row, size_t column, void *value);
ChemStatus chem_sparse_set(ChemSparse *, size_t block, size_t row, size_t column, const void *value);
size_t chem_sparse_diagonal_indices(const ChemSparse *, size_t block, size_t *indices, size_t capacity);
ChemStatus chem_sparse_add_diagonal(ChemSparse *, const void *value);
ChemStatus chem_sparse_block_view(const ChemSparse *, size_t row, size_t column, int read_only, ChemView *);
/* vector_index is a physical first-block index, as returned by sparse_index. */
ChemStatus chem_sparse_index_view(const ChemSparse *, size_t vector_index, int read_only, ChemView *);
size_t chem_sparse_format(const ChemSparse *, char *buffer, size_t capacity, int nonzeros_only);
#endif
