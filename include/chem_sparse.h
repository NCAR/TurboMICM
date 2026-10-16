#ifndef CHEM_SPARSE_H
#define CHEM_SPARSE_H
#include "chem_array.h"
typedef enum { CHEM_CSR, CHEM_CSC } chem_sparse_order;
typedef struct {
    size_t n, blocks, nnz;
    unsigned char *pattern;
    double initial;
} chem_sparse_builder;
typedef struct {
    size_t n, nnz;
    size_t *map, *starts, *indices;
    chem_sparse_order order;
    chem_array values;
} chem_sparse;
chem_status chem_sparse_builder_create(chem_sparse_builder *builder, size_t n);
chem_status chem_sparse_builder_add(chem_sparse_builder *builder, size_t row, size_t column);
void chem_sparse_builder_destroy(chem_sparse_builder *builder);
chem_status chem_sparse_create(chem_sparse *matrix, const chem_sparse_builder *builder,
                               chem_sparse_order ordering, size_t vector_length, chem_element_type type);
void chem_sparse_destroy(chem_sparse *matrix);
chem_status chem_sparse_clone(chem_sparse *destination, const chem_sparse *source);
chem_status chem_sparse_is_zero(const chem_sparse *matrix, size_t row, size_t column, int *zero);
/* CHEM_NONE means omitted block index, valid only for single-block matrices. */
chem_status chem_sparse_index(const chem_sparse *matrix, size_t block, size_t row, size_t column, size_t *index);
chem_status chem_sparse_get(const chem_sparse *matrix, size_t block, size_t row, size_t column, double *value);
chem_status chem_sparse_set(chem_sparse *matrix, size_t block, size_t row, size_t column, double value);
chem_status chem_sparse_diagonal_indices(const chem_sparse *matrix, size_t block,
                                         size_t *indices, size_t capacity, size_t *count);
void chem_sparse_add_diagonal(chem_sparse *matrix, double value);
chem_status chem_sparse_print(FILE *file, const chem_sparse *matrix, int nonzero_only);
/* index is a physical storage index in block zero; cache it before hot loops. */
chem_span chem_sparse_block_span(chem_array_group group, size_t index);
chem_span chem_sparse_const_block_span(chem_array_group group, size_t index);
chem_span chem_sparse_coordinate_span(const chem_sparse *matrix, chem_array_group group, size_t row, size_t column);
chem_raw_view chem_sparse_block(const chem_sparse *matrix, size_t index, int writable);
#endif
