#ifndef CHEM_SPARSE_MATRIX_H
#define CHEM_SPARSE_MATRIX_H

#include "chem_views.h"

#include <stddef.h>
#include <stdio.h>

typedef enum {
  CHEM_SPARSE_ROW_COMPRESSED = 0,
  CHEM_SPARSE_COLUMN_COMPRESSED = 1
} chem_sparse_order;

typedef enum {
  CHEM_SPARSE_STANDARD_BLOCKS = 0,
  CHEM_SPARSE_GROUPED_BLOCKS = 1
} chem_sparse_layout;

typedef struct {
  size_t dimension;
  size_t nonzero_count;
  size_t block_count;
  size_t lanes;
  size_t storage_length;
  chem_sparse_order order;
  chem_sparse_layout layout;
  size_t *rows;
  size_t *columns;
  double *values;
} chem_sparse_matrix;

/* Initialize an empty instance with {0}. Coordinates define the shared sparse
   pattern for all blocks. Standard storage requires lanes=1; grouped storage
   stores each lane group contiguously and pads its final partial group.
   Return -1 for invalid patterns/dimensions and -2 on allocation failure. */
int chem_sparse_matrix_init(
    chem_sparse_matrix *matrix, size_t dimension,
    const size_t *rows, const size_t *columns, size_t nonzero_count,
    size_t block_count, size_t lanes, chem_sparse_order order,
    chem_sparse_layout layout, double initial_value);
void chem_sparse_matrix_destroy(chem_sparse_matrix *matrix);
int chem_sparse_matrix_vector_index(
    const chem_sparse_matrix *matrix, size_t row, size_t column, size_t *index);
int chem_sparse_matrix_is_zero(
    const chem_sparse_matrix *matrix, size_t row, size_t column, int *is_zero);
int chem_sparse_matrix_get(
    const chem_sparse_matrix *matrix, size_t block, size_t row, size_t column, double *value);
int chem_sparse_matrix_set(
    chem_sparse_matrix *matrix, size_t block, size_t row, size_t column, double value);
int chem_sparse_matrix_block_group_view(
    chem_sparse_matrix *matrix, size_t group, size_t entry, chem_group_view *view);
int chem_sparse_matrix_const_block_group_view(
    const chem_sparse_matrix *matrix, size_t group, size_t entry, chem_const_group_view *view);
int chem_sparse_matrix_fill(chem_sparse_matrix *matrix, double value);
int chem_sparse_matrix_copy(chem_sparse_matrix *destination, const chem_sparse_matrix *source);
int chem_sparse_matrix_add_to_diagonal(chem_sparse_matrix *matrix, double value);
int chem_sparse_matrix_print(const chem_sparse_matrix *matrix, FILE *stream);
int chem_sparse_matrix_print_nonzero(const chem_sparse_matrix *matrix, FILE *stream);

#endif
