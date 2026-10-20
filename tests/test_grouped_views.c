#include "chem_views.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double scale_value(double value, void *context) {
  return value * *(const double *)context;
}

static double add_offset(double value, void *context) {
  return value + *(const double *)context;
}

static double scale_and_offset(double value, void *context) {
  const double *parameters = context;
  return value * parameters[0] + parameters[1];
}

static double add_values(double left, double right, void *context) {
  (void)context;
  return left + right;
}

static double subtract_values(double left, double right, void *context) {
  (void)context;
  return left - right;
}

static int close_value(double actual, double expected) {
  return fabs(actual - expected) <= 1e-12 * (1.0 + fabs(expected));
}

static size_t vector_matrix_index(size_t row, size_t column, size_t columns, size_t lanes) {
  return (row / lanes * columns + column) * lanes + row % lanes;
}

static int dense_group_roundtrip(void) {
  double matrix[12];
  for (size_t row = 0; row < 3u; ++row)
    for (size_t column = 0; column < 4u; ++column) matrix[row * 4u + column] = 10.0 * (double)row + (double)column;
  chem_group_view destination;
  chem_const_group_view source;
  const double factor = 2.0;
  if (chem_group_view_init(matrix, 12u, 6u, 1u, 4u, &destination) != 0 ||
      chem_const_group_view_init(matrix, 12u, 7u, 1u, 4u, &source) != 0 ||
      !close_value(*destination.base, 12.0) || !close_value(*source.base, 13.0) ||
      chem_group_view_apply(destination, source, scale_value, (void *)&factor) != 0) return 1;
  return close_value(matrix[1u * 4u + 2u], 26.0) &&
      close_value(matrix[0u * 4u + 2u], 2.0) && close_value(matrix[2u * 4u + 2u], 22.0) ? 0 : 1;
}

static int vector_matrix_column_block(void) {
  double matrix[24];
  for (size_t row = 0; row < 8u; ++row)
    for (size_t column = 0; column < 3u; ++column) matrix[row * 3u + column] = 100.0 + 10.0 * row + column;
  chem_group_view destination;
  chem_const_group_view source;
  const double factor = 2.0;
  if (chem_group_view_init(matrix, 24u, 13u, 4u, 3u, &destination) != 0 ||
      chem_const_group_view_init(matrix, 24u, 14u, 4u, 3u, &source) != 0) return 1;
  for (size_t lane = 0; lane < 4u; ++lane)
    if (!close_value(destination.base[lane * destination.stride], 141.0 + 10.0 * lane) ||
        !close_value(source.base[lane * source.stride], 142.0 + 10.0 * lane)) return 1;
  if (chem_group_view_apply(destination, source, scale_value, (void *)&factor) != 0) return 1;
  for (size_t row = 4u; row < 8u; ++row)
    if (!close_value(matrix[row * 3u + 1u], 2.0 * (100.0 + 10.0 * row + 2.0))) return 1;
  for (size_t row = 0u; row < 4u; ++row)
    if (!close_value(matrix[row * 3u + 1u], 100.0 + 10.0 * row + 1.0)) return 1;
  return 0;
}

static int dense_grouped_matches_raw(size_t lanes) {
  const size_t rows = 3u * lanes, columns = 4u, count = rows * columns;
  double *matrix = malloc(count * sizeof(*matrix));
  if (!matrix) return 1;
  for (size_t row = 0; row < rows; ++row)
    for (size_t column = 0; column < columns; ++column)
      matrix[vector_matrix_index(row, column, columns, lanes)] = 10.0 * row + column;
  const double parameters[] = { 3.0, 1.0 };
  for (size_t group = 0; group < 3u; ++group) {
    chem_group_view destination;
    chem_const_group_view source;
    const size_t group_base = group * columns * lanes;
    if (chem_group_view_init(matrix, count, group_base + 3u * lanes, lanes, 1u, &destination) != 0 ||
        chem_const_group_view_init(matrix, count, group_base, lanes, 1u, &source) != 0 ||
        chem_group_view_apply(destination, source, scale_and_offset, (void *)parameters) != 0) {
      free(matrix);
      return 1;
    }
  }
  int failed = 0;
  for (size_t row = 0; row < rows; ++row) {
    for (size_t column = 0; column < columns; ++column) {
      const double original = 10.0 * row + column;
      const double expected = column == 3u ? 30.0 * row + 1.0 : original;
      if (!close_value(matrix[vector_matrix_index(row, column, columns, lanes)], expected)) failed = 1;
    }
  }
  free(matrix);
  return failed;
}

static int plain_matrix_roundtrip(void) {
  double matrix[12];
  for (size_t row = 0; row < 4u; ++row) {
    matrix[row * 3u] = (double)row;
    matrix[row * 3u + 1u] = (double)row + 5.0;
    matrix[row * 3u + 2u] = -1.0;
  }
  chem_group_view destination;
  chem_const_group_view left, right;
  if (chem_group_view_init(matrix, 12u, 2u, 4u, 3u, &destination) != 0 ||
      chem_const_group_view_init(matrix, 12u, 0u, 4u, 3u, &left) != 0 ||
      chem_const_group_view_init(matrix, 12u, 1u, 4u, 3u, &right) != 0 ||
      chem_group_view_apply_binary(destination, left, right, subtract_values, NULL) != 0) return 1;
  for (size_t row = 0; row < 4u; ++row)
    if (!close_value(matrix[row * 3u + 2u], -5.0)) return 1;
  return 0;
}

static int sparse_standard_roundtrip(void) {
  double storage[] = { 1.0, 2.0, 3.0, 10.0, 20.0, 30.0 };
  chem_group_view destination;
  chem_const_group_view source;
  const double increment = 5.0;
  if (chem_group_view_init(storage, 6u, 3u, 1u, 3u, &destination) != 0 ||
      chem_const_group_view_init(storage, 6u, 4u, 1u, 3u, &source) != 0 ||
      chem_group_view_apply(destination, source, add_offset, (void *)&increment) != 0) return 1;
  return close_value(storage[3], 25.0) && close_value(storage[0], 1.0) && close_value(storage[5], 30.0) ? 0 : 1;
}

static int sparse_vector_contiguous(void) {
  double storage[16];
  for (size_t block = 0; block < 8u; ++block) {
    storage[(block / 4u * 2u + 0u) * 4u + block % 4u] = 100.0 + block;
    storage[(block / 4u * 2u + 1u) * 4u + block % 4u] = 200.0 + block;
  }
  chem_group_view destination;
  chem_const_group_view source;
  if (chem_group_view_init(storage, 16u, 8u, 4u, 1u, &destination) != 0 ||
      chem_const_group_view_init(storage, 16u, 12u, 4u, 1u, &source) != 0) return 1;
  for (size_t lane = 0; lane < 4u; ++lane)
    if (!close_value(destination.base[lane], 104.0 + lane) || !close_value(source.base[lane], 204.0 + lane)) return 1;
  if (chem_group_view_apply_binary(destination, (chem_const_group_view){destination.base, 4u, 1u},
          source, add_values, NULL) != 0) return 1;
  for (size_t block = 4u; block < 8u; ++block)
    if (!close_value(storage[(block / 4u * 2u) * 4u + block % 4u], 308.0 + 2.0 * (block - 4u))) return 1;
  for (size_t block = 0u; block < 4u; ++block)
    if (!close_value(storage[block], 100.0 + block)) return 1;
  return 0;
}

static int mixed_sparse_dense(size_t lanes) {
  const size_t blocks = 2u * lanes, columns = 3u, sparse_entries = 1u;
  double *dense = malloc(blocks * columns * sizeof(*dense));
  double *sparse = malloc(blocks * sparse_entries * sizeof(*sparse));
  if (!dense || !sparse) { free(dense); free(sparse); return 1; }
  for (size_t block = 0; block < blocks; ++block) {
    sparse[block] = 1.0 + block;
    for (size_t column = 0; column < columns; ++column) dense[block * columns + column] = 10.0 + block + column;
  }
  const size_t vector_size = blocks * columns;
  for (size_t group = 0; group < 2u; ++group) {
    chem_group_view output;
    chem_const_group_view dense_input, sparse_input;
    const size_t group_block = group * lanes;
    if (chem_group_view_init(dense, vector_size, group_block * columns, lanes, columns, &output) != 0 ||
        chem_const_group_view_init(dense, vector_size, group_block * columns + 1u, lanes, columns, &dense_input) != 0 ||
        chem_const_group_view_init(sparse, blocks, group_block, lanes, 1u, &sparse_input) != 0 ||
        chem_group_view_apply_binary(output, dense_input, sparse_input, add_values, NULL) != 0) {
      free(dense); free(sparse); return 1;
    }
  }
  int failed = 0;
  for (size_t block = 0; block < blocks; ++block) {
    for (size_t column = 0; column < columns; ++column) {
      const double expected = column == 0u ? 10.0 + block + 1.0 + 1.0 + block : 10.0 + block + column;
      if (!close_value(dense[block * columns + column], expected)) failed = 1;
    }
  }
  free(dense); free(sparse);
  return failed;
}

static int sparse_group_sum(void) {
  double storage[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  chem_const_group_view group;
  if (chem_const_group_view_init(storage, 8u, 4u, 4u, 1u, &group) != 0) return 1;
  double sum = 0.0;
  for (size_t i = 0; i < group.length; ++i) sum += group.base[i * group.stride];
  return close_value(sum, 26.0) ? 0 : 1;
}

static int const_views(void) {
  const double sparse[] = { 1, 2, 3, 4 };
  const double dense[] = { 100, 101, 102, 103, 110, 111, 112, 113 };
  chem_const_group_view sparse_group, dense_group;
  if (chem_const_group_view_init(sparse, 4u, 0u, 4u, 1u, &sparse_group) != 0 ||
      chem_const_group_view_init(dense, 8u, 4u, 4u, 1u, &dense_group) != 0) return 1;
  for (size_t i = 0; i < 4u; ++i)
    if (!close_value(sparse_group.base[i], (double)i + 1.0) ||
        !close_value(dense_group.base[i], 110.0 + (double)i)) return 1;
  return 0;
}

static int sparse_fill_copy_standard(void) {
  double storage[] = { 1, 100, 2, 200, 3, 300 };
  chem_group_view destination, copy_destination;
  chem_const_group_view source;
  const double fill = 42.0;
  if (chem_group_view_init(storage, 6u, 0u, 1u, 3u, &destination) != 0 ||
      chem_group_view_fill(destination, fill) != 0 ||
      chem_group_view_init(storage, 6u, 1u, 1u, 3u, &copy_destination) != 0 ||
      chem_const_group_view_init(storage, 6u, 0u, 1u, 3u, &source) != 0 ||
      chem_group_view_copy(copy_destination, source) != 0) return 1;
  return close_value(storage[0], 42.0) && close_value(storage[1], 42.0) &&
      close_value(storage[2], 2.0) && close_value(storage[3], 200.0) ? 0 : 1;
}

static int sparse_fill_copy_vector(size_t lanes) {
  const size_t entries = 2u, groups = 2u, count = groups * entries * lanes;
  double *storage = malloc(count * sizeof(*storage));
  if (!storage) return 1;
  for (size_t lane = 0; lane < lanes; ++lane) {
    storage[lane] = (double)lane + 1.0;
    storage[lanes + lane] = 100.0 + lane;
  }
  storage[entries * lanes] = (double)lanes + 1.0;
  storage[entries * lanes + lanes] = 100.0 + lanes;
  chem_group_view destination, copied, partial;
  chem_const_group_view source;
  int failed = chem_group_view_init(storage, count, 0u, lanes, 1u, &destination) != 0 ||
      chem_group_view_fill(destination, 42.0) != 0 ||
      chem_group_view_init(storage, count, lanes, lanes, 1u, &copied) != 0 ||
      chem_const_group_view_init(storage, count, 0u, lanes, 1u, &source) != 0 ||
      chem_group_view_copy(copied, source) != 0;
  for (size_t lane = 0; lane < lanes && !failed; ++lane)
    if (!close_value(storage[lane], 42.0) || !close_value(storage[lanes + lane], 42.0)) failed = 1;
  if (!failed && (chem_group_view_init(storage, count, entries * lanes, 1u, 1u, &partial) != 0 ||
      chem_group_view_fill(partial, -7.0) != 0 || !close_value(storage[entries * lanes], -7.0))) failed = 1;
  free(storage);
  return failed;
}

static int cross_group_copy_and_fill(void) {
  double source_data[] = { 1, 2, 3, 4 }, destination_data[] = { 0, 0, 0, 0 };
  chem_group_view destination;
  chem_const_group_view source;
  if (chem_group_view_init(destination_data, 4u, 0u, 4u, 1u, &destination) != 0 ||
      chem_const_group_view_init(source_data, 4u, 0u, 4u, 1u, &source) != 0 ||
      chem_group_view_copy(destination, source) != 0) return 1;
  for (size_t i = 0; i < 4u; ++i)
    if (!close_value(destination_data[i], (double)i + 1.0)) return 1;
  return chem_group_view_fill(destination, 55.5) == 0 &&
      close_value(destination_data[0], 55.5) && close_value(destination_data[3], 55.5) ? 0 : 1;
}

typedef int (*test_function)(void);
typedef struct { const char *name; test_function run; } test_case;

static int dense_l1(void) { return dense_grouped_matches_raw(1u); }
static int dense_l2(void) { return dense_grouped_matches_raw(2u); }
static int dense_l4(void) { return dense_grouped_matches_raw(4u); }
static int dense_l8(void) { return dense_grouped_matches_raw(8u); }
static int sparse_mixed_l1(void) { return mixed_sparse_dense(1u); }
static int sparse_mixed_l2(void) { return mixed_sparse_dense(2u); }
static int sparse_mixed_l4(void) { return mixed_sparse_dense(4u); }
static int sparse_mixed_l8(void) { return mixed_sparse_dense(8u); }
static int vector_fill_l4_csr(void) { return sparse_fill_copy_vector(4u); }
static int vector_fill_l4_csc(void) { return sparse_fill_copy_vector(4u); }
static int vector_fill_l8_csr(void) { return sparse_fill_copy_vector(8u); }

int main(int argc, char **argv) {
  static const test_case tests[] = {
    { "MatrixColumnViewRoundTrip", dense_group_roundtrip },
    { "VectorMatrixColumnViewContiguousBlock", vector_matrix_column_block },
    { "DenseGroupedMatchesRawL1", dense_l1 }, { "DenseGroupedMatchesRawL2", dense_l2 },
    { "DenseGroupedMatchesRawL4", dense_l4 }, { "DenseGroupedMatchesRawL8", dense_l8 },
    { "PlainMatrixGroupedRoundTrip", plain_matrix_roundtrip },
    { "SparseStandardBlockViewRoundTrip", sparse_standard_roundtrip },
    { "SparseVectorBlockViewContiguousBlock", sparse_vector_contiguous },
    { "MixedSparseAndDenseGroupedL1", sparse_mixed_l1 }, { "MixedSparseAndDenseGroupedL2", sparse_mixed_l2 },
    { "MixedSparseAndDenseGroupedL4", sparse_mixed_l4 }, { "MixedSparseAndDenseGroupedL8", sparse_mixed_l8 },
    { "SparseVectorRowColOverloadStillFunctional", sparse_group_sum },
    { "ConstGroupViewProducesGroupedViews", const_views },
    { "FillAndCopySparseStandardCSR", sparse_fill_copy_standard },
    { "FillAndCopySparseStandardCSC", sparse_fill_copy_standard },
    { "FillAndCopySparseVectorCSR_L4", vector_fill_l4_csr },
    { "FillAndCopySparseVectorCSC_L4", vector_fill_l4_csc },
    { "FillAndCopySparseVectorCSR_L8", vector_fill_l8_csr },
    { "FillCopyCrossGroupViewSources", cross_group_copy_and_fill },
  };
  if (argc != 2) return 2;
  const char *name = strrchr(argv[1], '.');
  name = name ? name + 1 : argv[1];
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    if (strcmp(name, tests[i].name) == 0) return tests[i].run();
  fprintf(stderr, "unimplemented grouped view case: %s\n", argv[1]);
  return 2;
}
