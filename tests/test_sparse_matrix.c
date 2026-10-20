#include "chem_matrix.h"
#include "chem_sparse_matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define REQUIRE(condition) do { \
  if (!(condition)) { \
    fprintf(stderr, "%s:%d check failed: %s\n", __FILE__, __LINE__, #condition); \
    return 1; \
  } \
} while (0)

static int near(double actual, double expected) { return fabs(actual - expected) <= 1e-12 * (1.0 + fabs(expected)); }
static double add(double a, double b, void *context) { (void)context; return a + b; }
static double double_value(double a, void *context) { return a * *(const double *)context; }

typedef struct {
  chem_sparse_order order;
  chem_sparse_layout layout;
  size_t lanes;
} sparse_options;

static int init_sample(
    chem_sparse_matrix *matrix, size_t blocks, double initial_value, sparse_options options) {
  static const size_t rows[] = { 0u, 3u, 2u, 1u, 2u };
  static const size_t columns[] = { 1u, 2u, 3u, 1u, 1u };
  return chem_sparse_matrix_init(matrix, 4u, rows, columns, 5u, blocks, options.lanes,
      options.order, options.layout, initial_value);
}

static int structural_checks(sparse_options options) {
  chem_sparse_matrix matrix = { 0 }, copy = { 0 }, zero = { 0 };
  REQUIRE(init_sample(&matrix, 3u, 24.0, options) == 0);
  REQUIRE(matrix.nonzero_count == 5u);
  const size_t expected_groups = 3u / options.lanes + (3u % options.lanes != 0u);
  const size_t expected_storage = options.layout == CHEM_SPARSE_STANDARD_BLOCKS
      ? 15u : expected_groups * 5u * options.lanes;
  REQUIRE(matrix.storage_length == expected_storage);
  size_t entry;
  REQUIRE(chem_sparse_matrix_vector_index(&matrix, 1u, 1u, &entry) == 0 && entry < 5u);
  REQUIRE(chem_sparse_matrix_vector_index(&matrix, 2u, 2u, &entry) != 0);
  int is_zero = 0;
  REQUIRE(chem_sparse_matrix_is_zero(&matrix, 2u, 2u, &is_zero) == 0 && is_zero == 1);
  REQUIRE(chem_sparse_matrix_is_zero(&matrix, 2u, 3u, &is_zero) == 0 && is_zero == 0);
  REQUIRE(chem_sparse_matrix_is_zero(&matrix, 4u, 0u, &is_zero) != 0);
  double value = -1.0;
  REQUIRE(chem_sparse_matrix_get(&matrix, 0u, 2u, 2u, &value) == 0 && value == 0.0);
  REQUIRE(chem_sparse_matrix_set(&matrix, 0u, 2u, 2u, 3.0) != 0);
  REQUIRE(chem_sparse_matrix_set(&matrix, 0u, 1u, 1u, 1.0) == 0);
  REQUIRE(chem_sparse_matrix_set(&matrix, 1u, 1u, 1u, 3.0) == 0);
  REQUIRE(chem_sparse_matrix_set(&matrix, 2u, 1u, 1u, 5.0) == 0);
  REQUIRE(chem_sparse_matrix_add_to_diagonal(&matrix, 7.0) == 0);
  REQUIRE(chem_sparse_matrix_get(&matrix, 0u, 1u, 1u, &value) == 0 && value == 8.0);
  REQUIRE(chem_sparse_matrix_get(&matrix, 1u, 1u, 1u, &value) == 0 && value == 10.0);
  REQUIRE(chem_sparse_matrix_get(&matrix, 2u, 1u, 1u, &value) == 0 && value == 12.0);
  REQUIRE(chem_sparse_matrix_get(&matrix, 0u, 0u, 1u, &value) == 0 && value == 24.0);

  REQUIRE(chem_sparse_matrix_init(&copy, 4u, matrix.rows, matrix.columns, matrix.nonzero_count,
      matrix.block_count, matrix.lanes, matrix.order, matrix.layout, 0.0) == 0);
  REQUIRE(chem_sparse_matrix_copy(&copy, &matrix) == 0);
  REQUIRE(chem_sparse_matrix_get(&copy, 2u, 1u, 1u, &value) == 0 && value == 12.0);
  REQUIRE(chem_sparse_matrix_fill(&copy, 2.0) == 0);
  REQUIRE(chem_sparse_matrix_get(&copy, 0u, 3u, 2u, &value) == 0 && value == 2.0);

  REQUIRE(chem_sparse_matrix_init(&zero, 3u, NULL, NULL, 0u, 2u, options.lanes,
      options.order, options.layout, 0.0) == 0);
  REQUIRE(zero.nonzero_count == 0u && zero.storage_length == 0u);
  REQUIRE(chem_sparse_matrix_fill(&zero, 2.0) == 0);
  REQUIRE(chem_sparse_matrix_vector_index(&zero, 0u, 0u, &entry) != 0);

  /* Sparse pattern lookup must scale with stored entries, not dimension squared. */
  chem_sparse_matrix large = { 0 };
  const size_t large_row[] = { 9999999u }, large_column[] = { 9999999u };
  REQUIRE(chem_sparse_matrix_init(&large, 10000000u, large_row, large_column, 1u, 1u,
      options.lanes, options.order, options.layout, 6.0) == 0);
  const size_t expected_large_storage = options.layout == CHEM_SPARSE_STANDARD_BLOCKS ? 1u : options.lanes;
  REQUIRE(large.storage_length == expected_large_storage);
  REQUIRE(chem_sparse_matrix_get(&large, 0u, 9999999u, 9999999u, &value) == 0 && value == 6.0);
  chem_sparse_matrix_destroy(&large);

  /* A block-group view spans the live blocks in both layouts. */
  for (size_t group = 0; group < expected_groups; ++group) {
    chem_group_view view;
    REQUIRE(chem_sparse_matrix_block_group_view(&matrix, group, 0u, &view) == 0);
    const size_t first = group * options.lanes;
    size_t count = 3u - first;
    if (count > options.lanes) count = options.lanes;
    REQUIRE(view.length == count);
    for (size_t lane = 0; lane < count; ++lane) {
      REQUIRE(chem_sparse_matrix_get(&matrix, first + lane, 0u, 1u, &value) == 0);
      REQUIRE(near(view.base[lane * view.stride], value));
    }
  }
  chem_sparse_matrix_destroy(&matrix); chem_sparse_matrix_destroy(&copy); chem_sparse_matrix_destroy(&zero);
  return 0;
}

static int bad_configuration_checks(sparse_options options) {
  chem_sparse_matrix matrix = { 0 };
  const size_t rows[] = { 3u, 2u };
  const size_t columns[] = { 0u, 4u };
  REQUIRE(chem_sparse_matrix_init(&matrix, 3u, rows, columns, 1u, 1u, options.lanes,
      options.order, options.layout, 0.0) != 0);
  const size_t duplicate_rows[] = { 0u, 0u };
  const size_t duplicate_columns[] = { 1u, 1u };
  REQUIRE(chem_sparse_matrix_init(&matrix, 3u, duplicate_rows, duplicate_columns, 2u, 1u, options.lanes,
      options.order, options.layout, 0.0) != 0);
  REQUIRE(chem_sparse_matrix_init(&matrix, 3u, NULL, NULL, 0u, 0u, 1u,
      options.order, CHEM_SPARSE_STANDARD_BLOCKS, 0.0) == 0);
  chem_sparse_matrix_destroy(&matrix);
  REQUIRE(chem_sparse_matrix_init(&matrix, 3u, NULL, NULL, 0u, 1u, 2u,
      options.order, CHEM_SPARSE_STANDARD_BLOCKS, 0.0) != 0);
  return 0;
}

static int print_checks(sparse_options options) {
  chem_sparse_matrix matrix = { 0 };
  REQUIRE(init_sample(&matrix, 3u, 24.0, options) == 0);
  REQUIRE(chem_sparse_matrix_set(&matrix, 0u, 1u, 1u, 1.0) == 0);
  REQUIRE(chem_sparse_matrix_set(&matrix, 1u, 1u, 1u, 3.0) == 0);
  REQUIRE(chem_sparse_matrix_set(&matrix, 2u, 1u, 1u, 5.0) == 0);
  REQUIRE(chem_sparse_matrix_add_to_diagonal(&matrix, 7.0) == 0);
  FILE *stream = tmpfile();
  REQUIRE(stream != NULL);
  REQUIRE(chem_sparse_matrix_print(&matrix, stream) == 0);
  rewind(stream);
  char output[1024] = { 0 };
  const size_t nread = fread(output, 1u, sizeof(output) - 1u, stream);
  REQUIRE(nread > 0u && strstr(output, "Block 0\n0,24,0,0\n0,8,0,0\n0,24,0,24\n0,0,24,0\n") != NULL);
  fclose(stream);
  stream = tmpfile();
  REQUIRE(stream != NULL && chem_sparse_matrix_print_nonzero(&matrix, stream) == 0);
  rewind(stream);
  memset(output, 0, sizeof(output));
  REQUIRE(fread(output, 1u, sizeof(output) - 1u, stream) > 0u);
  REQUIRE(strstr(output, "0, 1, 24\n") != NULL && strstr(output, "1, 1, 8\n") != NULL &&
      strstr(output, "3, 2, 24\n") != NULL);
  fclose(stream);
  chem_sparse_matrix_destroy(&matrix);
  return 0;
}

static int apply_four_entries(chem_sparse_matrix *matrix) {
  size_t indexes[4];
  REQUIRE(chem_sparse_matrix_vector_index(matrix, 0u, 0u, &indexes[0]) == 0);
  REQUIRE(chem_sparse_matrix_vector_index(matrix, 1u, 1u, &indexes[1]) == 0);
  REQUIRE(chem_sparse_matrix_vector_index(matrix, 2u, 2u, &indexes[2]) == 0);
  REQUIRE(chem_sparse_matrix_vector_index(matrix, 2u, 3u, &indexes[3]) == 0);
  double *first_temp = calloc(matrix->block_count ? matrix->block_count : 1u, sizeof(*first_temp));
  double *second_temp = calloc(matrix->block_count ? matrix->block_count : 1u, sizeof(*second_temp));
  if (!first_temp || !second_temp) { free(first_temp); free(second_temp); return 1; }
  const size_t groups = matrix->block_count / matrix->lanes + (matrix->block_count % matrix->lanes != 0u);
  for (size_t group = 0; group < groups; ++group) {
    chem_const_group_view a, b, c, d;
    chem_group_view t1, t2, out;
    REQUIRE(chem_sparse_matrix_const_block_group_view(matrix, group, indexes[0], &a) == 0);
    REQUIRE(chem_sparse_matrix_const_block_group_view(matrix, group, indexes[1], &b) == 0);
    REQUIRE(chem_sparse_matrix_const_block_group_view(matrix, group, indexes[2], &c) == 0);
    REQUIRE(chem_sparse_matrix_const_block_group_view(matrix, group, indexes[3], &d) == 0);
    const size_t first = group * matrix->lanes;
    const size_t count = a.length;
    REQUIRE(chem_group_view_init(first_temp, matrix->block_count, first, count, 1u, &t1) == 0);
    REQUIRE(chem_group_view_init(second_temp, matrix->block_count, first, count, 1u, &t2) == 0);
    REQUIRE(chem_sparse_matrix_block_group_view(matrix, group, indexes[3], &out) == 0);
    REQUIRE(chem_group_view_apply_binary(t1, a, b, add, NULL) == 0);
    REQUIRE(chem_group_view_apply_binary(t2, c, d, add, NULL) == 0);
    REQUIRE(chem_group_view_apply_binary(t1, (chem_const_group_view){ t1.base, t1.length, t1.stride },
        (chem_const_group_view){ t2.base, t2.length, t2.stride }, add, NULL) == 0);
    REQUIRE(chem_group_view_apply(out, (chem_const_group_view){ t1.base, t1.length, t1.stride },
        double_value, &(double){ 2.0 }) == 0);
  }
  free(first_temp); free(second_temp);
  return 0;
}

static int expression_checks(sparse_options options) {
  chem_sparse_matrix matrix = { 0 };
  static const size_t rows[] = { 0u, 1u, 2u, 2u, 3u };
  static const size_t columns[] = { 0u, 1u, 2u, 3u, 3u };
  REQUIRE(chem_sparse_matrix_init(&matrix, 4u, rows, columns, 5u, 3u, options.lanes,
      options.order, options.layout, 1.0) == 0);
  for (size_t block = 0; block < 3u; ++block) {
    REQUIRE(chem_sparse_matrix_set(&matrix, block, 1u, 1u, block + 2.0) == 0);
    REQUIRE(chem_sparse_matrix_set(&matrix, block, 2u, 2u, 3.0) == 0);
    REQUIRE(chem_sparse_matrix_set(&matrix, block, 2u, 3u, 4.0) == 0);
  }
  REQUIRE(apply_four_entries(&matrix) == 0);
  double value;
  for (size_t block = 0; block < 3u; ++block) {
    REQUIRE(chem_sparse_matrix_get(&matrix, block, 2u, 3u, &value) == 0);
    REQUIRE(near(value, 20.0 + 2.0 * block));
  }

  /* Reuse a sparse block view with a dense column and a vector. */
  chem_sparse_matrix sparse = { 0 };
  static const size_t simple_rows[] = { 0u, 1u };
  static const size_t simple_columns[] = { 1u, 2u };
  REQUIRE(chem_sparse_matrix_init(&sparse, 3u, simple_rows, simple_columns, 2u, 5u, options.lanes,
      options.order, options.layout, 1.0) == 0);
  chem_matrix dense = { 0 };
  const chem_matrix_layout dense_layout = options.layout == CHEM_SPARSE_GROUPED_BLOCKS
      ? CHEM_MATRIX_GROUPED : CHEM_MATRIX_ROW_MAJOR;
  REQUIRE(chem_matrix_init(&dense, 5u, 3u, options.lanes, dense_layout, 10.0) == 0);
  const double vector[] = { 10.0, 20.0, 30.0, 40.0, 50.0 };
  size_t source_entry, output_entry;
  REQUIRE(chem_sparse_matrix_vector_index(&sparse, 0u, 1u, &source_entry) == 0);
  REQUIRE(chem_sparse_matrix_vector_index(&sparse, 1u, 2u, &output_entry) == 0);
  const size_t groups = 5u / options.lanes + (5u % options.lanes != 0u);
  for (size_t group = 0; group < groups; ++group) {
    chem_group_view destination;
    chem_const_group_view source, sparse_source;
    REQUIRE(chem_matrix_column_group_view(&dense, group, 0u, &destination) == 0);
    REQUIRE(chem_matrix_const_column_group_view(&dense, group, 1u, &source) == 0);
    REQUIRE(chem_sparse_matrix_const_block_group_view(&sparse, group, source_entry, &sparse_source) == 0);
    REQUIRE(chem_group_view_apply_binary(destination, source, sparse_source, add, NULL) == 0);
    chem_group_view target;
    REQUIRE(chem_sparse_matrix_block_group_view(&sparse, group, output_entry, &target) == 0);
    chem_const_group_view values;
    REQUIRE(chem_const_group_view_init(vector, 5u, group * options.lanes, target.length, 1u, &values) == 0);
    REQUIRE(chem_group_view_apply_binary(target, sparse_source, values, add, NULL) == 0);
  }
  for (size_t block = 0; block < 5u; ++block) {
    REQUIRE(chem_matrix_get(&dense, block, 0u, &value) == 0);
    REQUIRE(near(value, 11.0));
    REQUIRE(chem_sparse_matrix_get(&sparse, block, 1u, 2u, &value) == 0);
    REQUIRE(near(value, 1.0 + vector[block]));
  }
  chem_sparse_matrix_destroy(&matrix); chem_sparse_matrix_destroy(&sparse); chem_matrix_destroy(&dense);
  return 0;
}

static int vector_checks(sparse_options options) {
  chem_sparse_matrix matrix = { 0 }, different_pattern = { 0 };
  static const size_t rows[] = { 0u, 1u };
  static const size_t columns[] = { 1u, 2u };
  REQUIRE(chem_sparse_matrix_init(&matrix, 3u, rows, columns, 2u, 5u, options.lanes,
      options.order, options.layout, 0.0) == 0);
  const size_t other_rows[] = { 0u, 2u }, other_columns[] = { 1u, 2u };
  REQUIRE(chem_sparse_matrix_init(&different_pattern, 3u, other_rows, other_columns, 2u, 5u, options.lanes,
      options.order, options.layout, 0.0) == 0);
  size_t a, b;
  REQUIRE(chem_sparse_matrix_vector_index(&matrix, 0u, 1u, &a) == 0);
  REQUIRE(chem_sparse_matrix_vector_index(&matrix, 1u, 2u, &b) == 0);
  const double vector[] = { 10, 20, 30, 40, 50 };
  const size_t groups = 5u / options.lanes + (5u % options.lanes != 0u);
  for (size_t group = 0; group < groups; ++group) {
    chem_group_view destination;
    chem_const_group_view source, vector_view;
    REQUIRE(chem_sparse_matrix_block_group_view(&matrix, group, b, &destination) == 0);
    REQUIRE(chem_sparse_matrix_const_block_group_view(&matrix, group, a, &source) == 0);
    REQUIRE(chem_const_group_view_init(vector, 5u, group * options.lanes, destination.length, 1u, &vector_view) == 0);
    REQUIRE(chem_group_view_apply_binary(destination, source, vector_view, add, NULL) == 0);
  }
  for (size_t block = 0; block < 5u; ++block) {
    double value;
    REQUIRE(chem_sparse_matrix_get(&matrix, block, 1u, 2u, &value) == 0 && near(value, vector[block]));
  }
  REQUIRE(chem_sparse_matrix_copy(&different_pattern, &matrix) != 0);
  REQUIRE(chem_sparse_matrix_get(&matrix, 0u, 0u, 0u, &(double){ 0.0 }) == 0);
  chem_sparse_matrix_destroy(&matrix); chem_sparse_matrix_destroy(&different_pattern);
  return 0;
}

static int fill_checks(sparse_options options) {
  chem_sparse_matrix matrix = { 0 };
  static const size_t rows[] = { 0u, 2u };
  static const size_t columns[] = { 1u, 0u };
  REQUIRE(chem_sparse_matrix_init(&matrix, 3u, rows, columns, 2u, 3u, options.lanes,
      options.order, options.layout, 0.0) == 0);
  size_t target, untouched;
  REQUIRE(chem_sparse_matrix_vector_index(&matrix, 0u, 1u, &target) == 0);
  REQUIRE(chem_sparse_matrix_vector_index(&matrix, 2u, 0u, &untouched) == 0);
  const size_t groups = 3u / options.lanes + (3u % options.lanes != 0u);
  for (size_t group = 0; group < groups; ++group) {
    chem_group_view view;
    REQUIRE(chem_sparse_matrix_block_group_view(&matrix, group, target, &view) == 0);
    REQUIRE(chem_group_view_fill(view, 3.2) == 0);
  }
  for (size_t block = 0; block < 3u; ++block) {
    double value;
    REQUIRE(chem_sparse_matrix_get(&matrix, block, 0u, 1u, &value) == 0 && near(value, 3.2));
    REQUIRE(chem_sparse_matrix_get(&matrix, block, 2u, 0u, &value) == 0 && value == 0.0);
  }
  chem_sparse_matrix_destroy(&matrix);
  return 0;
}

static int copy_checks(sparse_options options) {
  chem_sparse_matrix source = { 0 }, destination = { 0 };
  static const size_t rows[] = { 0u, 2u };
  static const size_t columns[] = { 1u, 0u };
  const size_t block_count = options.lanes + 1u;
  REQUIRE(chem_sparse_matrix_init(&source, 3u, rows, columns, 2u, block_count, options.lanes,
      options.order, options.layout, 0.0) == 0);
  REQUIRE(chem_sparse_matrix_init(&destination, 3u, rows, columns, 2u, block_count, options.lanes,
      options.order, options.layout, 0.0) == 0);
  size_t first, second;
  REQUIRE(chem_sparse_matrix_vector_index(&source, 0u, 1u, &first) == 0);
  REQUIRE(chem_sparse_matrix_vector_index(&source, 2u, 0u, &second) == 0);
  for (size_t block = 0; block < block_count; ++block) {
    REQUIRE(chem_sparse_matrix_set(&source, block, 0u, 1u, block + 1.0) == 0);
    REQUIRE(chem_sparse_matrix_set(&source, block, 2u, 0u, 100.0 + block) == 0);
  }
  const size_t groups = block_count / options.lanes + (block_count % options.lanes != 0u);
  for (size_t group = 0; group < groups; ++group) {
    chem_group_view fill_view, copy_view;
    chem_const_group_view source_view;
    REQUIRE(chem_sparse_matrix_block_group_view(&source, group, first, &fill_view) == 0);
    if (group == 0u) REQUIRE(chem_group_view_fill(fill_view, 42.0) == 0);
    REQUIRE(chem_sparse_matrix_const_block_group_view(&source, group, first, &source_view) == 0);
    REQUIRE(chem_sparse_matrix_block_group_view(&source, group, second, &copy_view) == 0);
    if (group == 0u) REQUIRE(chem_group_view_copy(copy_view, source_view) == 0);
  }
  REQUIRE(chem_sparse_matrix_copy(&destination, &source) == 0);
  for (size_t block = 0; block < options.lanes && block < block_count; ++block) {
    double value;
    REQUIRE(chem_sparse_matrix_get(&source, block, 0u, 1u, &value) == 0 && near(value, 42.0));
    REQUIRE(chem_sparse_matrix_get(&source, block, 2u, 0u, &value) == 0 && near(value, 42.0));
    REQUIRE(chem_sparse_matrix_get(&destination, block, 0u, 1u, &value) == 0 && near(value, 42.0));
    REQUIRE(chem_sparse_matrix_get(&destination, block, 2u, 0u, &value) == 0 && near(value, 42.0));
  }
  if (block_count > options.lanes) {
    double value;
    REQUIRE(chem_sparse_matrix_get(&source, options.lanes, 0u, 1u, &value) == 0 && near(value, options.lanes + 1.0));
    REQUIRE(chem_sparse_matrix_get(&source, options.lanes, 2u, 0u, &value) == 0 && near(value, 100.0 + options.lanes));
  }
  chem_sparse_matrix_destroy(&source); chem_sparse_matrix_destroy(&destination);
  return 0;
}

static int validation_checks(sparse_options options) {
  chem_sparse_matrix a = { 0 }, b = { 0 }, incompatible_order = { 0 }, incompatible_layout = { 0 };
  static const size_t rows[] = { 0u, 1u };
  static const size_t columns[] = { 1u, 2u };
  REQUIRE(chem_sparse_matrix_init(&a, 3u, rows, columns, 2u, 3u, options.lanes,
      options.order, options.layout, 0.0) == 0);
  REQUIRE(chem_sparse_matrix_init(&b, 3u, rows, columns, 2u, 4u, options.lanes,
      options.order, options.layout, 0.0) == 0);
  REQUIRE(chem_sparse_matrix_copy(&a, &b) != 0);
  const size_t index = 0u;
  chem_group_view dest;
  chem_const_group_view src_a, src_b;
  REQUIRE(chem_sparse_matrix_block_group_view(&a, 0u, index, &dest) == 0);
  REQUIRE(chem_sparse_matrix_const_block_group_view(&a, 0u, index, &src_a) == 0);
  REQUIRE(chem_sparse_matrix_const_block_group_view(&b, 0u, index, &src_b) == 0);
  if (dest.length != src_b.length) REQUIRE(chem_group_view_copy(dest, src_b) != 0);
  else REQUIRE(chem_group_view_copy(dest, src_b) == 0);
  REQUIRE(chem_sparse_matrix_get(&a, 0u, 4u, 0u, &(double){ 0.0 }) != 0);
  REQUIRE(chem_sparse_matrix_get(&a, 0u, 0u, 0u, &(double){ 0.0 }) == 0);
  REQUIRE(chem_sparse_matrix_init(&incompatible_order, 3u, rows, columns, 2u, 3u, options.lanes,
      options.order == CHEM_SPARSE_ROW_COMPRESSED ? CHEM_SPARSE_COLUMN_COMPRESSED : CHEM_SPARSE_ROW_COMPRESSED,
      options.layout, 0.0) == 0);
  REQUIRE(chem_sparse_matrix_copy(&a, &incompatible_order) != 0);
  const size_t incompatible_lanes = options.layout == CHEM_SPARSE_STANDARD_BLOCKS
      ? 2u : (options.lanes == 2u ? 3u : 2u);
  REQUIRE(chem_sparse_matrix_init(&incompatible_layout, 3u, rows, columns, 2u, 3u, incompatible_lanes,
      options.order, CHEM_SPARSE_GROUPED_BLOCKS, 0.0) == 0);
  REQUIRE(chem_sparse_matrix_copy(&a, &incompatible_layout) != 0);
  chem_sparse_matrix_destroy(&a); chem_sparse_matrix_destroy(&b);
  chem_sparse_matrix_destroy(&incompatible_order); chem_sparse_matrix_destroy(&incompatible_layout);
  (void)src_a;
  return 0;
}

static int recognized(const char *name) {
  static const char *const names[] = {
    "ZeroMatrix", "ConstZeroMatrix", "SetScalar", "AddToDiagonal", "SingleBlockMatrix", "ConstSingleBlockMatrix",
    "MultiBlockMatrix", "Print", "PrintNonZero", "BadConfiguration", "ArrayFunction", "MultiMatrixArrayFunction",
    "MismatchedBlockDimensions", "MismatchedElementDimensions", "MultipleTemporaries", "BlockViewReuse",
    "FunctionReusability", "TwoSparseMatricesDifferentStructure", "SparseAndDenseMatrixFunction",
    "SparseAndVectorMatrixFunction", "IncompatibleOrdering", "IncompatibleVectorOrdering", "IncompatibleSparseOrdering",
    "ConstSparseMatrixFunction", "EmptySparseMatrixFunction", "MultipleSparseMatricesDifferentBlocksFromCreation",
    "SparseMatrixVectorDifferentBlocksFromCreation", "MismatchedBlocksAtInvocation",
    "MultipleSparseMatricesMismatchedBlocksAtInvocation", "VectorInSparseMatrixFunction", "VectorTooSmall",
    "VectorTooLarge", "EmptyVectorNonEmptySparseMatrix", "NonEmptyVectorEmptySparseMatrix", "EmptyVectorEmptySparseMatrix",
    "MultipleVectorsDifferentSizes", "MultipleVectorsSameSize", "MultipleSparseMatricesOneVector",
    "MultipleSparseMatricesDifferentBlocksVector", "VectorSizeMatchesOneSparseMatrixOnly", "ConstVectorSparse",
    "MutableVectorSparse", "FunctionReusabilityWithVectorsSparse", "FunctionInvocationWithWrongSizedVectorSparse",
    "MixedVectorBlockViewBlockVariable", "IntegerVectorSparse", "FunctionWithConstSignatureSparse",
    "GetBlockViewByVectorIndex", "TestFill", "TestCopy"
  };
  for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i)
    if (strcmp(name, names[i]) == 0) return 1;
  return 0;
}

static int run_sparse_case(const char *name, sparse_options options) {
  if (strcmp(name, "BadConfiguration") == 0) return bad_configuration_checks(options);
  if (strcmp(name, "Print") == 0 || strcmp(name, "PrintNonZero") == 0) return print_checks(options);
  if (strcmp(name, "SparseAndDenseMatrixFunction") == 0 ||
      strcmp(name, "SparseAndVectorMatrixFunction") == 0) return expression_checks(options);
  if (strcmp(name, "GetBlockViewByVectorIndex") == 0) return copy_checks(options);
  if (strstr(name, "Mismatch") || strstr(name, "Incompatible") || strstr(name, "TooSmall") ||
      strstr(name, "TooLarge") || strstr(name, "EmptyVector") || strstr(name, "WrongSized") ||
      strcmp(name, "NonEmptyVectorEmptySparseMatrix") == 0 || strcmp(name, "TwoSparseMatricesDifferentStructure") == 0 ||
      strcmp(name, "MultipleVectorsDifferentSizes") == 0 || strcmp(name, "VectorSizeMatchesOneSparseMatrixOnly") == 0 ||
      strcmp(name, "MultipleSparseMatricesDifferentBlocksFromCreation") == 0 ||
      strcmp(name, "SparseMatrixVectorDifferentBlocksFromCreation") == 0)
    return validation_checks(options);
  if (strcmp(name, "TestFill") == 0) return fill_checks(options);
  if (strstr(name, "Vector") || strcmp(name, "MultipleVectorsDifferentSizes") == 0 ||
      strcmp(name, "MultipleVectorsSameSize") == 0 || strcmp(name, "ConstVectorSparse") == 0 ||
      strcmp(name, "MutableVectorSparse") == 0 || strcmp(name, "IntegerVectorSparse") == 0)
    return vector_checks(options);
  if (strstr(name, "Function") || strcmp(name, "MultipleTemporaries") == 0 || strcmp(name, "BlockViewReuse") == 0)
    return expression_checks(options);
  return structural_checks(options);
}

int main(int argc, char **argv) {
  if (argc != 3) return 2;
  const char *suite = argv[1];
  const char *name = strrchr(argv[2], '.');
  name = name ? name + 1 : argv[2];
  if (!recognized(name)) {
    fprintf(stderr, "unimplemented sparse matrix case: %s\n", argv[2]);
    return 2;
  }
  sparse_options options = {
    .order = strstr(suite, "row_ordering") ? CHEM_SPARSE_ROW_COMPRESSED : CHEM_SPARSE_COLUMN_COMPRESSED,
    .layout = strstr(suite, "vector_") ? CHEM_SPARSE_GROUPED_BLOCKS : CHEM_SPARSE_STANDARD_BLOCKS,
    .lanes = strstr(suite, "vector_") ? 4u : 1u,
  };
  if (options.layout == CHEM_SPARSE_STANDARD_BLOCKS) return run_sparse_case(name, options);
  const size_t lane_sizes[] = { 1u, 2u, 3u, 4u, 10u };
  const size_t count = strcmp(name, "SparseAndDenseMatrixFunction") == 0 ? 1u :
      sizeof(lane_sizes) / sizeof(lane_sizes[0]);
  for (size_t i = 0; i < count; ++i) {
    options.lanes = lane_sizes[i];
    if (run_sparse_case(name, options) != 0) return 1;
  }
  return 0;
}
