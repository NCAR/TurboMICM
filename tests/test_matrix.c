#include "chem_matrix.h"

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

static int near(double actual, double expected) {
  return fabs(actual - expected) <= 1e-10 * (1.0 + fabs(expected));
}

static double add(double a, double b, void *context) { (void)context; return a + b; }
static double subtract(double a, double b, void *context) { (void)context; return a - b; }
static double multiply(double a, double b, void *context) { (void)context; return a * b; }
static double double_value(double a, void *context) { return a * *(double *)context; }

static int storage_checks(chem_matrix_layout layout, size_t lanes) {
  chem_matrix matrix = { 0 }, initialized = { 0 }, empty_rows = { 0 }, empty_columns = { 0 };
  REQUIRE(chem_matrix_init(&matrix, 3u, 5u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_set(&matrix, 1u, 3u, 64.7) == 0);
  REQUIRE(chem_matrix_set(&matrix, 0u, 0u, 41.2) == 0);
  REQUIRE(chem_matrix_set(&matrix, 2u, 4u, 102.3) == 0);
  double value = 0.0;
  REQUIRE(chem_matrix_get(&matrix, 1u, 3u, &value) == 0 && near(value, 64.7));
  REQUIRE(chem_matrix_get(&matrix, 0u, 0u, &value) == 0 && near(value, 41.2));
  REQUIRE(chem_matrix_get(&matrix, 2u, 4u, &value) == 0 && near(value, 102.3));
  if (layout == CHEM_MATRIX_ROW_MAJOR) {
    REQUIRE(matrix.storage_length == 15u && near(matrix.data[0], 41.2));
    REQUIRE(near(matrix.data[8], 64.7) && near(matrix.data[14], 102.3));
  } else {
    const size_t i00 = (0u * 5u + 0u) * lanes;
    const size_t i13 = (1u / lanes * 5u + 3u) * lanes + 1u % lanes;
    const size_t i24 = (2u / lanes * 5u + 4u) * lanes + 2u % lanes;
    const size_t groups = 3u / lanes + (3u % lanes != 0u);
    REQUIRE(matrix.storage_length == groups * 5u * lanes);
    REQUIRE(near(matrix.data[i00], 41.2) && near(matrix.data[i13], 64.7) && near(matrix.data[i24], 102.3));
  }
  chem_group_view row_view;
  chem_const_group_view const_row_view;
  REQUIRE(chem_matrix_row_view(&matrix, 1u, &row_view) == 0);
  REQUIRE(chem_matrix_const_row_view(&matrix, 1u, &const_row_view) == 0);
  REQUIRE(row_view.length == 5u && row_view.stride == (layout == CHEM_MATRIX_ROW_MAJOR ? 1u : lanes));
  REQUIRE(near(row_view.base[3u * row_view.stride], 64.7));
  REQUIRE(near(const_row_view.base[3u * const_row_view.stride], 64.7));

  const double row_values[] = { 12.3, 15.1, 24.3, 88.0, -0.3 };
  REQUIRE(chem_matrix_assign_row(&matrix, 2u, row_values, 5u) == 0);
  REQUIRE(chem_matrix_assign_row(&matrix, 2u, row_values, 2u) != 0);
  double copied_row[5];
  REQUIRE(chem_matrix_read_row(&matrix, 2u, copied_row, 5u) == 0);
  REQUIRE(near(copied_row[0], 12.3) && near(copied_row[2], 24.3));
  chem_group_view column_group;
  const size_t group_count = 3u / lanes + (3u % lanes != 0u);
  REQUIRE(chem_matrix_column_group_view(&matrix, group_count - 1u, 3u, &column_group) == 0);
  REQUIRE(column_group.length == 3u - (group_count - 1u) * lanes);

  REQUIRE(chem_matrix_init(&initialized, 2u, 3u, lanes, layout, 12.4) == 0);
  REQUIRE(chem_matrix_get(&initialized, 0u, 0u, &value) == 0 && near(value, 12.4));
  REQUIRE(chem_matrix_get(&initialized, 1u, 2u, &value) == 0 && near(value, 12.4));
  REQUIRE(chem_matrix_init(&empty_rows, 0u, 3u, lanes, layout, 1.0) == 0);
  REQUIRE(chem_matrix_init(&empty_columns, 3u, 0u, lanes, layout, 1.0) == 0);
  REQUIRE(empty_rows.rows == 0u && empty_rows.columns == 3u && empty_columns.columns == 0u);
  REQUIRE(chem_matrix_fill(&empty_rows, 2.0) == 0 && chem_matrix_fill(&empty_columns, 2.0) == 0);
  chem_matrix_destroy(&matrix); chem_matrix_destroy(&initialized);
  chem_matrix_destroy(&empty_rows); chem_matrix_destroy(&empty_columns);
  return 0;
}

static int arithmetic_checks(chem_matrix_layout layout, size_t lanes) {
  chem_matrix matrix = { 0 }, other = { 0 }, destination = { 0 }, source = { 0 }, second = { 0 }, temp_matrix = { 0 };
  REQUIRE(chem_matrix_init(&matrix, 4u, 3u, lanes, layout, 100.0) == 0);
  REQUIRE(chem_matrix_init(&other, 4u, 3u, lanes, layout, 200.0) == 0);
  double expected_sum = 0.0;
  for (size_t row = 0; row < 4u; ++row) {
    for (size_t column = 0; column < 3u; ++column) {
      const double y = row * 10.3 + column * 100.5;
      const double x = row * 1.7 + column * 10.2;
      REQUIRE(chem_matrix_set(&matrix, row, column, y) == 0);
      REQUIRE(chem_matrix_set(&other, row, column, x) == 0);
      expected_sum += y + 1.39 * x;
    }
  }
  REQUIRE(chem_matrix_axpy(&matrix, 1.39, &other) == 0);
  double actual_sum = 0.0;
  for (size_t row = 0; row < 4u; ++row)
    for (size_t column = 0; column < 3u; ++column) {
      double value;
      REQUIRE(chem_matrix_get(&matrix, row, column, &value) == 0);
      actual_sum += value;
    }
  REQUIRE(near(actual_sum, expected_sum));

  REQUIRE(chem_matrix_init(&destination, 5u, 3u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_init(&source, 5u, 2u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_init(&second, 5u, 4u, lanes, layout, 0.0) == 0);
  for (size_t row = 0; row < 5u; ++row) {
    REQUIRE(chem_matrix_set(&source, row, 0u, (double)row - 2.0) == 0);
    REQUIRE(chem_matrix_set(&source, row, 1u, (double)row + 8.0) == 0);
    REQUIRE(chem_matrix_set(&second, row, 0u, (double)row * 2.0) == 0);
    REQUIRE(chem_matrix_set(&second, row, 3u, (double)row * 4.0) == 0);
  }
  REQUIRE(chem_matrix_transform_columns(&destination, 2u, &source, 0u, &second, 3u, add, NULL) == 0);
  for (size_t row = 0; row < 5u; ++row) {
    double value;
    REQUIRE(chem_matrix_get(&destination, row, 2u, &value) == 0);
    REQUIRE(near(value, 5.0 * row - 2.0));
  }
  REQUIRE(chem_matrix_transform_columns(&destination, 0u, &source, 0u, &second, 3u, subtract, NULL) == 0);
  REQUIRE(chem_matrix_transform_columns(&destination, 1u, &source, 0u, &second, 3u, multiply, NULL) == 0);

  chem_matrix temporaries[2] = { { 0 }, { 0 } };
  REQUIRE(chem_matrix_init(&temporaries[0], 4u, 1u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_init(&temporaries[1], 4u, 1u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_init(&temp_matrix, 4u, 5u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_fill(&temp_matrix, 0.0) == 0);
  for (size_t row = 0; row < 4u; ++row) {
    REQUIRE(chem_matrix_set(&temp_matrix, row, 0u, row + 1.0) == 0);
    REQUIRE(chem_matrix_set(&temp_matrix, row, 1u, 10.0 * (row + 1.0)) == 0);
  }
  REQUIRE(chem_matrix_transform_columns(&temporaries[0], 0u, &temp_matrix, 0u, &temp_matrix, 1u, multiply, NULL) == 0);
  REQUIRE(chem_matrix_transform_columns(&temporaries[1], 0u, &temp_matrix, 0u, &temp_matrix, 1u, add, NULL) == 0);
  REQUIRE(chem_matrix_transform_columns(&temp_matrix, 2u, &temporaries[0], 0u, &temporaries[1], 0u, add, NULL) == 0);
  REQUIRE(chem_matrix_transform_columns(&temp_matrix, 3u, &temporaries[0], 0u, &temporaries[1], 0u, subtract, NULL) == 0);
  REQUIRE(chem_matrix_transform_columns(&temp_matrix, 4u, &temporaries[0], 0u, &temporaries[1], 0u, multiply, NULL) == 0);
  double value;
  REQUIRE(chem_matrix_get(&temp_matrix, 0u, 2u, &value) == 0 && near(value, 21.0));
  REQUIRE(chem_matrix_get(&temp_matrix, 3u, 3u, &value) == 0 && near(value, 116.0));
  REQUIRE(chem_matrix_get(&temp_matrix, 3u, 4u, &value) == 0 && near(value, 7040.0));

  chem_matrix clone = { 0 };
  REQUIRE(chem_matrix_init(&clone, 5u, 3u, lanes, layout, -1.0) == 0);
  REQUIRE(chem_matrix_transform_unary(&clone, &destination, double_value, &(double){ 2.0 }) == 0);
  REQUIRE(chem_matrix_get(&clone, 4u, 2u, &value) == 0 && near(value, 36.0));
  REQUIRE(chem_matrix_copy(&destination, &clone) == 0);
  REQUIRE(chem_matrix_transform_binary(&destination, &destination, &clone, add, NULL) == 0);
  REQUIRE(chem_matrix_get(&destination, 4u, 2u, &value) == 0 && near(value, 72.0));
  REQUIRE(chem_matrix_fill(&clone, 1.0) == 0);
  REQUIRE(chem_matrix_max_scalar(&clone, 2.0) == 0);
  REQUIRE(chem_matrix_get(&clone, 0u, 0u, &value) == 0 && near(value, 2.0));
  REQUIRE(chem_matrix_set(&clone, 1u, 1u, 3.0) == 0 && chem_matrix_min_scalar(&clone, 2.0) == 0);
  REQUIRE(chem_matrix_get(&clone, 1u, 1u, &value) == 0 && near(value, 2.0));
  double minimum, maximum;
  REQUIRE(chem_matrix_fill(&destination, 0.0) == 0);
  REQUIRE(chem_matrix_set(&destination, 0u, 0u, -3.0) == 0);
  REQUIRE(chem_matrix_set(&destination, 1u, 1u, 5.0) == 0);
  REQUIRE(chem_matrix_min_max(&destination, &minimum, &maximum) == 0);
  REQUIRE(near(minimum, -3.0) && near(maximum, 5.0));

  /* Matrix columns can be filled and copied through strided group views. */
  for (size_t group = 0; group < 5u / lanes + (5u % lanes != 0u); ++group) {
    chem_group_view target;
    chem_const_group_view from;
    REQUIRE(chem_matrix_column_group_view(&destination, group, 0u, &target) == 0);
    REQUIRE(chem_matrix_const_column_group_view(&destination, group, 2u, &from) == 0);
    REQUIRE(chem_group_view_copy(target, from) == 0);
    REQUIRE(chem_group_view_fill(target, 3.2) == 0);
  }
  for (size_t row = 0; row < 5u; ++row) {
    REQUIRE(chem_matrix_get(&destination, row, 0u, &value) == 0 && near(value, 3.2));
  }

  FILE *stream = tmpfile();
  REQUIRE(stream != NULL);
  chem_matrix printable = { 0 };
  REQUIRE(chem_matrix_init(&printable, 2u, 3u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_set(&printable, 1u, 1u, 3.0) == 0 && chem_matrix_print(&printable, stream) == 0);
  rewind(stream);
  char printed[64] = { 0 };
  REQUIRE(fread(printed, 1u, sizeof(printed) - 1u, stream) == 12u);
  REQUIRE(strcmp(printed, "0,0,0\n0,3,0\n") == 0);
  fclose(stream);

  chem_matrix_destroy(&matrix); chem_matrix_destroy(&other); chem_matrix_destroy(&destination);
  chem_matrix_destroy(&source); chem_matrix_destroy(&second); chem_matrix_destroy(&temporaries[0]);
  chem_matrix_destroy(&temporaries[1]); chem_matrix_destroy(&clone); chem_matrix_destroy(&printable);
  chem_matrix_destroy(&temp_matrix);
  return 0;
}

static void reducer_sum_square(double value, double *accumulator, void *context) {
  (void)context; *accumulator += value * value;
}
static void reducer_max(double value, double *accumulator, void *context) {
  (void)context; if (value > *accumulator) *accumulator = value;
}
static void reducer_or(double value, double *accumulator, void *context) {
  *accumulator = (*accumulator != 0.0 || value > *(double *)context) ? 1.0 : 0.0;
}
static void reducer_finite(double value, double *accumulator, void *context) {
  (void)context; *accumulator = (*accumulator != 0.0 && isfinite(value)) ? 1.0 : 0.0;
}
static void reducer_count(double value, double *accumulator, void *context) {
  (void)value; (void)context; *accumulator += 1.0;
}

static int reducer_checks(chem_matrix_layout layout, size_t lanes) {
  chem_matrix matrix = { 0 };
  REQUIRE(chem_matrix_init(&matrix, 3u, 2u, lanes, layout, 0.0) == 0);
  const double values[] = { 1, 2, 3, 4, 5, 6 };
  for (size_t row = 0; row < 3u; ++row)
    for (size_t column = 0; column < 2u; ++column)
      REQUIRE(chem_matrix_set(&matrix, row, column, values[row * 2u + column]) == 0);
  double sum = 0.0, maximum = -INFINITY, any_large = 0.0, all_finite = 1.0, count = 0.0, threshold = 10.0;
  REQUIRE(chem_matrix_reduce(&matrix, &sum, reducer_sum_square, NULL) == 0 && near(sum, 91.0));
  const double max_values[] = { 1, 8, 5, 2, 3, 4 };
  for (size_t row = 0; row < 3u; ++row)
    for (size_t column = 0; column < 2u; ++column)
      REQUIRE(chem_matrix_set(&matrix, row, column, max_values[row * 2u + column]) == 0);
  REQUIRE(chem_matrix_reduce(&matrix, &maximum, reducer_max, NULL) == 0 && near(maximum, 8.0));
  REQUIRE(chem_matrix_reduce(&matrix, &any_large, reducer_or, &threshold) == 0 && any_large == 0.0);
  REQUIRE(chem_matrix_set(&matrix, 2u, 1u, 42.0) == 0);
  REQUIRE(chem_matrix_reduce(&matrix, &any_large, reducer_or, &threshold) == 0 && any_large == 1.0);
  REQUIRE(chem_matrix_fill(&matrix, 1.0) == 0);
  REQUIRE(chem_matrix_reduce(&matrix, &all_finite, reducer_finite, NULL) == 0 && all_finite == 1.0);
  REQUIRE(chem_matrix_set(&matrix, 1u, 0u, NAN) == 0);
  REQUIRE(chem_matrix_reduce(&matrix, &all_finite, reducer_finite, NULL) == 0 && all_finite == 0.0);
  REQUIRE(chem_matrix_reduce(&matrix, &count, reducer_count, NULL) == 0 && count == 6.0);
  chem_matrix strict_matrix = { 0 };
  REQUIRE(chem_matrix_init(&strict_matrix, 3u, 1u, lanes, layout, 0.0) == 0);
  count = 0.0;
  REQUIRE(chem_matrix_reduce(&strict_matrix, &count, reducer_count, NULL) == 0 && count == 3.0);
  chem_matrix_destroy(&strict_matrix);
  chem_matrix_destroy(&matrix);
  return 0;
}

static int vector_checks(chem_matrix_layout layout, size_t lanes) {
  chem_matrix matrix = { 0 }, second = { 0 }, single_column = { 0 };
  REQUIRE(chem_matrix_init(&matrix, 5u, 3u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_init(&second, 5u, 3u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_init(&single_column, 5u, 1u, lanes, layout, 0.0) == 0);
  double vector[] = { 100, 200, 300, 400, 500 };
  for (size_t row = 0; row < 5u; ++row) {
    REQUIRE(chem_matrix_set(&matrix, row, 0u, (double)row - 2.0) == 0);
    REQUIRE(chem_matrix_set(&matrix, row, 1u, (double)row + 8.0) == 0);
    REQUIRE(chem_matrix_set(&matrix, row, 2u, (double)row + 18.0) == 0);
    REQUIRE(chem_matrix_set(&second, row, 0u, (double)row * 2.0) == 0);
  }
  REQUIRE(chem_matrix_transform_column_vector(&matrix, 2u, &matrix, 2u, vector, 5u, multiply, NULL) == 0);
  double value;
  REQUIRE(chem_matrix_get(&matrix, 0u, 2u, &value) == 0 && near(value, 1800.0));
  REQUIRE(chem_matrix_get(&matrix, 4u, 2u, &value) == 0 && near(value, 11000.0));
  REQUIRE(chem_matrix_transform_column_vector(&matrix, 2u, &matrix, 0u, vector, 3u, add, NULL) != 0);
  REQUIRE(chem_matrix_transform_column_vector(&single_column, 0u, &matrix, 0u, NULL, 0u, add, NULL) != 0);
  REQUIRE(chem_matrix_transform_column_vector(&matrix, 1u, &second, 0u, vector, 5u, add, NULL) == 0);
  REQUIRE(chem_matrix_transform_columns(&single_column, 0u, &matrix, 0u, &second, 0u, add, NULL) == 0);
  REQUIRE(chem_matrix_transform_column_vector(&second, 1u, &second, 0u, vector, 5u, add, NULL) == 0);

  chem_matrix vector_a = { 0 }, vector_b = { 0 }, vector_output = { 0 };
  REQUIRE(chem_matrix_init(&vector_a, 5u, 1u, lanes, layout, 2.0) == 0);
  REQUIRE(chem_matrix_init(&vector_b, 5u, 1u, lanes, layout, 3.0) == 0);
  REQUIRE(chem_matrix_init(&vector_output, 5u, 1u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_transform_binary(&vector_output, &vector_a, &vector_b, add, NULL) == 0);
  chem_matrix short_vector = { 0 };
  REQUIRE(chem_matrix_init(&short_vector, 3u, 1u, lanes, layout, 3.0) == 0);
  REQUIRE(chem_matrix_transform_columns(&vector_output, 0u, &vector_a, 0u, &short_vector, 0u, add, NULL) != 0);

  /* A C callback can be reused with matrices having different row counts. */
  chem_matrix small = { 0 }, short_source = { 0 };
  REQUIRE(chem_matrix_init(&small, 2u, 3u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_init(&short_source, 2u, 3u, lanes, layout, 4.0) == 0);
  REQUIRE(chem_matrix_transform_columns(&small, 1u, &small, 0u, &short_source, 2u, add, NULL) == 0);
  REQUIRE(chem_matrix_transform_columns(&matrix, 1u, &matrix, 0u, &short_source, 2u, add, NULL) != 0);
  REQUIRE(chem_matrix_transform_column_vector(&small, 1u, &small, 0u, vector, 5u, add, NULL) != 0);
  int index_vector[] = { 2, 0, 1 };
  int integer_values[] = { 10, 20, 30 };
  int capture_sum = 0;
  for (size_t i = 0; i < 3u; ++i) capture_sum += integer_values[index_vector[i]];
  REQUIRE(capture_sum == 60 && integer_values[index_vector[1]] == 10);
  chem_matrix_destroy(&matrix); chem_matrix_destroy(&second); chem_matrix_destroy(&single_column);
  chem_matrix_destroy(&small); chem_matrix_destroy(&short_source); chem_matrix_destroy(&vector_a);
  chem_matrix_destroy(&vector_b); chem_matrix_destroy(&vector_output); chem_matrix_destroy(&short_vector);
  return 0;
}

static int validation_checks(chem_matrix_layout layout, size_t lanes) {
  chem_matrix a = { 0 }, b = { 0 }, output = { 0 }, empty = { 0 };
  REQUIRE(chem_matrix_init(&a, 3u, 3u, lanes, layout, 1.0) == 0);
  REQUIRE(chem_matrix_init(&b, 4u, 3u, lanes, layout, 2.0) == 0);
  REQUIRE(chem_matrix_init(&output, 3u, 4u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_transform_columns(&output, 1u, &a, 0u, &b, 0u, add, NULL) != 0);
  chem_matrix wrong_columns = { 0 };
  REQUIRE(chem_matrix_init(&wrong_columns, 3u, 4u, lanes, layout, 0.0) == 0);
  REQUIRE(chem_matrix_transform_binary(&a, &a, &wrong_columns, add, NULL) != 0);
  REQUIRE(chem_matrix_transform_columns(&output, 9u, &a, 0u, &a, 1u, add, NULL) != 0);
  REQUIRE(chem_matrix_transform_columns(&output, 1u, &a, 5u, &a, 1u, add, NULL) != 0);
  REQUIRE(chem_matrix_init(&empty, 0u, 3u, lanes, layout, 0.0) == 0);
  double *no_vector = NULL;
  REQUIRE(chem_matrix_transform_column_vector(&empty, 0u, &empty, 0u, no_vector, 0u, add, NULL) == 0);
  REQUIRE(chem_matrix_transform_column_vector(&a, 1u, &a, 0u, no_vector, 0u, add, NULL) != 0);
  REQUIRE(chem_matrix_transform_binary(&a, &a, &b, add, NULL) != 0);
  REQUIRE(chem_matrix_row_view(&a, 9u, &(chem_group_view){ 0 }) != 0);
  REQUIRE(chem_matrix_column_group_view(&a, 10u, 0u, &(chem_group_view){ 0 }) != 0);
  chem_matrix_destroy(&a); chem_matrix_destroy(&b); chem_matrix_destroy(&output); chem_matrix_destroy(&empty);
  chem_matrix_destroy(&wrong_columns);
  return 0;
}

static int fill_copy_checks(chem_matrix_layout layout, size_t lanes) {
  chem_matrix matrix = { 0 };
  REQUIRE(chem_matrix_init(&matrix, 3u, 2u, lanes, layout, 0.0) == 0);
  for (size_t row = 0; row < 3u; ++row) {
    REQUIRE(chem_matrix_set(&matrix, row, 0u, (double)row + 1.0) == 0);
    REQUIRE(chem_matrix_set(&matrix, row, 1u, (double)row + 4.0) == 0);
  }
  const size_t groups = 3u / lanes + (3u % lanes != 0u);
  for (size_t group = 0; group < groups; ++group) {
    chem_group_view destination;
    chem_const_group_view source;
    REQUIRE(chem_matrix_column_group_view(&matrix, group, 0u, &destination) == 0);
    REQUIRE(chem_matrix_const_column_group_view(&matrix, group, 1u, &source) == 0);
    REQUIRE(chem_group_view_copy(destination, source) == 0);
    REQUIRE(chem_group_view_fill(destination, 42.0) == 0);
  }
  for (size_t row = 0; row < 3u; ++row) {
    double value;
    REQUIRE(chem_matrix_get(&matrix, row, 0u, &value) == 0 && near(value, 42.0));
    REQUIRE(chem_matrix_get(&matrix, row, 1u, &value) == 0 && near(value, (double)row + 4.0));
  }
  chem_matrix_destroy(&matrix);
  return 0;
}

static int recognized_case(const char *name) {
  static const char *const names[] = {
    "SmallMatrix", "SmallConstMatrix", "InitializeMatrix", "InitializeConstMatrix", "LoopOverMatrix",
    "LoopOverConstMatrix", "IterateOverMatrix", "Strides", "ConversionToVector", "ConstConversionToVector",
    "ConversionFromVector", "AssignmentFromVector", "Axpy", "ForEach", "SetScaler", "Max", "Min", "Print",
    "ArrayFunction", "MultiMatrixArrayFunction", "MismatchedRowDimensions", "MismatchedColumnDimensions",
    "WrongMatrixDimensions", "MultipleTemporaries", "ColumnViewReuse", "FunctionReusability", "ConstMatrixFunction",
    "EmptyMatrixFunction", "MultiMatrixDifferentRowsFromCreation", "MatrixVectorDifferentRowsFromCreation",
    "MismatchedRowsAtInvocation", "MultipleMatricesMismatchedRowsAtInvocation", "WrongColumnCountAtInvocation",
    "VectorInMatrixFunction", "VectorTooSmall", "VectorTooLarge", "EmptyVectorNonEmptyMatrix",
    "NonEmptyVectorEmptyMatrix", "EmptyVectorEmptyMatrix", "MultipleVectorsDifferentSizes", "MultipleVectorsSameSize",
    "MultipleMatricesOneVector", "MultipleMatricesDifferentRowsVector", "VectorSizeMatchesOneMatrixOnly", "ConstVector",
    "MutableVector", "FunctionReusabilityWithVectors", "FunctionInvocationWithWrongSizedVector",
    "MixedVectorColumnViewRowVariable", "IntegerVector", "FunctionWithConstSignature", "TestFill", "TestCopy",
    "ReduceSum", "ReduceMax", "ReduceLOr", "ReduceLAnd", "ReduceStrict", "VectorCapture"
  };
  for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i)
    if (strcmp(name, names[i]) == 0) return 1;
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 3 || (strcmp(argv[1], "matrix") != 0 && strcmp(argv[1], "vector_matrix") != 0)) return 2;
  const int vector_layout = strcmp(argv[1], "vector_matrix") == 0;
  const char *name = strrchr(argv[2], '.');
  name = name ? name + 1 : argv[2];
  if (vector_layout) {
    if (strcmp(name, "SmallVectorMatrix") == 0) name = "SmallMatrix";
    else if (strcmp(name, "SmallConstVectorMatrix") == 0) name = "SmallConstMatrix";
    else if (strcmp(name, "InitializeVectorMatrix") == 0) name = "InitializeMatrix";
    else if (strcmp(name, "InitializeConstVectorMatrix") == 0) name = "InitializeConstMatrix";
    else if (strcmp(name, "LoopOverVectorMatrix") == 0) name = "LoopOverMatrix";
    else if (strcmp(name, "LoopOverConstVectorMatrix") == 0) name = "LoopOverConstMatrix";
  }
  if (!recognized_case(name)) {
    fprintf(stderr, "unimplemented matrix case: %s\n", argv[2]);
    return 2;
  }
  const int basic = strstr(name, "Small") == name || strstr(name, "Initialize") == name ||
      strstr(name, "LoopOver") == name || strcmp(name, "IterateOverMatrix") == 0 ||
      strstr(name, "Conversion") == name || strcmp(name, "Strides") == 0 ||
      strcmp(name, "AssignmentFromVector") == 0 || strcmp(name, "Print") == 0;
  const int validation = strstr(name, "Mismatch") != NULL || strstr(name, "Wrong") != NULL ||
      strstr(name, "TooSmall") != NULL || strstr(name, "TooLarge") != NULL ||
      strstr(name, "EmptyVector") != NULL || strstr(name, "VectorSizeMatches") != NULL ||
      strstr(name, "InvocationWithWrong") != NULL || strcmp(name, "NonEmptyVectorEmptyMatrix") == 0;
  const int vectors = strstr(name, "Vector") != NULL || strcmp(name, "MultipleVectorsDifferentSizes") == 0 ||
      strcmp(name, "MultipleVectorsSameSize") == 0 || strcmp(name, "MultipleMatricesOneVector") == 0 ||
      strcmp(name, "VectorCapture") == 0;
  const int reducers = strstr(name, "Reduce") != NULL;
  const int fill_copy = strcmp(name, "TestFill") == 0 || strcmp(name, "TestCopy") == 0;
  const size_t lane_counts[] = { 1u, 2u, 3u, 4u };
  const size_t run_count = vector_layout ? sizeof(lane_counts) / sizeof(lane_counts[0]) : 1u;
  for (size_t run = 0; run < run_count; ++run) {
    const chem_matrix_layout layout = vector_layout ? CHEM_MATRIX_GROUPED : CHEM_MATRIX_ROW_MAJOR;
    const size_t lanes = vector_layout ? lane_counts[run] : 1u;
    int status = 0;
    if (validation) status = validation_checks(layout, lanes);
    else if (basic) status = storage_checks(layout, lanes);
    else if (vectors) status = vector_checks(layout, lanes);
    else if (reducers) status = reducer_checks(layout, lanes);
    else if (fill_copy) status = fill_copy_checks(layout, lanes);
    else status = arithmetic_checks(layout, lanes);
    if (status != 0) return status;
  }
  return 0;
}
