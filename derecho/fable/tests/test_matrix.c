/* MICM test/unit/util/test_matrix.cpp: MICM's Matrix is tmc_dense with L = 1.
 * The shared policy checks are in dense_policy.c; the layout checks MICM makes
 * on the returned matrices are here. */
#include "dense_policy.h"

#include "tmc_dense_ext.h"

/* a case that runs one policy check on a Matrix */
#define MATRIX_CASE(name, check)                                                              \
  TMC_TEST(name)                                                                              \
  {                                                                                           \
    return check(1);                                                                          \
  }

TMC_TEST(small_matrix)
{
  tmc_dense *m = dense_policy_small_matrix(1);
  TMC_EXPECT_EQ(m->size, (size_t)(3 * 5));
  TMC_EXPECT_EQ(m->data[0], 41.2);
  TMC_EXPECT_EQ(m->data[3 * 5 - 1], 102.3);
  TMC_EXPECT_EQ(m->data[1 * 5 + 3], 64.7);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(small_const_matrix)
{
  tmc_dense *storage = dense_policy_small_const_matrix(1);
  const tmc_dense *m = storage;
  TMC_EXPECT_EQ(m->size, (size_t)(3 * 5));
  TMC_EXPECT_EQ(m->data[0], 41.2);
  TMC_EXPECT_EQ(m->data[3 * 5 - 1], 102.3);
  TMC_EXPECT_EQ(m->data[1 * 5 + 3], 64.7);
  tmc_dense_free(storage);
  return TMC_TEST_PASS;
}

MATRIX_CASE(initialize_matrix, dense_policy_initialize_matrix)
MATRIX_CASE(initialize_const_matrix, dense_policy_initialize_const_matrix)
MATRIX_CASE(loop_over_matrix, dense_policy_loop_over_matrix)
MATRIX_CASE(loop_over_const_matrix, dense_policy_loop_over_const_matrix)

/* Iterating over a row's elements: they start at the row's first element and
 * follow one another at the column stride. MICM's Matrix<int> holds the same
 * integers as doubles here. */
TMC_TEST(iterate_over_matrix)
{
  tmc_dense *m = tmc_dense_create(3, 4, 1, 0);
  double *row = m->data + tmc_dense_index(m, 1, 0);
  size_t stride = tmc_dense_column_stride(m);
  int value = 42;
  for (int j = 0; j < m->n_cols; ++j)
    row[(size_t)j * stride] = value++;

  tmc_dense *copy = tmc_dense_copy(m);

  value = 42;
  for (int j = 0; j < m->n_cols; ++j)
    TMC_EXPECT_EQ(row[(size_t)j * stride], value++);
  TMC_EXPECT_EQ(tmc_dense_get(m, 0, 3), 0);
  TMC_EXPECT_EQ(tmc_dense_get(m, 1, 0), 42);
  TMC_EXPECT_EQ(tmc_dense_get(m, 1, 1), 43);
  TMC_EXPECT_EQ(tmc_dense_get(m, 1, 3), 45);
  TMC_EXPECT_EQ(tmc_dense_get(m, 2, 0), 0);
  for (int j = 0; j < 4; ++j)
    TMC_EXPECT_EQ(tmc_dense_get(copy, 1, j), 42 + j);

  tmc_dense_free(m);
  tmc_dense_free(copy);
  return TMC_TEST_PASS;
}

TMC_TEST(strides)
{
  tmc_dense *m = dense_policy_strides(1);
  TMC_EXPECT_EQ(tmc_dense_row_stride(m), (size_t)4);
  TMC_EXPECT_EQ(tmc_dense_column_stride(m), (size_t)1);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

MATRIX_CASE(conversion_to_vector, dense_policy_conversion_to_vector)
MATRIX_CASE(const_conversion_to_vector, dense_policy_const_conversion_to_vector)
MATRIX_CASE(conversion_from_vector, dense_policy_conversion_from_vector)
MATRIX_CASE(assignment_from_vector, dense_policy_assignment_from_vector)
MATRIX_CASE(axpy, dense_policy_axpy)
MATRIX_CASE(for_each, dense_policy_for_each)
MATRIX_CASE(set_scalar, dense_policy_set_scalar)
MATRIX_CASE(max, dense_policy_max)
MATRIX_CASE(min, dense_policy_min)
MATRIX_CASE(array_function, dense_policy_array_function)
MATRIX_CASE(multi_matrix_array_function, dense_policy_multi_matrix_array_function)
MATRIX_CASE(mismatched_row_dimensions, dense_policy_mismatched_row_dimensions)
MATRIX_CASE(mismatched_column_dimensions, dense_policy_mismatched_column_dimensions)
MATRIX_CASE(wrong_matrix_dimensions, dense_policy_wrong_matrix_dimensions)
MATRIX_CASE(multiple_temporaries, dense_policy_multiple_temporaries)
MATRIX_CASE(column_view_reuse, dense_policy_column_view_reuse)
MATRIX_CASE(function_reusability, dense_policy_function_reusability)
MATRIX_CASE(const_matrix_function, dense_policy_const_matrix_function)
MATRIX_CASE(empty_matrix_function, dense_policy_empty_matrix_function)
MATRIX_CASE(multi_matrix_different_rows_from_creation, dense_policy_multi_matrix_different_rows_from_creation)
MATRIX_CASE(matrix_vector_different_rows_from_creation, dense_policy_matrix_vector_different_rows_from_creation)
MATRIX_CASE(mismatched_rows_at_invocation, dense_policy_mismatched_rows_at_invocation)
MATRIX_CASE(multiple_matrices_mismatched_rows_at_invocation, dense_policy_multiple_matrices_mismatched_rows_at_invocation)
MATRIX_CASE(wrong_column_count_at_invocation, dense_policy_wrong_column_count_at_invocation)
MATRIX_CASE(print, dense_policy_print)
MATRIX_CASE(vector_in_matrix_function, dense_policy_vector_in_matrix_function)
MATRIX_CASE(vector_too_small, dense_policy_vector_too_small)
MATRIX_CASE(vector_too_large, dense_policy_vector_too_large)
MATRIX_CASE(empty_vector_non_empty_matrix, dense_policy_empty_vector_non_empty_matrix)
MATRIX_CASE(non_empty_vector_empty_matrix, dense_policy_non_empty_vector_empty_matrix)
MATRIX_CASE(empty_vector_empty_matrix, dense_policy_empty_vector_empty_matrix)
MATRIX_CASE(multiple_vectors_different_sizes, dense_policy_multiple_vectors_different_sizes)
MATRIX_CASE(multiple_vectors_same_size, dense_policy_multiple_vectors_same_size)
MATRIX_CASE(multiple_matrices_one_vector, dense_policy_multiple_matrices_one_vector)
MATRIX_CASE(multiple_matrices_different_rows_vector, dense_policy_multiple_matrices_different_rows_vector)
MATRIX_CASE(vector_length_matches_one_matrix_only, dense_policy_vector_length_matches_one_matrix_only)
MATRIX_CASE(const_vector, dense_policy_const_vector)
MATRIX_CASE(mutable_vector, dense_policy_mutable_vector)
MATRIX_CASE(function_reusability_with_vectors, dense_policy_function_reusability_with_vectors)
MATRIX_CASE(function_invocation_with_wrong_sized_vector, dense_policy_function_invocation_with_wrong_sized_vector)
MATRIX_CASE(mixed_vector_column_view_row_variable, dense_policy_mixed_vector_column_view_row_variable)
MATRIX_CASE(integer_vector, dense_policy_integer_vector)
MATRIX_CASE(function_with_const_signature, dense_policy_function_with_const_signature)
MATRIX_CASE(fill, dense_policy_fill)
MATRIX_CASE(copy, dense_policy_copy)
MATRIX_CASE(reduce_sum, dense_policy_reduce_sum)
MATRIX_CASE(reduce_max, dense_policy_reduce_max)
MATRIX_CASE(reduce_lor, dense_policy_reduce_lor)
MATRIX_CASE(reduce_land, dense_policy_reduce_land)
MATRIX_CASE(reduce_strict, dense_policy_reduce_strict)
MATRIX_CASE(vector_capture, dense_policy_vector_capture)

static const tmc_test_case cases[] = {
  TMC_CASE("matrix:Matrix.SmallMatrix", small_matrix),
  TMC_CASE("matrix:Matrix.SmallConstMatrix", small_const_matrix),
  TMC_CASE("matrix:Matrix.InitializeMatrix", initialize_matrix),
  TMC_CASE("matrix:Matrix.InitializeConstMatrix", initialize_const_matrix),
  TMC_CASE("matrix:Matrix.LoopOverMatrix", loop_over_matrix),
  TMC_CASE("matrix:Matrix.LoopOverConstMatrix", loop_over_const_matrix),
  TMC_CASE("matrix:Matrix.IterateOverMatrix", iterate_over_matrix),
  TMC_CASE("matrix:Matrix.Strides", strides),
  TMC_CASE("matrix:Matrix.ConversionToVector", conversion_to_vector),
  TMC_CASE("matrix:Matrix.ConstConversionToVector", const_conversion_to_vector),
  TMC_CASE("matrix:Matrix.ConversionFromVector", conversion_from_vector),
  TMC_CASE("matrix:Matrix.AssignmentFromVector", assignment_from_vector),
  TMC_CASE("matrix:Matrix.Axpy", axpy),
  TMC_CASE("matrix:Matrix.ForEach", for_each),
  TMC_CASE("matrix:Matrix.SetScaler", set_scalar),
  TMC_CASE("matrix:Matrix.Max", max),
  TMC_CASE("matrix:Matrix.Min", min),
  TMC_CASE("matrix:Matrix.ArrayFunction", array_function),
  TMC_CASE("matrix:Matrix.MultiMatrixArrayFunction", multi_matrix_array_function),
  TMC_CASE("matrix:Matrix.MismatchedRowDimensions", mismatched_row_dimensions),
  TMC_CASE("matrix:Matrix.MismatchedColumnDimensions", mismatched_column_dimensions),
  TMC_CASE("matrix:Matrix.WrongMatrixDimensions", wrong_matrix_dimensions),
  TMC_CASE("matrix:Matrix.MultipleTemporaries", multiple_temporaries),
  TMC_CASE("matrix:Matrix.ColumnViewReuse", column_view_reuse),
  TMC_CASE("matrix:Matrix.FunctionReusability", function_reusability),
  TMC_CASE("matrix:Matrix.ConstMatrixFunction", const_matrix_function),
  TMC_CASE("matrix:Matrix.EmptyMatrixFunction", empty_matrix_function),
  TMC_CASE("matrix:Matrix.MultiMatrixDifferentRowsFromCreation", multi_matrix_different_rows_from_creation),
  TMC_CASE("matrix:Matrix.MatrixVectorDifferentRowsFromCreation", matrix_vector_different_rows_from_creation),
  TMC_CASE("matrix:Matrix.MismatchedRowsAtInvocation", mismatched_rows_at_invocation),
  TMC_CASE("matrix:Matrix.MultipleMatricesMismatchedRowsAtInvocation", multiple_matrices_mismatched_rows_at_invocation),
  TMC_CASE("matrix:Matrix.WrongColumnCountAtInvocation", wrong_column_count_at_invocation),
  TMC_CASE("matrix:Matrix.Print", print),
  TMC_CASE("matrix:Matrix.VectorInMatrixFunction", vector_in_matrix_function),
  TMC_CASE("matrix:Matrix.VectorTooSmall", vector_too_small),
  TMC_CASE("matrix:Matrix.VectorTooLarge", vector_too_large),
  TMC_CASE("matrix:Matrix.EmptyVectorNonEmptyMatrix", empty_vector_non_empty_matrix),
  TMC_CASE("matrix:Matrix.NonEmptyVectorEmptyMatrix", non_empty_vector_empty_matrix),
  TMC_CASE("matrix:Matrix.EmptyVectorEmptyMatrix", empty_vector_empty_matrix),
  TMC_CASE("matrix:Matrix.MultipleVectorsDifferentSizes", multiple_vectors_different_sizes),
  TMC_CASE("matrix:Matrix.MultipleVectorsSameSize", multiple_vectors_same_size),
  TMC_CASE("matrix:Matrix.MultipleMatricesOneVector", multiple_matrices_one_vector),
  TMC_CASE("matrix:Matrix.MultipleMatricesDifferentRowsVector", multiple_matrices_different_rows_vector),
  TMC_CASE("matrix:Matrix.VectorSizeMatchesOneMatrixOnly", vector_length_matches_one_matrix_only),
  TMC_CASE("matrix:Matrix.ConstVector", const_vector),
  TMC_CASE("matrix:Matrix.MutableVector", mutable_vector),
  TMC_CASE("matrix:Matrix.FunctionReusabilityWithVectors", function_reusability_with_vectors),
  TMC_CASE("matrix:Matrix.FunctionInvocationWithWrongSizedVector", function_invocation_with_wrong_sized_vector),
  TMC_CASE("matrix:Matrix.MixedVectorColumnViewRowVariable", mixed_vector_column_view_row_variable),
  TMC_CASE("matrix:Matrix.IntegerVector", integer_vector),
  TMC_CASE("matrix:Matrix.FunctionWithConstSignature", function_with_const_signature),
  TMC_CASE("matrix:Matrix.TestFill", fill),
  TMC_CASE("matrix:Matrix.TestCopy", copy),
  TMC_CASE("matrix:Matrix.ReduceSum", reduce_sum),
  TMC_CASE("matrix:Matrix.ReduceMax", reduce_max),
  TMC_CASE("matrix:Matrix.ReduceLOr", reduce_lor),
  TMC_CASE("matrix:Matrix.ReduceLAnd", reduce_land),
  TMC_CASE("matrix:Matrix.ReduceStrict", reduce_strict),
  TMC_CASE("matrix:Matrix.VectorCapture", vector_capture),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_matrix(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
