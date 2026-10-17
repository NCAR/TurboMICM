/* MICM test/unit/util/test_vector_matrix.cpp: MICM's VectorMatrix<L> is
 * tmc_dense with width L. A case MICM instantiates for VectorMatrix<1> to <4>
 * runs the policy check (dense_policy.c) at each of those widths; the others
 * run at the width MICM uses. The layout checks MICM makes on the returned
 * matrices are here. */
#include "dense_policy.h"

#include "tmc_dense_ext.h"

static int each_width(int (*check)(int L))
{
  int rc = TMC_TEST_PASS;
  for (int L = 1; L <= 4; ++L)
    if (check(L) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
  return rc;
}

/* a case MICM runs for VectorMatrix<1>, <2>, <3>, and <4> */
#define WIDTHS_1_TO_4_CASE(name, check)                                                       \
  TMC_TEST(name)                                                                              \
  {                                                                                           \
    return each_width(check);                                                                 \
  }

/* a case MICM runs for one width */
#define WIDTH_CASE(name, check, width)                                                        \
  TMC_TEST(name)                                                                              \
  {                                                                                           \
    return check(width);                                                                      \
  }

TMC_TEST(small_vector_matrix)
{
  tmc_dense *m = dense_policy_small_matrix(2);
  TMC_EXPECT_EQ(m->size, (size_t)(4 * 5));
  TMC_EXPECT_EQ(tmc_dense_group_size(m), 2 * 5);
  TMC_EXPECT_EQ(tmc_dense_number_of_groups(m), 2);
  TMC_EXPECT_EQ(tmc_dense_group_width(m), 2);
  TMC_EXPECT_EQ(m->data[0], 41.2);
  TMC_EXPECT_EQ(m->data[2 * 5 + 0 + 2 * 4], 102.3);
  TMC_EXPECT_EQ(m->data[1 + 2 * 3], 64.7);
  tmc_dense_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(small_const_vector_matrix)
{
  tmc_dense *storage = dense_policy_small_const_matrix(4);
  const tmc_dense *m = storage;
  TMC_EXPECT_EQ(m->size, (size_t)(4 * 5));
  TMC_EXPECT_EQ(tmc_dense_group_size(m), 4 * 5);
  TMC_EXPECT_EQ(tmc_dense_number_of_groups(m), 1);
  TMC_EXPECT_EQ(tmc_dense_group_width(m), 4);
  TMC_EXPECT_EQ(m->data[0], 41.2);
  TMC_EXPECT_EQ(m->data[2 + 4 * 4], 102.3);
  TMC_EXPECT_EQ(m->data[1 + 4 * 3], 64.7);
  tmc_dense_free(storage);
  return TMC_TEST_PASS;
}

WIDTH_CASE(initialize_vector_matrix, dense_policy_initialize_matrix, 1)
WIDTH_CASE(initialize_const_vector_matrix, dense_policy_initialize_const_matrix, 2)
WIDTH_CASE(loop_over_vector_matrix, dense_policy_loop_over_matrix, 2)
WIDTH_CASE(loop_over_const_vector_matrix, dense_policy_loop_over_const_matrix, 1)

TMC_TEST(strides)
{
  tmc_dense *matrix3vec = dense_policy_strides(3);
  TMC_EXPECT_EQ(tmc_dense_row_stride(matrix3vec), (size_t)1);
  TMC_EXPECT_EQ(tmc_dense_column_stride(matrix3vec), (size_t)3);
  tmc_dense *matrix4vec = dense_policy_strides(4);
  TMC_EXPECT_EQ(tmc_dense_row_stride(matrix4vec), (size_t)1);
  TMC_EXPECT_EQ(tmc_dense_column_stride(matrix4vec), (size_t)4);
  tmc_dense_free(matrix3vec);
  tmc_dense_free(matrix4vec);
  return TMC_TEST_PASS;
}

WIDTH_CASE(conversion_to_vector, dense_policy_conversion_to_vector, 3)
WIDTH_CASE(const_conversion_to_vector, dense_policy_const_conversion_to_vector, 1)
WIDTH_CASE(conversion_from_vector, dense_policy_conversion_from_vector, 2)
WIDTH_CASE(assignment_from_vector, dense_policy_assignment_from_vector, 2)
WIDTHS_1_TO_4_CASE(axpy, dense_policy_axpy)
WIDTHS_1_TO_4_CASE(for_each, dense_policy_for_each)
WIDTHS_1_TO_4_CASE(set_scalar, dense_policy_set_scalar)
WIDTHS_1_TO_4_CASE(max, dense_policy_max)
WIDTHS_1_TO_4_CASE(min, dense_policy_min)
WIDTHS_1_TO_4_CASE(print, dense_policy_print)
WIDTHS_1_TO_4_CASE(array_function, dense_policy_array_function)
WIDTHS_1_TO_4_CASE(multi_matrix_array_function, dense_policy_multi_matrix_array_function)
WIDTHS_1_TO_4_CASE(mismatched_row_dimensions, dense_policy_mismatched_row_dimensions)
WIDTHS_1_TO_4_CASE(mismatched_column_dimensions, dense_policy_mismatched_column_dimensions)
WIDTHS_1_TO_4_CASE(wrong_matrix_dimensions, dense_policy_wrong_matrix_dimensions)
WIDTHS_1_TO_4_CASE(multiple_temporaries, dense_policy_multiple_temporaries)
WIDTHS_1_TO_4_CASE(column_view_reuse, dense_policy_column_view_reuse)
WIDTHS_1_TO_4_CASE(function_reusability, dense_policy_function_reusability)
WIDTHS_1_TO_4_CASE(const_matrix_function, dense_policy_const_matrix_function)
WIDTHS_1_TO_4_CASE(empty_matrix_function, dense_policy_empty_matrix_function)
WIDTHS_1_TO_4_CASE(multi_matrix_different_rows_from_creation, dense_policy_multi_matrix_different_rows_from_creation)
WIDTHS_1_TO_4_CASE(matrix_vector_different_rows_from_creation, dense_policy_matrix_vector_different_rows_from_creation)
WIDTHS_1_TO_4_CASE(mismatched_rows_at_invocation, dense_policy_mismatched_rows_at_invocation)
WIDTHS_1_TO_4_CASE(multiple_matrices_mismatched_rows_at_invocation, dense_policy_multiple_matrices_mismatched_rows_at_invocation)
WIDTHS_1_TO_4_CASE(wrong_column_count_at_invocation, dense_policy_wrong_column_count_at_invocation)
WIDTHS_1_TO_4_CASE(vector_in_matrix_function, dense_policy_vector_in_matrix_function)
WIDTHS_1_TO_4_CASE(vector_too_small, dense_policy_vector_too_small)
WIDTHS_1_TO_4_CASE(vector_too_large, dense_policy_vector_too_large)
WIDTHS_1_TO_4_CASE(empty_vector_non_empty_matrix, dense_policy_empty_vector_non_empty_matrix)
WIDTHS_1_TO_4_CASE(non_empty_vector_empty_matrix, dense_policy_non_empty_vector_empty_matrix)
WIDTHS_1_TO_4_CASE(empty_vector_empty_matrix, dense_policy_empty_vector_empty_matrix)
WIDTHS_1_TO_4_CASE(multiple_vectors_different_sizes, dense_policy_multiple_vectors_different_sizes)
WIDTHS_1_TO_4_CASE(multiple_vectors_same_size, dense_policy_multiple_vectors_same_size)
WIDTHS_1_TO_4_CASE(multiple_matrices_one_vector, dense_policy_multiple_matrices_one_vector)
WIDTHS_1_TO_4_CASE(multiple_matrices_different_rows_vector, dense_policy_multiple_matrices_different_rows_vector)
WIDTHS_1_TO_4_CASE(vector_length_matches_one_matrix_only, dense_policy_vector_length_matches_one_matrix_only)
WIDTHS_1_TO_4_CASE(const_vector, dense_policy_const_vector)
WIDTHS_1_TO_4_CASE(mutable_vector, dense_policy_mutable_vector)
WIDTHS_1_TO_4_CASE(function_reusability_with_vectors, dense_policy_function_reusability_with_vectors)
WIDTHS_1_TO_4_CASE(function_invocation_with_wrong_sized_vector, dense_policy_function_invocation_with_wrong_sized_vector)
WIDTHS_1_TO_4_CASE(mixed_vector_column_view_row_variable, dense_policy_mixed_vector_column_view_row_variable)
WIDTHS_1_TO_4_CASE(integer_vector, dense_policy_integer_vector)
WIDTHS_1_TO_4_CASE(function_with_const_signature, dense_policy_function_with_const_signature)
WIDTHS_1_TO_4_CASE(fill, dense_policy_fill)
WIDTHS_1_TO_4_CASE(copy, dense_policy_copy)
WIDTHS_1_TO_4_CASE(reduce_sum, dense_policy_reduce_sum)
WIDTHS_1_TO_4_CASE(reduce_max, dense_policy_reduce_max)
WIDTHS_1_TO_4_CASE(reduce_lor, dense_policy_reduce_lor)
WIDTHS_1_TO_4_CASE(reduce_land, dense_policy_reduce_land)
WIDTHS_1_TO_4_CASE(reduce_strict, dense_policy_reduce_strict)
WIDTHS_1_TO_4_CASE(vector_capture, dense_policy_vector_capture)

static const tmc_test_case cases[] = {
  TMC_CASE("vector_matrix:VectorMatrix.SmallVectorMatrix", small_vector_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.SmallConstVectorMatrix", small_const_vector_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.InitializeVectorMatrix", initialize_vector_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.InitializeConstVectorMatrix", initialize_const_vector_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.LoopOverVectorMatrix", loop_over_vector_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.LoopOverConstVectorMatrix", loop_over_const_vector_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.Strides", strides),
  TMC_CASE("vector_matrix:VectorMatrix.ConversionToVector", conversion_to_vector),
  TMC_CASE("vector_matrix:VectorMatrix.ConstConversionToVector", const_conversion_to_vector),
  TMC_CASE("vector_matrix:VectorMatrix.ConversionFromVector", conversion_from_vector),
  TMC_CASE("vector_matrix:VectorMatrix.AssignmentFromVector", assignment_from_vector),
  TMC_CASE("vector_matrix:VectorMatrix.Axpy", axpy),
  TMC_CASE("vector_matrix:VectorMatrix.ForEach", for_each),
  TMC_CASE("vector_matrix:VectorMatrix.SetScaler", set_scalar),
  TMC_CASE("vector_matrix:VectorMatrix.Max", max),
  TMC_CASE("vector_matrix:VectorMatrix.Min", min),
  TMC_CASE("vector_matrix:VectorMatrix.Print", print),
  TMC_CASE("vector_matrix:VectorMatrix.ArrayFunction", array_function),
  TMC_CASE("vector_matrix:VectorMatrix.MultiMatrixArrayFunction", multi_matrix_array_function),
  TMC_CASE("vector_matrix:VectorMatrix.MismatchedRowDimensions", mismatched_row_dimensions),
  TMC_CASE("vector_matrix:VectorMatrix.MismatchedColumnDimensions", mismatched_column_dimensions),
  TMC_CASE("vector_matrix:VectorMatrix.WrongMatrixDimensions", wrong_matrix_dimensions),
  TMC_CASE("vector_matrix:VectorMatrix.MultipleTemporaries", multiple_temporaries),
  TMC_CASE("vector_matrix:VectorMatrix.ColumnViewReuse", column_view_reuse),
  TMC_CASE("vector_matrix:VectorMatrix.FunctionReusability", function_reusability),
  TMC_CASE("vector_matrix:VectorMatrix.ConstMatrixFunction", const_matrix_function),
  TMC_CASE("vector_matrix:VectorMatrix.EmptyMatrixFunction", empty_matrix_function),
  TMC_CASE("vector_matrix:VectorMatrix.MultiMatrixDifferentRowsFromCreation", multi_matrix_different_rows_from_creation),
  TMC_CASE("vector_matrix:VectorMatrix.MatrixVectorDifferentRowsFromCreation", matrix_vector_different_rows_from_creation),
  TMC_CASE("vector_matrix:VectorMatrix.MismatchedRowsAtInvocation", mismatched_rows_at_invocation),
  TMC_CASE("vector_matrix:VectorMatrix.MultipleMatricesMismatchedRowsAtInvocation", multiple_matrices_mismatched_rows_at_invocation),
  TMC_CASE("vector_matrix:VectorMatrix.WrongColumnCountAtInvocation", wrong_column_count_at_invocation),
  TMC_CASE("vector_matrix:VectorMatrix.VectorInMatrixFunction", vector_in_matrix_function),
  TMC_CASE("vector_matrix:VectorMatrix.VectorTooSmall", vector_too_small),
  TMC_CASE("vector_matrix:VectorMatrix.VectorTooLarge", vector_too_large),
  TMC_CASE("vector_matrix:VectorMatrix.EmptyVectorNonEmptyMatrix", empty_vector_non_empty_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.NonEmptyVectorEmptyMatrix", non_empty_vector_empty_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.EmptyVectorEmptyMatrix", empty_vector_empty_matrix),
  TMC_CASE("vector_matrix:VectorMatrix.MultipleVectorsDifferentSizes", multiple_vectors_different_sizes),
  TMC_CASE("vector_matrix:VectorMatrix.MultipleVectorsSameSize", multiple_vectors_same_size),
  TMC_CASE("vector_matrix:VectorMatrix.MultipleMatricesOneVector", multiple_matrices_one_vector),
  TMC_CASE("vector_matrix:VectorMatrix.MultipleMatricesDifferentRowsVector", multiple_matrices_different_rows_vector),
  TMC_CASE("vector_matrix:VectorMatrix.VectorSizeMatchesOneMatrixOnly", vector_length_matches_one_matrix_only),
  TMC_CASE("vector_matrix:VectorMatrix.ConstVector", const_vector),
  TMC_CASE("vector_matrix:VectorMatrix.MutableVector", mutable_vector),
  TMC_CASE("vector_matrix:VectorMatrix.FunctionReusabilityWithVectors", function_reusability_with_vectors),
  TMC_CASE("vector_matrix:VectorMatrix.FunctionInvocationWithWrongSizedVector", function_invocation_with_wrong_sized_vector),
  TMC_CASE("vector_matrix:VectorMatrix.MixedVectorColumnViewRowVariable", mixed_vector_column_view_row_variable),
  TMC_CASE("vector_matrix:VectorMatrix.IntegerVector", integer_vector),
  TMC_CASE("vector_matrix:VectorMatrix.FunctionWithConstSignature", function_with_const_signature),
  TMC_CASE("vector_matrix:VectorMatrix.TestFill", fill),
  TMC_CASE("vector_matrix:VectorMatrix.TestCopy", copy),
  TMC_CASE("vector_matrix:VectorMatrix.ReduceSum", reduce_sum),
  TMC_CASE("vector_matrix:VectorMatrix.ReduceMax", reduce_max),
  TMC_CASE("vector_matrix:VectorMatrix.ReduceLOr", reduce_lor),
  TMC_CASE("vector_matrix:VectorMatrix.ReduceLAnd", reduce_land),
  TMC_CASE("vector_matrix:VectorMatrix.ReduceStrict", reduce_strict),
  TMC_CASE("vector_matrix:VectorMatrix.VectorCapture", vector_capture),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_vector_matrix(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
