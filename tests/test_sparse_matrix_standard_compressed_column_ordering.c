/* MICM test/unit/util/test_sparse_matrix_standard_ordering_column.cpp: MICM's
 * SparseMatrix<T, SparseMatrixStandardOrderingCompressedSparseColumn> is
 * tmc_sparse with TMC_SPARSE_CSC and L = 1. The shared policy checks are in
 * sparse_policy.c; the data-order checks MICM makes on the returned matrices
 * (elements numbered column by column, block after block) are here. */
#include "sparse_policy.h"

#include "tmc_sparse_ext.h"

#define ORDERING TMC_SPARSE_CSC

/* a case that runs one policy check on a standard-ordered matrix */
#define STANDARD_CASE(name, check)                                                            \
  TMC_TEST(name)                                                                              \
  {                                                                                           \
    return check(ORDERING, 1);                                                                \
  }

/* MICM: elem = matrix.VectorIndex(row, col); EXPECT_EQ(elem, expected);
 * matrix.AsVector()[elem] = value; EXPECT_EQ(matrix.AsVector()[expected], value) */
static void expect_one_block_index(tmc_sparse *m, int row, int col, size_t expected, double value)
{
  tmc_error err;
  tmc_error_clear(&err);
  size_t e = tmc_sparse_vector_index_one_block(m, row, col, &err);
  TMC_EXPECT_EQ(err.category, TMC_OK);
  TMC_EXPECT_EQ(e, expected);
  if (e < m->size)
  {
    m->data[e] = value;
    TMC_EXPECT_EQ(m->data[expected], value);
    TMC_EXPECT_EQ(tmc_sparse_get(m, 0, row, col), value);
  }
}

static void expect_const_one_block_index(const tmc_sparse *m, int row, int col, size_t expected, double value)
{
  tmc_error err;
  tmc_error_clear(&err);
  TMC_EXPECT_EQ(tmc_sparse_vector_index_one_block(m, row, col, &err), expected);
  TMC_EXPECT_EQ(err.category, TMC_OK);
  if (expected < m->size)
    TMC_EXPECT_EQ(m->data[expected], value);
}

static void expect_block_index(tmc_sparse *m, int block, int row, int col, size_t expected, double value)
{
  tmc_error err;
  tmc_error_clear(&err);
  size_t e = tmc_sparse_vector_index(m, block, row, col, &err);
  TMC_EXPECT_EQ(err.category, TMC_OK);
  TMC_EXPECT_EQ(e, expected);
  if (e < m->size)
  {
    m->data[e] = value;
    TMC_EXPECT_EQ(m->data[expected], value);
    TMC_EXPECT_EQ(tmc_sparse_get(m, block, row, col), value);
  }
}

TMC_TEST(zero_matrix)
{
  tmc_sparse *m = sparse_policy_zero_matrix(ORDERING, 1);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(const_zero_matrix)
{
  tmc_sparse *m = sparse_policy_const_zero_matrix(ORDERING, 1);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

STANDARD_CASE(set_scalar, sparse_policy_set_scalar)
STANDARD_CASE(add_to_diagonal, sparse_policy_add_to_diagonal)

TMC_TEST(single_block_matrix)
{
  tmc_sparse *m = sparse_policy_single_block_matrix(ORDERING, 1);
  TMC_EXPECT_EQ(m->L, 1);
  TMC_EXPECT_EQ(m->size, (size_t)5);
  expect_one_block_index(m, 3, 2, 3, 42);
  expect_one_block_index(m, 2, 3, 4, 21);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(const_single_block_matrix)
{
  tmc_sparse *storage = sparse_policy_const_single_block_matrix(ORDERING, 1);
  const tmc_sparse *m = storage;
  expect_const_one_block_index(m, 3, 2, 3, 42);
  expect_const_one_block_index(m, 2, 3, 4, 21);
  tmc_sparse_free(storage);
  return TMC_TEST_PASS;
}

TMC_TEST(multi_block_matrix)
{
  tmc_sparse *m = sparse_policy_multi_block_matrix(ORDERING, 1);
  TMC_EXPECT_EQ(m->size, (size_t)(53 * 5));
  expect_block_index(m, 0, 2, 3, 4, 21);
  expect_block_index(m, 2, 2, 1, 12, 31);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

STANDARD_CASE(print, sparse_policy_print)
STANDARD_CASE(print_non_zero, sparse_policy_print_non_zero)

TMC_TEST(bad_configuration)
{
  return sparse_policy_builder_bad_configuration();
}

STANDARD_CASE(array_function, sparse_policy_array_function)
STANDARD_CASE(multi_matrix_array_function, sparse_policy_multi_matrix_array_function)
STANDARD_CASE(mismatched_block_dimensions, sparse_policy_mismatched_block_dimensions)
STANDARD_CASE(mismatched_element_dimensions, sparse_policy_mismatched_element_dimensions)
STANDARD_CASE(multiple_temporaries, sparse_policy_multiple_temporaries)
STANDARD_CASE(block_view_reuse, sparse_policy_block_view_reuse)
STANDARD_CASE(function_reusability, sparse_policy_function_reusability)
STANDARD_CASE(two_sparse_matrices_different_structure, sparse_policy_two_sparse_matrices_different_structure)

TMC_TEST(sparse_and_dense_matrix_function)
{
  /* with MICM's Matrix, L = 1 */
  return sparse_policy_sparse_and_dense_matrix_function(ORDERING, 1, 1);
}

STANDARD_CASE(const_sparse_matrix_function, sparse_policy_const_sparse_matrix_function)
STANDARD_CASE(empty_sparse_matrix_function, sparse_policy_empty_sparse_matrix_function)
STANDARD_CASE(multiple_sparse_matrices_different_blocks_from_creation, sparse_policy_multiple_sparse_matrices_different_blocks_from_creation)
STANDARD_CASE(sparse_matrix_vector_different_blocks_from_creation, sparse_policy_sparse_matrix_vector_different_blocks_from_creation)
STANDARD_CASE(mismatched_blocks_at_invocation, sparse_policy_mismatched_blocks_at_invocation)
STANDARD_CASE(multiple_sparse_matrices_mismatched_blocks_at_invocation, sparse_policy_multiple_sparse_matrices_mismatched_blocks_at_invocation)
STANDARD_CASE(vector_in_sparse_matrix_function, sparse_policy_vector_in_sparse_matrix_function)
STANDARD_CASE(vector_too_small, sparse_policy_vector_too_small)
STANDARD_CASE(vector_too_large, sparse_policy_vector_too_large)
STANDARD_CASE(empty_vector_non_empty_sparse_matrix, sparse_policy_empty_vector_non_empty_sparse_matrix)
STANDARD_CASE(non_empty_vector_empty_sparse_matrix, sparse_policy_non_empty_vector_empty_sparse_matrix)
STANDARD_CASE(empty_vector_empty_sparse_matrix, sparse_policy_empty_vector_empty_sparse_matrix)
STANDARD_CASE(multiple_vectors_different_sizes, sparse_policy_multiple_vectors_different_sizes)
STANDARD_CASE(multiple_vectors_same_size, sparse_policy_multiple_vectors_same_size)
STANDARD_CASE(multiple_sparse_matrices_one_vector, sparse_policy_multiple_sparse_matrices_one_vector)
STANDARD_CASE(multiple_sparse_matrices_different_blocks_vector, sparse_policy_multiple_sparse_matrices_different_blocks_vector)
STANDARD_CASE(vector_size_matches_one_sparse_matrix_only, sparse_policy_vector_size_matches_one_sparse_matrix_only)
STANDARD_CASE(const_vector_sparse, sparse_policy_const_vector_sparse)
STANDARD_CASE(mutable_vector_sparse, sparse_policy_mutable_vector_sparse)
STANDARD_CASE(function_reusability_with_vectors_sparse, sparse_policy_function_reusability_with_vectors_sparse)
STANDARD_CASE(function_invocation_with_wrong_sized_vector_sparse, sparse_policy_function_invocation_with_wrong_sized_vector_sparse)
STANDARD_CASE(mixed_vector_block_view_block_variable, sparse_policy_mixed_vector_block_view_block_variable)
STANDARD_CASE(integer_vector_sparse, sparse_policy_integer_vector_sparse)
STANDARD_CASE(function_with_const_signature_sparse, sparse_policy_function_with_const_signature_sparse)
STANDARD_CASE(get_block_view_by_vector_index, sparse_policy_get_block_view_by_vector_index)
STANDARD_CASE(test_fill, sparse_policy_fill)
STANDARD_CASE(test_copy, sparse_policy_copy)

#define ID(test) "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix." test

static const tmc_test_case cases[] = {
  TMC_CASE(ID("ZeroMatrix"), zero_matrix),
  TMC_CASE(ID("ConstZeroMatrix"), const_zero_matrix),
  TMC_CASE(ID("SetScalar"), set_scalar),
  TMC_CASE(ID("AddToDiagonal"), add_to_diagonal),
  TMC_CASE(ID("SingleBlockMatrix"), single_block_matrix),
  TMC_CASE(ID("ConstSingleBlockMatrix"), const_single_block_matrix),
  TMC_CASE(ID("MultiBlockMatrix"), multi_block_matrix),
  TMC_CASE(ID("Print"), print),
  TMC_CASE(ID("PrintNonZero"), print_non_zero),
  TMC_CASE(ID("ArrayFunction"), array_function),
  TMC_CASE(ID("MultiMatrixArrayFunction"), multi_matrix_array_function),
  TMC_CASE(ID("MismatchedBlockDimensions"), mismatched_block_dimensions),
  TMC_CASE(ID("MismatchedElementDimensions"), mismatched_element_dimensions),
  TMC_CASE(ID("MultipleTemporaries"), multiple_temporaries),
  TMC_CASE(ID("BlockViewReuse"), block_view_reuse),
  TMC_CASE(ID("FunctionReusability"), function_reusability),
  TMC_CASE(ID("TwoSparseMatricesDifferentStructure"), two_sparse_matrices_different_structure),
  TMC_CASE(ID("SparseAndDenseMatrixFunction"), sparse_and_dense_matrix_function),
  TMC_CASE(ID("ConstSparseMatrixFunction"), const_sparse_matrix_function),
  TMC_CASE(ID("EmptySparseMatrixFunction"), empty_sparse_matrix_function),
  TMC_CASE(ID("MultipleSparseMatricesDifferentBlocksFromCreation"), multiple_sparse_matrices_different_blocks_from_creation),
  TMC_CASE(ID("SparseMatrixVectorDifferentBlocksFromCreation"), sparse_matrix_vector_different_blocks_from_creation),
  TMC_CASE(ID("MismatchedBlocksAtInvocation"), mismatched_blocks_at_invocation),
  TMC_CASE(ID("MultipleSparseMatricesMismatchedBlocksAtInvocation"), multiple_sparse_matrices_mismatched_blocks_at_invocation),
  TMC_CASE(ID("VectorInSparseMatrixFunction"), vector_in_sparse_matrix_function),
  TMC_CASE(ID("VectorTooSmall"), vector_too_small),
  TMC_CASE(ID("VectorTooLarge"), vector_too_large),
  TMC_CASE(ID("EmptyVectorNonEmptySparseMatrix"), empty_vector_non_empty_sparse_matrix),
  TMC_CASE(ID("NonEmptyVectorEmptySparseMatrix"), non_empty_vector_empty_sparse_matrix),
  TMC_CASE(ID("EmptyVectorEmptySparseMatrix"), empty_vector_empty_sparse_matrix),
  TMC_CASE(ID("MultipleVectorsDifferentSizes"), multiple_vectors_different_sizes),
  TMC_CASE(ID("MultipleVectorsSameSize"), multiple_vectors_same_size),
  TMC_CASE(ID("MultipleSparseMatricesOneVector"), multiple_sparse_matrices_one_vector),
  TMC_CASE(ID("MultipleSparseMatricesDifferentBlocksVector"), multiple_sparse_matrices_different_blocks_vector),
  TMC_CASE(ID("VectorSizeMatchesOneSparseMatrixOnly"), vector_size_matches_one_sparse_matrix_only),
  TMC_CASE(ID("ConstVectorSparse"), const_vector_sparse),
  TMC_CASE(ID("MutableVectorSparse"), mutable_vector_sparse),
  TMC_CASE(ID("FunctionReusabilityWithVectorsSparse"), function_reusability_with_vectors_sparse),
  TMC_CASE(ID("FunctionInvocationWithWrongSizedVectorSparse"), function_invocation_with_wrong_sized_vector_sparse),
  TMC_CASE(ID("MixedVectorBlockViewBlockVariable"), mixed_vector_block_view_block_variable),
  TMC_CASE(ID("IntegerVectorSparse"), integer_vector_sparse),
  TMC_CASE(ID("FunctionWithConstSignatureSparse"), function_with_const_signature_sparse),
  TMC_CASE(ID("GetBlockViewByVectorIndex"), get_block_view_by_vector_index),
  TMC_CASE(ID("TestFill"), test_fill),
  TMC_CASE(ID("TestCopy"), test_copy),
  TMC_CASE("sparse_matrix_standard_compressed_column_ordering:SparseMatrixBuilder.BadConfiguration", bad_configuration),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_sparse_matrix_standard_compressed_column_ordering(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
