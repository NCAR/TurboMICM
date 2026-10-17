/* MICM test/unit/util/test_sparse_matrix_vector_ordering_row.cpp: MICM's
 * SparseMatrix<T, SparseMatrixVectorOrderingCompressedSparseRow<L>> is
 * tmc_sparse with TMC_SPARSE_CSR and width L. A case MICM instantiates for
 * several widths runs the policy check (sparse_policy.c) at each of them; the
 * others run at the width MICM uses. The data-order checks MICM makes on the
 * returned matrices (elements numbered row by row, L blocks interleaved,
 * group after group) are here. */
#include "sparse_policy.h"

#include "tmc_sparse_ext.h"

#define ORDERING TMC_SPARSE_CSR

static int each_width(int (*check)(tmc_sparse_ordering, int))
{
  int rc = TMC_TEST_PASS;
  for (int L = 1; L <= 4; ++L)
    if (check(ORDERING, L) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
  return rc;
}

/* a case MICM runs for widths 1, 2, 3, and 4 */
#define WIDTHS_1_TO_4_CASE(name, check)                                                       \
  TMC_TEST(name)                                                                              \
  {                                                                                           \
    return each_width(check);                                                                 \
  }

/* a case MICM runs for one width */
#define WIDTH_CASE(name, check, width)                                                        \
  TMC_TEST(name)                                                                              \
  {                                                                                           \
    return check(ORDERING, width);                                                            \
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

/* MICM's GroupVectorSize, GroupSize, and NumberOfGroups(blocks) */
static void expect_groups(const tmc_sparse *m, int width, int group_size, int blocks, int groups)
{
  TMC_EXPECT_EQ(tmc_sparse_group_width(m), width);
  TMC_EXPECT_EQ(tmc_sparse_group_size(m), group_size);
  TMC_EXPECT_EQ(tmc_sparse_number_of_groups(m, blocks), groups);
}

TMC_TEST(zero_matrix)
{
  tmc_sparse *m = sparse_policy_zero_matrix(ORDERING, 2);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(const_zero_matrix)
{
  tmc_sparse *m = sparse_policy_const_zero_matrix(ORDERING, 3);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

WIDTH_CASE(set_scalar, sparse_policy_set_scalar, 3)

TMC_TEST(add_to_diagonal)
{
  int rc = TMC_TEST_PASS;
  const int widths[] = { 1, 3, 4 };
  for (int k = 0; k < 3; ++k)
    if (sparse_policy_add_to_diagonal(ORDERING, widths[k]) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
  return rc;
}

TMC_TEST(single_block_matrix)
{
  tmc_sparse *m = sparse_policy_single_block_matrix(ORDERING, 4);
  TMC_EXPECT_EQ(m->size, (size_t)(4 * 5));
  expect_one_block_index(m, 3, 2, 16, 42);
  expect_one_block_index(m, 2, 3, 12, 21);
  expect_groups(m, 4, 5 * 4, 1, 1);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

TMC_TEST(const_single_block_matrix)
{
  tmc_sparse *storage = sparse_policy_const_single_block_matrix(ORDERING, 2);
  const tmc_sparse *m = storage;
  TMC_EXPECT_EQ(m->size, (size_t)(2 * 5));
  expect_const_one_block_index(m, 3, 2, 8, 42);
  expect_const_one_block_index(m, 2, 3, 6, 21);
  expect_groups(m, 2, 2 * 5, 1, 1);
  tmc_sparse_free(storage);
  return TMC_TEST_PASS;
}

TMC_TEST(multi_block_matrix)
{
  tmc_sparse *m = sparse_policy_multi_block_matrix(ORDERING, 2);
  TMC_EXPECT_EQ(m->size, (size_t)(54 * 5)); /* 53 blocks padded to 27 groups of 2 */
  expect_block_index(m, 0, 2, 3, 6, 21);
  expect_block_index(m, 2, 2, 1, 14, 31);
  expect_groups(m, 2, 2 * 5, 4, 2);
  TMC_EXPECT_EQ(tmc_sparse_number_of_groups(m, 53), 27);
  tmc_sparse_free(m);
  return TMC_TEST_PASS;
}

WIDTHS_1_TO_4_CASE(print, sparse_policy_print)
WIDTHS_1_TO_4_CASE(print_non_zero, sparse_policy_print_non_zero)
WIDTHS_1_TO_4_CASE(array_function, sparse_policy_array_function)
WIDTHS_1_TO_4_CASE(multi_matrix_array_function, sparse_policy_multi_matrix_array_function)
WIDTHS_1_TO_4_CASE(multiple_temporaries, sparse_policy_multiple_temporaries)
WIDTHS_1_TO_4_CASE(block_view_reuse, sparse_policy_block_view_reuse)
WIDTHS_1_TO_4_CASE(function_reusability, sparse_policy_function_reusability)
WIDTHS_1_TO_4_CASE(two_sparse_matrices_different_structure, sparse_policy_two_sparse_matrices_different_structure)

TMC_TEST(sparse_and_dense_matrix_function)
{
  /* only L = 1 goes with MICM's standard-ordered Matrix */
  return sparse_policy_sparse_and_dense_matrix_function(ORDERING, 1, 1);
}

TMC_TEST(sparse_and_vector_matrix_function)
{
  /* a vector-ordered sparse matrix with a VectorMatrix of the same L */
  int rc = TMC_TEST_PASS;
  const int widths[] = { 1, 2, 3, 4, 10 };
  for (int k = 0; k < 5; ++k)
    if (sparse_policy_sparse_and_dense_matrix_function(ORDERING, widths[k], widths[k]) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
  return rc;
}

TMC_TEST(incompatible_ordering)
{
  /* L > 1 with a standard-ordered (L = 1) dense matrix */
  int rc = TMC_TEST_PASS;
  for (int L = 2; L <= 4; ++L)
    if (sparse_policy_incompatible_ordering(ORDERING, L) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
  return rc;
}

TMC_TEST(incompatible_vector_ordering)
{
  /* a VectorMatrix of another L */
  int rc = TMC_TEST_PASS;
  const int pairs[4][2] = { { 2, 1 }, { 2, 3 }, { 4, 2 }, { 4, 10 } };
  for (int k = 0; k < 4; ++k)
    if (sparse_policy_incompatible_vector_ordering(ORDERING, pairs[k][0], pairs[k][1]) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
  return rc;
}

TMC_TEST(incompatible_sparse_ordering)
{
  /* two sparse matrices of different L */
  int rc = TMC_TEST_PASS;
  const int pairs[4][2] = { { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 10 } };
  for (int k = 0; k < 4; ++k)
    if (sparse_policy_incompatible_sparse_ordering(ORDERING, pairs[k][0], pairs[k][1]) != TMC_TEST_PASS)
      rc = TMC_TEST_FAIL;
  return rc;
}

WIDTHS_1_TO_4_CASE(mismatched_block_dimensions, sparse_policy_mismatched_block_dimensions)
WIDTHS_1_TO_4_CASE(mismatched_element_dimensions, sparse_policy_mismatched_element_dimensions)
WIDTHS_1_TO_4_CASE(const_sparse_matrix_function, sparse_policy_const_sparse_matrix_function)
WIDTHS_1_TO_4_CASE(empty_sparse_matrix_function, sparse_policy_empty_sparse_matrix_function)
WIDTHS_1_TO_4_CASE(multiple_sparse_matrices_different_blocks_from_creation, sparse_policy_multiple_sparse_matrices_different_blocks_from_creation)
WIDTHS_1_TO_4_CASE(sparse_matrix_vector_different_blocks_from_creation, sparse_policy_sparse_matrix_vector_different_blocks_from_creation)
WIDTHS_1_TO_4_CASE(mismatched_blocks_at_invocation, sparse_policy_mismatched_blocks_at_invocation)
WIDTHS_1_TO_4_CASE(multiple_sparse_matrices_mismatched_blocks_at_invocation, sparse_policy_multiple_sparse_matrices_mismatched_blocks_at_invocation)

/* the vector cases, which MICM runs at L = 1 */
WIDTH_CASE(vector_in_sparse_matrix_function, sparse_policy_vector_in_sparse_matrix_function, 1)
WIDTH_CASE(vector_too_small, sparse_policy_vector_too_small, 1)
WIDTH_CASE(vector_too_large, sparse_policy_vector_too_large, 1)
WIDTH_CASE(empty_vector_non_empty_sparse_matrix, sparse_policy_empty_vector_non_empty_sparse_matrix, 1)
WIDTH_CASE(non_empty_vector_empty_sparse_matrix, sparse_policy_non_empty_vector_empty_sparse_matrix, 1)
WIDTH_CASE(empty_vector_empty_sparse_matrix, sparse_policy_empty_vector_empty_sparse_matrix, 1)
WIDTH_CASE(multiple_vectors_different_sizes, sparse_policy_multiple_vectors_different_sizes, 1)
WIDTH_CASE(multiple_vectors_same_size, sparse_policy_multiple_vectors_same_size, 1)
WIDTH_CASE(multiple_sparse_matrices_one_vector, sparse_policy_multiple_sparse_matrices_one_vector, 1)
WIDTH_CASE(multiple_sparse_matrices_different_blocks_vector, sparse_policy_multiple_sparse_matrices_different_blocks_vector, 1)
WIDTH_CASE(vector_size_matches_one_sparse_matrix_only, sparse_policy_vector_size_matches_one_sparse_matrix_only, 1)
WIDTH_CASE(const_vector_sparse, sparse_policy_const_vector_sparse, 1)
WIDTH_CASE(mutable_vector_sparse, sparse_policy_mutable_vector_sparse, 1)
WIDTH_CASE(function_reusability_with_vectors_sparse, sparse_policy_function_reusability_with_vectors_sparse, 1)
WIDTH_CASE(function_invocation_with_wrong_sized_vector_sparse, sparse_policy_function_invocation_with_wrong_sized_vector_sparse, 1)
WIDTH_CASE(mixed_vector_block_view_block_variable, sparse_policy_mixed_vector_block_view_block_variable, 1)
WIDTH_CASE(integer_vector_sparse, sparse_policy_integer_vector_sparse, 1)
WIDTHS_1_TO_4_CASE(function_with_const_signature_sparse, sparse_policy_function_with_const_signature_sparse)
WIDTHS_1_TO_4_CASE(get_block_view_by_vector_index, sparse_policy_get_block_view_by_vector_index)
WIDTHS_1_TO_4_CASE(test_fill, sparse_policy_fill)
WIDTHS_1_TO_4_CASE(test_copy, sparse_policy_copy)

#define ID(test) "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix." test

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
  TMC_CASE(ID("MultipleTemporaries"), multiple_temporaries),
  TMC_CASE(ID("BlockViewReuse"), block_view_reuse),
  TMC_CASE(ID("FunctionReusability"), function_reusability),
  TMC_CASE(ID("TwoSparseMatricesDifferentStructure"), two_sparse_matrices_different_structure),
  TMC_CASE(ID("SparseAndDenseMatrixFunction"), sparse_and_dense_matrix_function),
  TMC_CASE(ID("SparseAndVectorMatrixFunction"), sparse_and_vector_matrix_function),
  TMC_CASE(ID("IncompatibleOrdering"), incompatible_ordering),
  TMC_CASE(ID("IncompatibleVectorOrdering"), incompatible_vector_ordering),
  TMC_CASE(ID("IncompatibleSparseOrdering"), incompatible_sparse_ordering),
  TMC_CASE(ID("MismatchedBlockDimensions"), mismatched_block_dimensions),
  TMC_CASE(ID("MismatchedElementDimensions"), mismatched_element_dimensions),
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
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_sparse_matrix_vector_compressed_row_ordering(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
