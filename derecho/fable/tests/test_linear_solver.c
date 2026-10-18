/* MICM test/unit/solver/test_linear_solver.cpp with
 * test_linear_solver_policy.hpp: MICM's LinearSolver factors A into separate
 * L and U matrices with its default (Doolittle) decomposition and solves
 * A x = b; here tmc_lu_decompose and tmc_lu_solve. Each check compares A x
 * with b. MICM's standard ordering is L = 1, its vector orderings L = 1 to 4.
 * The random matrices are MICM's own (libcxx_random.h). */
#include "harness.h"
#include "libcxx_random.h"

/* MICM's CheckResults: every row of A x equals b within tol. */
static void check_results(const tmc_sparse *A, const tmc_dense *b, const tmc_dense *x, double tol)
{
  TMC_EXPECT_EQ(A->n_blocks, b->n_rows);
  TMC_EXPECT_EQ(A->n_blocks, x->n_rows);
  for (int blk = 0; blk < A->n_blocks; ++blk)
    for (int i = 0; i < A->n; ++i)
    {
      double result = 0.0;
      for (int j = 0; j < A->n; ++j)
        if (!tmc_sparse_is_zero(A, i, j))
          result += tmc_sparse_get(A, blk, i, j) * tmc_dense_get(x, blk, j);
      TMC_EXPECT_NEAR(tmc_dense_get(b, blk, i), result, tol);
    }
}

/* Factor A into L and U (initialized to initial_value), solve for x = b, and check. */
static void factor_solve_and_check(const tmc_sparse *A, const tmc_dense *b, double initial_value, double tol)
{
  tmc_error err;
  tmc_error_clear(&err);
  tmc_lu *lu = tmc_lu_create(A, TMC_LU_DOOLITTLE, &err);
  TMC_EXPECT_TRUE(lu != NULL);
  if (!lu)
    return;
  tmc_sparse *Lm = tmc_lu_lower_matrix(lu, A->n_blocks, initial_value);
  tmc_sparse *U = tmc_lu_upper_matrix(lu, A->n_blocks, initial_value);
  tmc_dense *x = tmc_dense_copy(b);
  tmc_lu_decompose(lu, A, Lm, U);
  tmc_lu_solve(lu, Lm, U, x);
  check_results(A, b, x, tol);
  tmc_dense_free(x);
  tmc_sparse_free(Lm);
  tmc_sparse_free(U);
  tmc_lu_free(lu);
}

static tmc_sparse *build_matrix(tmc_sparse_builder *builder, int L)
{
  tmc_sparse *A = tmc_sparse_build(builder, TMC_SPARSE_CSR, L, false);
  tmc_sparse_builder_free(builder);
  return A;
}

/* A 3 x 3 dense system, one block. */
static void test_dense_matrix(int L)
{
  tmc_sparse_builder *builder = tmc_sparse_builder_create(3);
  tmc_sparse_builder_set_initial_value(builder, 0.0);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      tmc_sparse_builder_with_element(builder, i, j, NULL);
  tmc_sparse *A = build_matrix(builder, L);
  const double values[3][3] = { { 2, -1, -2 }, { -4, 6, 3 }, { -4, -2, 8 } };
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      tmc_sparse_set(A, 0, i, j, values[i][j]);
  tmc_dense *b = tmc_dense_create(1, 3, L, 0.0);
  tmc_dense_set(b, 0, 0, 23);
  tmc_dense_set(b, 0, 1, 42);
  tmc_dense_set(b, 0, 2, 9);
  factor_solve_and_check(A, b, 0.0, 1.0e-5);
  tmc_dense_free(b);
  tmc_sparse_free(A);
}

/* A random 10 x 10 pattern (diagonal always present) of lognormal values,
 * different in every block, and a lognormal right-hand side. */
static void test_random_matrix(int L, int n_blocks)
{
  tmc_bool_stream gen_bool = tmc_bool_stream_create();
  tmc_lognormal_stream get_double = tmc_lognormal_stream_create(-2.0, 2.0);
  tmc_sparse_builder *builder = tmc_sparse_builder_create(10);
  tmc_sparse_builder_set_number_of_blocks(builder, n_blocks);
  tmc_sparse_builder_set_initial_value(builder, 0.0);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      if (i == j || tmc_bool_stream_next(&gen_bool))
        tmc_sparse_builder_with_element(builder, i, j, NULL);
  tmc_sparse *A = build_matrix(builder, L);
  tmc_dense *b = tmc_dense_create(n_blocks, 10, L, 0.0);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      if (!tmc_sparse_is_zero(A, i, j))
        for (int blk = 0; blk < n_blocks; ++blk)
          tmc_sparse_set(A, blk, i, j, tmc_lognormal_stream_next(&get_double));
  for (int i = 0; i < 10; ++i)
    for (int blk = 0; blk < n_blocks; ++blk)
      tmc_dense_set(b, blk, i, tmc_lognormal_stream_next(&get_double));
  factor_solve_and_check(A, b, 0.0, 1.0e-6);
  tmc_dense_free(b);
  tmc_sparse_free(A);
}

/* A random 30 x 30 system from mt19937(12345) with L and U initialized to an
 * extreme value. */
static void test_extreme_initial_value(int L, int n_blocks, double initial_value)
{
  tmc_mt19937 bool_engine = tmc_mt19937_create(12345u);
  tmc_mt19937 double_engine = tmc_mt19937_create(12345u);
  tmc_lognormal get_double = tmc_lognormal_create(-2.0, 2.0);
  const int size = 30;
  tmc_sparse_builder *builder = tmc_sparse_builder_create(size);
  tmc_sparse_builder_set_number_of_blocks(builder, n_blocks);
  tmc_sparse_builder_set_initial_value(builder, 0.0);
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
      if (i == j || tmc_mt19937_bernoulli(&bool_engine, 0.5))
        tmc_sparse_builder_with_element(builder, i, j, NULL);
  tmc_sparse *A = build_matrix(builder, L);
  tmc_dense *b = tmc_dense_create(n_blocks, size, L, 0.0);
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
      if (!tmc_sparse_is_zero(A, i, j))
        for (int blk = 0; blk < n_blocks; ++blk)
          tmc_sparse_set(A, blk, i, j, tmc_mt19937_lognormal(&get_double, &double_engine));
  for (int i = 0; i < size; ++i)
    for (int blk = 0; blk < n_blocks; ++blk)
      tmc_dense_set(b, blk, i, tmc_mt19937_lognormal(&get_double, &double_engine));
  factor_solve_and_check(A, b, initial_value, 2.0e-6);
  tmc_dense_free(b);
  tmc_sparse_free(A);
}

/* A 6 x 6 diagonal system of lognormal(-2, 4) values. */
static void test_diagonal_matrix(int L, int n_blocks)
{
  tmc_lognormal_stream get_double = tmc_lognormal_stream_create(-2.0, 4.0);
  tmc_sparse_builder *builder = tmc_sparse_builder_create(6);
  tmc_sparse_builder_set_number_of_blocks(builder, n_blocks);
  tmc_sparse_builder_set_initial_value(builder, 0.0);
  for (int i = 0; i < 6; ++i)
    tmc_sparse_builder_with_element(builder, i, i, NULL);
  tmc_sparse *A = build_matrix(builder, L);
  tmc_dense *b = tmc_dense_create(n_blocks, 6, L, 0.0);
  for (int i = 0; i < 6; ++i)
    for (int blk = 0; blk < n_blocks; ++blk)
      tmc_sparse_set(A, blk, i, i, tmc_lognormal_stream_next(&get_double));
  for (int i = 0; i < 6; ++i)
    for (int blk = 0; blk < n_blocks; ++blk)
      tmc_dense_set(b, blk, i, tmc_lognormal_stream_next(&get_double));
  factor_solve_and_check(A, b, 0.0, 1.0e-5);
  tmc_dense_free(b);
  tmc_sparse_free(A);
}

/* MICM's TestMarkowitzReordering: the diagonal Markowitz order of a random
 * 50 x 50 pattern is a permutation, and factoring the reordered pattern
 * needs fewer L and U elements than the original. */
static void test_markowitz_reordering(int L)
{
  enum { order = 50 };
  tmc_bool_stream gen_bool = tmc_bool_stream_create();
  unsigned char orig[order * order];
  for (int i = 0; i < order; ++i)
    for (int j = 0; j < order; ++j)
      orig[i * order + j] = (i == j || tmc_bool_stream_next(&gen_bool)) ? 1 : 0;
  int reorder_map[order];
  tmc_markowitz_reorder(order, orig, reorder_map);

  tmc_sparse_builder *orig_builder = tmc_sparse_builder_create(order);
  tmc_sparse_builder *reordered_builder = tmc_sparse_builder_create(order);
  for (int i = 0; i < order; ++i)
    for (int j = 0; j < order; ++j)
    {
      if (orig[i * order + j] != 0)
        tmc_sparse_builder_with_element(orig_builder, i, j, NULL);
      if (orig[reorder_map[i] * order + reorder_map[j]] != 0)
        tmc_sparse_builder_with_element(reordered_builder, i, j, NULL);
    }
  tmc_sparse *orig_jac = build_matrix(orig_builder, L);
  tmc_sparse *reordered_jac = build_matrix(reordered_builder, L);

  tmc_error err;
  tmc_error_clear(&err);
  tmc_lu *orig_lu = tmc_lu_create(orig_jac, TMC_LU_DOOLITTLE, &err);
  tmc_lu *reordered_lu = tmc_lu_create(reordered_jac, TMC_LU_DOOLITTLE, &err);
  TMC_EXPECT_TRUE(orig_lu != NULL && reordered_lu != NULL);
  if (orig_lu && reordered_lu)
  {
    tmc_sparse *orig_L = tmc_lu_lower_matrix(orig_lu, 1, 0.0), *orig_U = tmc_lu_upper_matrix(orig_lu, 1, 0.0);
    tmc_sparse *reordered_L = tmc_lu_lower_matrix(reordered_lu, 1, 0.0), *reordered_U = tmc_lu_upper_matrix(reordered_lu, 1, 0.0);
    long sum_orig = 0, sum_reordered = 0;
    for (int i = 0; i < order; ++i)
    {
      sum_orig += i;
      sum_reordered += reorder_map[i];
    }
    TMC_EXPECT_EQ(sum_orig, sum_reordered);
    /* every index appears once (MICM checks only the sum) */
    bool seen[order] = { false };
    for (int i = 0; i < order; ++i)
      if (reorder_map[i] >= 0 && reorder_map[i] < order)
        seen[reorder_map[i]] = true;
    for (int i = 0; i < order; ++i)
      TMC_EXPECT_TRUE(seen[i]);
    TMC_EXPECT_TRUE(tmc_sparse_flat_block_size(orig_L) + tmc_sparse_flat_block_size(orig_U) >
                    tmc_sparse_flat_block_size(reordered_L) + tmc_sparse_flat_block_size(reordered_U));
    tmc_sparse_free(orig_L);
    tmc_sparse_free(orig_U);
    tmc_sparse_free(reordered_L);
    tmc_sparse_free(reordered_U);
  }
  tmc_lu_free(orig_lu);
  tmc_lu_free(reordered_lu);
  tmc_sparse_free(orig_jac);
  tmc_sparse_free(reordered_jac);
}

static const double initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };

TMC_TEST(dense_matrix_standard_ordering)
{
  test_dense_matrix(1);
  return TMC_TEST_PASS;
}

TMC_TEST(random_matrix_standard_ordering)
{
  test_random_matrix(1, 5);
  return TMC_TEST_PASS;
}

TMC_TEST(diagonal_matrix_standard_ordering)
{
  test_diagonal_matrix(1, 5);
  return TMC_TEST_PASS;
}

TMC_TEST(diagonal_markowitz_reorder)
{
  test_markowitz_reordering(1);
  return TMC_TEST_PASS;
}

TMC_TEST(standard_ordering_agnostic_to_initial_value)
{
  for (int v = 0; v < 5; ++v)
    test_extreme_initial_value(1, 5, initial_values[v]);
  return TMC_TEST_PASS;
}

TMC_TEST(dense_matrix_vector_ordering)
{
  for (int L = 1; L <= 4; ++L)
    test_dense_matrix(L);
  return TMC_TEST_PASS;
}

TMC_TEST(random_matrix_vector_ordering)
{
  for (int L = 1; L <= 4; ++L)
    test_random_matrix(L, 5);
  return TMC_TEST_PASS;
}

TMC_TEST(vector_ordering_agnostic_to_initial_value)
{
  /* MICM's block counts per width: 1, 2, 5, 5 */
  static const int n_blocks[4] = { 1, 2, 5, 5 };
  for (int v = 0; v < 5; ++v)
    for (int L = 1; L <= 4; ++L)
      test_extreme_initial_value(L, n_blocks[L - 1], initial_values[v]);
  return TMC_TEST_PASS;
}

TMC_TEST(diagonal_matrix_vector_ordering)
{
  for (int L = 1; L <= 4; ++L)
    test_diagonal_matrix(L, 5);
  return TMC_TEST_PASS;
}

TMC_TEST(vector_diagonal_markowitz_reordering)
{
  for (int L = 1; L <= 4; ++L)
    test_markowitz_reordering(L);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("linear_solver:LinearSolver.DenseMatrixStandardOrdering", dense_matrix_standard_ordering),
  TMC_CASE("linear_solver:LinearSolver.RandomMatrixStandardOrdering", random_matrix_standard_ordering),
  TMC_CASE("linear_solver:LinearSolver.DiagonalMatrixStandardOrdering", diagonal_matrix_standard_ordering),
  TMC_CASE("linear_solver:LinearSolver.DiagonalMarkowitzReorder", diagonal_markowitz_reorder),
  TMC_CASE("linear_solver:LinearSolver.StandardOrderingAgnosticToInitialValue", standard_ordering_agnostic_to_initial_value),
  TMC_CASE("linear_solver:LinearSolver.DenseMatrixVectorOrdering", dense_matrix_vector_ordering),
  TMC_CASE("linear_solver:LinearSolver.RandomMatrixVectorOrdering", random_matrix_vector_ordering),
  TMC_CASE("linear_solver:LinearSolver.VectorOrderingAgnosticToInitialValue", vector_ordering_agnostic_to_initial_value),
  TMC_CASE("linear_solver:LinearSolver.DiagonalMatrixVectorOrdering", diagonal_matrix_vector_ordering),
  TMC_CASE("linear_solver:LinearSolver.VectorDiagonalMarkowitzReordering", vector_diagonal_markowitz_reordering),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_linear_solver(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
