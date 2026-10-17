/* MICM test/unit/solver/test_lu_decomposition_policy.hpp: decompositions into
 * separate L and U matrices (tmc_lu_decompose), shared by the Doolittle and
 * Mozart test files. Each check factors a matrix and compares L U with A
 * element by element. The random matrices are MICM's own (libcxx_random.h). */
#ifndef TMC_TEST_LU_DECOMPOSITION_POLICY_H
#define TMC_TEST_LU_DECOMPOSITION_POLICY_H

#include "harness.h"
#include "libcxx_random.h"

/* MICM's CheckResults: A = L U within tol, L lower and U upper triangular. */
static void lu_check_results(const tmc_sparse *A, const tmc_sparse *Lm, const tmc_sparse *U, double tol)
{
  TMC_EXPECT_EQ(A->n_blocks, Lm->n_blocks);
  TMC_EXPECT_EQ(A->n_blocks, U->n_blocks);
  const int n = A->n;
  for (int b = 0; b < A->n_blocks; ++b)
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
      {
        double result = 0.0;
        for (int k = 0; k < n; ++k)
          if (!(tmc_sparse_is_zero(Lm, i, k) || tmc_sparse_is_zero(U, k, j)))
            result += tmc_sparse_get(Lm, b, i, k) * tmc_sparse_get(U, b, k, j);
        TMC_EXPECT_TRUE(i >= j || tmc_sparse_is_zero(Lm, i, j));
        TMC_EXPECT_TRUE(j >= i || tmc_sparse_is_zero(U, i, j));
        double expected = tmc_sparse_is_zero(A, i, j) ? 0.0 : tmc_sparse_get(A, b, i, j);
        TMC_EXPECT_NEAR(result, expected, tol);
      }
}

/* Create the decomposition of A's pattern, L and U initialized to
 * initial_value, decompose, and check. */
static void lu_decompose_and_check(const tmc_sparse *A, tmc_lu_algorithm algorithm, double initial_value, double tol)
{
  tmc_error err;
  tmc_error_clear(&err);
  tmc_lu *lu = tmc_lu_create(A, algorithm, &err);
  TMC_EXPECT_TRUE(lu != NULL);
  if (!lu)
    return;
  TMC_EXPECT_TRUE(tmc_lu_algorithm_of(lu) == algorithm);
  tmc_sparse *Lm = tmc_lu_lower_matrix(lu, A->n_blocks, initial_value);
  tmc_sparse *U = tmc_lu_upper_matrix(lu, A->n_blocks, initial_value);
  tmc_lu_decompose(lu, A, Lm, U);
  lu_check_results(A, Lm, U, tol);
  tmc_sparse_free(Lm);
  tmc_sparse_free(U);
  tmc_lu_free(lu);
}

/* A 3 x 3 dense example (geeksforgeeks' Doolittle example), one block. */
static void lu_test_dense_matrix(tmc_lu_algorithm algorithm, int L)
{
  tmc_sparse_builder *builder = tmc_sparse_builder_create(3);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      tmc_sparse_builder_with_element(builder, i, j, NULL);
  tmc_sparse_builder_set_initial_value(builder, 0.0);
  tmc_sparse *A = tmc_sparse_build(builder, TMC_SPARSE_CSR, L, false);
  tmc_sparse_builder_free(builder);
  const double values[3][3] = { { 2, -1, -2 }, { -4, 6, 3 }, { -4, -2, 8 } };
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      tmc_sparse_set(A, 0, i, j, values[i][j]);
  lu_decompose_and_check(A, algorithm, 0.0, 1.0e-10);
  tmc_sparse_free(A);
}

/* A random 10 x 10 pattern (diagonal always present) of lognormal values,
 * different in every block. */
static void lu_test_random_matrix(tmc_lu_algorithm algorithm, int L, int n_blocks)
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
  tmc_sparse *A = tmc_sparse_build(builder, TMC_SPARSE_CSR, L, false);
  tmc_sparse_builder_free(builder);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      if (!tmc_sparse_is_zero(A, i, j))
        for (int b = 0; b < n_blocks; ++b)
          tmc_sparse_set(A, b, i, j, tmc_lognormal_stream_next(&get_double));
  lu_decompose_and_check(A, algorithm, 0.0, 1.0e-9);
  tmc_sparse_free(A);
}

/* The random 10 x 10 matrix with A, L, and U initialized to an extreme value
 * (one block's values copied to every block): the decomposition must
 * overwrite every element of L and U. */
static void lu_test_extreme_value_initialization(tmc_lu_algorithm algorithm, int L, int n_blocks, double initial_value)
{
  tmc_bool_stream gen_bool = tmc_bool_stream_create();
  tmc_lognormal_stream get_double = tmc_lognormal_stream_create(-2.0, 2.0);
  const int size = 10;
  tmc_sparse_builder *builder = tmc_sparse_builder_create(size);
  tmc_sparse_builder_set_number_of_blocks(builder, n_blocks);
  tmc_sparse_builder_set_initial_value(builder, initial_value);
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
      if (i == j || tmc_bool_stream_next(&gen_bool))
        tmc_sparse_builder_with_element(builder, i, j, NULL);
  tmc_sparse *A = tmc_sparse_build(builder, TMC_SPARSE_CSR, L, false);
  tmc_sparse_builder_free(builder);
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
      if (!tmc_sparse_is_zero(A, i, j))
      {
        double v = tmc_lognormal_stream_next(&get_double);
        for (int b = 0; b < n_blocks; ++b)
          tmc_sparse_set(A, b, i, j, v);
      }
  lu_decompose_and_check(A, algorithm, initial_value, 1.0e-9);
  tmc_sparse_free(A);
}

/* A 6 x 6 diagonal matrix of lognormal(-2, 4) values. */
static void lu_test_diagonal_matrix(tmc_lu_algorithm algorithm, int L, int n_blocks)
{
  tmc_lognormal_stream get_double = tmc_lognormal_stream_create(-2.0, 4.0);
  tmc_sparse_builder *builder = tmc_sparse_builder_create(6);
  tmc_sparse_builder_set_number_of_blocks(builder, n_blocks);
  tmc_sparse_builder_set_initial_value(builder, 0.0);
  for (int i = 0; i < 6; ++i)
    tmc_sparse_builder_with_element(builder, i, i, NULL);
  tmc_sparse *A = tmc_sparse_build(builder, TMC_SPARSE_CSR, L, false);
  tmc_sparse_builder_free(builder);
  for (int i = 0; i < 6; ++i)
    for (int b = 0; b < n_blocks; ++b)
      tmc_sparse_set(A, b, i, i, tmc_lognormal_stream_next(&get_double));
  lu_decompose_and_check(A, algorithm, 0.0, 1.0e-10);
  tmc_sparse_free(A);
}

/* The cases of test_lu_decomposition_{doolittle,mozart}.cpp. MICM's standard
 * ordering is L = 1; its vector orderings are L = 1 to 4. */
static void lu_dense_matrix_standard_ordering(tmc_lu_algorithm algorithm)
{
  lu_test_dense_matrix(algorithm, 1);
}

static void lu_random_matrix_standard_ordering(tmc_lu_algorithm algorithm)
{
  lu_test_random_matrix(algorithm, 1, 1);
  lu_test_random_matrix(algorithm, 1, 5);
}

static void lu_diagonal_matrix_standard_ordering(tmc_lu_algorithm algorithm)
{
  lu_test_diagonal_matrix(algorithm, 1, 5);
}

static const double lu_initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };

static void lu_agnostic_to_initial_value_standard_ordering(tmc_lu_algorithm algorithm)
{
  for (int v = 0; v < 5; ++v)
    lu_test_extreme_value_initialization(algorithm, 1, 5, lu_initial_values[v]);
}

static void lu_dense_matrix_vector_ordering(tmc_lu_algorithm algorithm)
{
  for (int L = 1; L <= 4; ++L)
    lu_test_dense_matrix(algorithm, L);
}

static void lu_random_matrix_vector_ordering(tmc_lu_algorithm algorithm)
{
  for (int L = 1; L <= 4; ++L)
    lu_test_random_matrix(algorithm, L, 5);
}

static void lu_diagonal_matrix_vector_ordering(tmc_lu_algorithm algorithm)
{
  for (int L = 1; L <= 4; ++L)
    lu_test_diagonal_matrix(algorithm, L, 5);
}

static void lu_agnostic_to_initial_value_vector_ordering(tmc_lu_algorithm algorithm)
{
  for (int v = 0; v < 5; ++v)
    for (int L = 1; L <= 4; ++L)
      lu_test_extreme_value_initialization(algorithm, L, 5, lu_initial_values[v]);
}

#endif
