/* MICM test/unit/solver/test_lu_decomposition_in_place_policy.hpp: in-place
 * decompositions (tmc_lu_decompose_in_place), shared by the Doolittle and
 * Mozart in-place test files. Each check copies A into the combined LU
 * pattern, factors it in place, and compares L U with A element by element.
 * The random matrices are MICM's own (libcxx_random.h). The checks also
 * verify the pivot reciprocals the decomposition writes. */
#ifndef TMC_TEST_LU_DECOMPOSITION_IN_PLACE_POLICY_H
#define TMC_TEST_LU_DECOMPOSITION_IN_PLACE_POLICY_H

#include "harness.h"
#include "libcxx_random.h"

/* MICM's CheckResults: A = L U within tol, with L's unit diagonal implicit. */
static void lu_in_place_check_results(const tmc_sparse *A, const tmc_sparse *LU, double tol)
{
  TMC_EXPECT_EQ(A->n_blocks, LU->n_blocks);
  const int n = A->n;
  for (int b = 0; b < A->n_blocks; ++b)
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
      {
        double result = 0.0;
        for (int k = 0; k <= j; ++k)
        {
          if (i == k && !tmc_sparse_is_zero(LU, k, j))
            result += tmc_sparse_get(LU, b, k, j);
          else if (k < i && !(tmc_sparse_is_zero(LU, i, k) || tmc_sparse_is_zero(LU, k, j)))
            result += tmc_sparse_get(LU, b, i, k) * tmc_sparse_get(LU, b, k, j);
        }
        double expected = tmc_sparse_is_zero(A, i, j) ? 0.0 : tmc_sparse_get(A, b, i, j);
        TMC_EXPECT_NEAR(result, expected, tol);
      }
}

/* The reciprocal of every pivot of every block, as the decomposition stores them. */
static void lu_in_place_check_reciprocals(const tmc_sparse *LU, const double *reciprocals)
{
  const int n = LU->n, L = LU->L;
  for (int b = 0; b < LU->n_blocks; ++b)
    for (int i = 0; i < n; ++i)
    {
      double r = reciprocals[((size_t)(b / L) * (size_t)n + (size_t)i) * (size_t)L + (size_t)(b % L)];
      TMC_EXPECT_DOUBLE_EQ(r, 1.0 / tmc_sparse_get(LU, b, i, i));
    }
}

/* Copy A into an LU matrix initialized to initial_value (zeroing the fill
 * positions when zero_fill is set), factor it in place, and check. */
static void lu_in_place_decompose_and_check(const tmc_sparse *A, tmc_lu_algorithm algorithm, double initial_value, bool zero_fill, double tol)
{
  tmc_error err;
  tmc_error_clear(&err);
  tmc_lu *lu = tmc_lu_create(A, algorithm, &err);
  TMC_EXPECT_TRUE(lu != NULL);
  if (!lu)
    return;
  TMC_EXPECT_TRUE(tmc_lu_algorithm_of(lu) == algorithm);
  tmc_sparse *ALU = tmc_lu_matrix(lu, A->n_blocks, initial_value);
  const int n = A->n;
  for (int b = 0; b < A->n_blocks; ++b)
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
      {
        if (!tmc_sparse_is_zero(A, i, j))
          tmc_sparse_set(ALU, b, i, j, tmc_sparse_get(A, b, i, j));
        else if (zero_fill && !tmc_sparse_is_zero(ALU, i, j))
          tmc_sparse_set(ALU, b, i, j, 0.0);
      }
  double *reciprocals = calloc((size_t)ALU->n_groups * (size_t)ALU->L * (size_t)n, sizeof(double));
  tmc_lu_decompose_in_place(lu, ALU, reciprocals);
  lu_in_place_check_results(A, ALU, tol);
  lu_in_place_check_reciprocals(ALU, reciprocals);
  free(reciprocals);
  tmc_sparse_free(ALU);
  tmc_lu_free(lu);
}

/* A 3 x 3 dense example (geeksforgeeks' Doolittle example), one block. */
static void lu_in_place_test_dense_matrix(tmc_lu_algorithm algorithm, int L)
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
  lu_in_place_decompose_and_check(A, algorithm, 0.0, false, 1.0e-10);
  tmc_sparse_free(A);
}

/* A random pattern (diagonal always present) of lognormal values, one block's
 * values copied to every block, A built with initial_value. */
static tmc_sparse *lu_in_place_random_matrix(int size, int L, int n_blocks, double initial_value)
{
  tmc_bool_stream gen_bool = tmc_bool_stream_create();
  tmc_lognormal_stream get_double = tmc_lognormal_stream_create(-2.0, 2.0);
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
  return A;
}

static void lu_in_place_test_random_matrix(tmc_lu_algorithm algorithm, int L, int n_blocks)
{
  tmc_sparse *A = lu_in_place_random_matrix(10, L, n_blocks, 0.0);
  lu_in_place_decompose_and_check(A, algorithm, 0.0, false, 1.0e-12);
  tmc_sparse_free(A);
}

/* A random 5 x 5 matrix with A and the LU matrix initialized to an extreme
 * value; the fill positions are zeroed before the decomposition. */
static void lu_in_place_test_extreme_value_initialization(tmc_lu_algorithm algorithm, int L, int n_blocks, double initial_value)
{
  tmc_sparse *A = lu_in_place_random_matrix(5, L, n_blocks, initial_value);
  lu_in_place_decompose_and_check(A, algorithm, initial_value, true, 1.0e-12);
  tmc_sparse_free(A);
}

/* A 6 x 6 diagonal matrix of lognormal(-2, 4) values (the LU matrix zeroed first). */
static void lu_in_place_test_diagonal_matrix(tmc_lu_algorithm algorithm, int L, int n_blocks)
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
  lu_in_place_decompose_and_check(A, algorithm, 0.0, false, 1.0e-10);
  tmc_sparse_free(A);
}

/* The cases of test_lu_decomposition_{doolittle,mozart}_in_place.cpp. MICM's
 * standard ordering is L = 1; its vector orderings are L = 1 to 4. */
static void lu_in_place_dense_matrix_standard_ordering(tmc_lu_algorithm algorithm)
{
  lu_in_place_test_dense_matrix(algorithm, 1);
}

static void lu_in_place_random_matrix_standard_ordering(tmc_lu_algorithm algorithm)
{
  lu_in_place_test_random_matrix(algorithm, 1, 1);
  lu_in_place_test_random_matrix(algorithm, 1, 5);
}

static void lu_in_place_diagonal_matrix_standard_ordering(tmc_lu_algorithm algorithm)
{
  lu_in_place_test_diagonal_matrix(algorithm, 1, 5);
}

static const double lu_in_place_initial_values[5] = { -INFINITY, -1.0, 0.0, 1.0, INFINITY };

static void lu_in_place_agnostic_to_initial_value_standard_ordering(tmc_lu_algorithm algorithm)
{
  for (int v = 0; v < 5; ++v)
    lu_in_place_test_extreme_value_initialization(algorithm, 1, 5, lu_in_place_initial_values[v]);
}

static void lu_in_place_dense_matrix_vector_ordering(tmc_lu_algorithm algorithm)
{
  for (int L = 1; L <= 4; ++L)
    lu_in_place_test_dense_matrix(algorithm, L);
}

static void lu_in_place_random_matrix_vector_ordering(tmc_lu_algorithm algorithm)
{
  for (int L = 1; L <= 4; ++L)
    lu_in_place_test_random_matrix(algorithm, L, 5);
}

static void lu_in_place_diagonal_matrix_vector_ordering(tmc_lu_algorithm algorithm)
{
  for (int L = 1; L <= 4; ++L)
    lu_in_place_test_diagonal_matrix(algorithm, L, 5);
}

static void lu_in_place_agnostic_to_initial_value_vector_ordering(tmc_lu_algorithm algorithm)
{
  for (int v = 0; v < 5; ++v)
    for (int L = 1; L <= 4; ++L)
      lu_in_place_test_extreme_value_initialization(algorithm, L, 5, lu_in_place_initial_values[v]);
}

#endif
