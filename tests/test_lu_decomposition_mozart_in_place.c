/* MICM test/unit/solver/test_lu_decomposition_mozart_in_place.cpp:
 * in-place Mozart decomposition (lu_decomposition_in_place_policy.h). */
#include "lu_decomposition_in_place_policy.h"

#define ALGORITHM TMC_LU_MOZART

TMC_TEST(dense_matrix_standard_ordering)
{
  lu_in_place_dense_matrix_standard_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(random_matrix_standard_ordering)
{
  lu_in_place_random_matrix_standard_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(diagonal_matrix_standard_ordering)
{
  lu_in_place_diagonal_matrix_standard_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(agnostic_to_initial_value_standard_ordering)
{
  lu_in_place_agnostic_to_initial_value_standard_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(dense_matrix_vector_ordering)
{
  lu_in_place_dense_matrix_vector_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(random_matrix_vector_ordering)
{
  lu_in_place_random_matrix_vector_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(diagonal_matrix_vector_ordering)
{
  lu_in_place_diagonal_matrix_vector_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(agnostic_to_initial_value_vector_ordering)
{
  lu_in_place_agnostic_to_initial_value_vector_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.DenseMatrixStandardOrdering", dense_matrix_standard_ordering),
  TMC_CASE("lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.RandomMatrixStandardOrdering", random_matrix_standard_ordering),
  TMC_CASE("lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.DiagonalMatrixStandardOrdering", diagonal_matrix_standard_ordering),
  TMC_CASE("lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.AgnosticToInitialValueStandardOrdering", agnostic_to_initial_value_standard_ordering),
  TMC_CASE("lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.DenseMatrixVectorOrdering", dense_matrix_vector_ordering),
  TMC_CASE("lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.RandomMatrixVectorOrdering", random_matrix_vector_ordering),
  TMC_CASE("lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.DiagonalMatrixVectorOrdering", diagonal_matrix_vector_ordering),
  TMC_CASE("lu_decomposition_mozart_in_place:LuDecompositionMozartInPlace.AgnosticToInitialValueVectorOrdering", agnostic_to_initial_value_vector_ordering),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_lu_decomposition_mozart_in_place(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
