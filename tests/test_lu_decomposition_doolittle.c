/* MICM test/unit/solver/test_lu_decomposition_doolittle.cpp: Doolittle
 * decomposition into separate L and U matrices (lu_decomposition_policy.h). */
#include "lu_decomposition_policy.h"

#define ALGORITHM TMC_LU_DOOLITTLE

TMC_TEST(dense_matrix_standard_ordering)
{
  lu_dense_matrix_standard_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(random_matrix_standard_ordering)
{
  lu_random_matrix_standard_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(diagonal_matrix_standard_ordering)
{
  lu_diagonal_matrix_standard_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(agnostic_to_initial_value_standard_ordering)
{
  lu_agnostic_to_initial_value_standard_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(dense_matrix_vector_ordering)
{
  lu_dense_matrix_vector_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(random_matrix_vector_ordering)
{
  lu_random_matrix_vector_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(diagonal_matrix_vector_ordering)
{
  lu_diagonal_matrix_vector_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

TMC_TEST(agnostic_to_initial_value_vector_ordering)
{
  lu_agnostic_to_initial_value_vector_ordering(ALGORITHM);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("lu_decomposition_doolittle:LuDecompositionDoolittle.DenseMatrixStandardOrdering", dense_matrix_standard_ordering),
  TMC_CASE("lu_decomposition_doolittle:LuDecompositionDoolittle.RandomMatrixStandardOrdering", random_matrix_standard_ordering),
  TMC_CASE("lu_decomposition_doolittle:LuDecompositionDoolittle.DiagonalMatrixStandardOrdering", diagonal_matrix_standard_ordering),
  TMC_CASE("lu_decomposition_doolittle:LuDecompositionDoolittle.AgnosticToInitialValueStandardOrdering", agnostic_to_initial_value_standard_ordering),
  TMC_CASE("lu_decomposition_doolittle:LuDecompositionDoolittle.DenseMatrixVectorOrdering", dense_matrix_vector_ordering),
  TMC_CASE("lu_decomposition_doolittle:LuDecompositionDoolittle.RandomMatrixVectorOrdering", random_matrix_vector_ordering),
  TMC_CASE("lu_decomposition_doolittle:LuDecompositionDoolittle.DiagonalMatrixVectorOrdering", diagonal_matrix_vector_ordering),
  TMC_CASE("lu_decomposition_doolittle:LuDecompositionDoolittle.AgnosticToInitialValueVectorOrdering", agnostic_to_initial_value_vector_ordering),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_lu_decomposition_doolittle(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
