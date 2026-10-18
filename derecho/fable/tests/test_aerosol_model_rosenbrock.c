/* MICM test/integration/test_aerosol_model_rosenbrock.cpp: the aerosol model
 * checks (aerosol_model_policy.c) for each of the test's 27 Rosenbrock
 * builders. CpuSolverBuilder's defaults are Matrix, SparseMatrixStandardOrdering
 * (L = 1) and LuDecomposition, which is Doolittle. */
#include "harness.h"

#include "aerosol_model_policy.h"

#include <string.h>

#define ROS(name, solver, L, lu, ordering) { name, solver, L, lu, ordering }
#define ROS3(name, L, lu, ordering) ROS(name, AEROSOL_ROSENBROCK_3STAGE, L, lu, ordering)

static const aerosol_policy policies[] = {
  ROS("rosenbrock_2stage", AEROSOL_ROSENBROCK_2STAGE, 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS("rosenbrock_3stage", AEROSOL_ROSENBROCK_3STAGE, 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS("rosenbrock_4stage", AEROSOL_ROSENBROCK_4STAGE, 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS("rosenbrock_4stage_da", AEROSOL_ROSENBROCK_4STAGE_DAE, 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS("rosenbrock_6stage_da", AEROSOL_ROSENBROCK_6STAGE_DAE, 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_1", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_2", 2, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_3", 3, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_4", 4, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_standard_doolittle", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_doolittle_1", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_doolittle_2", 2, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_doolittle_3", 3, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_doolittle_4", 4, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  ROS3("rosenbrock_standard_mozart", 1, TMC_LU_MOZART, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_mozart_1", 1, TMC_LU_MOZART, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_mozart_2", 2, TMC_LU_MOZART, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_mozart_3", 3, TMC_LU_MOZART, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_mozart_4", 4, TMC_LU_MOZART, TMC_SPARSE_CSR),
  ROS3("rosenbrock_vector_doolittle_csc_1", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSC),
  ROS3("rosenbrock_vector_doolittle_csc_2", 2, TMC_LU_DOOLITTLE, TMC_SPARSE_CSC),
  ROS3("rosenbrock_vector_doolittle_csc_3", 3, TMC_LU_DOOLITTLE, TMC_SPARSE_CSC),
  ROS3("rosenbrock_vector_doolittle_csc_4", 4, TMC_LU_DOOLITTLE, TMC_SPARSE_CSC),
  ROS3("rosenbrock_vector_mozart_csc_1", 1, TMC_LU_MOZART, TMC_SPARSE_CSC),
  ROS3("rosenbrock_vector_mozart_csc_2", 2, TMC_LU_MOZART, TMC_SPARSE_CSC),
  ROS3("rosenbrock_vector_mozart_csc_3", 3, TMC_LU_MOZART, TMC_SPARSE_CSC),
  ROS3("rosenbrock_vector_mozart_csc_4", 4, TMC_LU_MOZART, TMC_SPARSE_CSC),
};
enum
{
  n_policies = (int)(sizeof policies / sizeof policies[0])
};

/* MICM's default base relative tolerance, and the two-stage builder's in the
 * two-model solves */
static double two_model_tolerance(const aerosol_policy *p)
{
  return strcmp(p->name, "rosenbrock_2stage") == 0 ? 5e-4 : 5e-5;
}

TMC_TEST(state_includes_stub_aerosol_model)
{
  for (int i = 0; i < n_policies; ++i)
    test_state_includes_stub_aerosol_model(&policies[i]);
  return TMC_TEST_PASS;
}

TMC_TEST(can_update_state_with_stub_aerosol_model)
{
  for (int i = 0; i < n_policies; ++i)
    test_update_state_with_stub_aerosol_model(&policies[i]);
  return TMC_TEST_PASS;
}

TMC_TEST(can_update_multi_cell_state_with_stub_aerosol_model)
{
  for (int i = 0; i < n_policies; ++i)
    test_update_multi_cell_state_with_stub_aerosol_model(&policies[i]);
  return TMC_TEST_PASS;
}

TMC_TEST(single_grid_cell_forcing)
{
  for (int i = 0; i < n_policies; ++i)
    test_single_cell_forcing_with_stub_aerosol_model(&policies[i]);
  return TMC_TEST_PASS;
}

TMC_TEST(single_grid_cell_jacobian)
{
  for (int i = 0; i < n_policies; ++i)
    test_single_cell_jacobian_with_stub_aerosol_model(&policies[i]);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_single_grid_cell_model_1)
{
  for (int i = 0; i < n_policies; ++i)
    test_solve_with_stub_aerosol_model_1(&policies[i], 5e-5);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_single_grid_cell_two_models)
{
  for (int i = 0; i < n_policies; ++i)
    test_solve_with_two_stub_aerosol_models(&policies[i], two_model_tolerance(&policies[i]));
  return TMC_TEST_PASS;
}

TMC_TEST(solve_multi_grid_cell_model_1)
{
  for (int i = 0; i < n_policies; ++i)
    test_solve_with_stub_aerosol_model_1_multi_cell(&policies[i], 5e-5);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_multi_grid_cell_two_models)
{
  for (int i = 0; i < n_policies; ++i)
    test_solve_with_two_stub_aerosol_models_multi_cell(&policies[i], two_model_tolerance(&policies[i]));
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.StateIncludesStubAerosolModel", state_includes_stub_aerosol_model),
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.CanUpdateStateWithStubAerosolModel", can_update_state_with_stub_aerosol_model),
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.CanUpdateMultiCellStateWithStubAerosolModel", can_update_multi_cell_state_with_stub_aerosol_model),
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.CanCalculateSingleGridCellForcingWithStubAerosolModel", single_grid_cell_forcing),
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.CanCalculateSingleGridCellJacobianWithStubAerosolModel", single_grid_cell_jacobian),
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.CanSolveSingleGridCellWithStubAerosolModel1", solve_single_grid_cell_model_1),
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.CanSolveSingleGridCellWithTwoStubAerosolModels", solve_single_grid_cell_two_models),
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.CanSolveMultiGridCellWithStubAerosolModel1", solve_multi_grid_cell_model_1),
  TMC_CASE("aerosol_model_rosenbrock:AerosolModelIntegration.CanSolveMultiGridCellWithTwoStubAerosolModels", solve_multi_grid_cell_two_models),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_aerosol_model_rosenbrock(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
