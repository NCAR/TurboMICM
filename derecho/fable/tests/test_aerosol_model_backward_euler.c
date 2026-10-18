/* MICM test/integration/test_aerosol_model_backward_euler.cpp: the aerosol
 * model checks (aerosol_model_policy.c) for each of the test's 29 backward
 * Euler builders. CpuSolverBuilder's defaults are Matrix,
 * SparseMatrixStandardOrdering (L = 1) and LuDecomposition, which is
 * Doolittle; the library's LU is always in place, so the in-place builders
 * differ from the others only in name. */
#include "harness.h"

#include "aerosol_model_policy.h"

#define BE(name, L, lu, ordering) { name, AEROSOL_BACKWARD_EULER, L, lu, ordering }

static const aerosol_policy policies[] = {
  BE("backward_euler", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_1", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_2", 2, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_3", 3, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_4", 4, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_doolittle_1", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_doolittle_2", 2, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_doolittle_3", 3, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_doolittle_4", 4, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_mozart_1", 1, TMC_LU_MOZART, TMC_SPARSE_CSR),
  BE("backward_euler_vector_mozart_2", 2, TMC_LU_MOZART, TMC_SPARSE_CSR),
  BE("backward_euler_vector_mozart_3", 3, TMC_LU_MOZART, TMC_SPARSE_CSR),
  BE("backward_euler_vector_mozart_4", 4, TMC_LU_MOZART, TMC_SPARSE_CSR),
  BE("backward_euler_vector_doolittle_csc_1", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSC),
  BE("backward_euler_vector_doolittle_csc_2", 2, TMC_LU_DOOLITTLE, TMC_SPARSE_CSC),
  BE("backward_euler_vector_doolittle_csc_3", 3, TMC_LU_DOOLITTLE, TMC_SPARSE_CSC),
  BE("backward_euler_vector_doolittle_csc_4", 4, TMC_LU_DOOLITTLE, TMC_SPARSE_CSC),
  BE("backward_euler_vector_mozart_csc_1", 1, TMC_LU_MOZART, TMC_SPARSE_CSC),
  BE("backward_euler_vector_mozart_csc_2", 2, TMC_LU_MOZART, TMC_SPARSE_CSC),
  BE("backward_euler_vector_mozart_csc_3", 3, TMC_LU_MOZART, TMC_SPARSE_CSC),
  BE("backward_euler_vector_mozart_csc_4", 4, TMC_LU_MOZART, TMC_SPARSE_CSC),
  BE("backward_euler_vector_doolittle_in_place_1", 1, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_doolittle_in_place_2", 2, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_doolittle_in_place_3", 3, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_doolittle_in_place_4", 4, TMC_LU_DOOLITTLE, TMC_SPARSE_CSR),
  BE("backward_euler_vector_mozart_in_place_1", 1, TMC_LU_MOZART, TMC_SPARSE_CSR),
  BE("backward_euler_vector_mozart_in_place_2", 2, TMC_LU_MOZART, TMC_SPARSE_CSR),
  BE("backward_euler_vector_mozart_in_place_3", 3, TMC_LU_MOZART, TMC_SPARSE_CSR),
  BE("backward_euler_vector_mozart_in_place_4", 4, TMC_LU_MOZART, TMC_SPARSE_CSR),
};
enum
{
  n_policies = (int)(sizeof policies / sizeof policies[0])
};

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
    test_solve_with_stub_aerosol_model_1(&policies[i], 2e-4);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_single_grid_cell_two_models)
{
  for (int i = 0; i < n_policies; ++i)
    test_solve_with_two_stub_aerosol_models(&policies[i], 1e-1);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_multi_grid_cell_model_1)
{
  for (int i = 0; i < n_policies; ++i)
    test_solve_with_stub_aerosol_model_1_multi_cell(&policies[i], 2e-4);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_multi_grid_cell_two_models)
{
  for (int i = 0; i < n_policies; ++i)
    test_solve_with_two_stub_aerosol_models_multi_cell(&policies[i], 1e-1);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.StateIncludesStubAerosolModel", state_includes_stub_aerosol_model),
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.CanUpdateStateWithStubAerosolModel", can_update_state_with_stub_aerosol_model),
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.CanUpdateMultiCellStateWithStubAerosolModel", can_update_multi_cell_state_with_stub_aerosol_model),
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.CanCalculateSingleGridCellForcingWithStubAerosolModel", single_grid_cell_forcing),
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.CanCalculateSingleGridCellJacobianWithStubAerosolModel", single_grid_cell_jacobian),
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.CanSolveSingleGridCellWithStubAerosolModel1", solve_single_grid_cell_model_1),
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.CanSolveSingleGridCellWithTwoStubAerosolModels", solve_single_grid_cell_two_models),
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.CanSolveMultiGridCellWithStubAerosolModel1", solve_multi_grid_cell_model_1),
  TMC_CASE("aerosol_model_backward_euler:AerosolModelIntegration.CanSolveMultiGridCellWithTwoStubAerosolModels", solve_multi_grid_cell_two_models),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_aerosol_model_backward_euler(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
