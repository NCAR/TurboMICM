/* MICM test/integration/aerosol_model_policy.hpp: the checks the
 * aerosol_model_rosenbrock and aerosol_model_backward_euler cases share, each
 * run for one of MICM's builders. Failures are recorded with the builder's
 * name (tmc_test_fail); a check MICM asserts ends that builder's run. */
#ifndef TMC_TEST_AEROSOL_MODEL_POLICY_H
#define TMC_TEST_AEROSOL_MODEL_POLICY_H

#include "tmc.h"

typedef enum aerosol_solver
{
  AEROSOL_ROSENBROCK_2STAGE,
  AEROSOL_ROSENBROCK_3STAGE,
  AEROSOL_ROSENBROCK_4STAGE,
  AEROSOL_ROSENBROCK_4STAGE_DAE,
  AEROSOL_ROSENBROCK_6STAGE_DAE,
  AEROSOL_BACKWARD_EULER
} aerosol_solver;

/* One of MICM's builders: its solver parameters, the vector width L (1 for
 * Matrix and SparseMatrixStandardOrdering), the LU algorithm (MICM's default
 * LuDecomposition is Doolittle; the in-place variants are the library's only
 * kind), and the element numbering of its sparse matrices (CSC for the
 * CompressedSparseColumn policies). The library's builder always numbers the
 * solver's Jacobian row by row; the numbering changes no value it computes,
 * and the Jacobian check hands the model a matrix in the policy's numbering. */
typedef struct aerosol_policy
{
  const char *name;
  aerosol_solver solver;
  int width;
  tmc_lu_algorithm lu;
  tmc_sparse_ordering ordering;
} aerosol_policy;

void test_state_includes_stub_aerosol_model(const aerosol_policy *p);
void test_update_state_with_stub_aerosol_model(const aerosol_policy *p);
void test_update_multi_cell_state_with_stub_aerosol_model(const aerosol_policy *p);
void test_single_cell_forcing_with_stub_aerosol_model(const aerosol_policy *p);
void test_single_cell_jacobian_with_stub_aerosol_model(const aerosol_policy *p);
void test_solve_with_stub_aerosol_model_1(const aerosol_policy *p, double base_relative_tolerance);
void test_solve_with_two_stub_aerosol_models(const aerosol_policy *p, double base_relative_tolerance);
void test_solve_with_stub_aerosol_model_1_multi_cell(const aerosol_policy *p, double base_relative_tolerance);
void test_solve_with_two_stub_aerosol_models_multi_cell(const aerosol_policy *p, double base_relative_tolerance);

#endif
