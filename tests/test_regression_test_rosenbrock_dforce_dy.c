/* MICM test/regression/RosenbrockChapman/regression_test_dforce_dy.cpp with
 * regression_test_dforce_dy_policy.hpp: the general solver's Jacobian for the
 * Chapman mechanism against the hard-coded solver's 23 elements. The general
 * solver stores -dF/dy in its sparse matrix, which has four more elements
 * (species with net-zero stoichiometry in a reaction), always zero. */
#include "chapman_fixed.h"
#include "libcxx_random.h"

static tmc_solver *chapman_solver(int L)
{
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_width(builder, L);
  return cf_chapman_solver(builder);
}

/* TestJacobian (std::mt19937 seeded 739 and lognormal(-2, 4); MICM draws the
 * concentrations twice, then the rate constants). */
static void test_jacobian(tmc_solver *solver)
{
  if (!solver)
    return;
  tmc_mt19937 engine = tmc_mt19937_create(739);
  tmc_lognormal dist = tmc_lognormal_create(-2.0, 4.0);
  tmc_state *state = tmc_solver_create_state(solver, 3);
  tmc_dense *variables = tmc_state_variables(state);
  tmc_dense *rate_constants = tmc_state_rate_constants(state);
  for (size_t k = 0; k < variables->size; ++k)
    variables->data[k] = tmc_mt19937_lognormal(&dist, &engine);
  for (size_t k = 0; k < variables->size; ++k)
    variables->data[k] = tmc_mt19937_lognormal(&dist, &engine);
  for (size_t k = 0; k < rate_constants->size; ++k)
    rate_constants->data[k] = tmc_mt19937_lognormal(&dist, &engine);

  tmc_sparse *jacobian = tmc_state_jacobian(state);
  tmc_sparse_fill(jacobian, 0.0);
  tmc_solver_subtract_jacobian_terms(solver, state, variables, jacobian);

  for (int i = 0; i < 3; ++i)
  {
    double cell_rate_constants[CF_N_RATE_CONSTANTS], cell_variables[CF_N_SPECIES], fixed_jacobian[CF_N_ELEMENTS];
    tmc_dense_get_row(rate_constants, i, cell_rate_constants);
    for (int j = 0; j < CF_N_SPECIES; ++j)
      cell_variables[j] = tmc_state_get(state, i, tmc_state_variable_index(state, cf_species_names[j]));
    cf_dforce_dy(cell_rate_constants, cell_variables, fixed_jacobian);
    for (int j = 0; j < CF_N_ELEMENTS; ++j)
    {
      const int modern_row = tmc_state_variable_index(state, cf_species_names[cf_elements[j][0]]);
      const int modern_col = tmc_state_variable_index(state, cf_species_names[cf_elements[j][1]]);
      if (tmc_sparse_is_zero(jacobian, modern_row, modern_col))
      {
        tmc_test_fail(__FILE__, __LINE__, "L = %d: element (%s, %s) is not in the Jacobian", tmc_state_width(state), cf_species_names[cf_elements[j][0]], cf_species_names[cf_elements[j][1]]);
        continue;
      }
      /* SubtractJacobianTerms stores -dF/dy; DforceDy returns +dF/dy */
      const double a = tmc_sparse_get(jacobian, i, modern_row, modern_col);
      const double b = fixed_jacobian[j];
      if (!(fabs(a - -b) <= (fabs(a) + fabs(b)) * 1.0e-8 + 1.0e-12))
        tmc_test_fail(__FILE__, __LINE__, "L = %d, cell %d, element %d: %.17g vs %.17g", tmc_state_width(state), i, j, a, -b);
    }
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
}

TMC_TEST(jacobian)
{
  test_jacobian(chapman_solver(1));
  return TMC_TEST_PASS;
}

TMC_TEST(vector_jacobian)
{
  for (int L = 1; L <= 4; ++L)
    test_jacobian(chapman_solver(L));
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.Jacobian", jacobian),
  TMC_CASE("regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.VectorJacobian", vector_jacobian),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_regression_test_rosenbrock_dforce_dy(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
