/* MICM test/regression/RosenbrockChapman/regression_test_solve.cpp with
 * regression_test_solve_policy.hpp: 500 s of the Chapman mechanism from
 * lognormal concentrations in three cells, the general solver against the
 * hard-coded solver, to a relative tolerance of 1e-4 (1e-2 for the two-stage
 * method). MICM's CpuSolverBuilder is width L = 1; its VectorBuilder<L> is L. */
#include "chapman_fixed.h"
#include "libcxx_random.h"

static tmc_solver *chapman_solver(tmc_rosenbrock_parameters parameters, int L)
{
  tmc_builder *builder = tmc_builder_create_rosenbrock(parameters);
  tmc_builder_set_width(builder, L);
  return cf_chapman_solver(builder);
}

static void test_solve(tmc_solver *solver, double relative_tolerance)
{
  if (!solver)
    return;
  /* std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine()) */
  tmc_lognormal_stream get_double = tmc_lognormal_stream_create(-2.0, 2.0);
  tmc_state *state = tmc_solver_create_state(solver, 3);

  static const double photo_rates[3][3] = { { 1.0e-4, 1.0e-5, 1.0e-6 }, { 3.2e-4, 7.3e-5, 3.2e-6 }, { 5.2e-4, 8.2e-5, 4.6e-6 } };
  static const double temperatures[3] = { 284.19, 215.02, 299.31 }; /* [K] */
  static const double pressures[3] = { 101245.0, 100789.2, 101398.0 }; /* [Pa] */
  tmc_conditions *conditions = tmc_state_conditions(state);
  for (int i = 0; i < 3; ++i)
  {
    tmc_dense_set_row(tmc_state_custom_rate_parameters(state), i, photo_rates[i]);
    conditions[i].temperature = temperatures[i];
    conditions[i].pressure = pressures[i];
  }
  double variables[3][CF_N_SPECIES];
  double abs_tol = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < CF_N_SPECIES; ++j)
    {
      variables[i][j] = tmc_lognormal_stream_next(&get_double);
      abs_tol = fmax(abs_tol, variables[i][j]);
    }
  for (int i = 0; i < 3; ++i)
    tmc_dense_set_row(tmc_state_variables(state), i, variables[i]);
  abs_tol *= 1.0e-12;

  /* run solvers */
  tmc_solver_update_state_parameters(solver, state);
  tmc_result results = tmc_solver_solve(solver, 500.0, state);
  (void)results;
  cf_result fixed_results[3];
  for (int i = 0; i < 3; ++i)
  {
    double fixed_rate_constants[CF_N_RATE_CONSTANTS], fixed_variables[CF_N_SPECIES];
    cf_update_state(conditions[i].temperature, conditions[i].pressure, photo_rates[i], fixed_rate_constants);
    for (int j = 0; j < CF_N_SPECIES; ++j)
      fixed_variables[j] = variables[i][tmc_state_variable_index(state, cf_species_names[j])];
    fixed_results[i] = cf_solve(0.0, 500.0, fixed_variables, fixed_rate_constants);
  }

  /* compare results */
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < CF_N_SPECIES; ++j)
    {
      const double a = tmc_state_get(state, i, tmc_state_variable_index(state, cf_species_names[j]));
      const double b = fixed_results[i].y[j];
      if (!(fabs(a - b) <= (fabs(a) + fabs(b)) * relative_tolerance + abs_tol))
        tmc_test_fail(__FILE__, __LINE__, "%d stages, L = %d, cell %d, %s: %.17g vs %.17g (hard-coded solver reached t = %g)", tmc_solver_rosenbrock_parameters(solver)->stages, tmc_state_width(state), i, cf_species_names[j], a, b, fixed_results[i].final_time);
    }
  tmc_state_free(state);
  tmc_solver_free(solver);
}

TMC_TEST(two_stage_solve)
{
  test_solve(chapman_solver(tmc_rosenbrock_two_stage(), 1), 1.0e-2);
  return TMC_TEST_PASS;
}

TMC_TEST(three_stage_solve)
{
  test_solve(chapman_solver(tmc_rosenbrock_three_stage(), 1), 1.0e-4);
  return TMC_TEST_PASS;
}

TMC_TEST(four_stage_solve)
{
  test_solve(chapman_solver(tmc_rosenbrock_four_stage(), 1), 1.0e-4);
  return TMC_TEST_PASS;
}

TMC_TEST(four_stage_da_solve)
{
  test_solve(chapman_solver(tmc_rosenbrock_four_stage_dae(), 1), 1.0e-4);
  return TMC_TEST_PASS;
}

TMC_TEST(six_stage_da_solve)
{
  test_solve(chapman_solver(tmc_rosenbrock_six_stage_dae(), 1), 1.0e-4);
  return TMC_TEST_PASS;
}

TMC_TEST(vector_solve)
{
  for (int L = 1; L <= 4; ++L)
    test_solve(chapman_solver(tmc_rosenbrock_three_stage(), L), 1.0e-4);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("regression_test_solve:RegressionRosenbrock.TwoStageSolve", two_stage_solve),
  TMC_CASE("regression_test_solve:RegressionRosenbrock.ThreeStageSolve", three_stage_solve),
  TMC_CASE("regression_test_solve:RegressionRosenbrock.FourStageSolve", four_stage_solve),
  TMC_CASE("regression_test_solve:RegressionRosenbrock.FourStageDASolve", four_stage_da_solve),
  TMC_CASE("regression_test_solve:RegressionRosenbrock.SixStageDASolve", six_stage_da_solve),
  TMC_CASE("regression_test_solve:RegressionRosenbrock.VectorSolve", vector_solve),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_regression_test_solve(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
