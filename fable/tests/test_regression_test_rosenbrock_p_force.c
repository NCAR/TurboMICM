/* MICM test/regression/RosenbrockChapman/regression_test_p_force.cpp with
 * regression_test_p_force_policy.hpp: the general solver's rate constants and
 * forcing for the Chapman mechanism against the hard-coded solver's. MICM's
 * CpuSolverBuilder is width L = 1; its VectorBuilder<L> is L. */
#include "chapman_fixed.h"
#include "libcxx_random.h"

static tmc_solver *chapman_solver(int L)
{
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_width(builder, L);
  return cf_chapman_solver(builder);
}

static const double photo_rates[3][3] = { { 1.0e-4, 1.0e-5, 1.0e-6 }, { 3.2e-4, 7.3e-5, 3.2e-6 }, { 5.2e-4, 8.2e-5, 4.6e-6 } };
static const double temperatures[3] = { 284.19, 215.02, 299.31 }; /* [K] */
static const double pressures[3] = { 101245.0, 100789.2, 101398.0 }; /* [Pa] */

/* TestRateConstants: every rate constant of three cells equals the hard-coded solver's */
static void test_rate_constants(tmc_solver *solver)
{
  if (!solver)
    return;
  tmc_state *state = tmc_solver_create_state(solver, 3);
  /* the photolysis rates are assigned by position, in the state's parameter order */
  static const char *const labels[3] = { "jO2", "jO3a", "jO3b" };
  TMC_EXPECT_EQ(tmc_solver_num_custom_parameters(solver), 3);
  for (int j = 0; j < 3; ++j)
    TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(solver, j), labels[j]);
  tmc_conditions *conditions = tmc_state_conditions(state);
  for (int i = 0; i < 3; ++i)
  {
    tmc_dense_set_row(tmc_state_custom_rate_parameters(state), i, photo_rates[i]);
    conditions[i].temperature = temperatures[i];
    conditions[i].pressure = pressures[i];
  }
  tmc_solver_update_state_parameters(solver, state);

  const tmc_dense *rate_constants = tmc_state_rate_constants(state);
  TMC_EXPECT_EQ(rate_constants->n_cols, CF_N_RATE_CONSTANTS);
  for (int i = 0; i < 3; ++i)
  {
    double fixed_rate_constants[CF_N_RATE_CONSTANTS];
    cf_update_state(conditions[i].temperature, conditions[i].pressure, photo_rates[i], fixed_rate_constants);
    for (int j = 0; j < CF_N_RATE_CONSTANTS && j < rate_constants->n_cols; ++j)
      if (!(tmc_dense_get(rate_constants, i, j) == fixed_rate_constants[j]))
        tmc_test_fail(__FILE__, __LINE__, "L = %d, cell %d, rate constant %d: %.17g vs %.17g", tmc_state_width(state), i, j, tmc_dense_get(rate_constants, i, j), fixed_rate_constants[j]);
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
}

/* TestForcing: the forcing of three cells with lognormal concentrations and
 * rate constants (std::mt19937 seeded 739, as MICM without MICM_TEST_SEED). */
static void test_forcing(tmc_solver *solver)
{
  if (!solver)
    return;
  tmc_mt19937 engine = tmc_mt19937_create(739);
  tmc_lognormal dist = tmc_lognormal_create(-2.0, 2.0);
  tmc_state *state = tmc_solver_create_state(solver, 3);
  tmc_dense *variables = tmc_state_variables(state);
  tmc_dense *rate_constants = tmc_state_rate_constants(state);
  for (size_t k = 0; k < variables->size; ++k)
    variables->data[k] = tmc_mt19937_lognormal(&dist, &engine);
  for (size_t k = 0; k < rate_constants->size; ++k)
    rate_constants->data[k] = tmc_mt19937_lognormal(&dist, &engine);

  tmc_dense *forcing = tmc_dense_create(3, CF_N_SPECIES, tmc_state_width(state), 0.0);
  tmc_solver_add_forcing_terms(solver, state, variables, forcing);
  TMC_EXPECT_EQ(forcing->n_cols, CF_N_SPECIES);
  for (int i = 0; i < 3; ++i)
  {
    double cell_rate_constants[CF_N_RATE_CONSTANTS], cell_variables[CF_N_SPECIES], fixed_forcing[CF_N_SPECIES];
    tmc_dense_get_row(rate_constants, i, cell_rate_constants);
    for (int j = 0; j < CF_N_SPECIES; ++j)
      cell_variables[j] = tmc_state_get(state, i, tmc_state_variable_index(state, cf_species_names[j]));
    cf_force(cell_rate_constants, cell_variables, fixed_forcing);
    for (int j = 0; j < CF_N_SPECIES; ++j)
    {
      const double a = tmc_dense_get(forcing, i, tmc_state_variable_index(state, cf_species_names[j]));
      const double b = fixed_forcing[j];
      if (!(fabs(a - b) <= (fabs(a) + fabs(b)) * 1.0e-8 + 1.0e-12))
        tmc_test_fail(__FILE__, __LINE__, "L = %d, cell %d, species %s: forcing %.17g vs %.17g", tmc_state_width(state), i, cf_species_names[j], a, b);
    }
  }
  tmc_dense_free(forcing);
  tmc_state_free(state);
  tmc_solver_free(solver);
}

TMC_TEST(rate_constants)
{
  test_rate_constants(chapman_solver(1));
  return TMC_TEST_PASS;
}

TMC_TEST(vector_rate_constants)
{
  for (int L = 1; L <= 4; ++L)
    test_rate_constants(chapman_solver(L));
  return TMC_TEST_PASS;
}

TMC_TEST(forcing)
{
  test_forcing(chapman_solver(1));
  return TMC_TEST_PASS;
}

TMC_TEST(vector_forcing)
{
  for (int L = 1; L <= 4; ++L)
    test_forcing(chapman_solver(L));
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("regression_test_rosenbrock_p_force:RegressionRosenbrock.RateConstants", rate_constants),
  TMC_CASE("regression_test_rosenbrock_p_force:RegressionRosenbrock.VectorRateConstants", vector_rate_constants),
  TMC_CASE("regression_test_rosenbrock_p_force:RegressionRosenbrock.Forcing", forcing),
  TMC_CASE("regression_test_rosenbrock_p_force:RegressionRosenbrock.VectorForcing", vector_forcing),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_regression_test_rosenbrock_p_force(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
