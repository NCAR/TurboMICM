/* MICM test/unit/solver/test_solver_lambda_rate_constant.cpp
 *
 * MICM's Solver::GetLambdaRateConstantByName returns the solver's own lambda
 * parameters for a label, so a caller can replace the function after the
 * build; here that is tmc_support_solver_lambda_rate_constant. */
#include "harness.h"
#include "tmc_support.h"

static double rate_1e3(const tmc_conditions *conditions, void *context)
{
  (void)context;
  return 1.0e-3 * conditions->temperature;
}

static double rate_2e3(const tmc_conditions *conditions, void *context)
{
  (void)context;
  return 2.0e-3 * conditions->temperature;
}

TMC_TEST(get_lambda_rate_constant_by_name_can_override_lambda)
{
  tmc_species a = tmc_species_create("A"), b = tmc_species_create("B");
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &a);
  tmc_phase_add_species(&gas_phase, &b);
  tmc_reaction reaction = tmc_reaction_create();
  tmc_reaction_add_reactant(&reaction, &a);
  tmc_reaction_add_product(&reaction, &b, 1);
  tmc_reaction_set_rate_constant(&reaction, tmc_rate_lambda("lambda_rc", rate_1e3, NULL));
  tmc_reaction_set_phase(&reaction, &gas_phase);
  tmc_system system = tmc_system_create(&gas_phase);

  tmc_builder *builder = tmc_builder_create_backward_euler(tmc_backward_euler_defaults());
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, &reaction, 1);
  tmc_error err;
  tmc_solver *solver = tmc_builder_build(builder, &err);
  tmc_builder_free(builder);
  TMC_ASSERT_TRUE(solver != NULL);

  tmc_lambda_parameters *lambda_ref = tmc_support_solver_lambda_rate_constant(solver, "lambda_rc", &err);
  TMC_ASSERT_TRUE(lambda_ref != NULL);
  TMC_EXPECT_STREQ(lambda_ref->label, "lambda_rc");
  lambda_ref->function = rate_2e3;

  /* an unknown label is MICM's RATE_CONSTANT_NOT_FOUND */
  TMC_EXPECT_TRUE(tmc_support_solver_lambda_rate_constant(solver, "no such lambda", &err) == NULL);
  TMC_EXPECT_EQ(err.category, TMC_ERROR_SOLVER);
  TMC_EXPECT_EQ(err.code, TMC_SOLVER_RATE_CONSTANT_NOT_FOUND);

  tmc_state *state = tmc_solver_create_state(solver, 1);
  const double tolerances[2] = { 1.0e-6, 1.0e-6 };
  TMC_EXPECT_EQ(tmc_state_set_absolute_tolerances(state, tolerances, 2, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_concentration(state, "A", 1.0, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_concentration(state, "B", 0.0, &err), TMC_OK);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 300.0;
  conditions[0].pressure = 101253.3;
  conditions[0].air_density = 1.0e6;

  tmc_solver_update_state_parameters(solver, state);
  TMC_EXPECT_NEAR(tmc_dense_get(tmc_state_rate_constants(state), 0, 0), 0.6, 1.0e-12);

  tmc_result result = tmc_solver_solve(solver, 1.0, state);
  TMC_EXPECT_EQ(result.stats.final_time, 1.0);

  tmc_state_free(state);
  tmc_solver_free(solver);
  tmc_system_free(&system);
  tmc_reaction_free(&reaction);
  tmc_phase_free(&gas_phase);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("solver_lambda_rate_constant:Solver.GetLambdaRateConstantByNameCanOverrideLambda", get_lambda_rate_constant_by_name_can_override_lambda),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_solver_lambda_rate_constant(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
