/* MICM test/integration/test_integrated_reaction_rates.cpp: two product-only
 * species, irr_1 and irr_2, accumulate the integrated rates of a -> b and
 * b -> c and must never decrease, including after b -> c is switched off. */
#include "harness.h"

TMC_TEST(can_build_chapman_system)
{
  tmc_species a = tmc_species_create("a"), b = tmc_species_create("b"), c = tmc_species_create("c");
  tmc_species irr_1 = tmc_species_create("irr_1"), irr_2 = tmc_species_create("irr_2");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_phase_add_species(&gas, &c);
  tmc_phase_add_species(&gas, &irr_1);
  tmc_phase_add_species(&gas, &irr_2);

  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = 2.15e-11;
  p.B = 0;
  p.C = 110;
  tmc_reaction r[2];
  r[0] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[0], &a);
  tmc_reaction_add_product(&r[0], &b, 1);
  tmc_reaction_add_product(&r[0], &irr_1, 1);
  tmc_reaction_set_rate_constant(&r[0], tmc_rate_arrhenius(p));
  tmc_reaction_set_phase(&r[0], &gas);
  r[1] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[1], &b);
  tmc_reaction_add_product(&r[1], &c, 1);
  tmc_reaction_add_product(&r[1], &irr_2, 1);
  tmc_reaction_set_rate_constant(&r[1], tmc_rate_user_defined("r2", 1.0));
  tmc_reaction_set_phase(&r[1], &gas);

  tmc_system system = tmc_system_create(&gas);
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, r, 2);
  tmc_builder_set_width(builder, 1);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  tmc_builder_free(builder);
  tmc_reaction_free(&r[0]);
  tmc_reaction_free(&r[1]);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  if (!solver)
    tmc_test_fail(__FILE__, __LINE__, "build failed: %s", err.message);
  TMC_ASSERT_TRUE(solver != NULL);

  tmc_state *state = tmc_solver_create_state(solver, 1);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "r2", 1.0, &err), TMC_OK);
  const int a_idx = tmc_state_variable_index(state, "a");
  const int irr1_idx = tmc_state_variable_index(state, "irr_1"), irr2_idx = tmc_state_variable_index(state, "irr_2");
  TMC_ASSERT_TRUE(a_idx >= 0 && irr1_idx >= 0 && irr2_idx >= 0);
  TMC_EXPECT_TRUE(tmc_state_variable_index(state, "b") >= 0 && tmc_state_variable_index(state, "c") >= 0);
  tmc_dense_fill(tmc_state_variables(state), 0.0);
  tmc_state_set(state, 0, a_idx, 1.0);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 273;
  conditions[0].pressure = 1000;

  for (int t = 0; t < 100; ++t)
  {
    if (t > 50)
      TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "r2", 0.0, &err), TMC_OK);
    tmc_solver_update_state_parameters(solver, state);
    const double irr1 = tmc_state_get(state, 0, irr1_idx), irr2 = tmc_state_get(state, 0, irr2_idx);
    tmc_result result = tmc_solver_solve(solver, 30.0, state);
    (void)result;
    if (!(tmc_state_get(state, 0, irr1_idx) >= irr1))
      tmc_test_fail(__FILE__, __LINE__, "step %d: irr_1 decreased from %.17g to %.17g", t, irr1, tmc_state_get(state, 0, irr1_idx));
    if (!(tmc_state_get(state, 0, irr2_idx) >= irr2))
      tmc_test_fail(__FILE__, __LINE__, "step %d: irr_2 decreased from %.17g to %.17g", t, irr2, tmc_state_get(state, 0, irr2_idx));
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("integrated_reaction_rates:ChapmanIntegration.CanBuildChapmanSystem", can_build_chapman_system),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_integrated_reaction_rates(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
