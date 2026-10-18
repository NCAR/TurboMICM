/* MICM test/integration/test_chapman_integration.cpp with chapman_policy.hpp:
 * the Chapman ozone mechanism integrated for 100 steps of 30 s, checking the
 * invariants it must obey: no NaN or negative value, oxygen atoms conserved,
 * and species with no net reaction unchanged. */
#include "harness.h"
#include "tmc_dense_ext.h"

static tmc_rate_constant arrhenius(double A, double B, double C)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  p.B = B;
  p.C = C;
  return tmc_rate_arrhenius(p);
}

/* Oxygen-atom count, weighting each species by the oxygen atoms it carries */
static double oxygen_atoms(const tmc_state *state)
{
  static const char *const names[4] = { "O", "O1D", "O2", "O3" };
  static const double weights[4] = { 1.0, 1.0, 2.0, 3.0 };
  double total = 0.0;
  for (int i = 0; i < 4; ++i)
    total += weights[i] * tmc_state_get(state, 0, tmc_state_variable_index(state, names[i]));
  return total;
}

/* |actual - reference| <= rel |reference| (MICM's EXPECT_REAL_REL) */
#define EXPECT_REL_TO(actual, reference, rel) TMC_EXPECT_NEAR((actual), (reference), (rel) * fabs(reference))

TMC_TEST(can_build_chapman_system)
{
  static const char *const names[9] = { "O", "O1D", "O2", "O3", "M", "Ar", "N2", "H2O", "CO2" };
  tmc_species s[9];
  tmc_phase gas = tmc_phase_create("gas");
  for (int i = 0; i < 9; ++i)
  {
    s[i] = tmc_species_create(names[i]);
    tmc_phase_add_species(&gas, &s[i]);
  }
  enum
  {
    O,
    O1D,
    O2,
    O3,
    M,
    Ar,
    N2,
    H2O,
    CO2
  };
  tmc_reaction r[7];
  for (int i = 0; i < 7; ++i)
  {
    r[i] = tmc_reaction_create();
    tmc_reaction_set_phase(&r[i], &gas);
  }
  /* O1D + N2 -> O + N2 */
  tmc_reaction_add_reactant(&r[0], &s[O1D]);
  tmc_reaction_add_reactant(&r[0], &s[N2]);
  tmc_reaction_add_product(&r[0], &s[O], 1);
  tmc_reaction_add_product(&r[0], &s[N2], 1);
  tmc_reaction_set_rate_constant(&r[0], arrhenius(2.15e-11, 0, 110));
  /* O1D + O2 -> O + O2 */
  tmc_reaction_add_reactant(&r[1], &s[O1D]);
  tmc_reaction_add_reactant(&r[1], &s[O2]);
  tmc_reaction_add_product(&r[1], &s[O], 1);
  tmc_reaction_add_product(&r[1], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[1], arrhenius(3.3e-11, 0, 55));
  /* O + O3 -> 2 O2 */
  tmc_reaction_add_reactant(&r[2], &s[O]);
  tmc_reaction_add_reactant(&r[2], &s[O3]);
  tmc_reaction_add_product(&r[2], &s[O2], 2);
  tmc_reaction_set_rate_constant(&r[2], arrhenius(8e-12, 0, -2060));
  /* O + O2 + M -> O3 + M */
  tmc_reaction_add_reactant(&r[3], &s[O]);
  tmc_reaction_add_reactant(&r[3], &s[O2]);
  tmc_reaction_add_reactant(&r[3], &s[M]);
  tmc_reaction_add_product(&r[3], &s[O3], 1);
  tmc_reaction_add_product(&r[3], &s[M], 1);
  tmc_reaction_set_rate_constant(&r[3], arrhenius(6.0e-34, 0, 2.4));
  /* O2 -> 2 O */
  tmc_reaction_add_reactant(&r[4], &s[O2]);
  tmc_reaction_add_product(&r[4], &s[O], 2);
  tmc_reaction_set_rate_constant(&r[4], tmc_rate_user_defined("jO2", 1.0));
  /* O3 -> O1D + O2 */
  tmc_reaction_add_reactant(&r[5], &s[O3]);
  tmc_reaction_add_product(&r[5], &s[O1D], 1);
  tmc_reaction_add_product(&r[5], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[5], tmc_rate_user_defined("jO3a", 1.0));
  /* O3 -> O + O2 */
  tmc_reaction_add_reactant(&r[6], &s[O3]);
  tmc_reaction_add_product(&r[6], &s[O], 1);
  tmc_reaction_add_product(&r[6], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[6], tmc_rate_user_defined("jO3b", 1.0));

  tmc_system system = tmc_system_create(&gas);
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_three_stage());
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, r, 7);
  tmc_builder_set_ignore_unused_species(builder, true);
  tmc_builder_set_width(builder, 1);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  tmc_builder_free(builder);
  for (int i = 0; i < 7; ++i)
    tmc_reaction_free(&r[i]);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  if (!solver)
    tmc_test_fail(__FILE__, __LINE__, "build failed: %s", err.message);
  TMC_ASSERT_TRUE(solver != NULL);

  tmc_state *state = tmc_solver_create_state(solver, 1);
  /* MICM assigns the 9 values to the row by position */
  static const double concentrations[9] = { 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
  TMC_EXPECT_EQ(tmc_dense_assign_row(tmc_state_variables(state), 0, concentrations, 9, &err), TMC_OK);
  static const double photo_rates[3] = { 0.1, 0.2, 0.3 };
  TMC_EXPECT_EQ(tmc_state_num_variables(state), 9); /* ignored unused species stay in the state, as in MICM */
  TMC_EXPECT_EQ(tmc_solver_num_custom_parameters(solver), 3);
  TMC_EXPECT_EQ(tmc_dense_assign_row(tmc_state_custom_rate_parameters(state), 0, photo_rates, 3, &err), TMC_OK);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 2;
  conditions[0].pressure = 3;

  const double oxygen_initial = oxygen_atoms(state);
  const double o_initial = tmc_state_get(state, 0, tmc_state_variable_index(state, "O"));
  const double o2_initial = tmc_state_get(state, 0, tmc_state_variable_index(state, "O2"));

  /* Species that no reaction can change, among those the solver kept in the state */
  static const char *const spectator_names[5] = { "M", "Ar", "N2", "H2O", "CO2" };
  int spectator_index[5], n_spectators = 0;
  double spectator_value[5];
  for (int i = 0; i < 5; ++i)
  {
    const int index = tmc_state_variable_index(state, spectator_names[i]);
    if (index >= 0)
    {
      spectator_index[n_spectators] = index;
      spectator_value[n_spectators++] = tmc_state_get(state, 0, index);
    }
  }
  TMC_EXPECT_EQ(n_spectators, 5);

  const int n_variables = tmc_state_num_variables(state);
  for (int t = 0; t < 100; ++t)
  {
    TMC_EXPECT_EQ(tmc_dense_assign_row(tmc_state_custom_rate_parameters(state), 0, photo_rates, 3, &err), TMC_OK);
    tmc_solver_update_state_parameters(solver, state);
    tmc_result result = tmc_solver_solve(solver, 30.0, state);
    if (result.state != TMC_CONVERGED)
    {
      tmc_test_fail(__FILE__, __LINE__, "step %d did not converge: %s", t, tmc_solver_state_name(result.state));
      break;
    }
    bool finite = true;
    for (int i = 0; i < n_variables; ++i)
    {
      const double value = tmc_state_get(state, 0, i);
      if (!isfinite(value))
      {
        tmc_test_fail(__FILE__, __LINE__, "%s not finite at step %d", tmc_state_variable_name(state, i), t);
        finite = false;
      }
      if (!(value >= 0.0))
        tmc_test_fail(__FILE__, __LINE__, "%s went negative at step %d", tmc_state_variable_name(state, i), t);
    }
    if (!finite)
      break;
    EXPECT_REL_TO(oxygen_atoms(state), oxygen_initial, 1.0e-8);
    for (int i = 0; i < n_spectators; ++i)
      EXPECT_REL_TO(tmc_state_get(state, 0, spectator_index[i]), spectator_value[i], 1.0e-8);
  }

  /* O2 is destroyed by jO2 = 0.1 /s over 100 steps of 30 s; essentially all of
   * the oxygen ends up as atomic O. */
  const double o_final = tmc_state_get(state, 0, tmc_state_variable_index(state, "O"));
  const double o2_final = tmc_state_get(state, 0, tmc_state_variable_index(state, "O2"));
  TMC_EXPECT_TRUE(o2_final < 0.5 * o2_initial);
  TMC_EXPECT_TRUE(o_final > o_initial);
  TMC_EXPECT_TRUE(o_final > 0.5 * oxygen_initial);

  tmc_state_free(state);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("chapman_integration:ChapmanIntegration.CanBuildChapmanSystem", can_build_chapman_system),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_chapman_integration(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
