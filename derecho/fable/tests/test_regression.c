/* MICM test/regression/test_regression.cpp with regression_policy.hpp: a flow
 * tube with wall loss (SOA1 and SOA2 lost at fixed rates, formed from O3 +
 * a-pinene), 3600 steps of 1 s, against MICM's recorded results in
 * expected_results/{rosenbrock,backward_euler}_wall_loss.csv, copied to tests/data/regression. MICM
 * compares the first 10 rows to 1e-10 absolute. */
#include "harness.h"

enum
{
  N_STEPS = 3600,
  N_SPECIES = 4,
  N_COLUMNS = N_SPECIES + 1 /* time first */
};

static const char *const header[N_COLUMNS] = { "time", "SOA1", "SOA2", "a-pinene", "O3" };

static double air_density_mol_m3(double pressure, double temperature)
{
  const double GAS_CONSTANT = 1.380649e-23 * 6.02214076e23;
  return pressure / (GAS_CONSTANT * temperature);
}

/* regression_policy.hpp's ReadCsv: the header names and the rows of values.
 * Returns the number of data rows (all are counted, the first `max_rows` kept). */
static int read_csv(const char *path, char (*names)[32], int *n_names, double (*rows)[N_COLUMNS], int max_rows)
{
  FILE *file = fopen(path, "r");
  *n_names = 0;
  if (!file)
  {
    tmc_test_fail(__FILE__, __LINE__, "cannot open %s", path);
    return -1;
  }
  char line[1024];
  int n_rows = 0;
  if (fgets(line, sizeof line, file))
  {
    line[strcspn(line, "\r\n")] = '\0';
    for (char *item = strtok(line, ","); item && *n_names < N_COLUMNS; item = strtok(NULL, ","))
      snprintf(names[(*n_names)++], 32, "%s", item);
  }
  while (fgets(line, sizeof line, file))
  {
    if (line[strspn(line, " \t\r\n")] == '\0')
      continue;
    if (n_rows < max_rows)
    {
      char *p = line;
      for (int j = 0; j < N_COLUMNS; ++j)
      {
        char *end;
        rows[n_rows][j] = strtod(p, &end);
        if (end == p)
          tmc_test_fail(__FILE__, __LINE__, "%s: bad number in row %d", path, n_rows + 1);
        p = end;
        if (*p == ',')
          ++p;
      }
    }
    ++n_rows;
  }
  fclose(file);
  return n_rows;
}

/* TestFlowTube: the builder is freed */
static void test_flow_tube(tmc_builder *builder, const char *expected_results_path)
{
  /*
   * SOA1 ->                                , k1 = 0.01
   * SOA2 ->                                , k2 = 0.05
   * O3 + a-pinene -> 0.18 SOA1 + 0.09 SOA2 , k3 = Arrhenius(A=8.8e-17)
   */
  tmc_species soa1 = tmc_species_create("SOA1"), soa2 = tmc_species_create("SOA2");
  tmc_species apinene = tmc_species_create("a-pinene"), o3 = tmc_species_create("O3");
  const double MOLES_M3_TO_MOLECULES_CM3 = 1.0e-6 * 6.02214076e23;
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &soa1);
  tmc_phase_add_species(&gas, &soa2);
  tmc_phase_add_species(&gas, &apinene);
  tmc_phase_add_species(&gas, &o3);

  tmc_reaction r[3];
  r[0] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[0], &soa1);
  tmc_reaction_set_rate_constant(&r[0], tmc_rate_user_defined("r1", 1.0));
  tmc_reaction_set_phase(&r[0], &gas);
  r[1] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[1], &soa2);
  tmc_reaction_set_rate_constant(&r[1], tmc_rate_user_defined("r2", 1.0));
  tmc_reaction_set_phase(&r[1], &gas);
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = 8.8e-17 * MOLES_M3_TO_MOLECULES_CM3;
  r[2] = tmc_reaction_create();
  tmc_reaction_add_reactant(&r[2], &apinene);
  tmc_reaction_add_reactant(&r[2], &o3);
  tmc_reaction_add_product(&r[2], &soa1, 0.18);
  tmc_reaction_add_product(&r[2], &soa2, 0.09);
  tmc_reaction_set_rate_constant(&r[2], tmc_rate_arrhenius(p));
  tmc_reaction_set_phase(&r[2], &gas);

  tmc_system system = tmc_system_create(&gas);
  tmc_builder_set_reorder_state(builder, false);
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, r, 3);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  tmc_builder_free(builder);
  for (int i = 0; i < 3; ++i)
    tmc_reaction_free(&r[i]);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  if (!solver)
  {
    tmc_test_fail(__FILE__, __LINE__, "build failed: %s", err.message);
    return;
  }

  tmc_state *state = tmc_solver_create_state(solver, 1);
  TMC_EXPECT_EQ(tmc_state_num_variables(state), N_SPECIES);
  for (int j = 0; j < N_SPECIES && j < tmc_state_num_variables(state); ++j)
    TMC_EXPECT_STREQ(tmc_state_variable_name(state, j), header[j + 1]);
  TMC_EXPECT_EQ(tmc_state_set_concentration(state, "SOA1", 0, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_concentration(state, "SOA2", 0, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_concentration(state, "a-pinene", 8e-8, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_concentration(state, "O3", 0.00002, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "r1", 0.01, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "r2", 0.05, &err), TMC_OK);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 298.15;
  conditions[0].pressure = 101325;
  conditions[0].air_density = air_density_mol_m3(101325, 298.15);
  tmc_solver_update_state_parameters(solver, state);

  static double model[N_STEPS + 1][N_SPECIES];
  tmc_dense_get_row(tmc_state_variables(state), 0, model[0]);
  const double time_step = 1;
  for (int i_time = 0; i_time < N_STEPS; ++i_time)
  {
    double actual_solve = 0;
    while (actual_solve < time_step)
    {
      tmc_result result = tmc_solver_solve(solver, time_step - actual_solve, state);
      if (!(result.stats.final_time > 0))
      {
        tmc_test_fail(__FILE__, __LINE__, "step %d: solve advanced no time (%s)", i_time, tmc_solver_state_name(result.state));
        break;
      }
      actual_solve += result.stats.final_time;
    }
    tmc_dense_get_row(tmc_state_variables(state), 0, model[i_time + 1]);
  }

  /* Check that the results we got match those we expect */
  static double data[N_STEPS + 1][N_COLUMNS];
  char names[N_COLUMNS][32];
  int n_names = 0;
  const int n_rows = read_csv(tmc_test_data_path(expected_results_path), names, &n_names, data, N_STEPS + 1);
  TMC_EXPECT_EQ(n_names, N_COLUMNS);
  for (int j = 0; j < n_names && j < N_COLUMNS; ++j)
    TMC_EXPECT_STREQ(names[j], header[j]);
  TMC_EXPECT_EQ(n_rows, N_STEPS + 1);
  if (n_rows >= 10)
    for (int i = 0; i < 10; ++i)
      for (int j = 0; j < N_SPECIES; ++j)
        if (!(fabs(model[i][j] - data[i][j + 1]) <= 1e-10))
          tmc_test_fail(__FILE__, __LINE__, "row %d, %s: %.17g, expected %.17g", i, header[j + 1], model[i][j], data[i][j + 1]);
  tmc_state_free(state);
  tmc_solver_free(solver);
}

TMC_TEST(rosenbrock_flow_tube)
{
  tmc_builder *builder = tmc_builder_create_rosenbrock(tmc_rosenbrock_two_stage());
  tmc_builder_set_width(builder, 1);
  test_flow_tube(builder, "regression/rosenbrock_wall_loss.csv");
  return TMC_TEST_PASS;
}

TMC_TEST(backward_euler_flow_tube)
{
  tmc_builder *builder = tmc_builder_create_backward_euler(tmc_backward_euler_defaults());
  tmc_builder_set_width(builder, 1);
  test_flow_tube(builder, "regression/backward_euler_wall_loss.csv");
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("regression:RosenbrockRegression.FlowTube", rosenbrock_flow_tube),
  TMC_CASE("regression:BackwardEulerRegression.FlowTube", backward_euler_flow_tube),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_regression(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
