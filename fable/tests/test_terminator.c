/* MICM test/integration/test_terminator.cpp with terminator.hpp: the
 * "Terminator" toy mechanism, Cl2 -> 2 Cl by a user-defined photolysis rate
 * that follows the sun and Cl + Cl -> Cl2, against its analytical solution.
 * MICM's CpuSolverBuilder is width L = 1; its VectorBuilder<L> is L. */
#include "harness.h"
#include "libcxx_random.h"

#define PI 3.14159265358979323846

static void test_terminator(int L, int number_of_grid_cells)
{
  tmc_species cl2 = tmc_species_create_with_property("Cl2", TMC_PROPERTY_ABSOLUTE_TOLERANCE, 1.0e-20);
  tmc_species cl = tmc_species_create_with_property("Cl", TMC_PROPERTY_ABSOLUTE_TOLERANCE, 1.0e-20);
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &cl2);
  tmc_phase_add_species(&gas, &cl);

  tmc_reaction toy_r1 = tmc_reaction_create();
  tmc_reaction_add_reactant(&toy_r1, &cl2);
  tmc_reaction_add_product(&toy_r1, &cl, 2.0);
  tmc_reaction_set_phase(&toy_r1, &gas);
  tmc_reaction_set_rate_constant(&toy_r1, tmc_rate_user_defined("toy_k1", 1.0));

  const double k2 = 1.0;
  tmc_arrhenius_parameters k2_parameters = tmc_arrhenius_defaults();
  k2_parameters.A = k2;
  tmc_reaction toy_r2 = tmc_reaction_create();
  tmc_reaction_add_reactant(&toy_r2, &cl);
  tmc_reaction_add_reactant(&toy_r2, &cl);
  tmc_reaction_add_product(&toy_r2, &cl2, 1.0);
  tmc_reaction_set_phase(&toy_r2, &gas);
  tmc_reaction_set_rate_constant(&toy_r2, tmc_rate_arrhenius(k2_parameters));
  tmc_reaction reactions[2] = { toy_r1, toy_r2 };

  tmc_rosenbrock_parameters parameters = tmc_rosenbrock_three_stage();
  parameters.max_number_of_steps = 100000;
  tmc_builder *builder = tmc_builder_create_rosenbrock(parameters);
  tmc_builder_set_ignore_unused_species(builder, true);
  tmc_builder_set_width(builder, L);
  tmc_system system = tmc_system_create(&gas);
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, reactions, 2);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  tmc_builder_free(builder);
  tmc_reaction_free(&toy_r1);
  tmc_reaction_free(&toy_r2);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  if (!solver)
  {
    tmc_test_fail(__FILE__, __LINE__, "L = %d: build failed: %s", L, err.message);
    return;
  }

  tmc_state *state = tmc_solver_create_state(solver, number_of_grid_cells);
  tmc_state_set_relative_tolerance(state, 1.0e-8);

  /* std::bind(std::lognormal_distribution(-2.0, 2.0), std::default_random_engine()) */
  tmc_lognormal_stream get_double = tmc_lognormal_stream_create(-2.0, 2.0);
  double cl2_initial[4], cl_initial[4];
  for (int i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
  {
    cl2_initial[i_cell] = tmc_lognormal_stream_next(&get_double) * 1.0e-6;
    cl_initial[i_cell] = tmc_lognormal_stream_next(&get_double) * 1.0e-10;
  }
  TMC_EXPECT_EQ(tmc_state_set_concentrations(state, "Cl2", cl2_initial, number_of_grid_cells, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_concentrations(state, "Cl", cl_initial, number_of_grid_cells, &err), TMC_OK);

  tmc_dense *orig_state_vars = tmc_dense_copy(tmc_state_variables(state));
  const int cl_index = tmc_state_variable_index(state, "Cl"), cl2_index = tmc_state_variable_index(state, "Cl2");
  const int steps = (int)floor(2.0 * PI / 0.3);
  tmc_conditions *conditions = tmc_state_conditions(state);
  for (int step = 0; step < steps; ++step)
  {
    const double lon = step * 0.3;
    tmc_dense_copy_into(tmc_state_variables(state), orig_state_vars);
    double k1_values[4];
    for (int i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    {
      const double k1_lat_center = PI * 20.0 / 180.0;
      const double k1_lon_center = PI * 300.0 / 180.0;
      const double lat = PI / 180.0 * (i_cell * (90.0 / number_of_grid_cells));
      const double k1 = fmax(0.0, sin(lat) * sin(k1_lat_center) + cos(lat) * cos(k1_lat_center) * cos(lon - k1_lon_center));
      k1_values[i_cell] = k1;
      conditions[i_cell].temperature = 298.0;  /* K */
      conditions[i_cell].pressure = 101300.0;  /* Pa */
      conditions[i_cell].air_density = 42.0;   /* mol m-3 */
    }
    TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameters(state, "toy_k1", k1_values, number_of_grid_cells, &err), TMC_OK);
    tmc_solver_update_state_parameters(solver, state);

    const double dt = 30.0;
    tmc_result result = tmc_solver_solve(solver, dt, state);
    if (result.state != TMC_CONVERGED)
      tmc_test_fail(__FILE__, __LINE__, "L = %d, %d cells, step %d: solver state %s", L, number_of_grid_cells, step, tmc_solver_state_name(result.state));

    for (int i_cell = 0; i_cell < number_of_grid_cells; ++i_cell)
    {
      const double r = k1_values[i_cell] / (4.0 * k2);
      const double cl_i = cl_initial[i_cell];
      const double cl2_i = cl2_initial[i_cell];
      const double cly = cl_i + 2.0 * cl2_i;
      const double det = sqrt(r * r + 2.0 * r * cly);
      const double e = exp(-4.0 * k2 * det * dt);
      const double l = (det * k2 * dt) > 1.0e-16 ? (1.0 - e) / det / dt : 4.0 * k2;
      const double cl_f = -l * (cl_i - det + r) * (cl_i + det + r) / (1.0 + e + dt * l * (cl_i + r));
      const double cl2_f = -cl_f / 2.0;
      const double cl_ref = cl_i + dt * cl_f;
      const double cl2_ref = cl2_i + dt * cl2_f;
      const double cl_value = tmc_state_get(state, i_cell, cl_index), cl2_value = tmc_state_get(state, i_cell, cl2_index);
      if (!(fabs(cl_value - cl_ref) <= fabs(cl_ref) * 1.0e-6 + 1.0e-14))
        tmc_test_fail(__FILE__, __LINE__, "L = %d, %d cells, step %d, cell %d: Cl %.17g vs %.17g", L, number_of_grid_cells, step, i_cell, cl_value, cl_ref);
      if (!(fabs(cl2_value - cl2_ref) <= fabs(cl2_ref) * 1.0e-6 + 1.0e-14))
        tmc_test_fail(__FILE__, __LINE__, "L = %d, %d cells, step %d, cell %d: Cl2 %.17g vs %.17g", L, number_of_grid_cells, step, i_cell, cl2_value, cl2_ref);
    }
  }
  tmc_dense_free(orig_state_vars);
  tmc_state_free(state);
  tmc_solver_free(solver);
}

TMC_TEST(terminator)
{
  static const int cells[8] = { 1, 2, 3, 4, 1, 2, 3, 4 };
  for (int i = 0; i < 8; ++i)
    test_terminator(1, cells[i]);
  return TMC_TEST_PASS;
}

TMC_TEST(vector_terminator)
{
  for (int L = 1; L <= 4; ++L)
    test_terminator(L, L);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("terminator:RosenbrockSolver.Terminator", terminator),
  TMC_CASE("terminator:RosenbrockSolver.VectorTerminator", vector_terminator),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_terminator(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
