/* MICM test/unit/solver/test_backward_euler.cpp. The convergence test MICM
 * calls (AbstractBackwardEuler::IsConverged) is tmc_backward_euler_is_converged
 * in tmc_solver_ext.h; MICM runs it on its standard matrix and VectorMatrix
 * widths 1 to 5. */
#include "harness.h"
#include "tmc_solver_ext.h"

static tmc_rate_constant arrhenius(double A, double B, double C)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  p.B = B;
  p.C = C;
  return tmc_rate_arrhenius(p);
}

TMC_TEST(can_call_solve)
{
  tmc_species a = tmc_species_create("A"), b = tmc_species_create("B"), c = tmc_species_create("C");
  tmc_phase gas = tmc_phase_create("gas");
  tmc_phase_add_species(&gas, &a);
  tmc_phase_add_species(&gas, &b);
  tmc_phase_add_species(&gas, &c);
  tmc_system the_system = tmc_system_create(&gas);
  tmc_reaction reactions[2];
  reactions[0] = tmc_reaction_create();
  tmc_reaction_add_reactant(&reactions[0], &a);
  tmc_reaction_add_product(&reactions[0], &b, 1);
  tmc_reaction_set_rate_constant(&reactions[0], arrhenius(2.15e-11, 0, 110));
  tmc_reaction_set_phase(&reactions[0], &gas);
  reactions[1] = tmc_reaction_create();
  tmc_reaction_add_reactant(&reactions[1], &b);
  tmc_reaction_add_product(&reactions[1], &c, 1);
  tmc_reaction_set_rate_constant(&reactions[1], arrhenius(3.3e-11, 0, 55));
  tmc_reaction_set_phase(&reactions[1], &gas);

  tmc_builder *builder = tmc_builder_create_backward_euler(tmc_backward_euler_defaults());
  tmc_builder_set_system(builder, &the_system);
  tmc_builder_set_reactions(builder, reactions, 2);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *be = tmc_builder_build(builder, &err);
  tmc_builder_free(builder);
  tmc_reaction_free(&reactions[0]);
  tmc_reaction_free(&reactions[1]);
  tmc_system_free(&the_system);
  tmc_phase_free(&gas);
  TMC_ASSERT_TRUE(be != NULL);
  TMC_EXPECT_EQ(tmc_solver_method(be), TMC_METHOD_BACKWARD_EULER);

  const double time_step = 1.0;
  tmc_state *state = tmc_solver_create_state(be, 1);
  const double atol[3] = { 1e-6, 1e-6, 1e-6 };
  TMC_EXPECT_EQ(tmc_state_set_absolute_tolerances(state, atol, 3, &err), TMC_OK);
  tmc_state_set(state, 0, 0, 1.0);
  tmc_state_set(state, 0, 1, 0.0);
  tmc_state_set(state, 0, 2, 0.0);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 272.5;
  conditions[0].pressure = 101253.3;
  conditions[0].air_density = 1e6;
  tmc_solver_update_state_parameters(be, state);

  /* MICM checks that Solve returns without throwing; the C call returns a result */
  tmc_result result = tmc_solver_solve(be, time_step, state);
  TMC_EXPECT_EQ(result.state, TMC_CONVERGED);
  TMC_EXPECT_DOUBLE_EQ(result.stats.final_time, time_step);
  tmc_state_free(state);
  tmc_solver_free(be);
  return TMC_TEST_PASS;
}

/* MICM's CheckIsConverged for one layout: residual and Yn1 of 4 cells by 3
 * variables (with L = 3 and 5 the last block is padded). */
static void check_is_converged(int L)
{
  tmc_backward_euler_parameters parameters = tmc_backward_euler_defaults();
  tmc_dense *residual = tmc_dense_create(4, 3, L, 0.0);
  tmc_dense *Yn1 = tmc_dense_create(4, 3, L, 0.0);
  parameters.small = 1e-6;
  const double relative_tolerance = 1e-3;
  double absolute_tolerance[3] = { 1e-6, 1e-6, 1e-6 };
#define IS_CONVERGED() tmc_backward_euler_is_converged(&parameters, residual, Yn1, absolute_tolerance, relative_tolerance)

  TMC_EXPECT_TRUE(IS_CONVERGED());
  tmc_dense_set(residual, 0, 1, 1e-5);
  TMC_EXPECT_FALSE(IS_CONVERGED());
  parameters.small = 1e-4;
  TMC_EXPECT_TRUE(IS_CONVERGED());
  tmc_dense_set(residual, 3, 2, 1e-3);
  TMC_EXPECT_FALSE(IS_CONVERGED());
  tmc_dense_set(Yn1, 3, 2, 10.0);
  TMC_EXPECT_TRUE(IS_CONVERGED());
  tmc_dense_set(residual, 3, 2, 1e-1);
  TMC_EXPECT_FALSE(IS_CONVERGED());
  absolute_tolerance[2] = 1.0;
  TMC_EXPECT_TRUE(IS_CONVERGED());

#undef IS_CONVERGED
  tmc_dense_free(residual);
  tmc_dense_free(Yn1);
}

TMC_TEST(is_converged)
{
  /* Matrix, then VectorMatrix<1> to <5> */
  check_is_converged(1);
  for (int L = 1; L <= 5; ++L)
    check_is_converged(L);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("backward_euler:BackwardEuler.CanCallSolve", can_call_solve),
  TMC_CASE("backward_euler:BackwardEuler.IsConverged", is_converged),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_backward_euler(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
