/* MICM test/integration/test_equilibrium_constraint.cpp with
 * equilibrium_constraint_policy.hpp.
 *
 * The builder is CpuSolverBuilder with VectorMatrix<MICM_DEFAULT_VECTOR_SIZE>
 * and SparseMatrixVectorOrdering<MICM_DEFAULT_VECTOR_SIZE>; MICM's default
 * vector size is 4, so the solvers here have vector width 4 and the builder's
 * Doolittle LU (cst_build). state.constraint_size_ and
 * state.upper_left_identity_diagonal_ are the solver's
 * tmc_solver_num_constraints and tmc_solver_mass_matrix_diagonal. */
#include "constraint_support.h"

#include <float.h>

#define WIDTH 4 /* MICM_DEFAULT_VECTOR_SIZE */

static tmc_species species_named(const char *name)
{
  return tmc_species_create(name);
}

/* The equilibrium constraint `name`: reactant <-> product (unit stoichiometry), product algebraic. */
static tmc_constraint pair_equilibrium(const char *name, const tmc_species *reactant, const tmc_species *product, double K_HLC_ref, double delta_H)
{
  return cst_equilibrium(name, product, (const tmc_stoich_species[]){ { *reactant, 1.0 } }, 1, (const tmc_stoich_species[]){ { *product, 1.0 } }, 1, K_HLC_ref, delta_H);
}

static int index_of(tmc_state *st, const char *name)
{
  return tmc_state_variable_index(st, name);
}

static double parameter(tmc_state *st, int index)
{
  return tmc_dense_get(tmc_state_custom_rate_parameters(st), 0, index);
}

static void set_conditions(tmc_state *st, double temperature)
{
  tmc_state_conditions(st)[0].temperature = temperature;
  tmc_state_conditions(st)[0].pressure = 101325.0;
}

/* The SolverBuilder accepts constraints, and the solver can be built. */
TMC_TEST(set_constraints_api_works)
{
  const tmc_species A = species_named("A"), B = species_named("B"), C = species_named("C");
  tmc_system system = cst_system((const tmc_species[]){ A, B, C }, 3);
  tmc_reaction r = cst_reaction(&A, 1, (const tmc_stoich_species[]){ { B, 1.0 } }, 1, cst_arrhenius(0.1));
  /* C is an algebraic variable, in no kinetic reaction */
  const double K_eq = 0.034;
  tmc_constraint c = pair_equilibrium("B_C_eq", &B, &C, K_eq, -2400.0);
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(tmc_rosenbrock_four_stage_dae(), &system, &r, 1, &c, 1, WIDTH, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  /* constraint metadata */
  CST_REQUIRE(tmc_state_num_variables(st) == 3);
  CST_REQUIRE(tmc_solver_num_constraints(s) == 1);
  CST_REQUIRE(index_of(st, "A") >= 0);
  CST_REQUIRE(index_of(st, "B") >= 0);
  CST_REQUIRE(index_of(st, "C") >= 0);
  CST_REQUIRE(tmc_state_custom_parameter_index(st, "B_C_eq") >= 0);
  /* the mass-matrix diagonal: species only, algebraic rows 0 */
  const double *M = tmc_solver_mass_matrix_diagonal(s);
  CST_REQUIRE(tmc_solver_num_variables(s) == 3);
  TMC_EXPECT_EQ(M[index_of(st, "A")], 1.0);
  TMC_EXPECT_EQ(M[index_of(st, "B")], 1.0);
  TMC_EXPECT_EQ(M[index_of(st, "C")], 0.0);
  const int B_C_eq = tmc_state_custom_parameter_index(st, "B_C_eq");
  tmc_state_set(st, 0, index_of(st, "A"), 1.0);
  tmc_state_set(st, 0, index_of(st, "B"), 0.1);
  tmc_state_set(st, 0, index_of(st, "C"), K_eq * 0.1); /* C = K_eq B */
  set_conditions(st, 300.0);
  tmc_solver_update_state_parameters(s, st);
  const double expected_K_eq = cst_equilibrium_constant(K_eq, -2400.0, 300.0);
  TMC_EXPECT_NEAR(parameter(st, B_C_eq), expected_K_eq, fabs(expected_K_eq) * 1e2 * DBL_EPSILON);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  return TMC_TEST_PASS;
}

/* Several constraints through the same API. */
TMC_TEST(set_constraints_api_multiple_constraints)
{
  const tmc_species A = species_named("A"), B = species_named("B"), C = species_named("C"), D = species_named("D"), E = species_named("E"), F = species_named("F");
  tmc_system system = cst_system((const tmc_species[]){ A, B, C, D, E, F }, 6);
  tmc_reaction r[2] = {
    cst_reaction(&A, 1, (const tmc_stoich_species[]){ { B, 1.0 } }, 1, cst_arrhenius(0.5)),
    cst_reaction(&D, 1, (const tmc_stoich_species[]){ { E, 1.0 } }, 1, cst_arrhenius(0.2)),
  };
  const double K_eq1 = 0.034, K_eq2 = 0.012, delta_H1 = -2400.0, delta_H2 = -2000.0;
  tmc_constraint c[2] = { pair_equilibrium("B_C_eq", &B, &C, K_eq1, delta_H1), pair_equilibrium("E_F_eq", &E, &F, K_eq2, delta_H2) };
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(tmc_rosenbrock_three_stage(), &system, r, 2, c, 2, WIDTH, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  CST_REQUIRE(tmc_state_num_variables(st) == 6);
  CST_REQUIRE(tmc_solver_num_constraints(s) == 2);
  const double *M = tmc_solver_mass_matrix_diagonal(s);
  const char *const names[] = { "A", "B", "C", "D", "E", "F" };
  const double expected_M[] = { 1.0, 1.0, 0.0, 1.0, 1.0, 0.0 };
  for (int i = 0; i < 6; ++i)
  {
    CST_REQUIRE(index_of(st, names[i]) >= 0);
    TMC_EXPECT_EQ(M[index_of(st, names[i])], expected_M[i]);
  }
  tmc_state_set(st, 0, index_of(st, "A"), 1.0);
  tmc_state_set(st, 0, index_of(st, "B"), 0.1);
  tmc_state_set(st, 0, index_of(st, "C"), K_eq1 * 0.1);
  tmc_state_set(st, 0, index_of(st, "D"), 0.5);
  tmc_state_set(st, 0, index_of(st, "E"), 0.05);
  tmc_state_set(st, 0, index_of(st, "F"), K_eq2 * 0.05);
  const double current_temp = 310.0;
  set_conditions(st, current_temp);
  tmc_solver_update_state_parameters(s, st);
  const double expected_K_eq1 = cst_equilibrium_constant(K_eq1, delta_H1, current_temp);
  const double expected_K_eq2 = cst_equilibrium_constant(K_eq2, delta_H2, current_temp);
  const int k1 = tmc_state_custom_parameter_index(st, "B_C_eq"), k2 = tmc_state_custom_parameter_index(st, "E_F_eq");
  CST_REQUIRE(k1 >= 0 && k2 >= 0);
  TMC_EXPECT_NEAR(parameter(st, k1), expected_K_eq1, 1e-10);
  TMC_EXPECT_NEAR(parameter(st, k2), expected_K_eq2, fabs(expected_K_eq2) * 1e2 * DBL_EPSILON);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c[0]);
  tmc_constraint_free(&c[1]);
  tmc_reaction_free(&r[0]);
  tmc_reaction_free(&r[1]);
  tmc_system_free(&system);
  return TMC_TEST_PASS;
}

/* A -> B (k = 1) with K_eq B - C = 0: integrate 0.1 s in steps of 0.001 s. */
static void solve_a_to_b_with_constraint(bool reorder_state, double temperature, bool exact_K_tolerance)
{
  const tmc_species A = species_named("A"), B = species_named("B"), C = species_named("C");
  tmc_system system = cst_system((const tmc_species[]){ A, B, C }, 3);
  tmc_reaction r = cst_reaction(&A, 1, (const tmc_stoich_species[]){ { B, 1.0 } }, 1, cst_arrhenius(1.0));
  const double K_eq = 2.0, delta_H = -2400.0;
  tmc_constraint c = pair_equilibrium("B_C_eq", &B, &C, K_eq, delta_H);
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(tmc_rosenbrock_three_stage(), &system, &r, 1, &c, 1, WIDTH, reorder_state, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  const int A_idx = index_of(st, "A"), B_idx = index_of(st, "B"), C_idx = index_of(st, "C");
  const int B_C_eq = tmc_state_custom_parameter_index(st, "B_C_eq");
  CST_REQUIRE(A_idx >= 0 && B_idx >= 0 && C_idx >= 0 && B_C_eq >= 0);
  tmc_state_set(st, 0, A_idx, 1.0);
  tmc_state_set(st, 0, B_idx, 0.0);
  tmc_state_set(st, 0, C_idx, 0.0);
  set_conditions(st, temperature);
  tmc_solver_update_state_parameters(s, st);
  const double expected_K_eq = cst_equilibrium_constant(K_eq, delta_H, temperature);
  TMC_EXPECT_NEAR(parameter(st, B_C_eq), expected_K_eq, exact_K_tolerance ? fabs(expected_K_eq) * 1e2 * DBL_EPSILON : 1e-10);
  const double dt = 0.001, total_time = 0.1;
  double time = 0.0;
  int steps = 0;
  while (time < total_time)
  {
    tmc_solver_update_state_parameters(s, st);
    tmc_result result = tmc_solver_solve(s, dt, st);
    if (result.state != TMC_CONVERGED)
    {
      tmc_test_fail(__FILE__, __LINE__, "DAE solve did not converge at step %d, time=%g: %s", steps, time, tmc_solver_state_name(result.state));
      goto done;
    }
    const double residual = parameter(st, B_C_eq) * tmc_state_get(st, 0, B_idx) - tmc_state_get(st, 0, C_idx);
    TMC_EXPECT_NEAR(residual, 0.0, 1.0e-6);
    time += dt;
    steps++;
  }
  if (!reorder_state)
  {
    const double expected_C = parameter(st, B_C_eq) * tmc_state_get(st, 0, B_idx);
    TMC_EXPECT_NEAR(tmc_state_get(st, 0, C_idx), expected_C, 1.0e-6);
    TMC_EXPECT_NEAR(parameter(st, B_C_eq) * tmc_state_get(st, 0, B_idx) - tmc_state_get(st, 0, C_idx), 0.0, 1.0e-6);
  }
  /* A + B is conserved by A -> B */
  TMC_EXPECT_NEAR(tmc_state_get(st, 0, A_idx) + tmc_state_get(st, 0, B_idx), 1.0, 0.01);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
}

TMC_TEST(dae_solve_with_constraint)
{
  solve_a_to_b_with_constraint(false, 270.0, true);
  return TMC_TEST_PASS;
}

TMC_TEST(dae_solve_with_constraint_and_reorder_state)
{
  solve_a_to_b_with_constraint(true, 400.0, false);
  return TMC_TEST_PASS;
}

/* A -> B (kinetic), B <-> C and B <-> D (constraints sharing B). */
TMC_TEST(dae_solve_with_two_coupled_constraints)
{
  const tmc_species A = species_named("A"), B = species_named("B"), C = species_named("C"), D = species_named("D");
  tmc_system system = cst_system((const tmc_species[]){ A, B, C, D }, 4);
  tmc_reaction r = cst_reaction(&A, 1, (const tmc_stoich_species[]){ { B, 1.0 } }, 1, cst_arrhenius(1.0));
  const double K_eq1 = 3.0, K_eq2 = 5.0;
  tmc_constraint c[2] = { pair_equilibrium("B_C_eq", &B, &C, K_eq1, -2400.0), pair_equilibrium("B_D_eq", &B, &D, K_eq2, -2400.0) };
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(tmc_rosenbrock_three_stage(), &system, &r, 1, c, 2, WIDTH, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  const int A_idx = index_of(st, "A"), B_idx = index_of(st, "B"), C_idx = index_of(st, "C"), D_idx = index_of(st, "D");
  const int B_C_eq = tmc_state_custom_parameter_index(st, "B_C_eq"), B_D_eq = tmc_state_custom_parameter_index(st, "B_D_eq");
  CST_REQUIRE(A_idx >= 0 && B_idx >= 0 && C_idx >= 0 && D_idx >= 0 && B_C_eq >= 0 && B_D_eq >= 0);
  tmc_state_set(st, 0, A_idx, 1.0);
  tmc_state_set(st, 0, B_idx, 0.0);
  tmc_state_set(st, 0, C_idx, 0.0);
  tmc_state_set(st, 0, D_idx, 0.0);
  set_conditions(st, 300.0);
  tmc_solver_update_state_parameters(s, st);
  TMC_EXPECT_NEAR(parameter(st, B_C_eq), cst_equilibrium_constant(K_eq1, -2400.0, 300.0), 1e-10);
  TMC_EXPECT_NEAR(parameter(st, B_D_eq), cst_equilibrium_constant(K_eq2, -2400.0, 300.0), 1e-10);
  const double dt = 0.001, total_time = 0.1;
  double time = 0.0;
  while (time < total_time)
  {
    tmc_solver_update_state_parameters(s, st);
    tmc_result result = tmc_solver_solve(s, dt, st);
    if (result.state != TMC_CONVERGED)
    {
      tmc_test_fail(__FILE__, __LINE__, "Coupled constraints did not converge at time=%g: %s", time, tmc_solver_state_name(result.state));
      goto done;
    }
    const double B_val = tmc_state_get(st, 0, B_idx);
    TMC_EXPECT_NEAR(parameter(st, B_C_eq) * B_val - tmc_state_get(st, 0, C_idx), 0.0, 1.0e-6);
    TMC_EXPECT_NEAR(parameter(st, B_D_eq) * B_val - tmc_state_get(st, 0, D_idx), 0.0, 1.0e-6);
    time += dt;
  }
  TMC_EXPECT_NEAR(parameter(st, B_C_eq) * tmc_state_get(st, 0, B_idx), tmc_state_get(st, 0, C_idx), 1.0e-6);
  TMC_EXPECT_NEAR(parameter(st, B_D_eq) * tmc_state_get(st, 0, B_idx), tmc_state_get(st, 0, D_idx), 1.0e-6);
  TMC_EXPECT_NEAR(tmc_state_get(st, 0, A_idx) + tmc_state_get(st, 0, B_idx), 1.0, 0.01);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c[0]);
  tmc_constraint_free(&c[1]);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  return TMC_TEST_PASS;
}

/* 2A <-> B: K_eq [A]^2 - [B] = 0, with C -> A feeding A. */
TMC_TEST(dae_solve_with_non_unit_stoichiometry)
{
  const tmc_species A = species_named("A"), B = species_named("B"), C = species_named("C");
  tmc_system system = cst_system((const tmc_species[]){ A, B, C }, 3);
  tmc_reaction r = cst_reaction(&C, 1, (const tmc_stoich_species[]){ { A, 1.0 } }, 1, cst_arrhenius(0.5));
  const double K_eq = 10.0;
  tmc_constraint c = cst_equilibrium("A2_B_eq", &B, (const tmc_stoich_species[]){ { A, 2.0 } }, 1, (const tmc_stoich_species[]){ { B, 1.0 } }, 1, K_eq, -2400.0);
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(tmc_rosenbrock_three_stage(), &system, &r, 1, &c, 1, WIDTH, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  const int A_idx = index_of(st, "A"), B_idx = index_of(st, "B"), C_idx = index_of(st, "C");
  const int A2_B_eq = tmc_state_custom_parameter_index(st, "A2_B_eq");
  CST_REQUIRE(A_idx >= 0 && B_idx >= 0 && C_idx >= 0 && A2_B_eq >= 0);
  tmc_state_set(st, 0, A_idx, 0.1);
  tmc_state_set(st, 0, B_idx, K_eq * 0.1 * 0.1); /* B = K_eq A^2 */
  tmc_state_set(st, 0, C_idx, 1.0);
  set_conditions(st, 298.0);
  tmc_solver_update_state_parameters(s, st);
  TMC_EXPECT_NEAR(parameter(st, A2_B_eq), cst_equilibrium_constant(K_eq, -2400.0, 298.0), 1e-10);
  const double dt = 0.001, total_time = 0.05;
  double time = 0.0;
  while (time < total_time)
  {
    tmc_solver_update_state_parameters(s, st);
    tmc_result result = tmc_solver_solve(s, dt, st);
    if (result.state != TMC_CONVERGED)
    {
      tmc_test_fail(__FILE__, __LINE__, "NonUnit stoich did not converge at time=%g: %s", time, tmc_solver_state_name(result.state));
      goto done;
    }
    const double A_val = tmc_state_get(st, 0, A_idx), B_val = tmc_state_get(st, 0, B_idx);
    TMC_EXPECT_NEAR(parameter(st, A2_B_eq) * A_val * A_val - B_val, 0.0, 1.0e-5);
    time += dt;
  }
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("equilibrium_constraint:EquilibriumIntegration.SetConstraintsAPIWorks", set_constraints_api_works),
  TMC_CASE("equilibrium_constraint:EquilibriumIntegration.SetConstraintsAPIMultipleConstraints", set_constraints_api_multiple_constraints),
  TMC_CASE("equilibrium_constraint:EquilibriumIntegration.DAESolveWithConstraint", dae_solve_with_constraint),
  TMC_CASE("equilibrium_constraint:EquilibriumIntegration.DAESolveWithConstraintAndReorderState", dae_solve_with_constraint_and_reorder_state),
  TMC_CASE("equilibrium_constraint:EquilibriumIntegration.DAESolveWithTwoCoupledConstraints", dae_solve_with_two_coupled_constraints),
  TMC_CASE("equilibrium_constraint:EquilibriumIntegration.DAESolveWithNonUnitStoichiometry", dae_solve_with_non_unit_stoichiometry),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_equilibrium_constraint(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
