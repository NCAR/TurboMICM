/* MICM test/integration/test_dae_constraint_overshoot.cpp with
 * dae_constraint_overshoot_policy.hpp: a regression test for overshoot of an
 * algebraic balance variable in the Rosenbrock DAE solver. With fast A -> B
 * kinetics and the conservation constraint A + B + C = C_total (C algebraic),
 * C must not go negative.
 *
 * The builder is CpuSolverBuilder with VectorMatrix<MICM_DEFAULT_VECTOR_SIZE>:
 * vector width 4 and the builder's Doolittle LU (cst_build). */
#include "constraint_support.h"

#define WIDTH 4          /* MICM_DEFAULT_VECTOR_SIZE */
#define ABS_TOL 1.0e-12  /* kAbsTol in double */
#define C_TOTAL 1.0e-6

/* A -> B (k = 1e4), A + B + C = C_TOTAL with C algebraic, integrated over 30 s
 * from A = 0.9e-6, C = 0.1e-6; checks conservation and C >= 0 after every
 * solve, and the end state. */
static void algebraic_variable_stays_non_negative(tmc_rosenbrock_parameters parameters, const char *name)
{
  const tmc_species A = tmc_species_create("A"), B = tmc_species_create("B"), C = tmc_species_create("C");
  tmc_system system = cst_system((const tmc_species[]){ A, B, C }, 3);
  tmc_reaction r = cst_reaction(&A, 1, (const tmc_stoich_species[]){ { B, 1.0 } }, 1, cst_arrhenius(1.0e4));
  tmc_constraint c = cst_linear("mass_conservation", &C, (const tmc_stoich_species[]){ { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, 3, C_TOTAL);
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(parameters, &system, &r, 1, &c, 1, WIDTH, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  tmc_state_set_relative_tolerance(st, 1.0e-6);
  CST_REQUIRE(tmc_state_set_absolute_tolerances(st, (const double[]){ ABS_TOL, ABS_TOL, ABS_TOL }, 3, &err) == TMC_OK);
  const int A_idx = tmc_state_variable_index(st, "A"), B_idx = tmc_state_variable_index(st, "B"), C_idx = tmc_state_variable_index(st, "C");
  tmc_state_set(st, 0, A_idx, 0.9e-6);
  tmc_state_set(st, 0, B_idx, 0.0);
  tmc_state_set(st, 0, C_idx, 0.1e-6); /* C_total - A */
  tmc_state_conditions(st)[0].temperature = 298.0;
  tmc_state_conditions(st)[0].pressure = 101325.0;
  tmc_solver_update_state_parameters(s, st);
  const double dt = 30.0;
  double advanced = 0.0;
  while (advanced < dt)
  {
    tmc_result result = tmc_solver_solve(s, dt - advanced, st);
    if (result.state != TMC_CONVERGED)
    {
      tmc_test_fail(__FILE__, __LINE__, "Solver did not converge for %s at t=%g: %s", name, advanced, tmc_solver_state_name(result.state));
      goto done;
    }
    advanced += result.stats.final_time;
    const double a = tmc_state_get(st, 0, A_idx), b = tmc_state_get(st, 0, B_idx), cc = tmc_state_get(st, 0, C_idx);
    TMC_EXPECT_NEAR(a + b + cc, C_TOTAL, 1.0e-12);
    if (!(cc >= -1.0e-18))
      tmc_test_fail(__FILE__, __LINE__, "Algebraic variable C went negative for %s (%g) at t=%g; A=%g, B=%g", name, cc, advanced, a, b);
  }
  /* after 30 s with k = 1e4, A is gone and B holds A's initial amount; C keeps its own */
  TMC_EXPECT_TRUE(tmc_state_get(st, 0, A_idx) < 1.0e-12);
  TMC_EXPECT_NEAR(tmc_state_get(st, 0, B_idx), 0.9e-6, 1.0e-10);
  TMC_EXPECT_NEAR(tmc_state_get(st, 0, C_idx), 0.1e-6, 1.0e-10);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
}

TMC_TEST(algebraic_variable_stays_non_negative_test)
{
  algebraic_variable_stays_non_negative(tmc_rosenbrock_four_stage_dae(), "four-stage-dae");
  return TMC_TEST_PASS;
}

/* A_gas <-> A_aq (equilibrium, A_aq algebraic), A_aq -> P (fast), and
 * A_gas + A_aq + P = C_total (A_gas algebraic). */
TMC_TEST(equilibrium_plus_conservation)
{
  const tmc_species A_gas = tmc_species_create("A_gas"), A_aq = tmc_species_create("A_aq"), P = tmc_species_create("P");
  tmc_system system = cst_system((const tmc_species[]){ A_gas, A_aq, P }, 3);
  tmc_reaction r = cst_reaction(&A_aq, 1, (const tmc_stoich_species[]){ { P, 1.0 } }, 1, cst_arrhenius(1.0e3));
  const double K_eq = 10.0;
  tmc_constraint c[2] = {
    cst_equilibrium("gas_aq_eq", &A_aq, (const tmc_stoich_species[]){ { A_gas, 1.0 } }, 1, (const tmc_stoich_species[]){ { A_aq, 1.0 } }, 1, K_eq, 0.0),
    cst_linear("mass_conservation", &A_gas, (const tmc_stoich_species[]){ { A_aq, 1.0 }, { P, 1.0 }, { A_gas, 1.0 } }, 3, C_TOTAL),
  };
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(tmc_rosenbrock_four_stage_dae(), &system, &r, 1, c, 2, WIDTH, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  tmc_state_set_relative_tolerance(st, 1.0e-6);
  const int A_gas_idx = tmc_state_variable_index(st, "A_gas"), A_aq_idx = tmc_state_variable_index(st, "A_aq"), P_idx = tmc_state_variable_index(st, "P");
  /* tight for the differential P, moderate for the algebraic A_gas and A_aq */
  double atols[3] = { 1.0e-12, 1.0e-12, 1.0e-12 };
  atols[A_gas_idx] = 1.0e-8;
  atols[A_aq_idx] = 1.0e-8;
  CST_REQUIRE(tmc_state_set_absolute_tolerances(st, atols, 3, &err) == TMC_OK);
  /* equilibrium satisfied, no product yet */
  const double A_gas_init = C_TOTAL / (1.0 + K_eq);
  const double A_aq_init = K_eq * A_gas_init;
  tmc_state_set(st, 0, A_gas_idx, A_gas_init);
  tmc_state_set(st, 0, A_aq_idx, A_aq_init);
  tmc_state_set(st, 0, P_idx, 0.0);
  tmc_state_conditions(st)[0].temperature = 298.0;
  tmc_state_conditions(st)[0].pressure = 101325.0;
  tmc_solver_update_state_parameters(s, st);
  const double dt = 30.0;
  double advanced = 0.0;
  while (advanced < dt)
  {
    tmc_result result = tmc_solver_solve(s, dt - advanced, st);
    if (result.state != TMC_CONVERGED)
    {
      tmc_test_fail(__FILE__, __LINE__, "Solver did not converge at t=%g: %s", advanced, tmc_solver_state_name(result.state));
      goto done;
    }
    advanced += result.stats.final_time;
    const double g = tmc_state_get(st, 0, A_gas_idx), aq = tmc_state_get(st, 0, A_aq_idx), p = tmc_state_get(st, 0, P_idx);
    TMC_EXPECT_NEAR(g + aq + p, C_TOTAL, 1.0e-12);
    if (!(g >= -1.0e-18))
      tmc_test_fail(__FILE__, __LINE__, "A_gas went negative (%g) at t=%g", g, advanced);
    if (!(aq >= -1.0e-18))
      tmc_test_fail(__FILE__, __LINE__, "A_aq went negative (%g) at t=%g", aq, advanced);
  }
  /* after 30 s with k = 1e3, nearly everything is in P */
  TMC_EXPECT_NEAR(tmc_state_get(st, 0, P_idx), C_TOTAL, 1.0e-8);
  TMC_EXPECT_TRUE(tmc_state_get(st, 0, A_gas_idx) >= -1.0e-18);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c[0]);
  tmc_constraint_free(&c[1]);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  return TMC_TEST_PASS;
}

/* Every Rosenbrock parameter set but the two-stage one, which lacks the
 * stability for k = 1e4 (MICM covers it elsewhere). */
TMC_TEST(all_rosenbrock_orders_constrained)
{
  algebraic_variable_stays_non_negative(tmc_rosenbrock_three_stage(), "three-stage");
  algebraic_variable_stays_non_negative(tmc_rosenbrock_four_stage(), "four-stage");
  algebraic_variable_stays_non_negative(tmc_rosenbrock_four_stage_dae(), "four-stage-dae");
  algebraic_variable_stays_non_negative(tmc_rosenbrock_six_stage_dae(), "six-stage-dae");
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("dae_constraint_overshoot:DAEConstraintOvershoot.AlgebraicVariableStaysNonNegative", algebraic_variable_stays_non_negative_test),
  TMC_CASE("dae_constraint_overshoot:DAEConstraintOvershoot.EquilibriumPlusConservation", equilibrium_plus_conservation),
  TMC_CASE("dae_constraint_overshoot:DAEConstraintOvershoot.AllRosenbrockOrdersConstrained", all_rosenbrock_orders_constrained),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_dae_constraint_overshoot(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
