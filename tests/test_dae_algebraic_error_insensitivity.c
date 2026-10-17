/* MICM test/integration/test_dae_algebraic_error_insensitivity.cpp with
 * dae_algebraic_error_insensitivity_policy.hpp: the Rosenbrock DAE solver's
 * step-change error estimate for algebraic variables makes it sensitive to
 * their tolerances and prevents overshoot.
 *
 * System: A_gas <-> A_aq (K1, A_aq algebraic), A_aq <-> B_aq (K2, B_aq
 * algebraic), B_aq -> P (k), and A_aq + B_aq + P + A_gas = C_total (A_gas the
 * algebraic balance). The builder is CpuSolverBuilder with its default Matrix
 * and standard sparse ordering: vector width 1 and the builder's Doolittle LU
 * (cst_build), with the four-stage DAE parameters. */
#include "constraint_support.h"

typedef struct cascade_result
{
  double A_gas_final, P_final, min_A_gas;
  long accepted, rejected;
} cascade_result;

/* RunCascadeSystem: 30 s from equilibrium with P = 0. */
static cascade_result run_cascade_system(double balance_atol, double k, double K1, double K2)
{
  cascade_result out = { NAN, NAN, 1e99, 0, 0 };
  const tmc_species A_gas = tmc_species_create("A_gas"), A_aq = tmc_species_create("A_aq"), B_aq = tmc_species_create("B_aq"), P = tmc_species_create("P");
  tmc_system system = cst_system((const tmc_species[]){ A_gas, A_aq, B_aq, P }, 4);
  const double C_total = 1.0e-6;
  tmc_reaction r = cst_reaction(&B_aq, 1, (const tmc_stoich_species[]){ { P, 1.0 } }, 1, cst_arrhenius(k));
  tmc_constraint c[3] = {
    cst_equilibrium("eq1", &A_aq, (const tmc_stoich_species[]){ { A_gas, 1.0 } }, 1, (const tmc_stoich_species[]){ { A_aq, 1.0 } }, 1, K1, 0.0),
    cst_equilibrium("eq2", &B_aq, (const tmc_stoich_species[]){ { A_aq, 1.0 } }, 1, (const tmc_stoich_species[]){ { B_aq, 1.0 } }, 1, K2, 0.0),
    cst_linear("mass", &A_gas, (const tmc_stoich_species[]){ { A_aq, 1.0 }, { B_aq, 1.0 }, { P, 1.0 }, { A_gas, 1.0 } }, 4, C_total),
  };
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(tmc_rosenbrock_four_stage_dae(), &system, &r, 1, c, 3, 1, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  tmc_state_set_relative_tolerance(st, 1.0e-6);
  const int gi = tmc_state_variable_index(st, "A_gas"), ai = tmc_state_variable_index(st, "A_aq");
  const int bi = tmc_state_variable_index(st, "B_aq"), pi = tmc_state_variable_index(st, "P");
  /* the balance variable under test, moderate for the equilibrium variables,
   * tight for P, the only differential variable */
  double atols[4] = { 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12 };
  atols[gi] = balance_atol;
  atols[ai] = 1.0e-8;
  atols[bi] = 1.0e-8;
  atols[pi] = 1.0e-12;
  CST_REQUIRE(tmc_state_set_absolute_tolerances(st, atols, 4, &err) == TMC_OK);
  const double denom = 1.0 + K1 + K1 * K2;
  const double Ag0 = C_total / denom;
  tmc_state_set(st, 0, gi, Ag0);
  tmc_state_set(st, 0, ai, K1 * Ag0);
  tmc_state_set(st, 0, bi, K2 * K1 * Ag0);
  tmc_state_set(st, 0, pi, 0.0);
  tmc_state_conditions(st)[0].temperature = 298.0;
  tmc_state_conditions(st)[0].pressure = 101325.0;
  tmc_solver_update_state_parameters(s, st);
  const double dt = 30.0;
  double advanced = 0.0;
  while (advanced < dt)
  {
    tmc_result result = tmc_solver_solve(s, dt - advanced, st);
    TMC_EXPECT_EQ(result.state, TMC_CONVERGED);
    if (result.state != TMC_CONVERGED)
      break;
    advanced += result.stats.final_time;
    out.accepted += result.stats.accepted;
    out.rejected += result.stats.rejected;
    out.min_A_gas = fmin(out.min_A_gas, tmc_state_get(st, 0, gi));
  }
  out.A_gas_final = tmc_state_get(st, 0, gi);
  out.P_final = tmc_state_get(st, 0, pi);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  for (int i = 0; i < 3; ++i)
    tmc_constraint_free(&c[i]);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  return out;
}

/* A tighter tolerance on the balance variable takes more steps during the
 * transient (moderate stiffness: k = 100, K1 = K2 = 2, time constant ~0.018 s). */
TMC_TEST(error_sensitive_to_balance_atol)
{
  cascade_result loose = run_cascade_system(1e-3, 100.0, 2.0, 2.0);
  cascade_result tight = run_cascade_system(1e-8, 100.0, 2.0, 2.0);
  if (!(tight.accepted > loose.accepted))
    tmc_test_fail(__FILE__, __LINE__, "Tight atol should require more steps than loose atol. loose=%ld tight=%ld", loose.accepted, tight.accepted);
  return TMC_TEST_PASS;
}

/* Ultra-stiff (k = 1e4, K1 = 100, K2 = 50): the balance variable does not
 * overshoot deeply negative, and the system converts fully. */
TMC_TEST(algebraic_variable_does_not_overshoot_deeply)
{
  cascade_result r = run_cascade_system(1e-8, 1.0e4, 100.0, 50.0);
  if (!(r.min_A_gas >= -1.0e-8))
    tmc_test_fail(__FILE__, __LINE__, "A_gas overshot deeply negative: min_A_gas=%g", r.min_A_gas);
  TMC_EXPECT_NEAR(r.P_final, 1.0e-6, 1.0e-8);
  TMC_EXPECT_NEAR(r.A_gas_final, 0.0, 1.0e-10);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("dae_algebraic_error_insensitivity:DAEAlgebraicError.ErrorSensitiveToBalanceAtol", error_sensitive_to_balance_atol),
  TMC_CASE("dae_algebraic_error_insensitivity:DAEAlgebraicError.AlgebraicVariableDoesNotOvershootDeeply", algebraic_variable_does_not_overshoot_deeply),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_dae_algebraic_error_insensitivity(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
