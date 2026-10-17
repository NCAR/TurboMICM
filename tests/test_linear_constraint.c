/* MICM test/integration/test_linear_constraint.cpp with linear_constraint_policy.hpp.
 *
 * The Terminator and Robertson mechanisms side by side, Robertson closed by
 * the DAE constraint A + B + C = 1 with C algebraic. The builder is
 * CpuSolverBuilder with VectorMatrix<MICM_DEFAULT_VECTOR_SIZE> and the
 * four-stage DAE parameters: vector width 4 and the builder's Doolittle LU
 * (cst_build). MICM's precision matchers reduce to EXPECT_NEAR in double:
 * EXPECT_REAL_CLOSE(a, b, tol) is |a - b| <= tol, EXPECT_REAL_SOLVE_REL(a, b,
 * rel) is |a - b| <= rel |b|. */
#include "constraint_support.h"

#define WIDTH 4 /* MICM_DEFAULT_VECTOR_SIZE */

static bool finite_and_non_negative(const char *name, double value, int step)
{
  if (!isfinite(value))
  {
    tmc_test_fail(__FILE__, __LINE__, "%s not finite at step %d", name, step);
    return false;
  }
  if (!(value >= 0.0))
  {
    tmc_test_fail(__FILE__, __LINE__, "%s went negative at step %d: %g", name, step, value);
    return false;
  }
  return true;
}

TMC_TEST(terminator_and_robertson)
{
  tmc_species Cl2 = tmc_species_create_with_property("Cl2", TMC_PROPERTY_ABSOLUTE_TOLERANCE, 1.0e-20);
  tmc_species Cl = tmc_species_create_with_property("Cl", TMC_PROPERTY_ABSOLUTE_TOLERANCE, 1.0e-20);
  const tmc_species A = tmc_species_create("A"), B = tmc_species_create("B"), C = tmc_species_create("C");
  tmc_system system = cst_system((const tmc_species[]){ Cl2, Cl, A, B, C }, 5);
  tmc_reaction reactions[5] = {
    /* terminator: Cl2 -> 2 Cl (terminator_k1, left at zero), Cl + Cl -> Cl2 (k = 1) */
    cst_reaction(&Cl2, 1, (const tmc_stoich_species[]){ { Cl, 2.0 } }, 1, tmc_rate_user_defined("terminator_k1", 1.0)),
    cst_reaction((const tmc_species[]){ Cl, Cl }, 2, (const tmc_stoich_species[]){ { Cl2, 1.0 } }, 1, cst_arrhenius(1.0)),
    /* Robertson: A -> B, B + B -> B + C, B + C -> A + C */
    cst_reaction(&A, 1, (const tmc_stoich_species[]){ { B, 1.0 } }, 1, tmc_rate_user_defined("robertson_r1", 1.0)),
    cst_reaction((const tmc_species[]){ B, B }, 2, (const tmc_stoich_species[]){ { B, 1.0 }, { C, 1.0 } }, 2, tmc_rate_user_defined("robertson_r2", 1.0)),
    cst_reaction((const tmc_species[]){ B, C }, 2, (const tmc_stoich_species[]){ { A, 1.0 }, { C, 1.0 } }, 2, tmc_rate_user_defined("robertson_r3", 1.0)),
  };
  const double sum_initial_conc = 1.0;
  tmc_constraint c = cst_linear("mass_conservation", &C, (const tmc_stoich_species[]){ { A, 1.0 }, { B, 1.0 }, { C, 1.0 } }, 3, sum_initial_conc);
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = cst_build(tmc_rosenbrock_four_stage_dae(), &system, reactions, 5, &c, 1, WIDTH, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  tmc_state_set_relative_tolerance(st, 1.0e-8);
  CST_REQUIRE(tmc_state_set_custom_rate_parameter(st, "robertson_r1", 0.04, &err) == TMC_OK);
  CST_REQUIRE(tmc_state_set_custom_rate_parameter(st, "robertson_r2", 3.0e7, &err) == TMC_OK);
  CST_REQUIRE(tmc_state_set_custom_rate_parameter(st, "robertson_r3", 1.0e4, &err) == TMC_OK);
  CST_REQUIRE(tmc_state_set_concentration(st, "Cl", 1.2e-6, &err) == TMC_OK);
  CST_REQUIRE(tmc_state_set_concentration(st, "Cl2", 1.8e-10, &err) == TMC_OK);
  CST_REQUIRE(tmc_state_set_concentration(st, "A", 1.0, &err) == TMC_OK);
  CST_REQUIRE(tmc_state_set_concentration(st, "B", 0.0, &err) == TMC_OK);
  CST_REQUIRE(tmc_state_set_concentration(st, "C", 0.0, &err) == TMC_OK);
  tmc_state_conditions(st)[0].temperature = 298.0;
  tmc_state_conditions(st)[0].pressure = 101300.0;
  tmc_state_conditions(st)[0].air_density = 42.0;
  /* total chlorine 2 [Cl2] + [Cl] is invariant under both terminator reactions */
  const double total_chlorine_initial = 2.0 * 1.8e-10 + 1.2e-6;
  const int iCl2 = tmc_state_variable_index(st, "Cl2"), iCl = tmc_state_variable_index(st, "Cl");
  const int iA = tmc_state_variable_index(st, "A"), iB = tmc_state_variable_index(st, "B"), iC = tmc_state_variable_index(st, "C");
  tmc_solver_update_state_parameters(s, st);

  const int N = 12;
  double time_step = 1.0;
  for (int i = 0; i < N; ++i)
  {
    double advanced = 0.0;
    while (advanced < time_step)
    {
      tmc_result result = tmc_solver_solve(s, time_step - advanced, st);
      if (!(result.state == TMC_CONVERGED || result.state == TMC_CONVERGENCE_EXCEEDED_MAX_STEPS))
      {
        tmc_test_fail(__FILE__, __LINE__, "solver returned %s at step %d, t = %g", tmc_solver_state_name(result.state), i, advanced);
        goto done;
      }
      const double advanced_to = advanced + result.stats.final_time;
      if (!(advanced_to > advanced))
      {
        tmc_test_fail(__FILE__, __LINE__, "solve advanced only %g s at step %d", result.stats.final_time, i);
        goto done;
      }
      advanced = advanced_to;
    }
    const double cl2 = tmc_state_get(st, 0, iCl2), cl = tmc_state_get(st, 0, iCl);
    const double a = tmc_state_get(st, 0, iA), b = tmc_state_get(st, 0, iB), cc = tmc_state_get(st, 0, iC);
    /* 1. mass conservation, enforced by the DAE constraint */
    TMC_EXPECT_NEAR(a + b + cc, sum_initial_conc, 1e-10);
    /* 2. total chlorine, conserved by the kinetics to the integration's accuracy */
    TMC_EXPECT_NEAR(2.0 * cl2 + cl, total_chlorine_initial, 1.0e-6 * fabs(total_chlorine_initial));
    /* 3. nothing negative or non-finite */
    finite_and_non_negative("Cl2", cl2, i);
    finite_and_non_negative("Cl", cl, i);
    finite_and_non_negative("A", a, i);
    finite_and_non_negative("B", b, i);
    finite_and_non_negative("C", cc, i);
    /* 4. the Robertson system moved: after 1 s, A ~ 0.966, B ~ 3.1e-5, C ~ 0.033 */
    if (i == 0)
    {
      TMC_EXPECT_TRUE(a < 0.99);
      TMC_EXPECT_TRUE(b > 0.0);
      TMC_EXPECT_TRUE(cc > 1.0e-3);
    }
    time_step *= 10.0;
  }
  /* 5. at t ~ 1e11 s Robertson has run to completion and the terminator
   *    has turned nearly all Cl into Cl2 */
  TMC_EXPECT_TRUE(tmc_state_get(st, 0, iC) > 0.95);
  TMC_EXPECT_TRUE(tmc_state_get(st, 0, iA) < 0.05);
  TMC_EXPECT_TRUE(tmc_state_get(st, 0, iCl2) > 0.4 * total_chlorine_initial);
  TMC_EXPECT_TRUE(tmc_state_get(st, 0, iCl) < 0.1 * total_chlorine_initial);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  for (int i = 0; i < 5; ++i)
    tmc_reaction_free(&reactions[i]);
  tmc_system_free(&system);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("linear_constraint:DAESolveWithConstraint.TerminatorAndRobertson", terminator_and_robertson),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_linear_constraint(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
