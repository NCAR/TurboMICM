/* MICM test/unit/solver/test_constraint_initialization.cpp.
 *
 * MICM's StandardBuilder (CpuSolverBuilder with Matrix and the standard
 * sparse ordering) is a solver of vector width 1 with the builder's Doolittle
 * LU (cst_build). Where MICM calls solver.solver_.InitializeConstraints(state,
 * parameters, stats) directly, the port calls
 * tmc_solver_initialize_constraints, the same projection without a step. */
#include "constraint_support.h"

#define K_RATE 0.5 /* SimpleConstrainedSystem::K */
#define K_EQ 2.0   /* SimpleConstrainedSystem::K_EQ */

/* A -> B (k = 0.5) with the equilibrium constraint K_EQ B - C = 0, C algebraic,
 * no temperature dependence, no reordering. */
static tmc_solver *simple_constrained_system(tmc_rosenbrock_parameters parameters, tmc_error *err)
{
  const tmc_species species[] = { tmc_species_create("A"), tmc_species_create("B"), tmc_species_create("C") };
  tmc_system system = cst_system(species, 3);
  tmc_reaction r = cst_reaction(&species[0], 1, (const tmc_stoich_species[]){ { species[1], 1.0 } }, 1, cst_arrhenius(K_RATE));
  tmc_constraint c = cst_equilibrium("B_C_eq", &species[2], (const tmc_stoich_species[]){ { species[1], 1.0 } }, 1, (const tmc_stoich_species[]){ { species[2], 1.0 } }, 1, K_EQ, 0.0);
  tmc_solver *s = cst_build(parameters, &system, &r, 1, &c, 1, 1, false, err);
  tmc_constraint_free(&c);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  return s;
}

/* Projection-only: Z is slaved to X by row_scale (X - Z) = 0, so Z = X for every scale. */
static tmc_solver *scaled_linear_constraint_solver(tmc_rosenbrock_parameters parameters, double row_scale, tmc_error *err)
{
  const tmc_species species[] = { tmc_species_create("X"), tmc_species_create("Z") };
  tmc_system system = cst_system(species, 2);
  tmc_constraint c = cst_linear("scaled_copy", &species[1], (const tmc_stoich_species[]){ { species[0], row_scale }, { species[1], -row_scale } }, 2, 0.0);
  tmc_solver *s = cst_build(parameters, &system, NULL, 0, &c, 1, 1, false, err);
  tmc_constraint_free(&c);
  tmc_system_free(&system);
  return s;
}

/* Projection-only and nonlinear in the algebraic variable: 2 [B] - [C]^2 = 0,
 * so C = 1 when B = 0.5. */
static tmc_solver *squared_constraint_solver(tmc_rosenbrock_parameters parameters, tmc_error *err)
{
  const tmc_species species[] = { tmc_species_create("B"), tmc_species_create("C") };
  tmc_system system = cst_system(species, 2);
  tmc_constraint c = cst_equilibrium("squared", &species[1], (const tmc_stoich_species[]){ { species[0], 1.0 } }, 1, (const tmc_stoich_species[]){ { species[1], 2.0 } }, 1, 2.0, 0.0);
  tmc_solver *s = cst_build(parameters, &system, NULL, 0, &c, 1, 1, false, err);
  tmc_constraint_free(&c);
  tmc_system_free(&system);
  return s;
}

static void set_conditions(tmc_state *st, int n_cells)
{
  for (int i = 0; i < n_cells; ++i)
  {
    tmc_state_conditions(st)[i].temperature = 298.15;
    tmc_state_conditions(st)[i].pressure = 101325.0;
  }
}

static double parameter(tmc_state *st, int cell, const char *label)
{
  return tmc_dense_get(tmc_state_custom_rate_parameters(st), cell, tmc_state_custom_parameter_index(st, label));
}

static void set_tolerances(tmc_state *st, double rtol, double atol)
{
  double atols[8];
  const int n = tmc_state_num_variables(st);
  for (int i = 0; i < n; ++i)
    atols[i] = atol;
  tmc_state_set_relative_tolerance(st, rtol);
  tmc_error err;
  if (tmc_state_set_absolute_tolerances(st, atols, n, &err) != TMC_OK)
    tmc_test_fail(__FILE__, __LINE__, "%s", err.message);
}

/* Consistent initial conditions need at most one initialization iteration. */
TMC_TEST(consistent_ics_unchanged)
{
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = simple_constrained_system(tmc_rosenbrock_three_stage(), &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  const int A = tmc_state_variable_index(st, "A"), B = tmc_state_variable_index(st, "B"), C = tmc_state_variable_index(st, "C");
  tmc_state_set(st, 0, A, 1.0);
  tmc_state_set(st, 0, B, 0.5);
  tmc_state_set(st, 0, C, K_EQ * 0.5); /* C = K_eq B */
  set_conditions(st, 1);
  tmc_solver_update_state_parameters(s, st);
  tmc_result r = tmc_solver_solve(s, 0.001, st);
  TMC_EXPECT_EQ(r.state, TMC_CONVERGED);
  TMC_EXPECT_TRUE(r.stats.constraint_init_iterations <= 1);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  return TMC_TEST_PASS;
}

/* Solve from C off the manifold; the constraint holds afterwards. MICM also
 * checks that the initialization ran for the mildly inconsistent case. */
static void inconsistent_ics(double C0, bool check_iterations)
{
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = simple_constrained_system(tmc_rosenbrock_three_stage(), &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  const int A = tmc_state_variable_index(st, "A"), B = tmc_state_variable_index(st, "B"), C = tmc_state_variable_index(st, "C");
  tmc_state_set(st, 0, A, 1.0);
  tmc_state_set(st, 0, B, 0.5);
  tmc_state_set(st, 0, C, C0); /* K_eq B = 1 */
  set_conditions(st, 1);
  tmc_solver_update_state_parameters(s, st);
  tmc_result r = tmc_solver_solve(s, 0.001, st);
  TMC_EXPECT_EQ(r.state, TMC_CONVERGED);
  if (check_iterations)
    TMC_EXPECT_TRUE(r.stats.constraint_init_iterations > 0);
  const double K_eq_actual = parameter(st, 0, "B_C_eq");
  TMC_EXPECT_NEAR(K_eq_actual * tmc_state_get(st, 0, B) - tmc_state_get(st, 0, C), 0.0, 1.0e-6);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
}

TMC_TEST(mildly_inconsistent_ics_corrected)
{
  inconsistent_ics(5.0, true);
  return TMC_TEST_PASS;
}

TMC_TEST(severely_inconsistent_ics_converge)
{
  inconsistent_ics(1000.0, false);
  return TMC_TEST_PASS;
}

#define TOTAL 1.3
#define A0 0.7
#define B0 0.15

/* A -> B with the conservation row A + B + C = TOTAL multiplied by scale;
 * returns the solver state and C after one solve from C = 5. */
static tmc_solver_state solve_with_row_scale(double scale, double *C_after)
{
  tmc_solver_state result = TMC_NOT_YET_CALLED;
  *C_after = NAN;
  const tmc_species species[] = { tmc_species_create("A"), tmc_species_create("B"), tmc_species_create("C") };
  tmc_system system = cst_system(species, 3);
  tmc_reaction r = cst_reaction(&species[0], 1, (const tmc_stoich_species[]){ { species[1], 1.0 } }, 1, cst_arrhenius(0.5));
  tmc_constraint c = cst_linear("mass", &species[2], (const tmc_stoich_species[]){ { species[0], scale }, { species[1], scale }, { species[2], scale } }, 3, scale * TOTAL);
  tmc_error err;
  tmc_state *st = NULL;
  tmc_solver *s = cst_build(tmc_rosenbrock_four_stage_dae(), &system, &r, 1, &c, 1, 1, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  tmc_state_set(st, 0, tmc_state_variable_index(st, "A"), A0);
  tmc_state_set(st, 0, tmc_state_variable_index(st, "B"), B0);
  tmc_state_set(st, 0, tmc_state_variable_index(st, "C"), 5.0); /* far off the manifold */
  set_conditions(st, 1);
  tmc_solver_update_state_parameters(s, st);
  result = tmc_solver_solve(s, 0.001, st).state;
  *C_after = tmc_state_get(st, 0, tmc_state_variable_index(st, "C"));
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  return result;
}

/* Row-scale invariance through the public solve, not just the projection. */
TMC_TEST(scaled_constraint_rows_solve_identically)
{
  const double C_CONSISTENT = TOTAL - A0 - B0;
  double reference_C, scaled_C;
  TMC_EXPECT_EQ(solve_with_row_scale(1.0, &reference_C), TMC_CONVERGED);
  TMC_EXPECT_NEAR(reference_C, C_CONSISTENT, 1.0e-5 * C_CONSISTENT);
  const double scales[] = { 1.0e-12, 1.0e12 };
  for (int i = 0; i < 2; ++i)
  {
    TMC_EXPECT_EQ(solve_with_row_scale(scales[i], &scaled_C), TMC_CONVERGED);
    TMC_EXPECT_NEAR(scaled_C, reference_C, 1.0e-5 * C_CONSISTENT);
  }
  return TMC_TEST_PASS;
}

/* A pure ODE system takes no constraint initialization. */
TMC_TEST(pure_ode_system_unaffected)
{
  const tmc_species species[] = { tmc_species_create("A"), tmc_species_create("B") };
  tmc_system system = cst_system(species, 2);
  tmc_reaction r = cst_reaction(&species[0], 1, (const tmc_stoich_species[]){ { species[1], 1.0 } }, 1, cst_arrhenius(0.5));
  tmc_error err;
  tmc_state *st = NULL;
  tmc_solver *s = cst_build(tmc_rosenbrock_three_stage(), &system, &r, 1, NULL, 0, 1, false, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  tmc_state_set(st, 0, tmc_state_variable_index(st, "A"), 1.0);
  tmc_state_set(st, 0, tmc_state_variable_index(st, "B"), 0.0);
  set_conditions(st, 1);
  tmc_result result = tmc_solver_solve(s, 0.01, st);
  TMC_EXPECT_EQ(result.state, TMC_CONVERGED);
  TMC_EXPECT_EQ(result.stats.constraint_init_iterations, 0);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  tmc_reaction_free(&r);
  tmc_system_free(&system);
  return TMC_TEST_PASS;
}

/* Three cells: consistent, mildly and severely inconsistent. */
TMC_TEST(multi_cell_systems)
{
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = simple_constrained_system(tmc_rosenbrock_three_stage(), &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 3);
  const int A = tmc_state_variable_index(st, "A"), B = tmc_state_variable_index(st, "B"), C = tmc_state_variable_index(st, "C");
  const double values[3][3] = { { 1.0, 0.5, K_EQ * 0.5 }, { 2.0, 0.3, 5.0 }, { 0.5, 0.1, 100.0 } };
  for (int i = 0; i < 3; ++i)
  {
    tmc_state_set(st, i, A, values[i][0]);
    tmc_state_set(st, i, B, values[i][1]);
    tmc_state_set(st, i, C, values[i][2]);
  }
  set_conditions(st, 3);
  tmc_solver_update_state_parameters(s, st);
  tmc_result r = tmc_solver_solve(s, 0.001, st);
  TMC_EXPECT_EQ(r.state, TMC_CONVERGED);
  const double K_eq_actual = parameter(st, 0, "B_C_eq");
  for (int i = 0; i < 3; ++i)
    TMC_EXPECT_NEAR(K_eq_actual * tmc_state_get(st, i, B) - tmc_state_get(st, i, C), 0.0, 1.0e-5);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  return TMC_TEST_PASS;
}

/* A state moved off the manifold between solves is projected again. */
TMC_TEST(subsequent_solve_calls_reinitialize)
{
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = simple_constrained_system(tmc_rosenbrock_three_stage(), &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  const int A = tmc_state_variable_index(st, "A"), B = tmc_state_variable_index(st, "B"), C = tmc_state_variable_index(st, "C");
  tmc_state_set(st, 0, A, 1.0);
  tmc_state_set(st, 0, B, 0.5);
  tmc_state_set(st, 0, C, K_EQ * 0.5);
  set_conditions(st, 1);
  tmc_solver_update_state_parameters(s, st);
  TMC_EXPECT_EQ(tmc_solver_solve(s, 0.001, st).state, TMC_CONVERGED);
  tmc_state_set(st, 0, C, 999.0); /* an external event */
  tmc_solver_update_state_parameters(s, st);
  TMC_EXPECT_EQ(tmc_solver_solve(s, 0.001, st).state, TMC_CONVERGED);
  const double K_eq_actual = parameter(st, 0, "B_C_eq");
  TMC_EXPECT_NEAR(K_eq_actual * tmc_state_get(st, 0, B) - tmc_state_get(st, 0, C), 0.0, 1.0e-6);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  return TMC_TEST_PASS;
}

TMC_TEST(solver_state_to_string_new_value)
{
  TMC_EXPECT_STREQ(tmc_solver_state_name(TMC_CONSTRAINT_INITIALIZATION_FAILED), "Constraint Initialization Failed");
  return TMC_TEST_PASS;
}

/* A failed projection leaves the caller's state exactly as passed in. */
TMC_TEST(failed_initialization_restores_caller_state)
{
  tmc_rosenbrock_parameters parameters = tmc_rosenbrock_three_stage();
  parameters.constraint_init_max_iterations = 1; /* too few updates to reach C = 1 from C = 100 */
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = squared_constraint_solver(parameters, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  const int B = tmc_state_variable_index(st, "B"), C = tmc_state_variable_index(st, "C");
  const double B_in = 0.5, C_in = 100.0;
  tmc_state_set(st, 0, B, B_in);
  tmc_state_set(st, 0, C, C_in);
  set_conditions(st, 1);
  tmc_solver_update_state_parameters(s, st);
  tmc_stats stats;
  memset(&stats, 0, sizeof stats);
  const tmc_solver_state status = tmc_solver_initialize_constraints(s, st, &parameters, &stats);
  CST_REQUIRE(status == TMC_CONSTRAINT_INITIALIZATION_FAILED);
  TMC_EXPECT_EQ(tmc_state_get(st, 0, C), C_in);
  TMC_EXPECT_EQ(tmc_state_get(st, 0, B), B_in);
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  return TMC_TEST_PASS;
}

/* A cold start needs damping to converge within the default Newton budget. */
TMC_TEST(successful_initialization_keeps_the_correction)
{
  tmc_rosenbrock_parameters parameters = tmc_rosenbrock_three_stage();
  tmc_state *st = NULL;
  tmc_error err;
  tmc_solver *s = squared_constraint_solver(parameters, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  const int B = tmc_state_variable_index(st, "B"), C = tmc_state_variable_index(st, "C");
  tmc_state_set(st, 0, B, 0.5);
  tmc_state_set(st, 0, C, 1.0e-6);
  set_conditions(st, 1);
  tmc_solver_update_state_parameters(s, st);
  tmc_stats stats;
  memset(&stats, 0, sizeof stats);
  CST_REQUIRE(tmc_solver_initialize_constraints(s, st, &parameters, &stats) == TMC_CONVERGED);
  /* the remaining correction is a tenth of the state-error scale: within ~1e-4 with atol 1e-3 */
  TMC_EXPECT_NEAR(tmc_state_get(st, 0, C), 1.0, 1.0e-3);
  const int backtracks[] = { 0, 2 };
  for (int i = 0; i < 2; ++i)
  {
    parameters.constraint_init_max_backtracks = backtracks[i];
    tmc_state_set(st, 0, C, 1.0e-6);
    memset(&stats, 0, sizeof stats);
    TMC_EXPECT_EQ(tmc_solver_initialize_constraints(s, st, &parameters, &stats), TMC_CONSTRAINT_INITIALIZATION_FAILED);
    TMC_EXPECT_EQ(tmc_state_get(st, 0, B), 0.5);
    TMC_EXPECT_EQ(tmc_state_get(st, 0, C), 1.0e-6);
  }
done:
  tmc_state_free(st);
  tmc_solver_free(s);
  return TMC_TEST_PASS;
}

/* The weighted Newton-correction measure is invariant to complete constraint-row scaling. */
TMC_TEST(weighted_correction_is_invariant_to_row_scaling)
{
  const double row_scales[] = { 1.0e-12, 1.0, 1.0e12 };
  for (int i = 0; i < 3; ++i)
  {
    const tmc_rosenbrock_parameters parameters = tmc_rosenbrock_three_stage();
    const double rtol = 1.0e-8, atol = 1.0e-12;
    tmc_error err;
    tmc_solver *s = scaled_linear_constraint_solver(parameters, row_scales[i], &err);
    if (!s)
    {
      tmc_test_fail(__FILE__, __LINE__, "build failed: %s", err.message);
      continue;
    }
    tmc_state *st = tmc_solver_create_state(s, 1);
    set_tolerances(st, rtol, atol);
    const int X = tmc_state_variable_index(st, "X"), Z = tmc_state_variable_index(st, "Z");
    tmc_state_set(st, 0, X, 1.0);
    tmc_state_set(st, 0, Z, 4.0);
    tmc_stats stats;
    memset(&stats, 0, sizeof stats);
    const tmc_solver_state status = tmc_solver_initialize_constraints(s, st, &parameters, &stats);
    if (status != TMC_CONVERGED || fabs(tmc_state_get(st, 0, Z) - 1.0) > atol)
      tmc_test_fail(__FILE__, __LINE__, "row scale=%g: %s, Z=%.17g", row_scales[i], tmc_solver_state_name(status), tmc_state_get(st, 0, Z));
    tmc_state_free(st);
    tmc_solver_free(s);
  }
  return TMC_TEST_PASS;
}

/* constraint_init_tolerance is a fraction of the state-error scale: the same
 * off-manifold state is converged under loose tolerances and needs a Newton
 * update under tight ones. The iteration count is one more than the updates. */
static void initialize_with_tolerances(double rtol, double atol, double deviation, tmc_solver_state *status, double *Z_before, double *Z_after, long *iterations)
{
  const tmc_rosenbrock_parameters parameters = tmc_rosenbrock_three_stage();
  tmc_state *st = NULL;
  tmc_error err;
  *status = TMC_NOT_YET_CALLED;
  tmc_solver *s = scaled_linear_constraint_solver(parameters, 1.0, &err);
  CST_REQUIRE_BUILT(s, err);
  st = tmc_solver_create_state(s, 1);
  set_tolerances(st, rtol, atol);
  const int X = tmc_state_variable_index(st, "X"), Z = tmc_state_variable_index(st, "Z");
  tmc_state_set(st, 0, X, 1.0);
  tmc_state_set(st, 0, Z, 1.0 + deviation);
  *Z_before = tmc_state_get(st, 0, Z);
  tmc_stats stats;
  memset(&stats, 0, sizeof stats);
  *status = tmc_solver_initialize_constraints(s, st, &parameters, &stats);
  *Z_after = tmc_state_get(st, 0, Z);
  *iterations = stats.constraint_init_iterations;
done:
  tmc_state_free(st);
  tmc_solver_free(s);
}

TMC_TEST(weighted_correction_uses_state_tolerances)
{
  TMC_EXPECT_EQ(tmc_rosenbrock_three_stage().constraint_init_tolerance, 0.1);
  const double deviation = 5.0e-8;
  tmc_solver_state status;
  double Z_before = NAN, Z_after = NAN;
  long iterations = -1;
  /* loose: rtol 1e-6, atol 1e-8; already converged, no update applied */
  initialize_with_tolerances(1.0e-6, 1.0e-8, deviation, &status, &Z_before, &Z_after, &iterations);
  TMC_EXPECT_EQ(status, TMC_CONVERGED);
  TMC_EXPECT_EQ(Z_after, Z_before);
  TMC_EXPECT_EQ(iterations, 1);
  /* tight: rtol 1e-8, atol 1e-12; one update applied */
  initialize_with_tolerances(1.0e-8, 1.0e-12, deviation, &status, &Z_before, &Z_after, &iterations);
  TMC_EXPECT_EQ(status, TMC_CONVERGED);
  TMC_EXPECT_NEAR(Z_after, 1.0, 1.0e-12);
  TMC_EXPECT_EQ(iterations, 2);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("constraint_initialization:ConstraintInitialization.ConsistentICsUnchanged", consistent_ics_unchanged),
  TMC_CASE("constraint_initialization:ConstraintInitialization.MildlyInconsistentICsCorrected", mildly_inconsistent_ics_corrected),
  TMC_CASE("constraint_initialization:ConstraintInitialization.SeverelyInconsistentICsConverge", severely_inconsistent_ics_converge),
  TMC_CASE("constraint_initialization:ConstraintInitialization.ScaledConstraintRowsSolveIdenticallyThroughSolve", scaled_constraint_rows_solve_identically),
  TMC_CASE("constraint_initialization:ConstraintInitialization.PureODESystemUnaffected", pure_ode_system_unaffected),
  TMC_CASE("constraint_initialization:ConstraintInitialization.MultiCellSystems", multi_cell_systems),
  TMC_CASE("constraint_initialization:ConstraintInitialization.SubsequentSolveCallsReinitialize", subsequent_solve_calls_reinitialize),
  TMC_CASE("constraint_initialization:ConstraintInitialization.SolverStateToStringNewValue", solver_state_to_string_new_value),
  TMC_CASE("constraint_initialization:ConstraintInitialization.FailedInitializationRestoresCallerState", failed_initialization_restores_caller_state),
  TMC_CASE("constraint_initialization:ConstraintInitialization.SuccessfulInitializationKeepsTheCorrection", successful_initialization_keeps_the_correction),
  TMC_CASE("constraint_initialization:ConstraintInitialization.WeightedCorrectionIsInvariantToConstraintRowScaling", weighted_correction_is_invariant_to_row_scaling),
  TMC_CASE("constraint_initialization:ConstraintInitialization.WeightedCorrectionUsesStateTolerances", weighted_correction_uses_state_tolerances),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_constraint_initialization(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
