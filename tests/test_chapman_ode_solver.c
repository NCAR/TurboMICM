/* MICM test/unit/solver/test_chapman_ode_solver.cpp.
 *
 * MICM's cases exercise ChapmanODESolver, a hard-coded Chapman solver kept
 * under test/regression/RosenbrockChapman, through its parts: forcing,
 * Jacobian, factored alpha I - J, triangular solves, and Solve. Here each
 * check runs on the library's generic equivalent for the same mechanism,
 * built as the hard-coded solver orders it (species M, Ar, CO2, H2O, N2, O1D,
 * O, O2, O3, with M a state variable; no reordering; Ar, CO2, and H2O take no
 * part in reactions):
 *   - the hard-coded 23-element matrix numbering is column by column; the
 *     library's LU pattern holds those 23 elements plus (M, O), (M, O2),
 *     (N2, O1D), and (O2, O1D), where a species is both reactant and product
 *     and the derivative is identically zero;
 *   - the hard-coded factorization stores pivot reciprocals on the diagonal;
 *     the library stores the pivots and writes the reciprocals separately;
 *   - Force and DforceDy of the hard-coded solver are kept here as
 *     references for the library's forcing and Jacobian;
 *   - Solve uses the hard-coded solver's settings: the three-stage method,
 *     relative tolerance 1e-4, absolute 1e-3, first step 1e-6 s, at most 100
 *     steps. */
#include "harness.h"
#include "libcxx_random.h"
#include "tmc_solver_ext.h"

static const double ABSOLUTE_TOLERANCE = 1e-4;

enum
{
  M,
  Ar,
  CO2,
  H2O,
  N2,
  O1D,
  O,
  O2,
  O3,
  N_SPECIES
};
static const char *const species_names[N_SPECIES] = { "M", "Ar", "CO2", "H2O", "N2", "O1D", "O", "O2", "O3" };
static const char *const photolysis_names[3] = { "O2_1", "O3_1", "O3_2" };

/* The hard-coded solver's matrix elements (row, column), numbered column by column. */
enum
{
  N_ELEMENTS = 23
};
static const int elements[N_ELEMENTS][2] = { { M, M },     { O, M },     { O2, M },   { O3, M },   { Ar, Ar },  { CO2, CO2 },
                                             { H2O, H2O }, { N2, N2 },   { O1D, N2 }, { O, N2 },   { O1D, O1D }, { O, O1D },
                                             { O, O },     { O2, O },    { O3, O },   { O1D, O2 }, { O, O2 },    { O2, O2 },
                                             { O3, O2 },   { O1D, O3 },  { O, O3 },   { O2, O3 },  { O3, O3 } };

/* ---- the mechanism ---------------------------------------------------------- */

static tmc_rate_constant arrhenius(double A, double B, double C)
{
  tmc_arrhenius_parameters p = tmc_arrhenius_defaults();
  p.A = A;
  p.B = B;
  p.C = C;
  return tmc_rate_arrhenius(p);
}

/* The Chapman mechanism. The library orders reactions by rate-constant type
 * (stably), so the four Arrhenius reactions come first, as the hard-coded
 * solver's rate constants 0-3, then the photolysis rates 4-6. */
static tmc_solver *chapman_solver(int L)
{
  tmc_species s[N_SPECIES];
  tmc_phase gas = tmc_phase_create("gas");
  for (int i = 0; i < N_SPECIES; ++i)
  {
    s[i] = tmc_species_create(species_names[i]);
    tmc_phase_add_species(&gas, &s[i]);
  }
  tmc_system system = tmc_system_create(&gas);
  tmc_reaction r[7];
  for (int i = 0; i < 7; ++i)
  {
    r[i] = tmc_reaction_create();
    tmc_reaction_set_phase(&r[i], &gas);
  }
  /* N2_O1D_1: N2 + O1D -> O + N2 */
  tmc_reaction_add_reactant(&r[0], &s[N2]);
  tmc_reaction_add_reactant(&r[0], &s[O1D]);
  tmc_reaction_add_product(&r[0], &s[O], 1);
  tmc_reaction_add_product(&r[0], &s[N2], 1);
  tmc_reaction_set_rate_constant(&r[0], arrhenius(2.15e-11, 0, 110));
  /* O1D_O2_1: O1D + O2 -> O + O2 */
  tmc_reaction_add_reactant(&r[1], &s[O1D]);
  tmc_reaction_add_reactant(&r[1], &s[O2]);
  tmc_reaction_add_product(&r[1], &s[O], 1);
  tmc_reaction_add_product(&r[1], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[1], arrhenius(3.3e-11, 0, 55));
  /* O_O3_1: O + O3 -> 2 O2 */
  tmc_reaction_add_reactant(&r[2], &s[O]);
  tmc_reaction_add_reactant(&r[2], &s[O3]);
  tmc_reaction_add_product(&r[2], &s[O2], 2);
  tmc_reaction_set_rate_constant(&r[2], arrhenius(8e-12, 0, -2060));
  /* M_O_O2_1: M + O + O2 -> O3 + M */
  tmc_reaction_add_reactant(&r[3], &s[M]);
  tmc_reaction_add_reactant(&r[3], &s[O]);
  tmc_reaction_add_reactant(&r[3], &s[O2]);
  tmc_reaction_add_product(&r[3], &s[O3], 1);
  tmc_reaction_add_product(&r[3], &s[M], 1);
  tmc_reaction_set_rate_constant(&r[3], arrhenius(6e-34, 2.4, 0));
  /* O2_1: O2 -> 2 O */
  tmc_reaction_add_reactant(&r[4], &s[O2]);
  tmc_reaction_add_product(&r[4], &s[O], 2);
  tmc_reaction_set_rate_constant(&r[4], tmc_rate_user_defined(photolysis_names[0], 1.0));
  /* O3_1: O3 -> O1D + O2 */
  tmc_reaction_add_reactant(&r[5], &s[O3]);
  tmc_reaction_add_product(&r[5], &s[O1D], 1);
  tmc_reaction_add_product(&r[5], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[5], tmc_rate_user_defined(photolysis_names[1], 1.0));
  /* O3_2: O3 -> O + O2 */
  tmc_reaction_add_reactant(&r[6], &s[O3]);
  tmc_reaction_add_product(&r[6], &s[O], 1);
  tmc_reaction_add_product(&r[6], &s[O2], 1);
  tmc_reaction_set_rate_constant(&r[6], tmc_rate_user_defined(photolysis_names[2], 1.0));

  tmc_rosenbrock_parameters parameters = tmc_rosenbrock_three_stage();
  parameters.max_number_of_steps = 100;
  parameters.h_start = 1.0e-6;
  tmc_builder *builder = tmc_builder_create_rosenbrock(parameters);
  tmc_builder_set_system(builder, &system);
  tmc_builder_set_reactions(builder, r, 7);
  tmc_builder_set_reorder_state(builder, false);
  tmc_builder_set_ignore_unused_species(builder, true);
  tmc_builder_set_width(builder, L);
  tmc_error err;
  tmc_error_clear(&err);
  tmc_solver *solver = tmc_builder_build(builder, &err);
  if (!solver)
    fprintf(stderr, "  build failed: %s\n", err.message);
  tmc_builder_free(builder);
  for (int i = 0; i < 7; ++i)
    tmc_reaction_free(&r[i]);
  tmc_system_free(&system);
  tmc_phase_free(&gas);
  return solver;
}

/* The library's reaction r is the hard-coded solver's rate constant r. */
static void expect_reaction_order(const tmc_solver *solver)
{
  int n = 0;
  const tmc_reaction *reactions = tmc_solver_reactions(solver, &n);
  TMC_EXPECT_EQ(n, 7);
  if (n != 7)
    return;
  static const double arrhenius_A[4] = { 2.15e-11, 3.3e-11, 8e-12, 6e-34 };
  for (int r = 0; r < 4; ++r)
  {
    TMC_EXPECT_EQ(reactions[r].rate_constant.type, TMC_RATE_ARRHENIUS);
    TMC_EXPECT_EQ(reactions[r].rate_constant.p.arrhenius.A, arrhenius_A[r]);
  }
  for (int r = 4; r < 7; ++r)
  {
    TMC_EXPECT_EQ(reactions[r].rate_constant.type, TMC_RATE_USER_DEFINED);
    TMC_EXPECT_STREQ(reactions[r].rate_constant.p.user_defined.label, photolysis_names[r - 4]);
  }
}

/* ---- the hard-coded solver's forcing and Jacobian, as references ------------ */

/* ChapmanODESolver::Force */
static void reference_force(const double *rate_constants, const double *number_densities, double *force)
{
  for (int i = 0; i < N_SPECIES; ++i)
    force[i] = 0.0;
  force[5] = force[5] + rate_constants[5] * number_densities[8];
  force[5] = force[5] - rate_constants[0] * number_densities[4] * number_densities[5];
  force[5] = force[5] - rate_constants[1] * number_densities[5] * number_densities[7];
  force[6] = force[6] + 2 * rate_constants[4] * number_densities[7];
  force[6] = force[6] + rate_constants[6] * number_densities[8];
  force[6] = force[6] + rate_constants[0] * number_densities[4] * number_densities[5];
  force[6] = force[6] + rate_constants[1] * number_densities[5] * number_densities[7];
  force[6] = force[6] - rate_constants[2] * number_densities[6] * number_densities[8];
  force[6] = force[6] - rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
  force[7] = force[7] - rate_constants[4] * number_densities[7];
  force[7] = force[7] + rate_constants[5] * number_densities[8];
  force[7] = force[7] + rate_constants[6] * number_densities[8];
  force[7] = force[7] + 2 * rate_constants[2] * number_densities[6] * number_densities[8];
  force[7] = force[7] - rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
  force[8] = force[8] - rate_constants[5] * number_densities[8];
  force[8] = force[8] - rate_constants[6] * number_densities[8];
  force[8] = force[8] - rate_constants[2] * number_densities[6] * number_densities[8];
  force[8] = force[8] + rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
}

/* ChapmanODESolver::DforceDy: J = d force / d y in the 23-element numbering */
static void reference_dforce_dy(const double *rate_constants, const double *number_densities, double *jacobian)
{
  for (int e = 0; e < N_ELEMENTS; ++e)
    jacobian[e] = 0.0;
  jacobian[1] = jacobian[1] - rate_constants[3] * number_densities[6] * number_densities[7];
  jacobian[2] = jacobian[2] - rate_constants[3] * number_densities[6] * number_densities[7];
  jacobian[3] = jacobian[3] + rate_constants[3] * number_densities[6] * number_densities[7];
  jacobian[8] = jacobian[8] - rate_constants[0] * number_densities[5];
  jacobian[9] = jacobian[9] + rate_constants[0] * number_densities[5];
  jacobian[10] = jacobian[10] - rate_constants[0] * number_densities[4];
  jacobian[10] = jacobian[10] - rate_constants[1] * number_densities[7];
  jacobian[11] = jacobian[11] + rate_constants[0] * number_densities[4];
  jacobian[11] = jacobian[11] + rate_constants[1] * number_densities[7];
  jacobian[12] = jacobian[12] - rate_constants[2] * number_densities[8];
  jacobian[12] = jacobian[12] - rate_constants[3] * number_densities[0] * number_densities[7];
  jacobian[13] = jacobian[13] + 2 * rate_constants[2] * number_densities[8];
  jacobian[13] = jacobian[13] - rate_constants[3] * number_densities[0] * number_densities[7];
  jacobian[14] = jacobian[14] - rate_constants[2] * number_densities[8];
  jacobian[14] = jacobian[14] + rate_constants[3] * number_densities[0] * number_densities[7];
  jacobian[15] = jacobian[15] - rate_constants[1] * number_densities[5];
  jacobian[16] = jacobian[16] + 2 * rate_constants[4];
  jacobian[16] = jacobian[16] + rate_constants[1] * number_densities[5];
  jacobian[16] = jacobian[16] - rate_constants[3] * number_densities[0] * number_densities[6];
  jacobian[17] = jacobian[17] - rate_constants[4];
  jacobian[17] = jacobian[17] - rate_constants[3] * number_densities[0] * number_densities[6];
  jacobian[18] = jacobian[18] + rate_constants[3] * number_densities[0] * number_densities[6];
  jacobian[19] = jacobian[19] + rate_constants[5];
  jacobian[20] = jacobian[20] + rate_constants[6];
  jacobian[20] = jacobian[20] - rate_constants[2] * number_densities[6];
  jacobian[21] = jacobian[21] + rate_constants[5];
  jacobian[21] = jacobian[21] + rate_constants[6];
  jacobian[21] = jacobian[21] + 2 * rate_constants[2] * number_densities[6];
  jacobian[22] = jacobian[22] - rate_constants[5];
  jacobian[22] = jacobian[22] - rate_constants[6];
  jacobian[22] = jacobian[22] - rate_constants[2] * number_densities[6];
}

/* ---- helpers ------------------------------------------------------------------- */

/* The library's forcing for one cell with the given rate constants and number densities. */
static void library_force(const tmc_solver *solver, const double *rate_constants, const double *number_densities, double *force)
{
  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_dense_set_row(tmc_state_rate_constants(state), 0, rate_constants);
  tmc_dense *variables = tmc_dense_create(1, N_SPECIES, tmc_state_width(state), 0.0);
  tmc_dense_set_row(variables, 0, number_densities);
  tmc_dense *forcing = tmc_dense_create(1, N_SPECIES, tmc_state_width(state), 0.0);
  tmc_solver_add_forcing_terms(solver, state, variables, forcing);
  tmc_dense_get_row(forcing, 0, force);
  tmc_dense_free(forcing);
  tmc_dense_free(variables);
  tmc_state_free(state);
}

/* Sets one cell's state to the regression tests' initial conditions, each
 * concentration multiplied by `multiply` and divided by `divide` (MICM's
 * elem *= 10 and elem /= 10). */
static void set_regression_state(const tmc_solver *solver, tmc_state *state, double multiply, double divide)
{
  static const double initial[N_SPECIES] = { 1, 3.92e-1, 1.69e-2, 0, 3.29e1, 0, 0, 8.84, 0 };
  for (int i = 0; i < N_SPECIES; ++i)
    tmc_state_set(state, 0, i, initial[i] * multiply / divide);
  tmc_conditions *conditions = tmc_state_conditions(state);
  conditions[0].temperature = 273.15;
  conditions[0].pressure = 1000 * 100; /* 1000 hPa */
  conditions[0].air_density = 2.7e19;
  tmc_error err;
  tmc_error_clear(&err);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "O2_1", 1.0e-4, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "O3_1", 1.0e-5, &err), TMC_OK);
  TMC_EXPECT_EQ(tmc_state_set_custom_rate_parameter(state, "O3_2", 1.0e-6, &err), TMC_OK);
  tmc_state_set_relative_tolerance(state, 1.0e-4);
  tmc_solver_update_state_parameters(solver, state);
}

static tmc_result solve_one_second(tmc_solver *solver, tmc_state *state)
{
  const double time_start = 0, time_end = 1;
  return tmc_solver_solve(solver, time_end - time_start, state);
}

/* ---- cases ------------------------------------------------------------------- */

TMC_TEST(default_constructor)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  const tmc_rosenbrock_parameters *p = tmc_solver_rosenbrock_parameters(solver);
  TMC_EXPECT_EQ(p->stages, 3);
  /* the hard-coded solver's ThreeStageRosenbrock coefficients */
  TMC_EXPECT_EQ(p->a[0], 1.0);
  TMC_EXPECT_EQ(p->a[1], 1.0);
  TMC_EXPECT_EQ(p->a[2], 0.0);
  TMC_EXPECT_EQ(p->c[0], -0.10156171083877702091975600115545e+01);
  TMC_EXPECT_EQ(p->c[1], 0.40759956452537699824805835358067e+01);
  TMC_EXPECT_EQ(p->c[2], 0.92076794298330791242156818474003e+01);
  TMC_EXPECT_TRUE(p->new_function_evaluation[0] && p->new_function_evaluation[1] && !p->new_function_evaluation[2]);
  TMC_EXPECT_EQ(p->m[0], 0.1e+01);
  TMC_EXPECT_EQ(p->m[1], 0.61697947043828245592553615689730e+01);
  TMC_EXPECT_EQ(p->m[2], -0.42772256543218573326238373806514);
  TMC_EXPECT_EQ(p->e[0], 0.5);
  TMC_EXPECT_EQ(p->e[1], -0.29079558716805469821718236208017e+01);
  TMC_EXPECT_EQ(p->e[2], 0.22354069897811569627360909276199);
  TMC_EXPECT_EQ(p->estimator_of_local_order, 3.0);
  TMC_EXPECT_EQ(p->alpha[0], 0.0);
  TMC_EXPECT_EQ(p->alpha[1], 0.43586652150845899941601945119356);
  TMC_EXPECT_EQ(p->alpha[2], 0.43586652150845899941601945119356);
  TMC_EXPECT_EQ(p->gamma[0], 0.43586652150845899941601945119356);
  TMC_EXPECT_EQ(p->gamma[1], 0.24291996454816804366592249683314);
  TMC_EXPECT_EQ(p->gamma[2], 0.21851380027664058511513169485832e+01);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

/* The library's LU matrix of the Chapman mechanism with the given values in
 * the 23 hard-coded elements and zero in the four others. */
static void set_chapman_elements(tmc_sparse *matrix, const double *values)
{
  tmc_sparse_fill(matrix, 0.0);
  for (int e = 0; e < N_ELEMENTS; ++e)
  {
    TMC_EXPECT_FALSE(tmc_sparse_is_zero(matrix, elements[e][0], elements[e][1]));
    tmc_sparse_set(matrix, 0, elements[e][0], elements[e][1], values[e]);
  }
}

TMC_TEST(lin_solve)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  const tmc_lu *lu = tmc_solver_lu(solver);
  TMC_EXPECT_EQ(tmc_solver_lu_pattern(solver)->nnz, N_ELEMENTS + 4);
  /* a factored matrix of ones (unit pivots, so the reciprocal convention does not matter) */
  double ones[N_ELEMENTS];
  for (int e = 0; e < N_ELEMENTS; ++e)
    ones[e] = 1.0;
  tmc_sparse *jacobian = tmc_lu_matrix(lu, 1, 0.0);
  set_chapman_elements(jacobian, ones);
  double reciprocals[N_SPECIES];
  for (int i = 0; i < N_SPECIES; ++i)
    reciprocals[i] = 1.0;
  static const double expected[N_SPECIES] = { 0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -1, 0.5, 0 };
  for (int pass = 0; pass < 2; ++pass)
  {
    tmc_dense *x = tmc_dense_create(1, N_SPECIES, 1, 0.5);
    tmc_lu_solve_in_place(lu, jacobian, pass == 0 ? NULL : reciprocals, x);
    for (int i = 0; i < N_SPECIES; ++i)
      TMC_EXPECT_EQ(tmc_dense_get(x, 0, i), expected[i]);
    tmc_dense_free(x);
  }
  tmc_sparse_free(jacobian);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

TMC_TEST(reaction_names)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  TMC_ASSERT_TRUE(tmc_solver_num_reactions(solver) == 7);
  expect_reaction_order(solver);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

TMC_TEST(photolysis_names_case)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  TMC_ASSERT_TRUE(tmc_solver_num_custom_parameters(solver) == 3);
  for (int i = 0; i < 3; ++i)
    TMC_EXPECT_STREQ(tmc_solver_custom_parameter_label(solver, i), photolysis_names[i]);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

TMC_TEST(species_names_case)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  TMC_ASSERT_TRUE(tmc_solver_num_variables(solver) == N_SPECIES);
  for (int i = 0; i < N_SPECIES; ++i)
    TMC_EXPECT_STREQ(tmc_solver_variable_name(solver, i), species_names[i]);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

TMC_TEST(simple_force)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  double rate_constants[9], number_densities[N_SPECIES];
  for (int i = 0; i < 9; ++i)
    rate_constants[i] = 1;
  for (int i = 0; i < N_SPECIES; ++i)
    number_densities[i] = 1;
  double forcing[N_SPECIES];
  library_force(solver, rate_constants, number_densities, forcing);
  /* the truth values were calculated in Fortran with old MICM */
  static const double expected[N_SPECIES] = { 0, 0, 0, 0, 0, -1, 3, 2, -2 };
  for (int i = 0; i < N_SPECIES; ++i)
    TMC_EXPECT_EQ(forcing[i], expected[i]);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

TMC_TEST(smaller_force)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  double rate_constants[9], number_densities[N_SPECIES];
  for (int i = 0; i < 9; ++i)
    rate_constants[i] = 3e-8;
  for (int i = 0; i < N_SPECIES; ++i)
    number_densities[i] = 5e-6;
  double forcing[N_SPECIES], reference[N_SPECIES];
  library_force(solver, rate_constants, number_densities, forcing);
  /* the truth values were calculated in Fortran with old MICM */
  for (int i = 0; i < 5; ++i)
    TMC_EXPECT_EQ(forcing[i], 0.0);
  TMC_EXPECT_NEAR(forcing[5], 1.49e-13, 0.01);
  TMC_EXPECT_NEAR(forcing[6], 4.55e-13, 0.01);
  TMC_EXPECT_NEAR(forcing[7], 1.5e-13, 0.01);
  TMC_EXPECT_NEAR(forcing[8], -3e-13, 0.01);
  /* and, tighter, the hard-coded solver's forcing */
  reference_force(rate_constants, number_densities, reference);
  for (int i = 0; i < N_SPECIES; ++i)
    TMC_EXPECT_REL(forcing[i], reference[i], 1e-14);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

TMC_TEST(factored_alpha_minus_jac)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  tmc_state *state = tmc_solver_create_state(solver, 1);
  tmc_sparse *jacobian = tmc_state_jacobian(state);
  /* dforce_dy = 1 in every element; the library's matrix holds -J */
  double minus_dforce_dy[N_ELEMENTS];
  for (int e = 0; e < N_ELEMENTS; ++e)
    minus_dforce_dy[e] = -1.0;
  set_chapman_elements(jacobian, minus_dforce_dy);
  const double alpha = 2;
  tmc_solver_alpha_minus_jacobian(solver, state, alpha);
  double reciprocals[N_SPECIES];
  tmc_lu_decompose_in_place(tmc_solver_lu(solver), jacobian, reciprocals);

  /* the truth values were calculated in Fortran with old MICM; diagonal
   * elements hold pivot reciprocals */
  static const double expected[N_ELEMENTS] = { 1.000,  -1.000, -1.000, -1.000, 1.000,  1.000,  1.000,  1.000,
                                               -1.000, -1.000, 1.000,  -1.000, 1.000,  -1.000, -1.000, -1.000,
                                               -2.000, -1.000, 3.000,  -1.000, -2.000, -3.000, 0.125 };
  for (int e = 0; e < N_ELEMENTS; ++e)
  {
    const int row = elements[e][0], col = elements[e][1];
    double value = tmc_sparse_get(jacobian, 0, row, col);
    if (row == col)
    {
      TMC_EXPECT_EQ(reciprocals[row], 1.0 / value);
      value = reciprocals[row];
    }
    TMC_EXPECT_NEAR(value, expected[e], 1e-14); /* MICM: within 0.01 */
  }
  tmc_state_free(state);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

TMC_TEST(dforce_dy_times_vector)
{
  /* MICM multiplies a Jacobian of ones by a vector of halves with the
   * hard-coded solver's product, which takes the derivatives that are not
   * identically zero. The library's Jacobian must have those non-zeros and no
   * others, and the hard-coded solver's values. */
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  expect_reaction_order(solver);
  tmc_state *state = tmc_solver_create_state(solver, 1);
  set_regression_state(solver, state, 1.0, 1.0);
  /* every species present, so no derivative vanishes by accident */
  static const double y[N_SPECIES] = { 1.0, 0.392, 0.0169, 0.5, 32.9, 2.0e-3, 1.5e-3, 8.84, 3.0e-2 };
  tmc_dense *variables = tmc_dense_create(1, N_SPECIES, 1, 0.0);
  tmc_dense_set_row(variables, 0, y);
  tmc_sparse *minus_jacobian = tmc_state_jacobian(state);
  tmc_sparse_fill(minus_jacobian, 0.0);
  tmc_solver_subtract_jacobian_terms(solver, state, variables, minus_jacobian);

  double rate_constants[7], reference[N_ELEMENTS];
  tmc_dense_get_row(tmc_state_rate_constants(state), 0, rate_constants);
  reference_dforce_dy(rate_constants, y, reference);
  bool hard_coded[N_SPECIES][N_SPECIES] = { { false } };
  for (int e = 0; e < N_ELEMENTS; ++e)
  {
    hard_coded[elements[e][0]][elements[e][1]] = true;
    TMC_EXPECT_REL(-tmc_sparse_get(minus_jacobian, 0, elements[e][0], elements[e][1]), reference[e], 1e-14);
  }

  /* J times a vector of halves, with J one wherever the library's derivative is non-zero */
  double product[N_ELEMENTS] = { 0 };
  for (int row = 0; row < N_SPECIES; ++row)
    for (int col = 0; col < N_SPECIES; ++col)
    {
      if (tmc_sparse_is_zero(minus_jacobian, row, col))
        continue;
      const double value = tmc_sparse_get(minus_jacobian, 0, row, col);
      if (!hard_coded[row][col])
        TMC_EXPECT_EQ(value, 0.0); /* the elements the hard-coded solver leaves out are zero */
      if (value != 0.0)
        product[row] += 1.0 * 0.5;
    }
  /* the truth values were calculated in Fortran with old MICM */
  static const double expected[N_ELEMENTS] = { 0, 0, 0, 0, 0, 2, 3, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  for (int i = 0; i < N_ELEMENTS; ++i)
    TMC_EXPECT_NEAR(product[i], expected[i], 0.01);
  tmc_dense_free(variables);
  tmc_state_free(state);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

static void check_solve(double multiply, double divide, const double *expected)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_EXPECT_TRUE(solver != NULL);
  if (!solver)
    return;
  expect_reaction_order(solver);
  tmc_state *state = tmc_solver_create_state(solver, 1);
  set_regression_state(solver, state, multiply, divide);
  tmc_result results = solve_one_second(solver, state);
  /* MICM checks the state in Solve only; the hard-coded solver reports Converged unconditionally */
  TMC_EXPECT_EQ(results.state, TMC_CONVERGED);
  TMC_EXPECT_DOUBLE_EQ(results.stats.final_time, 1.0);
  for (int i = 0; i < N_SPECIES; ++i)
    TMC_EXPECT_NEAR(tmc_state_get(state, 0, i), expected[i], ABSOLUTE_TOLERANCE);
  tmc_state_free(state);
  tmc_solver_free(solver);
}

TMC_TEST(solve)
{
  static const double expected[N_SPECIES] = { 1, 0.392, 0.0169, 0, 32.9, 1.8039e-41, 0.00176789, 8.83912, 4.5031e-36 };
  check_solve(1.0, 1.0, expected);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_10_times_larger)
{
  static const double expected[N_SPECIES] = { 10, 3.92, 0.169, 0, 329, 1.8039e-38, 0.0176789, 88.3912, 4.5031e-33 };
  check_solve(10.0, 1.0, expected);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_10_times_smaller)
{
  static const double expected[N_SPECIES] = { 0.1, 0.0392, 0.00169, 0, 3.29, 1.8039e-44, 0.000176789, 0.883912, 4.5031e-39 };
  check_solve(1.0, 10.0, expected);
  return TMC_TEST_PASS;
}

TMC_TEST(solve_with_random_number_densities)
{
  tmc_solver *solver = chapman_solver(1);
  TMC_ASSERT_TRUE(solver != NULL);
  tmc_state *state = tmc_solver_create_state(solver, 1);
  set_regression_state(solver, state, 1.0, 1.0);
  /* std::default_random_engine with uniform_real_distribution<float>(0, 1) */
  tmc_minstd generator = tmc_minstd_create();
  for (int i = 0; i < N_SPECIES; ++i)
    tmc_state_set(state, 0, i, (double)tmc_minstd_uniform_float(&generator));
  tmc_solver_update_state_parameters(solver, state);
  tmc_result results = solve_one_second(solver, state);
  TMC_EXPECT_EQ(results.state, TMC_CONVERGED);
  tmc_state_free(state);
  tmc_solver_free(solver);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.DefaultConstructor", default_constructor),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.lin_solve", lin_solve),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.ReactionNames", reaction_names),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.PhotolysisNames", photolysis_names_case),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.SpeciesNames", species_names_case),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.simple_force", simple_force),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.smaller_force", smaller_force),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.factored_alpha_minus_jac", factored_alpha_minus_jac),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.dforce_dy_times_vector", dforce_dy_times_vector),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.Solve", solve),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.solve_10_times_larger", solve_10_times_larger),
  TMC_CASE("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.solve_with_random_number_densities", solve_with_random_number_densities),
  TMC_CASE("chapman_ode_solver:RegressionChapmanODESolver.solve_10_times_smaller", solve_10_times_smaller),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_chapman_ode_solver(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
