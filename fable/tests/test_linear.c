/* MICM test/unit/constraint/type/test_linear.cpp.
 *
 * Construction uses tmc_linear_constraint. The tests that MICM runs "through
 * a ConstraintSet" build the C equivalent (a solver over the variable map
 * with the constraint, cst_constraint_set; linear constraints have no state
 * parameters) and evaluate its residual and Jacobian with the set's
 * evaluators. The finite-difference checks are MICM's jacobian_verification.hpp
 * with its double-precision defaults. */
#include "constraint_support.h"

/* sum(coefficient species) = constant over the named terms */
static tmc_constraint linear(const char *name, const char *algebraic, const char *const *species, const double *coefficients, int n, double constant)
{
  tmc_stoich_species terms[4];
  for (int i = 0; i < n; ++i)
    terms[i] = cst_stoich(species[i], coefficients[i]);
  tmc_species a = tmc_species_create(algebraic);
  return cst_linear(name, &a, terms, n, constant);
}

static void expect_dependencies(const tmc_constraint *c, const char *const *expected, int n_expected)
{
  char names[8][TMC_NAME_MAX];
  int n = tmc_constraint_dependencies(c, names, 8);
  TMC_EXPECT_EQ(n, n_expected);
  for (int i = 0; i < n && i < n_expected; ++i)
    TMC_EXPECT_STREQ(names[i], expected[i]);
}

static const char *const AB[] = { "A", "B" };
static const char *const ABC[] = { "A", "B", "C" };

TMC_TEST(construction)
{
  /* A + B = 1: G = [A] + [B] - 1 */
  tmc_constraint c = linear("A_B_conservation", "B", AB, (const double[]){ 1.0, 1.0 }, 2, 1.0);
  TMC_EXPECT_STREQ(c.name, "A_B_conservation");
  TMC_EXPECT_EQ(c.constant, 1.0);
  expect_dependencies(&c, AB, 2);
  TMC_EXPECT_EQ(c.n_terms, 2);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(algebraic_species)
{
  tmc_constraint c = linear("A_B_C_conservation", "C", ABC, (const double[]){ 1.0, 1.0, 1.0 }, 3, 10.0);
  TMC_EXPECT_STREQ(c.algebraic_species.name, "C");
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(weighted_terms)
{
  /* 2A + 3B - C = 5 */
  tmc_constraint c = linear("weighted_sum", "C", ABC, (const double[]){ 2.0, 3.0, -1.0 }, 3, 5.0);
  TMC_EXPECT_STREQ(c.name, "weighted_sum");
  TMC_EXPECT_EQ(tmc_constraint_dependencies(&c, NULL, 0), 3);
  TMC_EXPECT_EQ(c.constant, 5.0);
  TMC_EXPECT_EQ(c.n_terms, 3);
  TMC_EXPECT_EQ(c.terms[0].coefficient, 2.0);
  TMC_EXPECT_EQ(c.terms[1].coefficient, 3.0);
  TMC_EXPECT_EQ(c.terms[2].coefficient, -1.0);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(zero_constant)
{
  /* A - B = 0 */
  tmc_constraint c = linear("A_equals_B", "B", AB, (const double[]){ 1.0, -1.0 }, 2, 0.0);
  TMC_EXPECT_STREQ(c.name, "A_equals_B");
  TMC_EXPECT_EQ(c.constant, 0.0);
  TMC_EXPECT_EQ(tmc_constraint_dependencies(&c, NULL, 0), 2);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(fractional_coefficients)
{
  /* 0.5 A + 1.5 B = 2 */
  tmc_constraint c = linear("fractional_conservation", "B", AB, (const double[]){ 0.5, 1.5 }, 2, 2.0);
  TMC_EXPECT_STREQ(c.name, "fractional_conservation");
  TMC_EXPECT_EQ(c.constant, 2.0);
  TMC_EXPECT_EQ(c.terms[0].coefficient, 0.5);
  TMC_EXPECT_EQ(c.terms[1].coefficient, 1.5);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

/* ---- through a constraint set -------------------------------------------- */

/* The constraint set over the species, its Jacobian (one block), and a
 * parameter matrix with no columns (linear constraints have none). */
typedef struct linear_set
{
  tmc_solver *s;
  tmc_sparse *jacobian;
  tmc_dense *state, *forcing, *parameters;
} linear_set;

static bool linear_set_create(linear_set *ls, const tmc_constraint *c, const char *const *species, int n_species, int n_blocks)
{
  memset(ls, 0, sizeof *ls);
  tmc_error err;
  ls->s = cst_constraint_set(species, n_species, c, 1, 1, &err);
  if (!ls->s)
  {
    tmc_test_fail(__FILE__, __LINE__, "build failed: %s", err.message);
    return false;
  }
  TMC_EXPECT_EQ(tmc_solver_num_custom_parameters(ls->s), 0);
  ls->jacobian = cst_jacobian(ls->s, n_blocks);
  ls->state = tmc_dense_create(n_blocks, n_species, 1, 0.0);
  ls->forcing = tmc_dense_create(n_blocks, n_species, 1, 0.0);
  ls->parameters = tmc_dense_create(n_blocks, 0, 1, 0.0);
  return true;
}

static void linear_set_free(linear_set *ls)
{
  tmc_dense_free(ls->state);
  tmc_dense_free(ls->forcing);
  tmc_dense_free(ls->parameters);
  tmc_sparse_free(ls->jacobian);
  tmc_solver_free(ls->s);
}

/* the residual of the given state in cell 0, algebraic row `row` */
static double residual(linear_set *ls, const double *values, int row)
{
  tmc_dense_set_row(ls->state, 0, values);
  tmc_dense_fill(ls->forcing, 0.0);
  cst_add_forcing_terms(ls->s, ls->state, ls->parameters, ls->forcing);
  return tmc_dense_get(ls->forcing, 0, row);
}

TMC_TEST(residual_computation_through_constraint_set)
{
  tmc_constraint c = linear("A_B_conservation", "B", AB, (const double[]){ 1.0, 1.0 }, 2, 1.0);
  linear_set ls;
  if (linear_set_create(&ls, &c, AB, 2, 1))
  {
    TMC_EXPECT_NEAR(residual(&ls, (const double[]){ 0.3, 0.7 }, 1), 0.0, 1e-10); /* satisfied */
    TMC_EXPECT_NEAR(residual(&ls, (const double[]){ 0.5, 0.6 }, 1), 0.1, 1e-10); /* 0.5 + 0.6 - 1 */
  }
  linear_set_free(&ls);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(jacobian_computation_through_constraint_set)
{
  tmc_constraint c = linear("A_B_conservation", "B", AB, (const double[]){ 1.0, 1.0 }, 2, 1.0);
  linear_set ls;
  if (linear_set_create(&ls, &c, AB, 2, 1))
  {
    tmc_dense_set_row(ls.state, 0, (const double[]){ 0.3, 0.7 });
    cst_subtract_jacobian_terms(ls.s, ls.state, ls.parameters, ls.jacobian);
    TMC_EXPECT_NEAR(CST_JAC(ls.jacobian, 0, 1, 0), -1.0, 1e-10);
    TMC_EXPECT_NEAR(CST_JAC(ls.jacobian, 0, 1, 1), -1.0, 1e-10);
  }
  linear_set_free(&ls);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(weighted_sum_residual_and_jacobian)
{
  tmc_constraint c = linear("weighted_sum", "C", ABC, (const double[]){ 2.0, 3.0, -1.0 }, 3, 5.0);
  linear_set ls;
  if (linear_set_create(&ls, &c, ABC, 3, 1))
  {
    /* 2 + 6 - 3 - 5 = 0 */
    TMC_EXPECT_NEAR(residual(&ls, (const double[]){ 1.0, 2.0, 3.0 }, 2), 0.0, 1e-10);
    tmc_sparse_fill(ls.jacobian, 0.0);
    cst_subtract_jacobian_terms(ls.s, ls.state, ls.parameters, ls.jacobian);
    TMC_EXPECT_NEAR(CST_JAC(ls.jacobian, 0, 2, 0), -2.0, 1e-10);
    TMC_EXPECT_NEAR(CST_JAC(ls.jacobian, 0, 2, 1), -3.0, 1e-10);
    TMC_EXPECT_NEAR(CST_JAC(ls.jacobian, 0, 2, 2), 1.0, 1e-10);
  }
  linear_set_free(&ls);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(three_species_conservation_residual)
{
  tmc_constraint c = linear("ABC_total", "C", ABC, (const double[]){ 1.0, 1.0, 1.0 }, 3, 10.0);
  linear_set ls;
  if (linear_set_create(&ls, &c, ABC, 3, 1))
  {
    TMC_EXPECT_NEAR(residual(&ls, (const double[]){ 2.0, 3.0, 5.0 }, 2), 0.0, 1e-10);
    TMC_EXPECT_NEAR(residual(&ls, (const double[]){ 2.0, 3.0, 6.0 }, 2), 1.0, 1e-10);
  }
  linear_set_free(&ls);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(zero_constant_residual)
{
  tmc_constraint c = linear("A_equals_B", "B", AB, (const double[]){ 1.0, -1.0 }, 2, 0.0);
  linear_set ls;
  if (linear_set_create(&ls, &c, AB, 2, 1))
  {
    TMC_EXPECT_NEAR(residual(&ls, (const double[]){ 0.5, 0.5 }, 1), 0.0, 1e-10);
    TMC_EXPECT_NEAR(residual(&ls, (const double[]){ 0.6, 0.5 }, 1), 0.1, 1e-10);
  }
  linear_set_free(&ls);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(fractional_coefficients_residual_and_jacobian)
{
  tmc_constraint c = linear("fractional_conservation", "B", AB, (const double[]){ 0.5, 1.5 }, 2, 2.0);
  linear_set ls;
  if (linear_set_create(&ls, &c, AB, 2, 1))
  {
    TMC_EXPECT_NEAR(residual(&ls, (const double[]){ 1.0, 1.0 }, 1), 0.0, 1e-10);
    cst_subtract_jacobian_terms(ls.s, ls.state, ls.parameters, ls.jacobian);
    TMC_EXPECT_NEAR(CST_JAC(ls.jacobian, 0, 1, 0), -0.5, 1e-10);
    TMC_EXPECT_NEAR(CST_JAC(ls.jacobian, 0, 1, 1), -1.5, 1e-10);
  }
  linear_set_free(&ls);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(jacobian_independent_of_concentrations)
{
  tmc_constraint c = linear("A_B_sum", "B", AB, (const double[]){ 2.0, 3.0 }, 2, 1.0);
  linear_set ls;
  if (linear_set_create(&ls, &c, AB, 2, 1))
  {
    tmc_dense_set_row(ls.state, 0, (const double[]){ 0.1, 0.2 });
    cst_subtract_jacobian_terms(ls.s, ls.state, ls.parameters, ls.jacobian);
    const double jac1_dA = CST_JAC(ls.jacobian, 0, 1, 0), jac1_dB = CST_JAC(ls.jacobian, 0, 1, 1);
    tmc_sparse_fill(ls.jacobian, 0.0);
    tmc_dense_set_row(ls.state, 0, (const double[]){ 10.0, 20.0 });
    cst_subtract_jacobian_terms(ls.s, ls.state, ls.parameters, ls.jacobian);
    const double jac2_dA = CST_JAC(ls.jacobian, 0, 1, 0), jac2_dB = CST_JAC(ls.jacobian, 0, 1, 1);
    TMC_EXPECT_NEAR(jac1_dA, jac2_dA, 1e-10);
    TMC_EXPECT_NEAR(jac1_dB, jac2_dB, 1e-10);
    TMC_EXPECT_NEAR(jac1_dA, -2.0, 1e-10);
    TMC_EXPECT_NEAR(jac1_dB, -3.0, 1e-10);
  }
  linear_set_free(&ls);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

/* The analytical Jacobian against finite differences of the residual; with
 * check_sparsity, also no finite-difference entry outside the pattern. */
static void check_against_finite_differences(const tmc_constraint *c, const char *const *species, int n_species, const double *rows, int n_blocks, bool check_sparsity)
{
  tmc_dense *variables = NULL, *parameters = NULL, *fd = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set(species, n_species, c, 1, 1, &err);
  CST_REQUIRE_BUILT(s, err);
  jacobian = cst_jacobian(s, n_blocks);
  variables = tmc_dense_create(n_blocks, n_species, 1, 0.0);
  for (int b = 0; b < n_blocks; ++b)
    tmc_dense_set_row(variables, b, rows + b * n_species);
  parameters = tmc_dense_create(n_blocks, 0, 1, 0.0);
  cst_subtract_jacobian_terms(s, variables, parameters, jacobian);
  cst_forcing_context context = { s, parameters };
  fd = cst_finite_difference_jacobian(cst_constraint_forcing, &context, variables, n_species, CST_FD_PERTURBATION);
  cst_jacobian_comparison comparison = cst_compare_jacobian_to_finite_difference(jacobian, fd, n_species, CST_FD_ATOL, CST_FD_RTOL);
  if (!comparison.passed)
    tmc_test_fail(__FILE__, __LINE__, "Linear constraint Jacobian mismatch: block=%d row=%d col=%d analytical=%.17g fd=%.17g", comparison.worst_block, comparison.worst_row, comparison.worst_col, comparison.worst_analytical, comparison.worst_fd);
  if (check_sparsity)
  {
    cst_jacobian_comparison sparsity = cst_check_jacobian_sparsity_completeness(jacobian, fd, n_species, CST_FD_SPARSITY_THRESHOLD);
    if (!sparsity.passed)
      tmc_test_fail(__FILE__, __LINE__, "Missing sparsity at block=%d row=%d col=%d fd_value=%.17g", sparsity.worst_block, sparsity.worst_row, sparsity.worst_col, sparsity.worst_fd);
  }
done:
  tmc_dense_free(variables);
  tmc_dense_free(parameters);
  tmc_dense_free(fd);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
}

TMC_TEST(finite_difference_jacobian_simple_conservation)
{
  /* A + B = 1, B algebraic: dG/dA = dG/dB = 1 */
  tmc_constraint c = linear("conservation", "B", AB, (const double[]){ 1.0, 1.0 }, 2, 1.0);
  check_against_finite_differences(&c, AB, 2, (const double[]){ 0.3, 0.7, 0.8, 0.2 }, 2, true);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(finite_difference_jacobian_weighted_sum)
{
  /* 2A + 3B - C = 5, C algebraic (MICM checks the values only here) */
  tmc_constraint c = linear("weighted", "C", ABC, (const double[]){ 2.0, 3.0, -1.0 }, 3, 5.0);
  check_against_finite_differences(&c, ABC, 3, (const double[]){ 1.5, 0.8, 2.0 }, 1, false);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("linear:LinearConstraint.Construction", construction),
  TMC_CASE("linear:LinearConstraint.AlgebraicSpecies", algebraic_species),
  TMC_CASE("linear:LinearConstraint.WeightedTerms", weighted_terms),
  TMC_CASE("linear:LinearConstraint.ZeroConstant", zero_constant),
  TMC_CASE("linear:LinearConstraint.FractionalCoefficients", fractional_coefficients),
  TMC_CASE("linear:LinearConstraint.ResidualComputationThroughConstraintSet", residual_computation_through_constraint_set),
  TMC_CASE("linear:LinearConstraint.JacobianComputationThroughConstraintSet", jacobian_computation_through_constraint_set),
  TMC_CASE("linear:LinearConstraint.WeightedSumResidualAndJacobian", weighted_sum_residual_and_jacobian),
  TMC_CASE("linear:LinearConstraint.ThreeSpeciesConservationResidual", three_species_conservation_residual),
  TMC_CASE("linear:LinearConstraint.ZeroConstantResidual", zero_constant_residual),
  TMC_CASE("linear:LinearConstraint.FractionalCoefficientsResidualAndJacobian", fractional_coefficients_residual_and_jacobian),
  TMC_CASE("linear:LinearConstraint.JacobianIndependentOfConcentrations", jacobian_independent_of_concentrations),
  TMC_CASE("linear:LinearConstraint.FiniteDifferenceJacobianSimpleConservation", finite_difference_jacobian_simple_conservation),
  TMC_CASE("linear:LinearConstraint.FiniteDifferenceJacobianWeightedSum", finite_difference_jacobian_weighted_sum),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_linear(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
