/* MICM test/unit/constraint/type/test_equilibrium.cpp.
 *
 * Construction and validation use tmc_equilibrium_constraint, whose failures
 * return MICM's error category and code where MICM throws. The tests that
 * MICM runs "through a ConstraintSet" build the C equivalent (a solver over
 * the variable map with the constraint, cst_constraint_set) and evaluate its
 * residual and Jacobian with the set's evaluators, the state parameters
 * holding K_eq as in MICM. The finite-difference checks are MICM's
 * jacobian_verification.hpp with its double-precision defaults. */
#include "constraint_support.h"

/* A + B <-> AB (AB algebraic) */
static tmc_constraint a_b_ab(const char *name, double K_HLC_ref, double delta_H)
{
  tmc_species AB = tmc_species_create("AB");
  const tmc_stoich_species r[] = { cst_stoich("A", 1.0), cst_stoich("B", 1.0) };
  const tmc_stoich_species p[] = { { AB, 1.0 } };
  return cst_equilibrium(name, &AB, r, 2, p, 1, K_HLC_ref, delta_H);
}

/* 2A <-> B + C (B algebraic) */
static tmc_constraint dissociation(double K_HLC_ref, double delta_H)
{
  tmc_species B = tmc_species_create("B");
  const tmc_stoich_species r[] = { cst_stoich("A", 2.0) };
  const tmc_stoich_species p[] = { { B, 1.0 }, cst_stoich("C", 1.0) };
  return cst_equilibrium("dissociation", &B, r, 1, p, 2, K_HLC_ref, delta_H);
}

/* The species the constraint depends on, compared with the expected names. */
static void expect_dependencies(const tmc_constraint *c, const char *const *expected, int n_expected)
{
  char names[8][TMC_NAME_MAX];
  int n = tmc_constraint_dependencies(c, names, 8);
  TMC_EXPECT_EQ(n, n_expected);
  for (int i = 0; i < n && i < n_expected; ++i)
    TMC_EXPECT_STREQ(names[i], expected[i]);
}

TMC_TEST(construction)
{
  /* A + B <-> AB with K_eq = 1000: G = K_eq [A][B] - [AB] */
  tmc_constraint c = a_b_ab("A_B_equilibrium", 1000.0, -2400.0);
  TMC_EXPECT_STREQ(c.name, "A_B_equilibrium");
  expect_dependencies(&c, (const char *[]){ "A", "B", "AB" }, 3);
  TMC_EXPECT_EQ(c.n_reactants, 2);
  TMC_EXPECT_EQ(c.n_products, 1);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(algebraic_species)
{
  tmc_constraint c = a_b_ab("A_B_equilibrium", 1000.0, -2400.0);
  TMC_EXPECT_STREQ(c.algebraic_species.name, "AB");
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(single_reactant_single_product)
{
  /* A <-> B with K_eq = 10: G = K_eq [A] - [B] */
  tmc_species B = tmc_species_create("B");
  const tmc_stoich_species r[] = { cst_stoich("A", 1.0) }, p[] = { { B, 1.0 } };
  tmc_constraint c = cst_equilibrium("A_B_simple", &B, r, 1, p, 1, 10.0, -2400.0);
  TMC_EXPECT_STREQ(c.name, "A_B_simple");
  expect_dependencies(&c, (const char *[]){ "A", "B" }, 2);
  TMC_EXPECT_STREQ(c.algebraic_species.name, "B");
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(multiple_reactants_and_products)
{
  /* 2A <-> B + C with K_eq = 100: G = K_eq [A]^2 - [B][C] */
  tmc_constraint c = dissociation(100.0, -2400.0);
  TMC_EXPECT_STREQ(c.name, "dissociation");
  expect_dependencies(&c, (const char *[]){ "A", "B", "C" }, 3);
  TMC_EXPECT_EQ(c.n_reactants, 1);
  TMC_EXPECT_EQ(c.reactants[0].coefficient, 2.0);
  TMC_EXPECT_EQ(c.n_products, 2);
  TMC_EXPECT_STREQ(c.algebraic_species.name, "B");
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

/* tmc_equilibrium_constraint with one reactant and one product of the given
 * stoichiometries (or none, when a count is 0) must fail with the code. */
static void expect_invalid(double K_HLC_ref, int n_reactants, double reactant_stoich, int n_products, double product_stoich, int code)
{
  tmc_species B = tmc_species_create("B");
  const tmc_stoich_species r[] = { cst_stoich("A", reactant_stoich) }, p[] = { { B, product_stoich } };
  tmc_constraint c;
  tmc_error err;
  int rc = tmc_equilibrium_constraint(&c, "invalid", &B, r, n_reactants, p, n_products, (tmc_vant_hoff){ K_HLC_ref, -2400.0, 0.0, 0.0 }, &err);
  CST_EXPECT_ERROR(rc, err, TMC_ERROR_CONSTRAINT, code);
  if (rc == TMC_OK)
    tmc_constraint_free(&c);
}

TMC_TEST(invalid_equilibrium_constant)
{
  expect_invalid(-1.0, 1, 1.0, 1, 1.0, TMC_CONSTRAINT_INVALID_EQUILIBRIUM_CONSTANT);
  expect_invalid(0.0, 1, 1.0, 1, 1.0, TMC_CONSTRAINT_INVALID_EQUILIBRIUM_CONSTANT);
  return TMC_TEST_PASS;
}

TMC_TEST(empty_reactants_throws)
{
  expect_invalid(1.0, 0, 1.0, 1, 1.0, TMC_CONSTRAINT_EMPTY_REACTANTS);
  return TMC_TEST_PASS;
}

TMC_TEST(empty_products_throws)
{
  expect_invalid(1.0, 1, 1.0, 0, 1.0, TMC_CONSTRAINT_EMPTY_PRODUCTS);
  return TMC_TEST_PASS;
}

TMC_TEST(invalid_stoichiometry_throws)
{
  expect_invalid(1.0, 1, 0.0, 1, 1.0, TMC_CONSTRAINT_INVALID_STOICHIOMETRY);  /* zero, reactant */
  expect_invalid(1.0, 1, -1.0, 1, 1.0, TMC_CONSTRAINT_INVALID_STOICHIOMETRY); /* negative, reactant */
  expect_invalid(1.0, 1, 1.0, 1, 0.0, TMC_CONSTRAINT_INVALID_STOICHIOMETRY);  /* zero, product */
  expect_invalid(1.0, 1, 1.0, 1, -2.0, TMC_CONSTRAINT_INVALID_STOICHIOMETRY); /* negative, product */
  return TMC_TEST_PASS;
}

/* ---- through a constraint set -------------------------------------------- */

TMC_TEST(residual_computation_through_constraint_set)
{
  tmc_constraint c = a_b_ab("A_B_equilibrium", 1000.0, -2400.0);
  tmc_dense *state = NULL, *forcing = NULL, *parameters = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set((const char *[]){ "A", "B", "AB" }, 3, &c, 1, 1, &err);
  CST_REQUIRE_BUILT(s, err);
  jacobian = cst_jacobian(s, 1);
  CST_REQUIRE(cst_parameter_index(s, "A_B_equilibrium") >= 0);
  state = tmc_dense_create(1, 3, 1, 0.0);
  forcing = tmc_dense_create(1, 3, 1, 0.0);
  parameters = cst_parameters(s, 1, 1000.0);
  /* at equilibrium: K_eq [A][B] = 1000 x 0.001 x 0.001 = 0.001 = [AB] */
  tmc_dense_set_row(state, 0, (const double[]){ 0.001, 0.001, 0.001 });
  tmc_dense_fill(forcing, 0.0);
  cst_add_forcing_terms(s, state, parameters, forcing);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 2), 0.0, 1e-10);
  /* away from it: 1000 x 0.02 x 0.03 - 0.05 = 0.55 */
  tmc_dense_set_row(state, 0, (const double[]){ 0.02, 0.03, 0.05 });
  tmc_dense_fill(forcing, 0.0);
  cst_add_forcing_terms(s, state, parameters, forcing);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 2), 0.55, 1e-10);
done:
  tmc_dense_free(state);
  tmc_dense_free(forcing);
  tmc_dense_free(parameters);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(jacobian_computation_through_constraint_set)
{
  tmc_constraint c = a_b_ab("A_B_equilibrium", 1000.0, -2400.0);
  tmc_dense *state = NULL, *parameters = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set((const char *[]){ "A", "B", "AB" }, 3, &c, 1, 1, &err);
  CST_REQUIRE_BUILT(s, err);
  jacobian = cst_jacobian(s, 1);
  CST_REQUIRE(cst_parameter_index(s, "A_B_equilibrium") >= 0);
  state = tmc_dense_create(1, 3, 1, 0.0);
  tmc_dense_set_row(state, 0, (const double[]){ 0.01, 0.02, 0.05 });
  parameters = cst_parameters(s, 1, 1000.0);
  cst_subtract_jacobian_terms(s, state, parameters, jacobian);
  const double dG_dA = 1000.0 * 0.02, dG_dB = 1000.0 * 0.01, dG_dAB = -1.0;
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 0), -dG_dA, 1e-8);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 1), -dG_dB, 1e-8);
  TMC_EXPECT_NEAR(CST_JAC(jacobian, 0, 2, 2), -dG_dAB, 1e-8);
done:
  tmc_dense_free(state);
  tmc_dense_free(parameters);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(complex_stoichiometry_residual)
{
  tmc_constraint c = dissociation(100.0, -2400.0);
  tmc_dense *state = NULL, *forcing = NULL, *parameters = NULL;
  tmc_sparse *jacobian = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set((const char *[]){ "A", "B", "C" }, 3, &c, 1, 1, &err);
  CST_REQUIRE_BUILT(s, err);
  jacobian = cst_jacobian(s, 1);
  CST_REQUIRE(cst_parameter_index(s, "dissociation") >= 0);
  state = tmc_dense_create(1, 3, 1, 0.0);
  forcing = tmc_dense_create(1, 3, 1, 0.0);
  /* at equilibrium: K_eq [A]^2 = 100 x 0.01 = 1 = 0.5 x 2 = [B][C] */
  tmc_dense_set_row(state, 0, (const double[]){ 0.1, 0.5, 2.0 });
  parameters = cst_parameters(s, 1, 100.0);
  cst_add_forcing_terms(s, state, parameters, forcing);
  TMC_EXPECT_NEAR(tmc_dense_get(forcing, 0, 1), 0.0, 1e-10);
done:
  tmc_dense_free(state);
  tmc_dense_free(forcing);
  tmc_dense_free(parameters);
  tmc_sparse_free(jacobian);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

/* The analytical Jacobian against central differences of the residual, and
 * no finite-difference entry outside the sparsity pattern. */
static void check_against_finite_differences(const tmc_solver *s, const tmc_dense *variables, const tmc_dense *parameters, int n_species)
{
  tmc_sparse *jacobian = cst_jacobian(s, variables->n_rows);
  cst_subtract_jacobian_terms(s, variables, parameters, jacobian);
  cst_forcing_context context = { s, parameters };
  tmc_dense *fd = cst_finite_difference_jacobian(cst_constraint_forcing, &context, variables, n_species, CST_FD_PERTURBATION);
  cst_jacobian_comparison comparison = cst_compare_jacobian_to_finite_difference(jacobian, fd, n_species, CST_FD_ATOL, CST_FD_RTOL);
  if (!comparison.passed)
    tmc_test_fail(__FILE__, __LINE__, "Jacobian mismatch: block=%d row=%d col=%d analytical=%.17g fd=%.17g", comparison.worst_block, comparison.worst_row, comparison.worst_col, comparison.worst_analytical, comparison.worst_fd);
  cst_jacobian_comparison sparsity = cst_check_jacobian_sparsity_completeness(jacobian, fd, n_species, CST_FD_SPARSITY_THRESHOLD);
  if (!sparsity.passed)
    tmc_test_fail(__FILE__, __LINE__, "Missing sparsity at block=%d row=%d col=%d fd_value=%.17g", sparsity.worst_block, sparsity.worst_row, sparsity.worst_col, sparsity.worst_fd);
  tmc_dense_free(fd);
  tmc_sparse_free(jacobian);
}

TMC_TEST(finite_difference_jacobian_simple)
{
  /* A + B <-> AB, K_eq = 1000 (delta_H = 0) */
  tmc_constraint c = a_b_ab("eq", 1000.0, 0.0);
  tmc_dense *variables = NULL, *parameters = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set((const char *[]){ "A", "B", "AB" }, 3, &c, 1, 1, &err);
  CST_REQUIRE_BUILT(s, err);
  CST_REQUIRE(cst_parameter_index(s, "eq") >= 0);
  variables = tmc_dense_create(2, 3, 1, 0.0);
  tmc_dense_set_row(variables, 0, (const double[]){ 0.02, 0.03, 0.05 });
  tmc_dense_set_row(variables, 1, (const double[]){ 0.1, 0.2, 5.0 });
  parameters = cst_parameters(s, 2, 1000.0);
  check_against_finite_differences(s, variables, parameters, 3);
done:
  tmc_dense_free(variables);
  tmc_dense_free(parameters);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

TMC_TEST(finite_difference_jacobian_complex_stoichiometry)
{
  /* 2A <-> B + C, K_eq = 100: dG/dA = 2 K_eq [A], dG/dB = -[C], dG/dC = -[B] */
  tmc_constraint c = dissociation(100.0, 0.0);
  tmc_dense *variables = NULL, *parameters = NULL;
  tmc_error err;
  tmc_solver *s = cst_constraint_set((const char *[]){ "A", "B", "C" }, 3, &c, 1, 1, &err);
  CST_REQUIRE_BUILT(s, err);
  CST_REQUIRE(cst_parameter_index(s, "dissociation") >= 0);
  variables = tmc_dense_create(1, 3, 1, 0.0);
  tmc_dense_set_row(variables, 0, (const double[]){ 0.15, 0.8, 1.5 });
  parameters = cst_parameters(s, 1, 100.0);
  check_against_finite_differences(s, variables, parameters, 3);
done:
  tmc_dense_free(variables);
  tmc_dense_free(parameters);
  tmc_solver_free(s);
  tmc_constraint_free(&c);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("equilibrium:EquilibriumConstraint.Construction", construction),
  TMC_CASE("equilibrium:EquilibriumConstraint.AlgebraicSpecies", algebraic_species),
  TMC_CASE("equilibrium:EquilibriumConstraint.SingleReactantSingleProduct", single_reactant_single_product),
  TMC_CASE("equilibrium:EquilibriumConstraint.MultipleReactantsAndProducts", multiple_reactants_and_products),
  TMC_CASE("equilibrium:EquilibriumConstraint.InvalidEquilibriumConstant", invalid_equilibrium_constant),
  TMC_CASE("equilibrium:EquilibriumConstraint.EmptyReactantsThrows", empty_reactants_throws),
  TMC_CASE("equilibrium:EquilibriumConstraint.EmptyProductsThrows", empty_products_throws),
  TMC_CASE("equilibrium:EquilibriumConstraint.InvalidStoichiometryThrows", invalid_stoichiometry_throws),
  TMC_CASE("equilibrium:EquilibriumConstraint.ResidualComputationThroughConstraintSet", residual_computation_through_constraint_set),
  TMC_CASE("equilibrium:EquilibriumConstraint.JacobianComputationThroughConstraintSet", jacobian_computation_through_constraint_set),
  TMC_CASE("equilibrium:EquilibriumConstraint.ComplexStoichiometryResidual", complex_stoichiometry_residual),
  TMC_CASE("equilibrium:EquilibriumConstraint.FiniteDifferenceJacobianSimple", finite_difference_jacobian_simple),
  TMC_CASE("equilibrium:EquilibriumConstraint.FiniteDifferenceJacobianComplexStoichiometry", finite_difference_jacobian_complex_stoichiometry),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_equilibrium(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
