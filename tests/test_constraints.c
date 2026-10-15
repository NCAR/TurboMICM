/* Constraint quantities and bounds ported from the reference's unit cases. */
#include "test.h"
static int constraint_numerics(ck_constraint *constraint, size_t nspecies, const double *values,
                               double expected, double tolerance, int difference)
{
  const ck_species species[] = {{.name = "A"}, {.name = "B"}, {.name = "C"}};
  ck_mechanism mechanism = {.name = "constraints",
                            .nspecies = nspecies,
                            .species = species,
                            .temperature = 298.15,
                            .pressure = 101325,
                            .nconstraints = 1,
                            .constraints = constraint};
  ck_model *model;
  ck_state *state;
  CHECK(ck_model_create_with_order(&mechanism, 0, &model) == CK_OK);
  CHECK(ck_state_create(model, 1, &state) == CK_OK);
  CHECK(ck_update_rates(state) == CK_OK);
  CHECK(ck_model_nconstraints(model) == 1);
  CHECK(ck_is_algebraic(model, constraint->algebraic_species));
  memcpy(state->y, values, nspecies * sizeof(double));
  double f[3] = {0}, base[3] = {0}, plus[3] = {0}, minus[3] = {0};
  size_t nnz = ck_jacobian_size(model);
  double *jac = calloc(nnz, sizeof(double));
  CHECK(jac);
  ck_constraint_forcing(state, state->y, f);
  NEAR(f[constraint->algebraic_species], expected, tolerance);
  ck_constraint_jacobian(state, state->y, jac);
  if (difference)
  {
    memcpy(base, state->y, nspecies * sizeof(double));
    for (size_t j = 0; j < nspecies; ++j)
    {
      double h = 1e-6 * fmax(1, fabs(base[j]));
      state->y[j] = base[j] + h;
      ck_constraint_forcing(state, state->y, plus);
      state->y[j] = base[j] - h;
      ck_constraint_forcing(state, state->y, minus);
      state->y[j] = base[j];
      double derivative =
          -(plus[constraint->algebraic_species] - minus[constraint->algebraic_species]) / (2 * h);
      size_t off = ck_jacobian_index(model, constraint->algebraic_species, j);
      double analytic = off == SIZE_MAX ? 0 : jac[off];
      NEAR(analytic, derivative, 1e-7 + 1e-7 * fabs(derivative));
    }
  }
  else
  {
    for (size_t j = 0; j < constraint->nreactants; ++j)
    {
      ck_product term = constraint->reactants[j];
      double deriv = term.yield;
      if (constraint->type == CK_EQUILIBRIUM_CONSTRAINT)
      {
        deriv = constraint->equilibrium_ref * term.yield * pow(values[term.species], term.yield - 1);
        for (size_t k = 0; k < constraint->nreactants; ++k)
          if (k != j)
            deriv *= pow(values[constraint->reactants[k].species], constraint->reactants[k].yield);
      }
      NEAR(jac[ck_jacobian_index(model, constraint->algebraic_species, term.species)], -deriv, 1e-8);
    }
    for (size_t j = 0; j < constraint->nproducts; ++j)
    {
      ck_product term = constraint->products[j];
      double deriv = term.yield * pow(values[term.species], term.yield - 1);
      for (size_t k = 0; k < constraint->nproducts; ++k)
        if (k != j)
          deriv *= pow(values[constraint->products[k].species], constraint->products[k].yield);
      NEAR(jac[ck_jacobian_index(model, constraint->algebraic_species, term.species)], deriv, 1e-8);
    }
  }
  free(jac);
  ck_state_free(state);
  ck_model_free(model);
  return 0;
}
int test_constraint(const char *suite, const char *test)
{
  int eq = !strcmp(suite, "EquilibriumConstraint"), linear = !strcmp(suite, "LinearConstraint");
  if (!eq && !linear)
    return 2;
  ck_product reactants[] = {{0, 1}, {1, 1}, {2, 1}}, products[] = {{2, 1}, {1, 1}};
  ck_constraint p = {.type = eq ? CK_EQUILIBRIUM_CONSTRAINT : CK_LINEAR_CONSTRAINT,
                     .name = eq ? "A_B_equilibrium" : "A_B_conservation",
                     .algebraic_species = eq ? 2 : 1,
                     .nreactants = 2,
                     .reactants = reactants,
                     .nproducts = eq ? 1 : 0,
                     .products = products,
                     .constant = 1,
                     .equilibrium_ref = 1000,
                     .enthalpy = -2400};
  if (eq && (IS("InvalidEquilibriumConstant") || IS("EmptyReactantsThrows") || IS("EmptyProductsThrows") ||
             IS("InvalidStoichiometryThrows")))
  {
    if (IS("InvalidEquilibriumConstant"))
    {
      p.equilibrium_ref = -1;
      CHECK(ck_constraint_validate(&p) == CK_INVALID);
      p.equilibrium_ref = 0;
      CHECK(ck_constraint_validate(&p) == CK_INVALID);
    }
    if (IS("EmptyReactantsThrows"))
    {
      p.nreactants = 0;
      CHECK(ck_constraint_validate(&p) == CK_INVALID);
    }
    if (IS("EmptyProductsThrows"))
    {
      p.nproducts = 0;
      CHECK(ck_constraint_validate(&p) == CK_INVALID);
    }
    if (IS("InvalidStoichiometryThrows"))
    {
      reactants[0].yield = 0;
      CHECK(ck_constraint_validate(&p) == CK_INVALID);
      reactants[0].yield = -1;
      CHECK(ck_constraint_validate(&p) == CK_INVALID);
      reactants[0].yield = 1;
      products[0].yield = 0;
      CHECK(ck_constraint_validate(&p) == CK_INVALID);
      products[0].yield = -2;
      CHECK(ck_constraint_validate(&p) == CK_INVALID);
    }
    return 0;
  }
  if (IS("Construction"))
  {
    CHECK(ck_constraint_validate(&p) == CK_OK);
    CHECK(p.nreactants == 2);
    CHECK(p.reactants[0].species == 0 && p.reactants[1].species == 1);
    CHECK(p.nproducts == (eq ? 1u : 0u));
    CHECK(!strcmp(p.name, eq ? "A_B_equilibrium" : "A_B_conservation"));
    return 0;
  }
  if (IS("AlgebraicSpecies"))
  {
    p.algebraic_species = 2;
    CHECK(p.algebraic_species == 2);
    return 0;
  }
  if (eq && IS("SingleReactantSingleProduct"))
  {
    p.name = "A_B_simple";
    p.nreactants = 1;
    p.algebraic_species = 1;
    products[0].species = 1;
    CHECK(ck_constraint_validate(&p) == CK_OK);
    CHECK(p.nreactants + p.nproducts == 2);
    CHECK(p.reactants[0].species == 0 && p.products[0].species == 1);
    return 0;
  }
  int weighted = IS("WeightedTerms") || IS("WeightedSumResidualAndJacobian") ||
                 IS("FiniteDifferenceJacobianWeightedSum");
  int fractional = IS("FractionalCoefficients") || IS("FractionalCoefficientsResidualAndJacobian");
  int zero = IS("ZeroConstant") || IS("ZeroConstantResidual");
  int three = IS("ThreeSpeciesConservationResidual");
  int complex = eq && (IS("MultipleReactantsAndProducts") || IS("ComplexStoichiometryResidual") ||
                       IS("FiniteDifferenceJacobianComplexStoichiometry"));
  size_t n = eq ? 3 : 2;
  double values[] = {0.3, 0.7, 0}, expected = 0;
  if (weighted)
  {
    n = 3;
    p.algebraic_species = 2;
    p.nreactants = 3;
    p.constant = 5;
    reactants[0].yield = 2;
    reactants[1].yield = 3;
    reactants[2].yield = -1;
    values[0] = 1;
    values[1] = 2;
    values[2] = 3;
  }
  if (fractional)
  {
    reactants[0].yield = 0.5;
    reactants[1].yield = 1.5;
    p.constant = 2;
    values[0] = values[1] = 1;
  }
  if (zero)
  {
    p.constant = 0;
    reactants[1].yield = -1;
    values[0] = values[1] = 0.5;
  }
  if (three)
  {
    n = 3;
    p.algebraic_species = 2;
    p.nreactants = 3;
    p.constant = 10;
    values[0] = 2;
    values[1] = 3;
    values[2] = 5;
  }
  if (complex)
  {
    p.nreactants = 1;
    reactants[0].yield = 2;
    p.nproducts = 2;
    products[0].species = 1;
    products[1].species = 2;
    p.algebraic_species = 1;
    p.equilibrium_ref = 100;
    values[0] = 0.1;
    values[1] = 0.5;
    values[2] = 2;
  }
  if (eq && !complex)
  {
    values[0] = 0.001;
    values[1] = 0.001;
    values[2] = 0.001;
  }
  if (IS("WeightedTerms") || IS("FractionalCoefficients") || IS("ZeroConstant") ||
      IS("MultipleReactantsAndProducts"))
  {
    CHECK(ck_constraint_validate(&p) == CK_OK);
    CHECK(p.nreactants == (weighted ? 3u : complex ? 1u : 2u));
    if (weighted)
    {
      CHECK(p.constant == 5 && p.reactants[0].yield == 2 && p.reactants[1].yield == 3 &&
            p.reactants[2].yield == -1);
    }
    if (fractional)
    {
      CHECK(p.constant == 2 && p.reactants[0].yield == 0.5 && p.reactants[1].yield == 1.5);
    }
    if (zero)
      CHECK(p.constant == 0);
    if (complex)
    {
      CHECK(p.reactants[0].yield == 2 && p.nproducts == 2 && p.algebraic_species == 1);
    }
    return 0;
  }
  int finite_difference = strstr(test, "FiniteDifference") != NULL;
  if (IS("JacobianComputationThroughConstraintSet") && eq)
  {
    values[0] = 0.01;
    values[1] = 0.02;
    values[2] = 0.05;
    expected = 0.15;
  }
  if (IS("JacobianIndependentOfConcentrations"))
  {
    reactants[0].yield = 2;
    reactants[1].yield = 3;
    expected = 1.7;
  }
  if (!weighted && !fractional && !zero && !three && !complex && !finite_difference &&
      !IS("ResidualComputationThroughConstraintSet") && !IS("JacobianComputationThroughConstraintSet") &&
      !IS("JacobianIndependentOfConcentrations"))
    return 2;
  CHECK(constraint_numerics(&p, n, values, expected, 1e-10, finite_difference) == 0);
  if (IS("ResidualComputationThroughConstraintSet"))
  {
    if (eq)
    {
      values[0] = 0.02;
      values[1] = 0.03;
      values[2] = 0.05;
      expected = 0.55;
    }
    else
    {
      values[0] = 0.5;
      values[1] = 0.6;
      expected = 0.1;
    }
    CHECK(constraint_numerics(&p, n, values, expected, 1e-10, 0) == 0);
  }
  if (three)
  {
    values[2] = 6;
    CHECK(constraint_numerics(&p, n, values, 1, 1e-10, 0) == 0);
  }
  if (zero)
  {
    values[0] = 0.6;
    CHECK(constraint_numerics(&p, n, values, 0.1, 1e-10, 0) == 0);
  }
  if (IS("JacobianIndependentOfConcentrations"))
  {
    values[0] = 10;
    values[1] = 20;
    CHECK(constraint_numerics(&p, n, values, 79, 1e-10, 0) == 0);
  }
  return 0;
}
