#include "chem_constraints.h"
#include "chem.h"
#include "jacobian_verification.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
  return 1; \
} } while (0)
#define NEAR(actual, expected, tolerance) CHECK(fabs((actual) - (expected)) <= (tolerance))

#define GAS_CONSTANT 8.31446261815324
#define REFERENCE_TEMPERATURE 298.15

static chem_constraint equilibrium(
    const char *name,
    unsigned algebraic,
    const chem_constraint_term *reactants,
    size_t reactant_count,
    const chem_constraint_term *products,
    size_t product_count,
    double keq,
    double delta_h,
    size_t parameter_index) {
  return (chem_constraint){
    .name = name, .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = algebraic,
    .reactants = reactants, .reactant_count = reactant_count,
    .products = products, .product_count = product_count,
    .equilibrium_constant_ref = keq, .delta_enthalpy = delta_h,
    .gas_constant = GAS_CONSTANT, .reference_temperature = REFERENCE_TEMPERATURE,
    .parameter_index = parameter_index
  };
}

static chem_constraint linear(
    const char *name,
    unsigned algebraic,
    const chem_constraint_term *terms,
    size_t term_count,
    double constant) {
  return (chem_constraint){
    .name = name, .kind = CHEM_CONSTRAINT_LINEAR, .algebraic_species = algebraic,
    .terms = terms, .term_count = term_count, .constant = constant
  };
}

static chem_constraint_set constraint_set(
    size_t species_count, const chem_constraint *constraints, size_t count, size_t parameter_count) {
  return (chem_constraint_set){
    .species_count = species_count, .constraint_count = count,
    .state_parameter_count = parameter_count, .constraints = constraints
  };
}

static int equilibrium_construction(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint value = equilibrium("A_B_equilibrium", 2u, reactants, 2u, products, 1u, 1000.0, -2400.0, 0u);
  double keq = 0.0;
  CHECK(chem_constraint_validate(&value, 3u) == 0);
  CHECK(strcmp(value.name, "A_B_equilibrium") == 0);
  CHECK(value.kind == CHEM_CONSTRAINT_EQUILIBRIUM && value.algebraic_species == 2u);
  CHECK(value.reactant_count == 2u && value.product_count == 1u);
  CHECK(value.reactants[0].species == 0u && value.reactants[1].species == 1u && value.products[0].species == 2u);
  CHECK(chem_constraint_equilibrium_parameter(&value, REFERENCE_TEMPERATURE, &keq) == 0);
  NEAR(keq, 1000.0, 1e-12);
  return 0;
}

static int equilibrium_algebraic_species(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint value = equilibrium("A_B_equilibrium", 2u, reactants, 2u, products, 1u, 1000.0, -2400.0, 0u);
  CHECK(value.algebraic_species == 2u);
  return 0;
}

static int equilibrium_single_reactant(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 } };
  const chem_constraint_term products[] = { { 1u, 1.0 } };
  const chem_constraint value = equilibrium("A_B_simple", 1u, reactants, 1u, products, 1u, 10.0, -2400.0, 0u);
  CHECK(chem_constraint_validate(&value, 2u) == 0);
  CHECK(strcmp(value.name, "A_B_simple") == 0);
  CHECK(value.reactants[0].species == 0u && value.products[0].species == 1u);
  return 0;
}

static int equilibrium_multiple_terms(void) {
  const chem_constraint_term reactants[] = { { 0u, 2.0 } };
  const chem_constraint_term products[] = { { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_constraint value = equilibrium("dissociation", 1u, reactants, 1u, products, 2u, 100.0, -2400.0, 0u);
  CHECK(chem_constraint_validate(&value, 3u) == 0);
  CHECK(value.reactant_count == 1u && value.reactants[0].coefficient == 2.0);
  CHECK(value.product_count == 2u && value.products[0].species == 1u && value.products[1].species == 2u);
  return 0;
}

static int equilibrium_invalid_constant(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 } };
  const chem_constraint_term products[] = { { 1u, 1.0 } };
  chem_constraint bad = equilibrium("invalid", 1u, reactants, 1u, products, 1u, -1.0, -2400.0, 0u);
  CHECK(chem_constraint_validate(&bad, 2u) != 0);
  bad.equilibrium_constant_ref = 0.0;
  CHECK(chem_constraint_validate(&bad, 2u) != 0);
  return 0;
}

static int equilibrium_empty_reactants(void) {
  const chem_constraint_term products[] = { { 1u, 1.0 } };
  const chem_constraint bad = equilibrium("invalid", 1u, NULL, 0u, products, 1u, 1.0, -2400.0, 0u);
  CHECK(chem_constraint_validate(&bad, 2u) != 0);
  return 0;
}

static int equilibrium_empty_products(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 } };
  const chem_constraint bad = equilibrium("invalid", 1u, reactants, 1u, NULL, 0u, 1.0, -2400.0, 0u);
  CHECK(chem_constraint_validate(&bad, 2u) != 0);
  return 0;
}

static int equilibrium_invalid_stoichiometry(void) {
  const chem_constraint_term valid_r[] = { { 0u, 1.0 } };
  const chem_constraint_term valid_p[] = { { 1u, 1.0 } };
  chem_constraint bad = equilibrium("invalid", 1u, valid_r, 1u, valid_p, 1u, 1.0, -2400.0, 0u);
  chem_constraint_term terms[] = { { 0u, 0.0 } };
  bad.reactants = terms;
  CHECK(chem_constraint_validate(&bad, 2u) != 0);
  terms[0].coefficient = -1.0;
  CHECK(chem_constraint_validate(&bad, 2u) != 0);
  bad.reactants = valid_r;
  terms[0] = (chem_constraint_term){ 1u, 0.0 };
  bad.products = terms;
  CHECK(chem_constraint_validate(&bad, 2u) != 0);
  terms[0].coefficient = -2.0;
  CHECK(chem_constraint_validate(&bad, 2u) != 0);
  return 0;
}

static int equilibrium_residual_computation(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint c = equilibrium("A_B_equilibrium", 2u, reactants, 2u, products, 1u, 1000.0, -2400.0, 0u);
  const chem_constraint set_items[] = { c };
  const chem_constraint_set set = constraint_set(3u, set_items, 1u, 1u);
  const double parameters[] = { 1000.0 };
  double state[] = { 0.001, 0.001, 0.001 };
  double forcing[] = { 9.0, 8.0, 7.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, parameters, forcing) == 0);
  NEAR(forcing[2], 0.0, 1e-14);
  state[0] = 0.02; state[1] = 0.03; state[2] = 0.05;
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, parameters, forcing) == 0);
  NEAR(forcing[2], 0.55, 1e-14);
  CHECK(forcing[0] == 9.0 && forcing[1] == 8.0);
  return 0;
}

static int equilibrium_jacobian_computation(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint c = equilibrium("A_B_equilibrium", 2u, reactants, 2u, products, 1u, 1000.0, -2400.0, 0u);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(3u, items, 1u, 1u);
  const double state[] = { 0.01, 0.02, 0.05 }, parameters[] = { 1000.0 };
  double jacobian[9] = { 0.0 };
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, parameters, jacobian) == 0);
  NEAR(jacobian[6], -20.0, 1e-12);
  NEAR(jacobian[7], -10.0, 1e-12);
  NEAR(jacobian[8], 1.0, 1e-12);
  return 0;
}

static int equilibrium_complex_residual(void) {
  const chem_constraint_term reactants[] = { { 0u, 2.0 } };
  const chem_constraint_term products[] = { { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_constraint c = equilibrium("dissociation", 1u, reactants, 1u, products, 2u, 100.0, -2400.0, 0u);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(3u, items, 1u, 1u);
  const double state[] = { 0.1, 0.5, 2.0 }, parameters[] = { 100.0 };
  double forcing[] = { 0.0, 1.0, 0.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, parameters, forcing) == 0);
  NEAR(forcing[1], 0.0, 1e-12);
  return 0;
}

typedef struct {
  const chem_constraint_set *set;
  const double *parameters;
} constraint_fd_context;

static int constraint_forcing_callback(
    size_t block, const double *state, size_t species_count, double *forcing, void *opaque) {
  const constraint_fd_context *context = (const constraint_fd_context *)opaque;
  if (species_count != context->set->species_count) return -1;
  const double *parameters = context->parameters == NULL ? NULL :
      context->parameters + block * context->set->state_parameter_count;
  return chem_constraint_set_apply_forcing(context->set, 1u, state, parameters, forcing);
}

static int equilibrium_fd_simple(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint c = equilibrium("eq", 2u, reactants, 2u, products, 1u, 1000.0, 0.0, 0u);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(3u, items, 1u, 1u);
  const double state[] = { 0.02, 0.03, 0.05, 0.1, 0.2, 5.0 };
  const double parameters[] = { 1000.0, 1000.0 };
  double analytical[18] = { 0.0 }, finite_difference[18];
  unsigned char pattern[9];
  chem_jacobian_comparison comparison;
  constraint_fd_context context = { &set, parameters };
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 2u, state, parameters, analytical) == 0);
  CHECK(chem_finite_difference_jacobian(
      constraint_forcing_callback, &context, state, 2u, 3u, 1e-6, finite_difference) == 0);
  CHECK(chem_compare_jacobian_to_finite_difference(
      pattern, analytical, finite_difference, 2u, 3u, 1e-7, 1e-7, &comparison) == 0);
  CHECK(comparison.passed);
  CHECK(chem_check_jacobian_sparsity(pattern, finite_difference, 2u, 3u, 1e-8, &comparison) == 0);
  CHECK(comparison.passed);
  return 0;
}

static int equilibrium_fd_complex(void) {
  const chem_constraint_term reactants[] = { { 0u, 2.0 } };
  const chem_constraint_term products[] = { { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_constraint c = equilibrium("dissociation", 1u, reactants, 1u, products, 2u, 100.0, 0.0, 0u);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(3u, items, 1u, 1u);
  const double state[] = { 0.15, 0.8, 1.5 }, parameters[] = { 100.0 };
  double analytical[9] = { 0.0 }, finite_difference[9];
  unsigned char pattern[9];
  chem_jacobian_comparison comparison;
  constraint_fd_context context = { &set, parameters };
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, parameters, analytical) == 0);
  CHECK(chem_finite_difference_jacobian(
      constraint_forcing_callback, &context, state, 1u, 3u, 1e-6, finite_difference) == 0);
  CHECK(chem_compare_jacobian_to_finite_difference(
      pattern, analytical, finite_difference, 1u, 3u, 1e-7, 1e-7, &comparison) == 0);
  CHECK(comparison.passed);
  CHECK(chem_check_jacobian_sparsity(pattern, finite_difference, 1u, 3u, 1e-8, &comparison) == 0);
  CHECK(comparison.passed);
  return 0;
}

static int linear_construction(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint c = linear("A_B_conservation", 1u, terms, 2u, 1.0);
  CHECK(chem_constraint_validate(&c, 2u) == 0);
  CHECK(strcmp(c.name, "A_B_conservation") == 0 && c.constant == 1.0);
  CHECK(c.term_count == 2u && c.terms[0].species == 0u && c.terms[1].species == 1u);
  return 0;
}

static int linear_algebraic_species(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_constraint c = linear("A_B_C_conservation", 2u, terms, 3u, 10.0);
  CHECK(c.algebraic_species == 2u);
  return 0;
}

static int linear_weighted_terms(void) {
  const chem_constraint_term terms[] = { { 0u, 2.0 }, { 1u, 3.0 }, { 2u, -1.0 } };
  const chem_constraint c = linear("weighted_sum", 2u, terms, 3u, 5.0);
  CHECK(chem_constraint_validate(&c, 3u) == 0);
  CHECK(strcmp(c.name, "weighted_sum") == 0 && c.constant == 5.0 && c.term_count == 3u);
  CHECK(c.terms[0].coefficient == 2.0 && c.terms[1].coefficient == 3.0 && c.terms[2].coefficient == -1.0);
  return 0;
}

static int linear_zero_constant(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, -1.0 } };
  const chem_constraint c = linear("A_equals_B", 1u, terms, 2u, 0.0);
  CHECK(chem_constraint_validate(&c, 2u) == 0 && c.constant == 0.0 && c.term_count == 2u);
  return 0;
}

static int linear_fractional_coefficients(void) {
  const chem_constraint_term terms[] = { { 0u, 0.5 }, { 1u, 1.5 } };
  const chem_constraint c = linear("fractional_conservation", 1u, terms, 2u, 2.0);
  CHECK(chem_constraint_validate(&c, 2u) == 0);
  CHECK(c.constant == 2.0 && c.terms[0].coefficient == 0.5 && c.terms[1].coefficient == 1.5);
  return 0;
}

static int linear_residual(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint c = linear("A_B_conservation", 1u, terms, 2u, 1.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(2u, items, 1u, 0u);
  const double state[] = { 0.3, 0.7 };
  double forcing[] = { 0.0, 0.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, NULL, forcing) == 0);
  NEAR(forcing[1], 0.0, 1e-14);
  double other_state[] = { 0.5, 0.6 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, other_state, NULL, forcing) == 0);
  NEAR(forcing[1], 0.1, 1e-14);
  return 0;
}

static int linear_jacobian(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint c = linear("A_B_conservation", 1u, terms, 2u, 1.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(2u, items, 1u, 0u);
  const double state[] = { 0.3, 0.7 };
  double jacobian[4] = { 0.0 };
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, NULL, jacobian) == 0);
  NEAR(jacobian[2], -1.0, 1e-14);
  NEAR(jacobian[3], -1.0, 1e-14);
  return 0;
}

static int linear_weighted_residual_jacobian(void) {
  const chem_constraint_term terms[] = { { 0u, 2.0 }, { 1u, 3.0 }, { 2u, -1.0 } };
  const chem_constraint c = linear("weighted_sum", 2u, terms, 3u, 5.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(3u, items, 1u, 0u);
  const double state[] = { 1.0, 2.0, 3.0 };
  double forcing[3] = { 0.0 }, jacobian[9] = { 0.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, NULL, forcing) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, NULL, jacobian) == 0);
  NEAR(forcing[2], 0.0, 1e-14);
  NEAR(jacobian[6], -2.0, 1e-14);
  NEAR(jacobian[7], -3.0, 1e-14);
  NEAR(jacobian[8], 1.0, 1e-14);
  return 0;
}

static int linear_conservation_residual(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_constraint c = linear("ABC_total", 2u, terms, 3u, 10.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(3u, items, 1u, 0u);
  const double state[] = { 2.0, 3.0, 5.0 };
  double forcing[3] = { 0.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, NULL, forcing) == 0);
  NEAR(forcing[2], 0.0, 1e-14);
  const double violated[] = { 2.0, 3.0, 6.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, violated, NULL, forcing) == 0);
  NEAR(forcing[2], 1.0, 1e-14);
  return 0;
}

static int linear_zero_constant_residual(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, -1.0 } };
  const chem_constraint c = linear("A_equals_B", 1u, terms, 2u, 0.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(2u, items, 1u, 0u);
  const double state[] = { 0.5, 0.5 };
  double forcing[2] = { 0.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, NULL, forcing) == 0);
  NEAR(forcing[1], 0.0, 1e-14);
  const double violated[] = { 0.6, 0.5 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, violated, NULL, forcing) == 0);
  NEAR(forcing[1], 0.1, 1e-14);
  return 0;
}

static int linear_fractional_residual_jacobian(void) {
  const chem_constraint_term terms[] = { { 0u, 0.5 }, { 1u, 1.5 } };
  const chem_constraint c = linear("fractional_conservation", 1u, terms, 2u, 2.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(2u, items, 1u, 0u);
  const double state[] = { 1.0, 1.0 };
  double forcing[2] = { 0.0 }, jacobian[4] = { 0.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, NULL, forcing) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, NULL, jacobian) == 0);
  NEAR(forcing[1], 0.0, 1e-14);
  NEAR(jacobian[2], -0.5, 1e-14);
  NEAR(jacobian[3], -1.5, 1e-14);
  return 0;
}

static int linear_jacobian_independent(void) {
  const chem_constraint_term terms[] = { { 0u, 2.0 }, { 1u, 3.0 } };
  const chem_constraint c = linear("A_B_sum", 1u, terms, 2u, 1.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(2u, items, 1u, 0u);
  const double state1[] = { 0.1, 0.2 }, state2[] = { 10.0, 20.0 };
  double jacobian1[4] = { 0.0 }, jacobian2[4] = { 0.0 };
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state1, NULL, jacobian1) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state2, NULL, jacobian2) == 0);
  CHECK(jacobian1[2] == jacobian2[2] && jacobian1[3] == jacobian2[3]);
  CHECK(jacobian1[2] == -2.0 && jacobian1[3] == -3.0);
  return 0;
}

static int linear_fd_simple(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint c = linear("conservation", 1u, terms, 2u, 1.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(2u, items, 1u, 0u);
  const double state[] = { 0.3, 0.7, 0.8, 0.2 };
  double analytical[8] = { 0.0 }, finite_difference[8];
  unsigned char pattern[4];
  chem_jacobian_comparison comparison;
  const constraint_fd_context context = { &set, NULL };
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 2u, state, NULL, analytical) == 0);
  CHECK(chem_finite_difference_jacobian(
      constraint_forcing_callback, (void *)&context, state, 2u, 2u, 1e-6, finite_difference) == 0);
  CHECK(chem_compare_jacobian_to_finite_difference(
      pattern, analytical, finite_difference, 2u, 2u, 1e-8, 1e-8, &comparison) == 0 && comparison.passed);
  CHECK(chem_check_jacobian_sparsity(pattern, finite_difference, 2u, 2u, 1e-8, &comparison) == 0 && comparison.passed);
  return 0;
}

static int linear_fd_weighted(void) {
  const chem_constraint_term terms[] = { { 0u, 2.0 }, { 1u, 3.0 }, { 2u, -1.0 } };
  const chem_constraint c = linear("weighted", 2u, terms, 3u, 5.0);
  const chem_constraint items[] = { c };
  const chem_constraint_set set = constraint_set(3u, items, 1u, 0u);
  const double state[] = { 1.5, 0.8, 2.0 };
  double analytical[9] = { 0.0 }, finite_difference[9];
  unsigned char pattern[9];
  chem_jacobian_comparison comparison;
  const constraint_fd_context context = { &set, NULL };
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, NULL, analytical) == 0);
  CHECK(chem_finite_difference_jacobian(
      constraint_forcing_callback, (void *)&context, state, 1u, 3u, 1e-6, finite_difference) == 0);
  CHECK(chem_compare_jacobian_to_finite_difference(
      pattern, analytical, finite_difference, 1u, 3u, 1e-8, 1e-8, &comparison) == 0 && comparison.passed);
  return 0;
}

static int set_construction(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint item = equilibrium("A_B_eq", 2u, reactants, 2u, products, 1u, 3.3e-2, -24000.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  CHECK(chem_constraint_set_validate(&set) == 0);
  CHECK(chem_constraint_set_size(&set) == 1u);
  CHECK(chem_constraint_set_parameter_count(&set) == 1u);
  const chem_constraint duplicate_items[] = { item, item };
  const chem_constraint_set duplicate_set = constraint_set(3u, duplicate_items, 2u, 1u);
  CHECK(chem_constraint_set_validate(&duplicate_set) != 0);
  return 0;
}

static int set_algebraic_rows(void) {
  const chem_constraint_term reactants[] = { { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint item = equilibrium("B_C_eq", 2u, reactants, 1u, products, 1u, 3.3e-2, -24000.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  unsigned char differential[3], pattern[9];
  CHECK(chem_constraint_set_species_mask(&set, differential) == 0);
  CHECK(differential[0] == 1u && differential[1] == 1u && differential[2] == 0u);
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(pattern[2u * 3u + 1u] == 1u && pattern[2u * 3u + 2u] == 1u);
  CHECK(pattern[2u * 3u + 0u] == 0u);
  return 0;
}

static int set_nonzero_pattern(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint item = equilibrium("A_B_eq", 2u, reactants, 2u, products, 1u, 3.3e-2, -24000.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  unsigned char pattern[9];
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(pattern[6] == 1u && pattern[7] == 1u && pattern[8] == 1u);
  CHECK(pattern[0] == 0u && pattern[1] == 0u && pattern[2] == 0u);
  return 0;
}

static int set_multiple_constraints(void) {
  const chem_constraint_term r1[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term p1[] = { { 2u, 1.0 } };
  const chem_constraint_term r2[] = { { 3u, 1.0 } };
  const chem_constraint_term p2[] = { { 4u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "A_B_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 2u,
      .reactants = r1, .reactant_count = 2u, .products = p1, .product_count = 1u,
      .equilibrium_constant_ref = 3.3e-2, .delta_enthalpy = -24000.0, .gas_constant = GAS_CONSTANT,
      .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 0u },
    { .name = "C_D_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 4u,
      .reactants = r2, .reactant_count = 1u, .products = p2, .product_count = 1u,
      .equilibrium_constant_ref = 3.3e-2, .delta_enthalpy = -24000.0, .gas_constant = GAS_CONSTANT,
      .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 1u }
  };
  const chem_constraint_set set = constraint_set(5u, items, 2u, 2u);
  unsigned char pattern[25];
  CHECK(chem_constraint_set_validate(&set) == 0 && chem_constraint_set_size(&set) == 2u);
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(pattern[10u + 0u] && pattern[10u + 1u] && pattern[10u + 2u]);
  CHECK(pattern[20u + 3u] && pattern[20u + 4u]);
  size_t nonzero = 0u;
  for (size_t i = 0; i < 25u; ++i) nonzero += pattern[i] != 0u;
  CHECK(nonzero == 5u);
  return 0;
}

static int set_add_forcing(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint item = equilibrium("A_B_eq", 2u, reactants, 2u, products, 1u, 3.3e-2, -24000.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double state[] = { 0.2, 0.4, 0.6, 0.1, 0.3, 0.7 };
  const double parameters[] = { 3.3e-2, 3.3e-2 };
  double forcing[6] = { 0.0 };
  CHECK(chem_constraint_set_apply_forcing(&set, 2u, state, parameters, forcing) == 0);
  NEAR(forcing[2], -0.59736, 1e-12);
  NEAR(forcing[5], -0.69901, 1e-12);
  return 0;
}

static int set_subtract_jacobian(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint item = equilibrium("A_B_eq", 2u, reactants, 2u, products, 1u, 3.3e-2, -24000.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double state[] = { 0.01, 0.02, 0.05 }, parameters[] = { 3.3e-2 };
  double jacobian[9] = { 0.0 };
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, parameters, jacobian) == 0);
  NEAR(jacobian[6], -0.00066, 1e-14);
  NEAR(jacobian[7], -0.00033, 1e-14);
  NEAR(jacobian[8], 1.0, 1e-14);
  return 0;
}

static int set_empty(void) {
  const chem_constraint_set set = constraint_set(3u, NULL, 0u, 0u);
  const double state[] = { 0.1, 0.2, 0.3 };
  double forcing[] = { 1.0, 1.0, 1.0 };
  double jacobian[9];
  for (size_t i = 0; i < 9u; ++i) jacobian[i] = 2.0;
  unsigned char pattern[9], differential[3];
  CHECK(chem_constraint_set_size(&set) == 0u && chem_constraint_set_validate(&set) == 0);
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  for (size_t i = 0; i < 9u; ++i) CHECK(pattern[i] == 0u);
  CHECK(chem_constraint_set_species_mask(&set, differential) == 0);
  CHECK(differential[0] && differential[1] && differential[2]);
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, NULL, forcing) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, NULL, jacobian) == 0);
  CHECK(forcing[0] == 1.0 && forcing[1] == 1.0 && forcing[2] == 1.0);
  for (size_t i = 0; i < 9u; ++i) CHECK(jacobian[i] == 2.0);
  return 0;
}

static int set_unknown_species(void) {
  static const char *const names[] = { "A", "B" };
  unsigned index = 0u;
  CHECK(chem_constraint_find_species("A", names, 2u, &index) == 0 && index == 0u);
  CHECK(chem_constraint_find_species("X", names, 2u, &index) != 0);
  const chem_constraint_term invalid_terms[] = { { 2u, 1.0 } };
  const chem_constraint bad = linear("invalid", 0u, invalid_terms, 1u, 0.0);
  const chem_constraint_set set = constraint_set(2u, &bad, 1u, 0u);
  CHECK(chem_constraint_set_validate(&set) != 0);
  return 0;
}

static int set_three_species_one_constraint(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 } };
  const chem_constraint_term products[] = { { 1u, 1.0 } };
  const chem_constraint item = equilibrium("X_Y_eq", 1u, reactants, 1u, products, 1u, 3.3e-2, -24000.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double state[] = { 10.0, 0.2, 0.5, 10.0, 0.33, 0.3 };
  const double parameters[] = { 3.3e-2, 3.3e-2 };
  double forcing[6] = { 0.0 }, jacobian[18] = { 0.0 };
  unsigned char pattern[9];
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(pattern[3] && pattern[4] && !pattern[5]);
  CHECK(chem_constraint_set_apply_forcing(&set, 2u, state, parameters, forcing) == 0);
  NEAR(forcing[1], 0.13, 1e-14);
  NEAR(forcing[4], 0.0, 1e-14);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 2u, state, parameters, jacobian) == 0);
  NEAR(jacobian[3], -3.3e-2, 1e-14);
  NEAR(jacobian[4], 1.0, 1e-14);
  NEAR(jacobian[12], -3.3e-2, 1e-14);
  NEAR(jacobian[13], 1.0, 1e-14);
  return 0;
}

static int set_four_species_two_constraints(void) {
  const chem_constraint_term ar[] = { { 0u, 1.0 } }, ap[] = { { 1u, 1.0 } };
  const chem_constraint_term br[] = { { 2u, 1.0 }, { 3u, 1.0 } }, bp[] = { { 0u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "A_B_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 1u,
      .reactants = ar, .reactant_count = 1u, .products = ap, .product_count = 1u,
      .equilibrium_constant_ref = 3.3e-2, .delta_enthalpy = -24000.0, .gas_constant = GAS_CONSTANT,
      .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 0u },
    { .name = "CD_A_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 0u,
      .reactants = br, .reactant_count = 2u, .products = bp, .product_count = 1u,
      .equilibrium_constant_ref = 3.3e-2, .delta_enthalpy = -24000.0, .gas_constant = GAS_CONSTANT,
      .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 1u }
  };
  const chem_constraint_set set = constraint_set(4u, items, 2u, 2u);
  const double state[] = {
    0.33, 0.01089, 10.0, 1.0,
    0.33, 0.01089, 5.0, 1.0,
    0.33, 0.02, 5.0, 1.0
  };
  const double parameters[] = { 3.3e-2, 3.3e-2, 3.3e-2, 3.3e-2, 3.3e-2, 3.3e-2 };
  double forcing[12] = { 0.0 }, jacobian[48] = { 0.0 };
  unsigned char pattern[16];
  CHECK(chem_constraint_set_validate(&set) == 0);
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(pattern[4] && pattern[5] && pattern[2] && pattern[3] && pattern[0]);
  CHECK(chem_constraint_set_apply_forcing(&set, 3u, state, parameters, forcing) == 0);
  NEAR(forcing[0], 0.0, 1e-14); NEAR(forcing[1], 0.0, 1e-14);
  NEAR(forcing[4], -0.165, 1e-14); NEAR(forcing[5], 0.0, 1e-14);
  NEAR(forcing[8], -0.165, 1e-14); NEAR(forcing[9], -0.00911, 1e-14);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 3u, state, parameters, jacobian) == 0);
  NEAR(jacobian[4], -3.3e-2, 1e-14); NEAR(jacobian[5], 1.0, 1e-14);
  NEAR(jacobian[2], -3.3e-2, 1e-14); NEAR(jacobian[3], -0.33, 1e-14); NEAR(jacobian[0], 1.0, 1e-14);
  NEAR(jacobian[20], -3.3e-2, 1e-14); NEAR(jacobian[21], 1.0, 1e-14);
  NEAR(jacobian[18], -3.3e-2, 1e-14); NEAR(jacobian[19], -0.165, 1e-14); NEAR(jacobian[16], 1.0, 1e-14);
  return 0;
}

static int set_coupled_constraints(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 } };
  const chem_constraint_term b_products[] = { { 1u, 1.0 } }, c_products[] = { { 2u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "A_B_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 1u,
      .reactants = reactants, .reactant_count = 1u, .products = b_products, .product_count = 1u,
      .equilibrium_constant_ref = 3.3e-2, .delta_enthalpy = -24000.0, .gas_constant = GAS_CONSTANT,
      .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 0u },
    { .name = "A_C_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 2u,
      .reactants = reactants, .reactant_count = 1u, .products = c_products, .product_count = 1u,
      .equilibrium_constant_ref = 3.3e-2, .delta_enthalpy = -24000.0, .gas_constant = GAS_CONSTANT,
      .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 1u }
  };
  const chem_constraint_set set = constraint_set(3u, items, 2u, 2u);
  const double state[] = { 0.1, 0.0033, 0.0033 }, parameters[] = { 3.3e-2, 3.3e-2 };
  double forcing[3] = { 0.0 }, jacobian[9] = { 0.0 };
  unsigned char pattern[9];
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(pattern[3] && pattern[4] && pattern[6] && pattern[8]);
  size_t nnz = 0u; for (size_t i = 0; i < 9u; ++i) nnz += pattern[i] != 0u;
  CHECK(nnz == 4u);
  CHECK(chem_constraint_set_apply_forcing(&set, 1u, state, parameters, forcing) == 0);
  NEAR(forcing[1], 0.0, 1e-14); NEAR(forcing[2], 0.0, 1e-14);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, 1u, state, parameters, jacobian) == 0);
  NEAR(jacobian[3], -3.3e-2, 1e-14); NEAR(jacobian[4], 1.0, 1e-14);
  NEAR(jacobian[6], -3.3e-2, 1e-14); NEAR(jacobian[8], 1.0, 1e-14);
  return 0;
}

static int constraint_vector_cells(size_t cell_count) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_constraint_term products[] = { { 2u, 1.0 } };
  const chem_constraint item = equilibrium("A_B_eq", 2u, reactants, 2u, products, 1u, 3.3e-2, -24000.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  double state[9] = { 0.01, 0.02, 0.05, 0.03, 0.01, 0.2, 0.001, 0.002, 0.004 };
  double parameters[3] = { 3.3e-2, 3.3e-2, 3.3e-2 };
  double forcing[9] = { 0.0 }, jacobian[27] = { 0.0 };
  unsigned char pattern[9];
  CHECK(cell_count <= 3u);
  CHECK(chem_constraint_set_nonzero_pattern(&set, pattern) == 0);
  CHECK(pattern[6] && pattern[7] && pattern[8]);
  CHECK(chem_constraint_set_apply_forcing(&set, cell_count, state, parameters, forcing) == 0);
  CHECK(chem_constraint_set_accumulate_negative_jacobian(&set, cell_count, state, parameters, jacobian) == 0);
  for (size_t cell = 0; cell < cell_count; ++cell) {
    const size_t offset = cell * 3u;
    const double expected = 3.3e-2 * state[offset] * state[offset + 1u] - state[offset + 2u];
    NEAR(forcing[offset + 2u], expected, 1e-14);
    NEAR(jacobian[cell * 9u + 6u], -3.3e-2 * state[offset + 1u], 1e-14);
    NEAR(jacobian[cell * 9u + 7u], -3.3e-2 * state[offset], 1e-14);
    NEAR(jacobian[cell * 9u + 8u], 1.0, 1e-14);
  }
  return 0;
}

static int set_vectorized_indexing(void) {
  return constraint_vector_cells(3u);
}

static int constraint_vector_matrix_1(void) { return constraint_vector_cells(1u); }
static int constraint_vector_matrix_2(void) { return constraint_vector_cells(2u); }
static int constraint_vector_matrix_3(void) { return constraint_vector_cells(3u); }

static int initialize_with(
    const chem_constraint_set *set,
    size_t cells,
    const double *parameters,
    double *state,
    const double *absolute_tolerances,
    double relative_tolerance,
    size_t max_updates,
    double correction_tolerance,
    size_t max_backtracks,
    size_t *iterations) {
  return chem_constraint_set_initialize(
      set, cells, parameters, state, absolute_tolerances, relative_tolerance,
      max_updates, correction_tolerance, max_backtracks, 0.5, 1e-4, iterations);
}

static chem_constraint simple_equilibrium_constraint(size_t parameter_index) {
  static const chem_constraint_term reactants[] = { { 1u, 1.0 } };
  static const chem_constraint_term products[] = { { 2u, 1.0 } };
  return equilibrium("B_C_eq", 2u, reactants, 1u, products, 1u, 2.0, 0.0, parameter_index);
}

static int initialization_consistent(void) {
  const chem_constraint item = simple_equilibrium_constraint(0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double parameters[] = { 2.0 }, tolerances[] = { 1e-3, 1e-3, 1e-3 };
  double state[] = { 1.0, 0.5, 1.0 };
  size_t iterations = 0u;
  CHECK(initialize_with(&set, 1u, parameters, state, tolerances, 1e-4, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  CHECK(iterations == 1u && state[0] == 1.0 && state[1] == 0.5 && state[2] == 1.0);
  return 0;
}

static int initialization_mild(void) {
  const chem_constraint item = simple_equilibrium_constraint(0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double parameters[] = { 2.0 }, tolerances[] = { 1e-3, 1e-3, 1e-3 };
  double state[] = { 1.0, 0.5, 5.0 };
  size_t iterations = 0u;
  CHECK(initialize_with(&set, 1u, parameters, state, tolerances, 1e-4, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  CHECK(iterations > 1u && state[0] == 1.0 && state[1] == 0.5);
  NEAR(state[2], 1.0, 1e-12);
  return 0;
}

static int initialization_severe(void) {
  const chem_constraint item = simple_equilibrium_constraint(0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double parameters[] = { 2.0 }, tolerances[] = { 1e-3, 1e-3, 1e-3 };
  double state[] = { 1.0, 0.5, 1000.0 };
  size_t iterations = 0u;
  CHECK(initialize_with(&set, 1u, parameters, state, tolerances, 1e-4, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  NEAR(state[2], 1.0, 1e-12);
  return 0;
}

static int initialization_scaled_solve(void) {
  const double scales[] = { 1.0, 1e-12, 1e12 };
  const double tolerances[] = { 1e-3, 1e-3, 1e-3 };
  for (size_t i = 0; i < 3u; ++i) {
    const double scale = scales[i];
    const chem_constraint_term terms[] = { { 0u, scale }, { 1u, scale }, { 2u, scale } };
    const chem_constraint item = linear("mass", 2u, terms, 3u, scale * 1.3);
    const chem_constraint_set set = constraint_set(3u, &item, 1u, 0u);
    double state[] = { 0.7, 0.15, 5.0 };
    size_t iterations = 0u;
    CHECK(initialize_with(&set, 1u, NULL, state, tolerances, 1e-4, 10u, 0.1, 24u, &iterations) ==
        CHEM_CONSTRAINT_INIT_CONVERGED);
    NEAR(state[2], 0.45, 1e-12);
  }
  return 0;
}

static int initialization_empty_set(void) {
  const chem_constraint_set set = constraint_set(2u, NULL, 0u, 0u);
  double state[] = { 1.0, 2.0 };
  size_t iterations = 9u;
  CHECK(initialize_with(&set, 1u, NULL, state, NULL, 1e-4, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  CHECK(iterations == 0u && state[0] == 1.0 && state[1] == 2.0);
  return 0;
}

static int initialization_multicell(void) {
  const chem_constraint item = simple_equilibrium_constraint(0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double parameters[] = { 2.0, 2.0, 2.0 }, tolerances[] = { 1e-3, 1e-3, 1e-3 };
  double state[] = { 1.0, 0.5, 1.0, 2.0, 0.3, 5.0, 0.5, 0.1, 100.0 };
  size_t iterations = 0u;
  CHECK(initialize_with(&set, 3u, parameters, state, tolerances, 1e-4, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  NEAR(state[2], 1.0, 1e-12); NEAR(state[5], 0.6, 1e-12); NEAR(state[8], 0.2, 1e-12);
  return 0;
}

static int initialization_repeated(void) {
  const chem_constraint item = simple_equilibrium_constraint(0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double parameters[] = { 2.0 }, tolerances[] = { 1e-3, 1e-3, 1e-3 };
  double state[] = { 1.0, 0.5, 1.0 };
  size_t iterations = 0u;
  CHECK(initialize_with(&set, 1u, parameters, state, tolerances, 1e-4, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  state[2] = 999.0;
  CHECK(initialize_with(&set, 1u, parameters, state, tolerances, 1e-4, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  NEAR(state[2], 1.0, 1e-12);
  return 0;
}

static int initialization_status_name(void) {
  CHECK(strcmp(chem_constraint_init_status_name(CHEM_CONSTRAINT_INIT_FAILED),
      "Constraint Initialization Failed") == 0);
  return 0;
}

static int initialization_failure_restores_state(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 } };
  const chem_constraint_term products[] = { { 1u, 2.0 } };
  const chem_constraint item = equilibrium("squared", 1u, reactants, 1u, products, 1u, 2.0, 0.0, 0u);
  const chem_constraint_set set = constraint_set(2u, &item, 1u, 1u);
  const double parameters[] = { 2.0 }, tolerances[] = { 1e-3, 1e-3 };
  double state[] = { 0.5, 100.0 };
  size_t iterations = 0u;
  CHECK(initialize_with(&set, 1u, parameters, state, tolerances, 1e-4, 1u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_FAILED);
  CHECK(state[0] == 0.5 && state[1] == 100.0);
  return 0;
}

static int initialization_damped_success_and_limits(void) {
  const chem_constraint_term reactants[] = { { 0u, 1.0 } };
  const chem_constraint_term products[] = { { 1u, 2.0 } };
  const chem_constraint item = equilibrium("squared", 1u, reactants, 1u, products, 1u, 2.0, 0.0, 0u);
  const chem_constraint_set set = constraint_set(2u, &item, 1u, 1u);
  const double parameters[] = { 2.0 }, tolerances[] = { 1e-3, 1e-3 };
  double state[] = { 0.5, 1e-6 };
  size_t iterations = 0u;
  CHECK(initialize_with(&set, 1u, parameters, state, tolerances, 1e-4, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  NEAR(state[1], 1.0, 1e-3);
  for (size_t backtracks = 0u; backtracks <= 2u; backtracks += 2u) {
    state[1] = 1e-6;
    CHECK(initialize_with(&set, 1u, parameters, state, tolerances, 1e-4, 10u, 0.1, backtracks, &iterations) ==
        CHEM_CONSTRAINT_INIT_FAILED);
    CHECK(state[0] == 0.5 && state[1] == 1e-6);
  }
  return 0;
}

static int initialization_weighted_row_scale(void) {
  const double scales[] = { 1e-12, 1.0, 1e12 };
  const double tolerances[] = { 1e-12, 1e-12 };
  for (size_t i = 0; i < 3u; ++i) {
    const double scale = scales[i];
    const chem_constraint_term terms[] = { { 0u, scale }, { 1u, -scale } };
    const chem_constraint item = linear("scaled_copy", 1u, terms, 2u, 0.0);
    const chem_constraint_set set = constraint_set(2u, &item, 1u, 0u);
    double state[] = { 1.0, 4.0 };
    size_t iterations = 0u;
    CHECK(initialize_with(&set, 1u, NULL, state, tolerances, 1e-8, 10u, 0.1, 24u, &iterations) ==
        CHEM_CONSTRAINT_INIT_CONVERGED);
    NEAR(state[1], 1.0, 1e-12);
  }
  return 0;
}

static int initialization_weighted_tolerances(void) {
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, -1.0 } };
  const chem_constraint item = linear("scaled_copy", 1u, terms, 2u, 0.0);
  const chem_constraint_set set = constraint_set(2u, &item, 1u, 0u);
  const double loose_tolerance[] = { 1e-8, 1e-8 };
  const double tight_tolerance[] = { 1e-12, 1e-12 };
  double state[] = { 1.0, 1.0 + 5e-8 };
  size_t iterations = 0u;
  CHECK(initialize_with(&set, 1u, NULL, state, loose_tolerance, 1e-6, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  CHECK(iterations == 1u && state[1] == 1.0 + 5e-8);
  state[1] = 1.0 + 5e-8;
  CHECK(initialize_with(&set, 1u, NULL, state, tight_tolerance, 1e-8, 10u, 0.1, 24u, &iterations) ==
      CHEM_CONSTRAINT_INIT_CONVERGED);
  NEAR(state[1], 1.0, 1e-12);
  return 0;
}

static int dae_equilibrium_simple(void) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "B", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-3, 0u },
  };
  const chem_term reactants[] = { { 0u, 1.0 } };
  const chem_term products[] = { { 1u, 1.0 } };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .label = "A_to_B", .p = { [0] = 1.0, [3] = 300.0 } } }
  };
  const chem_mechanism mechanism = {
    .name = "small DAE",
    .species_count = 3u,
    .reaction_count = 1u,
    .reactant_term_count = 1u,
    .product_term_count = 1u,
    .species = species,
    .reactions = reactions,
    .reactant_terms = reactants,
    .product_terms = products,
    .temperature = 300.0,
    .pressure = 101325.0,
    .relative_tolerance = 1e-4
  };
  const chem_constraint_term constraint_reactants[] = { { 1u, 1.0 } };
  const chem_constraint_term constraint_products[] = { { 2u, 1.0 } };
  const chem_constraint item = equilibrium(
      "B_C_eq", 2u, constraint_reactants, 1u, constraint_products, 1u, 2.0, -2400.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double constraint_parameters[] = { 2.0 }, temperature[] = { 300.0 }, pressure[] = { 101325.0 };
  double expected_keq = 0.0;
  CHECK(chem_constraint_equilibrium_parameter(&item, 270.0, &expected_keq) == 0);
  NEAR(expected_keq, 2.0 * exp((-2400.0 / GAS_CONSTANT) * (1.0 / 270.0 - 1.0 / REFERENCE_TEMPERATURE)), 1e-14);
  unsigned char differential[] = { 1u, 1u, 1u };
  CHECK(chem_constraint_set_species_mask(&set, differential) == 0);
  CHECK(differential[0] == 1u && differential[1] == 1u && differential[2] == 0u);
  for (int method = CHEM_LU_MOZART_IN_PLACE; method <= CHEM_LU_DOOLITTLE; ++method) {
    double state[] = { 1.0, 0.0, 0.0 };
    chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
    solver.h_start_fraction = 0.005;
    solver.h_max_fraction = 0.5;
    solver.lu_method = (chem_lu_method)method;
    CHECK(chem_integrate_group_with_constraints(
        &mechanism, &set, 1u, state, NULL, constraint_parameters, temperature, pressure,
        0.1, &solver, NULL) == 0);
    NEAR(state[0] + state[1], 1.0, 1e-7);
    NEAR(state[2], 2.0 * state[1], 1e-10);
    NEAR(state[0], exp(-0.1), 1e-5);
  }
  return 0;
}

static int dae_equilibrium_reordered(void) {
  const chem_species species[] = {
    { "C", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "B", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "A", 0.0, 0.0, 0.0, 1e-3, 0u },
  };
  const chem_term reactants[] = { { 2u, 1.0 } };
  const chem_term products[] = { { 1u, 1.0 } };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .label = "A_to_B", .p = { [0] = 1.0, [3] = 300.0 } } }
  };
  const chem_mechanism mechanism = {
    .name = "reordered DAE", .species_count = 3u, .reaction_count = 1u,
    .reactant_term_count = 1u, .product_term_count = 1u, .species = species,
    .reactions = reactions, .reactant_terms = reactants, .product_terms = products,
    .temperature = 300.0, .pressure = 101325.0, .relative_tolerance = 1e-4
  };
  const chem_constraint_term cr[] = { { 1u, 1.0 } }, cp[] = { { 0u, 1.0 } };
  const chem_constraint item = equilibrium("B_C_eq", 0u, cr, 1u, cp, 1u, 2.0, -2400.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  const double constraint_parameters[] = { 2.0 }, temperature[] = { 400.0 }, pressure[] = { 101325.0 };
  double state[] = { 0.0, 0.0, 1.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.h_start_fraction = 0.005;
  solver.h_max_fraction = 0.5;
  CHECK(chem_integrate_group_with_constraints(
      &mechanism, &set, 1u, state, NULL, constraint_parameters, temperature, pressure,
      0.1, &solver, NULL) == 0);
  NEAR(state[1] + state[2], 1.0, 0.01);
  NEAR(state[0], 2.0 * state[1], 1e-10);
  return 0;
}

static int dae_equilibrium_coupled(void) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-3, 0u }, { "B", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-3, 0u }, { "D", 0.0, 0.0, 0.0, 1e-3, 0u },
  };
  const chem_term reactants[] = { { 0u, 1.0 } }, products[] = { { 1u, 1.0 } };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .label = "A_to_B", .p = { [0] = 1.0, [3] = 300.0 } } }
  };
  const chem_mechanism mechanism = {
    .name = "coupled DAE", .species_count = 4u, .reaction_count = 1u,
    .reactant_term_count = 1u, .product_term_count = 1u, .species = species,
    .reactions = reactions, .reactant_terms = reactants, .product_terms = products,
    .temperature = 300.0, .pressure = 101325.0, .relative_tolerance = 1e-4
  };
  const chem_constraint_term b[] = { { 1u, 1.0 } }, c[] = { { 2u, 1.0 } }, d[] = { { 3u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "B_C_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 2u,
      .reactants = b, .reactant_count = 1u, .products = c, .product_count = 1u,
      .equilibrium_constant_ref = 3.0, .delta_enthalpy = -2400.0,
      .gas_constant = GAS_CONSTANT, .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 0u },
    { .name = "B_D_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 3u,
      .reactants = b, .reactant_count = 1u, .products = d, .product_count = 1u,
      .equilibrium_constant_ref = 5.0, .delta_enthalpy = -2400.0,
      .gas_constant = GAS_CONSTANT, .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 1u },
  };
  const chem_constraint_set set = constraint_set(4u, items, 2u, 2u);
  const double constraint_parameters[] = { 3.0, 5.0 }, temperature[] = { 300.0 }, pressure[] = { 101325.0 };
  double state[] = { 1.0, 0.0, 0.0, 0.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.h_start_fraction = 0.005;
  solver.h_max_fraction = 0.5;
  CHECK(chem_integrate_group_with_constraints(
      &mechanism, &set, 1u, state, NULL, constraint_parameters, temperature, pressure,
      0.1, &solver, NULL) == 0);
  NEAR(state[0] + state[1], 1.0, 0.01);
  NEAR(state[2], 3.0 * state[1], 1e-10);
  NEAR(state[3], 5.0 * state[1], 1e-10);
  return 0;
}

static int dae_equilibrium_nonunit(void) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-3, 0u }, { "B", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-3, 0u },
  };
  const chem_term reactants[] = { { 2u, 1.0 } }, products[] = { { 0u, 1.0 } };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .label = "C_to_A", .p = { [0] = 0.5, [3] = 300.0 } } }
  };
  const chem_mechanism mechanism = {
    .name = "nonunit DAE", .species_count = 3u, .reaction_count = 1u,
    .reactant_term_count = 1u, .product_term_count = 1u, .species = species,
    .reactions = reactions, .reactant_terms = reactants, .product_terms = products,
    .temperature = 298.0, .pressure = 101325.0, .relative_tolerance = 1e-4
  };
  const chem_constraint_term ar[] = { { 0u, 2.0 } }, bp[] = { { 1u, 1.0 } };
  const chem_constraint item = equilibrium("A2_B_eq", 1u, ar, 1u, bp, 1u, 10.0, -2400.0, 0u);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 1u);
  double equilibrium_constant = 0.0;
  CHECK(chem_constraint_equilibrium_parameter(&item, 298.0, &equilibrium_constant) == 0);
  const double constraint_parameters[] = { equilibrium_constant }, temperature[] = { 298.0 }, pressure[] = { 101325.0 };
  double state[] = { 0.1, equilibrium_constant * 0.01, 1.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.h_start_fraction = 0.005;
  solver.h_max_fraction = 0.5;
  for (size_t step = 0; step < 50u; ++step) {
    CHECK(chem_integrate_group_with_constraints(
        &mechanism, &set, 1u, state, NULL, constraint_parameters, temperature, pressure,
        0.001, &solver, NULL) == 0);
    NEAR(equilibrium_constant * state[0] * state[0], state[1], 1e-5);
  }
  NEAR(state[0] + state[2], 1.1, 0.01);
  return 0;
}

static int run_linear_budget(const chem_rosenbrock_parameters *solver_parameters) {
  const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-12, 0u }, { "B", 0.0, 0.0, 0.0, 1e-12, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-12, 0u },
  };
  const chem_term reactants[] = { { 0u, 1.0 } }, products[] = { { 1u, 1.0 } };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .label = "A_to_B_fast", .p = { [0] = 1e4, [3] = 300.0 } } }
  };
  const chem_mechanism mechanism = {
    .name = "linear budget DAE", .species_count = 3u, .reaction_count = 1u,
    .reactant_term_count = 1u, .product_term_count = 1u, .species = species,
    .reactions = reactions, .reactant_terms = reactants, .product_terms = products,
    .temperature = 298.0, .pressure = 101325.0, .relative_tolerance = 1e-6
  };
  const chem_constraint_term terms[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_constraint item = linear("mass", 2u, terms, 3u, 1e-6);
  const chem_constraint_set set = constraint_set(3u, &item, 1u, 0u);
  const double temperature[] = { 298.0 }, pressure[] = { 101325.0 };
  double state[] = { 0.9e-6, 0.0, 0.1e-6 };
  CHECK(chem_integrate_group_with_constraints(
      &mechanism, &set, 1u, state, NULL, NULL, temperature, pressure,
      30.0, solver_parameters, NULL) == 0);
  NEAR(state[0] + state[1] + state[2], 1e-6, 1e-12);
  CHECK(state[2] >= -1e-18);
  CHECK(state[0] < 1e-12);
  NEAR(state[1], 0.9e-6, 1e-10);
  NEAR(state[2], 0.1e-6, 1e-10);
  return 0;
}

static int dae_linear_budget(void) {
  const chem_rosenbrock_parameters solver = chem_rosenbrock_four_stage_da_parameters();
  return run_linear_budget(&solver);
}

static int dae_equilibrium_conservation(void) {
  const chem_species species[] = {
    { "A_gas", 0.0, 0.0, 0.0, 1e-8, 0u }, { "A_aq", 0.0, 0.0, 0.0, 1e-8, 0u },
    { "P", 0.0, 0.0, 0.0, 1e-12, 0u },
  };
  const chem_term reactants[] = { { 1u, 1.0 } }, products[] = { { 2u, 1.0 } };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .label = "A_aq_to_P", .p = { [0] = 1e3, [3] = 300.0 } } }
  };
  const chem_mechanism mechanism = {
    .name = "equilibrium and conservation DAE", .species_count = 3u, .reaction_count = 1u,
    .reactant_term_count = 1u, .product_term_count = 1u, .species = species,
    .reactions = reactions, .reactant_terms = reactants, .product_terms = products,
    .temperature = 298.0, .pressure = 101325.0, .relative_tolerance = 1e-6
  };
  const chem_constraint_term er[] = { { 0u, 1.0 } }, ep[] = { { 1u, 1.0 } };
  const chem_constraint_term lt[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "gas_aq_eq", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 1u,
      .reactants = er, .reactant_count = 1u, .products = ep, .product_count = 1u,
      .equilibrium_constant_ref = 10.0, .delta_enthalpy = 0.0,
      .gas_constant = GAS_CONSTANT, .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 0u },
    { .name = "mass", .kind = CHEM_CONSTRAINT_LINEAR, .algebraic_species = 0u,
      .terms = lt, .term_count = 3u, .constant = 1e-6 },
  };
  const chem_constraint_set set = constraint_set(3u, items, 2u, 1u);
  const double constraint_parameters[] = { 10.0 }, temperature[] = { 298.0 }, pressure[] = { 101325.0 };
  double state[] = { 1e-6 / 11.0, 10e-6 / 11.0, 0.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_four_stage_da_parameters();
  CHECK(chem_integrate_group_with_constraints(
      &mechanism, &set, 1u, state, NULL, constraint_parameters, temperature, pressure,
      30.0, &solver, NULL) == 0);
  NEAR(state[0] + state[1] + state[2], 1e-6, 1e-12);
  CHECK(state[0] >= -1e-18 && state[1] >= -1e-18);
  NEAR(state[2], 1e-6, 1e-8);
  return 0;
}

static int dae_all_rosenbrock_orders(void) {
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.h_start_fraction = 1e-6;
  CHECK(run_linear_budget(&solver) == 0);
  solver = chem_rosenbrock_four_stage_parameters();
  CHECK(run_linear_budget(&solver) == 0);
  solver = chem_rosenbrock_four_stage_da_parameters();
  CHECK(run_linear_budget(&solver) == 0);
  solver = chem_rosenbrock_six_stage_da_parameters();
  CHECK(run_linear_budget(&solver) == 0);
  return 0;
}

static int run_cascade_error_case(
    double balance_tolerance, double rate, double k1, double k2,
    double *balance_final, double *product_final, unsigned long *accepted) {
  chem_species species[] = {
    { "A_gas", 0.0, 0.0, 0.0, balance_tolerance, 0u },
    { "A_aq", 0.0, 0.0, 0.0, 1e-8, 0u },
    { "B_aq", 0.0, 0.0, 0.0, 1e-8, 0u },
    { "P", 0.0, 0.0, 0.0, 1e-12, 0u },
  };
  const chem_term reactants[] = { { 2u, 1.0 } }, products[] = { { 3u, 1.0 } };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .label = "B_to_P", .p = { [0] = rate, [3] = 300.0 } } }
  };
  const chem_mechanism mechanism = {
    .name = "cascade DAE", .species_count = 4u, .reaction_count = 1u,
    .reactant_term_count = 1u, .product_term_count = 1u, .species = species,
    .reactions = reactions, .reactant_terms = reactants, .product_terms = products,
    .temperature = 298.0, .pressure = 101325.0, .relative_tolerance = 1e-6
  };
  const chem_constraint_term r1[] = { { 0u, 1.0 } }, p1[] = { { 1u, 1.0 } };
  const chem_constraint_term r2[] = { { 1u, 1.0 } }, p2[] = { { 2u, 1.0 } };
  const chem_constraint_term mass[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 3u, 1.0 } };
  const chem_constraint items[] = {
    { .name = "eq1", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 1u,
      .reactants = r1, .reactant_count = 1u, .products = p1, .product_count = 1u,
      .equilibrium_constant_ref = k1, .delta_enthalpy = 0.0,
      .gas_constant = GAS_CONSTANT, .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 0u },
    { .name = "eq2", .kind = CHEM_CONSTRAINT_EQUILIBRIUM, .algebraic_species = 2u,
      .reactants = r2, .reactant_count = 1u, .products = p2, .product_count = 1u,
      .equilibrium_constant_ref = k2, .delta_enthalpy = 0.0,
      .gas_constant = GAS_CONSTANT, .reference_temperature = REFERENCE_TEMPERATURE, .parameter_index = 1u },
    { .name = "mass", .kind = CHEM_CONSTRAINT_LINEAR, .algebraic_species = 0u,
      .terms = mass, .term_count = 4u, .constant = 1e-6 },
  };
  const chem_constraint_set set = constraint_set(4u, items, 3u, 2u);
  const double constraint_parameters[] = { k1, k2 }, temperature[] = { 298.0 }, pressure[] = { 101325.0 };
  const double denominator = 1.0 + k1 + k1 * k2;
  double state[] = { 1e-6 / denominator, k1 * 1e-6 / denominator,
      k1 * k2 * 1e-6 / denominator, 0.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_four_stage_da_parameters();
  solver.h_start_fraction = 1e-6;
  chem_statistics stats = { 0 };
  if (chem_integrate_group_with_constraints(
          &mechanism, &set, 1u, state, NULL, constraint_parameters, temperature, pressure,
          30.0, &solver, &stats) != 0) return -1;
  *balance_final = state[0];
  *product_final = state[3];
  *accepted = stats.accepted;
  return 0;
}

static int dae_error_sensitive_to_balance_tolerance(void) {
  double loose_balance = 0.0, loose_product = 0.0, tight_balance = 0.0, tight_product = 0.0;
  unsigned long loose_steps = 0u, tight_steps = 0u;
  CHECK(run_cascade_error_case(1e-3, 100.0, 2.0, 2.0, &loose_balance, &loose_product, &loose_steps) == 0);
  CHECK(run_cascade_error_case(1e-8, 100.0, 2.0, 2.0, &tight_balance, &tight_product, &tight_steps) == 0);
  (void)loose_balance; (void)loose_product; (void)tight_balance; (void)tight_product;
  CHECK(tight_steps > loose_steps);
  return 0;
}

static int dae_error_prevents_overshoot(void) {
  double balance = 0.0, product = 0.0;
  unsigned long accepted = 0u;
  CHECK(run_cascade_error_case(1e-8, 1e4, 100.0, 50.0, &balance, &product, &accepted) == 0);
  (void)accepted;
  CHECK(balance >= -1e-8);
  NEAR(product, 1e-6, 1e-8);
  NEAR(balance, 0.0, 1e-10);
  return 0;
}

static int dae_terminator_and_robertson(void) {
  const chem_species species[] = {
    { "Cl2", 0.0, 0.0, 0.0, 1e-20, 0u },
    { "Cl", 0.0, 0.0, 0.0, 1e-20, 0u },
    { "A", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "B", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-3, 0u },
  };
  const chem_term reactants[] = {
    { 0u, 1.0 }, { 1u, 1.0 }, { 1u, 1.0 },
    { 2u, 1.0 }, { 3u, 1.0 }, { 3u, 1.0 }, { 3u, 1.0 }, { 4u, 1.0 },
  };
  const chem_term products[] = {
    { 1u, 2.0 }, { 0u, 1.0 },
    { 3u, 1.0 }, { 3u, 1.0 }, { 4u, 1.0 }, { 2u, 1.0 }, { 4u, 1.0 },
  };
  const chem_parameter parameters[] = {
    { "terminator_k1", 0.0 }, { "robertson_r1", 0.04 },
    { "robertson_r2", 3.0e7 }, { "robertson_r3", 1.0e4 },
  };
  const chem_reaction reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
    { 1u, 2u, 1u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 1.0, [3] = 298.0 } } },
    { 3u, 1u, 2u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 1u } },
    { 4u, 2u, 3u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 2u } },
    { 6u, 2u, 5u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 3u } },
  };
  const chem_mechanism mechanism = {
    .name = "Terminator and Robertson", .species_count = 5u, .reaction_count = 5u,
    .reactant_term_count = 8u, .product_term_count = 7u, .parameter_count = 4u,
    .species = species, .reactions = reactions, .reactant_terms = reactants,
    .product_terms = products, .parameters = parameters, .temperature = 298.0,
    .pressure = 101300.0, .relative_tolerance = 1e-8
  };
  const chem_constraint_term conservation[] = { { 2u, 1.0 }, { 3u, 1.0 }, { 4u, 1.0 } };
  const chem_constraint item = linear("mass_conservation", 4u, conservation, 3u, 1.0);
  const chem_constraint_set set = constraint_set(5u, &item, 1u, 0u);
  const double temperatures[] = { 298.0 }, pressures[] = { 101300.0 };
  double state[] = { 1.8e-10, 1.2e-6, 1.0, 0.0, 0.0 };
  double rate_parameters[] = { 0.0, 0.04, 3.0e7, 1.0e4 };
  const double total_chlorine = 2.0 * state[0] + state[1];
  chem_rosenbrock_parameters solver = chem_rosenbrock_four_stage_da_parameters();
  solver.max_internal_steps = 100000u;

  double interval = 1.0;
  for (size_t step = 0; step < 12u; ++step) {
    chem_statistics stats = { 0 };
    if (chem_integrate_group_with_constraints(
            &mechanism, &set, 1u, state, rate_parameters, NULL, temperatures, pressures,
            interval, &solver, &stats) != 0) {
      fprintf(stderr, "Terminator/Robertson integration failed at interval %.17g s after %lu accepted steps\n",
          interval, stats.accepted);
      return 1;
    }

    NEAR(state[2] + state[3] + state[4], 1.0, 1e-10);
    if (fabs(2.0 * state[0] + state[1] - total_chlorine) > 1e-6 * total_chlorine) {
      fprintf(stderr, "chlorine drift at interval %.17g: Cl2=%.17g Cl=%.17g total=%.17g expected=%.17g\n",
          interval, state[0], state[1], 2.0 * state[0] + state[1], total_chlorine);
      return 1;
    }
    for (size_t species_index = 0; species_index < 5u; ++species_index) {
      CHECK(isfinite(state[species_index]));
      CHECK(state[species_index] >= -1e-18);
    }
    if (step == 0u) {
      CHECK(state[2] < 0.99);
      CHECK(state[3] > 0.0);
      CHECK(state[4] > 1e-3);
    }
    interval *= 10.0;
  }

  CHECK(state[4] > 0.95);
  CHECK(state[2] < 0.05);
  CHECK(state[0] > 0.4 * total_chlorine);
  CHECK(state[1] < 0.1 * total_chlorine);
  return 0;
}

typedef struct { const char *name; int (*run)(void); } test_case;

int main(int argc, char **argv) {
  static const test_case tests[] = {
    { "EquilibriumConstraint.Construction", equilibrium_construction },
    { "EquilibriumConstraint.AlgebraicSpecies", equilibrium_algebraic_species },
    { "EquilibriumConstraint.SingleReactantSingleProduct", equilibrium_single_reactant },
    { "EquilibriumConstraint.MultipleReactantsAndProducts", equilibrium_multiple_terms },
    { "EquilibriumConstraint.InvalidEquilibriumConstant", equilibrium_invalid_constant },
    { "EquilibriumConstraint.EmptyReactantsThrows", equilibrium_empty_reactants },
    { "EquilibriumConstraint.EmptyProductsThrows", equilibrium_empty_products },
    { "EquilibriumConstraint.InvalidStoichiometryThrows", equilibrium_invalid_stoichiometry },
    { "EquilibriumConstraint.ResidualComputationThroughConstraintSet", equilibrium_residual_computation },
    { "EquilibriumConstraint.JacobianComputationThroughConstraintSet", equilibrium_jacobian_computation },
    { "EquilibriumConstraint.ComplexStoichiometryResidual", equilibrium_complex_residual },
    { "EquilibriumConstraint.FiniteDifferenceJacobianSimple", equilibrium_fd_simple },
    { "EquilibriumConstraint.FiniteDifferenceJacobianComplexStoichiometry", equilibrium_fd_complex },
    { "LinearConstraint.Construction", linear_construction },
    { "LinearConstraint.AlgebraicSpecies", linear_algebraic_species },
    { "LinearConstraint.WeightedTerms", linear_weighted_terms },
    { "LinearConstraint.ZeroConstant", linear_zero_constant },
    { "LinearConstraint.FractionalCoefficients", linear_fractional_coefficients },
    { "LinearConstraint.ResidualComputationThroughConstraintSet", linear_residual },
    { "LinearConstraint.JacobianComputationThroughConstraintSet", linear_jacobian },
    { "LinearConstraint.WeightedSumResidualAndJacobian", linear_weighted_residual_jacobian },
    { "LinearConstraint.ThreeSpeciesConservationResidual", linear_conservation_residual },
    { "LinearConstraint.ZeroConstantResidual", linear_zero_constant_residual },
    { "LinearConstraint.FractionalCoefficientsResidualAndJacobian", linear_fractional_residual_jacobian },
    { "LinearConstraint.JacobianIndependentOfConcentrations", linear_jacobian_independent },
    { "LinearConstraint.FiniteDifferenceJacobianSimpleConservation", linear_fd_simple },
    { "LinearConstraint.FiniteDifferenceJacobianWeightedSum", linear_fd_weighted },
    { "ConstraintSet.Construction", set_construction },
    { "ConstraintSet.ReplaceStateRowsMapsToAlgebraicSpecies", set_algebraic_rows },
    { "ConstraintSet.NonZeroJacobianElements", set_nonzero_pattern },
    { "ConstraintSet.MultipleConstraints", set_multiple_constraints },
    { "ConstraintSet.AddForcingTerms", set_add_forcing },
    { "ConstraintSet.SubtractJacobianTerms", set_subtract_jacobian },
    { "ConstraintSet.EmptyConstraintSet", set_empty },
    { "ConstraintSet.UnknownSpeciesThrows", set_unknown_species },
    { "ConstraintSet.ThreeDStateOneConstraint", set_three_species_one_constraint },
    { "ConstraintSet.FourDStateTwoConstraints", set_four_species_two_constraints },
    { "ConstraintSet.CoupledConstraintsSharedSpecies", set_coupled_constraints },
    { "ConstraintSet.VectorizedMatricesRespectGridCellIndexing", set_vectorized_indexing },
    { "ConstraintSet.VectorMatrix1", constraint_vector_matrix_1 },
    { "ConstraintSet.VectorMatrix2", constraint_vector_matrix_2 },
    { "ConstraintSet.VectorMatrix3", constraint_vector_matrix_3 },
    { "ConstraintInitialization.ConsistentICsUnchanged", initialization_consistent },
    { "ConstraintInitialization.MildlyInconsistentICsCorrected", initialization_mild },
    { "ConstraintInitialization.SeverelyInconsistentICsConverge", initialization_severe },
    { "ConstraintInitialization.ScaledConstraintRowsSolveIdenticallyThroughSolve", initialization_scaled_solve },
    { "ConstraintInitialization.PureODESystemUnaffected", initialization_empty_set },
    { "ConstraintInitialization.MultiCellSystems", initialization_multicell },
    { "ConstraintInitialization.SubsequentSolveCallsReinitialize", initialization_repeated },
    { "ConstraintInitialization.SolverStateToStringNewValue", initialization_status_name },
    { "ConstraintInitialization.FailedInitializationRestoresCallerState", initialization_failure_restores_state },
    { "ConstraintInitialization.SuccessfulInitializationKeepsTheCorrection", initialization_damped_success_and_limits },
    { "ConstraintInitialization.WeightedCorrectionIsInvariantToConstraintRowScaling", initialization_weighted_row_scale },
    { "ConstraintInitialization.WeightedCorrectionUsesStateTolerances", initialization_weighted_tolerances },
    { "DAEIntegration.SimpleEquilibrium", dae_equilibrium_simple },
    { "EquilibriumIntegration.DAESolveWithConstraint", dae_equilibrium_simple },
    { "EquilibriumIntegration.SetConstraintsAPIWorks", dae_equilibrium_simple },
    { "EquilibriumIntegration.DAESolveWithConstraintAndReorderState", dae_equilibrium_reordered },
    { "EquilibriumIntegration.SetConstraintsAPIMultipleConstraints", dae_equilibrium_coupled },
    { "EquilibriumIntegration.DAESolveWithTwoCoupledConstraints", dae_equilibrium_coupled },
    { "EquilibriumIntegration.DAESolveWithNonUnitStoichiometry", dae_equilibrium_nonunit },
    { "DAEConstraintOvershoot.AlgebraicVariableStaysNonNegative", dae_linear_budget },
    { "DAEConstraintOvershoot.EquilibriumPlusConservation", dae_equilibrium_conservation },
    { "DAEConstraintOvershoot.AllRosenbrockOrdersConstrained", dae_all_rosenbrock_orders },
    { "DAEAlgebraicError.ErrorSensitiveToBalanceAtol", dae_error_sensitive_to_balance_tolerance },
    { "DAEAlgebraicError.AlgebraicVariableDoesNotOvershootDeeply", dae_error_prevents_overshoot },
    { "DAESolveWithConstraint.TerminatorAndRobertson", dae_terminator_and_robertson },
  };
  if (argc != 2) return 2;
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i) {
    if (strcmp(argv[1], tests[i].name) == 0) {
      return tests[i].run();
    }
  }
  fprintf(stderr, "unknown constraint test: %s\n", argv[1]);
  return 2;
}
