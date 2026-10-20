#include "chem_constraints.h"
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
