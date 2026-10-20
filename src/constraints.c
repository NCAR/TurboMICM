#include "chemistry.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static int valid_terms(const chem_constraint_term *terms, size_t count,
                       size_t species_count, int require_positive) {
  if (count && terms == NULL) return 0;
  for (size_t i = 0; i < count; ++i)
    if (terms[i].species >= species_count || !isfinite(terms[i].coefficient) ||
        (require_positive ? terms[i].coefficient <= 0.0 : terms[i].coefficient == 0.0))
      return 0;
  return 1;
}

int chem_constraint_set_validate(const chem_constraint_set *set) {
  if (set == NULL || set->species_count == 0 ||
      (set->constraint_count && set->constraints == NULL)) return 1;
  for (size_t i = 0; i < set->constraint_count; ++i) {
    const chem_constraint *c = &set->constraints[i];
    if (c->algebraic_species >= set->species_count) return 1;
    if (c->kind == CHEM_CONSTRAINT_LINEAR) {
      if (c->term_count == 0 || !isfinite(c->constant) ||
          !valid_terms(c->terms, c->term_count, set->species_count, 0)) return 1;
    } else if (c->kind == CHEM_CONSTRAINT_EQUILIBRIUM) {
      if (c->reactant_count == 0 || c->product_count == 0 ||
          !isfinite(c->equilibrium_constant_ref) || c->equilibrium_constant_ref <= 0.0 ||
          !isfinite(c->enthalpy_change) ||
          !isfinite(c->reference_temperature) || c->reference_temperature <= 0.0 ||
          !valid_terms(c->reactants, c->reactant_count, set->species_count, 1) ||
          !valid_terms(c->products, c->product_count, set->species_count, 1)) return 1;
    } else return 1;
  }
  return 0;
}

int chem_constraint_set_algebraic_mask(const chem_constraint_set *set,
                                       unsigned char *mask) {
  if (chem_constraint_set_validate(set) || (set->species_count && mask == NULL)) return 1;
  memset(mask, 0, set->species_count * sizeof(*mask));
  for (size_t i = 0; i < set->constraint_count; ++i)
    mask[set->constraints[i].algebraic_species] = 1;
  return 0;
}

int chem_constraint_set_jacobian_pattern(const chem_constraint_set *set,
                                         unsigned char *pattern) {
  if (chem_constraint_set_validate(set) ||
      (set->species_count && set->species_count > (size_t)-1/set->species_count) ||
      (set->species_count && pattern == NULL)) return 1;
  const size_t n = set->species_count;
  memset(pattern, 0, n*n*sizeof(*pattern));
  for (size_t i = 0; i < set->constraint_count; ++i) {
    const chem_constraint *c = &set->constraints[i];
    const size_t row = c->algebraic_species;
    pattern[row*n + row] = 1;
    const chem_constraint_term *terms = c->kind == CHEM_CONSTRAINT_LINEAR
        ? c->terms : c->reactants;
    const size_t count = c->kind == CHEM_CONSTRAINT_LINEAR
        ? c->term_count : c->reactant_count;
    for (size_t j = 0; j < count; ++j) pattern[row*n + terms[j].species] = 1;
    if (c->kind == CHEM_CONSTRAINT_EQUILIBRIUM)
      for (size_t j = 0; j < c->product_count; ++j)
        pattern[row*n + c->products[j].species] = 1;
  }
  return 0;
}

static double equilibrium_constant(const chem_constraint *c, double temperature) {
  const double gas_constant = 8.31446261815324;
  return c->equilibrium_constant_ref *
      exp((c->enthalpy_change/gas_constant) *
          (1.0/temperature - 1.0/c->reference_temperature));
}

static double concentration_product(const chem_constraint_term *terms,
                                    size_t count, const double *state) {
  double value = 1.0;
  for (size_t i = 0; i < count; ++i)
    value *= pow(state[terms[i].species], terms[i].coefficient);
  return value;
}

static double constraint_residual(const chem_constraint *c, const double *state,
                                  double temperature) {
  if (c->kind == CHEM_CONSTRAINT_LINEAR) {
    double residual = -c->constant;
    for (size_t i = 0; i < c->term_count; ++i)
      residual += c->terms[i].coefficient * state[c->terms[i].species];
    return residual;
  }
  return equilibrium_constant(c, temperature) *
             concentration_product(c->reactants, c->reactant_count, state) -
         concentration_product(c->products, c->product_count, state);
}

static double product_derivative(const chem_constraint_term *terms, size_t count,
                                 size_t q, const double *state) {
  double value = terms[q].coefficient;
  for (size_t i = 0; i < count; ++i) {
    const double exponent = terms[i].coefficient - (i == q ? 1.0 : 0.0);
    value *= pow(state[terms[i].species], exponent);
  }
  return value;
}

static double residual_derivative(const chem_constraint *c, const double *state,
                                  double temperature, size_t column) {
  double derivative = 0.0;
  if (c->kind == CHEM_CONSTRAINT_LINEAR) {
    for (size_t i = 0; i < c->term_count; ++i)
      if (c->terms[i].species == column) derivative += c->terms[i].coefficient;
    return derivative;
  }
  const double k = equilibrium_constant(c, temperature);
  for (size_t i = 0; i < c->reactant_count; ++i)
    if (c->reactants[i].species == column)
      derivative += k * product_derivative(c->reactants, c->reactant_count, i, state);
  for (size_t i = 0; i < c->product_count; ++i)
    if (c->products[i].species == column)
      derivative -= product_derivative(c->products, c->product_count, i, state);
  return derivative;
}

int chem_constraint_set_add_forcing(const chem_constraint_set *set,
                                    const double *states,
                                    const double *temperatures,
                                    size_t cell_count, double *forcing) {
  if (chem_constraint_set_validate(set) || (cell_count && states == NULL) ||
      (cell_count && forcing == NULL)) return 1;
  for (size_t cell = 0; cell < cell_count; ++cell) {
    const double *state = states + cell*set->species_count;
    for (size_t i = 0; i < set->constraint_count; ++i) {
      const chem_constraint *c = &set->constraints[i];
      const double temperature = temperatures ? temperatures[cell] : c->reference_temperature;
      if (c->kind == CHEM_CONSTRAINT_EQUILIBRIUM &&
          (!isfinite(temperature) || temperature <= 0.0)) return 1;
      forcing[cell*set->species_count + c->algebraic_species] +=
          constraint_residual(c, state, temperature);
    }
  }
  return 0;
}

int chem_constraint_set_add_negative_jacobian(
    const chem_constraint_set *set, const double *states,
    const double *temperatures, size_t cell_count,
    double *negative_jacobians) {
  if (chem_constraint_set_validate(set) || (cell_count && states == NULL) ||
      (cell_count && negative_jacobians == NULL) ||
      (set->species_count && set->species_count > (size_t)-1/set->species_count)) return 1;
  const size_t n = set->species_count;
  for (size_t cell = 0; cell < cell_count; ++cell) {
    const double *state = states + cell*n;
    double *jacobian = negative_jacobians + cell*n*n;
    for (size_t i = 0; i < set->constraint_count; ++i) {
      const chem_constraint *c = &set->constraints[i];
      const double temperature = temperatures ? temperatures[cell] : c->reference_temperature;
      if (c->kind == CHEM_CONSTRAINT_EQUILIBRIUM &&
          (!isfinite(temperature) || temperature <= 0.0)) return 1;
      const size_t row = c->algebraic_species;
      for (size_t column = 0; column < n; ++column)
        jacobian[row*n + column] -= residual_derivative(c, state, temperature, column);
    }
  }
  return 0;
}

int chem_forcing_add_with_algebraic_mask(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, double air_density,
    const unsigned char *algebraic_mask, double *forcing) {
  if (!mechanism || !algebraic_mask || !forcing || !mechanism->species_count)
    return 1;
  double *addition = calloc(mechanism->species_count, sizeof(*addition));
  if (!addition) return 1;
  const int status = chem_forcing_with_constants(mechanism, concentration,
                                                  rate_constants, air_density,
                                                  addition);
  if (!status)
    for (size_t i = 0; i < mechanism->species_count; ++i)
      if (!algebraic_mask[i]) forcing[i] += addition[i];
  free(addition);
  return status;
}

int chem_negative_jacobian_add_with_algebraic_mask(
    const chem_mechanism *mechanism, const double *concentration,
    const double *rate_constants, const unsigned char *algebraic_mask,
    double *negative_jacobian) {
  if (!mechanism || !algebraic_mask || !negative_jacobian ||
      !mechanism->species_count ||
      mechanism->species_count > (size_t)-1/mechanism->species_count) return 1;
  const size_t n = mechanism->species_count;
  double *addition = calloc(n*n, sizeof(*addition));
  if (!addition) return 1;
  const int status = chem_negative_jacobian_with_constants(
      mechanism, concentration, rate_constants, addition);
  if (!status)
    for (size_t row = 0; row < n; ++row)
      if (!algebraic_mask[row])
        for (size_t column = 0; column < n; ++column)
          negative_jacobian[row*n+column] += addition[row*n+column];
  free(addition);
  return status;
}
