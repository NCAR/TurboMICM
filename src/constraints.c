#include "chem_constraints.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static int term_list_valid(const chem_constraint_term *terms, size_t count, size_t species_count, int positive) {
  if (count == 0u || terms == NULL) return 0;
  for (size_t i = 0; i < count; ++i) {
    if (terms[i].species >= species_count || !isfinite(terms[i].coefficient) ||
        (positive && terms[i].coefficient <= 0.0)) return 0;
  }
  return 1;
}

int chem_constraint_find_species(
    const char *name, const char *const *species_names, size_t species_count, unsigned *species_index) {
  if (name == NULL || species_names == NULL || species_index == NULL) return -1;
  for (size_t i = 0; i < species_count; ++i) {
    if (species_names[i] != NULL && strcmp(name, species_names[i]) == 0) {
      if (i > (size_t)UINT_MAX) return -1;
      *species_index = (unsigned)i;
      return 0;
    }
  }
  return -1;
}

int chem_constraint_validate(const chem_constraint *constraint, size_t species_count) {
  if (constraint == NULL || constraint->name == NULL ||
      species_count == 0u || constraint->algebraic_species >= species_count) return -1;
  if (constraint->kind == CHEM_CONSTRAINT_LINEAR) {
    if (!term_list_valid(constraint->terms, constraint->term_count, species_count, 0) ||
        !isfinite(constraint->constant)) return -1;
    return 0;
  }
  if (constraint->kind != CHEM_CONSTRAINT_EQUILIBRIUM ||
      !term_list_valid(constraint->reactants, constraint->reactant_count, species_count, 1) ||
      !term_list_valid(constraint->products, constraint->product_count, species_count, 1) ||
      !isfinite(constraint->equilibrium_constant_ref) || constraint->equilibrium_constant_ref <= 0.0 ||
      !isfinite(constraint->delta_enthalpy) || !isfinite(constraint->gas_constant) ||
      constraint->gas_constant <= 0.0 || !isfinite(constraint->reference_temperature) ||
      constraint->reference_temperature <= 0.0) return -1;
  return 0;
}

int chem_constraint_set_validate(const chem_constraint_set *set) {
  if (set == NULL || (set->constraint_count != 0u && set->constraints == NULL)) return -1;
  if (set->species_count != 0u && set->species_count > SIZE_MAX / set->species_count) return -1;
  for (size_t i = 0; i < set->constraint_count; ++i) {
    const chem_constraint *constraint = &set->constraints[i];
    if (chem_constraint_validate(constraint, set->species_count) != 0) return -1;
    if (constraint->kind == CHEM_CONSTRAINT_EQUILIBRIUM &&
        constraint->parameter_index >= set->state_parameter_count) return -1;
    for (size_t j = 0; j < i; ++j)
      if (set->constraints[j].algebraic_species == constraint->algebraic_species) return -1;
  }
  return 0;
}

size_t chem_constraint_set_size(const chem_constraint_set *set) {
  return set == NULL ? 0u : set->constraint_count;
}

size_t chem_constraint_set_parameter_count(const chem_constraint_set *set) {
  return set == NULL ? 0u : set->state_parameter_count;
}

int chem_constraint_set_species_mask(const chem_constraint_set *set, unsigned char *differential_species) {
  if (chem_constraint_set_validate(set) != 0 || differential_species == NULL || set->species_count == 0u) return -1;
  memset(differential_species, 1, set->species_count * sizeof(*differential_species));
  for (size_t i = 0; i < set->constraint_count; ++i)
    differential_species[set->constraints[i].algebraic_species] = 0u;
  return 0;
}

int chem_constraint_set_nonzero_pattern(const chem_constraint_set *set, unsigned char *pattern) {
  if (chem_constraint_set_validate(set) != 0 || pattern == NULL || set->species_count == 0u) return -1;
  const size_t n = set->species_count;
  memset(pattern, 0, n * n * sizeof(*pattern));
  for (size_t i = 0; i < set->constraint_count; ++i) {
    const chem_constraint *constraint = &set->constraints[i];
    const size_t row = constraint->algebraic_species;
    pattern[row * n + row] = 1u;
    const chem_constraint_term *terms = constraint->kind == CHEM_CONSTRAINT_LINEAR
        ? constraint->terms : constraint->reactants;
    size_t count = constraint->kind == CHEM_CONSTRAINT_LINEAR
        ? constraint->term_count : constraint->reactant_count;
    for (size_t j = 0; j < count; ++j) pattern[row * n + terms[j].species] = 1u;
    if (constraint->kind == CHEM_CONSTRAINT_EQUILIBRIUM)
      for (size_t j = 0; j < constraint->product_count; ++j)
        pattern[row * n + constraint->products[j].species] = 1u;
  }
  return 0;
}

int chem_constraint_equilibrium_parameter(
    const chem_constraint *constraint, double temperature, double *value) {
  if (chem_constraint_validate(constraint, SIZE_MAX) != 0 || constraint->kind != CHEM_CONSTRAINT_EQUILIBRIUM ||
      !isfinite(temperature) || temperature <= 0.0 || value == NULL) return -1;
  *value = constraint->equilibrium_constant_ref * exp((constraint->delta_enthalpy / constraint->gas_constant) *
      (1.0 / temperature - 1.0 / constraint->reference_temperature));
  return isfinite(*value) && *value > 0.0 ? 0 : -1;
}

static double concentration_power(double value, double power) {
  return pow(value > 0.0 ? value : 0.0, power);
}

static double equilibrium_residual(
    const chem_constraint *constraint, const double *state, double equilibrium_constant) {
  double reactant_product = equilibrium_constant;
  double product_product = 1.0;
  for (size_t i = 0; i < constraint->reactant_count; ++i) {
    const chem_constraint_term *term = &constraint->reactants[i];
    reactant_product *= concentration_power(state[term->species], term->coefficient);
  }
  for (size_t i = 0; i < constraint->product_count; ++i) {
    const chem_constraint_term *term = &constraint->products[i];
    product_product *= concentration_power(state[term->species], term->coefficient);
  }
  return reactant_product - product_product;
}

static double equilibrium_derivative(
    const chem_constraint_term *terms,
    size_t term_count,
    size_t term_index,
    const double *state,
    double scale) {
  const chem_constraint_term *selected = &terms[term_index];
  const double selected_value = state[selected->species];
  if (selected_value <= 0.0) return 0.0;
  double derivative = scale * selected->coefficient * pow(selected_value, selected->coefficient - 1.0);
  for (size_t i = 0; i < term_count; ++i) {
    if (i == term_index) continue;
    const chem_constraint_term *other = &terms[i];
    derivative *= concentration_power(state[other->species], other->coefficient);
  }
  return derivative;
}

int chem_constraint_set_apply_forcing(
    const chem_constraint_set *set,
    size_t cell_count,
    const double *state,
    const double *state_parameters,
    double *forcing) {
  if (chem_constraint_set_validate(set) != 0 || set->species_count == 0u || cell_count == 0u ||
      state == NULL || forcing == NULL ||
      set->species_count > SIZE_MAX / cell_count / sizeof(double) ||
      set->state_parameter_count > SIZE_MAX / cell_count / sizeof(double) ||
      (set->state_parameter_count != 0u && state_parameters == NULL)) return -1;
  const size_t n = set->species_count;
  for (size_t cell = 0; cell < cell_count; ++cell) {
    const double *cell_state = state + cell * n;
    double *cell_forcing = forcing + cell * n;
    const double *cell_parameters = set->state_parameter_count == 0u ? NULL :
        state_parameters + cell * set->state_parameter_count;
    for (size_t i = 0; i < set->constraint_count; ++i) {
      const chem_constraint *constraint = &set->constraints[i];
      double residual = 0.0;
      if (constraint->kind == CHEM_CONSTRAINT_LINEAR) {
        residual = -constraint->constant;
        for (size_t j = 0; j < constraint->term_count; ++j)
          residual += constraint->terms[j].coefficient * cell_state[constraint->terms[j].species];
      } else {
        const double keq = cell_parameters[constraint->parameter_index];
        residual = equilibrium_residual(constraint, cell_state, keq);
      }
      cell_forcing[constraint->algebraic_species] = residual;
    }
  }
  return 0;
}

int chem_constraint_set_accumulate_negative_jacobian(
    const chem_constraint_set *set,
    size_t cell_count,
    const double *state,
    const double *state_parameters,
    double *negative_jacobian) {
  if (chem_constraint_set_validate(set) != 0 || set->species_count == 0u || cell_count == 0u ||
      state == NULL || negative_jacobian == NULL ||
      set->species_count > SIZE_MAX / set->species_count ||
      cell_count > SIZE_MAX / (set->species_count * set->species_count) / sizeof(double) ||
      set->state_parameter_count > SIZE_MAX / cell_count / sizeof(double) ||
      (set->state_parameter_count != 0u && state_parameters == NULL)) return -1;
  const size_t n = set->species_count;
  for (size_t cell = 0; cell < cell_count; ++cell) {
    const double *cell_state = state + cell * n;
    const double *cell_parameters = set->state_parameter_count == 0u ? NULL :
        state_parameters + cell * set->state_parameter_count;
    for (size_t i = 0; i < set->constraint_count; ++i) {
      const chem_constraint *constraint = &set->constraints[i];
      const size_t row_offset = (cell * n + constraint->algebraic_species) * n;
      if (constraint->kind == CHEM_CONSTRAINT_LINEAR) {
        for (size_t j = 0; j < constraint->term_count; ++j)
          negative_jacobian[row_offset + constraint->terms[j].species] -= constraint->terms[j].coefficient;
        continue;
      }
      const double keq = cell_parameters[constraint->parameter_index];
      for (size_t j = 0; j < constraint->reactant_count; ++j) {
        const double derivative = equilibrium_derivative(
            constraint->reactants, constraint->reactant_count, j, cell_state, keq);
        negative_jacobian[row_offset + constraint->reactants[j].species] -= derivative;
      }
      for (size_t j = 0; j < constraint->product_count; ++j) {
        const double derivative = equilibrium_derivative(
            constraint->products, constraint->product_count, j, cell_state, 1.0);
        negative_jacobian[row_offset + constraint->products[j].species] += derivative;
      }
    }
  }
  return 0;
}
