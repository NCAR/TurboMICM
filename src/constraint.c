#include "constraint.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define CHEM_DEFAULT_GAS_CONSTANT 8.31446261815324
#define CHEM_DEFAULT_REFERENCE_TEMPERATURE 298.15

static int valid_terms(const chem_constraint_term *terms, size_t count,
                       size_t species_count, int positive_coefficients)
{
  size_t i;
  if (!terms || !count)
    return 0;
  for (i = 0; i < count; ++i)
    if (terms[i].species_index >= species_count ||
        !isfinite(terms[i].coefficient) ||
        (positive_coefficients ? !(terms[i].coefficient > 0.0)
                               : terms[i].coefficient == 0.0))
      return 0;
  return 1;
}

int chem_constraint_linear_init(chem_constraint *constraint, const char *name,
                                size_t algebraic_species,
                                const chem_constraint_term *terms,
                                size_t term_count, double constant)
{
  if (!constraint || !name || !valid_terms(terms, term_count, SIZE_MAX, 0) ||
      !isfinite(constant))
    return -1;
  *constraint = (chem_constraint){
    .name=name, .kind=CHEM_CONSTRAINT_LINEAR,
    .algebraic_species=algebraic_species, .terms=terms,
    .term_count=term_count, .constant=constant
  };
  return 0;
}

int chem_constraint_equilibrium_init(
    chem_constraint *constraint, const char *name, size_t algebraic_species,
    const chem_constraint_term *reactants, size_t reactant_count,
    const chem_constraint_term *products, size_t product_count,
    double equilibrium_constant_reference, double enthalpy_change,
    double gas_constant, double reference_temperature)
{
  if (!constraint || !name || !valid_terms(reactants, reactant_count, SIZE_MAX, 1) ||
      !valid_terms(products, product_count, SIZE_MAX, 1) ||
      !(equilibrium_constant_reference > 0.0) ||
      !isfinite(equilibrium_constant_reference) || !isfinite(enthalpy_change) ||
      !(gas_constant > 0.0) || !isfinite(gas_constant) ||
      !(reference_temperature > 0.0) || !isfinite(reference_temperature))
    return -1;
  *constraint = (chem_constraint){
    .name=name, .kind=CHEM_CONSTRAINT_EQUILIBRIUM,
    .algebraic_species=algebraic_species, .reactants=reactants,
    .reactant_count=reactant_count, .products=products,
    .product_count=product_count,
    .equilibrium_constant_reference=equilibrium_constant_reference,
    .enthalpy_change=enthalpy_change, .gas_constant=gas_constant,
    .reference_temperature=reference_temperature
  };
  return 0;
}

static int validate_constraint(const chem_constraint *constraint,
                               size_t species_count)
{
  if (!constraint || !constraint->name ||
      constraint->algebraic_species >= species_count)
    return 0;
  if (constraint->kind == CHEM_CONSTRAINT_LINEAR)
    return valid_terms(constraint->terms, constraint->term_count,
                       species_count, 0);
  if (constraint->kind == CHEM_CONSTRAINT_EQUILIBRIUM)
    return constraint->equilibrium_constant_reference > 0.0 &&
           isfinite(constraint->equilibrium_constant_reference) &&
           constraint->gas_constant > 0.0 &&
           isfinite(constraint->gas_constant) &&
           constraint->reference_temperature > 0.0 &&
           isfinite(constraint->reference_temperature) &&
           isfinite(constraint->enthalpy_change) &&
           valid_terms(constraint->reactants, constraint->reactant_count,
                       species_count, 1) &&
           valid_terms(constraint->products, constraint->product_count,
                       species_count, 1);
  return 0;
}

int chem_constraint_set_init(chem_constraint_set *set,
                             const chem_constraint *constraints,
                             size_t constraint_count, size_t species_count)
{
  size_t i;
  if (!set || !species_count || (constraint_count && !constraints))
    return -1;
  for (i = 0; i < constraint_count; ++i)
    if (!validate_constraint(&constraints[i], species_count))
      return -1;
  *set = (chem_constraint_set){constraints, constraint_count, species_count};
  return 0;
}

size_t chem_constraint_set_size(const chem_constraint_set *set)
{
  return set ? set->constraint_count : 0;
}

int chem_constraint_set_algebraic_species(const chem_constraint_set *set,
                                          size_t *indices, size_t capacity)
{
  size_t i;
  if (!set || capacity < set->constraint_count ||
      (set->constraint_count && !indices))
    return -1;
  for (i = 0; i < set->constraint_count; ++i)
    indices[i] = set->constraints[i].algebraic_species;
  return 0;
}

static int add_jacobian_element(chem_jacobian_element *elements,
                                size_t capacity, size_t *count,
                                size_t row, size_t column)
{
  size_t i;
  for (i = 0; i < *count; ++i)
    if (elements[i].row == row && elements[i].column == column)
      return 0;
  if (*count >= capacity)
    return -1;
  elements[*count] = (chem_jacobian_element){row, column};
  ++*count;
  return 0;
}

int chem_constraint_set_nonzero_jacobian_elements(
    const chem_constraint_set *set, chem_jacobian_element *elements,
    size_t capacity, size_t *element_count)
{
  size_t i, j, count = 0;
  if (!set || !element_count ||
      (set->constraint_count && (!elements || !capacity)))
    return -1;
  for (i = 0; i < set->constraint_count; ++i) {
    const chem_constraint *constraint = &set->constraints[i];
    const size_t row = constraint->algebraic_species;
    const chem_constraint_term *terms = constraint->kind == CHEM_CONSTRAINT_LINEAR
                                            ? constraint->terms
                                            : constraint->reactants;
    const size_t term_count = constraint->kind == CHEM_CONSTRAINT_LINEAR
                                  ? constraint->term_count
                                  : constraint->reactant_count;
    for (j = 0; j < term_count; ++j)
      if (add_jacobian_element(elements, capacity, &count, row,
                               terms[j].species_index) != 0)
        return -1;
    if (constraint->kind == CHEM_CONSTRAINT_EQUILIBRIUM)
      for (j = 0; j < constraint->product_count; ++j)
        if (add_jacobian_element(elements, capacity, &count, row,
                                 constraint->products[j].species_index) != 0)
          return -1;
    if (add_jacobian_element(elements, capacity, &count, row, row) != 0)
      return -1;
  }
  *element_count = count;
  return 0;
}

static double effective_equilibrium_constant(const chem_constraint *constraint,
                                             double equilibrium_constant,
                                             double temperature)
{
  if (equilibrium_constant > 0.0)
    return equilibrium_constant;
  if (!(temperature > 0.0))
    return NAN;
  return constraint->equilibrium_constant_reference *
         exp((constraint->enthalpy_change / constraint->gas_constant) *
             (1.0 / temperature - 1.0 / constraint->reference_temperature));
}

static double mass_action_product(const chem_constraint_term *terms,
                                  size_t count, const double *concentrations)
{
  size_t i;
  double product = 1.0;
  for (i = 0; i < count; ++i)
    product *= pow(concentrations[terms[i].species_index], terms[i].coefficient);
  return product;
}

static double mass_action_derivative(const chem_constraint_term *terms,
                                     size_t count, const double *concentrations,
                                     size_t species)
{
  size_t i, j;
  double derivative = 0.0;
  for (i = 0; i < count; ++i) {
    double term = 1.0;
    if (terms[i].species_index != species)
      continue;
    term = terms[i].coefficient *
           pow(concentrations[species], terms[i].coefficient - 1.0);
    for (j = 0; j < count; ++j)
      if (j != i)
        term *= pow(concentrations[terms[j].species_index], terms[j].coefficient);
    derivative += term;
  }
  return derivative;
}

int chem_constraint_residual(const chem_constraint *constraint,
                             const double *concentrations, size_t species_count,
                             double equilibrium_constant, double temperature,
                             double *residual)
{
  size_t i;
  if (!residual || !validate_constraint(constraint, species_count) ||
      !concentrations)
    return -1;
  if (constraint->kind == CHEM_CONSTRAINT_LINEAR) {
    double value = -constraint->constant;
    for (i = 0; i < constraint->term_count; ++i)
      value += constraint->terms[i].coefficient *
               concentrations[constraint->terms[i].species_index];
    *residual = value;
  } else {
    const double k = effective_equilibrium_constant(constraint,
                                                    equilibrium_constant,
                                                    temperature);
    if (!(k > 0.0) || !isfinite(k))
      return -1;
    *residual = k * mass_action_product(constraint->reactants,
                                        constraint->reactant_count,
                                        concentrations) -
                mass_action_product(constraint->products,
                                    constraint->product_count,
                                    concentrations);
  }
  return isfinite(*residual) ? 0 : -1;
}

int chem_constraint_jacobian(const chem_constraint *constraint,
                             const double *concentrations, size_t species_count,
                             double equilibrium_constant, double temperature,
                             double *gradient)
{
  size_t species, i;
  if (!gradient || !validate_constraint(constraint, species_count) ||
      !concentrations)
    return -1;
  for (species = 0; species < species_count; ++species)
    gradient[species] = 0.0;
  if (constraint->kind == CHEM_CONSTRAINT_LINEAR) {
    for (i = 0; i < constraint->term_count; ++i)
      gradient[constraint->terms[i].species_index] +=
          constraint->terms[i].coefficient;
  } else {
    const double k = effective_equilibrium_constant(constraint,
                                                    equilibrium_constant,
                                                    temperature);
    if (!(k > 0.0) || !isfinite(k))
      return -1;
    for (species = 0; species < species_count; ++species)
      gradient[species] =
          k * mass_action_derivative(constraint->reactants,
                                    constraint->reactant_count,
                                    concentrations, species) -
          mass_action_derivative(constraint->products,
                                constraint->product_count,
                                concentrations, species);
  }
  for (species = 0; species < species_count; ++species)
    if (!isfinite(gradient[species]))
      return -1;
  return 0;
}

int chem_constraint_set_add_forcing(const chem_constraint_set *set,
                                    size_t cell_count,
                                    const double *concentrations,
                                    const double *equilibrium_constants,
                                    const double *temperatures,
                                    double *forcing)
{
  size_t cell, constraint_index;
  if (!set || (cell_count && (!concentrations || !forcing)) ||
      set->species_count > SIZE_MAX / (cell_count ? cell_count : 1) ||
      (cell_count && set->constraint_count > SIZE_MAX / cell_count))
    return -1;
  for (cell = 0; cell < cell_count; ++cell)
    for (constraint_index = 0; constraint_index < set->constraint_count;
         ++constraint_index) {
      const chem_constraint *constraint = &set->constraints[constraint_index];
      const double k = equilibrium_constants
                           ? equilibrium_constants[cell * set->constraint_count + constraint_index]
                           : 0.0;
      const double temperature = temperatures ? temperatures[cell]
                                               : CHEM_DEFAULT_REFERENCE_TEMPERATURE;
      if (chem_constraint_residual(constraint,
                                   concentrations + cell * set->species_count,
                                   set->species_count, k, temperature,
                                   &forcing[cell * set->species_count +
                                            constraint->algebraic_species]) != 0)
        return -1;
    }
  return 0;
}

int chem_constraint_set_subtract_jacobian(
    const chem_constraint_set *set, size_t cell_count,
    const double *concentrations, const double *equilibrium_constants,
    const double *temperatures, double *dense_jacobians)
{
  double *gradient;
  size_t cell, constraint_index, species, matrix_size;
  if (!set || (cell_count && (!concentrations || !dense_jacobians)) ||
      !set->species_count ||
      set->species_count > SIZE_MAX / set->species_count ||
      (cell_count && set->constraint_count > SIZE_MAX / cell_count))
    return -1;
  matrix_size = set->species_count * set->species_count;
  if (cell_count && matrix_size > SIZE_MAX / cell_count)
    return -1;
  gradient = (double *)malloc(set->species_count * sizeof(*gradient));
  if (!gradient)
    return -1;
  for (cell = 0; cell < cell_count; ++cell)
    for (constraint_index = 0; constraint_index < set->constraint_count;
         ++constraint_index) {
      const chem_constraint *constraint = &set->constraints[constraint_index];
      const double k = equilibrium_constants
                           ? equilibrium_constants[cell * set->constraint_count + constraint_index]
                           : 0.0;
      const double temperature = temperatures ? temperatures[cell]
                                               : CHEM_DEFAULT_REFERENCE_TEMPERATURE;
      double *matrix = dense_jacobians + cell * matrix_size;
      if (chem_constraint_jacobian(constraint,
                                   concentrations + cell * set->species_count,
                                   set->species_count, k, temperature,
                                   gradient) != 0) {
        free(gradient);
        return -1;
      }
      for (species = 0; species < set->species_count; ++species)
        matrix[constraint->algebraic_species * set->species_count + species] -=
            gradient[species];
    }
  free(gradient);
  return 0;
}
