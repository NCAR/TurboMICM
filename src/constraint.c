#include "constraint.h"
#include "matrix.h"

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

int chem_constraint_callback_init(
    chem_constraint *constraint, const char *name, size_t algebraic_species,
    const size_t *dependencies, size_t dependency_count,
    chem_constraint_residual_callback residual_callback,
    chem_constraint_jacobian_callback jacobian_callback, void *context)
{
  size_t i, j;
  if (!constraint || !name || !dependencies || !dependency_count ||
      !residual_callback || !jacobian_callback)
    return -1;
  for (i = 0; i < dependency_count; ++i)
    for (j = 0; j < i; ++j)
      if (dependencies[i] == dependencies[j])
        return -1;
  *constraint = (chem_constraint){
    .name=name, .kind=CHEM_CONSTRAINT_CALLBACK,
    .algebraic_species=algebraic_species,
    .dependencies=dependencies, .dependency_count=dependency_count,
    .residual_callback=residual_callback,
    .jacobian_callback=jacobian_callback, .callback_context=context
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
  if (constraint->kind == CHEM_CONSTRAINT_CALLBACK) {
    size_t i, j;
    if (!constraint->dependencies || !constraint->dependency_count ||
        !constraint->residual_callback || !constraint->jacobian_callback)
      return 0;
    for (i = 0; i < constraint->dependency_count; ++i) {
      if (constraint->dependencies[i] >= species_count)
        return 0;
      for (j = 0; j < i; ++j)
        if (constraint->dependencies[i] == constraint->dependencies[j])
          return 0;
    }
    return 1;
  }
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
    if (constraint->kind == CHEM_CONSTRAINT_CALLBACK) {
      for (j = 0; j < constraint->dependency_count; ++j)
        if (add_jacobian_element(elements, capacity, &count, row,
                                 constraint->dependencies[j]) != 0)
          return -1;
      if (add_jacobian_element(elements, capacity, &count, row, row) != 0)
        return -1;
      continue;
    }
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
  } else if (constraint->kind == CHEM_CONSTRAINT_EQUILIBRIUM) {
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
  } else if (constraint->residual_callback(
                 concentrations, species_count, equilibrium_constant,
                 temperature, residual, constraint->callback_context) != 0) {
    return -1;
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
  } else if (constraint->kind == CHEM_CONSTRAINT_EQUILIBRIUM) {
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
  } else if (constraint->jacobian_callback(
                 concentrations, species_count, equilibrium_constant,
                 temperature, gradient, constraint->callback_context) != 0) {
    return -1;
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

static int constraint_set_evaluate(const chem_constraint_set *set,
                                   const double *concentrations, size_t cell,
                                   const double *equilibrium_constants,
                                   const double *temperatures,
                                   double *residuals, double *jacobian,
                                   double *gradient)
{
  size_t i, j;
  const size_t m = set->constraint_count;
  const size_t n = set->species_count;
  const double temperature = temperatures ? temperatures[cell]
                                           : CHEM_DEFAULT_REFERENCE_TEMPERATURE;
  for (i = 0; i < m; ++i) {
    const chem_constraint *constraint = &set->constraints[i];
    const double k = equilibrium_constants ? equilibrium_constants[cell * m + i] : 0.0;
    if (chem_constraint_residual(constraint, concentrations, n, k, temperature,
                                 &residuals[i]) != 0 || !isfinite(residuals[i]))
      return -1;
    if (jacobian) {
      if (chem_constraint_jacobian(constraint, concentrations, n, k, temperature,
                                   gradient) != 0)
        return -1;
      for (j = 0; j < m; ++j) {
        const double value = gradient[set->constraints[j].algebraic_species];
        if (!isfinite(value))
          return -1;
        jacobian[i * m + j] = value;
      }
    }
  }
  return 0;
}

static double constraint_correction_norm(const chem_constraint_set *set,
                                         const double *concentrations,
                                         const double *correction,
                                         const double *absolute_tolerances,
                                         double relative_tolerance)
{
  size_t i;
  double maximum = 0.0;
  for (i = 0; i < set->constraint_count; ++i) {
    const size_t algebraic = set->constraints[i].algebraic_species;
    const double old_value = concentrations[algebraic];
    const double new_value = old_value + correction[i];
    const double absolute = absolute_tolerances
                                ? absolute_tolerances[algebraic] : 1.0e-3;
    const double scale = absolute + relative_tolerance *
                                      fmax(fabs(old_value), fabs(new_value));
    const double value = scale > 0.0 ? fabs(correction[i]) / scale
                                     : (correction[i] == 0.0 ? 0.0 : INFINITY);
    if (!isfinite(value))
      return INFINITY;
    if (value > maximum)
      maximum = value;
  }
  return maximum;
}

int chem_constraint_set_project(
    const chem_constraint_set *set, size_t cell_count, double *concentrations,
    const double *equilibrium_constants, const double *temperatures,
    const double *absolute_tolerances, double relative_tolerance,
    double correction_tolerance, size_t maximum_iterations,
    size_t maximum_backtracks, double backtrack_factor,
    double sufficient_decrease, size_t *iterations)
{
  size_t cell, i, j, total_iterations = 0, state_size;
  size_t m, n, matrix_size, work_count;
  double *work = NULL, *backup = NULL, *candidate = NULL;
  if (iterations) *iterations = 0;
  if (!set || (cell_count && set->species_count > SIZE_MAX / cell_count) ||
      (cell_count && set->constraint_count && !concentrations) ||
      !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
      !isfinite(correction_tolerance) || correction_tolerance < 0.0 ||
      !isfinite(backtrack_factor) || backtrack_factor <= 0.0 || backtrack_factor >= 1.0 ||
      !isfinite(sufficient_decrease) || sufficient_decrease < 0.0 || sufficient_decrease >= 1.0)
    return -1;
  m = set->constraint_count;
  n = set->species_count;
  if (!m) {
    if (iterations) *iterations = 0;
    return 0;
  }
  if (!cell_count) {
    if (iterations) *iterations = 0;
    return 0;
  }
  if (!set->constraints || !n || m > n || m > SIZE_MAX / m ||
      n > SIZE_MAX / sizeof(double) || cell_count > SIZE_MAX / n ||
      cell_count * n > SIZE_MAX / sizeof(double))
    return -1;
  state_size = cell_count * n;
  for (i = 0; i < m; ++i) {
    if (set->constraints[i].algebraic_species >= n)
      return -1;
    for (j = 0; j < i; ++j)
      if (set->constraints[i].algebraic_species == set->constraints[j].algebraic_species)
        return -1;
  }
  if (!maximum_iterations)
    maximum_iterations = 10;
  matrix_size = m * m;
  if (m > (SIZE_MAX - matrix_size) / 4 ||
      n > SIZE_MAX - matrix_size - 4 * m ||
      matrix_size + 4 * m + n > SIZE_MAX / sizeof(double))
    return -1;
  work_count = matrix_size + 4 * m + n;
  work = (double *)malloc(work_count * sizeof(*work));
  backup = (double *)malloc(state_size * sizeof(*backup));
  candidate = (double *)malloc(n * sizeof(*candidate));
  if (!work || !backup || !candidate) {
    free(work); free(backup); free(candidate);
    return -1;
  }
  memcpy(backup, concentrations, state_size * sizeof(*backup));
  for (cell = 0; cell < cell_count; ++cell) {
    double *y = concentrations + cell * n;
    double *matrix = work;
    double *rhs = matrix + matrix_size;
    double *residual = rhs + m;
    double *trial_residual = residual + m;
    double *trial_correction = trial_residual + m;
    double *gradient = trial_correction + m;
    size_t update;
    int cell_done = 0;
    for (update = 0; update <= maximum_iterations; ++update) {
      int residual_is_zero = 1;
      ++total_iterations;
      if (constraint_set_evaluate(set, y, cell, equilibrium_constants,
                                  temperatures, residual, matrix, gradient) != 0)
        break;
      for (i = 0; i < m; ++i) {
        if (residual[i] != 0.0)
          residual_is_zero = 0;
        rhs[i] = -residual[i];
      }
      if (residual_is_zero) {
        cell_done = 1;
        break;
      }
      if (chem_lu_factor_dense(matrix, m) != 0 ||
          chem_lu_solve_dense(matrix, m, rhs) != 0)
        break;
      for (i = 0; i < m; ++i)
        if (!isfinite(rhs[i]))
          break;
      if (i != m)
        break;
      {
        const double current_norm = constraint_correction_norm(
            set, y, rhs, absolute_tolerances, relative_tolerance);
        if (current_norm <= correction_tolerance) {
          cell_done = 1;
          break;
        }
        if (update == maximum_iterations)
          break;
        if (maximum_backtracks == 0) {
          for (i = 0; i < m; ++i)
            y[set->constraints[i].algebraic_species] += rhs[i];
          continue;
        }
        {
          double step = 1.0;
          size_t backtrack;
          int accepted = 0;
          for (backtrack = 0; backtrack <= maximum_backtracks; ++backtrack) {
            int trial_zero = 1;
            memcpy(candidate, y, n * sizeof(*candidate));
            for (i = 0; i < m; ++i)
              candidate[set->constraints[i].algebraic_species] += step * rhs[i];
            if (constraint_set_evaluate(set, candidate, cell,
                                        equilibrium_constants, temperatures,
                                        trial_residual, NULL, gradient) == 0) {
              for (i = 0; i < m; ++i) {
                if (trial_residual[i] != 0.0)
                  trial_zero = 0;
                trial_correction[i] = -trial_residual[i];
              }
              if (trial_zero) {
                accepted = 1;
                break;
              }
              if (chem_lu_solve_dense(matrix, m, trial_correction) == 0) {
                const double trial_norm = constraint_correction_norm(
                    set, candidate, trial_correction, absolute_tolerances,
                    relative_tolerance);
                if (trial_norm <= correction_tolerance ||
                    trial_norm < (1.0 - sufficient_decrease * step) * current_norm) {
                  accepted = 1;
                  break;
                }
              }
            }
            step *= backtrack_factor;
          }
          if (!accepted)
            break;
          memcpy(y, candidate, n * sizeof(*y));
        }
      }
    }
    if (!cell_done) {
      memcpy(concentrations, backup, state_size * sizeof(*backup));
      free(work); free(backup); free(candidate);
      if (iterations) *iterations = total_iterations;
      return -1;
    }
  }
  free(work); free(backup); free(candidate);
  if (iterations) *iterations = total_iterations;
  return 0;
}
