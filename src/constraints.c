#include "chem_constraints.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
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

static int constraint_size_product(size_t left, size_t right, size_t *product) {
  if (product == NULL || (right != 0u && left > SIZE_MAX / right)) return -1;
  *product = left * right;
  return 0;
}

static int constraint_lu_factor(double *matrix, size_t n, size_t *pivots) {
  for (size_t k = 0; k < n; ++k) {
    size_t pivot = k;
    double largest = fabs(matrix[k * n + k]);
    for (size_t row = k + 1u; row < n; ++row) {
      const double candidate = fabs(matrix[row * n + k]);
      if (candidate > largest) {
        largest = candidate;
        pivot = row;
      }
    }
    if (!(largest > 0.0) || !isfinite(largest)) return -1;
    pivots[k] = pivot;
    if (pivot != k) {
      for (size_t column = 0; column < n; ++column) {
        const double swap = matrix[k * n + column];
        matrix[k * n + column] = matrix[pivot * n + column];
        matrix[pivot * n + column] = swap;
      }
    }
    for (size_t row = k + 1u; row < n; ++row) {
      matrix[row * n + k] /= matrix[k * n + k];
      for (size_t column = k + 1u; column < n; ++column)
        matrix[row * n + column] -= matrix[row * n + k] * matrix[k * n + column];
    }
  }
  return 0;
}

static int constraint_lu_solve(
    const double *matrix, size_t n, const size_t *pivots, double *right_hand_side) {
  for (size_t k = 0; k < n; ++k) {
    if (pivots[k] != k) {
      const double swap = right_hand_side[k];
      right_hand_side[k] = right_hand_side[pivots[k]];
      right_hand_side[pivots[k]] = swap;
    }
  }
  for (size_t row = 0; row < n; ++row)
    for (size_t column = 0; column < row; ++column)
      right_hand_side[row] -= matrix[row * n + column] * right_hand_side[column];
  for (size_t row = n; row-- > 0u;) {
    for (size_t column = row + 1u; column < n; ++column)
      right_hand_side[row] -= matrix[row * n + column] * right_hand_side[column];
    right_hand_side[row] /= matrix[row * n + row];
    if (!isfinite(right_hand_side[row])) return -1;
  }
  return 0;
}

static double constraint_weighted_correction(
    const chem_constraint_set *set,
    size_t cell_count,
    const size_t *algebraic_species,
    const double *state,
    const double *delta,
    const double *absolute_tolerances,
    double relative_tolerance) {
  double maximum = 0.0;
  const size_t m = set->constraint_count;
  for (size_t cell = 0; cell < cell_count; ++cell) {
    for (size_t i = 0; i < m; ++i) {
      const size_t species = algebraic_species[i];
      const double value = state[cell * set->species_count + species];
      const double correction = delta[cell * m + i];
      const double corrected = value + correction;
      const double scale = absolute_tolerances[species] + relative_tolerance * fmax(fabs(value), fabs(corrected));
      const double ratio = scale > 0.0 ? fabs(correction) / scale :
          (correction == 0.0 ? 0.0 : INFINITY);
      if (!isfinite(ratio)) return ratio;
      maximum = fmax(maximum, ratio);
    }
  }
  return maximum;
}

static int constraint_residuals_finite(
    const chem_constraint_set *set, size_t cell_count, const size_t *algebraic_species, const double *forcing) {
  for (size_t cell = 0; cell < cell_count; ++cell)
    for (size_t i = 0; i < set->constraint_count; ++i)
      if (!isfinite(forcing[cell * set->species_count + algebraic_species[i]])) return 0;
  return 1;
}

static int constraint_state_finite(
    const chem_constraint_set *set, size_t cell_count, const size_t *algebraic_species, const double *state) {
  for (size_t cell = 0; cell < cell_count; ++cell)
    for (size_t i = 0; i < set->constraint_count; ++i)
      if (!isfinite(state[cell * set->species_count + algebraic_species[i]])) return 0;
  return 1;
}

chem_constraint_init_status chem_constraint_set_initialize(
    const chem_constraint_set *set,
    size_t cell_count,
    const double *state_parameters,
    double *state,
    const double *absolute_tolerances,
    double relative_tolerance,
    size_t max_updates,
    double correction_tolerance,
    size_t max_backtracks,
    double backtrack_factor,
    double sufficient_decrease,
    size_t *iterations) {
  if (iterations != NULL) *iterations = 0u;
  if (chem_constraint_set_validate(set) != 0 || cell_count == 0u || state == NULL || iterations == NULL ||
      set->species_count == 0u || max_updates == SIZE_MAX || !isfinite(relative_tolerance) ||
      relative_tolerance < 0.0 || !isfinite(correction_tolerance) || correction_tolerance < 0.0 ||
      !isfinite(backtrack_factor) || backtrack_factor <= 0.0 || backtrack_factor >= 1.0 ||
      !isfinite(sufficient_decrease) || sufficient_decrease < 0.0 || sufficient_decrease >= 1.0 ||
      max_backtracks == SIZE_MAX ||
      (set->state_parameter_count != 0u && state_parameters == NULL) ||
      (set->constraint_count != 0u && absolute_tolerances == NULL))
    return CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT;
  if (set->constraint_count == 0u) return CHEM_CONSTRAINT_INIT_CONVERGED;

  const size_t n = set->species_count;
  const size_t m = set->constraint_count;
  size_t state_size, dense_width, dense_size, small_width, small_size, delta_size, pivot_size;
  if (constraint_size_product(n, cell_count, &state_size) != 0 ||
      constraint_size_product(n, n, &dense_width) != 0 ||
      constraint_size_product(dense_width, cell_count, &dense_size) != 0 ||
      constraint_size_product(m, m, &small_width) != 0 ||
      constraint_size_product(small_width, cell_count, &small_size) != 0 ||
      constraint_size_product(m, cell_count, &delta_size) != 0 ||
      constraint_size_product(m, cell_count, &pivot_size) != 0 ||
      state_size > SIZE_MAX / sizeof(double) || dense_size > SIZE_MAX / sizeof(double) ||
      small_size > SIZE_MAX / sizeof(double) || delta_size > SIZE_MAX / sizeof(double) ||
      pivot_size > SIZE_MAX / sizeof(size_t) ||
      set->state_parameter_count > SIZE_MAX / cell_count / sizeof(double))
    return CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT;
  for (size_t i = 0; i < n; ++i)
    if (!isfinite(absolute_tolerances[i]) || absolute_tolerances[i] < 0.0)
      return CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT;

  size_t *algebraic_species = malloc(m * sizeof(*algebraic_species));
  size_t *pivots = malloc(pivot_size * sizeof(*pivots));
  double *original = malloc(state_size * sizeof(*original));
  double *forcing = calloc(state_size, sizeof(*forcing));
  double *trial_forcing = calloc(state_size, sizeof(*trial_forcing));
  double *negative_jacobian = calloc(dense_size, sizeof(*negative_jacobian));
  double *factor = malloc(small_size * sizeof(*factor));
  double *delta = malloc(delta_size * sizeof(*delta));
  double *trial_delta = malloc(delta_size * sizeof(*trial_delta));
  double *candidate = malloc(state_size * sizeof(*candidate));
  if (algebraic_species == NULL || pivots == NULL || original == NULL || forcing == NULL ||
      trial_forcing == NULL || negative_jacobian == NULL || factor == NULL || delta == NULL ||
      trial_delta == NULL || candidate == NULL) {
    free(algebraic_species); free(pivots); free(original); free(forcing); free(trial_forcing);
    free(negative_jacobian); free(factor); free(delta); free(trial_delta); free(candidate);
    return CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT;
  }
  for (size_t i = 0; i < m; ++i) algebraic_species[i] = set->constraints[i].algebraic_species;
  memcpy(original, state, state_size * sizeof(*state));

  chem_constraint_init_status status = CHEM_CONSTRAINT_INIT_FAILED;
  size_t updates = 0u;
  for (;;) {
    memset(forcing, 0, state_size * sizeof(*forcing));
    if (chem_constraint_set_apply_forcing(set, cell_count, state, state_parameters, forcing) != 0) {
      status = CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT;
      break;
    }
    if (!constraint_residuals_finite(set, cell_count, algebraic_species, forcing)) {
      status = CHEM_CONSTRAINT_INIT_INFINITY;
      for (size_t cell = 0; cell < cell_count; ++cell)
        for (size_t i = 0; i < m; ++i)
          if (isnan(forcing[cell * n + algebraic_species[i]])) status = CHEM_CONSTRAINT_INIT_NAN;
      break;
    }
    if (iterations != NULL) ++*iterations;

    int residual_is_zero = 1;
    for (size_t cell = 0; cell < cell_count; ++cell)
      for (size_t i = 0; i < m; ++i)
        residual_is_zero &= forcing[cell * n + algebraic_species[i]] == 0.0;
    if (residual_is_zero) {
      status = CHEM_CONSTRAINT_INIT_CONVERGED;
      break;
    }

    memset(negative_jacobian, 0, dense_size * sizeof(*negative_jacobian));
    if (chem_constraint_set_accumulate_negative_jacobian(
            set, cell_count, state, state_parameters, negative_jacobian) != 0) {
      status = CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT;
      break;
    }
    int factor_failed = 0;
    for (size_t cell = 0; cell < cell_count && !factor_failed; ++cell) {
      const double *cell_matrix = negative_jacobian + cell * dense_width;
      double *cell_factor = factor + cell * small_width;
      double *cell_delta = delta + cell * m;
      size_t *cell_pivots = pivots + cell * m;
      for (size_t row = 0; row < m; ++row) {
        const size_t source_row = algebraic_species[row];
        cell_delta[row] = forcing[cell * n + source_row];
        for (size_t column = 0; column < m; ++column)
          cell_factor[row * m + column] = cell_matrix[source_row * n + algebraic_species[column]];
      }
      if (constraint_lu_factor(cell_factor, m, cell_pivots) != 0 ||
          constraint_lu_solve(cell_factor, m, cell_pivots, cell_delta) != 0)
        factor_failed = 1;
    }
    if (factor_failed) {
      status = CHEM_CONSTRAINT_INIT_FAILED;
      break;
    }
    if (!constraint_state_finite(set, cell_count, algebraic_species, state)) {
      status = CHEM_CONSTRAINT_INIT_INFINITY;
      break;
    }
    const double current_norm = constraint_weighted_correction(
        set, cell_count, algebraic_species, state, delta, absolute_tolerances, relative_tolerance);
    if (isnan(current_norm)) {
      status = CHEM_CONSTRAINT_INIT_NAN;
      break;
    }
    if (isinf(current_norm)) {
      status = CHEM_CONSTRAINT_INIT_INFINITY;
      break;
    }
    if (current_norm <= correction_tolerance) {
      status = CHEM_CONSTRAINT_INIT_CONVERGED;
      break;
    }
    if (updates >= max_updates) break;

    if (max_backtracks == 0u) {
      for (size_t cell = 0; cell < cell_count; ++cell)
        for (size_t i = 0; i < m; ++i)
          state[cell * n + algebraic_species[i]] += delta[cell * m + i];
      if (!constraint_state_finite(set, cell_count, algebraic_species, state)) {
        status = CHEM_CONSTRAINT_INIT_INFINITY;
        break;
      }
      ++updates;
      continue;
    }

    int accepted = 0;
    double step = 1.0;
    for (size_t backtrack = 0; backtrack <= max_backtracks; ++backtrack, step *= backtrack_factor) {
      memcpy(candidate, state, state_size * sizeof(*state));
      for (size_t cell = 0; cell < cell_count; ++cell)
        for (size_t i = 0; i < m; ++i)
          candidate[cell * n + algebraic_species[i]] += step * delta[cell * m + i];
      if (!constraint_state_finite(set, cell_count, algebraic_species, candidate)) continue;

      memset(trial_forcing, 0, state_size * sizeof(*trial_forcing));
      if (chem_constraint_set_apply_forcing(set, cell_count, candidate, state_parameters, trial_forcing) != 0) {
        status = CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT;
        accepted = -1;
        break;
      }
      if (!constraint_residuals_finite(set, cell_count, algebraic_species, trial_forcing)) continue;
      int trial_solve_failed = 0;
      for (size_t cell = 0; cell < cell_count; ++cell) {
        double *cell_delta = trial_delta + cell * m;
        for (size_t i = 0; i < m; ++i)
          cell_delta[i] = trial_forcing[cell * n + algebraic_species[i]];
        if (constraint_lu_solve(factor + cell * small_width, m, pivots + cell * m, cell_delta) != 0) {
          trial_solve_failed = 1;
          break;
        }
      }
      if (accepted < 0) break;
      if (trial_solve_failed) continue;
      const double trial_norm = constraint_weighted_correction(
          set, cell_count, algebraic_species, candidate, trial_delta, absolute_tolerances, relative_tolerance);
      if (isnan(trial_norm) || isinf(trial_norm)) continue;
      if (trial_norm <= correction_tolerance ||
          trial_norm < (1.0 - sufficient_decrease * step) * current_norm) {
        memcpy(state, candidate, state_size * sizeof(*state));
        accepted = 1;
        break;
      }
    }
    if (accepted < 0) break;
    if (!accepted) {
      status = CHEM_CONSTRAINT_INIT_FAILED;
      break;
    }
    ++updates;
  }

  if (status != CHEM_CONSTRAINT_INIT_CONVERGED) memcpy(state, original, state_size * sizeof(*state));
  free(algebraic_species); free(pivots); free(original); free(forcing); free(trial_forcing);
  free(negative_jacobian); free(factor); free(delta); free(trial_delta); free(candidate);
  return status;
}

const char *chem_constraint_init_status_name(chem_constraint_init_status status) {
  switch (status) {
    case CHEM_CONSTRAINT_INIT_CONVERGED: return "Converged";
    case CHEM_CONSTRAINT_INIT_FAILED: return "Constraint Initialization Failed";
    case CHEM_CONSTRAINT_INIT_NAN: return "NaN Detected";
    case CHEM_CONSTRAINT_INIT_INFINITY: return "Inf Detected";
    case CHEM_CONSTRAINT_INIT_INVALID_ARGUMENT: return "Invalid Argument";
    default: return "Unknown";
  }
}
