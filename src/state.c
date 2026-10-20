#include "chem.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int multiply_size(size_t left, size_t right, size_t *result) {
  if (right != 0u && left > SIZE_MAX / right) return -1;
  *result = left * right;
  return 0;
}

static char *copy_name(const char *name) {
  const size_t length = strlen(name);
  if (length == SIZE_MAX) return NULL;
  char *copy = malloc(length + 1u);
  if (copy != NULL) memcpy(copy, name, length + 1u);
  return copy;
}

static void release_names(char **names, size_t count) {
  if (names == NULL) return;
  for (size_t i = 0; i < count; ++i) free(names[i]);
  free(names);
}

static int copy_names(char ***destination, const char *const *source, size_t count) {
  if (count == 0u) {
    *destination = NULL;
    return 0;
  }
  if (source == NULL || count > SIZE_MAX / sizeof(char *)) return -1;
  char **names = calloc(count, sizeof(*names));
  if (names == NULL) return -2;
  for (size_t i = 0; i < count; ++i) {
    if (source[i] == NULL || source[i][0] == '\0') {
      release_names(names, count);
      return -1;
    }
    names[i] = copy_name(source[i]);
    if (names[i] == NULL) {
      release_names(names, count);
      return -2;
    }
    for (size_t j = 0; j < i; ++j) {
      if (strcmp(names[i], names[j]) == 0) {
        release_names(names, count);
        return -1;
      }
    }
  }
  *destination = names;
  return 0;
}

static size_t find_name(char *const *names, size_t count, const char *name) {
  if (name == NULL) return count;
  for (size_t i = 0; i < count; ++i)
    if (strcmp(names[i], name) == 0) return i;
  return count;
}

void chem_state_destroy(chem_state *state) {
  if (state == NULL) return;
  release_names(state->species_names, state->species_count);
  release_names(state->parameter_names, state->parameter_count);
  free(state->concentrations);
  free(state->custom_parameters);
  free(state->rate_constants);
  free(state->absolute_tolerances);
  memset(state, 0, sizeof(*state));
}

int chem_state_init(
    chem_state *state,
    const char *const *species_names,
    size_t species_count,
    size_t rate_constant_count,
    const char *const *parameter_names,
    size_t parameter_count,
    size_t cell_count,
    double relative_tolerance,
    const double *absolute_tolerances) {
  if (state == NULL || (species_count != 0u && species_names == NULL) ||
      (parameter_count != 0u && parameter_names == NULL)) return -1;

  chem_state result = { 0 };
  result.cell_count = cell_count;
  result.species_count = species_count;
  result.parameter_count = parameter_count;
  result.rate_constant_count = rate_constant_count;
  result.relative_tolerance = relative_tolerance;

  int status = copy_names(&result.species_names, species_names, species_count);
  if (status != 0) goto fail;
  status = copy_names(&result.parameter_names, parameter_names, parameter_count);
  if (status != 0) goto fail;

  size_t cells_species, cells_parameters, cells_rates;
  if (multiply_size(cell_count, species_count, &cells_species) != 0 ||
      multiply_size(cell_count, parameter_count, &cells_parameters) != 0 ||
      multiply_size(cell_count, rate_constant_count, &cells_rates) != 0 ||
      cells_species > SIZE_MAX / sizeof(double) ||
      cells_parameters > SIZE_MAX / sizeof(double) ||
      cells_rates > SIZE_MAX / sizeof(double) ||
      species_count > SIZE_MAX / sizeof(double)) {
    status = -1;
    goto fail;
  }
  if (cells_species != 0u && (result.concentrations = calloc(cells_species, sizeof(double))) == NULL) {
    status = -2;
    goto fail;
  }
  if (cells_parameters != 0u && (result.custom_parameters = calloc(cells_parameters, sizeof(double))) == NULL) {
    status = -2;
    goto fail;
  }
  if (cells_rates != 0u && (result.rate_constants = calloc(cells_rates, sizeof(double))) == NULL) {
    status = -2;
    goto fail;
  }
  if (species_count != 0u) {
    result.absolute_tolerances = malloc(species_count * sizeof(double));
    if (result.absolute_tolerances == NULL) {
      status = -2;
      goto fail;
    }
    for (size_t i = 0; i < species_count; ++i)
      result.absolute_tolerances[i] = absolute_tolerances == NULL ? 0.0 : absolute_tolerances[i];
  }

  chem_state_destroy(state);
  *state = result;
  return 0;

fail:
  chem_state_destroy(&result);
  return status;
}

int chem_state_copy(chem_state *destination, const chem_state *source) {
  if (destination == NULL || source == NULL || destination == source) return -1;
  chem_state copy = { 0 };
  const int status = chem_state_init(&copy,
      (const char *const *)source->species_names, source->species_count,
      source->rate_constant_count, (const char *const *)source->parameter_names,
      source->parameter_count, source->cell_count, source->relative_tolerance,
      source->absolute_tolerances);
  if (status != 0) return status;
  size_t n = 0u;
  if (multiply_size(source->cell_count, source->species_count, &n) != 0) {
    chem_state_destroy(&copy);
    return -1;
  }
  if (n != 0u) memcpy(copy.concentrations, source->concentrations, n * sizeof(double));
  if (multiply_size(source->cell_count, source->parameter_count, &n) != 0) {
    chem_state_destroy(&copy);
    return -1;
  }
  if (n != 0u) memcpy(copy.custom_parameters, source->custom_parameters, n * sizeof(double));
  if (multiply_size(source->cell_count, source->rate_constant_count, &n) != 0) {
    chem_state_destroy(&copy);
    return -1;
  }
  if (n != 0u) memcpy(copy.rate_constants, source->rate_constants, n * sizeof(double));
  chem_state_destroy(destination);
  *destination = copy;
  return 0;
}

void chem_state_move(chem_state *destination, chem_state *source) {
  if (destination == NULL || source == NULL || destination == source) return;
  chem_state_destroy(destination);
  *destination = *source;
  memset(source, 0, sizeof(*source));
}

int chem_state_set_concentrations(chem_state *state, const char *species, const double *values, size_t value_count) {
  if (state == NULL || values == NULL || value_count != state->cell_count) return -1;
  const size_t index = find_name(state->species_names, state->species_count, species);
  if (index == state->species_count) return -1;
  for (size_t cell = 0; cell < state->cell_count; ++cell)
    state->concentrations[cell * state->species_count + index] = values[cell];
  return 0;
}

int chem_state_set_concentration(chem_state *state, const char *species, double value) {
  if (state == NULL || state->cell_count != 1u) return -1;
  return chem_state_set_concentrations(state, species, &value, 1u);
}

static int validate_named_values(
    char *const *names, size_t name_count, const chem_named_values *values, size_t count,
    size_t cell_count) {
  if ((count != 0u && values == NULL) || count > name_count) return -1;
  for (size_t i = 0; i < count; ++i) {
    if (values[i].values == NULL || values[i].value_count != cell_count ||
        find_name(names, name_count, values[i].name) == name_count) return -1;
    for (size_t j = 0; j < i; ++j)
      if (strcmp(values[i].name, values[j].name) == 0) return -1;
  }
  return 0;
}

int chem_state_set_all_concentrations(chem_state *state, const chem_named_values *values, size_t value_count) {
  if (state == NULL || validate_named_values(state->species_names, state->species_count,
          values, value_count, state->cell_count) != 0) return -1;
  for (size_t i = 0; i < value_count; ++i)
    if (chem_state_set_concentrations(state, values[i].name, values[i].values, values[i].value_count) != 0) return -1;
  return 0;
}

int chem_state_get_concentration(const chem_state *state, const char *species, size_t cell, double *value) {
  if (state == NULL || value == NULL || cell >= state->cell_count) return -1;
  const size_t index = find_name(state->species_names, state->species_count, species);
  if (index == state->species_count) return -1;
  *value = state->concentrations[cell * state->species_count + index];
  return 0;
}

int chem_state_set_custom_parameters(chem_state *state, const char *parameter, const double *values, size_t value_count) {
  if (state == NULL || values == NULL || value_count != state->cell_count) return -1;
  const size_t index = find_name(state->parameter_names, state->parameter_count, parameter);
  if (index == state->parameter_count) return -1;
  for (size_t cell = 0; cell < state->cell_count; ++cell)
    state->custom_parameters[cell * state->parameter_count + index] = values[cell];
  return 0;
}

int chem_state_set_custom_parameter(chem_state *state, const char *parameter, double value) {
  if (state == NULL || state->cell_count != 1u) return -1;
  return chem_state_set_custom_parameters(state, parameter, &value, 1u);
}

int chem_state_set_all_custom_parameters(chem_state *state, const chem_named_values *values, size_t value_count) {
  if (state == NULL || validate_named_values(state->parameter_names, state->parameter_count,
          values, value_count, state->cell_count) != 0) return -1;
  for (size_t i = 0; i < value_count; ++i)
    if (chem_state_set_custom_parameters(state, values[i].name, values[i].values, values[i].value_count) != 0) return -1;
  return 0;
}

int chem_state_set_custom_parameters_raw(chem_state *state, const double *values, size_t rows, size_t columns) {
  if (state == NULL || rows != state->cell_count || columns != state->parameter_count ||
      (rows != 0u && columns != 0u && values == NULL)) return -1;
  size_t count;
  if (multiply_size(rows, columns, &count) != 0) return -1;
  if (count != 0u) memcpy(state->custom_parameters, values, count * sizeof(double));
  return 0;
}
