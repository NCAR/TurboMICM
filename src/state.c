#include "chemistry.h"

#include <string.h>

static size_t chem_state_find(const char *const *labels, size_t count,
                              const char *label) {
  if (label == NULL) return (size_t)-1;
  for (size_t i = 0; i < count; ++i)
    if (labels[i] != NULL && strcmp(labels[i], label) == 0) return i;
  return (size_t)-1;
}

static int chem_state_values_valid(const chem_state_view *state,
                                   const double *values, size_t value_count) {
  return state != NULL && values != NULL && value_count == state->cell_count;
}

int chem_state_view_init(chem_state_view *state,
                         const char *const *species_names,
                         size_t species_count,
                         const char *const *custom_parameter_labels,
                         size_t custom_parameter_count, size_t cell_count,
                         double *concentrations, double *custom_parameters) {
  if (state == NULL || (species_count && species_names == NULL) ||
      (custom_parameter_count && custom_parameter_labels == NULL) ||
      (cell_count && species_count && concentrations == NULL) ||
      (cell_count && custom_parameter_count && custom_parameters == NULL))
    return 1;
  for (size_t i = 0; i < species_count; ++i)
    if (species_names[i] == NULL) return 1;
  for (size_t i = 0; i < custom_parameter_count; ++i)
    if (custom_parameter_labels[i] == NULL) return 1;
  *state = (chem_state_view){
      .species_names = species_names,
      .species_count = species_count,
      .custom_parameter_labels = custom_parameter_labels,
      .custom_parameter_count = custom_parameter_count,
      .cell_count = cell_count,
      .concentrations = concentrations,
      .custom_parameters = custom_parameters};
  return 0;
}

int chem_state_set_concentration(chem_state_view *state, const char *species,
                                 const double *values, size_t value_count) {
  if (state == NULL || !chem_state_values_valid(state, values, value_count) ||
      (state->cell_count && state->species_count && state->concentrations == NULL))
    return 1;
  const size_t column = chem_state_find(state->species_names,
                                        state->species_count, species);
  if (column == (size_t)-1) return 1;
  for (size_t cell = 0; cell < state->cell_count; ++cell)
    state->concentrations[cell * state->species_count + column] = values[cell];
  return 0;
}

int chem_state_get_concentration(const chem_state_view *state,
                                 const char *species, size_t cell,
                                 double *value) {
  if (state == NULL || value == NULL || cell >= state->cell_count ||
      state->concentrations == NULL) return 1;
  const size_t column = chem_state_find(state->species_names,
                                        state->species_count, species);
  if (column == (size_t)-1) return 1;
  *value = state->concentrations[cell * state->species_count + column];
  return 0;
}

int chem_state_set_concentrations(chem_state_view *state,
                                  const chem_state_series *series,
                                  size_t series_count) {
  if (state == NULL || (series_count && series == NULL)) return 1;
  for (size_t i = 0; i < series_count; ++i)
    if (chem_state_find(state->species_names, state->species_count,
                        series[i].label) == (size_t)-1 ||
        !chem_state_values_valid(state, series[i].values, series[i].value_count))
      return 1;
  for (size_t i = 0; i < series_count; ++i)
    if (chem_state_set_concentration(state, series[i].label, series[i].values,
                                     series[i].value_count)) return 1;
  return 0;
}

int chem_state_set_custom_parameter(chem_state_view *state, const char *label,
                                    const double *values,
                                    size_t value_count) {
  if (state == NULL || !chem_state_values_valid(state, values, value_count) ||
      (state->cell_count && state->custom_parameter_count &&
       state->custom_parameters == NULL)) return 1;
  const size_t column = chem_state_find(state->custom_parameter_labels,
                                        state->custom_parameter_count, label);
  if (column == (size_t)-1) return 1;
  for (size_t cell = 0; cell < state->cell_count; ++cell)
    state->custom_parameters[cell * state->custom_parameter_count + column] = values[cell];
  return 0;
}

int chem_state_set_custom_parameters(chem_state_view *state,
                                     const chem_state_series *series,
                                     size_t series_count) {
  if (state == NULL || (series_count && series == NULL)) return 1;
  for (size_t i = 0; i < series_count; ++i)
    if (chem_state_find(state->custom_parameter_labels,
                        state->custom_parameter_count, series[i].label) == (size_t)-1 ||
        !chem_state_values_valid(state, series[i].values, series[i].value_count))
      return 1;
  for (size_t i = 0; i < series_count; ++i)
    if (chem_state_set_custom_parameter(state, series[i].label, series[i].values,
                                        series[i].value_count)) return 1;
  return 0;
}

int chem_state_set_custom_parameter_matrix(chem_state_view *state,
                                           const double *values,
                                           size_t row_count,
                                           size_t column_count) {
  if (state == NULL || row_count != state->cell_count ||
      column_count != state->custom_parameter_count ||
      (row_count && column_count && (values == NULL || state->custom_parameters == NULL)))
    return 1;
  if (row_count && column_count)
    memcpy(state->custom_parameters, values,
           row_count * column_count * sizeof(*values));
  return 0;
}
