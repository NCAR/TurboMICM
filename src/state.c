#include "kinetics.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int chem_state_init(chem_state *state, const chem_mechanism *mechanism,
                    size_t cell_count)
{
  size_t concentration_count, custom_count, cell, species, parameter;
  if (!state || !mechanism || !mechanism->species || !mechanism->species_count)
    return -1;
  *state = (chem_state){0};
  if (cell_count && (mechanism->species_count > SIZE_MAX / cell_count ||
                     mechanism->custom_count > SIZE_MAX / cell_count ||
                     cell_count > SIZE_MAX / sizeof(double) ||
                     mechanism->species_count > SIZE_MAX / sizeof(double) / cell_count ||
                     mechanism->custom_count > SIZE_MAX / sizeof(double) / cell_count))
    return -1;
  concentration_count = cell_count * mechanism->species_count;
  custom_count = cell_count * mechanism->custom_count;
  state->mechanism = mechanism;
  state->cells = cell_count;
  state->concentrations = (double *)calloc(concentration_count ? concentration_count : 1,
                                           sizeof(*state->concentrations));
  state->temperature = (double *)malloc((cell_count ? cell_count : 1) * sizeof(*state->temperature));
  state->pressure = (double *)malloc((cell_count ? cell_count : 1) * sizeof(*state->pressure));
  state->custom_parameters = (double *)calloc(custom_count ? custom_count : 1,
                                               sizeof(*state->custom_parameters));
  if (!state->concentrations || !state->temperature || !state->pressure ||
      !state->custom_parameters) {
    chem_state_destroy(state);
    return -1;
  }
  for (cell = 0; cell < cell_count; ++cell) {
    state->temperature[cell] = mechanism->default_temperature;
    state->pressure[cell] = mechanism->default_pressure;
    for (species = 0; species < mechanism->species_count; ++species)
      state->concentrations[cell * mechanism->species_count + species] =
          mechanism->species[species].initial_value;
    for (parameter = 0; parameter < mechanism->custom_count; ++parameter)
      state->custom_parameters[cell * mechanism->custom_count + parameter] =
          mechanism->custom_parameters[parameter].default_value;
  }
  return 0;
}

int chem_state_copy(chem_state *destination, const chem_state *source)
{
  chem_state copied = {0};
  size_t concentration_count, custom_count;
  if (!destination || !source || !source->mechanism ||
      (source->cells && (!source->concentrations || !source->temperature ||
                         !source->pressure || !source->custom_parameters)))
    return -1;
  if (chem_state_init(&copied, source->mechanism, source->cells) != 0)
    return -1;
  concentration_count = source->cells * source->mechanism->species_count;
  custom_count = source->cells * source->mechanism->custom_count;
  if (concentration_count)
    memcpy(copied.concentrations, source->concentrations,
           concentration_count * sizeof(*source->concentrations));
  if (source->cells) {
    memcpy(copied.temperature, source->temperature,
           source->cells * sizeof(*source->temperature));
    memcpy(copied.pressure, source->pressure,
           source->cells * sizeof(*source->pressure));
  }
  if (custom_count)
    memcpy(copied.custom_parameters, source->custom_parameters,
           custom_count * sizeof(*source->custom_parameters));
  chem_state_destroy(destination);
  *destination = copied;
  return 0;
}

void chem_state_move(chem_state *destination, chem_state *source)
{
  if (!destination || !source || destination == source)
    return;
  chem_state_destroy(destination);
  *destination = *source;
  *source = (chem_state){0};
}

void chem_state_destroy(chem_state *state)
{
  if (!state)
    return;
  free(state->concentrations);
  free(state->temperature);
  free(state->pressure);
  free(state->custom_parameters);
  *state = (chem_state){0};
}

size_t chem_state_species_index(const chem_state *state, const char *name)
{
  size_t i;
  if (!state || !state->mechanism || !name)
    return SIZE_MAX;
  for (i = 0; i < state->mechanism->species_count; ++i)
    if (state->mechanism->species[i].name &&
        strcmp(state->mechanism->species[i].name, name) == 0)
      return i;
  return SIZE_MAX;
}

size_t chem_state_custom_parameter_index(const chem_state *state,
                                         const char *name)
{
  size_t i;
  if (!state || !state->mechanism || !name)
    return SIZE_MAX;
  for (i = 0; i < state->mechanism->custom_count; ++i)
    if (state->mechanism->custom_parameters &&
        state->mechanism->custom_parameters[i].name &&
        strcmp(state->mechanism->custom_parameters[i].name, name) == 0)
      return i;
  return SIZE_MAX;
}

int chem_state_set_concentration(chem_state *state, const char *species,
                                 size_t cell, double value)
{
  const size_t index = chem_state_species_index(state, species);
  if (index == SIZE_MAX || cell >= state->cells || !state->concentrations)
    return -1;
  state->concentrations[cell * state->mechanism->species_count + index] = value;
  return 0;
}

int chem_state_get_concentration(const chem_state *state, const char *species,
                                 size_t cell, double *value)
{
  const size_t index = chem_state_species_index(state, species);
  if (index == SIZE_MAX || cell >= state->cells || !state->concentrations || !value)
    return -1;
  *value = state->concentrations[cell * state->mechanism->species_count + index];
  return 0;
}

int chem_state_set_species_values(chem_state *state, const char *species,
                                 const double *values, size_t value_count)
{
  const size_t index = chem_state_species_index(state, species);
  size_t cell;
  if (index == SIZE_MAX || value_count != state->cells ||
      (value_count && !values))
    return -1;
  for (cell = 0; cell < state->cells; ++cell)
    state->concentrations[cell * state->mechanism->species_count + index] = values[cell];
  return 0;
}

int chem_state_set_cell_values(chem_state *state, size_t cell,
                               const double *values, size_t value_count)
{
  if (!state || !state->mechanism || cell >= state->cells || !values ||
      value_count != state->mechanism->species_count)
    return -1;
  memcpy(state->concentrations + cell * value_count, values,
         value_count * sizeof(*values));
  return 0;
}

int chem_state_set_custom_parameter(chem_state *state, const char *name,
                                    size_t cell, double value)
{
  const size_t index = chem_state_custom_parameter_index(state, name);
  if (!state || index == SIZE_MAX || cell >= state->cells ||
      !state->custom_parameters)
    return -1;
  state->custom_parameters[cell * state->mechanism->custom_count + index] = value;
  return 0;
}

int chem_state_set_custom_parameter_values(chem_state *state, const char *name,
                                           const double *values,
                                           size_t value_count)
{
  const size_t index = chem_state_custom_parameter_index(state, name);
  size_t cell;
  if (index == SIZE_MAX || value_count != state->cells ||
      (value_count && !values))
    return -1;
  for (cell = 0; cell < state->cells; ++cell)
    state->custom_parameters[cell * state->mechanism->custom_count + index] = values[cell];
  return 0;
}

int chem_state_set_custom_parameters(chem_state *state, const double *values,
                                     size_t cell_count,
                                     size_t parameter_count)
{
  size_t count;
  if (!state || !state->mechanism || cell_count != state->cells ||
      parameter_count != state->mechanism->custom_count ||
      (cell_count && parameter_count && !values))
    return -1;
  if (cell_count && parameter_count > SIZE_MAX / cell_count)
    return -1;
  count = cell_count * parameter_count;
  if (count)
    memcpy(state->custom_parameters, values, count * sizeof(*values));
  return 0;
}
