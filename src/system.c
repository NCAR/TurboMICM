#include "system.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *duplicate_string(const char *source)
{
  const size_t length = strlen(source);
  char *copy;
  if (length == SIZE_MAX)
    return NULL;
  copy = (char *)malloc(length + 1);
  if (copy)
    memcpy(copy, source, length + 1);
  return copy;
}

int chem_phase_species_init(chem_phase_species *phase_species,
                            const chem_species *species,
                            int has_diffusion_coefficient,
                            double diffusion_coefficient)
{
  if (!phase_species || !species || !species->name ||
      (has_diffusion_coefficient != 0 && has_diffusion_coefficient != 1))
    return -1;
  *phase_species = (chem_phase_species){0};
  if (chem_species_copy(&phase_species->species, species) != 0)
    return -1;
  phase_species->has_diffusion_coefficient = has_diffusion_coefficient;
  phase_species->diffusion_coefficient =
      has_diffusion_coefficient ? diffusion_coefficient : 0.0;
  return 0;
}

int chem_phase_species_copy(chem_phase_species *destination,
                            const chem_phase_species *source)
{
  chem_phase_species copied = {0};
  if (!destination || !source ||
      chem_species_copy(&copied.species, &source->species) != 0)
    return -1;
  copied.diffusion_coefficient = source->diffusion_coefficient;
  copied.has_diffusion_coefficient = source->has_diffusion_coefficient;
  chem_phase_species_destroy(destination);
  *destination = copied;
  return 0;
}

void chem_phase_species_destroy(chem_phase_species *phase_species)
{
  if (!phase_species)
    return;
  chem_species_destroy(&phase_species->species);
  *phase_species = (chem_phase_species){0};
}

int chem_phase_species_set_diffusion_coefficient(
    chem_phase_species *phase_species, double diffusion_coefficient)
{
  if (!phase_species)
    return -1;
  phase_species->diffusion_coefficient = diffusion_coefficient;
  phase_species->has_diffusion_coefficient = 1;
  return 0;
}

void chem_phase_destroy(chem_phase *phase)
{
  size_t i;
  if (!phase)
    return;
  for (i = 0; i < phase->species_count; ++i)
    chem_phase_species_destroy(&phase->species[i]);
  free(phase->species);
  free(phase->name);
  *phase = (chem_phase){0};
}

int chem_phase_init(chem_phase *phase, const char *name,
                    const chem_phase_species *species, size_t species_count)
{
  chem_phase created = {0};
  size_t i;
  if (!phase || !name || (species_count && !species) ||
      species_count > SIZE_MAX / sizeof(*created.species))
    return -1;
  created.name = duplicate_string(name);
  created.species = (chem_phase_species *)calloc(
      species_count ? species_count : 1, sizeof(*created.species));
  if (!created.name || !created.species) {
    free(created.name);
    free(created.species);
    return -1;
  }
  created.species_count = species_count;
  for (i = 0; i < species_count; ++i) {
    if (chem_phase_species_copy(&created.species[i], &species[i]) != 0) {
      chem_phase_destroy(&created);
      return -1;
    }
  }
  chem_phase_destroy(phase);
  *phase = created;
  return 0;
}

int chem_phase_copy(chem_phase *destination, const chem_phase *source)
{
  if (!destination || !source || !source->name)
    return -1;
  return chem_phase_init(destination, source->name, source->species,
                        source->species_count);
}

size_t chem_phase_state_size(const chem_phase *phase)
{
  size_t i, count = 0;
  if (!phase)
    return 0;
  for (i = 0; i < phase->species_count; ++i)
    count += chem_species_is_parameterized(&phase->species[i].species) == 0;
  return count;
}

int chem_phase_unique_name(const chem_phase *phase, size_t phase_species_index,
                           char *buffer, size_t buffer_size)
{
  const chem_species *species;
  const int written = phase && phase->name &&
                              phase_species_index < phase->species_count && buffer
                          ? (chem_species_is_parameterized(
                                 &phase->species[phase_species_index].species)
                                 ? 0
                                 : snprintf(buffer, buffer_size, "%s.%s", phase->name,
                                            phase->species[phase_species_index].species.name))
                          : -1;
  if (written == 0)
    return 0;
  if (written < 0 || (size_t)written >= buffer_size)
    return -1;
  species = &phase->species[phase_species_index].species;
  return species->name ? 1 : -1;
}

void chem_system_destroy(chem_system *system)
{
  if (!system)
    return;
  chem_phase_destroy(&system->gas_phase);
  *system = (chem_system){0};
}

int chem_system_init(chem_system *system, const chem_phase *gas_phase)
{
  chem_system created = {0};
  if (!system || !gas_phase || chem_phase_copy(&created.gas_phase, gas_phase) != 0)
    return -1;
  chem_system_destroy(system);
  *system = created;
  return 0;
}

int chem_system_copy(chem_system *destination, const chem_system *source)
{
  if (!destination || !source)
    return -1;
  return chem_system_init(destination, &source->gas_phase);
}

size_t chem_system_state_size(const chem_system *system)
{
  return system ? chem_phase_state_size(&system->gas_phase) : 0;
}

int chem_system_unique_name(const chem_system *system, size_t state_index,
                            const size_t *permutation, size_t permutation_size,
                            char *buffer, size_t buffer_size)
{
  const chem_phase *phase;
  size_t wanted, i, state = 0;
  if (!system || !buffer ||
      (permutation && permutation_size != chem_system_state_size(system)) ||
      (!permutation && permutation_size != 0) ||
      state_index >= chem_system_state_size(system))
    return -1;
  phase = &system->gas_phase;
  wanted = permutation ? permutation[state_index] : state_index;
  if (wanted >= chem_system_state_size(system))
    return -1;
  for (i = 0; i < phase->species_count; ++i) {
    const chem_species *species = &phase->species[i].species;
    int written;
    if (chem_species_is_parameterized(species))
      continue;
    if (state++ != wanted)
      continue;
    written = snprintf(buffer, buffer_size, "%s", species->name);
    return written >= 0 && (size_t)written < buffer_size ? 1 : -1;
  }
  return -1;
}
