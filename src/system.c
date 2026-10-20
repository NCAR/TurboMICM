#include "chemistry.h"

#include <math.h>
#include <stdio.h>

int chem_phase_species_init(chem_phase_species *phase_species,
                            const chem_species *species,
                            int has_diffusion_coefficient,
                            double diffusion_coefficient) {
  if (phase_species == NULL || species == NULL || species->name == NULL ||
      (has_diffusion_coefficient != 0 && has_diffusion_coefficient != 1) ||
      (has_diffusion_coefficient &&
       (!isfinite(diffusion_coefficient) || diffusion_coefficient < 0.0))) return 1;
  *phase_species = (chem_phase_species){species,has_diffusion_coefficient,
                                       diffusion_coefficient};
  return 0;
}

int chem_phase_species_set_diffusion_coefficient(
    chem_phase_species *phase_species, double diffusion_coefficient) {
  if (phase_species == NULL || phase_species->species == NULL ||
      !isfinite(diffusion_coefficient) || diffusion_coefficient < 0.0) return 1;
  phase_species->has_diffusion_coefficient = 1;
  phase_species->diffusion_coefficient = diffusion_coefficient;
  return 0;
}

static int valid_phase(const chem_phase *phase) {
  if (phase == NULL || phase->name == NULL ||
      (phase->species_count && phase->species == NULL)) return 0;
  for (size_t i = 0; i < phase->species_count; ++i)
    if (phase->species[i].species == NULL ||
        phase->species[i].species->name == NULL) return 0;
  return 1;
}

size_t chem_phase_state_size(const chem_phase *phase) {
  if (!valid_phase(phase)) return 0;
  size_t count = 0;
  for (size_t i = 0; i < phase->species_count; ++i)
    if (!chem_species_is_parameterized(phase->species[i].species)) ++count;
  return count;
}

int chem_phase_unique_name(const chem_phase *phase, size_t species_index,
                           char *buffer, size_t buffer_size) {
  if (!valid_phase(phase) || buffer == NULL || buffer_size == 0 ||
      species_index >= phase->species_count) return -1;
  const chem_species *species = phase->species[species_index].species;
  if (chem_species_is_parameterized(species)) return 1;
  const int size = snprintf(buffer,buffer_size,"%s.%s",phase->name,species->name);
  return size >= 0 && (size_t)size < buffer_size ? 0 : -1;
}

int chem_system_init(chem_system *system, const chem_phase *gas_phase) {
  if (system == NULL || !valid_phase(gas_phase)) return 1;
  system->gas_phase = gas_phase;
  return 0;
}

size_t chem_system_state_size(const chem_system *system) {
  return system && valid_phase(system->gas_phase)
      ? chem_phase_state_size(system->gas_phase) : 0;
}

int chem_system_unique_name(const chem_system *system, size_t species_index,
                            char *buffer, size_t buffer_size) {
  if (system == NULL || !valid_phase(system->gas_phase) || buffer == NULL ||
      buffer_size == 0) return -1;
  const chem_phase *phase = system->gas_phase;
  size_t unique_index = 0;
  for (size_t i = 0; i < phase->species_count; ++i) {
    const chem_species *species = phase->species[i].species;
    if (chem_species_is_parameterized(species)) continue;
    if (unique_index++ == species_index) {
      const int size = snprintf(buffer,buffer_size,"%s",species->name);
      return size >= 0 && (size_t)size < buffer_size ? 0 : -1;
    }
  }
  return -1;
}
