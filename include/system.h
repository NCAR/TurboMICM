#ifndef CHEM_SYSTEM_H
#define CHEM_SYSTEM_H

#include "kinetics.h"

typedef struct {
  chem_species species;
  double diffusion_coefficient;
  int has_diffusion_coefficient;
} chem_phase_species;

typedef struct {
  char *name;
  chem_phase_species *species;
  size_t species_count;
} chem_phase;

typedef struct {
  chem_phase gas_phase;
} chem_system;

int chem_phase_species_init(chem_phase_species *phase_species,
                            const chem_species *species,
                            int has_diffusion_coefficient,
                            double diffusion_coefficient);
int chem_phase_species_copy(chem_phase_species *destination,
                            const chem_phase_species *source);
void chem_phase_species_destroy(chem_phase_species *phase_species);
int chem_phase_species_set_diffusion_coefficient(
    chem_phase_species *phase_species, double diffusion_coefficient);

int chem_phase_init(chem_phase *phase, const char *name,
                    const chem_phase_species *species, size_t species_count);
int chem_phase_copy(chem_phase *destination, const chem_phase *source);
void chem_phase_destroy(chem_phase *phase);
size_t chem_phase_state_size(const chem_phase *phase);
int chem_phase_unique_name(const chem_phase *phase, size_t phase_species_index,
                           char *buffer, size_t buffer_size);

int chem_system_init(chem_system *system, const chem_phase *gas_phase);
int chem_system_copy(chem_system *destination, const chem_system *source);
void chem_system_destroy(chem_system *system);
size_t chem_system_state_size(const chem_system *system);
int chem_system_unique_name(const chem_system *system, size_t state_index,
                            const size_t *permutation, size_t permutation_size,
                            char *buffer, size_t buffer_size);

#endif
