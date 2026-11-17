#ifndef CHEM_SYSTEM_H
#define CHEM_SYSTEM_H

#include "chem.h"

/* Typed metadata used while assembling mechanisms. Strings and property arrays
 * belong to the species object. Initialize objects to zero before first use. */
typedef enum { CHEM_PROPERTY_STRING, CHEM_PROPERTY_REAL, CHEM_PROPERTY_INT, CHEM_PROPERTY_BOOL } ChemPropertyType;
typedef union { const char *string; double real; int integer; int boolean; } ChemPropertyValue;
typedef struct { char *name; ChemPropertyType type; ChemPropertyValue value; } ChemProperty;
typedef struct {
    char *name;
    ChemProperty *properties;
    size_t property_count, property_capacity;
    int parameterized;
    ChemParameterizedFunction parameterization;
} ChemSpeciesInfo;

int chem_species_info_init(ChemSpeciesInfo *species, const char *name);
void chem_species_info_destroy(ChemSpeciesInfo *species);
int chem_species_info_copy(ChemSpeciesInfo *destination, const ChemSpeciesInfo *source);
int chem_species_info_set(ChemSpeciesInfo *species, const char *name, ChemPropertyType type, ChemPropertyValue value);
int chem_species_info_get(const ChemSpeciesInfo *species, const char *name, ChemPropertyType type, ChemPropertyValue *value);
int chem_species_info_has(const ChemSpeciesInfo *species, const char *name);
void chem_species_info_set_parameterization(ChemSpeciesInfo *species, ChemParameterizedFunction function);
void chem_species_info_set_third_body(ChemSpeciesInfo *species);
double chem_parameterized_value(ChemParameterizedFunction function, double temperature, double pressure, double air_density);

/* Phase and system views borrow their species, arrays, and phase name. State
 * indices exclude parameterized species. A phase name includes the phase prefix;
 * a system name uses the unqualified gas-phase species name. */
typedef struct { const ChemSpeciesInfo *species; int has_diffusion; double diffusion; } ChemPhaseSpecies;
typedef struct { const char *name; const ChemPhaseSpecies *species; size_t species_count; } ChemPhase;
typedef struct { const ChemPhase *gas; } ChemSystem;
/* Construction-time association; the runtime reaction stores species indices. */
typedef struct { ChemReaction reaction; const ChemPhase *phase; } ChemProcessDefinition;
size_t chem_phase_state_size(const ChemPhase *phase);
int chem_phase_state_name(const ChemPhase *phase, size_t state_index, char *buffer, size_t size);
/* Materialize the non-parameterized phase species in state order. Names are
 * borrowed from the phase metadata. Optional standard properties are copied. */
int chem_phase_export_species(const ChemPhase *phase, ChemSpecies *species, size_t count);
size_t chem_system_state_size(const ChemSystem *system);
int chem_system_state_name(const ChemSystem *system, size_t state_index, char *buffer, size_t size);

#endif
