#ifndef CHEM_SPECIES_H
#define CHEM_SPECIES_H
#include "rates.h"
#include <stddef.h>

typedef enum { CHEM_PROPERTY_DOUBLE, CHEM_PROPERTY_STRING, CHEM_PROPERTY_INT, CHEM_PROPERTY_BOOL } ChemPropertyType;
typedef union { double real; const char *string; int integer; int boolean; } ChemPropertyValue;
typedef struct { char *name; ChemPropertyType type; ChemPropertyValue value; } ChemProperty;
typedef enum { CHEM_SPECIES_OK, CHEM_SPECIES_INVALID, CHEM_SPECIES_NO_MEMORY,
               CHEM_SPECIES_NOT_FOUND, CHEM_SPECIES_WRONG_TYPE, CHEM_SPECIES_RANGE } ChemSpeciesStatus;
typedef struct {
    char *name;
    size_t property_count;
    ChemProperty *properties;
    int is_parameterized;
    ChemParameterizedSpecies parameterization;
} ChemSpecies;
typedef struct { ChemSpecies species; int has_diffusion; double diffusion; } ChemPhaseSpecies;
typedef struct { char *name; size_t count; ChemPhaseSpecies *species; } ChemPhase;
typedef struct { ChemPhase gas; } ChemSystem;

/* Initialize owned objects to zero before init or clone. Replacement is
   transactional. Property strings and names are copied and owned. */
ChemSpeciesStatus chem_species_init(ChemSpecies *species, const char *name);
ChemSpeciesStatus chem_species_clone(ChemSpecies *destination, const ChemSpecies *source);
void chem_species_destroy(ChemSpecies *species);
ChemSpeciesStatus chem_species_set_property(ChemSpecies *species, const char *name,
                                            ChemPropertyType type, ChemPropertyValue value);
ChemSpeciesStatus chem_species_get_property(const ChemSpecies *species, const char *name,
                                            ChemPropertyType type, ChemPropertyValue *value);
void chem_species_set_parameterization(ChemSpecies *species, ChemParameterizedSpecies value);
void chem_species_set_third_body(ChemSpecies *species);
double chem_species_evaluate(const ChemSpecies *species, ChemConditions conditions);

ChemSpeciesStatus chem_phase_species_init(ChemPhaseSpecies *entry, const ChemSpecies *species);
void chem_phase_species_destroy(ChemPhaseSpecies *entry);
void chem_phase_species_set_diffusion(ChemPhaseSpecies *entry, double coefficient);
ChemSpeciesStatus chem_phase_init(ChemPhase *phase, const char *name,
                                  size_t count, const ChemPhaseSpecies *species);
ChemSpeciesStatus chem_phase_clone(ChemPhase *destination, const ChemPhase *source);
void chem_phase_destroy(ChemPhase *phase);
size_t chem_phase_state_size(const ChemPhase *phase);
/* State index excludes parameterized species. With qualified set, writes
   phase.species, otherwise species. Output includes a terminating null. */
ChemSpeciesStatus chem_phase_state_name(const ChemPhase *phase, size_t state_index, int qualified,
                                        char *buffer, size_t capacity);
ChemSpeciesStatus chem_system_init(ChemSystem *system, const ChemPhase *gas);
void chem_system_destroy(ChemSystem *system);
#endif
