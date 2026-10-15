#ifndef CHEM_SYSTEM_H
#define CHEM_SYSTEM_H
#include "chem.h"

/* Species metadata owns its name, properties, and string values. All getter
   pointers are borrowed. Metadata is assembled before mechanism construction. */
typedef enum { CHEM_PROPERTY_REAL, CHEM_PROPERTY_INTEGER, CHEM_PROPERTY_BOOLEAN, CHEM_PROPERTY_TEXT } ChemPropertyType;
typedef union { double real; int integer, boolean; const char *text; } ChemPropertyValue;
typedef struct { const char *name; ChemPropertyType type; ChemPropertyValue value; } ChemProperty;
typedef struct ChemSpecies ChemSpecies;
ChemStatus chem_species_create(const char *name, size_t count, const ChemProperty *, ChemSpecies **);
void chem_species_free(ChemSpecies *);
ChemStatus chem_species_clone(const ChemSpecies *, ChemSpecies **);
ChemStatus chem_species_assign(ChemSpecies **, const ChemSpecies *);
const char *chem_species_name(const ChemSpecies *);
size_t chem_species_property_count(const ChemSpecies *, ChemPropertyType);
ChemStatus chem_species_set_property(ChemSpecies *, ChemProperty);
ChemStatus chem_species_get_property(const ChemSpecies *, const char *name, ChemPropertyType, ChemPropertyValue *);
void chem_species_parameterize(ChemSpecies *, ChemParameterization);
ChemParameterization chem_species_parameterization(const ChemSpecies *);
void chem_species_third_body(ChemSpecies *);

typedef struct { ChemSpecies *species; double diffusion; int has_diffusion; } ChemPhaseSpecies;
typedef struct ChemPhase ChemPhase;
/* A phase owns independent copies of its members' metadata. A single gas phase
   is the system in this CPU API. System names omit the phase prefix. */
ChemStatus chem_phase_create(const char *name, size_t count, const ChemPhaseSpecies *, ChemPhase **);
void chem_phase_free(ChemPhase *);
const char *chem_phase_name(const ChemPhase *);
size_t chem_phase_count(const ChemPhase *);
const ChemPhaseSpecies *chem_phase_species(const ChemPhase *, size_t index);
size_t chem_phase_state_size(const ChemPhase *);
/* Variable indices exclude parameterized species. Return length excluding NUL,
   or SIZE_MAX for an invalid index. permutation is optional. */
size_t chem_phase_variable_name(const ChemPhase *, size_t variable, int prefix,
    const size_t *permutation, char *buffer, size_t capacity);
ChemStatus chem_surface_rate(const ChemPhaseSpecies *, double reaction_probability,
    int parameter_index, ChemRate *);
#endif
