#ifndef CHEM_SYSTEM_H
#define CHEM_SYSTEM_H
#include "chem.h"

typedef enum { CHEM_PROPERTY_DOUBLE, CHEM_PROPERTY_INT, CHEM_PROPERTY_BOOL, CHEM_PROPERTY_STRING } ChemPropertyType;
typedef struct ChemProperty {
  char *key;
  ChemPropertyType type;
  double number;
  int integer;
  char *text;
} ChemProperty;
typedef struct { char *name; size_t count; ChemSpecies *species; } ChemPhase;

/* Init writes to uninitialized storage; copy writes to initialized storage. */
ChemStatus chem_species_init(ChemSpecies *species, const char *name);
void chem_species_destroy(ChemSpecies *species);
ChemStatus chem_species_copy(ChemSpecies *destination, const ChemSpecies *source);
ChemStatus chem_species_set_property(ChemSpecies *, const char *, ChemPropertyType, const void *value);
ChemStatus chem_species_get_property(const ChemSpecies *, const char *, ChemPropertyType, void *value);
double chem_species_parameter(const ChemSpecies *, const ChemConditions *);
ChemStatus chem_phase_init(ChemPhase *, const char *name, const ChemSpecies *species, size_t count);
void chem_phase_destroy(ChemPhase *);
size_t chem_phase_state_size(const ChemPhase *);
ChemStatus chem_phase_names(const ChemPhase *, int qualified, const size_t *permutation, char ***names, size_t *count);
void chem_names_free(char **names, size_t count);
const char *chem_version(void);
unsigned chem_version_component(unsigned component);
#endif
