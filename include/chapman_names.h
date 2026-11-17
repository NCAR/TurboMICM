#ifndef CHEM_CHAPMAN_NAMES_H
#define CHEM_CHAPMAN_NAMES_H
#include <stddef.h>
/* Immutable metadata of the fixed Chapman reference mechanism. The returned
   arrays and strings have static lifetime and must not be modified. */
const char *const *chem_chapman_reaction_names(size_t *count);
const char *const *chem_chapman_photolysis_names(size_t *count);
const char *const *chem_chapman_species_names(size_t *count);
#endif
