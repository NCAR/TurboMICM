#ifndef CHEM_MECHANISM_METADATA_H
#define CHEM_MECHANISM_METADATA_H
#include "chem.h"

/* Owned names associated with a mechanism. Photolysis indices refer to the
   supplied reaction names, and may be in a different order. Initialize with
   {0}; initialization copies all strings and replaces the object on success. */
typedef struct {
    size_t reaction_count, photolysis_count, species_count;
    char **reaction_names, **photolysis_names, **species_names;
} ChemMechanismMetadata;
typedef enum { CHEM_METADATA_OK, CHEM_METADATA_INVALID, CHEM_METADATA_NO_MEMORY } ChemMetadataStatus;
ChemMetadataStatus chem_mechanism_metadata_init(ChemMechanismMetadata *metadata,
    const ChemMechanism *mechanism, const char *const *reaction_names,
    const char *const *species_names, const size_t *photolysis_indices, size_t photolysis_count);
void chem_mechanism_metadata_destroy(ChemMechanismMetadata *metadata);
#endif
