#ifndef CHEM_REORDER_H
#define CHEM_REORDER_H
#include "constraint.h"
#include "chem.h"
#include <stddef.h>

/* Symmetric Markowitz ordering of a row-major structural Jacobian pattern.
   old_to_new[old] and new_to_old[new] are inverse mappings. */
typedef struct { size_t size, *old_to_new, *new_to_old; } ChemPermutation;
typedef enum { CHEM_REORDER_OK, CHEM_REORDER_INVALID, CHEM_REORDER_NO_MEMORY } ChemReorderStatus;
ChemReorderStatus chem_permutation_init(ChemPermutation *permutation,size_t size,const unsigned char *pattern);
void chem_permutation_destroy(ChemPermutation *permutation);
/* Derive an ordering from the combined kinetic and constraint structure. */
ChemReorderStatus chem_permutation_for_system(ChemPermutation *permutation,
    const ChemMechanism *mechanism,const ChemConstraintSet *constraints);
/* Independently owned reaction tables, with species indices permuted. */
typedef struct {
    ChemMechanism mechanism;
    ChemReaction *reactions;
    size_t *reactants;
    ChemProduct *products;
} ChemReorderedMechanism;
ChemReorderStatus chem_reordered_mechanism_init(ChemReorderedMechanism *destination,
    const ChemMechanism *source,const ChemPermutation *permutation);
void chem_reordered_mechanism_destroy(ChemReorderedMechanism *mechanism);
/* A reordered constraint set owns its rows and terms; constraint names remain
   borrowed from the source. Custom callbacks require a caller-provided adapter
   and are rejected by this helper. Initialize destination with {0}. */
typedef struct { ChemConstraintSet set; ChemConstraint *rows; ChemConstraintTerm *terms; } ChemReorderedConstraints;
ChemReorderStatus chem_reordered_constraints_init(ChemReorderedConstraints *destination,
    const ChemConstraintSet *source,const ChemPermutation *permutation);
void chem_reordered_constraints_destroy(ChemReorderedConstraints *constraints);
#endif
