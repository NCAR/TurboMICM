#ifndef CHEM_RATE_STORE_H
#define CHEM_RATE_STORE_H
#include "chem.h"
/* Groups are in solver-builder order; entries retain their original reaction
 * index so that callers can map rates back to the mechanism. */
typedef enum {
    CHEM_RATE_ARRHENIUS, CHEM_RATE_TROE, CHEM_RATE_TERNARY, CHEM_RATE_TERNARY_JPL19,
    CHEM_RATE_BRANCHED, CHEM_RATE_TUNNELING, CHEM_RATE_TAYLOR, CHEM_RATE_REVERSIBLE,
    CHEM_RATE_USER, CHEM_RATE_SURFACE, CHEM_RATE_CALLBACK, CHEM_RATE_GROUPS
} ChemRateGroup;
typedef struct {
    ChemReaction reaction;
    size_t reaction_index;
    double branch_k0, branch_z, surface_speed_factor;
    ChemCallback callback;
} ChemRateEntry;
typedef struct { size_t count, offsets[CHEM_RATE_GROUPS+1]; ChemRateEntry *entries; } ChemRateStore;
/* Initialize the destination to zero before building; callbacks are borrowed. */
int chem_rate_store_build(ChemRateStore *, const ChemMechanism *, const ChemCallback *callbacks);
void chem_rate_store_free(ChemRateStore *);
#endif
