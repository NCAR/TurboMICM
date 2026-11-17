#ifndef CHEM_RATE_STORE_H
#define CHEM_RATE_STORE_H
#include "rates.h"

/* A surface descriptor distinguishes a missing property from a zero value. */
typedef struct {
    ChemRate rate;
    int has_diffusion_coefficient;
    int has_molecular_weight;
} ChemRateDescriptor;

typedef enum {
    CHEM_RATE_STORE_OK, CHEM_RATE_STORE_INVALID, CHEM_RATE_STORE_NO_MEMORY,
    CHEM_RATE_STORE_MISSING_DIFFUSION, CHEM_RATE_STORE_MISSING_MOLECULAR_WEIGHT
} ChemRateStoreStatus;

enum { CHEM_RATE_STORE_GROUPS = 11 };
typedef struct {
    ChemRate rate;
    size_t source_index;
    double branch_k0, branch_z, mean_free_speed_factor;
} ChemStoredRate;
typedef struct {
    size_t count, custom_count;
    /* Half-open groups in reference type order; offset[11] == count. */
    size_t offset[CHEM_RATE_STORE_GROUPS+1];
    ChemStoredRate *entries;
} ChemRateStore;

/* Stable-sort descriptors into rate-type groups and assign user and surface
   custom slots. The store owns its copied entries and must be destroyed.
   Initialize an empty store before use; on failure it remains empty. */
ChemRateStoreStatus chem_rate_store_init(ChemRateStore *store, size_t count,
                                         const ChemRateDescriptor *descriptors);
void chem_rate_store_destroy(ChemRateStore *store);
#endif
