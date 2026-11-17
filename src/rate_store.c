#include "rate_store.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int group(ChemRateType type)
{
    switch (type) {
    case CHEM_RATE_ARRHENIUS: return 0;
    case CHEM_RATE_TROE: return 1;
    case CHEM_RATE_TERNARY: return 2;
    case CHEM_RATE_TERNARY_JPL19: return 3;
    case CHEM_RATE_ALKOXY: case CHEM_RATE_NITRATE: return 4;
    case CHEM_RATE_TUNNELING: return 5;
    case CHEM_RATE_TAYLOR: return 6;
    case CHEM_RATE_REVERSIBLE: return 7;
    case CHEM_RATE_USER: return 8;
    case CHEM_RATE_SURFACE: return 9;
    case CHEM_RATE_CALLBACK: return 10;
    }
    return -1;
}

void chem_rate_store_destroy(ChemRateStore *store)
{
    if (store) { free(store->entries); memset(store, 0, sizeof(*store)); }
}

ChemRateStoreStatus chem_rate_store_init_diagnostic(ChemRateStore *store, size_t count,
                                                    const ChemRateDescriptor *descriptors,
                                                    size_t *invalid_index)
{
    if (invalid_index) *invalid_index=SIZE_MAX;
    if (!store || (count && !descriptors)) return CHEM_RATE_STORE_INVALID;
    memset(store, 0, sizeof(*store));
    if (count > SIZE_MAX / sizeof(*store->entries)) return CHEM_RATE_STORE_NO_MEMORY;
    for (size_t i=0;i<count;++i) {
        int g=group(descriptors[i].rate.type);
        if (g<0) {
            if (invalid_index) *invalid_index=i;
            memset(store,0,sizeof(*store)); return CHEM_RATE_STORE_INVALID;
        }
        if (g==9 && !descriptors[i].has_diffusion_coefficient) {
            if (invalid_index) *invalid_index=i;
            memset(store,0,sizeof(*store)); return CHEM_RATE_STORE_MISSING_DIFFUSION;
        }
        if (g==9 && !descriptors[i].has_molecular_weight) {
            if (invalid_index) *invalid_index=i;
            memset(store,0,sizeof(*store)); return CHEM_RATE_STORE_MISSING_MOLECULAR_WEIGHT;
        }
        ++store->offset[g+1];
    }
    for (size_t g=0;g<CHEM_RATE_STORE_GROUPS;++g) store->offset[g+1]+=store->offset[g];
    if (count) {
        store->entries=calloc(count, sizeof(*store->entries));
        if (!store->entries) { memset(store,0,sizeof(*store)); return CHEM_RATE_STORE_NO_MEMORY; }
    }
    size_t next[CHEM_RATE_STORE_GROUPS];
    memcpy(next,store->offset,sizeof(next));
    for (size_t i=0;i<count;++i) {
        const ChemRateDescriptor *d=descriptors+i;
        ChemStoredRate *e=store->entries+next[group(d->rate.type)]++;
        e->rate=d->rate; e->source_index=i;
    }
    const double avogadro=6.02214076e23, gas_constant=1.380649e-23*6.02214076e23;
    const double pi=3.14159265358979323846;
    for (size_t i=0;i<count;++i) {
        ChemStoredRate *e=store->entries+i;
        double *p=e->rate.p;
        if (e->rate.type==CHEM_RATE_ALKOXY || e->rate.type==CHEM_RATE_NITRATE) {
            e->branch_k0=2e-22*avogadro*1e-6*exp(p[3]);
            double a=e->branch_k0*(2.45e19/avogadro*1e6);
            double b=.43*pow(293.0/298.0,-8.0), ratio=a/b;
            double value=a/(1+ratio)*pow(.41,1/(1+pow(log10(ratio),2)));
            e->branch_z=value*(1-p[2])/p[2];
        } else if (e->rate.type==CHEM_RATE_USER) {
            e->rate.custom_index=store->custom_count++;
        } else if (e->rate.type==CHEM_RATE_SURFACE) {
            e->rate.custom_index=store->custom_count;
            store->custom_count+=2;
            e->mean_free_speed_factor=8*gas_constant/(pi*p[1]);
        }
    }
    store->count=count;
    return CHEM_RATE_STORE_OK;
}

ChemRateStoreStatus chem_rate_store_init(ChemRateStore *store, size_t count,
                                         const ChemRateDescriptor *descriptors)
{
    return chem_rate_store_init_diagnostic(store,count,descriptors,NULL);
}

int chem_rate_store_error_message(ChemRateStoreStatus status, const char *species,
                                  char *buffer, size_t buffer_size)
{
    if (!species) species="";
    switch (status) {
    case CHEM_RATE_STORE_MISSING_DIFFUSION:
        return snprintf(buffer,buffer_size,"Diffusion coefficient for species '%s' is not defined",species);
    case CHEM_RATE_STORE_MISSING_MOLECULAR_WEIGHT:
        return snprintf(buffer,buffer_size,"Species: '%s' Property: 'molecular weight [kg mol-1]'",species);
    case CHEM_RATE_STORE_OK: return snprintf(buffer,buffer_size,"No error");
    case CHEM_RATE_STORE_NO_MEMORY: return snprintf(buffer,buffer_size,"Out of memory");
    case CHEM_RATE_STORE_INVALID: return snprintf(buffer,buffer_size,"Invalid rate descriptor");
    }
    return snprintf(buffer,buffer_size,"Unknown rate-store error");
}
