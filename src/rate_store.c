#include "rate_store.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int group(ChemRateType type)
{
    switch(type) {
    case CHEM_ARRHENIUS:return CHEM_RATE_ARRHENIUS;
    case CHEM_TROE:return CHEM_RATE_TROE;
    case CHEM_TERNARY:return CHEM_RATE_TERNARY;
    case CHEM_TERNARY_JPL19:return CHEM_RATE_TERNARY_JPL19;
    case CHEM_BRANCH_ALKOXY:case CHEM_BRANCH_NITRATE:return CHEM_RATE_BRANCHED;
    case CHEM_TUNNELING:return CHEM_RATE_TUNNELING;
    case CHEM_TAYLOR:return CHEM_RATE_TAYLOR;
    case CHEM_REVERSIBLE:return CHEM_RATE_REVERSIBLE;
    case CHEM_USER:return CHEM_RATE_USER;
    case CHEM_SURFACE:return CHEM_RATE_SURFACE;
    case CHEM_CALLBACK:return CHEM_RATE_CALLBACK;
    }
    return -1;
}
int chem_rate_store_build(ChemRateStore *store,const ChemMechanism *m,const ChemCallback *callbacks)
{
    if(!store || !m || chem_validate_mechanism(m,NULL,0)!=CHEM_VALID || m->reaction_count>SIZE_MAX/sizeof(ChemRateEntry)) return 1;
    ChemRateStore result={0};result.count=m->reaction_count;
    result.entries=result.count?calloc(result.count,sizeof(*result.entries)):NULL;
    if(result.count && !result.entries)return 1;
    for(size_t i=0;i<result.count;++i) {int g=group(m->reactions[i].type);if(g<0){chem_rate_store_free(&result);return 1;}++result.offsets[g+1];}
    for(size_t g=0;g<CHEM_RATE_GROUPS;++g)result.offsets[g+1]+=result.offsets[g];
    size_t cursor[CHEM_RATE_GROUPS];memcpy(cursor,result.offsets,sizeof cursor);
    for(size_t i=0;i<result.count;++i) {
        const ChemReaction *r=&m->reactions[i];int g=group(r->type);ChemRateEntry *entry=&result.entries[cursor[g]++];
        entry->reaction=*r;entry->reaction_index=i;
        if(g==CHEM_RATE_BRANCHED) {
            double na=6.02214076e23,air=2.45e19/na*1e6;
            entry->branch_k0=2e-22*na*1e-6*exp(r->p[3]);
            double a=entry->branch_k0*air,b=.43*pow(293.0/298.0,-8),q=a/b,logq=log10(q);
            entry->branch_z=a/(1+q)*pow(.41,1/(1+logq*logq))*(1-r->p[2])/r->p[2];
        } else if(g==CHEM_RATE_SURFACE) entry->surface_speed_factor=8*(1.380649e-23*6.02214076e23)/(3.14159265358979323846*r->p[1]);
        else if(g==CHEM_RATE_CALLBACK && callbacks && r->parameter>=0 && (size_t)r->parameter<m->parameter_count) entry->callback=callbacks[r->parameter];
    }
    *store=result;return 0;
}
void chem_rate_store_free(ChemRateStore *store)
{if(!store)return;free(store->entries);*store=(ChemRateStore){0};}
