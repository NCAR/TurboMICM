#ifndef CHEM_EXTERNAL_INTERNAL_H
#define CHEM_EXTERNAL_INTERNAL_H
#include "chem_external.h"
chem_status chem_external_prepare(chem_model *,const chem_mechanism **);
void chem_external_destroy(chem_model *);
void chem_external_pattern(const chem_model *,unsigned char *);
chem_status chem_external_finalize(chem_model *);
void chem_external_update(chem_state *);
int chem_external_species_used(const chem_model *,size_t);
#endif
