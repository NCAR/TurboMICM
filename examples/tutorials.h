#ifndef CHEM_TUTORIALS_H
#define CHEM_TUTORIALS_H
#include "chem.h"
#include <stdio.h>
typedef enum { EXAMPLE_README, EXAMPLE_ROBERTSON, EXAMPLE_RATES, EXAMPLE_USER_RATES } example_kind;
typedef struct {
    chem_species species[7];
    chem_reaction reactions[10];
    size_t reactants[10][2];
    chem_product products[10][2];
    chem_mechanism mechanism;
    chem_model model;
} example_chemistry;
chem_status example_create(example_chemistry *,example_kind);
chem_status example_initialize(chem_state *,example_kind,int from_configuration);
int example_run(const char *name,FILE *output);
#endif
