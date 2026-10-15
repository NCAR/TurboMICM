#ifndef CHEM_INTERNAL_H
#define CHEM_INTERNAL_H
#include "chem.h"
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    size_t reactant_count, product_count, parameter_count;
    size_t *reactants, *parameters;
    ChemProduct *products;
    ChemRate rate;
} ChemProcess;
struct ChemMechanism {
    size_t species_count, variable_count, reaction_count, custom_count;
    ChemSpecies *species;
    size_t *variable_species;
    ChemProcess *processes;
    char **custom_names;
    size_t offsets[CHEM_RATE_TYPE_COUNT + 1];
};
unsigned char *chem_mechanism_pattern(const ChemMechanism *);
static inline void *chem_alloc(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size) return NULL;
    return calloc(n ? n : 1, size);
}
#endif
