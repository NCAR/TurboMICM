#ifndef CHEM_INTERNAL_H
#define CHEM_INTERNAL_H
#include "chem.h"
#include "chem_system.h"
#include "chem_constraints.h"
#include "chem_solver.h"
#include "chem_external.h"
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    size_t reactant_count, product_count, parameter_count;
    size_t *reactants, *parameters;
    ChemProduct *products;
    ChemRate rate;
} ChemProcess;
typedef struct {
    ChemConstraint definition;
    size_t row, parameter;
    size_t *reactants, *products;
} ChemBoundConstraint;
typedef struct {
    ChemModelBinding binding;
    ChemExternalModel callbacks;
    void *context;
} ChemBoundModel;
struct ChemMechanism {
    ChemSystem *system;
    size_t species_count, variable_count, reaction_count, custom_count;
    ChemSpecies *species;
    size_t *variable_species;
    ChemProcess *processes;
    char **custom_names;
    size_t offsets[CHEM_RATE_TYPE_COUNT + 1];
    size_t constraint_count, base_custom_count;
    ChemBoundConstraint *constraints;
    unsigned char *algebraic;
    size_t model_count, external_constraint_count;
    ChemBoundModel *models;
};
void chem_models_free(ChemMechanism *);
void chem_models_reorder(ChemMechanism *,const size_t *);
void chem_models_update(ChemState *,int constraints);
void chem_models_residual(const ChemState *,const double *,double *);
void chem_models_jacobian(const ChemState *,const double *,double *,const ChemSparse *,double,int constraints);
void chem_clear_algebraic_jacobian(const ChemState *,double *,const ChemSparse *);
void chem_constraints_free(ChemBoundConstraint *, size_t);
void chem_constraints_jacobian_add(const ChemState *, const double *, double *, const ChemSparse *, double);
unsigned char *chem_mechanism_pattern(const ChemMechanism *);
static inline void *chem_alloc(size_t n, size_t size)
{
    if (size && n > SIZE_MAX / size) return NULL;
    return calloc(n ? n : 1, size);
}
#endif
