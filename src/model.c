#include "chem.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int chem_validate_reaction_species(const ChemMechanism *m,int ignore_unused)
{
    if (!m || !m->species_count || !m->species || (m->reaction_count && !m->reactions) ||
        (m->parameter_count && (!m->parameter_names || !m->parameter_defaults))) return CHEM_INVALID;
    for (size_t i=0;i<m->species_count;++i) if (!m->species[i].name) return CHEM_INVALID;
    unsigned char *used=ignore_unused ? NULL : calloc(m->species_count,1);
    if (!ignore_unused && !used) return CHEM_MEMORY;
    for (size_t i=0;i<m->reaction_count;++i) {
        const ChemReaction *r=&m->reactions[i];
        if ((r->reactant_count && !m->reactants) || (r->product_count && !m->products) ||
            r->reactant_count>SIZE_MAX-r->reactant_start || r->product_count>SIZE_MAX-r->product_start) {
            free(used); return CHEM_INVALID;
        }
        for (size_t j=0;j<r->reactant_count;++j) {
            int index=m->reactants[r->reactant_start+j];
            if (index<0 || (size_t)index>=m->species_count) { free(used); return CHEM_INVALID; }
            if (used) used[index]=1;
        }
        for (size_t j=0;j<r->product_count;++j) {
            int index=m->products[r->product_start+j].species;
            if (index<0 || (size_t)index>=m->species_count) { free(used); return CHEM_INVALID; }
            if (used) used[index]=1;
        }
    }
    if (used) for (size_t i=0;i<m->species_count;++i) if (!used[i]) { free(used); return CHEM_UNUSED_SPECIES; }
    free(used); return CHEM_OK;
}

int chem_state_build(ChemState **state,const ChemMechanism *m,size_t cells,int ignore_unused)
{
    if (!state || !cells) return CHEM_INVALID;
    int status=chem_validate_reaction_species(m,ignore_unused);
    if (status) return status;
    ChemState *created=chem_state_create(m,cells);
    if (!created) return CHEM_INVALID;
    *state=created; return CHEM_OK;
}

void chem_mechanism_composition_destroy(ChemMechanismComposition *c)
{
    if (!c) return;
    free(c->species_storage);
    free(c->parameter_name_storage);
    free(c->parameter_default_storage);
    memset(c,0,sizeof(*c));
}

int chem_mechanism_compose(ChemMechanismComposition *result, const ChemMechanism *base,
                           const ChemExternalModel *models, size_t model_count)
{
    if (!result || !base || (model_count && !models) ||
        (base->species_count && !base->species) ||
        (base->parameter_count && (!base->parameter_names || !base->parameter_defaults))) return CHEM_INVALID;
    size_t ns=base->species_count, np=base->parameter_count;
    for (size_t i=0;i<model_count;++i) {
        const ChemExternalModel *m=&models[i];
        if ((m->species_count && !m->species) ||
            (m->parameter_count && (!m->parameter_names || !m->parameter_defaults)) ||
            m->species_count>SIZE_MAX-ns || m->parameter_count>SIZE_MAX-np) return CHEM_INVALID;
        ns+=m->species_count; np+=m->parameter_count;
    }
    if (ns>SIZE_MAX/sizeof(ChemSpecies) || np>SIZE_MAX/sizeof(const char *) ||
        np>SIZE_MAX/sizeof(double)) return CHEM_INVALID;
    ChemMechanismComposition c={0};
    c.species_storage=calloc(ns ? ns : 1,sizeof(*c.species_storage));
    c.parameter_name_storage=calloc(np ? np : 1,sizeof(*c.parameter_name_storage));
    c.parameter_default_storage=calloc(np ? np : 1,sizeof(*c.parameter_default_storage));
    if (!c.species_storage || !c.parameter_name_storage || !c.parameter_default_storage) {
        chem_mechanism_composition_destroy(&c); return CHEM_MEMORY;
    }
    c.mechanism=*base;
    c.mechanism.species=c.species_storage; c.mechanism.species_count=ns;
    c.mechanism.parameter_names=c.parameter_name_storage;
    c.mechanism.parameter_defaults=c.parameter_default_storage; c.mechanism.parameter_count=np;
    size_t next_s=0, next_p=0;
    for (size_t source=0;source<=model_count;++source) {
        const ChemSpecies *species=source ? models[source-1].species : base->species;
        const char *const *names=source ? models[source-1].parameter_names : base->parameter_names;
        const double *defaults=source ? models[source-1].parameter_defaults : base->parameter_defaults;
        size_t n_species=source ? models[source-1].species_count : base->species_count;
        size_t n_parameters=source ? models[source-1].parameter_count : base->parameter_count;
        for (size_t i=0;i<n_species;++i) {
            if (!species[i].name) { chem_mechanism_composition_destroy(&c); return CHEM_INVALID; }
            for (size_t j=0;j<next_s;++j) if (!strcmp(species[i].name,c.species_storage[j].name)) {
                chem_mechanism_composition_destroy(&c); return CHEM_INVALID;
            }
            c.species_storage[next_s++]=species[i];
        }
        for (size_t i=0;i<n_parameters;++i) {
            if (!names[i]) { chem_mechanism_composition_destroy(&c); return CHEM_INVALID; }
            for (size_t j=0;j<next_p;++j) if (!strcmp(names[i],c.parameter_name_storage[j])) {
                chem_mechanism_composition_destroy(&c); return CHEM_INVALID;
            }
            c.parameter_name_storage[next_p]=names[i]; c.parameter_default_storage[next_p++]=defaults[i];
        }
    }
    *result=c;
    return CHEM_OK;
}
