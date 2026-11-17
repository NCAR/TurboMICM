#include "process.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void chem_process_destroy(ChemProcess *p)
{
    if (!p) return;
    for (size_t i=0;i<p->reactant_count;++i) chem_species_destroy(p->reactants+i);
    for (size_t i=0;i<p->product_count;++i) chem_species_destroy(&p->products[i].species);
    free(p->reactants); free(p->products); chem_phase_destroy(&p->phase); *p=(ChemProcess){0};
}
ChemSpeciesStatus chem_process_init(ChemProcess *p,const ChemProcessConfig *c)
{
    if (!p || !c || !c->phase || (c->reactant_count && !c->reactants) ||
        (c->product_count && (!c->products || !c->yields))) return CHEM_SPECIES_INVALID;
    if (c->reactant_count>SIZE_MAX/sizeof(*p->reactants) || c->product_count>SIZE_MAX/sizeof(*p->products))
        return CHEM_SPECIES_NO_MEMORY;
    ChemProcess copy={0};
    copy.reactants=c->reactant_count ? calloc(c->reactant_count,sizeof(*copy.reactants)) : NULL;
    copy.products=c->product_count ? calloc(c->product_count,sizeof(*copy.products)) : NULL;
    if ((c->reactant_count && !copy.reactants) || (c->product_count && !copy.products)) {
        chem_process_destroy(&copy); return CHEM_SPECIES_NO_MEMORY;
    }
    ChemSpeciesStatus status=chem_phase_clone(&copy.phase,c->phase);
    if (status!=CHEM_SPECIES_OK) { chem_process_destroy(&copy); return status; }
    copy.rate=c->rate;
    for (size_t i=0;i<c->reactant_count;++i) {
        status=chem_species_clone(copy.reactants+i,c->reactants[i]);
        if (status!=CHEM_SPECIES_OK) { chem_process_destroy(&copy); return status; }
        ++copy.reactant_count;
    }
    for (size_t i=0;i<c->product_count;++i) {
        status=chem_species_clone(&copy.products[i].species,c->products[i]);
        if (status!=CHEM_SPECIES_OK) { chem_process_destroy(&copy); return status; }
        copy.products[i].yield=c->yields[i]; ++copy.product_count;
    }
    chem_process_destroy(p); *p=copy; return CHEM_SPECIES_OK;
}
ChemSpeciesStatus chem_process_clone(ChemProcess *destination,const ChemProcess *source)
{
    if (!destination || !source || !source->phase.name) return CHEM_SPECIES_INVALID;
    if (destination==source) return CHEM_SPECIES_OK;
    if (source->reactant_count>SIZE_MAX/sizeof(ChemSpecies *) || source->product_count>SIZE_MAX/sizeof(ChemSpecies *) ||
        source->product_count>SIZE_MAX/sizeof(double)) return CHEM_SPECIES_NO_MEMORY;
    const ChemSpecies **reactants=source->reactant_count ? malloc(source->reactant_count*sizeof(*reactants)) : NULL;
    const ChemSpecies **products=source->product_count ? malloc(source->product_count*sizeof(*products)) : NULL;
    double *yields=source->product_count ? malloc(source->product_count*sizeof(*yields)) : NULL;
    if ((source->reactant_count && !reactants) || (source->product_count && (!products || !yields))) {
        free(reactants); free(products); free(yields); return CHEM_SPECIES_NO_MEMORY;
    }
    for (size_t i=0;i<source->reactant_count;++i) reactants[i]=source->reactants+i;
    for (size_t i=0;i<source->product_count;++i) { products[i]=&source->products[i].species; yields[i]=source->products[i].yield; }
    ChemProcessConfig config={source->reactant_count,source->product_count,reactants,products,yields,&source->phase,source->rate};
    ChemSpeciesStatus status=chem_process_init(destination,&config);
    free(reactants); free(products); free(yields); return status;
}
void chem_built_mechanism_destroy(ChemBuiltMechanism *b)
{
    if (!b) return;
    chem_system_destroy(&b->system); free(b->reactions); free(b->reactants); free(b->products);
    free(b->rates); free(b->multipliers); free(b->functions); free(b->state_names);
    free(b->absolute_tolerances); *b=(ChemBuiltMechanism){0};
}
static size_t find_species(const ChemPhase *phase,const char *name)
{
    for (size_t i=0;i<phase->count;++i) if (!strcmp(phase->species[i].species.name,name)) return i;
    return phase->count;
}
ChemBuildStatus chem_build_mechanism(ChemBuiltMechanism *out,const ChemSystem *system,
                                     size_t count,const ChemProcess *processes,int ignore_unused)
{
    if (!out || (count && !processes)) return CHEM_BUILD_INVALID;
    if (!system || !system->gas.name) return CHEM_BUILD_MISSING_SYSTEM;
    size_t nr=0,np=0;
    for (size_t i=0;i<count;++i) {
        if (processes[i].reactant_count>SIZE_MAX-nr || processes[i].product_count>SIZE_MAX-np) return CHEM_BUILD_NO_MEMORY;
        nr+=processes[i].reactant_count; np+=processes[i].product_count;
    }
    size_t n=chem_phase_state_size(&system->gas), all=system->gas.count;
    if (count>SIZE_MAX/sizeof(ChemReaction) || count>SIZE_MAX/sizeof(ChemRate) ||
        count>SIZE_MAX/sizeof(ChemRateMultiplier) || nr>SIZE_MAX/sizeof(size_t) ||
        nr>SIZE_MAX/sizeof(ChemParameterizedSpecies) || np>SIZE_MAX/sizeof(ChemProduct) ||
        n>SIZE_MAX/sizeof(char *) || n>SIZE_MAX/sizeof(double) || all>SIZE_MAX/sizeof(size_t)) return CHEM_BUILD_NO_MEMORY;
    ChemBuiltMechanism b={0};
    if (chem_system_init(&b.system,&system->gas)!=CHEM_SPECIES_OK) return CHEM_BUILD_NO_MEMORY;
    b.reactions=count ? calloc(count,sizeof(*b.reactions)) : NULL;
    b.rates=count ? malloc(count*sizeof(*b.rates)) : NULL;
    b.multipliers=count ? calloc(count,sizeof(*b.multipliers)) : NULL;
    b.reactants=nr ? malloc(nr*sizeof(*b.reactants)) : NULL;
    b.functions=nr ? malloc(nr*sizeof(*b.functions)) : NULL;
    b.products=np ? malloc(np*sizeof(*b.products)) : NULL;
    b.state_names=n ? malloc(n*sizeof(*b.state_names)) : NULL;
    b.absolute_tolerances=n ? malloc(n*sizeof(*b.absolute_tolerances)) : NULL;
    size_t *map=all ? malloc(all*sizeof(*map)) : NULL;
    unsigned char *used=all ? calloc(all,1) : NULL;
    ChemBuildStatus result=CHEM_BUILD_OK;
    if ((count && (!b.reactions || !b.rates || !b.multipliers)) || (nr && (!b.reactants || !b.functions)) ||
        (np && !b.products) || (n && (!b.state_names || !b.absolute_tolerances)) || (all && (!map || !used))) {
        result=CHEM_BUILD_NO_MEMORY; goto done;
    }
    size_t state_index=0;
    for (size_t i=0;i<all;++i) {
        ChemSpecies *s=&b.system.gas.species[i].species;
        if (s->is_parameterized) { map[i]=n; continue; }
        map[i]=state_index; b.state_names[state_index]=s->name;
        ChemPropertyValue value;
        ChemSpeciesStatus status=chem_species_get_property(s,"absolute tolerance",CHEM_PROPERTY_DOUBLE,&value);
        if (status==CHEM_SPECIES_WRONG_TYPE || (status==CHEM_SPECIES_OK && (!isfinite(value.real) || value.real<=0))) {
            result=CHEM_BUILD_BAD_TOLERANCE; goto done;
        }
        b.absolute_tolerances[state_index++]=status==CHEM_SPECIES_OK ? value.real : .001;
    }
    size_t reactant_pos=0,product_pos=0,function_pos=0;
    for (size_t i=0;i<count;++i) {
        const ChemProcess *p=processes+i;
        ChemReaction *r=b.reactions+i; r->reactant_begin=reactant_pos; r->product_begin=product_pos;
        b.multipliers[i].begin=function_pos; b.rates[i]=p->rate;
        for (size_t j=0;j<p->reactant_count;++j) {
            size_t source=find_species(&b.system.gas,p->reactants[j].name);
            if (source==all) { result=CHEM_BUILD_UNKNOWN_SPECIES; goto done; }
            used[source]=1;
            const ChemSpecies *species=&b.system.gas.species[source].species;
            if (species->is_parameterized) { b.functions[function_pos++]=species->parameterization; ++b.multipliers[i].count; }
            else { b.reactants[reactant_pos++]=map[source]; ++r->reactant_count; }
        }
        for (size_t j=0;j<p->product_count;++j) {
            size_t source=find_species(&b.system.gas,p->products[j].species.name);
            if (source==all) { result=CHEM_BUILD_UNKNOWN_SPECIES; goto done; }
            used[source]=1;
            if (map[source]!=n) { b.products[product_pos++]=(ChemProduct){map[source],p->products[j].yield}; ++r->product_count; }
        }
    }
    if (!ignore_unused) for (size_t i=0;i<all;++i) if (!used[i]) { result=CHEM_BUILD_UNUSED_SPECIES; goto done; }
    b.mechanism=(ChemMechanism){n,count,b.reactions,b.reactants,b.products};
    chem_built_mechanism_destroy(out); *out=b; b=(ChemBuiltMechanism){0};
done:
    free(map); free(used); chem_built_mechanism_destroy(&b); return result;
}
