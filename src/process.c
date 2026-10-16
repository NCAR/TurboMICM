#include "chem_process.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *string_copy(const char *s)
{ size_t n=strlen(s)+1; char *copy=malloc(n); if (copy) memcpy(copy,s,n); return copy; }
static size_t phase_index(const chem_phase *p,const char *name)
{
    if (!name) return CHEM_NONE;
    for (size_t i=0;i<p->n_species;++i) if (!strcmp(p->species[i].info.species.name,name)) return i;
    return CHEM_NONE;
}
static void reaction_destroy(chem_reaction *r)
{ free((void *)r->reactants); free((void *)r->products); free((void *)r->rate.label); *r=(chem_reaction){0}; }
void chem_process_destroy(chem_process *p)
{ if (p) { reaction_destroy(&p->reaction); chem_phase_destroy(&p->phase); } }
chem_status chem_process_create(chem_process *p,const chem_phase *phase,size_t nr,const char *const *reactants,
    size_t np,const chem_named_product *products,chem_rate rate)
{
    if (!p || !phase || (nr && !reactants) || (np && !products)) return CHEM_INVALID;
    *p=(chem_process){0}; chem_status status=chem_phase_create(&p->phase,phase->name,phase->n_species,phase->species);
    if (status!=CHEM_OK) return status;
    size_t *rr=calloc(nr ? nr : 1,sizeof(*rr)); chem_product *pp=calloc(np ? np : 1,sizeof(*pp));
    p->reaction=(chem_reaction){nr,np,rr,pp,rate}; p->reaction.rate.label=NULL;
    if (rate.label) p->reaction.rate.label=string_copy(rate.label);
    if (!rr || !pp || (rate.label && !p->reaction.rate.label)) { status=CHEM_MEMORY; goto fail; }
    for (size_t i=0;i<nr;++i) if ((rr[i]=phase_index(phase,reactants[i]))==CHEM_NONE) { status=CHEM_NOT_FOUND; goto fail; }
    for (size_t i=0;i<np;++i) {
        pp[i]=(chem_product){phase_index(phase,products[i].name),products[i].yield};
        if (pp[i].species==CHEM_NONE) { status=CHEM_NOT_FOUND; goto fail; }
    }
    return CHEM_OK;
fail:
    chem_process_destroy(p); return status;
}
chem_status chem_process_clone(chem_process *dest,const chem_process *src)
{
    if (!dest || !src || dest==src) return CHEM_INVALID;
    size_t nr=src->reaction.n_reactants,np=src->reaction.n_products;
    const char **rr=calloc(nr ? nr : 1,sizeof(*rr)); chem_named_product *pp=calloc(np ? np : 1,sizeof(*pp));
    if (!rr || !pp) { free(rr); free(pp); return CHEM_MEMORY; }
    for (size_t i=0;i<nr;++i) rr[i]=src->phase.species[src->reaction.reactants[i]].info.species.name;
    for (size_t i=0;i<np;++i) pp[i]=(chem_named_product){src->phase.species[src->reaction.products[i].species].info.species.name,src->reaction.products[i].yield};
    chem_status status=chem_process_create(dest,&src->phase,nr,rr,np,pp,src->reaction.rate);
    free(rr); free(pp); return status;
}
chem_status chem_process_replace(chem_process *dest,const chem_process *src)
{
    if (!dest || !src) return CHEM_INVALID;
    if (dest==src) return CHEM_OK;
    chem_process copy={0}; chem_status status=chem_process_clone(&copy,src);
    if (status==CHEM_OK) { chem_process_destroy(dest); *dest=copy; }
    return status;
}
chem_status chem_surface_rate(chem_rate *rate,const char *label,const chem_phase_species *species,double probability)
{
    if (!rate || !label || !species) return CHEM_INVALID;
    *rate=chem_rate_default(CHEM_SURFACE); rate->label=label;
    rate->p.surface=(chem_surface){species->has_diffusion ? species->diffusion : 0,species->info.species.molecular_weight,probability,0};
    return chem_rate_prepare(rate);
}
int chem_rate_type_order(chem_rate_type type)
{
    if (type==CHEM_BRANCHED) return 4;
    if (type==CHEM_TUNNELING) return 5;
    return (int)type;
}
void chem_reactions_sort(chem_reaction *reactions,size_t count)
{
    for (size_t i=1;i<count;++i) {
        chem_reaction r=reactions[i]; size_t j=i;
        while (j && chem_rate_type_order(reactions[j-1].rate.type)>chem_rate_type_order(r.rate.type)) {
            reactions[j]=reactions[j-1]; --j;
        }
        reactions[j]=r;
    }
}
size_t chem_model_rate_count(const chem_model *m,chem_rate_type type)
{ size_t n=0; for (size_t i=0;i<m->n_reactions;++i) n+=m->rates[i].type==type; return n; }
size_t chem_model_rate_offset(const chem_model *m,chem_rate_type type)
{ size_t n=0; for (size_t i=0;i<m->n_reactions;++i) n+=chem_rate_type_order(m->rates[i].type)<chem_rate_type_order(type); return n; }
static int has_parameter(const chem_model *m,size_t r)
{
    const chem_reaction *rx=m->mechanism->reactions+r;
    for (size_t j=0;j<rx->n_reactants;++j) if (m->mechanism->species[rx->reactants[j]].type!=CHEM_SPECIES_VARIABLE) return 1;
    return 0;
}
size_t chem_model_parameterized_count(const chem_model *m)
{ size_t n=0; for (size_t i=0;i<m->n_reactions;++i) n+=(size_t)has_parameter(m,i); return n; }
size_t chem_model_parameterized_index(const chem_model *m,size_t index)
{ for (size_t i=0;i<m->n_reactions;++i) if (has_parameter(m,i) && index--==0) return i; return CHEM_NONE; }
double chem_model_parameterized_multiplier(const chem_model *m,size_t r,const chem_conditions *c)
{
    double result=1; const chem_reaction *rx=m->mechanism->reactions+r;
    for (size_t j=0;j<rx->n_reactants;++j) {
        const chem_species *sp=m->mechanism->species+rx->reactants[j];
        if (sp->type!=CHEM_SPECIES_VARIABLE) result*=chem_species_evaluate(sp,c);
    }
    return result;
}
int chem_model_species_used(const chem_model *m,size_t species)
{
    if (!m || species>=m->mechanism->n_species) return 0;
    for (size_t r=0;r<m->n_reactions;++r) {
        const chem_reaction *rx=m->mechanism->reactions+r;
        for (size_t j=0;j<rx->n_reactants;++j) if (rx->reactants[j]==species) return 1;
        for (size_t j=0;j<rx->n_products;++j) if (rx->products[j].species==species) return 1;
    }
    return 0;
}
chem_status chem_model_kinetic_pattern(const chem_model *m,unsigned char *pattern,size_t size)
{
    if (!m || !pattern) return CHEM_INVALID;
    size_t n=m->n_species; if (size!=n*n) return CHEM_SIZE;
    memset(pattern,0,size);
    for (size_t r=0;r<m->n_reactions;++r) {
        const chem_reaction *rx=m->mechanism->reactions+r;
        for (size_t j=0;j<rx->n_reactants;++j) {
            size_t col=m->spec_to_var[rx->reactants[j]]; if (col==CHEM_NONE) continue;
            for (size_t k=0;k<rx->n_reactants;++k) {
                size_t row=m->spec_to_var[rx->reactants[k]]; if (row!=CHEM_NONE) pattern[row*n+col]=1;
            }
            for (size_t k=0;k<rx->n_products;++k) {
                size_t row=m->spec_to_var[rx->products[k].species]; if (row!=CHEM_NONE) pattern[row*n+col]=1;
            }
        }
    }
    return CHEM_OK;
}
void chem_process_model_destroy(chem_process_model *p)
{
    if (!p) return;
    chem_model_destroy(&p->model);
    for (size_t i=0;i<p->mechanism.n_reactions;++i) reaction_destroy(p->reactions+i);
    free(p->reactions); free(p->species); chem_system_destroy(&p->system); *p=(chem_process_model){0};
}
chem_status chem_process_model_create(chem_process_model *p,const chem_system *system,size_t count,const chem_process *processes,
    int preserve_order,int reject_unused,chem_lu_algorithm algorithm,chem_error *error)
{
    if (error) *error=(chem_error){0};
    if (!p || !system || !system->gas.n_species || (count && !processes)) {
        if (error) { error->code=CHEM_INVALID; snprintf(error->message,sizeof(error->message),"Missing chemical system or processes"); }
        return CHEM_INVALID;
    }
    *p=(chem_process_model){0}; chem_status status=chem_system_create(&p->system,&system->gas);
    if (status!=CHEM_OK) goto fail;
    size_t n=system->gas.n_species;
    p->species=calloc(n,sizeof(*p->species)); p->reactions=calloc(count ? count : 1,sizeof(*p->reactions));
    unsigned char *used=calloc(n,1);
    if (!p->species || !p->reactions || !used) { free(used); status=CHEM_MEMORY; goto fail; }
    p->mechanism=(chem_mechanism){.n_species=n,.species=p->species,.reactions=p->reactions,.preserve_order=preserve_order};
    for (size_t i=0;i<n;++i) p->species[i]=p->system.gas.species[i].info.species;
    for (size_t i=0;i<count;++i) {
        const chem_process *src=processes+i; chem_process copy={0};
        status=chem_process_clone(&copy,src);
        if (status!=CHEM_OK) { free(used); goto fail; }
        chem_reaction *r=p->reactions+i; *r=copy.reaction; copy.reaction=(chem_reaction){0};
        chem_process_destroy(&copy); ++p->mechanism.n_reactions;
        for (size_t j=0;j<r->n_reactants;++j) {
            size_t index=phase_index(&p->system.gas,src->phase.species[r->reactants[j]].info.species.name);
            if (index==CHEM_NONE) { free(used); status=CHEM_NOT_FOUND; goto fail; }
            ((size_t *)r->reactants)[j]=index; used[index]=1;
        }
        for (size_t j=0;j<r->n_products;++j) {
            size_t index=phase_index(&p->system.gas,src->phase.species[r->products[j].species].info.species.name);
            if (index==CHEM_NONE) { free(used); status=CHEM_NOT_FOUND; goto fail; }
            ((chem_product *)r->products)[j].species=index; used[index]=1;
        }
    }
    if (reject_unused) for (size_t i=0;i<n;++i) if (!used[i]) {
        if (error) snprintf(error->message,sizeof(error->message),"Unused species: %s",p->species[i].name);
        free(used); status=CHEM_INVALID; goto fail;
    }
    free(used); chem_reactions_sort(p->reactions,count);
    status=chem_model_create_checked(&p->model,&p->mechanism,algorithm,error);
    if (status==CHEM_OK) return status;
fail:
    if (error) { error->code=status; if (!error->message[0]) snprintf(error->message,sizeof(error->message),"%s",chem_status_string(status)); }
    chem_process_model_destroy(p); return status;
}
