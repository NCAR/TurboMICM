/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem_system.h"
#include "internal.h"
#include <string.h>

struct ChemSpeciesInfo { char *name;size_t count;ChemProperty *properties;ChemParameterization parameter; };
struct ChemPhase { char *name;size_t count;ChemPhaseSpecies *species; };
struct ChemSystem { ChemPhase *gas_phase; };
struct ChemReactionInfo { ChemPhase *phase;ChemReaction reaction; };
static char *duplicate(const char *s)
{
    size_t n=strlen(s);char *copy=malloc(n+1);if(copy)memcpy(copy,s,n+1);return copy;
}
static int property_type(ChemPropertyType t)
{ return t>=CHEM_PROPERTY_REAL && t<=CHEM_PROPERTY_BOOL; }
void chem_species_info_free(ChemSpeciesInfo *s)
{
    if(!s)return;
    for(size_t i=0;i<s->count;++i){
        free((void *)s->properties[i].name);
        if(s->properties[i].type==CHEM_PROPERTY_STRING)free((void *)s->properties[i].value.string);
    }
    free(s->properties);free(s->name);free(s);
}
ChemStatus chem_species_set_property(ChemSpeciesInfo *s,ChemProperty p)
{
    if(!s||!p.name||!property_type(p.type)||(p.type==CHEM_PROPERTY_STRING&&!p.value.string))return CHEM_INVALID;
    size_t i=0;while(i<s->count&&(s->properties[i].type!=p.type||strcmp(s->properties[i].name,p.name)))++i;
    char *value=p.type==CHEM_PROPERTY_STRING?duplicate(p.value.string):NULL;
    if(p.type==CHEM_PROPERTY_STRING&&!value)return CHEM_MEMORY;
    if(i==s->count){
        char *name=duplicate(p.name);
        if(!name||s->count>=SIZE_MAX/sizeof(ChemProperty)){free(name);free(value);return CHEM_MEMORY;}
        ChemProperty *props=realloc(s->properties,(s->count+1)*sizeof(*props));
        if(!props){free(name);free(value);return CHEM_MEMORY;}
        s->properties=props;p.name=name;++s->count;
    }else{
        p.name=s->properties[i].name;
        if(p.type==CHEM_PROPERTY_STRING)free((void *)s->properties[i].value.string);
    }
    if(p.type==CHEM_PROPERTY_STRING)p.value.string=value;
    s->properties[i]=p;return CHEM_OK;
}
ChemStatus chem_species_info_create(const char *name,size_t n,const ChemProperty *props,ChemSpeciesInfo **out)
{
    if(!name||!out||(n&&!props))return CHEM_INVALID;
    ChemSpeciesInfo *s=calloc(1,sizeof(*s));if(!s)return CHEM_MEMORY;
    s->name=duplicate(name);if(!s->name){free(s);return CHEM_MEMORY;}
    for(size_t i=0;i<n;++i){ChemStatus status=chem_species_set_property(s,props[i]);if(status!=CHEM_OK){chem_species_info_free(s);return status;}}
    *out=s;return CHEM_OK;
}
ChemStatus chem_species_info_copy(const ChemSpeciesInfo *s,ChemSpeciesInfo **out)
{
    if(!s||!out)return CHEM_INVALID;
    ChemStatus status=chem_species_info_create(s->name,s->count,s->properties,out);
    if(status==CHEM_OK)(*out)->parameter=s->parameter;
    return status;
}
ChemStatus chem_species_info_assign(ChemSpeciesInfo **dst,const ChemSpeciesInfo *src)
{
    if(!dst||!src)return CHEM_INVALID;
    if(*dst==src)return CHEM_OK;
    ChemSpeciesInfo *copy=NULL;ChemStatus status=chem_species_info_copy(src,&copy);
    if(status==CHEM_OK){chem_species_info_free(*dst);*dst=copy;}return status;
}
const char *chem_species_info_name(const ChemSpeciesInfo *s){return s?s->name:NULL;}
size_t chem_species_property_count(const ChemSpeciesInfo *s,ChemPropertyType t)
{size_t n=0;if(s)for(size_t i=0;i<s->count;++i)n+=s->properties[i].type==t;return n;}
int chem_species_has_property(const ChemSpeciesInfo *s,const char *name)
{if(s&&name)for(size_t i=0;i<s->count;++i)if(!strcmp(s->properties[i].name,name))return 1;return 0;}
ChemStatus chem_species_get_property(const ChemSpeciesInfo *s,const char *name,ChemPropertyType t,ChemPropertyValue *out)
{
    if(!s||!name||!out||!property_type(t))return CHEM_INVALID;
    for(size_t i=0;i<s->count;++i)if(s->properties[i].type==t&&!strcmp(s->properties[i].name,name)){*out=s->properties[i].value;return CHEM_OK;}
    return CHEM_NOT_FOUND;
}
void chem_species_parameterize(ChemSpeciesInfo *s,ChemParameterization p){if(s)s->parameter=p;}
void chem_species_set_third_body(ChemSpeciesInfo *s)
{chem_species_parameterize(s,(ChemParameterization){.air_density=1,.enabled=1});}
int chem_species_is_parameterized(const ChemSpeciesInfo *s){return s&&s->parameter.enabled;}
double chem_species_parameter_value(const ChemSpeciesInfo *s,const ChemConditions *c)
{
    const ChemParameterization *p=&s->parameter;
    return p->constant+p->temperature*c->temperature+p->pressure*c->pressure+p->air_density*c->air_density;
}
ChemStatus chem_phase_species_create(const ChemSpeciesInfo *s,ChemPhaseSpecies *out)
{
    if(!out)return CHEM_INVALID;
    ChemPhaseSpecies p={0};ChemStatus status=chem_species_info_copy(s,&p.species);
    if(status==CHEM_OK)*out=p;
    return status;
}
void chem_phase_species_free(ChemPhaseSpecies *s)
{if(s){chem_species_info_free(s->species);*s=(ChemPhaseSpecies){0};}}
void chem_phase_species_diffusion(ChemPhaseSpecies *s,double v){if(s){s->diffusion=v;s->has_diffusion=1;}}
void chem_phase_species_density(ChemPhaseSpecies *s,double v){if(s){s->density=v;s->has_density=1;}}
void chem_phase_free(ChemPhase *p)
{if(p){if(p->species)for(size_t i=0;i<p->count;++i)chem_phase_species_free(p->species+i);free(p->species);free(p->name);free(p);}}
ChemStatus chem_phase_create(const char *name,size_t n,const ChemPhaseSpecies *species,ChemPhase **out)
{
    if(!name||!out||(n&&!species))return CHEM_INVALID;
    ChemPhase *p=calloc(1,sizeof(*p));if(!p)return CHEM_MEMORY;
    p->name=duplicate(name);p->count=n;p->species=chem_alloc(n,sizeof(*p->species));
    if(!p->name||!p->species){chem_phase_free(p);return CHEM_MEMORY;}
    for(size_t i=0;i<n;++i){
        ChemStatus status=chem_phase_species_create(species[i].species,p->species+i);
        if(status!=CHEM_OK){chem_phase_free(p);return status;}
        p->species[i].diffusion=species[i].diffusion;p->species[i].has_diffusion=species[i].has_diffusion;
        p->species[i].density=species[i].density;p->species[i].has_density=species[i].has_density;
    }
    *out=p;return CHEM_OK;
}
ChemStatus chem_phase_copy(const ChemPhase *p,ChemPhase **out)
{return p?chem_phase_create(p->name,p->count,p->species,out):CHEM_INVALID;}
const char *chem_phase_name(const ChemPhase *p){return p?p->name:NULL;}
size_t chem_phase_count(const ChemPhase *p){return p?p->count:0;}
size_t chem_phase_state_size(const ChemPhase *p)
{size_t n=0;if(p)for(size_t i=0;i<p->count;++i)n+=!chem_species_is_parameterized(p->species[i].species);return n;}
const ChemPhaseSpecies *chem_phase_species(const ChemPhase *p,size_t i){return p&&i<p->count?p->species+i:NULL;}
void chem_names_free(ChemNames *names)
{if(names){if(names->values)for(size_t i=0;i<names->count;++i)free(names->values[i]);free(names->values);*names=(ChemNames){0};}}
ChemStatus chem_phase_names(const ChemPhase *p,int prefix,ChemNames *out)
{
    if(!p||!out)return CHEM_INVALID;
    ChemNames names={chem_phase_state_size(p),NULL};names.values=chem_alloc(names.count,sizeof(char *));
    if(!names.values)return CHEM_MEMORY;
    size_t j=0;
    for(size_t i=0;i<p->count;++i){
        if(chem_species_is_parameterized(p->species[i].species))continue;
        const char *name=p->species[i].species->name;size_t n=strlen(name),len=prefix?strlen(p->name)+1:0;
        if(n>SIZE_MAX-len-1){chem_names_free(&names);return CHEM_MEMORY;}
        char *s=malloc(len+n+1);if(!s){chem_names_free(&names);return CHEM_MEMORY;}
        if(prefix){memcpy(s,p->name,len-1);s[len-1]='.';}memcpy(s+len,name,n+1);names.values[j++]=s;
    }
    *out=names;return CHEM_OK;
}
ChemStatus chem_system_create(const ChemPhase *p,ChemSystem **out)
{
    if(!p||!out)return CHEM_INVALID;
    ChemSystem *s=calloc(1,sizeof(*s));if(!s)return CHEM_MEMORY;
    ChemStatus status=chem_phase_copy(p,&s->gas_phase);if(status!=CHEM_OK){free(s);return status;}
    *out=s;return CHEM_OK;
}
void chem_system_free(ChemSystem *s){if(s){chem_phase_free(s->gas_phase);free(s);}}
const ChemPhase *chem_system_gas_phase(const ChemSystem *s){return s?s->gas_phase:NULL;}
size_t chem_system_state_size(const ChemSystem *s){return s?chem_phase_state_size(s->gas_phase):0;}
ChemStatus chem_system_names(const ChemSystem *s,ChemNameOrder order,void *ctx,ChemNames *out)
{
    if(!s||!out)return CHEM_INVALID;
    ChemNames names={0};ChemStatus status=chem_phase_names(s->gas_phase,0,&names);
    if(status!=CHEM_OK)return status;
    if(!order){*out=names;return CHEM_OK;}
    ChemNames sorted={names.count,chem_alloc(names.count,sizeof(char *))};
    if(!sorted.values){chem_names_free(&names);return CHEM_MEMORY;}
    for(size_t i=0;i<names.count;++i){
        const char *name=order((const char *const *)names.values,names.count,i,ctx);
        if(!name){status=CHEM_INVALID;break;}
        sorted.values[i]=duplicate(name);if(!sorted.values[i]){status=CHEM_MEMORY;break;}
    }
    chem_names_free(&names);if(status!=CHEM_OK){chem_names_free(&sorted);return status;}
    *out=sorted;return CHEM_OK;
}
static ChemStatus phase_index(const ChemPhase *p,const char *name,size_t *index)
{
    if(!name)return CHEM_INVALID;
    for(size_t i=0;i<p->count;++i)if(!strcmp(p->species[i].species->name,name)){*index=i;return CHEM_OK;}
    return CHEM_NOT_FOUND;
}
void chem_reaction_info_free(ChemReactionInfo *r)
{if(r){chem_phase_free(r->phase);free((void *)r->reaction.reactants);free((void *)r->reaction.products);free((void *)r->reaction.rate.label);free(r);}}
ChemStatus chem_reaction_info_create(const ChemPhase *p,size_t nr,const char *const *reactants,
    size_t np,const ChemNamedProduct *products,ChemRate rate,ChemReactionInfo **out)
{
    if(!p||!out||(nr&&!reactants)||(np&&!products))return CHEM_INVALID;
    ChemReactionInfo *r=calloc(1,sizeof(*r));if(!r)return CHEM_MEMORY;
    ChemStatus status=chem_phase_copy(p,&r->phase);if(status!=CHEM_OK){free(r);return status;}
    size_t *rs=chem_alloc(nr,sizeof(*rs));ChemProduct *ps=chem_alloc(np,sizeof(*ps));
    r->reaction=(ChemReaction){nr,np,rs,ps,rate};r->reaction.rate.label=duplicate(rate.label?rate.label:"");
    if(!rs||!ps||!r->reaction.rate.label){chem_reaction_info_free(r);return CHEM_MEMORY;}
    for(size_t i=0;i<nr;++i){status=phase_index(p,reactants[i],rs+i);if(status!=CHEM_OK)goto fail;}
    for(size_t i=0;i<np;++i){status=phase_index(p,products[i].name,&ps[i].species);ps[i].yield=products[i].yield;if(status!=CHEM_OK)goto fail;}
    *out=r;return CHEM_OK;
fail:chem_reaction_info_free(r);return status;
}
const ChemPhase *chem_reaction_phase(const ChemReactionInfo *r){return r?r->phase:NULL;}
const ChemReaction *chem_reaction_definition(const ChemReactionInfo *r){return r?&r->reaction:NULL;}
const char *chem_reaction_reactant(const ChemReactionInfo *r,size_t i)
{return r&&i<r->reaction.reactant_count?r->phase->species[r->reaction.reactants[i]].species->name:NULL;}
const char *chem_reaction_product(const ChemReactionInfo *r,size_t i)
{return r&&i<r->reaction.product_count?r->phase->species[r->reaction.products[i].species].species->name:NULL;}
ChemStatus chem_reaction_info_copy(const ChemReactionInfo *r,ChemReactionInfo **out)
{
    if(!r||!out)return CHEM_INVALID;
    size_t nr=r->reaction.reactant_count,np=r->reaction.product_count;
    const char **rs=chem_alloc(nr,sizeof(*rs));ChemNamedProduct *ps=chem_alloc(np,sizeof(*ps));
    if(!rs||!ps){free(rs);free(ps);return CHEM_MEMORY;}
    for(size_t i=0;i<nr;++i)rs[i]=chem_reaction_reactant(r,i);
    for(size_t i=0;i<np;++i)ps[i]=(ChemNamedProduct){chem_reaction_product(r,i),r->reaction.products[i].yield};
    ChemStatus status=chem_reaction_info_create(r->phase,nr,rs,np,ps,r->reaction.rate,out);free(rs);free(ps);return status;
}
ChemStatus chem_reaction_info_assign(ChemReactionInfo **dst,const ChemReactionInfo *src)
{
    if(!dst||!src)return CHEM_INVALID;
    if(*dst==src)return CHEM_OK;
    ChemReactionInfo *copy=NULL;ChemStatus s=chem_reaction_info_copy(src,&copy);
    if(s==CHEM_OK){chem_reaction_info_free(*dst);*dst=copy;}return s;
}
static double parameter(const ChemConditions *c,void *ctx)
{return chem_species_parameter_value(ctx,c);}
ChemStatus chem_system_mechanism(const ChemSystem *s,size_t n,const ChemReactionInfo *const *reactions,ChemMechanism **out)
{
    if(!s||!out||(n&&!reactions))return CHEM_INVALID;
    ChemSystem *owned=NULL;ChemStatus status=chem_system_create(s->gas_phase,&owned);if(status!=CHEM_OK)return status;
    ChemPhase *phase=owned->gas_phase;
    ChemSpecies *species=chem_alloc(phase->count,sizeof(*species));ChemReaction *rs=chem_alloc(n,sizeof(*rs));
    if(!species||!rs){free(species);free(rs);chem_system_free(owned);return CHEM_MEMORY;}
    for(size_t i=0;i<phase->count;++i){
        ChemSpeciesInfo *info=phase->species[i].species;ChemPropertyValue v;
        species[i].name=info->name;species[i].absolute_tolerance=1e-3;
        if(chem_species_get_property(info,"absolute tolerance",CHEM_PROPERTY_REAL,&v)==CHEM_OK)species[i].absolute_tolerance=v.real;
        if(chem_species_get_property(info,"molecular weight [kg mol-1]",CHEM_PROPERTY_REAL,&v)==CHEM_OK)species[i].molecular_weight=v.real;
        if(phase->species[i].has_diffusion)species[i].diffusion=phase->species[i].diffusion;
        if(chem_species_is_parameterized(info)){species[i].parameter=parameter;species[i].parameter_context=info;}
    }
    for(size_t i=0;i<n;++i){
        if(!reactions[i]){status=CHEM_INVALID;goto done;}
        const ChemReaction *r=&reactions[i]->reaction;size_t *reactants=chem_alloc(r->reactant_count,sizeof(size_t));
        ChemProduct *products=chem_alloc(r->product_count,sizeof(ChemProduct));rs[i]=*r;rs[i].reactants=reactants;rs[i].products=products;
        if(!reactants||!products){status=CHEM_MEMORY;goto done;}
        for(size_t j=0;j<r->reactant_count;++j){status=phase_index(phase,chem_reaction_reactant(reactions[i],j),reactants+j);if(status!=CHEM_OK)goto done;}
        for(size_t j=0;j<r->product_count;++j){status=phase_index(phase,chem_reaction_product(reactions[i],j),&products[j].species);products[j].yield=r->products[j].yield;if(status!=CHEM_OK)goto done;}
    }
    status=chem_mechanism_create(phase->count,species,n,rs,out);
    if(status==CHEM_OK){(*out)->system=owned;owned=NULL;}
done:
    for(size_t i=0;i<n;++i){free((void *)rs[i].reactants);free((void *)rs[i].products);}
    free(rs);free(species);chem_system_free(owned);return status;
}
const ChemSystem *chem_mechanism_system(const ChemMechanism *m){return m?m->system:NULL;}
