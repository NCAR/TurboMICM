#include "chem_system.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *copy_string(const char *s)
{
    if (!s) return NULL;
    size_t size=strlen(s)+1;
    char *copy=malloc(size);
    if (copy) memcpy(copy,s,size);
    return copy;
}
static int valid_type(ChemPropertyType type)
{ return type>=CHEM_PROPERTY_STRING && type<=CHEM_PROPERTY_BOOL; }

int chem_species_info_init(ChemSpeciesInfo *species,const char *name)
{
    if (!species || !name || species->name || species->properties || species->property_count) return CHEM_INVALID;
    char *copy=copy_string(name);
    if (!copy) return CHEM_MEMORY;
    *species=(ChemSpeciesInfo){.name=copy};
    return CHEM_OK;
}
void chem_species_info_destroy(ChemSpeciesInfo *species)
{
    if (!species) return;
    for (size_t i=0;i<species->property_count;++i) {
        free(species->properties[i].name);
        if (species->properties[i].type==CHEM_PROPERTY_STRING) free((void *)species->properties[i].value.string);
    }
    free(species->properties); free(species->name); memset(species,0,sizeof(*species));
}
int chem_species_info_set(ChemSpeciesInfo *species,const char *name,ChemPropertyType type,ChemPropertyValue value)
{
    if (!species || !species->name || !name) return CHEM_INVALID;
    if (!valid_type(type)) return CHEM_INVALID_PROPERTY_TYPE;
    if (type==CHEM_PROPERTY_STRING && !value.string) return CHEM_INVALID;
    char *string_value=NULL;
    if (type==CHEM_PROPERTY_STRING && !(string_value=copy_string(value.string))) return CHEM_MEMORY;
    for (size_t i=0;i<species->property_count;++i) {
        ChemProperty *property=&species->properties[i];
        if (property->type==type && !strcmp(property->name,name)) {
            if (type==CHEM_PROPERTY_STRING) {
                free((void *)property->value.string); value.string=string_value;
            }
            property->value=value; return CHEM_OK;
        }
    }
    char *name_copy=copy_string(name);
    if (!name_copy) { free(string_value); return CHEM_MEMORY; }
    if (species->property_count==species->property_capacity) {
        size_t capacity=species->property_capacity ? species->property_capacity*2 : 4;
        if (capacity<species->property_capacity || capacity>SIZE_MAX/sizeof(ChemProperty)) {
            free(name_copy); free(string_value); return CHEM_MEMORY;
        }
        ChemProperty *properties=realloc(species->properties,capacity*sizeof(*properties));
        if (!properties) { free(name_copy); free(string_value); return CHEM_MEMORY; }
        species->properties=properties; species->property_capacity=capacity;
    }
    if (type==CHEM_PROPERTY_STRING) value.string=string_value;
    species->properties[species->property_count++]=(ChemProperty){name_copy,type,value};
    return CHEM_OK;
}
int chem_species_info_get(const ChemSpeciesInfo *species,const char *name,ChemPropertyType type,ChemPropertyValue *value)
{
    if (!species || !species->name || !name || !value) return CHEM_INVALID;
    if (!valid_type(type)) return CHEM_INVALID_PROPERTY_TYPE;
    for (size_t i=0;i<species->property_count;++i)
        if (species->properties[i].type==type && !strcmp(species->properties[i].name,name)) {
            *value=species->properties[i].value; return CHEM_OK;
        }
    return CHEM_PROPERTY_NOT_FOUND;
}
int chem_species_info_has(const ChemSpeciesInfo *species,const char *name)
{
    if (!species || !name) return 0;
    for (size_t i=0;i<species->property_count;++i) if (!strcmp(species->properties[i].name,name)) return 1;
    return 0;
}
int chem_species_info_copy(ChemSpeciesInfo *destination,const ChemSpeciesInfo *source)
{
    if (!destination || !source || !source->name) return CHEM_INVALID;
    if (destination==source) return CHEM_OK;
    ChemSpeciesInfo copy={0};
    int status=chem_species_info_init(&copy,source->name);
    if (status) return status;
    copy.parameterized=source->parameterized; copy.parameterization=source->parameterization;
    for (size_t i=0;i<source->property_count;++i) {
        ChemProperty p=source->properties[i];
        status=chem_species_info_set(&copy,p.name,p.type,p.value);
        if (status) { chem_species_info_destroy(&copy); return status; }
    }
    chem_species_info_destroy(destination); *destination=copy;
    return CHEM_OK;
}
void chem_species_info_set_parameterization(ChemSpeciesInfo *species,ChemParameterizedFunction function)
{ if (species) { species->parameterized=1; species->parameterization=function; } }
void chem_species_info_set_third_body(ChemSpeciesInfo *species)
{ chem_species_info_set_parameterization(species,(ChemParameterizedFunction){.air_density=1}); }
double chem_parameterized_value(ChemParameterizedFunction f,double t,double p,double density)
{ return f.constant+f.temperature*t+f.pressure*p+f.air_density*density; }

size_t chem_phase_state_size(const ChemPhase *phase)
{
    if (!phase || (phase->species_count && !phase->species)) return 0;
    size_t count=0;
    for (size_t i=0;i<phase->species_count;++i)
        if (phase->species[i].species && !phase->species[i].species->parameterized) ++count;
    return count;
}
static const char *phase_name_at(const ChemPhase *phase,size_t index)
{
    if (!phase || (phase->species_count && !phase->species)) return NULL;
    for (size_t i=0;i<phase->species_count;++i) {
        const ChemSpeciesInfo *s=phase->species[i].species;
        if (s && !s->parameterized && index--==0) return s->name;
    }
    return NULL;
}
int chem_phase_state_name(const ChemPhase *phase,size_t state_index,char *buffer,size_t size)
{
    const char *name=phase_name_at(phase,state_index);
    if (!name || !phase->name || !buffer || !size) return CHEM_INVALID;
    int written=snprintf(buffer,size,"%s.%s",phase->name,name);
    return written<0 || (size_t)written>=size ? CHEM_INVALID : CHEM_OK;
}
int chem_phase_export_species(const ChemPhase *phase,ChemSpecies *output,size_t count)
{
    if (!phase || (phase->species_count && !phase->species) || count!=chem_phase_state_size(phase) ||
        (count && !output)) return CHEM_INVALID;
    size_t at=0;
    for (size_t i=0;i<phase->species_count;++i) {
        const ChemPhaseSpecies *item=&phase->species[i]; const ChemSpeciesInfo *info=item->species;
        if (!info || !info->name) return CHEM_INVALID;
        if (info->parameterized) continue;
        ChemSpecies s={info->name,0,1e-3,0,item->has_diffusion ? item->diffusion : 0}; ChemPropertyValue v;
        if (chem_species_info_get(info,"absolute tolerance",CHEM_PROPERTY_REAL,&v)==CHEM_OK) s.absolute_tolerance=v.real;
        if (chem_species_info_get(info,"molecular weight [kg mol-1]",CHEM_PROPERTY_REAL,&v)==CHEM_OK) s.molecular_weight=v.real;
        if (chem_species_info_get(info,"initial concentration [mol m-3]",CHEM_PROPERTY_REAL,&v)==CHEM_OK) s.initial=v.real;
        output[at++]=s;
    }
    return CHEM_OK;
}
size_t chem_system_state_size(const ChemSystem *system)
{ return system ? chem_phase_state_size(system->gas) : 0; }
int chem_system_state_name(const ChemSystem *system,size_t state_index,char *buffer,size_t size)
{
    const char *name=system ? phase_name_at(system->gas,state_index) : NULL;
    if (!name || !buffer || !size) return CHEM_INVALID;
    int written=snprintf(buffer,size,"%s",name);
    return written<0 || (size_t)written>=size ? CHEM_INVALID : CHEM_OK;
}
