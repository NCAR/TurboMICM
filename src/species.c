#include "species.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static char *copy_string(const char *s)
{
    if (!s) return NULL;
    size_t n=strlen(s);
    char *out=malloc(n+1);
    if (out) memcpy(out,s,n+1);
    return out;
}
void chem_species_destroy(ChemSpecies *s)
{
    if (!s) return;
    for (size_t i=0;i<s->property_count;++i) {
        free(s->properties[i].name);
        if (s->properties[i].type==CHEM_PROPERTY_STRING) free((void *)s->properties[i].value.string);
    }
    free(s->properties); free(s->name); *s=(ChemSpecies){0};
}
ChemSpeciesStatus chem_species_init(ChemSpecies *s,const char *name)
{
    if (!s || !name) return CHEM_SPECIES_INVALID;
    char *copy=copy_string(name); if (!copy) return CHEM_SPECIES_NO_MEMORY;
    chem_species_destroy(s); s->name=copy; return CHEM_SPECIES_OK;
}
ChemSpeciesStatus chem_species_set_property(ChemSpecies *s,const char *name,
                                            ChemPropertyType type,ChemPropertyValue value)
{
    if (!s || !s->name || !name || type<CHEM_PROPERTY_DOUBLE || type>CHEM_PROPERTY_BOOL ||
        (type==CHEM_PROPERTY_STRING && !value.string)) return CHEM_SPECIES_INVALID;
    char *key=copy_string(name), *text=type==CHEM_PROPERTY_STRING ? copy_string(value.string) : NULL;
    if (!key || (type==CHEM_PROPERTY_STRING && !text)) { free(key); free(text); return CHEM_SPECIES_NO_MEMORY; }
    size_t index=0;
    for (;index<s->property_count;++index) if (!strcmp(s->properties[index].name,name)) break;
    if (index==s->property_count) {
        if (index>=SIZE_MAX/sizeof(*s->properties)) { free(key); free(text); return CHEM_SPECIES_NO_MEMORY; }
        ChemProperty *updated=realloc(s->properties,(index+1)*sizeof(*updated));
        if (!updated) { free(key); free(text); return CHEM_SPECIES_NO_MEMORY; }
        s->properties=updated; ++s->property_count;
    } else {
        free(s->properties[index].name);
        if (s->properties[index].type==CHEM_PROPERTY_STRING) free((void *)s->properties[index].value.string);
    }
    if (type==CHEM_PROPERTY_STRING) value.string=text;
    s->properties[index]=(ChemProperty){key,type,value};
    return CHEM_SPECIES_OK;
}
ChemSpeciesStatus chem_species_get_property(const ChemSpecies *s,const char *name,
                                            ChemPropertyType type,ChemPropertyValue *value)
{
    if (!s || !name || !value) return CHEM_SPECIES_INVALID;
    for (size_t i=0;i<s->property_count;++i) if (!strcmp(name,s->properties[i].name)) {
        if (s->properties[i].type!=type) return CHEM_SPECIES_WRONG_TYPE;
        *value=s->properties[i].value; return CHEM_SPECIES_OK;
    }
    return CHEM_SPECIES_NOT_FOUND;
}
ChemSpeciesStatus chem_species_clone(ChemSpecies *destination,const ChemSpecies *source)
{
    if (!destination || !source || !source->name) return CHEM_SPECIES_INVALID;
    if (destination==source) return CHEM_SPECIES_OK;
    ChemSpecies copy={0}; ChemSpeciesStatus status=chem_species_init(&copy,source->name);
    if (status!=CHEM_SPECIES_OK) return status;
    copy.is_parameterized=source->is_parameterized; copy.parameterization=source->parameterization;
    for (size_t i=0;i<source->property_count;++i) {
        const ChemProperty *p=source->properties+i;
        status=chem_species_set_property(&copy,p->name,p->type,p->value);
        if (status!=CHEM_SPECIES_OK) { chem_species_destroy(&copy); return status; }
    }
    chem_species_destroy(destination); *destination=copy; return CHEM_SPECIES_OK;
}
void chem_species_set_parameterization(ChemSpecies *s,ChemParameterizedSpecies p)
{ if (s) { s->parameterization=p; s->is_parameterized=1; } }
void chem_species_set_third_body(ChemSpecies *s)
{ chem_species_set_parameterization(s,(ChemParameterizedSpecies){0,0,0,1}); }
double chem_species_evaluate(const ChemSpecies *s,ChemConditions c)
{
    const ChemParameterizedSpecies *p=&s->parameterization;
    return p->constant+p->temperature*c.temperature+p->pressure*c.pressure+p->air_density*c.air_density;
}
void chem_phase_species_destroy(ChemPhaseSpecies *entry)
{ if (entry) { chem_species_destroy(&entry->species); *entry=(ChemPhaseSpecies){0}; } }
ChemSpeciesStatus chem_phase_species_init(ChemPhaseSpecies *entry,const ChemSpecies *species)
{
    if (!entry || !species) return CHEM_SPECIES_INVALID;
    ChemPhaseSpecies copy={0}; ChemSpeciesStatus status=chem_species_clone(&copy.species,species);
    if (status!=CHEM_SPECIES_OK) return status;
    chem_phase_species_destroy(entry); *entry=copy; return CHEM_SPECIES_OK;
}
void chem_phase_species_set_diffusion(ChemPhaseSpecies *entry,double coefficient)
{ if (entry) { entry->has_diffusion=1; entry->diffusion=coefficient; } }
void chem_phase_destroy(ChemPhase *phase)
{
    if (!phase) return;
    for (size_t i=0;i<phase->count;++i) chem_phase_species_destroy(phase->species+i);
    free(phase->species); free(phase->name); *phase=(ChemPhase){0};
}
ChemSpeciesStatus chem_phase_init(ChemPhase *phase,const char *name,size_t count,const ChemPhaseSpecies *species)
{
    if (!phase || !name || (count && !species)) return CHEM_SPECIES_INVALID;
    if (count>SIZE_MAX/sizeof(*phase->species)) return CHEM_SPECIES_NO_MEMORY;
    ChemPhase copy={0}; copy.name=copy_string(name);
    copy.species=count ? calloc(count,sizeof(*copy.species)) : NULL;
    if (!copy.name || (count && !copy.species)) { chem_phase_destroy(&copy); return CHEM_SPECIES_NO_MEMORY; }
    for (size_t i=0;i<count;++i) {
        ChemSpeciesStatus status=chem_phase_species_init(copy.species+i,&species[i].species);
        if (status!=CHEM_SPECIES_OK) { chem_phase_destroy(&copy); return status; }
        copy.species[i].has_diffusion=species[i].has_diffusion; copy.species[i].diffusion=species[i].diffusion;
        ++copy.count;
    }
    chem_phase_destroy(phase); *phase=copy; return CHEM_SPECIES_OK;
}
ChemSpeciesStatus chem_phase_clone(ChemPhase *destination,const ChemPhase *source)
{
    if (!destination || !source || !source->name) return CHEM_SPECIES_INVALID;
    if (destination==source) return CHEM_SPECIES_OK;
    return chem_phase_init(destination,source->name,source->count,source->species);
}
size_t chem_phase_state_size(const ChemPhase *phase)
{
    size_t count=0;
    if (phase) for (size_t i=0;i<phase->count;++i) count+=!phase->species[i].species.is_parameterized;
    return count;
}
ChemSpeciesStatus chem_phase_state_name(const ChemPhase *phase,size_t index,int qualified,char *buffer,size_t capacity)
{
    if (!phase || !buffer) return CHEM_SPECIES_INVALID;
    for (size_t i=0;i<phase->count;++i) if (!phase->species[i].species.is_parameterized) {
        if (index--) continue;
        const char *name=phase->species[i].species.name;
        size_t prefix=qualified ? strlen(phase->name)+1 : 0, length=strlen(name);
        if (prefix>SIZE_MAX-length-1 || capacity<prefix+length+1) return CHEM_SPECIES_RANGE;
        if (qualified) { memcpy(buffer,phase->name,prefix-1); buffer[prefix-1]='.'; }
        memcpy(buffer+prefix,name,length+1); return CHEM_SPECIES_OK;
    }
    return CHEM_SPECIES_RANGE;
}
ChemSpeciesStatus chem_system_init(ChemSystem *system,const ChemPhase *gas)
{ return system ? chem_phase_clone(&system->gas,gas) : CHEM_SPECIES_INVALID; }
void chem_system_destroy(ChemSystem *system) { if (system) chem_phase_destroy(&system->gas); }
