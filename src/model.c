#include "chem.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void chem_species_init(ChemSpecies *s, const char *name)
{
    *s = (ChemSpecies){0};
    s->name = name;
}

void chem_species_free(ChemSpecies *s)
{
    if (!s) return;
    free(s->properties);
    *s = (ChemSpecies){0};
}

ChemPropertyStatus chem_species_copy(ChemSpecies *destination, const ChemSpecies *source)
{
    if (!destination || !source) return CHEM_PROPERTY_NOT_FOUND;
    if (destination == source) return CHEM_PROPERTY_OK;
    ChemSpecies copy = *source;
    copy.properties = NULL;
    copy.property_capacity = copy.property_count;
    if (copy.property_count) {
        if (copy.property_count > SIZE_MAX/sizeof(*copy.properties)) return CHEM_PROPERTY_NO_MEMORY;
        copy.properties = malloc(copy.property_count*sizeof(*copy.properties));
        if (!copy.properties) return CHEM_PROPERTY_NO_MEMORY;
        memcpy(copy.properties,source->properties,copy.property_count*sizeof(*copy.properties));
    }
    free(destination->properties);
    *destination=copy;
    return CHEM_PROPERTY_OK;
}

ChemPropertyStatus chem_species_set(ChemSpecies *s, const char *name, ChemPropertyType type, ChemPropertyValue value)
{
    if (type < CHEM_PROPERTY_REAL || type > CHEM_PROPERTY_BOOL) return CHEM_PROPERTY_INVALID_TYPE;
    if (!s || !name) return CHEM_PROPERTY_NOT_FOUND;
    for (size_t i=0;i<s->property_count;++i) if (s->properties[i].type==type && !strcmp(s->properties[i].name,name)) {
        s->properties[i].value=value; return CHEM_PROPERTY_OK;
    }
    if (s->property_count==s->property_capacity) {
        size_t capacity=s->property_capacity ? s->property_capacity*2 : 4;
        if (capacity < s->property_capacity || capacity > SIZE_MAX/sizeof(*s->properties)) return CHEM_PROPERTY_NO_MEMORY;
        ChemProperty *properties=realloc(s->properties,capacity*sizeof(*properties));
        if (!properties) return CHEM_PROPERTY_NO_MEMORY;
        s->properties=properties; s->property_capacity=capacity;
    }
    s->properties[s->property_count++]=(ChemProperty){name,type,value};
    return CHEM_PROPERTY_OK;
}

ChemPropertyStatus chem_species_get(const ChemSpecies *s, const char *name, ChemPropertyType type, ChemPropertyValue *value)
{
    if (type < CHEM_PROPERTY_REAL || type > CHEM_PROPERTY_BOOL) return CHEM_PROPERTY_INVALID_TYPE;
    if (!s || !name || !value) return CHEM_PROPERTY_NOT_FOUND;
    for (size_t i=0;i<s->property_count;++i) if (s->properties[i].type==type && !strcmp(s->properties[i].name,name)) {
        *value=s->properties[i].value; return CHEM_PROPERTY_OK;
    }
    return CHEM_PROPERTY_NOT_FOUND;
}

void chem_species_set_third_body(ChemSpecies *s)
{
    s->is_parameterized=1;
    s->parameterization=(ChemParameterization){0,0,0,1};
}

void chem_phase_species_set_diffusion(ChemPhaseSpecies *species, double diffusion)
{
    species->has_diffusion=1;
    species->diffusion=diffusion;
}

size_t chem_phase_state_size(ChemPhase phase)
{
    size_t count=0;
    for (size_t i=0;i<phase.count;++i) if (!phase.species[i].species->is_parameterized) ++count;
    return count;
}

int chem_phase_state_name(ChemPhase phase, size_t index, int qualified, char *output, size_t capacity)
{
    if (!output || !capacity) return 1;
    for (size_t i=0;i<phase.count;++i) {
        const ChemSpecies *s=phase.species[i].species;
        if (s->is_parameterized) continue;
        if (index-- == 0) {
            int written=qualified ? snprintf(output,capacity,"%s.%s",phase.name,s->name) : snprintf(output,capacity,"%s",s->name);
            return written < 0 || (size_t)written >= capacity;
        }
    }
    return 1;
}

size_t chem_system_state_size(ChemSystem system) { return chem_phase_state_size(system.gas_phase); }

int chem_system_state_name(ChemSystem system, size_t index, const size_t *order, char *output, size_t capacity)
{
    if (index >= chem_system_state_size(system)) return 1;
    return chem_phase_state_name(system.gas_phase,order ? order[index] : index,0,output,capacity);
}
