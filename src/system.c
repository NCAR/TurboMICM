#include "chem_system.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct ChemSpecies {
    char *name;
    size_t count;
    ChemProperty *properties;
    ChemParameterization parameterization;
};
struct ChemPhase { char *name; size_t count; ChemPhaseSpecies *members; };

static char *copy_string(const char *source)
{
    size_t length = strlen(source) + 1;
    char *result = malloc(length);
    if (result) memcpy(result, source, length);
    return result;
}
static void free_property(ChemProperty *property)
{
    free((char *)property->name);
    if (property->type == CHEM_PROPERTY_TEXT) free((char *)property->value.text);
}
ChemStatus chem_species_set_property(ChemSpecies *s, ChemProperty property)
{
    if (!s || !property.name || property.type < CHEM_PROPERTY_REAL || property.type > CHEM_PROPERTY_TEXT
        || (property.type == CHEM_PROPERTY_TEXT && !property.value.text)) return CHEM_INVALID;
    size_t position = 0;
    for (; position < s->count; ++position)
        if (s->properties[position].type == property.type && !strcmp(s->properties[position].name, property.name)) break;
    ChemProperty copy = property;
    copy.name = copy_string(property.name);
    if (!copy.name) return CHEM_MEMORY;
    if (property.type == CHEM_PROPERTY_TEXT) {
        copy.value.text = copy_string(property.value.text);
        if (!copy.value.text) { free((char *)copy.name); return CHEM_MEMORY; }
    }
    if (position == s->count) {
        if (s->count == SIZE_MAX / sizeof(ChemProperty)) { free_property(&copy); return CHEM_MEMORY; }
        ChemProperty *properties = realloc(s->properties, (s->count + 1) * sizeof(*properties));
        if (!properties) { free_property(&copy); return CHEM_MEMORY; }
        s->properties = properties;
        ++s->count;
    } else free_property(&s->properties[position]);
    s->properties[position] = copy;
    return CHEM_OK;
}
ChemStatus chem_species_create(const char *name, size_t count, const ChemProperty *properties, ChemSpecies **out)
{
    if (!out) return CHEM_INVALID;
    *out = NULL;
    if (!name || (count && !properties)) return CHEM_INVALID;
    ChemSpecies *s = calloc(1, sizeof(*s));
    if (!s) return CHEM_MEMORY;
    s->name = copy_string(name);
    if (!s->name) { free(s); return CHEM_MEMORY; }
    for (size_t i = 0; i < count; ++i) {
        ChemStatus status = chem_species_set_property(s, properties[i]);
        if (status != CHEM_OK) { chem_species_free(s); return status; }
    }
    *out = s;
    return CHEM_OK;
}
void chem_species_free(ChemSpecies *s)
{
    if (!s) return;
    for (size_t i = 0; i < s->count; ++i) free_property(&s->properties[i]);
    free(s->properties); free(s->name); free(s);
}
ChemStatus chem_species_clone(const ChemSpecies *s, ChemSpecies **out)
{
    if (!s) return CHEM_INVALID;
    ChemStatus status = chem_species_create(s->name, s->count, s->properties, out);
    if (status == CHEM_OK) (*out)->parameterization = s->parameterization;
    return status;
}
ChemStatus chem_species_assign(ChemSpecies **out, const ChemSpecies *source)
{
    if (!out) return CHEM_INVALID;
    ChemSpecies *copy = NULL;
    ChemStatus status = chem_species_clone(source, &copy);
    if (status == CHEM_OK) { chem_species_free(*out); *out = copy; }
    return status;
}
const char *chem_species_name(const ChemSpecies *s) { return s ? s->name : NULL; }
size_t chem_species_property_count(const ChemSpecies *s, ChemPropertyType type)
{
    size_t count = 0;
    if (s) for (size_t i = 0; i < s->count; ++i) count += s->properties[i].type == type;
    return count;
}
ChemStatus chem_species_get_property(const ChemSpecies *s, const char *name, ChemPropertyType type, ChemPropertyValue *out)
{
    if (!s || !name || !out) return CHEM_INVALID;
    if (type < CHEM_PROPERTY_REAL || type > CHEM_PROPERTY_TEXT) return CHEM_WRONG_TYPE;
    for (size_t i = 0; i < s->count; ++i)
        if (s->properties[i].type == type && !strcmp(s->properties[i].name, name)) {
            *out = s->properties[i].value;
            return CHEM_OK;
        }
    return CHEM_NOT_FOUND;
}
void chem_species_parameterize(ChemSpecies *s, ChemParameterization p) { s->parameterization = p; }
ChemParameterization chem_species_parameterization(const ChemSpecies *s) { return s->parameterization; }
void chem_species_third_body(ChemSpecies *s) { s->parameterization = (ChemParameterization){.air_density = 1, .has_value = 1}; }

ChemStatus chem_phase_create(const char *name, size_t count, const ChemPhaseSpecies *members, ChemPhase **out)
{
    if (!out) return CHEM_INVALID;
    *out = NULL;
    if (!name || (count && !members) || count > SIZE_MAX / sizeof(*members)) return CHEM_INVALID;
    ChemPhase *p = calloc(1, sizeof(*p));
    if (!p) return CHEM_MEMORY;
    p->name = copy_string(name);
    p->members = calloc(count + (count == 0), sizeof(*members));
    if (!p->name || !p->members) { chem_phase_free(p); return CHEM_MEMORY; }
    for (size_t i = 0; i < count; ++i) {
        ChemSpecies *copy = NULL;
        ChemStatus status = chem_species_clone(members[i].species, &copy);
        if (status != CHEM_OK) { chem_phase_free(p); return status; }
        p->members[p->count++] = (ChemPhaseSpecies){copy, members[i].diffusion, members[i].has_diffusion};
    }
    *out = p;
    return CHEM_OK;
}
void chem_phase_free(ChemPhase *p)
{
    if (!p) return;
    for (size_t i = 0; i < p->count; ++i) chem_species_free(p->members[i].species);
    free(p->members); free(p->name); free(p);
}
const char *chem_phase_name(const ChemPhase *p) { return p ? p->name : NULL; }
size_t chem_phase_count(const ChemPhase *p) { return p ? p->count : 0; }
const ChemPhaseSpecies *chem_phase_species(const ChemPhase *p, size_t i)
{ return p && i < p->count ? &p->members[i] : NULL; }
size_t chem_phase_state_size(const ChemPhase *p)
{
    size_t count = 0;
    if (p) for (size_t i = 0; i < p->count; ++i) count += !p->members[i].species->parameterization.has_value;
    return count;
}
size_t chem_phase_variable_name(const ChemPhase *p, size_t variable, int prefix,
    const size_t *permutation, char *buffer, size_t capacity)
{
    size_t size = chem_phase_state_size(p);
    if (variable >= size || (capacity && !buffer)) return SIZE_MAX;
    if (permutation) variable = permutation[variable];
    if (variable >= size) return SIZE_MAX;
    for (size_t i = 0; i < p->count; ++i) if (!p->members[i].species->parameterization.has_value) {
        if (variable-- == 0) {
            int length = prefix ? snprintf(buffer, capacity, "%s.%s", p->name, p->members[i].species->name)
                                : snprintf(buffer, capacity, "%s", p->members[i].species->name);
            return length < 0 ? SIZE_MAX : (size_t)length;
        }
    }
    return SIZE_MAX;
}
ChemStatus chem_surface_rate(const ChemPhaseSpecies *species, double probability, int parameter, ChemRate *rate)
{
    if (!species || !species->species || !rate || !species->has_diffusion || !(species->diffusion > 0)
        || parameter < 0) return CHEM_INVALID;
    ChemPropertyValue molecular_weight;
    ChemStatus status = chem_species_get_property(species->species, "molecular weight [kg mol-1]", CHEM_PROPERTY_REAL, &molecular_weight);
    if (status != CHEM_OK) return status;
    if (!(molecular_weight.real > 0)) return CHEM_INVALID;
    *rate = chem_rate_default(CHEM_SURFACE);
    rate->parameter = parameter;
    rate->p[0] = species->diffusion; rate->p[1] = molecular_weight.real; rate->p[2] = probability;
    return CHEM_OK;
}
