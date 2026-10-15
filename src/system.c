#include "internal.h"
#include "chem_system.h"
#include <stdio.h>

ChemStatus chem_species_init(ChemSpecies *s, const char *name) {
  if (!s || !name) return CHEM_INVALID;
  memset(s, 0, sizeof *s); s->name = chem_strdup(name);
  s->absolute_tolerance = 1e-3; s->diffusion = s->molecular_weight = NAN; s->variable = SIZE_MAX;
  return s->name ? CHEM_OK : CHEM_MEMORY;
}
void chem_species_destroy(ChemSpecies *s) {
  if (!s) return;
  free(s->name);
  for (size_t i = 0; i < s->nproperties; ++i) { free(s->properties[i].key); free(s->properties[i].text); }
  free(s->properties); memset(s, 0, sizeof *s);
}
ChemStatus chem_species_set_property(ChemSpecies *s, const char *key, ChemPropertyType type, const void *value) {
  if (!s || !key || !value || (unsigned)type > CHEM_PROPERTY_STRING) return CHEM_INVALID;
  size_t index = s->nproperties;
  for (size_t i = 0; i < s->nproperties; ++i) if (s->properties[i].type == type && !strcmp(s->properties[i].key,key)) { index = i; break; }
  ChemProperty p = {0}; p.key = chem_strdup(key); p.type = type;
  if (type == CHEM_PROPERTY_STRING) p.text = chem_strdup(value);
  else if (type == CHEM_PROPERTY_DOUBLE) p.number = *(const double *)value;
  else p.integer = *(const int *)value;
  if (!p.key || (type == CHEM_PROPERTY_STRING && !p.text) ||
      (index == s->nproperties && !chem_resize((void **)&s->properties,index+1,sizeof p))) {
    free(p.key); free(p.text); return CHEM_MEMORY;
  }
  if (index < s->nproperties) { free(s->properties[index].key); free(s->properties[index].text); }
  else ++s->nproperties;
  s->properties[index] = p;
  if (type == CHEM_PROPERTY_DOUBLE) {
    if (!strcmp(key,"molecular weight [kg mol-1]")) s->molecular_weight = p.number;
    if (!strcmp(key,"diffusion coefficient [m2 s-1]")) s->diffusion = p.number;
    if (!strcmp(key,"absolute tolerance")) s->absolute_tolerance = p.number;
  }
  return CHEM_OK;
}
ChemStatus chem_species_get_property(const ChemSpecies *s, const char *key, ChemPropertyType type, void *value) {
  if (!s || !key || !value || (unsigned)type > CHEM_PROPERTY_STRING) return CHEM_INVALID;
  for (size_t i = 0; i < s->nproperties; ++i) {
    const ChemProperty *p = &s->properties[i];
    if (p->type == type && !strcmp(p->key,key)) {
      if (type == CHEM_PROPERTY_STRING) *(const char **)value = p->text;
      else if (type == CHEM_PROPERTY_DOUBLE) *(double *)value = p->number;
      else *(int *)value = p->integer;
      return CHEM_OK;
    }
  }
  return CHEM_UNKNOWN;
}
ChemStatus chem_species_copy(ChemSpecies *destination, const ChemSpecies *source) {
  if (!destination || !source) return CHEM_INVALID;
  if (destination == source) return CHEM_OK;
  ChemSpecies s = *source; s.name = chem_strdup(source->name); s.properties = NULL; s.nproperties = 0;
  if (!s.name) return CHEM_MEMORY;
  for (size_t i = 0; i < source->nproperties; ++i) {
    ChemProperty *p = &source->properties[i];
    const void *v = p->type == CHEM_PROPERTY_STRING ? (const void *)p->text :
      (p->type == CHEM_PROPERTY_DOUBLE ? (const void *)&p->number : (const void *)&p->integer);
    ChemStatus code = chem_species_set_property(&s,p->key,p->type,v);
    if (code) { chem_species_destroy(&s); return code; }
  }
  chem_species_destroy(destination); *destination = s; return CHEM_OK;
}
double chem_species_parameter(const ChemSpecies *s, const ChemConditions *c) {
  if (!s || !c) return NAN;
  if (s->parameterized == 1) return c->air_density;
  if (s->parameterized == 2 && s->parameterization) return s->parameterization(c,NULL,s->context);
  return NAN;
}
ChemStatus chem_phase_init(ChemPhase *p, const char *name, const ChemSpecies *species, size_t count) {
  if (!p || !name || (count && !species)) return CHEM_INVALID;
  memset(p,0,sizeof *p); p->name = chem_strdup(name); p->species = calloc(count ? count : 1,sizeof(ChemSpecies));
  if (!p->name || !p->species) { chem_phase_destroy(p); return CHEM_MEMORY; }
  for (size_t i = 0; i < count; ++i) {
    ChemStatus code = chem_species_copy(&p->species[i],&species[i]);
    if (code) { chem_phase_destroy(p); return code; }
    ++p->count;
  }
  return CHEM_OK;
}
void chem_phase_destroy(ChemPhase *p) {
  if (!p) return;
  for (size_t i = 0; i < p->count; ++i) chem_species_destroy(&p->species[i]);
  free(p->species); free(p->name); memset(p,0,sizeof *p);
}
size_t chem_phase_state_size(const ChemPhase *p) {
  size_t n = 0;
  if (p) for (size_t i = 0; i < p->count; ++i) n += !p->species[i].parameterized;
  return n;
}
void chem_names_free(char **names, size_t count) {
  for (size_t i = 0; i < count; ++i) free(names[i]);
  free(names);
}
ChemStatus chem_phase_names(const ChemPhase *p, int qualified, const size_t *permutation, char ***names, size_t *count) {
  if (!p || !names || !count) return CHEM_INVALID;
  size_t n = chem_phase_state_size(p);
  char **result = calloc(n ? n : 1,sizeof(char *));
  const char **raw = malloc((n ? n : 1)*sizeof(char *));
  if (!result || !raw) { free(result); free(raw); return CHEM_MEMORY; }
  size_t index = 0;
  for (size_t i = 0; i < p->count; ++i) if (!p->species[i].parameterized) raw[index++] = p->species[i].name;
  for (size_t i = 0; i < n; ++i) {
    size_t j = permutation ? permutation[i] : i;
    if (j >= n) { chem_names_free(result,n); free(raw); return CHEM_INVALID; }
    size_t length = strlen(raw[j])+1+(qualified ? strlen(p->name)+1 : 0);
    result[i] = malloc(length);
    if (!result[i]) { chem_names_free(result,n); free(raw); return CHEM_MEMORY; }
    if (qualified) snprintf(result[i],length,"%s.%s",p->name,raw[j]); else memcpy(result[i],raw[j],length);
  }
  free(raw); *names = result; *count = n; return CHEM_OK;
}
const char *chem_version(void) { return "1.0.0.0"; }
unsigned chem_version_component(unsigned component) { return component == 0 ? 1 : 0; }
