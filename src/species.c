/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

ch_species_spec *ch_species_create(const char *name)
{
  if (!name || !*name) return NULL;
  ch_species_spec *s = calloc(1,sizeof(*s));
  if (!s) return NULL;
  s->name = ch_string_copy(name);
  if (!s->name) { free(s); return NULL; }
  return s;
}
void ch_species_free(ch_species_spec *s)
{
  if (!s) return;
  for (size_t i = 0; i < s->nproperties; ++i) {
    free(s->properties[i].key);
    if (s->properties[i].type == CH_PROPERTY_STRING) free(s->properties[i].value.string);
  }
  free(s->properties); free(s->name); free(s);
}
ch_species_spec *ch_species_clone(const ch_species_spec *s)
{
  if (!s) return NULL;
  ch_species_spec *copy = ch_species_create(s->name);
  if (!copy) return NULL;
  copy->parameterized = s->parameterized;
  memcpy(copy->coefficients,s->coefficients,sizeof(copy->coefficients));
  for (size_t i = 0; i < s->nproperties; ++i) {
    const ch_property *p = s->properties+i;
    const void *v = p->type == CH_PROPERTY_STRING ? (void *)p->value.string : (const void *)&p->value;
    if (ch_species_set_property(copy,p->key,p->type,v)) { ch_species_free(copy); return NULL; }
  }
  return copy;
}
const char *ch_species_name(const ch_species_spec *s) { return s ? s->name : NULL; }
size_t ch_species_property_count(const ch_species_spec *s, int type)
{
  size_t n = 0;
  if (s) for (size_t i = 0; i < s->nproperties; ++i) n += s->properties[i].type == type;
  return n;
}
int ch_species_set_property(ch_species_spec *s, const char *key, int type, const void *value)
{
  if (!s || !key || !value || type < CH_PROPERTY_REAL || type > CH_PROPERTY_BOOL) return CH_INVALID;
  size_t i = 0;
  for (; i < s->nproperties; ++i) if (s->properties[i].type == type && !strcmp(s->properties[i].key,key)) break;
  char *string = type == CH_PROPERTY_STRING ? ch_string_copy(value) : NULL;
  if (type == CH_PROPERTY_STRING && !string) return CH_MEMORY;
  if (i == s->nproperties) {
    ch_property *p = realloc(s->properties,(i+1)*sizeof(*p));
    if (!p) { free(string); return CH_MEMORY; }
    s->properties = p;
    p[i] = (ch_property){0}; p[i].key = ch_string_copy(key); p[i].type = type;
    if (!p[i].key) { free(string); return CH_MEMORY; }
    ++s->nproperties;
  }
  ch_property *p = s->properties+i;
  if (type == CH_PROPERTY_REAL) p->value.real = *(const double *)value;
  else if (type == CH_PROPERTY_STRING) { free(p->value.string); p->value.string = string; }
  else p->value.integer = type == CH_PROPERTY_BOOL ? *(const int *)value != 0 : *(const int *)value;
  return CH_OK;
}
int ch_species_get_property(const ch_species_spec *s, const char *key, int type, void *value)
{
  if (!s || !key || !value || type < CH_PROPERTY_REAL || type > CH_PROPERTY_BOOL) return CH_INVALID;
  for (size_t i = 0; i < s->nproperties; ++i) {
    const ch_property *p = s->properties+i;
    if (p->type != type || strcmp(p->key,key)) continue;
    if (type == CH_PROPERTY_REAL) *(double *)value = p->value.real;
    else if (type == CH_PROPERTY_STRING) *(const char **)value = p->value.string;
    else *(int *)value = p->value.integer;
    return CH_OK;
  }
  return CH_NOT_FOUND;
}
int ch_species_parameterize(ch_species_spec *s, int enabled, double a, double t, double p, double d)
{
  if (!s || !isfinite(a) || !isfinite(t) || !isfinite(p) || !isfinite(d)) return CH_INVALID;
  s->parameterized = enabled != 0;
  s->coefficients[0] = a; s->coefficients[1] = t; s->coefficients[2] = p; s->coefficients[3] = d;
  return CH_OK;
}
int ch_species_is_parameterized(const ch_species_spec *s) { return s && s->parameterized; }
double ch_species_parameter_value(const ch_species_spec *s, const ch_conditions *c)
{
  if (!s || !c) return NAN;
  return s->coefficients[0] + s->coefficients[1]*c->temperature + s->coefficients[2]*c->pressure + s->coefficients[3]*c->air_density;
}
ch_phase *ch_phase_create(const char *name)
{
  if (!name) return NULL;
  ch_phase *p = calloc(1,sizeof(*p)); if (!p) return NULL;
  p->name = ch_string_copy(name);
  if (!p->name) { free(p); return NULL; }
  return p;
}
void ch_phase_free(ch_phase *p)
{
  if (!p) return;
  for (size_t i = 0; i < p->count; ++i) { ch_species_free(p->species[i]); free(p->qualified[i]); }
  free(p->name); free(p->species); free(p->qualified); free(p);
}
const char *ch_phase_name(const ch_phase *p) { return p ? p->name : NULL; }
int ch_phase_add(ch_phase *p, const ch_species_spec *s)
{
  if (!p || !s) return CH_INVALID;
  for (size_t i = 0; i < p->count; ++i) if (!strcmp(s->name,p->species[i]->name)) return CH_INVALID;
  ch_species_spec *copy = ch_species_clone(s);
  char *qualified = malloc(strlen(p->name)+strlen(s->name)+2);
  if (!copy || !qualified) { ch_species_free(copy); free(qualified); return CH_MEMORY; }
  sprintf(qualified,"%s.%s",p->name,s->name);
  ch_species_spec **list = realloc(p->species,(p->count+1)*sizeof(*list));
  if (!list) { ch_species_free(copy); free(qualified); return CH_MEMORY; }
  p->species = list;
  char **names = realloc(p->qualified,(p->count+1)*sizeof(*names));
  if (!names) { ch_species_free(copy); free(qualified); return CH_MEMORY; }
  p->qualified = names; p->species[p->count] = copy; p->qualified[p->count] = qualified;
  return (int)p->count++;
}
size_t ch_phase_species_count(const ch_phase *p) { return p ? p->count : 0; }
size_t ch_phase_state_size(const ch_phase *p)
{
  size_t n = 0; if (p) for (size_t i = 0; i < p->count; ++i) n += !p->species[i]->parameterized;
  return n;
}
const ch_species_spec *ch_phase_species(const ch_phase *p, size_t i) { return p && i < p->count ? p->species[i] : NULL; }
const char *ch_phase_state_name(const ch_phase *p, size_t index, int qualified)
{
  if (!p) return NULL;
  for (size_t i = 0; i < p->count; ++i) if (!p->species[i]->parameterized) {
    if (!index--) return qualified ? p->qualified[i] : p->species[i]->name;
  }
  return NULL;
}
const ch_species_spec *ch_model_species_info(const ch_model *m, size_t i)
{ return m && i < m->ns ? m->species[i].spec : NULL; }
size_t ch_model_builder_species_count(const ch_model *m) { return m ? m->ns : 0; }
int ch_model_add_species_spec(ch_model *m, const ch_species_spec *s)
{
  if (!m || !s) return CH_INVALID;
  double atol = 1e-3;
  (void)ch_species_get_property(s,"absolute tolerance",CH_PROPERTY_REAL,&atol);
  ch_species_spec *copy = ch_species_clone(s); if (!copy) return CH_MEMORY;
  int id = ch_model_add_species(m,s->name,atol,s->parameterized);
  if (id < 0) { ch_species_free(copy); return id; }
  ch_species *dest = m->species+id;
  ch_species_free(dest->spec); dest->spec = copy;
  dest->has_weight = ch_species_get_property(s,"molecular weight [kg mol-1]",CH_PROPERTY_REAL,&dest->molecular_weight) == CH_OK;
  dest->has_diffusion = ch_species_get_property(s,"diffusion coefficient [m2 s-1]",CH_PROPERTY_REAL,&dest->diffusion) == CH_OK;
  return id;
}
int ch_model_add_phase(ch_model *m, const ch_phase *p, int gas)
{
  if (!m || !p) return CH_INVALID;
  if (m->finalized) return CH_FROZEN;
  for (size_t i = 0; i < p->count; ++i) {
    ch_species_spec *s = ch_species_clone(p->species[i]); if (!s) return CH_MEMORY;
    if (!gas) {
      free(s->name); s->name = ch_string_copy(p->qualified[i]);
      if (!s->name) { ch_species_free(s); return CH_MEMORY; }
    }
    int status = ch_model_add_species_spec(m,s); ch_species_free(s);
    if (status < 0) return status;
  }
  return CH_OK;
}
ch_model *ch_model_clone_builder(const ch_model *m)
{
  if (!m || m->finalized) return NULL;
  ch_model *copy = ch_model_create(); if (!copy) return NULL;
  for (size_t i = 0; i < m->ns; ++i) if (ch_model_add_species_spec(copy,m->species[i].spec) < 0) goto error;
  for (size_t i = 0; i < m->nparam; ++i) if (ch_model_add_parameter(copy,m->parameters[i].name,m->parameters[i].initial) < 0) goto error;
  for (size_t i = 0; i < m->nr; ++i) {
    const ch_reaction *r = m->reactions+i;
    if (ch_model_add_reaction(copy,r->nr,r->reactants,r->np,r->products,r->yields,&r->rate) < 0) goto error;
  }
  for (size_t i = 0; i < m->nconstraint; ++i) {
    const ch_constraint *c = m->constraints+i;
    int status = c->kind == 0 ? ch_model_add_linear(copy,c->name,c->row,c->nr,c->ids,c->coefficients,c->constant)
      : ch_model_add_equilibrium(copy,c->name,c->row,c->nr,c->ids,c->coefficients,c->np,c->ids+c->nr,c->coefficients+c->nr,c->constant,c->enthalpy);
    if (status < 0) goto error;
  }
  for (size_t i = 0; i < m->nexternal; ++i) if (ch_model_add_external(copy,m->externals+i) < 0) goto error;
  return copy;
error:
  ch_model_free(copy); return NULL;
}
