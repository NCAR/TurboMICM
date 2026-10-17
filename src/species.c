/* Species, phases, and systems. */
#include "internal.h"

#include <stdio.h>

#define TMC_GAS_CONSTANT 8.31446261815324 /* J mol-1 K-1 */

void tmc_conditions_calculate_ideal_air_density(tmc_conditions *c)
{
  c->air_density = c->pressure / (TMC_GAS_CONSTANT * c->temperature);
}

tmc_species tmc_species_create(const char *name)
{
  tmc_species s;
  memset(&s, 0, sizeof s);
  tmc_strlcpy(s.name, name, sizeof s.name);
  return s;
}

tmc_species tmc_species_create_with_property(const char *name, const char *key, double value)
{
  tmc_species s = tmc_species_create(name);
  tmc_species_set_double(&s, key, value);
  return s;
}

void tmc_species_set_third_body(tmc_species *s)
{
  s->parameterized = true;
  s->c_rho = 1.0;
}

bool tmc_species_is_parameterized(const tmc_species *s)
{
  return s->parameterized;
}

double tmc_species_parameterized_value(const tmc_species *s, const tmc_conditions *c)
{
  return s->c_rho * c->air_density;
}

static tmc_property *find_property(const tmc_species *s, const char *key)
{
  for (int i = 0; i < s->n_properties; ++i)
    if (strcmp(s->properties[i].key, key) == 0)
      return (tmc_property *)&s->properties[i];
  return NULL;
}

bool tmc_species_has_property(const tmc_species *s, const char *key)
{
  return find_property(s, key) != NULL;
}

static tmc_property *property_slot(tmc_species *s, const char *key, tmc_property_type type)
{
  tmc_property *p = find_property(s, key);
  if (!p)
  {
    if (s->n_properties == TMC_MAX_PROPERTIES)
      return NULL;
    p = &s->properties[s->n_properties++];
    memset(p, 0, sizeof *p);
    tmc_strlcpy(p->key, key, sizeof p->key);
  }
  p->type = type;
  return p;
}

int tmc_species_set_double(tmc_species *s, const char *key, double value)
{
  tmc_property *p = property_slot(s, key, TMC_PROPERTY_DOUBLE);
  if (!p)
    return TMC_ERROR_SPECIES;
  p->value_double = value;
  return TMC_OK;
}

int tmc_species_set_string(tmc_species *s, const char *key, const char *value)
{
  tmc_property *p = property_slot(s, key, TMC_PROPERTY_STRING);
  if (!p)
    return TMC_ERROR_SPECIES;
  tmc_strlcpy(p->value_string, value, sizeof p->value_string);
  return TMC_OK;
}

int tmc_species_set_int(tmc_species *s, const char *key, int value)
{
  tmc_property *p = property_slot(s, key, TMC_PROPERTY_INT);
  if (!p)
    return TMC_ERROR_SPECIES;
  p->value_int = value;
  return TMC_OK;
}

int tmc_species_set_bool(tmc_species *s, const char *key, bool value)
{
  tmc_property *p = property_slot(s, key, TMC_PROPERTY_BOOL);
  if (!p)
    return TMC_ERROR_SPECIES;
  p->value_bool = value;
  return TMC_OK;
}

static int get_property(const tmc_species *s, const char *key, tmc_property_type type, const tmc_property **out, tmc_error *err)
{
  const tmc_property *p = find_property(s, key);
  if (!p)
    return tmc_error_set(err, TMC_ERROR_SPECIES, TMC_SPECIES_PROPERTY_NOT_FOUND, "Species: '%s' Property: '%s'", s->name, key);
  if (p->type != type)
    return tmc_error_set(err, TMC_ERROR_SPECIES, TMC_SPECIES_INVALID_TYPE_FOR_PROPERTY, "Species: '%s' Property: '%s'", s->name, key);
  *out = p;
  return TMC_OK;
}

int tmc_species_get_double(const tmc_species *s, const char *key, double *value, tmc_error *err)
{
  const tmc_property *p;
  int rc = get_property(s, key, TMC_PROPERTY_DOUBLE, &p, err);
  if (rc == TMC_OK)
    *value = p->value_double;
  return rc;
}

int tmc_species_get_string(const tmc_species *s, const char *key, const char **value, tmc_error *err)
{
  const tmc_property *p;
  int rc = get_property(s, key, TMC_PROPERTY_STRING, &p, err);
  if (rc == TMC_OK)
    *value = p->value_string;
  return rc;
}

int tmc_species_get_int(const tmc_species *s, const char *key, int *value, tmc_error *err)
{
  const tmc_property *p;
  int rc = get_property(s, key, TMC_PROPERTY_INT, &p, err);
  if (rc == TMC_OK)
    *value = p->value_int;
  return rc;
}

int tmc_species_get_bool(const tmc_species *s, const char *key, bool *value, tmc_error *err)
{
  const tmc_property *p;
  int rc = get_property(s, key, TMC_PROPERTY_BOOL, &p, err);
  if (rc == TMC_OK)
    *value = p->value_bool;
  return rc;
}

/* ---- phases ---- */

tmc_phase tmc_phase_create(const char *name)
{
  tmc_phase p;
  memset(&p, 0, sizeof p);
  tmc_strlcpy(p.name, name, sizeof p.name);
  return p;
}

tmc_phase tmc_phase_copy(const tmc_phase *src)
{
  tmc_phase p = *src;
  p.species = NULL;
  if (src->n_species)
  {
    p.species = tmc_malloc((size_t)src->n_species * sizeof *p.species);
    memcpy(p.species, src->species, (size_t)src->n_species * sizeof *p.species);
  }
  return p;
}

void tmc_phase_free(tmc_phase *p)
{
  free(p->species);
  p->species = NULL;
  p->n_species = 0;
}

static tmc_phase_species *phase_grow(tmc_phase *p)
{
  p->species = tmc_realloc(p->species, (size_t)(p->n_species + 1) * sizeof *p->species);
  tmc_phase_species *ps = &p->species[p->n_species++];
  memset(ps, 0, sizeof *ps);
  return ps;
}

void tmc_phase_add_species(tmc_phase *p, const tmc_species *s)
{
  tmc_phase_species *ps = phase_grow(p);
  ps->species = *s;
  double d;
  if (tmc_species_get_double(s, TMC_PROPERTY_DIFFUSION_COEFFICIENT, &d, NULL) == TMC_OK)
  {
    ps->has_diffusion_coefficient = true;
    ps->diffusion_coefficient = d;
  }
}

void tmc_phase_add_species_with_diffusion(tmc_phase *p, const tmc_species *s, double diffusion_coefficient)
{
  tmc_phase_species *ps = phase_grow(p);
  ps->species = *s;
  ps->has_diffusion_coefficient = true;
  ps->diffusion_coefficient = diffusion_coefficient;
}

int tmc_phase_state_size(const tmc_phase *p)
{
  int n = 0;
  for (int i = 0; i < p->n_species; ++i)
    if (!p->species[i].species.parameterized)
      ++n;
  return n;
}

const tmc_phase_species *tmc_phase_find(const tmc_phase *p, const char *name)
{
  for (int i = 0; i < p->n_species; ++i)
    if (strcmp(p->species[i].species.name, name) == 0)
      return &p->species[i];
  return NULL;
}

int tmc_phase_species_names(const tmc_phase *p, bool qualified, char (*names)[TMC_NAME_MAX], int max_names)
{
  int n = 0;
  for (int i = 0; i < p->n_species; ++i)
  {
    if (p->species[i].species.parameterized)
      continue;
    if (names && n < max_names)
    {
      if (qualified)
        snprintf(names[n], TMC_NAME_MAX, "%s.%s", p->name, p->species[i].species.name);
      else
        tmc_strlcpy(names[n], p->species[i].species.name, TMC_NAME_MAX);
    }
    ++n;
  }
  return n;
}

/* ---- systems ---- */

tmc_system tmc_system_create(const tmc_phase *gas_phase)
{
  tmc_system s;
  memset(&s, 0, sizeof s);
  s.gas_phase = tmc_phase_copy(gas_phase);
  return s;
}

void tmc_system_add_phase(tmc_system *s, const tmc_phase *p)
{
  s->phases = tmc_realloc(s->phases, (size_t)(s->n_phases + 1) * sizeof *s->phases);
  s->phases[s->n_phases++] = tmc_phase_copy(p);
}

tmc_system tmc_system_copy(const tmc_system *src)
{
  tmc_system s;
  memset(&s, 0, sizeof s);
  s.gas_phase = tmc_phase_copy(&src->gas_phase);
  for (int i = 0; i < src->n_phases; ++i)
    tmc_system_add_phase(&s, &src->phases[i]);
  return s;
}

void tmc_system_free(tmc_system *s)
{
  tmc_phase_free(&s->gas_phase);
  for (int i = 0; i < s->n_phases; ++i)
    tmc_phase_free(&s->phases[i]);
  free(s->phases);
  s->phases = NULL;
  s->n_phases = 0;
}

int tmc_system_state_size(const tmc_system *s)
{
  return tmc_phase_state_size(&s->gas_phase);
}

int tmc_system_unique_names(const tmc_system *s, char (*names)[TMC_NAME_MAX], int max_names)
{
  return tmc_phase_species_names(&s->gas_phase, false, names, max_names);
}
