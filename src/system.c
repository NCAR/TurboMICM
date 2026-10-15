#include "internal.h"
#include <stdio.h>
ck_status ck_species_property(const ck_species *s, const char *name, ck_property_type type, void *out)
{
  if (!s || !name || !out)
    return CK_INVALID;
  for (size_t j = 0; j < s->nproperties; ++j)
    if (!strcmp(s->properties[j].name, name))
    {
      const ck_property *p = &s->properties[j];
      if (p->type != type)
        return CK_INVALID;
      switch (type)
      {
      case CK_PROPERTY_REAL:
        *(double *)out = p->value.real;
        break;
      case CK_PROPERTY_STRING:
        *(const char **)out = p->value.string;
        break;
      case CK_PROPERTY_INT:
        *(int *)out = p->value.integer;
        break;
      case CK_PROPERTY_BOOL:
        *(int *)out = p->value.boolean;
        break;
      default:
        return CK_INVALID;
      }
      return CK_OK;
    }
  return CK_UNKNOWN;
}
double ck_species_parameter(const ck_species *s, const ck_conditions *c)
{
  if (!s || !c)
    return NAN;
  if (s->parameterized == 1)
    return c->air_density;
  return s->parameterized == 2 && s->function ? s->function(c, NULL, s->context) : NAN;
}
size_t ck_phase_state_size(const ck_phase *phase)
{
  size_t count = 0;
  if (phase)
    for (size_t j = 0; j < phase->nspecies; ++j)
      count += !phase->species[j].parameterized;
  return count;
}
ck_status ck_phase_unique_name(const ck_phase *phase, size_t variable, char *buffer, size_t capacity)
{
  if (!phase || !phase->name || !buffer)
    return CK_INVALID;
  size_t index = 0;
  for (size_t j = 0; j < phase->nspecies; ++j)
    if (!phase->species[j].parameterized)
    {
      if (index++ != variable)
        continue;
      size_t required = strlen(phase->name) + strlen(phase->species[j].name) + 2;
      if (capacity < required)
        return CK_INVALID;
      snprintf(buffer, capacity, "%s.%s", phase->name, phase->species[j].name);
      return CK_OK;
    }
  return CK_UNKNOWN;
}
ck_status ck_model_validate_species_used(const ck_mechanism *d)
{
  if (!d || !d->nspecies || !d->species)
    return CK_INVALID;
  unsigned char *used = ck_allocate(d->nspecies, 1);
  if (!used)
    return CK_MEMORY;
  ck_status status = CK_OK;
  for (size_t r = 0; r < d->nreactions; ++r)
  {
    const ck_reaction *rx = &d->reactions[r];
    for (size_t j = 0; j < rx->nreactants + rx->nproducts; ++j)
    {
      size_t s = j < rx->nreactants ? rx->reactants[j] : rx->products[j - rx->nreactants].species;
      if (s >= d->nspecies)
      {
        status = CK_UNKNOWN;
        goto done;
      }
      used[s] = 1;
    }
  }
  for (size_t k = 0; k < d->nconstraints; ++k)
  {
    const ck_constraint *c = &d->constraints[k];
    if (c->algebraic_species >= d->nspecies)
    {
      status = CK_UNKNOWN;
      goto done;
    }
    used[c->algebraic_species] = 1;
    for (size_t j = 0; j < c->nreactants + c->nproducts; ++j)
    {
      size_t s = j < c->nreactants ? c->reactants[j].species : c->products[j - c->nreactants].species;
      if (s >= d->nspecies)
      {
        status = CK_UNKNOWN;
        goto done;
      }
      used[s] = 1;
    }
  }
  for (size_t j = 0; j < d->nspecies; ++j)
    if (!used[j])
    {
      status = CK_INVALID;
      break;
    }
done:
  free(used);
  return status;
}
ck_status ck_model_set_lambda(ck_model *m, const char *name, ck_rate_function fn, void *context)
{
  if (!m || !name || !fn)
    return CK_INVALID;
  for (size_t r = 0; r < m->nr; ++r)
    if (m->rates[r].type == CK_LAMBDA && m->rates[r].label && !strcmp(m->rates[r].label, name))
    {
      m->rates[r].function = fn;
      m->rates[r].context = context;
      return CK_OK;
    }
  return CK_UNKNOWN;
}
const char *ck_version(void)
{
  return "0.1.0";
}
int ck_version_major(void)
{
  return 0;
}
int ck_version_minor(void)
{
  return 1;
}
int ck_version_patch(void)
{
  return 0;
}
int ck_version_tweak(void)
{
  return 0;
}
