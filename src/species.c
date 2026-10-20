#include "kinetics.h"

#include <stdlib.h>
#include <string.h>

static char *duplicate(const char *text)
{
  const size_t length = strlen(text);
  char *copy;
  if (length == SIZE_MAX)
    return NULL;
  copy = (char *)malloc(length + 1);
  if (copy)
    memcpy(copy, text, length + 1);
  return copy;
}

static void property_destroy(chem_species_property *property)
{
  free(property->name);
  if (property->kind == CHEM_PROPERTY_STRING)
    free((char *)property->value.string_value);
  *property = (chem_species_property){0};
}

static int property_kind_valid(chem_species_property_kind kind)
{
  return kind >= CHEM_PROPERTY_DOUBLE && kind <= CHEM_PROPERTY_BOOLEAN;
}

static int property_copy(chem_species_property *destination,
                         const chem_species_property *source)
{
  if (!source->name || !property_kind_valid(source->kind) ||
      (source->kind == CHEM_PROPERTY_STRING && !source->value.string_value))
    return -1;
  destination->name = duplicate(source->name);
  if (!destination->name)
    return -1;
  destination->kind = source->kind;
  destination->value = source->value;
  if (source->kind == CHEM_PROPERTY_STRING) {
    destination->value.string_value = duplicate(source->value.string_value);
    if (!destination->value.string_value) {
      free(destination->name);
      destination->name = NULL;
      return -1;
    }
  }
  return 0;
}

static int ensure_capacity(chem_species *species, size_t required)
{
  chem_species_property *properties;
  size_t capacity, i;
  if (!species->owns_properties) {
    capacity = required > 4 ? required : 4;
    if (capacity > SIZE_MAX / sizeof(*properties))
      return -1;
    properties = (chem_species_property *)calloc(capacity, sizeof(*properties));
    if (!properties)
      return -1;
    for (i = 0; i < species->property_count; ++i) {
      if (property_copy(&properties[i], &species->properties[i]) != 0) {
        while (i)
          property_destroy(&properties[--i]);
        free(properties);
        return -1;
      }
    }
    species->properties = properties;
    species->property_capacity = capacity;
    species->owns_properties = 1;
    return 0;
  }
  if (required <= species->property_capacity)
    return 0;
  capacity = species->property_capacity ? species->property_capacity : 4;
  while (capacity < required) {
    if (capacity > SIZE_MAX / 2) {
      capacity = required;
      break;
    }
    capacity *= 2;
  }
  if (capacity > SIZE_MAX / sizeof(*properties))
    return -1;
  properties = (chem_species_property *)realloc(
      species->properties, capacity * sizeof(*properties));
  if (!properties)
    return -1;
  memset(properties + species->property_capacity, 0,
         (capacity - species->property_capacity) * sizeof(*properties));
  species->properties = properties;
  species->property_capacity = capacity;
  return 0;
}

int chem_species_init(chem_species *species, const char *name)
{
  if (!species || !name)
    return -1;
  *species = (chem_species){0};
  species->name = name;
  return 0;
}

void chem_species_destroy(chem_species *species)
{
  size_t i;
  if (!species)
    return;
  if (species->owns_properties) {
    for (i = 0; i < species->property_count; ++i)
      property_destroy(&species->properties[i]);
    free(species->properties);
  }
  *species = (chem_species){0};
}

int chem_species_copy(chem_species *destination, const chem_species *source)
{
  chem_species copied;
  size_t i;
  if (!destination || !source ||
      (source->property_count && !source->properties))
    return -1;
  copied = *source;
  copied.properties = NULL;
  copied.property_capacity = 0;
  copied.owns_properties = 0;
  if (source->property_count) {
    copied.properties = (chem_species_property *)calloc(
        source->property_count, sizeof(*copied.properties));
    if (!copied.properties)
      return -1;
    copied.property_capacity = source->property_count;
    copied.owns_properties = 1;
    for (i = 0; i < source->property_count; ++i) {
      if (property_copy(&copied.properties[i], &source->properties[i]) != 0) {
        copied.property_count = i;
        chem_species_destroy(&copied);
        return -1;
      }
    }
  }
  chem_species_destroy(destination);
  *destination = copied;
  return 0;
}

int chem_species_set_property(chem_species *species, const char *name,
                              chem_species_property_kind kind,
                              chem_species_property_value value)
{
  char *property_name, *string_value = NULL;
  size_t i;
  if (!species || !name || !property_kind_valid(kind) ||
      (kind == CHEM_PROPERTY_STRING && !value.string_value))
    return -1;
  property_name = duplicate(name);
  if (!property_name)
    return -1;
  if (kind == CHEM_PROPERTY_STRING) {
    string_value = duplicate(value.string_value);
    if (!string_value) {
      free(property_name);
      return -1;
    }
  }
  for (i = 0; i < species->property_count; ++i)
    if (species->properties[i].name &&
        strcmp(species->properties[i].name, name) == 0)
      break;
  if (ensure_capacity(species, i == species->property_count
                                   ? species->property_count + 1
                                   : species->property_count) != 0) {
    free(property_name);
    free(string_value);
    return -1;
  }
  if (i == species->property_count)
    ++species->property_count;
  else
    property_destroy(&species->properties[i]);
  species->properties[i].name = property_name;
  species->properties[i].kind = kind;
  species->properties[i].value = value;
  if (kind == CHEM_PROPERTY_STRING)
    species->properties[i].value.string_value = string_value;
  return 0;
}

int chem_species_get_property(const chem_species *species, const char *name,
                              chem_species_property_kind kind, void *value)
{
  size_t i;
  if (!species || !name || !value || !property_kind_valid(kind))
    return -1;
  for (i = 0; i < species->property_count; ++i) {
    const chem_species_property *property = &species->properties[i];
    if (!property->name || strcmp(property->name, name) != 0)
      continue;
    if (property->kind != kind)
      return -3;
    switch (kind) {
    case CHEM_PROPERTY_DOUBLE: *(double *)value = property->value.real_value; break;
    case CHEM_PROPERTY_STRING: *(const char **)value = property->value.string_value; break;
    case CHEM_PROPERTY_INTEGER: *(int *)value = property->value.integer_value; break;
    case CHEM_PROPERTY_BOOLEAN: *(int *)value = property->value.boolean_value; break;
    }
    return 0;
  }
  return -2;
}

int chem_species_set_parameterization(chem_species *species,
                                      chem_species_parameter_kind kind,
                                      double c0, double c_temperature,
                                      double c_air_density)
{
  if (!species || kind < CHEM_SPECIES_UNPARAMETERIZED ||
      kind > CHEM_SPECIES_AIR_DENSITY)
    return -1;
  species->parameter_kind = kind;
  species->parameter_c0 = c0;
  species->parameter_c_temperature = c_temperature;
  species->parameter_c_air_density = c_air_density;
  return 0;
}

int chem_species_set_third_body(chem_species *species)
{
  return chem_species_set_parameterization(species, CHEM_SPECIES_AIR_DENSITY,
                                           0.0, 0.0, 1.0);
}

int chem_species_is_parameterized(const chem_species *species)
{
  return species && species->parameter_kind != CHEM_SPECIES_UNPARAMETERIZED;
}

double chem_species_parameterize(const chem_species *species,
                                 double temperature, double air_density)
{
  if (!species)
    return 0.0;
  switch (species->parameter_kind) {
  case CHEM_SPECIES_CONSTANT: return species->parameter_c0;
  case CHEM_SPECIES_TEMPERATURE: return species->parameter_c_temperature * temperature;
  case CHEM_SPECIES_AIR_DENSITY: return species->parameter_c_air_density * air_density;
  case CHEM_SPECIES_UNPARAMETERIZED: break;
  }
  return 0.0;
}
