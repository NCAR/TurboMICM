#include "chem.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static char *copy_string(const char *value) {
  if (value == NULL) return NULL;
  const size_t length = strlen(value);
  if (length == SIZE_MAX) return NULL;
  char *copy = (char *)malloc(length + 1u);
  if (copy != NULL) memcpy(copy, value, length + 1u);
  return copy;
}

static int reserve_properties(chem_species_info *species, size_t count) {
  if (count > SIZE_MAX / sizeof(*species->properties)) return -1;
  chem_species_property *properties = (chem_species_property *)realloc(
      species->properties, count * sizeof(*species->properties));
  if (properties == NULL && count != 0u) return -3;
  species->properties = properties;
  return 0;
}

static size_t property_index(
    const chem_species_info *species, const char *name, chem_property_type type) {
  for (size_t i = 0; i < species->property_count; ++i)
    if (species->properties[i].type == type && strcmp(species->properties[i].name, name) == 0) return i;
  return SIZE_MAX;
}

static int set_property(
    chem_species_info *species,
    const char *name,
    chem_property_type type,
    const void *value) {
  if (species == NULL || species->name == NULL || name == NULL || value == NULL ||
      type < CHEM_PROPERTY_STRING || type > CHEM_PROPERTY_BOOLEAN) return -1;
  char *name_copy = NULL;
  char *string_copy = NULL;
  if (type == CHEM_PROPERTY_STRING) {
    string_copy = copy_string(*(const char *const *)value);
    if (string_copy == NULL) return -3;
  }

  size_t index = property_index(species, name, type);
  if (index == SIZE_MAX) {
    if (species->property_count == SIZE_MAX) {
      free(string_copy);
      return -1;
    }
    name_copy = copy_string(name);
    if (name_copy == NULL) {
      free(string_copy);
      return -3;
    }
    const int status = reserve_properties(species, species->property_count + 1u);
    if (status != 0) {
      free(name_copy);
      free(string_copy);
      return status;
    }
    species->properties[species->property_count].name = name_copy;
    species->properties[species->property_count].type = type;
    index = species->property_count++;
  } else if (type == CHEM_PROPERTY_STRING) {
    free(species->properties[index].value.string_value);
  }

  switch (type) {
    case CHEM_PROPERTY_STRING: species->properties[index].value.string_value = string_copy; break;
    case CHEM_PROPERTY_DOUBLE: memcpy(&species->properties[index].value.double_value, value, sizeof(double)); break;
    case CHEM_PROPERTY_INTEGER: memcpy(&species->properties[index].value.integer_value, value, sizeof(int)); break;
    case CHEM_PROPERTY_BOOLEAN: {
      int boolean_value = 0;
      memcpy(&boolean_value, value, sizeof(int));
      species->properties[index].value.boolean_value = boolean_value != 0;
      break;
    }
    default: return -1;
  }
  return 0;
}

int chem_species_info_init(chem_species_info *species, const char *name) {
  if (species == NULL || name == NULL) return -1;
  char *name_copy = copy_string(name);
  if (name_copy == NULL) return -3;
  chem_species_info result = { 0 };
  result.name = name_copy;
  chem_species_info_destroy(species);
  *species = result;
  return 0;
}

void chem_species_info_destroy(chem_species_info *species) {
  if (species == NULL) return;
  free(species->name);
  for (size_t i = 0; i < species->property_count; ++i) {
    free(species->properties[i].name);
    if (species->properties[i].type == CHEM_PROPERTY_STRING)
      free(species->properties[i].value.string_value);
  }
  free(species->properties);
  memset(species, 0, sizeof(*species));
}

int chem_species_info_copy(chem_species_info *destination, const chem_species_info *source) {
  if (destination == NULL || source == NULL || destination == source || source->name == NULL) return -1;
  chem_species_info copy = { 0 };
  copy.name = copy_string(source->name);
  if (copy.name == NULL) return -3;
  copy.parameterization = source->parameterization;
  if (source->property_count != 0u) {
    copy.properties = (chem_species_property *)calloc(source->property_count, sizeof(*copy.properties));
    if (copy.properties == NULL) {
      chem_species_info_destroy(&copy);
      return -3;
    }
    for (size_t i = 0; i < source->property_count; ++i) {
      const chem_species_property *from = &source->properties[i];
      chem_species_property *to = &copy.properties[i];
      to->name = copy_string(from->name);
      if (to->name == NULL) {
        copy.property_count = i;
        chem_species_info_destroy(&copy);
        return -3;
      }
      to->type = from->type;
      if (from->type == CHEM_PROPERTY_STRING) {
        to->value.string_value = copy_string(from->value.string_value);
        if (to->value.string_value == NULL) {
          copy.property_count = i + 1u;
          chem_species_info_destroy(&copy);
          return -3;
        }
      } else {
        to->value = from->value;
      }
      copy.property_count = i + 1u;
    }
  }
  chem_species_info_destroy(destination);
  *destination = copy;
  return 0;
}

int chem_species_info_set_string(chem_species_info *species, const char *name, const char *value) {
  if (value == NULL) return -1;
  return set_property(species, name, CHEM_PROPERTY_STRING, &value);
}

int chem_species_info_set_double(chem_species_info *species, const char *name, double value) {
  return set_property(species, name, CHEM_PROPERTY_DOUBLE, &value);
}

int chem_species_info_set_integer(chem_species_info *species, const char *name, int value) {
  return set_property(species, name, CHEM_PROPERTY_INTEGER, &value);
}

int chem_species_info_set_boolean(chem_species_info *species, const char *name, int value) {
  return set_property(species, name, CHEM_PROPERTY_BOOLEAN, &value);
}

int chem_species_info_get_property(
    const chem_species_info *species, const char *name, chem_property_type type, void *value) {
  if (species == NULL || species->name == NULL || name == NULL || value == NULL) return -1;
  if (type < CHEM_PROPERTY_STRING || type > CHEM_PROPERTY_BOOLEAN) return -2;
  const size_t index = property_index(species, name, type);
  if (index == SIZE_MAX) return -1;
  const chem_species_property *property = &species->properties[index];
  switch (type) {
    case CHEM_PROPERTY_STRING: *(const char **)value = property->value.string_value; break;
    case CHEM_PROPERTY_DOUBLE: memcpy(value, &property->value.double_value, sizeof(double)); break;
    case CHEM_PROPERTY_INTEGER: memcpy(value, &property->value.integer_value, sizeof(int)); break;
    case CHEM_PROPERTY_BOOLEAN: memcpy(value, &property->value.boolean_value, sizeof(int)); break;
    default: return -2;
  }
  return 0;
}

int chem_species_info_set_parameterization(
    chem_species_info *species, chem_species_parameterization parameterization) {
  if (species == NULL || species->name == NULL) return -1;
  species->parameterization = parameterization;
  species->parameterization.enabled = parameterization.enabled != 0;
  return 0;
}

void chem_species_info_set_third_body(chem_species_info *species) {
  if (species == NULL || species->name == NULL) return;
  species->parameterization = (chem_species_parameterization){ .air_density = 1.0, .enabled = 1 };
}

int chem_species_info_is_parameterized(const chem_species_info *species) {
  return species != NULL && species->parameterization.enabled != 0;
}

double chem_species_info_parameterized_value(
    const chem_species_info *species, double temperature, double pressure, double air_density) {
  if (!chem_species_info_is_parameterized(species)) return 0.0;
  const chem_species_parameterization *p = &species->parameterization;
  return p->constant + p->temperature * temperature + p->pressure * pressure + p->air_density * air_density;
}

int chem_phase_species_init(
    chem_phase_species *phase_species,
    const chem_species_info *species,
    const double *diffusion_coefficient,
    const double *density) {
  if (phase_species == NULL || species == NULL || species->name == NULL) return -1;
  chem_phase_species result = { 0 };
  int status = chem_species_info_copy(&result.species, species);
  if (status != 0) return status;
  if (diffusion_coefficient != NULL) {
    result.diffusion_coefficient = *diffusion_coefficient;
    result.has_diffusion_coefficient = 1;
  }
  if (density != NULL) {
    result.density = *density;
    result.has_density = 1;
  }
  chem_phase_species_destroy(phase_species);
  *phase_species = result;
  return 0;
}

void chem_phase_species_destroy(chem_phase_species *phase_species) {
  if (phase_species == NULL) return;
  chem_species_info_destroy(&phase_species->species);
  memset(phase_species, 0, sizeof(*phase_species));
}

int chem_phase_species_set_diffusion_coefficient(chem_phase_species *phase_species, double value) {
  if (phase_species == NULL || phase_species->species.name == NULL) return -1;
  phase_species->diffusion_coefficient = value;
  phase_species->has_diffusion_coefficient = 1;
  return 0;
}

int chem_phase_species_set_density(chem_phase_species *phase_species, double value) {
  if (phase_species == NULL || phase_species->species.name == NULL) return -1;
  phase_species->density = value;
  phase_species->has_density = 1;
  return 0;
}

int chem_phase_init(chem_phase *phase, const char *name, const chem_phase_species *species, size_t species_count) {
  if (phase == NULL || name == NULL || (species_count != 0u && species == NULL) ||
      species_count > SIZE_MAX / sizeof(*species)) return -1;
  chem_phase result = { 0 };
  result.name = copy_string(name);
  if (result.name == NULL) return -3;
  if (species_count != 0u) {
    result.species = (chem_phase_species *)calloc(species_count, sizeof(*result.species));
    if (result.species == NULL) {
      chem_phase_destroy(&result);
      return -3;
    }
    for (size_t i = 0; i < species_count; ++i) {
      const chem_phase_species *from = &species[i];
      int status = chem_phase_species_init(
          &result.species[i], &from->species,
          from->has_diffusion_coefficient ? &from->diffusion_coefficient : NULL,
          from->has_density ? &from->density : NULL);
      if (status != 0) {
        result.species_count = i;
        chem_phase_destroy(&result);
        return status;
      }
      ++result.species_count;
    }
  }
  chem_phase_destroy(phase);
  *phase = result;
  return 0;
}

void chem_phase_destroy(chem_phase *phase) {
  if (phase == NULL) return;
  free(phase->name);
  for (size_t i = 0; i < phase->species_count; ++i) chem_phase_species_destroy(&phase->species[i]);
  free(phase->species);
  memset(phase, 0, sizeof(*phase));
}

size_t chem_phase_state_size(const chem_phase *phase) {
  if (phase == NULL) return 0u;
  size_t count = 0u;
  for (size_t i = 0; i < phase->species_count; ++i)
    if (!chem_species_info_is_parameterized(&phase->species[i].species)) ++count;
  return count;
}

void chem_string_list_destroy(chem_string_list *list) {
  if (list == NULL) return;
  for (size_t i = 0; i < list->count; ++i) free(list->items[i]);
  free(list->items);
  memset(list, 0, sizeof(*list));
}

static int make_name_list(
    const chem_phase *phase, int qualified, chem_name_reorder reorder, void *context, chem_string_list *names) {
  if (phase == NULL || phase->name == NULL || names == NULL) return -1;
  const size_t count = chem_phase_state_size(phase);
  if (count > SIZE_MAX / sizeof(const char *) || count > SIZE_MAX / sizeof(char *)) return -1;
  const char **raw = count == 0u ? NULL : (const char **)malloc(count * sizeof(*raw));
  if (count != 0u && raw == NULL) return -3;
  size_t position = 0u;
  for (size_t i = 0; i < phase->species_count; ++i)
    if (!chem_species_info_is_parameterized(&phase->species[i].species)) raw[position++] = phase->species[i].species.name;

  chem_string_list result = { 0 };
  if (count != 0u) {
    result.items = (char **)calloc(count, sizeof(*result.items));
    if (result.items == NULL) {
      free(raw);
      return -3;
    }
    for (size_t i = 0; i < count; ++i) {
      size_t source_index = reorder == NULL ? i : reorder(raw, count, i, context);
      if (source_index >= count) {
        result.count = i;
        chem_string_list_destroy(&result);
        free(raw);
        return -1;
      }
      const char *base = raw[source_index];
      if (qualified) {
        size_t phase_length = strlen(phase->name), species_length = strlen(base);
        if (phase_length > SIZE_MAX - 2u || species_length > SIZE_MAX - phase_length - 2u) {
          result.count = i;
          chem_string_list_destroy(&result);
          free(raw);
          return -1;
        }
        result.items[i] = (char *)malloc(phase_length + species_length + 2u);
        if (result.items[i] != NULL) {
          memcpy(result.items[i], phase->name, phase_length);
          result.items[i][phase_length] = '.';
          memcpy(result.items[i] + phase_length + 1u, base, species_length + 1u);
        }
      } else {
        result.items[i] = copy_string(base);
      }
      if (result.items[i] == NULL) {
        result.count = i;
        chem_string_list_destroy(&result);
        free(raw);
        return -3;
      }
      result.count = i + 1u;
    }
  }
  free(raw);
  chem_string_list_destroy(names);
  *names = result;
  return 0;
}

int chem_phase_unique_names(const chem_phase *phase, chem_string_list *names) {
  return make_name_list(phase, 1, NULL, NULL, names);
}

int chem_system_init(chem_system *system, const chem_phase *gas_phase) {
  if (system == NULL || gas_phase == NULL || gas_phase->name == NULL) return -1;
  chem_system result = { 0 };
  int status = chem_phase_init(&result.gas_phase, gas_phase->name, gas_phase->species, gas_phase->species_count);
  if (status != 0) return status;
  chem_system_destroy(system);
  *system = result;
  return 0;
}

void chem_system_destroy(chem_system *system) {
  if (system == NULL) return;
  chem_phase_destroy(&system->gas_phase);
  memset(system, 0, sizeof(*system));
}

size_t chem_system_state_size(const chem_system *system) {
  return system == NULL ? 0u : chem_phase_state_size(&system->gas_phase);
}

int chem_system_unique_names(
    const chem_system *system, chem_name_reorder reorder, void *context, chem_string_list *names) {
  if (system == NULL) return -1;
  return make_name_list(&system->gas_phase, 0, reorder, context, names);
}
