#include "chem.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
  return 1; \
} } while (0)

static int species_string_constructor(void) {
  chem_species_info species = { 0 };
  CHECK(chem_species_info_init(&species, "thing") == 0);
  CHECK(strcmp(species.name, "thing") == 0);
  CHECK(!chem_species_info_is_parameterized(&species));
  CHECK(chem_species_info_set_parameterization(&species, (chem_species_parameterization){
      .temperature = 100.0, .enabled = 1 }) == 0);
  CHECK(chem_species_info_is_parameterized(&species));
  CHECK(chem_species_info_parameterized_value(&species, 12.5, 0.0, 0.0) == 1250.0);
  chem_species_info_destroy(&species);
  return 0;
}

static int species_properties(void) {
  chem_species_info species = { 0 };
  const char *string_value = NULL;
  double double_value = 0.0;
  int integer_value = 0, boolean_value = 0;
  CHECK(chem_species_info_init(&species, "thing") == 0);
  CHECK(chem_species_info_set_double(&species, "name [units]", 1.0) == 0);
  CHECK(chem_species_info_set_double(&species, "name2 [units2]", 2.0) == 0);
  CHECK(species.property_count == 2u);
  CHECK(chem_species_info_get_property(&species, "name [units]", CHEM_PROPERTY_DOUBLE, &double_value) == 0);
  CHECK(double_value == 1.0);
  CHECK(chem_species_info_get_property(&species, "name2 [units2]", CHEM_PROPERTY_DOUBLE, &double_value) == 0);
  CHECK(double_value == 2.0);
  CHECK(chem_species_info_set_string(&species, "foo", "bar") == 0);
  CHECK(chem_species_info_set_integer(&species, "baz", 42) == 0);
  CHECK(chem_species_info_set_boolean(&species, "qux", 1) == 0);
  CHECK(chem_species_info_get_property(&species, "foo", CHEM_PROPERTY_STRING, &string_value) == 0);
  CHECK(strcmp(string_value, "bar") == 0);
  CHECK(chem_species_info_get_property(&species, "baz", CHEM_PROPERTY_INTEGER, &integer_value) == 0);
  CHECK(integer_value == 42);
  CHECK(chem_species_info_get_property(&species, "qux", CHEM_PROPERTY_BOOLEAN, &boolean_value) == 0);
  CHECK(boolean_value == 1);
  CHECK(chem_species_info_get_property(&species, "not there", CHEM_PROPERTY_STRING, &string_value) == -1);
  CHECK(chem_species_info_get_property(&species, "not there", CHEM_PROPERTY_DOUBLE, &double_value) == -1);
  CHECK(chem_species_info_get_property(&species, "not there", CHEM_PROPERTY_INTEGER, &integer_value) == -1);
  CHECK(chem_species_info_get_property(&species, "not there", CHEM_PROPERTY_BOOLEAN, &boolean_value) == -1);
  CHECK(chem_species_info_get_property(&species, "name [units]", (chem_property_type)99, &double_value) == -2);
  chem_species_info_destroy(&species);
  return 0;
}

static int species_third_body(void) {
  chem_species_info species = { 0 };
  CHECK(chem_species_info_init(&species, "M") == 0);
  chem_species_info_set_third_body(&species);
  CHECK(strcmp(species.name, "M") == 0);
  CHECK(chem_species_info_is_parameterized(&species));
  CHECK(chem_species_info_parameterized_value(&species, 250.0, 100000.0, 42.4) == 42.4);
  chem_species_info_destroy(&species);
  return 0;
}

static int species_copy(void) {
  for (int assignment = 0; assignment < 2; ++assignment) {
    chem_species_info original = { 0 }, copy = { 0 };
    const char *str = NULL;
    int integer = 0, boolean = 0;
    double real = 0.0;
    CHECK(chem_species_info_init(&original, "thing") == 0);
    CHECK(chem_species_info_set_double(&original, "name [units]", 1.0) == 0);
    CHECK(chem_species_info_set_double(&original, "name2 [units2]", 2.0) == 0);
    CHECK(chem_species_info_set_string(&original, "foo", "bar") == 0);
    CHECK(chem_species_info_set_integer(&original, "baz", 42) == 0);
    CHECK(chem_species_info_set_boolean(&original, "qux", 1) == 0);
    if (assignment != 0)
      CHECK(chem_species_info_set_parameterization(&original, (chem_species_parameterization){
          .constant = 15.4, .enabled = 1 }) == 0);
    CHECK(chem_species_info_copy(&copy, &original) == 0);
    CHECK(strcmp(copy.name, "thing") == 0 && copy.name != original.name);
    CHECK(copy.property_count == 5u && copy.properties != original.properties);
    CHECK(chem_species_info_get_property(&copy, "name [units]", CHEM_PROPERTY_DOUBLE, &real) == 0 && real == 1.0);
    CHECK(chem_species_info_get_property(&copy, "name2 [units2]", CHEM_PROPERTY_DOUBLE, &real) == 0 && real == 2.0);
    CHECK(chem_species_info_get_property(&copy, "foo", CHEM_PROPERTY_STRING, &str) == 0 && strcmp(str, "bar") == 0);
    CHECK(chem_species_info_get_property(&copy, "baz", CHEM_PROPERTY_INTEGER, &integer) == 0 && integer == 42);
    CHECK(chem_species_info_get_property(&copy, "qux", CHEM_PROPERTY_BOOLEAN, &boolean) == 0 && boolean == 1);
    CHECK(chem_species_info_is_parameterized(&copy) == (assignment != 0));
    if (assignment != 0) CHECK(chem_species_info_parameterized_value(&copy, 0.0, 0.0, 0.0) == 15.4);
    chem_species_info_destroy(&original);
    CHECK(strcmp(copy.name, "thing") == 0);
    chem_species_info_destroy(&copy);
  }
  return 0;
}

static int phase_species_no_diffusion(void) {
  chem_species_info species = { 0 };
  chem_phase_species phase_species = { 0 };
  CHECK(chem_species_info_init(&species, "CO2") == 0);
  CHECK(chem_phase_species_init(&phase_species, &species, NULL, NULL) == 0);
  CHECK(strcmp(phase_species.species.name, "CO2") == 0);
  CHECK(!phase_species.has_diffusion_coefficient);
  chem_phase_species_destroy(&phase_species);
  chem_species_info_destroy(&species);
  return 0;
}

static int phase_species_with_diffusion(void) {
  chem_species_info species = { 0 };
  chem_phase_species phase_species = { 0 };
  const double diffusion = 1.23e-5;
  CHECK(chem_species_info_init(&species, "CO2") == 0);
  CHECK(chem_phase_species_init(&phase_species, &species, &diffusion, NULL) == 0);
  CHECK(strcmp(phase_species.species.name, "CO2") == 0);
  CHECK(phase_species.has_diffusion_coefficient && phase_species.diffusion_coefficient == diffusion);
  chem_phase_species_destroy(&phase_species);
  chem_species_info_destroy(&species);
  return 0;
}

static int phase_species_set_diffusion(void) {
  chem_species_info species = { 0 };
  chem_phase_species phase_species = { 0 };
  const double diffusion = 2.5e-6;
  CHECK(chem_species_info_init(&species, "CO2") == 0);
  CHECK(chem_phase_species_init(&phase_species, &species, NULL, NULL) == 0);
  CHECK(!phase_species.has_diffusion_coefficient);
  CHECK(chem_phase_species_set_diffusion_coefficient(&phase_species, diffusion) == 0);
  CHECK(phase_species.has_diffusion_coefficient && phase_species.diffusion_coefficient == diffusion);
  chem_phase_species_destroy(&phase_species);
  chem_species_info_destroy(&species);
  return 0;
}

static int init_phase_named(chem_phase *phase, const char *const *names, size_t count, size_t parameterized_index) {
  chem_species_info species[3] = { { 0 } };
  chem_phase_species members[3] = { 0 };
  int result = 1;
  if (count > 3u) return 1;
  for (size_t i = 0; i < count; ++i) {
    if (chem_species_info_init(&species[i], names[i]) != 0) goto done;
    if (i == parameterized_index && chem_species_info_set_parameterization(
            &species[i], (chem_species_parameterization){ .constant = 42.0, .enabled = 1 }) != 0) goto done;
    if (chem_phase_species_init(&members[i], &species[i], NULL, NULL) != 0) goto done;
  }
  if (chem_phase_init(phase, "gas", members, count) != 0) goto done;
  result = 0;
done:
  for (size_t i = 0; i < count; ++i) {
    chem_phase_species_destroy(&members[i]);
    chem_species_info_destroy(&species[i]);
  }
  return result;
}

static int init_basic_phase(chem_phase *phase, int parameterized_middle) {
  static const char *const names[] = { "foo", "bar", "baz" };
  return init_phase_named(phase, names, 3u, parameterized_middle ? 1u : SIZE_MAX);
}

static int phase_constructor(void) {
  static const char *const member_names[] = { "species1", "species2" };
  chem_phase phase = { 0 };
  chem_string_list names = { 0 };
  CHECK(init_phase_named(&phase, member_names, 2u, SIZE_MAX) == 0);
  CHECK(strcmp(phase.name, "gas") == 0 && phase.species_count == 2u);
  CHECK(chem_phase_state_size(&phase) == 2u);
  CHECK(chem_phase_unique_names(&phase, &names) == 0);
  CHECK(names.count == 2u);
  CHECK(strcmp(names.items[0], "gas.species1") == 0);
  CHECK(strcmp(names.items[1], "gas.species2") == 0);
  chem_string_list_destroy(&names);
  chem_phase_destroy(&phase);
  return 0;
}

static int phase_parameterized_species(void) {
  chem_phase phase = { 0 };
  CHECK(init_basic_phase(&phase, 1) == 0);
  CHECK(phase.species_count == 3u);
  CHECK(chem_phase_state_size(&phase) == 2u);
  chem_phase_destroy(&phase);
  return 0;
}

static int phase_unique_parameterized_names(void) {
  chem_phase phase = { 0 };
  chem_string_list names = { 0 };
  CHECK(init_basic_phase(&phase, 1) == 0);
  CHECK(chem_phase_unique_names(&phase, &names) == 0);
  CHECK(names.count == 2u);
  CHECK(strcmp(names.items[0], "gas.foo") == 0);
  CHECK(strcmp(names.items[1], "gas.baz") == 0);
  chem_string_list_destroy(&names);
  chem_phase_destroy(&phase);
  return 0;
}

static size_t reverse_names(const char *const *names, size_t count, size_t output_index, void *context) {
  (void)names;
  (void)context;
  return count - output_index - 1u;
}

static int system_constructor(void) {
  static const char *const species_names[] = { "foo", "bar" };
  chem_phase phase = { 0 };
  chem_system system = { 0 };
  chem_string_list names = { 0 };
  CHECK(init_phase_named(&phase, species_names, 2u, SIZE_MAX) == 0);
  CHECK(chem_system_init(&system, &phase) == 0);
  CHECK(system.gas_phase.species_count == 2u);
  CHECK(chem_system_state_size(&system) == 2u);
  CHECK(chem_system_unique_names(&system, NULL, NULL, &names) == 0);
  CHECK(names.count == 2u);
  CHECK(strcmp(names.items[0], "foo") == 0);
  CHECK(strcmp(names.items[1], "bar") == 0);
  chem_string_list_destroy(&names);
  chem_system_destroy(&system);
  chem_phase_destroy(&phase);
  return 0;
}

static int system_parameterized_species(void) {
  static const char *const species_names[] = { "foo", "bar", "paramSpecies" };
  chem_phase phase = { 0 };
  chem_system system = { 0 };
  chem_string_list names = { 0 };
  CHECK(init_phase_named(&phase, species_names, 3u, 2u) == 0);
  CHECK(chem_system_init(&system, &phase) == 0);
  CHECK(system.gas_phase.species_count == 3u);
  CHECK(chem_system_state_size(&system) == 2u);
  CHECK(chem_system_unique_names(&system, NULL, NULL, &names) == 0);
  CHECK(names.count == 2u && strcmp(names.items[0], "foo") == 0 && strcmp(names.items[1], "bar") == 0);
  chem_string_list_destroy(&names);
  CHECK(chem_system_unique_names(&system, reverse_names, NULL, &names) == 0);
  CHECK(names.count == 2u && strcmp(names.items[0], "bar") == 0 && strcmp(names.items[1], "foo") == 0);
  chem_string_list_destroy(&names);
  chem_system_destroy(&system);
  chem_phase_destroy(&phase);
  return 0;
}

typedef struct { const char *name; int (*run)(void); } test_case;

int main(int argc, char **argv) {
  static const test_case tests[] = {
    { "Species.StringConstructor", species_string_constructor },
    { "Species.StringAndVectorConstructor", species_properties },
    { "Species.GetProperty", species_properties },
    { "Species.SetThirdBody", species_third_body },
    { "Species.CopyConstructor", species_copy },
    { "Species.CopyAssignment", species_copy },
    { "PhaseSpecies.ConstructsWithoutDiffusionCoefficient", phase_species_no_diffusion },
    { "PhaseSpecies.ConstructsWithDiffusionCoefficient", phase_species_with_diffusion },
    { "PhaseSpecies.SetDiffusionCoefficient", phase_species_set_diffusion },
    { "Phase.Constructor", phase_constructor },
    { "Phase.ConstructorWithParameterizedSpecies", phase_parameterized_species },
    { "Phase.UniqueNamesWithParameterizedSpecies", phase_unique_parameterized_names },
    { "System.ConstructorWithAllParameters", system_constructor },
    { "System.ConstructorWithParameterizedSpecies", system_parameterized_species },
  };
  if (argc != 2) return 2;
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    if (strcmp(argv[1], tests[i].name) == 0) return tests[i].run();
  fprintf(stderr, "unknown object test: %s\n", argv[1]);
  return 2;
}
