/* MICM test/unit/system/: test_species.cpp, test_phase.cpp, test_system.cpp
 * (the species, phase, and system executables).
 *
 * MICM's Species, PhaseSpecies, Phase, and System are tmc_species,
 * tmc_phase_species, tmc_phase, and tmc_system: value types whose copies are
 * plain assignments (species) or tmc_*_copy (phases and systems, which own an
 * array). A species parameterization (MICM's parameterize_) is set with
 * tmc_species_set_parameterization and evaluated with
 * tmc_species_parameterized_value. */
#include "harness.h"
#include "tmc_support.h"

static int count_properties(const tmc_species *s, tmc_property_type type)
{
  int n = 0;
  for (int i = 0; i < s->n_properties; ++i)
    if (s->properties[i].type == type)
      ++n;
  return n;
}

/* MICM: Species("thing", { { "name [units]", 1.0 }, { "name2 [units2]", 2.0 } }) */
static tmc_species thing_with_two_doubles(void)
{
  tmc_species species = tmc_species_create_with_property("thing", "name [units]", 1.0);
  tmc_species_set_double(&species, "name2 [units2]", 2.0);
  return species;
}

static void expect_double_property(const tmc_species *s, const char *key, double expected)
{
  double value = -1.0;
  tmc_error err;
  TMC_EXPECT_EQ(tmc_species_get_double(s, key, &value, &err), TMC_OK);
  TMC_EXPECT_EQ(value, expected);
}

/* ---- species ---- */

TMC_TEST(species_string_constructor)
{
  tmc_species species = tmc_species_create("thing");
  TMC_EXPECT_STREQ(species.name, "thing");
  TMC_EXPECT_FALSE(tmc_species_is_parameterized(&species));
  /* MICM: species.parameterize_ = { .c_T_ = 100.0, .has_value_ = true } */
  tmc_species_set_parameterization(&species, 0.0, 100.0, 0.0, 0.0);
  TMC_EXPECT_TRUE(tmc_species_is_parameterized(&species));
  tmc_conditions conditions = { .temperature = 12.5 };
  TMC_EXPECT_EQ(tmc_species_parameterized_value(&species, &conditions), 12.5 * 100.0);
  return TMC_TEST_PASS;
}

TMC_TEST(species_string_and_vector_constructor)
{
  tmc_species species = thing_with_two_doubles();
  TMC_EXPECT_STREQ(species.name, "thing");
  TMC_EXPECT_EQ(count_properties(&species, TMC_PROPERTY_DOUBLE), 2);
  expect_double_property(&species, "name [units]", 1.0);
  expect_double_property(&species, "name2 [units2]", 2.0);
  return TMC_TEST_PASS;
}

TMC_TEST(species_get_property)
{
  tmc_species species = thing_with_two_doubles();
  expect_double_property(&species, "name [units]", 1.0);
  expect_double_property(&species, "name2 [units2]", 2.0);

  /* a missing property, requested as each of the four property types */
  tmc_error err;
  const char *string_value = NULL;
  double double_value = 0.0;
  int int_value = 0;
  bool bool_value = false;
  TMC_EXPECT_EQ(tmc_species_get_string(&species, "not there", &string_value, &err), TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_PROPERTY_NOT_FOUND);
  TMC_EXPECT_EQ(tmc_species_get_double(&species, "not there", &double_value, &err), TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_PROPERTY_NOT_FOUND);
  TMC_EXPECT_EQ(tmc_species_get_int(&species, "not there", &int_value, &err), TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_PROPERTY_NOT_FOUND);
  TMC_EXPECT_EQ(tmc_species_get_bool(&species, "not there", &bool_value, &err), TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_PROPERTY_NOT_FOUND);

  /* MICM requests "name [units]" as long double, a type no property can hold,
   * and gets INVALID_TYPE_FOR_PROPERTY. The C getters are typed, so the same
   * request is the double property asked for as each of the other types. */
  TMC_EXPECT_EQ(tmc_species_get_int(&species, "name [units]", &int_value, &err), TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_INVALID_TYPE_FOR_PROPERTY);
  TMC_EXPECT_EQ(tmc_species_get_string(&species, "name [units]", &string_value, &err), TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_INVALID_TYPE_FOR_PROPERTY);
  TMC_EXPECT_EQ(tmc_species_get_bool(&species, "name [units]", &bool_value, &err), TMC_ERROR_SPECIES);
  TMC_EXPECT_EQ(err.code, TMC_SPECIES_INVALID_TYPE_FOR_PROPERTY);
  return TMC_TEST_PASS;
}

TMC_TEST(species_set_third_body)
{
  tmc_species species = tmc_species_create("M");
  tmc_species_set_third_body(&species);
  TMC_EXPECT_STREQ(species.name, "M");
  TMC_EXPECT_TRUE(tmc_species_is_parameterized(&species));
  tmc_conditions conditions = { .air_density = 42.4 };
  TMC_EXPECT_NEAR(tmc_species_parameterized_value(&species, &conditions), 42.4, 0.0);
  return TMC_TEST_PASS;
}

static tmc_species thing_with_all_property_types(bool parameterized)
{
  tmc_species species = thing_with_two_doubles();
  if (parameterized) /* MICM: species.parameterize_ = { .c0_ = 15.4, .has_value_ = true } */
    tmc_species_set_parameterization(&species, 15.4, 0.0, 0.0, 0.0);
  tmc_species_set_string(&species, "foo", "bar");
  tmc_species_set_int(&species, "baz", 42);
  tmc_species_set_bool(&species, "qux", true);
  return species;
}

static void expect_copied_thing(const tmc_species *species2, bool parameterized)
{
  tmc_error err;
  TMC_EXPECT_STREQ(species2->name, "thing");
  TMC_EXPECT_EQ(count_properties(species2, TMC_PROPERTY_DOUBLE), 2);
  expect_double_property(species2, "name [units]", 1.0);
  expect_double_property(species2, "name2 [units2]", 2.0);
  const char *s = NULL;
  int i = 0;
  bool b = false;
  TMC_EXPECT_EQ(tmc_species_get_string(species2, "foo", &s, &err), TMC_OK);
  TMC_EXPECT_TRUE(s != NULL && strcmp(s, "bar") == 0);
  TMC_EXPECT_EQ(tmc_species_get_int(species2, "baz", &i, &err), TMC_OK);
  TMC_EXPECT_EQ(i, 42);
  TMC_EXPECT_EQ(tmc_species_get_bool(species2, "qux", &b, &err), TMC_OK);
  TMC_EXPECT_EQ(b, true);
  TMC_EXPECT_EQ(tmc_species_is_parameterized(species2), parameterized);
}

/* the copy is independent of the original */
static void expect_independent(tmc_species *original, const tmc_species *copy)
{
  tmc_species_set_double(original, "name [units]", -1.0);
  tmc_species_set_string(original, "foo", "changed");
  original->name[0] = 'X';
  expect_double_property(copy, "name [units]", 1.0);
  const char *s = NULL;
  tmc_error err;
  TMC_EXPECT_EQ(tmc_species_get_string(copy, "foo", &s, &err), TMC_OK);
  TMC_EXPECT_TRUE(s != NULL && strcmp(s, "bar") == 0);
  TMC_EXPECT_STREQ(copy->name, "thing");
}

TMC_TEST(species_copy_constructor)
{
  {
    tmc_species species = thing_with_all_property_types(false);
    tmc_species species2 = species;
    expect_copied_thing(&species2, false);
    expect_independent(&species, &species2);
  }
  {
    tmc_species species = thing_with_all_property_types(true);
    tmc_species species2 = species;
    expect_copied_thing(&species2, true);
    tmc_conditions conditions = { 0 };
    TMC_EXPECT_NEAR(tmc_species_parameterized_value(&species, &conditions), 15.4, 0.0);
    TMC_EXPECT_NEAR(tmc_species_parameterized_value(&species2, &conditions), 15.4, 0.0);
    expect_independent(&species, &species2);
  }
  return TMC_TEST_PASS;
}

TMC_TEST(species_copy_assignment)
{
  {
    tmc_species species = thing_with_all_property_types(false);
    tmc_species species2 = tmc_species_create_with_property("other", "other [units]", 3.0);
    species2 = species;
    expect_copied_thing(&species2, false);
    TMC_EXPECT_FALSE(tmc_species_has_property(&species2, "other [units]"));
    expect_independent(&species, &species2);
  }
  {
    tmc_species species = thing_with_all_property_types(true);
    tmc_species species2 = tmc_species_create("other");
    species2 = species;
    expect_copied_thing(&species2, true);
    tmc_conditions conditions = { 0 };
    TMC_EXPECT_NEAR(tmc_species_parameterized_value(&species, &conditions), 15.4, 0.0);
    TMC_EXPECT_NEAR(tmc_species_parameterized_value(&species2, &conditions), 15.4, 0.0);
    expect_independent(&species, &species2);
  }
  return TMC_TEST_PASS;
}

/* ---- phase ---- */

TMC_TEST(phase_species_without_diffusion_coefficient)
{
  tmc_species CO2 = tmc_species_create("CO2");
  tmc_phase phase = tmc_phase_create("gas");
  tmc_phase_add_species(&phase, &CO2); /* MICM: PhaseSpecies gas_CO2(CO2) */
  TMC_ASSERT_TRUE(phase.n_species == 1);
  const tmc_phase_species *gas_CO2 = &phase.species[0];
  TMC_EXPECT_STREQ(gas_CO2->species.name, "CO2");
  TMC_EXPECT_FALSE(gas_CO2->has_diffusion_coefficient);
  tmc_phase_free(&phase);
  return TMC_TEST_PASS;
}

TMC_TEST(phase_species_with_diffusion_coefficient)
{
  tmc_species CO2 = tmc_species_create("CO2");
  double diff_coeff = 1.23e-5;
  tmc_phase phase = tmc_phase_create("gas");
  tmc_phase_add_species_with_diffusion(&phase, &CO2, diff_coeff); /* MICM: PhaseSpecies gas_CO2(CO2, diff_coeff) */
  TMC_ASSERT_TRUE(phase.n_species == 1);
  const tmc_phase_species *gas_CO2 = &phase.species[0];
  TMC_EXPECT_STREQ(gas_CO2->species.name, "CO2");
  TMC_ASSERT_TRUE(gas_CO2->has_diffusion_coefficient);
  TMC_EXPECT_DOUBLE_EQ(gas_CO2->diffusion_coefficient, diff_coeff);
  tmc_phase_free(&phase);
  return TMC_TEST_PASS;
}

TMC_TEST(phase_species_set_diffusion_coefficient)
{
  tmc_species CO2 = tmc_species_create("CO2");
  tmc_phase phase = tmc_phase_create("gas");
  tmc_phase_add_species(&phase, &CO2);
  TMC_ASSERT_TRUE(phase.n_species == 1);
  tmc_phase_species *gas_CO2 = &phase.species[0];
  TMC_EXPECT_FALSE(gas_CO2->has_diffusion_coefficient);

  /* MICM's SetDiffusionCoefficient; the C phase species is a plain struct
   * whose optional coefficient is a flag and a value */
  double diff_coeff = 2.5e-6;
  gas_CO2->has_diffusion_coefficient = true;
  gas_CO2->diffusion_coefficient = diff_coeff;
  TMC_ASSERT_TRUE(gas_CO2->has_diffusion_coefficient);
  TMC_EXPECT_DOUBLE_EQ(gas_CO2->diffusion_coefficient, diff_coeff);
  /* the phase's lookup sees the change */
  const tmc_phase_species *found = tmc_phase_find(&phase, "CO2");
  TMC_ASSERT_TRUE(found != NULL);
  TMC_EXPECT_TRUE(found->has_diffusion_coefficient);
  TMC_EXPECT_DOUBLE_EQ(found->diffusion_coefficient, diff_coeff);
  tmc_phase_free(&phase);
  return TMC_TEST_PASS;
}

TMC_TEST(phase_constructor)
{
  tmc_species species1 = tmc_species_create("species1"), species2 = tmc_species_create("species2");
  tmc_phase phase = tmc_phase_create("test_phase");
  tmc_phase_add_species(&phase, &species1);
  tmc_phase_add_species(&phase, &species2);
  TMC_EXPECT_STREQ(phase.name, "test_phase");
  TMC_EXPECT_EQ(phase.n_species, 2);
  TMC_EXPECT_EQ(tmc_phase_state_size(&phase), 2);

  char names[4][TMC_NAME_MAX];
  TMC_EXPECT_EQ(tmc_phase_species_names(&phase, true, names, 4), 2);
  TMC_EXPECT_STREQ(names[0], "test_phase.species1");
  TMC_EXPECT_STREQ(names[1], "test_phase.species2");
  tmc_phase_free(&phase);
  return TMC_TEST_PASS;
}

/* foo, bar, baz with bar parameterized (MICM: bar.parameterize_ = { .c0_ = 42.0 }) */
static tmc_phase gas_with_parameterized_bar(void)
{
  tmc_species foo = tmc_species_create("foo"), bar = tmc_species_create("bar"), baz = tmc_species_create("baz");
  tmc_species_set_parameterization(&bar, 42.0, 0.0, 0.0, 0.0);
  tmc_phase phase = tmc_phase_create("gas");
  tmc_phase_add_species(&phase, &foo);
  tmc_phase_add_species(&phase, &bar);
  tmc_phase_add_species(&phase, &baz);
  return phase;
}

TMC_TEST(phase_constructor_with_parameterized_species)
{
  tmc_phase phase = gas_with_parameterized_bar();
  TMC_EXPECT_EQ(phase.n_species, 3);
  TMC_EXPECT_EQ(tmc_phase_state_size(&phase), 2);
  tmc_phase_free(&phase);
  return TMC_TEST_PASS;
}

TMC_TEST(phase_unique_names_with_parameterized_species)
{
  tmc_phase phase = gas_with_parameterized_bar();
  char names[4][TMC_NAME_MAX];
  int n = tmc_phase_species_names(&phase, true, names, 4);
  TMC_EXPECT_EQ(n, 2);
  if (n == 2)
  {
    TMC_EXPECT_STREQ(names[0], "gas.foo");
    TMC_EXPECT_STREQ(names[1], "gas.baz");
  }
  tmc_phase_free(&phase);
  return TMC_TEST_PASS;
}

/* ---- system ---- */

/* MICM compares the names as multisets */
static bool same_names(char (*names)[TMC_NAME_MAX], int n, const char *const *expected, int n_expected)
{
  if (n != n_expected)
    return false;
  bool used[16] = { false };
  for (int i = 0; i < n; ++i)
  {
    bool matched = false;
    for (int j = 0; j < n_expected && !matched; ++j)
      if (!used[j] && strcmp(names[i], expected[j]) == 0)
        used[j] = matched = true;
    if (!matched)
      return false;
  }
  return true;
}

TMC_TEST(system_constructor_with_all_parameters)
{
  tmc_species foo = tmc_species_create("foo"), bar = tmc_species_create("bar");
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &foo);
  tmc_phase_add_species(&gas_phase, &bar);
  tmc_system system = tmc_system_create(&gas_phase);

  TMC_EXPECT_EQ(system.gas_phase.n_species, 2);
  char names[4][TMC_NAME_MAX];
  int n = tmc_system_unique_names(&system, names, 4);
  const char *expected[] = { "foo", "bar" };
  TMC_EXPECT_EQ(n, 2);
  TMC_EXPECT_TRUE(same_names(names, n, expected, 2));

  tmc_system_free(&system);
  tmc_phase_free(&gas_phase);
  return TMC_TEST_PASS;
}

static const char *reorder_names(const char (*names)[TMC_NAME_MAX], int n_names, int i, void *context)
{
  const int *reorder = context;
  (void)n_names;
  return names[reorder[i]];
}

TMC_TEST(system_constructor_with_parameterized_species)
{
  tmc_species foo = tmc_species_create("foo"), bar = tmc_species_create("bar");
  tmc_species param_species = tmc_species_create("paramSpecies");
  tmc_species_set_parameterization(&param_species, 64.2, 0.0, 0.0, 0.0); /* MICM: { .c0_ = 64.2 } */
  tmc_phase gas_phase = tmc_phase_create("gas");
  tmc_phase_add_species(&gas_phase, &foo);
  tmc_phase_add_species(&gas_phase, &bar);
  tmc_phase_add_species(&gas_phase, &param_species);
  tmc_system system = tmc_system_create(&gas_phase);

  TMC_EXPECT_EQ(system.gas_phase.n_species, 3);
  TMC_EXPECT_EQ(tmc_system_state_size(&system), 3 - 1); /* one parameterized species */

  char names[4][TMC_NAME_MAX];
  int n = tmc_system_unique_names(&system, names, 4);
  const char *expected[] = { "foo", "bar" };
  TMC_EXPECT_EQ(n, 2);
  TMC_EXPECT_TRUE(same_names(names, n, expected, 2));

  /* MICM's UniqueNames(f) with f(variables, i) = variables[reorder[i]] */
  int reorder[2] = { 1, 0 };
  char reordered_names[4][TMC_NAME_MAX];
  int n_reordered = tmc_support_system_unique_names_ordered(&system, reorder_names, reorder, reordered_names, 4);
  TMC_EXPECT_EQ(n_reordered, 2);
  if (n == 2 && n_reordered == 2)
  {
    TMC_EXPECT_STREQ(reordered_names[0], names[1]);
    TMC_EXPECT_STREQ(reordered_names[1], names[0]);
  }

  tmc_system_free(&system);
  tmc_phase_free(&gas_phase);
  return TMC_TEST_PASS;
}

static const tmc_test_case cases[] = {
  TMC_CASE("species:Species.StringConstructor", species_string_constructor),
  TMC_CASE("species:Species.StringAndVectorConstructor", species_string_and_vector_constructor),
  TMC_CASE("species:Species.GetProperty", species_get_property),
  TMC_CASE("species:Species.SetThirdBody", species_set_third_body),
  TMC_CASE("species:Species.CopyConstructor", species_copy_constructor),
  TMC_CASE("species:Species.CopyAssignment", species_copy_assignment),
  TMC_CASE("phase:PhaseSpecies.ConstructsWithoutDiffusionCoefficient", phase_species_without_diffusion_coefficient),
  TMC_CASE("phase:PhaseSpecies.ConstructsWithDiffusionCoefficient", phase_species_with_diffusion_coefficient),
  TMC_CASE("phase:PhaseSpecies.SetDiffusionCoefficient", phase_species_set_diffusion_coefficient),
  TMC_CASE("phase:Phase.Constructor", phase_constructor),
  TMC_CASE("phase:Phase.ConstructorWithParameterizedSpecies", phase_constructor_with_parameterized_species),
  TMC_CASE("phase:Phase.UniqueNamesWithParameterizedSpecies", phase_unique_names_with_parameterized_species),
  TMC_CASE("system:System.ConstructorWithAllParameters", system_constructor_with_all_parameters),
  TMC_CASE("system:System.ConstructorWithParameterizedSpecies", system_constructor_with_parameterized_species),
  { NULL, NULL }
};

const tmc_test_case *tmc_tests_system(int *n)
{
  *n = (int)(sizeof cases / sizeof cases[0]) - 1;
  return cases;
}
