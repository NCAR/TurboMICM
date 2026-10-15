#include "test.h"
#include "chem_matrix.h"
static double parameter(const ck_conditions *c, const double *custom, void *context)
{
  (void)custom;
  return context ? *(double *)context : 100 * c->temperature;
}
int test_system(const char *suite, const char *test)
{
  if (!strcmp(suite, "Species"))
  {
    ck_property properties[] = {{.name = "name [units]", .type = CK_PROPERTY_REAL, .value.real = 1},
                                {.name = "name2 [units2]", .type = CK_PROPERTY_REAL, .value.real = 2},
                                {.name = "foo", .type = CK_PROPERTY_STRING, .value.string = "bar"},
                                {.name = "baz", .type = CK_PROPERTY_INT, .value.integer = 42},
                                {.name = "qux", .type = CK_PROPERTY_BOOL, .value.boolean = 1}};
    ck_species species = {.name = "thing", .nproperties = 2, .properties = properties};
    if (IS("StringConstructor"))
    {
      species.nproperties = 0;
      CHECK(!strcmp(species.name, "thing") && !species.parameterized);
      species.parameterized = 2;
      species.function = parameter;
      ck_conditions c = {.temperature = 12.5};
      CHECK(species.parameterized);
      NEAR(ck_species_parameter(&species, &c), 1250, 0);
      return 0;
    }
    if (IS("SetThirdBody"))
    {
      species.name = "M";
      species.parameterized = 1;
      ck_conditions c = {.air_density = 42.4};
      CHECK(!strcmp(species.name, "M") && species.parameterized);
      NEAR(ck_species_parameter(&species, &c), 42.4, 0);
      return 0;
    }
    if (!IS("StringAndVectorConstructor") && !IS("GetProperty") && !IS("CopyConstructor") &&
        !IS("CopyAssignment"))
      return 2;
    size_t runs = IS("CopyConstructor") || IS("CopyAssignment") ? 2 : 1;
    for (size_t run = 0; run < runs; ++run)
    {
      double constant = 15.4;
      if (runs == 2)
      {
        species.nproperties = 5;
        if (run)
        {
          species.parameterized = 2;
          species.function = parameter;
          species.context = &constant;
        }
      }
      ck_species copy = species;
      CHECK(!strcmp(copy.name, "thing"));
      CHECK(copy.nproperties == (runs == 2 ? 5u : 2u));
      double value;
      CHECK(ck_species_property(&copy, "name [units]", CK_PROPERTY_REAL, &value) == CK_OK);
      NEAR(value, 1, 0);
      CHECK(ck_species_property(&copy, "name2 [units2]", CK_PROPERTY_REAL, &value) == CK_OK);
      NEAR(value, 2, 0);
      if (runs == 2)
      {
        const char *text;
        int number;
        CHECK(ck_species_property(&copy, "foo", CK_PROPERTY_STRING, &text) == CK_OK && !strcmp(text, "bar"));
        CHECK(ck_species_property(&copy, "baz", CK_PROPERTY_INT, &number) == CK_OK && number == 42);
        CHECK(ck_species_property(&copy, "qux", CK_PROPERTY_BOOL, &number) == CK_OK && number == 1);
        CHECK(copy.parameterized == (run ? 2 : 0));
        if (run)
        {
          ck_conditions c = {0};
          NEAR(ck_species_parameter(&copy, &c), 15.4, 0);
        }
      }
      if (IS("GetProperty"))
      {
        for (int type = CK_PROPERTY_REAL; type <= CK_PROPERTY_BOOL; ++type)
          CHECK(ck_species_property(&species, "not there", (ck_property_type)type, &value) == CK_UNKNOWN);
        CHECK(ck_species_property(&species, "name [units]", (ck_property_type)99, &value) == CK_INVALID);
      }
    }
    return 0;
  }
  if (!strcmp(suite, "PhaseSpecies"))
  {
    ck_species species = {.name = "CO2"};
    CHECK(!strcmp(species.name, "CO2"));
    if (IS("ConstructsWithoutDiffusionCoefficient"))
    {
      NEAR(species.diffusion, 0, 0);
      return 0;
    }
    if (IS("ConstructsWithDiffusionCoefficient"))
    {
      species.diffusion = 1.23e-5;
      NEAR(species.diffusion, 1.23e-5, 0);
      return 0;
    }
    if (IS("SetDiffusionCoefficient"))
    {
      NEAR(species.diffusion, 0, 0);
      species.diffusion = 2.5e-6;
      NEAR(species.diffusion, 2.5e-6, 0);
      return 0;
    }
    return 2;
  }
  if (!strcmp(suite, "Phase") || !strcmp(suite, "System"))
  {
    int system = !strcmp(suite, "System"), param = strstr(test, "Parameterized") != NULL;
    if (!IS("Constructor") && !IS("ConstructorWithAllParameters") &&
        !IS("ConstructorWithParameterizedSpecies") && !IS("UniqueNamesWithParameterizedSpecies"))
      return 2;
    double val = system ? 64.2 : 42;
    ck_species species[] = {{.name = system || param ? "foo" : "species1"},
                            {.name = system || param ? "bar" : "species2"},
                            {.name = system ? "paramSpecies" : "baz"}};
    if (param)
    {
      size_t index = system ? 2 : 1;
      species[index].parameterized = 2;
      species[index].function = parameter;
      species[index].context = &val;
    }
    ck_phase phase = {
        .name = system || param ? "gas" : "test_phase", .nspecies = param ? 3 : 2, .species = species};
    CHECK(phase.nspecies == (param ? 3u : 2u));
    CHECK(ck_phase_state_size(&phase) == 2);
    char name[64];
    CHECK(ck_phase_unique_name(&phase, 0, name, sizeof name) == CK_OK);
    CHECK(!strcmp(name, system || param ? "gas.foo" : "test_phase.species1"));
    CHECK(ck_phase_unique_name(&phase, 1, name, sizeof name) == CK_OK);
    CHECK(!strcmp(name, system ? "gas.bar" : param ? "gas.baz" : "test_phase.species2"));
    ck_mechanism desc = {.name = phase.name, .nspecies = phase.nspecies, .species = species};
    ck_model *m;
    CHECK(ck_model_create_with_order(&desc, 0, &m) == CK_OK);
    CHECK(ck_model_nspecies(m) == 2);
    if (system)
    {
      CHECK(!strcmp(ck_species_name(m, 0), "foo") && !strcmp(ck_species_name(m, 1), "bar"));
      size_t reorder[] = {1, 0};
      CHECK(!strcmp(ck_species_name(m, reorder[0]), "bar"));
      CHECK(!strcmp(ck_species_name(m, reorder[1]), "foo"));
    }
    ck_model_free(m);
    return 0;
  }
  if (!strcmp(suite, "Version"))
  {
    if (IS("FullVersion"))
    {
      CHECK(ck_version() != NULL);
      return 0;
    }
    if (IS("VersionMajor"))
    {
      CHECK(ck_version_major() >= 0);
      return 0;
    }
    if (IS("VersionMinor"))
    {
      CHECK(ck_version_minor() >= 0);
      return 0;
    }
    if (IS("VersionPatch"))
    {
      CHECK(ck_version_patch() >= 0);
      return 0;
    }
    if (IS("VersionTweak"))
    {
      CHECK(ck_version_tweak() >= 0);
      return 0;
    }
    return 2;
  }
  if (!strcmp(suite, "ViewCategorySFINAE"))
  {
    if (!IS("RuntimeVectorUse"))
      return 3;
    double values[] = {1, 2, 3, 4, 5};
    ck_view view;
    CHECK(ck_vector_view(values, 5, sizeof(double), 0, &view) == CK_OK);
    double sum = 0;
    for (size_t j = 0; j < view.rows; ++j)
      sum += *(const double *)ck_view_const_at(&view, j);
    NEAR(sum, 15, 0);
    return 0;
  }
  return 2;
}
