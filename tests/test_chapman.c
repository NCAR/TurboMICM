#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static const chem_species chapman_species[] = {
  { "M", 0.0, 0.0, 0.0, 1e-3, 0u }, { "Ar", 0.0, 0.0, 0.0, 1e-3, 0u },
  { "CO2", 0.0, 0.0, 0.0, 1e-3, 0u }, { "H2O", 0.0, 0.0, 0.0, 1e-3, 0u },
  { "N2", 0.0, 0.0, 0.0, 1e-3, 0u }, { "O1D", 0.0, 0.0, 0.0, 1e-3, 0u },
  { "O", 0.0, 0.0, 0.0, 1e-3, 0u }, { "O2", 0.0, 0.0, 0.0, 1e-3, 0u },
  { "O3", 0.0, 0.0, 0.0, 1e-3, 0u },
};
static const chem_parameter photolysis[] = {
  { "O2_1", 0.0 }, { "O3_1", 0.0 }, { "O3_2", 0.0 },
};
static const chem_term reactants[] = {
  { 4u, 1.0 }, { 5u, 1.0 }, { 5u, 1.0 }, { 7u, 1.0 }, { 6u, 1.0 }, { 8u, 1.0 },
  { 0u, 1.0 }, { 6u, 1.0 }, { 7u, 1.0 }, { 7u, 1.0 }, { 8u, 1.0 }, { 8u, 1.0 },
  { 6u, 1.0 }, { 8u, 1.0 }, { 0u, 1.0 }, { 6u, 1.0 }, { 7u, 1.0 },
};
static const chem_term products[] = {
  { 6u, 1.0 }, { 4u, 1.0 }, { 6u, 1.0 }, { 7u, 1.0 }, { 7u, 2.0 },
  { 6u, 2.0 }, { 5u, 1.0 }, { 7u, 1.0 }, { 6u, 1.0 }, { 7u, 1.0 },
  { 8u, 1.0 }, { 0u, 1.0 },
};
static const chem_reaction reactions[] = {
  { 0u, 2u, 0u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .label = "N2_O1D_1", .p = { [0] = 2.15e-11, [2] = -110.0, [3] = 300.0 } } },
  { 2u, 2u, 2u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .label = "O1D_O2_1", .p = { [0] = 3.3e-11, [2] = -55.0, [3] = 300.0 } } },
  { 4u, 2u, 4u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .label = "O_O3_1", .p = { [0] = 8.0e-12, [2] = 2060.0, [3] = 300.0 } } },
  { 6u, 3u, 10u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .label = "M_O_O2_1", .p = { [0] = 6.0e-34, [1] = 2.4, [3] = 300.0 } } },
  { 9u, 1u, 5u, 1u, { .kind = CHEM_RATE_USER, .label = "O2_1", .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 10u, 1u, 6u, 2u, { .kind = CHEM_RATE_USER, .label = "O3_1", .p = { [1] = 1.0 }, .custom_parameter = 1u } },
  { 11u, 1u, 8u, 2u, { .kind = CHEM_RATE_USER, .label = "O3_2", .p = { [1] = 1.0 }, .custom_parameter = 2u } },
};
static const chem_mechanism chapman = {
  .name = "Chapman mechanism test descriptor",
  .species_count = 9u,
  .reaction_count = 7u,
  .reactant_term_count = sizeof(reactants) / sizeof(reactants[0]),
  .product_term_count = sizeof(products) / sizeof(products[0]),
  .parameter_count = 3u,
  .species = chapman_species,
  .reactions = reactions,
  .reactant_terms = reactants,
  .product_terms = products,
  .parameters = photolysis,
  .temperature = 273.15,
  .pressure = 100000.0,
  .relative_tolerance = 1e-4,
};

static int default_constructor(void) {
  const chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  return solver.stage_count == 3u ? 0 : 1;
}

static int reaction_names(void) {
  static const char *const expected[] = {
    "O2_1", "O3_1", "O3_2", "N2_O1D_1", "O1D_O2_1", "O_O3_1", "M_O_O2_1"
  };
  for (size_t i = 0; i < sizeof(expected) / sizeof(expected[0]); ++i) {
    int found = 0;
    for (size_t r = 0; r < chapman.reaction_count; ++r)
      if (strcmp(expected[i], chapman.reactions[r].rate.label) == 0) found = 1;
    if (!found) return 1;
  }
  return chapman.reaction_count == 7u ? 0 : 1;
}

static int photolysis_names(void) {
  static const char *const expected[] = { "O2_1", "O3_1", "O3_2" };
  if (chapman.parameter_count != sizeof(expected) / sizeof(expected[0])) return 1;
  for (size_t i = 0; i < chapman.parameter_count; ++i)
    if (strcmp(photolysis[i].name, expected[i]) != 0) return 1;
  return 0;
}

static int species_names(void) {
  static const char *const expected[] = { "M", "Ar", "CO2", "H2O", "N2", "O1D", "O", "O2", "O3" };
  if (chapman.species_count != sizeof(expected) / sizeof(expected[0])) return 1;
  for (size_t i = 0; i < chapman.species_count; ++i)
    if (strcmp(chapman.species[i].name, expected[i]) != 0) return 1;
  return 0;
}

static int force_case(double concentration, double rate, const double *expected, double tolerance) {
  double y[9], rates[7], forcing[9];
  for (size_t i = 0; i < 9u; ++i) y[i] = concentration;
  for (size_t i = 0; i < 7u; ++i) rates[i] = rate;
  if (chem_evaluate_forcing(&chapman, 1u, y, rates, forcing) != 0) return 0;
  for (size_t i = 0; i < 9u; ++i) {
    if (!isfinite(forcing[i]) || fabs(forcing[i] - expected[i]) > tolerance) {
      fprintf(stderr, "Chapman forcing species %zu: %.17g expected %.17g\n", i, forcing[i], expected[i]);
      return 0;
    }
  }
  return 1;
}

static int simple_force(void) {
  const double expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 3.0, 2.0, -2.0 };
  return force_case(1.0, 1.0, expected, 0.0) ? 0 : 1;
}

static int smaller_force(void) {
  const double expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.49e-13, 4.55e-13, 1.5e-13, -3e-13 };
  return force_case(5e-6, 3e-8, expected, 0.01) ? 0 : 1;
}

static int integrate_case(const double *initial, const double *expected) {
  double y[9];
  double custom[] = { 1e-4, 1e-5, 1e-6 };
  const double temperature[] = { 273.15 }, pressure[] = { 100000.0 };
  memcpy(y, initial, sizeof(y));

  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.h_start_fraction = 0.005;
  solver.h_max_fraction = 0.5;
  solver.max_internal_steps = 1000u;
  if (chem_integrate_group_with_parameters(&chapman, 1u, y, custom, temperature,
          pressure, 1.0, &solver, NULL) != 0) return 0;
  for (size_t i = 0; i < 9u; ++i) {
    if (!isfinite(y[i]) || fabs(y[i] - expected[i]) > 1e-4) {
      fprintf(stderr, "Chapman solve species %zu: %.17g expected %.17g\n", i, y[i], expected[i]);
      return 0;
    }
  }
  return 1;
}

static int solve_one(void) {
  const double initial[] = { 1.0, 0.392, 0.0169, 0.0, 32.9, 0.0, 0.0, 8.84, 0.0 };
  const double expected[] = { 1.0, 0.392, 0.0169, 0.0, 32.9, 1.8039e-41, 0.00176789, 8.83912, 4.5031e-36 };
  return integrate_case(initial, expected) ? 0 : 1;
}

static int solve_large(void) {
  const double initial[] = { 10.0, 3.92, 0.169, 0.0, 329.0, 0.0, 0.0, 88.4, 0.0 };
  const double expected[] = { 10.0, 3.92, 0.169, 0.0, 329.0, 1.8039e-38, 0.0176789, 88.3912, 4.5031e-33 };
  return integrate_case(initial, expected) ? 0 : 1;
}

static int solve_small(void) {
  const double initial[] = { 0.1, 0.0392, 0.00169, 0.0, 3.29, 0.0, 0.0, 0.884, 0.0 };
  const double expected[] = { 0.1, 0.0392, 0.00169, 0.0, 3.29, 1.8039e-44, 0.000176789, 0.883912, 4.5031e-39 };
  return integrate_case(initial, expected) ? 0 : 1;
}

static int solve_random(void) {
  const double initial[] = { 0.23, 0.81, 0.14, 0.99, 0.52, 0.61, 0.37, 0.72, 0.44 };
  double y[9];
  double custom[] = { 1e-4, 1e-5, 1e-6 };
  const double temperature[] = { 273.15 }, pressure[] = { 100000.0 };
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.h_start_fraction = 0.005;
  solver.h_max_fraction = 0.5;
  solver.max_internal_steps = 1000u;
  memcpy(y, initial, sizeof(y));
  if (chem_integrate_group_with_parameters(&chapman, 1u, y, custom, temperature,
          pressure, 1.0, &solver, NULL) != 0) return 1;
  for (size_t i = 0; i < 9u; ++i)
    if (!isfinite(y[i]) || y[i] < 0.0) return 1;
  return 0;
}

typedef struct { const char *name; int (*run)(void); } test_case;

int main(int argc, char **argv) {
  static const test_case tests[] = {
    { "DefaultConstructor", default_constructor }, { "ReactionNames", reaction_names },
    { "PhotolysisNames", photolysis_names }, { "SpeciesNames", species_names },
    { "simple_force", simple_force }, { "smaller_force", smaller_force },
    { "Solve", solve_one }, { "solve_10_times_larger", solve_large },
    { "solve_10_times_smaller", solve_small }, { "solve_with_random_number_densities", solve_random },
  };
  if (argc != 2) return 2;
  const char *name = strrchr(argv[1], '.');
  name = name == NULL ? argv[1] : name + 1;
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    if (strcmp(name, tests[i].name) == 0) return tests[i].run();
  fprintf(stderr, "unimplemented Chapman case: %s\n", argv[1]);
  return 2;
}
