#include "chem.h"
#include "jacobian_verification.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #condition); \
  return 1; \
} } while (0)

#define CHEM_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CHEM_PI 3.14159265358979323846264338327950288

static int near_relative(double actual, double expected, double tolerance) {
  return fabs(actual - expected) <= tolerance * fmax(fabs(expected), DBL_MIN);
}

static int builds_chemical_reaction(void) {
  const chem_species species[] = {
    { "O3", 0.0, 0.0, 0.0, 1e-3, 0u }, { "NO", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "NO2", 0.0, 0.0, 0.0, 1e-3, 0u }, { "O2", 0.0, 0.0, 0.0, 1e-3, 0u }
  };
  const chem_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_term products[] = { { 2u, 1.0 }, { 3u, 1.0 } };
  const chem_reaction reaction = {
    .reactant_offset = 0u, .reactant_count = 2u,
    .product_offset = 0u, .product_count = 2u,
    .rate = { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 1.0, [3] = 300.0 } }
  };
  const chem_mechanism mechanism = {
    .name = "chemical reaction builder equivalent", .species_count = 4u, .reaction_count = 1u,
    .reactant_term_count = 2u, .product_term_count = 2u,
    .species = species, .reactions = &reaction, .reactant_terms = reactants, .product_terms = products,
    .temperature = 298.15, .pressure = 101325.0, .relative_tolerance = 1e-4
  };
  CHECK(strcmp(mechanism.name, "chemical reaction builder equivalent") == 0);
  CHECK(mechanism.reaction_count == 1u && mechanism.species_count == 4u);
  CHECK(reaction.reactant_count == 2u && reaction.product_count == 2u);
  CHECK(strcmp(species[reactants[0].species].name, "O3") == 0);
  CHECK(strcmp(species[reactants[1].species].name, "NO") == 0);
  CHECK(strcmp(species[products[0].species].name, "NO2") == 0);
  CHECK(strcmp(species[products[1].species].name, "O2") == 0);
  return 0;
}

static int chemical_reaction_copy(void) {
  const chem_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
  const chem_term products[] = { { 2u, 1.0 } };
  chem_reaction original = {
    .reactant_offset = 0u, .reactant_count = 2u,
    .product_offset = 0u, .product_count = 1u,
    .rate = { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 12.2, [2] = 300.0 } }
  };
  chem_reaction copy = original;
  CHECK(copy.reactant_offset == original.reactant_offset);
  CHECK(copy.reactant_count == original.reactant_count);
  CHECK(copy.product_offset == original.product_offset);
  CHECK(copy.product_count == original.product_count);
  CHECK(copy.rate.kind == original.rate.kind && copy.rate.p[0] == original.rate.p[0]);
  CHECK(copy.rate.p[2] == original.rate.p[2]);
  CHECK(reactants[copy.reactant_offset].species == 0u);
  CHECK(products[copy.product_offset].species == 2u);
  return 0;
}

static int process_rate_constants(void) {
  const chem_species species[] = {
    { "foo", 0.025, 2.3e2, 0.0, 1e-3, 0u }, { "bar", 0.0, 0.0, 0.0, 1e-3, CHEM_SPECIES_THIRD_BODY }
  };
  const chem_parameter parameters[] = { { "bar_user", 0.0 }, { "foo_surf.effective radius [m]", 0.0 },
    { "foo_surf.particle number concentration [# m-3]", 0.0 } };
  const chem_reaction reactions[] = {
    { .rate = { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 12.2 * 0.82, [2] = 300.0 }, .air_density_power = 1u } },
    { .rate = { .kind = CHEM_RATE_SURFACE, .p = { [0] = 1.0 }, .custom_parameter = 1u, .surface_species = 0u } },
    { .rate = { .kind = CHEM_RATE_USER, .p = { [1] = 0.82 }, .custom_parameter = 0u, .air_density_power = 1u } }
  };
  const chem_mechanism mechanism = {
    .name = "process rates", .species_count = 2u, .reaction_count = 3u, .parameter_count = 3u,
    .species = species, .reactions = reactions, .parameters = parameters,
    .temperature = 285.0, .pressure = 101100.0, .relative_tolerance = 1e-4
  };
  const double temperatures[] = { 280.0, 285.0, 290.0, 295.0, 300.0 };
  const double pressures[] = { 90000.0, 95000.0, 100000.0, 105000.0, 110000.0 };
  for (size_t cell = 0; cell < 5u; ++cell) {
    const double air_density = pressures[cell] / (CHEM_GAS_CONSTANT * temperatures[cell]);
    const double user = 0.001 * (double)(cell + 1u);
    const double radius = 1e-8 * (double)(cell + 1u);
    const double number = 1e5 * (double)(cell + 1u);
    const double cell_parameters[] = { user, radius, number };
    double rate = 0.0;
    CHECK(chem_evaluate_rate(&mechanism, 0u, cell_parameters, temperatures[cell], pressures[cell], &rate) == 0);
    CHECK(near_relative(rate, 12.2 * exp(300.0 / temperatures[cell]) * air_density * 0.82, 1e-12));
    CHECK(chem_evaluate_rate(&mechanism, 1u, cell_parameters, temperatures[cell], pressures[cell], &rate) == 0);
    const double speed = sqrt(8.0 * CHEM_GAS_CONSTANT * temperatures[cell] / (CHEM_PI * 0.025));
    const double expected_surface = 4.0 * number * CHEM_PI * radius * radius /
                                    (radius / 2.3e2 + 4.0 / speed);
    CHECK(near_relative(rate, expected_surface, 1e-10));
    CHECK(chem_evaluate_rate(&mechanism, 2u, cell_parameters, temperatures[cell], pressures[cell], &rate) == 0);
    CHECK(near_relative(rate, user * air_density * 0.82, 1e-12));
  }
  return 0;
}

static const chem_species process_set_species[] = {
  { "foo", 0.0, 0.0, 0.0, 1e-3, 0u }, { "bar", 0.0, 0.0, 0.0, 1e-3, 0u },
  { "baz", 0.0, 0.0, 0.0, 1e-3, 0u }, { "quz", 0.0, 0.0, 0.0, 1e-3, 0u },
  { "quuz", 0.0, 0.0, 0.0, 1e-3, 0u }
};
static const chem_term process_set_reactants[] = {
  { 0u, 1.0 }, { 2u, 1.0 }, { 1u, 1.0 }, { 3u, 1.0 }, { 2u, 1.0 }
};
static const chem_term process_set_products[] = {
  { 1u, 1.0 }, { 4u, 2.4 }, { 0u, 1.0 }, { 3u, 1.4 }, { 1u, 1.0 }, { 3u, 2.5 }
};
static const chem_reaction process_set_reactions[] = {
  { .reactant_offset = 0u, .reactant_count = 2u, .product_offset = 0u, .product_count = 2u },
  { .reactant_offset = 2u, .reactant_count = 1u, .product_offset = 2u, .product_count = 2u },
  { .reactant_offset = 3u, .reactant_count = 1u, .product_offset = 4u, .product_count = 0u },
  { .reactant_offset = 4u, .reactant_count = 1u, .product_offset = 4u, .product_count = 2u }
};
static const chem_mechanism process_set_mechanism = {
  .name = "process set", .species_count = 5u, .reaction_count = 4u,
  .reactant_term_count = 5u, .product_term_count = 6u,
  .species = process_set_species, .reactions = process_set_reactions,
  .reactant_terms = process_set_reactants, .product_terms = process_set_products,
  .temperature = 298.15, .pressure = 101325.0, .relative_tolerance = 1e-4
};

static int process_set_forcing(void) {
  const double y[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 1.1, 1.2, 1.3, 1.4, 1.5 };
  const double rates[] = { 10.0, 20.0 * 70.0 * 0.72, 30.0, 40.0 * 70.0 * 0.72,
    110.0, 120.0 * 80.0 * 0.72, 130.0, 140.0 * 80.0 * 0.72 };
  double forcing[10];
  for (size_t i = 0; i < 10u; ++i) forcing[i] = 1000.0;
  CHECK(chem_accumulate_forcing(&process_set_mechanism, 2u, y, rates, NULL, forcing) == 0);
  const double expected[] = {
    1000.0 - 10.0 * 0.1 * 0.3 + 20.0 * 70.0 * 0.72 * 0.2,
    1000.0 + 10.0 * 0.1 * 0.3 - 20.0 * 0.2 * 70.0 * 0.72 + 40.0 * 70.0 * 0.72 * 0.3,
    1000.0 - 10.0 * 0.1 * 0.3 - 40.0 * 70.0 * 0.72 * 0.3,
    1000.0 + 20.0 * 70.0 * 0.72 * 0.2 * 1.4 - 30.0 * 0.4 + 40.0 * 70.0 * 0.72 * 2.5 * 0.3,
    1000.0 + 10.0 * 0.1 * 0.3 * 2.4,
    1000.0 - 110.0 * 1.1 * 1.3 + 120.0 * 80.0 * 0.72 * 1.2,
    1000.0 + 110.0 * 1.1 * 1.3 - 120.0 * 1.2 * 80.0 * 0.72 + 140.0 * 80.0 * 0.72 * 1.3,
    1000.0 - 110.0 * 1.1 * 1.3 - 140.0 * 80.0 * 0.72 * 1.3,
    1000.0 + 120.0 * 80.0 * 0.72 * 1.2 * 1.4 - 130.0 * 1.4 + 140.0 * 80.0 * 0.72 * 2.5 * 1.3,
    1000.0 + 110.0 * 1.1 * 1.3 * 2.4
  };
  for (size_t i = 0; i < 10u; ++i) CHECK(near_relative(forcing[i], expected[i], 1e-12));
  return 0;
}

static const chem_species mask_species[] = {
  { "A", 0.0, 0.0, 0.0, 1e-3, 0u }, { "B", 0.0, 0.0, 0.0, 1e-3, 0u },
  { "C", 0.0, 0.0, 0.0, 1e-3, 0u }, { "D", 0.0, 0.0, 0.0, 1e-3, 0u }
};
static const chem_term mask_reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
static const chem_term mask_products[] = { { 2u, 1.0 }, { 3u, 1.0 } };
static const chem_reaction mask_reactions[] = {
  { .reactant_offset = 0u, .reactant_count = 2u, .product_offset = 0u, .product_count = 2u }
};
static const chem_mechanism mask_mechanism = {
  .name = "masked process set", .species_count = 4u, .reaction_count = 1u,
  .reactant_term_count = 2u, .product_term_count = 2u, .species = mask_species,
  .reactions = mask_reactions, .reactant_terms = mask_reactants, .product_terms = mask_products,
  .temperature = 298.15, .pressure = 101325.0, .relative_tolerance = 1e-4
};

static int process_set_algebraic_mask(void) {
  const unsigned char differential[] = { 1u, 0u, 1u, 0u };
  const double y[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
  const double rates[] = { 10.0, 20.0 };
  double forcing[8], negative_jacobian[32];
  for (size_t i = 0; i < 8u; ++i) forcing[i] = 1000.0;
  for (size_t i = 0; i < 32u; ++i) negative_jacobian[i] = 500.0;
  CHECK(chem_accumulate_forcing(&mask_mechanism, 2u, y, rates, differential, forcing) == 0);
  CHECK(forcing[0] == 980.0 && forcing[1] == 1000.0 && forcing[2] == 1020.0 && forcing[3] == 1000.0);
  CHECK(forcing[4] == 400.0 && forcing[5] == 1000.0 && forcing[6] == 1600.0 && forcing[7] == 1000.0);
  CHECK(chem_accumulate_negative_jacobian(
      &mask_mechanism, 2u, y, rates, differential, negative_jacobian) == 0);
  CHECK(negative_jacobian[0] == 520.0 && negative_jacobian[1] == 510.0);
  CHECK(negative_jacobian[4] == 500.0 && negative_jacobian[5] == 500.0);
  CHECK(negative_jacobian[8] == 480.0 && negative_jacobian[9] == 490.0);
  CHECK(negative_jacobian[12] == 500.0 && negative_jacobian[13] == 500.0);
  CHECK(negative_jacobian[16] == 620.0 && negative_jacobian[17] == 600.0);
  CHECK(negative_jacobian[24] == 380.0 && negative_jacobian[25] == 400.0);
  return 0;
}

typedef struct {
  const chem_mechanism *mechanism;
  const double *rate_constants;
} finite_difference_context;

static int forcing_callback(size_t block, const double *state, size_t species_count, double *forcing, void *opaque) {
  finite_difference_context *context = (finite_difference_context *)opaque;
  if (species_count != context->mechanism->species_count) return -1;
  const double *rates = context->rate_constants + block * context->mechanism->reaction_count;
  return chem_evaluate_forcing(context->mechanism, 1u, state, rates, forcing);
}

static int process_set_finite_difference_jacobian(void) {
  static const chem_species species[] = {
    { "A", 0.0, 0.0, 0.0, 1e-3, 0u }, { "B", 0.0, 0.0, 0.0, 1e-3, 0u },
    { "C", 0.0, 0.0, 0.0, 1e-3, 0u }
  };
  static const chem_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 } };
  static const chem_term products[] = { { 2u, 2.0 }, { 0u, 1.0 } };
  static const chem_reaction reactions[] = {
    { .reactant_offset = 0u, .reactant_count = 2u, .product_offset = 0u, .product_count = 1u },
    { .reactant_offset = 2u, .reactant_count = 1u, .product_offset = 1u, .product_count = 1u }
  };
  const chem_mechanism mechanism = {
    .name = "finite difference reaction set", .species_count = 3u, .reaction_count = 2u,
    .reactant_term_count = 3u, .product_term_count = 2u,
    .species = species, .reactions = reactions, .reactant_terms = reactants, .product_terms = products,
    .temperature = 298.15, .pressure = 101325.0, .relative_tolerance = 1e-4
  };
  const double state[] = { 1.0, 2.0, 3.0, 0.5, 4.0, 1.5 };
  const double rates[] = { 5.0, 10.0, 8.0, 12.0 };
  double negative_jacobian[18] = { 0.0 }, fd_jacobian[18];
  const unsigned char pattern[] = { 1u, 1u, 1u, 1u, 1u, 0u, 1u, 1u, 1u };
  chem_jacobian_comparison comparison;
  finite_difference_context context = { &mechanism, rates };
  CHECK(chem_accumulate_negative_jacobian(&mechanism, 2u, state, rates, NULL, negative_jacobian) == 0);
  CHECK(chem_finite_difference_jacobian(forcing_callback, &context, state, 2u, 3u, 1e-6, fd_jacobian) == 0);
  CHECK(chem_compare_jacobian_to_finite_difference(
      pattern, negative_jacobian, fd_jacobian, 2u, 3u, 1e-6, 1e-7, &comparison) == 0);
  CHECK(comparison.passed);
  CHECK(chem_check_jacobian_sparsity(pattern, fd_jacobian, 2u, 3u, 1e-8, &comparison) == 0);
  CHECK(comparison.passed);
  return 0;
}

static unsigned next_random(unsigned *state) {
  *state = *state * 1664525u + 1013904223u;
  return *state;
}

static int random_process_set_case(size_t cell_count, size_t reaction_count, size_t species_count, unsigned seed) {
  int result = 1;
  chem_species *species = (chem_species *)calloc(species_count, sizeof(*species));
  chem_reaction *reactions = (chem_reaction *)calloc(reaction_count, sizeof(*reactions));
  chem_term *reactants = (chem_term *)calloc(reaction_count * 3u, sizeof(*reactants));
  chem_term *products = (chem_term *)calloc(reaction_count * 10u, sizeof(*products));
  double *y = (double *)calloc(cell_count * species_count, sizeof(*y));
  double *rates = (double *)calloc(cell_count * reaction_count, sizeof(*rates));
  double *forcing = (double *)calloc(cell_count * species_count, sizeof(*forcing));
  if (!species || !reactions || !reactants || !products || !y || !rates || !forcing) goto done;
  char **names = (char **)calloc(species_count, sizeof(*names));
  if (names == NULL) goto done;
  for (size_t s = 0; s < species_count; ++s) {
    names[s] = (char *)malloc(24u);
    if (names[s] == NULL) {
      for (size_t i = 0; i < s; ++i) free(names[i]);
      free(names);
      goto done;
    }
    (void)snprintf(names[s], 24u, "species_%zu", s);
    species[s].name = names[s];
    species[s].absolute_tolerance = 1e-3;
  }
  size_t reactant_count = 0u, product_count = 0u;
  unsigned rng = seed;
  for (size_t r = 0; r < reaction_count; ++r) {
    size_t reactant_n = next_random(&rng) % 4u;
    size_t product_n = next_random(&rng) % 11u;
    reactions[r].reactant_offset = (unsigned)reactant_count;
    reactions[r].reactant_count = (unsigned)reactant_n;
    reactions[r].product_offset = (unsigned)product_count;
    reactions[r].product_count = (unsigned)product_n;
    for (size_t q = 0; q < reactant_n; ++q)
      reactants[reactant_count++] = (chem_term){ next_random(&rng) % (unsigned)species_count, 1.0 };
    for (size_t q = 0; q < product_n; ++q)
      products[product_count++] = (chem_term){ next_random(&rng) % (unsigned)species_count, 1.2 };
  }
  for (size_t i = 0; i < cell_count * species_count; ++i) y[i] = 0.01 + (next_random(&rng) % 10000u) / 1000.0;
  for (size_t i = 0; i < cell_count * reaction_count; ++i) rates[i] = 0.01 + (next_random(&rng) % 10000u) / 100.0;
  const chem_mechanism mechanism = {
    .name = "deterministic random process set", .species_count = species_count, .reaction_count = reaction_count,
    .reactant_term_count = reactant_count, .product_term_count = product_count,
    .species = species, .reactions = reactions, .reactant_terms = reactants, .product_terms = products,
    .temperature = 298.15, .pressure = 101325.0, .relative_tolerance = 1e-4
  };
  if (chem_evaluate_forcing(&mechanism, cell_count, y, rates, forcing) != 0) {
    for (size_t s = 0; s < species_count; ++s) free(names[s]);
    free(names);
    goto done;
  }
  for (size_t i = 0; i < cell_count * species_count; ++i) if (!isfinite(forcing[i])) {
    for (size_t s = 0; s < species_count; ++s) free(names[s]);
    free(names);
    goto done;
  }
  for (size_t s = 0; s < species_count; ++s) free(names[s]);
  free(names);
  result = 0;
done:
  free(species); free(reactions); free(reactants); free(products); free(y); free(rates); free(forcing);
  return result;
}

static int random_process_set(void) {
  CHECK(random_process_set_case(200u, 50u, 40u, 1u) == 0);
  CHECK(random_process_set_case(300u, 30u, 20u, 2u) == 0);
  CHECK(random_process_set_case(400u, 100u, 80u, 3u) == 0);
  return 0;
}

typedef struct { const char *name; int (*run)(void); } test_case;

int main(int argc, char **argv) {
  static const test_case tests[] = {
    { "Process.BuildsChemicalReaction", builds_chemical_reaction },
    { "Process.ChemicalReactionCopyAssignmentSucceeds", chemical_reaction_copy },
    { "Process.Matrix", process_rate_constants },
    { "Process.VectorMatrix", process_rate_constants },
    { "ProcessSet.Matrix", process_set_forcing },
    { "ProcessSet.VectorMatrix", process_set_forcing },
    { "ProcessSetAlgebraicVariables.CudaMatrix", process_set_algebraic_mask },
    { "ProcessSetFiniteDifferenceJacobian.Matrix", process_set_finite_difference_jacobian },
    { "RandomProcessSet.Matrix", random_process_set }
  };
  if (argc != 2) return 2;
  for (size_t i = 0; i < sizeof(tests) / sizeof(tests[0]); ++i)
    if (strcmp(argv[1], tests[i].name) == 0) return tests[i].run();
  fprintf(stderr, "unknown process test: %s\n", argv[1]);
  return 2;
}
