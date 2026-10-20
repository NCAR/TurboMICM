#include "kinetics.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

extern const chem_mechanism chem_model;

#define MOLECULES_PER_MOL_M3 (1.0e-6 * 6.02214076e23)

#if EXAMPLE_USER_CONFIG
static const chem_species expected_species[] = {
  {"A", 0.0, 0.0, 0.0, 0.0}, {"B", 0.0, 0.0, 0.0, 0.0},
  {"D", 0.0, 0.0, 0.0, 0.0}, {"F", 0.0, 0.0, 0.0, 0.0},
  {"C", 0.025, 230.0, 0.0, 0.0}, {"E", 0.0, 0.0, 0.0, 0.0},
  {"G", 0.0, 0.0, 0.0, 0.0}
};
static const uint32_t expected_reactants[] = {0, 1, 1, 4, 2, 5, 5, 3, 4, 1};
static const uint32_t expected_products[] = {1, 4, 2, 5, 3, 6, 6, 6, 0};
static const double expected_yields[] = {1, 1, 1, 1, 2, 1, 1, 1, 1};
static const chem_reaction expected_reactions[] = {
  {.kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
   .product_offset=0, .product_count=1, .p={2.15e-4, 0, 110, 300}},
  {.kind=CHEM_BRANCHED, .reactant_offset=1, .reactant_count=1,
   .product_offset=1, .product_count=1, .p={204.3, 1.2, 1e-3, 2, 0}},
  {.kind=CHEM_BRANCHED, .reactant_offset=2, .reactant_count=1,
   .product_offset=2, .product_count=1, .p={204.3, 1.2, 1e-3, 2, 1}},
  {.kind=CHEM_SURFACE, .reactant_offset=3, .reactant_count=1,
   .product_offset=3, .product_count=1, .custom0=0, .custom1=1,
   .surface_species=4, .p={0.90}},
  {.kind=CHEM_TERNARY, .reactant_offset=4, .reactant_count=1,
   .product_offset=4, .product_count=1,
   .p={1.2, 2.3, 302.3, 2.6 / MOLECULES_PER_MOL_M3, -3.1, 402.1, 0.9, 1.2}},
  {.kind=CHEM_TROE, .reactant_offset=5, .reactant_count=2,
   .product_offset=5, .product_count=1,
   .p={1.8e-30 * MOLECULES_PER_MOL_M3 * MOLECULES_PER_MOL_M3,
       -3.0, 0.0, 2.8e-11 * MOLECULES_PER_MOL_M3, 0.0, 0.0, 0.6, 1.0}},
  {.kind=CHEM_TUNNELING, .reactant_offset=7, .reactant_count=1,
   .product_offset=6, .product_count=1, .p={1.2, 2.3, 302.3}},
  {.kind=CHEM_USER_DEFINED, .reactant_offset=8, .reactant_count=1,
   .product_offset=7, .product_count=1, .custom0=2, .p={1.0}},
  {.kind=CHEM_USER_DEFINED, .reactant_offset=9, .reactant_count=1,
   .product_offset=8, .product_count=0, .custom0=3, .p={1.0}},
  {.kind=CHEM_USER_DEFINED, .reactant_offset=10, .reactant_count=0,
   .product_offset=8, .product_count=1, .custom0=4, .p={1.0}}
};
static const chem_custom_parameter expected_custom[] = {
  {"C.effective radius [m]", 0.0},
  {"C.particle number concentration [# m-3]", 0.0},
  {"my photolysis rate", 0.0},
  {"my loss rate", 0.0},
  {"my emission rate", 0.0}
};
static const chem_mechanism expected_model = {
  .name="reaction rates no user defined", .species_count=7, .reaction_count=10,
  .reactant_count=10, .product_count=9, .custom_count=5,
  .default_temperature=298.15, .default_pressure=101325.0,
  .species=expected_species, .reactions=expected_reactions,
  .reactants=expected_reactants, .products=expected_products,
  .yields=expected_yields, .custom_parameters=expected_custom
};
#else
static const chem_species expected_species[] = {
  {"A", 0.0, 0.0, 0.0, 0.0}, {"B", 0.0, 0.0, 0.0, 0.0},
  {"C", 0.025, 230.0, 0.0, 0.0}, {"D", 0.0, 0.0, 0.0, 0.0},
  {"E", 0.0, 0.0, 0.0, 0.0}, {"F", 0.0, 0.0, 0.0, 0.0},
  {"G", 0.0, 0.0, 0.0, 0.0}
};
static const uint32_t expected_reactants[] = {0, 1, 1, 2, 3, 4, 4, 5};
static const uint32_t expected_products[] = {1, 2, 3, 4, 5, 6, 6};
static const double expected_yields[] = {1, 1, 1, 1, 2, 1, 1};
static const chem_reaction expected_reactions[] = {
  {.kind=CHEM_ARRHENIUS, .reactant_offset=0, .reactant_count=1,
   .product_offset=0, .product_count=1, .p={2.15e-4, 0, 110, 300}},
  {.kind=CHEM_BRANCHED, .reactant_offset=1, .reactant_count=1,
   .product_offset=1, .product_count=1, .p={204.3, 1.2, 1e-3, 2, 0}},
  {.kind=CHEM_BRANCHED, .reactant_offset=2, .reactant_count=1,
   .product_offset=2, .product_count=1, .p={204.3, 1.2, 1e-3, 2, 1}},
  {.kind=CHEM_SURFACE, .reactant_offset=3, .reactant_count=1,
   .product_offset=3, .product_count=1, .custom0=0, .custom1=1,
   .surface_species=2, .p={0.90}},
  {.kind=CHEM_TERNARY, .reactant_offset=4, .reactant_count=1,
   .product_offset=4, .product_count=1,
   .p={1.2, 2.3, 302.3, 2.6 / MOLECULES_PER_MOL_M3, -3.1, 402.1, 0.9, 1.2}},
  {.kind=CHEM_TROE, .reactant_offset=5, .reactant_count=2,
   .product_offset=5, .product_count=1,
   .p={1.8e-30 * MOLECULES_PER_MOL_M3 * MOLECULES_PER_MOL_M3,
       -3.0, 0.0, 2.8e-11 * MOLECULES_PER_MOL_M3, 0.0, 0.0, 0.6, 1.0}},
  {.kind=CHEM_TUNNELING, .reactant_offset=7, .reactant_count=1,
   .product_offset=6, .product_count=1, .p={1.2, 2.3, 302.3}}
};
static const chem_custom_parameter expected_custom[] = {
  {"C.effective radius [m]", 0.0},
  {"C.particle number concentration [# m-3]", 0.0}
};
static const chem_mechanism expected_model = {
  .name="reaction rates no user defined", .species_count=7, .reaction_count=7,
  .reactant_count=8, .product_count=7, .custom_count=2,
  .default_temperature=298.15, .default_pressure=101325.0,
  .species=expected_species, .reactions=expected_reactions,
  .reactants=expected_reactants, .products=expected_products,
  .yields=expected_yields, .custom_parameters=expected_custom
};
#endif

static int nearly_equal(double a, double b)
{
  return isfinite(a) && isfinite(b) &&
         fabs(a - b) <= 1.0e-30 + 1.0e-12 * fabs(b);
}

static int find_species(const chem_mechanism *m, const char *name)
{
  size_t i;
  for (i = 0; i < m->species_count; ++i)
    if (strcmp(m->species[i].name, name) == 0)
      return (int)i;
  return -1;
}

static int find_custom(const chem_mechanism *m, const char *name)
{
  size_t i;
  for (i = 0; i < m->custom_count; ++i)
    if (strcmp(m->custom_parameters[i].name, name) == 0)
      return (int)i;
  return -1;
}

static int same_custom(const chem_mechanism *a, int ai,
                       const chem_mechanism *b, int bi)
{
  if (ai < 0 || bi < 0)
    return ai == bi;
  if ((size_t)ai >= a->custom_count || (size_t)bi >= b->custom_count)
    return 0;
  return strcmp(a->custom_parameters[ai].name, b->custom_parameters[bi].name) == 0;
}

static int mechanisms_match(const chem_mechanism *actual,
                            const chem_mechanism *expected)
{
  size_t i, j;
  if (actual->species_count != expected->species_count ||
      actual->reaction_count != expected->reaction_count ||
      actual->custom_count != expected->custom_count ||
      actual->reactant_count != expected->reactant_count ||
      actual->product_count != expected->product_count ||
      !nearly_equal(actual->default_temperature, expected->default_temperature) ||
      !nearly_equal(actual->default_pressure, expected->default_pressure)) {
    fprintf(stderr, "table sizes/defaults differ (species %zu/%zu reactions %zu/%zu custom %zu/%zu reactants %zu/%zu products %zu/%zu)\n",
            actual->species_count, expected->species_count,
            actual->reaction_count, expected->reaction_count,
            actual->custom_count, expected->custom_count,
            actual->reactant_count, expected->reactant_count,
            actual->product_count, expected->product_count);
    return 0;
  }
  for (i = 0; i < expected->species_count; ++i) {
    int index = find_species(actual, expected->species[i].name);
    const chem_species *a;
    const chem_species *b = &expected->species[i];
    if (index < 0) {
      fprintf(stderr, "missing expected species %s\n", b->name);
      return 0;
    }
    a = &actual->species[index];
    if (!nearly_equal(a->molecular_weight, b->molecular_weight) ||
        !nearly_equal(a->diffusion_coefficient, b->diffusion_coefficient) ||
        !nearly_equal(a->initial_value, b->initial_value) ||
        !nearly_equal(a->absolute_tolerance, b->absolute_tolerance))
      { fprintf(stderr, "species properties differ for %s\n", b->name); return 0; }
  }
  for (i = 0; i < expected->custom_count; ++i) {
    int index = find_custom(actual, expected->custom_parameters[i].name);
    if (index < 0 || !nearly_equal(actual->custom_parameters[index].default_value,
                                   expected->custom_parameters[i].default_value)) {
      fprintf(stderr, "missing custom parameter %s\n", expected->custom_parameters[i].name);
      return 0;
    }
  }
  for (i = 0; i < expected->reaction_count; ++i) {
    const chem_reaction *a = &actual->reactions[i];
    const chem_reaction *b = &expected->reactions[i];
    if (a->kind != b->kind || a->reactant_count != b->reactant_count ||
        a->product_count != b->product_count || a->body_order != b->body_order ||
        a->coefficient_count != b->coefficient_count)
      { fprintf(stderr, "reaction %zu header differs: kind %d/%d reactants %u/%u products %u/%u custom %d/%d %d/%d\n", i, a->kind, b->kind, a->reactant_count, b->reactant_count, a->product_count, b->product_count, a->custom0, b->custom0, a->custom1, b->custom1); return 0; }
    if ((a->kind == CHEM_USER_DEFINED || a->kind == CHEM_SURFACE) &&
        !same_custom(actual, a->custom0, expected, b->custom0))
      { fprintf(stderr, "reaction %zu custom parameter 0 differs\n", i); return 0; }
    if (a->kind == CHEM_SURFACE &&
        !same_custom(actual, a->custom1, expected, b->custom1))
      { fprintf(stderr, "reaction %zu custom parameter 1 differs\n", i); return 0; }
    if (a->kind == CHEM_SURFACE) {
      if (a->surface_species < 0 || b->surface_species < 0 ||
          strcmp(actual->species[a->surface_species].name,
                 expected->species[b->surface_species].name) != 0)
        { fprintf(stderr, "reaction %zu surface species differs\n", i); return 0; }
    }
    for (j = 0; j < 16; ++j)
      if (!nearly_equal(a->p[j], b->p[j]))
        { fprintf(stderr, "reaction %zu parameter %zu differs: %.17g vs %.17g\n", i, j, a->p[j], b->p[j]); return 0; }
    for (j = 0; j < a->reactant_count; ++j) {
      const char *actual_name = actual->species[actual->reactants[a->reactant_offset + j]].name;
      const char *expected_name = expected->species[expected->reactants[b->reactant_offset + j]].name;
      if (strcmp(actual_name, expected_name) != 0)
        { fprintf(stderr, "reaction %zu reactant %zu differs\n", i, j); return 0; }
    }
    for (j = 0; j < a->product_count; ++j) {
      const char *actual_name = actual->species[actual->products[a->product_offset + j]].name;
      const char *expected_name = expected->species[expected->products[b->product_offset + j]].name;
      if (strcmp(actual_name, expected_name) != 0 ||
          !nearly_equal(actual->yields[a->product_offset + j],
                        expected->yields[b->product_offset + j]))
        { fprintf(stderr, "reaction %zu product %zu differs\n", i, j); return 0; }
    }
  }
  return 1;
}

static void set_custom_values(const chem_mechanism *m, double *values,
                              double photolysis_rate)
{
  size_t i;
  for (i = 0; i < m->custom_count; ++i)
    values[i] = m->custom_parameters[i].default_value;
  values[find_custom(m, "C.effective radius [m]")] = 1.0e-7;
  values[find_custom(m, "C.particle number concentration [# m-3]")] = 2.5e6;
#if EXAMPLE_USER_CONFIG
  values[find_custom(m, "my photolysis rate")] = photolysis_rate;
  values[find_custom(m, "my loss rate")] = 1.0e-23;
  values[find_custom(m, "my emission rate")] = 1.0e-20;
#else
  (void)photolysis_rate;
#endif
}

static int initialize_state(const chem_mechanism *m, double *y)
{
  size_t i;
  int a = find_species(m, "A");
  if (a < 0)
    return -1;
  for (i = 0; i < m->species_count; ++i)
    y[i] = m->species[i].initial_value;
  y[a] = 1.0;
  return 0;
}

static int compare_states(const chem_mechanism *a, const double *ya,
                          const chem_mechanism *b, const double *yb)
{
  size_t i;
  for (i = 0; i < a->species_count; ++i) {
    int other = find_species(b, a->species[i].name);
    if (other < 0 || !nearly_equal(ya[i], yb[other]))
      return 0;
  }
  return 1;
}

int main(void)
{
  chem_rosenbrock_parameters parameters;
  double generated_y[16] = {0}, expected_y[16] = {0};
  double generated_custom[16] = {0}, expected_custom_values[16] = {0};
  double temperature = 287.45, pressure = 101319.9, photo_rate = 1.0e-10;
  size_t step;
  if (!mechanisms_match(&chem_model, &expected_model) ||
      initialize_state(&chem_model, generated_y) != 0 ||
      initialize_state(&expected_model, expected_y) != 0) {
    puts("FAIL generated mechanism differs from the hand-built mechanism tables");
    return 1;
  }
  chem_rosenbrock_default(&parameters);
  for (step = 0; step < 10; ++step) {
    chem_solver *generated_solver;
    chem_solver *expected_solver;
    chem_solver_stats generated_stats, expected_stats;
    int generated_result, expected_result;
    set_custom_values(&chem_model, generated_custom, photo_rate);
    set_custom_values(&expected_model, expected_custom_values, photo_rate);
    generated_solver = chem_solver_create(&chem_model, 1, &temperature,
                                         &pressure, generated_custom,
                                         &parameters, 1);
    expected_solver = chem_solver_create(&expected_model, 1, &temperature,
                                         &pressure, expected_custom_values,
                                         &parameters, 1);
    if (!generated_solver || !expected_solver) {
      chem_solver_destroy(generated_solver);
      chem_solver_destroy(expected_solver);
      puts("FAIL could not construct both example solvers");
      return 1;
    }
    generated_result = chem_solver_step(generated_solver, generated_y, 500.0,
                                        &generated_stats);
    expected_result = chem_solver_step(expected_solver, expected_y, 500.0,
                                       &expected_stats);
    chem_solver_destroy(generated_solver);
    chem_solver_destroy(expected_solver);
    if (generated_result != 0 || expected_result != 0 ||
        !compare_states(&chem_model, generated_y, &expected_model, expected_y)) {
      printf("FAIL state comparison at example interval %zu\n", step + 1);
      return 1;
    }
    photo_rate *= 1.5;
  }
  puts("PASS configuration mechanism and ten interval states match hand-built C mechanism");
  return 0;
}
