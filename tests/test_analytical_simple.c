#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
  int model;
  unsigned rate;
  unsigned fast_rate;
} analytical_rate_context;

static const chem_species simple_species[] = {
  { "A", 0.0, 0.0, 0.0, 1.0e-14, 0u },
  { "B", 0.0, 0.0, 0.0, 1.0e-14, 0u },
  { "C", 0.0, 0.0, 0.0, 1.0e-14, 0u },
};
static const chem_species stiff_species[] = {
  { "A1", 0.0, 0.0, 0.0, 1.0e-14, 0u },
  { "A2", 0.0, 0.0, 0.0, 1.0e-14, 0u },
  { "B", 0.0, 0.0, 0.0, 1.0e-14, 0u },
  { "C", 0.0, 0.0, 0.0, 1.0e-14, 0u },
};
static const chem_term reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
static const chem_term products[] = { { 1u, 1.0 }, { 2u, 1.0 } };
static const chem_term stiff_reactants[] = {
  { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 0u, 1.0 }, { 1u, 1.0 },
};
static const chem_term stiff_products[] = {
  { 2u, 1.0 }, { 2u, 1.0 }, { 3u, 1.0 }, { 1u, 1.0 }, { 0u, 1.0 },
};
static analytical_rate_context contexts[5];

static double troe_rate(double temperature, double air_density, unsigned reaction) {
  const double k0 = reaction == 0u ? 4.0e-11 : 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 1.6);
  const double kinf = reaction == 0u ? 1.0 : 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  const double fc = reaction == 0u ? 0.6 : 0.9;
  const double n = reaction == 0u ? 1.0 : 0.8;
  const double ratio = k0 * air_density / kinf;
  return k0 * air_density / (1.0 + ratio) * pow(fc, n / (n + pow(log10(ratio), 2.0)));
}

static double ternary_rate(double temperature, double air_density, unsigned reaction) {
  const double k0 = reaction == 0u ? 4.0e-5 : 1.2e-3 * exp(3.0 / temperature) * pow(temperature / 300.0, 1.6);
  const double kinf = reaction == 0u ? 1.0 : 136.0 * exp(24.0 / temperature) * pow(temperature / 300.0, 5.0);
  const double fc = reaction == 0u ? 0.6 : 0.9;
  const double n = reaction == 0u ? 1.0 : 0.8;
  const double ratio = k0 * air_density / kinf;
  return k0 / (1.0 + ratio) * pow(fc, n / (n + pow(log10(ratio), 2.0)));
}

static double branched_rate(double temperature, double air_density, unsigned reaction) {
  const double air_number_density = air_density * 6.02214076e23 * 1.0e-6;
  const double a_ref = 2.0e-22 * exp(2.0) * 2.45e19;
  const double b_ref = 0.43 * pow(293.0 / 298.0, -8.0);
  const double reference = a_ref / (1.0 + a_ref / b_ref) *
      pow(0.41, 1.0 / (1.0 + pow(log10(a_ref / b_ref), 2.0)));
  const double z = reference * (1.0 - 1.0e-3) / 1.0e-3;
  const double a = 2.0e-22 * exp(2.0) * air_number_density;
  const double b = 0.43 * pow(temperature / 298.0, -8.0);
  const double value = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
  return 1.0e-4 * exp(-204.3 / temperature) * (reaction == 0u ? z : value) /
      (z + value);
}

static double analytical_rate(double temperature, double pressure, double air_density, void *opaque) {
  const analytical_rate_context *context = opaque;
  if (context->fast_rate == 1u) return 4.0e10;
  if (context->fast_rate == 2u) return 3.6e10;
  const unsigned rate = context->rate;
  if (context->model == 0) return troe_rate(temperature, air_density, rate);
  if (context->model == 1) return rate == 0u ? 2.0e-3 : 3.0e-3;
  if (context->model == 2) return ternary_rate(temperature, air_density, rate);
  if (context->model == 3) {
    return rate == 0u ? 4.0e-3 :
        1.2e-4 * exp(-167.0 / temperature + 1.0e8 / (temperature * temperature * temperature));
  }
  if (context->model == 4) {
    return rate == 0u ? 4.0e-3 * exp(50.0 / temperature) :
        1.2e-4 * exp(75.0 / temperature) * pow(temperature / 50.0, 7.0) * (1.0 + 0.5 * pressure);
  }
  return branched_rate(temperature, air_density, rate);
}

static int model_and_tolerance(const char *case_name, int stiff, int backward_euler, int *model, double *tolerance) {
  if (strcmp(case_name, "Troe") == 0) { *model = 0; *tolerance = stiff ? 1.0e-5 : 1.0e-6; }
  else if (strcmp(case_name, "Photolysis") == 0) { *model = 1; *tolerance = stiff ? (backward_euler ? 1.0e-3 : 2.0e-5) : 1.0e-3; }
  else if (strcmp(case_name, "TernaryChemicalActivation") == 0) {
    *model = 2; *tolerance = stiff ? (backward_euler ? 1.0e-2 : 2.0e-3) : 1.0e-5;
  }
  else if (strcmp(case_name, "Tunneling") == 0) { *model = 3; *tolerance = stiff ? (backward_euler ? 1.0e-3 : 1.0e-4) : 1.0e-3; }
  else if (strcmp(case_name, "Arrhenius") == 0) { *model = 4; *tolerance = stiff ? (backward_euler ? 1.0e-3 : 1.0e-4) : 1.0e-3; }
  else if (strcmp(case_name, "Branched") == 0) { *model = 5; *tolerance = stiff ? (backward_euler ? 1.0e-2 : 2.0e-3) : 1.0e-5; }
  else return -1;
  return 0;
}

static int run_analytical(const char *solver_name, const char *case_name, int stiff) {
  int model = -1;
  double absolute_tolerance = 0.0;
  const int backward_euler = strcmp(solver_name, "backward_euler") == 0;
  if (model_and_tolerance(case_name, stiff, backward_euler, &model, &absolute_tolerance) != 0) return 2;
  if (!backward_euler && strcmp(solver_name, "rosenbrock") != 0) return 2;
  const size_t reaction_count = stiff ? 5u : 2u;
  for (size_t i = 0; i < reaction_count; ++i) {
    const unsigned rate = stiff ? (i == 2u ? 1u : 0u) : (unsigned)i;
    const unsigned fast_rate = stiff && i >= 3u ? (unsigned)(i - 2u) : 0u;
    contexts[i] = (analytical_rate_context){ model, rate, fast_rate };
  }
  const chem_reaction simple_reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = analytical_rate, .context = &contexts[0] } },
    { 1u, 1u, 1u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = analytical_rate, .context = &contexts[1] } },
  };
  const chem_reaction stiff_reactions[] = {
    { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = analytical_rate, .context = &contexts[0] } },
    { 1u, 1u, 0u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = analytical_rate, .context = &contexts[1] } },
    { 2u, 1u, 2u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = analytical_rate, .context = &contexts[2] } },
    { 3u, 1u, 3u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = analytical_rate, .context = &contexts[3] } },
    { 4u, 1u, 4u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = analytical_rate, .context = &contexts[4] } },
  };
  const chem_reaction *reactions = stiff ? stiff_reactions : simple_reactions;
  const chem_mechanism mechanism = {
    .name = stiff ? "analytical stiff A to B to C" : "analytical A to B to C",
    .species_count = stiff ? 4u : 3u,
    .reaction_count = reaction_count,
    .reactant_term_count = stiff ? 5u : 2u,
    .product_term_count = stiff ? 5u : 2u,
    .species = stiff ? stiff_species : simple_species,
    .reactions = reactions,
    .reactant_terms = stiff ? stiff_reactants : reactants,
    .product_terms = stiff ? stiff_products : products,
    .temperature = 298.15,
    .pressure = 101325.0,
    .relative_tolerance = 1.0e-6,
  };
  static const double temperature[] = { 272.5, 254.7, 312.6 };
  static const double pressure[] = { 101253.3, 100672.5, 101319.8 };
  double rates[15];
  const size_t species_count = stiff ? 4u : 3u;
  double state[12] = { 0 };
  double initial_a[] = { 1.0, 2.0 / 3.0, 1.0 / 3.0 };
  for (size_t cell = 0; cell < 3u; ++cell) {
    if (stiff) {
      state[cell * species_count] = initial_a[cell] / 2.0;
      state[cell * species_count + 1u] = initial_a[cell] / 2.0;
    } else {
      state[cell * species_count] = initial_a[cell];
    }
    for (size_t reaction = 0; reaction < reaction_count; ++reaction) {
      if (chem_evaluate_rate(&mechanism, reaction, NULL, temperature[cell], pressure[cell],
              &rates[cell * reaction_count + reaction]) != 0) return 1;
    }
  }
  for (size_t step = 1u; step < 1000u; ++step) {
    const int status = backward_euler
        ? chem_integrate_group_backward_euler(&mechanism, 3u, state, NULL, temperature, pressure, 1.0, NULL)
        : chem_integrate_group(&mechanism, 3u, state, NULL, temperature, pressure, 1.0, NULL);
    if (status != 0) return 1;
    const double time = (double)step;
    for (size_t cell = 0; cell < 3u; ++cell) {
      const double k1 = rates[cell * reaction_count], k2 = rates[cell * reaction_count + (stiff ? 2u : 1u)];
      const double exp1 = exp(-k1 * time), exp2 = exp(-k2 * time);
      const double expected[] = {
        initial_a[cell] * exp1,
        initial_a[cell] * (k1 / (k2 - k1)) * (exp1 - exp2),
        initial_a[cell] * (1.0 + (k1 * exp2 - k2 * exp1) / (k2 - k1)),
      };
      for (size_t s = 0; s < 3u; ++s) {
        const double actual = stiff
            ? (s == 0u ? state[cell * species_count] + state[cell * species_count + 1u]
                       : state[cell * species_count + s + 1u])
            : state[cell * species_count + s];
        if (!isfinite(actual) || fabs(actual - expected[s]) > absolute_tolerance) {
          fprintf(stderr, "%s %s%s step %zu cell %zu species %zu: %.17g vs %.17g (tol %.3g)\n",
              solver_name, case_name, stiff ? "SuperStiffButAnalytical" : "", step, cell, s,
              actual, expected[s], absolute_tolerance);
          return 1;
        }
      }
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 3) return 2;
  const char *case_name = strrchr(argv[2], '.');
  case_name = case_name ? case_name + 1 : argv[2];
  const size_t name_length = strlen(case_name);
  static const char suffix[] = "SuperStiffButAnalytical";
  if (name_length > sizeof(suffix) - 1u &&
      strcmp(case_name + name_length - (sizeof(suffix) - 1u), suffix) == 0) {
    char base_name[64];
    const size_t base_length = name_length - (sizeof(suffix) - 1u);
    if (base_length >= sizeof(base_name)) return 2;
    memcpy(base_name, case_name, base_length);
    base_name[base_length] = '\0';
    return run_analytical(argv[1], base_name, 1);
  }
  return run_analytical(argv[1], case_name, 0);
}
