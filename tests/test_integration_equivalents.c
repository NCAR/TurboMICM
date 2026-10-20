#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static const chem_species terminator_species[] = {
  { "Cl2", 0.0, 0.0, 0.0, 1.0e-20, 0u },
  { "Cl", 0.0, 0.0, 0.0, 1.0e-20, 0u },
};
static const chem_parameter terminator_parameters[] = { { "toy_k1", 0.0 } };
static const chem_term terminator_reactants[] = { { 0u, 1.0 }, { 1u, 1.0 }, { 1u, 1.0 } };
static const chem_term terminator_products[] = { { 1u, 2.0 }, { 0u, 1.0 } };
static const chem_reaction terminator_reactions[] = {
  { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 1u, 2u, 1u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 1.0, [3] = 300.0 } } },
};
static const chem_mechanism terminator_mechanism = {
  .name = "chlorine terminator",
  .species_count = 2,
  .reaction_count = 2,
  .reactant_term_count = 3,
  .product_term_count = 2,
  .parameter_count = 1,
  .species = terminator_species,
  .reactions = terminator_reactions,
  .reactant_terms = terminator_reactants,
  .product_terms = terminator_products,
  .parameters = terminator_parameters,
  .temperature = 298.0,
  .pressure = 101300.0,
  .relative_tolerance = 1.0e-8,
};

static double terminator_photolysis(double latitude, double longitude) {
  const double lat_center = 3.14159265358979323846 * 20.0 / 180.0;
  const double lon_center = 3.14159265358979323846 * 300.0 / 180.0;
  const double value = sin(latitude) * sin(lat_center) +
      cos(latitude) * cos(lat_center) * cos(longitude - lon_center);
  return fmax(0.0, value);
}

static int terminator_case(void) {
  const double pi = 3.14159265358979323846;
  chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  solver.max_internal_steps = 100000u;
  const double dt = 30.0;
  const size_t steps = (size_t)floor(2.0 * pi / 0.3);
  int passed = 1;
  for (size_t cells = 1; cells <= 4; ++cells) {
    double original[8] = { 0.0 };
    double state[8] = { 0.0 };
    double custom[4] = { 0.0 };
    double temperature[4], pressure[4];
    for (size_t cell = 0; cell < cells; ++cell) {
      const double factor = (double)(cell + 1u);
      original[2u * cell] = factor * 1.0e-6;
      original[2u * cell + 1u] = (1.0 + 0.1 * (double)cell) * 1.0e-10;
      temperature[cell] = 298.0;
      pressure[cell] = 101300.0;
    }
    for (size_t step = 0; step < steps; ++step) {
      memcpy(state, original, 2u * cells * sizeof(double));
      const double longitude = (double)step * 0.3;
      for (size_t cell = 0; cell < cells; ++cell) {
        const double latitude = pi / 180.0 * ((double)cell * (90.0 / (double)cells));
        custom[cell] = terminator_photolysis(latitude, longitude);
      }
      if (chem_integrate_group_with_parameters(
              &terminator_mechanism, cells, state, custom, temperature, pressure, dt, &solver, NULL) != 0) {
        passed = 0;
        break;
      }
      for (size_t cell = 0; cell < cells; ++cell) {
        const double k1 = custom[cell];
        const double k2 = 1.0;
        const double cl_initial = original[2u * cell + 1u];
        const double cl2_initial = original[2u * cell];
        const double r = k1 / (4.0 * k2);
        const double cly = cl_initial + 2.0 * cl2_initial;
        const double det = sqrt(r * r + 2.0 * r * cly);
        const double exponent = exp(-4.0 * k2 * det * dt);
        const double l = det * k2 * dt > 1.0e-16 ? (1.0 - exponent) / det / dt : 4.0 * k2;
        const double cl_flux = -l * (cl_initial - det + r) * (cl_initial + det + r) /
            (1.0 + exponent + dt * l * (cl_initial + r));
        const double expected_cl = cl_initial + dt * cl_flux;
        const double expected_cl2 = cl2_initial - 0.5 * dt * cl_flux;
        const double tolerance_cl = fabs(expected_cl) * 1.0e-6 + 1.0e-14;
        const double tolerance_cl2 = fabs(expected_cl2) * 1.0e-6 + 1.0e-14;
        if (!isfinite(state[2u * cell]) || !isfinite(state[2u * cell + 1u]) ||
            fabs(state[2u * cell + 1u] - expected_cl) > tolerance_cl ||
            fabs(state[2u * cell] - expected_cl2) > tolerance_cl2) {
          fprintf(stderr, "terminator step %zu cell %zu: Cl2 %.17g/%.17g Cl %.17g/%.17g\n",
                  step, cell, state[2u * cell], expected_cl2, state[2u * cell + 1u], expected_cl);
          passed = 0;
        }
      }
    }
  }
  return passed;
}

static const chem_species integrated_species[] = {
  { "a", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "b", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "c", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "irr_1", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "irr_2", 0.0, 0.0, 0.0, 1.0e-3, 0u },
};
static const chem_parameter integrated_parameters[] = { { "r2", 1.0 } };
static const chem_term integrated_reactants[] = { { 0u, 1.0 }, { 1u, 1.0 } };
static const chem_term integrated_products[] = { { 1u, 1.0 }, { 3u, 1.0 }, { 2u, 1.0 }, { 4u, 1.0 } };
static const chem_reaction integrated_reactions[] = {
  { 0u, 1u, 0u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 2.15e-11, [3] = 300.0, [2] = 110.0 } } },
  { 1u, 1u, 2u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
};
static const chem_mechanism integrated_mechanism = {
  .name = "integrated reaction rates example",
  .species_count = 5,
  .reaction_count = 2,
  .reactant_term_count = 2,
  .product_term_count = 4,
  .parameter_count = 1,
  .species = integrated_species,
  .reactions = integrated_reactions,
  .reactant_terms = integrated_reactants,
  .product_terms = integrated_products,
  .parameters = integrated_parameters,
  .temperature = 273.0,
  .pressure = 1000.0,
  .relative_tolerance = 1.0e-6,
};

static int integrated_reaction_rates_case(void) {
  double state[] = { 1.0, 0.0, 0.0, 0.0, 0.0 };
  const double temperature[] = { 273.0 }, pressure[] = { 1000.0 };
  double custom[] = { 1.0 };
  for (size_t step = 0; step < 100; ++step) {
    if (step > 50u) custom[0] = 0.0;
    const double irr1_before = state[3], irr2_before = state[4];
    if (chem_integrate_group(&integrated_mechanism, 1u, state, custom, temperature, pressure, 30.0, NULL) != 0 ||
        !isfinite(state[3]) || !isfinite(state[4]) || state[3] < irr1_before || state[4] < irr2_before)
      return 0;
  }
  return 1;
}

static const chem_species chapman_species[] = {
  { "O", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "O1D", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "O2", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "O3", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "M", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "Ar", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "N2", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "H2O", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "CO2", 0.0, 0.0, 0.0, 1.0e-3, 0u },
};
static const chem_parameter chapman_parameters[] = { { "jO2", 0.1 }, { "jO3a", 0.2 }, { "jO3b", 0.3 } };
static const chem_term chapman_reactants[] = {
  { 1u, 1.0 }, { 6u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 0u, 1.0 }, { 3u, 1.0 },
  { 0u, 1.0 }, { 2u, 1.0 }, { 4u, 1.0 }, { 2u, 1.0 }, { 3u, 1.0 }, { 3u, 1.0 },
};
static const chem_term chapman_products[] = {
  { 0u, 1.0 }, { 6u, 1.0 }, { 0u, 1.0 }, { 2u, 1.0 }, { 2u, 2.0 }, { 3u, 1.0 },
  { 4u, 1.0 }, { 0u, 2.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 0u, 1.0 }, { 2u, 1.0 },
};
static const chem_reaction chapman_reactions[] = {
  { 0u, 2u, 0u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 2.15e-11, [2] = 110.0, [3] = 300.0 } } },
  { 2u, 2u, 2u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 3.3e-11, [2] = 55.0, [3] = 300.0 } } },
  { 4u, 2u, 4u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 8.0e-12, [2] = -2060.0, [3] = 300.0 } } },
  { 6u, 3u, 5u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 6.0e-34, [2] = 2.4, [3] = 300.0 } } },
  { 9u, 1u, 7u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 10u, 1u, 8u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 1u } },
  { 11u, 1u, 10u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 2u } },
};
static const chem_mechanism chapman_mechanism = {
  .name = "Chapman integration invariant example",
  .species_count = 9,
  .reaction_count = 7,
  .reactant_term_count = 12,
  .product_term_count = 12,
  .parameter_count = 3,
  .species = chapman_species,
  .reactions = chapman_reactions,
  .reactant_terms = chapman_reactants,
  .product_terms = chapman_products,
  .parameters = chapman_parameters,
  .temperature = 2.0,
  .pressure = 3.0,
  .relative_tolerance = 1.0e-6,
};

static int chapman_invariants_case(void) {
  double state[] = { 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3 };
  double custom[] = { 0.1, 0.2, 0.3 };
  const double temperature[] = { 2.0 }, pressure[] = { 3.0 };
  const double oxygen_initial = state[0] + state[1] + 2.0 * state[2] + 3.0 * state[3];
  const double spectators[] = { state[4], state[5], state[6], state[7], state[8] };
  for (size_t step = 0; step < 100; ++step) {
    if (chem_integrate_group(&chapman_mechanism, 1u, state, custom, temperature, pressure, 30.0, NULL) != 0)
      return 0;
    for (size_t species = 0; species < 9; ++species)
      if (!isfinite(state[species]) || state[species] < 0.0) return 0;
    const double oxygen = state[0] + state[1] + 2.0 * state[2] + 3.0 * state[3];
    if (fabs(oxygen - oxygen_initial) > 1.0e-8 * fmax(1.0, fabs(oxygen_initial))) return 0;
    for (size_t i = 0; i < 5; ++i)
      if (fabs(state[i + 4u] - spectators[i]) > 1.0e-8 * fmax(1.0, fabs(spectators[i]))) return 0;
  }
  return state[2] < 0.05 && state[0] > 0.05 && state[0] > 0.5 * oxygen_initial;
}

static int runtime_vector_use_case(void) {
  const double values[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  double sum = 0.0;
  for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i)
    sum += values[i];
  return sum == 15.0;
}

static double callback_scale = 2.0;
static double temperature_rate_callback(double temperature, double pressure, double air_density, void *context) {
  (void)pressure;
  (void)air_density;
  return *(double *)context * temperature;
}
static const chem_species callback_species[] = {
  { "A", 0.0, 0.0, 0.0, 1.0e-12, 0u },
  { "B", 0.0, 0.0, 0.0, 1.0e-12, 0u },
};
static const chem_term callback_reactants[] = { { 0u, 1.0 } };
static const chem_term callback_products[] = { { 1u, 1.0 } };
static const chem_reaction callback_reactions[] = {
  { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_CALLBACK, .callback = temperature_rate_callback, .context = &callback_scale } },
};
static const chem_mechanism callback_mechanism = {
  .name = "temperature callback rate",
  .species_count = 2,
  .reaction_count = 1,
  .reactant_term_count = 1,
  .product_term_count = 1,
  .parameter_count = 0,
  .species = callback_species,
  .reactions = callback_reactions,
  .reactant_terms = callback_reactants,
  .product_terms = callback_products,
  .temperature = 300.0,
  .pressure = 101253.3,
  .relative_tolerance = 1.0e-9,
};

static int callback_rate_case(int solver_override) {
  double rate = 0.0;
  const double temperature[] = { 300.0 }, pressure[] = { 101253.3 };
  if (solver_override) {
    callback_scale = 1.0e-3;
    if (chem_evaluate_rate(&callback_mechanism, 0u, NULL, temperature[0], pressure[0], &rate) != 0 ||
        fabs(rate - 0.3) > 1.0e-14)
      return 0;
    callback_scale = 2.0e-3;
    if (chem_evaluate_rate(&callback_mechanism, 0u, NULL, temperature[0], pressure[0], &rate) != 0 ||
        fabs(rate - 0.6) > 1.0e-14)
      return 0;
    double state[] = { 1.0, 0.0 };
    if (chem_integrate_group(&callback_mechanism, 1u, state, NULL, temperature, pressure, 1.0, NULL) != 0)
      return 0;
    const double expected_a = exp(-0.6);
    return fabs(state[0] - expected_a) < 1.0e-7 && fabs(state[1] - (1.0 - expected_a)) < 1.0e-7;
  }
  callback_scale = 2.0;
  return chem_evaluate_rate(&callback_mechanism, 0u, NULL, temperature[0], pressure[0], &rate) == 0 &&
      fabs(rate - 600.0) <= 1.0e-13 * 600.0;
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  if (strcmp(argv[1], "terminator") == 0) return terminator_case() ? 0 : 1;
  if (strcmp(argv[1], "integrated_reaction_rates") == 0)
    return integrated_reaction_rates_case() ? 0 : 1;
  if (strcmp(argv[1], "chapman_integration") == 0) return chapman_invariants_case() ? 0 : 1;
  if (strcmp(argv[1], "runtime_vector_use") == 0) return runtime_vector_use_case() ? 0 : 1;
  if (strcmp(argv[1], "lambda_rate_constant") == 0) return callback_rate_case(0) ? 0 : 1;
  if (strcmp(argv[1], "solver_lambda_rate_constant") == 0) return callback_rate_case(1) ? 0 : 1;
  return 2;
}
