#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define COUNT(items) (sizeof(items) / sizeof((items)[0]))

static const chem_species flowtube_species[] = {
  { "SOA1", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "SOA2", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "a-pinene", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "O3", 0.0, 0.0, 0.0, 1.0e-3, 0u },
};
static const chem_parameter flowtube_parameters[] = {
  { "r1", 0.01 }, { "r2", 0.05 },
};
static const chem_term flowtube_reactants[] = {
  { 0u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 3u, 1.0 },
};
static const chem_term flowtube_products[] = {
  { 0u, 0.18 }, { 1u, 0.09 },
};
static const chem_reaction flowtube_reactions[] = {
  { 0u, 1u, 0u, 0u, { .kind = CHEM_RATE_USER, .custom_parameter = 0u, .p = { [1] = 1.0 } } },
  { 1u, 1u, 0u, 0u, { .kind = CHEM_RATE_USER, .custom_parameter = 1u, .p = { [1] = 1.0 } } },
  { 2u, 2u, 0u, 2u, { .kind = CHEM_RATE_ARRHENIUS,
      .p = { [0] = 8.8e-17 * 1.0e-6 * 6.02214076e23, [3] = 298.15 } } },
};
static const chem_mechanism flowtube_mechanism = {
  .name = "FlowTube regression", .species_count = COUNT(flowtube_species),
  .reaction_count = COUNT(flowtube_reactions), .reactant_term_count = COUNT(flowtube_reactants),
  .product_term_count = COUNT(flowtube_products), .parameter_count = COUNT(flowtube_parameters),
  .species = flowtube_species, .reactions = flowtube_reactions,
  .reactant_terms = flowtube_reactants, .product_terms = flowtube_products,
  .parameters = flowtube_parameters, .temperature = 298.15, .pressure = 101325.0,
  .relative_tolerance = 1.0e-6,
};

static const double flowtube_expected[2][10][4] = {
  {
    { 0, 0, 8.00000e-8, 2.00000e-5 },
    { 1.51776e-11, 7.43233e-12, 7.99153e-8, 1.99999e-5 },
    { 3.01881e-11, 1.44946e-11, 7.98306e-8, 1.99998e-5 },
    { 4.50332e-11, 2.12049e-11, 7.97460e-8, 1.99997e-5 },
    { 5.97144e-11, 2.75804e-11, 7.96616e-8, 1.99997e-5 },
    { 7.42334e-11, 3.36373e-11, 7.95772e-8, 1.99996e-5 },
    { 8.85919e-11, 3.93913e-11, 7.94929e-8, 1.99995e-5 },
    { 1.02791e-10, 4.48570e-11, 7.94087e-8, 1.99994e-5 },
    { 1.16834e-10, 5.00486e-11, 7.93245e-8, 1.99993e-5 },
    { 1.30720e-10, 5.49794e-11, 7.92405e-8, 1.99992e-5 },
  },
  {
    { 0, 0, 8.00000e-8, 2.00000e-5 },
    { 1.50953e-11, 7.26014e-12, 7.99153e-8, 1.99999e-5 },
    { 3.00252e-11, 1.41668e-11, 7.98307e-8, 1.99998e-5 },
    { 4.47911e-11, 2.07369e-11, 7.97462e-8, 1.99997e-5 },
    { 5.93949e-11, 2.69865e-11, 7.96617e-8, 1.99997e-5 },
    { 7.38381e-11, 3.29307e-11, 7.95774e-8, 1.99996e-5 },
    { 8.81223e-11, 3.85842e-11, 7.94931e-8, 1.99995e-5 },
    { 1.02249e-10, 4.39608e-11, 7.94090e-8, 1.99994e-5 },
    { 1.16220e-10, 4.90738e-11, 7.93249e-8, 1.99993e-5 },
    { 1.30037e-10, 5.39355e-11, 7.92409e-8, 1.99992e-5 },
  },
};

static int flowtube(int use_backward_euler) {
  double state[] = { 0.0, 0.0, 8.0e-8, 2.0e-5 };
  double parameters[] = { 0.01, 0.05 };
  const double temperature[] = { 298.15 }, pressure[] = { 101325.0 };
  const double tolerance = 1.0e-10;
  for (size_t output = 0; output < COUNT(flowtube_expected[0]); ++output) {
    if (output != 0u) {
      const int result = use_backward_euler
          ? chem_integrate_group_backward_euler(&flowtube_mechanism, 1u, state, parameters,
                temperature, pressure, 1.0, NULL)
          : chem_integrate_group(&flowtube_mechanism, 1u, state, parameters,
                temperature, pressure, 1.0, NULL);
      if (result != 0) return 1;
    }
    for (size_t species = 0; species < COUNT(state); ++species) {
      const double expected = flowtube_expected[use_backward_euler][output][species];
      if (!isfinite(state[species]) || fabs(state[species] - expected) > tolerance) {
        fprintf(stderr, "FlowTube %s output %zu species %zu: %.17g vs %.17g\n",
            use_backward_euler ? "backward_euler" : "rosenbrock", output, species,
            state[species], expected);
        return 1;
      }
    }
  }
  return 0;
}

static const chem_species chapman_species[] = {
  { "M", 0.0, 0.0, 0.0, 1.0e-3, 0u }, { "Ar", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "CO2", 0.0, 0.0, 0.0, 1.0e-3, 0u }, { "H2O", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "N2", 0.0, 0.0, 0.0, 1.0e-3, 0u }, { "O1D", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "O", 0.0, 0.0, 0.0, 1.0e-3, 0u }, { "O2", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "O3", 0.0, 0.0, 0.0, 1.0e-3, 0u },
};
static const chem_term chapman_reactants[] = {
  { 4u, 1.0 }, { 5u, 1.0 }, { 5u, 1.0 }, { 7u, 1.0 }, { 6u, 1.0 }, { 8u, 1.0 },
  { 0u, 1.0 }, { 6u, 1.0 }, { 7u, 1.0 }, { 7u, 1.0 }, { 8u, 1.0 }, { 8u, 1.0 },
  { 6u, 1.0 }, { 8u, 1.0 }, { 0u, 1.0 }, { 6u, 1.0 }, { 7u, 1.0 },
};
static const chem_term chapman_products[] = {
  { 6u, 1.0 }, { 4u, 1.0 }, { 6u, 1.0 }, { 7u, 1.0 }, { 7u, 2.0 }, { 6u, 2.0 },
  { 5u, 1.0 }, { 7u, 1.0 }, { 6u, 1.0 }, { 7u, 1.0 }, { 8u, 1.0 }, { 0u, 1.0 },
};
static const chem_reaction chapman_reactions[] = {
  { 0u, 2u, 0u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 2.15e-11, [2] = 110.0, [3] = 300.0 } } },
  { 2u, 2u, 2u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 3.3e-11, [2] = 55.0, [3] = 300.0 } } },
  { 4u, 2u, 4u, 1u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 8.0e-12, [2] = -2060.0, [3] = 300.0 } } },
  { 6u, 3u, 10u, 2u, { .kind = CHEM_RATE_ARRHENIUS, .p = { [0] = 6.0e-34, [1] = 2.4, [3] = 300.0 } } },
  { 9u, 1u, 5u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 10u, 1u, 6u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 1u } },
  { 11u, 1u, 8u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 2u } },
};
static const chem_parameter chapman_parameters[] = {
  { "O2_1", 0.0 }, { "O3_1", 0.0 }, { "O3_2", 0.0 },
};
static const chem_mechanism chapman_mechanism = {
  .name = "Chapman regression", .species_count = COUNT(chapman_species),
  .reaction_count = COUNT(chapman_reactions), .reactant_term_count = COUNT(chapman_reactants),
  .product_term_count = COUNT(chapman_products), .parameter_count = COUNT(chapman_parameters),
  .species = chapman_species, .reactions = chapman_reactions,
  .reactant_terms = chapman_reactants, .product_terms = chapman_products,
  .parameters = chapman_parameters, .temperature = 273.15, .pressure = 100000.0,
  .relative_tolerance = 1.0e-4,
};

static chem_rosenbrock_parameters solver_parameters(const char *method) {
  if (strcmp(method, "TwoStageSolve") == 0) return chem_rosenbrock_two_stage_parameters();
  if (strcmp(method, "FourStageSolve") == 0) return chem_rosenbrock_four_stage_parameters();
  if (strcmp(method, "FourStageDASolve") == 0) return chem_rosenbrock_four_stage_da_parameters();
  if (strcmp(method, "SixStageDASolve") == 0) return chem_rosenbrock_six_stage_da_parameters();
  return chem_rosenbrock_default_parameters();
}

static int chapman_solve(const char *test_name) {
  const int vector_case = strcmp(test_name, "VectorSolve") == 0;
  const size_t first_cells = vector_case ? 1u : 3u;
  const size_t last_cells = vector_case ? 4u : 3u;
  const double temperatures[] = { 284.19, 215.02, 299.31, 267.50 };
  const double pressures[] = { 101245.0, 100789.2, 101398.0, 99950.0 };
  const double photo_rates[4][3] = {
    { 1.0e-4, 1.0e-5, 1.0e-6 }, { 3.2e-4, 7.3e-5, 3.2e-6 },
    { 5.2e-4, 8.2e-5, 4.6e-6 }, { 2.0e-4, 4.0e-5, 2.0e-6 },
  };
  double state[4u * 9u], baseline[4u * 9u], parameters[4u * 3u];
  for (size_t cell = 0; cell < last_cells; ++cell) {
    for (size_t species = 0; species < 9u; ++species)
      state[cell * 9u + species] = 0.04 + 0.07 * (double)(species + 1u) + 0.025 * (double)cell;
    for (size_t parameter = 0; parameter < 3u; ++parameter)
      parameters[cell * 3u + parameter] = photo_rates[cell][parameter];
  }
  memcpy(baseline, state, last_cells * 9u * sizeof(double));
  chem_rosenbrock_parameters reference = chem_rosenbrock_default_parameters();
  reference.h_start_fraction = 0.005;
  reference.h_max_fraction = 0.5;
  reference.max_internal_steps = 1000u;
  for (size_t cell = 0; cell < last_cells; ++cell) {
    double p[3];
    memcpy(p, parameters + cell * 3u, sizeof(p));
    const double t[] = { temperatures[cell] }, pressure[] = { pressures[cell] };
    if (chem_integrate_group_with_parameters(&chapman_mechanism, 1u, baseline + cell * 9u,
            p, t, pressure, 500.0, &reference, NULL) != 0) return 1;
  }

  chem_rosenbrock_parameters solver = solver_parameters(test_name);
  solver.max_internal_steps = 1000u;
  const int two_stage = strcmp(test_name, "TwoStageSolve") == 0;
  const double relative_tolerance = two_stage ? 1.0e-2 : 1.0e-4;
  for (size_t cells = first_cells; cells <= last_cells; ++cells) {
    double actual[4u * 9u];
    memcpy(actual, state, cells * 9u * sizeof(double));
    if (chem_integrate_group_with_parameters(&chapman_mechanism, cells, actual,
            parameters, temperatures, pressures, 500.0, &solver, NULL) != 0) return 1;
    double maximum = 0.0;
    for (size_t i = 0; i < cells * 9u; ++i) maximum = fmax(maximum, state[i]);
    const double absolute_tolerance = maximum * 1.0e-12;
    for (size_t cell = 0; cell < cells; ++cell)
      for (size_t species = 0; species < 9u; ++species) {
        const size_t i = cell * 9u + species;
        const double a = actual[i], b = baseline[i];
        if (!isfinite(a) || fabs(a - b) > (fabs(a) + fabs(b)) * relative_tolerance + absolute_tolerance) {
          fprintf(stderr, "%s cell %zu species %zu: %.17g vs %.17g\n",
              test_name, cell, species, a, b);
          return 1;
        }
      }
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  if (strcmp(argv[1], "RosenbrockRegression.FlowTube") == 0) return flowtube(0);
  if (strcmp(argv[1], "BackwardEulerRegression.FlowTube") == 0) return flowtube(1);
  const char *name = strrchr(argv[1], '.');
  name = name ? name + 1 : argv[1];
  if (strcmp(name, "TwoStageSolve") == 0 || strcmp(name, "ThreeStageSolve") == 0 ||
      strcmp(name, "FourStageSolve") == 0 || strcmp(name, "FourStageDASolve") == 0 ||
      strcmp(name, "SixStageDASolve") == 0 || strcmp(name, "VectorSolve") == 0)
    return chapman_solve(name);
  return 2;
}
