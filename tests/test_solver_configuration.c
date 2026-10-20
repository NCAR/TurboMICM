#include "chem.h"

#include <math.h>
#include <stdio.h>

static const chem_species species[] = {
  { "A", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "B", 0.0, 0.0, 0.0, 1.0e-3, 0u },
  { "C", 0.0, 0.0, 0.0, 1.0e-3, 0u },
};
static const chem_parameter parameters[] = {
  { "r1", 0.04 }, { "r2", 3.0e7 }, { "r3", 1.0e4 },
};
static const chem_term reactants[] = {
  { 0u, 1.0 }, { 1u, 1.0 }, { 1u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 },
};
static const chem_term products[] = {
  { 1u, 1.0 }, { 1u, 1.0 }, { 2u, 1.0 }, { 0u, 1.0 }, { 2u, 1.0 },
};
static const chem_reaction reactions[] = {
  { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
  { 1u, 2u, 1u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 1u } },
  { 3u, 2u, 3u, 2u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 2u } },
};
static const chem_mechanism mechanism = {
  .name = "solver configuration example",
  .species_count = 3,
  .reaction_count = 3,
  .reactant_term_count = 5,
  .product_term_count = 5,
  .parameter_count = 3,
  .species = species,
  .reactions = reactions,
  .reactant_terms = reactants,
  .product_terms = products,
  .parameters = parameters,
  .temperature = 272.5,
  .pressure = 101253.3,
  .relative_tolerance = 1.0e-6,
};

static const double expected[][11][3] = {
  {
    { 1.0, 0.0, 0.0 }, { 0.537, 4.50e-6, 0.463 }, { 0.451, 3.23e-6, 0.549 },
    { 0.401, 2.64e-6, 0.599 }, { 0.365, 2.28e-6, 0.635 }, { 0.338, 2.02e-6, 0.662 },
    { 0.316, 1.83e-6, 0.684 }, { 0.297, 1.68e-6, 0.703 }, { 0.282, 1.56e-6, 0.718 },
    { 0.268, 1.46e-6, 0.732 }, { 0.256, 1.37e-6, 0.744 },
  },
  {
    { 1.0, 0.0, 0.0 }, { 0.535, 4.48e-6, 0.465 }, { 0.450, 3.23e-6, 0.550 },
    { 0.400, 2.63e-6, 0.600 }, { 0.364, 2.27e-6, 0.636 }, { 0.337, 2.01e-6, 0.663 },
    { 0.315, 1.82e-6, 0.685 }, { 0.296, 1.67e-6, 0.704 }, { 0.281, 1.55e-6, 0.719 },
    { 0.267, 1.45e-6, 0.733 }, { 0.255, 1.37e-6, 0.745 },
  },
  {
    { 1.0, 0.0, 0.0 }, { 0.536, 4.49e-6, 0.464 }, { 0.450, 3.20e-6, 0.550 },
    { 0.400, 2.62e-6, 0.600 }, { 0.364, 2.26e-6, 0.636 }, { 0.337, 2.01e-6, 0.663 },
    { 0.315, 1.82e-6, 0.685 }, { 0.297, 1.67e-6, 0.703 }, { 0.281, 1.55e-6, 0.719 },
    { 0.267, 1.45e-6, 0.733 }, { 0.256, 1.37e-6, 0.744 },
  },
  {
    { 1.0, 0.0, 0.0 }, { 0.536, 4.49e-6, 0.464 }, { 0.451, 3.23e-6, 0.549 },
    { 0.400, 2.63e-6, 0.600 }, { 0.364, 2.27e-6, 0.636 }, { 0.337, 2.01e-6, 0.663 },
    { 0.315, 1.82e-6, 0.685 }, { 0.297, 1.68e-6, 0.704 }, { 0.281, 1.55e-6, 0.719 },
    { 0.268, 1.45e-6, 0.732 }, { 0.256, 1.37e-6, 0.744 },
  },
  {
    { 1.0, 0.0, 0.0 }, { 0.536, 4.49e-6, 0.464 }, { 0.451, 3.22e-6, 0.549 },
    { 0.400, 2.63e-6, 0.600 }, { 0.364, 2.27e-6, 0.636 }, { 0.337, 2.01e-6, 0.663 },
    { 0.315, 1.82e-6, 0.685 }, { 0.297, 1.67e-6, 0.704 }, { 0.281, 1.55e-6, 0.719 },
    { 0.267, 1.45e-6, 0.733 }, { 0.256, 1.37e-6, 0.744 },
  },
};

static const char *method_names[] = { "two-stage", "three-stage", "four-stage", "four-stage-DA", "six-stage-DA" };
static chem_rosenbrock_parameters (*const method_builders[])(void) = {
  chem_rosenbrock_two_stage_parameters,
  chem_rosenbrock_default_parameters,
  chem_rosenbrock_four_stage_parameters,
  chem_rosenbrock_four_stage_da_parameters,
  chem_rosenbrock_six_stage_da_parameters,
};
static const chem_statistics expected_statistics[] = {
  { 149, 149, 0, 298, 149, 149, 298 },
  { 98, 98, 0, 196, 98, 98, 294 },
  { 97, 97, 0, 291, 97, 97, 388 },
  { 99, 99, 0, 297, 99, 99, 396 },
  { 97, 97, 0, 582, 97, 97, 582 },
};

static int near_printed(double actual, double printed) {
  const double unit = pow(10.0, floor(log10(fabs(printed))) - 2.0);
  return isfinite(actual) && fabs(actual - printed) <= 0.51 * unit + 1.0e-30;
}

static int run_method(size_t method) {
  int passed = 1;
  chem_rosenbrock_parameters solver = method_builders[method]();
  solver.lu_method = CHEM_LU_DOOLITTLE;
  double state[] = { 1.0, 0.0, 0.0 };
  double custom[] = { 0.04, 3.0e7, 1.0e4 };
  const double temperature[] = { 272.5 }, pressure[] = { 101253.3 };
  chem_statistics stats = { 0 };
  for (size_t step = 0; step < 10; ++step) {
    if (chem_integrate_group_with_parameters(&mechanism, 1, state, custom, temperature, pressure,
                                             200.0, &solver, &stats) != 0)
      return 0;
    for (size_t species_index = 0; species_index < 3; ++species_index)
      if (!near_printed(state[species_index], expected[method][step + 1][species_index])) {
        fprintf(stderr, "%s step %zu species %zu: %.17g expected %.17g\n", method_names[method],
                step + 1, species_index, state[species_index], expected[method][step + 1][species_index]);
        passed = 0;
      }
  }
  if (stats.internal_steps != expected_statistics[method].internal_steps ||
      stats.accepted != expected_statistics[method].accepted || stats.rejected != expected_statistics[method].rejected ||
      stats.function_calls != expected_statistics[method].function_calls ||
      stats.jacobian_updates != expected_statistics[method].jacobian_updates ||
      stats.decompositions != expected_statistics[method].decompositions ||
      stats.solves != expected_statistics[method].solves) {
    fprintf(stderr, "%s stats steps=%lu accepted=%lu functions=%lu jacobians=%lu solves=%lu\n",
            method_names[method], stats.internal_steps, stats.accepted, stats.function_calls,
            stats.jacobian_updates, stats.solves);
    passed = 0;
  }
  return passed;
}

int main(void) {
  int passed = 1;
  for (size_t method = 0; method < sizeof(method_builders) / sizeof(method_builders[0]); ++method)
    if (!run_method(method)) passed = 0;
  return passed ? 0 : 1;
}
