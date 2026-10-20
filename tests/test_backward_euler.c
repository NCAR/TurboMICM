#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define REQUIRE(condition) do { \
  if (!(condition)) { \
    fprintf(stderr, "%s:%d check failed: %s\n", __FILE__, __LINE__, #condition); \
    return 1; \
  } \
} while (0)

static const chem_species species[] = {
  { "A", 0.0, 0.0, 0.0, 1.0e-12, 0u },
  { "B", 0.0, 0.0, 0.0, 1.0e-12, 0u },
};
static const chem_parameter parameters[] = { { "k", 2.0 } };
static const chem_term reactants[] = { { 0u, 1.0 } };
static const chem_term products[] = { { 1u, 1.0 } };
static const chem_reaction reactions[] = {
  { 0u, 1u, 0u, 1u, { .kind = CHEM_RATE_USER, .p = { [1] = 1.0 }, .custom_parameter = 0u } },
};
static const chem_mechanism mechanism = {
  .name = "Backward Euler transfer",
  .species_count = 2u,
  .reaction_count = 1u,
  .reactant_term_count = 1u,
  .product_term_count = 1u,
  .parameter_count = 1u,
  .species = species,
  .reactions = reactions,
  .reactant_terms = reactants,
  .product_terms = products,
  .parameters = parameters,
  .temperature = 298.15,
  .pressure = 101325.0,
  .relative_tolerance = 1.0e-6,
};

static int solve_checks(void) {
  const double y0[] = { 1.0, 0.0, 2.0, 0.0, 3.0, 0.0 };
  const double rate_parameters[] = { 2.0, 3.0, 4.0 };
  const double temperature[] = { 298.15, 298.15, 298.15 };
  const double pressure[] = { 101325.0, 101325.0, 101325.0 };
  double grouped[6];
  memcpy(grouped, y0, sizeof(grouped));
  chem_statistics statistics = { 0 };
  REQUIRE(chem_integrate_group_backward_euler(&mechanism, 3u, grouped, rate_parameters,
      temperature, pressure, 1.0, &statistics) == 0);
  REQUIRE(statistics.internal_steps == 2u && statistics.accepted == 1u && statistics.rejected == 0u);
  REQUIRE(statistics.function_calls == 2u && statistics.jacobian_updates == 2u &&
      statistics.decompositions == 2u && statistics.solves == 2u);
  for (size_t cell = 0; cell < 3u; ++cell) {
    double separate[2] = { y0[2u * cell], y0[2u * cell + 1u] };
    const double *cell_parameters = rate_parameters + cell;
    REQUIRE(chem_integrate_group_backward_euler(&mechanism, 1u, separate, cell_parameters,
        temperature + cell, pressure + cell, 1.0, NULL) == 0);
    const double expected_a = y0[2u * cell] / (1.0 + rate_parameters[cell]);
    REQUIRE(fabs(grouped[2u * cell] - expected_a) < 1.0e-12);
    REQUIRE(fabs(grouped[2u * cell + 1u] - (y0[2u * cell] - expected_a)) < 1.0e-12);
    REQUIRE(fabs(grouped[2u * cell] - separate[0]) < 1.0e-14);
    REQUIRE(fabs(grouped[2u * cell + 1u] - separate[1]) < 1.0e-14);
  }
  double explicit_parameters_state[] = { 1.0, 0.0 };
  const double explicit_rate_parameter[] = { 2.0 };
  const double one_temperature[] = { 298.15 }, one_pressure[] = { 101325.0 };
  chem_backward_euler_parameters custom = chem_backward_euler_default_parameters();
  REQUIRE(chem_integrate_group_backward_euler_with_parameters(&mechanism, 1u, explicit_parameters_state,
      explicit_rate_parameter, one_temperature, one_pressure, 1.0, &custom, NULL) == 0);
  REQUIRE(fabs(explicit_parameters_state[0] - 1.0 / 3.0) < 1.0e-12);
  return 0;
}

static int missing_mechanism_check(void) {
  double state[] = { 1.0, 0.0 };
  const double temperature[] = { 298.15 }, pressure[] = { 101325.0 };
  REQUIRE(chem_integrate_group_backward_euler(NULL, 1u, state, NULL,
      temperature, pressure, 1.0, NULL) != 0);
  return 0;
}

static int convergence_checks(void) {
  chem_backward_euler_parameters parameters = chem_backward_euler_default_parameters();
  double residual[12] = { 0.0 };
  double next[12] = { 0.0 };
  double absolute_tolerance[] = { 1.0e-6, 1.0e-6, 1.0e-6 };
  int converged = 0;
  parameters.small = 1.0e-6;
  REQUIRE(chem_backward_euler_is_converged(&parameters, 4u, 3u, residual, next,
      absolute_tolerance, 1.0e-3, &converged) == 0 && converged);
  residual[1] = 1.0e-5;
  REQUIRE(chem_backward_euler_is_converged(&parameters, 4u, 3u, residual, next,
      absolute_tolerance, 1.0e-3, &converged) == 0 && !converged);
  parameters.small = 1.0e-4;
  REQUIRE(chem_backward_euler_is_converged(&parameters, 4u, 3u, residual, next,
      absolute_tolerance, 1.0e-3, &converged) == 0 && converged);
  residual[11] = 1.0e-3;
  REQUIRE(chem_backward_euler_is_converged(&parameters, 4u, 3u, residual, next,
      absolute_tolerance, 1.0e-3, &converged) == 0 && !converged);
  next[11] = 10.0;
  REQUIRE(chem_backward_euler_is_converged(&parameters, 4u, 3u, residual, next,
      absolute_tolerance, 1.0e-3, &converged) == 0 && converged);
  residual[11] = 1.0e-1;
  REQUIRE(chem_backward_euler_is_converged(&parameters, 4u, 3u, residual, next,
      absolute_tolerance, 1.0e-3, &converged) == 0 && !converged);
  absolute_tolerance[2] = 1.0;
  REQUIRE(chem_backward_euler_is_converged(&parameters, 4u, 3u, residual, next,
      absolute_tolerance, 1.0e-3, &converged) == 0 && converged);
  residual[0] = INFINITY;
  REQUIRE(chem_backward_euler_is_converged(&parameters, 4u, 3u, residual, next,
      absolute_tolerance, 1.0e-3, &converged) == 0 && !converged);
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  if (strcmp(argv[1], "CanCallSolve") == 0 || strcmp(argv[1], "CanBuildBackwardEuler") == 0 ||
      strcmp(argv[1], "CanBuildBackwardEulerOverloadedSolverMethod") == 0)
    return solve_checks();
  if (strcmp(argv[1], "IsConverged") == 0) return convergence_checks();
  if (strcmp(argv[1], "ThrowsMissingSystem") == 0) return missing_mechanism_check();
  return 2;
}
