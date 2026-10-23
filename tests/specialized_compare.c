#include "chemistry.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern const chem_mechanism chem_generated_mechanism;

static int close_arrays(const double *a, const double *b, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    const double scale = fmax(fabs(a[i]), fabs(b[i]));
    if (fabs(a[i] - b[i]) > 1.0e-13 * scale + 1.0e-300) {
      fprintf(stderr, "entry %zu differs: %.17g versus %.17g\n", i, a[i], b[i]);
      return 0;
    }
  }
  return 1;
}

int main(void) {
  const chem_mechanism *specialized = &chem_generated_mechanism;
  chem_mechanism generic = *specialized;
  generic.specialized_ops = NULL;
  const size_t n = specialized->species_count;
  const size_t nr = specialized->reaction_count;
  const size_t np = specialized->custom_parameter_count;
  double *y = calloc(n, sizeof(*y));
  double *parameters = calloc(np ? np : 1, sizeof(*parameters));
  double *rates_a = calloc(nr, sizeof(*rates_a));
  double *rates_b = calloc(nr, sizeof(*rates_b));
  double *forcing_a = calloc(n, sizeof(*forcing_a));
  double *forcing_b = calloc(n, sizeof(*forcing_b));
  double *jacobian_a = calloc(n*n, sizeof(*jacobian_a));
  double *jacobian_b = calloc(n*n, sizeof(*jacobian_b));
  double *state_a = calloc(n, sizeof(*state_a));
  double *state_b = calloc(n, sizeof(*state_b));
  if (!y || !parameters || !rates_a || !rates_b || !forcing_a || !forcing_b ||
      !jacobian_a || !jacobian_b || !state_a || !state_b) return 2;
  for (size_t i = 0; i < n; ++i)
    y[i] = 1.0e-10 * (1.0 + 0.03125 * (double)i);
  if (np) {
    for (size_t i = 0; i < np; ++i)
      parameters[i] = specialized->default_custom_parameters[i];
  }
  const double temperature = 287.45, pressure = 101319.9;
  const double density = pressure / (8.31446261815324 * temperature);
  if (chem_rate_constants(specialized, temperature, pressure, density,
                          parameters, rates_a) ||
      chem_rate_constants(&generic, temperature, pressure, density,
                          parameters, rates_b) ||
      !close_arrays(rates_a, rates_b, nr) ||
      specialized->specialized_ops->forcing_with_constants(
          specialized, y, rates_a, density, forcing_a) ||
      chem_forcing_with_constants(&generic, y, rates_b, density, forcing_b) ||
      !close_arrays(forcing_a, forcing_b, n) ||
      specialized->specialized_ops->negative_jacobian_with_pattern(
          specialized, y, rates_a, NULL, NULL, jacobian_a) ||
      chem_negative_jacobian_with_constants(&generic, y, rates_b, jacobian_b) ||
      !close_arrays(jacobian_a, jacobian_b, n*n)) return 1;
  for (size_t i = 0; i < n; ++i)
    state_a[i] = state_b[i] = 1.0e-10 * (double)(i + 1);
  chem_solver_workspace workspace_a = {0}, workspace_b = {0};
  chem_solver_options options = chem_solver_default_options();
  chem_solver_stats stats_a, stats_b;
  if (chem_solver_workspace_init(&workspace_a, specialized) ||
      chem_solver_workspace_init(&workspace_b, &generic) ||
      chem_rosenbrock3_solve_with_workspace(
          specialized, &workspace_a, state_a, temperature, pressure, density,
          parameters, 1.0e-6, &options, &stats_a) ||
      chem_rosenbrock3_solve_with_workspace(
          &generic, &workspace_b, state_b, temperature, pressure, density,
          parameters, 1.0e-6, &options, &stats_b) ||
      !close_arrays(state_a, state_b, n)) return 1;
  chem_solver_workspace_destroy(&workspace_a);
  chem_solver_workspace_destroy(&workspace_b);
  free(y); free(parameters); free(rates_a); free(rates_b);
  free(forcing_a); free(forcing_b); free(jacobian_a); free(jacobian_b);
  free(state_a); free(state_b);
  return 0;
}
