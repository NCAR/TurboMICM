#include "chem_specialized.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double parameter_value(const char *name, double initial) {
  if (strstr(name, "effective radius")) return 1.0e-7;
  if (strstr(name, "particle number")) return 2.0e8;
  if (strcmp(name, "my photolysis rate") == 0) return 1.0e-8;
  if (strcmp(name, "my loss rate") == 0) return 1.0e-10;
  if (strcmp(name, "my emission rate") == 0) return 1.0e-14;
  return initial;
}

int main(void) {
  const chem_mechanism *source = chem_generated_mechanism();
  if (!source || !source->species_count || !source->reaction_count) {
    fputs("generated mechanism is empty\n", stderr);
    return 1;
  }
  const size_t cells = 3, n = source->species_count, pcount = source->parameter_count;
  chem_mechanism generic = *source;
  chem_species *species = malloc(n * sizeof(*species));
  double *actual = malloc(cells * n * sizeof(*actual));
  double *expected = malloc(cells * n * sizeof(*expected));
  double *parameters = calloc(cells * (pcount ? pcount : 1u), sizeof(*parameters));
  double temperature[3] = { 287.45, 287.45, 287.45 };
  double pressure[3] = { 101319.9, 101319.9, 101319.9 };
  if (!species || !actual || !expected || !parameters) {
    fputs("allocation failed\n", stderr);
    free(species); free(actual); free(expected); free(parameters);
    return 1;
  }
  memcpy(species, source->species, n * sizeof(*species));
  for (size_t s = 0; s < n; ++s) species[s].absolute_tolerance = 1.0e-15;
  generic.species = species;
  generic.relative_tolerance = 1.0e-6;
  for (size_t c = 0; c < cells; ++c) {
    for (size_t s = 0; s < n; ++s) {
      double value = source->species[s].initial_value;
      actual[c * n + s] = expected[c * n + s] = value;
    }
    for (size_t p = 0; p < pcount; ++p)
      parameters[c * pcount + p] = parameter_value(source->parameters[p].name, source->parameters[p].initial_value);
  }
  chem_statistics specialized_stats = { 0 }, generic_stats = { 0 };
  const int specialized_status = chem_generated_integrate_group_steps(
      cells, actual, parameters, temperature, pressure, 1.0e-3, 1u, 1.0e-6, 1.0e-15, &specialized_stats);
  const int generic_status = chem_integrate_group_steps(
      &generic, cells, expected, parameters, temperature, pressure, 1.0e-3, 1u, &generic_stats);
  if (specialized_status != 0 || generic_status != 0) {
    fprintf(stderr, "solver failed: specialized=%d generic=%d\n", specialized_status, generic_status);
    free(species); free(actual); free(expected); free(parameters);
    return 1;
  }
  for (size_t i = 0; i < cells * n; ++i) {
    const double difference = fabs(actual[i] - expected[i]);
    const double limit = 1.0e-30 + 1.0e-10 * fabs(expected[i]);
    if (!isfinite(actual[i]) || difference > limit) {
      fprintf(stderr, "specialized mismatch at cell %zu species %s: %.17g vs %.17g (error %.3g)\n",
              i / n, source->species[i % n].name, actual[i], expected[i], difference);
      free(species); free(actual); free(expected); free(parameters);
      return 1;
    }
  }
  printf("PASS specialized equivalence: %s (%zu species, %zu reactions, %lu internal steps)\n",
         source->name, n, source->reaction_count, specialized_stats.internal_steps);
  free(species); free(actual); free(expected); free(parameters);
  return 0;
}
