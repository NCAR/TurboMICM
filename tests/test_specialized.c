#include "kinetics.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

extern const chem_mechanism chem_model;

static int close_enough(double actual, double reference)
{
  return isfinite(actual) && isfinite(reference) &&
         fabs(actual - reference) <= 1.0e-30 + 1.0e-10 * fabs(reference);
}

static void initialize_custom(const chem_mechanism *mechanism, size_t cells,
                              double *custom)
{
  size_t cell, parameter;
  for (cell = 0; cell < cells; ++cell) {
    for (parameter = 0; parameter < mechanism->custom_count; ++parameter) {
      const char *name = mechanism->custom_parameters[parameter].name;
      double value = mechanism->custom_parameters[parameter].default_value;
      if (strstr(name, "effective radius"))
        value = 1.0e-6;
      else if (strstr(name, "particle number"))
        value = 1.0e8;
      else if (name[0] == 'j' && !(value > 0.0))
        value = 1.0e-5;
      else if (!(value > 0.0))
        value = 1.0e-4;
      custom[cell * mechanism->custom_count + parameter] = value;
    }
  }
}

int main(void)
{
  enum { cells = 4 };
  chem_mechanism generic = chem_model;
  chem_rosenbrock_parameters parameters;
  chem_solver *specialized_solver = NULL;
  chem_solver *generic_solver = NULL;
  chem_solver_stats stats;
  double temperature[cells], pressure[cells];
  double *custom = NULL, *specialized_state = NULL, *generic_state = NULL;
  size_t cell, species, parameter, step;
  int status = 1;

#ifdef _OPENMP
  omp_set_num_threads(1);
#endif
  if (!chem_model.species_count || !chem_model.species ||
      !chem_model.reactions || !chem_model.kernels) {
    fprintf(stderr, "generated model has no specialized kernels\n");
    return 1;
  }
  generic.kernels = NULL;
  custom = (double *)malloc((chem_model.custom_count != 0
                                 ? cells * chem_model.custom_count : 1) * sizeof(*custom));
  specialized_state = (double *)malloc(cells * chem_model.species_count * sizeof(*specialized_state));
  generic_state = (double *)malloc(cells * chem_model.species_count * sizeof(*generic_state));
  if (!custom || !specialized_state || !generic_state)
    goto done;
  initialize_custom(&chem_model, cells, custom);
  for (cell = 0; cell < cells; ++cell) {
    temperature[cell] = 278.0 + 8.0 * (double)cell;
    pressure[cell] = 97000.0 + 1500.0 * (double)cell;
    for (species = 0; species < chem_model.species_count; ++species) {
      double value = chem_model.species[species].initial_value;
      if (value != 0.0)
        value *= 1.0 + 0.005 * (double)cell;
      specialized_state[cell * chem_model.species_count + species] = value;
      generic_state[cell * chem_model.species_count + species] = value;
    }
  }
  chem_rosenbrock_default(&parameters);
  specialized_solver = chem_solver_create(&chem_model, cells, temperature, pressure,
                                          custom, &parameters, 2);
  generic_solver = chem_solver_create(&generic, cells, temperature, pressure,
                                      custom, &parameters, 2);
  if (!specialized_solver || !generic_solver) {
    fprintf(stderr, "could not create test solvers for %s\n", chem_model.name);
    goto done;
  }
  for (step = 0; step < 3; ++step) {
    for (cell = 0; cell < cells; ++cell) {
      for (parameter = 0; parameter < chem_model.custom_count; ++parameter) {
        const char *name = chem_model.custom_parameters[parameter].name;
        if (!strstr(name, "effective radius") && !strstr(name, "particle number"))
          custom[cell * chem_model.custom_count + parameter] *= 0.75;
      }
    }
    if (chem_solver_update_custom_rates(specialized_solver, temperature, pressure,
                                        custom) != 0 ||
        chem_solver_update_custom_rates(generic_solver, temperature, pressure,
                                        custom) != 0 ||
        chem_solver_step(specialized_solver, specialized_state, 30.0, &stats) != 0 ||
        chem_solver_step(generic_solver, generic_state, 30.0, &stats) != 0) {
      fprintf(stderr, "solver step failed for %s at interval %zu\n",
              chem_model.name, step + 1);
      goto done;
    }
  }
  for (cell = 0; cell < cells; ++cell) {
    for (species = 0; species < chem_model.species_count; ++species) {
      const size_t index = cell * chem_model.species_count + species;
      if (!close_enough(specialized_state[index], generic_state[index])) {
        fprintf(stderr, "%s: cell %zu species %s differs: %.17g vs %.17g\n",
                chem_model.name, cell, chem_model.species[species].name,
                specialized_state[index], generic_state[index]);
        goto done;
      }
    }
  }
  printf("PASS specialized solver agrees with generic solver: %s\n", chem_model.name);
  status = 0;

done:
  chem_solver_destroy(specialized_solver);
  chem_solver_destroy(generic_solver);
  free(custom);
  free(specialized_state);
  free(generic_state);
  return status;
}
