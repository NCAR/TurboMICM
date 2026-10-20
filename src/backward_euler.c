#include "chem.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static double *allocate_values(size_t count, size_t width) {
  if (width && count > SIZE_MAX / width / sizeof(double)) return NULL;
  const size_t total = count * width;
  return calloc(total ? total : 1u, sizeof(double));
}

chem_backward_euler_parameters chem_backward_euler_default_parameters(void) {
  chem_backward_euler_parameters parameters = { 0 };
  parameters.small = 1.0e-40;
  parameters.max_iterations = 11u;
  parameters.reduction_count = 5u;
  parameters.time_step_reductions[0] = 0.5;
  parameters.time_step_reductions[1] = 0.5;
  parameters.time_step_reductions[2] = 0.5;
  parameters.time_step_reductions[3] = 0.5;
  parameters.time_step_reductions[4] = 0.1;
  return parameters;
}

static int valid_backward_euler_parameters(const chem_backward_euler_parameters *parameters) {
  if (!parameters || !isfinite(parameters->small) || parameters->small < 0.0 ||
      !isfinite(parameters->h_start) || parameters->h_start < 0.0 ||
      !parameters->max_iterations || parameters->reduction_count > 5u) return 0;
  for (size_t i = 0; i < parameters->reduction_count; ++i)
    if (!isfinite(parameters->time_step_reductions[i]) ||
        parameters->time_step_reductions[i] <= 0.0 || parameters->time_step_reductions[i] >= 1.0)
      return 0;
  return 1;
}

int chem_backward_euler_is_converged(
    const chem_backward_euler_parameters *solver_parameters,
    size_t cell_count, size_t species_count,
    const double *residual, const double *next_state,
    const double *absolute_tolerances, double relative_tolerance,
    int *is_converged) {
  if (!valid_backward_euler_parameters(solver_parameters) || !cell_count || !species_count ||
      !residual || !next_state || !absolute_tolerances || !is_converged ||
      !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
      cell_count > SIZE_MAX / species_count) return -1;
  *is_converged = 1;
  for (size_t cell = 0; cell < cell_count; ++cell) {
    for (size_t species = 0; species < species_count; ++species) {
      const size_t i = cell * species_count + species;
      const double correction = fabs(residual[i]);
      const double value = fabs(next_state[i]);
      if (!isfinite(correction) || !isfinite(value) || !isfinite(absolute_tolerances[species]) ||
          absolute_tolerances[species] < 0.0 ||
          (correction > solver_parameters->small &&
           correction > absolute_tolerances[species] &&
           correction > relative_tolerance * value)) {
        *is_converged = 0;
        return 0;
      }
    }
  }
  return 0;
}

static int calculate_rates(
    const chem_mechanism *mechanism, size_t cell_count, const double *parameters,
    const double *temperature, const double *pressure, double *rates) {
  const size_t parameter_count = mechanism->parameter_count;
  for (size_t cell = 0; cell < cell_count; ++cell) {
    const double *cell_parameters = parameter_count
        ? parameters + cell * parameter_count : NULL;
    for (size_t reaction = 0; reaction < mechanism->reaction_count; ++reaction) {
      if (chem_evaluate_rate(mechanism, reaction, cell_parameters, temperature[cell], pressure[cell],
              &rates[cell * mechanism->reaction_count + reaction]) != 0)
        return -1;
    }
  }
  return 0;
}

int chem_integrate_group_backward_euler_with_parameters(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    const double *restrict parameters_external,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    const chem_backward_euler_parameters *solver_parameters,
    chem_statistics *statistics) {
  if (!mechanism || !cell_count || !y || !temperature || !pressure || !isfinite(interval) || interval <= 0.0 ||
      !mechanism->species_count || !mechanism->reaction_count ||
      !mechanism->species || !mechanism->reactions ||
      (mechanism->parameter_count && !parameters_external && !mechanism->parameters) ||
      !isfinite(mechanism->relative_tolerance) || mechanism->relative_tolerance < 0.0 ||
      !valid_backward_euler_parameters(solver_parameters) ||
      cell_count > SIZE_MAX / mechanism->species_count ||
      mechanism->species_count > SIZE_MAX / mechanism->species_count ||
      cell_count > SIZE_MAX / mechanism->reaction_count ||
      (mechanism->parameter_count && cell_count > SIZE_MAX / mechanism->parameter_count)) return -1;

  const size_t n = mechanism->species_count;
  const size_t cells = cell_count;
  const size_t state_length = cells * n;
  const size_t matrix_length = n * n;
  const size_t rate_length = cells * mechanism->reaction_count;
  double *initial = allocate_values(state_length, 1u);
  double *next = allocate_values(state_length, 1u);
  double *forcing = allocate_values(state_length, 1u);
  double *residual = allocate_values(state_length, 1u);
  double *rates = allocate_values(rate_length, 1u);
  double *jacobian = allocate_values(cells, matrix_length);
  double *owned_parameters = NULL;
  double *absolute_tolerances = allocate_values(n, 1u);
  if (mechanism->parameter_count && !parameters_external)
    owned_parameters = allocate_values(cells, mechanism->parameter_count);
  if (!initial || !next || !forcing || !residual || !rates || !jacobian || !absolute_tolerances ||
      (mechanism->parameter_count && !parameters_external && !owned_parameters)) {
    free(initial); free(next); free(forcing); free(residual); free(rates); free(jacobian);
    free(owned_parameters); free(absolute_tolerances);
    return -2;
  }
  const double *parameters = parameters_external ? parameters_external : owned_parameters;
  if (owned_parameters) {
    for (size_t cell = 0; cell < cells; ++cell)
      for (size_t p = 0; p < mechanism->parameter_count; ++p)
        owned_parameters[cell * mechanism->parameter_count + p] = mechanism->parameters[p].initial_value;
  }
  for (size_t species = 0; species < n; ++species) {
    if (!isfinite(mechanism->species[species].absolute_tolerance) ||
        mechanism->species[species].absolute_tolerance < 0.0) {
      free(initial); free(next); free(forcing); free(residual); free(rates); free(jacobian);
      free(owned_parameters); free(absolute_tolerances);
      return -1;
    }
    absolute_tolerances[species] = mechanism->species[species].absolute_tolerance;
  }
  for (size_t cell = 0; cell < cells; ++cell) {
    if (!isfinite(temperature[cell]) || temperature[cell] <= 0.0 ||
        !isfinite(pressure[cell]) || pressure[cell] < 0.0) {
      free(initial); free(next); free(forcing); free(residual); free(rates); free(jacobian);
      free(owned_parameters); free(absolute_tolerances);
      return -1;
    }
  }
  if (chem_evaluate_forcing(mechanism, cells, y, rates, forcing) != 0) {
    free(initial); free(next); free(forcing); free(residual); free(rates); free(jacobian);
    free(owned_parameters); free(absolute_tolerances);
    return -1;
  }
  if (calculate_rates(mechanism, cells, parameters, temperature, pressure, rates) != 0) {
    free(initial); free(next); free(forcing); free(residual); free(rates); free(jacobian);
    free(owned_parameters); free(absolute_tolerances);
    return -1;
  }

  memcpy(initial, y, state_length * sizeof(double));
  memcpy(next, y, state_length * sizeof(double));
  double step = solver_parameters->h_start == 0.0 ? interval : solver_parameters->h_start;
  double present_time = 0.0;
  size_t successful_steps = 0u, failed_steps = 0u;
  int status = 0;

  while (present_time < interval) {
    int converged = 0;
    size_t iterations = 0u;
    do {
      memset(forcing, 0, state_length * sizeof(double));
      memset(jacobian, 0, cells * matrix_length * sizeof(double));
      if (chem_evaluate_forcing(mechanism, cells, next, rates, forcing) != 0 ||
          chem_accumulate_negative_jacobian(mechanism, cells, next, rates, NULL, jacobian) != 0) {
        status = -1;
        goto cleanup;
      }
      if (statistics) {
        ++statistics->internal_steps;
        ++statistics->function_calls;
        ++statistics->jacobian_updates;
        ++statistics->decompositions;
        ++statistics->solves;
      }

      for (size_t cell = 0; cell < cells; ++cell) {
        double *cell_jacobian = jacobian + cell * matrix_length;
        double *cell_residual = residual + cell * n;
        const double *cell_forcing = forcing + cell * n;
        const double *cell_next = next + cell * n;
        const double *cell_initial = initial + cell * n;
        for (size_t species = 0; species < n; ++species) {
          cell_jacobian[species * n + species] += 1.0 / step;
          cell_residual[species] = cell_forcing[species] - (cell_next[species] - cell_initial[species]) / step;
        }
        if (chem_lu_factor_dense_in_place(cell_jacobian, n, CHEM_LU_DOOLITTLE) != 0 ||
            chem_lu_solve_dense_in_place(cell_jacobian, n, cell_residual) != 0) {
          status = -1;
          goto cleanup;
        }
        double *cell_state = next + cell * n;
        for (size_t species = 0; species < n; ++species) {
          const double updated = cell_state[species] + cell_residual[species];
          cell_state[species] = updated < 0.0 ? 0.0 : updated;
        }
      }

      if (iterations++ != 0u && chem_backward_euler_is_converged(
              solver_parameters, cells, n, residual, next, absolute_tolerances,
              mechanism->relative_tolerance, &converged) != 0) {
        status = -1;
        goto cleanup;
      }
    } while (!converged && iterations < solver_parameters->max_iterations);

    if (!converged) {
      if (statistics) ++statistics->rejected;
      successful_steps = 0u;
      if (failed_steps >= solver_parameters->reduction_count) {
        if (!isfinite(step) || step <= 0.0) {
          status = -1;
          goto cleanup;
        }
        present_time += step;
        for (size_t i = 0; i < state_length; ++i)
          if (!isfinite(next[i])) {
            status = -1;
            goto cleanup;
          }
      } else {
        memcpy(next, initial, state_length * sizeof(double));
        step *= solver_parameters->time_step_reductions[failed_steps++];
      }
    } else {
      if (statistics) ++statistics->accepted;
      present_time += step;
      memcpy(initial, next, state_length * sizeof(double));
      if (++successful_steps >= 2u) {
        successful_steps = 0u;
        step *= 2.0;
      }
    }
    if (present_time < interval && step > interval - present_time) step = interval - present_time;
  }

  memcpy(y, next, state_length * sizeof(double));

cleanup:
  free(initial); free(next); free(forcing); free(residual); free(rates); free(jacobian);
  free(owned_parameters); free(absolute_tolerances);
  return status;
}

int chem_integrate_group_backward_euler(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    const double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    chem_statistics *statistics) {
  const chem_backward_euler_parameters solver_parameters = chem_backward_euler_default_parameters();
  return chem_integrate_group_backward_euler_with_parameters(
      mechanism, cell_count, y, parameters, temperature, pressure, interval, &solver_parameters, statistics);
}
