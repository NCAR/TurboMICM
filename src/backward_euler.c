#include "kinetics.h"
#include "matrix.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

int chem_backward_euler_is_converged(
    const chem_backward_euler_parameters *parameters,
    const double *residuals, const double *new_values,
    const double *absolute_tolerances, double relative_tolerance,
    size_t cell_count, size_t species_count)
{
  size_t cell, species;
  if (!parameters || (cell_count && species_count && (!residuals || !new_values)) ||
      !isfinite(parameters->small) || parameters->small < 0.0 ||
      !isfinite(relative_tolerance) || relative_tolerance < 0.0)
    return 0;
  for (cell = 0; cell < cell_count; ++cell)
    for (species = 0; species < species_count; ++species) {
      const size_t index = cell * species_count + species;
      const double residual = residuals[index];
      const double absolute = absolute_tolerances
                                  ? absolute_tolerances[species] : 0.0;
      if (!isfinite(residual) || !isfinite(new_values[index]))
        return 0;
      if (fabs(residual) > parameters->small &&
          fabs(residual) > absolute &&
          fabs(residual) > relative_tolerance * fabs(new_values[index]))
        return 0;
    }
  return 1;
}

int chem_backward_euler_step(
    const chem_mechanism *mechanism, size_t cell_count,
    const double *rate_constants, double *concentrations, double time_step,
    const chem_backward_euler_parameters *parameters,
    double relative_tolerance, const double *absolute_tolerances,
    chem_solver_stats *stats)
{
  const size_t n = mechanism ? mechanism->species_count : 0;
  const size_t reactions = mechanism ? mechanism->reaction_count : 0;
  const size_t max_iterations = parameters && parameters->max_iterations
                                    ? parameters->max_iterations : 50;
  chem_backward_euler_parameters defaults = {1.0e-6, max_iterations};
  const chem_backward_euler_parameters *p = parameters ? parameters : &defaults;
  double *work;
  double *old_values, *residual, *forcing, *jacobian, *rhs;
  size_t cell, i, j;
  if (!mechanism || !cell_count || !n || !concentrations || !stats ||
      (reactions && !rate_constants) || !(time_step > 0.0) ||
      !isfinite(time_step) || !(p->small >= 0.0) || !isfinite(p->small) ||
      !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
      n > SIZE_MAX / n || cell_count > SIZE_MAX / n)
    return -1;
  if (n > SIZE_MAX / (4 * sizeof(double)) ||
      n * n > SIZE_MAX / sizeof(double))
    return -1;
  work = (double *)malloc((4 * n + n * n) * sizeof(*work));
  if (!work)
    return -1;
  old_values = work;
  residual = old_values + n;
  forcing = residual + n;
  rhs = forcing + n;
  jacobian = rhs + n;
  memset(stats, 0, sizeof(*stats));
  for (cell = 0; cell < cell_count; ++cell) {
    const double *rates = rate_constants ? rate_constants + cell * reactions : NULL;
    double *y = concentrations + cell * n;
    size_t iteration;
    memcpy(old_values, y, n * sizeof(*old_values));
    for (iteration = 0; iteration < max_iterations; ++iteration) {
      if (chem_mechanism_forcing(mechanism, y, rates, forcing) != 0 ||
          chem_mechanism_jacobian_dense(mechanism, y, rates, jacobian) != 0) {
        free(work);
        return -1;
      }
      ++stats->function_calls;
      ++stats->jacobian_updates;
      for (i = 0; i < n; ++i) {
        residual[i] = y[i] - old_values[i] - time_step * forcing[i];
        rhs[i] = -residual[i];
      }
      if (chem_backward_euler_is_converged(p, residual, y,
                                           absolute_tolerances,
                                           relative_tolerance, 1, n)) {
        ++stats->accepted;
        ++stats->internal_steps;
        break;
      }
      for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
          jacobian[i * n + j] = (i == j ? 1.0 : 0.0) -
                                time_step * jacobian[i * n + j];
      if (chem_lu_factor_dense(jacobian, n) != 0 ||
          chem_lu_solve_dense(jacobian, n, rhs) != 0) {
        memcpy(y, old_values, n * sizeof(*y));
        free(work);
        return -1;
      }
      ++stats->decompositions;
      ++stats->solves;
      for (i = 0; i < n; ++i)
        y[i] += rhs[i];
      if (iteration + 1 == max_iterations) {
        memcpy(y, old_values, n * sizeof(*y));
        ++stats->rejected;
        free(work);
        return -1;
      }
    }
  }
  free(work);
  return 0;
}
