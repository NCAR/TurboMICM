#include "chem.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHEM_GAS_CONSTANT (1.380649e-23 * 6.02214076e23)
#define CHEM_PI 3.14159265358979323846264338327950288

typedef struct {
  size_t reaction;
  size_t independent;
  size_t dependent_offset;
  size_t dependent_count;
} chem_jacobian_job;

typedef struct {
  size_t species_count;
  size_t reaction_count;
  size_t job_count;
  chem_jacobian_job *jobs;
  size_t *dependent_species;
  unsigned char *pattern;
  size_t *position;
  size_t nonzero_count;
  size_t *lower;
  size_t *upper;
  size_t *lower_count;
  size_t *upper_count;
} chem_symbolic;

static double chem_falloff(
    const chem_rate *rate,
    double temperature,
    double air_density,
    int multiply_density) {
  const double k0 = rate->p[5] * exp(rate->p[7] / temperature) * pow(temperature / 300.0, rate->p[6]);
  const double kinf = rate->p[8] * exp(rate->p[10] / temperature) * pow(temperature / 300.0, rate->p[9]);
  const double ratio = k0 * air_density / kinf;
  const double numerator = multiply_density ? air_density : 1.0;
  return k0 * numerator / (1.0 + ratio) *
         pow(rate->p[11], rate->p[12] / (rate->p[12] + pow(log10(ratio), 2.0)));
}

static double chem_rate_value(
    const chem_mechanism *mechanism,
    const chem_reaction *reaction,
    const double *parameters,
    size_t cell,
    size_t cell_count,
    double temperature,
    double pressure,
    double air_density) {
  const chem_rate *rate = &reaction->rate;
  double value = 0.0;
  switch (rate->kind) {
    case CHEM_RATE_ARRHENIUS:
      value = rate->p[0] * exp(rate->p[2] / temperature) * pow(temperature / rate->p[3], rate->p[1]) *
              (1.0 + rate->p[4] * pressure);
      break;
    case CHEM_RATE_TROE:
      value = chem_falloff(rate, temperature, air_density, 1);
      break;
    case CHEM_RATE_TERNARY:
      value = chem_falloff(rate, temperature, air_density, 0);
      break;
    case CHEM_RATE_TERNARY_JPL19: {
      const double k0 = rate->p[5] * exp(rate->p[7] / temperature) * pow(temperature / rate->p[18], rate->p[6]);
      const double kinf = rate->p[8] * exp(rate->p[10] / temperature) * pow(temperature / rate->p[19], rate->p[9]);
      const double kint = rate->p[20] * exp(rate->p[22] / temperature) * pow(temperature / rate->p[23], rate->p[21]);
      const double ratio = k0 * air_density / kinf;
      const double kf = k0 * air_density / (1.0 + ratio) *
                        pow(rate->p[11], rate->p[12] / (rate->p[12] + pow(log10(ratio), 2.0)));
      value = kf + kint * (1.0 - kf / kinf);
      break;
    }
    case CHEM_RATE_BRANCHED: {
      const double k0 = 2.0e-22 * 6.02214076e23 * 1.0e-6 * exp(rate->p[16]);
      const double air_ref = 2.45e19 / 6.02214076e23 * 1.0e6;
      const double a_ref = k0 * air_ref;
      const double b_ref = 0.43 * pow(293.0 / 298.0, -8.0);
      const double aval_ref = a_ref / (1.0 + a_ref / b_ref) *
                              pow(0.41, 1.0 / (1.0 + pow(log10(a_ref / b_ref), 2.0)));
      const double z = aval_ref * (1.0 - rate->p[15]) / rate->p[15];
      const double a = k0 * air_density;
      const double b = 0.43 * pow(temperature / 298.0, -8.0);
      const double aval = a / (1.0 + a / b) * pow(0.41, 1.0 / (1.0 + pow(log10(a / b), 2.0)));
      const double pre = rate->p[13] * exp(-rate->p[14] / temperature);
      value = rate->branch ? pre * (z / (z + aval)) : pre * (aval / (aval + z));
      break;
    }
    case CHEM_RATE_TUNNELING:
      value = rate->p[0] * exp(-rate->p[1] / temperature + rate->p[2] / (temperature * temperature * temperature));
      break;
    case CHEM_RATE_TAYLOR: {
      double polynomial = 0.0;
      double power = 1.0;
      for (size_t i = 0; i < rate->coefficient_count; ++i, power *= temperature)
        polynomial += rate->coefficients[i] * power;
      value = polynomial * rate->p[0] * exp(rate->p[2] / temperature) *
              pow(temperature / rate->p[3], rate->p[1]) * (1.0 + rate->p[4] * pressure);
      break;
    }
    case CHEM_RATE_REVERSIBLE:
      value = rate->p[0] * exp(rate->p[2] / temperature) * rate->p[17];
      break;
    case CHEM_RATE_USER:
      value = parameters[rate->custom_parameter * cell_count + cell] * rate->p[1];
      break;
    case CHEM_RATE_SURFACE: {
      const chem_species *species = &mechanism->species[rate->surface_species];
      const double radius = parameters[rate->custom_parameter * cell_count + cell];
      const double number = parameters[(rate->custom_parameter + 1u) * cell_count + cell];
      const double speed = sqrt(8.0 * CHEM_GAS_CONSTANT * temperature /
                                (CHEM_PI * species->molecular_weight));
      value = 4.0 * number * CHEM_PI * radius * radius /
              (radius / species->diffusion_coefficient + 4.0 / (speed * rate->p[0]));
      break;
    }
    default:
      return NAN;
  }
  for (unsigned i = 0; i < reaction->rate.air_density_power; ++i)
    value *= air_density;
  return value;
}

static void chem_symbolic_destroy(chem_symbolic *symbolic) {
  free(symbolic->jobs);
  free(symbolic->dependent_species);
  free(symbolic->pattern);
  free(symbolic->position);
  free(symbolic->lower);
  free(symbolic->upper);
  free(symbolic->lower_count);
  free(symbolic->upper_count);
  memset(symbolic, 0, sizeof(*symbolic));
}

static int chem_symbolic_build(const chem_mechanism *mechanism, chem_symbolic *symbolic) {
  const size_t n = mechanism->species_count;
  size_t jobs = 0, dependent_count = 0;
  memset(symbolic, 0, sizeof(*symbolic));
  symbolic->species_count = n;
  symbolic->reaction_count = mechanism->reaction_count;
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    jobs += mechanism->reactions[r].reactant_count;
    dependent_count += mechanism->reactions[r].reactant_count * mechanism->reactions[r].reactant_count;
  }
  symbolic->jobs = calloc(jobs ? jobs : 1, sizeof(*symbolic->jobs));
  symbolic->dependent_species = malloc((dependent_count ? dependent_count : 1) * sizeof(size_t));
  symbolic->pattern = calloc(n * n, sizeof(unsigned char));
  symbolic->lower = malloc(n * n * sizeof(size_t));
  symbolic->upper = malloc(n * n * sizeof(size_t));
  symbolic->lower_count = calloc(n, sizeof(size_t));
  symbolic->upper_count = calloc(n, sizeof(size_t));
  if (!symbolic->jobs || !symbolic->dependent_species || !symbolic->pattern || !symbolic->lower ||
      !symbolic->upper || !symbolic->lower_count || !symbolic->upper_count) {
    chem_symbolic_destroy(symbolic);
    return -1;
  }
  symbolic->job_count = jobs;
  for (size_t i = 0; i < n; ++i)
    symbolic->pattern[i * n + i] = 1;

  size_t job_pos = 0, dep_pos = 0;
  for (size_t independent = 0; independent < n; ++independent) {
    for (size_t r = 0; r < mechanism->reaction_count; ++r) {
      const chem_reaction *reaction = &mechanism->reactions[r];
      for (size_t q = 0; q < reaction->reactant_count; ++q) {
        const size_t column = mechanism->reactant_terms[reaction->reactant_offset + q].species;
        if (column != independent)
          continue;
        size_t found = 0;
        chem_jacobian_job *job = &symbolic->jobs[job_pos++];
        job->reaction = r;
        job->independent = independent;
        job->dependent_offset = dep_pos;
        for (size_t d = 0; d < reaction->reactant_count; ++d) {
          const size_t row = mechanism->reactant_terms[reaction->reactant_offset + d].species;
          if (row == independent && !found) {
            found = 1;
            continue;
          }
          symbolic->dependent_species[dep_pos++] = row;
          ++job->dependent_count;
          symbolic->pattern[row * n + independent] = 1;
        }
        for (size_t d = 0; d < reaction->product_count; ++d) {
          const size_t row = mechanism->product_terms[reaction->product_offset + d].species;
          symbolic->pattern[row * n + independent] = 1;
        }
      }
    }
  }

  /* This is the symbolic fill used by the in-place Mozart LU: no row pivoting. */
  for (size_t k = 0; k < n; ++k) {
    for (size_t i = k + 1; i < n; ++i) {
      if (!symbolic->pattern[i * n + k])
        continue;
      for (size_t j = k + 1; j < n; ++j)
        if (symbolic->pattern[k * n + j])
          symbolic->pattern[i * n + j] = 1;
    }
  }
  symbolic->position = malloc(n * n * sizeof(size_t));
  if (!symbolic->position) {
    chem_symbolic_destroy(symbolic);
    return -1;
  }
  for (size_t i = 0; i < n * n; ++i)
    symbolic->position[i] = SIZE_MAX;
  for (size_t i = 0; i < n * n; ++i)
    if (symbolic->pattern[i])
      symbolic->position[i] = symbolic->nonzero_count++;
  for (size_t k = 0; k < n; ++k) {
    for (size_t i = k + 1; i < n; ++i)
      if (symbolic->pattern[i * n + k])
        symbolic->lower[k * n + symbolic->lower_count[k]++] = i;
    for (size_t j = k + 1; j < n; ++j)
      if (symbolic->pattern[k * n + j])
        symbolic->upper[k * n + symbolic->upper_count[k]++] = j;
  }
  return 0;
}

static void chem_calculate_rates(
    const chem_mechanism *mechanism,
    size_t cell_count,
    const double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double *restrict rates) {
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    double *rate_row = rates + r * cell_count;
    for (size_t c = 0; c < cell_count; ++c) {
      const double air_density = pressure[c] / (CHEM_GAS_CONSTANT * temperature[c]);
      rate_row[c] = chem_rate_value(mechanism, reaction, parameters, c, cell_count,
                                    temperature[c], pressure[c], air_density);
    }
  }
}

static void chem_add_forcing(
    const chem_mechanism *mechanism,
    size_t cell_count,
    const double *restrict y,
    const double *restrict rates,
    double *restrict forcing,
    double *restrict flux) {
  memset(forcing, 0, mechanism->species_count * cell_count * sizeof(double));
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    const double *rate_row = rates + r * cell_count;
    for (size_t c = 0; c < cell_count; ++c) {
      double rate = rate_row[c];
      for (size_t q = 0; q < reaction->reactant_count; ++q) {
        const size_t species = mechanism->reactant_terms[reaction->reactant_offset + q].species;
        rate *= y[species * cell_count + c];
      }
      flux[c] = rate;
    }
    for (size_t q = 0; q < reaction->reactant_count; ++q) {
      const size_t species = mechanism->reactant_terms[reaction->reactant_offset + q].species;
      double *forcing_row = forcing + species * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        forcing_row[c] -= flux[c];
    }
    for (size_t q = 0; q < reaction->product_count; ++q) {
      const chem_term *product = &mechanism->product_terms[reaction->product_offset + q];
      double *forcing_row = forcing + product->species * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        forcing_row[c] += product->coefficient * flux[c];
    }
  }
}

static void chem_build_jacobian(
    const chem_mechanism *mechanism,
    const chem_symbolic *symbolic,
    size_t cell_count,
    const double *restrict y,
    const double *restrict rates,
    double *restrict jacobian,
    double *restrict scratch) {
  const size_t n = mechanism->species_count;
  memset(jacobian, 0, cell_count * symbolic->nonzero_count * sizeof(double));
  for (size_t j = 0; j < symbolic->job_count; ++j) {
    const chem_jacobian_job *job = &symbolic->jobs[j];
    const chem_reaction *reaction = &mechanism->reactions[job->reaction];
    const double *rate_row = rates + job->reaction * cell_count;
    for (size_t c = 0; c < cell_count; ++c) {
      double derivative = rate_row[c];
      for (size_t d = 0; d < job->dependent_count; ++d) {
        const size_t species = symbolic->dependent_species[job->dependent_offset + d];
        derivative *= y[species * cell_count + c];
      }
      scratch[c] = derivative;
    }
    for (size_t d = 0; d < job->dependent_count; ++d) {
      const size_t row = symbolic->dependent_species[job->dependent_offset + d];
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        jacobian[symbolic->position[row * n + job->independent] * cell_count + c] += scratch[c];
    }
    #pragma omp simd
    for (size_t c = 0; c < cell_count; ++c)
      jacobian[symbolic->position[job->independent * n + job->independent] * cell_count + c] += scratch[c];
    for (size_t q = 0; q < reaction->product_count; ++q) {
      const chem_term *product = &mechanism->product_terms[reaction->product_offset + q];
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        jacobian[symbolic->position[product->species * n + job->independent] * cell_count + c] -=
            product->coefficient * scratch[c];
    }
  }
}

static int chem_factor(
    const chem_symbolic *symbolic,
    size_t cell_count,
    double *matrix,
    double *restrict inverse,
    size_t *bad_pivot,
    size_t *bad_cell) {
  const size_t n = symbolic->species_count;
  for (size_t k = 0; k < n; ++k) {
    const size_t diagonal_position = symbolic->position[k * n + k];
    double *diagonal = matrix + diagonal_position * cell_count;
    for (size_t c = 0; c < cell_count; ++c) {
      if (diagonal[c] == 0.0 || !isfinite(diagonal[c])) {
        *bad_pivot = k;
        *bad_cell = c;
        return -1;
      }
      inverse[c] = 1.0 / diagonal[c];
    }
    const size_t *lower = symbolic->lower + k * n;
    const size_t *upper = symbolic->upper + k * n;
    for (size_t i = 0; i < symbolic->lower_count[k]; ++i) {
      double *lower_entry = matrix + symbolic->position[lower[i] * n + k] * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        lower_entry[c] *= inverse[c];
    }
    for (size_t j = 0; j < symbolic->upper_count[k]; ++j) {
      const size_t col = upper[j];
      const double *pivot_row_value = matrix + symbolic->position[k * n + col] * cell_count;
      for (size_t i = 0; i < symbolic->lower_count[k]; ++i) {
        const size_t row = lower[i];
        if (symbolic->pattern[row * n + col]) {
          double *entry = matrix + symbolic->position[row * n + col] * cell_count;
          const double *lower_entry = matrix + symbolic->position[row * n + k] * cell_count;
          #pragma omp simd
          for (size_t c = 0; c < cell_count; ++c)
            entry[c] -= lower_entry[c] * pivot_row_value[c];
        }
      }
    }
  }
  return 0;
}

static void chem_solve(const chem_symbolic *symbolic, size_t cell_count, const double *matrix, double *rhs) {
  const size_t n = symbolic->species_count;
  for (size_t i = 0; i < n; ++i) {
    double *rhs_i = rhs + i * cell_count;
    for (size_t j = 0; j < i; ++j) {
      if (symbolic->pattern[i * n + j]) {
        const double *entry = matrix + symbolic->position[i * n + j] * cell_count;
        const double *rhs_j = rhs + j * cell_count;
        #pragma omp simd
        for (size_t c = 0; c < cell_count; ++c)
          rhs_i[c] -= entry[c] * rhs_j[c];
      }
    }
  }
  for (size_t ii = n; ii-- > 0;) {
    double *rhs_i = rhs + ii * cell_count;
    for (size_t j = ii + 1; j < n; ++j) {
      if (symbolic->pattern[ii * n + j]) {
        const double *entry = matrix + symbolic->position[ii * n + j] * cell_count;
        const double *rhs_j = rhs + j * cell_count;
        #pragma omp simd
        for (size_t c = 0; c < cell_count; ++c)
          rhs_i[c] -= entry[c] * rhs_j[c];
      }
    }
    const double *diagonal = matrix + symbolic->position[ii * n + ii] * cell_count;
    #pragma omp simd
    for (size_t c = 0; c < cell_count; ++c)
      rhs_i[c] /= diagonal[c];
  }
}

static double chem_normalized_error(
    const chem_mechanism *mechanism,
    size_t cell_count,
    const double *restrict y,
    const double *restrict ynew,
    const double *restrict error) {
  const size_t n = mechanism->species_count;
  double sum = 0.0;
  for (size_t species = 0; species < n; ++species) {
    const chem_species *item = &mechanism->species[species];
    const double atol = item->absolute_tolerance;
    const double *yrow = y + species * cell_count;
    const double *newrow = ynew + species * cell_count;
    const double *erow = error + species * cell_count;
    for (size_t c = 0; c < cell_count; ++c) {
      const double ymax = fmax(fabs(yrow[c]), fabs(newrow[c]));
      const double scaled = erow[c] / (atol + mechanism->relative_tolerance * ymax);
      sum += scaled * scaled;
    }
  }
  const double result = sqrt(sum / (double)(cell_count * mechanism->species_count));
  return fmax(result, 1.0e-10);
}

static int chem_allocate(size_t count, size_t width, double **out) {
  if (width && count > SIZE_MAX / width / sizeof(double))
    return -1;
  const size_t total = count * width;
  *out = calloc(total ? total : 1, sizeof(double));
  return *out ? 0 : -1;
}

int chem_integrate_group_steps(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y_external,
    double *restrict parameters_external,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    size_t interval_count,
    chem_statistics *statistics) {
  if (!mechanism || !cell_count || !interval_count || !y_external || !temperature || !pressure || interval <= 0.0 ||
      !mechanism->species_count || !mechanism->reaction_count)
    return -1;
  const size_t n = mechanism->species_count;
  const size_t nr = mechanism->reaction_count;
  chem_symbolic symbolic;
  if (chem_symbolic_build(mechanism, &symbolic) != 0)
    return -1;

  double *y = NULL, *ynew = NULL, *initial_forcing = NULL, *error = NULL;
  double *k[3] = { NULL, NULL, NULL };
  double *parameters = NULL, *rates = NULL, *jacobian = NULL, *scratch = NULL;
  int failed = 0;
  if (chem_allocate(n, cell_count, &y) || chem_allocate(n, cell_count, &ynew) ||
      chem_allocate(n, cell_count, &initial_forcing) || chem_allocate(n, cell_count, &error) ||
      chem_allocate(n, cell_count, &k[0]) || chem_allocate(n, cell_count, &k[1]) ||
      chem_allocate(n, cell_count, &k[2]) || chem_allocate(mechanism->parameter_count, cell_count, &parameters) ||
      chem_allocate(nr, cell_count, &rates) || chem_allocate(symbolic.nonzero_count, cell_count, &jacobian) ||
      chem_allocate(1, cell_count, &scratch)) {
    failed = 1;
    goto cleanup;
  }
  for (size_t c = 0; c < cell_count; ++c)
    for (size_t s = 0; s < n; ++s)
      y[s * cell_count + c] = y_external[c * n + s];
  if (parameters_external)
    for (size_t c = 0; c < cell_count; ++c)
      for (size_t p = 0; p < mechanism->parameter_count; ++p)
        parameters[p * cell_count + c] = parameters_external[c * mechanism->parameter_count + p];
  else
    for (size_t c = 0; c < cell_count; ++c)
      for (size_t p = 0; p < mechanism->parameter_count; ++p)
        parameters[p * cell_count + c] = mechanism->parameters[p].initial_value;

  /* Conditions and rate parameters are fixed during an ODE solve, so all
     reaction constants are invariant across its internal stages and intervals. */
  chem_calculate_rates(mechanism, cell_count, parameters, temperature, pressure, rates);

  const double gamma0 = 0.435866521508459;
  const double c10 = -1.0156171083877702091975600115545;
  const double c20 = 4.0759956452537699824805835358067;
  const double c21 = 9.2076794298330791242156818474003;
  const double a10 = 1.0;
  const double m[3] = { 1.0, 6.1697947043828245592553615689730, -0.42772256543218573326238373806514 };
  const double e[3] = { 0.5, -2.907955871680469821793236208017, 0.22354069897811569627360909276199 };
  const double h_min = 1.0e-15 * interval;
  const double h_max = interval;
  const double round_off = DBL_EPSILON;

  for (size_t interval_index = 0; interval_index < interval_count; ++interval_index) {
    double h = fmin(fmax(h_min, 1.0e-6 * interval), h_max);
    double present_time = 0.0;
    int reject_last_h = 0, reject_more_h = 0;
    size_t step_attempts = 0, accepted_this_interval = 0;
    while ((present_time - interval + round_off) <= 0.0) {
    if (step_attempts > 1000) {
      fprintf(stderr, "solver exceeded 1000 internal steps\n");
      failed = 1;
      goto cleanup;
    }
    if ((present_time + 0.1 * h) == present_time || h <= fmin(round_off, h_min)) {
      fprintf(stderr, "solver step size too small: t=%.17g h=%.17g\n", present_time, h);
      failed = 1;
      goto cleanup;
    }
    h = fmin(h, fabs(interval - present_time));
    chem_add_forcing(mechanism, cell_count, y, rates, initial_forcing, scratch);
    if (statistics) ++statistics->function_calls;
    chem_build_jacobian(mechanism, &symbolic, cell_count, y, rates, jacobian, scratch);
    if (statistics) ++statistics->jacobian_updates;

    int accepted = 0;
    while (!accepted) {
      const double alpha = 1.0 / (h * gamma0);
      for (size_t s = 0; s < n; ++s) {
        double *diagonal = jacobian + symbolic.position[s * n + s] * cell_count;
        #pragma omp simd
        for (size_t c = 0; c < cell_count; ++c)
          diagonal[c] += alpha;
      }
      size_t bad_pivot = 0, bad_cell = 0;
      if (chem_factor(&symbolic, cell_count, jacobian, scratch, &bad_pivot, &bad_cell) != 0) {
        const double bad_value = jacobian[symbolic.position[bad_pivot * n + bad_pivot] * cell_count + bad_cell];
        fprintf(stderr, "linear factorization failed at cell %zu time %.17g h %.17g pivot %zu value %.17g\n",
                bad_cell, present_time, h, bad_pivot, bad_value);
        failed = 1;
        goto cleanup;
      }
      if (statistics) ++statistics->decompositions;

      memcpy(k[0], initial_forcing, n * cell_count * sizeof(double));
      chem_solve(&symbolic, cell_count, jacobian, k[0]);
      if (statistics) ++statistics->solves;
      memcpy(ynew, y, n * cell_count * sizeof(double));
      for (size_t i = 0; i < n * cell_count; ++i)
        ynew[i] += a10 * k[0][i];
      chem_add_forcing(mechanism, cell_count, ynew, rates, k[1], scratch);
      if (statistics) ++statistics->function_calls;
      memcpy(k[2], k[1], n * cell_count * sizeof(double));
      const double c0 = c10 / h;
      for (size_t i = 0; i < n * cell_count; ++i)
        k[1][i] += c0 * k[0][i];
      chem_solve(&symbolic, cell_count, jacobian, k[1]);
      if (statistics) ++statistics->solves;

      const double c20_over_h = c20 / h;
      const double c21_over_h = c21 / h;
      for (size_t i = 0; i < n * cell_count; ++i)
        k[2][i] += c20_over_h * k[0][i] + c21_over_h * k[1][i];
      chem_solve(&symbolic, cell_count, jacobian, k[2]);
      if (statistics) ++statistics->solves;

      for (size_t i = 0; i < n * cell_count; ++i) {
        ynew[i] = y[i] + m[0] * k[0][i] + m[1] * k[1][i] + m[2] * k[2][i];
        error[i] = e[0] * k[0][i] + e[1] * k[1][i] + e[2] * k[2][i];
      }
      const double normalized_error = chem_normalized_error(mechanism, cell_count, y, ynew, error);
      double factor = 0.9 / pow(normalized_error, 1.0 / 3.0);
      factor = fmin(6.0, fmax(0.2, factor));
      double hnew = h * factor;
      if (statistics) ++statistics->internal_steps;
      ++step_attempts;
      if (!isfinite(normalized_error)) {
        fprintf(stderr, "non-finite normalized error %.17g at time %.17g h %.17g\n", normalized_error, present_time, h);
        failed = 1;
        goto cleanup;
      }
      if (normalized_error < 1.0 || h < h_min) {
        ++accepted_this_interval;
        if (statistics) ++statistics->accepted;
        present_time += h;
        memcpy(y, ynew, n * cell_count * sizeof(double));
        hnew = fmax(h_min, fmin(hnew, h_max));
        if (reject_last_h)
          hnew = fmin(hnew, h);
        reject_last_h = 0;
        reject_more_h = 0;
        h = hnew;
        accepted = 1;
      } else {
        if (reject_more_h)
          hnew = h * 0.1;
        reject_more_h = reject_last_h;
        reject_last_h = 1;
        h = hnew;
        if (statistics && accepted_this_interval >= 1)
          ++statistics->rejected;
        chem_build_jacobian(mechanism, &symbolic, cell_count, y, rates, jacobian, scratch);
        if (statistics) ++statistics->jacobian_updates;
      }
    }
  }
  }

  for (size_t c = 0; c < cell_count; ++c)
    for (size_t s = 0; s < n; ++s)
      y_external[c * n + s] = y[s * cell_count + c];

cleanup:
  free(y);
  free(ynew);
  free(initial_forcing);
  free(error);
  free(k[0]);
  free(k[1]);
  free(k[2]);
  free(parameters);
  free(rates);
  free(jacobian);
  free(scratch);
  chem_symbolic_destroy(&symbolic);
  return failed ? -1 : 0;
}

int chem_integrate_group(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    chem_statistics *statistics) {
  return chem_integrate_group_steps(
      mechanism, cell_count, y, parameters, temperature, pressure, interval, 1, statistics);
}
