#include "chem.h"
#include "jacobian_verification.h"

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
  size_t target_position;
  size_t lower_position;
  size_t upper_position;
} chem_lu_update;

typedef struct {
  size_t species_count;
  size_t reaction_count;
  size_t job_count;
  chem_jacobian_job *jobs;
  size_t *dependent_species;
  unsigned char *pattern;
  size_t *position;
  size_t nonzero_count;
  size_t *permutation;
  size_t *inverse_permutation;
  size_t *lower;
  size_t *upper;
  size_t *lower_count;
  size_t *upper_count;
  size_t *row_lower;
  size_t *row_upper;
  size_t *row_lower_count;
  size_t *row_upper_count;
  size_t *update_offset;
  size_t *update_count;
  chem_lu_update *updates;
  size_t total_updates;
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
    case CHEM_RATE_CALLBACK:
      value = rate->callback
          ? rate->callback(temperature, pressure, air_density, rate->context)
          : NAN;
      break;
    default:
      return NAN;
  }
  for (unsigned i = 0; i < reaction->rate.air_density_power; ++i)
    value *= air_density;
  return value;
}

int chem_evaluate_rate(
    const chem_mechanism *mechanism,
    size_t reaction_index,
    const double *parameters,
    double temperature,
    double pressure,
    double *rate_constant) {
  if (!mechanism || !rate_constant || reaction_index >= mechanism->reaction_count ||
      temperature <= 0.0 || pressure < 0.0)
    return -1;
  const chem_reaction *reaction = &mechanism->reactions[reaction_index];
  if (reaction->rate.kind == CHEM_RATE_USER &&
      (!parameters || reaction->rate.custom_parameter >= mechanism->parameter_count))
    return -1;
  if (reaction->rate.kind == CHEM_RATE_CALLBACK && !reaction->rate.callback)
    return -1;
  if (reaction->rate.kind == CHEM_RATE_SURFACE) {
    if (reaction->rate.surface_species >= mechanism->species_count || !parameters ||
        reaction->rate.custom_parameter >= mechanism->parameter_count ||
        reaction->rate.custom_parameter + 1u >= mechanism->parameter_count)
      return -1;
    const chem_species *species = &mechanism->species[reaction->rate.surface_species];
    if (species->molecular_weight <= 0.0 || species->diffusion_coefficient <= 0.0)
      return -2;
  }
  const double air_density = pressure / (CHEM_GAS_CONSTANT * temperature);
  *rate_constant = chem_rate_value(mechanism, reaction, parameters, 0, 1,
                                   temperature, pressure, air_density);
  return isfinite(*rate_constant) ? 0 : -1;
}

chem_rosenbrock_parameters chem_rosenbrock_default_parameters(void) {
  chem_rosenbrock_parameters parameters = { 0 };
  parameters.stage_count = 3;
  parameters.new_function_evaluation[0] = 1u;
  parameters.new_function_evaluation[1] = 1u;
  parameters.a[0] = 1.0;
  parameters.c[0] = -1.0156171083877702091975600115545;
  parameters.c[1] = 4.0759956452537699824805835358067;
  parameters.c[2] = 9.2076794298330791242156818474003;
  parameters.m[0] = 1.0;
  parameters.m[1] = 6.1697947043828245592553615689730;
  parameters.m[2] = -0.42772256543218573326238373806514;
  parameters.e[0] = 0.5;
  parameters.e[1] = -2.907955871680469821793236208017;
  parameters.e[2] = 0.22354069897811569627360909276199;
  parameters.gamma0 = 0.43586652150845899941601945119356;
  parameters.estimator_order = 3.0;
  parameters.round_off = DBL_EPSILON;
  parameters.factor_min = 0.2;
  parameters.factor_max = 6.0;
  parameters.rejection_factor_decrease = 0.1;
  parameters.safety_factor = 0.9;
  parameters.h_start_fraction = 1.0e-6;
  parameters.h_min_fraction = 1.0e-15;
  parameters.h_max_fraction = 1.0;
  parameters.max_internal_steps = 1000u;
  return parameters;
}

static chem_rosenbrock_parameters chem_rosenbrock_empty_method(size_t stages) {
  chem_rosenbrock_parameters parameters = chem_rosenbrock_default_parameters();
  parameters.stage_count = stages;
  memset(parameters.new_function_evaluation, 0, sizeof(parameters.new_function_evaluation));
  memset(parameters.a, 0, sizeof(parameters.a));
  memset(parameters.c, 0, sizeof(parameters.c));
  memset(parameters.m, 0, sizeof(parameters.m));
  memset(parameters.e, 0, sizeof(parameters.e));
  return parameters;
}

chem_rosenbrock_parameters chem_rosenbrock_two_stage_parameters(void) {
  chem_rosenbrock_parameters parameters = chem_rosenbrock_empty_method(2u);
  const double gamma = 1.0 + 1.0 / sqrt(2.0);
  parameters.new_function_evaluation[0] = 1u;
  parameters.new_function_evaluation[1] = 1u;
  parameters.a[0] = 1.0 / gamma;
  parameters.c[0] = -2.0 / gamma;
  parameters.m[0] = 3.0 / (2.0 * gamma);
  parameters.m[1] = 1.0 / (2.0 * gamma);
  parameters.e[0] = 1.0 / (2.0 * gamma);
  parameters.e[1] = 1.0 / (2.0 * gamma);
  parameters.gamma0 = gamma;
  parameters.estimator_order = 2.0;
  return parameters;
}

chem_rosenbrock_parameters chem_rosenbrock_four_stage_parameters(void) {
  chem_rosenbrock_parameters parameters = chem_rosenbrock_empty_method(4u);
  const double a[] = {
    2.0, 1.867943637803922, 0.2344449711399156, 1.867943637803922,
    0.2344449711399156, 0.0,
  };
  const double c[] = {
    -7.137615036412310, 2.580708087951457, 0.6515950076447975,
    -2.137148994382534, -0.3214669691237626, -0.6949742501781779,
  };
  memcpy(parameters.a, a, sizeof(a));
  memcpy(parameters.c, c, sizeof(c));
  parameters.new_function_evaluation[0] = 1u;
  parameters.new_function_evaluation[1] = 1u;
  parameters.new_function_evaluation[2] = 1u;
  parameters.m[0] = 2.255570073418735;
  parameters.m[1] = 0.2870493262186792;
  parameters.m[2] = 0.4353179431840180;
  parameters.m[3] = 1.093502252409163;
  parameters.e[0] = -0.2815431932141155;
  parameters.e[1] = -0.07276199124938920;
  parameters.e[2] = -0.1082196201495311;
  parameters.e[3] = -1.093502252409163;
  parameters.gamma0 = 0.5728200000000000;
  parameters.estimator_order = 4.0;
  return parameters;
}

chem_rosenbrock_parameters chem_rosenbrock_four_stage_da_parameters(void) {
  chem_rosenbrock_parameters parameters = chem_rosenbrock_empty_method(4u);
  const double a[] = { 0.0, 2.0, 0.0, 2.0, 0.0, 1.0 };
  const double c[] = { 4.0, 1.0, -1.0, 1.0, -1.0, -8.0 / 3.0 };
  memcpy(parameters.a, a, sizeof(a));
  memcpy(parameters.c, c, sizeof(c));
  parameters.new_function_evaluation[0] = 1u;
  parameters.new_function_evaluation[2] = 1u;
  parameters.new_function_evaluation[3] = 1u;
  parameters.m[0] = 2.0;
  parameters.m[2] = 1.0;
  parameters.m[3] = 1.0;
  parameters.e[3] = 1.0;
  parameters.gamma0 = 0.5;
  parameters.estimator_order = 3.0;
  return parameters;
}

chem_rosenbrock_parameters chem_rosenbrock_six_stage_da_parameters(void) {
  chem_rosenbrock_parameters parameters = chem_rosenbrock_empty_method(6u);
  const double a[] = {
    1.544, 0.9466785280815826, 0.2557011698983284, 3.314825187068521,
    2.896124015972201, 0.9986419139977817, 1.221224509226641,
    6.019134481288629, 12.53708332932087, -0.6878860361058950,
    1.221224509226641, 6.019134481288629, 12.53708332932087,
    -0.6878860361058950, 1.0,
  };
  const double c[] = {
    -5.6688, -2.430093356833875, -0.2063599157091915,
    -0.1073529058151375, -9.594562251023355, -20.47028614809616,
    7.496443313967647, -10.24680431464352, -33.99990352819905,
    11.70890893206160, 8.083246795921522, -7.981132988064893,
    -31.52159432874371, 16.31930543123136, -6.058818238834054,
  };
  memcpy(parameters.a, a, sizeof(a));
  memcpy(parameters.c, c, sizeof(c));
  for (size_t stage = 0; stage < 6; ++stage)
    parameters.new_function_evaluation[stage] = 1u;
  parameters.m[0] = a[6];
  parameters.m[1] = a[7];
  parameters.m[2] = a[8];
  parameters.m[3] = a[9];
  parameters.m[4] = 1.0;
  parameters.m[5] = 1.0;
  parameters.e[5] = 1.0;
  parameters.gamma0 = 0.25;
  parameters.estimator_order = 4.0;
  return parameters;
}

static int chem_rosenbrock_parameters_valid(const chem_rosenbrock_parameters *parameters) {
  if (!parameters || !parameters->stage_count || parameters->stage_count > CHEM_ROSENBROCK_MAX_STAGES ||
      !parameters->new_function_evaluation[0] || !isfinite(parameters->gamma0) || parameters->gamma0 <= 0.0 ||
      (parameters->lu_method != CHEM_LU_MOZART_IN_PLACE && parameters->lu_method != CHEM_LU_DOOLITTLE) ||
      !isfinite(parameters->estimator_order) || parameters->estimator_order <= 0.0 ||
      !isfinite(parameters->round_off) || parameters->round_off < 0.0 ||
      !isfinite(parameters->factor_min) || parameters->factor_min <= 0.0 ||
      !isfinite(parameters->factor_max) || parameters->factor_max < parameters->factor_min ||
      !isfinite(parameters->rejection_factor_decrease) || parameters->rejection_factor_decrease <= 0.0 ||
      !isfinite(parameters->safety_factor) || parameters->safety_factor <= 0.0 ||
      !isfinite(parameters->h_start_fraction) || parameters->h_start_fraction <= 0.0 ||
      !isfinite(parameters->h_min_fraction) || parameters->h_min_fraction <= 0.0 ||
      !isfinite(parameters->h_max_fraction) || parameters->h_max_fraction < parameters->h_min_fraction ||
      parameters->h_max_fraction <= 0.0 ||
      !parameters->max_internal_steps)
    return 0;
  const size_t coefficient_count = parameters->stage_count * (parameters->stage_count - 1u) / 2u;
  for (size_t i = 0; i < parameters->stage_count; ++i) {
    if (parameters->new_function_evaluation[i] > 1u || !isfinite(parameters->m[i]) || !isfinite(parameters->e[i]))
      return 0;
  }
  for (size_t i = 0; i < coefficient_count; ++i)
    if (!isfinite(parameters->a[i]) || !isfinite(parameters->c[i])) return 0;
  return 1;
}

static void chem_symbolic_destroy(chem_symbolic *symbolic) {
  free(symbolic->jobs);
  free(symbolic->dependent_species);
  free(symbolic->pattern);
  free(symbolic->position);
  free(symbolic->permutation);
  free(symbolic->inverse_permutation);
  free(symbolic->lower);
  free(symbolic->upper);
  free(symbolic->lower_count);
  free(symbolic->upper_count);
  free(symbolic->row_lower);
  free(symbolic->row_upper);
  free(symbolic->row_lower_count);
  free(symbolic->row_upper_count);
  free(symbolic->update_offset);
  free(symbolic->update_count);
  free(symbolic->updates);
  memset(symbolic, 0, sizeof(*symbolic));
}

static int chem_markowitz_order(const unsigned char *pattern, size_t n, size_t *permutation) {
  unsigned char *graph = malloc(n * n * sizeof(unsigned char));
  size_t *row_degree = calloc(n, sizeof(size_t));
  size_t *column_degree = calloc(n, sizeof(size_t));
  unsigned char *alive = malloc(n * sizeof(unsigned char));
  size_t *columns = malloc(n * sizeof(size_t));
  size_t *incoming = malloc(n * sizeof(size_t));
  if (!graph || !row_degree || !column_degree || !alive || !columns || !incoming) {
    free(graph); free(row_degree); free(column_degree); free(alive); free(columns); free(incoming);
    return -1;
  }
  memcpy(graph, pattern, n * n * sizeof(unsigned char));
  memset(alive, 1, n * sizeof(unsigned char));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (graph[i * n + j]) {
        ++row_degree[i];
        ++column_degree[j];
      }
    }
  }
  for (size_t step = 0; step < n; ++step) {
    size_t pivot = n;
    size_t best_cost = SIZE_MAX;
    for (size_t v = 0; v < n; ++v) {
      if (!alive[v])
        continue;
      const size_t cost = (row_degree[v] - 1u) * (column_degree[v] - 1u);
      if (pivot == n || cost < best_cost) {
        best_cost = cost;
        pivot = v;
      }
    }
    permutation[step] = pivot;
    alive[pivot] = 0;
    size_t column_count = 0, incoming_count = 0;
    for (size_t c = 0; c < n; ++c)
      if (c != pivot && alive[c] && graph[pivot * n + c])
        columns[column_count++] = c;
    for (size_t i = 0; i < n; ++i)
      if (i != pivot && alive[i] && graph[i * n + pivot])
        incoming[incoming_count++] = i;
    for (size_t i = 0; i < incoming_count; ++i) {
      const size_t row = incoming[i];
      for (size_t j = 0; j < column_count; ++j) {
        const size_t col = columns[j];
        if (!graph[row * n + col]) {
          graph[row * n + col] = 1;
          ++row_degree[row];
          ++column_degree[col];
        }
      }
    }
    for (size_t j = 0; j < column_count; ++j) {
      const size_t col = columns[j];
      graph[pivot * n + col] = 0;
      --column_degree[col];
    }
    for (size_t i = 0; i < incoming_count; ++i) {
      const size_t row = incoming[i];
      graph[row * n + pivot] = 0;
      --row_degree[row];
    }
  }
  free(graph); free(row_degree); free(column_degree); free(alive); free(columns); free(incoming);
  return 0;
}

static int chem_symbolic_build(
    const chem_mechanism *mechanism, chem_lu_method lu_method, chem_symbolic *symbolic) {
  const size_t n = mechanism->species_count;
  size_t jobs = 0, dependent_count = 0;
  memset(symbolic, 0, sizeof(*symbolic));
  if (n == 0 || n > SIZE_MAX / n || n * n > SIZE_MAX / sizeof(size_t))
    return -1;
  symbolic->species_count = n;
  symbolic->reaction_count = mechanism->reaction_count;
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const size_t reactants = mechanism->reactions[r].reactant_count;
    if (jobs > SIZE_MAX - reactants || (reactants && reactants > SIZE_MAX / reactants))
      return -1;
    const size_t dependencies = reactants * reactants;
    if (dependent_count > SIZE_MAX - dependencies)
      return -1;
    jobs += reactants;
    dependent_count += dependencies;
  }
  if (jobs > SIZE_MAX / sizeof(*symbolic->jobs) || dependent_count > SIZE_MAX / sizeof(size_t))
    return -1;
  symbolic->jobs = calloc(jobs ? jobs : 1, sizeof(*symbolic->jobs));
  symbolic->dependent_species = malloc((dependent_count ? dependent_count : 1) * sizeof(size_t));
  symbolic->pattern = calloc(n * n, sizeof(unsigned char));
  symbolic->permutation = malloc(n * sizeof(size_t));
  symbolic->inverse_permutation = malloc(n * sizeof(size_t));
  symbolic->lower = malloc(n * n * sizeof(size_t));
  symbolic->upper = malloc(n * n * sizeof(size_t));
  symbolic->lower_count = calloc(n, sizeof(size_t));
  symbolic->upper_count = calloc(n, sizeof(size_t));
  symbolic->row_lower = malloc(n * n * sizeof(size_t));
  symbolic->row_upper = malloc(n * n * sizeof(size_t));
  symbolic->row_lower_count = calloc(n, sizeof(size_t));
  symbolic->row_upper_count = calloc(n, sizeof(size_t));
  symbolic->update_offset = calloc(n, sizeof(size_t));
  symbolic->update_count = calloc(n, sizeof(size_t));
  if (!symbolic->jobs || !symbolic->dependent_species || !symbolic->pattern || !symbolic->permutation ||
      !symbolic->inverse_permutation || !symbolic->lower ||
      !symbolic->upper || !symbolic->lower_count || !symbolic->upper_count || !symbolic->row_lower ||
      !symbolic->row_upper || !symbolic->row_lower_count || !symbolic->row_upper_count ||
      !symbolic->update_offset || !symbolic->update_count) {
    chem_symbolic_destroy(symbolic);
    return -1;
  }
  symbolic->job_count = jobs;
  for (size_t i = 0; i < n; ++i)
    symbolic->pattern[i * n + i] = 1;
  for (size_t r = 0; r < mechanism->reaction_count; ++r) {
    const chem_reaction *reaction = &mechanism->reactions[r];
    for (size_t q = 0; q < reaction->reactant_count; ++q) {
      const size_t column = mechanism->reactant_terms[reaction->reactant_offset + q].species;
      for (size_t d = 0; d < reaction->reactant_count; ++d)
        symbolic->pattern[mechanism->reactant_terms[reaction->reactant_offset + d].species * n + column] = 1;
      for (size_t d = 0; d < reaction->product_count; ++d)
        symbolic->pattern[mechanism->product_terms[reaction->product_offset + d].species * n + column] = 1;
    }
  }

  if (lu_method == CHEM_LU_DOOLITTLE) {
    for (size_t i = 0; i < n; ++i)
      symbolic->permutation[i] = i;
  } else if (chem_markowitz_order(symbolic->pattern, n, symbolic->permutation) != 0) {
    chem_symbolic_destroy(symbolic);
    return -1;
  }
  for (size_t i = 0; i < n; ++i)
    symbolic->inverse_permutation[symbolic->permutation[i]] = i;
  unsigned char *permuted_pattern = calloc(n * n, sizeof(unsigned char));
  if (!permuted_pattern) {
    chem_symbolic_destroy(symbolic);
    return -1;
  }
  for (size_t row = 0; row < n; ++row)
    for (size_t col = 0; col < n; ++col)
      if (symbolic->pattern[symbolic->permutation[row] * n + symbolic->permutation[col]])
        permuted_pattern[row * n + col] = 1;
  memcpy(symbolic->pattern, permuted_pattern, n * n * sizeof(unsigned char));
  free(permuted_pattern);

  size_t job_pos = 0, dep_pos = 0;
  for (size_t independent = 0; independent < n; ++independent) {
    const size_t original_independent = symbolic->permutation[independent];
    for (size_t r = 0; r < mechanism->reaction_count; ++r) {
      const chem_reaction *reaction = &mechanism->reactions[r];
      for (size_t q = 0; q < reaction->reactant_count; ++q) {
        const size_t column = mechanism->reactant_terms[reaction->reactant_offset + q].species;
        if (column != original_independent)
          continue;
        size_t found = 0;
        chem_jacobian_job *job = &symbolic->jobs[job_pos++];
        job->reaction = r;
        job->independent = independent;
        job->dependent_offset = dep_pos;
        for (size_t d = 0; d < reaction->reactant_count; ++d) {
          const size_t original_row = mechanism->reactant_terms[reaction->reactant_offset + d].species;
          if (original_row == original_independent && !found) {
            found = 1;
            continue;
          }
          symbolic->dependent_species[dep_pos++] = symbolic->inverse_permutation[original_row];
          ++job->dependent_count;
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
  for (size_t row = 0; row < n; ++row) {
    for (size_t col = 0; col < row; ++col)
      if (symbolic->pattern[row * n + col])
        symbolic->row_lower[row * n + symbolic->row_lower_count[row]++] = col;
    for (size_t col = row + 1; col < n; ++col)
      if (symbolic->pattern[row * n + col])
        symbolic->row_upper[row * n + symbolic->row_upper_count[row]++] = col;
  }
  for (size_t k = 0; k < n; ++k) {
    symbolic->update_offset[k] = symbolic->total_updates;
    for (size_t uj = 0; uj < symbolic->upper_count[k]; ++uj) {
      const size_t col = symbolic->upper[k * n + uj];
      for (size_t li = 0; li < symbolic->lower_count[k]; ++li) {
        const size_t row = symbolic->lower[k * n + li];
        if (symbolic->pattern[row * n + col])
          ++symbolic->total_updates;
      }
    }
    symbolic->update_count[k] = symbolic->total_updates - symbolic->update_offset[k];
  }
  symbolic->updates = malloc((symbolic->total_updates ? symbolic->total_updates : 1) * sizeof(*symbolic->updates));
  if (!symbolic->updates) {
    chem_symbolic_destroy(symbolic);
    return -1;
  }
  size_t update_pos = 0;
  for (size_t k = 0; k < n; ++k) {
    for (size_t uj = 0; uj < symbolic->upper_count[k]; ++uj) {
      const size_t col = symbolic->upper[k * n + uj];
      for (size_t li = 0; li < symbolic->lower_count[k]; ++li) {
        const size_t row = symbolic->lower[k * n + li];
        if (symbolic->pattern[row * n + col]) {
          symbolic->updates[update_pos++] = (chem_lu_update){
            symbolic->position[row * n + col],
            symbolic->position[row * n + k],
            symbolic->position[k * n + col],
          };
        }
      }
    }
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
    const chem_symbolic *symbolic,
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
        const size_t species = symbolic->inverse_permutation[
            mechanism->reactant_terms[reaction->reactant_offset + q].species];
        rate *= y[species * cell_count + c];
      }
      flux[c] = rate;
    }
    for (size_t q = 0; q < reaction->reactant_count; ++q) {
      const size_t species = symbolic->inverse_permutation[
          mechanism->reactant_terms[reaction->reactant_offset + q].species];
      double *forcing_row = forcing + species * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        forcing_row[c] -= flux[c];
    }
    for (size_t q = 0; q < reaction->product_count; ++q) {
      const chem_term *product = &mechanism->product_terms[reaction->product_offset + q];
      double *forcing_row = forcing + symbolic->inverse_permutation[product->species] * cell_count;
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
      jacobian[symbolic->position[symbolic->inverse_permutation[product->species] * n + job->independent] * cell_count + c] -=
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
    for (size_t i = 0; i < symbolic->lower_count[k]; ++i) {
      double *lower_entry = matrix + symbolic->position[lower[i] * n + k] * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        lower_entry[c] *= inverse[c];
    }
    for (size_t i = 0; i < symbolic->update_count[k]; ++i) {
      const chem_lu_update *update = &symbolic->updates[symbolic->update_offset[k] + i];
      double *target = matrix + update->target_position * cell_count;
      const double *lower_entry = matrix + update->lower_position * cell_count;
      const double *upper_entry = matrix + update->upper_position * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        target[c] -= lower_entry[c] * upper_entry[c];
    }
  }
  return 0;
}

static int chem_factor_doolittle(
    const chem_symbolic *symbolic,
    size_t cell_count,
    double *matrix,
    size_t *bad_pivot,
    size_t *bad_cell) {
  const size_t n = symbolic->species_count;
  for (size_t row = 0; row < n; ++row) {
    const size_t diagonal_position = symbolic->position[row * n + row];
    const size_t upper_count = symbolic->row_upper_count[row];
    for (size_t upper_index = 0; upper_index <= upper_count; ++upper_index) {
      const size_t col = upper_index == 0 ? row : symbolic->row_upper[row * n + upper_index - 1u];
      const size_t target_position = symbolic->position[row * n + col];
      double *target = matrix + target_position * cell_count;
      for (size_t lower_index = 0; lower_index < symbolic->row_lower_count[row]; ++lower_index) {
        const size_t inner = symbolic->row_lower[row * n + lower_index];
        const size_t upper_position = symbolic->position[inner * n + col];
        if (upper_position == SIZE_MAX)
          continue;
        const double *lower_value = matrix + symbolic->position[row * n + inner] * cell_count;
        const double *upper_value = matrix + upper_position * cell_count;
        #pragma omp simd
        for (size_t cell = 0; cell < cell_count; ++cell)
          target[cell] -= lower_value[cell] * upper_value[cell];
      }
    }

    const double *diagonal = matrix + diagonal_position * cell_count;
    for (size_t cell = 0; cell < cell_count; ++cell)
      if (diagonal[cell] == 0.0 || !isfinite(diagonal[cell])) {
        *bad_pivot = row;
        *bad_cell = cell;
        return -1;
      }

    for (size_t lower_index = 0; lower_index < symbolic->lower_count[row]; ++lower_index) {
      const size_t lower_row = symbolic->lower[row * n + lower_index];
      const size_t target_position = symbolic->position[lower_row * n + row];
      double *target = matrix + target_position * cell_count;
      for (size_t inner_index = 0; inner_index < symbolic->row_lower_count[lower_row]; ++inner_index) {
        const size_t inner = symbolic->row_lower[lower_row * n + inner_index];
        if (inner >= row)
          break;
        const size_t upper_position = symbolic->position[inner * n + row];
        if (upper_position == SIZE_MAX)
          continue;
        const double *lower_value = matrix + symbolic->position[lower_row * n + inner] * cell_count;
        const double *upper_value = matrix + upper_position * cell_count;
        #pragma omp simd
        for (size_t cell = 0; cell < cell_count; ++cell)
          target[cell] -= lower_value[cell] * upper_value[cell];
      }
      #pragma omp simd
      for (size_t cell = 0; cell < cell_count; ++cell)
        target[cell] /= diagonal[cell];
    }
  }
  return 0;
}

static void chem_solve(const chem_symbolic *symbolic, size_t cell_count, const double *matrix, double *rhs) {
  const size_t n = symbolic->species_count;
  for (size_t i = 0; i < n; ++i) {
    double *rhs_i = rhs + i * cell_count;
    for (size_t q = 0; q < symbolic->row_lower_count[i]; ++q) {
      const size_t j = symbolic->row_lower[i * n + q];
      const double *entry = matrix + symbolic->position[i * n + j] * cell_count;
      const double *rhs_j = rhs + j * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        rhs_i[c] -= entry[c] * rhs_j[c];
    }
  }
  for (size_t ii = n; ii-- > 0;) {
    double *rhs_i = rhs + ii * cell_count;
    for (size_t q = 0; q < symbolic->row_upper_count[ii]; ++q) {
      const size_t j = symbolic->row_upper[ii * n + q];
      const double *entry = matrix + symbolic->position[ii * n + j] * cell_count;
      const double *rhs_j = rhs + j * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        rhs_i[c] -= entry[c] * rhs_j[c];
    }
    const double *diagonal = matrix + symbolic->position[ii * n + ii] * cell_count;
    #pragma omp simd
    for (size_t c = 0; c < cell_count; ++c)
      rhs_i[c] /= diagonal[c];
  }
}

static double chem_normalized_error(
    const chem_mechanism *mechanism,
    const chem_symbolic *symbolic,
    size_t cell_count,
    const double *restrict y,
    const double *restrict ynew,
    const double *restrict error) {
  const size_t n = mechanism->species_count;
  double sum = 0.0;
  for (size_t species = 0; species < n; ++species) {
    const chem_species *item = &mechanism->species[symbolic->permutation[species]];
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

static double chem_step_size_factor(const chem_rosenbrock_parameters *parameters, double normalized_error) {
  const double factor = parameters->safety_factor / pow(normalized_error, 1.0 / parameters->estimator_order);
  return fmin(parameters->factor_max, fmax(parameters->factor_min, factor));
}

static int chem_allocate(size_t count, size_t width, double **out) {
  if (width && count > SIZE_MAX / width / sizeof(double))
    return -1;
  const size_t total = count * width;
  *out = calloc(total ? total : 1, sizeof(double));
  return *out ? 0 : -1;
}

int chem_integrate_group_steps_with_parameters(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y_external,
    double *restrict parameters_external,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    size_t interval_count,
    const chem_rosenbrock_parameters *solver_parameters,
    chem_statistics *statistics) {
  if (!mechanism || !cell_count || !interval_count || !y_external || !temperature || !pressure || interval <= 0.0 ||
      !mechanism->species_count || !mechanism->reaction_count ||
      !chem_rosenbrock_parameters_valid(solver_parameters))
    return -1;
  const chem_rosenbrock_parameters *solver = solver_parameters;
  const size_t n = mechanism->species_count;
  const size_t nr = mechanism->reaction_count;
  chem_symbolic symbolic;
  if (chem_symbolic_build(mechanism, solver->lu_method, &symbolic) != 0)
    return -1;
  double *y = NULL, *ynew = NULL, *initial_forcing = NULL, *stage_forcing = NULL, *error = NULL;
  double *k[CHEM_ROSENBROCK_MAX_STAGES] = { NULL };
  double *parameters = NULL, *rates = NULL, *jacobian = NULL, *scratch = NULL;
  int failed = 0;
  if (chem_allocate(n, cell_count, &y) || chem_allocate(n, cell_count, &ynew) ||
      chem_allocate(n, cell_count, &initial_forcing) || chem_allocate(n, cell_count, &stage_forcing) ||
      chem_allocate(n, cell_count, &error) || chem_allocate(mechanism->parameter_count, cell_count, &parameters) ||
      chem_allocate(nr, cell_count, &rates) || chem_allocate(symbolic.nonzero_count, cell_count, &jacobian) ||
      chem_allocate(1, cell_count, &scratch)) {
    failed = 1;
    goto cleanup;
  }
  for (size_t stage = 0; stage < solver->stage_count; ++stage)
    if (chem_allocate(n, cell_count, &k[stage]) != 0) {
      failed = 1;
      goto cleanup;
    }
  for (size_t c = 0; c < cell_count; ++c)
    for (size_t s = 0; s < n; ++s)
      y[s * cell_count + c] = y_external[c * n + symbolic.permutation[s]];
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

  const double h_min = solver->h_min_fraction * interval;
  const double h_max = fmin(solver->h_max_fraction * interval, interval);
  const double h_start = fmin(solver->h_start_fraction * interval, h_max);
  const double round_off = solver->round_off;

  for (size_t interval_index = 0; interval_index < interval_count; ++interval_index) {
    double h = fmin(fmax(h_min, h_start), h_max);
    double present_time = 0.0;
    int reject_last_h = 0, reject_more_h = 0;
    size_t step_attempts = 0, accepted_this_interval = 0;
    while ((present_time - interval + round_off) <= 0.0) {
      if (step_attempts > solver->max_internal_steps) {
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
      chem_add_forcing(mechanism, &symbolic, cell_count, y, rates, initial_forcing, scratch);
      if (statistics) ++statistics->function_calls;
      chem_build_jacobian(mechanism, &symbolic, cell_count, y, rates, jacobian, scratch);
      if (statistics) ++statistics->jacobian_updates;

      int accepted = 0;
      while (!accepted) {
        const double alpha = 1.0 / (h * solver->gamma0);
        for (size_t s = 0; s < n; ++s) {
          double *diagonal = jacobian + symbolic.position[s * n + s] * cell_count;
          #pragma omp simd
          for (size_t c = 0; c < cell_count; ++c)
            diagonal[c] += alpha;
        }
        size_t bad_pivot = 0, bad_cell = 0;
        const int factor_status = solver->lu_method == CHEM_LU_DOOLITTLE
            ? chem_factor_doolittle(&symbolic, cell_count, jacobian, &bad_pivot, &bad_cell)
            : chem_factor(&symbolic, cell_count, jacobian, scratch, &bad_pivot, &bad_cell);
        if (factor_status != 0) {
          const double bad_value = jacobian[symbolic.position[bad_pivot * n + bad_pivot] * cell_count + bad_cell];
          fprintf(stderr, "linear factorization failed at cell %zu time %.17g h %.17g pivot %zu value %.17g\n",
                  bad_cell, present_time, h, bad_pivot, bad_value);
          failed = 1;
          goto cleanup;
        }
        if (statistics) ++statistics->decompositions;

        const size_t state_size = n * cell_count;
        memcpy(stage_forcing, initial_forcing, state_size * sizeof(double));
        for (size_t stage = 0; stage < solver->stage_count; ++stage) {
          if (stage == 0) {
            memcpy(k[stage], initial_forcing, state_size * sizeof(double));
          } else {
            if (solver->new_function_evaluation[stage]) {
              const size_t coefficient_offset = stage * (stage - 1u) / 2u;
              memcpy(ynew, y, state_size * sizeof(double));
              for (size_t previous = 0; previous < stage; ++previous) {
                const double coefficient = solver->a[coefficient_offset + previous];
                for (size_t i = 0; i < state_size; ++i)
                  ynew[i] += coefficient * k[previous][i];
              }
              chem_add_forcing(mechanism, &symbolic, cell_count, ynew, rates, stage_forcing, scratch);
              if (statistics) ++statistics->function_calls;
            }
            memcpy(k[stage], stage_forcing, state_size * sizeof(double));
            const size_t coefficient_offset = stage * (stage - 1u) / 2u;
            for (size_t previous = 0; previous < stage; ++previous) {
              const double coefficient = solver->c[coefficient_offset + previous] / h;
              for (size_t i = 0; i < state_size; ++i)
                k[stage][i] += coefficient * k[previous][i];
            }
          }
          chem_solve(&symbolic, cell_count, jacobian, k[stage]);
          if (statistics) ++statistics->solves;
        }

        memcpy(ynew, y, state_size * sizeof(double));
        memset(error, 0, state_size * sizeof(double));
        for (size_t stage = 0; stage < solver->stage_count; ++stage) {
          for (size_t i = 0; i < state_size; ++i) {
            ynew[i] += solver->m[stage] * k[stage][i];
            error[i] += solver->e[stage] * k[stage][i];
          }
        }
        const double normalized_error = chem_normalized_error(mechanism, &symbolic, cell_count, y, ynew, error);
        const double factor = chem_step_size_factor(solver, normalized_error);
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
            hnew = h * solver->rejection_factor_decrease;
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
    /* Remove accepted roundoff undershoots within the configured tolerance at
       interval boundaries, so they cannot become artificial sources in the
       next interval through very fast reactions. */
    for (size_t species = 0; species < n; ++species) {
      const double atol = mechanism->species[symbolic.permutation[species]].absolute_tolerance;
      double *row = y + species * cell_count;
      #pragma omp simd
      for (size_t c = 0; c < cell_count; ++c)
        if (row[c] < 0.0 && row[c] >= -atol)
          row[c] = 0.0;
    }
  }

  for (size_t c = 0; c < cell_count; ++c)
    for (size_t s = 0; s < n; ++s)
      y_external[c * n + symbolic.permutation[s]] = y[s * cell_count + c];

cleanup:
  free(y);
  free(ynew);
  free(initial_forcing);
  free(stage_forcing);
  free(error);
  for (size_t stage = 0; stage < CHEM_ROSENBROCK_MAX_STAGES; ++stage)
    free(k[stage]);
  free(parameters);
  free(rates);
  free(jacobian);
  free(scratch);
  chem_symbolic_destroy(&symbolic);
  return failed ? -1 : 0;
}

int chem_integrate_group_steps(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    size_t interval_count,
    chem_statistics *statistics) {
  const chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  return chem_integrate_group_steps_with_parameters(
      mechanism, cell_count, y, parameters, temperature, pressure, interval, interval_count, &solver, statistics);
}

int chem_integrate_group_with_parameters(
    const chem_mechanism *mechanism,
    size_t cell_count,
    double *restrict y,
    double *restrict parameters,
    const double *restrict temperature,
    const double *restrict pressure,
    double interval,
    const chem_rosenbrock_parameters *solver_parameters,
    chem_statistics *statistics) {
  return chem_integrate_group_steps_with_parameters(
      mechanism, cell_count, y, parameters, temperature, pressure, interval, 1, solver_parameters, statistics);
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
  const chem_rosenbrock_parameters solver = chem_rosenbrock_default_parameters();
  return chem_integrate_group_with_parameters(
      mechanism, cell_count, y, parameters, temperature, pressure, interval, &solver, statistics);
}

int chem_finite_difference_jacobian(
    chem_forcing_callback forcing,
    void *context,
    const double *state,
    size_t block_count,
    size_t species_count,
    double perturbation,
    double *jacobian) {
  if (!forcing || !state || !jacobian || !block_count || !species_count ||
      !isfinite(perturbation) || perturbation <= 0.0 || species_count > SIZE_MAX / species_count ||
      block_count > SIZE_MAX / (species_count * species_count) || species_count > SIZE_MAX / (5u * sizeof(double)))
    return -1;
  const size_t width = species_count * species_count;
  double *workspace = malloc(5u * species_count * sizeof(double));
  if (!workspace)
    return -1;
  double *plus_state = workspace;
  double *minus_state = plus_state + species_count;
  double *plus_forcing = minus_state + species_count;
  double *minus_forcing = plus_forcing + species_count;
  double *base_forcing = minus_forcing + species_count;
  for (size_t block = 0; block < block_count; ++block) {
    const double *base_state = state + block * species_count;
    double *block_jacobian = jacobian + block * width;
    for (size_t column = 0; column < species_count; ++column) {
      memcpy(plus_state, base_state, species_count * sizeof(double));
      memcpy(minus_state, base_state, species_count * sizeof(double));
      const double x = base_state[column];
      const double h = perturbation * fmax(1.0, fabs(x));
      if (!isfinite(x) || !isfinite(h) || h <= 0.0) {
        free(workspace);
        return -1;
      }
      const int use_central = x - h >= 0.0;
      plus_state[column] = x + h;
      if (use_central)
        minus_state[column] = x - h;
      memset(plus_forcing, 0, species_count * sizeof(double));
      if (forcing(block, plus_state, species_count, plus_forcing, context) != 0) {
        free(workspace);
        return -1;
      }
      for (size_t row = 0; row < species_count; ++row) {
        if (!isfinite(plus_forcing[row])) {
          free(workspace);
          return -1;
        }
      }
      if (use_central) {
        memset(minus_forcing, 0, species_count * sizeof(double));
        if (forcing(block, minus_state, species_count, minus_forcing, context) != 0) {
          free(workspace);
          return -1;
        }
        for (size_t row = 0; row < species_count; ++row)
          if (!isfinite(minus_forcing[row])) {
            free(workspace);
            return -1;
          }
      } else {
        memset(base_forcing, 0, species_count * sizeof(double));
        if (forcing(block, base_state, species_count, base_forcing, context) != 0) {
          free(workspace);
          return -1;
        }
      }
      for (size_t row = 0; row < species_count; ++row) {
        if (!use_central && !isfinite(base_forcing[row])) {
          free(workspace);
          return -1;
        }
        block_jacobian[row * species_count + column] = use_central
            ? (plus_forcing[row] - minus_forcing[row]) / (2.0 * h)
            : (plus_forcing[row] - base_forcing[row]) / h;
        if (!isfinite(block_jacobian[row * species_count + column])) {
          free(workspace);
          return -1;
        }
      }
    }
  }
  free(workspace);
  return 0;
}

int chem_compare_jacobian_to_finite_difference(
    const unsigned char *pattern,
    const double *negative_jacobian,
    const double *finite_difference_jacobian,
    size_t block_count,
    size_t species_count,
    double absolute_tolerance,
    double relative_tolerance,
    chem_jacobian_comparison *result) {
  if (!pattern || !negative_jacobian || !finite_difference_jacobian || !result || !block_count ||
      !species_count || species_count > SIZE_MAX / species_count ||
      block_count > SIZE_MAX / (species_count * species_count) || !isfinite(absolute_tolerance) ||
      !isfinite(relative_tolerance) || absolute_tolerance < 0.0 || relative_tolerance < 0.0)
    return -1;
  memset(result, 0, sizeof(*result));
  result->passed = 1;
  const size_t width = species_count * species_count;
  for (size_t block = 0; block < block_count; ++block) {
    for (size_t row = 0; row < species_count; ++row) {
      for (size_t column = 0; column < species_count; ++column) {
        const size_t index = row * species_count + column;
        if (!pattern[index]) continue;
        const size_t offset = block * width + index;
        const double analytical = -negative_jacobian[offset];
        const double fd = finite_difference_jacobian[offset];
        if (!isfinite(analytical) || !isfinite(fd)) {
          result->passed = 0;
          result->max_absolute_error = INFINITY;
          result->max_relative_error = INFINITY;
          result->worst_block = block;
          result->worst_row = row;
          result->worst_column = column;
          result->worst_analytical = analytical;
          result->worst_finite_difference = fd;
          continue;
        }
        const double absolute_error = fabs(analytical - fd);
        const double scale = fmax(fabs(analytical), fabs(fd));
        const double relative_error = scale > 0.0 ? absolute_error / scale : 0.0;
        if (absolute_error > result->max_absolute_error) {
          result->max_absolute_error = absolute_error;
          result->worst_block = block;
          result->worst_row = row;
          result->worst_column = column;
          result->worst_analytical = analytical;
          result->worst_finite_difference = fd;
        }
        if (relative_error > result->max_relative_error)
          result->max_relative_error = relative_error;
        if (absolute_error > absolute_tolerance + relative_tolerance * scale)
          result->passed = 0;
      }
    }
  }
  return 0;
}

int chem_check_jacobian_sparsity(
    const unsigned char *pattern,
    const double *finite_difference_jacobian,
    size_t block_count,
    size_t species_count,
    double threshold,
    chem_jacobian_comparison *result) {
  if (!pattern || !finite_difference_jacobian || !result || !block_count || !species_count ||
      species_count > SIZE_MAX / species_count || block_count > SIZE_MAX / (species_count * species_count) ||
      !isfinite(threshold) || threshold < 0.0)
    return -1;
  memset(result, 0, sizeof(*result));
  result->passed = 1;
  const size_t width = species_count * species_count;
  for (size_t block = 0; block < block_count; ++block) {
    for (size_t row = 0; row < species_count; ++row) {
      for (size_t column = 0; column < species_count; ++column) {
        const size_t index = row * species_count + column;
        if (pattern[index]) continue;
        const double fd = fabs(finite_difference_jacobian[block * width + index]);
        if (!isfinite(fd)) {
          result->passed = 0;
          result->max_absolute_error = INFINITY;
          result->worst_block = block;
          result->worst_row = row;
          result->worst_column = column;
          result->worst_analytical = 0.0;
          result->worst_finite_difference = fd;
          continue;
        }
        if (fd > threshold && fd > result->max_absolute_error) {
          result->passed = 0;
          result->max_absolute_error = fd;
          result->worst_block = block;
          result->worst_row = row;
          result->worst_column = column;
          result->worst_analytical = 0.0;
          result->worst_finite_difference = fd;
        }
      }
    }
  }
  return 0;
}
