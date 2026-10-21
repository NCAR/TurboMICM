#include "chemistry.h"
#include "mechanism_internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { SOLVER_OK = CHEM_SOLVER_STATUS_OK,
       SOLVER_BAD_ARGUMENT = CHEM_SOLVER_STATUS_BAD_ARGUMENT,
       SOLVER_SINGULAR = CHEM_SOLVER_STATUS_SINGULAR,
       SOLVER_STEP_FAILURE = CHEM_SOLVER_STATUS_STEP_FAILURE,
       SOLVER_NONFINITE = CHEM_SOLVER_STATUS_NONFINITE };

const char *chem_solver_status_string(chem_solver_status status) {
  switch (status) {
  case CHEM_SOLVER_STATUS_OK: return "Converged";
  case CHEM_SOLVER_STATUS_BAD_ARGUMENT: return "Bad Argument";
  case CHEM_SOLVER_STATUS_SINGULAR: return "Singular Matrix";
  case CHEM_SOLVER_STATUS_STEP_FAILURE: return "Step Failure";
  case CHEM_SOLVER_STATUS_NONFINITE: return "Non-finite Result";
  case CHEM_SOLVER_STATUS_CONSTRAINT_INITIALIZATION_FAILED:
    return "Constraint Initialization Failed";
  case CHEM_SOLVER_STATUS_ACCEPTING_UNCONVERGED:
    return "Accepting Unconverged Integration";
  }
  return "Unknown";
}

typedef struct {
  const chem_mechanism *mechanism;
  size_t n, capacity;
  double *jac, *a, *forcing, *stage, *ynew, *yerror, *stage_y, *reused_forcing;
  double *rate_constants, *cached_parameters;
  double cached_temperature, cached_pressure, cached_density;
  int rate_constants_valid;
  unsigned char *pattern;
  size_t *row_counts, *col_counts, *row_cols, *col_rows;
  double *group_ynew, *group_yerror, *group_rate_constants, *group_cached_parameters;
  double *group_temperature, *group_pressure, *group_density;
  unsigned char *group_rate_constants_valid;
  chem_jacobian_plan jacobian_plan;
} solver_workspace_data;

static int make_pattern(const chem_mechanism *m, unsigned char *pattern,
                        const chem_constraint_set *constraints, size_t n,
                        size_t *row_counts, size_t *col_counts,
                        size_t *row_cols, size_t *col_rows);

static void free_workspace_data(solver_workspace_data *w) {
  if (!w) return;
  free(w->jac); free(w->a); free(w->forcing); free(w->stage); free(w->ynew);
  free(w->yerror); free(w->stage_y); free(w->reused_forcing); free(w->pattern);
  free(w->rate_constants); free(w->cached_parameters);
  free(w->row_counts); free(w->col_counts); free(w->row_cols); free(w->col_rows);
  free(w->group_ynew); free(w->group_yerror); free(w->group_rate_constants);
  free(w->group_cached_parameters); free(w->group_temperature);
  free(w->group_pressure); free(w->group_density);
  free(w->group_rate_constants_valid);
  chem_jacobian_plan_destroy(&w->jacobian_plan);
  free(w);
}

int chem_solver_workspace_init(chem_solver_workspace *workspace,
                               const chem_mechanism *mechanism) {
  return chem_solver_workspace_init_group(workspace, mechanism, 1);
}

int chem_solver_workspace_init_group(chem_solver_workspace *workspace,
                                     const chem_mechanism *mechanism,
                                     size_t cell_capacity) {
  if (!workspace || !mechanism || !mechanism->species || !mechanism->species_count)
    return SOLVER_BAD_ARGUMENT;
  if (!cell_capacity || cell_capacity > (size_t)-1 / mechanism->species_count ||
      cell_capacity > (size_t)-1 / (mechanism->reaction_count ? mechanism->reaction_count : 1) ||
      cell_capacity > (size_t)-1 / (mechanism->custom_parameter_count ? mechanism->custom_parameter_count : 1) ||
      cell_capacity*mechanism->species_count > (size_t)-1 / sizeof(double) ||
      cell_capacity*(mechanism->reaction_count ? mechanism->reaction_count : 1) > (size_t)-1 / sizeof(double) ||
      cell_capacity*(mechanism->custom_parameter_count ? mechanism->custom_parameter_count : 1) > (size_t)-1 / sizeof(double))
    return SOLVER_BAD_ARGUMENT;
  workspace->implementation = NULL;
  const size_t n = mechanism->species_count;
  if (n > (size_t)-1 / n || n*n > (size_t)-1 / sizeof(size_t))
    return SOLVER_BAD_ARGUMENT;
  const size_t nn = n*n;
  solver_workspace_data *w = calloc(1, sizeof(*w));
  if (!w) return SOLVER_BAD_ARGUMENT;
  w->mechanism = mechanism;
  w->n = n;
  w->capacity = cell_capacity;
  w->jac = calloc(nn, sizeof(double));
  w->a = calloc(nn, sizeof(double));
  w->forcing = calloc(n, sizeof(double));
  w->stage = calloc(6*n, sizeof(double));
  w->ynew = calloc(n, sizeof(double));
  w->yerror = calloc(n, sizeof(double));
  w->stage_y = calloc(n, sizeof(double));
  w->reused_forcing = calloc(n, sizeof(double));
  w->rate_constants = calloc(mechanism->reaction_count ? mechanism->reaction_count : 1,
                             sizeof(double));
  w->cached_parameters = calloc(mechanism->custom_parameter_count ? mechanism->custom_parameter_count : 1,
                                sizeof(double));
  w->pattern = calloc(nn, sizeof(unsigned char));
  w->row_counts = calloc(n, sizeof(size_t));
  w->col_counts = calloc(n, sizeof(size_t));
  w->row_cols = calloc(nn, sizeof(size_t));
  w->col_rows = calloc(nn, sizeof(size_t));
  w->group_ynew = calloc(cell_capacity*n, sizeof(double));
  w->group_yerror = calloc(cell_capacity*n, sizeof(double));
  w->group_rate_constants = calloc(cell_capacity*(mechanism->reaction_count ? mechanism->reaction_count : 1), sizeof(double));
  w->group_cached_parameters = calloc(cell_capacity*(mechanism->custom_parameter_count ? mechanism->custom_parameter_count : 1), sizeof(double));
  w->group_temperature = calloc(cell_capacity, sizeof(double));
  w->group_pressure = calloc(cell_capacity, sizeof(double));
  w->group_density = calloc(cell_capacity, sizeof(double));
  w->group_rate_constants_valid = calloc(cell_capacity, sizeof(unsigned char));
  if (!w->jac || !w->a || !w->forcing || !w->stage || !w->ynew || !w->yerror ||
      !w->stage_y || !w->reused_forcing || !w->pattern || !w->row_counts ||
      !w->col_counts || !w->row_cols || !w->col_rows || !w->rate_constants ||
      !w->cached_parameters || !w->group_ynew || !w->group_yerror || !w->group_rate_constants ||
      !w->group_cached_parameters || !w->group_temperature || !w->group_pressure ||
      !w->group_density || !w->group_rate_constants_valid) {
    free_workspace_data(w);
    return SOLVER_BAD_ARGUMENT;
  }
  if (chem_jacobian_plan_init(&w->jacobian_plan, mechanism)) {
    free_workspace_data(w);
    return SOLVER_BAD_ARGUMENT;
  }
  make_pattern(mechanism, w->pattern, NULL, n, w->row_counts, w->col_counts,
               w->row_cols, w->col_rows);
  workspace->implementation = w;
  return SOLVER_OK;
}

void chem_solver_workspace_destroy(chem_solver_workspace *workspace) {
  if (!workspace) return;
  free_workspace_data((solver_workspace_data *)workspace->implementation);
  workspace->implementation = NULL;
}

chem_solver_options chem_solver_default_options(void) {
  chem_solver_options result = {
      .relative_tolerance=1.0e-6, .absolute_tolerance=1.0e-3,
      .absolute_tolerances=NULL, .absolute_tolerance_count=0, .max_steps=1000,
      .initial_step=0.0};
  return result;
}

static int solver_tolerances_valid(const chem_solver_options *options,
                                   size_t species_count) {
  if (options->relative_tolerance <= 0.0 ||
      options->absolute_tolerance <= 0.0 || options->max_steps == 0 ||
      options->initial_step < 0.0 || !isfinite(options->initial_step))
    return 0;
  if ((options->absolute_tolerances &&
       options->absolute_tolerance_count != species_count) ||
      (!options->absolute_tolerances && options->absolute_tolerance_count != 0))
    return 0;
  for (size_t i = 0; options->absolute_tolerances && i < species_count; ++i)
    if (!(options->absolute_tolerances[i] > 0.0) ||
        !isfinite(options->absolute_tolerances[i])) return 0;
  return 1;
}

static double solver_absolute_tolerance(const chem_solver_options *options,
                                        size_t species) {
  return options->absolute_tolerances
      ? options->absolute_tolerances[species] : options->absolute_tolerance;
}

chem_rosenbrock_parameters chem_rosenbrock_builtin_parameters(
    chem_rosenbrock_parameter_set parameter_set) {
  chem_rosenbrock_parameters p = {0};
  if (parameter_set == CHEM_ROSENBROCK_TWO_STAGE) {
    const double g = 1.0 + 1.0 / sqrt(2.0);
    p.stages = 2;
    p.a[0] = 1.0/g;
    p.c[0] = -2.0/g;
    p.m[0] = 3.0/(2.0*g); p.m[1] = 1.0/(2.0*g);
    p.e[0] = 1.0/(2.0*g); p.e[1] = 1.0/(2.0*g);
    p.gamma[0] = g; p.gamma[1] = -g;
    p.new_function_evaluation[0] = 1; p.new_function_evaluation[1] = 1;
    p.estimator_order = 2.0;
  } else if (parameter_set == CHEM_ROSENBROCK_FOUR_STAGE) {
    p.stages = 4;
    p.a[0] = 2.0; p.a[1] = 1.867943637803922;
    p.a[2] = 0.2344449711399156; p.a[3] = p.a[1];
    p.a[4] = p.a[2];
    p.c[0] = -7.137615036412310; p.c[1] = 2.580708087951457;
    p.c[2] = 0.6515950076447975; p.c[3] = -2.137148994382534;
    p.c[4] = -0.3214669691237626; p.c[5] = -0.6949742501781779;
    p.m[0] = 2.255570073418735; p.m[1] = 0.2870493262186792;
    p.m[2] = 0.4353179431840180; p.m[3] = 1.093502252409163;
    p.e[0] = -0.2815431932141155; p.e[1] = -0.07276199124938920;
    p.e[2] = -0.1082196201495311; p.e[3] = -1.093502252409163;
    p.gamma[0] = 0.572820; p.gamma[1] = -1.769193891319233;
    p.gamma[2] = 0.7592633437920482; p.gamma[3] = -0.1049021087100450;
    p.new_function_evaluation[0] = 1; p.new_function_evaluation[1] = 1;
    p.new_function_evaluation[2] = 1; p.new_function_evaluation[3] = 0;
    p.estimator_order = 4.0;
  } else if (parameter_set == CHEM_ROSENBROCK_FOUR_STAGE_DAE) {
    p.stages = 4;
    p.a[0] = 0.0; p.a[1] = 2.0; p.a[2] = 0.0;
    p.a[3] = 2.0; p.a[4] = 0.0; p.a[5] = 1.0;
    p.c[0] = 4.0; p.c[1] = 1.0; p.c[2] = -1.0;
    p.c[3] = 1.0; p.c[4] = -1.0; p.c[5] = -8.0/3.0;
    p.m[0] = 2.0; p.m[2] = 1.0; p.m[3] = 1.0;
    p.e[3] = 1.0;
    p.gamma[0] = 0.5; p.gamma[1] = 1.5;
    p.new_function_evaluation[0] = 1; p.new_function_evaluation[1] = 0;
    p.new_function_evaluation[2] = 1; p.new_function_evaluation[3] = 1;
    p.estimator_order = 3.0;
  } else if (parameter_set == CHEM_ROSENBROCK_SIX_STAGE_DAE) {
    p.stages = 6;
    p.a[0] = 1.544; p.a[1] = 0.9466785280815826;
    p.a[2] = 0.2557011698983284; p.a[3] = 3.314825187068521;
    p.a[4] = 2.896124015972201; p.a[5] = 0.9986419139977817;
    p.a[6] = 1.221224509226641; p.a[7] = 6.019134481288629;
    p.a[8] = 12.53708332932087; p.a[9] = -0.6878860361058950;
    p.a[10] = p.a[6]; p.a[11] = p.a[7]; p.a[12] = p.a[8];
    p.a[13] = p.a[9]; p.a[14] = 1.0;
    p.c[0] = -5.6688; p.c[1] = -2.430093356833875;
    p.c[2] = -0.2063599157091915; p.c[3] = -0.1073529058151375;
    p.c[4] = -9.594562251023355; p.c[5] = -20.47028614809616;
    p.c[6] = 7.496443313967647; p.c[7] = -10.24680431464352;
    p.c[8] = -33.99990352819905; p.c[9] = 11.70890893206160;
    p.c[10] = 8.083246795921522; p.c[11] = -7.981132988064893;
    p.c[12] = -31.52159432874371; p.c[13] = 16.31930543123136;
    p.c[14] = -6.058818238834054;
    p.m[0] = p.a[6]; p.m[1] = p.a[7]; p.m[2] = p.a[8];
    p.m[3] = p.a[9]; p.m[4] = 1.0; p.m[5] = 1.0;
    p.e[5] = 1.0;
    p.gamma[0] = 0.25; p.gamma[1] = -0.1043; p.gamma[2] = 0.1035;
    p.gamma[3] = -0.0362;
    for (size_t i = 0; i < 6; ++i) p.new_function_evaluation[i] = 1;
    p.estimator_order = 4.0;
  } else {
    p.stages = 3;
    p.a[0] = 1.0; p.a[1] = 1.0; p.a[2] = 0.0;
    p.c[0] = -1.0156171083877702091975600115545;
    p.c[1] = 4.0759956452537699824805835358067;
    p.c[2] = 9.2076794298330791242156818474003;
    p.m[0] = 1.0; p.m[1] = 6.1697947043828245592553615689730;
    p.m[2] = -0.42772256543218573326238373806514;
    p.e[0] = 0.5; p.e[1] = -2.9079558716805469821718236208017;
    p.e[2] = 0.22354069897811569627360909276199;
    p.gamma[0] = 0.43586652150845899941601945119356;
    p.new_function_evaluation[0] = 1; p.new_function_evaluation[1] = 1;
    p.new_function_evaluation[2] = 0;
    p.estimator_order = 3.0;
  }
  return p;
}

static int make_pattern(const chem_mechanism *m, unsigned char *pattern,
                        const chem_constraint_set *constraints, size_t n,
                        size_t *row_counts, size_t *col_counts,
                        size_t *row_cols, size_t *col_rows) {
  memset(pattern, 0, n * n);
  for (size_t i = 0; i < n; ++i) pattern[i*n + i] = 1;
  for (size_t r = 0; r < m->reaction_count; ++r) {
    const chem_reaction *rx = &m->reactions[r];
    for (size_t i = 0; i < rx->reactant_count; ++i) {
      size_t col = rx->reactants[i].species;
      if (col >= n || m->species[col].third_body ||
          m->species[col].density_parameterized) continue;
      for (size_t j = 0; j < rx->reactant_count; ++j) {
        size_t row = rx->reactants[j].species;
        if (row < n && !m->species[row].third_body) pattern[row*n + col] = 1;
      }
      for (size_t j = 0; j < rx->product_count; ++j) {
        size_t row = rx->products[j].species;
        if (row < n && !m->species[row].third_body) pattern[row*n + col] = 1;
      }
    }
  }
  if (constraints) {
    unsigned char *constraint_pattern = calloc(n*n, sizeof(*constraint_pattern));
    if (constraint_pattern == NULL) return SOLVER_BAD_ARGUMENT;
    if (chem_constraint_set_jacobian_pattern(constraints, constraint_pattern)) {
      free(constraint_pattern);
      return SOLVER_BAD_ARGUMENT;
    }
    for (size_t i = 0; i < n*n; ++i)
      pattern[i] = (unsigned char)(pattern[i] || constraint_pattern[i]);
    free(constraint_pattern);
  }
  /* Symbolic fill for the same natural-order, no-pivot factorization used below. */
  for (size_t k = 0; k < n; ++k)
    for (size_t i = k + 1; i < n; ++i)
      if (pattern[i*n + k])
        for (size_t j = k + 1; j < n; ++j)
          if (pattern[k*n + j]) pattern[i*n + j] = 1;

  memset(row_counts, 0, n * sizeof(*row_counts));
  memset(col_counts, 0, n * sizeof(*col_counts));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      if (pattern[i*n + j]) {
        row_cols[i*n + row_counts[i]++] = j;
        col_rows[j*n + col_counts[j]++] = i;
      }
  return 0;
}

static int factor(double *a, size_t n, const size_t *row_counts,
                  const size_t *col_counts, const size_t *row_cols,
                  const size_t *col_rows) {
  for (size_t k = 0; k < n; ++k) {
    double diagonal = a[k*n + k];
    if (diagonal == 0.0 || !isfinite(diagonal)) return SOLVER_SINGULAR;
    const double inverse = 1.0 / diagonal;
    for (size_t q = 0; q < col_counts[k]; ++q) {
      size_t i = col_rows[k*n + q];
      if (i > k) a[i*n + k] *= inverse;
    }
    for (size_t p = 0; p < col_counts[k]; ++p) {
      const size_t i = col_rows[k*n + p];
      if (i <= k) continue;
      const double aik = a[i*n + k];
      #pragma omp simd
      for (size_t q = 0; q < row_counts[k]; ++q) {
        const size_t j = row_cols[k*n + q];
        if (j > k) a[i*n + j] -= aik * a[k*n + j];
      }
    }
  }
  return SOLVER_OK;
}

static void solve_factored(const double *a, size_t n, const size_t *row_counts,
                           const size_t *row_cols, double *b) {
  for (size_t i = 0; i < n; ++i)
    for (size_t q = 0; q < row_counts[i]; ++q) {
      size_t j = row_cols[i*n + q];
      if (j < i) b[i] -= a[i*n + j] * b[j];
    }
  for (size_t ii = n; ii-- > 0;) {
    for (size_t q = 0; q < row_counts[ii]; ++q) {
      size_t j = row_cols[ii*n + q];
      if (j > ii) b[ii] -= a[ii*n + j] * b[j];
    }
    b[ii] /= a[ii*n + ii];
  }
}

static int row_is_algebraic(const chem_constraint_set *constraints, size_t row) {
  if (constraints == NULL) return 0;
  for (size_t i = 0; i < constraints->constraint_count; ++i)
    if (constraints->constraints[i].algebraic_species == row) return 1;
  return 0;
}

static int project_constraints(const chem_constraint_set *constraints,
                               double *state, double temperature,
                               double *residual, double *negative_jacobian,
                               double relative_tolerance,
                               double absolute_tolerance, size_t *iterations) {
  if (constraints == NULL || constraints->constraint_count == 0) return 0;
  if (iterations) *iterations = 0;
  const size_t n = constraints->species_count;
  const size_t count = constraints->constraint_count;
  if (count > n || n > (size_t)-1/n || count > (size_t)-1/count) return 1;
  size_t *variables = malloc(count*sizeof(*variables));
  double *original = malloc(count*sizeof(*original));
  double *matrix = malloc(count*count*sizeof(*matrix));
  double *row_norms = malloc(count*sizeof(*row_norms));
  double *rhs = malloc(count*sizeof(*rhs));
  double *trial = malloc(n*sizeof(*trial));
  double *trial_residual = malloc(n*sizeof(*trial_residual));
  if (!variables || !original || !matrix || !row_norms || !rhs ||
      !trial || !trial_residual) {
    free(variables); free(original); free(matrix); free(row_norms); free(rhs);
    free(trial); free(trial_residual);
    return 1;
  }
  for (size_t i = 0; i < count; ++i) {
    variables[i] = constraints->constraints[i].algebraic_species;
    for (size_t j = 0; j < i; ++j)
      if (variables[j] == variables[i]) {
        free(variables); free(original); free(matrix); free(row_norms); free(rhs);
        free(trial); free(trial_residual);
        return 1;
      }
    original[i] = state[variables[i]];
  }
  int result = 1;
  for (size_t iteration = 0; iteration < 30; ++iteration) {
    memset(residual, 0, n*sizeof(*residual));
    memset(negative_jacobian, 0, n*n*sizeof(*negative_jacobian));
    if (chem_constraint_set_add_forcing(constraints, state, &temperature, 1,
                                        residual) ||
        chem_constraint_set_add_negative_jacobian(
            constraints, state, &temperature, 1, negative_jacobian)) break;
    double maximum_residual = 0.0;
    for (size_t i = 0; i < count; ++i) {
      const size_t row = variables[i];
      if (fabs(residual[row]) > maximum_residual)
        maximum_residual = fabs(residual[row]);
      rhs[i] = -residual[row];
      row_norms[i] = 0.0;
      for (size_t j = 0; j < count; ++j) {
        matrix[i*count+j] = -negative_jacobian[row*n + variables[j]];
        if (fabs(matrix[i*count+j]) > row_norms[i])
          row_norms[i] = fabs(matrix[i*count+j]);
      }
    }
    if (maximum_residual == 0.0) {
      if (iterations) *iterations = iteration + 1;
      result = 0;
      break;
    }
    for (size_t k = 0; k < count; ++k) {
      size_t pivot = k;
      for (size_t i = k+1; i < count; ++i)
        if (fabs(matrix[i*count+k]) > fabs(matrix[pivot*count+k])) pivot = i;
      if (fabs(matrix[pivot*count+k]) < 1.0e-30) goto done;
      if (pivot != k) {
        for (size_t j = k; j < count; ++j) {
          const double tmp = matrix[k*count+j];
          matrix[k*count+j] = matrix[pivot*count+j];
          matrix[pivot*count+j] = tmp;
        }
        const double tmp = rhs[k]; rhs[k] = rhs[pivot]; rhs[pivot] = tmp;
      }
      for (size_t i = k+1; i < count; ++i) {
        const double factor = matrix[i*count+k]/matrix[k*count+k];
        for (size_t j = k+1; j < count; ++j)
          matrix[i*count+j] -= factor*matrix[k*count+j];
        rhs[i] -= factor*rhs[k];
      }
    }
    for (size_t ii = count; ii-- > 0;) {
      for (size_t j = ii+1; j < count; ++j) rhs[ii] -= matrix[ii*count+j]*rhs[j];
      rhs[ii] /= matrix[ii*count+ii];
    }
    int correction_is_small = 1;
    for (size_t i = 0; i < count; ++i) {
      if (!isfinite(rhs[i])) goto done;
      const double tolerance = 0.1 * (absolute_tolerance +
          relative_tolerance * fabs(state[variables[i]]));
      if (fabs(rhs[i]) > tolerance) correction_is_small = 0;
    }
    if (correction_is_small) {
      if (iterations) *iterations = iteration + 1;
      result = 0;
      break;
    }
    double current_error = 0.0;
    for (size_t i = 0; i < count; ++i) {
      if (row_norms[i] < 1.0e-30) goto done;
      const double normalized = fabs(residual[variables[i]])/row_norms[i];
      if (normalized > current_error) current_error = normalized;
    }
    double damping = 1.0;
    int accepted = 0;
    for (size_t backtrack = 0; backtrack < 40; ++backtrack, damping *= 0.5) {
      memcpy(trial, state, n*sizeof(*trial));
      int valid = 1;
      for (size_t i = 0; i < count; ++i) {
        trial[variables[i]] += damping*rhs[i];
        if (!isfinite(trial[variables[i]]) || trial[variables[i]] < 0.0)
          valid = 0;
      }
      if (!valid) continue;
      memset(trial_residual, 0, n*sizeof(*trial_residual));
      if (chem_constraint_set_add_forcing(constraints, trial, &temperature, 1,
                                          trial_residual)) continue;
      double trial_error = 0.0;
      for (size_t i = 0; i < count; ++i) {
        const double normalized = fabs(trial_residual[variables[i]])/row_norms[i];
        if (normalized > trial_error) trial_error = normalized;
      }
      if (trial_error < current_error) {
        memcpy(state, trial, n*sizeof(*state));
        accepted = 1;
        break;
      }
    }
    if (!accepted) goto done;
  }
done:
  if (result != 0)
    for (size_t i = 0; i < count; ++i) state[variables[i]] = original[i];
  free(variables); free(original); free(matrix); free(row_norms); free(rhs);
  free(trial); free(trial_residual);
  return result;
}

static int rosenbrock_solve_impl(
    const chem_mechanism *m, chem_solver_workspace *workspace, double *y,
    double temperature, double pressure, double air_density,
    const double *rate_parameters,
    double duration, const chem_rosenbrock_parameters *parameters,
    const chem_constraint_set *constraints,
    const chem_solver_options *requested,
    chem_solver_stats *statistics) {
  if (!m || !y || !isfinite(duration) || duration <= 0.0 ||
      !isfinite(temperature) || temperature <= 0.0 ||
      !isfinite(pressure) || pressure < 0.0 ||
      !isfinite(air_density) || air_density < 0.0 || m->species_count == 0)
    return SOLVER_BAD_ARGUMENT;
  if (constraints && (chem_constraint_set_validate(constraints) ||
                      constraints->species_count != m->species_count))
    return SOLVER_BAD_ARGUMENT;
  solver_workspace_data *w = workspace ? workspace->implementation : NULL;
  if (!w || w->mechanism != m || w->n != m->species_count)
    return SOLVER_BAD_ARGUMENT;
  const size_t n = m->species_count;
  const chem_solver_options defaults = chem_solver_default_options();
  const chem_solver_options *o = requested ? requested : &defaults;
  if (!solver_tolerances_valid(o, n)) return SOLVER_BAD_ARGUMENT;
  const chem_rosenbrock_parameters default_parameters =
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  const chem_rosenbrock_parameters *p = parameters ? parameters : &default_parameters;
  if (p->stages < 2 || p->stages > 6 || p->gamma[0] == 0.0 ||
      p->estimator_order <= 0.0) return SOLVER_BAD_ARGUMENT;
  chem_solver_stats local_stats = {0};
  if (statistics) *statistics = local_stats;
  double *jac = w->jac, *a = w->a, *forcing = w->forcing;
  double *stage = w->stage, *ynew = w->ynew, *yerror = w->yerror;
  double *stage_y = w->stage_y, *reused_forcing = w->reused_forcing;
  size_t *row_counts = w->row_counts, *col_counts = w->col_counts;
  size_t *row_cols = w->row_cols, *col_rows = w->col_rows;
  int status = SOLVER_OK;
  if (make_pattern(m, w->pattern, constraints, n, row_counts, col_counts,
                   row_cols, col_rows) != SOLVER_OK) return SOLVER_BAD_ARGUMENT;

  int same_parameters = m->custom_parameter_count == 0 ||
      (rate_parameters && memcmp(w->cached_parameters, rate_parameters,
          m->custom_parameter_count*sizeof(double)) == 0);
  if (!w->rate_constants_valid || w->cached_temperature != temperature ||
      w->cached_pressure != pressure ||
      w->cached_density != air_density || !same_parameters) {
    status = chem_rate_constants(m, temperature, pressure, air_density, rate_parameters,
                                 w->rate_constants);
    if (status) return SOLVER_BAD_ARGUMENT;
    if (m->custom_parameter_count && rate_parameters)
      memcpy(w->cached_parameters, rate_parameters,
             m->custom_parameter_count*sizeof(double));
    w->cached_temperature = temperature;
    w->cached_pressure = pressure;
    w->cached_density = air_density;
    w->rate_constants_valid = 1;
  }
  if (constraints && project_constraints(constraints, y, temperature,
                                         stage_y, jac, o->relative_tolerance,
                                         o->absolute_tolerance,
                                         &local_stats.constraint_init_iterations))
    return CHEM_SOLVER_STATUS_CONSTRAINT_INITIALIZATION_FAILED;

  const double eps = 2.2204460492503131e-16;
  const double hmin = 1.0e-15 * duration;
  const double hmax = duration;
  double h = o->initial_step > 0.0 ? fmin(o->initial_step,duration)
                                    : 1.0e-6 * duration;
  double time = 0.0;
  int rejected_previous = 0, rejected_twice = 0;
  while (time - duration + eps <= 0.0) {
    if (local_stats.internal_steps > o->max_steps) { status = SOLVER_STEP_FAILURE; break; }
    if (time + 0.1*h == time || h <= fmin(eps, hmin)) { status = SOLVER_STEP_FAILURE; break; }
    if (h > duration-time) h = duration-time;
    memset(forcing, 0, n*sizeof(*forcing));
    status = chem_forcing_with_constants(m, y, w->rate_constants,
                                         air_density, forcing);
    if (status) { status = SOLVER_BAD_ARGUMENT; break; }
    if (constraints && chem_constraint_set_add_forcing(
            constraints, y, &temperature, 1, forcing)) {
      status = SOLVER_BAD_ARGUMENT; break;
    }
    ++local_stats.function_calls;
    status = chem_negative_jacobian_with_plan_pattern(
        m, y, w->rate_constants, &w->jacobian_plan, row_counts, row_cols,
        jac);
    if (status) { status = SOLVER_BAD_ARGUMENT; break; }
    if (constraints && chem_constraint_set_add_negative_jacobian(
            constraints, y, &temperature, 1, jac)) {
      status = SOLVER_BAD_ARGUMENT; break;
    }
    ++local_stats.jacobian_updates;
    int accepted = 0;
    while (!accepted) {
      const double used_h = h;
      const double alpha = 1.0 / (h*p->gamma[0]);
      for (size_t row = 0; row < n; ++row)
        for (size_t q = 0; q < row_counts[row]; ++q) {
          const size_t column = row_cols[row*n + q];
          a[row*n + column] = jac[row*n + column];
        }
      for (size_t i = 0; i < n; ++i)
        if (!row_is_algebraic(constraints, i)) a[i*n + i] += alpha;
      status = factor(a, n, row_counts, col_counts, row_cols, col_rows);
      ++local_stats.decompositions;
      if (status) break;
      for (size_t i = 0; i < n; ++i) stage[i] = forcing[i];
      for (size_t s = 0; s < p->stages; ++s) {
        double *ks = stage + s*n;
        const size_t coefficient_offset = s*(s-1)/2;
        if (s > 0 && p->new_function_evaluation[s]) {
          for (size_t i = 0; i < n; ++i) {
            stage_y[i] = y[i];
            for (size_t j = 0; j < s; ++j)
              stage_y[i] += p->a[coefficient_offset+j]*stage[j*n+i];
          }
          memset(ks, 0, n*sizeof(*ks));
          status = chem_forcing_with_constants(m, stage_y, w->rate_constants,
                                               air_density, ks);
          ++local_stats.function_calls;
          if (status) { status = SOLVER_BAD_ARGUMENT; break; }
          if (constraints && chem_constraint_set_add_forcing(
                  constraints, stage_y, &temperature, 1, ks)) {
            status = SOLVER_BAD_ARGUMENT; break;
          }
        } else if (s > 0) {
          memcpy(ks, reused_forcing, n*sizeof(*ks));
        }
        if (s + 1 < p->stages && !p->new_function_evaluation[s+1])
          memcpy(reused_forcing, ks, n*sizeof(*reused_forcing));
        for (size_t j = 0; j < s; ++j) {
          const double scale = p->c[coefficient_offset+j] / h;
          for (size_t i = 0; i < n; ++i)
            if (!row_is_algebraic(constraints, i))
              ks[i] += scale*stage[j*n + i];
        }
        solve_factored(a, n, row_counts, row_cols, ks);
        ++local_stats.solves;
      }
      if (status) break;
      for (size_t i = 0; i < n; ++i) {
        ynew[i] = y[i];
        yerror[i] = 0.0;
        for (size_t s = 0; s < p->stages; ++s) {
          ynew[i] += p->m[s]*stage[s*n + i];
          yerror[i] += p->e[s]*stage[s*n + i];
        }
      }
      if (constraints && project_constraints(constraints, ynew, temperature,
                                             stage_y, jac, 1.0e-12,
                                             1.0e-12, NULL)) {
        status = SOLVER_STEP_FAILURE;
        break;
      }
      double sum = 0.0;
      size_t error_variables = 0;
      for (size_t i = 0; i < n; ++i) {
        if (m->species[i].third_body) continue;
        ++error_variables;
        const double scale = solver_absolute_tolerance(o, i) +
            o->relative_tolerance*fmax(fabs(y[i]), fabs(ynew[i]));
        const double estimate = row_is_algebraic(constraints, i)
            ? ynew[i] - y[i] : yerror[i];
        const double v = estimate/scale;
        sum += v*v;
      }
      double error = error_variables
          ? fmax(sqrt(sum/(double)error_variables), 1.0e-10) : 1.0e-10;
      double fac = fmin(6.0, fmax(0.2, 0.9/pow(error, 1.0/p->estimator_order)));
      double hnew = h*fac;
      ++local_stats.internal_steps;
      if (!isfinite(error)) { status = SOLVER_NONFINITE; break; }
      if (error < 1.0 || h < hmin) {
        ++local_stats.accepted;
        time += h;
        memcpy(y, ynew, n*sizeof(*y));
        double next_h = fmax(hmin, fmin(hnew, hmax));
        if (rejected_previous) next_h = fmin(next_h, used_h);
        h = next_h;
        rejected_previous = 0;
        rejected_twice = 0;
        accepted = 1;
      } else {
        if (rejected_twice) hnew = h*0.1;
        rejected_twice = rejected_previous;
        rejected_previous = 1;
        h = hnew;
        ++local_stats.rejected;
        ++local_stats.jacobian_updates;
      }
    }
    if (status) break;
  }
  for (size_t i = 0; i < n; ++i)
    if (!m->species[i].third_body && !row_is_algebraic(constraints, i) &&
        y[i] < 0.0) y[i] = 0.0;
  if (!status && constraints && project_constraints(
          constraints, y, temperature, stage_y, jac, 1.0e-12, 1.0e-12,
          NULL))
    status = SOLVER_STEP_FAILURE;
  if (statistics) *statistics = local_stats;
  return status;
}

int chem_rosenbrock_solve_with_parameters(
    const chem_mechanism *m, chem_solver_workspace *workspace, double *y,
    double temperature, double pressure, double air_density,
    const double *rate_parameters, double duration,
    const chem_rosenbrock_parameters *parameters,
    const chem_solver_options *requested, chem_solver_stats *statistics) {
  return rosenbrock_solve_impl(m, workspace, y, temperature, pressure,
                               air_density, rate_parameters, duration,
                               parameters, NULL, requested, statistics);
}

int chem_rosenbrock_solve_constrained_with_parameters(
    const chem_mechanism *m, chem_solver_workspace *workspace, double *y,
    double temperature, double pressure, double air_density,
    const double *rate_parameters, double duration,
    const chem_rosenbrock_parameters *parameters,
    const chem_constraint_set *constraints,
    const chem_solver_options *requested, chem_solver_stats *statistics) {
  if (constraints == NULL) return SOLVER_BAD_ARGUMENT;
  return rosenbrock_solve_impl(m, workspace, y, temperature, pressure,
                               air_density, rate_parameters, duration,
                               parameters, constraints, requested, statistics);
}

int chem_rosenbrock_solve_constrained_group_with_parameters(
    const chem_mechanism *m, chem_solver_workspace *workspace, double *states,
    const double *temperatures, const double *pressures,
    const double *air_densities, const double *rate_parameters,
    size_t cell_count, double duration,
    const chem_rosenbrock_parameters *parameters,
    const chem_constraint_set *constraints,
    const chem_solver_options *requested, chem_solver_stats *statistics) {
  solver_workspace_data *w = workspace ? workspace->implementation : NULL;
  if (!m || !w || !states || !temperatures || !pressures || !air_densities ||
      cell_count > w->capacity || w->mechanism != m ||
      (m->custom_parameter_count && rate_parameters == NULL))
    return SOLVER_BAD_ARGUMENT;
  chem_solver_stats total = {0};
  for (size_t cell = 0; cell < cell_count; ++cell) {
    chem_solver_stats one = {0};
    const double *parameters_for_cell = rate_parameters
        ? rate_parameters + cell*m->custom_parameter_count : NULL;
    const int status = rosenbrock_solve_impl(
        m, workspace, states + cell*m->species_count, temperatures[cell],
        pressures[cell], air_densities[cell], parameters_for_cell, duration,
        parameters, constraints, requested, &one);
    total.internal_steps += one.internal_steps;
    total.accepted += one.accepted;
    total.rejected += one.rejected;
    total.function_calls += one.function_calls;
    total.jacobian_updates += one.jacobian_updates;
    total.decompositions += one.decompositions;
    total.solves += one.solves;
    total.constraint_init_iterations += one.constraint_init_iterations;
    if (status) {
      if (statistics) *statistics = total;
      return status;
    }
  }
  if (statistics) *statistics = total;
  return SOLVER_OK;
}

int chem_rosenbrock3_solve_with_workspace(
    const chem_mechanism *m, chem_solver_workspace *workspace, double *y,
    double temperature, double pressure, double air_density,
    const double *rate_parameters, double duration,
    const chem_solver_options *requested, chem_solver_stats *statistics) {
  const chem_rosenbrock_parameters parameters =
      chem_rosenbrock_builtin_parameters(CHEM_ROSENBROCK_THREE_STAGE);
  return chem_rosenbrock_solve_with_parameters(
      m, workspace, y, temperature, pressure, air_density, rate_parameters,
      duration, &parameters, requested, statistics);
}

int chem_rosenbrock3_solve_group_with_workspace(
    const chem_mechanism *m, chem_solver_workspace *workspace, double *states,
    const double *temperatures, const double *pressures,
    const double *air_densities, const double *rate_parameters,
    size_t cell_count, double duration, const chem_solver_options *requested,
    chem_solver_stats *statistics) {
  if (!m || !workspace || !states || !temperatures || !pressures ||
      !air_densities || !isfinite(duration) || duration <= 0.0 ||
      cell_count == 0) return SOLVER_BAD_ARGUMENT;
  solver_workspace_data *w = workspace->implementation;
  if (!w || w->mechanism != m || w->n != m->species_count ||
      cell_count > w->capacity) return SOLVER_BAD_ARGUMENT;
  const size_t n = m->species_count, nr = m->reaction_count;
  const chem_solver_options defaults = chem_solver_default_options();
  const chem_solver_options *o = requested ? requested : &defaults;
  if (!solver_tolerances_valid(o, n)) return SOLVER_BAD_ARGUMENT;
  for (size_t c = 0; c < cell_count; ++c)
    if (!isfinite(temperatures[c]) || temperatures[c] <= 0.0 ||
        !isfinite(pressures[c]) || pressures[c] < 0.0 ||
        !isfinite(air_densities[c]) || air_densities[c] < 0.0)
      return SOLVER_BAD_ARGUMENT;

  chem_solver_stats local = {0};
  const double eps = 2.2204460492503131e-16;
  const double hmin = 1.0e-15 * duration, hmax = duration;
  double h = o->initial_step > 0.0 ? fmin(o->initial_step,duration)
                                    : 1.0e-6 * duration;
  double time = 0.0;
  int rejected_previous = 0, rejected_twice = 0, status = SOLVER_OK;
  const double aa[3] = {1.0, 1.0, 0.0};
  const double cc[3] = {-1.0156171083877702091975600115545,
                         4.0759956452537699824805835358067,
                         9.2076794298330791242156818474003};
  const double mm[3] = {1.0, 6.1697947043828245592553615689730,
                        -0.42772256543218573326238373806514};
  const double ee[3] = {0.5, -2.9079558716805469821718236208017,
                        0.22354069897811569627360909276199};
  const double gamma = 0.43586652150845899941601945119356;
  while (time - duration + eps <= 0.0) {
    if (local.internal_steps > o->max_steps) { status = SOLVER_STEP_FAILURE; break; }
    if (time + 0.1*h == time || h <= fmin(eps, hmin)) { status = SOLVER_STEP_FAILURE; break; }
    if (h > duration-time) h = duration-time;
    int accepted = 0;
    while (!accepted) {
      const double used_h = h;
      for (size_t c = 0; c < cell_count; ++c) {
        const double *params = rate_parameters
            ? rate_parameters + c*m->custom_parameter_count : NULL;
        double *cached_params = w->group_cached_parameters +
            c*(m->custom_parameter_count ? m->custom_parameter_count : 1);
        double *rates = w->group_rate_constants +
            c*(nr ? nr : 1);
        int same_params = m->custom_parameter_count == 0 ||
            (params && memcmp(cached_params, params,
                 m->custom_parameter_count*sizeof(double)) == 0);
        if (!w->group_rate_constants_valid[c] ||
            w->group_temperature[c] != temperatures[c] ||
            w->group_pressure[c] != pressures[c] ||
            w->group_density[c] != air_densities[c] || !same_params) {
          status = chem_rate_constants(m, temperatures[c], pressures[c],
                                       air_densities[c], params, rates);
          if (status) { status = SOLVER_BAD_ARGUMENT; break; }
          if (m->custom_parameter_count && params)
            memcpy(cached_params, params,
                   m->custom_parameter_count*sizeof(double));
          w->group_temperature[c] = temperatures[c];
          w->group_pressure[c] = pressures[c];
          w->group_density[c] = air_densities[c];
          w->group_rate_constants_valid[c] = 1;
        }
        const double *y = states + c*n;
        double *ynew = w->group_ynew + c*n;
        double *yerror = w->group_yerror + c*n;
        status = chem_forcing_with_constants(m, y, rates, air_densities[c],
                                             w->forcing);
        if (status) { status = SOLVER_BAD_ARGUMENT; break; }
        status = chem_negative_jacobian_with_plan_pattern(
            m, y, rates, &w->jacobian_plan, w->row_counts, w->row_cols,
            w->jac);
        if (status) { status = SOLVER_BAD_ARGUMENT; break; }
        double alpha = 1.0 / (h*gamma);
        for (size_t row = 0; row < n; ++row)
          for (size_t q = 0; q < w->row_counts[row]; ++q) {
            const size_t column = w->row_cols[row*n + q];
            w->a[row*n + column] = w->jac[row*n + column];
          }
        for (size_t i = 0; i < n; ++i) w->a[i*n + i] += alpha;
        status = factor(w->a, n, w->row_counts, w->col_counts,
                        w->row_cols, w->col_rows);
        if (status) break;
        double *stage = w->stage;
        memcpy(stage, w->forcing, n*sizeof(double));
        for (size_t s = 0; s < 3; ++s) {
          if (s == 1) {
            for (size_t i = 0; i < n; ++i)
              w->stage_y[i] = y[i] + aa[0]*stage[i];
            status = chem_forcing_with_constants(m, w->stage_y, rates,
                                                 air_densities[c], stage+n);
            if (status) { status = SOLVER_BAD_ARGUMENT; break; }
            memcpy(w->reused_forcing, stage+n, n*sizeof(double));
          } else if (s == 2) {
            memcpy(stage+2*n, w->reused_forcing, n*sizeof(double));
          }
          double *ks = stage + s*n;
          for (size_t j = 0; j < s; ++j) {
            const size_t ci = s == 1 ? 0 : (j == 0 ? 1 : 2);
            const double scale = cc[ci] / h;
            for (size_t i = 0; i < n; ++i)
              ks[i] += scale*stage[j*n + i];
          }
          solve_factored(w->a, n, w->row_counts, w->row_cols, ks);
        }
        if (status) break;
        for (size_t i = 0; i < n; ++i) {
          ynew[i] = y[i];
          yerror[i] = 0.0;
          for (size_t s = 0; s < 3; ++s) {
            ynew[i] += mm[s]*stage[s*n + i];
            yerror[i] += ee[s]*stage[s*n + i];
          }
        }
      }
      if (status) break;
      double error_sum = 0.0;
      size_t error_variables = 0;
      /* MICM reduces each species column across group rows before moving to
         the next species; preserve that summation order near controller limits. */
      for (size_t i = 0; i < n; ++i) {
        if (m->species[i].third_body) continue;
        for (size_t c = 0; c < cell_count; ++c) {
          const double *y = states + c*n;
          const double *ynew = w->group_ynew + c*n;
          const double *yerror = w->group_yerror + c*n;
          const double scale = solver_absolute_tolerance(o, i) + o->relative_tolerance *
              fmax(fabs(y[i]), fabs(ynew[i]));
          const double value = yerror[i] / scale;
          error_sum += value*value;
          ++error_variables;
        }
      }
      const double error = fmax(sqrt(error_sum/(double)error_variables), 1.0e-10);
      const double fac = fmin(6.0, fmax(0.2, 0.9/pow(error, 1.0/3.0)));
      double hnew = h*fac;
      ++local.internal_steps;
      ++local.function_calls;
      ++local.function_calls;
      ++local.jacobian_updates;
      ++local.decompositions;
      local.solves += 3;
      if (!isfinite(error)) { status = SOLVER_NONFINITE; break; }
      if (error < 1.0 || h < hmin) {
        ++local.accepted;
        time += h;
        for (size_t c = 0; c < cell_count; ++c)
          memcpy(states+c*n, w->group_ynew+c*n, n*sizeof(double));
        double next_h = fmax(hmin, fmin(hnew, hmax));
        if (rejected_previous) next_h = fmin(next_h, used_h);
        h = next_h;
        rejected_previous = 0;
        rejected_twice = 0;
        accepted = 1;
      } else {
        if (rejected_twice) hnew = h*0.1;
        rejected_twice = rejected_previous;
        rejected_previous = 1;
        h = hnew;
        ++local.rejected;
        ++local.jacobian_updates;
      }
    }
    if (status) break;
  }
  for (size_t c = 0; c < cell_count; ++c)
    for (size_t i = 0; i < n; ++i)
      if (!m->species[i].third_body && states[c*n+i] < 0.0)
        states[c*n+i] = 0.0;
  if (statistics) *statistics = local;
  return status;
}

int chem_rosenbrock3_solve(const chem_mechanism *m, double *y,
                           double temperature, double air_density,
                           const double *rate_parameters, double duration,
                           const chem_solver_options *requested,
                           chem_solver_stats *statistics) {
  chem_solver_workspace workspace = {0};
  int status = chem_solver_workspace_init(&workspace, m);
  if (status) return status;
  status = chem_rosenbrock3_solve_with_workspace(
      m, &workspace, y, temperature,
      air_density * 8.31446261815324 * temperature, air_density,
      rate_parameters, duration,
      requested, statistics);
  chem_solver_workspace_destroy(&workspace);
  return status;
}

chem_backward_euler_options chem_backward_euler_default_options(void) {
  return (chem_backward_euler_options){
      .small=1.0e-40,.initial_step=0.0,.relative_tolerance=1.0e-6,
      .absolute_tolerance=1.0e-3,.absolute_tolerances=NULL,
      .absolute_tolerance_count=0,.max_iterations=11,
      .time_step_reductions={0.5,0.5,0.5,0.5,0.1}};
}

int chem_backward_euler_solve_with_workspace(
    const chem_mechanism *m, chem_solver_workspace *workspace, double *y,
    double temperature, double pressure, double air_density,
    const double *rate_parameters, double duration,
    const chem_backward_euler_options *requested,
    chem_solver_stats *statistics) {
  solver_workspace_data *w = workspace ? workspace->implementation : NULL;
  if (!m || !w || w->mechanism != m || !y || !m->species_count ||
      !isfinite(temperature) || temperature <= 0.0 ||
      !isfinite(pressure) || pressure < 0.0 ||
      !isfinite(air_density) || air_density < 0.0 ||
      !isfinite(duration) || duration < 0.0 ||
      (m->custom_parameter_count && !rate_parameters))
    return SOLVER_BAD_ARGUMENT;
  const chem_backward_euler_options defaults =
      chem_backward_euler_default_options();
  const chem_backward_euler_options *o = requested ? requested : &defaults;
  if (!(o->small >= 0.0) || !isfinite(o->small) ||
      !(o->initial_step >= 0.0) || !isfinite(o->initial_step) ||
      !(o->relative_tolerance >= 0.0) || !isfinite(o->relative_tolerance) ||
      !(o->absolute_tolerance >= 0.0) || !isfinite(o->absolute_tolerance) ||
      !o->max_iterations) return SOLVER_BAD_ARGUMENT;
  if ((o->absolute_tolerances && o->absolute_tolerance_count != m->species_count) ||
      (!o->absolute_tolerances && o->absolute_tolerance_count != 0))
    return SOLVER_BAD_ARGUMENT;
  for (size_t i = 0; o->absolute_tolerances && i < m->species_count; ++i)
    if (!(o->absolute_tolerances[i] >= 0.0) ||
        !isfinite(o->absolute_tolerances[i])) return SOLVER_BAD_ARGUMENT;
  for (size_t i = 0; i < 5; ++i)
    if (!(o->time_step_reductions[i] > 0.0) ||
        !(o->time_step_reductions[i] < 1.0) ||
        !isfinite(o->time_step_reductions[i])) return SOLVER_BAD_ARGUMENT;
  chem_solver_stats local = {0};
  if (statistics) *statistics = local;
  if (duration == 0.0) return SOLVER_OK;
  const size_t n = m->species_count;
  int status = chem_rate_constants(m,temperature,pressure,air_density,
                                   rate_parameters,w->rate_constants);
  if (status) return SOLVER_BAD_ARGUMENT;
  double *old = w->stage_y;
  double *next = w->ynew;
  double *delta = w->forcing;
  double *negative_jacobian = w->jac;
  double *linear_matrix = w->a;
  memcpy(old,y,n*sizeof(*old));
  memcpy(next,y,n*sizeof(*next));
  double time = 0.0;
  double h = o->initial_step == 0.0 ? duration : fmin(o->initial_step,duration);
  size_t consecutive_successes = 0, failed_reductions = 0;
  while (time < duration) {
    h = fmin(h,duration-time);
    memcpy(next,old,n*sizeof(*next));
    int converged = 0;
    for (size_t iteration = 0; iteration < o->max_iterations; ++iteration) {
      status = chem_forcing_with_constants(m,next,w->rate_constants,
                                           air_density,delta);
      if (status) { status = SOLVER_BAD_ARGUMENT; break; }
      ++local.function_calls;
      status = chem_negative_jacobian_with_plan(
          m,next,w->rate_constants,&w->jacobian_plan,negative_jacobian);
      if (status) { status = SOLVER_BAD_ARGUMENT; break; }
      ++local.jacobian_updates;
      for (size_t row = 0; row < n; ++row)
        for (size_t column = 0; column < n; ++column)
          linear_matrix[row*n+column] =
              negative_jacobian[row*n+column] +
              (row == column ? 1.0/h : 0.0);
      for (size_t i = 0; i < n; ++i)
        delta[i] -= (next[i]-old[i])/h;
      if (chem_linear_solve_in_place(n,linear_matrix,delta)) {
        status = CHEM_SOLVER_STATUS_SINGULAR;
        break;
      }
      ++local.decompositions;
      ++local.solves;
      ++local.internal_steps;
      for (size_t i = 0; i < n; ++i) {
        next[i] = fmax(0.0,next[i]+delta[i]);
        if (!isfinite(next[i]) || !isfinite(delta[i])) {
          status = CHEM_SOLVER_STATUS_NONFINITE;
          break;
        }
      }
      if (status) break;
      if (iteration == 0) continue;
      converged = 1;
      for (size_t i = 0; i < n; ++i) {
        const double absolute_tolerance = o->absolute_tolerances
            ? o->absolute_tolerances[i] : o->absolute_tolerance;
        if (fabs(delta[i]) > o->small &&
            fabs(delta[i]) > absolute_tolerance &&
            fabs(delta[i]) > o->relative_tolerance*fabs(next[i])) {
          converged = 0;
          break;
        }
      }
      if (converged) break;
    }
    if (status) break;
    if (!converged) {
      ++local.rejected;
      consecutive_successes = 0;
      if (failed_reductions >= 5) {
        time += h;
        status = CHEM_SOLVER_STATUS_ACCEPTING_UNCONVERGED;
        break;
      }
      memcpy(next,old,n*sizeof(*next));
      h *= o->time_step_reductions[failed_reductions++];
      continue;
    }
    memcpy(y,next,n*sizeof(*y));
    memcpy(old,next,n*sizeof(*old));
    time += h;
    ++local.accepted;
    ++consecutive_successes;
    if (consecutive_successes >= 2) {
      consecutive_successes = 0;
      h *= 2.0;
    }
  }
  if (statistics) *statistics = local;
  return status;
}
