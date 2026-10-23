#include "kinetics.h"
#include "solver_internal.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

enum { K_WORK_BASE, YNEW_WORK = 6, ERROR_WORK, TEMP_WORK, WORK_COUNT };

#define CHEM_DEFAULT_REFERENCE_TEMPERATURE 298.15

struct chem_solver {
  const chem_mechanism *mechanism;
  size_t cells;
  size_t n;
  size_t nr;
  size_t ncustom;
  size_t groups;
  size_t group_size;
  int thread_count;
  chem_rosenbrock_parameters parameters;
  const chem_constraint_set *constraints;
  const double *equilibrium_constants;
  const double *constraint_temperatures;
  chem_sparse_pattern pattern;
  double *rates;       /* cell-major, reaction */
  double *work;        /* work-vector-major, cell-major, species */
  double *factor_work; /* one sparse LU workspace per solver group */
  int allocation_error;
};

const char *chem_solver_state_string(chem_solver_state state)
{
  switch (state) {
  case CHEM_SOLVER_CONVERGED: return "Converged";
  case CHEM_SOLVER_RUNNING: return "Running";
  case CHEM_SOLVER_FAILED: return "Failed";
  case CHEM_SOLVER_CONSTRAINT_INITIALIZATION_FAILED:
    return "Constraint Initialization Failed";
  }
  return "Unknown";
}

double chem_normalized_error(const double *old_values,
                             const double *new_values,
                             const double *estimated_errors,
                             const chem_species *species, size_t cell_count,
                             size_t species_count, double relative_tolerance,
                             double default_absolute_tolerance)
{
  size_t cell, variable;
  size_t count;
  double sum = 0.0;
  if (cell_count && species_count > SIZE_MAX / cell_count)
    return NAN;
  count = cell_count * species_count;
  if (count && (!old_values || !new_values || !estimated_errors || !species))
    return NAN;
  if (!(relative_tolerance >= 0.0) || !(default_absolute_tolerance > 0.0))
    return NAN;
  for (cell = 0; cell < cell_count; ++cell) {
    for (variable = 0; variable < species_count; ++variable) {
      const size_t index = cell * species_count + variable;
      double absolute_tolerance = species[variable].absolute_tolerance;
      const double scale = fmax(fabs(old_values[index]), fabs(new_values[index]));
      double normalized;
      if (!(absolute_tolerance > 0.0))
        absolute_tolerance = default_absolute_tolerance;
      normalized = estimated_errors[index] /
                   (absolute_tolerance + relative_tolerance * scale);
      sum += normalized * normalized;
    }
  }
  return fmax(sqrt(sum / (double)(count ? count : 1)), 1.0e-10);
}

static double *work_vector(chem_solver *s, size_t which)
{
  return s->work + which * s->cells * s->n;
}

static void set_pair(unsigned char *pattern, size_t n, size_t row, size_t column)
{
  pattern[row * n + column] = 1;
}

int chem_make_pattern(const chem_mechanism *m, chem_sparse_pattern *p)
{
  size_t n = m->species_count;
  unsigned char *matrix;
  size_t i, j, k, nnz = 0;
  memset(p, 0, sizeof(*p));
  p->n = n;
  if (!n)
    return 0;
  matrix = (unsigned char *)calloc(n * n, sizeof(*matrix));
  if (!matrix)
    return -1;
  for (i = 0; i < n; ++i)
    set_pair(matrix, n, i, i);
  for (i = 0; i < m->reaction_count; ++i) {
    const chem_reaction *r = &m->reactions[i];
    unsigned d, row;
    for (d = 0; d < r->reactant_count; ++d) {
      size_t column = m->reactants[r->reactant_offset + d];
      for (row = 0; row < r->reactant_count; ++row)
        set_pair(matrix, n, m->reactants[r->reactant_offset + row], column);
      for (row = 0; row < r->product_count; ++row)
        set_pair(matrix, n, m->products[r->product_offset + row], column);
    }
  }
  /* Symbolic fill for the no-pivot Mozart elimination ordering. */
  for (k = 0; k < n; ++k) {
    for (i = k + 1; i < n; ++i) {
      if (matrix[i * n + k]) {
        for (j = k + 1; j < n; ++j)
          if (matrix[k * n + j])
            matrix[i * n + j] = 1;
      }
    }
  }
  for (i = 0; i < n * n; ++i)
    nnz += matrix[i] != 0;
  p->nnz = nnz;
  p->row_offsets = (size_t *)calloc(n + 1, sizeof(*p->row_offsets));
  p->column_offsets = (size_t *)calloc(n + 1, sizeof(*p->column_offsets));
  p->columns = (uint32_t *)malloc(nnz * sizeof(*p->columns));
  p->rows = (uint32_t *)malloc(nnz * sizeof(*p->rows));
  p->diagonal = (size_t *)malloc(n * sizeof(*p->diagonal));
  p->position = (size_t *)malloc(n * n * sizeof(*p->position));
  p->mass_diagonal = (double *)malloc(n * sizeof(*p->mass_diagonal));
  if (!p->row_offsets || !p->column_offsets || !p->columns || !p->rows ||
      !p->diagonal || !p->position || !p->mass_diagonal) {
    free(matrix);
    chem_free_pattern(p);
    return -1;
  }
  for (i = 0; i < n * n; ++i)
    p->position[i] = SIZE_MAX;
  for (i = 0; i < n; ++i) {
    p->mass_diagonal[i] = 1.0;
    p->row_offsets[i] = p->row_offsets[i + 1] = 0;
    p->column_offsets[i] = p->column_offsets[i + 1] = 0;
  }
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      if (matrix[i * n + j]) {
        ++p->row_offsets[i + 1];
        ++p->column_offsets[j + 1];
      }
  for (i = 0; i < n; ++i) {
    p->row_offsets[i + 1] += p->row_offsets[i];
    p->column_offsets[i + 1] += p->column_offsets[i];
  }
  {
    size_t *row_next = (size_t *)malloc(n * sizeof(*row_next));
    size_t *column_next = (size_t *)malloc(n * sizeof(*column_next));
    if (!row_next || !column_next) {
      free(row_next);
      free(column_next);
      free(matrix);
      chem_free_pattern(p);
      return -1;
    }
    memcpy(row_next, p->row_offsets, n * sizeof(*row_next));
    memcpy(column_next, p->column_offsets, n * sizeof(*column_next));
    for (i = 0; i < n; ++i) {
      for (j = 0; j < n; ++j) {
        if (matrix[i * n + j]) {
          size_t pos = row_next[i]++;
          p->columns[pos] = (uint32_t)j;
          p->position[i * n + j] = pos;
          p->rows[column_next[j]++] = (uint32_t)i;
          if (i == j)
            p->diagonal[i] = pos;
        }
      }
    }
    free(row_next);
    free(column_next);
  }
  free(matrix);
  return 0;
}

void chem_free_pattern(chem_sparse_pattern *p)
{
  if (!p)
    return;
  free(p->row_offsets);
  free(p->columns);
  free(p->column_offsets);
  free(p->rows);
  free(p->diagonal);
  free(p->position);
  free(p->mass_diagonal);
  memset(p, 0, sizeof(*p));
}

static int chem_make_dense_pattern(size_t n, chem_sparse_pattern *pattern)
{
  size_t row, column, nnz;
  if (!n || n > UINT32_MAX || n > SIZE_MAX / n)
    return -1;
  nnz = n * n;
  if (nnz > SIZE_MAX / sizeof(uint32_t) ||
      n == SIZE_MAX || n + 1 > SIZE_MAX / sizeof(size_t))
    return -1;
  *pattern = (chem_sparse_pattern){0};
  pattern->n = n;
  pattern->nnz = nnz;
  pattern->row_offsets = (size_t *)malloc((n + 1) * sizeof(size_t));
  pattern->column_offsets = (size_t *)malloc((n + 1) * sizeof(size_t));
  pattern->columns = (uint32_t *)malloc(nnz * sizeof(uint32_t));
  pattern->rows = (uint32_t *)malloc(nnz * sizeof(uint32_t));
  pattern->diagonal = (size_t *)malloc(n * sizeof(size_t));
  pattern->position = (size_t *)malloc(nnz * sizeof(size_t));
  pattern->mass_diagonal = (double *)malloc(n * sizeof(*pattern->mass_diagonal));
  if (!pattern->row_offsets || !pattern->column_offsets || !pattern->columns ||
      !pattern->rows || !pattern->diagonal || !pattern->position ||
      !pattern->mass_diagonal) {
    chem_free_pattern(pattern);
    return -1;
  }
  for (row = 0; row < n; ++row) {
    pattern->mass_diagonal[row] = 1.0;
    pattern->row_offsets[row] = row * n;
    pattern->column_offsets[row] = row * n;
    pattern->diagonal[row] = row * n + row;
    for (column = 0; column < n; ++column) {
      pattern->columns[row * n + column] = (uint32_t)column;
      pattern->rows[column * n + row] = (uint32_t)row;
      pattern->position[row * n + column] = row * n + column;
    }
  }
  pattern->row_offsets[n] = nnz;
  pattern->column_offsets[n] = nnz;
  return 0;
}

static int evaluate_rhs(const chem_mechanism *m, const chem_sparse_pattern *pattern,
                        const double *y, const double *rates, double *forcing,
                        double *negative_jacobian,
                        const chem_constraint_set *constraints,
                        const double *equilibrium_constants,
                        double temperature, size_t cell)
{
  size_t i;
  if (!constraints && m->kernels && m->kernels->evaluate &&
      m->kernels->species_count == pattern->n &&
      m->kernels->pattern_nonzeros == pattern->nnz)
    return m->kernels->evaluate(y, rates, forcing, negative_jacobian, pattern);
  memset(forcing, 0, m->species_count * sizeof(*forcing));
  if (negative_jacobian)
    memset(negative_jacobian, 0, pattern->nnz * sizeof(*negative_jacobian));
  for (i = 0; i < m->reaction_count; ++i) {
    const chem_reaction *r = &m->reactions[i];
    double flux = rates[i];
    unsigned a, b;
    for (a = 0; a < r->reactant_count; ++a)
      flux *= y[m->reactants[r->reactant_offset + a]];
    for (a = 0; a < r->reactant_count; ++a)
      forcing[m->reactants[r->reactant_offset + a]] -= flux;
    for (a = 0; a < r->product_count; ++a)
      forcing[m->products[r->product_offset + a]] +=
          m->yields[r->product_offset + a] * flux;
    if (!negative_jacobian)
      continue;
    for (a = 0; a < r->reactant_count; ++a) {
      const uint32_t independent = m->reactants[r->reactant_offset + a];
      double derivative = rates[i];
      for (b = 0; b < r->reactant_count; ++b)
        if (b != a)
          derivative *= y[m->reactants[r->reactant_offset + b]];
      for (b = 0; b < r->reactant_count; ++b) {
        const uint32_t row = m->reactants[r->reactant_offset + b];
        negative_jacobian[pattern->position[(size_t)row * m->species_count + independent]] += derivative;
      }
      for (b = 0; b < r->product_count; ++b) {
        const uint32_t row = m->products[r->product_offset + b];
        negative_jacobian[pattern->position[(size_t)row * m->species_count + independent]] -=
            m->yields[r->product_offset + b] * derivative;
      }
    }
  }
  if (constraints) {
    double *gradient = NULL;
    size_t constraint_index, variable;
    if (negative_jacobian) {
      gradient = (double *)malloc(m->species_count * sizeof(*gradient));
      if (!gradient)
        return -1;
    }
    for (constraint_index = 0; constraint_index < constraints->constraint_count;
         ++constraint_index) {
      const chem_constraint *constraint = &constraints->constraints[constraint_index];
      const size_t row = constraint->algebraic_species;
      const double k = equilibrium_constants
                           ? equilibrium_constants[cell * constraints->constraint_count + constraint_index]
                           : 0.0;
      double residual;
      if (chem_constraint_residual(constraint, y, m->species_count, k,
                                   temperature, &residual) != 0) {
        free(gradient);
        return -1;
      }
      forcing[row] = residual;
      if (negative_jacobian) {
        if (chem_constraint_jacobian(constraint, y, m->species_count, k,
                                     temperature, gradient) != 0) {
          free(gradient);
          return -1;
        }
        for (variable = 0; variable < m->species_count; ++variable)
          negative_jacobian[pattern->position[row * m->species_count + variable]] =
              -gradient[variable];
      }
    }
    free(gradient);
  }
  return 0;
}

static int factor(const chem_sparse_pattern *p, double alpha,
                  const double *negative_jacobian, double *lu)
{
  size_t i, k;
  memcpy(lu, negative_jacobian, p->nnz * sizeof(*lu));
  for (i = 0; i < p->n; ++i)
    lu[p->diagonal[i]] += alpha * p->mass_diagonal[i];
  for (k = 0; k < p->n; ++k) {
    const double pivot = lu[p->diagonal[k]];
    if (pivot == 0.0 || !isfinite(pivot))
      return -1;
    for (i = p->column_offsets[k]; i < p->column_offsets[k + 1]; ++i) {
      const size_t row = p->rows[i];
      size_t q;
      double multiplier;
      if (row <= k)
        continue;
      multiplier = lu[p->position[row * p->n + k]] / pivot;
      lu[p->position[row * p->n + k]] = multiplier;
      for (q = p->row_offsets[k]; q < p->row_offsets[k + 1]; ++q) {
        const size_t column = p->columns[q];
        if (column > k)
          /* Preserve small Schur complements when stiff reversible paths make
             the product and updated matrix entry nearly cancel. */
          lu[p->position[row * p->n + column]] =
              fma(-multiplier, lu[q], lu[p->position[row * p->n + column]]);
      }
    }
  }
  return 0;
}

static void solve_lu(const chem_sparse_pattern *p, const double *lu,
                     double *right_hand_side)
{
  size_t i, j;
  /* Forward substitution through the unit-diagonal lower factor. */
  for (i = 0; i < p->n; ++i) {
    double value = right_hand_side[i];
    for (j = p->row_offsets[i]; j < p->row_offsets[i + 1]; ++j) {
      const size_t column = p->columns[j];
      if (column < i)
        value = fma(-lu[j], right_hand_side[column], value);
    }
    right_hand_side[i] = value;
  }
  /* Back substitution through the upper factor. */
  for (i = p->n; i-- > 0;) {
    double value = right_hand_side[i];
    for (j = p->row_offsets[i]; j < p->row_offsets[i + 1]; ++j) {
      const size_t column = p->columns[j];
      if (column > i)
        value = fma(-lu[j], right_hand_side[column], value);
    }
    right_hand_side[i] = value / lu[p->diagonal[i]];
  }
}

static int default_factor_callback(const chem_sparse_pattern *pattern,
                                   double alpha,
                                   const double *negative_jacobian,
                                   double *lu, void *context)
{
  (void)context;
  return factor(pattern, alpha, negative_jacobian, lu);
}

static void default_solve_callback(const chem_sparse_pattern *pattern,
                                   const double *lu, double *right_hand_side,
                                   void *context)
{
  (void)context;
  solve_lu(pattern, lu, right_hand_side);
}

static double default_step_factor(double error, double order,
                                  double safety_factor, void *context)
{
  (void)context;
  return safety_factor / pow(error, 1.0 / order);
}

void chem_rosenbrock_default(chem_rosenbrock_parameters *p)
{
  if (!p)
    return;
  memset(p, 0, sizeof(*p));
  p->stages = 3;
  p->new_function_evaluation[0] = 1;
  p->new_function_evaluation[1] = 1;
  p->new_function_evaluation[2] = 0;
  p->a[0] = 1.0;
  p->a[1] = 1.0;
  p->c[0] = -0.10156171083877702091975600115545e+01;
  p->c[1] = 0.40759956452537699824805835358067e+01;
  p->c[2] = 0.92076794298330791242156818474003e+01;
  p->m[0] = 1.0;
  p->m[1] = 0.61697947043828245592553615689730e+01;
  p->m[2] = -0.42772256543218573326238373806514;
  p->e[0] = 0.5;
  p->e[1] = -0.29079558716805469821718236208017e+01;
  p->e[2] = 0.22354069897811569627360909276199;
  p->gamma[0] = 0.43586652150845899941601945119356;
  p->gamma[1] = 0.24291996454816804366592249683314;
  p->gamma[2] = 0.21851380027664058511513169485832e+01;
  p->estimator_order = 3.0;
  p->relative_tolerance = 1.0e-6;
  p->default_absolute_tolerance = 1.0e-3;
  p->factor_min = 0.2;
  p->factor_max = 6.0;
  p->rejection_factor_decrease = 0.1;
  p->safety_factor = 0.9;
  p->h_min_fraction = 1.0e-15;
  p->h_start_fraction = 1.0e-6;
  p->round_off = DBL_EPSILON;
  p->max_number_of_steps = 1000;
  p->constraint_init_max_iterations = 10;
  p->constraint_init_tolerance = 0.1;
  p->constraint_init_max_backtracks = 24;
  p->constraint_init_backtrack_factor = 0.5;
  p->constraint_init_sufficient_decrease = 1.0e-4;
  p->controller = default_step_factor;
  p->factor = default_factor_callback;
  p->solve = default_solve_callback;
}

void chem_rosenbrock_two_stage(chem_rosenbrock_parameters *p)
{
  const double gamma = 1.0 + 1.0 / sqrt(2.0);
  if (!p)
    return;
  chem_rosenbrock_default(p);
  p->stages = 2;
  memset(p->a, 0, sizeof(p->a));
  memset(p->c, 0, sizeof(p->c));
  memset(p->m, 0, sizeof(p->m));
  memset(p->e, 0, sizeof(p->e));
  memset(p->gamma, 0, sizeof(p->gamma));
  memset(p->new_function_evaluation, 1, sizeof(p->new_function_evaluation));
  p->a[0] = 1.0 / gamma;
  p->c[0] = -2.0 / gamma;
  p->m[0] = 3.0 / (2.0 * gamma);
  p->m[1] = 1.0 / (2.0 * gamma);
  p->e[0] = 1.0 / (2.0 * gamma);
  p->e[1] = 1.0 / (2.0 * gamma);
  p->gamma[0] = gamma;
  p->gamma[1] = -gamma;
  p->estimator_order = 2.0;
}

void chem_rosenbrock_four_stage(chem_rosenbrock_parameters *p)
{
  if (!p)
    return;
  chem_rosenbrock_default(p);
  p->stages = 4;
  memset(p->a, 0, sizeof(p->a));
  memset(p->c, 0, sizeof(p->c));
  memset(p->m, 0, sizeof(p->m));
  memset(p->e, 0, sizeof(p->e));
  memset(p->gamma, 0, sizeof(p->gamma));
  memset(p->new_function_evaluation, 0, sizeof(p->new_function_evaluation));
  p->new_function_evaluation[0] = 1;
  p->new_function_evaluation[1] = 1;
  p->new_function_evaluation[2] = 1;
  p->a[0] = 2.0;
  p->a[1] = 1.867943637803922;
  p->a[2] = 0.2344449711399156;
  p->a[3] = p->a[1];
  p->a[4] = p->a[2];
  p->c[0] = -7.137615036412310;
  p->c[1] = 2.580708087951457;
  p->c[2] = 0.6515950076447975;
  p->c[3] = -2.137148994382534;
  p->c[4] = -0.3214669691237626;
  p->c[5] = -0.6949742501781779;
  p->m[0] = 2.255570073418735;
  p->m[1] = 0.2870493262186792;
  p->m[2] = 0.4353179431840180;
  p->m[3] = 1.093502252409163;
  p->e[0] = -0.2815431932141155;
  p->e[1] = -0.07276199124938920;
  p->e[2] = -0.1082196201495311;
  p->e[3] = -1.093502252409163;
  p->gamma[0] = 0.5728200000000000;
  p->gamma[1] = -1.769193891319233;
  p->gamma[2] = 0.7592633437920482;
  p->gamma[3] = -0.1049021087100450;
  p->estimator_order = 4.0;
}

void chem_rosenbrock_four_stage_dae(chem_rosenbrock_parameters *p)
{
  if (!p)
    return;
  chem_rosenbrock_default(p);
  p->stages = 4;
  memset(p->a, 0, sizeof(p->a));
  memset(p->c, 0, sizeof(p->c));
  memset(p->m, 0, sizeof(p->m));
  memset(p->e, 0, sizeof(p->e));
  memset(p->gamma, 0, sizeof(p->gamma));
  memset(p->new_function_evaluation, 0, sizeof(p->new_function_evaluation));
  p->new_function_evaluation[0] = 1;
  p->new_function_evaluation[2] = 1;
  p->new_function_evaluation[3] = 1;
  p->a[1] = 2.0;
  p->a[3] = 2.0;
  p->a[5] = 1.0;
  p->c[0] = 4.0;
  p->c[1] = 1.0;
  p->c[2] = -1.0;
  p->c[3] = 1.0;
  p->c[4] = -1.0;
  p->c[5] = -8.0 / 3.0;
  p->m[0] = 2.0;
  p->m[2] = 1.0;
  p->m[3] = 1.0;
  p->e[3] = 1.0;
  p->gamma[0] = 0.5;
  p->gamma[1] = 1.5;
  p->estimator_order = 3.0;
}

void chem_rosenbrock_six_stage_dae(chem_rosenbrock_parameters *p)
{
  if (!p)
    return;
  chem_rosenbrock_default(p);
  p->stages = 6;
  memset(p->a, 0, sizeof(p->a));
  memset(p->c, 0, sizeof(p->c));
  memset(p->m, 0, sizeof(p->m));
  memset(p->e, 0, sizeof(p->e));
  memset(p->gamma, 0, sizeof(p->gamma));
  memset(p->new_function_evaluation, 1, sizeof(p->new_function_evaluation));
  p->gamma[0] = 0.25;
  p->gamma[1] = -0.1043;
  p->gamma[2] = 0.1035;
  p->gamma[3] = -0.0362;
  p->a[0] = 1.544;
  p->a[1] = 0.9466785280815826;
  p->a[2] = 0.2557011698983284;
  p->a[3] = 3.314825187068521;
  p->a[4] = 2.896124015972201;
  p->a[5] = 0.9986419139977817;
  p->a[6] = 1.221224509226641;
  p->a[7] = 6.019134481288629;
  p->a[8] = 12.53708332932087;
  p->a[9] = -0.6878860361058950;
  p->a[10] = p->a[6];
  p->a[11] = p->a[7];
  p->a[12] = p->a[8];
  p->a[13] = p->a[9];
  p->a[14] = 1.0;
  p->c[0] = -5.6688;
  p->c[1] = -2.430093356833875;
  p->c[2] = -0.2063599157091915;
  p->c[3] = -0.1073529058151375;
  p->c[4] = -9.594562251023355;
  p->c[5] = -20.47028614809616;
  p->c[6] = 7.496443313967647;
  p->c[7] = -10.24680431464352;
  p->c[8] = -33.99990352819905;
  p->c[9] = 11.70890893206160;
  p->c[10] = 8.083246795921522;
  p->c[11] = -7.981132988064893;
  p->c[12] = -31.52159432874371;
  p->c[13] = 16.31930543123136;
  p->c[14] = -6.058818238834054;
  p->m[0] = p->a[6];
  p->m[1] = p->a[7];
  p->m[2] = p->a[8];
  p->m[3] = p->a[9];
  p->m[4] = 1.0;
  p->m[5] = 1.0;
  p->e[5] = 1.0;
  p->estimator_order = 4.0;
}

static int valid_mechanism(const chem_mechanism *m)
{
  size_t i, j;
  if (!m || !m->species || !m->reactions || m->species_count == 0 ||
      (m->reactant_count && !m->reactants) ||
      (m->product_count && (!m->products || !m->yields)))
    return 0;
  for (i = 0; i < m->reaction_count; ++i) {
    const chem_reaction *r = &m->reactions[i];
    if (r->kind < CHEM_ARRHENIUS || r->kind > CHEM_CALLBACK ||
        r->reactant_offset + r->reactant_count > m->reactant_count ||
        r->product_offset + r->product_count > m->product_count)
      return 0;
    for (j = 0; j < r->reactant_count; ++j)
      if (m->reactants[r->reactant_offset + j] >= m->species_count)
        return 0;
    for (j = 0; j < r->product_count; ++j)
      if (m->products[r->product_offset + j] >= m->species_count)
        return 0;
    if (r->kind == CHEM_USER_DEFINED &&
        (r->custom0 < 0 || (size_t)r->custom0 >= m->custom_count))
      return 0;
    if (r->kind == CHEM_SURFACE) {
      const chem_species *surface;
      if (r->surface_species < 0 || (size_t)r->surface_species >= m->species_count ||
          r->custom0 < 0 || r->custom1 < 0 ||
          (size_t)r->custom0 >= m->custom_count ||
          (size_t)r->custom1 >= m->custom_count)
        return 0;
      surface = &m->species[r->surface_species];
      if (!(surface->molecular_weight > 0.0) ||
          !(surface->diffusion_coefficient > 0.0))
        return 0;
    }
  }
  return 1;
}

static int solve_group(chem_solver *s, size_t thread,
                       size_t first, size_t last, double *y, double dt,
                       chem_solver_stats *stats)
{
  const chem_mechanism *m = s->mechanism;
  const chem_rosenbrock_parameters *p = &s->parameters;
  const size_t n = s->n;
  double *k_all[6];
  double *ynew_all = work_vector(s, YNEW_WORK);
  double *error_all = work_vector(s, ERROR_WORK);
  double *temp_all = work_vector(s, TEMP_WORK);
  double *lu = s->factor_work + thread * 2 * s->pattern.nnz;
  double *negative_jacobian = lu + s->pattern.nnz;
  chem_lu_factor_callback factor_callback = p->factor ? p->factor : default_factor_callback;
  chem_lu_solve_callback solve_callback = p->solve ? p->solve : default_solve_callback;
  const double h_min = p->h_min_fraction * dt;
  const double h_max = dt;
  double h = fmin(h_max, fmax(h_min, p->h_start_fraction * dt));
  double present = 0.0;
  int reject_last = 0, reject_more = 0;
  size_t stage;

  if (!p->stages || p->stages > 6 || !(p->gamma[0] != 0.0) ||
      !(p->estimator_order > 0.0))
    return -1;
  if (m->kernels && m->kernels->factor &&
      m->kernels->species_count == s->pattern.n &&
      m->kernels->pattern_nonzeros == s->pattern.nnz &&
      (!p->factor || p->factor == default_factor_callback))
    factor_callback = m->kernels->factor;
  if (m->kernels && m->kernels->solve &&
      m->kernels->species_count == s->pattern.n &&
      m->kernels->pattern_nonzeros == s->pattern.nnz &&
      (!p->solve || p->solve == default_solve_callback))
    solve_callback = m->kernels->solve;
  for (stage = 0; stage < p->stages; ++stage)
    k_all[stage] = work_vector(s, K_WORK_BASE + stage);

  while (present - dt + p->round_off <= 0.0) {
    int accepted = 0;
    if (stats->internal_steps > p->max_number_of_steps)
      return -1;
    if ((present + 0.1 * h == present) || h <= fmin(p->round_off, h_min))
      return -1;
    h = fmin(h, fabs(dt - present));
    while (!accepted) {
      const double alpha = 1.0 / (h * p->gamma[0]);
      double error;
      size_t cell, v;
      uint64_t function_evaluations = 1;
      for (stage = 1; stage < p->stages; ++stage)
        function_evaluations += p->new_function_evaluation[stage] != 0;
      for (cell = first; cell < last; ++cell) {
        const size_t offset = cell * n;
        double *ynew = ynew_all + offset;
        double *err = error_all + offset;
        double *temporary = temp_all + offset;
        const double *rates = s->rates + cell * s->nr;
        const size_t n_constraints = s->constraints
                                         ? s->constraints->constraint_count : 0;
        const double *equilibrium_constants = s->equilibrium_constants
            ? s->equilibrium_constants + cell * n_constraints : NULL;
        const double temperature = s->constraint_temperatures
            ? s->constraint_temperatures[cell] : CHEM_DEFAULT_REFERENCE_TEMPERATURE;
        if (evaluate_rhs(m, &s->pattern, y + offset, rates,
                         k_all[0] + offset, negative_jacobian,
                         s->constraints, equilibrium_constants, temperature,
                         cell) != 0)
          return -1;
        if (p->stages > 1 && !p->new_function_evaluation[1])
          memcpy(k_all[1] + offset, k_all[0] + offset, n * sizeof(double));
        if (factor_callback(&s->pattern, alpha, negative_jacobian, lu,
                            p->callback_context) != 0)
          return -1;
        solve_callback(&s->pattern, lu, k_all[0] + offset, p->callback_context);
        for (stage = 1; stage < p->stages; ++stage) {
          double *k = k_all[stage] + offset;
          size_t j;
          const size_t combination = stage * (stage - 1) / 2;
          if (p->new_function_evaluation[stage]) {
            for (v = 0; v < n; ++v)
              temporary[v] = y[offset + v];
            for (j = 0; j < stage; ++j) {
              const double a = p->a[combination + j];
              if (a != 0.0)
                for (v = 0; v < n; ++v)
                  temporary[v] += a * k_all[j][offset + v];
            }
            if (evaluate_rhs(m, &s->pattern, temporary, rates, k, NULL,
                             s->constraints, equilibrium_constants, temperature,
                             cell) != 0)
              return -1;
            if (stage + 1 < p->stages && !p->new_function_evaluation[stage + 1])
              memcpy(k_all[stage + 1] + offset, k, n * sizeof(double));
          }
          for (j = 0; j < stage; ++j) {
            const double coefficient = p->c[combination + j] / h;
            if (coefficient != 0.0) {
              for (v = 0; v < n; ++v)
                if (!s->constraints || s->pattern.mass_diagonal[v] != 0.0)
                  k[v] += coefficient * k_all[j][offset + v];
            }
          }
          solve_callback(&s->pattern, lu, k, p->callback_context);
        }
        for (v = 0; v < n; ++v) {
          double y1 = y[offset + v];
          double error = 0.0;
          for (stage = 0; stage < p->stages; ++stage) {
            y1 += p->m[stage] * k_all[stage][offset + v];
            error += p->e[stage] * k_all[stage][offset + v];
          }
          if (s->constraints && s->pattern.mass_diagonal[v] == 0.0)
            error = y1 - y[offset + v];
          err[v] = error;
          ynew[v] = y1;
        }
      }
      error = chem_normalized_error(
          y + first * n, ynew_all + first * n, error_all + first * n,
          m->species, last - first, n, p->relative_tolerance,
          p->default_absolute_tolerance);
      stats->function_calls += function_evaluations;
      ++stats->jacobian_updates;
      ++stats->decompositions;
      stats->solves += p->stages;
      {
        double factor = p->controller
                            ? p->controller(error, p->estimator_order, p->safety_factor,
                                            p->callback_context)
                            : default_step_factor(error, p->estimator_order,
                                                  p->safety_factor, NULL);
        double hnew;
        if (!isfinite(error) || !isfinite(factor))
          return -1;
        factor = fmin(p->factor_max, fmax(p->factor_min, factor));
        hnew = h * factor;
        ++stats->internal_steps;
        if (error < 1.0 || h < h_min) {
          ++stats->accepted;
          present += h;
          for (cell = first; cell < last; ++cell)
            memcpy(y + cell * n, ynew_all + cell * n, n * sizeof(*y));
          hnew = fmax(h_min, fmin(hnew, h_max));
          if (reject_last)
            hnew = fmin(hnew, h);
          reject_last = 0;
          reject_more = 0;
          h = hnew;
          accepted = 1;
        } else {
          if (reject_more)
            hnew = h * p->rejection_factor_decrease;
          reject_more = reject_last;
          reject_last = 1;
          h = hnew;
          if (stats->accepted >= 1)
            ++stats->rejected;
        }
      }
    }
  }
  return 0;
}

chem_solver *chem_solver_create(const chem_mechanism *m, size_t cells,
                                const double *temperature, const double *pressure,
                                const double *custom_parameters,
                                const chem_rosenbrock_parameters *parameters,
                                size_t cell_group_size)
{
  chem_solver *s = (chem_solver *)calloc(1, sizeof(*s));
  size_t cell, reaction;
  if (!s)
    return NULL;
  if (!valid_mechanism(m)) {
    free(s);
    return NULL;
  }
  s->mechanism = m;
  s->cells = cells;
  s->n = m->species_count;
  s->nr = m->reaction_count;
  s->ncustom = m->custom_count;
  if (parameters)
    s->parameters = *parameters;
  else
    chem_rosenbrock_default(&s->parameters);
#ifdef _OPENMP
  s->thread_count = omp_get_max_threads();
#else
  s->thread_count = 1;
#endif
  if (s->thread_count < 1)
    s->thread_count = 1;
  s->group_size = cell_group_size;
  s->groups = cell_group_size
                  ? (cells ? (cells + cell_group_size - 1) / cell_group_size : 1)
                  : (size_t)s->thread_count;
  if (chem_make_pattern(m, &s->pattern) != 0) {
    chem_solver_destroy(s);
    return NULL;
  }
  s->rates = (double *)malloc((cells != 0 && s->nr != 0 ? cells * s->nr : 1) * sizeof(*s->rates));
  s->work = (double *)malloc((cells != 0 && s->n != 0 ? WORK_COUNT * cells * s->n : 1) * sizeof(*s->work));
  s->factor_work = (double *)malloc((size_t)s->thread_count * 2 * s->pattern.nnz * sizeof(*s->factor_work));
  if (!s->rates || !s->work || !s->factor_work) {
    chem_solver_destroy(s);
    return NULL;
  }
  for (cell = 0; cell < cells; ++cell) {
    const double t = temperature ? temperature[cell] : m->default_temperature;
    const double p = pressure ? pressure[cell] : m->default_pressure;
    const double *custom = custom_parameters ? custom_parameters + cell * m->custom_count : NULL;
    if (m->kernels && m->kernels->rates) {
      m->kernels->rates(t, p, custom, s->rates + cell * s->nr);
    } else {
      for (reaction = 0; reaction < s->nr; ++reaction) {
        const chem_reaction *r = &m->reactions[reaction];
        s->rates[cell * s->nr + reaction] =
            chem_effective_rate_constant(m, r, t, p, custom);
      }
    }
  }
  return s;
}

void chem_solver_destroy(chem_solver *s)
{
  if (!s)
    return;
  chem_free_pattern(&s->pattern);
  free(s->rates);
  free(s->work);
  free(s->factor_work);
  free(s);
}

int chem_solver_set_constraints(chem_solver *s,
                                const chem_constraint_set *constraints,
                                const double *equilibrium_constants,
                                const double *constraint_temperatures)
{
  chem_sparse_pattern dense = {0};
  chem_constraint_set validated;
  double *factor_work = NULL;
  size_t constraint_index, prior;
  if (!s)
    return -1;
  if (constraints &&
      (constraints->species_count != s->n ||
       chem_constraint_set_init(&validated, constraints->constraints,
                                constraints->constraint_count,
                                constraints->species_count) != 0))
    return -1;
  if (constraints && constraints->constraint_count) {
    size_t work_values;
    if (chem_make_dense_pattern(s->n, &dense) != 0)
      return -1;
    for (constraint_index = 0; constraint_index < constraints->constraint_count;
         ++constraint_index) {
      const size_t algebraic = constraints->constraints[constraint_index].algebraic_species;
      for (prior = 0; prior < constraint_index; ++prior)
        if (constraints->constraints[prior].algebraic_species == algebraic) {
          chem_free_pattern(&dense);
          return -1;
        }
      dense.mass_diagonal[algebraic] = 0.0;
    }
    if (dense.nnz > SIZE_MAX / (2 * sizeof(*factor_work))) {
      chem_free_pattern(&dense);
      return -1;
    }
    work_values = 2 * dense.nnz;
    if (work_values > SIZE_MAX / sizeof(*factor_work) ||
        (size_t)s->thread_count > SIZE_MAX / (work_values * sizeof(*factor_work))) {
      chem_free_pattern(&dense);
      return -1;
    }
    factor_work = (double *)malloc((size_t)s->thread_count * work_values * sizeof(*factor_work));
    if (!factor_work) {
      chem_free_pattern(&dense);
      return -1;
    }
    chem_free_pattern(&s->pattern);
    s->pattern = dense;
    free(s->factor_work);
    s->factor_work = factor_work;
  }
  s->constraints = constraints;
  s->equilibrium_constants = equilibrium_constants;
  s->constraint_temperatures = constraint_temperatures;
  return 0;
}

double chem_solver_absolute_tolerance(const chem_solver *s, size_t species)
{
  double tolerance;
  if (!s || species >= s->n)
    return NAN;
  tolerance = s->mechanism->species[species].absolute_tolerance;
  return tolerance > 0.0 ? tolerance : s->parameters.default_absolute_tolerance;
}

double chem_solver_relative_tolerance(const chem_solver *s)
{
  return s ? s->parameters.relative_tolerance : NAN;
}

int chem_solver_step(chem_solver *s, double *concentrations, double dt,
                     chem_solver_stats *stats)
{
  int group;
  chem_solver_stats *group_stats;
  int *group_status;
  if (!s || !concentrations || !stats || dt <= 0.0)
    return -1;
  memset(stats, 0, sizeof(*stats));
  if (s->constraints && s->constraints->constraint_count) {
    size_t projection_iterations = 0;
    double *absolute_tolerances = (double *)malloc(s->n * sizeof(*absolute_tolerances));
    size_t species;
    int projection_status;
    if (!absolute_tolerances)
      return -1;
    for (species = 0; species < s->n; ++species)
      absolute_tolerances[species] = chem_solver_absolute_tolerance(s, species);
    projection_status = chem_constraint_set_project(
        s->constraints, s->cells, concentrations, s->equilibrium_constants,
        s->constraint_temperatures, absolute_tolerances,
        s->parameters.relative_tolerance,
        s->parameters.constraint_init_tolerance,
        s->parameters.constraint_init_max_iterations,
        s->parameters.constraint_init_max_backtracks,
        s->parameters.constraint_init_backtrack_factor,
        s->parameters.constraint_init_sufficient_decrease,
        &projection_iterations);
    free(absolute_tolerances);
    stats->constraint_init_iterations = projection_iterations;
    if (projection_status != 0)
      return (int)CHEM_SOLVER_CONSTRAINT_INITIALIZATION_FAILED;
  }
  group_stats = (chem_solver_stats *)calloc(s->groups, sizeof(*group_stats));
  group_status = (int *)calloc(s->groups, sizeof(*group_status));
  if (!group_stats || !group_status) {
    free(group_stats);
    free(group_status);
    return -1;
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(s->thread_count)
#endif
  for (group = 0; group < (int)s->groups; ++group) {
    const size_t first = s->group_size
                             ? (size_t)group * s->group_size
                             : s->cells * (size_t)group / s->groups;
    const size_t last = s->group_size
                            ? (first + s->group_size < s->cells ? first + s->group_size : s->cells)
                            : s->cells * (size_t)(group + 1) / s->groups;
    size_t thread = 0;
#ifdef _OPENMP
    thread = (size_t)omp_get_thread_num();
#endif
    if (first < last)
      group_status[group] = solve_group(s, thread, first, last,
                                        concentrations, dt, &group_stats[group]);
  }
  for (group = 0; group < (int)s->groups; ++group) {
    if (group_status[group] != 0) {
      free(group_stats);
      free(group_status);
      return -1;
    }
    stats->internal_steps += group_stats[group].internal_steps;
    stats->accepted += group_stats[group].accepted;
    stats->rejected += group_stats[group].rejected;
    stats->function_calls += group_stats[group].function_calls;
    stats->jacobian_updates += group_stats[group].jacobian_updates;
    stats->decompositions += group_stats[group].decompositions;
    stats->solves += group_stats[group].solves;
  }
  /* MICM clamps ODE state variables after each top-level solve. Keep
     algebraic variables untouched, as their values are set by constraints. */
  {
    size_t cell, species;
    for (cell = 0; cell < s->cells; ++cell) {
      for (species = 0; species < s->n; ++species) {
        double *value;
        if (s->pattern.mass_diagonal[species] == 0.0)
          continue;
        value = &concentrations[cell * s->n + species];
        if (!(*value > 0.0))
          *value = 0.0;
      }
    }
  }
  free(group_stats);
  free(group_status);
  return 0;
}

int chem_solver_update_custom_rates(chem_solver *s,
                                   const double *temperature,
                                   const double *pressure,
                                   const double *custom_parameters)
{
  size_t cell, reaction;
  if (!s || (s->ncustom && s->cells && !custom_parameters))
    return -1;
  for (cell = 0; cell < s->cells; ++cell) {
    const double t = temperature ? temperature[cell] : s->mechanism->default_temperature;
    const double p = pressure ? pressure[cell] : s->mechanism->default_pressure;
    const double *custom = custom_parameters
                               ? custom_parameters + cell * s->ncustom : NULL;
    double *rates = s->rates + cell * s->nr;
    if (s->mechanism->kernels && s->mechanism->kernels->custom_rates) {
      s->mechanism->kernels->custom_rates(t, p, custom, rates);
    } else {
      for (reaction = 0; reaction < s->nr; ++reaction) {
        const chem_reaction *r = &s->mechanism->reactions[reaction];
        if (r->kind == CHEM_USER_DEFINED)
          rates[reaction] = chem_effective_rate_constant(s->mechanism, r, t, p, custom);
      }
    }
  }
  return 0;
}
