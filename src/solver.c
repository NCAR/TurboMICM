#include "chemistry.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { SOLVER_OK = 0, SOLVER_BAD_ARGUMENT = 1, SOLVER_SINGULAR = 2,
       SOLVER_STEP_FAILURE = 3, SOLVER_NONFINITE = 4 };

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
} solver_workspace_data;

static int make_pattern(const chem_mechanism *m, unsigned char *pattern,
                        size_t n, size_t *row_counts, size_t *col_counts,
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
  w->stage = calloc(3*n, sizeof(double));
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
  make_pattern(mechanism, w->pattern, n, w->row_counts, w->col_counts,
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
  chem_solver_options result = {1.0e-6, 1.0e-3, 1000};
  return result;
}

static int make_pattern(const chem_mechanism *m, unsigned char *pattern,
                        size_t n, size_t *row_counts, size_t *col_counts,
                        size_t *row_cols, size_t *col_rows) {
  memset(pattern, 0, n * n);
  for (size_t i = 0; i < n; ++i) pattern[i*n + i] = 1;
  for (size_t r = 0; r < m->reaction_count; ++r) {
    const chem_reaction *rx = &m->reactions[r];
    for (size_t i = 0; i < rx->reactant_count; ++i) {
      size_t col = rx->reactants[i].species;
      if (col >= n || m->species[col].third_body) continue;
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

int chem_rosenbrock3_solve_with_workspace(
    const chem_mechanism *m, chem_solver_workspace *workspace, double *y,
    double temperature, double pressure, double air_density,
    const double *rate_parameters,
    double duration, const chem_solver_options *requested,
    chem_solver_stats *statistics) {
  if (!m || !y || !isfinite(duration) || duration <= 0.0 ||
      !isfinite(temperature) || temperature <= 0.0 ||
      !isfinite(pressure) || pressure < 0.0 ||
      !isfinite(air_density) || air_density < 0.0 || m->species_count == 0)
    return SOLVER_BAD_ARGUMENT;
  solver_workspace_data *w = workspace ? workspace->implementation : NULL;
  if (!w || w->mechanism != m || w->n != m->species_count)
    return SOLVER_BAD_ARGUMENT;
  const size_t n = m->species_count;
  const chem_solver_options defaults = chem_solver_default_options();
  const chem_solver_options *o = requested ? requested : &defaults;
  if (o->relative_tolerance <= 0.0 || o->absolute_tolerance <= 0.0 ||
      o->max_steps == 0) return SOLVER_BAD_ARGUMENT;
  chem_solver_stats local_stats = {0};
  if (statistics) *statistics = local_stats;
  const size_t nn = n*n;
  double *jac = w->jac, *a = w->a, *forcing = w->forcing;
  double *stage = w->stage, *ynew = w->ynew, *yerror = w->yerror;
  double *stage_y = w->stage_y, *reused_forcing = w->reused_forcing;
  size_t *row_counts = w->row_counts, *col_counts = w->col_counts;
  size_t *row_cols = w->row_cols, *col_rows = w->col_rows;
  int status = SOLVER_OK;

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

  const double eps = 2.2204460492503131e-16;
  const double hmin = 1.0e-15 * duration;
  const double hmax = duration;
  double h = 1.0e-6 * duration;
  double time = 0.0;
  int rejected_previous = 0, rejected_twice = 0;
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
    if (local_stats.internal_steps > o->max_steps) { status = SOLVER_STEP_FAILURE; break; }
    if (time + 0.1*h == time || h <= fmin(eps, hmin)) { status = SOLVER_STEP_FAILURE; break; }
    if (h > duration-time) h = duration-time;
    memset(forcing, 0, n*sizeof(*forcing));
    status = chem_forcing_with_constants(m, y, w->rate_constants,
                                         air_density, forcing);
    if (status) { status = SOLVER_BAD_ARGUMENT; break; }
    ++local_stats.function_calls;
    status = chem_negative_jacobian_with_constants(m, y, w->rate_constants,
                                                    jac);
    if (status) { status = SOLVER_BAD_ARGUMENT; break; }
    ++local_stats.jacobian_updates;
    int accepted = 0;
    while (!accepted) {
      const double used_h = h;
      const double alpha = 1.0 / (h*gamma);
      for (size_t i = 0; i < nn; ++i) a[i] = jac[i];
      for (size_t i = 0; i < n; ++i) a[i*n + i] += alpha;
      status = factor(a, n, row_counts, col_counts, row_cols, col_rows);
      ++local_stats.decompositions;
      if (status) break;
      for (size_t i = 0; i < n; ++i) stage[i] = forcing[i];
      for (size_t s = 0; s < 3; ++s) {
        if (s == 1) {
          for (size_t i = 0; i < n; ++i) stage_y[i] = y[i] + aa[0]*stage[i];
          memset(stage+n, 0, n*sizeof(*stage));
          status = chem_forcing_with_constants(m, stage_y, w->rate_constants,
                                               air_density, stage+n);
          ++local_stats.function_calls;
          if (status) { status = SOLVER_BAD_ARGUMENT; break; }
          memcpy(reused_forcing, stage+n, n*sizeof(*reused_forcing));
        } else if (s == 2) {
          memcpy(stage + 2*n, reused_forcing, n*sizeof(*stage));
        }
        double *ks = stage + s*n;
        for (size_t j = 0; j < s; ++j) {
          const size_t ci = s == 1 ? 0 : (j == 0 ? 1 : 2);
          const double scale = cc[ci] / h;
          for (size_t i = 0; i < n; ++i) ks[i] += scale*stage[j*n + i];
        }
        solve_factored(a, n, row_counts, row_cols, ks);
        ++local_stats.solves;
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
      double sum = 0.0;
      size_t error_variables = 0;
      for (size_t i = 0; i < n; ++i) {
        if (m->species[i].third_body) continue;
        ++error_variables;
        const double scale = o->absolute_tolerance + o->relative_tolerance*fmax(fabs(y[i]), fabs(ynew[i]));
        const double v = yerror[i]/scale;
        sum += v*v;
      }
      double error = fmax(sqrt(sum/(double)error_variables), 1.0e-10);
      double fac = fmin(6.0, fmax(0.2, 0.9/pow(error, 1.0/3.0)));
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
        if (local_stats.accepted >= 1) ++local_stats.rejected;
        ++local_stats.jacobian_updates;
      }
    }
    if (status) break;
  }
  if (statistics) *statistics = local_stats;
  return status;
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
  if (o->relative_tolerance <= 0.0 || o->absolute_tolerance <= 0.0 ||
      o->max_steps == 0) return SOLVER_BAD_ARGUMENT;
  for (size_t c = 0; c < cell_count; ++c)
    if (!isfinite(temperatures[c]) || temperatures[c] <= 0.0 ||
        !isfinite(pressures[c]) || pressures[c] < 0.0 ||
        !isfinite(air_densities[c]) || air_densities[c] < 0.0)
      return SOLVER_BAD_ARGUMENT;

  chem_solver_stats local = {0};
  const size_t nn = n*n;
  const double eps = 2.2204460492503131e-16;
  const double hmin = 1.0e-15 * duration, hmax = duration;
  double h = 1.0e-6 * duration, time = 0.0;
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
        status = chem_negative_jacobian_with_constants(m, y, rates, w->jac);
        if (status) { status = SOLVER_BAD_ARGUMENT; break; }
        double alpha = 1.0 / (h*gamma);
        for (size_t i = 0; i < nn; ++i) w->a[i] = w->jac[i];
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
          const double scale = o->absolute_tolerance + o->relative_tolerance *
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
        if (local.accepted >= 1) ++local.rejected;
        ++local.jacobian_updates;
      }
    }
    if (status) break;
  }
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
