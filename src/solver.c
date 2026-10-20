#include "chemistry.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { SOLVER_OK = 0, SOLVER_BAD_ARGUMENT = 1, SOLVER_SINGULAR = 2,
       SOLVER_STEP_FAILURE = 3, SOLVER_NONFINITE = 4 };

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
    for (size_t q = 0; q < row_counts[k]; ++q) {
      size_t j = row_cols[k*n + q];
      if (j <= k) continue;
      const double aik = a[k*n + j];
      for (size_t p = 0; p < col_counts[k]; ++p) {
        size_t i = col_rows[k*n + p];
        if (i > k) a[i*n + j] -= a[i*n + k] * aik;
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

int chem_rosenbrock3_solve(const chem_mechanism *m, double *y,
                           double temperature, double air_density,
                           const double *rate_parameters, double duration,
                           const chem_solver_options *requested,
                           chem_solver_stats *statistics) {
  if (!m || !y || !isfinite(duration) || duration <= 0.0 ||
      !isfinite(temperature) || temperature <= 0.0 ||
      !isfinite(air_density) || air_density < 0.0 || m->species_count == 0)
    return SOLVER_BAD_ARGUMENT;
  const size_t n = m->species_count;
  const chem_solver_options defaults = chem_solver_default_options();
  const chem_solver_options *o = requested ? requested : &defaults;
  if (o->relative_tolerance <= 0.0 || o->absolute_tolerance <= 0.0 ||
      o->max_steps == 0) return SOLVER_BAD_ARGUMENT;
  chem_solver_stats local_stats = {0};
  if (statistics) *statistics = local_stats;
  if (n > (size_t)-1 / n || n*n > (size_t)-1 / sizeof(double)) return SOLVER_BAD_ARGUMENT;
  const size_t nn = n*n;
  double *jac = calloc(nn, sizeof(*jac));
  double *a = calloc(nn, sizeof(*a));
  unsigned char *pattern_tmp = calloc(nn, sizeof(*pattern_tmp));
  size_t *row_counts = calloc(n, sizeof(*row_counts));
  size_t *col_counts = calloc(n, sizeof(*col_counts));
  size_t *row_cols = calloc(nn, sizeof(*row_cols));
  size_t *col_rows = calloc(nn, sizeof(*col_rows));
  double *forcing = calloc(n, sizeof(*forcing));
  double *stage = calloc(3*n, sizeof(*stage));
  double *ynew = calloc(n, sizeof(*ynew));
  double *yerror = calloc(n, sizeof(*yerror));
  double *stage_y = calloc(n, sizeof(*stage_y));
  double *reused_forcing = calloc(n, sizeof(*reused_forcing));
  if (!jac || !a || !pattern_tmp || !row_counts || !col_counts || !row_cols ||
      !col_rows || !forcing || !stage || !ynew || !yerror || !stage_y ||
      !reused_forcing) {
    free(jac); free(a); free(pattern_tmp); free(row_counts); free(col_counts);
    free(row_cols); free(col_rows); free(forcing); free(stage); free(ynew);
    free(yerror); free(stage_y); free(reused_forcing);
    return SOLVER_BAD_ARGUMENT;
  }
  make_pattern(m, pattern_tmp, n, row_counts, col_counts, row_cols, col_rows);

  const double eps = 2.2204460492503131e-16;
  const double hmin = 1.0e-15 * duration;
  const double hmax = duration;
  double h = 1.0e-6 * duration;
  double time = 0.0;
  int rejected_previous = 0, rejected_twice = 0;
  int status = SOLVER_OK;
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
    status = chem_forcing(m, y, temperature, air_density, rate_parameters, forcing);
    if (status) { status = SOLVER_BAD_ARGUMENT; break; }
    ++local_stats.function_calls;
    status = chem_jacobian(m, y, temperature, air_density, rate_parameters, jac);
    if (status) { status = SOLVER_BAD_ARGUMENT; break; }
    ++local_stats.jacobian_updates;
    int accepted = 0;
    while (!accepted) {
      const double used_h = h;
      const double alpha = 1.0 / (h*gamma);
      for (size_t i = 0; i < nn; ++i) a[i] = -jac[i];
      for (size_t i = 0; i < n; ++i) a[i*n + i] += alpha;
      status = factor(a, n, row_counts, col_counts, row_cols, col_rows);
      ++local_stats.decompositions;
      if (status) break;
      for (size_t i = 0; i < n; ++i) stage[i] = forcing[i];
      for (size_t s = 0; s < 3; ++s) {
        if (s == 1) {
          for (size_t i = 0; i < n; ++i) stage_y[i] = y[i] + aa[0]*stage[i];
          memset(stage+n, 0, n*sizeof(*stage));
          status = chem_forcing(m, stage_y, temperature, air_density,
                                rate_parameters, stage+n);
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
  free(jac); free(a); free(pattern_tmp); free(row_counts); free(col_counts);
  free(row_cols); free(col_rows); free(forcing); free(stage); free(ynew);
  free(yerror); free(stage_y); free(reused_forcing);
  return status;
}
