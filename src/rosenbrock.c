/* Translated numerical method: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

ChemSolver chem_solver(int method) {
  ChemSolver s = {0};
  s.parameters = chem_rosenbrock_parameters(method);
  s.controller = chem_default_controller;
  s.factor = chem_lu_mozart; s.solve = chem_linear_solve;
  return s;
}
double chem_default_controller(double error, const ChemRosenbrock *p, void *context) {
  (void)context;
  return fmin(p->factor_max, fmax(p->factor_min, p->safety/pow(error, 1.0/p->order)));
}
double chem_normalized_error(const ChemState *s, const double *y, const double *yn, const double *e) {
  double total = 0.0;
  for (size_t i = 0; i < s->model->nvariables; ++i) for (size_t c = 0; c < s->cells; ++c) {
    size_t at = i*s->stride+c;
    double ymax = fabs(y[at]) > fabs(yn[at]) ? fabs(y[at]) : fabs(yn[at]);
    double scaled = e[at]/(s->absolute_tolerance[i]+s->relative_tolerance*ymax);
    total += scaled*scaled;
  }
  /* Match std::max(value, floor), including propagation of a NaN value. */
  double error = sqrt(total/(s->cells*s->model->nvariables));
  return error < 1e-10 ? 1e-10 : error;
}
static void factor_shift(ChemSolver *solver, ChemState *s, double alpha, ChemStats *stats) {
  for (size_t i = 0; i < s->model->nvariables; ++i) if (!s->model->algebraic[i]) {
    double *d = s->jacobian+s->model->diagonal[i]*s->stride;
    #pragma omp simd
    for (size_t c = 0; c < s->cells; ++c) d[c] += alpha;
  }
  solver->factor(s->model, s->cells, s->stride, s->jacobian); ++stats->decompositions;
}
static ChemStatus finite_algebraic(const ChemState *s, const double *v, double *maximum) {
  *maximum = 0.0;
  int inf = 0;
  for (size_t i = 0; i < s->model->nvariables; ++i) if (s->model->algebraic[i])
    for (size_t c = 0; c < s->cells; ++c) {
      double val = v[i*s->stride+c];
      if (isnan(val)) return CHEM_NAN;
      if (isinf(val)) inf = 1;
      if (fabs(val) > *maximum) *maximum = fabs(val);
    }
  return inf ? CHEM_INF : CHEM_OK;
}
static double weighted_correction(const ChemState *s, const double *y, const double *d) {
  double maximum = 0.0;
  for (size_t i = 0; i < s->model->nvariables; ++i) if (s->model->algebraic[i])
    for (size_t c = 0; c < s->cells; ++c) {
      size_t at = i*s->stride+c;
      double corrected = y[at]+d[at];
      double zmax = fmax(fabs(y[at]), fabs(corrected));
      double scale = s->absolute_tolerance[i]+s->relative_tolerance*zmax;
      double q = scale > 0.0 ? fabs(d[at])/scale : (d[at] == 0.0 ? 0.0 : INFINITY);
      if (q > maximum) maximum = q;
    }
  return maximum;
}
ChemStatus chem_initialize_constraints(ChemSolver *solver, ChemState *s, ChemStats *stats) {
  const ChemModel *m = s->model; const ChemRosenbrock *p = &solver->parameters;
  if (!m->nconstraints) return CHEM_OK;
  size_t ny = m->nvariables*s->stride;
  double *delta = s->forcing, *original = s->stages[0], *candidate = s->new_y, *trial = s->error;
  int modified = 0;
  ChemStatus status = CHEM_INIT_FAILED;
  for (size_t update = 0; update <= p->init_iterations; ++update) {
    double maximum;
    memset(delta, 0, ny*sizeof(double)); chem_constraints(s, s->variables, delta, NULL);
    status = finite_algebraic(s, delta, &maximum); if (status) break;
    ++stats->constraint_init_iterations;
    if (maximum == 0.0) return CHEM_OK;
    memset(s->jacobian, 0, m->nnz*s->stride*sizeof(double));
    chem_constraints(s, s->variables, NULL, s->jacobian); ++stats->jacobian_updates;
    factor_shift(solver, s, 1.0, stats);
    solver->solve(m, s->cells, s->stride, s->jacobian, delta); ++stats->solves;
    status = finite_algebraic(s, delta, &maximum); if (status) break;
    double current = weighted_correction(s, s->variables, delta);
    if (current <= p->init_tolerance) return CHEM_OK;
    status = CHEM_INIT_FAILED;
    if (update == p->init_iterations) break;
    if (!modified) { memcpy(original, s->variables, ny*sizeof(double)); modified = 1; }
    if (p->init_backtracks) {
      double step = 1.0;
      int accepted = 0;
      for (size_t backtrack = 0; backtrack <= p->init_backtracks; ++backtrack, step *= p->init_backtrack_factor) {
        memcpy(candidate, s->variables, ny*sizeof(double));
        for (size_t i = 0; i < m->nvariables; ++i) if (m->algebraic[i])
          for (size_t c = 0; c < s->cells; ++c) candidate[i*s->stride+c] += step*delta[i*s->stride+c];
        if (finite_algebraic(s, candidate, &maximum)) continue;
        memset(trial, 0, ny*sizeof(double)); chem_constraints(s, candidate, trial, NULL);
        if (finite_algebraic(s, trial, &maximum)) continue;
        if (maximum == 0.0) { accepted = 1; break; }
        solver->solve(m, s->cells, s->stride, s->jacobian, trial); ++stats->solves;
        if (finite_algebraic(s, trial, &maximum)) continue;
        double q = weighted_correction(s, candidate, trial);
        if (q <= p->init_tolerance || q < (1.0-p->init_sufficient_decrease*step)*current) { accepted = 1; break; }
      }
      if (!accepted) break;
      memcpy(s->variables, candidate, ny*sizeof(double));
    } else {
      for (size_t i = 0; i < m->nvariables; ++i) if (m->algebraic[i])
        for (size_t c = 0; c < s->cells; ++c) s->variables[i*s->stride+c] += delta[i*s->stride+c];
    }
  }
  if (modified) memcpy(s->variables, original, ny*sizeof(double));
  return status;
}

ChemStatus chem_solve(ChemSolver *solver, double time_step, ChemState *s, ChemStats *stats) {
  if (!solver || !s || !stats || !isfinite(time_step) || time_step < 0.0 || solver->parameters.stages < 1 ||
      solver->parameters.stages > 6 || !solver->controller || !solver->factor || !solver->solve) return CHEM_INVALID;
  memset(stats, 0, sizeof *stats);
  const ChemModel *m = s->model; const ChemRosenbrock *p = &solver->parameters;
  size_t ny = m->nvariables*s->stride;
  double hmin = p->h_min == 0.0 ? 1e-15*time_step : p->h_min;
  double hmax = p->h_max == 0.0 ? time_step : fmin(time_step, p->h_max);
  double hstart = p->h_start == 0.0 ? 1e-6*time_step : fmin(hmax, p->h_start);
  double h = fmin(fmax(hmin, fabs(hstart)), fabs(hmax)), time = 0.0;
  ChemStatus status = chem_initialize_constraints(solver, s, stats);
  if (status) return status;
  int reject_last = 0, reject_more = 0;
  while (time-time_step+p->roundoff <= 0.0) {
    if (stats->internal_steps > p->max_steps) { status = CHEM_MAX_STEPS; break; }
    if (time+0.1*h == time || h <= fmin(p->roundoff, hmin)) { status = CHEM_SMALL_STEP; break; }
    h = fmin(h, fabs(time_step-time));
    chem_forcing(s, s->variables, s->forcing); ++stats->function_calls;
    chem_jacobian(s, s->variables, s->jacobian); ++stats->jacobian_updates;
    int accepted = 0;
    while (!accepted) {
      factor_shift(solver, s, 1.0/(h*p->gamma[0]), stats);
      for (int stage = 0; stage < p->stages; ++stage) {
        size_t offset = (size_t)(stage*(stage-1)/2);
        if (!stage) memcpy(s->stages[stage], s->forcing, ny*sizeof(double));
        else if (p->new_function[stage]) {
          memcpy(s->new_y, s->variables, ny*sizeof(double));
          for (int j = 0; j < stage; ++j) chem_axpy(ny, s->new_y, p->a[offset+j], s->stages[j]);
          chem_forcing(s, s->new_y, s->stages[stage]); ++stats->function_calls;
        }
        if (stage+1 < p->stages && !p->new_function[stage+1])
          memcpy(s->stages[stage+1], s->stages[stage], ny*sizeof(double));
        for (int j = 0; j < stage; ++j) {
          double a = p->c[offset+j]/h;
          if (!m->nconstraints) chem_axpy(ny, s->stages[stage], a, s->stages[j]);
          else for (size_t i = 0; i < m->nvariables; ++i) if (!m->algebraic[i])
            chem_axpy(s->cells, s->stages[stage]+i*s->stride, a, s->stages[j]+i*s->stride);
        }
        solver->solve(m, s->cells, s->stride, s->jacobian, s->stages[stage]); ++stats->solves;
      }
      memcpy(s->new_y, s->variables, ny*sizeof(double)); memset(s->error, 0, ny*sizeof(double));
      for (int stage = 0; stage < p->stages; ++stage) chem_axpy(ny, s->new_y, p->m[stage], s->stages[stage]);
      for (int stage = 0; stage < p->stages; ++stage) chem_axpy(ny, s->error, p->e[stage], s->stages[stage]);
      if (m->nconstraints) for (size_t i = 0; i < m->nvariables; ++i) if (m->algebraic[i])
        for (size_t c = 0; c < s->cells; ++c) s->error[i*s->stride+c] = s->new_y[i*s->stride+c]-s->variables[i*s->stride+c];
      double error = chem_normalized_error(s, s->variables, s->new_y, s->error);
      double hnew = h*solver->controller(error, p, solver->controller_context);
      ++stats->internal_steps;
      if (isnan(error)) { status = CHEM_NAN; goto finished; }
      if (isinf(error)) { status = CHEM_INF; goto finished; }
      if (error < 1.0 || h < hmin) {
        ++stats->accepted; time += h;
        /* Keep allocation ownership and user-visible pointers stable. */
        memcpy(s->variables, s->new_y, ny*sizeof(double));
        hnew = fmax(hmin, fmin(hnew, hmax));
        if (reject_last) hnew = fmin(hnew, h);
        reject_last = reject_more = 0; h = hnew; accepted = 1;
      } else {
        if (reject_more) hnew = h*p->rejection_factor;
        reject_more = reject_last; reject_last = 1; h = hnew;
        if (stats->accepted >= 1) ++stats->rejected;
        chem_jacobian(s, s->variables, s->jacobian); ++stats->jacobian_updates;
        if (!(h > 0.0) || time+0.1*h == time) { status = CHEM_SMALL_STEP; goto finished; }
      }
    }
  }
finished:
  stats->final_time = time;
  return status;
}
