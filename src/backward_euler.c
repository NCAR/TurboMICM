/* Translated numerical method: Copyright (C) 2023-2026 University Corporation for Atmospheric Research.
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

ChemEuler chem_euler_parameters(void) {
  ChemEuler p = {0};
  p.small = 1e-40; p.max_steps = 11; p.nreductions = 5;
  for (size_t i = 0; i < 4; ++i) p.reductions[i] = 0.5;
  p.reductions[4] = 0.1; p.factor = chem_lu_mozart; p.solve = chem_linear_solve;
  return p;
}
int chem_euler_converged(const ChemState *s, double small, const double *residual, const double *y) {
  for (size_t i = 0; i < s->model->nvariables; ++i) for (size_t c = 0; c < s->cells; ++c) {
    size_t at = i*s->stride+c; double v = fabs(residual[at]);
    if (!isfinite(v) || (v > small && v > s->absolute_tolerance[i] && v > s->relative_tolerance*fabs(y[at]))) return 0;
  }
  return 1;
}
ChemStatus chem_backward_euler(const ChemEuler *p, double dt, ChemState *s, ChemStats *stats) {
  if (!p || !s || !stats || !isfinite(dt) || dt < 0.0 || !p->factor || !p->solve || p->nreductions > 5) return CHEM_INVALID;
  memset(stats, 0, sizeof *stats);
  double h = p->h_start == 0.0 ? dt : p->h_start, time = 0.0;
  size_t successes = 0, failures = 0, ny = s->model->nvariables*s->stride;
  double *old = s->new_y, *residual = s->forcing;
  memcpy(old, s->variables, ny*sizeof(double));
  ChemStatus status = CHEM_OK;
  while (time < dt) {
    size_t iterations = 0; int converged = 0;
    do {
      ++stats->internal_steps;
      chem_forcing(s, s->variables, residual); ++stats->function_calls;
      chem_jacobian(s, s->variables, s->jacobian); ++stats->jacobian_updates;
      for (size_t i = 0; i < s->model->nvariables; ++i)
        for (size_t c = 0; c < s->cells; ++c) s->jacobian[s->model->diagonal[i]*s->stride+c] += 1.0/h;
      p->factor(s->model, s->cells, s->stride, s->jacobian); ++stats->decompositions;
      #pragma omp simd
      for (size_t i = 0; i < ny; ++i) residual[i] -= (s->variables[i]-old[i])/h;
      p->solve(s->model, s->cells, s->stride, s->jacobian, residual); ++stats->solves;
      #pragma omp simd
      for (size_t i = 0; i < ny; ++i) {
        double y = s->variables[i]+residual[i]; s->variables[i] = 0.0 > y ? 0.0 : y;
      }
      if (iterations++ == 0) continue;
      converged = chem_euler_converged(s, p->small, residual, s->variables);
    } while (!converged && iterations < p->max_steps);
    if (!converged) {
      ++stats->rejected; successes = 0;
      if (failures >= p->nreductions) {
        time += h; status = CHEM_UNCONVERGED;
        for (size_t i = 0; i < ny; ++i) if (isnan(s->variables[i])) { status = CHEM_NAN; goto done; }
        for (size_t i = 0; i < ny; ++i) if (isinf(s->variables[i])) { status = CHEM_INF; goto done; }
        break;
      }
      memcpy(s->variables, old, ny*sizeof(double)); h *= p->reductions[failures++];
    } else {
      ++stats->accepted; time += h; memcpy(old, s->variables, ny*sizeof(double));
      if (++successes >= 2) { successes = 0; h *= 2.0; }
    }
    h = fmin(h, dt-time);
    if (!(h > 0.0) && time < dt) { status = CHEM_SMALL_STEP; break; }
  }
done:
  stats->final_time = time; return status;
}
