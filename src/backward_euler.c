/* Fully implicit Euler with the reference Newton loop and retry policy.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <math.h>
#include <string.h>

ChemBackwardEuler chem_backward_euler_parameters(void)
{
    return (ChemBackwardEuler){.small=1e-40, .max_iterations=11, .reduction_count=5, .reductions={.5,.5,.5,.5,.1}};
}
int chem_backward_euler_converged(const ChemState *s, const ChemBackwardEuler *p, const double *residual, const double *y)
{
    size_t n = s->mechanism->data->species_count, nc = s->cells;
    for (size_t i = 0; i < n; ++i)
        for (size_t c = 0; c < nc; ++c) {
            size_t j = i * nc + c;
            double e = fabs(residual[j]);
            if (!isfinite(e) || (e > p->small && e > s->atol[i] && e > s->rtol * fabs(y[j]))) return 0;
        }
    return 1;
}
static ChemStatus integrate_backward_euler(ChemState *s, const ChemBackwardEuler *p, double dt, ChemStats *stats)
{
    if (!s || !p || !stats || !(dt > 0) || !isfinite(dt) || !p->max_iterations || p->reduction_count > 5) return CHEM_INVALID;
    *stats = (ChemStats){0};
    const ChemMechanism *m = s->mechanism;
    size_t nc = s->cells, n = m->data->species_count, size = n * nc;
    double *old = s->candidate, *f = s->forcing;
    memcpy(old, s->y, size * sizeof(double));
    double h = p->h_start == 0 ? dt : p->h_start, time = 0;
    if (!(h > 0) || !isfinite(h)) return CHEM_INVALID;
    size_t successes = 0, failures = 0;
    ChemStatus result = CHEM_OK;
    while (time < dt) {
        int converged = 0;
        size_t iteration = 0;
        do {
            ++stats->internal_steps;
            chem_forcing(s, s->y, f); ++stats->function_calls;
            chem_negative_jacobian(s, s->y, s->jacobian); ++stats->jacobian_updates;
            for (size_t i = 0; i < n; ++i)
                for (size_t c = 0; c < nc; ++c) s->jacobian[m->diagonal[i] * nc + c] += 1 / h;
            ChemStatus status = s->linear.factor(m, nc, s->jacobian);
            ++stats->decompositions;
            if (status != CHEM_OK) return status;
            for (size_t j = 0; j < size; ++j) f[j] -= (s->y[j] - old[j]) / h;
            s->linear.solve(m, nc, s->jacobian, f); ++stats->solves;
            for (size_t j = 0; j < size; ++j) {
                double value = s->y[j] + f[j];
                s->y[j] = value < 0 ? 0 : value;
            }
            if (iteration++ == 0) continue;
            converged = chem_backward_euler_converged(s, p, f, s->y);
        } while (!converged && iteration < p->max_iterations);
        if (!converged) {
            ++stats->rejected; successes = 0;
            if (failures >= p->reduction_count) {
                time += h;
                result = CHEM_UNCONVERGED;
                for (size_t j = 0; j < size; ++j) if (!isfinite(s->y[j])) result = CHEM_NONFINITE;
                break;
            }
            memcpy(s->y, old, size * sizeof(double));
            double reduction = p->reductions[failures++];
            if (!(reduction > 0) || !(reduction < 1)) return CHEM_INVALID;
            h *= reduction;
        } else {
            ++stats->accepted; time += h;
            memcpy(old, s->y, size * sizeof(double));
            if (++successes >= 2) { successes = 0; h *= 2; }
        }
        h = fmin(h, dt - time);
    }
    stats->final_time = time;
    return result;
}

ChemStatus chem_solve_backward_euler(ChemState *s, const ChemBackwardEuler *p, double dt, ChemStats *stats)
{
    ChemStatus status = integrate_backward_euler(s, p, dt, stats);
    if (s && status != CHEM_INVALID) chem_post_solve_clamp(s);
    return status;
}
