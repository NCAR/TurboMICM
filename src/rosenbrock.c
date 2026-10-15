/* Rosenbrock integration and error control ported from the reference.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <float.h>
#include <math.h>
#include <string.h>

static void axpy(size_t n, double a, const double *restrict x, double *restrict y)
{
    #pragma omp simd
    for (size_t i = 0; i < n; ++i) y[i] += a * x[i];
}
static double normalized_error(const ChemState *s, const double *ynew, const double *error)
{
    double sum = 0;
    size_t nc = s->cells, n = s->mechanism->data->species_count;
    /* Strict species, then cell order, matching the reference's error reduction. */
    for (size_t i = 0; i < n; ++i)
        for (size_t c = 0; c < nc; ++c) {
            size_t j = i * nc + c;
            double ymax = fmax(fabs(s->y[j]), fabs(ynew[j]));
            double e = error[j] / (s->atol[i] + s->rtol * ymax);
            sum += e * e;
        }
    return fmax(sqrt(sum / (n * nc)), 1e-10);
}
ChemStatus chem_solve(ChemState *s, const ChemRosenbrock *p, double dt, ChemStats *stats)
{
    if (!s || !p || !stats || !(dt > 0) || !isfinite(dt) || p->stages < 1 || p->stages > 6
        || !(p->gamma[0] > 0) || !(s->rtol >= 0)) return CHEM_INVALID;
    *stats = (ChemStats){0};
    const ChemMechanism *m = s->mechanism;
    size_t n = m->data->species_count, nc = s->cells, size = n * nc;
    for (size_t i = 0; i < n; ++i) if (!(s->atol[i] > 0)) return CHEM_INVALID;
    double hmin = p->h_min == 0 ? 1e-15 * dt : p->h_min;
    double hmax = p->h_max == 0 ? dt : fmin(dt, p->h_max);
    double hstart = p->h_start == 0 ? 1e-6 * dt : fmin(hmax, p->h_start);
    double h = fmin(fmax(hmin, fabs(hstart)), fabs(hmax)), time = 0;
    int reject_last = 0, reject_more = 0;
    while (time - dt + DBL_EPSILON <= 0) {
        stats->final_time = time;
        if (stats->internal_steps > p->max_steps) return CHEM_MAX_STEPS;
        if (time + .1 * h == time || h <= fmin(DBL_EPSILON, hmin)) return CHEM_SMALL_STEP;
        h = fmin(h, fabs(dt - time));
        chem_forcing(s, s->y, s->forcing);
        ++stats->function_calls;
        chem_negative_jacobian(s, s->y, s->jacobian);
        ++stats->jacobian_updates;
        for (;;) {
            double alpha = 1.0 / (h * p->gamma[0]);
            for (size_t i = 0; i < n; ++i) {
                double *restrict diag = s->jacobian + m->diagonal[i] * nc;
                #pragma omp simd
                for (size_t c = 0; c < nc; ++c) diag[c] += alpha;
            }
            ChemStatus status = chem_factor(m, nc, s->jacobian);
            ++stats->decompositions;
            if (status != CHEM_OK) return status;
            for (int stage = 0; stage < p->stages; ++stage) {
                double *k = s->stages + (size_t)stage * size;
                int offset = stage * (stage - 1) / 2;
                if (stage == 0) memcpy(k, s->forcing, size * sizeof(double));
                else if (p->new_function[stage]) {
                    memcpy(s->candidate, s->y, size * sizeof(double));
                    for (int j = 0; j < stage; ++j)
                        axpy(size, p->a[offset + j], s->stages + (size_t)j * size, s->candidate);
                    chem_forcing(s, s->candidate, k);
                    ++stats->function_calls;
                }
                if (stage + 1 < p->stages && !p->new_function[stage + 1])
                    memcpy(k + size, k, size * sizeof(double));
                for (int j = 0; j < stage; ++j)
                    axpy(size, p->c[offset + j] / h, s->stages + (size_t)j * size, k);
                chem_linear_solve(m, nc, s->jacobian, k);
                ++stats->solves;
            }
            memcpy(s->candidate, s->y, size * sizeof(double));
            memset(s->error, 0, size * sizeof(double));
            for (int j = 0; j < p->stages; ++j)
                axpy(size, p->m[j], s->stages + (size_t)j * size, s->candidate);
            for (int j = 0; j < p->stages; ++j)
                axpy(size, p->e[j], s->stages + (size_t)j * size, s->error);
            double error = normalized_error(s, s->candidate, s->error);
            ++stats->internal_steps;
            /* Inspect the state too: fmax follows C's NaN suppression rules. */
            if (!isfinite(error)) return CHEM_NONFINITE;
            for (size_t i = 0; i < size; ++i)
                if (!isfinite(s->candidate[i]) || !isfinite(s->error[i])) return CHEM_NONFINITE;
            double factor = fmin(p->factor_max, fmax(p->factor_min, p->safety / pow(error, 1.0 / p->order)));
            double hnew = h * factor;
            if (error < 1 || h < hmin) {
                ++stats->accepted;
                time += h;
                double *swap = s->y; s->y = s->candidate; s->candidate = swap;
                hnew = fmax(hmin, fmin(hnew, hmax));
                if (reject_last) hnew = fmin(hnew, h);
                reject_last = reject_more = 0;
                h = hnew;
                break;
            }
            if (reject_more) hnew = h * p->rejection_factor;
            reject_more = reject_last;
            reject_last = 1;
            h = hnew;
            if (stats->accepted >= 1) ++stats->rejected;
            chem_negative_jacobian(s, s->y, s->jacobian);
            ++stats->jacobian_updates;
            /* A rejection can otherwise underflow without returning to the outer loop. */
            if (time + .1 * h == time || h <= 0) return CHEM_SMALL_STEP;
        }
    }
    stats->final_time = time;
    return CHEM_OK;
}
