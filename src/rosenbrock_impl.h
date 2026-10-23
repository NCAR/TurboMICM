/* Shared Rosenbrock algorithm. Generated kernels bind these hooks at compile time.
 * Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <float.h>
#include <math.h>
#include <string.h>
#ifndef CHEM_NORMALIZED_ERROR
#define CHEM_NORMALIZED_ERROR chem_normalized_error
#endif
#ifndef CHEM_SOLVE
#define CHEM_SOLVE chem_solve
#endif
#ifndef CHEM_FORCING
#define CHEM_FORCING chem_forcing
#endif
#ifndef CHEM_JACOBIAN
#define CHEM_JACOBIAN chem_negative_jacobian
#endif
#ifndef CHEM_ALPHA
#define CHEM_ALPHA chem_alpha_minus_jacobian
#endif
#ifndef CHEM_FACTOR
#define CHEM_FACTOR(m,nc,a) s->linear.factor(m,nc,a)
#endif
#ifndef CHEM_BACKSOLVE
#define CHEM_BACKSOLVE(m,nc,a,x) s->linear.solve(m,nc,a,x)
#endif
#ifndef CHEM_SPECIES
#define CHEM_SPECIES(s) s->mechanism->data->species_count
#endif
#ifndef CHEM_CONSTRAINED
#define CHEM_CONSTRAINED(m) m->data->constraint_count != 0
#endif
#ifndef CHEM_VALID_STATE
#define CHEM_VALID_STATE(s) 1
#endif
static void axpy(size_t n, double a, const double *restrict x, double *restrict y)
{
    #pragma omp simd
    for (size_t i = 0; i < n; ++i) y[i] += a * x[i];
}
static ChemStatus integrate(ChemState *s, const ChemRosenbrock *p, double dt, ChemStats *stats)
{
    if (!s || !s->mechanism || !s->cells || !CHEM_VALID_STATE(s) || !p || !stats || !(dt > 0) || !isfinite(dt) || p->stages < 1 || p->stages > 6
        || !(p->gamma[0] > 0) || !(s->rtol >= 0)) return CHEM_INVALID;
    *stats = (ChemStats){0};
    const ChemMechanism *m = s->mechanism;
    size_t n = CHEM_SPECIES(s), nc = s->cells, size = n * nc;
    for (size_t i = 0; i < n; ++i) if (!(s->atol[i] > 0)) return CHEM_INVALID;
    ChemStatus initialization = chem_initialize_constraints(s, p, stats);
    if (initialization != CHEM_OK) return initialization;
    int constrained = CHEM_CONSTRAINED(m);
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
        CHEM_FORCING(s, s->y, s->forcing);
        if (constrained) chem_constraint_forcing(s, s->y, s->forcing);
        ++stats->function_calls;
        CHEM_JACOBIAN(s, s->y, s->jacobian);
        if (constrained) chem_constraint_negative_jacobian(s, s->y, s->jacobian);
        ++stats->jacobian_updates;
        for (;;) {
            double alpha = 1.0 / (h * p->gamma[0]);
            CHEM_ALPHA(m, nc, s->jacobian, alpha);
            ChemStatus status = CHEM_FACTOR(m, nc, s->jacobian);
            ++stats->decompositions;
            if (status != CHEM_OK) return status;
#ifdef CHEM_FUSED_STAGES
            if (p->stages == 3 && p->new_function[1] && !p->new_function[2]) {
                double *restrict k0 = s->stages;
                double *restrict k1 = s->stages + size;
                double *restrict k2 = s->stages + 2 * size;
                double *restrict candidate = s->candidate, *restrict estimator = s->error;
                const double *restrict old = s->y;
                memcpy(k0, s->forcing, size * sizeof(double));
                CHEM_BACKSOLVE(m, nc, s->jacobian, k0);
                #pragma omp simd
                for (size_t i = 0; i < size; ++i) {
                    double v = old[i];
                    v += p->a[0] * k0[i];
                    candidate[i] = v;
                }
                CHEM_FORCING(s, candidate, k1);
                ++stats->function_calls;
                double c0 = p->c[0] / h, c1 = p->c[1] / h, c2 = p->c[2] / h;
                #pragma omp simd
                for (size_t i = 0; i < size; ++i) {
                    k2[i] = k1[i];
                    k1[i] += c0 * k0[i];
                }
                CHEM_BACKSOLVE(m, nc, s->jacobian, k1);
                #pragma omp simd
                for (size_t i = 0; i < size; ++i) {
                    k2[i] += c1 * k0[i];
                    k2[i] += c2 * k1[i];
                }
                CHEM_BACKSOLVE(m, nc, s->jacobian, k2);
                stats->solves += 3;
                #pragma omp simd
                for (size_t i = 0; i < size; ++i) {
                    double v = old[i], error = 0;
                    v += p->m[0] * k0[i];
                    v += p->m[1] * k1[i];
                    v += p->m[2] * k2[i];
                    error += p->e[0] * k0[i];
                    error += p->e[1] * k1[i];
                    error += p->e[2] * k2[i];
                    candidate[i] = v;
                    estimator[i] = error;
                }
            } else
#endif
            {
            for (int stage = 0; stage < p->stages; ++stage) {
                double *k = s->stages + (size_t)stage * size;
                int offset = stage * (stage - 1) / 2;
                if (stage == 0) memcpy(k, s->forcing, size * sizeof(double));
                else if (p->new_function[stage]) {
                    memcpy(s->candidate, s->y, size * sizeof(double));
                    for (int j = 0; j < stage; ++j)
                        axpy(size, p->a[offset + j], s->stages + (size_t)j * size, s->candidate);
                    CHEM_FORCING(s, s->candidate, k);
                    if (constrained) chem_constraint_forcing(s, s->candidate, k);
                    ++stats->function_calls;
                }
                if (stage + 1 < p->stages && !p->new_function[stage + 1])
                    memcpy(k + size, k, size * sizeof(double));
                for (int j = 0; j < stage; ++j) {
                    double coefficient = p->c[offset + j] / h;
                    const double *previous = s->stages + (size_t)j * size;
                    if (!constrained) axpy(size, coefficient, previous, k);
                    else for (size_t i = 0; i < n; ++i) if (!m->algebraic[i])
                        axpy(nc, coefficient, previous + i * nc, k + i * nc);
                }
                CHEM_BACKSOLVE(m, nc, s->jacobian, k);
                ++stats->solves;
            }
            memcpy(s->candidate, s->y, size * sizeof(double));
            memset(s->error, 0, size * sizeof(double));
            for (int j = 0; j < p->stages; ++j)
                axpy(size, p->m[j], s->stages + (size_t)j * size, s->candidate);
            for (int j = 0; j < p->stages; ++j)
                axpy(size, p->e[j], s->stages + (size_t)j * size, s->error);
            }
            if (constrained) for (size_t i = 0; i < n; ++i) if (m->algebraic[i])
                for (size_t c = 0; c < nc; ++c) {
                    size_t j = i * nc + c;
                    s->error[j] = s->candidate[j] - s->y[j];
                }
            double error = CHEM_NORMALIZED_ERROR(s, s->y, s->candidate, s->error);
            ++stats->internal_steps;
            /* Inspect the state too: fmax follows C's NaN suppression rules. */
            if (!isfinite(error)) return CHEM_NONFINITE;
            for (size_t i = 0; i < size; ++i)
                if (!isfinite(s->candidate[i]) || !isfinite(s->error[i])) return CHEM_NONFINITE;
            double factor = p->controller
                ? p->controller(error, p->order, p->safety, p->factor_min, p->factor_max, p->controller_context)
                : fmin(p->factor_max, fmax(p->factor_min, p->safety / pow(error, 1.0 / p->order)));
            if (!(factor > 0) || !isfinite(factor)) return CHEM_INVALID;
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
            CHEM_JACOBIAN(s, s->y, s->jacobian);
        if (constrained) chem_constraint_negative_jacobian(s, s->y, s->jacobian);
            ++stats->jacobian_updates;
            /* A rejection can otherwise underflow without returning to the outer loop. */
            if (time + .1 * h == time || h <= 0) return CHEM_SMALL_STEP;
        }
    }
    stats->final_time = time;
    return CHEM_OK;
}

ChemStatus CHEM_SOLVE(ChemState *s, const ChemRosenbrock *p, double dt, ChemStats *stats)
{
    ChemStatus status = integrate(s, p, dt, stats);
    if (s && status != CHEM_INVALID) chem_post_solve_clamp(s);
    return status;
}

#undef CHEM_NORMALIZED_ERROR
#undef CHEM_FUSED_STAGES
#undef CHEM_SOLVE
#undef CHEM_FORCING
#undef CHEM_JACOBIAN
#undef CHEM_ALPHA
#undef CHEM_FACTOR
#undef CHEM_BACKSOLVE
#undef CHEM_SPECIES
#undef CHEM_CONSTRAINED
#undef CHEM_VALID_STATE
