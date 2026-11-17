#include "internal.h"

#include <float.h>
#include <math.h>
#include <string.h>

ChemRosenbrock chem_rosenbrock_three_stage(void)
{
    ChemRosenbrock p = {0};
    p.stages = 3;
    p.a[0] = 1.0; p.a[1] = 1.0;
    p.c[0] = -0.10156171083877702091975600115545e+01;
    p.c[1] =  0.40759956452537699824805835358067e+01;
    p.c[2] =  0.92076794298330791242156818474003e+01;
    p.m[0] = 1.0; p.m[1] = 0.61697947043828245592553615689730e+01;
    p.m[2] = -0.42772256543218573326238373806514;
    p.e[0] = 0.5; p.e[1] = -0.29079558716805469821718236208017e+01;
    p.e[2] = 0.22354069897811569627360909276199;
    p.gamma[0] = 0.43586652150845899941601945119356;
    p.new_function[0] = p.new_function[1] = 1;
    p.order = 3.0; p.relative_tolerance = 1e-6;
    p.factor_min = 0.2; p.factor_max = 6.0; p.safety = 0.9;
    p.rejection_factor = 0.1; p.max_steps = 1000;
    return p;
}

static void axpy(size_t n, double *restrict dst, const double *restrict src, double a)
{
    #pragma omp simd
    for (size_t i = 0; i < n; ++i) dst[i] += a*src[i];
}

static double normalized_error(const ChemState *s, const ChemRosenbrock *p)
{
    const ChemMechanism *m = s->mechanism;
    size_t cells = s->cells;
    double sum = 0;
    for (size_t i = 0; i < m->species_count; ++i)
        for (size_t c = 0; c < cells; ++c) {
            size_t at = i*cells+c;
            double scale = m->species[i].absolute_tolerance +
                p->relative_tolerance * fmax(fabs(s->values[at]), fabs(s->candidate[at]));
            double ratio = s->error[at] / scale;
            sum += ratio*ratio;
        }
    double error = sqrt(sum / (m->species_count*cells));
    return error < 1e-10 ? 1e-10 : error;
}

int chem_solve(ChemState *s, const ChemRosenbrock *p, double dt, ChemStats *result)
{
    if (!s || !p || !result || !isfinite(dt) || dt < 0 || p->stages < 1 || p->stages > 6 ||
        p->gamma[0] == 0 || p->order <= 0 || p->relative_tolerance < 0 ||
        p->factor_min <= 0 || p->factor_max < p->factor_min || p->safety <= 0 ||
        p->rejection_factor <= 0 || p->h_min < 0 || p->h_max < 0 || p->h_start < 0)
        return CHEM_INVALID;
    memset(result, 0, sizeof(*result));
    if (dt == 0) return CHEM_OK;
    size_t ny = s->mechanism->species_count * s->cells;
    double h_min = p->h_min == 0.0 ? 1e-15*dt : p->h_min;
    double h_max = p->h_max == 0.0 ? dt : fmin(dt, p->h_max);
    double h_start = p->h_start == 0.0 ? 1e-6*dt : fmin(h_max, p->h_start);
    double h = fmin(fmax(h_min, fabs(h_start)), fabs(h_max));
    double time = 0;
    int reject_last = 0, reject_more = 0;
    while (time - dt + DBL_EPSILON <= 0) {
        if (result->internal_steps > p->max_steps) return CHEM_TOO_MANY_STEPS;
        if (time + 0.1*h == time || h <= fmin(DBL_EPSILON, h_min)) return CHEM_STEP_TOO_SMALL;
        h = fmin(h, fabs(dt-time));
        chem_forcing(s, s->values, s->forcing, s->scratch);
        ++result->function_calls;
        chem_negative_jacobian(s, s->values);
        ++result->jacobian_updates;
        int accepted = 0;
        while (!accepted) {
            int status = chem_factor(s, 1.0 / (h*p->gamma[0]));
            ++result->decompositions;
            if (status != CHEM_OK) return status;
            for (size_t stage = 0; stage < p->stages; ++stage) {
                size_t first = stage ? stage*(stage-1)/2 : 0;
                double *k = s->stages + stage*ny;
                if (stage == 0) {
                    memcpy(k, s->forcing, ny*sizeof(*k));
                } else if (p->new_function[stage]) {
                    memcpy(s->candidate, s->values, ny*sizeof(*k));
                    for (size_t j = 0; j < stage; ++j)
                        axpy(ny, s->candidate, s->stages+j*ny, p->a[first+j]);
                    chem_forcing(s, s->candidate, k, s->scratch);
                    ++result->function_calls;
                }
                if (stage+1 < p->stages && !p->new_function[stage+1])
                    memcpy(s->stages+(stage+1)*ny, k, ny*sizeof(*k));
                for (size_t j = 0; j < stage; ++j) axpy(ny, k, s->stages+j*ny, p->c[first+j]/h);
                chem_substitute(s, k);
                ++result->solves;
            }
            memcpy(s->candidate, s->values, ny*sizeof(*s->candidate));
            memset(s->error, 0, ny*sizeof(*s->error));
            for (size_t stage = 0; stage < p->stages; ++stage) {
                axpy(ny, s->candidate, s->stages+stage*ny, p->m[stage]);
                axpy(ny, s->error, s->stages+stage*ny, p->e[stage]);
            }
            double error = normalized_error(s, p);
            if (!isfinite(error)) return CHEM_NONFINITE;
            double factor = fmin(p->factor_max, fmax(p->factor_min, p->safety / pow(error, 1.0/p->order)));
            double h_new = h*factor;
            ++result->internal_steps;
            if (error < 1 || h < h_min) {
                ++result->accepted;
                time += h;
                double *swap = s->values; s->values = s->candidate; s->candidate = swap;
                h_new = fmax(h_min, fmin(h_new, h_max));
                if (reject_last) h_new = fmin(h_new, h);
                reject_last = reject_more = 0;
                h = h_new;
                accepted = 1;
            } else {
                if (reject_more) h_new = h*p->rejection_factor;
                reject_more = reject_last;
                reject_last = 1;
                h = h_new;
                if (result->accepted >= 1) ++result->rejected;
                chem_negative_jacobian(s, s->values);
                ++result->jacobian_updates;
            }
        }
    }
    result->final_time = time;
    return CHEM_OK;
}
