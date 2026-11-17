#include "rosenbrock.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

ChemRosenbrockMethod chem_rosenbrock_three_stage(void)
{
    ChemRosenbrockMethod p = {0};
    p.stages = 3;
    p.a[0] = 1; p.a[1] = 1;
    p.c[0] = -0.10156171083877702091975600115545e+01;
    p.c[1] = 0.40759956452537699824805835358067e+01;
    p.c[2] = 0.92076794298330791242156818474003e+01;
    p.m[0] = 1; p.m[1] = 0.61697947043828245592553615689730e+01;
    p.m[2] = -0.42772256543218573326238373806514;
    p.e[0] = 0.5; p.e[1] = -0.29079558716805469821718236208017e+01;
    p.e[2] = 0.22354069897811569627360909276199;
    p.gamma0 = 0.43586652150845899941601945119356;
    p.error_order = 3;
    p.new_function[0] = p.new_function[1] = 1;
    return p;
}

ChemRosenbrockOptions chem_rosenbrock_defaults(const double *atol)
{
    ChemRosenbrockOptions p = {1e-6, atol, 0, 0, 0, 0.2, 6, 0.9, 0.1, 1000};
    return p;
}

/* Dense form of the no-pivot Mozart elimination: scale a column, then
   update each trailing column. The factor is stored in place. */
static int factor(double *a, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        if (a[i*n+i] == 0 || !isfinite(a[i*n+i])) return 0;
        double inverse = 1.0 / a[i*n+i];
        for (size_t j = i+1; j < n; ++j) a[j*n+i] *= inverse;
        for (size_t k = i+1; k < n; ++k)
            for (size_t j = i+1; j < n; ++j) a[j*n+k] -= a[j*n+i]*a[i*n+k];
    }
    return 1;
}

static void solve(const double *a, size_t n, double *b)
{
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < i; ++j) b[i] -= a[i*n+j]*b[j];
    for (size_t i = n; i-- > 0;) {
        for (size_t j = i+1; j < n; ++j) b[i] -= a[i*n+j]*b[j];
        b[i] /= a[i*n+i];
    }
}

ChemSolveStatus chem_rosenbrock_solve(size_t n, double *y, double duration,
    ChemOdeFunction rhs, ChemOdeFunction jacobian, void *context,
    const ChemRosenbrockMethod *p, const ChemRosenbrockOptions *o,
    ChemSolverStats *stats)
{
    ChemSolverStats local_stats = {0};
    if (!stats) stats = &local_stats;
    memset(stats, 0, sizeof(*stats));
    if (!n || !y || !rhs || !jacobian || !p || !o || !o->absolute_tolerance ||
        !isfinite(duration) || duration <= 0 || !p->stages || p->stages > 6 ||
        !isfinite(p->gamma0) || p->gamma0 == 0 || !isfinite(p->error_order) || p->error_order <= 0 ||
        !isfinite(o->relative_tolerance) || o->relative_tolerance < 0 ||
        !isfinite(o->factor_min) || !isfinite(o->factor_max) || !isfinite(o->safety) ||
        o->factor_min <= 0 || o->factor_max < o->factor_min || o->safety <= 0 ||
        !isfinite(o->rejection_factor) || o->rejection_factor <= 0 ||
        !isfinite(o->h_min) || !isfinite(o->h_max) || !isfinite(o->h_start) ||
        o->h_min < 0 || o->h_max < 0 || o->h_start < 0)
        return CHEM_SOLVE_INVALID;
    for (size_t i = 0; i < n; ++i)
        if (!isfinite(o->absolute_tolerance[i]) || !(o->absolute_tolerance[i] > 0)) return CHEM_SOLVE_INVALID;
    if (n > SIZE_MAX / n || n*n > SIZE_MAX/sizeof(double) ||
        n > (SIZE_MAX/sizeof(double) - n*n)/10) return CHEM_SOLVE_ALLOCATION;
    double *work = calloc(n*n + 10*n, sizeof(*work));
    if (!work) return CHEM_SOLVE_ALLOCATION;
    double *matrix = work, *initial = matrix+n*n, *next = initial+n;
    double *stage = next+n, *forcing = stage+6*n, *err = forcing+n;
    double hmin = o->h_min == 0 ? 1e-15*duration : o->h_min;
    double hmax = o->h_max == 0 ? duration : fmin(duration, o->h_max);
    double hstart = o->h_start == 0 ? 1e-6*duration : fmin(hmax, o->h_start);
    double h = fmin(fmax(hmin, fabs(hstart)), fabs(hmax));
    double t = 0;
    int rejected_last = 0, rejected_more = 0;
    ChemSolveStatus status = CHEM_SOLVE_OK;
    while ((t - duration + DBL_EPSILON) <= 0) {
        if (stats->steps > o->max_steps) { status = CHEM_SOLVE_MAX_STEPS; break; }
        if (t + 0.1*h == t || h <= fmin(DBL_EPSILON, hmin)) {
            status = CHEM_SOLVE_STEP_TOO_SMALL; break;
        }
        h = fmin(h, fabs(duration-t));
        rhs(context, y, initial); ++stats->function_calls;
        int accepted = 0;
        while (!accepted) {
            jacobian(context, y, matrix); ++stats->jacobian_updates;
            for (size_t i = 0; i < n*n; ++i) matrix[i] = -matrix[i];
            for (size_t i = 0; i < n; ++i) matrix[i*n+i] += 1.0/(h*p->gamma0);
            ++stats->decompositions;
            if (!factor(matrix, n)) { status = CHEM_SOLVE_SINGULAR; goto done; }
            for (size_t s = 0; s < p->stages; ++s) {
                double *ks = stage + s*n;
                size_t offset = s*(s-1)/2;
                if (s == 0) memcpy(ks, initial, n*sizeof(*ks));
                else if (p->new_function[s]) {
                    memcpy(next, y, n*sizeof(*next));
                    for (size_t j = 0; j < s; ++j)
                        for (size_t i = 0; i < n; ++i) next[i] += p->a[offset+j]*stage[j*n+i];
                    rhs(context, next, ks); ++stats->function_calls;
                }
                if (s+1 < p->stages && !p->new_function[s+1])
                    memcpy(stage+(s+1)*n, ks, n*sizeof(*ks));
                for (size_t j = 0; j < s; ++j)
                    for (size_t i = 0; i < n; ++i) ks[i] += (p->c[offset+j]/h)*stage[j*n+i];
                solve(matrix, n, ks); ++stats->solves;
            }
            memcpy(next, y, n*sizeof(*next));
            memset(err, 0, n*sizeof(*err));
            for (size_t s = 0; s < p->stages; ++s)
                for (size_t i = 0; i < n; ++i) {
                    next[i] += p->m[s]*stage[s*n+i];
                    err[i] += p->e[s]*stage[s*n+i];
                }
            double sum = 0;
            for (size_t i = 0; i < n; ++i) {
                double scale = o->absolute_tolerance[i] + o->relative_tolerance*fmax(fabs(y[i]), fabs(next[i]));
                double e = err[i]/scale;
                sum += e*e;
            }
            double error = fmax(sqrt(sum/n), 1e-10);
            ++stats->steps;
            if (!isfinite(sum)) { status = CHEM_SOLVE_NONFINITE; goto done; }
            double fac = fmin(o->factor_max, fmax(o->factor_min, o->safety/pow(error, 1.0/p->error_order)));
            double hnew = h*fac;
            if (error < 1 || h < hmin) {
                ++stats->accepted; t += h; memcpy(y, next, n*sizeof(*y));
                hnew = fmax(hmin, fmin(hnew, hmax));
                if (rejected_last) hnew = fmin(hnew, h);
                rejected_last = rejected_more = 0;
                h = hnew; accepted = 1;
            } else {
                if (rejected_more) hnew = h*o->rejection_factor;
                rejected_more = rejected_last; rejected_last = 1; h = hnew;
                if (stats->accepted >= 1) ++stats->rejected;
            }
        }
    }
done:
    /* The CPU solver clamps ODE state after each public solve call. */
    for (size_t i = 0; i < n; ++i) y[i] = y[i] > 0.0 ? y[i] : 0.0;
    free(work);
    return status;
}
