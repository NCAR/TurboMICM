#include "rosenbrock.h"
#include "linear.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

ChemRosenbrockMethod chem_rosenbrock_two_stage(void)
{
    ChemRosenbrockMethod p = {0};
    double g = 1.0 + 1.0/sqrt(2.0);
    p.stages = 2; p.gamma0 = g; p.error_order = 2;
    p.a[0] = 1.0/g; p.c[0] = -2.0/g;
    p.m[0] = 3.0/(2.0*g); p.m[1] = 1.0/(2.0*g);
    p.e[0] = p.e[1] = 1.0/(2.0*g);
    p.new_function[0] = p.new_function[1] = 1;
    return p;
}

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

ChemRosenbrockMethod chem_rosenbrock_four_stage(void)
{
    ChemRosenbrockMethod p = {0};
    p.stages = 4; p.gamma0 = 0.5728200000000000; p.error_order = 4;
    p.a[0]=0.2000000000000000E+01; p.a[1]=0.1867943637803922E+01;
    p.a[2]=0.2344449711399156; p.a[3]=p.a[1]; p.a[4]=p.a[2];
    p.c[0]=-0.7137615036412310E+01; p.c[1]=0.2580708087951457E+01;
    p.c[2]=0.6515950076447975; p.c[3]=-0.2137148994382534E+01;
    p.c[4]=-0.3214669691237626; p.c[5]=-0.6949742501781779;
    p.m[0]=0.2255570073418735E+01; p.m[1]=0.2870493262186792;
    p.m[2]=0.4353179431840180; p.m[3]=0.1093502252409163E+01;
    p.e[0]=-0.2815431932141155; p.e[1]=-0.7276199124938920E-01;
    p.e[2]=-0.1082196201495311; p.e[3]=-0.1093502252409163E+01;
    p.new_function[0]=p.new_function[1]=p.new_function[2]=1;
    return p;
}

ChemRosenbrockMethod chem_rosenbrock_four_stage_da(void)
{
    ChemRosenbrockMethod p = {0};
    p.stages=4; p.gamma0=.5; p.error_order=3;
    p.a[1]=2; p.a[3]=2; p.a[5]=1;
    p.c[0]=4; p.c[1]=1; p.c[2]=-1; p.c[3]=1; p.c[4]=-1; p.c[5]=-8.0/3.0;
    p.m[0]=2; p.m[2]=1; p.m[3]=1; p.e[3]=1;
    p.new_function[0]=p.new_function[2]=p.new_function[3]=1;
    return p;
}

ChemRosenbrockMethod chem_rosenbrock_six_stage_da(void)
{
    ChemRosenbrockMethod p = {0};
    p.stages=6; p.gamma0=.25; p.error_order=4;
    p.a[0]=0.1544000000000000E+01; p.a[1]=0.9466785280815826;
    p.a[2]=0.2557011698983284; p.a[3]=0.3314825187068521E+01;
    p.a[4]=0.2896124015972201E+01; p.a[5]=0.9986419139977817;
    p.a[6]=0.1221224509226641E+01; p.a[7]=0.6019134481288629E+01;
    p.a[8]=0.1253708332932087E+02; p.a[9]=-0.6878860361058950;
    p.a[10]=p.a[6]; p.a[11]=p.a[7]; p.a[12]=p.a[8]; p.a[13]=p.a[9]; p.a[14]=1;
    p.c[0]=-0.5668800000000000E+01; p.c[1]=-0.2430093356833875E+01;
    p.c[2]=-0.2063599157091915; p.c[3]=-0.1073529058151375;
    p.c[4]=-0.9594562251023355E+01; p.c[5]=-0.2047028614809616E+02;
    p.c[6]=0.7496443313967647E+01; p.c[7]=-0.1024680431464352E+02;
    p.c[8]=-0.3399990352819905E+02; p.c[9]=0.1170890893206160E+02;
    p.c[10]=0.8083246795921522E+01; p.c[11]=-0.7981132988064893E+01;
    p.c[12]=-0.3152159432874371E+02; p.c[13]=0.1631930543123136E+02;
    p.c[14]=-0.6058818238834054E+01;
    p.m[0]=p.a[6]; p.m[1]=p.a[7]; p.m[2]=p.a[8]; p.m[3]=p.a[9];
    p.m[4]=p.m[5]=1; p.e[5]=1;
    for (size_t i=0; i<6; ++i) p.new_function[i]=1;
    return p;
}

double chem_rosenbrock_normalized_error(size_t cells, size_t variables,
    const double *old_state, const double *new_state, const double *errors,
    const double *atol, double rtol)
{
    if (!cells || !variables || cells > SIZE_MAX/variables || !old_state ||
        !new_state || !errors || !atol) return NAN;
    double sum = 0;
    for (size_t i=0; i<variables; ++i)
        for (size_t cell=0; cell<cells; ++cell) {
            size_t k=cell*variables+i;
            double ymax=fmax(fabs(old_state[k]),fabs(new_state[k]));
            double scaled=errors[k]/(atol[i]+rtol*ymax);
            sum += scaled*scaled;
        }
    if (!isfinite(sum)) return NAN;
    return fmax(sqrt(sum/(cells*variables)),1e-10);
}

ChemRosenbrockOptions chem_rosenbrock_defaults(const double *atol)
{
    ChemRosenbrockOptions p = {1e-6, atol, 0, 0, 0, 0.2, 6, 0.9, 0.1, 1000, NULL};
    return p;
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
        o->h_min < 0 || o->h_max < 0 || o->h_start < 0 ||
        (o->sparse_plan && chem_sparse_plan_size(o->sparse_plan) != n))
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
            if (!(o->sparse_plan ? chem_sparse_factor(o->sparse_plan, matrix) : chem_dense_factor(matrix, n))) { status = CHEM_SOLVE_SINGULAR; goto done; }
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
                if (o->sparse_plan) chem_sparse_solve(o->sparse_plan, matrix, ks);
                else chem_dense_solve(matrix, n, ks);
                ++stats->solves;
            }
            memcpy(next, y, n*sizeof(*next));
            memset(err, 0, n*sizeof(*err));
            for (size_t s = 0; s < p->stages; ++s)
                for (size_t i = 0; i < n; ++i) {
                    next[i] += p->m[s]*stage[s*n+i];
                    err[i] += p->e[s]*stage[s*n+i];
                }
            double error = chem_rosenbrock_normalized_error(1,n,y,next,err,o->absolute_tolerance,o->relative_tolerance);
            ++stats->steps;
            if (!isfinite(error)) { status = CHEM_SOLVE_NONFINITE; goto done; }
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
    stats->final_time = t;
    /* The CPU solver clamps ODE state after each public solve call. */
    for (size_t i = 0; i < n; ++i) y[i] = y[i] > 0.0 ? y[i] : 0.0;
    free(work);
    return status;
}
