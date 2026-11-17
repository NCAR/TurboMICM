#include "internal.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static ChemRosenbrock defaults(void)
{
    ChemRosenbrock p = {0};
    p.relative_tolerance = 1e-6;
    p.factor_min = 0.2; p.factor_max = 6.0; p.safety = 0.9;
    p.rejection_factor = 0.1; p.max_steps = 1000;
    p.constraint_max_iterations=10; p.constraint_max_backtracks=24;
    p.constraint_tolerance=.1; p.constraint_backtrack_factor=.5; p.constraint_sufficient_decrease=1e-4;
    return p;
}

ChemRosenbrock chem_rosenbrock_two_stage(void)
{
    ChemRosenbrock p = defaults();
    double g = 1.0 + 1.0 / sqrt(2.0);
    p.stages=2; p.a[0]=1/g; p.c[0]=-2/g;
    p.m[0]=3/(2*g); p.m[1]=1/(2*g);
    p.e[0]=p.e[1]=1/(2*g); p.gamma[0]=g; p.gamma[1]=-g;
    p.new_function[0]=p.new_function[1]=1; p.order=2;
    return p;
}

ChemRosenbrock chem_rosenbrock_three_stage(void)
{
    ChemRosenbrock p = defaults();
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
    p.gamma[1] = 0.24291996454816804366592249683314;
    p.gamma[2] = 0.21851380027664058511513169485832e+01;
    p.new_function[0] = p.new_function[1] = 1;
    p.order = 3.0;
    return p;
}

ChemRosenbrock chem_rosenbrock_four_stage(void)
{
    ChemRosenbrock p=defaults();
    p.stages=4;
    const double a[]={2.0,1.867943637803922,.2344449711399156,1.867943637803922,.2344449711399156,0};
    const double c[]={-7.137615036412310,2.580708087951457,.6515950076447975,-2.137148994382534,-.3214669691237626,-.6949742501781779};
    const double m[]={2.255570073418735,.2870493262186792,.4353179431840180,1.093502252409163};
    const double e[]={-.2815431932141155,-.07276199124938920,-.1082196201495311,-1.093502252409163};
    memcpy(p.a,a,sizeof a); memcpy(p.c,c,sizeof c); memcpy(p.m,m,sizeof m); memcpy(p.e,e,sizeof e);
    p.new_function[0]=p.new_function[1]=p.new_function[2]=1;
    p.gamma[0]=.57282; p.gamma[1]=-1.769193891319233; p.gamma[2]=.7592633437920482; p.gamma[3]=-.1049021087100450;
    p.order=4;
    return p;
}

ChemRosenbrock chem_rosenbrock_four_stage_dae(void)
{
    ChemRosenbrock p=defaults();
    p.stages=4;
    const double a[]={0,2,0,2,0,1}, c[]={4,1,-1,1,-1,-8.0/3.0};
    const double m[]={2,0,1,1}, e[]={0,0,0,1};
    memcpy(p.a,a,sizeof a); memcpy(p.c,c,sizeof c); memcpy(p.m,m,sizeof m); memcpy(p.e,e,sizeof e);
    p.new_function[0]=p.new_function[2]=p.new_function[3]=1;
    p.gamma[0]=.5; p.gamma[1]=1.5; p.order=3;
    return p;
}

ChemRosenbrock chem_rosenbrock_six_stage_dae(void)
{
    ChemRosenbrock p=defaults();
    p.stages=6;
    const double a[]={1.544,.9466785280815826,.2557011698983284,3.314825187068521,2.896124015972201,
                      .9986419139977817,1.221224509226641,6.019134481288629,12.53708332932087,
                      -.6878860361058950,1.221224509226641,6.019134481288629,12.53708332932087,-.6878860361058950,1};
    const double c[]={-5.6688,-2.430093356833875,-.2063599157091915,-.1073529058151375,-9.594562251023355,
                      -20.47028614809616,7.496443313967647,-10.24680431464352,-33.99990352819905,
                      11.70890893206160,8.083246795921522,-7.981132988064893,-31.52159432874371,
                      16.31930543123136,-6.058818238834054};
    memcpy(p.a,a,sizeof a); memcpy(p.c,c,sizeof c);
    for (size_t i=0;i<4;++i) p.m[i]=a[6+i];
    p.m[4]=p.m[5]=1; p.e[5]=1;
    for (size_t i=0;i<6;++i) p.new_function[i]=1;
    p.gamma[0]=.25; p.gamma[1]=-.1043; p.gamma[2]=.1035; p.gamma[3]=-.03620000000000023;
    p.order=4;
    return p;
}

static void axpy(size_t n, double *restrict dst, const double *restrict src, double a)
{
    #pragma omp simd
    for (size_t i = 0; i < n; ++i) dst[i] += a*src[i];
}

int chem_rosenbrock_normalized_error(size_t species, size_t cells, const double *old_values,
                                     const double *new_values, const double *errors,
                                     const double *atol, double rtol, double *result)
{
    if (!species || !cells || species>SIZE_MAX/cells || !old_values || !new_values || !errors || !atol ||
        !result || !isfinite(rtol) || rtol<0) return CHEM_INVALID;
    double sum = 0;
    for (size_t i = 0; i < species; ++i)
        for (size_t c = 0; c < cells; ++c) {
            size_t at = i*cells+c;
            double scale = atol[i] + rtol * fmax(fabs(old_values[at]), fabs(new_values[at]));
            double ratio = errors[at] / scale;
            sum += ratio*ratio;
        }
    double error = sqrt(sum / (species*cells));
    *result=error < 1e-10 ? 1e-10 : error;
    return CHEM_OK;
}

static double normalized_error(const ChemState *s, const ChemRosenbrock *p)
{
    double result;
    int status=chem_rosenbrock_normalized_error(s->mechanism->species_count,s->cells,s->values,
        s->candidate,s->error,s->absolute_tolerances,p->relative_tolerance,&result);
    return status==CHEM_OK ? result : NAN;
}

/* Preserve both the partial solution and the time reached when integration
 * stops early, so callers can resume after a step-limit return. */
static int finish_solve(ChemState *s, ChemStats *result, size_t ny, double time, int status)
{
    result->final_time = time;
    for (size_t i = 0; i < ny; ++i)
        if (!s->algebraic[i/s->cells] && s->values[i] < 0) s->values[i] = 0;
    return status;
}

int chem_solve(ChemState *s, const ChemRosenbrock *p, double dt, ChemStats *result)
{
    if (!s || !p || !result || !isfinite(dt) || dt < 0 || p->stages < 1 || p->stages > 6 ||
        p->gamma[0] == 0 || p->order <= 0 || p->relative_tolerance < 0 ||
        p->factor_min <= 0 || p->factor_max < p->factor_min || p->safety <= 0 ||
        p->rejection_factor <= 0 || p->h_min < 0 || p->h_max < 0 || p->h_start < 0 ||
        p->constraint_tolerance < 0 || p->constraint_backtrack_factor <= 0 ||
        p->constraint_backtrack_factor >= 1 || p->constraint_sufficient_decrease < 0)
        return CHEM_INVALID;
    memset(result, 0, sizeof(*result));
    size_t ny = s->mechanism->species_count * s->cells;
    if (s->constraint_count || s->external_constraint_count) {
        int status=chem_project_constraints(s,p,result);
        if (status) return status;
    }
    if (dt == 0) {
        for (size_t i = 0; i < ny; ++i) if (!s->algebraic[i/s->cells] && s->values[i] < 0) s->values[i] = 0;
        return CHEM_OK;
    }
    double h_min = p->h_min == 0.0 ? 1e-15*dt : p->h_min;
    double h_max = p->h_max == 0.0 ? dt : fmin(dt, p->h_max);
    double h_start = p->h_start == 0.0 ? 1e-6*dt : fmin(h_max, p->h_start);
    double h = fmin(fmax(h_min, fabs(h_start)), fabs(h_max));
    double time = 0;
    int reject_last = 0, reject_more = 0;
    while (time - dt + DBL_EPSILON <= 0) {
        if (result->internal_steps > p->max_steps)
            return finish_solve(s, result, ny, time, CHEM_TOO_MANY_STEPS);
        if (time + 0.1*h == time || h <= fmin(DBL_EPSILON, h_min))
            return finish_solve(s, result, ny, time, CHEM_STEP_TOO_SMALL);
        h = fmin(h, fabs(dt-time));
        chem_forcing(s, s->values, s->forcing, s->scratch);
        ++result->function_calls;
        chem_negative_jacobian(s, s->values);
        ++result->jacobian_updates;
        int accepted = 0;
        while (!accepted) {
            int status = chem_linear_factor(s, 1.0 / (h*p->gamma[0]));
            ++result->decompositions;
            if (status != CHEM_OK) return finish_solve(s, result, ny, time, status);
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
                for (size_t j = 0; j < stage; ++j) {
                    double factor=p->c[first+j]/h;
                    if (!s->constraint_count && !s->external_constraint_count) axpy(ny,k,s->stages+j*ny,factor);
                    else for (size_t i=0; i<s->mechanism->species_count; ++i) if (!s->algebraic[i])
                        axpy(s->cells,k+i*s->cells,s->stages+j*ny+i*s->cells,factor);
                }
                status = chem_linear_solve(s, k);
                ++result->solves;
                if (status != CHEM_OK) return finish_solve(s, result, ny, time, status);
            }
            memcpy(s->candidate, s->values, ny*sizeof(*s->candidate));
            memset(s->error, 0, ny*sizeof(*s->error));
            for (size_t stage = 0; stage < p->stages; ++stage) {
                axpy(ny, s->candidate, s->stages+stage*ny, p->m[stage]);
                axpy(ny, s->error, s->stages+stage*ny, p->e[stage]);
            }
            for (size_t i=0; i<s->mechanism->species_count; ++i) if (s->algebraic[i])
                for (size_t c=0; c<s->cells; ++c) {
                    size_t at=i*s->cells+c;
                    s->error[at]=s->candidate[at]-s->values[at];
                }
            double error = normalized_error(s, p);
            if (!isfinite(error)) return finish_solve(s, result, ny, time, CHEM_NONFINITE);
            double factor = p->step_controller ? p->step_controller(error,h,p,p->step_controller_context) :
                fmin(p->factor_max, fmax(p->factor_min, p->safety / pow(error, 1.0/p->order)));
            if (!isfinite(factor) || factor <= 0) return finish_solve(s, result, ny, time, CHEM_INVALID);
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
    return finish_solve(s, result, ny, time, CHEM_OK);
}
