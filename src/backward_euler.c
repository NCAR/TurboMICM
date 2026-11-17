#include "internal.h"

#include <math.h>
#include <string.h>

ChemBackwardEuler chem_backward_euler_defaults(void)
{
    ChemBackwardEuler p = {1e-40, 0, 1e-6, 11, {.5, .5, .5, .5, .1}};
    return p;
}

int chem_backward_euler_converged(size_t n, size_t cells, const double *residual,
                                  const double *y, const double *atol, double rtol, double small)
{
    if (!n || !cells || !residual || !y || !atol || rtol < 0 || small < 0) return 0;
    for (size_t i = 0; i < n; ++i) for (size_t c = 0; c < cells; ++c) {
        size_t k = i*cells+c;
        double d = fabs(residual[k]);
        if (!isfinite(residual[k]) || (d > small && d > atol[i] && d > rtol*fabs(y[k]))) return 0;
    }
    return 1;
}

int chem_solve_backward_euler(ChemState *s, const ChemBackwardEuler *p, double dt, ChemStats *stats)
{
    if (!s || !p || !stats || s->constraint_count || !isfinite(dt) || dt < 0 || !isfinite(p->h_start) || p->h_start < 0 ||
        p->relative_tolerance < 0 || p->small < 0 || !p->max_iterations) return CHEM_INVALID;
    for (size_t i = 0; i < 5; ++i)
        if (p->time_step_reductions[i] <= 0 || p->time_step_reductions[i] >= 1) return CHEM_INVALID;
    memset(stats, 0, sizeof(*stats));
    if (dt == 0) return CHEM_OK;
    size_t n=s->mechanism->species_count, cells=s->cells, ny=n*cells;
    double h=p->h_start == 0 ? dt : p->h_start, time=0;
    unsigned successes=0, failures=0;
    memcpy(s->candidate,s->values,ny*sizeof(double));
    while (time < dt) {
        unsigned iterations=0;
        int converged=0;
        do {
            ++stats->internal_steps;
            chem_forcing(s,s->values,s->forcing,s->scratch); ++stats->function_calls;
            chem_negative_jacobian(s,s->values); ++stats->jacobian_updates;
            int status=chem_linear_factor(s,1/h); ++stats->decompositions;
            if (status) { stats->final_time=time; return status; }
            #pragma omp simd
            for (size_t i=0;i<ny;++i) s->forcing[i]-=(s->values[i]-s->candidate[i])/h;
            status=chem_linear_solve(s,s->forcing); ++stats->solves;
            if (status) { stats->final_time=time; return status; }
            #pragma omp simd
            for (size_t i=0;i<ny;++i) {
                double next=s->values[i]+s->forcing[i];
                s->values[i]=0.0 > next ? 0.0 : next;
            }
            ++iterations;
            if (iterations > 1)
                converged=chem_backward_euler_converged(n,cells,s->forcing,s->values,s->absolute_tolerances,
                                                        p->relative_tolerance,p->small);
        } while (!converged && iterations < p->max_iterations);
        if (!converged) {
            ++stats->rejected; successes=0;
            if (failures >= 5) {
                stats->final_time=time+h;
                for (size_t i=0;i<ny;++i) if (!isfinite(s->values[i])) return CHEM_NONFINITE;
                return CHEM_UNCONVERGED;
            }
            memcpy(s->values,s->candidate,ny*sizeof(double));
            h*=p->time_step_reductions[failures++];
        } else {
            ++stats->accepted;
            time+=h;
            memcpy(s->candidate,s->values,ny*sizeof(double));
            if (++successes >= 2) { successes=0; h*=2; }
        }
        h=fmin(h,dt-time);
        if (time < dt && (h <= 0 || time+h == time)) { stats->final_time=time; return CHEM_STEP_TOO_SMALL; }
    }
    stats->final_time=time;
    return CHEM_OK;
}
