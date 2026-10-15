/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Translation of the MICM backward Euler Newton iteration and step reductions. */
#include "solver_internal.h"
#include <math.h>
#include <string.h>

ChemBackwardEulerParameters chem_backward_euler_parameters(void)
{
    return (ChemBackwardEulerParameters){1e-40,0,11,{.5,.5,.5,.5,.1}};
}
int chem_backward_euler_converged(const ChemState *s,const double *residual,
                                 const double *values,const ChemBackwardEulerParameters *p)
{
    for (size_t v=0;v<s->mechanism->variable_count;++v) for (size_t c=0;c<s->cells;++c) {
        size_t i=v*s->cells+c;
        double r=fabs(residual[i]);
        if (!isfinite(r) || (r>p->small && r>s->absolute_tolerance[v] && r>s->relative_tolerance*fabs(values[i]))) return 0;
    }
    return 1;
}
ChemStatus chem_solve_backward_euler(ChemSolver *s,ChemState *state,double dt,
                                    const ChemBackwardEulerParameters *p,ChemSolverStats *stats)
{
    if (!s || !state || !p || !stats || dt<0 || !isfinite(dt) ||
        s->mechanism!=state->mechanism || s->cells!=state->cells || p->h_start<0) return CHEM_INVALID;
    memset(stats,0,sizeof(*stats));
    double h=p->h_start==0?dt:p->h_start,time=0;
    size_t successful=0,failures=0,n=s->count;
    double *previous=s->next,*forcing=s->initial;
    memcpy(previous,state->values,n*sizeof(double));
    ChemStatus status=CHEM_OK;
    while (time<dt) {
        int converged=0;
        size_t iterations=0;
        do {
            ++stats->internal_steps;
            chem_forcing(state,state->values,forcing);++stats->function_calls;
            chem_solver_negative_jacobian(s,state,state->values);++stats->jacobian_updates;
            for (size_t v=0;v<s->mechanism->variable_count;++v) {
                double *diag=s->matrix+chem_sparse_index(s->sparse,v,v)*s->cells;
                #pragma omp simd
                for (size_t c=0;c<s->cells;++c) diag[c]+=1/h;
            }
            status=s->linear.factor(s->sparse,s->cells,s->matrix,s->linear.context);++stats->decompositions;
            if (status) goto finish;
            #pragma omp simd
            for (size_t i=0;i<n;++i) forcing[i]-=(state->values[i]-previous[i])/h;
            s->linear.solve(s->sparse,s->cells,s->matrix,forcing,s->linear.context);++stats->solves;
            #pragma omp simd
            for (size_t i=0;i<n;++i) {
                double next=state->values[i]+forcing[i];
                state->values[i]=0>next?0:next;
            }
            if (iterations++==0) continue;
            converged=chem_backward_euler_converged(state,forcing,state->values,p);
        } while (!converged && iterations<p->max_iterations);
        if (!converged) {
            ++stats->rejected;successful=0;
            if (failures>=5) {
                time+=h;status=CHEM_UNCONVERGED;
                int nan=0,inf=0;
                for (size_t i=0;i<n;++i) { nan|=isnan(state->values[i]);inf|=isinf(state->values[i]); }
                if (nan) status=CHEM_NAN;
                else if (inf) status=CHEM_INF;
                break;
            }
            memcpy(state->values,previous,n*sizeof(double));
            h*=p->time_step_reductions[failures++];
        } else {
            ++stats->accepted;time+=h;memcpy(previous,state->values,n*sizeof(double));
            if (++successful>=2) { successful=0;h*=2; }
        }
        h=fmin(h,dt-time);
    }
finish:
    #pragma omp simd
    for (size_t i=0;i<n;++i) if (state->values[i]<0) state->values[i]=0;
    stats->final_time=time;
    return status;
}
