/* Backward Euler iteration and reductions ported from MICM/UCAR, Apache-2.0. */
#include "chem.h"
#include <math.h>
#include <string.h>

chem_be_parameters chem_be_parameters_default(void)
{ return (chem_be_parameters){1e-40,0,11,{.5,.5,.5,.5,.1}}; }

int chem_be_converged(const chem_state *s, const chem_be_parameters *p, const double *residual, const double *y)
{
    for (size_t i=0; i<s->model->n_species; ++i) for (size_t c=0; c<s->cells; ++c) {
        size_t k=i*s->cells+c;
        double r=fabs(residual[k]);
        if (!isfinite(r) || (r>p->small && r>s->atol[i] && r>s->rtol*fabs(y[k]))) return 0;
    }
    return 1;
}

chem_result chem_be_solve(chem_state *s, double dt, const chem_be_parameters *p)
{
    chem_result result={CHEM_OK,{0}};
    if (!s || !p || dt<0 || !isfinite(dt)) { result.status=CHEM_INVALID; return result; }
    const chem_model *m=s->model;
    size_t cells=s->cells, n=m->n_species*cells, successful=0, failures=0;
    double *yn=s->work, *f=s->work+n, *y=s->y;
    double h=p->h_start==0 ? dt : p->h_start, time=0;
    if (h<0 || !isfinite(h)) { result.status=CHEM_INVALID; return result; }
    memcpy(yn,y,n*sizeof(double));
    while (time<dt) {
        int converged=0;
        size_t iterations=0;
        do {
            ++result.stats.internal_steps;
            memset(f,0,n*sizeof(double));
            m->forcing(m,cells,s->rates,y,f); ++result.stats.function_calls;
            memset(s->jacobian,0,m->lu.nnz*cells*sizeof(double));
            m->jacobian(m,cells,s->rates,y,s->jacobian); ++result.stats.jacobian_updates;
            for (size_t i=0; i<m->n_species; ++i) {
                double *diagonal=s->jacobian+m->lu.diagonal[i]*cells;
                #pragma omp simd
                for (size_t c=0; c<cells; ++c) diagonal[c]+=1/h;
            }
            result.status=chem_state_factor(s); ++result.stats.decompositions;
            if (result.status!=CHEM_OK) goto done;
            #pragma omp simd
            for (size_t i=0; i<n; ++i) f[i]-=(y[i]-yn[i])/h;
            result.status=chem_state_solve_linear(s,f); ++result.stats.solves;
            if (result.status!=CHEM_OK) goto done;
            #pragma omp simd
            for (size_t i=0; i<n; ++i) y[i]=0>y[i]+f[i] ? 0 : y[i]+f[i];
            if (iterations++==0) continue;
            converged=chem_be_converged(s,p,f,y);
        } while (!converged && iterations<p->max_steps);
        if (!converged) {
            ++result.stats.rejected; successful=0;
            if (failures>=5) {
                time+=h; result.status=CHEM_UNCONVERGED;
                int has_nan=0,has_inf=0;
                for (size_t i=0; i<n; ++i) { has_nan|=isnan(y[i]); has_inf|=isinf(y[i]); }
                if (has_nan) result.status=CHEM_NAN;
                else if (has_inf) result.status=CHEM_INF;
                break;
            }
            memcpy(y,yn,n*sizeof(double)); h*=p->reductions[failures++];
        } else {
            ++result.stats.accepted; time+=h;
            memcpy(yn,y,n*sizeof(double));
            if (++successful>=2) { successful=0; h*=2; }
        }
        if (h>dt-time) h=dt-time;
    }
done:
    result.stats.final_time=time;
    return result;
}
