/* Rosenbrock integration and step control ported from MICM/UCAR, Apache-2.0. */
#include "chem.h"
#include <math.h>
#include <string.h>

static double minimum(double a,double b) { return b<a ? b : a; }
static double maximum(double a,double b) { return a<b ? b : a; }
static void axpy(size_t n, double a, const double *restrict x, double *restrict y)
{
    #pragma omp simd
    for (size_t i=0; i<n; ++i) y[i]+=a*x[i];
}

double chem_normalized_error(const chem_state *s, const double *y, const double *ynew, const double *error)
{
    double sum=0;
    for (size_t i=0; i<s->model->n_species; ++i) {
        for (size_t c=0; c<s->cells; ++c) {
            size_t k=i*s->cells+c;
            double ymax=maximum(fabs(y[k]),fabs(ynew[k]));
            double e=error[k]/(s->atol[i]+s->rtol*ymax);
            sum+=e*e;
        }
    }
    return maximum(sqrt(sum/(double)(s->cells*s->model->n_species)),1.0e-10);
}

double chem_ros_controller(double h, double error, const chem_ros_parameters *p, void *context)
{
    (void)context;
    return h*minimum(p->factor_max, maximum(p->factor_min,p->safety_factor/pow(error,1/p->order)));
}

chem_result chem_ros_solve(chem_state *s, double dt, const chem_ros_parameters *p,
                           chem_step_controller controller, void *context)
{
    chem_result result={CHEM_OK,{0}};
    if (!s || !p || dt<0 || !isfinite(dt) || p->stages<1 || p->stages>6) { result.status=CHEM_INVALID; return result; }
    const chem_model *m=s->model;
    size_t cells=s->cells, n=m->n_species*cells;
    double *y=s->y, *ynew=s->work, *initial=s->work+n, *error=s->work+2*n;
    double *k[6];
    for (size_t i=0; i<6; ++i) k[i]=s->work+(3+i)*n;
    double hmin=p->h_min==0 ? 1e-15*dt : p->h_min;
    double hmax=p->h_max==0 ? dt : minimum(dt,p->h_max);
    double hstart=p->h_start==0 ? 1e-6*dt : minimum(hmax,p->h_start);
    double h=minimum(maximum(hmin,fabs(hstart)),fabs(hmax)), time=0;
    int reject_last=0, reject_more=0;
    if (!controller) controller=chem_ros_controller;
    while ((time-dt+p->round_off)<=0 && result.status==CHEM_OK) {
        if (result.stats.internal_steps>p->max_steps) { result.status=CHEM_MAX_STEPS; break; }
        if ((time+0.1*h)==time || h<=minimum(p->round_off,hmin)) { result.status=CHEM_SMALL_STEP; break; }
        h=minimum(h,fabs(dt-time));
        memset(initial,0,n*sizeof(double));
        m->forcing(m,cells,s->rates,y,initial); ++result.stats.function_calls;
        memset(s->jacobian,0,m->lu.nnz*cells*sizeof(double));
        m->jacobian(m,cells,s->rates,y,s->jacobian); ++result.stats.jacobian_updates;
        int accepted=0;
        while (!accepted) {
            double alpha=1.0/(h*p->gamma[0]);
            for (size_t i=0; i<m->n_species; ++i) {
                double *diag=s->jacobian+m->lu.diagonal[i]*cells;
                #pragma omp simd
                for (size_t c=0; c<cells; ++c) diag[c]+=alpha;
            }
            chem_lu_factor(&m->lu,cells,s->jacobian,s->scratch); ++result.stats.decompositions;
            for (size_t stage=0; stage<p->stages; ++stage) {
                size_t base=stage ? stage*(stage-1)/2 : 0;
                if (stage==0) memcpy(k[stage],initial,n*sizeof(double));
                else if (p->new_function[stage]) {
                    memcpy(ynew,y,n*sizeof(double));
                    for (size_t j=0; j<stage; ++j) axpy(n,p->a[base+j],k[j],ynew);
                    memset(k[stage],0,n*sizeof(double));
                    m->forcing(m,cells,s->rates,ynew,k[stage]); ++result.stats.function_calls;
                }
                if (stage+1<p->stages && !p->new_function[stage+1]) memcpy(k[stage+1],k[stage],n*sizeof(double));
                for (size_t j=0; j<stage; ++j) axpy(n,p->c[base+j]/h,k[j],k[stage]);
                chem_lu_solve(&m->lu,cells,s->jacobian,k[stage]); ++result.stats.solves;
            }
            memcpy(ynew,y,n*sizeof(double));
            for (size_t stage=0; stage<p->stages; ++stage) axpy(n,p->m[stage],k[stage],ynew);
            memset(error,0,n*sizeof(double));
            for (size_t stage=0; stage<p->stages; ++stage) axpy(n,p->e[stage],k[stage],error);
            double norm=chem_normalized_error(s,y,ynew,error);
            double next=controller(h,norm,p,context);
            ++result.stats.internal_steps;
            if (isnan(norm)) { result.status=CHEM_NAN; break; }
            if (isinf(norm)) { result.status=CHEM_INF; break; }
            if (norm<1 || h<hmin) {
                ++result.stats.accepted; time+=h;
                memcpy(y,ynew,n*sizeof(double));
                next=maximum(hmin,minimum(next,hmax));
                if (reject_last) next=minimum(next,h);
                reject_last=0; reject_more=0; h=next; accepted=1;
            } else {
                if (reject_more) next=h*p->rejection_factor;
                reject_more=reject_last; reject_last=1; h=next;
                if (result.stats.accepted>=1) ++result.stats.rejected;
                memset(s->jacobian,0,m->lu.nnz*cells*sizeof(double));
                m->jacobian(m,cells,s->rates,y,s->jacobian); ++result.stats.jacobian_updates;
            }
        }
    }
    result.stats.final_time=time;
    return result;
}
