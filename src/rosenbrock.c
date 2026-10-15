/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Rosenbrock formulas and controller translated from MICM. */
#include "chem_solver.h"
#include "solver_internal.h"
#include <math.h>
#include <string.h>

static ChemStatus default_factor(const ChemSparse *a,size_t n,double *v,void *ctx)
{ (void)ctx; return chem_lu_factor(a,n,v,CHEM_LU_MOZART); }
static void default_solve(const ChemSparse *a,size_t n,const double *v,double *b,void *ctx)
{ (void)ctx; chem_lu_solve(a,n,v,b); }
void chem_solver_set_controller(ChemSolver *s,ChemStepController fn,void *ctx)
{ s->controller=fn?fn:chem_default_step_controller; s->controller_context=ctx; }
void chem_solver_set_linear_solver(ChemSolver *s,ChemLinearSolver linear)
{ if (linear.factor && linear.solve) s->linear=linear; }
void chem_solver_free(ChemSolver *s)
{
    if (!s) return;
    if (s->derivatives) for (size_t i=0;i<s->derivative_count;++i) {
        free(s->derivatives[i].dependent); free(s->derivatives[i].updates); free(s->derivatives[i].yields);
    }
    free(s->derivatives); chem_sparse_free(s->sparse); free(s->matrix);
    free(s->next); free(s->initial); free(s->error); free(s->stages); free(s->work); free(s);
}
ChemStatus chem_solver_create(const ChemMechanism *m,size_t cells,ChemSolver **out)
{
    if (!out || !m) return CHEM_INVALID;
    *out=NULL;
    ChemSolver *s=calloc(1,sizeof(*s));
    if (!s) return CHEM_MEMORY;
    s->mechanism=m; s->cells=cells;
    ChemStatus status=CHEM_MEMORY;
    unsigned char *pattern=chem_mechanism_pattern(m);
    if (!pattern) goto fail;
    status=chem_sparse_create(m->variable_count,pattern,&s->sparse);
    free(pattern);
    if (status) goto fail;
    status=CHEM_MEMORY;
    if ((m->variable_count && cells>SIZE_MAX/m->variable_count) ||
        (chem_sparse_elements(s->sparse) && cells>SIZE_MAX/chem_sparse_elements(s->sparse))) goto fail;
    s->count=cells*m->variable_count;
    s->matrix=chem_alloc(cells*chem_sparse_elements(s->sparse),sizeof(double));
    s->next=chem_alloc(s->count,sizeof(double)); s->initial=chem_alloc(s->count,sizeof(double));
    s->error=chem_alloc(s->count,sizeof(double)); s->stages=chem_alloc(s->count,6*sizeof(double));
    s->work=chem_alloc(cells,sizeof(double));
    if (!s->matrix || !s->next || !s->initial || !s->error || !s->stages || !s->work) goto fail;
    for (size_t r=0;r<m->reaction_count;++r) s->derivative_count+=m->processes[r].reactant_count;
    s->derivatives=chem_alloc(s->derivative_count,sizeof(Derivative));
    if (!s->derivatives) goto fail;
    size_t at=0;
    /* Column then reaction order matches the source Jacobian accumulation. */
    for (size_t col=0;col<m->variable_count;++col) for (size_t r=0;r<m->reaction_count;++r) {
        const ChemProcess *p=&m->processes[r];
        for (size_t occurrence=0;occurrence<p->reactant_count;++occurrence) {
            if (p->reactants[occurrence]!=col) continue;
            Derivative *d=&s->derivatives[at++];
            d->reaction=r; d->independent=col; d->dependent_count=p->reactant_count-1; d->product_count=p->product_count;
            d->dependent=chem_alloc(d->dependent_count,sizeof(size_t));
            d->updates=chem_alloc(p->reactant_count+p->product_count,sizeof(size_t));
            d->yields=chem_alloc(p->product_count,sizeof(double));
            if (!d->dependent || !d->updates || !d->yields) goto fail;
            size_t k=0; int skipped=0;
            for (size_t j=0;j<p->reactant_count;++j) {
                if (!skipped && p->reactants[j]==col) { skipped=1; continue; }
                d->dependent[k]=p->reactants[j];
                d->updates[k++]=m->algebraic[p->reactants[j]]?SIZE_MAX:chem_sparse_index(s->sparse,p->reactants[j],col);
            }
            d->updates[k++]=m->algebraic[col]?SIZE_MAX:chem_sparse_index(s->sparse,col,col);
            for (size_t j=0;j<p->product_count;++j) {
                d->updates[k++]=m->algebraic[p->products[j].species]?SIZE_MAX:chem_sparse_index(s->sparse,p->products[j].species,col);
                d->yields[j]=p->products[j].yield;
            }
        }
    }
    s->controller=chem_default_step_controller;
    s->linear=(ChemLinearSolver){default_factor,default_solve,NULL};
    *out=s; return CHEM_OK;
fail:
    chem_solver_free(s); return status;
}
void chem_solver_negative_jacobian(ChemSolver *s,const ChemState *state,const double *y)
{
    size_t cells=s->cells;
    double *restrict work=s->work;
    memset(s->matrix,0,chem_sparse_elements(s->sparse)*cells*sizeof(double));
    for (size_t i=0;i<s->derivative_count;++i) {
        const Derivative *d=&s->derivatives[i];
        memcpy(work,state->rates+d->reaction*cells,cells*sizeof(double));
        for (size_t k=0;k<d->dependent_count;++k) {
            const double *v=y+d->dependent[k]*cells;
            #pragma omp simd
            for (size_t c=0;c<cells;++c) work[c]*=v[c];
        }
        for (size_t k=0;k<=d->dependent_count;++k) {
            if (d->updates[k]==SIZE_MAX) continue;
            double *restrict dst=s->matrix+d->updates[k]*cells;
            #pragma omp simd
            for (size_t c=0;c<cells;++c) dst[c]+=work[c];
        }
        for (size_t k=0;k<d->product_count;++k) {
            if (d->updates[d->dependent_count+1+k]==SIZE_MAX) continue;
            double *restrict dst=s->matrix+d->updates[d->dependent_count+1+k]*cells;
            double yield=d->yields[k];
            #pragma omp simd
            for (size_t c=0;c<cells;++c) dst[c]-=yield*work[c];
        }
    }
    chem_constraints_jacobian_add(state,y,s->matrix,s->sparse,-1);
}
void chem_add_mass_diagonal(const ChemMechanism *m,const ChemSparse *pattern,size_t cells,double alpha,double *matrix)
{
    for (size_t v=0;v<m->variable_count;++v) if (!m->algebraic[v]) {
        double *diag=matrix+chem_sparse_index(pattern,v,v)*cells;
        #pragma omp simd
        for (size_t c=0;c<cells;++c) diag[c]+=alpha;
    }
}
static void axpy(size_t n,double *restrict y,double a,const double *restrict x)
{
    #pragma omp simd
    for (size_t i=0;i<n;++i) y[i]+=a*x[i];
}
double chem_normalized_error(const ChemState *s,const double *y,const double *next,const double *error)
{
    double sum=0;
    size_t n=s->mechanism->variable_count, cells=s->cells;
    for (size_t v=0;v<n;++v) for (size_t c=0;c<cells;++c) {
        size_t i=v*cells+c;
        double max=fabs(y[i])>fabs(next[i])?fabs(y[i]):fabs(next[i]);
        double e=error[i]/(s->absolute_tolerance[v]+s->relative_tolerance*max);
        sum+=e*e;
    }
    double result=sqrt(sum/((n && cells)?n*cells:1));
    return result<1e-10?1e-10:result;
}
double chem_default_step_controller(double error,const ChemRosenbrockParameters *p,void *ctx)
{
    (void)ctx;
    return fmin(p->factor_max,fmax(p->factor_min,p->safety_factor/pow(error,1/p->local_order)));
}
ChemStatus chem_solve(ChemSolver *s,ChemState *state,double dt,const ChemRosenbrockParameters *p,ChemSolverStats *stats)
{
    if (!s || !state || !p || !stats || dt<0 || !isfinite(dt) || p->stages<1 || p->stages>6 ||
        s->mechanism!=state->mechanism || s->cells!=state->cells) return CHEM_INVALID;
    memset(stats,0,sizeof(*stats));
    double hmin=p->h_min==0?1e-15*dt:p->h_min;
    double hmax=p->h_max==0?dt:fmin(dt,p->h_max);
    double start=p->h_start==0?1e-6*dt:fmin(hmax,p->h_start);
    double h=fmin(fmax(hmin,fabs(start)),fabs(hmax));
    double time=0;
    int reject_last=0,reject_more=0;
    size_t n=s->count;
    ChemStatus status=chem_initialize_constraints(s,state,p,stats);
    if (status) return status;
    while (time-dt+p->round_off<=0) {
        if (stats->internal_steps>p->max_steps) { status=CHEM_MAX_STEPS; break; }
        if (time+.1*h==time || h<=fmin(p->round_off,hmin)) { status=CHEM_STEP_TOO_SMALL; break; }
        h=fmin(h,fabs(dt-time));
        chem_forcing(state,state->values,s->initial); ++stats->function_calls;
        chem_solver_negative_jacobian(s,state,state->values); ++stats->jacobian_updates;
        int accepted=0;
        while (!accepted) {
            double alpha=1/(h*p->gamma[0]);
            chem_add_mass_diagonal(s->mechanism,s->sparse,s->cells,alpha,s->matrix);
            status=s->linear.factor(s->sparse,s->cells,s->matrix,s->linear.context);
            ++stats->decompositions;
            if (status) goto finish;
            for (size_t stage=0;stage<p->stages;++stage) {
                size_t offset=stage*(stage-1)/2;
                double *k=s->stages+stage*n;
                if (stage==0) memcpy(k,s->initial,n*sizeof(double));
                else if (p->new_function[stage]) {
                    memcpy(s->next,state->values,n*sizeof(double));
                    for (size_t j=0;j<stage;++j) axpy(n,s->next,p->a[offset+j],s->stages+j*n);
                    chem_forcing(state,s->next,k); ++stats->function_calls;
                }
                if (stage+1<p->stages && !p->new_function[stage+1]) memcpy(k+n,k,n*sizeof(double));
                for (size_t j=0;j<stage;++j) {
                    if (!s->mechanism->constraint_count) axpy(n,k,p->c[offset+j]/h,s->stages+j*n);
                    else for (size_t v=0;v<s->mechanism->variable_count;++v)
                        if (!s->mechanism->algebraic[v])
                            axpy(s->cells,k+v*s->cells,p->c[offset+j]/h,s->stages+j*n+v*s->cells);
                }
                s->linear.solve(s->sparse,s->cells,s->matrix,k,s->linear.context); ++stats->solves;
            }
            memcpy(s->next,state->values,n*sizeof(double)); memset(s->error,0,n*sizeof(double));
            for (size_t stage=0;stage<p->stages;++stage) axpy(n,s->next,p->m[stage],s->stages+stage*n);
            for (size_t stage=0;stage<p->stages;++stage) axpy(n,s->error,p->e[stage],s->stages+stage*n);
            for (size_t v=0;v<s->mechanism->variable_count;++v) if (s->mechanism->algebraic[v])
                for (size_t c=0;c<s->cells;++c) {
                    size_t i=v*s->cells+c;s->error[i]=s->next[i]-state->values[i];
                }
            double error=chem_normalized_error(state,state->values,s->next,s->error);
            double hnew=h*s->controller(error,p,s->controller_context);
            ++stats->internal_steps;
            if (!isfinite(error)) { status=isnan(error)?CHEM_NAN:CHEM_INF; goto finish; }
            if (error<1 || h<hmin) {
                ++stats->accepted; time+=h;
                double *swap=state->values; state->values=s->next; s->next=swap;
                hnew=fmax(hmin,fmin(hnew,hmax));
                if (reject_last) hnew=fmin(hnew,h);
                reject_last=reject_more=0; h=hnew; accepted=1;
            } else {
                if (reject_more) hnew=h*p->rejection_factor;
                reject_more=reject_last; reject_last=1; h=hnew;
                if (stats->accepted>=1) ++stats->rejected;
                chem_solver_negative_jacobian(s,state,state->values); ++stats->jacobian_updates;
            }
        }
    }
finish:
    /* Public MICM Solve clamps concentrations after each requested interval,
     * not after individual internal stages or accepted substeps. */
    for (size_t v=0;v<s->mechanism->variable_count;++v) if (!s->mechanism->algebraic[v]) {
        double *values=state->values+v*s->cells;
        #pragma omp simd
        for (size_t c=0;c<s->cells;++c) if (values[c]<0) values[c]=0;
    }
    stats->final_time=time; return status;
}
