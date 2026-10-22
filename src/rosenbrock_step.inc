/* Shared MICM stage equations and control; instantiated by generated kernels.
 * SPDX-License-Identifier: Apache-2.0 */
#ifndef CHEM_STAGE_COUNT
#define CHEM_STAGE_COUNT p->stages
#endif
#ifndef CHEM_ALGEBRAIC
#define CHEM_ALGEBRAIC(s,v) ((s)->mechanism->algebraic[v])
#endif
#ifndef CHEM_INITIAL
#define CHEM_INITIAL(s,state) do { \
    chem_forcing(state,(state)->values,(s)->initial); \
    chem_solver_negative_jacobian(s,state,(state)->values); \
} while (0)
#endif
static void axpy(size_t n,double *restrict y,double a,const double *restrict x)
{
    #pragma omp simd
    for (size_t i=0;i<n;++i) y[i]+=a*x[i];
}
ChemStatus CHEM_SOLVE(ChemSolver *s,ChemState *state,double dt,const ChemRosenbrockParameters *p,ChemSolverStats *stats)
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
        CHEM_INITIAL(s,state);
        ++stats->function_calls; ++stats->jacobian_updates;
        int accepted=0;
        while (!accepted) {
            double alpha=1/(h*p->gamma[0]);
            chem_add_mass_diagonal(s->mechanism,s->sparse,s->cells,alpha,s->matrix);
            status=s->linear.factor(s->sparse,s->cells,s->matrix,s->linear.context);
            ++stats->decompositions;
            if (status) goto finish;
            for (size_t stage=0;stage<CHEM_STAGE_COUNT;++stage) {
                size_t offset=stage*(stage-1)/2;
                double *k=s->stages+stage*n;
                if (stage==0) memcpy(k,s->initial,n*sizeof(double));
                else if (p->new_function[stage]) {
                    memcpy(s->next,state->values,n*sizeof(double));
                    for (size_t j=0;j<stage;++j) axpy(n,s->next,p->a[offset+j],s->stages+j*n);
                    chem_forcing(state,s->next,k); ++stats->function_calls;
                }
                if (stage+1<CHEM_STAGE_COUNT && !p->new_function[stage+1]) memcpy(k+n,k,n*sizeof(double));
                for (size_t j=0;j<stage;++j) {
                    if (!chem_constraint_count(s->mechanism)) axpy(n,k,p->c[offset+j]/h,s->stages+j*n);
                    else for (size_t v=0;v<s->mechanism->variable_count;++v)
                        if (!CHEM_ALGEBRAIC(s,v))
                            axpy(s->cells,k+v*s->cells,p->c[offset+j]/h,s->stages+j*n+v*s->cells);
                }
                s->linear.solve(s->sparse,s->cells,s->matrix,k,s->linear.context); ++stats->solves;
            }
#ifdef CHEM_FUSED_ESTIMATE
            #pragma omp simd
            for (size_t i=0;i<n;++i) {
                double next=state->values[i], error=0;
                for (size_t stage=0;stage<CHEM_STAGE_COUNT;++stage) next+=p->m[stage]*s->stages[stage*n+i];
                for (size_t stage=0;stage<CHEM_STAGE_COUNT;++stage) error+=p->e[stage]*s->stages[stage*n+i];
                s->next[i]=next; s->error[i]=error;
            }
#else
            memcpy(s->next,state->values,n*sizeof(double)); memset(s->error,0,n*sizeof(double));
            for (size_t stage=0;stage<CHEM_STAGE_COUNT;++stage) axpy(n,s->next,p->m[stage],s->stages+stage*n);
            for (size_t stage=0;stage<CHEM_STAGE_COUNT;++stage) axpy(n,s->error,p->e[stage],s->stages+stage*n);
#endif
            for (size_t v=0;v<s->mechanism->variable_count;++v) if (CHEM_ALGEBRAIC(s,v))
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
    for (size_t v=0;v<s->mechanism->variable_count;++v) if (!CHEM_ALGEBRAIC(s,v)) {
        double *values=state->values+v*s->cells;
        #pragma omp simd
        for (size_t c=0;c<s->cells;++c) if (values[c]<0) values[c]=0;
    }
    stats->final_time=time; return status;
}
