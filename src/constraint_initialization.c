/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * MICM's weighted Newton projection, with transactional backtracking. */
#include "solver_internal.h"
#include <math.h>
#include <string.h>

static ChemStatus inspect(const ChemState *state,const double *v,double *maximum)
{
    int inf=0,nan=0;
    *maximum=0;
    for (size_t a=0;a<state->mechanism->variable_count;++a) if (state->mechanism->algebraic[a])
        for (size_t c=0;c<state->cells;++c) {
            double value=v[a*state->cells+c];
            nan|=isnan(value);inf|=isinf(value);
            if (isfinite(value) && fabs(value)>*maximum) *maximum=fabs(value);
        }
    return nan?CHEM_NAN:inf?CHEM_INF:CHEM_OK;
}
static double weighted_correction(const ChemState *state,const double *y,const double *delta)
{
    double maximum=0;
    for (size_t a=0;a<state->mechanism->variable_count;++a) if (state->mechanism->algebraic[a])
        for (size_t c=0;c<state->cells;++c) {
            size_t i=a*state->cells+c;
            double corrected=y[i]+delta[i];
            double zmax=fabs(y[i])>fabs(corrected)?fabs(y[i]):fabs(corrected);
            double scale=state->absolute_tolerance[a]+state->relative_tolerance*zmax;
            double q=scale>0?fabs(delta[i])/scale:delta[i]==0?0:INFINITY;
            maximum=maximum>q?maximum:q;
        }
    return maximum;
}
static void update(const ChemState *state,double *y,const double *delta,double step)
{
    for (size_t a=0;a<state->mechanism->variable_count;++a) if (state->mechanism->algebraic[a])
        for (size_t c=0;c<state->cells;++c) {
            size_t i=a*state->cells+c;
            double d=step*delta[i];
            if (isfinite(d)) y[i]+=d;
        }
}
ChemStatus chem_initialize_constraints(ChemSolver *s,ChemState *state,
                                       const ChemRosenbrockParameters *p,ChemSolverStats *stats)
{
    if (!s || !state || !p || !stats || state->mechanism!=s->mechanism || state->cells!=s->cells)
        return CHEM_INVALID;
    if (!s->mechanism->constraint_count) return CHEM_OK;
    double *delta=s->initial,*original=s->stages,*candidate=s->next,*trial=s->error;
    size_t bytes=s->count*sizeof(double);
    int modified=0;
    ChemStatus status=CHEM_CONSTRAINT_FAILURE;
    for (size_t iteration=0;iteration<=p->constraint_max_iterations;++iteration) {
        memset(delta,0,bytes);
        chem_constraints_forcing(state,state->values,delta);
        double residual;
        status=inspect(state,delta,&residual);
        if (status) goto restore;
        ++stats->constraint_init_iterations;
        if (residual==0) return CHEM_OK;
        memset(s->matrix,0,chem_sparse_elements(s->sparse)*s->cells*sizeof(double));
        chem_constraints_jacobian_add(state,state->values,s->matrix,s->sparse,-1);
        ++stats->jacobian_updates;
        chem_add_mass_diagonal(s->mechanism,s->sparse,s->cells,1,s->matrix);
        status=s->linear.factor(s->sparse,s->cells,s->matrix,s->linear.context);
        ++stats->decompositions;
        if (status) goto restore;
        s->linear.solve(s->sparse,s->cells,s->matrix,delta,s->linear.context);
        ++stats->solves;
        status=inspect(state,delta,&residual);
        if (status) goto restore;
        double norm=weighted_correction(state,state->values,delta);
        if (norm<=p->constraint_tolerance) return CHEM_OK;
        if (iteration==p->constraint_max_iterations) break;
        if (!modified) { memcpy(original,state->values,bytes);modified=1; }
        if (p->constraint_max_backtracks) {
            double step=1;
            for (size_t backtrack=0;;++backtrack,step*=p->constraint_backtrack_factor) {
                if (backtrack>p->constraint_max_backtracks) {
                    status=CHEM_CONSTRAINT_FAILURE;goto restore;
                }
                memcpy(candidate,state->values,bytes);
                update(state,candidate,delta,step);
                if (inspect(state,candidate,&residual)) continue;
                memset(trial,0,bytes);
                chem_constraints_forcing(state,candidate,trial);
                if (inspect(state,trial,&residual)) continue;
                if (residual==0) break;
                s->linear.solve(s->sparse,s->cells,s->matrix,trial,s->linear.context);
                ++stats->solves;
                if (inspect(state,trial,&residual)) continue;
                double q=weighted_correction(state,candidate,trial);
                if (q<=p->constraint_tolerance || q<(1-p->constraint_sufficient_decrease*step)*norm) break;
            }
            memcpy(state->values,candidate,bytes);
        } else update(state,state->values,delta,1);
    }
    status=CHEM_CONSTRAINT_FAILURE;
restore:
    if (modified) memcpy(state->values,original,bytes);
    return status;
}
