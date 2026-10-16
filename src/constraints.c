/* Constraint formulas and projection ported from MICM/UCAR, Apache-2.0. */
#include "chem.h"
#include <math.h>
#include <string.h>

chem_status chem_constraint_validate(const chem_constraint *c, size_t n)
{
    if (!c || !c->name || c->algebraic_species>=n) return CHEM_NOT_FOUND;
    if (c->type==CHEM_CUSTOM_CONSTRAINT) {
        if (!c->residual || !c->jacobian || (c->n_dependencies && !c->dependencies)) return CHEM_INVALID;
        for (size_t i=0; i<c->n_dependencies; ++i) if (c->dependencies[i]>=n) return CHEM_NOT_FOUND;
        return CHEM_OK;
    }
    if (!c->n_reactants || !c->reactants) return CHEM_INVALID;
    if (c->type==CHEM_EQUILIBRIUM_CONSTRAINT && (!c->n_products || !c->products || c->equilibrium.K<=0)) return CHEM_INVALID;
    for (size_t i=0; i<c->n_reactants; ++i) {
        if (c->reactants[i].species>=n) return CHEM_NOT_FOUND;
        if (c->type==CHEM_EQUILIBRIUM_CONSTRAINT && c->reactants[i].yield<=0) return CHEM_INVALID;
    }
    for (size_t i=0; i<c->n_products; ++i) {
        if (c->products[i].species>=n) return CHEM_NOT_FOUND;
        if (c->type==CHEM_EQUILIBRIUM_CONSTRAINT && c->products[i].yield<=0) return CHEM_INVALID;
    }
    return CHEM_OK;
}

void chem_update_constraint_parameters(chem_state *s)
{
    const chem_model *m=s->model;
    for (size_t k=0; k<m->mechanism->n_constraints; ++k) {
        const chem_constraint *constraint=&m->mechanism->constraints[k];
        if (constraint->type!=CHEM_EQUILIBRIUM_CONSTRAINT) continue;
        const chem_vant_hoff *p=&constraint->equilibrium;
        double R=p->R==0 ? CHEM_GAS_CONSTANT : p->R;
        double ref=p->reference_temperature==0 ? 298.15 : p->reference_temperature;
        size_t param=m->constraint_offsets[k];
        for (size_t c=0; c<s->cells; ++c)
            s->custom[param*s->cells+c]=p->K*exp((p->delta_H/R)*(1.0/s->conditions[c].temperature-1.0/ref));
    }
}

static double term_product(const chem_model *m, size_t cells, size_t cell, const double *y,
                           size_t count, const chem_product *terms, size_t skip, double initial)
{
    for (size_t j=0; j<count; ++j) if (j!=skip) {
        double conc=y[m->spec_to_var[terms[j].species]*cells+cell];
        initial*=pow(conc>0 ? conc : 0,terms[j].yield);
    }
    return initial;
}

void chem_constraint_forcing(const chem_state *s, const double *y, double *f)
{
    const chem_model *m=s->model;
    for (size_t k=0; k<m->mechanism->n_constraints; ++k) {
        const chem_constraint *con=&m->mechanism->constraints[k];
        if (con->type==CHEM_CUSTOM_CONSTRAINT) { con->residual(s,y,f,con->context); continue; }
        size_t row=m->constraint_rows[k];
        for (size_t c=0; c<s->cells; ++c) {
            double value=0;
            if (con->type==CHEM_LINEAR_CONSTRAINT) {
                for (size_t j=0; j<con->n_reactants; ++j)
                    value+=con->reactants[j].yield*y[m->spec_to_var[con->reactants[j].species]*s->cells+c];
                value-=con->constant;
            } else {
                double K=s->custom[m->constraint_offsets[k]*s->cells+c];
                value=term_product(m,s->cells,c,y,con->n_reactants,con->reactants,CHEM_NONE,K)
                     -term_product(m,s->cells,c,y,con->n_products,con->products,CHEM_NONE,1);
            }
            f[row*s->cells+c]=value;
        }
    }
}

void chem_constraint_jacobian(const chem_state *s, const double *y, double *jac)
{
    const chem_model *m=s->model;
    for (size_t k=0; k<m->mechanism->n_constraints; ++k) {
        const chem_constraint *con=&m->mechanism->constraints[k];
        if (con->type==CHEM_CUSTOM_CONSTRAINT) { con->jacobian(s,y,jac,con->context); continue; }
        size_t row=m->constraint_rows[k];
        for (int side=0; side<(con->type==CHEM_LINEAR_CONSTRAINT ? 1 : 2); ++side) {
            size_t count=side ? con->n_products : con->n_reactants;
            const chem_product *terms=side ? con->products : con->reactants;
            for (size_t j=0; j<count; ++j) {
                size_t var=m->spec_to_var[terms[j].species], flat=m->lu.map[row*m->n_species+var];
                for (size_t c=0; c<s->cells; ++c) {
                    double derivative=terms[j].yield;
                    if (con->type==CHEM_EQUILIBRIUM_CONSTRAINT) {
                        double initial=side ? 1 : s->custom[m->constraint_offsets[k]*s->cells+c];
                        double prod=term_product(m,s->cells,c,y,count,terms,j,initial);
                        double conc=y[var*s->cells+c];
                        derivative=terms[j].yield==1 ? prod : conc>0 ? terms[j].yield*prod*pow(conc,terms[j].yield-1) : 0;
                    }
                    if (side) jac[flat*s->cells+c]+=derivative;
                    else jac[flat*s->cells+c]-=derivative;
                }
            }
        }
    }
}

void chem_alpha_minus_jacobian(chem_state *s, double alpha)
{
    const chem_model *m=s->model;
    for (size_t i=0; i<m->n_species; ++i) {
        double shift=m->algebraic[i] ? 0 : alpha;
        double *diag=s->jacobian+m->lu.diagonal[i]*s->cells;
        #pragma omp simd
        for (size_t c=0; c<s->cells; ++c) diag[c]+=shift;
    }
}

static chem_status algebraic_max(const chem_state *s, const double *x, double *maximum)
{
    *maximum=0;
    int has_nan=0,has_inf=0;
    for (size_t i=0; i<s->model->n_species; ++i) if (s->model->algebraic[i])
        for (size_t c=0; c<s->cells; ++c) {
            double v=fabs(x[i*s->cells+c]);
            has_nan|=isnan(v); has_inf|=isinf(v);
            if (isfinite(v) && v>*maximum) *maximum=v;
        }
    return has_nan ? CHEM_NAN : has_inf ? CHEM_INF : CHEM_OK;
}

static double correction_norm(const chem_state *s, const double *y, const double *delta)
{
    double norm=0;
    for (size_t i=0; i<s->model->n_species; ++i) if (s->model->algebraic[i])
        for (size_t c=0; c<s->cells; ++c) {
            size_t k=i*s->cells+c;
            double z=fabs(y[k]), corrected=fabs(y[k]+delta[k]);
            double scale=s->atol[i]+s->rtol*(z>corrected ? z : corrected);
            double q=scale>0 ? fabs(delta[k])/scale : delta[k]==0 ? 0 : INFINITY;
            norm=norm>q ? norm : q;
        }
    return norm;
}

chem_status chem_initialize_constraints(chem_state *s, const chem_ros_parameters *p, chem_stats *stats)
{
    const chem_model *m=s->model;
    if (!m->mechanism->n_constraints) return CHEM_OK;
    size_t n=m->n_species*s->cells;
    double *candidate=s->work, *delta=s->work+n, *trial=s->work+2*n, *original=s->work+3*n;
    int modified=0;
    chem_status status=CHEM_OK;
    for (size_t update=0; update<=p->init_max_iterations; ++update) {
        double maximum;
        memset(delta,0,n*sizeof(double)); chem_constraint_forcing(s,s->y,delta);
        status=algebraic_max(s,delta,&maximum);
        if (status!=CHEM_OK) goto restore;
        ++stats->constraint_init_iterations;
        if (maximum==0) return CHEM_OK;
        memset(s->jacobian,0,m->lu.nnz*s->cells*sizeof(double));
        chem_constraint_jacobian(s,s->y,s->jacobian); ++stats->jacobian_updates;
        chem_alpha_minus_jacobian(s,1);
        chem_lu_factor(&m->lu,s->cells,s->jacobian,s->scratch); ++stats->decompositions;
        chem_lu_solve(&m->lu,s->cells,s->jacobian,delta); ++stats->solves;
        status=algebraic_max(s,delta,&maximum);
        if (status!=CHEM_OK) goto restore;
        double norm=correction_norm(s,s->y,delta);
        if (norm<=p->init_tolerance) return CHEM_OK;
        if (update==p->init_max_iterations) break;
        if (!modified) { memcpy(original,s->y,n*sizeof(double)); modified=1; }
        if (p->init_max_backtracks) {
            double step=1;
            for (size_t backtrack=0;; ++backtrack,step*=p->init_backtrack_factor) {
                if (backtrack>p->init_max_backtracks) { status=CHEM_CONSTRAINT_INIT; goto restore; }
                memcpy(candidate,s->y,n*sizeof(double));
                for (size_t i=0; i<m->n_species; ++i) if (m->algebraic[i])
                    for (size_t c=0; c<s->cells; ++c) candidate[i*s->cells+c]+=step*delta[i*s->cells+c];
                if (algebraic_max(s,candidate,&maximum)!=CHEM_OK) continue;
                memset(trial,0,n*sizeof(double)); chem_constraint_forcing(s,candidate,trial);
                if (algebraic_max(s,trial,&maximum)!=CHEM_OK) continue;
                if (maximum==0) break;
                chem_lu_solve(&m->lu,s->cells,s->jacobian,trial); ++stats->solves;
                if (algebraic_max(s,trial,&maximum)!=CHEM_OK) continue;
                double next_norm=correction_norm(s,candidate,trial);
                if (next_norm<=p->init_tolerance || next_norm<(1-p->init_sufficient_decrease*step)*norm) break;
            }
            memcpy(s->y,candidate,n*sizeof(double));
        } else {
            for (size_t i=0; i<m->n_species; ++i) if (m->algebraic[i])
                for (size_t c=0; c<s->cells; ++c) s->y[i*s->cells+c]+=delta[i*s->cells+c];
        }
    }
    status=CHEM_CONSTRAINT_INIT;
restore:
    if (modified) memcpy(s->y,original,n*sizeof(double));
    return status;
}
