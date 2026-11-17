#include "internal.h"

#include <math.h>
#include <string.h>

static int finite_algebraic(const ChemState *s, const double *values, double *maximum)
{
    size_t n=s->mechanism->species_count, cells=s->cells;
    *maximum=0;
    for (size_t i=0; i<n; ++i) if (s->algebraic[i])
        for (size_t c=0; c<cells; ++c) {
            double v=values[i*cells+c];
            if (!isfinite(v)) return 0;
            if (fabs(v)>*maximum) *maximum=fabs(v);
        }
    return 1;
}

static double weighted_correction(const ChemState *s, const ChemRosenbrock *p,
                                  const double *y, const double *delta)
{
    size_t n=s->mechanism->species_count, cells=s->cells;
    double maximum=0;
    for (size_t i=0; i<n; ++i) if (s->algebraic[i])
        for (size_t c=0; c<cells; ++c) {
            size_t k=i*cells+c;
            double scale=s->absolute_tolerances[i]+p->relative_tolerance*fmax(fabs(y[k]),fabs(y[k]+delta[k]));
            double q=scale>0 ? fabs(delta[k])/scale : (delta[k]==0 ? 0 : INFINITY);
            if (q>maximum) maximum=q;
        }
    return maximum;
}

int chem_project_constraints(ChemState *s, const ChemRosenbrock *p, ChemStats *stats)
{
    size_t n=s->mechanism->species_count, cells=s->cells, ny=n*cells;
    double *delta=s->forcing, *original=s->stages, *trial=s->error, *candidate=s->candidate;
    int modified=0, result=CHEM_CONSTRAINT_INITIALIZATION_FAILED;
    for (unsigned update=0; update<=p->constraint_max_iterations; ++update) {
        double maximum;
        memset(delta,0,ny*sizeof(*delta));
        chem_constraint_residuals(s,s->values,delta);
        if (!finite_algebraic(s,delta,&maximum)) { result=CHEM_NONFINITE; break; }
        ++stats->constraint_init_iterations;
        if (maximum==0) return CHEM_OK;
        memset(s->matrix,0,s->nonzeros*cells*sizeof(*s->matrix));
        chem_constraint_jacobian(s,s->values); ++stats->jacobian_updates;
        int factor_status=chem_linear_factor(s,1); ++stats->decompositions;
        if (factor_status) break;
        int solve_status=chem_linear_solve(s,delta); ++stats->solves;
        if (solve_status) { result=solve_status; break; }
        if (!finite_algebraic(s,delta,&maximum)) { result=CHEM_NONFINITE; break; }
        double current=weighted_correction(s,p,s->values,delta);
        if (current <= p->constraint_tolerance) return CHEM_OK;
        if (update == p->constraint_max_iterations) break;
        if (!modified) { memcpy(original,s->values,ny*sizeof(*original)); modified=1; }
        if (p->constraint_max_backtracks) {
            double step=1;
            int accepted=0;
            for (unsigned backtrack=0; backtrack<=p->constraint_max_backtracks; ++backtrack) {
                memcpy(candidate,s->values,ny*sizeof(*candidate));
                for (size_t i=0; i<n; ++i) if (s->algebraic[i])
                    for (size_t c=0; c<cells; ++c) candidate[i*cells+c]+=step*delta[i*cells+c];
                if (!finite_algebraic(s,candidate,&maximum)) { step*=p->constraint_backtrack_factor; continue; }
                memset(trial,0,ny*sizeof(*trial));
                chem_constraint_residuals(s,candidate,trial);
                if (!finite_algebraic(s,trial,&maximum)) { step*=p->constraint_backtrack_factor; continue; }
                if (maximum==0) { accepted=1; break; }
                solve_status=chem_linear_solve(s,trial); ++stats->solves;
                if (solve_status) { result=solve_status; break; }
                if (!finite_algebraic(s,trial,&maximum)) { step*=p->constraint_backtrack_factor; continue; }
                double q=weighted_correction(s,p,candidate,trial);
                if (q<=p->constraint_tolerance || q<(1-p->constraint_sufficient_decrease*step)*current) {
                    accepted=1; break;
                }
                step*=p->constraint_backtrack_factor;
            }
            if (!accepted) break;
            memcpy(s->values,candidate,ny*sizeof(*candidate));
        } else {
            for (size_t i=0; i<n; ++i) if (s->algebraic[i])
                for (size_t c=0; c<cells; ++c) s->values[i*cells+c]+=delta[i*cells+c];
            if (!finite_algebraic(s,s->values,&maximum)) { result=CHEM_NONFINITE; break; }
        }
    }
    if (modified) memcpy(s->values,original,ny*sizeof(*original));
    return result;
}

int chem_initialize_constraints(ChemState *s, const ChemRosenbrock *p, ChemStats *stats)
{
    if (!s || !p || !stats || p->constraint_tolerance < 0 ||
        p->constraint_backtrack_factor <= 0 || p->constraint_backtrack_factor >= 1 ||
        p->constraint_sufficient_decrease < 0 || p->relative_tolerance < 0) return CHEM_INVALID;
    memset(stats,0,sizeof(*stats));
    return s->constraint_count ? chem_project_constraints(s,p,stats) : CHEM_OK;
}
