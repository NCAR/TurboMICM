#include "rosenbrock.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int factor(const ChemRosenbrockOptions *o, double *matrix, size_t n)
{
    if (o->sparse_plan)
        return o->lu_algorithm==CHEM_LU_DOOLITTLE ? chem_sparse_factor_doolittle(o->sparse_plan,matrix) : chem_sparse_factor(o->sparse_plan,matrix);
    return o->lu_algorithm==CHEM_LU_DOOLITTLE ? chem_dense_factor_doolittle(matrix,n) : chem_dense_factor(matrix,n);
}
static void solve(const ChemRosenbrockOptions *o, const double *matrix, size_t n, double *delta)
{
    if (o->sparse_plan) chem_sparse_solve(o->sparse_plan,matrix,delta);
    else chem_dense_solve(matrix,n,delta);
}
static double correction_norm(const ChemConstraintSet *set, const double *y, const double *delta,
                              const ChemRosenbrockOptions *o)
{
    double norm=0;
    for (size_t k=0;k<set->count;++k) {
        size_t i=set->constraints[k].row;
        double scale=o->absolute_tolerance[i]+o->relative_tolerance*fmax(fabs(y[i]),fabs(y[i]+delta[i]));
        double q=scale>0 ? fabs(delta[i])/scale : (delta[i]==0 ? 0 : INFINITY);
        if (!isfinite(delta[i]) || !isfinite(q)) return INFINITY;
        norm=fmax(norm,q);
    }
    return norm;
}
static double residual_norm(const ChemConstraintSet *set, const double *residual)
{
    double norm=0;
    for (size_t k=0;k<set->count;++k) {
        double value=residual[set->constraints[k].row];
        if (!isfinite(value)) return INFINITY;
        norm=fmax(norm,fabs(value));
    }
    return norm;
}
ChemSolveStatus chem_rosenbrock_initialize_constraints(size_t n, double *y,
    const ChemRosenbrockOptions *o, ChemSolverStats *stats)
{
    ChemSolverStats local={0}; if (!stats) stats=&local;
    if (!o || !y || !n || !o->absolute_tolerance) return CHEM_SOLVE_INVALID;
    const ChemConstraintSet *set=o->constraints;
    if (!set || !set->count) return CHEM_SOLVE_OK;
    if (set->variables!=n || chem_constraint_validate(set)!=CHEM_CONSTRAINT_OK ||
        (o->sparse_plan && chem_sparse_plan_size(o->sparse_plan)!=n) ||
        (o->lu_algorithm!=CHEM_LU_MOZART && o->lu_algorithm!=CHEM_LU_DOOLITTLE) ||
        !isfinite(o->constraint_temperature) || o->constraint_temperature<=0 ||
        !isfinite(o->constraint_init_tolerance) || o->constraint_init_tolerance<0 ||
        !isfinite(o->constraint_init_backtrack_factor) || o->constraint_init_backtrack_factor<=0 || o->constraint_init_backtrack_factor>=1 ||
        !isfinite(o->constraint_init_sufficient_decrease) || o->constraint_init_sufficient_decrease<0 ||
        !isfinite(o->relative_tolerance) || o->relative_tolerance<0) return CHEM_SOLVE_INVALID;
    for (size_t i=0;i<n;++i) if (!isfinite(o->absolute_tolerance[i]) || o->absolute_tolerance[i]<0) return CHEM_SOLVE_INVALID;
    if (n>SIZE_MAX/n || n*n>SIZE_MAX/sizeof(double) || n>(SIZE_MAX/sizeof(double)-n*n)/4) return CHEM_SOLVE_ALLOCATION;
    double *work=calloc(n*n+4*n,sizeof(double));
    unsigned char *algebraic=calloc(n,1);
    if (!work || !algebraic) { free(work); free(algebraic); return CHEM_SOLVE_ALLOCATION; }
    double *matrix=work,*delta=matrix+n*n,*original=delta+n,*candidate=original+n,*trial=candidate+n;
    for (size_t i=0;i<set->count;++i) algebraic[set->constraints[i].row]=1;
    memcpy(original,y,n*sizeof(double));
    ChemSolveStatus status=CHEM_SOLVE_CONSTRAINT_INIT;
    for (size_t update=0;update<=o->constraint_init_max_iterations;++update) {
        memset(delta,0,n*sizeof(double));
        chem_constraint_residuals(set,1,y,&o->constraint_temperature,o->constraint_parameters,delta);
        double residual=residual_norm(set,delta);
        if (!isfinite(residual)) { status=CHEM_SOLVE_NONFINITE; break; }
        ++stats->constraint_init_iterations;
        if (residual==0) { status=CHEM_SOLVE_OK; break; }
        memset(matrix,0,n*n*sizeof(double));
        chem_constraint_subtract_jacobian(set,1,y,&o->constraint_temperature,o->constraint_parameters,matrix);
        ++stats->jacobian_updates;
        for (size_t i=0;i<n;++i) if (!algebraic[i]) matrix[i*n+i]+=1;
        ++stats->decompositions;
        if (!factor(o,matrix,n)) { status=CHEM_SOLVE_SINGULAR; break; }
        solve(o,matrix,n,delta); ++stats->solves;
        double current=correction_norm(set,y,delta,o);
        if (!isfinite(current)) { status=CHEM_SOLVE_NONFINITE; break; }
        if (current<=o->constraint_init_tolerance) { status=CHEM_SOLVE_OK; break; }
        if (update==o->constraint_init_max_iterations) break;
        if (!o->constraint_init_max_backtracks) {
            for (size_t i=0;i<n;++i) if (algebraic[i]) y[i]+=delta[i];
            continue;
        }
        int accepted=0; double step=1;
        for (size_t backtrack=0;backtrack<=o->constraint_init_max_backtracks;++backtrack) {
            memcpy(candidate,y,n*sizeof(double));
            for (size_t i=0;i<n;++i) if (algebraic[i]) candidate[i]+=step*delta[i];
            int finite=1;
            for (size_t i=0;i<n;++i) if (algebraic[i] && !isfinite(candidate[i])) finite=0;
            if (finite) {
                memset(trial,0,n*sizeof(double));
                chem_constraint_residuals(set,1,candidate,&o->constraint_temperature,o->constraint_parameters,trial);
                double r=residual_norm(set,trial);
                if (r==0) { accepted=1; break; }
                if (isfinite(r)) {
                    solve(o,matrix,n,trial); ++stats->solves;
                    double q=correction_norm(set,candidate,trial,o);
                    if (q<=o->constraint_init_tolerance || q<(1-o->constraint_init_sufficient_decrease*step)*current) {
                        accepted=1; break;
                    }
                }
            }
            step*=o->constraint_init_backtrack_factor;
        }
        if (!accepted) break;
        memcpy(y,candidate,n*sizeof(double));
    }
    if (status!=CHEM_SOLVE_OK) memcpy(y,original,n*sizeof(double));
    free(algebraic); free(work); return status;
}
