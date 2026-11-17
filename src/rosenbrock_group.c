#include "rosenbrock.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    size_t cells, variables;
    ChemOdeFunction rhs, jacobian;
    void *const *contexts;
    double *cell_jacobian;
} Group;

typedef struct {
    size_t n, offset;
    const ChemConstraint *constraint;
    double temperature, *scratch;
} GroupConstraint;

static double grouped_residual(void *v,const double *y,double temperature,const double *parameter)
{
    GroupConstraint *g=v; (void)temperature;
    ChemConstraintSet one={g->n,1,g->constraint};
    chem_constraint_residuals(&one,1,y+g->offset,&g->temperature,parameter,g->scratch);
    return g->scratch[g->constraint->row];
}

static void grouped_subtract_jacobian(void *v,const double *y,double temperature,
                                      const double *parameter,double *row)
{
    GroupConstraint *g=v; (void)temperature;
    ChemConstraintSet one={g->n,1,g->constraint};
    memset(g->scratch,0,g->n*g->n*sizeof(*g->scratch));
    chem_constraint_subtract_jacobian(&one,1,y+g->offset,&g->temperature,parameter,g->scratch);
    const double *local=g->scratch+g->constraint->row*g->n;
    for (size_t i=0;i<g->n;++i) row[g->offset+i]+=local[i];
}

static void group_rhs(void *v, const double *y, double *out)
{
    Group *g=v;
    for (size_t c=0;c<g->cells;++c)
        g->rhs(g->contexts[c],y+c*g->variables,out+c*g->variables);
}

static void group_jacobian(void *v, const double *y, double *out)
{
    Group *g=v;
    size_t n=g->variables,total=g->cells*n;
    memset(out,0,total*total*sizeof(*out));
    for (size_t c=0;c<g->cells;++c) {
        g->jacobian(g->contexts[c],y+c*n,g->cell_jacobian);
        for (size_t i=0;i<n;++i)
            memcpy(out+(c*n+i)*total+c*n,g->cell_jacobian+i*n,n*sizeof(*out));
    }
}

ChemSolveStatus chem_rosenbrock_solve_grouped_temperatures(size_t cells, size_t n,
    double *y, double duration, ChemOdeFunction rhs, ChemOdeFunction jacobian,
    void *const *contexts, const unsigned char *cell_pattern, const double *temperatures,
    const ChemRosenbrockMethod *method, const ChemRosenbrockOptions *options,
    ChemSolverStats *stats)
{
    if (stats) memset(stats,0,sizeof(*stats));
    if (!cells || !n || !y || !rhs || !jacobian || !contexts || !method || !options ||
        !options->absolute_tolerance || options->sparse_plan) return CHEM_SOLVE_INVALID;
    if (cells>SIZE_MAX/n) return CHEM_SOLVE_ALLOCATION;
    size_t total=cells*n;
    if (total>SIZE_MAX/total || total*total>SIZE_MAX/sizeof(double) ||
        n>SIZE_MAX/n || n*n>SIZE_MAX/sizeof(double) || total>SIZE_MAX/sizeof(double))
        return CHEM_SOLVE_ALLOCATION;
    const ChemConstraintSet *original=options->constraints;
    size_t per_cell=original ? original->count : 0;
    if (per_cell && (original->variables!=n || chem_constraint_validate(original)!=CHEM_CONSTRAINT_OK))
        return CHEM_SOLVE_INVALID;
    if (per_cell && (cells>SIZE_MAX/per_cell || cells*per_cell>SIZE_MAX/n ||
        cells*per_cell*n>SIZE_MAX/sizeof(size_t) || cells*per_cell>SIZE_MAX/sizeof(GroupConstraint) ||
        cells*per_cell>SIZE_MAX/sizeof(ChemConstraint) || cells*per_cell>SIZE_MAX/sizeof(ChemCustomConstraint)))
        return CHEM_SOLVE_ALLOCATION;
    double *atol=malloc(total*sizeof(*atol));
    double *cell_jacobian=malloc(n*n*sizeof(*cell_jacobian));
    if (!atol || !cell_jacobian) { free(atol); free(cell_jacobian); return CHEM_SOLVE_ALLOCATION; }
    size_t count=cells*per_cell;
    GroupConstraint *constraint_context=count ? calloc(count,sizeof(*constraint_context)) : NULL;
    ChemCustomConstraint *custom=count ? calloc(count,sizeof(*custom)) : NULL;
    ChemConstraint *constraints=count ? calloc(count,sizeof(*constraints)) : NULL;
    size_t *dependencies=count ? malloc(count*n*sizeof(*dependencies)) : NULL;
    double *constraint_scratch=count ? malloc(n*n*sizeof(*constraint_scratch)) : NULL;
    if (count && (!constraint_context || !custom || !constraints || !dependencies || !constraint_scratch)) {
        free(atol); free(cell_jacobian); free(constraint_context); free(custom); free(constraints);
        free(dependencies); free(constraint_scratch); return CHEM_SOLVE_ALLOCATION;
    }
    ChemSolveStatus result=CHEM_SOLVE_OK;
    for (size_t c=0;c<cells && count;++c) {
        double temperature=temperatures ? temperatures[c] : options->constraint_temperature;
        if (!isfinite(temperature) || temperature<=0) { result=CHEM_SOLVE_INVALID; goto done; }
        for (size_t j=0;j<per_cell;++j) {
            size_t k=c*per_cell+j;
            constraint_context[k]=(GroupConstraint){n,c*n,original->constraints+j,temperature,constraint_scratch};
            for (size_t i=0;i<n;++i) dependencies[k*n+i]=c*n+i;
            custom[k]=(ChemCustomConstraint){constraint_context+k,dependencies+k*n,n,grouped_residual,grouped_subtract_jacobian};
            constraints[k]=(ChemConstraint){.name=original->constraints[j].name,.type=CHEM_CONSTRAINT_CUSTOM,
                .row=c*n+original->constraints[j].row,.custom=custom+k};
        }
    }
    ChemConstraintSet grouped_constraints={total,count,constraints};
    for (size_t c=0;c<cells;++c)
        memcpy(atol+c*n,options->absolute_tolerance,n*sizeof(*atol));
    ChemSparsePlan *plan=NULL;
    if (cell_pattern) {
        unsigned char *pattern=calloc(total*total,1);
        if (!pattern) { result=CHEM_SOLVE_ALLOCATION; goto done; }
        for (size_t c=0;c<cells;++c) for (size_t i=0;i<n;++i)
            memcpy(pattern+(c*n+i)*total+c*n,cell_pattern+i*n,n);
        if (count) chem_constraint_sparsity(&grouped_constraints,pattern);
        plan=chem_sparse_plan_from_pattern(total,pattern);
        free(pattern);
        if (!plan) { result=CHEM_SOLVE_ALLOCATION; goto done; }
    }
    Group group={cells,n,rhs,jacobian,contexts,cell_jacobian};
    ChemRosenbrockOptions copy=*options;
    copy.absolute_tolerance=atol; copy.sparse_plan=plan;
    copy.constraints=count ? &grouped_constraints : NULL;
    if (count && temperatures) copy.constraint_temperature=temperatures[0];
    result=chem_rosenbrock_solve(total,y,duration,group_rhs,group_jacobian,
                                 &group,method,&copy,stats);
    chem_sparse_plan_destroy(plan);
done:
    free(atol); free(cell_jacobian);
    free(constraint_context); free(custom); free(constraints); free(dependencies); free(constraint_scratch);
    return result;
}

ChemSolveStatus chem_rosenbrock_solve_grouped(size_t cells, size_t n,
    double *y, double duration, ChemOdeFunction rhs, ChemOdeFunction jacobian,
    void *const *contexts, const unsigned char *cell_pattern,
    const ChemRosenbrockMethod *method, const ChemRosenbrockOptions *options,
    ChemSolverStats *stats)
{
    return chem_rosenbrock_solve_grouped_temperatures(cells,n,y,duration,rhs,jacobian,
        contexts,cell_pattern,NULL,method,options,stats);
}
