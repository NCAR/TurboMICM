#include "rosenbrock.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    size_t cells, variables;
    ChemOdeFunction rhs, jacobian;
    void *const *contexts;
    double *cell_jacobian;
} Group;

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

ChemSolveStatus chem_rosenbrock_solve_grouped(size_t cells, size_t n,
    double *y, double duration, ChemOdeFunction rhs, ChemOdeFunction jacobian,
    void *const *contexts, const unsigned char *cell_pattern,
    const ChemRosenbrockMethod *method, const ChemRosenbrockOptions *options,
    ChemSolverStats *stats)
{
    if (stats) memset(stats,0,sizeof(*stats));
    if (!cells || !n || !y || !rhs || !jacobian || !contexts || !method || !options ||
        !options->absolute_tolerance || options->sparse_plan ||
        (options->constraints && options->constraints->count)) return CHEM_SOLVE_INVALID;
    if (cells>SIZE_MAX/n) return CHEM_SOLVE_ALLOCATION;
    size_t total=cells*n;
    if (total>SIZE_MAX/total || total*total>SIZE_MAX/sizeof(double) ||
        n>SIZE_MAX/n || n*n>SIZE_MAX/sizeof(double) || total>SIZE_MAX/sizeof(double))
        return CHEM_SOLVE_ALLOCATION;
    double *atol=malloc(total*sizeof(*atol));
    double *cell_jacobian=malloc(n*n*sizeof(*cell_jacobian));
    if (!atol || !cell_jacobian) { free(atol); free(cell_jacobian); return CHEM_SOLVE_ALLOCATION; }
    for (size_t c=0;c<cells;++c)
        memcpy(atol+c*n,options->absolute_tolerance,n*sizeof(*atol));
    ChemSparsePlan *plan=NULL;
    if (cell_pattern) {
        unsigned char *pattern=calloc(total*total,1);
        if (!pattern) { free(atol); free(cell_jacobian); return CHEM_SOLVE_ALLOCATION; }
        for (size_t c=0;c<cells;++c) for (size_t i=0;i<n;++i)
            memcpy(pattern+(c*n+i)*total+c*n,cell_pattern+i*n,n);
        plan=chem_sparse_plan_from_pattern(total,pattern);
        free(pattern);
        if (!plan) { free(atol); free(cell_jacobian); return CHEM_SOLVE_ALLOCATION; }
    }
    Group group={cells,n,rhs,jacobian,contexts,cell_jacobian};
    ChemRosenbrockOptions copy=*options;
    copy.absolute_tolerance=atol; copy.sparse_plan=plan;
    ChemSolveStatus result=chem_rosenbrock_solve(total,y,duration,group_rhs,group_jacobian,
                                                 &group,method,&copy,stats);
    chem_sparse_plan_destroy(plan);
    free(atol); free(cell_jacobian);
    return result;
}
