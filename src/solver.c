#include "solver.h"
#include <stdlib.h>
#include <string.h>

ChemSolverBuilder chem_solver_builder(ChemSolverMethod method)
{
    ChemSolverBuilder builder={0};builder.method=method;builder.rosenbrock=chem_rosenbrock_three_stage();
    builder.backward_euler=chem_backward_euler_defaults();return builder;
}
static int unused(const ChemSolverBuilder *b)
{
    size_t n=b->mechanism->species_count;if(n && n>SIZE_MAX/n)return -1;
    unsigned char *used=calloc(n,1);if(!used)return -1;
    for(size_t i=0;i<b->mechanism->reaction_count;++i) {
        const ChemReaction *r=&b->mechanism->reactions[i];
        for(size_t j=0;j<r->reactant_count;++j) {int s=b->mechanism->reactants[r->reactant_begin+j].species;if(s>=0 && (size_t)s<n)used[s]=1;}
        for(size_t j=0;j<r->product_count;++j) {int s=b->mechanism->products[r->product_begin+j].species;if(s>=0 && (size_t)s<n)used[s]=1;}
    }
    for(size_t i=0;i<b->constraint_count;++i) {
        const ChemConstraint *c=&b->constraints[i];
        if(c->row<n) used[c->row]=1;
        for(size_t j=0;j<c->reactant_count;++j) if(c->reactants[j].species>=0 && (size_t)c->reactants[j].species<n)used[c->reactants[j].species]=1;
        for(size_t j=0;j<c->product_count;++j) if(c->products[j].species>=0 && (size_t)c->products[j].species<n)used[c->products[j].species]=1;
    }
    for(size_t i=0;i<b->model_count;++i) if(b->models[i].jacobian_pattern) for(size_t j=0;j<n*n;++j) if(b->models[i].jacobian_pattern[j]) {used[j/n]=1;used[j%n]=1;}
    for(size_t i=0;i<b->model_count;++i) if(b->models[i].constraint_pattern) for(size_t j=0;j<n*n;++j) if(b->models[i].constraint_pattern[j]) {used[j/n]=1;used[j%n]=1;}
    int found=0;for(size_t i=0;i<n;++i)if(!used[i])found=1;free(used);return found;
}
ChemBuildStatus chem_solver_build(const ChemSolverBuilder *builder,ChemSolver *solver)
{
    if(!builder || !solver)return CHEM_BUILD_INVALID;
    if(!builder->system)return CHEM_BUILD_MISSING_SYSTEM;
    if(!builder->mechanism || chem_validate_mechanism(builder->mechanism,NULL,0)!=CHEM_VALID ||
       (builder->constraint_count && !builder->constraints) || (builder->model_count && !builder->models) ||
       (builder->method!=CHEM_SOLVER_ROSENBROCK && builder->method!=CHEM_SOLVER_BACKWARD_EULER))return CHEM_BUILD_INVALID;
    if(!builder->ignore_unused_species) {int status=unused(builder);if(status<0)return CHEM_BUILD_NO_MEMORY;if(status)return CHEM_BUILD_UNUSED_SPECIES;}
    ChemWorkspace *workspace=chem_workspace_create_external(builder->mechanism,builder->constraints,builder->constraint_count,builder->models,builder->model_count);
    if(!workspace)return CHEM_BUILD_INVALID;
    size_t n=builder->mechanism->species_count;unsigned char *pattern=malloc(n*n);if(!pattern){chem_workspace_free(workspace);return CHEM_BUILD_NO_MEMORY;}
    chem_workspace_jacobian_pattern(workspace,pattern,n*n);
    for(size_t i=0;i<builder->model_count;++i) {
        const ChemExternalModel *model=&builder->models[i];int status=0;
        if(model->finalize_dense) status=model->finalize_dense(model->context,builder->group_width?builder->group_width:1,pattern,n);
        else if(model->finalize_sparse) status=model->finalize_sparse(model->context,pattern,n);
        if(status){free(pattern);chem_workspace_free(workspace);return CHEM_BUILD_INVALID;}
    }
    free(pattern);*solver=(ChemSolver){*builder,workspace};return CHEM_BUILD_OK;
}
void chem_solver_free(ChemSolver *solver)
{if(!solver)return;chem_workspace_free(solver->workspace);*solver=(ChemSolver){0};}
size_t chem_solver_maximum_cells(const ChemSolver *solver)
{return !solver?0:solver->options.group_width?solver->options.group_width:SIZE_MAX;}
int chem_solver_solve(ChemSolver *solver,double dt,ChemState *state,const void *override,ChemStatistics *statistics)
{
    if(!solver || !solver->workspace || !state || !statistics || !state->cells || state->cells>chem_solver_maximum_cells(solver) ||
       state->parameters.variable_count!=solver->options.mechanism->species_count ||
       state->parameters.rate_count!=solver->options.mechanism->reaction_count)return CHEM_SOLVE_BAD_ARGUMENT;
    const ChemMechanism *m=solver->options.mechanism;size_t n=m->species_count,nr=m->reaction_count,cells=state->cells;
    if(n>SIZE_MAX-nr || cells>SIZE_MAX/(n+nr) || cells*(n+nr)>SIZE_MAX/sizeof(double) || cells>SIZE_MAX/sizeof(ChemWorkspace *))return CHEM_SOLVE_BAD_ARGUMENT;
    double *data=malloc(cells*(n+nr)*sizeof(double));if(!data)return CHEM_SOLVE_BAD_ARGUMENT;
    double *y=data,*rates=data+cells*n;
    for(size_t cell=0;cell<cells;++cell) {
        if(chem_matrix_get_row(&state->variables,cell,y+cell*n,n) || chem_matrix_get_row(&state->rate_constants,cell,rates+cell*nr,nr)){free(data);return CHEM_SOLVE_BAD_ARGUMENT;}
    }
    int status=0;ChemWorkspace **work=NULL;
    if(cells>1) {
        work=calloc(cells,sizeof(*work));
        if(!work){free(data);return CHEM_SOLVE_BAD_ARGUMENT;}
        work[0]=solver->workspace;
        for(size_t cell=1;cell<cells;++cell) {
            work[cell]=chem_workspace_create_external(m,solver->options.constraints,solver->options.constraint_count,solver->options.models,solver->options.model_count);
            if(!work[cell]){status=CHEM_SOLVE_BAD_ARGUMENT;break;}
        }
    }
    if(solver->options.method==CHEM_SOLVER_ROSENBROCK) {
        if(override)solver->options.rosenbrock=*(const ChemRosenbrockParameters *)override;
        ChemRosenbrockParameters p=solver->options.rosenbrock;
        p.absolute_tolerances=state->absolute_tolerances;p.relative_tolerance=state->parameters.relative_tolerance;
        if(!status)status=cells==1?chem_rosenbrock_solve(m,solver->workspace,&p,state->conditions[0],rates,dt,y,statistics):
                                  chem_rosenbrock_solve_group(m,work,cells,&p,state->conditions,rates,dt,y,statistics);
    } else {
        if(override)solver->options.backward_euler=*(const ChemBackwardEulerParameters *)override;
        ChemBackwardEulerParameters p=solver->options.backward_euler;
        p.absolute_tolerances=state->absolute_tolerances;p.relative_tolerance=state->parameters.relative_tolerance;
        if(!status)status=cells==1?chem_backward_euler_solve(m,solver->workspace,&p,state->conditions[0],rates,dt,y,statistics):
                                  chem_backward_euler_solve_group(m,work,cells,&p,state->conditions,rates,dt,y,statistics);
    }
    if(work){for(size_t cell=1;cell<cells;++cell)chem_workspace_free(work[cell]);free(work);}
    for(size_t cell=0;cell<cells;++cell)chem_matrix_set_row(&state->variables,cell,y+cell*n,n);
    free(data);return status;
}
