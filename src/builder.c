/* SPDX-License-Identifier: Apache-2.0 */
#include "chem_builder.h"
#include "solver_internal.h"

struct ChemSolverPlan {
    ChemMechanism *mechanism;
    ChemSparse *pattern;
    size_t maximum_cells;
    ChemSolverOptions options;
};
ChemSolverBuilder chem_solver_builder(ChemSolverMethod method)
{
    ChemSolverBuilder b={.ignore_unused_species=1,.reorder=1,.options={.method=method}};
    if (method==CHEM_SOLVER_BACKWARD_EULER) b.options.parameters.backward_euler=chem_backward_euler_parameters();
    else b.options.parameters.rosenbrock=chem_rosenbrock_parameters(CHEM_ROS3);
    return b;
}
void chem_solver_plan_free(ChemSolverPlan *p)
{ if (p) { chem_sparse_free(p->pattern);chem_mechanism_free(p->mechanism);free(p); } }
static ChemStatus unused_species(const ChemMechanism *m)
{
    unsigned char *used=chem_alloc(m->variable_count,1);
    if (!used) return CHEM_MEMORY;
    for (size_t i=0;i<m->reaction_count;++i) {
        const ChemProcess *r=m->processes+i;
        for (size_t j=0;j<r->reactant_count;++j) used[r->reactants[j]]=1;
        for (size_t j=0;j<r->product_count;++j) used[r->products[j].species]=1;
    }
    for (size_t i=0;i<m->model_count;++i) {
        const ChemModelBinding *b=&m->models[i].binding;
        for (size_t j=0;j<b->variable_count;++j) used[b->variables[j]]=1;
        for (size_t j=0;j<b->algebraic_count;++j) used[b->algebraic[j]]=1;
        for (size_t kind=0;kind<2;++kind) {
            size_t count=kind?b->constraint_element_count:b->process_element_count;
            const ChemJacobianElement *e=kind?b->constraint_elements:b->process_elements;
            for (size_t j=0;j<count;++j) used[e[j].row]=used[e[j].column]=1;
        }
    }
    for (size_t i=0;i<m->constraint_count;++i) {
        const ChemBoundConstraint *c=m->constraints+i;used[c->row]=1;
        for (size_t j=0;j<c->definition.reactant_count;++j) used[c->reactants[j]]=1;
        for (size_t j=0;j<c->definition.product_count;++j) used[c->products[j]]=1;
    }
    ChemStatus status=CHEM_OK;
    for (size_t i=0;i<m->variable_count;++i) if (!used[i]) { status=CHEM_UNUSED_SPECIES;break; }
    free(used);return status;
}
static int valid_method(ChemSolverMethod m)
{ return m==CHEM_SOLVER_ROSENBROCK || m==CHEM_SOLVER_BACKWARD_EULER; }
ChemStatus chem_solver_build(const ChemSolverBuilder *b,ChemSolverPlan **out)
{
    if (!out) return CHEM_INVALID;
    *out=NULL;
    if (!b) return CHEM_INVALID;
    if (!b->system) return CHEM_MISSING_SYSTEM;
    if (!valid_method(b->options.method)) return CHEM_INVALID;
    ChemSolverPlan *p=calloc(1,sizeof(*p));
    if (!p) return CHEM_MEMORY;
    p->options=b->options;p->maximum_cells=b->maximum_cells?b->maximum_cells:SIZE_MAX;
    ChemStatus status=chem_system_mechanism_with_models(b->system,b->reaction_count,b->reactions,b->model_count,b->models,&p->mechanism);
    if (status) goto fail;
    status=chem_mechanism_set_constraints(p->mechanism,b->constraint_count,b->constraints);
    if (status) goto fail;
    if (!b->ignore_unused_species) { status=unused_species(p->mechanism);if (status) goto fail; }
    if (b->reorder) { status=chem_mechanism_reorder(p->mechanism);if (status) goto fail; }
    unsigned char *pattern=chem_mechanism_pattern(p->mechanism);
    if (!pattern) { status=CHEM_MEMORY;goto fail; }
    status=chem_sparse_create(p->mechanism->variable_count,pattern,&p->pattern);free(pattern);
    if (status) goto fail;
    ChemDenseStorage dense={p->maximum_cells,1};ChemSparseStorage sparse={p->maximum_cells,p->pattern};
    for (size_t i=0;i<p->mechanism->model_count;++i) {
        ChemBoundModel *model=p->mechanism->models+i;
        if (model->callbacks.finalize) status=model->callbacks.finalize(&model->binding,model->context,&dense,&sparse);
        else if (model->callbacks.finalize_sparse) status=model->callbacks.finalize_sparse(&model->binding,model->context,&sparse);
        if (status) goto fail;
    }
    *out=p;return CHEM_OK;
fail:
    chem_solver_plan_free(p);return status;
}
const ChemMechanism *chem_solver_plan_mechanism(const ChemSolverPlan *p) { return p?p->mechanism:NULL; }
const ChemSystem *chem_solver_plan_system(const ChemSolverPlan *p) { return p?chem_mechanism_system(p->mechanism):NULL; }
size_t chem_solver_plan_maximum_cells(const ChemSolverPlan *p) { return p?p->maximum_cells:0; }
const ChemSolverOptions *chem_solver_plan_options(const ChemSolverPlan *p) { return p?&p->options:NULL; }
ChemStatus chem_solver_plan_state(const ChemSolverPlan *p,size_t cells,ChemState **out)
{
    if (!out) return CHEM_INVALID;
    *out=NULL;
    if (!p) return CHEM_INVALID;
    if (cells>p->maximum_cells) return CHEM_DIMENSION;
    return chem_state_create(p->mechanism,cells,out);
}
ChemStatus chem_solver_plan_workspace(const ChemSolverPlan *p,size_t cells,ChemSolver **out)
{
    if (!out) return CHEM_INVALID;
    *out=NULL;
    if (!p) return CHEM_INVALID;
    if (cells>p->maximum_cells) return CHEM_DIMENSION;
    return chem_solver_create(p->mechanism,cells,out);
}
ChemStatus chem_solver_plan_solve(ChemSolverPlan *p,ChemSolver *solver,ChemState *state,double dt,
                                  const ChemSolverOptions *override,ChemSolverStats *stats)
{
    if (!p || !solver || !state || solver->mechanism!=p->mechanism || state->mechanism!=p->mechanism ||
        state->cells>p->maximum_cells || (override && !valid_method(override->method))) return CHEM_INVALID;
    if (override) p->options=*override;
    if (p->options.method==CHEM_SOLVER_BACKWARD_EULER)
        return chem_solve_backward_euler(solver,state,dt,&p->options.parameters.backward_euler,stats);
    return chem_solve(solver,state,dt,&p->options.parameters.rosenbrock,stats);
}
