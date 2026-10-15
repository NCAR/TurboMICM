/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_builder.h"
#include <stdint.h>
#include <stdlib.h>

static ChemStatus build_extra(ChemSolverMethod method,size_t limit,int unused,int ignore,const ChemExternalModel *model,const ChemConstraint *constraint,ChemSolverPlan **out)
{
    ChemSpeciesInfo *species[4]={0};ChemPhaseSpecies entries[4]={{0}};const char *names[]={"A","B","C","D"};
    size_t count=unused?4:3;
    for (size_t i=0;i<count;++i) {
        CHECK(chem_species_info_create(names[i],0,NULL,&species[i])==CHEM_OK);
        CHECK(chem_phase_species_create(species[i],entries+i)==CHEM_OK);chem_species_info_free(species[i]);
    }
    ChemPhase *phase=NULL;ChemSystem *system=NULL;ChemReactionInfo *reactions[2]={0};
    CHECK(chem_phase_create("gas",count,entries,&phase)==CHEM_OK);
    for (size_t i=0;i<count;++i) chem_phase_species_free(entries+i);
    CHECK(chem_system_create(phase,&system)==CHEM_OK);
    for (size_t i=0;i<2;++i) {
        ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);
        rate.data.arrhenius.A=i?3.3e-11:2.15e-11;rate.data.arrhenius.C=i?55:110;
        ChemNamedProduct product={names[i+1],1};
        CHECK(chem_reaction_info_create(phase,1,names+i,1,&product,rate,&reactions[i])==CHEM_OK);
    }
    const ChemReactionInfo *inputs[]={reactions[0],reactions[1]};
    ChemSolverBuilder b=chem_solver_builder(method);b.system=system;b.reaction_count=2;b.reactions=inputs;
    b.maximum_cells=limit;b.ignore_unused_species=ignore;b.model_count=model?1:0;b.models=model;
    b.constraint_count=constraint?1:0;b.constraints=constraint;
    ChemStatus status=chem_solver_build(&b,out);
    for (size_t i=0;i<2;++i) chem_reaction_info_free(reactions[i]);
    chem_system_free(system);chem_phase_free(phase);return status;
}
static ChemStatus build(ChemSolverMethod method,size_t limit,int unused,int ignore,const ChemExternalModel *model,ChemSolverPlan **out)
{ return build_extra(method,limit,unused,ignore,model,NULL,out); }
static void missing(void)
{
    ChemSolverBuilder b=chem_solver_builder(CHEM_SOLVER_BACKWARD_EULER);ChemSolverPlan *p=NULL;
    CHECK(chem_solver_build(&b,&p)==CHEM_MISSING_SYSTEM);CHECK(p==NULL);
}
static void unused(void)
{
    ChemSolverPlan *p=NULL;CHECK(build(CHEM_SOLVER_ROSENBROCK,0,1,0,NULL,&p)==CHEM_UNUSED_SPECIES);CHECK(p==NULL);
    CHECK(build(CHEM_SOLVER_ROSENBROCK,0,1,1,NULL,&p)==CHEM_OK);chem_solver_plan_free(p);
    const char *names[]={"D"};ChemExternalModel model={.variable_count=1,.variables=names};
    CHECK(build(CHEM_SOLVER_ROSENBROCK,0,1,0,&model,&p)==CHEM_OK);chem_solver_plan_free(p);
    ChemConstraintTerm terms[]={{"D",1},{"A",-1}};
    ChemConstraint constraint={.type=CHEM_LINEAR_CONSTRAINT,.name="balance",.algebraic_species="D",.reactant_count=2,.reactants=terms};
    CHECK(build_extra(CHEM_SOLVER_ROSENBROCK,0,1,0,NULL,&constraint,&p)==CHEM_OK);
    CHECK(chem_constraint_count(chem_solver_plan_mechanism(p))==1);chem_solver_plan_free(p);
}
static void can_build(ChemSolverMethod method)
{
    for (size_t limit=0;limit<=4;limit+=4) {
        ChemSolverPlan *p=NULL;CHECK(build(method,limit,0,0,NULL,&p)==CHEM_OK);
        const ChemPhase *phase=chem_system_gas_phase(chem_solver_plan_system(p));
        CHECK(!strcmp(chem_phase_name(phase),"gas"));CHECK(chem_phase_count(phase)==3);
        if (limit) CHECK(chem_solver_plan_maximum_cells(p)==4);
        else CHECK(chem_solver_plan_maximum_cells(p)>100000000);
        ChemState *s=NULL;ChemSolver *solver=NULL;
        CHECK(chem_solver_plan_state(p,limit?4:7,&s)==CHEM_OK);
        CHECK(chem_solver_plan_workspace(p,limit?4:7,&solver)==CHEM_OK);
        chem_solver_free(solver);chem_state_free(s);
        if (limit) {
            CHECK(chem_solver_plan_state(p,5,&s)==CHEM_DIMENSION);CHECK(s==NULL);
            CHECK(chem_solver_plan_workspace(p,5,&solver)==CHEM_DIMENSION);CHECK(solver==NULL);
        }
        chem_solver_plan_free(p);
    }
}
static void statistics(const ChemSolverStats *s,double time,size_t functions,size_t jacobians,size_t steps,size_t solves)
{
    CHECK(s->final_time==time);CHECK(s->function_calls==functions);CHECK(s->jacobian_updates==jacobians);
    CHECK(s->internal_steps==steps);CHECK(s->solves==solves);
}
static void overrides(ChemSolverMethod method)
{
    ChemSolverPlan *p=NULL;ChemState *s=NULL;ChemSolver *solver=NULL;ChemSolverStats stats;
    CHECK(build(method,0,0,1,NULL,&p)==CHEM_OK);
    CHECK(chem_solver_plan_state(p,1,&s)==CHEM_OK);CHECK(chem_solver_plan_workspace(p,1,&solver)==CHEM_OK);
    ChemSolverOptions options=*chem_solver_plan_options(p);
    if (method==CHEM_SOLVER_BACKWARD_EULER) {
        CHECK(chem_solver_plan_solve(p,solver,s,5,&options,&stats)==CHEM_OK);statistics(&stats,5,2,2,2,2);
        options.parameters.backward_euler.small=1;options.parameters.backward_euler.max_iterations=1;
        CHECK(chem_solver_plan_solve(p,solver,s,5,&options,&stats)==CHEM_UNCONVERGED);statistics(&stats,.03125,6,6,6,6);
    } else {
        s->values[0]=1;CHECK(chem_solver_plan_solve(p,solver,s,5,NULL,&stats)==CHEM_OK);statistics(&stats,5,18,9,9,27);
        options.parameters.rosenbrock.h_min=15;options.parameters.rosenbrock.max_steps=6;
        s->values[0]=1;s->values[1]=s->values[2]=0;
        CHECK(chem_solver_plan_solve(p,solver,s,5,&options,&stats)==CHEM_OK);statistics(&stats,5,2,1,1,3);
        CHECK(chem_solver_plan_options(p)->parameters.rosenbrock.h_min==15);
        CHECK(chem_solver_plan_options(p)->parameters.rosenbrock.max_steps==6);
        CHECK(chem_solver_plan_solve(p,solver,s,5,NULL,&stats)==CHEM_OK);statistics(&stats,5,2,1,1,3);
    }
    chem_solver_free(solver);chem_state_free(s);chem_solver_plan_free(p);
}
typedef struct { size_t calls,dense_limit,sparse_limit,order;int dense_seen,sparse_seen; } Record;
typedef struct { Record *record; } Recorder;
static void *clone_record(const void *context)
{ Recorder *copy=malloc(sizeof(*copy));if (copy) *copy=*(const Recorder *)context;return copy; }
static ChemStatus record_sparse(const ChemModelBinding *binding,void *context,const ChemSparseStorage *sparse)
{
    Record *r=((Recorder *)context)->record;++r->calls;r->sparse_seen=1;
    r->sparse_limit=sparse->maximum_cells;r->order=chem_sparse_order(sparse->pattern);
    CHECK(binding->variable_count==1);CHECK(binding->variables[0]<r->order);
    CHECK(chem_sparse_index(sparse->pattern,binding->variables[0],binding->variables[0])!=SIZE_MAX);
    return CHEM_OK;
}
static ChemStatus record_dense(const ChemModelBinding *binding,void *context,const ChemDenseStorage *dense,const ChemSparseStorage *sparse)
{
    Record *r=((Recorder *)context)->record;r->dense_seen=dense->variable_major;r->dense_limit=dense->maximum_cells;
    return record_sparse(binding,context,sparse);
}
static void finalize(int legacy)
{
    for (size_t limit=legacy?4:0;limit<=4;limit+=4) {
        Record record={0};Recorder context={&record};const char *names[]={"A"};
        ChemExternalModel model={.variable_count=1,.variables=names,.context=&context,.clone=clone_record,.destroy=free};
        if (legacy) model.finalize_sparse=record_sparse;else model.finalize=record_dense;
        ChemSolverPlan *p=NULL;CHECK(build(CHEM_SOLVER_ROSENBROCK,limit,0,1,&model,&p)==CHEM_OK);
        CHECK(record.calls==1);CHECK(record.dense_seen==!legacy);CHECK(record.sparse_seen);
        CHECK(record.sparse_limit==(limit?limit:SIZE_MAX));CHECK(record.order==3);
        if (!legacy) CHECK(record.dense_limit==(limit?limit:SIZE_MAX));
        ChemSolver *solver=NULL;CHECK(chem_solver_plan_workspace(p,1,&solver)==CHEM_OK);
        CHECK(record.calls==1);chem_solver_free(solver);chem_solver_plan_free(p);
    }
}
void test_builder(void)
{
    RUN("solver_builder:SolverBuilder.ThrowsMissingSystem",missing());
    RUN("solver_builder:SolverBuilder.ThrowsUnusedSpecies",unused());
    RUN("solver_builder:SolverBuilder.CanBuildBackwardEuler",can_build(CHEM_SOLVER_BACKWARD_EULER));
    RUN("solver_builder:SolverBuilder.CanBuildRosenbrock",can_build(CHEM_SOLVER_ROSENBROCK));
    RUN("solver_builder:SolverBuilder.CanBuildBackwardEulerOverloadedSolverMethod",overrides(CHEM_SOLVER_BACKWARD_EULER));
    RUN("solver_builder:SolverBuilder.CanBuildRosenbrockOverloadedSolveMethod",overrides(CHEM_SOLVER_ROSENBROCK));
    RUN("solver_builder:SolverBuilder.FinalizeReceivesDenseMatrixPolicy",finalize(0));
    RUN("solver_builder:SolverBuilder.FinalizeSupportsSparseOnlySignature",finalize(1));
}
