#include "solver.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {int calls;size_t width,order;int sparse_only;} Record;
static int no_op(void *ctx,ChemConditions c,const double *y,double *out)
{(void)ctx;(void)c;(void)y;(void)out;return 0;}
static int dense_finalize(void *context,size_t width,const unsigned char *pattern,size_t n)
{Record *r=context;if(!pattern)return 1;r->calls++;r->width=width;r->order=n;return 0;}
static int sparse_finalize(void *context,const unsigned char *pattern,size_t n)
{Record *r=context;if(!pattern)return 1;r->calls++;r->order=n;r->sparse_only=1;return 0;}
static int run(int kind)
{
    const char *names[]={"A","B","C","D"};ChemSpecies species[4];ChemPhaseSpecies members[4];
    for(size_t i=0;i<4;++i){chem_species_init(&species[i],names[i]);members[i]=(ChemPhaseSpecies){&species[i],0,0};}
    ChemSystem system={{"gas",members,kind==1?4u:3u}};
    ChemTerm rs[]={{0,1},{1,1}},ps[]={{1,1},{2,1}};
    ChemReaction reactions[]={{.type=CHEM_ARRHENIUS,.reactant_count=1,.product_count=1,.p={2.15e-11,0,110,300}},
        {.type=CHEM_ARRHENIUS,.reactant_begin=1,.reactant_count=1,.product_begin=1,.product_count=1,.p={3.3e-11,0,55,300}}};
    ChemMechanism m={.species_count=kind==1?4u:3u,.reaction_count=2,.species_names=names,.reactions=reactions,.reactants=rs,.products=ps};
    ChemSolverBuilder b=chem_solver_builder(kind==2 || kind==4?CHEM_SOLVER_BACKWARD_EULER:CHEM_SOLVER_ROSENBROCK);
    b.system=kind==0?NULL:&system;b.mechanism=&m;
    int ok=1;ChemSolver solver={0};
    if(kind==0)ok=chem_solver_build(&b,&solver)==CHEM_BUILD_MISSING_SYSTEM;
    else if(kind==1)ok=chem_solver_build(&b,&solver)==CHEM_BUILD_UNUSED_SPECIES;
    else if(kind==2 || kind==3) {
        ok=chem_solver_build(&b,&solver)==CHEM_BUILD_OK && chem_solver_maximum_cells(&solver)>100000000;
        ChemSolver vector={0};b.group_width=4;
        if(chem_solver_build(&b,&vector)!=CHEM_BUILD_OK)ok=0;
        if(chem_solver_maximum_cells(&vector)!=4 || strcmp(vector.options.system->gas_phase.name,"gas") || vector.options.system->gas_phase.count!=3)ok=0;
        ChemState state={0};ChemStateParameters sp={2,3,0,names,NULL,1e-6,NULL};
        if(chem_state_init(&state,&sp,3,4))ok=0;
        for(size_t cell=0;cell<3;++cell) {
            *chem_matrix_at(&state.variables,cell,0)=1+cell;
            *chem_matrix_at(&state.rate_constants,cell,0)=.1;
            *chem_matrix_at(&state.rate_constants,cell,1)=.2;
        }
        ChemStatistics stats;
        if(chem_solver_solve(&vector,1,&state,NULL,&stats))ok=0;
        for(size_t cell=0;cell<3;++cell) {
            double a=*chem_matrix_at(&state.variables,cell,0), bb=*chem_matrix_at(&state.variables,cell,1), c=*chem_matrix_at(&state.variables,cell,2);
            if(!(a<1+cell && bb>0 && c>0) || fabs(a+bb+c-(1+cell))>1e-12)ok=0;
        }
        chem_state_free(&state);
        chem_solver_free(&vector);
    } else if(kind==4 || kind==5) {
        if(chem_solver_build(&b,&solver)!=CHEM_BUILD_OK)ok=0;
        ChemState state={0};ChemStateParameters p={2,3,0,names,NULL,1e-6,NULL};if(chem_state_init(&state,&p,1,1))ok=0;
        if(kind==5)*chem_state_variable(&state,0,"A")=1;
        ChemStatistics stats;int status=chem_solver_solve(&solver,5,&state,NULL,&stats);
        if(status || stats.final_time!=5 || stats.function_calls!=(unsigned long)(kind==4?2:18) ||
           stats.jacobian_updates!=(unsigned long)(kind==4?2:9) || stats.internal_steps!=(unsigned long)(kind==4?2:9) || stats.solves!=(unsigned long)(kind==4?2:27))ok=0;
        if(kind==4) {
            ChemBackwardEulerParameters options=chem_backward_euler_defaults();options.small=1;options.max_iterations=1;
            status=chem_solver_solve(&solver,5,&state,&options,&stats);
            if(status!=CHEM_SOLVE_ACCEPTING_UNCONVERGED || stats.final_time!=.03125 || stats.function_calls!=6 || stats.jacobian_updates!=6 || stats.internal_steps!=6 || stats.solves!=6)ok=0;
        } else {
            ChemRosenbrockParameters options=chem_rosenbrock_three_stage();options.h_min=15;options.max_steps=6;*chem_state_variable(&state,0,"A")=1;
            status=chem_solver_solve(&solver,5,&state,&options,&stats);
            if(status || stats.final_time!=5 || stats.function_calls!=2 || stats.jacobian_updates!=1 || stats.internal_steps!=1 || stats.solves!=3 || solver.options.rosenbrock.h_min!=15 || solver.options.rosenbrock.max_steps!=6)ok=0;
        }
        chem_state_free(&state);
    } else {
        unsigned char pattern[9]={0};Record record={0};
        ChemExternalModel model={.jacobian_pattern=pattern,.add_forcing=no_op,.add_jacobian=no_op,.context=&record};
        if(kind==6)model.finalize_dense=dense_finalize;else model.finalize_sparse=sparse_finalize;
        b.models=&model;b.model_count=1;b.group_width=kind==7?4:0;
        if(chem_solver_build(&b,&solver)!=CHEM_BUILD_OK || record.calls!=1 || record.order!=3)ok=0;
        if(kind==6) {
            if(record.width!=1)ok=0;
            chem_solver_free(&solver);memset(&record,0,sizeof record);b.group_width=4;
            if(chem_solver_build(&b,&solver)!=CHEM_BUILD_OK || record.calls!=1 || record.width!=4 || record.order!=3)ok=0;
        } else if(!record.sparse_only || record.width)ok=0;
    }
    chem_solver_free(&solver);
    for(size_t i=0;i<4;++i)chem_species_free(&species[i]);
    return ok;
}
int main(void)
{
    const char *names[]={"ThrowsMissingSystem","ThrowsUnusedSpecies","CanBuildBackwardEuler","CanBuildRosenbrock","CanBuildBackwardEulerOverloadedSolverMethod","CanBuildRosenbrockOverloadedSolveMethod","FinalizeReceivesDenseMatrixPolicy","FinalizeSupportsSparseOnlySignature"};
    for(int i=0;i<8;++i)printf("solver_builder:SolverBuilder.%s %d\n",names[i],run(i));
    return 0;
}
