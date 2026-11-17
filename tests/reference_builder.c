#include "chem.h"

#include <stdio.h>

static const ChemSpecies species[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0},{"D",0,1e-3,0,0}};
static const int reactants[]={0,1};
static const ChemProduct products[]={{1,1},{2,1}};
static const ChemReaction reactions[]={
    {0,1,0,1,0,{.kind=CHEM_ARRHENIUS,.a=2.15e-11,.c=110,.d=300}},
    {1,1,1,1,0,{.kind=CHEM_ARRHENIUS,.a=3.3e-11,.c=55,.d=300}}};
static const ChemMechanism mechanism={3,2,0,species,reactions,reactants,products,NULL,NULL,298,101325};

static int missing(void)
{
    ChemState *state=NULL;
    return chem_state_build(&state,NULL,1,0)==CHEM_INVALID && state==NULL;
}
static int unused(void)
{
    ChemMechanism m=mechanism; m.species_count=4; ChemState *state=NULL;
    int ok=chem_state_build(&state,&m,1,0)==CHEM_UNUSED_SPECIES && state==NULL;
    ok &= chem_state_build(&state,&m,1,1)==CHEM_OK && state && chem_cell_count(state)==1;
    chem_state_free(state); return ok;
}
static int build(int backward)
{
    int ok=1;
    for (size_t cells=1;cells<=7;cells+=3) {
        ChemState *state=NULL; ChemStats stats;
        ok &= chem_state_build(&state,&mechanism,cells,0)==CHEM_OK && state && chem_cell_count(state)==cells;
        if (!state) return 0;
        if (backward) { ChemBackwardEuler p=chem_backward_euler_defaults(); ok &= chem_solve_backward_euler(state,&p,5,&stats)==CHEM_OK; }
        else { ChemRosenbrock p=chem_rosenbrock_three_stage(); ok &= chem_solve(state,&p,5,&stats)==CHEM_OK; }
        ok &= stats.final_time==5;
        chem_state_free(state);
    }
    return ok;
}
static int backward_options(void)
{
    ChemState *s=NULL; if (chem_state_build(&s,&mechanism,1,0)!=CHEM_OK) return 0;
    ChemBackwardEuler p=chem_backward_euler_defaults(); ChemStats stats;
    int ok=chem_solve_backward_euler(s,&p,5,&stats)==CHEM_OK && stats.final_time==5 &&
        stats.function_calls==2 && stats.jacobian_updates==2 && stats.internal_steps==2 && stats.solves==2;
    p.small=1; p.max_iterations=1;
    ok &= chem_solve_backward_euler(s,&p,5,&stats)==CHEM_UNCONVERGED && stats.final_time==.03125 &&
          stats.function_calls==6 && stats.jacobian_updates==6 && stats.internal_steps==6 && stats.solves==6;
    chem_state_free(s); return ok;
}
static int rosenbrock_options(void)
{
    ChemState *s=NULL; if (chem_state_build(&s,&mechanism,1,0)!=CHEM_OK) return 0;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats stats; chem_species_values(s,0)[0]=1;
    int ok=chem_solve(s,&p,5,&stats)==CHEM_OK && stats.final_time==5 && stats.function_calls==18 &&
           stats.jacobian_updates==9 && stats.internal_steps==9 && stats.solves==27;
    p.h_min=15; p.max_steps=6; chem_species_values(s,0)[0]=1;
    ok &= chem_solve(s,&p,5,&stats)==CHEM_OK && stats.final_time==5 && stats.function_calls==2 &&
          stats.jacobian_updates==1 && stats.internal_steps==1 && stats.solves==3 && p.h_min==15 && p.max_steps==6;
    if (!ok) fprintf(stderr,"Rosenbrock options: final=%g calls=%lu jac=%lu steps=%lu solves=%lu\n",stats.final_time,
                     stats.function_calls,stats.jacobian_updates,stats.internal_steps,stats.solves);
    chem_state_free(s); return ok;
}
typedef struct { unsigned calls; size_t species,cells,nonzeros; int diagonal; } Record;
static int finalize(const ChemMechanism *m,const ChemLinearSystem *p,void *context)
{
    Record *record=context; ++record->calls; record->species=m->species_count; record->cells=p->cells;
    record->nonzeros=p->nonzeros; record->diagonal=p->lookup[0];
    return CHEM_OK;
}
static int finalize_case(int multiple)
{
    int ok=1;
    for (size_t cells=1;cells<=(multiple ? 4 : 1);cells+=3) {
        Record record={0}; ChemExternalModel model={.finalize=finalize,.context=&record};
        ChemState *state=chem_state_create_extended(&mechanism,NULL,0,&model,1,cells);
        ok &= state && record.calls==1 && record.species==3 && record.cells==cells && record.nonzeros>=3 && record.diagonal>=0;
        ChemState *second=chem_state_create_like(state,cells+1);
        ok &= second && record.calls==1; chem_state_free(state); chem_state_free(second);
    }
    return ok;
}
static int failures;
static void report(const char *name,int ok)
{
    printf("solver_builder:SolverBuilder.%s %s %s\n",name,ok ? "PASS" : "FAIL",
           ok ? "reference construction behavior matched" : "construction behavior differed"); failures+=!ok;
}
int main(void)
{
    report("ThrowsMissingSystem",missing());
    report("ThrowsUnusedSpecies",unused());
    report("CanBuildBackwardEuler",build(1));
    report("CanBuildRosenbrock",build(0));
    report("CanBuildBackwardEulerOverloadedSolverMethod",backward_options());
    report("CanBuildRosenbrockOverloadedSolveMethod",rosenbrock_options());
    report("FinalizeReceivesDenseMatrixPolicy",finalize_case(1));
    report("FinalizeSupportsSparseOnlySignature",finalize_case(0));
    return failures ? 1 : 0;
}
