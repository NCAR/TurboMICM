#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "cpu_solver.h"

int main(void)
{
    ChemSpecies species[3]={{0}}; ChemPhaseSpecies entries[3]={{0}};
    const char *names[]={"A","B","C"};
    for (size_t i=0;i<3;++i) {
        assert(chem_species_init(species+i,names[i])==CHEM_SPECIES_OK);
        assert(chem_phase_species_init(entries+i,species+i)==CHEM_SPECIES_OK);
    }
    ChemPhase gas={0}; ChemSystem system={0};
    assert(chem_phase_init(&gas,"gas",3,entries)==CHEM_SPECIES_OK);
    assert(chem_system_init(&system,&gas)==CHEM_SPECIES_OK);
    ChemProcess processes[2]={{0}}; double yield=1;
    for (size_t i=0;i<2;++i) {
        const ChemSpecies *reactants[]={species+i},*products[]={species+i+1};
        ChemProcessConfig config={1,1,reactants,products,&yield,&gas,
            {CHEM_RATE_ARRHENIUS,{i ? 3.3e-11 : 2.15e-11,0,i ? 55 : 110,300,0},0}};
        assert(chem_process_init(processes+i,&config)==CHEM_SPECIES_OK);
    }
    for (int kind=0;kind<2;++kind) {
        ChemCpuSolver scalar={0},vector={0};
        assert(chem_cpu_solver_init(&scalar,&system,2,processes,0,(ChemSolverKind)kind,0,0,NULL)==CHEM_BUILD_OK);
        assert(chem_cpu_solver_init(&vector,&system,2,processes,0,(ChemSolverKind)kind,4,0,NULL)==CHEM_BUILD_OK);
        assert(!strcmp(vector.built.system.gas.name,"gas") && vector.built.system.gas.count==3);
        assert(chem_cpu_solver_maximum_cells(&scalar)>100000000 && chem_cpu_solver_maximum_cells(&vector)==4);
        ChemState state={0}; assert(chem_cpu_solver_state(&vector,5,&state)==CHEM_STATE_INVALID);
        assert(chem_cpu_solver_state(&scalar,1,&state)==CHEM_STATE_OK);
        ChemSolverStats stats={0};
        if (kind==CHEM_SOLVER_BACKWARD_EULER) {
            assert(chem_cpu_solver_solve(&scalar,&state,5,NULL,NULL,&stats)==CHEM_SOLVE_OK);
            assert(stats.final_time==5 && stats.function_calls==2 && stats.jacobian_updates==2 && stats.steps==2 && stats.solves==2);
            ChemBackwardEulerOptions options=chem_backward_euler_defaults(NULL); options.small=1; options.max_iterations=1;
            assert(chem_cpu_solver_solve(&scalar,&state,5,NULL,&options,&stats)==CHEM_SOLVE_UNCONVERGED);
            assert(stats.final_time==.03125 && stats.function_calls==6 && stats.jacobian_updates==6 && stats.steps==6 && stats.solves==6);
            assert(scalar.backward_euler.small==1 && scalar.backward_euler.max_iterations==1);
            assert(chem_cpu_solver_solve(&scalar,&state,5,NULL,NULL,&stats)==CHEM_SOLVE_UNCONVERGED);
            assert(stats.steps==6 && stats.final_time==.03125);
        } else {
            state.concentrations[0]=1;
            assert(chem_cpu_solver_solve(&scalar,&state,5,NULL,NULL,&stats)==CHEM_SOLVE_OK);
            assert(stats.final_time==5 && stats.function_calls==18 && stats.jacobian_updates==9 && stats.steps==9 && stats.solves==27);
            ChemRosenbrockOptions options=chem_rosenbrock_defaults(NULL); options.h_min=15; options.max_steps=6;
            state.concentrations[0]=1; state.concentrations[1]=state.concentrations[2]=0;
            assert(chem_cpu_solver_solve(&scalar,&state,5,&options,NULL,&stats)==CHEM_SOLVE_OK);
            assert(stats.final_time==5 && stats.function_calls==2 && stats.jacobian_updates==1 && stats.steps==1 && stats.solves==3);
            assert(scalar.rosenbrock.h_min==15 && scalar.rosenbrock.max_steps==6);
            assert(chem_cpu_solver_solve(&scalar,&state,5,NULL,NULL,&stats)==CHEM_SOLVE_OK);
            assert(stats.steps==1 && stats.final_time==5);
        }
        chem_state_destroy(&state);
        assert(chem_cpu_solver_state(&vector,3,&state)==CHEM_STATE_OK);
        for (size_t cell=0;cell<3;++cell) { state.conditions[cell]=(ChemConditions){300,100000,1}; state.concentrations[3*cell]=1+cell; }
        assert(chem_cpu_solver_update_parameters(&vector,&state)==CHEM_STATE_OK);
        assert(state.rates[0]>0 && state.rates[1]>0 && state.rates[2]==state.rates[0]);
        for (size_t cell=0;cell<3;++cell) { state.rates[2*cell]=.1; state.rates[2*cell+1]=.2; }
        assert(chem_cpu_solver_solve(&vector,&state,5,NULL,NULL,&stats)==CHEM_SOLVE_OK);
        for (size_t cell=0;cell<3;++cell) { assert(state.concentrations[3*cell]<(double)(1+cell)); assert(fabs(state.concentrations[3*cell]+state.concentrations[3*cell+1]+state.concentrations[3*cell+2]-(1+cell))<1e-12); }
        chem_state_destroy(&state); chem_cpu_solver_destroy(&scalar); chem_cpu_solver_destroy(&vector);
    }
    chem_system_destroy(&system); chem_phase_destroy(&gas);
    for (size_t i=0;i<2;++i) chem_process_destroy(processes+i);
    for (size_t i=0;i<3;++i) { chem_phase_species_destroy(entries+i); chem_species_destroy(species+i); }
    return 0;
}
