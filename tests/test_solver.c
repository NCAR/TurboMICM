#include "chem.h"
#include "solver.h"
#include <math.h>
#include <stdio.h>

/* A quiescent cell must follow the nonlinear cell's convergence decisions,
 * including any shared step reductions, rather than completing independently. */
static int backward_group(void)
{
    const char *names[]={"A","B"};
    const ChemTerm reactants[]={{0,1},{0,1}},products[]={{1,1}};
    const ChemReaction reaction={CHEM_ARRHENIUS,0,2,0,1,-1,{10,0,0,300}};
    const ChemMechanism m={2,1,0,names,NULL,NULL,NULL,&reaction,reactants,products,0,NULL};
    const ChemConditions c[2]={{300,100000,0,0},{300,100000,0,0}};
    ChemWorkspace *work[2]={chem_workspace_create(&m),chem_workspace_create(&m)};
    ChemWorkspace *single=chem_workspace_create(&m);
    if(!work[0] || !work[1] || !single) return 0;
    ChemBackwardEulerParameters p=chem_backward_euler_defaults();p.relative_tolerance=1e-8;p.absolute_tolerance=1e-12;p.max_iterations=6;
    const double rates[]={10,10};double y[]={0,0,1,0},reference[]={1,0};ChemStatistics gs,ss;
    int ok=!chem_backward_euler_solve(&m,single,&p,c[1],rates+1,1,reference,&ss) &&
        !chem_backward_euler_solve_group(&m,work,2,&p,c,rates,1,y,&gs);
    if(y[0]!=0 || y[1]!=0 || y[2]!=reference[0] || y[3]!=reference[1] ||
       gs.internal_steps!=ss.internal_steps || gs.accepted!=ss.accepted || gs.rejected!=ss.rejected ||
       gs.function_calls!=ss.function_calls || gs.jacobian_updates!=ss.jacobian_updates ||
       gs.decompositions!=ss.decompositions || gs.solves!=ss.solves || gs.final_time!=ss.final_time ||
       !gs.rejected || gs.internal_steps<=2 || fabs(y[2]+2*y[3]-1)>1e-12)ok=0;
    ChemSpecies species[2];ChemPhaseSpecies members[2];
    for(size_t i=0;i<2;++i){chem_species_init(&species[i],names[i]);members[i]=(ChemPhaseSpecies){&species[i],0,0};}
    ChemSystem system={{"gas",members,2}};ChemSolverBuilder builder=chem_solver_builder(CHEM_SOLVER_BACKWARD_EULER);
    builder.system=&system;builder.mechanism=&m;builder.backward_euler=p;
    ChemSolver solver={0};const double atol[]={1e-12,1e-12};ChemState state={0};
    ChemStateParameters sp={1,2,0,names,NULL,1e-8,atol};
    if(chem_solver_build(&builder,&solver)!=CHEM_BUILD_OK || chem_state_init(&state,&sp,2,1))ok=0;
    else {
        for(size_t cell=0;cell<2;++cell){state.conditions[cell]=c[cell];*chem_matrix_at(&state.rate_constants,cell,0)=rates[cell];}
        *chem_matrix_at(&state.variables,1,0)=1;
        ChemStatistics built_stats;
        if(chem_solver_solve(&solver,1,&state,NULL,&built_stats) ||
           *chem_matrix_at(&state.variables,0,0)!=0 || *chem_matrix_at(&state.variables,0,1)!=0 ||
           *chem_matrix_at(&state.variables,1,0)!=reference[0] || *chem_matrix_at(&state.variables,1,1)!=reference[1] ||
           built_stats.internal_steps!=ss.internal_steps || built_stats.accepted!=ss.accepted ||
           built_stats.rejected!=ss.rejected || built_stats.final_time!=ss.final_time)ok=0;
    }
    chem_state_free(&state);chem_solver_free(&solver);
    for(size_t i=0;i<2;++i)chem_species_free(&species[i]);
    chem_workspace_free(work[0]);chem_workspace_free(work[1]);chem_workspace_free(single);return ok;
}
int main(void)
{
    /* A -> B, with k = 2/s; total concentration is conserved. */
    const char *names[] = {"A", "B"};
    const double init[] = {1, 0};
    const ChemTerm reactants[] = {{0, 1}}, products[] = {{1, 1}};
    const ChemReaction reactions[] = {{CHEM_ARRHENIUS, 0, 1, 0, 1, -1, {2, 0, 0, 300, 0}}};
    const ChemMechanism mechanism = {2, 1, 0, names, NULL, init, NULL, reactions, reactants, products, 0, NULL};
    ChemConditions conditions = {300, 100000, 0, 0};
    ChemWorkspace *w = chem_workspace_create(&mechanism);
    double rates[1];
    if (!w || chem_rate_constants(&mechanism, conditions, NULL, rates)) return 1;
    for (int algorithm = 0; algorithm < 2; ++algorithm) for (int method = 0; method < 3; ++method) {
        ChemRosenbrockParameters p = method == 0 ? chem_rosenbrock_two_stage() :
                                     method == 1 ? chem_rosenbrock_three_stage() : chem_rosenbrock_four_stage();
        p.absolute_tolerance = 1e-9; p.relative_tolerance = 1e-6;
        p.linear_solver = (ChemLUAlgorithm)algorithm;
        double y[] = {1, 0};
        ChemStatistics stats;
        if (chem_rosenbrock_solve(&mechanism, w, &p, conditions, rates, .1, y, &stats)) return 1;
        if (fabs(y[0] - exp(-.2)) > 1e-6 || fabs(y[0]+y[1]-1) > 1e-12 ||
            stats.accepted == 0 || stats.solves != (unsigned)p.stages*stats.internal_steps) {
            fprintf(stderr, "unexpected method %d solution %.17g %.17g\n", method, y[0], y[1]); return 2;
        }
        if (!chem_rosenbrock_solve(&mechanism, w, &p, conditions, rates, 0, y, &stats)) return 3;
    }
    /* A one-cell group must follow the same step sequence as the scalar API. */
    ChemWorkspace *other=chem_workspace_create(&mechanism);
    if (!other) return 5;
    ChemWorkspace *group[]={other}; ChemRosenbrockParameters rp=chem_rosenbrock_three_stage();
    double scalar[]={1,0}, grouped[]={1,0}; ChemStatistics ss,gs;
    if (chem_rosenbrock_solve(&mechanism,w,&rp,conditions,rates,.1,scalar,&ss) ||
        chem_rosenbrock_solve_group(&mechanism,group,1,&rp,&conditions,rates,.1,grouped,&gs) ||
        scalar[0]!=grouped[0] || scalar[1]!=grouped[1] || ss.internal_steps!=gs.internal_steps ||
        ss.accepted!=gs.accepted || ss.rejected!=gs.rejected || ss.solves!=gs.solves) return 6;
    chem_workspace_free(other);
    ChemBackwardEulerParameters be = chem_backward_euler_defaults();
    for (int algorithm = 0; algorithm < 2; ++algorithm) {
        be.linear_solver = (ChemLUAlgorithm)algorithm;
        double backward_y[] = {1, 0};
        ChemStatistics backward_stats;
        if (chem_backward_euler_solve(&mechanism, w, &be, conditions, rates, .1, backward_y, &backward_stats) ||
            fabs(backward_y[0] - 1.0/1.2) > 1e-14 || fabs(backward_y[0]+backward_y[1]-1) > 1e-14) return 4;
    }
    if (!backward_group()) return 7;
    chem_workspace_free(w);
    puts("Rosenbrock and backward Euler decay and conservation: PASS");
    return 0;
}
