#include "chem.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
    const char *names[] = {"A", "B"};
    const ChemTerm reactants[] = {{0,1}}, products[] = {{1,1}}, terms[] = {{0,1},{1,1}};
    const ChemReaction reaction[] = {{CHEM_ARRHENIUS,0,1,0,1,-1,{.1,0,0,300,0}}};
    const ChemMechanism mechanism = {2,1,0,names,NULL,NULL,NULL,reaction,reactants,products,0,NULL};
    const ChemConstraint constraint = {"total",CHEM_LINEAR_CONSTRAINT,1,terms,NULL,2,0,1,0,0,0};
    ChemWorkspace *work = chem_workspace_create_constrained(&mechanism,&constraint,1);
    ChemConditions conditions = {300,100000,0,0};
    double rates[1];
    if (!work || chem_rate_constants(&mechanism,conditions,NULL,rates)) return 1;
    for (int algorithm=0; algorithm<2; ++algorithm) for (int method=0; method<2; ++method) {
        ChemRosenbrockParameters p = method ? chem_rosenbrock_six_stage_dae() : chem_rosenbrock_four_stage_dae();
        p.linear_solver = (ChemLUAlgorithm)algorithm;
        double y[] = {.3,0};
        ChemStatistics stats;
        int status = chem_rosenbrock_solve(&mechanism,work,&p,conditions,rates,.01,y,&stats);
        if (status || fabs(y[0]-.3*exp(-.001))>1e-7 || fabs(y[0]+y[1]-1)>1e-10) {
            fprintf(stderr,"DAE method %d status %d state %.17g %.17g\n", method,status,y[0],y[1]); return 2;
        }
    }
    chem_workspace_free(work);
    const ChemTerm eqreact[] = {{0,1}}, eqprod[] = {{1,2}};
    const ChemConstraint equilibrium = {"nonlinear",CHEM_EQUILIBRIUM_CONSTRAINT,1,eqreact,eqprod,1,1,2,0,8.31446261815324,298.15};
    work = chem_workspace_create_constrained(&mechanism,&equilibrium,1);
    if (!work) return 3;
    for (int method=0; method<2; ++method) {
        ChemRosenbrockParameters p = method ? chem_rosenbrock_six_stage_dae() : chem_rosenbrock_four_stage_dae();
        double y[] = {.3,.1}; ChemStatistics stats;
        int status = chem_rosenbrock_solve(&mechanism,work,&p,conditions,rates,.001,y,&stats);
        if (status || fabs(2*y[0]-y[1]*y[1])>1e-6) {
            fprintf(stderr,"nonlinear DAE method %d status %d state %.17g %.17g\n",method,status,y[0],y[1]); return 4;
        }
    }
    chem_workspace_free(work);
    puts("DAE initialization and integration: PASS");
    return 0;
}
