/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_solver.h"
#include "chem_version.h"
#include <omp.h>

static double lambda_rate(const ChemConditions *c,void *ctx)
{return *(const double *)ctx*c->temperature;}
static void override_lambda(void)
{
    ChemSpecies species[]={{.name="A"},{.name="B"}};size_t reactant=0;ChemProduct product={1,1};double original=.001,replacement=.002;
    ChemReaction reaction={.reactant_count=1,.product_count=1,.reactants=&reactant,.products=&product,.rate=chem_rate_default(CHEM_CALLBACK)};
    reaction.rate.label="lambda_rc";reaction.rate.data.callback.function=lambda_rate;reaction.rate.data.callback.context=&original;
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(2,species,1,&reaction,&m)==CHEM_OK);
    CHECK(chem_set_rate_callback(m,"lambda_rc",lambda_rate,&replacement)==CHEM_OK);
    ChemState *state=NULL;ChemSolver *solver=NULL;CHECK(chem_state_create(m,1,&state)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    state->absolute_tolerance[0]=state->absolute_tolerance[1]=1e-6;state->values[0]=1;state->values[1]=0;
    state->conditions[0]=(ChemConditions){300,101253.3,1e6};chem_update_rates(state);NEAR(state->rates[0],.6,1e-12);
    ChemSolverStats stats;ChemBackwardEulerParameters p=chem_backward_euler_parameters();
    CHECK(chem_solve_backward_euler(solver,state,1,&p,&stats)==CHEM_OK);CHECK(stats.final_time==1);
    chem_solver_free(solver);chem_state_free(state);chem_mechanism_free(m);
}
static void threaded_states(void)
{
    ChemSpecies species[]={{.name="A"},{.name="B"},{.name="C"}};
    size_t lhs0[]={0},lhs1[]={1,1},lhs2[]={1,2};
    ChemProduct rhs0[]={{1,1}},rhs1[]={{1,1},{2,1}},rhs2[]={{0,1},{2,1}};
    ChemReaction reactions[]={{1,1,lhs0,rhs0,{0}},{2,2,lhs1,rhs1,{0}},{2,2,lhs2,rhs2,{0}}};
    const char *labels[]={"r1","r2","r3"};double rates[]={.04,3e7,1e4};
    for(size_t i=0;i<3;++i){reactions[i].rate=chem_rate_default(CHEM_USER);reactions[i].rate.label=labels[i];}
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(3,species,3,reactions,&m)==CHEM_OK);
    ChemState *states[8]={0};ChemSolver *work[8]={0};int failures[8]={0};double results[8][3]={{0}};
    for(size_t t=0;t<8;++t){
        CHECK(chem_state_create(m,1,&states[t])==CHEM_OK);CHECK(chem_solver_create(m,1,&work[t])==CHEM_OK);
        states[t]->values[0]=1;states[t]->conditions[0]=(ChemConditions){272.5,101253.3,1e6};
        for(size_t i=0;i<3;++i)CHECK(chem_state_set_custom(states[t],labels[i],rates+i,1)==CHEM_OK);
    }
    ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS3);
    /* Immutable mechanism shared; each state has a separate solver workspace. */
    #pragma omp parallel for num_threads(8) schedule(static,1)
    for(int t=0;t<8;++t){
        for(size_t step=0;step<10;++step){
            chem_update_rates(states[t]);double elapsed=0;
            while(elapsed<200){
                ChemSolverStats stats;ChemStatus status=chem_solve(work[t],states[t],200-elapsed,&p,&stats);
                if(status!=CHEM_OK||stats.final_time<=0){failures[t]=1;break;}
                elapsed=stats.final_time;
            }
        }
        for(size_t j=0;j<3;++j)results[t][j]=states[t]->values[j];
    }
    for(size_t t=0;t<8;++t){
        CHECK(!failures[t]);for(size_t j=0;j<3;++j)CHECK(results[0][j]==results[t][j]);
        chem_solver_free(work[t]);chem_state_free(states[t]);
    }
    chem_mechanism_free(m);
}
void test_runtime(void)
{
    RUN("solver_lambda_rate_constant:Solver.GetLambdaRateConstantByNameCanOverrideLambda",override_lambda());
    RUN("openmp:OpenMP.OneSolverManyStates",threaded_states());
    RUN("version:Version.FullVersion",CHECK(!strcmp(chem_micm_version(),"3.14.0")));
    RUN("version:Version.VersionMajor",CHECK(chem_micm_version_major()==3));
    RUN("version:Version.VersionMinor",CHECK(chem_micm_version_minor()==14));
    RUN("version:Version.VersionPatch",CHECK(chem_micm_version_patch()==0));
    RUN("version:Version.VersionTweak",CHECK(chem_micm_version_tweak()==0));
}
