/* SPDX-License-Identifier: Apache-2.0
 * MICM regression_policy.hpp: the original 3600-interval wall-loss fixture. */
#include "test.h"
#include "chem_solver.h"
#include <stdlib.h>
static void flowtube(int euler)
{
    ChemSpecies species[]={{.name="SOA1"},{.name="SOA2"},{.name="a-pinene"},{.name="O3"}};
    size_t r0=0,r1=1,r2[]={2,3};ChemProduct products[]={{0,.18},{1,.09}};
    ChemRate rates[]={chem_rate_default(CHEM_USER),chem_rate_default(CHEM_USER),chem_rate_default(CHEM_ARRHENIUS)};
    rates[0].label="r1";rates[1].label="r2";rates[2].data.arrhenius.A=8.8e-17*(1e-6*6.02214076e23);
    ChemReaction reactions[]={{1,0,&r0,NULL,rates[0]},{1,0,&r1,NULL,rates[1]},{2,2,r2,products,rates[2]}};
    ChemMechanism *m=NULL;ChemState *s=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_create(4,species,3,reactions,&m)==CHEM_OK);
    CHECK(chem_state_create(m,1,&s)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    s->values[2]=8e-8;s->values[3]=.00002;
    CHECK(chem_state_fill_custom(s,"r1",.01)==CHEM_OK);CHECK(chem_state_fill_custom(s,"r2",.05)==CHEM_OK);
    s->conditions[0]=(ChemConditions){298.15,101325,0};chem_ideal_air_density(s->conditions);chem_update_rates(s);
    double first[10][4];memcpy(first[0],s->values,sizeof(first[0]));
    ChemRosenbrockParameters ros=chem_rosenbrock_parameters(CHEM_ROS2);
    ChemBackwardEulerParameters be=chem_backward_euler_parameters();ChemSolverStats stats;
    for (size_t i=0;i<3600;++i) {
        double actual=0;
        while (actual<1) {
            ChemStatus status=euler?chem_solve_backward_euler(solver,s,1-actual,&be,&stats):chem_solve(solver,s,1-actual,&ros,&stats);
            CHECK(status==CHEM_OK);CHECK(stats.final_time>0);
            if (status || !(stats.final_time>0)) goto release;
            actual+=stats.final_time;
        }
        if (i<9) memcpy(first[i+1],s->values,sizeof(first[i+1]));
    }
    const char *path=euler?"tests/data/backward_euler_wall_loss.csv":"tests/data/rosenbrock_wall_loss.csv";
    FILE *file=fopen(path,"r");CHECK(file!=NULL);
    if (file) {
        char line[512];CHECK(fgets(line,sizeof(line),file)!=NULL);
        line[strcspn(line,"\r\n")]=0;CHECK(!strcmp(line,"time,SOA1,SOA2,a-pinene,O3"));size_t rows=0;
        while (fgets(line,sizeof(line),file)) {
            double time,expected[4];char extra;
            int fields=sscanf(line,"%lf,%lf,%lf,%lf,%lf %c",&time,expected,expected+1,expected+2,expected+3,&extra);
            CHECK(fields==5);
            if (fields==5 && rows<10) for (size_t v=0;v<4;++v) NEAR(first[rows][v],expected[v],1e-10);
            ++rows;
        }
        CHECK(rows==3601);CHECK(!ferror(file));fclose(file);
    }
release:
    chem_solver_free(solver);chem_state_free(s);chem_mechanism_free(m);
}
void test_flowtube(void)
{
    RUN("regression:RosenbrockRegression.FlowTube",flowtube(0));
    RUN("regression:BackwardEulerRegression.FlowTube",flowtube(1));
}
