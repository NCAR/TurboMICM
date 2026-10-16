#include "test.h"
#include "chem_process.h"
#include "test_backend.h"
static void flow_tube(int backward_euler)
{
    chem_species species[]={{.name="SOA1"},{.name="SOA2"},{.name="a-pinene"},{.name="O3"}};
    size_t reactants[]={0,1,2,3}; chem_product products[]={{0,.18},{1,.09}};
    chem_rate arr=chem_rate_default(CHEM_ARRHENIUS),u1=chem_rate_default(CHEM_USER),u2=u1;
    arr.p.arrhenius.A=8.8e-17*(1e-6*6.02214076e23); u1.label="r1"; u2.label="r2";
    chem_reaction reactions[]={{1,0,reactants,NULL,u1},{1,0,reactants+1,NULL,u2},{2,2,reactants+2,products,arr}};
    chem_reactions_sort(reactions,3);
    chem_mechanism mechanism={.n_species=4,.n_reactions=3,.species=species,.reactions=reactions,.preserve_order=1};
    chem_model model={0}; chem_state state={0}; integration_backend backend={0};
    REQUIRE(chem_model_create(&model,&mechanism,CHEM_LU_DOOLITTLE)==CHEM_OK);
    REQUIRE(chem_state_create(&state,&model,1)==CHEM_OK); REQUIRE(integration_backend_create(&backend,&state,0,CHEM_CSR,0)==CHEM_OK);
    memcpy(state.y,(double[]){0,0,8e-8,.00002},4*sizeof(double)); memcpy(state.custom,(double[]){.01,.05},2*sizeof(double));
    state.conditions[0]=(chem_conditions){298.15,101325,101325/(CHEM_GAS_CONSTANT*298.15)}; chem_calculate_rates(&state);
    FILE *file=fopen(backward_euler ? "tests/reference/backward_euler_wall_loss.csv" : "tests/reference/rosenbrock_wall_loss.csv","r"); REQUIRE(file);
    char line[256]; REQUIRE(fgets(line,sizeof(line),file)); CHECK(!strcmp(line,"time,SOA1,SOA2,a-pinene,O3\n"));
    chem_ros_parameters rp=chem_ros_parameters_default(CHEM_ROS2); chem_be_parameters bp=chem_be_parameters_default();
    for (size_t step=0;step<=3600;++step) {
        if (step) {
            double advanced=0;
            while (advanced<1) {
                chem_result r=backward_euler ? chem_be_solve(&state,1-advanced,&bp) : chem_ros_solve(&state,1-advanced,&rp,NULL,NULL);
                chem_state_clamp(&state); REQUIRE(advanced+r.stats.final_time>advanced); advanced+=r.stats.final_time;
            }
        }
        double time,expected[4]; REQUIRE(fgets(line,sizeof(line),file));
        CHECK(sscanf(line,"%lf,%lf,%lf,%lf,%lf",&time,expected,expected+1,expected+2,expected+3)==5); NEAR(time,step,0);
        /* The source checks its first ten rows. Check the entire supplied file
         * at that same absolute tolerance after executing all 3600 steps. */
        for (size_t j=0;j<4;++j) NEAR(state.y[j],expected[j],1e-10);
    }
    CHECK(fgets(line,sizeof(line),file)==NULL); fclose(file);
    integration_backend_destroy(&backend); chem_state_destroy(&state); chem_model_destroy(&model);
}
TEST(regression,RosenbrockRegression,FlowTube) { flow_tube(0); }
TEST(regression,BackwardEulerRegression,FlowTube) { flow_tube(1); }
