#include "test.h"
#include "chem_solver.h"
#include "test_backend.h"
#include "test_random.h"
#include <omp.h>

TEST(integrated_reaction_rates,ChapmanIntegration,CanBuildChapmanSystem)
{
    chem_species species[]={{.name="a"},{.name="b"},{.name="c"},{.name="irr_1"},{.name="irr_2"}};
    size_t reactants[]={0,1}; chem_product products[]={{1,1},{3,1},{2,1},{4,1}};
    chem_rate arr=chem_rate_default(CHEM_ARRHENIUS); arr.p.arrhenius.A=2.15e-11; arr.p.arrhenius.C=110;
    chem_reaction reactions[]={{1,2,reactants,products,arr},{1,2,reactants+1,products+2,{.type=CHEM_USER,.label="r2",.p.user={1}}}};
    chem_mechanism mechanism={.n_species=5,.n_reactions=2,.species=species,.reactions=reactions};
    chem_model model={0}; chem_state state={0}; integration_backend backend={0};
    REQUIRE(chem_model_create(&model,&mechanism,CHEM_LU_DOOLITTLE)==CHEM_OK);
    REQUIRE(chem_state_create(&state,&model,1)==CHEM_OK); REQUIRE(integration_backend_create(&backend,&state,0,CHEM_CSR,0)==CHEM_OK);
    REQUIRE(chem_state_set_custom(&state,"r2",(double[]){1},1)==CHEM_OK);
    state.y[chem_species_index(&model,"a")]=1; state.conditions[0]=(chem_conditions){273,1000,0};
    size_t i1=chem_species_index(&model,"irr_1"),i2=chem_species_index(&model,"irr_2");
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);
    for (size_t step=0;step<100;++step) {
        if (step>50) REQUIRE(chem_state_set_custom(&state,"r2",(double[]){0},1)==CHEM_OK);
        chem_calculate_rates(&state); double r1=state.y[i1],r2=state.y[i2];
        chem_result result=chem_ros_solve(&state,30,&p,NULL,NULL); chem_state_clamp(&state); CHECK(result.status==CHEM_OK);
        CHECK(state.y[i1]>=r1); CHECK(state.y[i2]>=r2);
    }
    integration_backend_destroy(&backend); chem_state_destroy(&state); chem_model_destroy(&model);
}

static chem_status threaded_solver_create(chem_solver *solver)
{
    const char *names[]={"A","B","C"}; chem_phase_species species[3]={0}; chem_phase phase={0}; chem_system system={0};
    chem_process processes[3]={0}; chem_status status=CHEM_OK;
    for (size_t i=0;i<3 && status==CHEM_OK;++i) status=chem_species_info_create(&species[i].info,names[i],0,NULL);
    if (status==CHEM_OK) status=chem_phase_create(&phase,"gas",3,species);
    if (status==CHEM_OK) status=chem_system_create(&system,&phase);
    const char *reactants[][2]={{"A",NULL},{"B","B"},{"B","C"}};
    chem_named_product products[][2]={{{"B",1},{NULL,0}},{{"B",1},{"C",1}},{{"A",1},{"C",1}}};
    const char *labels[]={"r1","r2","r3"};
    for (size_t i=0;i<3 && status==CHEM_OK;++i) {
        chem_rate rate=chem_rate_default(CHEM_USER); rate.label=labels[i];
        status=chem_process_create(processes+i,&phase,i ? 2 : 1,reactants[i],i ? 2 : 1,products[i],rate);
    }
    chem_solver_options options=chem_solver_options_default(CHEM_SOLVER_ROSENBROCK);
    if (status==CHEM_OK) status=chem_solver_create(solver,&system,3,processes,&options,NULL);
    for (size_t i=0;i<3;++i) { chem_phase_species_destroy(species+i); chem_process_destroy(processes+i); }
    chem_phase_destroy(&phase); chem_system_destroy(&system); return status;
}
TEST(openmp,OpenMP,OneSolverManyStates)
{
    chem_solver solver={0}; REQUIRE(threaded_solver_create(&solver)==CHEM_OK);
    double results[8][3]={{0}}; int errors[8]={0},seen[8]={0};
    #pragma omp parallel num_threads(8)
    {
        int thread=omp_get_thread_num(); seen[thread]=1; chem_state state={0};
        if (chem_solver_state(&solver,&state,1)!=CHEM_OK) errors[thread]=1;
        else {
            /* The reference initializes the physical row, after reordering. */
            state.y[0]=1;
            errors[thread]|=chem_state_set_custom(&state,"r1",(double[]){.04},1)!=CHEM_OK;
            errors[thread]|=chem_state_set_custom(&state,"r2",(double[]){3e7},1)!=CHEM_OK;
            errors[thread]|=chem_state_set_custom(&state,"r3",(double[]){1e4},1)!=CHEM_OK;
            state.conditions[0]=(chem_conditions){272.5,101253.3,1e6};
            for (size_t step=0;step<10 && !errors[thread];++step) {
                chem_calculate_rates(&state); double advanced=0;
                while (advanced<200) {
                    chem_result r=chem_solver_solve(&solver,&state,200-advanced);
                    if ((r.status!=CHEM_OK && r.status!=CHEM_MAX_STEPS) || advanced+r.stats.final_time<=advanced) { errors[thread]=1; break; }
                    advanced+=r.stats.final_time;
                }
            }
            memcpy(results[thread],state.y,sizeof(results[thread]));
        }
        chem_state_destroy(&state);
    }
    for (size_t thread=0;thread<8;++thread) {
        CHECK(seen[thread]); CHECK(!errors[thread]);
        for (size_t j=0;j<3;++j) { CHECK(isfinite(results[thread][j])); NEAR(results[thread][j],results[0][j],0); }
    }
    chem_solver_destroy(&solver);
}

static void terminator(size_t cells,size_t width)
{
    chem_species species[]={{.name="Cl2",.absolute_tolerance=1e-20},{.name="Cl",.absolute_tolerance=1e-20}};
    size_t reactants[]={0,1,1}; chem_product products[]={{1,2},{0,1}};
    chem_rate arr=chem_rate_default(CHEM_ARRHENIUS); arr.p.arrhenius.A=1;
    chem_rate photo=chem_rate_default(CHEM_USER); photo.label="toy_k1";
    chem_reaction reactions[]={{1,1,reactants,products,photo},{2,1,reactants+1,products+1,arr}};
    chem_reactions_sort(reactions,2);
    chem_mechanism mechanism={.n_species=2,.n_reactions=2,.species=species,.reactions=reactions};
    chem_model model={0}; chem_state state={0}; integration_backend backend={0};
    REQUIRE(chem_model_create(&model,&mechanism,CHEM_LU_DOOLITTLE)==CHEM_OK);
    REQUIRE(chem_state_create(&state,&model,cells)==CHEM_OK); REQUIRE(integration_backend_create(&backend,&state,width,CHEM_CSR,0)==CHEM_OK);
    state.rtol=1e-8; double initial[8],k1[4]; size_t cl2=chem_species_index(&model,"Cl2"),cl=chem_species_index(&model,"Cl");
    random_state rng; random_init(&rng,1,0);
    for (size_t c=0;c<cells;++c) { initial[cl2*cells+c]=lognormal(&rng,2)*1e-6; initial[cl*cells+c]=lognormal(&rng,2)*1e-10; }
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3); p.max_steps=100000;
    double pi=acos(-1),lat_center=pi*20/180,lon_center=pi*300/180;
    for (size_t step=0;step<(size_t)floor(2*pi/.3);++step) {
        memcpy(state.y,initial,2*cells*sizeof(double));
        for (size_t c=0;c<cells;++c) {
            double lat=pi/180*(c*(90.0/cells));
            k1[c]=fmax(0,sin(lat)*sin(lat_center)+cos(lat)*cos(lat_center)*cos(step*.3-lon_center));
            state.conditions[c]=(chem_conditions){298,101300,42};
        }
        REQUIRE(chem_state_assign_custom(&state,"toy_k1",k1,cells)==CHEM_OK); chem_calculate_rates(&state);
        chem_result result=chem_ros_solve(&state,30,&p,NULL,NULL); chem_state_clamp(&state); CHECK(result.status==CHEM_OK);
        for (size_t c=0;c<cells;++c) {
            double r=k1[c]/4,cl_i=initial[cl*cells+c],cl2_i=initial[cl2*cells+c],cly=cl_i+2*cl2_i;
            double det=sqrt(r*r+2*r*cly),e=exp(-4*det*30),l=det*30>1e-16 ? (1-e)/det/30 : 4;
            double cl_f=-l*(cl_i-det+r)*(cl_i+det+r)/(1+e+30*l*(cl_i+r));
            double cl_ref=cl_i+30*cl_f,cl2_ref=cl2_i+30*(-cl_f/2);
            NEAR(state.y[cl*cells+c],cl_ref,fabs(cl_ref)*1e-6+1e-14);
            NEAR(state.y[cl2*cells+c],cl2_ref,fabs(cl2_ref)*1e-6+1e-14);
        }
    }
    integration_backend_destroy(&backend); chem_state_destroy(&state); chem_model_destroy(&model);
}
TEST(terminator,RosenbrockSolver,Terminator)
{ for (size_t repeat=0;repeat<2;++repeat) for (size_t cells=1;cells<=4;++cells) terminator(cells,0); }
TEST(terminator,RosenbrockSolver,VectorTerminator)
{ for (size_t cells=1;cells<=4;++cells) terminator(cells,cells); }
