#include "test.h"
#include "chem_solver.h"
static chem_status builder_fixture(chem_solver *solver,const chem_solver_options *options,int unused,chem_error *error)
{
    chem_phase_species ps[4]={0}; const char *names[]={"A","B","C","D"}; size_t n=unused ? 4 : 3;
    chem_status status=CHEM_OK; chem_phase phase={0}; chem_system system={0}; chem_process processes[2]={0};
    for (size_t i=0;i<n && status==CHEM_OK;++i) status=chem_species_info_create(&ps[i].info,names[i],0,NULL);
    if (status==CHEM_OK) status=chem_phase_create(&phase,"gas",n,ps);
    if (status==CHEM_OK) status=chem_system_create(&system,&phase);
    for (size_t i=0;i<2 && status==CHEM_OK;++i) {
        chem_rate rate=chem_rate_default(CHEM_ARRHENIUS); rate.p.arrhenius.A=i ? 3.3e-11 : 2.15e-11; rate.p.arrhenius.C=i ? 55 : 110;
        status=chem_process_create(processes+i,&phase,1,names+i,1,&(chem_named_product){names[i+1],1},rate);
    }
    if (status==CHEM_OK) status=chem_solver_create(solver,&system,2,processes,options,error);
    for (size_t i=0;i<4;++i) chem_phase_species_destroy(ps+i);
    for (size_t i=0;i<2;++i) chem_process_destroy(processes+i);
    chem_phase_destroy(&phase); chem_system_destroy(&system); return status;
}
TEST(solver_builder,SolverBuilder,ThrowsMissingSystem)
{
    chem_solver s={0}; chem_solver_options options=chem_solver_options_default(CHEM_SOLVER_BACKWARD_EULER);
    CHECK(chem_solver_create(&s,NULL,0,NULL,&options,NULL)==CHEM_INVALID); chem_solver_destroy(&s);
}
TEST(solver_builder,SolverBuilder,ThrowsUnusedSpecies)
{
    chem_solver s={0}; chem_solver_options options=chem_solver_options_default(CHEM_SOLVER_ROSENBROCK); options.reject_unused=1;
    chem_error error={0}; CHECK(builder_fixture(&s,&options,1,&error)==CHEM_INVALID);
    CHECK(error.code==CHEM_INVALID); CHECK(strstr(error.message,"Unused species: D")!=NULL); chem_solver_destroy(&s);
}
static void builder_test(chem_solver_method method)
{
    for (size_t width=0;width<=4;width+=4) {
        chem_solver s={0}; chem_solver_options options=chem_solver_options_default(method); options.layout.vector_length=width;
        REQUIRE(builder_fixture(&s,&options,0,NULL)==CHEM_OK);
        CHECK(!strcmp(s.chemistry.system.gas.name,"gas")); CHECK(s.chemistry.system.gas.n_species==3);
        if (width) CHECK(chem_solver_max_cells(&s)==4); else CHECK(chem_solver_max_cells(&s)>100000000);
        chem_state state={0}; REQUIRE(chem_solver_state(&s,&state,4)==CHEM_OK); chem_state_destroy(&state);
        REQUIRE(chem_solver_state(&s,&state,5)==CHEM_OK); chem_state_destroy(&state);
        chem_solver_destroy(&s);
    }
}
TEST(solver_builder,SolverBuilder,CanBuildBackwardEuler) { builder_test(CHEM_SOLVER_BACKWARD_EULER); }
TEST(solver_builder,SolverBuilder,CanBuildRosenbrock) { builder_test(CHEM_SOLVER_ROSENBROCK); }
static void check_stats(chem_result r,double final,uint64_t functions,uint64_t jacobians,uint64_t steps,uint64_t solves)
{
    NEAR(r.stats.final_time,final,0); CHECK(r.stats.function_calls==functions); CHECK(r.stats.jacobian_updates==jacobians);
    CHECK(r.stats.internal_steps==steps); CHECK(r.stats.solves==solves);
}
TEST(solver_builder,SolverBuilder,CanBuildBackwardEulerOverloadedSolverMethod)
{
    chem_solver solver={0}; chem_solver_options options=chem_solver_options_default(CHEM_SOLVER_BACKWARD_EULER);
    REQUIRE(builder_fixture(&solver,&options,0,NULL)==CHEM_OK); chem_state state={0}; REQUIRE(chem_solver_state(&solver,&state,1)==CHEM_OK);
    chem_be_parameters p=chem_be_parameters_default(); check_stats(chem_solver_solve_be(&solver,&state,5,&p),5,2,2,2,2);
    p.small=1; p.max_steps=1; check_stats(chem_solver_solve_be(&solver,&state,5,&p),.03125,6,6,6,6);
    chem_state_destroy(&state); chem_solver_destroy(&solver);
}
TEST(solver_builder,SolverBuilder,CanBuildRosenbrockOverloadedSolveMethod)
{
    chem_solver solver={0}; chem_solver_options options=chem_solver_options_default(CHEM_SOLVER_ROSENBROCK);
    REQUIRE(builder_fixture(&solver,&options,0,NULL)==CHEM_OK); chem_state state={0}; REQUIRE(chem_solver_state(&solver,&state,1)==CHEM_OK);
    state.y[0]=1; check_stats(chem_solver_solve(&solver,&state,5),5,18,9,9,27);
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3); p.h_min=15; p.max_steps=6; state.y[0]=1; state.y[1]=state.y[2]=0;
    check_stats(chem_solver_solve_ros(&solver,&state,5,&p,NULL,NULL),5,2,1,1,3);
    NEAR(solver.options.rosenbrock.h_min,15,0); CHECK(solver.options.rosenbrock.max_steps==6);
    chem_state_destroy(&state); chem_solver_destroy(&solver);
}
typedef struct { int calls,dense_seen; size_t width; chem_sparse_order order; } finalize_record;
static chem_status finalize_dense(const chem_model *model,const chem_solver_layout *layout,void *context)
{
    finalize_record *record=context; ++record->calls; record->dense_seen=1; record->width=layout->vector_length; record->order=layout->sparse_order;
    return chem_species_index(model,"A")==CHEM_NONE ? CHEM_NOT_FOUND : CHEM_OK;
}
static chem_status finalize_sparse(const chem_model *model,chem_sparse_order order,size_t width,void *context)
{
    finalize_record *record=context; ++record->calls; record->width=width; record->order=order;
    return model->lu.nnz ? CHEM_OK : CHEM_INVALID;
}
TEST(solver_builder,SolverBuilder,FinalizeReceivesDenseMatrixPolicy)
{
    for (size_t width=0;width<=4;width+=4) {
        finalize_record record={0}; chem_setup_callback setup={.finalize=finalize_dense,.context=&record};
        chem_solver_options options=chem_solver_options_default(CHEM_SOLVER_ROSENBROCK); options.layout.vector_length=width; options.n_setups=1; options.setups=&setup;
        chem_solver solver={0}; REQUIRE(builder_fixture(&solver,&options,0,NULL)==CHEM_OK);
        CHECK(record.calls==1); CHECK(record.dense_seen==1); CHECK(record.width==width); CHECK(record.order==CHEM_CSR); chem_solver_destroy(&solver);
    }
}
TEST(solver_builder,SolverBuilder,FinalizeSupportsSparseOnlySignature)
{
    finalize_record record={0}; chem_setup_callback setup={.finalize_sparse=finalize_sparse,.context=&record};
    chem_solver_options options=chem_solver_options_default(CHEM_SOLVER_ROSENBROCK); options.layout.vector_length=4; options.n_setups=1; options.setups=&setup;
    chem_solver solver={0}; REQUIRE(builder_fixture(&solver,&options,0,NULL)==CHEM_OK);
    CHECK(record.calls==1); CHECK(record.dense_seen==0); CHECK(record.width==4); CHECK(record.order==CHEM_CSR); chem_solver_destroy(&solver);
}
