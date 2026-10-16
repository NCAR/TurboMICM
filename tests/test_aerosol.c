#include "test.h"
#include "chem_external.h"
#include "chem_solver.h"
#include "test_backend.h"

typedef struct {
    const char *reactants[2],*products[2],*parameters[2];
    double rates[2]; size_t r[2],p[2],param[2],rr[2],pr[2];
    int temperature_dependent,finalized;
} aerosol_context;
static chem_status aerosol_finalize(const chem_model *m,void *context)
{
    aerosol_context *a=context; ++a->finalized;
    for (size_t i=0;i<2;++i) {
        a->r[i]=chem_species_index(m,a->reactants[i]); a->p[i]=chem_species_index(m,a->products[i]);
        a->param[i]=a->parameters[i] ? chem_custom_index(m,a->parameters[i]) : CHEM_NONE;
        if (a->r[i]==CHEM_NONE || a->p[i]==CHEM_NONE || (a->parameters[i] && a->param[i]==CHEM_NONE)) return CHEM_NOT_FOUND;
        a->rr[i]=m->lu.map[a->r[i]*m->n_species+a->r[i]]; a->pr[i]=m->lu.map[a->p[i]*m->n_species+a->r[i]];
        if (a->rr[i]==CHEM_NONE || a->pr[i]==CHEM_NONE) return CHEM_NOT_FOUND;
    }
    return CHEM_OK;
}
static void aerosol_update(chem_state *s,void *context)
{
    const aerosol_context *a=context;
    if (a->temperature_dependent) for (size_t c=0;c<s->cells;++c) s->custom[a->param[1]*s->cells+c]=.005*s->conditions[c].temperature;
}
static void aerosol_forcing(const chem_state *s,const double *y,double *out,void *context)
{
    const aerosol_context *a=context;
    for (size_t i=0;i<2;++i) for (size_t c=0;c<s->cells;++c) {
        double k=a->param[i]==CHEM_NONE ? a->rates[i] : s->custom[a->param[i]*s->cells+c];
        double rate=k*y[a->r[i]*s->cells+c]; out[a->r[i]*s->cells+c]-=rate; out[a->p[i]*s->cells+c]+=rate;
    }
}
static void aerosol_jacobian(const chem_state *s,const double *y,double *out,void *context)
{
    (void)y; const aerosol_context *a=context;
    for (size_t i=0;i<2;++i) for (size_t c=0;c<s->cells;++c) {
        double k=a->param[i]==CHEM_NONE ? a->rates[i] : s->custom[a->param[i]*s->cells+c];
        out[a->rr[i]*s->cells+c]+=k; out[a->pr[i]*s->cells+c]-=k;
    }
}
static const char *const aerosol1_names[]={"STUB1.MODE1.QUUX.BAZ","STUB1.MODE1.QUUX.QUX","STUB1.MODE2.CORGE.BAZ","STUB1.MODE2.CORGE.FO2","STUB1.MODE2.CORGE.QUX","STUB1.MODE2.QUUX.BAZ","STUB1.MODE2.QUUX.QUX"};
static const char *const aerosol2_names[]={"STUB2.MODE1.NUMBER","STUB2.MODE1.QUUX.BAZ","STUB2.MODE1.QUUX.QUX","STUB2.MODE2.CORGE.BAZ","STUB2.MODE2.CORGE.FO2","STUB2.MODE2.CORGE.QUX","STUB2.MODE2.NUMBER","STUB2.MODE3.CORGE.BAZ","STUB2.MODE3.CORGE.FO2","STUB2.MODE3.CORGE.QUX","STUB2.MODE3.NUMBER","STUB2.MODE3.QUUX.BAZ","STUB2.MODE3.QUUX.QUX"};
static const char *const aerosol_parameters[]={"STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT","STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT"};
static chem_status aerosol_create(chem_solver *solver,const integration_variant *variant)
{
    chem_species species1[7]={0},species2[13]={0};
    for (size_t i=0;i<7;++i) species1[i].name=aerosol1_names[i];
    for (size_t i=0;i<13;++i) species2[i].name=aerosol2_names[i];
    aerosol_context contexts[]={
        {.reactants={"FO2","STUB1.MODE1.QUUX.BAZ"},.products={"STUB1.MODE2.CORGE.FO2","STUB1.MODE2.QUUX.BAZ"},.rates={1e-3,2e-3}},
        {.reactants={"STUB2.MODE2.CORGE.FO2","STUB2.MODE3.QUUX.BAZ"},.products={"STUB2.MODE2.CORGE.BAZ","STUB2.MODE3.QUUX.QUX"},.parameters={aerosol_parameters[0],aerosol_parameters[1]},.temperature_dependent=1}};
    chem_named_element elements[2][4];
    for (size_t i=0;i<2;++i) for (size_t j=0;j<2;++j) {
        elements[i][2*j]=(chem_named_element){contexts[i].reactants[j],contexts[i].reactants[j]};
        elements[i][2*j+1]=(chem_named_element){contexts[i].products[j],contexts[i].reactants[j]};
    }
    chem_named_value parameters[]={{aerosol_parameters[0],0},{aerosol_parameters[1],0}}; chem_external_model external[2]={0};
    for (size_t i=0;i<2;++i) external[i]=(chem_external_model){.n_species=i ? 13 : 7,.species=i ? species2 : species1,
        .n_parameters=i ? 2 : 0,.parameters=i ? parameters : NULL,.n_process_elements=4,.process_elements=elements[i],
        .forcing=aerosol_forcing,.jacobian=aerosol_jacobian,.update=aerosol_update,.finalize=aerosol_finalize,
        .context=contexts+i,.context_size=sizeof(aerosol_context)};
    chem_phase_species gas_species[2]={0}; chem_phase gas={0}; chem_system system={0};
    chem_status status=chem_species_info_create(&gas_species[0].info,"FO2",0,NULL);
    if (status==CHEM_OK) status=chem_species_info_create(&gas_species[1].info,"BAR",0,NULL);
    if (status==CHEM_OK) status=chem_phase_create(&gas,"GAS",2,gas_species);
    if (status==CHEM_OK) status=chem_system_create(&system,&gas);
    chem_solver_options options=chem_solver_options_default(variant->method<0 ? CHEM_SOLVER_BACKWARD_EULER : CHEM_SOLVER_ROSENBROCK);
    if (variant->method>=0) options.rosenbrock=chem_ros_parameters_default((chem_ros_method)variant->method);
    options.lu=variant->lu; options.layout=(chem_solver_layout){variant->width,variant->order}; options.n_external=2; options.external=external;
    if (status==CHEM_OK) status=chem_solver_create(solver,&system,0,NULL,&options,NULL);
    chem_system_destroy(&system); chem_phase_destroy(&gas);
    for (size_t i=0;i<2;++i) chem_phase_species_destroy(gas_species+i);
    return status;
}
static const char *const aerosol_test_names[]={"FO2","BAR","STUB1.MODE1.QUUX.BAZ","STUB1.MODE2.CORGE.FO2","STUB2.MODE3.CORGE.QUX","STUB2.MODE3.CORGE.BAZ","STUB2.MODE1.NUMBER","STUB2.MODE2.NUMBER","STUB1.MODE2.QUUX.BAZ","STUB2.MODE2.CORGE.FO2","STUB2.MODE2.CORGE.BAZ","STUB2.MODE3.QUUX.BAZ","STUB2.MODE3.QUUX.QUX"};
static void aerosol_case(integration_variant variant,int kind)
{
    size_t cells=kind==2 || kind>=7 ? 3 : 1; int both=kind==6 || kind==8;
    chem_solver solver={0}; REQUIRE(aerosol_create(&solver,&variant)==CHEM_OK);
    chem_state state={0}; REQUIRE(chem_solver_state(&solver,&state,cells)==CHEM_OK); const chem_model *m=state.model;
    CHECK(m->n_species==22); CHECK(m->n_custom==2);
    for (size_t i=0;i<7;++i) CHECK(chem_species_index(m,aerosol1_names[i])!=CHEM_NONE);
    for (size_t i=0;i<13;++i) CHECK(chem_species_index(m,aerosol2_names[i])!=CHEM_NONE);
    CHECK(chem_species_index(m,"FO2")!=CHEM_NONE); CHECK(chem_species_index(m,"BAR")!=CHEM_NONE);
    for (size_t i=0;i<2;++i) {
        const aerosol_context *context=chem_model_external_context(m,i); REQUIRE(context); CHECK(context->finalized==1);
    }
    if (kind==1 || kind==2) {
        const double values[8][3]={{1.34,1.35,1.36},{2.53,2.54,2.55},{.75,.76,.77},{1.23,1.24,1.25},{.42,.43,.44},{.33,.34,.35},{1000,1001,1002},{500,501,502}};
        for (size_t i=0;i<8;++i) REQUIRE(chem_state_assign_species(&state,aerosol_test_names[i],values[i],cells)==CHEM_OK);
        for (size_t i=0;i<8;++i) {
            size_t index=chem_species_index(m,aerosol_test_names[i]); const double *column=chem_state_species_column_const(&state,aerosol_test_names[i]); REQUIRE(column);
            for (size_t c=0;c<cells;++c) { NEAR(column[c],values[i][c],0); NEAR(state.y[index*cells+c],values[i][c],0); }
        }
    } else if (kind>=3) {
        const double values[13][3]={{1,1.5,2},{2,2.5,3},{.5,.7,.9},{.8,1,1.2},{.3,.4,.5},{.2,.3,.4},{1000,1100,1200},{500,550,600},{0,0,0},{.4,.6,.8},{.1,.15,.2},{.2,.3,.4},{.3,.4,.5}};
        for (size_t i=0;i<(both ? 13u : 8u);++i) {
            if (both && i==4) continue;
            const double *input=both && i==5 ? (double[]){.2,.25,.3} : values[i];
            REQUIRE(chem_state_assign_species(&state,aerosol_test_names[i],input,cells)==CHEM_OK);
        }
        if (both) {
            REQUIRE(chem_state_assign_custom(&state,aerosol_parameters[0],(double[]){.01,.015,.02},cells)==CHEM_OK);
            for (size_t c=0;c<cells;++c) state.conditions[c].temperature=(double[]){275,285,295}[c];
        }
        size_t ids[13]; for (size_t i=0;i<13;++i) ids[i]=chem_species_index(m,aerosol_test_names[i]);
        if (kind==3 || kind==4) {
            const aerosol_context *context=chem_model_external_context(m,0); double forcing[22]={0};
            if (kind==3) {
                aerosol_forcing(&state,state.y,forcing,(void *)context);
                NEAR(forcing[ids[0]],-.001,0); NEAR(forcing[ids[3]],.001,0); NEAR(forcing[ids[2]],-.001,0); NEAR(forcing[ids[8]],.001,0);
            } else {
                aerosol_jacobian(&state,state.y,state.jacobian,(void *)context);
                NEAR(state.jacobian[m->lu.map[ids[0]*22+ids[0]]],.001,0);
                NEAR(state.jacobian[m->lu.map[ids[3]*22+ids[0]]],-.001,0);
                NEAR(state.jacobian[m->lu.map[ids[2]*22+ids[2]]],.002,0);
                NEAR(state.jacobian[m->lu.map[ids[8]*22+ids[2]]],-.002,0);
            }
        } else {
            integration_backend backend={0}; REQUIRE(integration_backend_create(&backend,&state,variant.width,variant.order,variant.inplace)==CHEM_OK);
            chem_calculate_rates(&state); chem_result result=chem_solver_solve(&solver,&state,10); CHECK(result.status==CHEM_OK);
            for (size_t c=0;c<cells;++c) {
                double delta1=values[0][c]*(1-exp(-.001*10)),delta2=values[2][c]*(1-exp(-.002*10));
                NEAR(state.y[ids[0]*cells+c],values[0][c]-delta1,variant.tolerance*values[0][c]);
                NEAR(state.y[ids[2]*cells+c],values[2][c]-delta2,variant.tolerance*values[2][c]);
                NEAR(state.y[ids[8]*cells+c],delta2,variant.tolerance*values[2][c]);
                NEAR(state.y[ids[3]*cells+c],values[3][c]+delta1,variant.tolerance*values[0][c]);
                for (size_t i=0;i<8;++i) if (i==1 || i>=6 || (!both && (i==4 || i==5))) NEAR(state.y[ids[i]*cells+c],values[i][c],0);
                if (both) {
                    double d3=values[9][c]*(1-exp(-(double[]){.01,.015,.02}[c]*10));
                    double d4=values[11][c]*(1-exp(-state.conditions[c].temperature*.005*10));
                    NEAR(state.y[ids[9]*cells+c],values[9][c]-d3,variant.tolerance*values[9][c]);
                    NEAR(state.y[ids[10]*cells+c],values[10][c]+d3,variant.tolerance*values[9][c]);
                    NEAR(state.y[ids[11]*cells+c],values[11][c]-d4,variant.tolerance*values[11][c]);
                    NEAR(state.y[ids[12]*cells+c],values[12][c]+d4,variant.tolerance*values[11][c]);
                }
            }
            integration_backend_destroy(&backend);
        }
    }
    chem_state_destroy(&state); chem_solver_destroy(&solver);
}

static unsigned failed_clone_destroys;
static void *fail_clone(const void *context) { (void)context; return NULL; }
static void count_failed_clone_destroy(void *context) { (void)context; ++failed_clone_destroys; }
static void check_failed_clone(void)
{
    chem_species species={.name="A"};
    chem_external_model external={.clone=fail_clone,.destroy=count_failed_clone_destroy};
    chem_mechanism description={.n_species=1,.species=&species,.n_external=1,.external=&external};
    chem_model model={0}; failed_clone_destroys=0;
    CHECK(chem_model_create(&model,&description,CHEM_LU_MOZART)==CHEM_MEMORY);
    CHECK(failed_clone_destroys==0); CHECK(model.composition==NULL);
    chem_model_destroy(&model); CHECK(failed_clone_destroys==0);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,StateIncludesStubAerosolModel)
{
    check_failed_clone();
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],0);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,CanUpdateStateWithStubAerosolModel)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],1);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,CanUpdateMultiCellStateWithStubAerosolModel)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],2);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],3);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,5e-5},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,5e-5},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],4);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,CanSolveSingleGridCellWithStubAerosolModel1)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,2e-4},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,2e-4},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,2e-4},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,2e-4},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-4},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-4},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-4},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-4},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,2e-4},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,2e-4},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,2e-4},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,2e-4},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-4},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-4},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-4},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-4},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,2e-4},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,2e-4},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,2e-4},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,2e-4}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],5);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,CanSolveSingleGridCellWithTwoStubAerosolModels)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,1e-1},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,1e-1},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,1e-1},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,1e-1},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-1},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-1},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-1},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-1},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-1},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-1},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-1},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-1},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-1},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-1},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-1},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-1},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-1}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],6);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,CanSolveMultiGridCellWithStubAerosolModel1)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-4},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,2e-4},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,2e-4},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,2e-4},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,2e-4},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-4},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-4},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-4},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-4},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,2e-4},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,2e-4},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,2e-4},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,2e-4},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-4},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-4},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-4},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-4},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,2e-4},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,2e-4},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,2e-4},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,2e-4}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],7);
}

TEST(aerosol_model_backward_euler,AerosolModelIntegration,CanSolveMultiGridCellWithTwoStubAerosolModels)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-1},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,1e-1},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,1e-1},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,1e-1},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,1e-1},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-1},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-1},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-1},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-1},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-1},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-1},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-1},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-1},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-1},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-1},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-1},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-1},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-1}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],8);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,StateIncludesStubAerosolModel)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],0);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,CanUpdateStateWithStubAerosolModel)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],1);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,CanUpdateMultiCellStateWithStubAerosolModel)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],2);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,CanCalculateSingleGridCellForcingWithStubAerosolModel)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],3);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,CanCalculateSingleGridCellJacobianWithStubAerosolModel)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],4);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,CanSolveSingleGridCellWithStubAerosolModel1)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],5);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,CanSolveSingleGridCellWithTwoStubAerosolModels)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-4},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],6);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,CanSolveMultiGridCellWithStubAerosolModel1)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],7);
}

TEST(aerosol_model_rosenbrock,AerosolModelIntegration,CanSolveMultiGridCellWithTwoStubAerosolModels)
{
    const integration_variant variants[]={
        {0,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-4},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-5},
        {1,0,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-5},
        {1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,5e-5},
        {1,1,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,2,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,3,CHEM_LU_MOZART,CHEM_CSC,0,5e-5},
        {1,4,CHEM_LU_MOZART,CHEM_CSC,0,5e-5}
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) aerosol_case(variants[i],8);
}
