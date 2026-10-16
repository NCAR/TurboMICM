#include "test.h"
#include "test_backend.h"
#include "analytical_reference.h"
typedef struct {
    chem_species species[8]; chem_reaction reactions[9]; size_t reactants[9][2]; chem_product products[9][3];
    chem_mechanism mechanism; chem_model model; chem_state state; integration_backend backend;
} special_fixture;
static const char *special_names[]={"Y0","Y1","Y2","Y3","Y4","Y5","Y6","Y7"};
static const char *special_labels[]={"r1","r2","r3","r4","r5","r6","r7","r8","r9"};
static void special_init(special_fixture *f,size_t count)
{
    *f=(special_fixture){0};
    for (size_t i=0;i<count;++i) f->species[i].name=special_names[i];
    f->mechanism=(chem_mechanism){.n_species=count,.species=f->species,.reactions=f->reactions,.preserve_order=1};
}
static void add_reaction(special_fixture *f,size_t nr,const size_t *reactants,size_t np,const chem_product *products)
{
    size_t i=f->mechanism.n_reactions++;
    if (nr) memcpy(f->reactants[i],reactants,nr*sizeof(size_t));
    if (np) memcpy(f->products[i],products,np*sizeof(chem_product));
    chem_rate rate=chem_rate_default(CHEM_USER); rate.label=special_labels[i];
    f->reactions[i]=(chem_reaction){nr,np,f->reactants[i],f->products[i],rate};
}
static chem_status special_create(special_fixture *f,integration_variant v)
{
    chem_status status=chem_model_create(&f->model,&f->mechanism,v.lu);
    if (status==CHEM_OK) status=chem_state_create(&f->state,&f->model,1);
    if (status==CHEM_OK) status=integration_backend_create(&f->backend,&f->state,v.width,v.order,v.inplace);
    return status;
}
static void special_destroy(special_fixture *f)
{ integration_backend_destroy(&f->backend); chem_state_destroy(&f->state); chem_model_destroy(&f->model); }
static void special_advance(special_fixture *f,integration_variant v,double dt,size_t substeps,int check_status)
{
    chem_ros_parameters rp=chem_ros_parameters_default(v.method<0 ? CHEM_ROS3 : (chem_ros_method)v.method);
    chem_be_parameters bp=chem_be_parameters_default(); double done=0,maxstep=dt/(double)substeps;
    while (done<dt) {
        double step=fmin(maxstep,dt-done);
        chem_result result=v.method<0 ? chem_be_solve(&f->state,step,&bp) : chem_ros_solve(&f->state,step,&rp,NULL,NULL);
        chem_state_clamp(&f->state);
        if (check_status) REQUIRE(result.status==CHEM_OK || result.status==CHEM_MAX_STEPS);
        REQUIRE(isfinite(result.stats.final_time) && done+result.stats.final_time>done);
        done+=result.stats.final_time;
    }
}
static void surface_case(integration_variant v)
{
    special_fixture f; special_init(&f,3); f.mechanism.preserve_order=0;
    add_reaction(&f,1,(size_t[]){0},2,(chem_product[]){{1,1},{2,.4}});
    f.reactions[0].rate=chem_rate_default(CHEM_SURFACE); f.reactions[0].rate.label="foo";
    f.reactions[0].rate.p.surface=(chem_surface){.diffusion=.95e-5,.molecular_weight=.04607,.probability=2e-2};
    REQUIRE(special_create(&f,v)==CHEM_OK); chem_state *s=&f.state;
    double diameter=1e-6,sd=.1,t=272.5;
    double radius=diameter/2*exp(5*log(sd)*log(sd)/2);
    double number=6/(CHEM_PI*pow(diameter,3)*exp(9./2*log(sd)*log(sd)))*(2e-3/1000+3e-3/1000);
    s->conditions[0]=(chem_conditions){t,101253.3,0}; s->custom[0]=radius; s->custom[1]=number;
    s->y[f.model.spec_to_var[0]]=1;
    double speed=sqrt(8*CHEM_GAS_CONSTANT/(CHEM_PI*.04607)*t);
    double k=4*number*CHEM_PI*radius*radius/(radius/.95e-5+4/(speed*2e-2));
    double dt=.1/k;
    for (size_t step=1;step<=10;++step) {
        chem_calculate_rates(s); special_advance(&f,v,dt,1,1); NEAR(k,s->rates[0],1e-8);
        double foo=exp(-k*(step*dt));
        REL(s->y[f.model.spec_to_var[0]],foo,v.tolerance); REL(s->y[f.model.spec_to_var[1]],1-foo,v.tolerance);
        REL(s->y[f.model.spec_to_var[2]],.4*(1-foo),v.tolerance);
    }
    special_destroy(&f);
}
static void robertson_case(integration_variant v)
{
    special_fixture f; special_init(&f,3);
    add_reaction(&f,1,(size_t[]){0},1,(chem_product[]){{1,1}});
    add_reaction(&f,2,(size_t[]){1,1},2,(chem_product[]){{1,1},{2,1}});
    add_reaction(&f,2,(size_t[]){1,2},2,(chem_product[]){{0,1},{2,1}});
    REQUIRE(special_create(&f,v)==CHEM_OK); chem_state *s=&f.state; s->rtol=1e-10;
    for (size_t i=0;i<3;++i) s->atol[i]=s->rtol*1e-2;
    memcpy(s->custom,(double[]){.04,3e7,1e4},3*sizeof(double)); s->y[0]=1;
    s->conditions[0]=(chem_conditions){272.5,101253.3,1e6}; chem_calculate_rates(s);
    double target=1,current=0;
    for (size_t step=1;step<=12;++step) {
        special_advance(&f,v,target-current,1,0); current=target; target*=10;
        for (size_t i=0;i<3;++i) {
            double error=fabs(s->y[i]-robertson_reference[step][i]);
            CHECK(error<.3e-3 || error/fabs(robertson_reference[step][i])<v.tolerance);
        }
    }
    special_destroy(&f);
}
static void hires_case(integration_variant v)
{
    special_fixture f; special_init(&f,8);
    add_reaction(&f,1,(size_t[]){0},1,(chem_product[]){{1,1}});
    add_reaction(&f,1,(size_t[]){1},2,(chem_product[]){{0,.43/8.75},{3,8.32/8.75}});
    add_reaction(&f,1,(size_t[]){2},2,(chem_product[]){{0,8.32/10.03},{3,1.71/10.03}});
    add_reaction(&f,0,NULL,1,(chem_product[]){{0,1}});
    add_reaction(&f,1,(size_t[]){3},2,(chem_product[]){{2,.43/1.12},{5,.69/1.12}});
    add_reaction(&f,1,(size_t[]){4},2,(chem_product[]){{2,.035/1.745},{5,1.71/1.745}});
    add_reaction(&f,1,(size_t[]){5},1,(chem_product[]){{4,1}});
    add_reaction(&f,1,(size_t[]){6},3,(chem_product[]){{4,.43/1.81},{5,.69/1.81},{7,1}});
    add_reaction(&f,2,(size_t[]){5,7},1,(chem_product[]){{6,1}});
    REQUIRE(special_create(&f,v)==CHEM_OK); chem_state *s=&f.state; s->rtol=1e-6;
    for (size_t i=0;i<8;++i) s->atol[i]=s->rtol*1e-2;
    memcpy(s->custom,(double[]){1.71,8.75,10.03,.0007,1.12,1.745,.43,1.81,280},9*sizeof(double));
    s->y[0]=1; s->y[7]=.0057; chem_calculate_rates(s);
    double dt=321.8122;
    for (size_t step=1;step<=2;++step) {
        special_advance(&f,v,dt,1,0); dt+=100;
        for (size_t i=0;i<8;++i) NEAR(s->y[i],hires_reference[step][i],v.tolerance);
    }
    special_destroy(&f);
}
static void e5_case(integration_variant v)
{
    special_fixture f; special_init(&f,6);
    add_reaction(&f,1,(size_t[]){0},2,(chem_product[]){{1,1},{2,1}});
    add_reaction(&f,2,(size_t[]){1,2},1,(chem_product[]){{4,1}});
    add_reaction(&f,2,(size_t[]){0,2},1,(chem_product[]){{3,1}});
    add_reaction(&f,1,(size_t[]){3},2,(chem_product[]){{2,1},{5,1}});
    REQUIRE(special_create(&f,v)==CHEM_OK); chem_state *s=&f.state; s->rtol=1e-13;
    for (size_t i=0;i<6;++i) s->atol[i]=(i==0 || i>=4) ? 1e-7 : 1e-17;
    memcpy(s->custom,(double[]){7.89e-10,1.13e9,1.1e7,1.13e3},4*sizeof(double)); s->y[0]=1.76e-3; chem_calculate_rates(s);
    double target=10,current=0;
    for (size_t step=0;step<=7;++step) {
        if (step) { special_advance(&f,v,target-current,1,0); current=target; target*=100; }
        for (size_t i=0;i<4;++i) {
            double error=fabs(s->y[i]-e5_reference[step][i]);
            CHECK(error<(i==0 ? 1e-6 : 1e-13) || error/fabs(e5_reference[step][i])<v.tolerance);
        }
    }
    special_destroy(&f);
}
static void oregonator_case(integration_variant v)
{
    special_fixture f; special_init(&f,5);
    add_reaction(&f,1,(size_t[]){1},1,(chem_product[]){{0,1}});
    add_reaction(&f,2,(size_t[]){0,1},1,(chem_product[]){{3,1}});
    add_reaction(&f,1,(size_t[]){0},2,(chem_product[]){{2,1},{0,2}});
    add_reaction(&f,2,(size_t[]){0,0},1,(chem_product[]){{4,1}});
    add_reaction(&f,1,(size_t[]){2},1,(chem_product[]){{1,1}});
    REQUIRE(special_create(&f,v)==CHEM_OK); chem_state *s=&f.state;
    const double sc=77.27,q=8.375e-6,w=.161,k2=1.6e9,k5=1;
    double k1=k5/(sc*w),k3=sc*k5/w,k4=sc*sc*q*k2/2;
    double alpha=k5/(sc*w*k2),eta=sc*k5/(w*k2),rho=k5*k5/(w*w*k2),tau=w/k5;
    s->rtol=1e-9; double floor=s->rtol*1e-2;
    memcpy(s->atol,(double[]){alpha*floor,eta*floor,rho*floor,eta,eta},5*sizeof(double));
    memcpy(s->custom,(double[]){k1,k2,k3,k4,k5},5*sizeof(double));
    s->y[0]=alpha; s->y[1]=2*eta; s->y[2]=3*rho; chem_calculate_rates(s);
    for (size_t step=1;step<=12;++step) {
        special_advance(&f,v,30*tau,v.method<0 ? 18000 : 1,1);
        for (size_t i=0;i<3;++i) {
            double expected=oregonator_reference[step][i]*(i==0 ? alpha : i==1 ? eta : rho);
            CHECK(fabs(s->y[i]-expected)/fabs(expected)<v.tolerance);
        }
    }
    special_destroy(&f);
}

TEST(analytical_backward_euler,AnalyticalExamples,SurfaceRxn)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,0.05},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) surface_case(variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,HIRES)
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
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-1},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) hires_case(variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,Oregonator)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,5e-3},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,5e-3},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,5e-3},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,5e-3},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,5e-3},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,5e-3},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,5e-3},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,5e-3},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,5e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) oregonator_case(variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,SurfaceRxn)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) surface_case(variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,Robertson)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,0,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) robertson_case(variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,E5)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,0,CHEM_LU_MOZART,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-10},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-10},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-10},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-10},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-10},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-10},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-10},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-10},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-10},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-10},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-10},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-10},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-10},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-10},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-10},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-10},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) e5_case(variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,Oregonator)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,0,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) oregonator_case(variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,HIRES)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,0,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,0,1e-7},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-7},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-7},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-7},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-7},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-7},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-7},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-7},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-7},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-7},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-7},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-7},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-7},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-7},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-7},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-7},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-7},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) hires_case(variants[i]);
}
