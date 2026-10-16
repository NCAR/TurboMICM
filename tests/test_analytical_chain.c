/* Analytical reaction chains from the supplied MICM integration policies. */
#include "test.h"
#include "chem_process.h"
#include "test_backend.h"
enum { CHAIN_TROE,CHAIN_PHOTO,CHAIN_TERNARY,CHAIN_TUNNEL,CHAIN_ARRHENIUS,CHAIN_BRANCHED };
static chem_rate chain_rate(int kind,int second,int stiff)
{
    chem_rate r=chem_rate_default(kind==CHAIN_TROE ? CHEM_TROE : kind==CHAIN_PHOTO ? CHEM_USER :
        kind==CHAIN_TERNARY ? CHEM_TERNARY : kind==CHAIN_TUNNEL ? CHEM_TUNNELING : kind==CHAIN_ARRHENIUS ? CHEM_ARRHENIUS : CHEM_BRANCHED);
    if (kind==CHAIN_TROE || kind==CHAIN_TERNARY) {
        if (second) r.p.falloff=(chem_falloff){1.2e-3,1.6,3,136,5,24,.9,.8};
        else { r.p.falloff.k0_A=kind==CHAIN_TROE ? 4e-11 : 4e-5; r.p.falloff.kinf_A=1; }
    } else if (kind==CHAIN_PHOTO) r.label=second ? "photoB" : "photoA";
    else if (kind==CHAIN_TUNNEL) r.p.tunneling=(chem_tunneling){second ? 1.2e-4 : 4e-3,second ? 167 : 0,second ? 1e8 : 0};
    else if (kind==CHAIN_ARRHENIUS) {
        if (second) r.p.arrhenius=(chem_arrhenius){1.2e-4,stiff ? 1.6 : 7,75,50,.5};
        else { r.p.arrhenius.A=4e-3; r.p.arrhenius.C=50; }
    } else r.p.branched=(chem_branched){.nitrate=second,.X=1e-4,.Y=204.3,.a0=1e-3,.n=2};
    return r;
}
static double chain_expected_rate(int kind,int second,int stiff,const chem_conditions *c)
{
    double t=c->temperature,p=c->pressure,rho=c->air_density;
    if (kind==CHAIN_TROE || kind==CHAIN_TERNARY) {
        double k0=second ? 1.2e-3*exp(3/t)*pow(t/300,1.6) : kind==CHAIN_TROE ? 4e-11 : 4e-5;
        double ki=second ? 136*exp(24/t)*pow(t/300,5) : 1;
        double fc=second ? .9 : .6,N=second ? .8 : 1,logratio=log10(k0*rho/ki);
        double exponent=(kind==CHAIN_TROE && !stiff && second) ? N/(N+pow(logratio,2)) : 1/(1+(1/N)*pow(logratio,2));
        return k0*(kind==CHAIN_TROE ? rho : 1)/(1+k0*rho/ki)*pow(fc,exponent);
    }
    if (kind==CHAIN_PHOTO) return second ? 3e-3 : 2e-3;
    if (kind==CHAIN_TUNNEL) return second ? 1.2e-4*exp(-167/t+1e8/pow(t,3)) : 4e-3;
    if (kind==CHAIN_ARRHENIUS) return second ? 1.2e-4*exp(75/t)*pow(t/50,stiff ? 1.6 : 7)*(1+.5*p) : 4e-3*exp(50/t);
    double air=rho*6.02214076e23*1e-6,a=2e-22*exp(2)*2.45e19,b=.43*pow(293./298.,-8);
    double z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
    a=2e-22*exp(2)*air; b=.43*pow(t/298,-8);
    double A=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
    return 1e-4*exp(-204.3/t)*((second ? A : z)/(z+A));
}
static void analytical_chain(int kind,int stiff,integration_variant variant)
{
    chem_species species[]={{.name=stiff ? "A1" : "A"},{.name=stiff ? "A2" : "B"},{.name=stiff ? "B" : "C"},{.name="C"}};
    size_t rr[]={0,1,2,0,1}; chem_product pp[5];
    chem_reaction rx[5]; size_t count=stiff ? 5 : 2;
    for (size_t i=0;i<count;++i) {
        pp[i]=(chem_product){stiff ? (i<2 ? 2 : i==2 ? 3 : i==3 ? 1 : 0) : i+1,1};
        chem_rate rate=chain_rate(kind,stiff ? i==2 : i==1,stiff);
        if (stiff && i>=3) { rate=chem_rate_default(CHEM_ARRHENIUS); rate.p.arrhenius.A=i==3 ? 4e10 : .9*4e10; }
        rx[i]=(chem_reaction){1,1,rr+i,pp+i,rate};
    }
    chem_reactions_sort(rx,count);
    chem_mechanism mechanism={.n_species=stiff ? 4 : 3,.species=species,.n_reactions=count,.reactions=rx};
    chem_model m={0}; chem_state s={0}; REQUIRE(chem_model_create(&m,&mechanism,variant.lu)==CHEM_OK);
    REQUIRE(chem_state_create(&s,&m,3)==CHEM_OK); integration_backend backend;
    REQUIRE(integration_backend_create(&backend,&s,variant.width,variant.order,variant.inplace)==CHEM_OK);
    const double temperatures[]={272.5,254.7,312.6},pressures[]={101253.3,100672.5,101319.8};
    double k1[3],k2[3]; size_t a1=m.spec_to_var[0],a2=m.spec_to_var[1],b=m.spec_to_var[stiff ? 2 : 1],cidx=m.spec_to_var[stiff ? 3 : 2];
    for (size_t c=0;c<3;++c) {
        s.conditions[c]=(chem_conditions){temperatures[c],pressures[c],pressures[c]/(CHEM_GAS_CONSTANT*temperatures[c])};
        k1[c]=chain_expected_rate(kind,0,stiff,s.conditions+c); k2[c]=chain_expected_rate(kind,1,stiff,s.conditions+c);
        s.y[a1*3+c]=(stiff ? .5 : 1)*(1-(double)c/3);
        if (stiff) s.y[a2*3+c]=.5*(1-(double)c/3);
    }
    if (kind==CHAIN_PHOTO) {
        REQUIRE(chem_state_set_custom(&s,"photoA",(double[]){2e-3},1)==CHEM_OK);
        REQUIRE(chem_state_set_custom(&s,"photoB",(double[]){3e-3},1)==CHEM_OK);
    }
    chem_ros_parameters rp=chem_ros_parameters_default(variant.method<0 ? CHEM_ROS3 : (chem_ros_method)variant.method);
    chem_be_parameters bp=chem_be_parameters_default();
    for (size_t step=1;step<1000;++step) {
        chem_calculate_rates(&s);
        chem_result result=variant.method<0 ? chem_be_solve(&s,1,&bp) : chem_ros_solve(&s,1,&rp,NULL,NULL);
        REQUIRE(result.status==CHEM_OK); chem_state_clamp(&s);
        for (size_t c=0;c<3;++c) {
            if (!stiff) { REL(k1[c],s.rates[c],1e-9); REL(k2[c],s.rates[3+c],1e-9); }
            double initial=1-(double)c/3,t=(double)step;
            double ea=initial*exp(-k1[c]*t),eb=initial*(k1[c]/(k2[c]-k1[c]))*(exp(-k1[c]*t)-exp(-k2[c]*t));
            double ec=initial*(1+(k1[c]*exp(-k2[c]*t)-k2[c]*exp(-k1[c]*t))/(k2[c]-k1[c]));
            NEAR(s.y[a1*3+c]+(stiff ? s.y[a2*3+c] : 0),ea,variant.tolerance);
            NEAR(s.y[b*3+c],eb,variant.tolerance); NEAR(s.y[cidx*3+c],ec,variant.tolerance);
        }
    }
    integration_backend_destroy(&backend); chem_state_destroy(&s); chem_model_destroy(&m);
}

TEST(analytical_backward_euler,AnalyticalExamples,Troe)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-6},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-6},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,0,1e-6},
        {-1,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {-1,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {-1,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {-1,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-6},
        {-1,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {-1,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {-1,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
        {-1,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-6},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TROE,0,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,TroeSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TROE,1,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,Photolysis)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_PHOTO,0,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,PhotolysisSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_PHOTO,1,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,TernaryChemicalActivation)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TERNARY,0,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,TernaryChemicalActivationSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TERNARY,1,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,Tunneling)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TUNNEL,0,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,TunnelingSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TUNNEL,1,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,Arrhenius)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_ARRHENIUS,0,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,ArrheniusSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_ARRHENIUS,1,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,Branched)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_BRANCHED,0,variants[i]);
}
TEST(analytical_backward_euler,AnalyticalExamples,BranchedSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {-1,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {-1,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {-1,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {-1,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
        {-1,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-2},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_BRANCHED,1,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,Troe)
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
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
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
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TROE,0,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,TroeSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-5},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-5},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-5},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-5},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-5},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-5},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-5},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-5},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-5},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TROE,1,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,Photolysis)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-6},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,2e-6},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,2e-6},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,2e-6},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,2e-6},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-6},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-6},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-6},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-6},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,2e-6},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,2e-6},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,2e-6},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,2e-6},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_PHOTO,0,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,PhotolysisSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,2e-5},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_PHOTO,1,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,TernaryChemicalActivation)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-8},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TERNARY,0,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,TernaryChemicalActivationSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,2e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TERNARY,1,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,Tunneling)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-8},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-8},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-8},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-8},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TUNNEL,0,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,TunnelingSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-4},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-4},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-4},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-4},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-4},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-4},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-4},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-4},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-4},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-4},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-4},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-4},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-4},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-4},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-4},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-4},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_TUNNEL,1,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,Arrhenius)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,4e-6},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-9},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-9},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-9},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-9},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-9},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-9},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-9},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-9},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-9},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-9},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-9},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-9},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-9},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-9},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-9},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-9},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-9},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-9},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-9},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-9},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-9},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-9},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-9},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-9},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_ARRHENIUS,0,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,ArrheniusSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-4},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,2e-5},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,2e-5},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,2e-5},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_ARRHENIUS,1,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,Branched)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-10},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-13},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-13},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-13},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-13},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-13},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-13},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-13},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,1e-13},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-13},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-13},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-13},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,1e-13},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,1e-13},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,1e-13},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,1e-13},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,1e-13},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-13},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-13},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-13},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,1e-13},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,1e-13},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,1e-13},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,1e-13},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,1e-13},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_BRANCHED,0,variants[i]);
}
TEST(analytical_rosenbrock_integration,AnalyticalExamples,BranchedSuperStiffButAnalytical)
{
    const integration_variant variants[]={
        {CHEM_ROS2,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_RODA3,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_RODA4,0,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,0,2e-3},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSC,0,2e-3},
        {CHEM_ROS3,1,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,2,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,3,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,4,CHEM_LU_DOOLITTLE,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,1,CHEM_LU_MOZART,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,2,CHEM_LU_MOZART,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,3,CHEM_LU_MOZART,CHEM_CSR,1,2e-3},
        {CHEM_ROS3,4,CHEM_LU_MOZART,CHEM_CSR,1,2e-3},
    };
    for (size_t i=0;i<sizeof(variants)/sizeof(*variants);++i) analytical_chain(CHAIN_BRANCHED,1,variants[i]);
}
