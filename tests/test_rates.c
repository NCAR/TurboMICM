/* Cases ported from test/unit/process/rate_constant in MICM (Apache-2.0). */
#include "test.h"

static void arrhenius_case(void)
{
    chem_conditions c={.temperature=301.24};
    chem_rate r=chem_rate_default(CHEM_ARRHENIUS);
    REL(chem_rate_calculate(&r,&c,NULL),1.0,1e-13);
    r.p.arrhenius.A=1;
    REL(chem_rate_calculate(&r,&c,NULL),1.0,1e-13);
    r.p.arrhenius.A=2.2e-10;
    REL(chem_rate_calculate(&r,&c,NULL),2.2e-10,1e-13);
    r.p.arrhenius.A=3e-11; r.p.arrhenius.C=-200;
    REL(chem_rate_calculate(&r,&c,NULL),3e-11*exp(-200/301.24),1e-13);
    r.p.arrhenius.A=2.6e-12; r.p.arrhenius.C=-350;
    REL(chem_rate_calculate(&r,&c,NULL),2.6e-12*exp(-350/301.24),1e-13);
}
TEST(arrhenius_rate_constant,ArrheniusRateConstant,CalculateWithSystem) { arrhenius_case(); }
TEST(arrhenius_rate_constant,ArrheniusRateConstant,CalculateWithPrescribedArguments) { arrhenius_case(); }

static void branched_case(int nitrate)
{
    chem_conditions c={.temperature=301.24,.air_density=42.2};
    chem_rate r=chem_rate_default(CHEM_BRANCHED);
    r.p.branched=(chem_branched){.nitrate=nitrate,.X=1.2,.Y=204.3,.a0=1e-3,.n=2};
    REQUIRE(chem_rate_prepare(&r)==CHEM_OK);
    double air=42.2*6.02214076e23*1e-6;
    double a=2e-22*exp(2)*2.45e19, b=0.43*pow(293.0/298,-8);
    double z=a/(1+a/b)*pow(0.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
    a=2e-22*exp(2)*air; b=0.43*pow(c.temperature/298,-8);
    double A=a/(1+a/b)*pow(0.41,1/(1+pow(log10(a/b),2)));
    double expected=1.2*exp(-204.3/c.temperature)*(nitrate ? A/(z+A) : z/(z+A));
    REL(chem_rate_calculate(&r,&c,NULL),expected,1e-13);
}
TEST(branched_rate_constant,BranchedRateConstant,CalculateAlkoxyBranchWithAllArguments) { branched_case(0); }
TEST(branched_rate_constant,BranchedRateConstant,CalculateNitrateBranchWithAllArguments) { branched_case(1); }

static double lambda_rate(const chem_conditions *c,const double *custom,void *context)
{ (void)custom; (void)context; return 2*c->temperature; }
TEST(lambda_rate_constant,LambdaRateConstant,CalculateWithSystem)
{
    chem_rate r=chem_rate_default(CHEM_LAMBDA); r.label="test lambda"; r.p.lambda.function=lambda_rate;
    chem_conditions c={.temperature=300};
    REL(chem_rate_calculate(&r,&c,NULL),600.0,1e-13);
}

TEST(reversible_rate_constant,ReversibleRateConstant,DefaultConstructor)
{
    chem_rate r=chem_rate_default(CHEM_REVERSIBLE); chem_conditions c={.temperature=298.15};
    NEAR(chem_rate_calculate(&r,&c,NULL),0.0,1e-13);
}
static void reversible(double A,double C,double kr,double T)
{
    chem_rate r=chem_rate_default(CHEM_REVERSIBLE); r.p.reversible=(chem_reversible){A,C,kr};
    chem_conditions c={.temperature=T}; double expected=A*exp(C/T)*kr;
    REL(chem_rate_calculate(&r,&c,NULL),expected,1e-13);
}
TEST(reversible_rate_constant,ReversibleRateConstant,CalculateWithSystem) { reversible(1.14e-2,2300,.32,301.24); }
TEST(reversible_rate_constant,ReversibleRateConstant,CalculateWithPrescribedArguments) { reversible(1.14e-2,2300,.32,298.15); }
TEST(reversible_rate_constant,ReversibleRateConstant,TemperatureDependence)
{ reversible(1.5,2000,.5,273.15); reversible(1.5,2000,.5,298.15); reversible(1.5,2000,.5,323.15); }
TEST(reversible_rate_constant,ReversibleRateConstant,ZeroReverseRate) { reversible(5,1000,0,298.15); }
TEST(reversible_rate_constant,ReversibleRateConstant,EquilibriumConstantOnly) { reversible(3.5,1800,1,298.15); }

static void surface(double probability)
{
    chem_rate r=chem_rate_default(CHEM_SURFACE);
    r.p.surface.diffusion=2.3e2; r.p.surface.molecular_weight=.025; r.p.surface.probability=probability;
    REQUIRE(chem_rate_prepare(&r)==CHEM_OK);
    chem_conditions c={.temperature=273.65}; double custom[]={1e-7,2.5e6};
    double expected=4*2.5e6*CHEM_PI*pow(1e-7,2)/(1e-7/2.3e2+4/(sqrt(8*CHEM_GAS_CONSTANT*273.65/(CHEM_PI*.025))*probability));
    REL(chem_rate_calculate(&r,&c,custom),expected,1e-13);
}
TEST(surface_rate_constant,SurfaceRateConstant,CalculateDefaultProbability) { surface(1); }
TEST(surface_rate_constant,SurfaceRateConstant,CalculateSpecifiedProbability) { surface(.74); }
static void missing_surface_property(int diffusion_missing)
{
    chem_species species[]={{.name="foo"}};
    size_t reactants[]={0};
    chem_reaction reaction={.n_reactants=1,.reactants=reactants};
    reaction.rate=chem_rate_default(CHEM_SURFACE); reaction.rate.label="foo";
    if (diffusion_missing) reaction.rate.p.surface.molecular_weight=.025;
    else reaction.rate.p.surface.diffusion=2.3e2;
    chem_mechanism mechanism={.n_species=1,.species=species,.n_reactions=1,.reactions=&reaction};
    chem_model model; chem_error error;
    CHECK(chem_model_create_checked(&model,&mechanism,CHEM_LU_MOZART,&error)==CHEM_NOT_FOUND);
    CHECK(error.code==CHEM_NOT_FOUND);
    CHECK(strstr(error.message,"foo")!=NULL);
    if (diffusion_missing) CHECK(strstr(error.message,"Diffusion coefficient for species 'foo' is not defined")!=NULL);
    else CHECK(strstr(error.message,"molecular weight")!=NULL);
}
TEST(surface_rate_constant,SurfaceRateConstant,DiffusionCoefficientIsMissing) { missing_surface_property(1); }
TEST(surface_rate_constant,SurfaceRateConstant,MolecularWeightIsMissing) { missing_surface_property(0); }

static void taylor(void)
{
    chem_rate r=chem_rate_default(CHEM_TAYLOR); chem_conditions c={.temperature=301.24};
    chem_taylor *p=&r.p.taylor;
    REL(chem_rate_calculate(&r,&c,NULL),1,1e-13);
    p->arrhenius.A=1; REL(chem_rate_calculate(&r,&c,NULL),1,1e-13);
    p->arrhenius.A=2.2e-10; REL(chem_rate_calculate(&r,&c,NULL),2.2e-10,1e-13);
    p->arrhenius.A=3e-11; p->arrhenius.C=-200; p->coefficients[0]=12.5; p->coefficients[1]=1.3e-2; p->coefficients[2]=5.2e-4; p->n=3;
    double expected=3e-11*exp(-200/301.24)*(12.5+1.3e-2*301.24+5.2e-4*pow(301.24,2));
    REL(chem_rate_calculate(&r,&c,NULL),expected,1e-13);
    p->arrhenius.A=2.6e-12; p->arrhenius.C=-350; p->coefficients[0]=1; p->coefficients[1]=4.3e-1; p->coefficients[2]=7.3e-3;
    expected=2.6e-12*exp(-350/301.24)*(1+4.3e-1*301.24+7.3e-3*pow(301.24,2));
    REL(chem_rate_calculate(&r,&c,NULL),expected,1e-13);
}
TEST(taylor_series_rate_constant,TaylorSeriesRateConstant,CalculateWithSystem) { taylor(); }
TEST(taylor_series_rate_constant,TaylorSeriesRateConstant,CalculateWithPrescribedArguments) { taylor(); }

static void falloff(chem_rate_type type,int full)
{
    chem_rate r=chem_rate_default(type); chem_conditions c={.temperature=301.24,.air_density=42.2};
    double k0=1,kinf=1,expected;
    if (full) {
        r.p.falloff=(chem_falloff){1.2,2.3,302.3,2.6,-3.1,402.1,.9,1.2};
        k0=1.2*exp(302.3/301.24)*pow(301.24/300,2.3);
        kinf=2.6*exp(402.1/301.24)*pow(301.24/300,-3.1);
        expected=k0*(type==CHEM_TROE ? 42.2 : 1)/(1+42.2*k0/kinf)*pow(.9,1/(1+1/1.2*pow(log10(42.2*k0/kinf),2)));
    } else {
        expected=k0*(type==CHEM_TROE ? 42.2 : 1)/(1+42.2*k0/kinf)*pow(.6,1/(1+pow(log10(42.2*k0/kinf),2)));
    }
    NEAR(chem_rate_calculate(&r,&c,NULL),expected,0);
}
TEST(troe_rate_constant,TroeRateConstant,CalculateWithMinimalArguments) { falloff(CHEM_TROE,0); }
TEST(troe_rate_constant,TroeRateConstant,CalculateWithAllArguments) { falloff(CHEM_TROE,1); }
TEST(ternary_chemical_activation_rate_constant,TernaryChemicalActivationRateConstant,CalculateWithMinimalArguments) { falloff(CHEM_TERNARY,0); }
TEST(ternary_chemical_activation_rate_constant,TernaryChemicalActivationRateConstant,CalculateWithAllArguments) { falloff(CHEM_TERNARY,1); }
TEST(troe_rate_constant,TroeRateConstant,AnalyticalTroeExampleAB)
{
    chem_rate r=chem_rate_default(CHEM_TROE); chem_conditions c={.temperature=301.24,.air_density=42.2};
    r.p.falloff.k0_A=4e-18;
    double k0=4e-18, kinf=1;
    double expected=k0*42.2/(1+k0*42.2/kinf)*pow(.6,1/(1+1/1.0*pow(log10(k0*42.2/kinf),2)));
    NEAR(chem_rate_calculate(&r,&c,NULL),expected,0);
}
TEST(troe_rate_constant,TroeRateConstant,AnalyticalTroeExampleBC)
{
    chem_rate r=chem_rate_default(CHEM_TROE); chem_conditions c={.temperature=301.24,.air_density=42.2};
    r.p.falloff=(chem_falloff){1.2e-12,167,3,136,5,24,.9,.8};
    double k0=1.2e-12*exp(3/301.24)*pow(301.24/300,167), kinf=136*exp(24/301.24)*pow(301.24/300,5);
    double expected=k0*42.2/(1+k0*42.2/kinf)*pow(.9,1/(1+(1/.8)*pow(log10(k0*42.2/kinf),2)));
    double actual=chem_rate_calculate(&r,&c,NULL);
    NEAR(actual,expected,1e-14*fmax(fabs(actual),fabs(expected)));
}
TEST(ternary_chemical_activation_rate_constant,TernaryChemicalActivationJPL19RateConstant,CalculateWithAllArguments)
{
    chem_rate r=chem_rate_default(CHEM_TERNARY_JPL19); chem_conditions c={.temperature=301.24,.air_density=42.2};
    r.p.jpl19=(chem_jpl19){1.2,2.3,302.3,298,2.6,-3.1,402.1,310,.7,.5,120,290,.9,1.2};
    double T=c.temperature, k0=1.2*exp(302.3/T)*pow(T/298,2.3), kinf=2.6*exp(402.1/T)*pow(T/310,-3.1);
    double kint=.7*exp(120/T)*pow(T/290,.5), ratio=k0*42.2/kinf;
    double kf=k0*42.2/(1+ratio)*pow(.9,1.2/(1.2+pow(log10(ratio),2)));
    NEAR(chem_rate_calculate(&r,&c,NULL),kf+kint*(1-kf/kinf),0);
}
TEST(ternary_chemical_activation_rate_constant,TernaryChemicalActivationJPL19RateConstant,ZeroKintReducesToTroe)
{
    chem_rate r=chem_rate_default(CHEM_TERNARY_JPL19); chem_conditions c={.temperature=250,.air_density=30};
    r.p.jpl19.k0_A=3e-3; r.p.jpl19.k0_B=-2; r.p.jpl19.kinf_A=1.5; r.p.jpl19.kinf_B=-1; r.p.jpl19.kint_A=0;
    double k0=3e-3*pow(250.0/298,-2), kinf=1.5*pow(250.0/298,-1), ratio=k0*30/kinf;
    double expected=k0*30/(1+ratio)*pow(.6,1/(1+pow(log10(ratio),2)));
    NEAR(chem_rate_calculate(&r,&c,NULL),expected,1e-30);
}
TEST(tunneling_rate_constant,TunnelingRateConstant,CalculateWithMinimalArguments)
{
    chem_rate r=chem_rate_default(CHEM_TUNNELING); chem_conditions c={.temperature=301.24};
    REL(chem_rate_calculate(&r,&c,NULL),1.0,1e-13);
}
TEST(tunneling_rate_constant,TunnelingRateConstant,CalculateWithAllArguments)
{
    chem_rate r=chem_rate_default(CHEM_TUNNELING); chem_conditions c={.temperature=301.24};
    r.p.tunneling=(chem_tunneling){1.2,2.3,302.3};
    double expected=1.2*exp(-2.3/301.24)*exp(302.3/pow(301.24,3));
    REL(chem_rate_calculate(&r,&c,NULL),expected,1e-13);
}
static void user_rate(double value,double scaling,const char *label)
{
    chem_rate r=chem_rate_default(CHEM_USER); r.p.user.scaling=scaling; r.label=label;
    chem_conditions c={0}; double custom[]={value};
    NEAR(chem_rate_calculate(&r,&c,custom),value*scaling,0);
}
TEST(user_defined_rate_constant,UserDefinedRateConstant,CalculateWithSystem) { user_rate(.5,1,""); }
TEST(user_defined_rate_constant,UserDefinedRateConstant,ConstructorWithRate) { user_rate(1.1,1,""); }
TEST(user_defined_rate_constant,UserDefinedRateConstant,ConstructorWithRateAndName) { user_rate(1.1,1,"name"); }
TEST(user_defined_rate_constant,UserDefinedRateConstant,CustomScalingFactor) { user_rate(1.2,2,""); }

static double lambda_override(const chem_conditions *c,const double *custom,void *context)
{ (void)custom; return *(const double *)context*c->temperature; }
TEST(solver_lambda_rate_constant,Solver,GetLambdaRateConstantByNameCanOverrideLambda)
{
    chem_species species[]={{.name="A"},{.name="B"}};
    size_t reactants[]={0}; chem_product products[]={{1,1}};
    chem_reaction rx={1,1,reactants,products,{0}};
    double scaling=1e-3, replacement=2e-3;
    rx.rate=chem_rate_default(CHEM_LAMBDA); rx.rate.label="lambda_rc";
    rx.rate.p.lambda=(chem_lambda){lambda_override,&scaling};
    chem_mechanism mechanism={.n_species=2,.n_reactions=1,.species=species,.reactions=&rx};
    chem_model model; REQUIRE(chem_model_create(&model,&mechanism,CHEM_LU_MOZART)==CHEM_OK);
    CHECK(chem_model_set_lambda(&model,"lambda_rc",lambda_override,&replacement)==CHEM_OK);
    chem_state state; REQUIRE(chem_state_create(&state,&model,1)==CHEM_OK);
    state.atol[0]=state.atol[1]=1e-6;
    double initial=1; CHECK(chem_state_set_species(&state,"A",&initial,1)==CHEM_OK);
    state.conditions[0]=(chem_conditions){300,101253.3,1e6};
    chem_calculate_rates(&state); NEAR(state.rates[0],.6,1e-12);
    chem_be_parameters parameters=chem_be_parameters_default();
    chem_result result=chem_be_solve(&state,1,&parameters);
    CHECK(result.status==CHEM_OK); NEAR(result.stats.final_time,1,0);
    chem_state_destroy(&state); chem_model_destroy(&model);
}
