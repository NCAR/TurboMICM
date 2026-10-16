#include "test.h"
#include "chem_process.h"
#include "test_random.h"

static chem_status process_phase(chem_phase *phase,int simple)
{
    const char *names[]={simple ? "O3" : "foo",simple ? "NO" : "bar","NO2","O2"};
    size_t count=simple ? 4 : 2; chem_phase_species species[4]={0}; chem_status status=CHEM_OK;
    for (size_t i=0;i<count;++i) {
        status=chem_species_info_create(&species[i].info,names[i],0,NULL); if (status!=CHEM_OK) break;
        if (!simple && i==0) {
            status=chem_species_property_set(&species[i].info,"molecular weight [kg mol-1]",CHEM_PROPERTY_DOUBLE,(chem_property_value){.real=.025});
            chem_phase_species_set_diffusion(species+i,2.3e2);
        } else if (!simple && i==1) chem_species_parameterize(&species[i].info,(chem_parameterization){.air_density=.82});
    }
    if (status==CHEM_OK) status=chem_phase_create(phase,"gas",count,species);
    for (size_t i=0;i<count;++i) chem_phase_species_destroy(species+i);
    return status;
}
TEST(process,Process,BuildsChemicalReaction)
{
    chem_phase phase={0}; REQUIRE(process_phase(&phase,1)==CHEM_OK); chem_process process={0};
    REQUIRE(chem_process_create(&process,&phase,2,(const char *[]){"O3","NO"},2,
        (chem_named_product[]){{"NO2",1},{"O2",1}},chem_rate_default(CHEM_ARRHENIUS))==CHEM_OK);
    chem_phase_destroy(&phase); CHECK(process.reaction.n_reactants==2);
    CHECK(!strcmp(process.phase.species[process.reaction.products[0].species].info.species.name,"NO2"));
    CHECK(!strcmp(process.phase.name,"gas")); chem_process_destroy(&process);
}
TEST(process,Process,ChemicalReactionCopyAssignmentSucceeds)
{
    chem_phase phase={0}; REQUIRE(process_phase(&phase,0)==CHEM_OK); chem_process original={0},copy={0};
    chem_rate rate=chem_rate_default(CHEM_ARRHENIUS); rate.p.arrhenius.A=12.2; rate.p.arrhenius.C=300;
    REQUIRE(chem_process_create(&original,&phase,2,(const char *[]){"foo","bar"},0,NULL,rate)==CHEM_OK);
    REQUIRE(chem_process_replace(&copy,&original)==CHEM_OK);
    CHECK(!strcmp(copy.phase.species[copy.reaction.reactants[0]].info.species.name,original.phase.species[original.reaction.reactants[0]].info.species.name));
    CHECK(copy.reaction.n_products==original.reaction.n_products); CHECK(!strcmp(copy.phase.name,original.phase.name));
    CHECK(copy.reaction.rate.type==CHEM_ARRHENIUS); NEAR(copy.reaction.rate.p.arrhenius.A,original.reaction.rate.p.arrhenius.A,0);
    chem_process_destroy(&original); chem_phase_destroy(&phase); NEAR(copy.reaction.rate.p.arrhenius.A,12.2,0);
    chem_process_destroy(&copy);
}
static void process_update(size_t width)
{
    chem_phase phase={0}; REQUIRE(process_phase(&phase,0)==CHEM_OK); chem_process processes[3]={0};
    chem_rate arr=chem_rate_default(CHEM_ARRHENIUS),surf,user=chem_rate_default(CHEM_USER);
    arr.p.arrhenius.A=12.2; arr.p.arrhenius.C=300; user.label="bar_user";
    REQUIRE(chem_surface_rate(&surf,"foo_surf",phase.species,1)==CHEM_OK);
    REQUIRE(chem_process_create(processes,&phase,2,(const char *[]){"foo","bar"},0,NULL,arr)==CHEM_OK);
    REQUIRE(chem_process_create(processes+1,&phase,1,(const char *[]){"foo"},0,NULL,surf)==CHEM_OK);
    REQUIRE(chem_process_create(processes+2,&phase,1,(const char *[]){"bar"},0,NULL,user)==CHEM_OK);
    chem_system system={0}; REQUIRE(chem_system_create(&system,&phase)==CHEM_OK); chem_process_model model={0};
    REQUIRE(chem_process_model_create(&model,&system,3,processes,1,0,CHEM_LU_MOZART,NULL)==CHEM_OK);
    for (size_t i=0;i<3;++i) chem_process_destroy(processes+i);
    chem_system_destroy(&system); chem_phase_destroy(&phase);
    random_state random; random_init(&random,1,0);
    if (!width) width=5;
    for (size_t start=0;start<5;start+=width) {
        size_t cells=5-start<width ? 5-start : width; chem_state s={0}; REQUIRE(chem_state_create(&s,&model.model,cells)==CHEM_OK);
        for (size_t c=0;c<cells;++c) {
            s.conditions[c].temperature=lognormal_at(&random,0,.01)*285;
            s.conditions[c].pressure=lognormal_at(&random,0,.01)*101100;
            s.conditions[c].air_density=lognormal_at(&random,0,.01)*10;
            s.custom[c]=lognormal_at(&random,0,.01)*1e-2;
            s.custom[cells+c]=lognormal_at(&random,0,.01)*1e-8; s.custom[2*cells+c]=lognormal_at(&random,0,.01)*1e5;
        }
        chem_calculate_rates(&s);
        for (size_t c=0;c<cells;++c) {
            double t=s.conditions[c].temperature,rho=s.conditions[c].air_density,radius=s.custom[cells+c],number=s.custom[2*cells+c];
            double expected=12.2*exp(300/t)*(rho*.82); NEAR(s.rates[c],expected,1e-12*expected);
            expected=s.custom[c]*(rho*.82); NEAR(s.rates[cells+c],expected,1e-12*expected);
            double speed=sqrt(8*CHEM_GAS_CONSTANT/(CHEM_PI*.025)*t);
            expected=4*number*CHEM_PI*radius*radius/(radius/2.3e2+4/speed); NEAR(s.rates[2*cells+c],expected,1e-10*expected);
        }
        chem_state_destroy(&s);
    }
    chem_process_model_destroy(&model);
}
TEST(process,Process,Matrix) { process_update(0); }
TEST(process,Process,VectorMatrix) { for (size_t width=1;width<=4;++width) process_update(width); }

typedef struct { chem_species species[3]; chem_reaction reactions[8]; chem_mechanism mechanism; chem_model model; size_t reactant; chem_product product; } store_fixture;
static void store_init(store_fixture *f,size_t count)
{
    *f=(store_fixture){0}; f->species[0].name="a"; f->species[1].name="b"; f->species[2].name="c"; f->product=(chem_product){1,1};
    f->mechanism=(chem_mechanism){.n_species=3,.species=f->species,.n_reactions=count,.reactions=f->reactions,.preserve_order=1};
    for (size_t i=0;i<count;++i) f->reactions[i]=(chem_reaction){1,1,&f->reactant,&f->product,chem_rate_default(CHEM_ARRHENIUS)};
}
static chem_status store_build(store_fixture *f)
{ chem_reactions_sort(f->reactions,f->mechanism.n_reactions); return chem_model_create(&f->model,&f->mechanism,CHEM_LU_MOZART); }
static chem_rate surface_params(double diffusion,double probability)
{
    chem_rate r=chem_rate_default(CHEM_SURFACE); r.label="surf"; r.p.surface=(chem_surface){diffusion,.025,probability,0}; return r;
}
TEST(reaction_rate_store,ReactionRateConstantStore,OffsetsAreContiguousCumulativeSizes)
{
    store_fixture f; store_init(&f,6); f.reactions[2].rate=chem_rate_default(CHEM_TROE);
    f.reactions[2].rate.p.falloff.k0_A=1; f.reactions[2].rate.p.falloff.kinf_A=1;
    f.reactions[3].rate=chem_rate_default(CHEM_TUNNELING);
    f.reactions[4].rate=chem_rate_default(CHEM_USER); f.reactions[4].rate.label="ud1";
    f.reactions[5].rate=surface_params(1e-5,1); REQUIRE(store_build(&f)==CHEM_OK);
    const chem_rate_type types[]={CHEM_ARRHENIUS,CHEM_TROE,CHEM_TERNARY,CHEM_TERNARY_JPL19,CHEM_BRANCHED,CHEM_TUNNELING,CHEM_TAYLOR,CHEM_REVERSIBLE,CHEM_USER,CHEM_SURFACE,CHEM_LAMBDA};
    const size_t counts[]={2,1,0,0,0,1,0,0,1,1,0},offsets[]={0,2,3,3,3,3,4,4,4,5,6};
    for (size_t i=0;i<11;++i) { CHECK(chem_model_rate_count(&f.model,types[i])==counts[i]); CHECK(chem_model_rate_offset(&f.model,types[i])==offsets[i]); }
    CHECK(chem_model_rate_offset(&f.model,CHEM_LAMBDA)+chem_model_rate_count(&f.model,CHEM_LAMBDA)==6); chem_model_destroy(&f.model);
}
TEST(reaction_rate_store,ReactionRateConstantStore,ArrheniusParametersPreserved)
{
    store_fixture f; store_init(&f,1); f.reactions[0].rate.p.arrhenius=(chem_arrhenius){2.15e-4,1.2,110,300,0};
    REQUIRE(store_build(&f)==CHEM_OK); CHECK(chem_model_rate_count(&f.model,CHEM_ARRHENIUS)==1);
    chem_arrhenius p=f.model.rates[0].p.arrhenius; NEAR(p.A,2.15e-4,0); NEAR(p.B,1.2,0); NEAR(p.C,110,0); NEAR(p.D,300,0); chem_model_destroy(&f.model);
}
TEST(reaction_rate_store,ReactionRateConstantStore,BranchedDerivedFieldsComputed)
{
    store_fixture f; store_init(&f,1); f.reactions[0].rate=chem_rate_default(CHEM_BRANCHED);
    f.reactions[0].rate.p.branched=(chem_branched){.nitrate=0,.X=1,.Y=0,.a0=.5,.n=3}; REQUIRE(store_build(&f)==CHEM_OK);
    CHECK(chem_model_rate_count(&f.model,CHEM_BRANCHED)==1); double k0=2e-22*6.02214076e23*1e-6*exp(3);
    NEAR(f.model.rates[0].p.branched.k0,k0,1e-10*k0); double air=2.45e19/6.02214076e23*1e6,a=k0*air,b=.43*pow(293./298.,-8);
    double value=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2))),z=value*(1-.5)/.5;
    NEAR(f.model.rates[0].p.branched.z,z,1e-10*fabs(z)); chem_model_destroy(&f.model);
}
TEST(reaction_rate_store,ReactionRateConstantStore,UserDefinedCustomParamIndex)
{
    store_fixture f; store_init(&f,2);
    for (size_t i=0;i<2;++i) { f.reactions[i].rate=chem_rate_default(CHEM_USER); f.reactions[i].rate.label=i ? "rate2" : "rate1"; f.reactions[i].rate.p.user.scaling=i ? .5 : 2; }
    REQUIRE(store_build(&f)==CHEM_OK); CHECK(chem_model_rate_count(&f.model,CHEM_USER)==2);
    CHECK(f.model.custom_offset[0]==0); CHECK(f.model.custom_offset[1]==1); NEAR(f.model.rates[0].p.user.scaling,2,0); NEAR(f.model.rates[1].p.user.scaling,.5,0);
    chem_model_destroy(&f.model);
}
TEST(reaction_rate_store,ReactionRateConstantStore,SurfaceDataFieldsAndCustomParamIndex)
{
    store_fixture f; store_init(&f,1); f.reactions[0].rate=surface_params(2.3e2,.74); REQUIRE(store_build(&f)==CHEM_OK);
    CHECK(chem_model_rate_count(&f.model,CHEM_SURFACE)==1); chem_surface p=f.model.rates[0].p.surface;
    NEAR(p.diffusion,2.3e2,0); NEAR(p.speed_factor,8*CHEM_GAS_CONSTANT/(CHEM_PI*.025),1e-14); NEAR(p.probability,.74,0);
    CHECK(f.model.custom_offset[0]==0); chem_model_destroy(&f.model);
}
TEST(reaction_rate_store,ReactionRateConstantStore,SurfaceCustomParamIndexAfterUserDefined)
{
    store_fixture f; store_init(&f,2); f.reactions[0].rate=chem_rate_default(CHEM_USER); f.reactions[0].rate.label="ud";
    f.reactions[1].rate=surface_params(1e-5,1); REQUIRE(store_build(&f)==CHEM_OK);
    CHECK(chem_model_rate_count(&f.model,CHEM_USER)==1); CHECK(chem_model_rate_count(&f.model,CHEM_SURFACE)==1);
    CHECK(f.model.custom_offset[0]==0); CHECK(f.model.custom_offset[1]==1); CHECK(f.model.n_custom==3); chem_model_destroy(&f.model);
}
static double store_lambda(const chem_conditions *c,const double *custom,void *context)
{ (void)custom; (void)context; return c->temperature*1e-3; }
TEST(reaction_rate_store,ReactionRateConstantStore,LambdaEntriesRcIndex)
{
    store_fixture f; store_init(&f,2); f.reactions[1].rate=chem_rate_default(CHEM_LAMBDA);
    f.reactions[1].rate.label="lam"; f.reactions[1].rate.p.lambda.function=store_lambda; REQUIRE(store_build(&f)==CHEM_OK);
    CHECK(chem_model_rate_count(&f.model,CHEM_ARRHENIUS)==1); CHECK(chem_model_rate_count(&f.model,CHEM_LAMBDA)==1);
    CHECK(chem_model_rate_offset(&f.model,CHEM_LAMBDA)==1); CHECK(f.model.rates[1].p.lambda.function!=NULL);
    NEAR(chem_rate_calculate(f.model.rates+1,&(chem_conditions){.temperature=300},NULL),.3,1e-14); chem_model_destroy(&f.model);
}
TEST(reaction_rate_store,ReactionRateConstantStore,ParameterizedMultipliers)
{
    store_fixture f; store_init(&f,2); f.species[1].type=CHEM_SPECIES_PARAMETER; f.species[1].affine.air_density=2;
    size_t reactant=1; chem_product product={0,1}; f.reactions[1].reactants=&reactant; f.reactions[1].products=&product;
    REQUIRE(store_build(&f)==CHEM_OK); CHECK(chem_model_parameterized_count(&f.model)==1); CHECK(chem_model_parameterized_index(&f.model,0)==1);
    NEAR(chem_model_parameterized_multiplier(&f.model,1,&(chem_conditions){.air_density=5}),10,1e-14); chem_model_destroy(&f.model);
}
TEST(reaction_rate_store,ReactionRateConstantStore,TotalRateConstantsMatchesProcessCount)
{
    store_fixture f; store_init(&f,5); f.reactions[1].rate=chem_rate_default(CHEM_TROE);
    f.reactions[1].rate.p.falloff.k0_A=1; f.reactions[1].rate.p.falloff.kinf_A=1;
    f.reactions[2].rate=chem_rate_default(CHEM_USER); f.reactions[2].rate.label="ud";
    f.reactions[3].rate=surface_params(1e-5,1); f.reactions[4].rate=chem_rate_default(CHEM_LAMBDA);
    f.reactions[4].rate.label="lam"; f.reactions[4].rate.p.lambda.function=store_lambda;
    REQUIRE(store_build(&f)==CHEM_OK); CHECK(chem_model_rate_offset(&f.model,CHEM_LAMBDA)+chem_model_rate_count(&f.model,CHEM_LAMBDA)==5); chem_model_destroy(&f.model);
}
TEST(reaction_rate_store,ReactionRateConstantStore,SurfaceMissingDiffusionCoefficientThrows)
{
    store_fixture f; store_init(&f,1); f.reactions[0].rate=surface_params(0,1); CHECK(store_build(&f)==CHEM_INVALID);
}
TEST(reaction_rate_store,ReactionRateConstantStore,SurfaceMissingMolecularWeightThrows)
{
    store_fixture f; store_init(&f,1); f.reactions[0].rate=surface_params(1e-5,1); f.reactions[0].rate.p.surface.molecular_weight=0;
    CHECK(store_build(&f)==CHEM_INVALID);
}
