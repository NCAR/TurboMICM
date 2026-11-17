#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define S(name) {name,0,1e-3,0,0}
static const ChemSpecies gas[]={S("FO2"),S("BAR")};
static const ChemSpecies model1_species[]={
    S("STUB1.MODE1.QUUX.BAZ"),S("STUB1.MODE1.QUUX.QUX"),S("STUB1.MODE2.QUUX.BAZ"),
    S("STUB1.MODE2.QUUX.QUX"),S("STUB1.MODE2.CORGE.FO2"),S("STUB1.MODE2.CORGE.BAZ"),S("STUB1.MODE2.CORGE.QUX")};
static const ChemSpecies model2_species[]={
    S("STUB2.MODE1.NUMBER"),S("STUB2.MODE1.QUUX.BAZ"),S("STUB2.MODE1.QUUX.QUX"),
    S("STUB2.MODE2.NUMBER"),S("STUB2.MODE2.CORGE.FO2"),S("STUB2.MODE2.CORGE.BAZ"),
    S("STUB2.MODE2.CORGE.QUX"),S("STUB2.MODE3.NUMBER"),S("STUB2.MODE3.QUUX.BAZ"),
    S("STUB2.MODE3.QUUX.QUX"),S("STUB2.MODE3.CORGE.FO2"),S("STUB2.MODE3.CORGE.BAZ"),
    S("STUB2.MODE3.CORGE.QUX")};
static const char *const parameter_names[]={"STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT",
                                            "STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT"};
static const double parameter_defaults[]={0,0};
static const ChemJacobianEntry pattern1[]={{0,0},{6,0},{2,2},{4,2}};
static const ChemJacobianEntry pattern2[]={{13,13},{14,13},{17,17},{18,17}};
static const ChemMechanism base={2,0,0,gas,NULL,NULL,NULL,NULL,NULL,298.15,101325};

typedef struct { const char *source,*product; int src,dst,parameter,ss,ps; double rate; } Transfer;
typedef struct { Transfer transfers[2]; int temperature_parameter,finalized; } Model;
typedef struct { Model contexts[2]; ChemExternalModel models[2]; ChemMechanismComposition composition; ChemState *state; } Fixture;

static int finalize(const ChemMechanism *m, const ChemLinearSystem *j, void *v)
{
    Model *model=v;
    for (size_t i=0;i<2;++i) {
        Transfer *t=&model->transfers[i];
        t->src=chem_find_species(m,t->source); t->dst=chem_find_species(m,t->product);
        if (t->src<0 || t->dst<0) return CHEM_INVALID;
        t->ss=j->lookup[(size_t)t->src*j->species_count+(size_t)t->src];
        t->ps=j->lookup[(size_t)t->dst*j->species_count+(size_t)t->src];
        if (t->ss<0 || t->ps<0) return CHEM_INVALID;
    }
    ++model->finalized;
    return CHEM_OK;
}
static int update(ChemState *s, void *v)
{
    Model *model=v;
    if (model->temperature_parameter>=0) {
        double *p=chem_parameter_values(s,(size_t)model->temperature_parameter), *t=chem_temperatures(s);
        for (size_t c=0;c<chem_cell_count(s);++c) p[c]=.005*t[c];
    }
    return CHEM_OK;
}
static void forcing(const ChemState *s, const double *y, double *out, void *v)
{
    Model *model=v; size_t cells=chem_cell_count(s);
    for (size_t i=0;i<2;++i) {
        Transfer *t=&model->transfers[i];
        const double *k=t->parameter>=0 ? chem_parameter_values_const(s,(size_t)t->parameter) : NULL;
        for (size_t c=0;c<cells;++c) {
            double flux=(k ? k[c] : t->rate)*y[(size_t)t->src*cells+c];
            out[(size_t)t->src*cells+c]-=flux; out[(size_t)t->dst*cells+c]+=flux;
        }
    }
}
static void jacobian(const ChemState *s, const double *y, ChemLinearSystem *j, void *v)
{
    (void)y; Model *model=v;
    for (size_t i=0;i<2;++i) {
        Transfer *t=&model->transfers[i];
        const double *k=t->parameter>=0 ? chem_parameter_values_const(s,(size_t)t->parameter) : NULL;
        for (size_t c=0;c<j->cells;++c) {
            double rate=k ? k[c] : t->rate;
            j->values[(size_t)t->ss*j->cells+c]+=rate;
            j->values[(size_t)t->ps*j->cells+c]-=rate;
        }
    }
}
static int fixture_init(Fixture *f, size_t cells)
{
    memset(f,0,sizeof(*f));
    f->contexts[0]=(Model){.transfers={{"FO2","STUB1.MODE2.CORGE.FO2",0,0,-1,0,0,.001},
                                      {"STUB1.MODE1.QUUX.BAZ","STUB1.MODE2.QUUX.BAZ",0,0,-1,0,0,.002}},.temperature_parameter=-1};
    f->contexts[1]=(Model){.transfers={{"STUB2.MODE2.CORGE.FO2","STUB2.MODE2.CORGE.BAZ",0,0,0,0,0,0},
                                      {"STUB2.MODE3.QUUX.BAZ","STUB2.MODE3.QUUX.QUX",0,0,1,0,0,0}},.temperature_parameter=1};
    f->models[0]=(ChemExternalModel){.jacobian_pattern=pattern1,.jacobian_count=4,.finalize=finalize,.update=update,
        .add_forcing=forcing,.subtract_jacobian=jacobian,.context=&f->contexts[0],.species=model1_species,.species_count=7};
    f->models[1]=(ChemExternalModel){.jacobian_pattern=pattern2,.jacobian_count=4,.finalize=finalize,.update=update,
        .add_forcing=forcing,.subtract_jacobian=jacobian,.context=&f->contexts[1],.species=model2_species,.species_count=13,
        .parameter_names=parameter_names,.parameter_defaults=parameter_defaults,.parameter_count=2};
    if (chem_mechanism_compose(&f->composition,&base,f->models,2)!=CHEM_OK) return 0;
    f->state=chem_state_create_extended(&f->composition.mechanism,NULL,0,f->models,2,cells);
    if (!f->state) { chem_mechanism_composition_destroy(&f->composition); return 0; }
    return 1;
}
static void fixture_free(Fixture *f)
{ chem_state_free(f->state); chem_mechanism_composition_destroy(&f->composition); }
static int index_of(Fixture *f,const char *name) { return chem_find_species(&f->composition.mechanism,name); }
static double *value(Fixture *f,const char *name)
{ int i=index_of(f,name); return i<0 ? NULL : chem_species_values(f->state,(size_t)i); }
static int includes(void)
{
    Fixture f; if (!fixture_init(&f,1)) return 0;
    int ok=f.composition.mechanism.species_count==22 && f.composition.mechanism.parameter_count==2 &&
           f.contexts[0].finalized==1 && f.contexts[1].finalized==1;
    for (size_t i=0;i<2;++i) ok &= index_of(&f,gas[i].name)>=0;
    for (size_t i=0;i<7;++i) ok &= index_of(&f,model1_species[i].name)>=0;
    for (size_t i=0;i<13;++i) ok &= index_of(&f,model2_species[i].name)>=0;
    for (size_t i=0;i<2;++i) ok &= chem_find_parameter(&f.composition.mechanism,parameter_names[i])>=0;
    fixture_free(&f); return ok;
}
static int update_state(size_t cells)
{
    Fixture f; if (!fixture_init(&f,cells)) return 0;
    const char *names[]={"FO2","BAR","STUB1.MODE1.QUUX.BAZ","STUB1.MODE2.CORGE.FO2",
                         "STUB2.MODE3.CORGE.QUX","STUB2.MODE3.CORGE.BAZ","STUB2.MODE1.NUMBER","STUB2.MODE2.NUMBER"};
    const double initial[]={1.34,2.53,.75,1.23,.42,.33,1000,500};
    double values[8][3]; ChemNamedValues entries[8];
    for (size_t i=0;i<8;++i) {
        for (size_t c=0;c<cells;++c) values[i][c]=initial[i]+.01*c;
        entries[i]=(ChemNamedValues){names[i],values[i],cells};
    }
    int ok=chem_set_concentrations(f.state,entries,8)==CHEM_OK;
    for (size_t i=0;i<8;++i) for (size_t c=0;c<cells;++c) ok &= value(&f,names[i])[c]==values[i][c];
    fixture_free(&f); return ok;
}
static int diagnose(int jac)
{
    Fixture f; if (!fixture_init(&f,1)) return 0;
    value(&f,"FO2")[0]=1; value(&f,"STUB1.MODE1.QUUX.BAZ")[0]=.5;
    double force[22], matrix[22*22];
    int ok=chem_update_rates(f.state)==CHEM_OK && chem_evaluate(f.state,0,force,matrix)==CHEM_OK;
    if (jac) ok &= matrix[0*22+0]==.001 && matrix[6*22+0]==-.001 &&
                   matrix[2*22+2]==.002 && matrix[4*22+2]==-.002;
    else ok &= force[0]==-.001 && force[6]==.001 && force[2]==-.001 && force[4]==.001;
    fixture_free(&f); return ok;
}
static int solve_one(int backward, size_t cells, int two_models, int variant, int doolittle)
{
    Fixture f; if (!fixture_init(&f,cells)) return 0;
    ChemLinearSolver linear=chem_doolittle_linear_solver();
    int ok=!doolittle || chem_set_linear_solver(f.state,&linear,NULL)==CHEM_OK;
    for (size_t c=0;c<cells;++c) {
        value(&f,"FO2")[c]=1+.5*c; value(&f,"BAR")[c]=2+.5*c;
        value(&f,"STUB1.MODE1.QUUX.BAZ")[c]=.5+.2*c; value(&f,"STUB1.MODE2.CORGE.FO2")[c]=.8+.2*c;
        value(&f,"STUB2.MODE1.NUMBER")[c]=1000+100*c; value(&f,"STUB2.MODE2.NUMBER")[c]=500+50*c;
        value(&f,"STUB2.MODE3.CORGE.BAZ")[c]=.2+.05*c;
        if (two_models) {
            value(&f,"STUB2.MODE2.CORGE.FO2")[c]=.4+.2*c; value(&f,"STUB2.MODE2.CORGE.BAZ")[c]=.1+.05*c;
            value(&f,"STUB2.MODE3.QUUX.BAZ")[c]=.2+.1*c; value(&f,"STUB2.MODE3.QUUX.QUX")[c]=.3+.1*c;
            chem_parameter_values(f.state,0)[c]=.01+.005*c; chem_temperatures(f.state)[c]=275+10*c;
        }
    }
    ok &= chem_update_rates(f.state)==CHEM_OK;
    ChemStats stats; int solve_status;
    if (backward) {
        ChemBackwardEuler p=chem_backward_euler_defaults();
        solve_status=chem_solve_backward_euler(f.state,&p,10,&stats); ok &= solve_status==CHEM_OK;
    } else {
        ChemRosenbrock p=variant==0 ? chem_rosenbrock_two_stage() : variant==1 ? chem_rosenbrock_three_stage() :
                         variant==2 ? chem_rosenbrock_four_stage() : variant==3 ? chem_rosenbrock_four_stage_dae() :
                         chem_rosenbrock_six_stage_dae();
        solve_status=chem_solve(f.state,&p,10,&stats); ok &= solve_status==CHEM_OK;
    }
    double tolerance=backward ? (two_models ? .1 : .0002) : (two_models && variant==0 ? .0005 : .00005);
    for (size_t c=0;c<cells;++c) {
        double a=1+.5*c,b=.5+.2*c,x=.8+.2*c, da=a*(1-exp(-.01)), db=b*(1-exp(-.02));
        ok &= fabs(value(&f,"FO2")[c]-(a-da))<=tolerance*a && value(&f,"BAR")[c]==2+.5*c;
        ok &= fabs(value(&f,"STUB1.MODE1.QUUX.BAZ")[c]-(b-db))<=tolerance*b &&
              fabs(value(&f,"STUB1.MODE2.QUUX.BAZ")[c]-db)<=tolerance*b &&
              fabs(value(&f,"STUB1.MODE2.CORGE.FO2")[c]-(x+da))<=tolerance*a;
        ok &= value(&f,"STUB2.MODE1.NUMBER")[c]==1000+100*c && value(&f,"STUB2.MODE2.NUMBER")[c]==500+50*c;
        if (two_models) {
            double d=.4+.2*c,e=.1+.05*c,g=.2+.1*c,h=.3+.1*c;
            double de=d*(1-exp(-(.01+.005*c)*10)), dg=g*(1-exp(-(275.0+10*c)*.005*10));
            ok &= fabs(value(&f,"STUB2.MODE2.CORGE.FO2")[c]-(d-de))<=tolerance*d &&
                  fabs(value(&f,"STUB2.MODE2.CORGE.BAZ")[c]-(e+de))<=tolerance*d &&
                  fabs(value(&f,"STUB2.MODE3.QUUX.BAZ")[c]-(g-dg))<=tolerance*g &&
                  fabs(value(&f,"STUB2.MODE3.QUUX.QUX")[c]-(h+dg))<=tolerance*g;
        }
    }
    if (!ok) fprintf(stderr,"solve b=%d cells=%zu two=%d variant=%d lu=%d status=%d time=%g FO2=%g BAR=%g BAZ1=%g BAZ2=%g FO2p=%g FO2m=%g BAZm=%g BAZ3=%g QUX3=%g numbers=%g,%g\n",
                     backward,cells,two_models,variant,doolittle,solve_status,stats.final_time,value(&f,"FO2")[0],value(&f,"BAR")[0],
                     value(&f,"STUB1.MODE1.QUUX.BAZ")[0],value(&f,"STUB1.MODE2.QUUX.BAZ")[0],value(&f,"STUB1.MODE2.CORGE.FO2")[0],
                     value(&f,"STUB2.MODE2.CORGE.FO2")[0],value(&f,"STUB2.MODE2.CORGE.BAZ")[0],
                     value(&f,"STUB2.MODE3.QUUX.BAZ")[0],value(&f,"STUB2.MODE3.QUUX.QUX")[0],
                     value(&f,"STUB2.MODE1.NUMBER")[0],value(&f,"STUB2.MODE2.NUMBER")[0]);
    fixture_free(&f); return ok;
}
static int solve_variants(int backward,size_t cells,int two_models)
{
    int ok=1;
    for (int variant=0;variant<(backward ? 1 : 5);++variant)
        for (int linear=0;linear<2;++linear) ok &= solve_one(backward,cells,two_models,variant,linear);
    return ok;
}
static int failures;
static void report(const char *prefix,const char *name,int ok)
{
    printf("%s:AerosolModelIntegration.%s %s %s\n",prefix,name,ok ? "PASS" : "FAIL",
           ok ? "reference external model behavior matched" : "external model behavior differed"); failures+=!ok;
}
int main(void)
{
    const char *prefixes[]={"aerosol_model_rosenbrock","aerosol_model_backward_euler"};
    for (int b=0;b<2;++b) {
        const char *p=prefixes[b];
        report(p,"StateIncludesStubAerosolModel",includes());
        report(p,"CanUpdateStateWithStubAerosolModel",update_state(1));
        report(p,"CanUpdateMultiCellStateWithStubAerosolModel",update_state(3));
        report(p,"CanCalculateSingleGridCellForcingWithStubAerosolModel",diagnose(0));
        report(p,"CanCalculateSingleGridCellJacobianWithStubAerosolModel",diagnose(1));
        report(p,"CanSolveSingleGridCellWithStubAerosolModel1",solve_variants(b,1,0));
        report(p,"CanSolveSingleGridCellWithTwoStubAerosolModels",solve_variants(b,1,1));
        report(p,"CanSolveMultiGridCellWithStubAerosolModel1",solve_variants(b,3,0));
        report(p,"CanSolveMultiGridCellWithTwoStubAerosolModels",solve_variants(b,3,1));
    }
    return failures ? 1 : 0;
}
