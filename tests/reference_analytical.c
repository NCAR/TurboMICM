#include "chem.h"

#include <math.h>
#include <stdio.h>

/* The reference tests use three distinct cells and 999 one-second solves. */
enum { CELLS=3, STEPS=999, KINDS=6 };
static const double temperatures[]={272.5,254.7,312.6};
static const double pressures[]={101253.3,100672.5,101319.8};
static const char *const names[]={"Troe","Photolysis","TernaryChemicalActivation","Tunneling","Arrhenius","Branched"};
static int failures;

static void report(int backward, int kind, int stiff, int ok)
{
    printf("analytical_%s:AnalyticalExamples.%s%s %s %s\n",
           backward ? "backward_euler" : "rosenbrock_integration", names[kind],
           stiff ? "SuperStiffButAnalytical" : "", ok ? "PASS" : "FAIL",
           ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}

static void rate_descriptors(int kind, int stiff, ChemRate *first, ChemRate *second)
{
    ChemRate a, b;
    if (kind==0 || kind==2) {
        a=chem_rate_defaults(kind==0 ? CHEM_TROE : CHEM_TERNARY);
        b=a;
        a.k0_a=kind==0 ? 4e-11 : 4e-5;
        b.k0_a=1.2e-3; b.k0_b=1.6; b.k0_c=3;
        b.kinf_a=136; b.kinf_b=5; b.kinf_c=24; b.fc=.9; b.n=.8;
    } else if (kind==1) {
        a=chem_rate_defaults(CHEM_USER); b=a;
        a.parameter=0; b.parameter=stiff ? 2 : 1;
    } else if (kind==3) {
        a=chem_rate_defaults(CHEM_TUNNELING); b=a;
        a.a=4e-3; b.a=1.2e-4; b.b=167; b.c=1e8;
    } else if (kind==4) {
        a=chem_rate_defaults(CHEM_ARRHENIUS); b=a;
        a.a=4e-3; a.c=50;
        b.a=1.2e-4; b.b=stiff ? 1.6 : 7; b.c=75; b.d=50; b.e=.5;
    } else {
        a=chem_rate_defaults(CHEM_BRANCHED);
        a.branch_x=1e-4; a.branch_y=204.3; a.branch_a0=1e-3; a.branch_n=2;
        b=a; b.nitrate_branch=1;
    }
    *first=a; *second=b;
}

/* Closed forms are deliberately independent of the implementation's rate
 * evaluator. These are also the rate constants used by the reference ODE. */
static void expected_rates(int kind, int stiff, double temperature, double pressure, double *k1, double *k2)
{
    double density=pressure/((1.380649e-23*6.02214076e23)*temperature);
    if (kind==0 || kind==2) {
        double a=kind==0 ? 4e-11 : 4e-5, ratio=a*density;
        *k1=a*(kind==0 ? density : 1)/(1+ratio)*pow(.6,1/(1+pow(log10(ratio),2)));
        double low=1.2e-3*exp(3/temperature)*pow(temperature/300,1.6);
        double high=136*exp(24/temperature)*pow(temperature/300,5);
        ratio=low*density/high;
        *k2=low*(kind==0 ? density : 1)/(1+ratio)*pow(.9,1/(1+pow(log10(ratio),2)/.8));
    } else if (kind==1) {
        *k1=2e-3; *k2=3e-3;
    } else if (kind==3) {
        *k1=4e-3; *k2=1.2e-4*exp(-167/temperature+1e8/pow(temperature,3));
    } else if (kind==4) {
        *k1=4e-3*exp(50/temperature);
        *k2=1.2e-4*exp(75/temperature)*pow(temperature/50,stiff ? 1.6 : 7)*(1+.5*pressure);
    } else {
        double a=2e-22*exp(2)*2.45e19;
        double b=.43*pow(293.0/298,-8);
        double z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
        a=2e-22*exp(2)*density*6.02214076e23*1e-6;
        b=.43*pow(temperature/298,-8);
        double branching=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
        double total=1e-4*exp(-204.3/temperature);
        *k1=total*(z/(z+branching)); *k2=total*(branching/(z+branching));
    }
}

static double tolerance(int kind, int stiff, int backward, int method)
{
    static const double normal_r[]={1e-10,2e-6,1e-8,1e-8,1e-9,1e-13};
    static const double normal_b[]={1e-6,1e-3,1e-5,1e-3,1e-3,1e-5};
    static const double stiff_r[]={1e-5,2e-5,2e-3,1e-4,2e-5,2e-3};
    static const double stiff_b[]={1e-5,1e-3,1e-2,1e-3,1e-3,1e-2};
    if (backward) return stiff ? stiff_b[kind] : normal_b[kind];
    if (!stiff && method==0) {
        if (kind==3) return 2e-5;
        if (kind==4) return 4e-6;
        if (kind==5) return 1e-10;
    }
    if (stiff && kind==4 && method==0) return 1e-4;
    if (stiff && kind==4 && method==4) return 1e-5;
    return stiff ? stiff_r[kind] : normal_r[kind];
}

static int check_integration(int kind, int stiff, int backward, int method, int linear_solver, int width)
{
    ChemSpecies species[]={{"A1",0,1e-3,0,0},{"A2",0,1e-3,0,0},
                           {"B",0,1e-3,0,0},{"C",0,1e-3,0,0}};
    int reactants[]={0,1,2,0,1};
    ChemProduct products[]={{2,1},{2,1},{3,1},{1,1},{0,1}};
    ChemRate first, second; rate_descriptors(kind,stiff,&first,&second);
    ChemRate exchange=chem_rate_defaults(CHEM_ARRHENIUS); exchange.a=4e10;
    ChemRate back=exchange; back.a=.9*4e10;
    ChemReaction reactions[]={
        {0,1,0,1,0,first},{1,1,1,1,0,first},{2,1,2,1,0,second},
        {3,1,3,1,0,exchange},{4,1,4,1,0,back}};
    if (kind==1 && stiff) reactions[1].rate.parameter=1;
    if (!stiff) {
        species[1]=species[2]; species[2]=species[3];
        reactants[1]=1; products[0].species=1; products[1].species=2;
        reactions[1].rate=second;
    }
    const char *params[]={"photoA1B","photoA2B","photoB"};
    const double defaults_stiff[]={2e-3,2e-3,3e-3}, defaults_normal[]={2e-3,3e-3};
    ChemMechanism mechanism={(size_t)(stiff?4:3),(size_t)(stiff?5:2),(size_t)(kind==1 ? (stiff?3:2) : 0),
        species,reactions,reactants,products,params,stiff?defaults_stiff:defaults_normal,298,101325};
    ChemState *groups[CELLS]={0};
    int group_count=(CELLS+width-1)/width, ok=1;
    ChemLinearSolver alternate=chem_doolittle_linear_solver();
    ChemRosenbrock methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),chem_rosenbrock_four_stage(),
                              chem_rosenbrock_four_stage_dae(),chem_rosenbrock_six_stage_dae()};
    ChemBackwardEuler be=chem_backward_euler_defaults();
    for (int g=0;g<group_count && ok;++g) {
        int begin=g*width, count=CELLS-begin<width ? CELLS-begin : width;
        groups[g]=chem_state_create(&mechanism,(size_t)count);
        if (!groups[g]) { ok=0; break; }
        if (linear_solver) ok=chem_set_linear_solver(groups[g],&alternate,NULL)==CHEM_OK;
        for (int c=0;c<count;++c) {
            int cell=begin+c; double initial=1-(double)cell/CELLS;
            chem_temperatures(groups[g])[c]=temperatures[cell]; chem_pressures(groups[g])[c]=pressures[cell];
            chem_species_values(groups[g],0)[c]=initial*(stiff?.5:1);
            if (stiff) chem_species_values(groups[g],1)[c]=initial*.5;
        }
    }
    double tol=tolerance(kind,stiff,backward,method);
    for (int step=1;step<=STEPS && ok;++step) for (int g=0;g<group_count && ok;++g) {
        ChemState *s=groups[g]; int begin=g*width, count=(int)chem_cell_count(s);
        ok=chem_update_rates(s)==CHEM_OK;
        for (int c=0;c<count && ok && !stiff;++c) {
            int cell=begin+c; double k1,k2; expected_rates(kind,stiff,temperatures[cell],pressures[cell],&k1,&k2);
            ok=fabs(chem_rate_constants(s,0)[c]-k1)<=1e-9*fabs(k1) &&
               fabs(chem_rate_constants(s,1)[c]-k2)<=1e-9*fabs(k2);
        }
        ChemStats stats;
        int status=backward ? chem_solve_backward_euler(s,&be,1,&stats) : chem_solve(s,&methods[method],1,&stats);
        ok &= status==CHEM_OK;
        for (int c=0;c<count && ok;++c) {
            int cell=begin+c; double k1,k2; expected_rates(kind,stiff,temperatures[cell],pressures[cell],&k1,&k2);
            double a0=1-(double)cell/CELLS, e1=exp(-k1*step), e2=exp(-k2*step);
            double expected[]={a0*e1,a0*(k1/(k2-k1))*(e1-e2),a0*(1+(k1*e2-k2*e1)/(k2-k1))};
            double actual[]={chem_species_values(s,0)[c]+(stiff?chem_species_values(s,1)[c]:0),
                             chem_species_values(s,stiff?2:1)[c],chem_species_values(s,stiff?3:2)[c]};
            for (int i=0;i<3;++i) if (!isfinite(actual[i]) || fabs(actual[i]-expected[i])>tol) {
                fprintf(stderr,"analytical %s stiff=%d backward=%d method=%d linear=%d width=%d step=%d cell=%d species=%d got %.17g expected %.17g tolerance %g status %d\n",
                        names[kind],stiff,backward,method,linear_solver,width,step,cell,i,actual[i],expected[i],tol,status);
                fprintf(stderr,"rates: actual %.17g %.17g expected %.17g %.17g\n",chem_rate_constants(s,0)[c],chem_rate_constants(s,stiff?2:1)[c],k1,k2);
                ok=0; break;
            }
        }
    }
    for (int g=0;g<group_count;++g) chem_state_free(groups[g]);
    return ok;
}

static int check_surface(int backward, int method)
{
    const double pi=3.14159265358979323846, diameter=1e-6, dispersion=.1;
    const double radius=diameter/2*exp(2.5*log(dispersion)*log(dispersion));
    const double number=6/(pi*pow(diameter,3)*exp(4.5*log(dispersion)*log(dispersion)))*(.002/1000+.003/1000);
    const ChemSpecies species[]={{"foo",1,1e-3,.04607,.95e-5},{"bar",0,1e-3,0,0},{"baz",0,1e-3,0,0}};
    const int reactants[]={0}; const ChemProduct products[]={{1,1},{2,.4}};
    const ChemReaction reaction={0,1,0,2,0,{.kind=CHEM_SURFACE,.parameter=0,.molecular_weight=.04607,
                                           .diffusion=.95e-5,.probability=.02}};
    const char *parameters[]={"radius","number"}; const double defaults[]={radius,number};
    const ChemMechanism m={3,1,2,species,&reaction,reactants,products,parameters,defaults,272.5,101253.3};
    ChemState *s=chem_state_create(&m,1);
    if (!s) return 0;
    const double speed=sqrt(8*(1.380649e-23*6.02214076e23)/(pi*.04607)*272.5);
    const double k=4*number*pi*radius*radius/(radius/.95e-5+4/(speed*.02)), dt=.1/k;
    const double tolerances[]={1e-2,1e-5,1e-6,1e-5,1e-7};
    double tol=backward ? .05 : tolerances[method];
    ChemRosenbrock methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),chem_rosenbrock_four_stage(),
                              chem_rosenbrock_four_stage_dae(),chem_rosenbrock_six_stage_dae()};
    ChemBackwardEuler be=chem_backward_euler_defaults();
    int ok=1;
    for (int step=1;step<=10 && ok;++step) {
        ChemStats stats;
        ok=chem_update_rates(s)==CHEM_OK && fabs(chem_rate_constants(s,0)[0]-k)<=1e-8;
        int status=backward ? chem_solve_backward_euler(s,&be,dt,&stats) : chem_solve(s,&methods[method],dt,&stats);
        ok &= status==CHEM_OK;
        double foo=exp(-k*step*dt), expected[]={foo,1-foo,.4*(1-foo)};
        for (size_t i=0;i<3;++i)
            ok &= fabs(chem_species_values(s,i)[0]-expected[i])<=tol*fabs(expected[i]);
    }
    chem_state_free(s); return ok;
}

int main(void)
{
    for (int backward=0;backward<2;++backward) for (int kind=0;kind<KINDS;++kind) for (int stiff=0;stiff<2;++stiff) {
        int ok=1;
        for (int method=0;method<(backward?1:5);++method)
            ok &= check_integration(kind,stiff,backward,method,1,3);
        for (int linear=0;linear<2;++linear) for (int width=1;width<=4;++width)
            ok &= check_integration(kind,stiff,backward,backward?0:1,linear,width);
        report(backward,kind,stiff,ok);
    }
    int surface=1;
    for (int method=0;method<5;++method) surface &= check_surface(0,method);
    printf("analytical_rosenbrock_integration:AnalyticalExamples.SurfaceRxn %s %s\n",
           surface ? "PASS" : "FAIL",surface ? "reference quantities matched" : "reference quantity differed");
    failures+=!surface;
    surface=check_surface(1,0);
    printf("analytical_backward_euler:AnalyticalExamples.SurfaceRxn %s %s\n",
           surface ? "PASS" : "FAIL",surface ? "reference quantities matched" : "reference quantity differed");
    failures+=!surface;
    return failures ? 1 : 0;
}
