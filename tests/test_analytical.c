#include "chem.h"
#include "state.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static double falloff(double k0,double ki,double f,double n,double air,int troe)
{
    double q=k0*air/ki, logq=log10(q);
    return k0*(troe?air:1)/(1+q)*pow(f,n/(n+logq*logq));
}
static void setup_rate(int kind,ChemReaction *r,double *parameters)
{
    memset(r,0,2*sizeof(*r)); r[0].reactant_count=r[0].product_count=1;
    r[1].reactant_begin=r[1].product_begin=1; r[1].reactant_count=r[1].product_count=1;
    if(kind==0 || kind==2) {
        r[0].type=r[1].type=kind==0?CHEM_TROE:CHEM_TERNARY;
        double a[]={kind==0?4e-11:4e-5,0,0,1,0,0,.6,1};
        double b[]={1.2e-3,1.6,3,136,5,24,.9,.8};
        memcpy(r[0].p,a,sizeof a); memcpy(r[1].p,b,sizeof b);
    } else if(kind==1) {
        r[0].type=r[1].type=CHEM_USER; r[0].parameter=0;r[1].parameter=1;r[0].p[0]=r[1].p[0]=1;
        parameters[0]=2e-3;parameters[1]=3e-3;
    } else if(kind==3) {
        r[0].type=r[1].type=CHEM_TUNNELING;r[0].p[0]=4e-3;
        r[1].p[0]=1.2e-4;r[1].p[1]=167;r[1].p[2]=1e8;
    } else if(kind==4) {
        r[0].type=r[1].type=CHEM_ARRHENIUS;
        double a[]={4e-3,0,50,300,0},b[]={1.2e-4,7,75,50,.5};
        memcpy(r[0].p,a,sizeof a);memcpy(r[1].p,b,sizeof b);
    } else {
        r[0].type=CHEM_BRANCH_ALKOXY;r[1].type=CHEM_BRANCH_NITRATE;
        double p[]={1e-4,204.3,1e-3,2}; memcpy(r[0].p,p,sizeof p);memcpy(r[1].p,p,sizeof p);
    }
}
static void reference_rate(int kind,ChemConditions c,double *k)
{
    double t=c.temperature,air=c.air_density;
    if(kind==0 || kind==2) {
        k[0]=falloff(kind==0?4e-11:4e-5,1,.6,1,air,kind==0);
        k[1]=falloff(1.2e-3*exp(3/t)*pow(t/300,1.6),136*exp(24/t)*pow(t/300,5),.9,.8,air,kind==0);
    } else if(kind==1) {k[0]=2e-3;k[1]=3e-3;}
    else if(kind==3) {k[0]=4e-3;k[1]=1.2e-4*exp(-167/t+1e8/pow(t,3));}
    else if(kind==4) {k[0]=4e-3*exp(50/t);k[1]=1.2e-4*exp(75/t)*pow(t/50,7)*(1+.5*c.pressure);}
    else {
        double avogadro=6.02214076e23;
        double a=2e-22*exp(2)*2.45e19,b=.43*pow(293.0/298,-8);
        double z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
        a=2e-22*exp(2)*air*avogadro*1e-6;b=.43*pow(t/298,-8);
        double branch=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2))),pre=1e-4*exp(-204.3/t);
        k[0]=pre*z/(z+branch);k[1]=pre*branch/(z+branch);
    }
}
static int run_case(int kind,int backward,int stiff,int method,int algorithm,size_t width)
{
    const double temperature[]={272.5,254.7,312.6},pressure[]={101253.3,100672.5,101319.8};
    const double rtolerance[]={1e-10,2e-6,1e-8,1e-8,1e-9,1e-13};
    const double btolerance[]={1e-6,1e-3,1e-5,1e-3,1e-3,1e-5};
    const double two_stage_tolerance[]={1e-10,2e-6,1e-8,2e-5,4e-6,1e-10};
    double tolerance=backward?btolerance[kind]:method==0?two_stage_tolerance[kind]:rtolerance[kind];
    const double stiff_r[]={1e-5,2e-5,2e-3,1e-4,2e-5,2e-3};
    const double stiff_b[]={1e-5,1e-3,1e-2,1e-3,1e-3,1e-2};
    if(stiff) { tolerance=backward?stiff_b[kind]:stiff_r[kind]; if(!backward && kind==4 && method==0) tolerance=1e-4; if(!backward && kind==4 && method==4) tolerance=1e-5; }
    const char *simple_names[]={"A","B","C"},*stiff_names[]={"A1","A2","B","C"},*parameter_names[]={"photoA","photoB"};
    const char **names=stiff?stiff_names:simple_names;
    ChemReaction rx[5]={0},simple[2];double parameters[2]={0},rates[3][5],expected_rates[3][2];setup_rate(kind,simple,parameters);
    ChemTerm simple_rs[]={{0,1},{1,1}},simple_ps[]={{1,1},{2,1}};
    ChemTerm stiff_rs[]={{0,1},{1,1},{2,1},{0,1},{1,1}},stiff_ps[]={{2,1},{2,1},{3,1},{1,1},{0,1}};
    size_t n=stiff?4:3,nr=stiff?5:2;
    if(stiff) {
        rx[0]=simple[0]; rx[1]=simple[0]; rx[2]=simple[1];
        if(kind==4) rx[2].p[1]=1.6;
        rx[3]=(ChemReaction){.type=CHEM_ARRHENIUS,.p={4e10,0,0,300}};
        rx[4]=(ChemReaction){.type=CHEM_ARRHENIUS,.p={3.6e10,0,0,300}};
        for(size_t i=0;i<5;++i) {rx[i].reactant_begin=i;rx[i].product_begin=i;rx[i].reactant_count=rx[i].product_count=1;}
    } else {rx[0]=simple[0];rx[1]=simple[1];}
    ChemMechanism m={.species_count=n,.reaction_count=nr,.parameter_count=kind==1?2:0,.species_names=names,
        .parameter_names=parameter_names,.reactions=rx,.reactants=stiff?stiff_rs:simple_rs,.products=stiff?stiff_ps:simple_ps};
    ChemStateParameters sp={nr,n,kind==1?2u:0u,names,parameter_names,1e-6,NULL}; ChemState state={0};
    if(chem_state_init(&state,&sp,3,width)) return 0;
    ChemWorkspace *work=chem_workspace_create(&m);if(!work) {chem_state_free(&state);return 0;}
    for(size_t cell=0;cell<3;++cell) {
        double gas=1.380649e-23*6.02214076e23;
        state.conditions[cell]=(ChemConditions){temperature[cell],pressure[cell],pressure[cell]/(gas*temperature[cell]),1};
        *chem_matrix_at(&state.variables,cell,0)=(stiff?.5:1)*(1-(double)cell/3);
        if(stiff) *chem_matrix_at(&state.variables,cell,1)=.5*(1-(double)cell/3);
        if(kind==1) chem_matrix_set_row(&state.custom_parameters,cell,parameters,2);
        reference_rate(kind,state.conditions[cell],expected_rates[cell]);
        if(stiff && kind==4) expected_rates[cell][1]=1.2e-4*exp(75/temperature[cell])*pow(temperature[cell]/50,1.6)*(1+.5*pressure[cell]);
    }
    int ok=!chem_state_update_rate_constants(&state,&m,NULL);
    for(size_t cell=0;cell<3;++cell) for(size_t j=0;j<nr;++j) {
        rates[cell][j]=*chem_matrix_at(&state.rate_constants,cell,j);
        double expected=stiff?(j<2?expected_rates[cell][0]:j==2?expected_rates[cell][1]:j==3?4e10:3.6e10):expected_rates[cell][j];
        if(fabs(rates[cell][j]-expected)>1e-9*fabs(expected)) ok=0;
    }
    ChemRosenbrockParameters rp=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
    ChemBackwardEulerParameters bp=chem_backward_euler_defaults();rp.linear_solver=bp.linear_solver=(ChemLUAlgorithm)algorithm;
    for(int step=1;step<1000 && ok;++step) for(size_t cell=0;cell<3 && ok;++cell) {
        double y[4]; ChemStatistics stats;chem_matrix_get_row(&state.variables,cell,y,n);
        int status=backward?chem_backward_euler_solve(&m,work,&bp,state.conditions[cell],rates[cell],1,y,&stats):chem_rosenbrock_solve(&m,work,&rp,state.conditions[cell],rates[cell],1,y,&stats);
        if(status) {fprintf(stderr,"kind=%d method=%d algorithm=%d width=%zu step=%d cell=%zu status=%d\n",kind,method,algorithm,width,step,cell,status);ok=0;break;} chem_matrix_set_row(&state.variables,cell,y,n);
        double k1=expected_rates[cell][0],k2=expected_rates[cell][1],t=step,initial=1-(double)cell/3;
        double e1=exp(-k1*t),e2=exp(-k2*t);
        double exact[]={initial*e1,initial*(k1/(k2-k1))*(e1-e2),initial*(1+(k1*e2-k2*e1)/(k2-k1))};
        double actual[]={stiff?y[0]+y[1]:y[0],stiff?y[2]:y[1],stiff?y[3]:y[2]};
        for(size_t i=0;i<3;++i) if(!isfinite(actual[i]) || fabs(actual[i]-exact[i])>tolerance) {
            fprintf(stderr,"kind=%d method=%d algorithm=%d width=%zu step=%d cell=%zu species=%zu error=%.4g tolerance=%.4g\n",kind,method,algorithm,width,step,cell,i,fabs(actual[i]-exact[i]),tolerance);
            ok=0;break;
        }
    }
    chem_workspace_free(work);chem_state_free(&state);return ok;
}
int main(void)
{
    const char *names[]={"Troe","Photolysis","TernaryChemicalActivation","Tunneling","Arrhenius","Branched"};
    for(int backward=0;backward<2;++backward) for(int stiff=0;stiff<2;++stiff) for(int kind=0;kind<6;++kind) {
        int ok=1;
        for(int method=0;method<(backward?1:5) && ok;++method) for(int algorithm=0;algorithm<2 && ok;++algorithm)
            for(size_t width=1;width<=4 && ok;++width) ok=run_case(kind,backward,stiff,method,algorithm,width);
        printf("%s:AnalyticalExamples.%s%s %d\n",backward?"analytical_backward_euler":"analytical_rosenbrock_integration",names[kind],stiff?"SuperStiffButAnalytical":"",ok);
        fflush(stdout);
    }
    return 0;
}
