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
    ChemReaction rx[5]={0},simple[2];double parameters[2]={0},rates[15],expected_rates[3][2];setup_rate(kind,simple,parameters);
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
    ChemWorkspace *work[3]={0};
    for(size_t cell=0;cell<3;++cell)if(!(work[cell]=chem_workspace_create(&m))) {
        for(size_t i=0;i<3;++i)chem_workspace_free(work[i]);
        chem_state_free(&state);return 0;
    }
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
        rates[cell*nr+j]=*chem_matrix_at(&state.rate_constants,cell,j);
        double expected=stiff?(j<2?expected_rates[cell][0]:j==2?expected_rates[cell][1]:j==3?4e10:3.6e10):expected_rates[cell][j];
        if(fabs(rates[cell*nr+j]-expected)>1e-9*fabs(expected)) ok=0;
    }
    ChemRosenbrockParameters rp=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
    ChemBackwardEulerParameters bp=chem_backward_euler_defaults();rp.linear_solver=bp.linear_solver=(ChemLUAlgorithm)algorithm;
    for(int step=1;step<1000 && ok;++step) {
        double y[12];ChemStatistics stats;
        for(size_t cell=0;cell<3;++cell)chem_matrix_get_row(&state.variables,cell,y+cell*n,n);
        int status=backward?chem_backward_euler_solve_group(&m,work,3,&bp,state.conditions,rates,1,y,&stats):
                            chem_rosenbrock_solve_group(&m,work,3,&rp,state.conditions,rates,1,y,&stats);
        if(status) {fprintf(stderr,"kind=%d method=%d algorithm=%d width=%zu step=%d status=%d\n",kind,method,algorithm,width,step,status);ok=0;break;}
        for(size_t cell=0;cell<3 && ok;++cell) {
        double *cell_y=y+cell*n;chem_matrix_set_row(&state.variables,cell,cell_y,n);
        double k1=expected_rates[cell][0],k2=expected_rates[cell][1],t=step,initial=1-(double)cell/3;
        double e1=exp(-k1*t),e2=exp(-k2*t);
        double exact[]={initial*e1,initial*(k1/(k2-k1))*(e1-e2),initial*(1+(k1*e2-k2*e1)/(k2-k1))};
        double actual[]={stiff?cell_y[0]+cell_y[1]:cell_y[0],stiff?cell_y[2]:cell_y[1],stiff?cell_y[3]:cell_y[2]};
        for(size_t i=0;i<3;++i) if(!isfinite(actual[i]) || fabs(actual[i]-exact[i])>tolerance) {
            fprintf(stderr,"kind=%d method=%d algorithm=%d width=%zu step=%d cell=%zu species=%zu error=%.4g tolerance=%.4g\n",kind,method,algorithm,width,step,cell,i,fabs(actual[i]-exact[i]),tolerance);
            ok=0;break;
        }
        }
    }
    for(size_t cell=0;cell<3;++cell)chem_workspace_free(work[cell]);
    chem_state_free(&state);return ok;
}
static int robertson_case(int method,int algorithm)
{
    const char *names[]={"A","B","C"},*parameters[]={"r1","r2","r3"};
    const ChemTerm r[]={{0,1},{1,1},{1,1},{1,1},{2,1}},q[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
    const ChemReaction rx[]={{CHEM_USER,0,1,0,1,0,{1}},{CHEM_USER,1,2,1,2,1,{1}},{CHEM_USER,3,2,3,2,2,{1}}};
    const ChemMechanism m={3,3,3,names,parameters,NULL,NULL,rx,r,q,0,NULL};
    const double expected[12][3]={
        {.9664597373330035,.3074626578578675e-4,.3350951640121071e-1},
        {.8413699238414729,.1623390937990473e-4,.1586138422491472},
        {.6172348823960878,.6153591274639123e-5,.3827589640126376},
        {.3368745306607069,.2013702318261393e-5,.6631234556369748},
        {.1073004285378040,.4800166972571660e-6,.8926990914454987},
        {.1786592114209946e-1,.7274751468436319e-7,.9821340061103859},
        {.2031483924973415e-2,.8142277783356159e-8,.9979685079327488},
        {.2076093439016395e-3,.8306077485067610e-9,.9997923898254906},
        {.2082417512179460e-4,.8329841429908955e-10,.9999791757415798},
        {.2083229471647004e-5,.8332935037760723e-11,.9999979167621954},
        {.2083328471883087e-6,.8333315602809495e-12,.9999997916663195},
        {.2083340149701284e-7,.8333360770334744e-13,.9999999791665152}
    };
    ChemWorkspace *w=chem_workspace_create(&m);if(!w)return 0;
    ChemRosenbrockParameters p=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
    p.relative_tolerance=1e-10;p.absolute_tolerance=1e-12;p.linear_solver=(ChemLUAlgorithm)algorithm;
    ChemConditions c={272.5,101253.3,1e6,1};double rates[]={.04,3e7,1e4},y[]={1,0,0},target=1,now=0;int ok=1;
    for(int i=0;i<12 && ok;++i){double duration=target-now,advanced=0;
        while(advanced<duration){ChemStatistics stats;int status=chem_rosenbrock_solve(&m,w,&p,c,rates,duration-advanced,y,&stats);
            if((status && status!=2) || !(stats.final_time>0)){fprintf(stderr,"Robertson method %d algorithm %d step %d status %d\n",method,algorithm,i,status);ok=0;break;}advanced+=stats.final_time;}
        for(int j=0;j<3;++j)if(fabs(y[j]-expected[i][j])>=.3e-3 && fabs(y[j]-expected[i][j])>=1e-6*fabs(expected[i][j])){fprintf(stderr,"Robertson method %d step %d species %d got %.8g expected %.8g\n",method,i,j,y[j],expected[i][j]);ok=0;}
        now=target;target*=10;
    }
    chem_workspace_free(w);return ok;
}
static int surface_case(int backward,int method)
{
    const char *names[]={"foo","bar","baz"},*parameters[]={"radius","number"};
    const ChemTerm reactants[]={{0,1}},products[]={{1,1},{2,.4}};
    const ChemReaction rx[]={{CHEM_SURFACE,0,1,0,2,0,{.95e-5,.04607,.02}}};
    const ChemMechanism m={3,1,2,names,parameters,NULL,NULL,rx,reactants,products,0,NULL};
    const double pi=acos(-1),gmd=1e-6,gsd=.1,temperature=272.5;
    double custom[]={gmd/2*exp(2.5*pow(log(gsd),2)),6/(pi*pow(gmd,3)*exp(4.5*pow(log(gsd),2)))*(.002/1000+.003/1000)};
    const double speed=sqrt(8*(1.380649e-23*6.02214076e23)*temperature/(pi*.04607));
    const double reference_rate=4*custom[1]*pi*custom[0]*custom[0]/(custom[0]/.95e-5+4/(speed*.02));
    const double tolerances[]={1e-2,1e-5,1e-6,1e-5,1e-7};double tol=backward?.05:tolerances[method],rate,y[]={1,0,0};
    ChemConditions c={temperature,101253.3,0,0};ChemWorkspace *w=chem_workspace_create(&m);if(!w)return 0;
    ChemRosenbrockParameters rp=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
    ChemBackwardEulerParameters bp=chem_backward_euler_defaults();int ok=!chem_rate_constants(&m,c,custom,&rate) && fabs(rate-reference_rate)<=1e-8;
    for(int i=1;i<=10 && ok;++i){ChemStatistics stats;int status=backward?chem_backward_euler_solve(&m,w,&bp,c,&rate,.1/reference_rate,y,&stats):chem_rosenbrock_solve(&m,w,&rp,c,&rate,.1/reference_rate,y,&stats);
        double foo=exp(-.1*i),exact[]={foo,1-foo,.4*(1-foo)};if(status)ok=0;
        for(int j=0;j<3;++j)if(fabs(y[j]-exact[j])>tol*fabs(exact[j]))ok=0;
    }
    chem_workspace_free(w);return ok;
}
int main(void)
{
    int robertson=1;for(int method=0;method<5;++method)for(int algorithm=0;algorithm<2;++algorithm)robertson &= robertson_case(method,algorithm);
    printf("analytical_rosenbrock_integration:AnalyticalExamples.Robertson %d\n",robertson);
    int surface=1;for(int method=0;method<5;++method)surface &= surface_case(0,method);
    printf("analytical_rosenbrock_integration:AnalyticalExamples.SurfaceRxn %d\n",surface);
    printf("analytical_backward_euler:AnalyticalExamples.SurfaceRxn %d\n",surface_case(1,0));
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
