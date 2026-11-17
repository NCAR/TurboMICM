#include "chem.h"

#include <math.h>
#include <stdio.h>

static int failures;
static void report(const char *solver, const char *name, int ok)
{
    printf("analytical_%s:AnalyticalExamples.%s %s %s\n",solver,name,ok ? "PASS" : "FAIL",
           ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static ChemRosenbrock method(int i)
{
    switch (i) {
    case 0: return chem_rosenbrock_two_stage();
    case 1: return chem_rosenbrock_three_stage();
    case 2: return chem_rosenbrock_four_stage();
    case 3: return chem_rosenbrock_four_stage_dae();
    default: return chem_rosenbrock_six_stage_dae();
    }
}

/* The reference resumes the same state after an exhausted step budget. */
static int advance(ChemState *s, const ChemRosenbrock *p, const ChemBackwardEuler *be, double duration, double max_step)
{
    double elapsed=0;
    while (elapsed<duration) {
        double dt=fmin(duration-elapsed,max_step);
        ChemStats stats;
        int status=be ? chem_solve_backward_euler(s,be,dt,&stats) : chem_solve(s,p,dt,&stats);
        if ((status!=CHEM_OK && status!=CHEM_TOO_MANY_STEPS) || stats.final_time<=0 || elapsed+stats.final_time<=elapsed) {
            fprintf(stderr,"integration failed: status=%d dt=%.17g elapsed=%.17g advanced=%.17g\n",status,dt,elapsed,stats.final_time);
            return 0;
        }
        elapsed+=stats.final_time;
    }
    return 1;
}

static int robertson(int m, int alternative)
{
    const ChemSpecies species[]={{"A",1,1e-12,0,0},{"B",0,1e-12,0,0},{"C",0,1e-12,0,0}};
    const int reactants[]={0,1,1,1,2};
    const ChemProduct products[]={{1,1},{1,1},{2,1},{0,1},{2,1}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_USER,.a=1,.parameter=0}},
        {1,2,1,2,0,{.kind=CHEM_USER,.a=1,.parameter=1}},
        {3,2,3,2,0,{.kind=CHEM_USER,.a=1,.parameter=2}}};
    const char *names[]={"r1","r2","r3"}; const double rates[]={.04,3e7,1e4};
    const ChemMechanism mechanism={3,3,3,species,reactions,reactants,products,names,rates,272.5,101253.3};
    const double reference[12][3]={
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
        {.2083340149701284e-7,.8333360770334744e-13,.9999999791665152}};
    ChemState *s=chem_state_create(&mechanism,1); if (!s) return 0;
    ChemLinearSolver linear=chem_doolittle_linear_solver();
    int ok=(!alternative || chem_set_linear_solver(s,&linear,NULL)==CHEM_OK) &&
           chem_set_air_density(s,0,1e6)==CHEM_OK && chem_update_rates(s)==CHEM_OK;
    ChemRosenbrock p=method(m); p.relative_tolerance=1e-10;
    double target=1,current=0;
    for (int step=0;step<12 && ok;++step) {
        ok=advance(s,&p,NULL,target-current,target-current);
        for (size_t i=0;i<3 && ok;++i) {
            double actual=chem_species_values(s,i)[0], expected=reference[step][i], error=fabs(actual-expected);
            ok=error<.3e-3 || error/fabs(expected)<1e-6;
        }
        current=target; target*=10;
    }
    chem_state_free(s); return ok;
}

static int e5(int m, int alternative)
{
    const ChemSpecies species[]={{"A1",1.76e-3,1e-7,0,0},{"A2",0,1e-17,0,0},{"A3",0,1e-17,0,0},
                                 {"A4",0,1e-17,0,0},{"A5",0,1e-7,0,0},{"A6",0,1e-7,0,0}};
    const int reactants[]={0,1,2,0,2,3};
    const ChemProduct products[]={{1,1},{2,1},{4,1},{3,1},{2,1},{5,1}};
    const ChemReaction reactions[]={
        {0,1,0,2,0,{.kind=CHEM_USER,.a=1,.parameter=0}},
        {1,2,2,1,0,{.kind=CHEM_USER,.a=1,.parameter=1}},
        {3,2,3,1,0,{.kind=CHEM_USER,.a=1,.parameter=2}},
        {5,1,4,2,0,{.kind=CHEM_USER,.a=1,.parameter=3}}};
    const char *names[]={"r1","r2","r3","r4"}; const double rates[]={7.89e-10,1.13e9,1.1e7,1.13e3};
    const ChemMechanism mechanism={6,4,4,species,reactions,reactants,products,names,rates,298,101325};
    const double reference[7][4]={
        {1.7599259497677897058e-3,1.3846281519376516449e-11,7.6370038530073911180e-13,1.3082581134075777338e-11},
        {1.6180769999072942552e-3,1.3822370304983735443e-10,8.2515735006838336088e-12,1.2997212954915352082e-10},
        {7.4813208224292220114e-6,2.3734781561205975019e-12,2.2123586689581663654e-12,1.6111948716243113653e-13},
        {4.7150333630401632232e-10,1.8188895860807021729e-14,1.8188812376786725407e-14,8.3484020296321693074e-20},
        {3.1317148329356996037e-14,1.4840957952870064294e-16,1.4840957948345691466e-16,4.5243728279782625194e-26},
        {3.8139035189787091771e-49,1.0192582567660293322e-20,1.0192582567660293322e-20,3.7844935507486221171e-65},
        {0,8.8612334976263783420e-23,8.8612334976263783421e-23,0}};
    ChemState *s=chem_state_create(&mechanism,1); if (!s) return 0;
    ChemLinearSolver linear=chem_doolittle_linear_solver();
    int ok=(!alternative || chem_set_linear_solver(s,&linear,NULL)==CHEM_OK) && chem_update_rates(s)==CHEM_OK;
    ChemRosenbrock p=method(m); p.relative_tolerance=1e-13;
    double target=10,current=0;
    for (int step=0;step<7 && ok;++step) {
        ok=advance(s,&p,NULL,target-current,target-current);
        for (size_t i=0;i<4 && ok;++i) {
            double actual=chem_species_values(s,i)[0], expected=reference[step][i], error=fabs(actual-expected);
            ok=error<(i==0 ? 1e-6 : 1e-13) || (expected!=0 && error/fabs(expected)<1e-10);
        }
        current=target; target*=100;
    }
    chem_state_free(s); return ok;
}

static int hires(int backward, int m, int alternative)
{
    const ChemSpecies species[]={{"Y0",1,1e-8,0,0},{"Y1",0,1e-8,0,0},{"Y2",0,1e-8,0,0},
        {"Y3",0,1e-8,0,0},{"Y4",0,1e-8,0,0},{"Y5",0,1e-8,0,0},{"Y6",0,1e-8,0,0},{"Y7",.0057,1e-8,0,0}};
    const int reactants[]={0,1,2,3,4,5,6,5,7};
    const ChemProduct products[]={{1,1},{0,.43/8.75},{3,8.32/8.75},{0,8.32/10.03},{3,1.71/10.03},
        {0,1},{2,.43/1.12},{5,.69/1.12},{2,.035/1.745},{5,1.71/1.745},{4,1},
        {4,.43/1.81},{5,.69/1.81},{7,1},{6,1}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_USER,.a=1,.parameter=0}}, {1,1,1,2,0,{.kind=CHEM_USER,.a=1,.parameter=1}},
        {2,1,3,2,0,{.kind=CHEM_USER,.a=1,.parameter=2}}, {3,0,5,1,0,{.kind=CHEM_USER,.a=1,.parameter=3}},
        {3,1,6,2,0,{.kind=CHEM_USER,.a=1,.parameter=4}}, {4,1,8,2,0,{.kind=CHEM_USER,.a=1,.parameter=5}},
        {5,1,10,1,0,{.kind=CHEM_USER,.a=1,.parameter=6}}, {6,1,11,3,0,{.kind=CHEM_USER,.a=1,.parameter=7}},
        {7,2,14,1,0,{.kind=CHEM_USER,.a=1,.parameter=8}}};
    const char *names[]={"r1","r2","r3","r4","r5","r6","r7","r8","r9"};
    const double rates[]={1.71,8.75,10.03,.0007,1.12,1.745,.43,1.81,280};
    const ChemMechanism mechanism={8,9,9,species,reactions,reactants,products,names,rates,298,101325};
    const double reference[2][8]={
        {.000737131257332567,.000144248572631618,.000058887297409676,.001175651343283149,
         .002386356198831330,.006238968252742796,.002849998395185769,.002850001604814231},
        {.000670305503581864,.000130996846986347,.000046862231597733,.001044668020551705,
         .000594883830951485,.001399628833942774,.001014492757718480,.004685507242281520}};
    ChemState *s=chem_state_create(&mechanism,1); if (!s) return 0;
    ChemLinearSolver linear=chem_doolittle_linear_solver();
    int ok=(!alternative || chem_set_linear_solver(s,&linear,NULL)==CHEM_OK) && chem_update_rates(s)==CHEM_OK;
    ChemRosenbrock p=method(m); p.relative_tolerance=1e-6;
    ChemBackwardEuler be=chem_backward_euler_defaults();
    double tol=backward ? .1 : (m==0 || m>=3 ? 1e-6 : 1e-7), dt=321.8122;
    for (int step=0;step<2 && ok;++step) {
        ok=advance(s,&p,backward?&be:NULL,dt,dt);
        for (size_t i=0;i<8 && ok;++i)
            ok=fabs(chem_species_values(s,i)[0]-reference[step][i])<=tol;
        dt+=100;
    }
    chem_state_free(s); return ok;
}

static int oregonator(int backward, int m, int alternative)
{
    const double sc=77.27, qc=8.375e-6, wc=.161, k2=1.6e9, k5=1;
    const double alpha=k5/(sc*wc*k2), eta=sc*k5/(wc*k2), rho=k5*k5/(wc*wc*k2);
    const ChemSpecies species[]={{"X",alpha,alpha*1e-11,0,0},{"Y",2*eta,eta*1e-11,0,0},
                                 {"Z",3*rho,rho*1e-11,0,0},{"P",0,eta,0,0},{"Q",0,eta,0,0}};
    const int reactants[]={1,0,1,0,0,0,2};
    const ChemProduct products[]={{0,1},{3,1},{2,1},{0,2},{4,1},{1,1}};
    const ChemReaction reactions[]={
        {0,1,0,1,0,{.kind=CHEM_USER,.a=1,.parameter=0}}, {1,2,1,1,0,{.kind=CHEM_USER,.a=1,.parameter=1}},
        {3,1,2,2,0,{.kind=CHEM_USER,.a=1,.parameter=2}}, {4,2,4,1,0,{.kind=CHEM_USER,.a=1,.parameter=3}},
        {6,1,5,1,0,{.kind=CHEM_USER,.a=1,.parameter=4}}};
    const char *names[]={"r1","r2","r3","r4","r5"};
    const double rates[]={k5/(sc*wc),k2,sc*k5/wc,sc*sc*qc*k2/2,k5};
    const ChemMechanism mechanism={5,5,5,species,reactions,reactants,products,names,rates,298,101325};
    const double reference[12][3]={
        {1.000661467180497,1512.778937348249,10358.54312767229},
        {1.000874625199626,1144.336972384497,83.72149966624639},
        {1.001890368438751,529.9926232295553,1.662279579042420},
        {1.004118022612645,243.8326079910346,1.008822224048647},
        {1.008995416634061,112.1664388662539,1.007783229065319},
        {1.019763472537298,51.59761322947535,1.016985778956374},
        {1.043985088527474,23.73442027531524,1.037691843544522},
        {1.100849071667922,10.91533805469020,1.085831969810860},
        {1.249102130020572,5.013945178605446,1.208326626237875},
        {1.779724751937019,2.281852385542403,1.613754023671725},
        {1.000889326903503,1125.438585746596,16410.49483777168},
        {1.000814870318523,1228.178521549889,132.0554942846513}};
    const double scales[]={alpha,eta,rho}, dt=30*wc;
    ChemState *s=chem_state_create(&mechanism,1); if (!s) return 0;
    ChemLinearSolver linear=chem_doolittle_linear_solver();
    int ok=(!alternative || chem_set_linear_solver(s,&linear,NULL)==CHEM_OK) && chem_update_rates(s)==CHEM_OK;
    ChemRosenbrock p=method(m); p.relative_tolerance=1e-9;
    ChemBackwardEuler be=chem_backward_euler_defaults(); be.relative_tolerance=1e-9;
    double tol=backward ? 5e-3 : 1e-6;
    for (int step=0;step<12 && ok;++step) {
        ok=advance(s,&p,backward?&be:NULL,dt,backward ? dt/18000 : dt);
        for (size_t i=0;i<3 && ok;++i) {
            double actual=chem_species_values(s,i)[0], expected=reference[step][i]*scales[i];
            if (fabs(actual-expected)>=tol*fabs(expected)) {
                fprintf(stderr,"oregonator step=%d species=%zu actual=%.17g expected=%.17g tolerance=%g method=%d backward=%d\n",
                        step,i,actual,expected,tol,m,backward);
                ok=0;
            }
        }
    }
    chem_state_free(s); return ok;
}

int main(void)
{
    int r=1, e=1, h=1, o=1;
    for (int m=0;m<5;++m) {
        r &= robertson(m,1); e &= e5(m,1); h &= hires(0,m,1); o &= oregonator(0,m,1);
    }
    r &= robertson(1,0); e &= e5(1,0); h &= hires(0,1,0); o &= oregonator(0,1,0);
    report("rosenbrock_integration","Robertson",r); report("rosenbrock_integration","E5",e);
    report("rosenbrock_integration","HIRES",h); report("rosenbrock_integration","Oregonator",o);
    h=hires(1,1,0) && hires(1,1,1); o=oregonator(1,1,0) && oregonator(1,1,1);
    report("backward_euler","HIRES",h); report("backward_euler","Oregonator",o);
    return failures ? 1 : 0;
}
