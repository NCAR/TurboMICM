#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "backward_euler.h"
#include "rates.h"

extern const ChemMechanism generated_mechanism;
typedef struct { const ChemMechanism *m; const double *k; } System;
static void rhs(void *v,const double *y,double *out) { System *s=v; chem_rhs(s->m,y,s->k,out); }
static void jac(void *v,const double *y,double *out) { System *s=v; chem_jacobian(s->m,y,s->k,out); }
static void advance(size_t n,double *y,double dt,System *s,ChemRosenbrockMethod *method,ChemRosenbrockOptions *options)
{
    double elapsed=0;
    for (size_t i=0;elapsed<dt && i<100000;++i) {
        ChemSolverStats stats;
        ChemSolveStatus status=chem_rosenbrock_solve(n,y,dt-elapsed,rhs,jac,s,method,options,&stats);
        assert((status==CHEM_SOLVE_OK || status==CHEM_SOLVE_MAX_STEPS) && stats.final_time>0);
        elapsed+=stats.final_time;
    }
    assert(elapsed>=dt);
}
static void robertson(void)
{
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
    double k[]={.04,3e7,1e4}, atol[]={1e-12,1e-12,1e-12};
    System sys={&generated_mechanism,k};
    ChemSparsePlan *plan=chem_sparse_plan_create(sys.m); assert(plan);
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    for (size_t m=0;m<5;++m) for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu) {
        double y[]={1,0,0},target=1,current=0;
        ChemRosenbrockOptions opts=chem_rosenbrock_defaults(atol);
        opts.relative_tolerance=1e-10; opts.sparse_plan=sparse ? plan : NULL; opts.lu_algorithm=(ChemLuAlgorithm)lu;
        for (size_t t=0;t<12;++t) {
            advance(3,y,target-current,&sys,&methods[m],&opts);
            for (size_t i=0;i<3;++i) {
                double error=fabs(y[i]-reference[t][i]);
                assert(error<.3e-3 || error < 1e-6*fabs(reference[t][i]));
            }
            current=target; target*=10;
        }
    }
    chem_sparse_plan_destroy(plan);
}
static void hires(void)
{
    const size_t react[]={0,1,2,3,4,5,6,5,7};
    const ChemProduct prod[]={
        {1,1},{0,.43/8.75},{3,8.32/8.75},{0,8.32/10.03},{3,1.71/10.03},{0,1},
        {2,.43/1.12},{5,.69/1.12},{2,.035/1.745},{5,1.71/1.745},{4,1},
        {4,.43/1.81},{5,.69/1.81},{7,1},{6,1}};
    const ChemReaction reactions[]={
        {0,1,0,1},{1,1,1,2},{2,1,3,2},{3,0,5,1},{3,1,6,2},
        {4,1,8,2},{5,1,10,1},{6,1,11,3},{7,2,14,1}};
    const ChemMechanism mech={8,9,reactions,react,prod};
    const double k[]={1.71,8.75,10.03,.0007,1.12,1.745,.43,1.81,280};
    const double target[2][8]={
        {.000737131257332567,.000144248572631618,.000058887297409676,.001175651343283149,
         .002386356198831330,.006238968252742796,.002849998395185769,.002850001604814231},
        {.000670305503581864,.000130996846986347,.000046862231597733,.001044668020551705,
         .000594883830951485,.001399628833942774,.001014492757718480,.004685507242281520}};
    const double atol[]={1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8};
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    System sys={&mech,k}; ChemSparsePlan *plan=chem_sparse_plan_create(&mech); assert(plan);
    for (size_t m=0;m<6;++m) for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu) {
        double y[]={1,0,0,0,0,0,0,.0057};
        ChemRosenbrockOptions ro=chem_rosenbrock_defaults(atol);
        ChemBackwardEulerOptions be=chem_backward_euler_defaults(atol);
        ro.sparse_plan=be.sparse_plan=sparse ? plan : NULL;
        ro.lu_algorithm=be.lu_algorithm=(ChemLuAlgorithm)lu;
        for (size_t t=0;t<2;++t) {
            double dt=321.8122+100*t;
            if (m<5) advance(8,y,dt,&sys,&methods[m],&ro);
            else assert(chem_backward_euler_solve(8,y,dt,rhs,jac,&sys,&be,NULL)==CHEM_SOLVE_OK);
            double tolerance=m==5 ? .1 : (m==1 || m==2 ? 1e-7 : 1e-6);
            for (size_t i=0;i<8;++i) {
                if (fabs(y[i]-target[t][i])>tolerance)
                    fprintf(stderr,"hires method=%zu sparse=%zu lu=%zu t=%zu i=%zu value=%.17g target=%.17g\n",m,sparse,lu,t,i,y[i],target[t][i]);
                assert(fabs(y[i]-target[t][i])<=tolerance);
            }
        }
    }
    chem_sparse_plan_destroy(plan);
}
static void e5(void)
{
    const size_t reactants[]={0,1,2,0,2,3};
    const ChemProduct products[]={{1,1},{2,1},{4,1},{3,1},{2,1},{5,1}};
    const ChemReaction reactions[]={{0,1,0,2},{1,2,2,1},{3,2,3,1},{5,1,4,2}};
    const ChemMechanism mechanism={6,4,reactions,reactants,products};
    const double rates[]={7.89e-10,1.13e9,1.1e7,1.13e3};
    const double atol[]={1e-7,1e-17,1e-17,1e-17,1e-7,1e-7};
    const double reference[7][4]={
        {1.7599259497677897058e-3,1.3846281519376516449e-11,7.6370038530073911180e-13,1.3082581134075777338e-11},
        {1.6180769999072942552e-3,1.3822370304983735443e-10,8.2515735006838336088e-12,1.2997212954915352082e-10},
        {7.4813208224292220114e-6,2.3734781561205975019e-12,2.2123586689581663654e-12,1.6111948716243113653e-13},
        {4.7150333630401632232e-10,1.8188895860807021729e-14,1.8188812376786725407e-14,8.3484020296321693074e-20},
        {3.1317148329356996037e-14,1.4840957952870064294e-16,1.4840957948345691466e-16,4.5243728279782625194e-26},
        {3.8139035189787091771e-49,1.0192582567660293322e-20,1.0192582567660293322e-20,3.7844935507486221171e-65},
        {0,8.8612334976263783420e-23,8.8612334976263783421e-23,0}};
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    System sys={&mechanism,rates}; ChemSparsePlan *plan=chem_sparse_plan_create(&mechanism); assert(plan);
    for (size_t m=0;m<5;++m) for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu) {
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
        options.relative_tolerance=1e-13; options.sparse_plan=sparse ? plan : NULL; options.lu_algorithm=(ChemLuAlgorithm)lu;
        double y[]={1.76e-3,0,0,0,0,0}, current=0,target=10;
        for (size_t t=0;t<7;++t) {
            advance(6,y,target-current,&sys,&methods[m],&options);
            for (size_t i=0;i<4;++i) {
                double error=fabs(y[i]-reference[t][i]);
                if (!(error<(i==0 ? 1e-6 : 1e-13) || error<1e-10*fabs(reference[t][i])))
                    fprintf(stderr,"E5 method=%zu sparse=%zu lu=%zu t=%zu i=%zu value=%.17g target=%.17g\n",m,sparse,lu,t,i,y[i],reference[t][i]);
                assert(error<(i==0 ? 1e-6 : 1e-13) || error<1e-10*fabs(reference[t][i]));
            }
            current=target; target*=100;
        }
    }
    chem_sparse_plan_destroy(plan);
}
static void oregonator(void)
{
    const double sc=77.27,q=8.375e-6,w=.161,k2=1.6e9,k5=1;
    const double alpha=k5/(sc*w*k2),eta=sc*k5/(w*k2),rho=k5*k5/(w*w*k2),tau=w/k5;
    const double rates[]={k5/(sc*w),k2,sc*k5/w,sc*sc*q*k2/2,k5};
    const double atol[]={alpha*1e-11,eta*1e-11,rho*1e-11,eta,eta};
    const size_t reactants[]={1,0,1,0,0,0,2};
    const ChemProduct products[]={{0,1},{3,1},{2,1},{0,2},{4,1},{1,1}};
    const ChemReaction reactions[]={{0,1,0,1},{1,2,1,1},{3,1,2,2},{4,2,4,1},{6,1,5,1}};
    const ChemMechanism mechanism={5,5,reactions,reactants,products};
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
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    System sys={&mechanism,rates}; ChemSparsePlan *plan=chem_sparse_plan_create(&mechanism); assert(plan);
    for (size_t m=0;m<6;++m) for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu) {
        ChemRosenbrockOptions ro=chem_rosenbrock_defaults(atol); ro.relative_tolerance=1e-9;
        ChemBackwardEulerOptions be=chem_backward_euler_defaults(atol); be.relative_tolerance=1e-9;
        ro.sparse_plan=be.sparse_plan=sparse ? plan : NULL;
        ro.lu_algorithm=be.lu_algorithm=(ChemLuAlgorithm)lu;
        double y[]={alpha,2*eta,3*rho,0,0};
        for (size_t t=0;t<12;++t) {
            if (m<5) advance(5,y,30*tau,&sys,&methods[m],&ro);
            else for (size_t step=0;step<18000;++step)
                assert(chem_backward_euler_solve(5,y,30*tau/18000,rhs,jac,&sys,&be,NULL)==CHEM_SOLVE_OK);
            const double scales[]={alpha,eta,rho};
            for (size_t i=0;i<3;++i) {
                double expected=reference[t][i]*scales[i], relative=fabs(y[i]-expected)/fabs(expected);
                if (!(relative<(m==5 ? .005 : 1e-6)))
                    fprintf(stderr,"oregonator method=%zu sparse=%zu lu=%zu t=%zu i=%zu rel=%g\n",m,sparse,lu,t,i,relative);
                assert(relative<(m==5 ? .005 : 1e-6));
            }
        }
    }
    chem_sparse_plan_destroy(plan);
}
static void surface(void)
{
    const double pi=3.14159265358979323846, gmd=1e-6, gsd=.1, mw=.04607, diffusion=.95e-5, gamma=.02;
    const double radius=gmd/2*exp(2.5*log(gsd)*log(gsd));
    const double number=6/(pi*pow(gmd,3)*exp(4.5*log(gsd)*log(gsd)))*(.002/1000+.003/1000);
    const double temperature=272.5, gas=1.380649e-23*6.02214076e23;
    const double speed=sqrt(8*gas*temperature/(pi*mw));
    const double reference_rate=4*number*pi*radius*radius/(radius/diffusion+4/(speed*gamma));
    const double custom[]={radius,number};
    const ChemRate rate={CHEM_RATE_SURFACE,{diffusion,mw,gamma},0};
    const ChemConditions conditions={temperature,101253.3,0};
    const double k[]={chem_rate(&rate,conditions,custom)};
    assert(fabs(k[0]-reference_rate)<=1e-8);
    const size_t reactants[]={0};
    const ChemProduct products[]={{1,1},{2,.4}};
    const ChemReaction reactions[]={{0,1,0,2}};
    const ChemMechanism mechanism={3,1,reactions,reactants,products};
    const double atol[]={1e-3,1e-3,1e-3};
    const double tolerance[]={1e-2,1e-5,1e-6,1e-5,1e-7,.05};
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    System sys={&mechanism,k}; ChemSparsePlan *plan=chem_sparse_plan_create(&mechanism); assert(plan);
    for (size_t m=0;m<6;++m) for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu) {
        ChemRosenbrockOptions ro=chem_rosenbrock_defaults(atol);
        ChemBackwardEulerOptions be=chem_backward_euler_defaults(atol);
        ro.sparse_plan=be.sparse_plan=sparse ? plan : NULL;
        ro.lu_algorithm=be.lu_algorithm=(ChemLuAlgorithm)lu;
        double y[]={1,0,0},dt=.1/reference_rate;
        for (size_t t=1;t<=10;++t) {
            ChemSolveStatus status=m<5 ? chem_rosenbrock_solve(3,y,dt,rhs,jac,&sys,&methods[m],&ro,NULL) :
                chem_backward_euler_solve(3,y,dt,rhs,jac,&sys,&be,NULL);
            assert(status==CHEM_SOLVE_OK);
            double first=exp(-reference_rate*t*dt), expected[]={first,1-first,.4*(1-first)};
            for (size_t i=0;i<3;++i) {
                double rel=fabs(y[i]-expected[i])/fabs(expected[i]);
                if (!(rel<=tolerance[m])) fprintf(stderr,"surface method=%zu step=%zu i=%zu rel=%g\n",m,t,i,rel);
                assert(rel<=tolerance[m]);
            }
        }
    }
    chem_sparse_plan_destroy(plan);
}
static void terminator(void)
{
    const double pi=3.14159265358979323846,atol[]={1e-20,1e-20};
    const size_t reactants[]={0,1,1};
    const ChemProduct products[]={{1,2},{0,1}};
    const ChemReaction reactions[]={{0,1,0,1},{1,2,1,1}};
    const ChemMechanism mechanism={2,2,reactions,reactants,products};
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    ChemSparsePlan *plan=chem_sparse_plan_create(&mechanism); assert(plan);
    for (size_t cells=1;cells<=4;++cells) for (size_t step=0;step<20;++step)
        for (size_t cell=0;cell<cells;++cell) {
            double lat=pi/180*(cell*(90.0/cells)),lon=step*.3;
            double k1=fmax(0,sin(lat)*sin(pi*20/180)+cos(lat)*cos(pi*20/180)*cos(lon-pi*300/180));
            double k[]={k1,1}, cl2_i=(.13+cell*.17)*1e-6, cl_i=(.1+cell*.08)*1e-10;
            double y[]={cl2_i,cl_i}; System sys={&mechanism,k};
            ChemRosenbrockOptions opts=chem_rosenbrock_defaults(atol); opts.relative_tolerance=1e-8;
            opts.max_steps=100000; opts.sparse_plan=plan;
            assert(chem_rosenbrock_solve(2,y,30,rhs,jac,&sys,&method,&opts,NULL)==CHEM_SOLVE_OK);
            double r=k1/4,cly=cl_i+2*cl2_i,det=sqrt(r*r+2*r*cly),e=exp(-4*det*30);
            double l=det*30>1e-16 ? (1-e)/det/30 : 4;
            double cl_f=-l*(cl_i-det+r)*(cl_i+det+r)/(1+e+30*l*(cl_i+r));
            double cl=cl_i+30*cl_f, cl2=cl2_i-15*cl_f;
            assert(fabs(y[1]-cl)<=fabs(cl)*1e-6+1e-14);
            assert(fabs(y[0]-cl2)<=fabs(cl2)*1e-6+1e-14);
        }
    chem_sparse_plan_destroy(plan);
}
static void integrated_rates(void)
{
    const size_t reactants[]={0,1};
    const ChemProduct products[]={{1,1},{3,1},{2,1},{4,1}};
    const ChemReaction reactions[]={{0,1,0,2},{1,1,2,2}};
    const ChemMechanism mechanism={5,2,reactions,reactants,products};
    const ChemRate rate={CHEM_RATE_ARRHENIUS,{2.15e-11,0,110,300,0},0};
    const ChemConditions conditions={273,1000,0};
    double k[]={chem_rate(&rate,conditions,NULL),1};
    const double atol[]={1e-3,1e-3,1e-3,1e-3,1e-3};
    double y[]={1,0,0,0,0}; System sys={&mechanism,k};
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    ChemRosenbrockOptions opts=chem_rosenbrock_defaults(atol);
    for (size_t t=0;t<100;++t) {
        if (t>50) k[1]=0;
        double prev1=y[3],prev2=y[4];
        assert(chem_rosenbrock_solve(5,y,30,rhs,jac,&sys,&method,&opts,NULL)==CHEM_SOLVE_OK);
        assert(y[3]>=prev1 && y[4]>=prev2);
    }
}
int main(void) { robertson(); hires(); e5(); oregonator(); surface(); terminator(); integrated_rates(); puts("nonlinear benchmark trajectories PASS"); return 0; }
