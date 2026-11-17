#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "backward_euler.h"
#include "rates.h"

typedef struct { const ChemMechanism *mechanism; double rates[5]; } Context;
static void forcing(void *v, const double *y, double *out)
{ Context *c=v; chem_rhs(c->mechanism,y,c->rates,out); }
static void jacobian(void *v, const double *y, double *out)
{ Context *c=v; chem_jacobian(c->mechanism,y,c->rates,out); }
static void check_close(double a, double b, double tolerance, int family, int method, size_t cell, size_t step)
{
    if (!(fabs(a-b)<=tolerance)) {
        fprintf(stderr,"family=%d method=%d cell=%zu step=%zu actual=%.17g expected=%.17g tolerance=%g\n",family,method,cell,step,a,b,tolerance);
        assert(0);
    }
}
static void rates_for_family(int family, ChemConditions c, double *actual, double *expected)
{
    ChemRate r[2] = {0};
    double t=c.temperature, m=c.air_density;
    switch(family) {
    case 0: /* Troe */
    case 1: { /* Ternary chemical activation */
        double a=family ? 4e-5 : 4e-11;
        ChemRateType type=family ? CHEM_RATE_TERNARY : CHEM_RATE_TROE;
        r[0]=(ChemRate){type,{a,0,0,1,0,0,.6,1},0};
        r[1]=(ChemRate){type,{1.2e-3,1.6,3,136,5,24,.9,.8},0};
        double low[]={a,1.2e-3*exp(3/t)*pow(t/300,1.6)};
        double high[]={1,136*exp(24/t)*pow(t/300,5)};
        for (int i=0;i<2;++i) {
            double q=low[i]*m/high[i];
            expected[i]=low[i]*(family ? 1 : m)/(1+q)*pow(i ? .9 : .6,(i ? .8 : 1)/((i ? .8 : 1)+pow(log10(q),2)));
        }
        break;
    }
    case 2:
        r[0]=(ChemRate){CHEM_RATE_TUNNELING,{.004,0,0},0};
        r[1]=(ChemRate){CHEM_RATE_TUNNELING,{.00012,167,1e8},0};
        expected[0]=.004; expected[1]=.00012*exp(-167/t+1e8/pow(t,3)); break;
    case 3:
        r[0]=(ChemRate){CHEM_RATE_ARRHENIUS,{.004,0,50,300,0},0};
        r[1]=(ChemRate){CHEM_RATE_ARRHENIUS,{.00012,7,75,50,.5},0};
        expected[0]=.004*exp(50/t);
        expected[1]=.00012*exp(75/t)*pow(t/50,7)*(1+.5*c.pressure); break;
    case 4: {
        r[0]=(ChemRate){CHEM_RATE_ALKOXY,{1e-4,204.3,1e-3,2},0};
        r[1]=(ChemRate){CHEM_RATE_NITRATE,{1e-4,204.3,1e-3,2},0};
        double z0=2e-22*exp(2)*2.45e19, z1=.43*pow(293.0/298,-8);
        double z=z0/(1+z0/z1)*pow(.41,1/(1+pow(log10(z0/z1),2)))*999;
        double a=2e-22*exp(2)*m*6.02214076e23*1e-6, b=.43*pow(t/298,-8);
        double x=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
        expected[0]=1e-4*exp(-204.3/t)*z/(z+x);
        expected[1]=1e-4*exp(-204.3/t)*x/(z+x); break;
    }
    default: {
        double custom[]={.002,.003};
        r[0]=(ChemRate){CHEM_RATE_USER,{1},0}; r[1]=(ChemRate){CHEM_RATE_USER,{1},1};
        expected[0]=.002; expected[1]=.003;
        chem_rates(2,r,c,custom,actual); return;
    }}
    chem_rates(2,r,c,NULL,actual);
}
int main(void)
{
    const size_t reactants[]={0,1}; const ChemProduct products[]={{1,1},{2,1}};
    const ChemReaction reactions[]={{0,1,0,1},{1,1,1,1}};
    const ChemMechanism mechanism={3,2,reactions,reactants,products};
    const double temp[]={272.5,254.7,312.6}, pressure[]={101253.3,100672.5,101319.8};
    const double tolerances_ros[]={1e-10,1e-8,1e-8,1e-9,1e-13,2e-6};
    const double tolerances_two_stage[]={1e-10,1e-8,2e-5,4e-6,1e-10,2e-6};
    const double tolerances_be[]={1e-6,1e-5,1e-3,1e-3,1e-5,1e-3};
    const double atol[]={1e-3,1e-3,1e-3};
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),
        chem_rosenbrock_four_stage(),chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    ChemSparsePlan *plan=chem_sparse_plan_create(&mechanism); assert(plan);
    for (int family=0;family<6;++family) for (size_t cell=0;cell<3;++cell) {
        ChemConditions c={temp[cell],pressure[cell],pressure[cell]/(1.380649e-23*6.02214076e23*temp[cell])};
        Context context={&mechanism,{0,0}}; double expected_rates[2];
        rates_for_family(family,c,context.rates,expected_rates);
        for (int i=0;i<2;++i) assert(fabs(context.rates[i]-expected_rates[i]) <= 1e-9*fabs(expected_rates[i]));
        for (int method=0;method<6;++method) for (int sparse=0;sparse<2;++sparse) for (int algorithm=0;algorithm<2;++algorithm) {
            ChemRosenbrockOptions ro=chem_rosenbrock_defaults(atol);
            ChemBackwardEulerOptions be=chem_backward_euler_defaults(atol);
            ro.sparse_plan=be.sparse_plan=sparse ? plan : NULL;
            ro.lu_algorithm=be.lu_algorithm=(ChemLuAlgorithm)algorithm;
            double a0=1-(double)cell/3, y[]={a0,0,0};
            for (size_t step=1;step<1000;++step) {
                ChemSolveStatus status=method<5 ? chem_rosenbrock_solve(3,y,1,forcing,jacobian,&context,&methods[method],&ro,NULL) :
                    chem_backward_euler_solve(3,y,1,forcing,jacobian,&context,&be,NULL);
                assert(status==CHEM_SOLVE_OK);
                double k1=expected_rates[0],k2=expected_rates[1], t=(double)step;
                double e1=exp(-k1*t),e2=exp(-k2*t);
                double target[]={a0*e1,a0*(k1/(k2-k1))*(e1-e2),a0*(1+(k1*e2-k2*e1)/(k2-k1))};
                for (int i=0;i<3;++i) check_close(y[i],target[i],method==0 ? tolerances_two_stage[family] : (method<5 ? tolerances_ros[family] : tolerances_be[family]),family,method,cell,step);
            }
        }
    }
    chem_sparse_plan_destroy(plan);
    const size_t stiff_reactants[]={0,1,2,0,1};
    const ChemProduct stiff_products[]={{2,1},{2,1},{3,1},{1,1},{0,1}};
    const ChemReaction stiff_reactions[]={{0,1,0,1},{1,1,1,1},{2,1,2,1},{3,1,3,1},{4,1,4,1}};
    const ChemMechanism stiff={4,5,stiff_reactions,stiff_reactants,stiff_products};
    const double stiff_atol[]={1e-3,1e-3,1e-3,1e-3};
    const double stiff_ros[]={1e-5,2e-3,1e-4,2e-5,2e-3,2e-5};
    const double stiff_be[]={1e-5,1e-2,1e-3,1e-3,1e-2,1e-3};
    plan=chem_sparse_plan_create(&stiff); assert(plan);
    for (int family=0;family<6;++family) for (size_t cell=0;cell<3;++cell) {
        ChemConditions c={temp[cell],pressure[cell],pressure[cell]/(1.380649e-23*6.02214076e23*temp[cell])};
        double rates[2], expected_rates[2]; rates_for_family(family,c,rates,expected_rates);
        if (family==3) {
            ChemRate rate={CHEM_RATE_ARRHENIUS,{1.2e-4,1.6,75,50,.5},0};
            rates[1]=chem_rate(&rate,c,NULL);
            expected_rates[1]=1.2e-4*exp(75/c.temperature)*pow(c.temperature/50,1.6)*(1+.5*c.pressure);
            assert(fabs(rates[1]-expected_rates[1]) <= 1e-9*fabs(expected_rates[1]));
        }
        Context context={&stiff,{rates[0],rates[0],rates[1],4e10,3.6e10}};
        for (int method=0;method<6;++method) for (int sparse=0;sparse<2;++sparse) for (int algorithm=0;algorithm<2;++algorithm) {
            ChemRosenbrockOptions ro=chem_rosenbrock_defaults(stiff_atol);
            ChemBackwardEulerOptions be=chem_backward_euler_defaults(stiff_atol);
            ro.sparse_plan=be.sparse_plan=sparse ? plan : NULL;
            ro.lu_algorithm=be.lu_algorithm=(ChemLuAlgorithm)algorithm;
            double a0=1-(double)cell/3, y[]={a0*.5,a0*.5,0,0};
            for (size_t step=1;step<1000;++step) {
                ChemSolveStatus status=method<5 ? chem_rosenbrock_solve(4,y,1,forcing,jacobian,&context,&methods[method],&ro,NULL) :
                    chem_backward_euler_solve(4,y,1,forcing,jacobian,&context,&be,NULL);
                assert(status==CHEM_SOLVE_OK);
                double k1=expected_rates[0],k2=expected_rates[1], t=(double)step;
                double e1=exp(-k1*t),e2=exp(-k2*t);
                double target[]={a0*e1,a0*(k1/(k2-k1))*(e1-e2),a0*(1+(k1*e2-k2*e1)/(k2-k1))};
                double observed[]={y[0]+y[1],y[2],y[3]};
                double tolerance=method<5 ? stiff_ros[family] : stiff_be[family];
                if (family==3 && method==0) tolerance=1e-4;
                if (family==3 && method==4) tolerance=1e-5;
                for (int i=0;i<3;++i) check_close(observed[i],target[i],tolerance,family,method,cell,step);
            }
        }
    }
    chem_sparse_plan_destroy(plan);
    puts("analytical rate-family trajectories PASS");
    return 0;
}
