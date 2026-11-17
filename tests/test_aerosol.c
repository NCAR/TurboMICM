#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "backward_euler.h"
#include "composite.h"
enum { F, BAR, BAZ1, BAZ2, FO2_AQ, F2, B2, BAZ3, QUX3, QUX_C, BAZ_C, N1, N2,
       UNUSED1, UNUSED2, UNUSED3, UNUSED4, UNUSED5, UNUSED6, UNUSED7, UNUSED8, UNUSED9, NVAR };
typedef struct { size_t from[2],to[2]; double k[2]; int second; } Transfers;
static void update(void *v,ChemConditions conditions,double *parameters)
{
    Transfers *t=v;
    if (t->second) { t->k[0]=parameters[0]; parameters[1]=.005*conditions.temperature; t->k[1]=parameters[1]; }
}
static void add_forcing(void *v,const double *y,double *f)
{
    Transfers *t=v;
    for (size_t i=0;i<2;++i) { double rate=t->k[i]*y[t->from[i]]; f[t->from[i]]-=rate; f[t->to[i]]+=rate; }
}
static void add_jacobian(void *v,const double *y,double *j)
{
    Transfers *t=v; (void)y;
    for (size_t i=0;i<2;++i) { j[t->from[i]*NVAR+t->from[i]]-=t->k[i]; j[t->to[i]*NVAR+t->from[i]]+=t->k[i]; }
}
int main(void)
{
    double atol[NVAR];
    /* This standalone composition test solves cells independently. */
    for (size_t i=0;i<NVAR;++i) atol[i]=.0001;
    ChemRosenbrockMethod methods[]={chem_rosenbrock_two_stage(),chem_rosenbrock_three_stage(),chem_rosenbrock_four_stage(),
        chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    unsigned char pattern[NVAR*NVAR]={0};
    const size_t sources[]={F,BAZ1,F2,BAZ3},products[]={FO2_AQ,BAZ2,B2,QUX3};
    for (size_t i=0;i<4;++i) { pattern[sources[i]*NVAR+sources[i]]=1; pattern[products[i]*NVAR+sources[i]]=1; }
    ChemSparsePlan *plan=chem_sparse_plan_from_pattern(NVAR,pattern); assert(plan);
    for (size_t both=0;both<2;++both) for (size_t method=0;method<6;++method)
      for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu) for (size_t cell=0;cell<3;++cell) {
        Transfers first={{F,BAZ1},{FO2_AQ,BAZ2},{.001,.002},0};
        Transfers second={{F2,BAZ3},{B2,QUX3},{0,0},1};
        ChemOdeComponent components[]={ {&first,update,add_forcing,add_jacobian},{&second,update,add_forcing,add_jacobian} };
        ChemCompositeOde system={NVAR,2,components};
        double parameters[]={both ? .01+.005*cell : 0,0};
        ChemConditions conditions={both ? 275+10*cell : 0,101325,0};
        chem_composite_update(&system,conditions,parameters);
        assert(parameters[1]==.005*conditions.temperature);
        double fo2=1+.5*cell, baz=.5+.2*cell, fo2aq=.8+.2*cell;
        double f2=both ? .4+.2*cell : 0, b2=both ? .1+.05*cell : 0;
        double baz3=both ? .2+.1*cell : 0,qux3=both ? .3+.1*cell : 0;
        double y[NVAR]={fo2,2+.5*cell,baz,0,fo2aq,f2,b2,baz3,qux3,.3+.1*cell,.2+.1*cell,1000+100*cell,500+50*cell};
        if (!both && !cell) {
            double f[NVAR],j[NVAR*NVAR];
            chem_composite_rhs(&system,y,f); chem_composite_jacobian(&system,y,j);
            assert(f[F]==-.001*fo2 && f[FO2_AQ]==.001*fo2 && f[BAZ1]==-.002*baz && f[BAZ2]==.002*baz);
            assert(j[F*NVAR+F]==-.001 && j[FO2_AQ*NVAR+F]==.001 && j[BAZ1*NVAR+BAZ1]==-.002 && j[BAZ2*NVAR+BAZ1]==.002);
        }
        ChemRosenbrockOptions ro=chem_rosenbrock_defaults(atol); ChemBackwardEulerOptions be=chem_backward_euler_defaults(atol);
        ro.sparse_plan=be.sparse_plan=sparse ? plan : NULL; ro.lu_algorithm=be.lu_algorithm=(ChemLuAlgorithm)lu;
        ChemSolveStatus status=method<5 ? chem_rosenbrock_solve(NVAR,y,10,chem_composite_rhs,chem_composite_jacobian,&system,&methods[method],&ro,NULL) :
            chem_backward_euler_solve(NVAR,y,10,chem_composite_rhs,chem_composite_jacobian,&system,&be,NULL);
        assert(status==CHEM_SOLVE_OK);
        double d1=fo2*(1-exp(-.001*10)),d2=baz*(1-exp(-.002*10));
        double d3=f2*(1-exp(-parameters[0]*10)),d4=baz3*(1-exp(-parameters[1]*10));
        double tolerance=method==5 ? (both ? .1 : 2e-4) : (method==0 && both ? 5e-4 : 5e-5);
        assert(fabs(y[F]-(fo2-d1))<=tolerance*fo2 && fabs(y[FO2_AQ]-(fo2aq+d1))<=tolerance*fo2);
        assert(fabs(y[BAZ1]-(baz-d2))<=tolerance*baz && fabs(y[BAZ2]-d2)<=tolerance*baz);
        if (both) {
            assert(fabs(y[F2]-(f2-d3))<=tolerance*f2 && fabs(y[B2]-(b2+d3))<=tolerance*f2);
            if (!(fabs(y[BAZ3]-(baz3-d4))<=tolerance*baz3 && fabs(y[QUX3]-(qux3+d4))<=tolerance*baz3))
                fprintf(stderr,"method=%zu sparse=%zu lu=%zu cell=%zu baz=%.17g expected=%.17g tol=%g\n",method,sparse,lu,cell,y[BAZ3],baz3-d4,tolerance*baz3);
            assert(fabs(y[BAZ3]-(baz3-d4))<=tolerance*baz3 && fabs(y[QUX3]-(qux3+d4))<=tolerance*baz3);
        }
        assert(y[BAR]==2+.5*cell && y[QUX_C]==.3+.1*cell && y[BAZ_C]==.2+.1*cell);
        assert(y[N1]==1000+100*cell && y[N2]==500+50*cell);
    }
    chem_sparse_plan_destroy(plan);
    return 0;
}
