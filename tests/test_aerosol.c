#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "backward_euler.h"
#include "composite.h"
#include "state.h"
enum { F, BAR, BAZ1, BAZ2, FO2_AQ, F2, B2, BAZ3, QUX3, QUX_C, BAZ_C, N1, N2,
       UNUSED1, UNUSED2, UNUSED3, UNUSED4, UNUSED5, UNUSED6, UNUSED7, UNUSED8, UNUSED9, NVAR };
static const char *const names[NVAR]={"FO2","BAR","STUB1.MODE1.QUUX.BAZ","STUB1.MODE2.QUUX.BAZ",
    "STUB1.MODE2.CORGE.FO2","STUB2.MODE2.CORGE.FO2","STUB2.MODE2.CORGE.BAZ",
    "STUB2.MODE3.QUUX.BAZ","STUB2.MODE3.QUUX.QUX","STUB2.MODE3.CORGE.QUX",
    "STUB2.MODE3.CORGE.BAZ","STUB2.MODE1.NUMBER","STUB2.MODE2.NUMBER",
    "STUB1.MODE1.QUUX.QUX","STUB1.MODE2.QUUX.QUX","STUB1.MODE2.CORGE.BAZ",
    "STUB1.MODE2.CORGE.QUX","STUB2.MODE1.QUUX.BAZ","STUB2.MODE1.QUUX.QUX",
    "STUB2.MODE2.CORGE.QUX","STUB2.MODE3.NUMBER","STUB2.MODE3.CORGE.FO2"};
static const char *const parameter_names[]={"STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT",
    "STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT"};
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
static void grouped_solves(const unsigned char *pattern, const ChemRosenbrockMethod *methods)
{
    double atol[NVAR]; for (size_t i=0;i<NVAR;++i) atol[i]=.001;
    for (size_t cells=1;cells<=3;cells+=2) for (size_t both=0;both<2;++both)
      for (size_t method=0;method<6;++method) for (size_t sparse=0;sparse<2;++sparse)
      for (size_t lu=0;lu<2;++lu) {
        Transfers first[3],second[3]; ChemOdeComponent components[3][2];
        ChemCompositeOde systems[3]; void *contexts[3];
        ChemState state={0};
        ChemStateConfig config={cells,NVAR,2,0,names,parameter_names,1e-6,atol};
        assert(chem_state_init(&state,&config)==CHEM_STATE_OK);
        double *y=state.concentrations, custom[3];
        for (size_t c=0;c<cells;++c) custom[c]=both ? .01+.005*c : 0;
        assert(chem_state_set_parameter(&state,parameter_names[0],custom,cells)==CHEM_STATE_OK);
        for (size_t c=0;c<cells;++c) {
            first[c]=(Transfers){{F,BAZ1},{FO2_AQ,BAZ2},{.001,.002},0};
            second[c]=(Transfers){{F2,BAZ3},{B2,QUX3},{0,0},1};
            components[c][0]=(ChemOdeComponent){&first[c],update,add_forcing,add_jacobian};
            components[c][1]=(ChemOdeComponent){&second[c],update,add_forcing,add_jacobian};
            systems[c]=(ChemCompositeOde){NVAR,2,components[c]}; contexts[c]=&systems[c];
            state.conditions[c]=(ChemConditions){both ? 275+10*c : 0,101325,0};
            chem_composite_update(&systems[c],state.conditions[c],state.parameters+2*c);
            double *x=y+c*NVAR;
            x[F]=1+.5*c; x[BAR]=2+.5*c; x[BAZ1]=.5+.2*c; x[FO2_AQ]=.8+.2*c;
            x[F2]=both ? .4+.2*c : 0; x[B2]=both ? .1+.05*c : 0;
            x[BAZ3]=both ? .2+.1*c : 0; x[QUX3]=both ? .3+.1*c : 0;
            x[QUX_C]=.3+.1*c; x[BAZ_C]=.2+.1*c; x[N1]=1000+100*c; x[N2]=500+50*c;
        }
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
        options.lu_algorithm=(ChemLuAlgorithm)lu;
        ChemBackwardEulerOptions be=chem_backward_euler_defaults(atol); be.lu_algorithm=(ChemLuAlgorithm)lu;
        ChemSolveStatus status=method<5 ? chem_rosenbrock_solve_grouped(cells,NVAR,y,10,chem_composite_rhs,
               chem_composite_jacobian,contexts,sparse ? pattern : NULL,&methods[method],&options,NULL) :
            chem_backward_euler_solve_grouped(cells,NVAR,y,10,chem_composite_rhs,
               chem_composite_jacobian,contexts,sparse ? pattern : NULL,&be,NULL);
        assert(status==CHEM_SOLVE_OK);
        for (size_t c=0;c<cells;++c) {
            double *x=y+c*NVAR, fo2=1+.5*c,baz=.5+.2*c,fo2aq=.8+.2*c;
            assert(chem_state_concentration(&state,c,"FO2")==x+F);
            assert(chem_state_concentration(&state,c,"STUB2.MODE3.QUUX.BAZ")==x+BAZ3);
            double delta1=fo2*(1-exp(-.01)),delta2=baz*(1-exp(-.02));
            double tolerance=method==5 ? (both ? .1 : 2e-4) : (method==0 && both ? 5e-4 : 5e-5);
            assert(fabs(x[F]-(fo2-delta1))<=tolerance*fo2 && fabs(x[FO2_AQ]-(fo2aq+delta1))<=tolerance*fo2);
            assert(fabs(x[BAZ1]-(baz-delta2))<=tolerance*baz && fabs(x[BAZ2]-delta2)<=tolerance*baz);
            assert(x[BAR]==2+.5*c && x[N1]==1000+100*c && x[N2]==500+50*c);
            if (both) {
                double f2=.4+.2*c,b2=.1+.05*c,baz3=.2+.1*c,qux3=.3+.1*c;
                double delta3=f2*(1-exp(-(.01+.005*c)*10));
                double delta4=baz3*(1-exp(-(275.0+10*c)*.005*10));
                assert(fabs(x[F2]-(f2-delta3))<=tolerance*f2 && fabs(x[B2]-(b2+delta3))<=tolerance*f2);
                if (fabs(x[BAZ3]-(baz3-delta4))>tolerance*baz3)
                    fprintf(stderr,"group cells=%zu method=%zu cell=%zu error=%g limit=%g\n",cells,method,c,fabs(x[BAZ3]-(baz3-delta4)),tolerance*baz3);
                assert(fabs(x[BAZ3]-(baz3-delta4))<=tolerance*baz3 && fabs(x[QUX3]-(qux3+delta4))<=tolerance*baz3);
            } else assert(x[QUX_C]==.3+.1*c && x[BAZ_C]==.2+.1*c);
        }
        chem_state_destroy(&state);
    }
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
    grouped_solves(pattern,methods);
    ChemSparsePlan *plan=chem_sparse_plan_from_pattern(NVAR,pattern); assert(plan);
    {
        ChemRosenbrockOptions invalid=chem_rosenbrock_defaults(atol);
        double y[NVAR]={0}; void *contexts[]={NULL}; ChemSolverStats stats={.steps=9};
        invalid.sparse_plan=plan;
        assert(chem_rosenbrock_solve_grouped(1,NVAR,y,10,chem_composite_rhs,chem_composite_jacobian,
               contexts,NULL,methods,&invalid,&stats)==CHEM_SOLVE_INVALID && stats.steps==0);
        invalid.sparse_plan=NULL;
        assert(chem_rosenbrock_solve_grouped(0,NVAR,y,10,chem_composite_rhs,chem_composite_jacobian,
               contexts,NULL,methods,&invalid,NULL)==CHEM_SOLVE_INVALID);
    }
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
