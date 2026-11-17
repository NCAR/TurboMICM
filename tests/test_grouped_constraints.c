#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "rosenbrock.h"

typedef struct { double rate; } System;
static void rhs(void *v,const double *y,double *f)
{
    double rate=((System *)v)->rate*y[0];
    f[0]=-rate; f[1]=rate; f[2]=0;
}
static void jac(void *v,const double *y,double *j)
{
    (void)y; memset(j,0,9*sizeof(*j));
    double rate=((System *)v)->rate; j[0]=-rate; j[3]=rate;
}
static double equilibrium(void *v,const double *y,double t,const double *parameter)
{
    (void)v; double k=parameter ? *parameter : 2*exp(3000*(1/298.15-1/t));
    return k*y[1]-y[2];
}
static void equilibrium_jac(void *v,const double *y,double t,const double *parameter,double *row)
{
    (void)v; (void)y; double k=parameter ? *parameter : 2*exp(3000*(1/298.15-1/t));
    row[1]-=k; row[2]+=1;
}
int main(void)
{
    const ChemConstraintTerm reactants[]={{1,1}},products[]={{2,1}};
    const size_t dependencies[]={1,2};
    ChemCustomConstraint custom={NULL,dependencies,2,equilibrium,equilibrium_jac};
    ChemConstraint constraint={.name="balance",.type=CHEM_CONSTRAINT_EQUILIBRIUM,.row=2,
        .reactants=reactants,.products=products,.reactant_count=1,.product_count=1,
        .constant=2,.enthalpy=-3000*8.31446261815324,.gas_constant=8.31446261815324,
        .reference_temperature=298.15};
    ChemConstraintSet set={3,1,&constraint};
    const double atol[]={1e-3,1e-3,1e-3};
    const double temperatures[]={298.15,320,350}, overrides[]={2,3,4};
    const unsigned char pattern[]={1,0,0,1,0,0,0,0,0};
    System systems[]={{.1},{.2},{.3}}; void *contexts[]={systems,systems+1,systems+2};
    ChemRosenbrockMethod methods[]={chem_rosenbrock_four_stage_da(),chem_rosenbrock_six_stage_da()};
    for (size_t mode=0;mode<2;++mode) for (size_t override=0;override<2;++override)
      for (size_t m=0;m<2;++m) for (size_t sparse=0;sparse<2;++sparse) for (size_t lu=0;lu<2;++lu) {
        constraint.type=mode ? CHEM_CONSTRAINT_CUSTOM : CHEM_CONSTRAINT_EQUILIBRIUM;
        constraint.custom=mode ? &custom : NULL;
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
        options.max_steps=20000;
        options.constraints=&set; options.constraint_parameters=override ? overrides : NULL;
        options.lu_algorithm=(ChemLuAlgorithm)lu;
        double y[]={1,.1,99, 2,.1,99, 3,.1,99};
        ChemSolveStatus status=chem_rosenbrock_solve_grouped_temperatures(3,3,y,3,rhs,jac,contexts,
            sparse ? pattern : NULL,temperatures,methods+m,&options,NULL);
        if (status!=CHEM_SOLVE_OK) fprintf(stderr,"mode=%zu override=%zu method=%zu sparse=%zu lu=%zu status=%s\n",
            mode,override,m,sparse,lu,chem_solve_status_string(status));
        assert(status==CHEM_SOLVE_OK);
        for (size_t c=0;c<3;++c) {
            double a=(c+1)*exp(-systems[c].rate*3), b=.1+(c+1)-a;
            double k=override ? overrides[c] : 2*exp(3000*(1/298.15-1/temperatures[c]));
            assert(fabs(y[3*c]-a)<1e-5 && fabs(y[3*c+1]-b)<1e-5);
            assert(fabs(y[3*c+2]-k*y[3*c+1])<1e-7);
        }
    }
    {
        const ChemConstraintTerm terms[]={{1,1},{2,1}};
        ChemConstraint balance={.name="total",.type=CHEM_CONSTRAINT_LINEAR,.row=2,
            .reactants=terms,.reactant_count=2,.constant=1};
        ChemConstraintSet linear={3,1,&balance};
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
        options.constraints=&linear; options.max_steps=20000;
        double y[]={.2,.1,99, .3,.2,99, .4,.3,99};
        assert(chem_rosenbrock_solve_grouped(3,3,y,1,rhs,jac,contexts,pattern,
                                            methods,&options,NULL)==CHEM_SOLVE_OK);
        for (size_t c=0;c<3;++c) assert(fabs(y[c*3+1]+y[c*3+2]-1)<1e-8);
    }
    ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol); options.constraints=&set;
    double y[]={1,.1,99,2,.1,99,3,.1,99}, original[9]; memcpy(original,y,sizeof(y));
    const double bad_temperatures[]={298.15,-1,350};
    ChemSolverStats stats={.steps=9};
    assert(chem_rosenbrock_solve_grouped_temperatures(3,3,y,3,rhs,jac,contexts,pattern,
           bad_temperatures,methods,&options,&stats)==CHEM_SOLVE_INVALID);
    assert(!memcmp(y,original,sizeof(y)) && stats.steps==0);
    return 0;
}
