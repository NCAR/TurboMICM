#include "rosenbrock.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct { size_t n, source, product; double rate; } Reaction;
static void rhs(void *context,const double *y,double *f)
{
    Reaction *r=context; memset(f,0,r->n*sizeof(double));
    f[r->source]-=r->rate*y[r->source]; f[r->product]+=r->rate*y[r->source];
}
static void jac(void *context,const double *y,double *j)
{
    Reaction *r=context; (void)y; memset(j,0,r->n*r->n*sizeof(double));
    j[r->source*r->n+r->source]-=r->rate; j[r->product*r->n+r->source]+=r->rate;
}
static void run(ChemConstraintSet *set, Reaction *reaction, double *initial, size_t steps, double dt, double temperature)
{
    size_t n=set->variables;
    unsigned char pattern[16]={0};
    pattern[reaction->source*n+reaction->source]=1;
    pattern[reaction->product*n+reaction->source]=1;
    chem_constraint_sparsity(set,pattern);
    ChemSparsePlan *plan=chem_sparse_plan_from_pattern(n,pattern);
    assert(plan);
    const double atol[]={1e-3,1e-3,1e-3,1e-3};
    ChemRosenbrockMethod method=chem_rosenbrock_three_stage();
    for (size_t sparse=0;sparse<2;++sparse)
      for (size_t algorithm=0;algorithm<2;++algorithm) {
        ChemRosenbrockOptions options=chem_rosenbrock_defaults(atol);
        options.constraints=set; options.constraint_temperature=temperature;
        options.sparse_plan=sparse ? plan : NULL; options.lu_algorithm=(ChemLuAlgorithm)algorithm;
        double y[4]; memcpy(y,initial,n*sizeof(double));
        for (size_t step=0;step<steps;++step) {
            assert(chem_rosenbrock_solve(n,y,dt,rhs,jac,reaction,&method,&options,NULL)==CHEM_SOLVE_OK);
            for (size_t i=0;i<set->count;++i) {
                const ChemConstraint *c=set->constraints+i;
                double k=c->constant*exp((c->enthalpy/c->gas_constant)*(1/temperature-1/c->reference_temperature));
                double left=k,right=1;
                for (size_t j=0;j<c->reactant_count;++j) left*=pow(y[c->reactants[j].species],c->reactants[j].coefficient);
                for (size_t j=0;j<c->product_count;++j) right*=pow(y[c->products[j].species],c->products[j].coefficient);
                assert(fabs(left-right)<1e-6);
            }
        }
        assert(fabs(y[reaction->source]+y[reaction->product]-initial[reaction->source]-initial[reaction->product])<.01);
      }
    chem_sparse_plan_destroy(plan);
}
int main(void)
{
    const ChemConstraintTerm b[]={{1,1}}, c[]={{2,1}}, d[]={{3,1}};
    ChemConstraint constraints[]={
        {"B_C",CHEM_CONSTRAINT_EQUILIBRIUM,2,b,c,1,1,2,-2400,8.31446261815324,298.15,NULL},
        {"B_D",CHEM_CONSTRAINT_EQUILIBRIUM,3,b,d,1,1,5,-2400,8.31446261815324,298.15,NULL}
    };
    ChemConstraintSet set={3,1,constraints}; Reaction reaction={3,0,1,1};
    double y[]={1,0,0}; run(&set,&reaction,y,100,.001,270);
    const ChemConstraintTerm reordered_b[]={{2,1}}, reordered_c[]={{0,1}};
    constraints[0]=(ChemConstraint){"B_C",CHEM_CONSTRAINT_EQUILIBRIUM,0,reordered_b,reordered_c,1,1,2,-2400,8.31446261815324,298.15,NULL};
    reaction=(Reaction){3,1,2,1}; double reordered[]={0,1,0}; run(&set,&reaction,reordered,100,.001,400);
    constraints[0]=(ChemConstraint){"B_C",CHEM_CONSTRAINT_EQUILIBRIUM,2,b,c,1,1,3,-2400,8.31446261815324,298.15,NULL};
    reaction=(Reaction){3,0,1,1}; set.variables=4; set.count=2; reaction.n=4;
    double y4[]={1,0,0,0}; run(&set,&reaction,y4,100,.001,300);
    const ChemConstraintTerm a2[]={{0,2}}, bproduct[]={{1,1}};
    constraints[0]=(ChemConstraint){"A2_B",CHEM_CONSTRAINT_EQUILIBRIUM,1,a2,bproduct,1,1,10,-2400,8.31446261815324,298.15,NULL};
    set=(ChemConstraintSet){3,1,constraints}; reaction=(Reaction){3,2,0,.5};
    double yn[]={.1,.1,1}; run(&set,&reaction,yn,50,.001,298);
    puts("equilibrium DAE integration PASS"); return 0;
}
