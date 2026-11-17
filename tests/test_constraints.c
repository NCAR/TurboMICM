#include "constraint.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static void finite_difference(const ChemConstraintSet *set, const double *y, const double *parameters)
{
    size_t n=set->variables;
    double jac[16]={0}, fp[4]={0}, fm[4]={0}, point[4];
    chem_constraint_subtract_jacobian(set,1,y,NULL,parameters,jac);
    for (size_t col=0;col<n;++col) {
        memcpy(point,y,n*sizeof(double)); point[col]+=1e-6;
        chem_constraint_residuals(set,1,point,NULL,parameters,fp);
        point[col]-=2e-6;
        chem_constraint_residuals(set,1,point,NULL,parameters,fm);
        for (size_t i=0;i<set->count;++i) {
            size_t row=set->constraints[i].row;
            double fd=(fp[row]-fm[row])/2e-6;
            assert(fabs(-jac[row*n+col]-fd) <= 1e-7+1e-7*fabs(fd));
        }
    }
}
int main(void)
{
    const ChemConstraintTerm ab[] = {{0,1},{1,1}};
    const ChemConstraintTerm weighted[] = {{0,2},{1,3},{2,-1}};
    const ChemConstraintTerm fractional[] = {{0,.5},{1,1.5}};
    ChemConstraint c = {"balance",CHEM_CONSTRAINT_LINEAR,1,ab,NULL,2,0,1,0,0,0};
    ChemConstraintSet set = {2,1,&c};
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_OK);
    assert(!strcmp(c.name,"balance") && c.row==1 && c.reactant_count==2 && c.constant==1);
    double state[] = {.3,.7, .5,.6}, forcing[] = {10,20,30,40}, jac[8]={0};
    chem_constraint_residuals(&set,2,state,NULL,NULL,forcing);
    assert(fabs(forcing[1])<1e-10 && fabs(forcing[3]-.1)<1e-10);
    assert(forcing[0]==10 && forcing[2]==30);
    chem_constraint_subtract_jacobian(&set,2,state,NULL,NULL,jac);
    assert(jac[2]==-1 && jac[3]==-1 && jac[6]==-1 && jac[7]==-1);
    finite_difference(&set,state,NULL);
    double y3[] = {1,2,3}, f3[3]={0}, j3[9]={0};
    c=(ChemConstraint){"weighted",CHEM_CONSTRAINT_LINEAR,2,weighted,NULL,3,0,5,0,0,0};
    set=(ChemConstraintSet){3,1,&c};
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_OK);
    assert(c.reactants[0].coefficient==2 && c.reactants[1].coefficient==3 && c.reactants[2].coefficient==-1);
    chem_constraint_residuals(&set,1,y3,NULL,NULL,f3);
    chem_constraint_subtract_jacobian(&set,1,y3,NULL,NULL,j3);
    assert(f3[2]==0 && j3[6]==-2 && j3[7]==-3 && j3[8]==1);
    finite_difference(&set,y3,NULL);
    const ChemConstraintTerm abc[] = {{0,1},{1,1},{2,1}};
    c=(ChemConstraint){"total",CHEM_CONSTRAINT_LINEAR,2,abc,NULL,3,0,10,0,0,0};
    double total[] = {2,3,5};
    chem_constraint_residuals(&set,1,total,NULL,NULL,f3); assert(f3[2]==0);
    total[2]=6; chem_constraint_residuals(&set,1,total,NULL,NULL,f3); assert(f3[2]==1);
    const ChemConstraintTerm diff[]={{0,1},{1,-1}};
    c=(ChemConstraint){"difference",CHEM_CONSTRAINT_LINEAR,1,diff,NULL,2,0,0,0,0,0};
    set=(ChemConstraintSet){2,1,&c};
    double pair[] = {.5,.5};
    chem_constraint_residuals(&set,1,pair,NULL,NULL,forcing); assert(forcing[1]==0);
    pair[0]=.6; chem_constraint_residuals(&set,1,pair,NULL,NULL,forcing); assert(fabs(forcing[1]-.1)<1e-10);
    c=(ChemConstraint){"fractional",CHEM_CONSTRAINT_LINEAR,1,fractional,NULL,2,0,2,0,0,0};
    pair[0]=pair[1]=1; memset(jac,0,sizeof(jac));
    chem_constraint_residuals(&set,1,pair,NULL,NULL,forcing);
    chem_constraint_subtract_jacobian(&set,1,pair,NULL,NULL,jac);
    assert(forcing[1]==0 && jac[2]==-.5 && jac[3]==-1.5);
    pair[0]=23; pair[1]=41; memset(jac,0,sizeof(jac));
    chem_constraint_subtract_jacobian(&set,1,pair,NULL,NULL,jac);
    assert(jac[2]==-.5 && jac[3]==-1.5);

    const ChemConstraintTerm independent[]={{0,2},{1,3}};
    c=(ChemConstraint){"A_B_sum",CHEM_CONSTRAINT_LINEAR,1,independent,NULL,2,0,1,0,0,0};
    pair[0]=.1; pair[1]=.2; memset(jac,0,sizeof(jac));
    chem_constraint_subtract_jacobian(&set,1,pair,NULL,NULL,jac);
    assert(jac[2]==-2 && jac[3]==-3);
    pair[0]=10; pair[1]=20; memset(jac,0,sizeof(jac));
    chem_constraint_subtract_jacobian(&set,1,pair,NULL,NULL,jac);
    assert(jac[2]==-2 && jac[3]==-3);
    const ChemConstraintTerm inputs[]={{0,1},{1,1}}, outputs[]={{2,1}};
    c=(ChemConstraint){"eq",CHEM_CONSTRAINT_EQUILIBRIUM,2,inputs,outputs,2,1,1000,-2400,8.31446261815324,298.15};
    set=(ChemConstraintSet){3,1,&c};
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_OK);
    assert(c.row==2 && c.reactant_count==2 && c.product_count==1 && !strcmp(c.name,"eq"));
    const ChemConstraintTerm single_in[]={{0,1}}, single_out[]={{1,1}};
    ChemConstraint single={"A_B_simple",CHEM_CONSTRAINT_EQUILIBRIUM,1,single_in,single_out,1,1,10,-2400,8.31446261815324,298.15};
    ChemConstraintSet single_set={2,1,&single};
    assert(chem_constraint_validate(&single_set)==CHEM_CONSTRAINT_OK && single.row==1);
    assert(!strcmp(single.name,"A_B_simple") && single.reactants[0].species==0 && single.products[0].species==1);
    double param[]={1000}, eq[] = {.001,.001,.001};
    chem_constraint_residuals(&set,1,eq,NULL,param,f3); assert(fabs(f3[2])<1e-10);
    eq[0]=.02; eq[1]=.03; eq[2]=.05;
    chem_constraint_residuals(&set,1,eq,NULL,param,f3); assert(fabs(f3[2]-.55)<1e-10);
    memset(j3,0,sizeof(j3)); eq[0]=.01; eq[1]=.02;
    chem_constraint_subtract_jacobian(&set,1,eq,NULL,param,j3);
    assert(fabs(j3[6]+20)<1e-8 && fabs(j3[7]+10)<1e-8 && fabs(j3[8]-1)<1e-8);
    finite_difference(&set,eq,param);
    double eq_other[]={.1,.2,5}; finite_difference(&set,eq_other,param);
    double temp[]={298.15};
    chem_constraint_residuals(&set,1,eq,temp,NULL,f3);
    assert(fabs(f3[2]-(1000*.01*.02-.05))<1e-10);
    const ChemConstraintTerm in2[]={{0,2}}, out2[]={{1,1},{2,1}};
    c=(ChemConstraint){"dissociation",CHEM_CONSTRAINT_EQUILIBRIUM,1,in2,out2,1,2,100,0,8.31446261815324,298.15};
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_OK);
    assert(c.row==1 && c.reactant_count==1 && c.product_count==2 && c.reactants[0].coefficient==2);
    double eq2[]={.1,.5,2}, p2[]={100};
    chem_constraint_residuals(&set,1,eq2,NULL,p2,f3); assert(fabs(f3[1])<1e-10);
    eq2[0]=.15; eq2[1]=.8; eq2[2]=1.5; finite_difference(&set,eq2,p2);
    c.constant=0; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_BAD_EQUILIBRIUM);
    c.constant=-1; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_BAD_EQUILIBRIUM);
    c.constant=1; c.reactant_count=0; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_EMPTY_REACTANTS);
    c.reactant_count=1; c.product_count=0; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_EMPTY_PRODUCTS);
    c.product_count=2; ChemConstraintTerm invalid[]={{0,0}}; c.reactants=invalid;
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_BAD_STOICHIOMETRY);
    invalid[0].coefficient=-1; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_BAD_STOICHIOMETRY);
    c.reactants=in2; c.products=invalid; c.product_count=1;
    invalid[0].coefficient=0; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_BAD_STOICHIOMETRY);
    invalid[0].coefficient=-2; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_BAD_STOICHIOMETRY);
    c.products=out2; c.product_count=2;
    unsigned char pattern[9]={0}; chem_constraint_sparsity(&set,pattern);
    assert(pattern[3] && pattern[4] && pattern[5]);
    puts("linear and equilibrium constraints PASS");
    return 0;
}
