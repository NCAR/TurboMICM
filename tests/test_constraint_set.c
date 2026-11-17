#include "constraint.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define CLOSE(a,b) assert(fabs((a)-(b)) < 1e-10)
#define EQ(row,rin,rnum,pout,pnum) {"eq",CHEM_CONSTRAINT_EQUILIBRIUM,row,rin,pout,rnum,pnum,.033,-24000,8.31446261815324,298.15,NULL}
int main(void)
{
    const ChemConstraintTerm a[]={{0,1}}, b[]={{1,1}}, c[]={{2,1}}, d[]={{3,1}};
    const ChemConstraintTerm ab[]={{0,1},{1,1}}, cd[]={{2,1},{3,1}};
    ChemConstraint constraints[]={EQ(2,ab,2,c,1),EQ(4,d,1,((const ChemConstraintTerm[]){{4,1}}),1)};
    ChemConstraintSet set={5,2,constraints};
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_OK && set.count==2);
    unsigned char pattern[25]={0}; chem_constraint_sparsity(&set,pattern);
    assert(pattern[10] && pattern[11] && pattern[12] && pattern[23] && pattern[24]);
    int count=0; for (size_t i=0;i<25;++i) count+=!!pattern[i]; assert(count==5);
    constraints[1].row=2; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_DUPLICATE_ROW);
    constraints[1].row=5; assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_UNKNOWN_SPECIES);
    constraints[1].row=4;
    ChemConstraintSet empty={3,0,NULL}; double untouched[]={1,1,1,1}, dummy[]={.1,.2,.3};
    assert(chem_constraint_validate(&empty)==CHEM_CONSTRAINT_OK);
    chem_constraint_residuals(&empty,1,dummy,NULL,NULL,untouched);
    for (size_t i=0;i<4;++i) assert(untouched[i]==1);
    ChemConstraintTerm unknown[]={{5,1}}; constraints[0].reactants=unknown; constraints[0].reactant_count=1;
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_UNKNOWN_SPECIES);

    constraints[0]=(ChemConstraint)EQ(2,ab,2,c,1); set=(ChemConstraintSet){3,1,constraints};
    double test_state[]={.2,.4,.6, .1,.3,.7}, test_params[]={.033,.033}, test_f[6]={0};
    chem_constraint_residuals(&set,2,test_state,NULL,test_params,test_f);
    CLOSE(test_f[2],-.59736); CLOSE(test_f[5],-.69901);
    double state[]={.01,.02,.05, .03,.01,.2, .001,.002,.004};
    double parameters[]={.033,.033,.033}, forcing[9]={0}, jac[27]={0};
    chem_constraint_residuals(&set,3,state,NULL,parameters,forcing);
    chem_constraint_subtract_jacobian(&set,3,state,NULL,parameters,jac);
    for (size_t i=0;i<3;++i) {
        CLOSE(forcing[3*i+2],.033*state[3*i]*state[3*i+1]-state[3*i+2]);
        CLOSE(jac[9*i+6],-.033*state[3*i+1]);
        CLOSE(jac[9*i+7],-.033*state[3*i]); CLOSE(jac[9*i+8],1);
    }
    constraints[0]=(ChemConstraint)EQ(2,b,1,c,1);
    memset(pattern,0,sizeof(pattern)); chem_constraint_sparsity(&set,pattern);
    assert(pattern[7] && pattern[8] && !pattern[6]);
    double params2[]={.033,.033}, xy[]={10,.2,.5, 10,.33,.3}, fxy[6]={0}, jxy[18]={0};
    constraints[0]=(ChemConstraint)EQ(1,a,1,b,1);
    chem_constraint_residuals(&set,2,xy,NULL,params2,fxy);
    chem_constraint_subtract_jacobian(&set,2,xy,NULL,params2,jxy);
    CLOSE(fxy[1],.13); CLOSE(fxy[4],0); CLOSE(jxy[3],-.033); CLOSE(jxy[4],1);
    CLOSE(jxy[12],-.033); CLOSE(jxy[13],1); assert(jxy[8]==0 && jxy[17]==0);
    memset(pattern,0,sizeof(pattern)); chem_constraint_sparsity(&set,pattern);
    assert(pattern[3] && pattern[4] && !pattern[5]);

    constraints[0]=(ChemConstraint)EQ(1,a,1,b,1);
    constraints[1]=(ChemConstraint)EQ(0,cd,2,a,1);
    set=(ChemConstraintSet){4,2,constraints};
    assert(chem_constraint_validate(&set)==CHEM_CONSTRAINT_OK);
    double state4[]={.33,.01089,10,1, .33,.01089,5,1, .33,.02,5,1};
    double params4[]={.033,.033,.033,.033,.033,.033};
    double f4[12]={0}, j4[48]={0};
    chem_constraint_residuals(&set,3,state4,NULL,params4,f4);
    chem_constraint_subtract_jacobian(&set,3,state4,NULL,params4,j4);
    CLOSE(f4[0],0); CLOSE(f4[1],0); CLOSE(f4[4],-.165); CLOSE(f4[5],0);
    CLOSE(f4[8],-.165); CLOSE(f4[9],-.00911);
    for (size_t i=0;i<3;++i) {
        CLOSE(j4[16*i+4],-.033); CLOSE(j4[16*i+5],1);
        CLOSE(j4[16*i],1); CLOSE(j4[16*i+2],-.033*state4[4*i+3]);
        CLOSE(j4[16*i+3],-.033*state4[4*i+2]);
    }
    memset(pattern,0,sizeof(pattern)); chem_constraint_sparsity(&set,pattern);
    assert(pattern[0] && pattern[2] && pattern[3] && pattern[4] && pattern[5]);
    count=0; for(size_t i=0;i<16;++i) count+=!!pattern[i]; assert(count==5);
    double errors[12]={0}, next[12];
    for(size_t i=0;i<12;++i) next[i]=state4[i]+.01*(i+1);
    chem_constraint_step_errors(&set,3,state4,next,errors);
    for(size_t i=0;i<3;++i) { CLOSE(errors[4*i],next[4*i]-state4[4*i]); CLOSE(errors[4*i+1],next[4*i+1]-state4[4*i+1]); assert(errors[4*i+2]==0); }

    constraints[1]=(ChemConstraint)EQ(2,a,1,c,1); set=(ChemConstraintSet){3,2,constraints};
    double shared[]={.1,.0033,.0033}, shared_params[]={.033,.033}, shared_f[3]={0}, shared_j[9]={0};
    chem_constraint_residuals(&set,1,shared,NULL,shared_params,shared_f);
    chem_constraint_subtract_jacobian(&set,1,shared,NULL,shared_params,shared_j);
    CLOSE(shared_f[1],0); CLOSE(shared_f[2],0); CLOSE(shared_j[3],-.033); CLOSE(shared_j[4],1);
    CLOSE(shared_j[6],-.033); CLOSE(shared_j[8],1);
    puts("constraint set checks PASS");
    return 0;
}
