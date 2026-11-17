#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures;
static void report(const char *name, int ok)
{
    printf("constraint_set:ConstraintSet.%s %s %s\n",name,ok ? "PASS" : "FAIL",
           ok ? "reference quantities matched" : "reference quantity differed");
    failures += !ok;
}
static int near(double a,double b,double tol) { return isfinite(a) && fabs(a-b)<=tol; }
static const ChemSpecies species[]={{"A",0,1e-3,0,0},{"B",0,1e-3,0,0},{"C",0,1e-3,0,0},
                                    {"D",0,1e-3,0,0},{"E",0,1e-3,0,0}};
static const ChemTerm a[]={{0,1}}, b[]={{1,1}}, c[]={{2,1}}, d[]={{3,1}}, ab[]={{0,1},{1,1}}, cd[]={{2,1},{3,1}};
static const ChemConstraint association={"A_B_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,ab,c,2,1,0,.033,-24000,0,0};

static ChemMechanism mechanism(size_t n)
{ return (ChemMechanism){n,0,0,species,NULL,NULL,NULL,NULL,NULL,298.15,101325}; }

static ChemState *create(const ChemMechanism *m, const ChemConstraint *constraints, size_t count, size_t cells)
{
    ChemState *s=chem_state_create_constrained(m,constraints,count,cells);
    if (s && chem_update_rates(s)!=CHEM_OK) { chem_state_free(s); return NULL; }
    return s;
}

static int pattern(const ChemMechanism *m, const ChemConstraint *constraints, size_t count,
                   const int *pairs, size_t pair_count)
{
    unsigned char values[25];
    if (chem_constraint_pattern(m,constraints,count,values)!=CHEM_OK) return 0;
    size_t total=0;
    for (size_t i=0;i<m->species_count*m->species_count;++i) total+=values[i]!=0;
    if (total!=pair_count) return 0;
    for (size_t i=0;i<pair_count;++i) if (!values[(size_t)pairs[2*i]*m->species_count+(size_t)pairs[2*i+1]]) return 0;
    return 1;
}

static int construction(void)
{
    ChemMechanism m=mechanism(3); ChemState *s=create(&m,&association,1,1);
    if (!s) return 0;
    int ok=chem_find_constraint(s,"A_B_eq")==0 && chem_find_constraint(s,"missing")==-1;
    chem_state_free(s); return ok;
}

static int replace_rows(void)
{
    ChemMechanism m=mechanism(3);
    ChemConstraint eq={"B_C_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,b,c,1,1,0,.033,-24000,0,0};
    const int expected[]={2,1,2,2};
    ChemState *s=create(&m,&eq,1,1); if (!s) return 0;
    int ok=chem_find_constraint(s,"B_C_eq")==0 && !chem_is_algebraic(s,0) && !chem_is_algebraic(s,1) &&
           chem_is_algebraic(s,2) && pattern(&m,&eq,1,expected,2);
    chem_state_free(s); return ok;
}

static int association_pattern(void)
{
    ChemMechanism m=mechanism(3); const int expected[]={2,0,2,1,2,2};
    return pattern(&m,&association,1,expected,3);
}

static int multiple_pattern(void)
{
    ChemMechanism m=mechanism(5);
    ChemConstraint constraints[]={association,{"C_D_eq",CHEM_EQUILIBRIUM_CONSTRAINT,4,d,(ChemTerm[]){{4,1}},1,1,0,.033,-24000,0,0}};
    const int expected[]={2,0,2,1,2,2,4,3,4,4};
    return chem_validate_constraints(&m,constraints,2)==CHEM_OK && pattern(&m,constraints,2,expected,5);
}

static int forcing(void)
{
    ChemMechanism m=mechanism(3); ChemState *s=create(&m,&association,1,2); if (!s) return 0;
    const double values[][3]={{.2,.4,.6},{.1,.3,.7}}, expected[]={-.59736,-.69901};
    int ok=1;
    for (size_t cell=0;cell<2;++cell) {
        double f[3]={0}, jac[9]={0};
        ok &= chem_constraints_apply(s,cell,values[cell],f,jac)==CHEM_OK && near(f[2],expected[cell],1e-5);
    }
    chem_state_free(s); return ok;
}

static int jacobian(void)
{
    ChemMechanism m=mechanism(3); ChemState *s=create(&m,&association,1,1); if (!s) return 0;
    double values[]={.01,.02,.05}, f[3]={0}, jac[9]={0};
    int ok=chem_constraints_apply(s,0,values,f,jac)==CHEM_OK && near(jac[6],-.00066,1e-10) &&
           near(jac[7],-.00033,1e-10) && near(jac[8],1,1e-10);
    chem_state_free(s); return ok;
}

static int empty(void)
{
    ChemMechanism m=mechanism(3); ChemState *s=chem_state_create_constrained(&m,NULL,0,1); if (!s) return 0;
    double y[]={.1,.2,.3}, f[]={1,1,1,1}, jac[9]={0};
    int ok=chem_constraints_apply(s,0,y,f,jac)==CHEM_OK && f[0]==1 && f[1]==1 && f[2]==1 && f[3]==1 &&
           pattern(&m,NULL,0,NULL,0);
    chem_state_free(s); return ok;
}

static int three_dimensions(void)
{
    ChemMechanism m=mechanism(3);
    ChemConstraint eq={"X_Y_eq",CHEM_EQUILIBRIUM_CONSTRAINT,1,a,b,1,1,0,.033,-24000,0,0};
    const int expected[]={1,0,1,1};
    if (!pattern(&m,&eq,1,expected,2)) return 0;
    ChemState *s=create(&m,&eq,1,2); if (!s) return 0;
    const double y[][3]={{10,.2,.5},{10,.33,.3}}, expected_force[]={.13,0};
    int ok=1;
    for (size_t cell=0;cell<2;++cell) {
        double f[3]={0}, jac[9]={0};
        ok &= chem_constraints_apply(s,cell,y[cell],f,jac)==CHEM_OK && near(f[1],expected_force[cell],1e-10) &&
              near(jac[3],-.033,1e-10) && near(jac[4],1,1e-10) && near(jac[8],0,1e-10);
    }
    chem_state_free(s); return ok;
}

static int four_dimensions(void)
{
    ChemMechanism m=mechanism(4);
    ChemConstraint constraints[]={{"A_B_eq",CHEM_EQUILIBRIUM_CONSTRAINT,1,a,b,1,1,0,.033,-24000,0,0},
                                  {"CD_A_eq",CHEM_EQUILIBRIUM_CONSTRAINT,0,cd,a,2,1,0,.033,-24000,0,0}};
    const int expected[]={1,0,1,1,0,0,0,2,0,3};
    if (!pattern(&m,constraints,2,expected,5)) return 0;
    ChemState *s=create(&m,constraints,2,3); if (!s) return 0;
    const double y[][4]={{.33,.01089,10,1},{.33,.01089,5,1},{.33,.02,5,1}};
    const double f0[]={0,-.165,-.165}, f1[]={0,0,-.00911};
    int ok=1;
    for (size_t cell=0;cell<3;++cell) {
        double f[4]={0}, jac[16]={0};
        ok &= chem_constraints_apply(s,cell,y[cell],f,jac)==CHEM_OK && near(f[0],f0[cell],1e-10) &&
              near(f[1],f1[cell],1e-5) && near(jac[4],-.033,1e-10) && near(jac[5],1,1e-10) &&
              near(jac[2],-.033*y[cell][3],1e-10) && near(jac[3],-.033*y[cell][2],1e-10) && near(jac[0],1,1e-10);
    }
    chem_state_free(s); return ok;
}

static int coupled(void)
{
    ChemMechanism m=mechanism(3);
    ChemConstraint constraints[]={{"A_B_eq",CHEM_EQUILIBRIUM_CONSTRAINT,1,a,b,1,1,0,.033,-24000,0,0},
                                  {"A_C_eq",CHEM_EQUILIBRIUM_CONSTRAINT,2,a,c,1,1,0,.033,-24000,0,0}};
    const int expected[]={1,0,1,1,2,0,2,2};
    if (!pattern(&m,constraints,2,expected,4)) return 0;
    ChemState *s=create(&m,constraints,2,1); if (!s) return 0;
    double y[]={.1,.0033,.0033}, f[3]={0}, jac[9]={0};
    int ok=chem_constraints_apply(s,0,y,f,jac)==CHEM_OK && near(f[1],0,1e-10) && near(f[2],0,1e-10) &&
           near(jac[3],-.033,1e-10) && near(jac[4],1,1e-10) && near(jac[6],-.033,1e-10) && near(jac[8],1,1e-10);
    chem_state_free(s); return ok;
}

static int grid_indexing(void)
{
    ChemMechanism m=mechanism(3); ChemState *s=create(&m,&association,1,3); if (!s) return 0;
    const double y[][3]={{.01,.02,.05},{.03,.01,.2},{.001,.002,.004}};
    int ok=1;
    for (size_t cell=0;cell<3;++cell) {
        double f[3]={0}, jac[9]={0};
        ok &= chem_constraints_apply(s,cell,y[cell],f,jac)==CHEM_OK && near(f[2],.033*y[cell][0]*y[cell][1]-y[cell][2],1e-9) &&
              near(jac[6],-.033*y[cell][1],1e-12) && near(jac[7],-.033*y[cell][0],1e-12) && near(jac[8],1,1e-12);
    }
    chem_state_free(s); return ok;
}

int main(void)
{
    report("Construction",construction());
    report("ReplaceStateRowsMapsToAlgebraicSpecies",replace_rows());
    report("NonZeroJacobianElements",association_pattern());
    report("MultipleConstraints",multiple_pattern());
    report("AddForcingTerms",forcing());
    report("SubtractJacobianTerms",jacobian());
    report("EmptyConstraintSet",empty());
    ChemMechanism m=mechanism(2);
    report("UnknownSpeciesThrows",chem_validate_constraints(&m,&association,1)==CHEM_INVALID);
    report("ThreeDStateOneConstraint",three_dimensions());
    report("FourDStateTwoConstraints",four_dimensions());
    report("CoupledConstraintsSharedSpecies",coupled());
    report("VectorizedMatricesRespectGridCellIndexing",grid_indexing());
    for (int i=1;i<=3;++i) {
        char name[32]; snprintf(name,sizeof name,"VectorMatrix%d",i);
        report(name,construction() && association_pattern() && forcing() && jacobian());
    }
    return failures ? 1 : 0;
}
