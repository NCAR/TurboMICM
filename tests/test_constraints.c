/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Ports of MICM test/unit/constraint/type and constraint_set_policy. */
#include "test.h"
#include "chem_constraints.h"
#include "chem_jacobian.h"
#include <stdlib.h>
#include <stdint.h>

typedef struct { ChemMechanism *m;ChemState *s; } Fixture;
static Fixture fixture(size_t n,const char *const *names,size_t cells,size_t count,const ChemConstraint *constraints)
{
    Fixture f={0};ChemSpecies *sp=calloc(n,sizeof(*sp));CHECK(sp!=NULL);
    for (size_t i=0;i<n;++i) sp[i].name=names[i];
    CHECK(chem_mechanism_create(n,sp,0,NULL,&f.m)==CHEM_OK);free(sp);
    CHECK(chem_mechanism_set_constraints(f.m,count,constraints)==CHEM_OK);
    CHECK(chem_state_create(f.m,cells,&f.s)==CHEM_OK);
    for (size_t c=0;c<cells;++c) f.s->conditions[c].temperature=298.15;
    chem_update_rates(f.s);
    return f;
}
static void release(Fixture f) { chem_state_free(f.s);chem_mechanism_free(f.m); }
static ChemConstraint linear(const char *name,const char *row,size_t n,const ChemConstraintTerm *terms,double constant)
{ return (ChemConstraint){.type=CHEM_LINEAR_CONSTRAINT,.name=name,.algebraic_species=row,.reactant_count=n,.reactants=terms,.constant=constant}; }
static ChemConstraint equilibrium(const char *name,const char *row,size_t nr,const ChemConstraintTerm *r,
                                  size_t np,const ChemConstraintTerm *p,double k,double h)
{ return (ChemConstraint){.type=CHEM_EQUILIBRIUM_CONSTRAINT,.name=name,.algebraic_species=row,.reactant_count=nr,.reactants=r,
                         .product_count=np,.products=p,.equilibrium=chem_vant_hoff(k,h)}; }
static const char *const abc[]={"A","B","C"};
static const char *const ab_ab[]={"A","B","AB"};
static void constraint_forcing(const double *y,double *f,void *ctx) { chem_constraints_forcing(ctx,y,f); }
static void verify_fd(Fixture f,int sparsity)
{
    size_t n=chem_species_count(f.m),cells=f.s->cells;
    unsigned char *pattern=calloc(n*n,1);
    double *analytic=calloc(n*n*cells,sizeof(double)),*fd=calloc(n*n*cells,sizeof(double));
    CHECK(pattern && analytic && fd);
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) pattern[i*n+j]=i==j || chem_constraint_has_element(f.m,i,j);
    ChemSparse *sparse=NULL;CHECK(chem_sparse_create(n,pattern,&sparse)==CHEM_OK);
    chem_constraints_subtract_jacobian(f.s,f.s->values,analytic);
    double *packed=calloc(chem_sparse_elements(sparse)*cells,sizeof(double));CHECK(packed!=NULL);
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
        size_t index=chem_sparse_index(sparse,i,j);
        if (index!=SIZE_MAX) for (size_t c=0;c<cells;++c) packed[index*cells+c]=analytic[(i*n+j)*cells+c];
    }
    CHECK(chem_finite_difference(n,cells,constraint_forcing,f.s,f.s->values,1e-8,fd)==CHEM_OK);
    CHECK(chem_compare_jacobian(sparse,cells,packed,fd,1e-7,1e-7).passed);
    if (sparsity) CHECK(chem_check_sparsity(sparse,cells,fd,1e-6).passed);
    free(pattern);free(analytic);free(packed);free(fd);chem_sparse_free(sparse);
}
static void linear_metadata(int which)
{
    ChemConstraintTerm terms[]={{"A",1},{"B",1},{"C",1}};
    const char *names[]={"A_B_conservation","A_B_C_conservation","weighted_sum","A_equals_B","fractional_conservation"};
    size_t n=(which==1 || which==2)?3:2;double constant=which==1?10:which==2?5:which==3?0:which==4?2:1;
    if (which==2) { terms[0].coefficient=2;terms[1].coefficient=3;terms[2].coefficient=-1; }
    if (which==3) terms[1].coefficient=-1;
    if (which==4) { terms[0].coefficient=.5;terms[1].coefficient=1.5; }
    ChemConstraint c=linear(names[which],abc[n-1],n,terms,constant);
    Fixture f=fixture(n,abc,1,1,&c);const ChemConstraint *stored=chem_constraint_definition(f.m,0);
    CHECK(strcmp(stored->name,names[which])==0);CHECK(stored->constant==constant);
    CHECK(stored->reactant_count==n);CHECK(!strcmp(stored->algebraic_species,abc[n-1]));
    for (size_t i=0;i<n;++i) { CHECK(!strcmp(stored->reactants[i].species,abc[i]));CHECK(stored->reactants[i].coefficient==terms[i].coefficient); }
    release(f);
}
static void linear_values(int which)
{
    ChemConstraintTerm terms[]={{"A",1},{"B",1},{"C",1}};
    size_t n=(which==2 || which==3 || which==8)?3:2,cells=which==7?2:1;
    double constant=which==2 || which==8?5:which==3?10:which==4?0:which==5?2:1;
    if (which==2 || which==8) { terms[0].coefficient=2;terms[1].coefficient=3;terms[2].coefficient=-1; }
    if (which==4) terms[1].coefficient=-1;
    if (which==5) { terms[0].coefficient=.5;terms[1].coefficient=1.5; }
    if (which==6) { terms[0].coefficient=2;terms[1].coefficient=3; }
    ChemConstraint c=linear("linear",abc[n-1],n,terms,constant);Fixture f=fixture(n,abc,cells,1,&c);
    double *y=f.s->values,force[6]={0},jac[18]={0};
    if (which==2) { y[0]=1;y[1]=2;y[2]=3; }
    else if (which==3) { y[0]=2;y[1]=3;y[2]=5; }
    else if (which==4) { y[0]=.5;y[1]=.5; }
    else if (which==5) { y[0]=1;y[1]=1; }
    else if (which==6) { y[0]=.1;y[1]=.2; }
    else if (which==7) { y[0]=.3;y[1]=.8;y[2]=.7;y[3]=.2; }
    else if (which==8) { y[0]=1.5;y[1]=.8;y[2]=2; }
    else { y[0]=.3;y[1]=.7; }
    if (which>=7) { verify_fd(f,which==7);release(f);return; }
    chem_constraints_forcing(f.s,y,force);
    if (which!=6) NEAR(force[n-1],0,1e-10);
    chem_constraints_subtract_jacobian(f.s,y,jac);
    for (size_t i=0;i<n;++i) NEAR(jac[(n-1)*n+i],-terms[i].coefficient,1e-10);
    if (which==0) { y[0]=.5;y[1]=.6;chem_constraints_forcing(f.s,y,force);NEAR(force[1],.1,1e-10); }
    if (which==3) { y[2]=6;chem_constraints_forcing(f.s,y,force);NEAR(force[2],1,1e-10); }
    if (which==4) { y[0]=.6;chem_constraints_forcing(f.s,y,force);NEAR(force[1],.1,1e-10); }
    if (which==6) {
        double old[2]={jac[2],jac[3]};y[0]=10;y[1]=20;memset(jac,0,sizeof(jac));
        chem_constraints_subtract_jacobian(f.s,y,jac);
        NEAR(jac[2],old[0],1e-10);NEAR(jac[3],old[1],1e-10);
    }
    release(f);
}
static void equilibrium_metadata(int which)
{
    ChemConstraintTerm r[]={{"A",1},{"B",1}},p[]={{"AB",1},{"C",1}};
    size_t nr=2,np=1,n=3;const char *name="A_B_equilibrium",*row="AB";double k=1000;
    const char *const *names=ab_ab;
    if (which==2) { nr=1;n=2;p[0].species="B";names=abc;row="B";name="A_B_simple";k=10; }
    if (which==3) { nr=1;np=2;r[0].coefficient=2;p[0].species="B";names=abc;row="B";name="dissociation";k=100; }
    ChemConstraint c=equilibrium(name,row,nr,r,np,p,k,-2400);Fixture f=fixture(n,names,1,1,&c);
    const ChemConstraint *stored=chem_constraint_definition(f.m,0);
    CHECK(!strcmp(stored->name,name));CHECK(!strcmp(stored->algebraic_species,row));
    CHECK(stored->reactant_count==nr);CHECK(stored->product_count==np);
    for (size_t i=0;i<nr;++i) { CHECK(!strcmp(stored->reactants[i].species,r[i].species));CHECK(stored->reactants[i].coefficient==r[i].coefficient); }
    for (size_t i=0;i<np;++i) CHECK(!strcmp(stored->products[i].species,p[i].species));
    release(f);
}
static void equilibrium_invalid(int which)
{
    ChemConstraintTerm r={"A",1},p={"B",1};ChemConstraint c=equilibrium("invalid","B",1,&r,1,&p,1,-2400);
    if (which==0) {
        c.equilibrium.reference_constant=-1;CHECK(chem_constraint_validate(&c)==CHEM_BAD_EQUILIBRIUM);
        c.equilibrium.reference_constant=0;CHECK(chem_constraint_validate(&c)==CHEM_BAD_EQUILIBRIUM);
    } else if (which==1) { c.reactant_count=0;CHECK(chem_constraint_validate(&c)==CHEM_EMPTY_REACTANTS); }
    else if (which==2) { c.product_count=0;CHECK(chem_constraint_validate(&c)==CHEM_EMPTY_PRODUCTS); }
    else {
        r.coefficient=0;CHECK(chem_constraint_validate(&c)==CHEM_BAD_STOICHIOMETRY);
        r.coefficient=-1;CHECK(chem_constraint_validate(&c)==CHEM_BAD_STOICHIOMETRY);
        r.coefficient=1;p.coefficient=0;CHECK(chem_constraint_validate(&c)==CHEM_BAD_STOICHIOMETRY);
        p.coefficient=-2;CHECK(chem_constraint_validate(&c)==CHEM_BAD_STOICHIOMETRY);
    }
}
static void equilibrium_values(int which)
{
    ChemConstraintTerm r[]={{"A",1},{"B",1}},p[]={{"AB",1},{"C",1}};
    int complex=which==2 || which==4;size_t cells=which==3?2:1;
    if (complex) { r[0].coefficient=2;p[0].species="B"; }
    ChemConstraint c=equilibrium("eq",complex?"B":"AB",complex?1:2,r,complex?2:1,p,complex?100:1000,which>=3?0:-2400);
    Fixture f=fixture(3,complex?abc:ab_ab,cells,1,&c);double *y=f.s->values;
    double force[6]={0},jac[9]={0};
    if (which==0) {
        y[0]=y[1]=y[2]=.001;chem_constraints_forcing(f.s,y,force);NEAR(force[2],0,1e-10);
        y[0]=.02;y[1]=.03;y[2]=.05;chem_constraints_forcing(f.s,y,force);NEAR(force[2],.55,1e-10);
    } else if (which==1) {
        y[0]=.01;y[1]=.02;y[2]=.05;chem_constraints_subtract_jacobian(f.s,y,jac);
        NEAR(jac[6],-20,1e-8);NEAR(jac[7],-10,1e-8);NEAR(jac[8],1,1e-8);
    } else if (which==2) {
        y[0]=.1;y[1]=.5;y[2]=2;chem_constraints_forcing(f.s,y,force);NEAR(force[1],0,1e-10);
    } else if (which==3) {
        y[0]=.02;y[1]=.1;y[2]=.03;y[3]=.2;y[4]=.05;y[5]=5;verify_fd(f,1);
    } else { y[0]=.15;y[1]=.8;y[2]=1.5;verify_fd(f,1); }
    release(f);
}

static void check_pattern(Fixture f,size_t elements,const unsigned char *expected)
{
    size_t n=chem_species_count(f.m),count=0;
    for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
        int present=chem_constraint_has_element(f.m,i,j);count+=present;
        if (expected) CHECK(present==(expected[i*n+j]!=0));
    }
    CHECK(count==elements);
}
static void set_basic(size_t group)
{
    ChemConstraintTerm r[]={{"A",1},{"B",1}},p={"AB",1};
    ChemConstraint c=equilibrium("A_B_eq","AB",2,r,1,&p,.033,-24000);
    Fixture f=fixture(3,ab_ab,1,1,&c);
    CHECK(chem_constraint_count(f.m)==1);
    unsigned char pattern[]={0,0,0,0,0,0,1,1,1};check_pattern(f,3,pattern);
    f.s->values[0]=.01;f.s->values[1]=.02;f.s->values[2]=.05;
    double jac[9]={0};chem_constraints_subtract_jacobian(f.s,f.s->values,jac);
    NEAR(jac[6],-.00066,1e-10);NEAR(jac[7],-.00033,1e-10);NEAR(jac[8],1,1e-10);release(f);
    double input[2][3]={{.2,.4,.6},{.1,.3,.7}},expected[]={-.59736,-.69901};
    for (size_t start=0;start<2;start+=group) {
        size_t cells=group<2-start?group:2-start;f=fixture(3,ab_ab,cells,1,&c);
        for (size_t cell=0;cell<cells;++cell) for (size_t v=0;v<3;++v) f.s->values[v*cells+cell]=input[start+cell][v];
        double forcing[6]={0};chem_constraints_forcing(f.s,f.s->values,forcing);
        for (size_t cell=0;cell<cells;++cell) NEAR(forcing[2*cells+cell],expected[start+cell],1e-5);
        release(f);
    }
}
static void set_replace(void)
{
    ChemConstraintTerm r={"B",1},p={"C",1};ChemConstraint c=equilibrium("B_C_eq","C",1,&r,1,&p,.033,-24000);
    Fixture f=fixture(3,abc,1,1,&c);CHECK(chem_constraint_count(f.m)==1);
    CHECK(chem_constraint_row(f.m,0)==2);CHECK(chem_is_algebraic(f.m,2));
    CHECK(!chem_is_algebraic(f.m,0));CHECK(!chem_is_algebraic(f.m,1));
    unsigned char pattern[]={0,0,0,0,0,0,0,1,1};check_pattern(f,2,pattern);release(f);
}
static void set_multiple(void)
{
    const char *names[]={"A","B","AB","C","D"};
    ChemConstraintTerm r[]={{"A",1},{"B",1}},p={"AB",1},r2={"C",1},p2={"D",1};
    ChemConstraint c[]={equilibrium("A_B_eq","AB",2,r,1,&p,.033,-24000),equilibrium("C_D_eq","D",1,&r2,1,&p2,.033,-24000)};
    Fixture f=fixture(5,names,1,2,c);CHECK(chem_constraint_count(f.m)==2);
    unsigned char pattern[25]={0};pattern[10]=pattern[11]=pattern[12]=pattern[23]=pattern[24]=1;
    check_pattern(f,5,pattern);release(f);
}
static void set_empty(void)
{
    Fixture f=fixture(3,abc,1,0,NULL);CHECK(chem_constraint_count(f.m)==0);check_pattern(f,0,NULL);
    f.s->values[0]=.1;f.s->values[1]=.2;f.s->values[2]=.3;double forcing[]={1,1,1,1};
    chem_constraints_forcing(f.s,f.s->values,forcing);for (size_t i=0;i<4;++i) CHECK(forcing[i]==1);
    release(f);
}
static void set_unknown(void)
{
    Fixture f=fixture(2,abc,1,0,NULL);
    ChemConstraintTerm r[]={{"X",1},{"Y",1}},p={"XY",1};ChemConstraint c=equilibrium("invalid","XY",2,r,1,&p,.033,-24000);
    CHECK(chem_mechanism_set_constraints(f.m,1,&c)==CHEM_NOT_FOUND);release(f);
}
static void set_three(void)
{
    const char *names[]={"X","Y","Z"};ChemConstraintTerm r={"X",1},p={"Y",1};
    ChemConstraint c=equilibrium("X_Y_eq","Y",1,&r,1,&p,.033,-24000);Fixture f=fixture(3,names,2,1,&c);
    CHECK(chem_constraint_count(f.m)==1);unsigned char pattern[]={0,0,0,1,1,0,0,0,0};check_pattern(f,2,pattern);
    double y[]={10,10,.2,.33,.5,.3};memcpy(f.s->values,y,sizeof(y));
    double forcing[6]={0},jac[18]={0};chem_constraints_forcing(f.s,y,forcing);chem_constraints_subtract_jacobian(f.s,y,jac);
    NEAR(forcing[2],.13,1e-10);NEAR(forcing[3],0,1e-10);
    for (size_t cell=0;cell<2;++cell) { NEAR(jac[6+cell],-.033,1e-10);NEAR(jac[8+cell],1,1e-10);NEAR(jac[16+cell],0,1e-10); }
    release(f);
}
static void set_four(void)
{
    const char *names[]={"A","B","C","D"};ChemConstraintTerm a={"A",1},b={"B",1},cd[]={{"C",1},{"D",1}};
    ChemConstraint constraints[]={equilibrium("A_B_eq","B",1,&a,1,&b,.033,-24000),equilibrium("CD_A_eq","A",2,cd,1,&a,.033,-24000)};
    Fixture f=fixture(4,names,3,2,constraints);CHECK(chem_constraint_count(f.m)==2);
    unsigned char pattern[]={1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0};check_pattern(f,5,pattern);
    double y[]={.33,.33,.33,.01089,.01089,.02,10,5,5,1,1,1};memcpy(f.s->values,y,sizeof(y));
    double forcing[12]={0},jac[48]={0};chem_constraints_forcing(f.s,y,forcing);chem_constraints_subtract_jacobian(f.s,y,jac);
    for (size_t c=0;c<3;++c) {
        NEAR(forcing[3+c],c==2?-.00911:0,1e-5);NEAR(forcing[c],c?-.165:0,1e-10);
        NEAR(jac[12+c],-.033,1e-10);NEAR(jac[15+c],1,1e-10);
        NEAR(jac[6+c],-.033*y[9+c],1e-10);NEAR(jac[9+c],-.033*y[6+c],1e-10);NEAR(jac[c],1,1e-10);
    }
    release(f);
}
static void set_coupled(void)
{
    ChemConstraintTerm a={"A",1},b={"B",1},c={"C",1};
    ChemConstraint constraints[]={equilibrium("A_B_eq","B",1,&a,1,&b,.033,-24000),equilibrium("A_C_eq","C",1,&a,1,&c,.033,-24000)};
    Fixture f=fixture(3,abc,1,2,constraints);CHECK(chem_constraint_count(f.m)==2);
    unsigned char pattern[]={0,0,0,1,1,0,1,0,1};check_pattern(f,4,pattern);
    double y[]={.1,.0033,.0033},forcing[3]={0},jac[9]={0};memcpy(f.s->values,y,sizeof(y));
    chem_constraints_forcing(f.s,y,forcing);chem_constraints_subtract_jacobian(f.s,y,jac);
    NEAR(forcing[1],0,1e-10);NEAR(forcing[2],0,1e-10);
    NEAR(jac[3],-.033,1e-10);NEAR(jac[4],1,1e-10);NEAR(jac[6],-.033,1e-10);NEAR(jac[8],1,1e-10);release(f);
}
static void set_indexing(void)
{
    ChemConstraintTerm r[]={{"A",1},{"B",1}},p={"AB",1};ChemConstraint constraint=equilibrium("A_B_eq","AB",2,r,1,&p,.033,-24000);
    Fixture f=fixture(3,ab_ab,3,1,&constraint);double y[]={.01,.03,.001,.02,.01,.002,.05,.2,.004};
    memcpy(f.s->values,y,sizeof(y));double forcing[9]={0},jac[27]={0};
    chem_constraints_forcing(f.s,y,forcing);chem_constraints_subtract_jacobian(f.s,y,jac);
    for (size_t c=0;c<3;++c) {
        NEAR(forcing[6+c],.033*y[c]*y[3+c]-y[6+c],1e-9);
        NEAR(jac[18+c],-.033*y[3+c],1e-12);NEAR(jac[21+c],-.033*y[c],1e-12);NEAR(jac[24+c],1,1e-12);
    }
    release(f);
}
void test_constraints(void)
{
#define L(name,call) RUN("linear:LinearConstraint." name,call)
    L("Construction",linear_metadata(0));L("AlgebraicSpecies",linear_metadata(1));
    L("WeightedTerms",linear_metadata(2));L("ZeroConstant",linear_metadata(3));L("FractionalCoefficients",linear_metadata(4));
    L("ResidualComputationThroughConstraintSet",linear_values(0));L("JacobianComputationThroughConstraintSet",linear_values(1));
    L("WeightedSumResidualAndJacobian",linear_values(2));L("ThreeSpeciesConservationResidual",linear_values(3));
    L("ZeroConstantResidual",linear_values(4));L("FractionalCoefficientsResidualAndJacobian",linear_values(5));
    L("JacobianIndependentOfConcentrations",linear_values(6));L("FiniteDifferenceJacobianSimpleConservation",linear_values(7));
    L("FiniteDifferenceJacobianWeightedSum",linear_values(8));
#undef L
#define E(name,call) RUN("equilibrium:EquilibriumConstraint." name,call)
    E("Construction",equilibrium_metadata(0));E("AlgebraicSpecies",equilibrium_metadata(1));
    E("SingleReactantSingleProduct",equilibrium_metadata(2));E("MultipleReactantsAndProducts",equilibrium_metadata(3));
    E("InvalidEquilibriumConstant",equilibrium_invalid(0));E("EmptyReactantsThrows",equilibrium_invalid(1));
    E("EmptyProductsThrows",equilibrium_invalid(2));E("InvalidStoichiometryThrows",equilibrium_invalid(3));
    E("ResidualComputationThroughConstraintSet",equilibrium_values(0));E("JacobianComputationThroughConstraintSet",equilibrium_values(1));
    E("ComplexStoichiometryResidual",equilibrium_values(2));E("FiniteDifferenceJacobianSimple",equilibrium_values(3));
    E("FiniteDifferenceJacobianComplexStoichiometry",equilibrium_values(4));
#undef E
#define S(name,call) RUN("constraint_set:ConstraintSet." name,call)
    S("Construction",set_basic(2));S("ReplaceStateRowsMapsToAlgebraicSpecies",set_replace());
    S("NonZeroJacobianElements",set_basic(2));S("MultipleConstraints",set_multiple());
    S("AddForcingTerms",set_basic(2));S("SubtractJacobianTerms",set_basic(2));
    S("EmptyConstraintSet",set_empty());S("UnknownSpeciesThrows",set_unknown());
    S("ThreeDStateOneConstraint",set_three());S("FourDStateTwoConstraints",set_four());
    S("CoupledConstraintsSharedSpecies",set_coupled());S("VectorizedMatricesRespectGridCellIndexing",set_indexing());
    S("VectorMatrix1",set_basic(1));S("VectorMatrix2",set_basic(2));S("VectorMatrix3",set_basic(3));
#undef S
}
