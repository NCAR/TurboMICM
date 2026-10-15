/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Equilibrium integration and DAE constraint overshoot regression ports. */
#include "test.h"
#include "chem_constraints.h"
#include "chem_solver.h"
#include <float.h>

static ChemConstraint eq(const char *name,const char *row,const ChemConstraintTerm *r,const ChemConstraintTerm *p,double k,double h)
{
    return (ChemConstraint){.type=CHEM_EQUILIBRIUM_CONSTRAINT,.name=name,.algebraic_species=row,
        .reactant_count=1,.reactants=r,.product_count=1,.products=p,.equilibrium=chem_vant_hoff(k,h)};
}
static ChemReaction reaction(const size_t *r,const ChemProduct *p,double rate)
{
    ChemRate k=chem_rate_default(CHEM_ARRHENIUS);k.data.arrhenius.A=rate;
    return (ChemReaction){1,1,r,p,k};
}
static void equilibrium_integration(int which)
{
    const char *names[]={"A","B","C","D","E","F"};size_t n=which==1?6:which==4?4:3;
    ChemSpecies species[6]={{0}};for (size_t i=0;i<n;++i) species[i].name=names[i];
    size_t r[]={0,3};ChemProduct products[]={{1,1},{4,1}};
    ChemReaction reactions[]={reaction(r,products,which==0?.1:which==1?.5:1),reaction(r+1,products+1,.2)};
    ChemConstraintTerm b={"B",1},c={"C",1},e={"E",1},f={"F",1},d={"D",1};
    double k1=which<2?.034:which==4?3:2,k2=which==1?.012:5;
    ChemConstraint constraints[]={eq("B_C_eq","C",&b,&c,k1,-2400),
        which==1?eq("E_F_eq","F",&e,&f,k2,-2000):eq("B_D_eq","D",&b,&d,k2,-2400)};
    size_t nc=which==1 || which==4?2:1;
    ChemMechanism *m=NULL;ChemState *state=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_create(n,species,which==1?2:1,reactions,&m)==CHEM_OK);
    CHECK(chem_mechanism_set_constraints(m,nc,constraints)==CHEM_OK);
    if (which==3) CHECK(chem_mechanism_reorder(m)==CHEM_OK);
    CHECK(chem_state_create(m,1,&state)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    CHECK(chem_species_count(m)==n);CHECK(chem_constraint_count(m)==nc);
    size_t index[6]={0},parameter[2]={0};
    for (size_t i=0;i<n;++i) CHECK(chem_species_index(m,names[i],index+i)==CHEM_OK);
    CHECK(chem_custom_index(m,"B_C_eq",parameter)==CHEM_OK);
    if (nc==2) CHECK(chem_custom_index(m,which==1?"E_F_eq":"B_D_eq",parameter+1)==CHEM_OK);
    if (which<2) for (size_t i=0;i<n;++i) CHECK(chem_is_algebraic(m,index[i])==(i==2 || i==5));
    state->values[index[0]]=1;
    if (which<2) { state->values[index[1]]=.1;state->values[index[2]]=k1*.1; }
    if (which==1) { state->values[index[3]]=.5;state->values[index[4]]=.05;state->values[index[5]]=k2*.05; }
    double temperature=which==1?310:which==2?270:which==3?400:300;
    state->conditions[0]=(ChemConditions){temperature,101325,0};chem_update_rates(state);
    double expected=k1*exp((-2400/CHEM_GAS_CONSTANT)*(1/temperature-1/298.15));
    NEAR(state->custom[parameter[0]],expected,which==1 || which==3 || which==4?1e-10:fabs(expected)*100*DBL_EPSILON);
    if (nc==2) {
        expected=k2*exp(((which==1?-2000:-2400)/CHEM_GAS_CONSTANT)*(1/temperature-1/298.15));
        NEAR(state->custom[parameter[1]],expected,which==1?fabs(expected)*100*DBL_EPSILON:1e-10);
    }
    if (which>=2) {
        ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS3);ChemSolverStats stats;
        for (double time=0;time<.1;time+=.001) {
            chem_update_rates(state);ChemStatus status=chem_solve(solver,state,.001,&p,&stats);CHECK(status==CHEM_OK);
            if (status) break;
            NEAR(state->custom[parameter[0]]*state->values[index[1]]-state->values[index[2]],0,1e-6);
            if (which==4) NEAR(state->custom[parameter[1]]*state->values[index[1]]-state->values[index[3]],0,1e-6);
        }
        NEAR(state->values[index[2]],state->custom[parameter[0]]*state->values[index[1]],1e-6);
        if (which==4) NEAR(state->values[index[3]],state->custom[parameter[1]]*state->values[index[1]],1e-6);
        NEAR(state->values[index[0]]+state->values[index[1]],1,.01);
    }
    chem_solver_free(solver);chem_state_free(state);chem_mechanism_free(m);
}
static void nonunit_stoichiometry(void)
{
    ChemSpecies species[]={{.name="A"},{.name="B"},{.name="C"}};
    size_t reactant=2;ChemProduct product={0,1};ChemReaction rxn=reaction(&reactant,&product,.5);
    ChemConstraintTerm a={"A",2},b={"B",1};ChemConstraint constraint=eq("A2_B_eq","B",&a,&b,10,-2400);
    ChemMechanism *m=NULL;ChemState *state=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_create(3,species,1,&rxn,&m)==CHEM_OK);
    CHECK(chem_mechanism_set_constraints(m,1,&constraint)==CHEM_OK);
    CHECK(chem_state_create(m,1,&state)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    state->values[0]=.1;state->values[1]=10*.1*.1;state->values[2]=1;
    state->conditions[0]=(ChemConditions){298,101325,0};chem_update_rates(state);
    double expected=10*exp((-2400/CHEM_GAS_CONSTANT)*(1/298.0-1/298.15));NEAR(state->custom[0],expected,1e-10);
    ChemRosenbrockParameters options=chem_rosenbrock_parameters(CHEM_ROS3);ChemSolverStats stats;
    for (double time=0;time<.05;time+=.001) {
        chem_update_rates(state);ChemStatus status=chem_solve(solver,state,.001,&options,&stats);CHECK(status==CHEM_OK);
        if (status) break;
        NEAR(state->custom[0]*state->values[0]*state->values[0]-state->values[1],0,1e-5);
    }
    chem_solver_free(solver);chem_state_free(state);chem_mechanism_free(m);
}
static void overshoot(ChemRosenbrockMethod method,int equilibrium)
{
    ChemSpecies species[]={{.name=equilibrium?"A_gas":"A"},{.name=equilibrium?"A_aq":"B"},{.name=equilibrium?"P":"C"}};
    size_t r=equilibrium?1:0;ChemProduct product={equilibrium?2:1,1};ChemReaction rxn=reaction(&r,&product,equilibrium?1e3:1e4);
    ChemConstraintTerm gas={"A_gas",1},aq={"A_aq",1};
    ChemConstraintTerm terms[3]={{equilibrium?"A_aq":"A",1},{equilibrium?"P":"B",1},{equilibrium?"A_gas":"C",1}};
    ChemConstraint conservation={.type=CHEM_LINEAR_CONSTRAINT,.name="mass_conservation",
        .algebraic_species=equilibrium?"A_gas":"C",.reactant_count=3,.reactants=terms,.constant=1e-6};
    ChemConstraint constraints[]={equilibrium?eq("gas_aq_eq","A_aq",&gas,&aq,10,0):conservation,conservation};
    ChemMechanism *m=NULL;ChemState *state=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_create(3,species,1,&rxn,&m)==CHEM_OK);
    CHECK(chem_mechanism_set_constraints(m,equilibrium?2:1,constraints)==CHEM_OK);
    CHECK(chem_state_create(m,1,&state)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    state->relative_tolerance=1e-6;
    for (size_t i=0;i<3;++i) state->absolute_tolerance[i]=equilibrium && i<2?1e-8:1e-12;
    state->values[0]=equilibrium?1e-6/11:.9e-6;state->values[1]=equilibrium?10*state->values[0]:0;
    state->values[2]=equilibrium?0:.1e-6;state->conditions[0]=(ChemConditions){298,101325,0};chem_update_rates(state);
    ChemRosenbrockParameters p=chem_rosenbrock_parameters(method);ChemSolverStats stats;
    double advanced=0;
    while (advanced<30) {
        ChemStatus status=chem_solve(solver,state,30-advanced,&p,&stats);CHECK(status==CHEM_OK);
        if (status || stats.final_time<=0) break;
        advanced+=stats.final_time;
        NEAR(state->values[0]+state->values[1]+state->values[2],1e-6,1e-12);
        CHECK(state->values[equilibrium?0:2]>=-1e-18);
        if (equilibrium) CHECK(state->values[1]>=-1e-18);
    }
    if (equilibrium) { NEAR(state->values[2],1e-6,1e-8);CHECK(state->values[0]>=-1e-18); }
    else { CHECK(state->values[0]<1e-12);NEAR(state->values[1],.9e-6,1e-10);NEAR(state->values[2],.1e-6,1e-10); }
    chem_solver_free(solver);chem_state_free(state);chem_mechanism_free(m);
}
static void all_orders(void)
{ for (int method=CHEM_ROS3;method<=CHEM_ROS6_DA;++method) overshoot((ChemRosenbrockMethod)method,0); }

static void terminator_robertson(void)
{
    ChemSpecies species[]={{.name="Cl2",.absolute_tolerance=1e-20},{.name="Cl",.absolute_tolerance=1e-20},
        {.name="A"},{.name="B"},{.name="C"}};
    size_t r0=0,r1[]={1,1},r2=2,r3[]={3,3},r4[]={3,4};
    ChemProduct p0={1,2},p1={0,1},p2={3,1},p3[]={{3,1},{4,1}},p4[]={{2,1},{4,1}};
    ChemRate rates[5];const char *labels[]={"terminator_k1",NULL,"robertson_r1","robertson_r2","robertson_r3"};
    for (size_t i=0;i<5;++i) { rates[i]=chem_rate_default(i==1?CHEM_ARRHENIUS:CHEM_USER);rates[i].label=labels[i]; }
    rates[1].data.arrhenius.A=1;
    ChemReaction reactions[]={{1,1,&r0,&p0,rates[0]},{2,1,r1,&p1,rates[1]},{1,1,&r2,&p2,rates[2]},
        {2,2,r3,p3,rates[3]},{2,2,r4,p4,rates[4]}};
    ChemConstraintTerm terms[]={{"A",1},{"B",1},{"C",1}};
    ChemConstraint constraint={.type=CHEM_LINEAR_CONSTRAINT,.name="mass_conservation",.algebraic_species="C",
        .reactant_count=3,.reactants=terms,.constant=1};
    ChemMechanism *m=NULL;ChemState *state=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_create(5,species,5,reactions,&m)==CHEM_OK);
    CHECK(chem_mechanism_set_constraints(m,1,&constraint)==CHEM_OK);
    CHECK(chem_state_create(m,1,&state)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    state->relative_tolerance=1e-8;double k[]={.04,3e7,1e4};
    for (size_t i=0;i<3;++i) CHECK(chem_state_set_custom(state,labels[i+2],k+i,1)==CHEM_OK);
    double initial[]={1.8e-10,1.2e-6,1,0,0};memcpy(state->values,initial,sizeof(initial));
    state->conditions[0]=(ChemConditions){298,101300,42};chem_update_rates(state);
    const double chlorine=2*1.8e-10+1.2e-6;
    ChemRosenbrockParameters options=chem_rosenbrock_parameters(CHEM_ROS4_DA);ChemSolverStats stats;
    double dt=1;
    for (size_t i=0;i<12;++i,dt*=10) {
        double advanced=0;
        while (advanced<dt) {
            ChemStatus status=chem_solve(solver,state,dt-advanced,&options,&stats);
            CHECK(status==CHEM_OK || status==CHEM_MAX_STEPS);CHECK(advanced+stats.final_time>advanced);
            if ((status!=CHEM_OK && status!=CHEM_MAX_STEPS) || advanced+stats.final_time<=advanced) goto finish;
            advanced+=stats.final_time;
        }
        double *y=state->values;
        NEAR(y[2]+y[3]+y[4],1,1e-10);REL(2*y[0]+y[1],chlorine,1e-6);
        for (size_t j=0;j<5;++j) { CHECK(isfinite(y[j]));CHECK(y[j]>=0); }
        if (!i) { CHECK(y[2]<.99);CHECK(y[3]>0);CHECK(y[4]>1e-3); }
    }
    CHECK(state->values[4]>.95);CHECK(state->values[2]<.05);
    CHECK(state->values[0]>.4*chlorine);CHECK(state->values[1]<.1*chlorine);
finish:
    chem_solver_free(solver);chem_state_free(state);chem_mechanism_free(m);
}
void test_dae(void)
{
#define E(name,call) RUN("equilibrium_constraint:EquilibriumIntegration." name,call)
    E("SetConstraintsAPIWorks",equilibrium_integration(0));E("SetConstraintsAPIMultipleConstraints",equilibrium_integration(1));
    E("DAESolveWithConstraint",equilibrium_integration(2));E("DAESolveWithConstraintAndReorderState",equilibrium_integration(3));
    E("DAESolveWithTwoCoupledConstraints",equilibrium_integration(4));
    E("DAESolveWithNonUnitStoichiometry",nonunit_stoichiometry());
#undef E
#define O(name,call) RUN("dae_constraint_overshoot:DAEConstraintOvershoot." name,call)
    O("AlgebraicVariableStaysNonNegative",overshoot(CHEM_ROS4_DA,0));
    O("EquilibriumPlusConservation",overshoot(CHEM_ROS4_DA,1));O("AllRosenbrockOrdersConstrained",all_orders());
#undef O
    RUN("linear_constraint:DAESolveWithConstraint.TerminatorAndRobertson",terminator_robertson());
}
