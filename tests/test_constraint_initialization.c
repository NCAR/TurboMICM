/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Port of test/unit/solver/test_constraint_initialization.cpp. */
#include "test.h"
#include "chem_constraints.h"
#include "chem_solver.h"

typedef struct { ChemMechanism *m;ChemState *state;ChemSolver *solver; } System;
static System system_create(size_t n,const char *const *names,size_t cells,int reaction,size_t count,const ChemConstraint *constraints)
{
    System s={0};ChemSpecies species[3]={{0}};
    for (size_t i=0;i<n;++i) species[i].name=names[i];
    size_t reactant=0;ChemProduct product={1,1};ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);rate.data.arrhenius.A=.5;
    ChemReaction rxn={1,1,&reactant,&product,rate};
    CHECK(chem_mechanism_create(n,species,reaction?1:0,&rxn,&s.m)==CHEM_OK);
    CHECK(chem_mechanism_set_constraints(s.m,count,constraints)==CHEM_OK);
    CHECK(chem_state_create(s.m,cells,&s.state)==CHEM_OK);CHECK(chem_solver_create(s.m,cells,&s.solver)==CHEM_OK);
    for (size_t c=0;c<cells;++c) { s.state->conditions[c].temperature=298.15;s.state->conditions[c].pressure=101325; }
    return s;
}
static void destroy(System s) { chem_solver_free(s.solver);chem_state_free(s.state);chem_mechanism_free(s.m); }
static void simple(int which)
{
    const char *names[]={"A","B","C"};ChemConstraintTerm r={"B",1},p={"C",1};
    ChemConstraint constraint={.type=CHEM_EQUILIBRIUM_CONSTRAINT,.name="B_C_eq",.algebraic_species="C",
        .reactant_count=1,.reactants=&r,.product_count=1,.products=&p,.equilibrium=chem_vant_hoff(2,0)};
    size_t cells=which==4?3:1;System s=system_create(which==3?2:3,names,cells,1,which==3?0:1,&constraint);
    ChemRosenbrockParameters options=chem_rosenbrock_parameters(CHEM_ROS3);ChemSolverStats stats;
    s.state->values[0]=1;s.state->values[cells]=which==3?0:.5;
    if (which!=3) s.state->values[2*cells]=which==1?5:which==2?1000:1;
    if (which==4) {
        s.state->values[1]=2;s.state->values[2]=.5;s.state->values[4]=.3;s.state->values[5]=.1;
        s.state->values[7]=5;s.state->values[8]=100;
    }
    if (which!=3) chem_update_rates(s.state);
    CHECK(chem_solve(s.solver,s.state,which==3?.01:.001,&options,&stats)==CHEM_OK);
    if (which==0) CHECK(stats.constraint_init_iterations<=1);
    if (which==1) CHECK(stats.constraint_init_iterations>0);
    if (which==3) CHECK(stats.constraint_init_iterations==0);
    if (which==5) {
        s.state->values[2]=999;chem_update_rates(s.state);
        CHECK(chem_solve(s.solver,s.state,.001,&options,&stats)==CHEM_OK);
    }
    if (which!=3) for (size_t c=0;c<cells;++c)
        NEAR(s.state->custom[c]*s.state->values[cells+c]-s.state->values[2*cells+c],0,which==4?1e-5:1e-6);
    destroy(s);
}
static void scaled_solve(void)
{
    double scales[]={1,1e-12,1e12},reference=0,consistent=1.3-.7-.15;
    const char *names[]={"A","B","C"};
    for (size_t i=0;i<3;++i) {
        double scale=scales[i];ChemConstraintTerm terms[]={{"A",scale},{"B",scale},{"C",scale}};
        ChemConstraint c={.type=CHEM_LINEAR_CONSTRAINT,.name="mass",.algebraic_species="C",
            .reactant_count=3,.reactants=terms,.constant=scale*1.3};
        System s=system_create(3,names,1,1,1,&c);s.state->values[0]=.7;s.state->values[1]=.15;s.state->values[2]=5;
        chem_update_rates(s.state);ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS4_DA);ChemSolverStats stats;
        CHECK(chem_solve(s.solver,s.state,.001,&p,&stats)==CHEM_OK);
        if (!i) { reference=s.state->values[2];NEAR(reference,consistent,1e-5*consistent); }
        else NEAR(s.state->values[2],reference,1e-5*consistent);
        destroy(s);
    }
}
static void squared(int success)
{
    const char *names[]={"B","C"};ChemConstraintTerm r={"B",1},product={"C",2};
    ChemConstraint c={.type=CHEM_EQUILIBRIUM_CONSTRAINT,.name="squared",.algebraic_species="C",
        .reactant_count=1,.reactants=&r,.product_count=1,.products=&product,.equilibrium=chem_vant_hoff(2,0)};
    System s=system_create(2,names,1,0,1,&c);s.state->values[0]=.5;s.state->values[1]=success?1e-6:100;
    chem_update_rates(s.state);ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS3);ChemSolverStats stats={0};
    if (!success) p.constraint_max_iterations=1;
    CHECK(chem_initialize_constraints(s.solver,s.state,&p,&stats)==(success?CHEM_OK:CHEM_CONSTRAINT_FAILURE));
    if (!success) { CHECK(s.state->values[0]==.5);CHECK(s.state->values[1]==100); }
    else {
        NEAR(s.state->values[1],1,1e-3);
        for (size_t backtracks=0;backtracks<=2;backtracks+=2) {
            p.constraint_max_backtracks=backtracks;s.state->values[1]=1e-6;memset(&stats,0,sizeof(stats));
            CHECK(chem_initialize_constraints(s.solver,s.state,&p,&stats)==CHEM_CONSTRAINT_FAILURE);
            CHECK(s.state->values[0]==.5);CHECK(s.state->values[1]==1e-6);
        }
    }
    destroy(s);
}
static void weighted(int state_tolerances)
{
    const char *names[]={"X","Z"};double scales[]={1e-12,1,1e12};
    for (size_t i=0;i<(state_tolerances?2:3);++i) {
        double scale=state_tolerances?1:scales[i];ChemConstraintTerm terms[]={{"X",scale},{"Z",-scale}};
        ChemConstraint c={.type=CHEM_LINEAR_CONSTRAINT,.name="scaled_copy",.algebraic_species="Z",
            .reactant_count=2,.reactants=terms,.constant=0};
        System s=system_create(2,names,1,0,1,&c);
        double atol=state_tolerances && !i?1e-8:1e-12,rtol=state_tolerances && !i?1e-6:1e-8;
        s.state->absolute_tolerance[0]=s.state->absolute_tolerance[1]=atol;s.state->relative_tolerance=rtol;
        s.state->values[0]=1;double initial=s.state->values[1]=state_tolerances?1+5e-8:4;
        ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS3);ChemSolverStats stats={0};
        if (state_tolerances) CHECK(p.constraint_tolerance==.1);
        CHECK(chem_initialize_constraints(s.solver,s.state,&p,&stats)==CHEM_OK);
        if (state_tolerances && !i) { CHECK(s.state->values[1]==initial);CHECK(stats.constraint_init_iterations==1); }
        else { NEAR(s.state->values[1],1,atol);if (state_tolerances) CHECK(stats.constraint_init_iterations==2); }
        destroy(s);
    }
}
void test_constraint_initialization(void)
{
#define I(name,call) RUN("constraint_initialization:ConstraintInitialization." name,call)
    I("ConsistentICsUnchanged",simple(0));I("MildlyInconsistentICsCorrected",simple(1));I("SeverelyInconsistentICsConverge",simple(2));
    I("ScaledConstraintRowsSolveIdenticallyThroughSolve",scaled_solve());I("PureODESystemUnaffected",simple(3));
    I("MultiCellSystems",simple(4));I("SubsequentSolveCallsReinitialize",simple(5));
    I("SolverStateToStringNewValue",CHECK(!strcmp(chem_status_string(CHEM_CONSTRAINT_FAILURE),"Constraint Initialization Failed")));
    I("FailedInitializationRestoresCallerState",squared(0));I("SuccessfulInitializationKeepsTheCorrection",squared(1));
    I("WeightedCorrectionIsInvariantToConstraintRowScaling",weighted(0));I("WeightedCorrectionUsesStateTolerances",weighted(1));
#undef I
}
