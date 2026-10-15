/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_solver.h"
static void can_solve(void)
{
    ChemSpecies sp[]={{.name="A",.absolute_tolerance=1e-6},{.name="B",.absolute_tolerance=1e-6},{.name="C",.absolute_tolerance=1e-6}};
    const size_t lhs[]={0,1};const ChemProduct rhs[]={{1,1},{2,1}};
    ChemReaction reactions[2]={0};
    for (size_t i=0;i<2;++i) {
        reactions[i].reactant_count=reactions[i].product_count=1;
        reactions[i].reactants=lhs+i;reactions[i].products=rhs+i;
        reactions[i].rate=chem_rate_default(CHEM_ARRHENIUS);
    }
    reactions[0].rate.data.arrhenius.A=2.15e-11;reactions[0].rate.data.arrhenius.C=110;
    reactions[1].rate.data.arrhenius.A=3.3e-11;reactions[1].rate.data.arrhenius.C=55;
    ChemMechanism *m=NULL;ChemState *state=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_create(3,sp,2,reactions,&m)==CHEM_OK);
    if (!m) return;
    CHECK(chem_state_create(m,1,&state)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    if (state && solver) {
        state->values[0]=1;state->conditions[0]=(ChemConditions){272.5,101253.3,1e6};
        chem_update_rates(state);
        ChemBackwardEulerParameters p=chem_backward_euler_parameters();ChemSolverStats stats;
        CHECK(chem_solve_backward_euler(solver,state,1,&p,&stats)==CHEM_OK);
        /* Additional analytical one-step implicit-Euler check. */
        REL(state->values[0],1/(1+state->rates[0]),1e-13);
        NEAR(state->values[0]+state->values[1]+state->values[2],1,1e-14);
    }
    chem_solver_free(solver);chem_state_free(state);chem_mechanism_free(m);
}
static void is_converged(void)
{
    ChemSpecies sp[]={{.name="A"},{.name="B"},{.name="C"}};
    ChemMechanism *m=NULL;ChemState *s=NULL;
    CHECK(chem_mechanism_create(3,sp,0,NULL,&m)==CHEM_OK);
    if (!m) return;
    CHECK(chem_state_create(m,4,&s)==CHEM_OK);
    if (s) {
        double residual[12]={0};
        for (size_t v=0;v<3;++v) s->absolute_tolerance[v]=1e-6;
        s->relative_tolerance=1e-3;
        ChemBackwardEulerParameters p=chem_backward_euler_parameters();p.small=1e-6;
        CHECK(chem_backward_euler_converged(s,residual,s->values,&p));
        residual[4]=1e-5;CHECK(!chem_backward_euler_converged(s,residual,s->values,&p));
        p.small=1e-4;CHECK(chem_backward_euler_converged(s,residual,s->values,&p));
        residual[11]=1e-3;CHECK(!chem_backward_euler_converged(s,residual,s->values,&p));
        s->values[11]=10;CHECK(chem_backward_euler_converged(s,residual,s->values,&p));
        residual[11]=1e-1;CHECK(!chem_backward_euler_converged(s,residual,s->values,&p));
        s->absolute_tolerance[2]=1;CHECK(chem_backward_euler_converged(s,residual,s->values,&p));
        residual[11]=INFINITY;CHECK(!chem_backward_euler_converged(s,residual,s->values,&p));
        residual[11]=NAN;CHECK(!chem_backward_euler_converged(s,residual,s->values,&p));
    }
    chem_state_free(s);chem_mechanism_free(m);
}
void test_backward_euler(void)
{
    RUN("backward_euler:BackwardEuler.CanCallSolve",can_solve());
    RUN("backward_euler:BackwardEuler.IsConverged",is_converged());
}
