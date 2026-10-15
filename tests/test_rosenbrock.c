/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_solver.h"
#include "chem_constraints.h"
#include <stdlib.h>
static ChemMechanism *mechanism(void)
{
    ChemSpecies species[]={{.name="foo"},{.name="bar"},{.name="baz"},{.name="quz"},{.name="quuz"}};
    size_t r1[]={0,2},r2=1,r3=3;ChemProduct p1[]={{1,1},{4,2.4}},p2[]={{0,1},{3,1.4}};
    ChemRate rates[3];for (size_t i=0;i<3;++i) rates[i]=chem_rate_default(CHEM_ARRHENIUS);
    rates[0].data.arrhenius.A=2e-11;rates[0].data.arrhenius.C=110;rates[1].data.arrhenius.A=1e-6;rates[2].data.arrhenius.A=3.5e-6;
    ChemReaction reactions[]={{2,2,r1,p1,rates[0]},{1,2,&r2,p2,rates[1]},{1,0,&r3,NULL,rates[2]}};
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(5,species,3,reactions,&m)==CHEM_OK);return m;
}
static void alpha(size_t cells)
{
    ChemMechanism *m=mechanism();
    size_t rows[]={0,0,0,1,1,1,2,2,3,3,4,4,4},cols[]={0,1,2,0,1,2,0,2,1,3,0,2,4};
    unsigned char pattern[25]={0};for (size_t i=0;i<13;++i) pattern[rows[i]*5+cols[i]]=1;
    ChemSparse *sparse=NULL;CHECK(chem_sparse_create(5,pattern,&sparse)==CHEM_OK);CHECK(chem_sparse_order(sparse)==5);
    size_t count=chem_sparse_elements(sparse)*cells;CHECK(count>=13*cells);
    double *matrix=calloc(count,sizeof(double));CHECK(matrix!=NULL);
    for (size_t c=0;c<cells;++c) {
        double values[]={-12.2,-24.3*(c+2),-42.3,-.43,-23.4,-83.4/(c+3),-4.74,-6.91,-59.1,-83.4,-78.5,-53.6,-1};
        for (size_t i=0;i<13;++i) matrix[chem_sparse_index(sparse,rows[i],cols[i])*cells+c]=values[i];
    }
    chem_add_mass_diagonal(m,sparse,cells,42.042,matrix);
    for (size_t c=0;c<cells;++c) {
        double values[]={42.042-12.2,-24.3*(c+2),-42.3,-.43,42.042-23.4,-83.4/(c+3),-4.74,42.042-6.91,-59.1,42.042-83.4,-78.5,-53.6,42.042-1};
        for (size_t i=0;i<13;++i) CHECK(matrix[chem_sparse_index(sparse,rows[i],cols[i])*cells+c]==values[i]);
    }
    free(matrix);chem_sparse_free(sparse);chem_mechanism_free(m);
}
static void alpha_suite(int vector)
{ const size_t counts[]={1,4,3,2};for (size_t i=0;i<4;++i) alpha(vector?counts[i]:i+1); }
static void tolerances(void)
{
    ChemSpecies sp[]={{.name="foo",.absolute_tolerance=1e-7},{.name="bar",.absolute_tolerance=1e-8}};
    size_t reactant=0;ChemProduct product={1,1};ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);
    rate.data.arrhenius.A=2e-11;rate.data.arrhenius.C=110;ChemReaction r={1,1,&reactant,&product,rate};
    for (size_t cells=1;cells<=10;++cells) {
        ChemMechanism *m=NULL;ChemState *s=NULL;CHECK(chem_mechanism_create(2,sp,1,&r,&m)==CHEM_OK);
        CHECK(chem_mechanism_reorder(m)==CHEM_OK);CHECK(chem_state_create(m,cells,&s)==CHEM_OK);
        CHECK(chem_species_count(m)==2);CHECK(s->absolute_tolerance[0]==1e-7);CHECK(s->absolute_tolerance[1]==1e-8);
        chem_state_free(s);chem_mechanism_free(m);
    }
}
static void normalized(size_t cells,int constraints,ChemRosenbrockMethod method)
{
    ChemMechanism *m=NULL;ChemState *s=NULL;ChemSolver *solver=NULL;
    if (!constraints) m=mechanism();
    else {
        ChemSpecies sp[]={{.name="A"},{.name="B"},{.name="C"}};size_t reactant=0;ChemProduct product={1,1};
        ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);rate.data.arrhenius.A=.5;ChemReaction r={1,1,&reactant,&product,rate};
        CHECK(chem_mechanism_create(3,sp,1,&r,&m)==CHEM_OK);
        ChemConstraintTerm b={"B",1},c={"C",1};ChemConstraint constraint={.type=CHEM_EQUILIBRIUM_CONSTRAINT,
            .name="B_C_eq",.algebraic_species="C",.reactant_count=1,.reactants=&b,.product_count=1,.products=&c,.equilibrium=chem_vant_hoff(10,-2400)};
        CHECK(chem_mechanism_set_constraints(m,1,&constraint)==CHEM_OK);CHECK(chem_constraint_count(m)==1);
    }
    CHECK(chem_state_create(m,cells,&s)==CHEM_OK);CHECK(chem_solver_create(m,cells,&solver)==CHEM_OK);
    /* Parameter sets share the same error norm; construct every original set. */
    ChemRosenbrockParameters parameters=chem_rosenbrock_parameters(method);CHECK(parameters.stages>=2);
    size_t n=chem_species_count(m),count=n*cells;double *data=calloc(count*3,sizeof(double));CHECK(data!=NULL);
    double *old=data,*next=data+count,*error=data+count*2,sum=0;
    for (size_t i=0;i<cells;++i) for (size_t j=0;j<n;++j) {
        size_t index=j*cells+i;
        old[index]=constraints?1+i+.1*j:7.7*i+j;next[index]=constraints?.8+.5*i+.2*j:-13.9/(j+1)-i;
        error[index]=constraints?.01*(1+i+j):81.57/(i+7)/(j+3);
        double scale=s->absolute_tolerance[j]+s->relative_tolerance*fmax(fabs(old[index]),fabs(next[index]));
        sum+=error[index]*error[index]/(scale*scale);
    }
    double expected=fmax(sqrt(sum/count),1e-10),actual=chem_normalized_error(s,old,next,error);
    if (constraints) NEAR(actual,expected,1e-12);
    else CHECK(fabs(actual-expected)/fmax(fabs(actual),fabs(expected))<=1e-11);
    free(data);chem_solver_free(solver);chem_state_free(s);chem_mechanism_free(m);
}
static void normalized_suite(int vector,int constraints)
{
    if (constraints) {
        for (int method=CHEM_ROS2;method<=CHEM_ROS6_DA;++method) normalized(vector?3:2,1,(ChemRosenbrockMethod)method);
        normalized(5,1,CHEM_ROS3);
    } else {
        for (size_t i=1;i<=4;++i) normalized(i,0,CHEM_ROS3);
        if (vector) { const size_t cells[]={1,2,3,5,3};for (size_t i=0;i<5;++i) normalized(cells[i],0,CHEM_ROS3); }
    }
}
void test_rosenbrock(void)
{
#define R(name,call) RUN("rosenbrock:RosenbrockSolver." name,call)
    R("StandardAlphaMinusJacobian",alpha_suite(0));R("VectorAlphaMinusJacobian",alpha_suite(1));R("CanSetTolerances",tolerances());
    R("StandardNormalizedError",normalized_suite(0,0));R("VectorNormalizedError",normalized_suite(1,0));
    R("StandardNormalizedErrorWithConstraints",normalized_suite(0,1));R("VectorNormalizedErrorWithConstraints",normalized_suite(1,1));
#undef R
}
