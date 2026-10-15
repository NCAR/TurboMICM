/* Copyright (C) 2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Ports of aerosol_model_policy.hpp and the two external stub models. */
#include "test.h"
#include "chem_external.h"
#include <stdint.h>
#include <stdlib.h>

static const char *const stub1[]={
    "STUB1.MODE1.QUUX.BAZ","STUB1.MODE1.QUUX.QUX",
    "STUB1.MODE2.QUUX.BAZ","STUB1.MODE2.QUUX.QUX",
    "STUB1.MODE2.CORGE.FO2","STUB1.MODE2.CORGE.BAZ","STUB1.MODE2.CORGE.QUX"};
static const char *const stub2[]={
    "STUB2.MODE1.NUMBER","STUB2.MODE1.QUUX.BAZ","STUB2.MODE1.QUUX.QUX",
    "STUB2.MODE2.NUMBER","STUB2.MODE2.CORGE.FO2","STUB2.MODE2.CORGE.BAZ","STUB2.MODE2.CORGE.QUX",
    "STUB2.MODE3.NUMBER","STUB2.MODE3.QUUX.BAZ","STUB2.MODE3.QUUX.QUX",
    "STUB2.MODE3.CORGE.FO2","STUB2.MODE3.CORGE.BAZ","STUB2.MODE3.CORGE.QUX"};
static const char *const parameter_names[]={
    "STUB2.PARAM.MODE2.CORGE.FO2_TO_BAZ_RATE_CONSTANT","STUB2.PARAM.MODE3.QUUX.BAZ_TO_QUX_RATE_CONSTANT"};
static void update(const ChemModelBinding *b,const void *context,ChemState *s)
{
    (void)context;
    for (size_t c=0;c<s->cells;++c) s->custom[b->parameters[1]*s->cells+c]=.005*s->conditions[c].temperature;
}
static void forcing(const ChemModelBinding *b,const void *context,const ChemState *s,const double *y,double *f)
{
    (void)context;size_t n=s->cells;
    for (size_t r=0;r<2;++r) for (size_t c=0;c<n;++c) {
        double k=b->parameter_count?s->custom[b->parameters[r]*n+c]:r?2e-3:1e-3;
        double rate=k*y[b->variables[2*r]*n+c];
        f[b->variables[2*r]*n+c]-=rate;f[b->variables[2*r+1]*n+c]+=rate;
    }
}
static void derivative(const ChemModelBinding *b,const void *context,const ChemState *s,const double *y,ChemModelJacobian *j)
{
    (void)context;(void)y;
    for (size_t r=0;r<2;++r) for (size_t c=0;c<s->cells;++c) {
        double k=b->parameter_count?s->custom[b->parameters[r]*s->cells+c]:r?2e-3:1e-3;
        CHECK(chem_model_jacobian_add(j,b->variables[2*r],b->variables[2*r],c,-k)==CHEM_OK);
        CHECK(chem_model_jacobian_add(j,b->variables[2*r+1],b->variables[2*r],c,k)==CHEM_OK);
    }
}
static ChemMechanism *mechanism(void)
{
    ChemSpecies gas[]={{.name="FO2"},{.name="BAR"}},added[2][13]={0};
    const char *const *names[]={stub1,stub2};size_t counts[]={7,13};
    const char *vars[2][4]={{"FO2",stub1[4],stub1[0],stub1[2]},{stub2[4],stub2[5],stub2[8],stub2[9]}};
    ChemNamedJacobianElement elements[2][4];ChemExternalModel models[2];
    for (size_t a=0;a<2;++a) {
        for (size_t i=0;i<counts[a];++i) added[a][i].name=names[a][i];
        for (size_t r=0;r<2;++r) {
            elements[a][2*r]=(ChemNamedJacobianElement){vars[a][2*r],vars[a][2*r]};
            elements[a][2*r+1]=(ChemNamedJacobianElement){vars[a][2*r+1],vars[a][2*r]};
        }
        models[a]=(ChemExternalModel){.additional_species_count=counts[a],.additional_species=added[a],
            .variable_count=4,.variables=vars[a],.parameter_count=a?2:0,.parameters=a?parameter_names:NULL,
            .process_element_count=4,.process_elements=elements[a],.add_forcing=forcing,.process_jacobian=derivative,
            .update_process_parameters=a?update:NULL};
    }
    CHECK(models[0].additional_species_count==7);CHECK(models[1].additional_species_count==13);
    CHECK(models[0].parameter_count==0);CHECK(models[1].parameter_count==2);
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create_with_models(2,gas,0,NULL,2,models,&m)==CHEM_OK);
    CHECK(chem_mechanism_reorder(m)==CHEM_OK);return m;
}
static size_t index_of(ChemState *s,const char *name)
{ size_t index=SIZE_MAX;CHECK(chem_species_index(s->mechanism,name,&index)==CHEM_OK);return index; }
static double value(ChemState *s,const char *name,size_t c) { return s->values[index_of(s,name)*s->cells+c]; }
static void set_values(ChemState *s,const char *name,const double *values,size_t offset)
{ CHECK(chem_state_set_concentration(s,name,values+offset,s->cells)==CHEM_OK); }
static void metadata(ChemState *s)
{
    CHECK(chem_species_count(s->mechanism)==22);CHECK(chem_custom_count(s->mechanism)==2);
    CHECK(chem_external_model_count(s->mechanism)==2);
    for (size_t a=0;a<2;++a) for (size_t i=0;i<(a?13:7);++i) CHECK(index_of(s,(a?stub2:stub1)[i])<22);
    CHECK(index_of(s,"FO2")<22);CHECK(index_of(s,"BAR")<22);
}
static const char *const update_names[]={"FO2","BAR","STUB1.MODE1.QUUX.BAZ","STUB1.MODE2.CORGE.FO2",
    "STUB2.MODE3.CORGE.QUX","STUB2.MODE3.CORGE.BAZ","STUB2.MODE1.NUMBER","STUB2.MODE2.NUMBER"};
static void update_state(ChemState *s,size_t offset)
{
    const double values[][3]={{1.34,1.35,1.36},{2.53,2.54,2.55},{.75,.76,.77},{1.23,1.24,1.25},
        {.42,.43,.44},{.33,.34,.35},{1000,1001,1002},{500,501,502}};
    for (size_t i=0;i<8;++i) {
        set_values(s,update_names[i],values[i],offset);
        const double *named=chem_state_concentration_const_data(s,update_names[i]);CHECK(named!=NULL);
        for (size_t c=0;c<s->cells;++c) {
            CHECK(value(s,update_names[i],c)==values[i][offset+c]);CHECK(named[c]==values[i][offset+c]);
        }
    }
    /* The C equivalent of an aerosol Species/Number object is its unique name. */
    char name[80];snprintf(name,sizeof(name),"%s.MODE%d.%s.%s","STUB2",3,"CORGE","QUX");
    CHECK(!strcmp(name,stub2[12]));CHECK(value(s,name,0)==values[4][offset]);
    snprintf(name,sizeof(name),"%s.MODE%d.NUMBER","STUB2",2);
    CHECK(value(s,name,0)==values[7][offset]);
}
static const double initial[][3]={
    {1,1.5,2},{2,2.5,3},{.5,.7,.9},{.8,1,1.2},{.3,.4,.5},{.2,.3,.4},{1000,1100,1200},{500,550,600},
    {.4,.6,.8},{.1,.15,.2},{.2,.3,.4},{.3,.4,.5},{.01,.015,.02},{275,285,295}};
static void initialize(ChemState *s,size_t offset,int two)
{
    for (size_t i=0;i<8;++i) set_values(s,update_names[i],initial[i],offset);
    if (two) {
        const char *names[]={stub2[4],stub2[5],stub2[8],stub2[9]};
        for (size_t i=0;i<4;++i) set_values(s,names[i],initial[8+i],offset);
        const double inert[]={.2,.25,.3};set_values(s,stub2[11],inert,offset);
        CHECK(chem_state_set_custom(s,parameter_names[0],initial[12]+offset,s->cells)==CHEM_OK);
        for (size_t c=0;c<s->cells;++c) s->conditions[c].temperature=initial[13][offset+c];
    }
}
static void evaluate(ChemState *s,int jacobian)
{
    size_t n=chem_species_count(s->mechanism);double *f=calloc(n*n*s->cells,sizeof(double));CHECK(f!=NULL);
    if (jacobian) chem_external_subtract_jacobian(s,s->values,f);
    else chem_external_add_forcing(s,s->values,f);
    const char *names[]={"FO2",stub1[4],stub1[0],stub1[2]};
    for (size_t r=0;r<2;++r) {
        size_t from=index_of(s,names[r*2]),to=index_of(s,names[r*2+1]);double k=r?.002:.001;
        if (jacobian) { CHECK(f[from*n+from]==k);CHECK(f[to*n+from]==-k); }
        else { CHECK(f[from]==-k*s->values[from]);CHECK(f[to]==k*s->values[from]); }
    }
    free(f);
}
typedef struct { ChemLuAlgorithm algorithm;int copy;double *factors;size_t capacity; } Linear;
static ChemStatus factor(const ChemSparse *p,size_t cells,double *a,void *context)
{
    Linear *l=context;
    if (!l->copy) return chem_lu_factor(p,cells,a,l->algorithm);
    size_t count=chem_sparse_elements(p)*cells;
    if (count>l->capacity) {
        double *v=realloc(l->factors,count*sizeof(double));if (!v) return CHEM_MEMORY;
        l->factors=v;l->capacity=count;
    }
    return chem_lu_factor_copy(p,cells,a,l->factors,l->algorithm);
}
static void solve(const ChemSparse *p,size_t cells,const double *a,double *b,void *context)
{ Linear *l=context;chem_lu_solve(p,cells,l->copy?l->factors:a,b); }
static void integrate(ChemState *s,size_t offset,int two,int euler,ChemRosenbrockMethod method,Linear *linear)
{
    ChemSolver *solver=NULL;CHECK(chem_solver_create(s->mechanism,s->cells,&solver)==CHEM_OK);
    chem_solver_set_linear_solver(solver,(ChemLinearSolver){factor,solve,linear});
    ChemSolverStats stats;chem_update_rates(s);
    if (euler) {
        ChemBackwardEulerParameters p=chem_backward_euler_parameters();
        CHECK(chem_solve_backward_euler(solver,s,10,&p,&stats)==CHEM_OK);
    } else {
        ChemRosenbrockParameters p=chem_rosenbrock_parameters(method);
        CHECK(chem_solve(solver,s,10,&p,&stats)==CHEM_OK);
    }
    double tolerance=euler?(two?.1:2e-4):two && method==CHEM_ROS2?5e-4:5e-5;
    for (size_t c=0;c<s->cells;++c) {
        size_t i=offset+c;
        double d1=initial[0][i]*(1-exp(-.001*10)),d2=initial[2][i]*(1-exp(-.002*10));
        NEAR(value(s,"FO2",c),initial[0][i]-d1,tolerance*initial[0][i]);
        CHECK(value(s,"BAR",c)==initial[1][i]);
        NEAR(value(s,stub1[0],c),initial[2][i]-d2,tolerance*initial[2][i]);
        NEAR(value(s,stub1[2],c),d2,tolerance*initial[2][i]);
        NEAR(value(s,stub1[4],c),initial[3][i]+d1,tolerance*initial[0][i]);
        CHECK(value(s,stub2[0],c)==initial[6][i]);CHECK(value(s,stub2[3],c)==initial[7][i]);
        if (two) {
            double d3=initial[8][i]*(1-exp(-initial[12][i]*10)),d4=initial[10][i]*(1-exp(-initial[13][i]*.005*10));
            NEAR(value(s,stub2[4],c),initial[8][i]-d3,tolerance*initial[8][i]);
            NEAR(value(s,stub2[5],c),initial[9][i]+d3,tolerance*initial[8][i]);
            NEAR(value(s,stub2[8],c),initial[10][i]-d4,tolerance*initial[10][i]);
            NEAR(value(s,stub2[9],c),initial[11][i]+d4,tolerance*initial[10][i]);
        } else { CHECK(value(s,stub2[12],c)==initial[4][i]);CHECK(value(s,stub2[11],c)==initial[5][i]); }
    }
    chem_solver_free(solver);
}
static void all_policies(int which,int euler)
{
    /* Template layouts map to runtime batches (unbounded, 1, 2, 3, 4), both
     * original LU algorithms, and separate/in-place factor ownership. Sparse
     * ordering is canonical in the integrator; CSR/CSC view layouts have their
     * own utility cases. As in the source, vector policies use ROS3; the other
     * Rosenbrock methods solve the entire fixture in a single state. */
    ChemMechanism *m=mechanism();size_t total=which==2 || which>=7?3:1;
    for (int method=0;method<(euler?1:5);++method) for (int algorithm=0;algorithm<2;++algorithm)
        for (int copy=0;copy<2;++copy) for (size_t group=0;group<=(euler || method==CHEM_ROS3?4u:0u);++group) {
            Linear linear={(ChemLuAlgorithm)algorithm,copy,NULL,0};size_t width=group?group:total;
            for (size_t offset=0;offset<total;offset+=width) {
                size_t cells=total-offset<width?total-offset:width;ChemState *s=NULL;
                CHECK(chem_state_create(m,cells,&s)==CHEM_OK);
                if (which==0) metadata(s);
                else if (which==1 || which==2) update_state(s,offset);
                else {
                    int two=which==6 || which==8;initialize(s,offset,two);
                    if (which==3 || which==4) evaluate(s,which==4);
                    else {
                        int before=test_failures;
                        integrate(s,offset,two,euler,(ChemRosenbrockMethod)method,&linear);
                        if (test_failures!=before) fprintf(stderr,"aerosol case=%d euler=%d method=%d algorithm=%d copy=%d group=%zu offset=%zu cells=%zu\n",which,euler,method,algorithm,copy,group,offset,cells);
                    }
                }
                chem_state_free(s);
            }
            free(linear.factors);
        }
    chem_mechanism_free(m);
}
void test_aerosol(void)
{
    const char *names[]={"StateIncludesStubAerosolModel","CanUpdateStateWithStubAerosolModel",
        "CanUpdateMultiCellStateWithStubAerosolModel","CanCalculateSingleGridCellForcingWithStubAerosolModel",
        "CanCalculateSingleGridCellJacobianWithStubAerosolModel","CanSolveSingleGridCellWithStubAerosolModel1",
        "CanSolveSingleGridCellWithTwoStubAerosolModels","CanSolveMultiGridCellWithStubAerosolModel1",
        "CanSolveMultiGridCellWithTwoStubAerosolModels"};
    for (int euler=0;euler<2;++euler) for (int i=0;i<9;++i) {
        char id[160];snprintf(id,sizeof(id),"aerosol_model_%s:AerosolModelIntegration.%s",euler?"backward_euler":"rosenbrock",names[i]);
        RUN(id,all_policies(i,euler));
    }
}
