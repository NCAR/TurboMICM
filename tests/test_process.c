/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Process and process-set policy ports, including accumulated forcing/Jacobians. */
#include "test.h"
#include "chem_constraints.h"
#include "chem_jacobian.h"
#include <stdint.h>
#include <stdlib.h>
#define EQ(a,b) do { double expected_eq_=(b);NEAR(a,expected_eq_,4*fabs(nextafter(expected_eq_,INFINITY)-expected_eq_)); } while(0)
static uint64_t random_state=1;
static double uniform(void) { random_state=random_state*48271%2147483647;return (double)random_state/2147483647; }
static double draw(double mean,double sigma) { return exp(mean+sigma*sqrt(-2*log(uniform()))*cos(2*CHEM_PI*uniform())); }
static double parameter(const ChemConditions *c,void *context) { return c->air_density*(*(double *)context); }
static void process_rates(size_t group)
{
    double coefficient=.82;ChemSpecies species[]={{.name="foo",.molecular_weight=.025,.diffusion=2.3e2},
        {.name="bar",.parameter=parameter,.parameter_context=&coefficient}};
    ChemRate rates[]={chem_rate_default(CHEM_ARRHENIUS),chem_rate_default(CHEM_SURFACE),chem_rate_default(CHEM_USER)};
    rates[0].data.arrhenius.A=12.2;rates[0].data.arrhenius.C=300;
    rates[1].label="foo_surf";rates[1].data.surface.diffusion=2.3e2;rates[1].data.surface.molecular_weight=.025;
    rates[2].label="bar_user";
    size_t both[]={0,1},foo=0,bar=1;ChemReaction reactions[]={{2,0,both,NULL,rates[0]},
        {1,0,&foo,NULL,rates[1]},{1,0,&bar,NULL,rates[2]}};
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(2,species,3,reactions,&m)==CHEM_OK);random_state=1;
    for (size_t start=0;start<5;start+=group) {
        size_t cells=group<5-start?group:5-start;ChemState *s=NULL;CHECK(chem_state_create(m,cells,&s)==CHEM_OK);
        size_t user,radius,number;CHECK(chem_custom_index(m,"bar_user",&user)==CHEM_OK);
        CHECK(chem_custom_index(m,"foo_surf.effective radius [m]",&radius)==CHEM_OK);
        CHECK(chem_custom_index(m,"foo_surf.particle number concentration [# m-3]",&number)==CHEM_OK);
        for (size_t c=0;c<cells;++c) {
            s->conditions[c]=(ChemConditions){draw(0,.01)*285,draw(0,.01)*101100,draw(0,.01)*10};
            s->custom[user*cells+c]=draw(0,.01)*1e-2;s->custom[radius*cells+c]=draw(0,.01)*1e-8;
            s->custom[number*cells+c]=draw(0,.01)*1e5;
        }
        chem_update_rates(s);
        for (size_t c=0;c<cells;++c) {
            double t=s->conditions[c].temperature,rho=s->conditions[c].air_density;
            double expected=12.2*exp(300/t)*(rho*.82);REL(s->rates[c],expected,1e-12);
            expected=s->custom[user*cells+c]*(rho*.82);REL(s->rates[chem_rate_offset(m,CHEM_USER)*cells+c],expected,1e-12);
            double speed=sqrt(8*CHEM_GAS_CONSTANT/(CHEM_PI*.025)*t),r=s->custom[radius*cells+c],num=s->custom[number*cells+c];
            expected=4*num*CHEM_PI*r*r/(r/2.3e2+4/speed);REL(s->rates[chem_rate_offset(m,CHEM_SURFACE)*cells+c],expected,1e-10);
        }
        chem_state_free(s);
    }
    chem_mechanism_free(m);
}
static void process_set(size_t group)
{
    double parameter_value=.72;const char *names[]={"foo","bar","baz","quz","quuz","corge","qux"};
    ChemSpecies species[7]={{0}};for (size_t i=0;i<7;++i) species[i].name=names[i];
    species[6].parameter=parameter;species[6].parameter_context=&parameter_value;
    size_t r0[]={0,2},r1[]={1,6},r2=3,r3[]={2,6};ChemProduct p0[]={{1,1},{4,2.4}},p1[]={{0,1},{3,1.4}},p3[]={{1,1},{3,2.5}};
    ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);rate.data.arrhenius.A=12.2;rate.data.arrhenius.C=300;
    ChemReaction reactions[]={{2,2,r0,p0,rate},{2,2,r1,p1,rate},{1,0,&r2,NULL,rate},{2,2,r3,p3,rate}};
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(7,species,4,reactions,&m)==CHEM_OK);
    size_t used=0;for (size_t i=0;i<7;++i) { CHECK(chem_reactions_use_species(m,names[i])==(i!=5));used+=chem_reactions_use_species(m,names[i]); }CHECK(used==6);
    size_t rows[]={0,0,0,1,1,1,2,2,3,3,3,4,4},cols[]={0,1,2,0,1,2,0,2,1,2,3,0,2};
    size_t count=0;for (size_t i=0;i<6;++i) for (size_t j=0;j<6;++j) if (chem_kinetic_has_element(m,i,j)) {
        CHECK(count<13);if (count<13) { CHECK(rows[count]==i);CHECK(cols[count]==j); }++count;
    }CHECK(count==13);
    for (size_t start=0;start<2;start+=group) {
        size_t cells=group<2-start?group:2-start;ChemState *s=NULL;CHECK(chem_state_create(m,cells,&s)==CHEM_OK);
        CHECK(s->cells==cells);CHECK(chem_species_count(m)==6);double forcing[12],jac[72];
        for (size_t i=0;i<12;++i) forcing[i]=1000;
        for (size_t i=0;i<72;++i) jac[i]=100;
        for (size_t c=0;c<cells;++c) {
            double cell=start+c,rho=70+10*cell;
            for (size_t v=0;v<5;++v) s->values[v*cells+c]=.1*(v+1)+cell;
            s->conditions[c].air_density=rho;s->rates[c]=10+100*cell;s->rates[cells+c]=(20+100*cell)*rho*.72;
            s->rates[2*cells+c]=30+100*cell;s->rates[3*cells+c]=(40+100*cell)*rho*.72;
        }
        chem_add_kinetic_forcing(s,s->values,forcing);chem_subtract_kinetic_jacobian(s,s->values,jac);
        for (size_t c=0;c<cells;++c) {
            double cell=start+c,rho=70+10*cell,a=.1+cell,b=.2+cell,z=.3+cell,q=.4+cell;
            double k0=10+100*cell,k1=20+100*cell,k2=30+100*cell,k3=40+100*cell;
            double expected[]={1000-k0*a*z+k1*rho*.72*b,1000+k0*a*z-k1*b*rho*.72+k3*rho*.72*z,
                1000-k0*a*z-k3*rho*.72*z,1000+k1*rho*.72*b*1.4-k2*q+k3*rho*.72*2.5*z,1000+k0*a*z*2.4};
            for (size_t v=0;v<5;++v) {
                if (v==1 && start+c==1) REL(forcing[v*cells+c],expected[v],1e-12);
                else EQ(forcing[v*cells+c],expected[v]);
            }
            double derivative[]={100+k0*z,100-k1*rho*.72,100+k0*a,100-k0*z,100+k1*rho*.72,
                100-k0*a-k3*rho*.72,100+k0*z,100+k0*a+k3*rho*.72,100-1.4*k1*rho*.72,
                100-2.5*k3*rho*.72,100+k2,100-2.4*k0*z,100-2.4*k0*a};
            for (size_t i=0;i<13;++i) EQ(jac[(rows[i]*6+cols[i])*cells+c],derivative[i]);
        }
        chem_state_free(s);
    }
    chem_mechanism_free(m);
}
static void masking(size_t group)
{
    ChemSpecies species[]={{.name="A"},{.name="B"},{.name="C"},{.name="D"}};
    size_t reactants[]={0,1};ChemProduct products[]={{2,1},{3,1}};ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);
    rate.data.arrhenius.A=12.2;rate.data.arrhenius.C=300;ChemReaction reaction={2,2,reactants,products,rate};
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(4,species,1,&reaction,&m)==CHEM_OK);
    ChemConstraintTerm terms[]={{"B",1},{"D",1}};
    ChemConstraint constraints[]={{.type=CHEM_LINEAR_CONSTRAINT,.name="B",.algebraic_species="B",.reactant_count=1,.reactants=terms},
        {.type=CHEM_LINEAR_CONSTRAINT,.name="D",.algebraic_species="D",.reactant_count=1,.reactants=terms+1}};
    CHECK(chem_mechanism_set_constraints(m,2,constraints)==CHEM_OK);
    for (size_t start=0;start<2;start+=group) {
        size_t cells=group<2-start?group:2-start;ChemState *s=NULL;CHECK(chem_state_create(m,cells,&s)==CHEM_OK);
        double force[8],jac[32];for (size_t i=0;i<8;++i) force[i]=1000;for (size_t i=0;i<32;++i) jac[i]=500;
        for (size_t c=0;c<cells;++c) { for (size_t v=0;v<4;++v) s->values[v*cells+c]=1+v+4*(start+c);s->rates[c]=10+10*(start+c); }
        chem_add_kinetic_forcing(s,s->values,force);chem_subtract_kinetic_jacobian(s,s->values,jac);
        for (size_t c=0;c<cells;++c) {
            double flux=start+c?600:20,da=start+c?120:20,db=start+c?100:10;
            for (size_t row=0;row<4;++row) {
                double sign=row==0?1:row==2?-1:0;
                EQ(force[row*cells+c],1000-sign*flux);
                EQ(jac[(row*4)*cells+c],500+sign*da);EQ(jac[(row*4+1)*cells+c],500+sign*db);
            }
        }
        chem_state_free(s);
    }
    chem_mechanism_free(m);
}
static void forcing_callback(const double *values,double *f,void *ctx) { chem_add_kinetic_forcing(ctx,values,f); }
static void fd_test(void)
{
    ChemSpecies species[]={{.name="A"},{.name="B"},{.name="C"}};size_t r0[]={0,1},r1=2;
    ChemProduct p0={2,2},p1={0,1};ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);
    ChemReaction reactions[]={{2,1,r0,&p0,rate},{1,1,&r1,&p1,rate}};ChemMechanism *m=NULL;ChemState *s=NULL;
    CHECK(chem_mechanism_create(3,species,2,reactions,&m)==CHEM_OK);CHECK(chem_state_create(m,2,&s)==CHEM_OK);
    double y[]={1,.5,2,4,3,1.5},k[]={5,8,10,12};memcpy(s->values,y,sizeof(y));memcpy(s->rates,k,sizeof(k));
    unsigned char pattern[9];for (size_t i=0;i<3;++i) for (size_t j=0;j<3;++j) pattern[i*3+j]=chem_kinetic_has_element(m,i,j);
    ChemSparse *sparse=NULL;CHECK(chem_sparse_create(3,pattern,&sparse)==CHEM_OK);
    double jac[18]={0},fd[18],packed[18]={0};chem_subtract_kinetic_jacobian(s,y,jac);
    for (size_t i=0;i<3;++i) for (size_t j=0;j<3;++j) {
        size_t index=chem_sparse_index(sparse,i,j);if (index==SIZE_MAX) continue;
        for (size_t c=0;c<2;++c) packed[index*2+c]=jac[(i*3+j)*2+c];
    }
    CHECK(chem_finite_difference(3,2,forcing_callback,s,y,1e-8,fd)==CHEM_OK);
    CHECK(chem_compare_jacobian(sparse,2,packed,fd,1e-7,1e-7).passed);CHECK(chem_check_sparsity(sparse,2,fd,1e-6).passed);
    chem_sparse_free(sparse);chem_state_free(s);chem_mechanism_free(m);
}
static void random_system(size_t cells,size_t nr,size_t ns)
{
    ChemSpecies *species=calloc(ns,sizeof(*species));char (*names)[32]=calloc(ns,sizeof(*names));
    ChemReaction *reactions=calloc(nr,sizeof(*reactions));size_t (*reactants)[3]=calloc(nr,sizeof(*reactants));
    ChemProduct (*products)[10]=calloc(nr,sizeof(*products));CHECK(species && names && reactions && reactants && products);
    random_state=1;
    for (size_t i=0;i<ns;++i) { snprintf(names[i],32,"%zu",i);species[i].name=names[i]; }
    for (size_t i=0;i<nr;++i) {
        size_t nreact=(size_t)(uniform()*4),nprod=(size_t)(uniform()*11);
        for (size_t j=0;j<nreact;++j) reactants[i][j]=(size_t)(uniform()*ns);
        for (size_t j=0;j<nprod;++j) products[i][j]=(ChemProduct){(size_t)(uniform()*ns),1.2};
        ChemRate rate=chem_rate_default(CHEM_ARRHENIUS);rate.data.arrhenius.A=12.2;rate.data.arrhenius.C=300;
        reactions[i]=(ChemReaction){nreact,nprod,reactants[i],products[i],rate};
    }
    ChemMechanism *m=NULL;ChemState *s=NULL;CHECK(chem_mechanism_create(ns,species,nr,reactions,&m)==CHEM_OK);
    CHECK(chem_state_create(m,cells,&s)==CHEM_OK);double *forcing=malloc(cells*ns*sizeof(double));CHECK(forcing!=NULL);
    for (size_t i=0;i<cells*ns;++i) { s->values[i]=draw(-2,4);forcing[i]=1000; }
    for (size_t i=0;i<cells*nr;++i) s->rates[i]=draw(-2,4);
    chem_add_kinetic_forcing(s,s->values,forcing);
    /* The original random-system case exercises construction and execution only. */
    free(forcing);chem_state_free(s);chem_mechanism_free(m);free(species);free(names);free(reactions);free(reactants);free(products);
}
static void random_suite(void) { random_system(200,50,40);random_system(300,30,20);random_system(400,100,80); }
void test_process(void)
{
    RUN("process:Process.Matrix",process_rates(5));
    RUN("process:Process.VectorMatrix",for (size_t g=1;g<=4;++g) process_rates(g));
    RUN("process_set:ProcessSet.Matrix",process_set(2));
    RUN("process_set:ProcessSet.VectorMatrix",for (size_t g=1;g<=4;++g) process_set(g));
    RUN("process_set:RandomProcessSet.Matrix",random_suite());
    RUN("process_set:ProcessSetAlgebraicVariables.CudaMatrix",for (size_t g=1;g<=4;++g) masking(g));
    RUN("process_set:ProcessSetFiniteDifferenceJacobian.Matrix",fd_test());
}
