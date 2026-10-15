/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Ports of analytical_policy.hpp and analytical_surface_rxn_policy.hpp. */
#include "test.h"
#include "chem_solver.h"
#include <stdlib.h>
#include <stdint.h>

typedef struct { ChemLuAlgorithm algorithm;int separate;double *factors;size_t capacity; } LinearPolicy;
static ChemStatus factor(const ChemSparse *p,size_t cells,double *a,void *context)
{
    LinearPolicy *l=context;
    if (!l->separate) return chem_lu_factor(p,cells,a,l->algorithm);
    size_t count=chem_sparse_elements(p)*cells;
    if (count>l->capacity) {
        double *v=realloc(l->factors,count*sizeof(double));if (!v) return CHEM_MEMORY;
        l->factors=v;l->capacity=count;
    }
    return chem_lu_factor_copy(p,cells,a,l->factors,l->algorithm);
}
static void linear_solve(const ChemSparse *p,size_t cells,const double *a,double *b,void *context)
{ LinearPolicy *l=context;chem_lu_solve(p,cells,l->separate?l->factors:a,b); }
static ChemStatus solve(ChemSolver *solver,ChemState *s,double dt,int euler,ChemRosenbrockMethod method)
{
    ChemSolverStats stats;
    if (euler) { ChemBackwardEulerParameters p=chem_backward_euler_parameters();return chem_solve_backward_euler(solver,s,dt,&p,&stats); }
    ChemRosenbrockParameters p=chem_rosenbrock_parameters(method);return chem_solve(solver,s,dt,&p,&stats);
}
static ChemRate rate(int kind,int second,int stiff)
{
    const ChemRateType types[]={CHEM_TROE,CHEM_USER,CHEM_TERNARY,CHEM_TUNNELING,CHEM_ARRHENIUS,CHEM_BRANCHED};
    ChemRate r=chem_rate_default(types[kind]);
    if (kind==0 || kind==2) {
        if (second) r.data.falloff=(ChemFalloff){1.2e-3,1.6,3,136,5,24,.9,.8};
        else r.data.falloff.k0_A=kind==0?4e-11:4e-5;
    } else if (kind==1) r.label=second?"photoB":"photoA";
    else if (kind==3) r.data.tunneling=second?(ChemTunneling){1.2e-4,167,1e8}:(ChemTunneling){4e-3,0,0};
    else if (kind==4) r.data.arrhenius=second?(ChemArrhenius){1.2e-4,stiff?1.6:7,75,50,.5}:(ChemArrhenius){4e-3,0,50,300,0};
    else r.data.branched=(ChemBranched){second,1e-4,204.3,1e-3,2,0,0};
    return r;
}
static void analytical_rates(int kind,int stiff,double t,double p,double density,double *k1,double *k2)
{
    if (kind==0 || kind==2) {
        double k0=kind==0?4e-11:4e-5,ki=1,scale=kind==0?density:1;
        *k1=k0*scale/(1+k0*density/ki)*pow(.6,1/(1+pow(log10(k0*density/ki),2)));
        k0=1.2e-3*exp(3/t)*pow(t/300,1.6);ki=136*exp(24/t)*pow(t/300,5);
        *k2=k0*scale/(1+k0*density/ki)*pow(.9,.8/(.8+pow(log10(k0*density/ki),2)));
    } else if (kind==1) { *k1=2e-3;*k2=3e-3; }
    else if (kind==3) { *k1=4e-3;*k2=1.2e-4*exp(-167/t+1e8/pow(t,3)); }
    else if (kind==4) { *k1=4e-3*exp(50/t);*k2=1.2e-4*exp(75/t)*pow(t/50,stiff?1.6:7)*(1+.5*p); }
    else {
        double a=2e-22*exp(2)*2.45e19,b=.43*pow(293.0/298,-8);
        double z=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)))*(1-1e-3)/1e-3;
        a=2e-22*exp(2)*(density*CHEM_AVOGADRO*1e-6);b=.43*pow(t/298,-8);
        double branch=a/(1+a/b)*pow(.41,1/(1+pow(log10(a/b),2)));
        *k1=1e-4*exp(-204.3/t)*(z/(z+branch));*k2=1e-4*exp(-204.3/t)*(branch/(z+branch));
    }
}
static double chain_tolerance(int kind,int stiff,int euler,ChemRosenbrockMethod method)
{
    const double simple[]={1e-10,2e-6,1e-8,1e-8,1e-9,1e-13};
    const double fast[]={1e-5,2e-5,2e-3,1e-4,2e-5,2e-3};
    const double be_simple[]={1e-6,1e-3,1e-5,1e-3,1e-3,1e-5};
    const double be_fast[]={1e-5,1e-3,1e-2,1e-3,1e-3,1e-2};
    if (euler) return stiff?be_fast[kind]:be_simple[kind];
    if (stiff) return kind==4?(method==CHEM_ROS2?1e-4:method==CHEM_ROS6_DA?1e-5:2e-5):fast[kind];
    if (method==CHEM_ROS2) { if (kind==3) return 2e-5;if (kind==4) return 4e-6;if (kind==5) return 1e-10; }
    return simple[kind];
}
static ChemMechanism *chain_mechanism(int kind,int stiff)
{
    ChemSpecies species[]={{.name=stiff?"A1":"A"},{.name=stiff?"A2":"B"},{.name=stiff?"B":"C"},{.name="C"}};
    size_t reactants[]={0,stiff?1:1,2,0,1};ChemProduct products[]={{stiff?2:1,1},{stiff?2:2,1},{3,1},{1,1},{0,1}};
    ChemReaction reactions[5];
    for (size_t i=0;i<(stiff?5u:2u);++i) {
        ChemRate r=rate(kind,stiff?i==2:i==1,stiff);
        if (kind==1 && stiff && i<2) r.label=i?"photoA2B":"photoA1B";
        if (i>=3) { r=chem_rate_default(CHEM_ARRHENIUS);r.data.arrhenius.A=i==3?4e10:.9*4e10; }
        reactions[i]=(ChemReaction){1,1,reactants+i,products+i,r};
    }
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(stiff?4:3,species,stiff?5:2,reactions,&m)==CHEM_OK);
    CHECK(chem_mechanism_reorder(m)==CHEM_OK);return m;
}
static void chain(int kind,int stiff,int euler,ChemRosenbrockMethod method,size_t width,LinearPolicy *linear)
{
    ChemMechanism *m=chain_mechanism(kind,stiff);size_t a,b,c,a2=0;
    CHECK(chem_species_index(m,stiff?"A1":"A",&a)==CHEM_OK);CHECK(chem_species_index(m,"B",&b)==CHEM_OK);
    CHECK(chem_species_index(m,"C",&c)==CHEM_OK);if (stiff) CHECK(chem_species_index(m,"A2",&a2)==CHEM_OK);
    const double temperatures[]={272.5,254.7,312.6},pressures[]={101253.3,100672.5,101319.8};
    double tolerance=chain_tolerance(kind,stiff,euler,method);
    for (size_t offset=0;offset<3;offset+=width) {
        size_t cells=3-offset<width?3-offset:width;ChemState *s=NULL;ChemSolver *solver=NULL;
        CHECK(chem_state_create(m,cells,&s)==CHEM_OK);CHECK(chem_solver_create(m,cells,&solver)==CHEM_OK);
        chem_solver_set_linear_solver(solver,(ChemLinearSolver){factor,linear_solve,linear});double k1[3],k2[3];
        for (size_t i=0;i<cells;++i) {
            double initial=1-(double)(offset+i)/3;
            s->values[a*cells+i]=stiff?.5*initial:initial;if (stiff) s->values[a2*cells+i]=.5*initial;
            s->conditions[i]=(ChemConditions){temperatures[offset+i],pressures[offset+i],0};chem_ideal_air_density(s->conditions+i);
            analytical_rates(kind,stiff,s->conditions[i].temperature,s->conditions[i].pressure,s->conditions[i].air_density,k1+i,k2+i);
        }
        if (kind==1) {
            CHECK(chem_state_fill_custom(s,stiff?"photoA1B":"photoA",2e-3)==CHEM_OK);
            if (stiff) CHECK(chem_state_fill_custom(s,"photoA2B",2e-3)==CHEM_OK);
            CHECK(chem_state_fill_custom(s,"photoB",3e-3)==CHEM_OK);
        }
        for (size_t step=1;step<1000;++step) {
            chem_update_rates(s);ChemStatus status=solve(solver,s,1,euler,method);CHECK(status==CHEM_OK);
            if (status) break;
            for (size_t i=0;i<cells;++i) {
                if (!stiff) { REL(k1[i],s->rates[i],1e-9);REL(k2[i],s->rates[cells+i],1e-9); }
                double initial=1-(double)(offset+i)/3,t=(double)step;
                double ya=initial*exp(-k1[i]*t);
                double yb=initial*(k1[i]/(k2[i]-k1[i]))*(exp(-k1[i]*t)-exp(-k2[i]*t));
                double yc=initial*(1+(k1[i]*exp(-k2[i]*t)-k2[i]*exp(-k1[i]*t))/(k2[i]-k1[i]));
                NEAR(s->values[a*cells+i]+(stiff?s->values[a2*cells+i]:0),ya,tolerance);
                NEAR(s->values[b*cells+i],yb,tolerance);NEAR(s->values[c*cells+i],yc,tolerance);
            }
        }
        chem_solver_free(solver);chem_state_free(s);
    }
    chem_mechanism_free(m);
}
static void surface(int euler,ChemRosenbrockMethod method,LinearPolicy *linear)
{
    double radius=1e-6/2*exp(5*log(.1)*log(.1)/2);
    double number=6/(CHEM_PI*pow(1e-6,3)*exp(9.0/2*log(.1)*log(.1)))*(.002/1000+.003/1000);
    double speed=sqrt(8*CHEM_GAS_CONSTANT/(CHEM_PI*.04607)*272.5);
    double k=4*number*CHEM_PI*radius*radius/(radius/.95e-5+4/(speed*.02)),dt=.1/k;
    ChemSpecies sp[]={{.name="foo",.molecular_weight=.04607,.diffusion=.95e-5},{.name="bar"},{.name="baz"}};
    ChemRate r=chem_rate_default(CHEM_SURFACE);r.label="foo";r.data.surface=(ChemSurface){.95e-5,.04607,.02,0};
    size_t reactant=0;ChemProduct products[]={{1,1},{2,.4}};ChemReaction reaction={1,2,&reactant,products,r};
    ChemMechanism *m=NULL;ChemState *s=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_create(3,sp,1,&reaction,&m)==CHEM_OK);CHECK(chem_mechanism_reorder(m)==CHEM_OK);
    CHECK(chem_state_create(m,1,&s)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    chem_solver_set_linear_solver(solver,(ChemLinearSolver){factor,linear_solve,linear});
    s->conditions[0]=(ChemConditions){272.5,101253.3,0};
    CHECK(chem_state_fill_custom(s,"foo.effective radius [m]",radius)==CHEM_OK);
    CHECK(chem_state_fill_custom(s,"foo.particle number concentration [# m-3]",number)==CHEM_OK);
    CHECK(chem_state_fill_concentration(s,"foo",1)==CHEM_OK);
    size_t a,b,c;CHECK(chem_species_index(m,"foo",&a)==CHEM_OK);CHECK(chem_species_index(m,"bar",&b)==CHEM_OK);
    CHECK(chem_species_index(m,"baz",&c)==CHEM_OK);const double tolerances[]={1e-2,1e-5,1e-6,1e-5,1e-7};
    for (size_t step=1;step<=10;++step) {
        chem_update_rates(s);CHECK(solve(solver,s,dt,euler,method)==CHEM_OK);NEAR(k,s->rates[0],1e-8);
        double ya=exp(-k*(step*dt)),yb=1-ya,yc=.4*(1-ya),tol=euler?.05:tolerances[method];
        REL(s->values[a],ya,tol);REL(s->values[b],yb,tol);REL(s->values[c],yc,tol);
    }
    chem_solver_free(solver);chem_state_free(s);chem_mechanism_free(m);
}
static void policies(int kind,int stiff,int euler)
{
    for (int method=0;method<(euler?1:5);++method) for (int algorithm=0;algorithm<2;++algorithm)
        for (int separate=0;separate<2;++separate) {
            LinearPolicy linear={(ChemLuAlgorithm)algorithm,separate,NULL,0};int before=test_failures;
            if (kind==6) surface(euler,(ChemRosenbrockMethod)method,&linear);
            else chain(kind,stiff,euler,(ChemRosenbrockMethod)method,3,&linear);
            if (test_failures!=before) fprintf(stderr,"analytical kind=%d stiff=%d euler=%d method=%d lu=%d separate=%d\n",kind,stiff,euler,method,algorithm,separate);
            free(linear.factors);
        }
    /* The vector policies use ROS3. Runtime batches exercise the same logical
     * three-cell fixture at widths 1, 2, 3 and 4, including partial final groups. */
    if (kind!=6) for (size_t width=1;width<=4;++width) {
        LinearPolicy linear={CHEM_LU_MOZART,0,NULL,0};int before=test_failures;
        chain(kind,stiff,euler,CHEM_ROS3,width,&linear);
        if (test_failures!=before) fprintf(stderr,"analytical kind=%d stiff=%d euler=%d width=%zu\n",kind,stiff,euler,width);
        free(linear.factors);
    }
}
/* Fixed-reference stiff systems retain the source reaction order and state
 * order. In particular, output targets for Robertson and E5 are absolute. */
#include "analytical_reference.h"
typedef struct {
    size_t reactant_count,product_count,reactants[2];
    ChemProduct products[3];
    double rate;
} ClassicReaction;
static int advance(ChemSolver *solver,ChemState *s,double interval,double max_step,int euler,ChemRosenbrockMethod method)
{
    double time=0;
    while (time<interval) {
        ChemSolverStats stats;double dt=fmin(max_step,interval-time);ChemStatus status;
        if (euler) { ChemBackwardEulerParameters p=chem_backward_euler_parameters();status=chem_solve_backward_euler(solver,s,dt,&p,&stats); }
        else { ChemRosenbrockParameters p=chem_rosenbrock_parameters(method);status=chem_solve(solver,s,dt,&p,&stats); }
        CHECK(status==CHEM_OK || status==CHEM_MAX_STEPS);
        CHECK(time+stats.final_time>time);
        if ((status!=CHEM_OK && status!=CHEM_MAX_STEPS) || !(time+stats.final_time>time)) return 0;
        time+=stats.final_time;
    }
    return 1;
}
static void classic(int kind,int euler,ChemRosenbrockMethod method,LinearPolicy *linear)
{
    const double sc=77.27,qc=8.375e-6,wc=.161,k2=1.6e9,k5=1;
    const double scales[]={k5/(sc*wc*k2),sc*k5/(wc*k2),k5*k5/(wc*wc*k2)};
    ClassicReaction robertson[]={
        {1,1,{0},{{1,1}},.04},
        {2,2,{1,1},{{1,1},{2,1}},3e7},
        {2,2,{1,2},{{0,1},{2,1}},1e4}};
    ClassicReaction oregonator[]={
        {1,1,{1},{{0,1}},k5/(sc*wc)},
        {2,1,{0,1},{{3,1}},k2},
        {1,2,{0},{{2,1},{0,2}},sc*k5/wc},
        {2,1,{0,0},{{4,1}},sc*sc*qc*k2/2},
        {1,1,{2},{{1,1}},k5}};
    ClassicReaction hires[]={
        {1,1,{0},{{1,1}},1.71},
        {1,2,{1},{{0,.43/8.75},{3,8.32/8.75}},8.75},
        {1,2,{2},{{0,8.32/10.03},{3,1.71/10.03}},10.03},
        {0,1,{0},{{0,1}},.0007},
        {1,2,{3},{{2,.43/1.12},{5,.69/1.12}},1.12},
        {1,2,{4},{{2,.035/1.745},{5,1.71/1.745}},1.745},
        {1,1,{5},{{4,1}},.43},
        {1,3,{6},{{4,.43/1.81},{5,.69/1.81},{7,1}},1.81},
        {2,1,{5,7},{{6,1}},280}};
    ClassicReaction e5[]={
        {1,2,{0},{{1,1},{2,1}},7.89e-10},
        {2,1,{1,2},{{4,1}},1.13e9},
        {2,1,{0,2},{{3,1}},1.1e7},
        {1,2,{3},{{2,1},{5,1}},1.13e3}};
    ClassicReaction *definitions[]={robertson,oregonator,hires,e5};
    const size_t species_count[]={3,5,8,6},reaction_count[]={3,5,9,4},outputs[]={12,12,2,7};
    const double *references[]={&reference_robertson[0][0],&reference_oregonator[0][0],&reference_hires[0][0],&reference_e5[0][0]};
    const double rtols[]={1e-10,1e-9,1e-6,1e-13},atols[]={1e-12,1e-11,1e-8,1e-17};
    size_t n=species_count[kind],nr=reaction_count[kind];ChemSpecies species[8];ChemReaction reactions[9];
    char names[8][8],labels[9][8];
    for (size_t i=0;i<n;++i) { snprintf(names[i],sizeof(names[i]),"Y%zu",i);species[i]=(ChemSpecies){.name=names[i]}; }
    for (size_t i=0;i<nr;++i) {
        ClassicReaction *r=&definitions[kind][i];snprintf(labels[i],sizeof(labels[i]),"r%zu",i+1);
        ChemRate rate=chem_rate_default(CHEM_USER);rate.label=labels[i];
        reactions[i]=(ChemReaction){r->reactant_count,r->product_count,r->reactants,r->products,rate};
    }
    ChemMechanism *m=NULL;ChemState *s=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_create(n,species,nr,reactions,&m)==CHEM_OK);
    CHECK(chem_state_create(m,1,&s)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    chem_solver_set_linear_solver(solver,(ChemLinearSolver){factor,linear_solve,linear});
    s->relative_tolerance=rtols[kind];
    for (size_t i=0;i<n;++i) { s->absolute_tolerance[i]=atols[kind];s->values[i]=references[kind][i]; }
    for (size_t i=0;i<nr;++i) CHECK(chem_state_fill_custom(s,labels[i],definitions[kind][i].rate)==CHEM_OK);
    if (kind==0) s->conditions[0]=(ChemConditions){272.5,101253.3,1e6};
    if (kind==1) {
        for (size_t i=0;i<3;++i) { s->values[i]*=scales[i];s->absolute_tolerance[i]*=scales[i]; }
        s->absolute_tolerance[3]=s->absolute_tolerance[4]=scales[1];
    }
    if (kind==3) s->absolute_tolerance[0]=s->absolute_tolerance[4]=s->absolute_tolerance[5]=1e-7;
    chem_update_rates(s);double target=kind==0?1:10,current=0;
    for (size_t output=1;output<=outputs[kind];++output) {
        double dt=kind==1?30*wc:kind==2?321.8122+100*(output-1):target-current;
        if (!advance(solver,s,dt,kind==1 && euler?dt/18000:dt,euler,method)) break;
        size_t checked=kind==1?3:kind==3?4:n;
        for (size_t i=0;i<checked;++i) {
            double expected=references[kind][output*n+i]*(kind==1?scales[i]:1);
            if (kind==1) REL(s->values[i],expected,euler?5e-3:1e-6);
            else if (kind==2) NEAR(s->values[i],expected,euler?.1:(method==CHEM_ROS3 || method==CHEM_ROS4?1e-7:1e-6));
            else {
                double err=fabs(s->values[i]-expected),atol=kind==0?3e-4:i==0?1e-6:1e-13,rtol=kind==0?1e-6:1e-10;
                CHECK(err<atol || err<fabs(expected)*rtol);
            }
        }
        current=target;target*=kind==0?10:100;
    }
    chem_solver_free(solver);chem_state_free(s);chem_mechanism_free(m);
}
static void classic_policies(int kind,int euler)
{
    for (int method=0;method<(euler?1:5);++method) for (int algorithm=0;algorithm<2;++algorithm)
        for (int separate=0;separate<2;++separate) {
            LinearPolicy linear={(ChemLuAlgorithm)algorithm,separate,NULL,0};int before=test_failures;
            classic(kind,euler,(ChemRosenbrockMethod)method,&linear);
            if (test_failures!=before) fprintf(stderr,"classic kind=%d euler=%d method=%d lu=%d separate=%d\n",kind,euler,method,algorithm,separate);
            free(linear.factors);
        }
}

void test_analytical(void)
{
    const char *names[]={"Troe","Photolysis","TernaryChemicalActivation","Tunneling","Arrhenius","Branched","SurfaceRxn"};
    for (int euler=0;euler<2;++euler) for (int kind=0;kind<7;++kind) for (int stiff=0;stiff<(kind==6?1:2);++stiff) {
        char id[180];snprintf(id,sizeof(id),"%s:AnalyticalExamples.%s%s",euler?"analytical_backward_euler":"analytical_rosenbrock_integration",
            names[kind],stiff?"SuperStiffButAnalytical":"");
        RUN(id,policies(kind,stiff,euler));
    }
    const char *classic_names[]={"Robertson","Oregonator","HIRES","E5"};
    for (int euler=0;euler<2;++euler) for (int kind=0;kind<4;++kind) {
        if (euler && (kind==0 || kind==3)) continue;
        char id[160];snprintf(id,sizeof(id),"%s:AnalyticalExamples.%s",euler?"analytical_backward_euler":"analytical_rosenbrock_integration",classic_names[kind]);
        RUN(id,classic_policies(kind,euler));
    }

}
