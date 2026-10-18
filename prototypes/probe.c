/* Small generic nonlinear conservative mechanism for independent comparisons. */
#include "chem_adaptive.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

static void vector(const char *name,const double *v,size_t n)
{
    printf("\"%s\":[",name);
    for (size_t i=0;i<n;++i) printf("%s%.17g",i ? "," : "",v[i]);
    printf("],");
}
int main(int argc,char **argv)
{
    if (argc<6 || argc>8) return 2;
    int method=atoi(argv[1]);
    size_t age=(size_t)atoi(argv[2]); int controller=atoi(argv[3]);
    double window=strtod(argv[4],NULL),stiff=strtod(argv[5],NULL);
    double fixed=argc>=7 ? strtod(argv[6],NULL) : 0;
    size_t cells=argc==8 ? (size_t)atoi(argv[7]) : 1;
    chem_species sp[3]={
        {.name="A",.absolute_tolerance=1e-10,.initial=.8},
        {.name="B",.absolute_tolerance=1e-10,.initial=.15},
        {.name="C",.absolute_tolerance=1e-10,.initial=.05}};
    size_t reactants[4][2]={{0,0},{1,0},{2,0},{0,1}};
    chem_product products[4]={{1,1},{2,1},{0,1},{1,2}};
    chem_reaction reactions[4]={0};
    double constants[]={2*stiff,.7,.1,.3};
    if (fixed>0) constants[3]=0;
    for (size_t i=0;i<4;++i) {
        reactions[i].n_reactants=i==3 ? 2 : 1; reactions[i].reactants=reactants[i];
        reactions[i].n_products=1; reactions[i].products=products+i;
        reactions[i].rate=chem_rate_default(CHEM_ARRHENIUS);
        reactions[i].rate.p.arrhenius.A=constants[i];
    }
    chem_mechanism mechanism={.n_species=3,.n_reactions=4,.species=sp,.reactions=reactions,
        .conditions={300,101325,40},.preserve_order=1};
    chem_model m={0}; chem_state s={0}; chem_adaptive_workspace w={0};
    if (chem_model_create(&m,&mechanism,CHEM_LU_MOZART)!=CHEM_OK ||
        chem_state_create(&s,&m,cells)!=CHEM_OK || chem_adaptive_create(&w,&s)!=CHEM_OK) return 1;
    chem_calculate_rates(&s); s.rtol=1e-5;
    chem_adaptive_parameters p=chem_adaptive_defaults();
    p.ros=(chem_ros_method)(method<5 ? method : 0);
    p.method=method<5 ? CHEM_ADAPT_ROS : (method<8 ? CHEM_ADAPT_EXTRAP : CHEM_ADAPT_SDIRK);
    p.extrap_order=method<5 ? 3 : (size_t)method-(method<8 ? 3 : 6);
    if (method==10) { p.method=CHEM_ADAPT_EXTRAP; p.extrap_order=4; p.adaptive_order=1; }
    p.h_start=.1; p.jacobian_age=age; p.controller=controller; p.reuse_window=window;
    if (fixed>0) { p.h_start=p.h_max=fixed; p.growth=1; s.rtol=1e6; }
    chem_ros_parameters ros=chem_ros_parameters_default(p.ros);
    printf("{"); vector("a",ros.a,15); vector("c",ros.c,15); vector("m",ros.m,6); vector("e",ros.e,6);
    vector("gamma",ros.gamma,6); printf("\"order\":%.17g,\"stages\":%zu,\"new_function\":[",ros.order,ros.stages);
    for (size_t i=0;i<6;++i) printf("%s%d",i ? "," : "",ros.new_function[i]);
    printf("],");
    chem_result r=chem_adaptive_solve(&s,3,&p,&w);
    vector("y",s.y,3*cells);
    printf("\"status\":%d,\"counts\":[%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64
        ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64
        ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 "]}\n",r.status,r.stats.internal_steps,
        r.stats.accepted,r.stats.rejected,r.stats.function_calls,r.stats.jacobian_updates,
        r.stats.decompositions,r.stats.solves,w.jacobian_reuses,w.factor_reuses,w.nonlinear_iterations,
        w.extrap_orders[0],w.extrap_orders[1],w.extrap_orders[2]);
    chem_adaptive_destroy(&w); chem_state_destroy(&s); chem_model_destroy(&m);
    return r.status!=CHEM_OK;
}
