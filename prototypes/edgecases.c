#include "chem_adaptive.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"Failed line %d: %s\n",__LINE__,#x); exit(1); } } while (0)

static chem_status bad_factor(chem_state *s,void *context)
{ (void)s; (void)context; return CHEM_INVALID; }
static chem_status unused_solve(const chem_state *s,double *v,void *context)
{ (void)s; (void)v; (void)context; return CHEM_INVALID; }

int main(void)
{
    chem_state empty={0}; chem_adaptive_workspace empty_work={0};
    chem_adaptive_parameters defaults=chem_adaptive_defaults();
    CHECK(chem_adaptive_create(&empty_work,&empty)==CHEM_INVALID);
    CHECK(chem_adaptive_solve(&empty,0,&defaults,&empty_work).status==CHEM_INVALID);
    chem_species sp[2]={{.name="reactant",.initial=1,.absolute_tolerance=1e-18},
                        {.name="product",.initial=0,.absolute_tolerance=1e-18}};
    size_t reactant=0; chem_product product={1,1};
    chem_reaction reaction={.n_reactants=1,.reactants=&reactant,.n_products=1,.products=&product};
    reaction.rate=chem_rate_default(CHEM_ARRHENIUS); reaction.rate.p.arrhenius.A=1e6;
    chem_mechanism mechanism={.n_species=2,.n_reactions=1,.species=sp,.reactions=&reaction,
        .conditions={300,101325,40},.preserve_order=1};
    chem_model model={0}; CHECK(chem_model_create(&model,&mechanism,CHEM_LU_MOZART)==CHEM_OK);
    for (int method=0;method<11;++method) {
        chem_state s={0}; chem_adaptive_workspace w={0};
        CHECK(chem_state_create(&s,&model,3)==CHEM_OK);
        CHECK(chem_adaptive_create(&w,&s)==CHEM_OK); s.rtol=1e-4;
        chem_calculate_rates(&s);
        chem_adaptive_parameters p=chem_adaptive_defaults();
        p.ros=(chem_ros_method)(method<5 ? method : 0);
        p.method=method<5 ? CHEM_ADAPT_ROS : (method<8 ? CHEM_ADAPT_EXTRAP : CHEM_ADAPT_SDIRK);
        p.extrap_order=method<5 ? 3 : (size_t)method-(method<8 ? 3 : 6);
        if (method==10) { p.method=CHEM_ADAPT_EXTRAP; p.extrap_order=4; p.adaptive_order=1; }
        CHECK(chem_adaptive_solve(&s,NAN,&p,&w).status==CHEM_INVALID);
        CHECK(chem_adaptive_solve(&s,-1,&p,&w).status==CHEM_INVALID);
        chem_result r=chem_adaptive_solve(&s,0,&p,&w);
        CHECK(r.status==CHEM_OK && r.stats.internal_steps==0 && s.y[0]==1);
        CHECK(chem_state_set_linear_backend(&s,bad_factor,unused_solve,NULL)==CHEM_OK);
        r=chem_adaptive_solve(&s,1,&p,&w);
        CHECK(r.status==CHEM_INVALID && r.stats.internal_steps==1 && s.y[0]==1);
        CHECK(chem_state_set_linear_backend(&s,NULL,NULL,NULL)==CHEM_OK);
        w.next_h=0;
        r=chem_adaptive_solve(&s,100,&p,&w);
        CHECK(r.status==CHEM_OK && fabs(r.stats.final_time-100)<1e-12);
        printf("method=%d reactant=%.17g product=%.17g clipped_mass=%.17g\n",method,s.y[0],s.y[3],w.clipped_mass);
        double total=0;
        for (size_t c=0;c<3;++c) {
            CHECK(isfinite(s.y[c]) && s.y[c]>=0 && s.y[c]<1e-12);
            CHECK(isfinite(s.y[3+c]) && fabs(s.y[3+c]-1)<1e-3);
            total+=s.y[c]+s.y[3+c];
        }
        /* Projection changes mass. Verify that every bit of that change is
         * accounted for, in addition to testing the requested 0.1% accuracy. */
        CHECK(fabs(total-3-w.clipped_mass)<1e-10);
        chem_adaptive_destroy(&w); chem_state_destroy(&s);
    }
    chem_model_destroy(&model);
    puts("11 methods: zero duration, invalid inputs, backend errors, and 1e8 stiffness interval pass");
    return 0;
}
