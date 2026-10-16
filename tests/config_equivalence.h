#ifndef CHEM_CONFIG_EQUIVALENCE_H
#define CHEM_CONFIG_EQUIVALENCE_H
#include "tutorials.h"
extern const chem_mechanism robertson_mechanism, rates_mechanism, user_rates_mechanism;
static int same_reaction(const chem_reaction *a,const chem_reaction *b)
{
    if (a->n_reactants!=b->n_reactants || a->n_products!=b->n_products || a->rate.type!=b->rate.type) return 0;
    for (size_t i=0;i<a->n_reactants;++i) if (a->reactants[i]!=b->reactants[i]) return 0;
    for (size_t i=0;i<a->n_products;++i) if (a->products[i].species!=b->products[i].species || a->products[i].yield!=b->products[i].yield) return 0;
    return 1;
}
static void configuration_equivalence(example_kind kind,const chem_mechanism *description)
{
    example_chemistry hand; REQUIRE(example_create(&hand,kind)==CHEM_OK);
    for (size_t algorithm=0;algorithm<2;++algorithm) for (size_t kernels=0;kernels<2;++kernels) {
        chem_model model={0}; chem_state a={0},b={0};
        REQUIRE(chem_model_create(&model,description,(chem_lu_algorithm)algorithm)==CHEM_OK);
        if (!kernels) { model.forcing=chem_forcing; model.jacobian=chem_jacobian; }
        REQUIRE(chem_state_create(&a,&hand.model,3)==CHEM_OK); REQUIRE(chem_state_create(&b,&model,3)==CHEM_OK);
        REQUIRE(example_initialize(&a,kind,0)==CHEM_OK); REQUIRE(example_initialize(&b,kind,1)==CHEM_OK);
        CHECK(a.model->n_species==b.model->n_species); CHECK(a.model->n_reactions==b.model->n_reactions);
        CHECK(a.model->n_custom==b.model->n_custom);
        size_t map[7],reaction_map[10];
        for (size_t i=0;i<a.model->n_species;++i) {
            const char *name=a.model->mechanism->species[a.model->var_to_spec[i]].name;
            map[i]=chem_species_index(b.model,name); REQUIRE(map[i]!=CHEM_NONE);
            for (size_t c=0;c<3;++c) { double y=.001*(i+1)*(c+1); a.y[i*3+c]=y; b.y[map[i]*3+c]=y; }
        }
        for (size_t i=0;i<a.model->n_reactions;++i) {
            reaction_map[i]=CHEM_NONE;
            for (size_t j=0;j<b.model->n_reactions;++j) if (same_reaction(hand.reactions+i,description->reactions+j)) {
                CHECK(reaction_map[i]==CHEM_NONE); reaction_map[i]=j;
            }
            REQUIRE(reaction_map[i]!=CHEM_NONE);
        }
        for (size_t c=0;c<3;++c) {
            a.conditions[c].temperature+=c*20; a.conditions[c].pressure-=c*10000;
            chem_ideal_air_density(a.conditions+c); b.conditions[c]=a.conditions[c];
        }
        chem_calculate_rates(&a); chem_calculate_rates(&b);
        for (size_t i=0;i<a.model->n_reactions;++i) for (size_t c=0;c<3;++c) {
            double expected=a.rates[i*3+c]; NEAR(b.rates[reaction_map[i]*3+c],expected,1e-30+2e-14*fabs(expected));
        }
        double af[21]={0},bf[21]={0},aj[147]={0},bj[147]={0};
        chem_state_forcing(&a,a.y,af); chem_state_forcing(&b,b.y,bf);
        chem_state_jacobian(&a,a.y,aj); chem_state_jacobian(&b,b.y,bj);
        for (size_t i=0;i<a.model->n_species;++i) for (size_t c=0;c<3;++c) NEAR(bf[map[i]*3+c],af[i*3+c],1e-20+1e-13*fabs(af[i*3+c]));
        for (size_t i=0;i<a.model->n_species;++i) for (size_t j=0;j<a.model->n_species;++j) {
            size_t ia=a.model->lu.map[i*a.model->n_species+j],ib=b.model->lu.map[map[i]*b.model->n_species+map[j]];
            for (size_t c=0;c<3;++c) {
                double va=ia==CHEM_NONE ? 0 : aj[ia*3+c],vb=ib==CHEM_NONE ? 0 : bj[ib*3+c];
                NEAR(vb,va,1e-20+1e-13*fabs(va));
            }
        }
        /* Repeat the hand-built tutorial trajectories, with independent cells. */
        memset(a.y,0,3*a.model->n_species*sizeof(double)); memset(b.y,0,3*b.model->n_species*sizeof(double));
        REQUIRE(example_initialize(&a,kind,0)==CHEM_OK); REQUIRE(example_initialize(&b,kind,1)==CHEM_OK);
        double initial[]={1,2,.5}; REQUIRE(chem_state_assign_species(&a,"A",initial,3)==CHEM_OK); REQUIRE(chem_state_assign_species(&b,"A",initial,3)==CHEM_OK);
        chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3); double photo=1e-10;
        for (size_t step=0;step<10;++step) {
            if (kind==EXAMPLE_USER_RATES) {
                REQUIRE(chem_state_set_custom(&a,"my photolysis rate",&photo,1)==CHEM_OK);
                REQUIRE(chem_state_set_custom(&b,"my photolysis rate",&photo,1)==CHEM_OK); photo*=1.5;
            }
            chem_calculate_rates(&a); chem_calculate_rates(&b);
            double dt=kind==EXAMPLE_ROBERTSON ? 200 : 500;
            chem_result ra=chem_ros_solve(&a,dt,&p,NULL,NULL),rb=chem_ros_solve(&b,dt,&p,NULL,NULL);
            REQUIRE(ra.status==CHEM_OK && rb.status==CHEM_OK); chem_state_clamp(&a); chem_state_clamp(&b);
            CHECK(ra.stats.accepted==rb.stats.accepted); CHECK(ra.stats.rejected==rb.stats.rejected);
            for (size_t i=0;i<a.model->n_species;++i) for (size_t c=0;c<3;++c)
                NEAR(b.y[map[i]*3+c],a.y[i*3+c],1e-30+1e-10*fabs(a.y[i*3+c]));
        }
        chem_state_destroy(&a); chem_state_destroy(&b); chem_model_destroy(&model);
    }
    chem_model_destroy(&hand.model);
}
#endif
