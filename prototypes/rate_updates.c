#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"Failed line %d: %s\n",__LINE__,#x); exit(1); } } while (0)

static double callback(const chem_conditions *c,const double *custom,void *context)
{ (void)custom; (void)context; return c->air_density; }

int main(void)
{
    chem_species species[3]={{.name="A",.initial=1},{.name="B"},
        {.name="M",.type=CHEM_SPECIES_THIRD_BODY}};
    size_t reactants[4][2]={{0,2},{1,0},{0,0},{1,0}};
    chem_product products[4]={{1,1},{0,1},{1,1},{0,1}};
    chem_reaction reactions[4]={0};
    chem_rate_type types[4]={CHEM_USER,CHEM_ARRHENIUS,CHEM_SURFACE,CHEM_USER};
    for (size_t i=0;i<4;++i) {
        reactions[i].n_reactants=i ? 1 : 2; reactions[i].reactants=reactants[i];
        reactions[i].n_products=1; reactions[i].products=products+i;
        reactions[i].rate=chem_rate_default(types[i]);
    }
    reactions[0].rate.label="light"; reactions[0].rate.p.user.scaling=2;
    reactions[2].rate.label="surface";
    reactions[2].rate.p.surface.diffusion=1e-5;
    reactions[2].rate.p.surface.molecular_weight=.028;
    reactions[2].rate.p.surface.probability=.2;
    reactions[3].rate.label="fixed";
    chem_mechanism mechanism={.n_species=3,.n_reactions=4,.species=species,.reactions=reactions,
        .conditions={300,101325,40},.preserve_order=1};
    for (int parameter_callback=0;parameter_callback<2;++parameter_callback) {
        if (parameter_callback) { species[2].type=CHEM_SPECIES_PARAMETER; species[2].parameter=callback; }
        chem_model m={0}; chem_state initial={0};
        CHECK(chem_model_create(&m,&mechanism,CHEM_LU_MOZART)==CHEM_OK);
        CHECK(chem_state_create(&initial,&m,3)==CHEM_OK);
        CHECK(m.n_custom==4);
        for (size_t j=0;j<m.n_custom;++j) for (size_t c=0;c<3;++c)
            initial.custom[j*3+c]=.1*(1+c);
        for (size_t c=0;c<3;++c) {
            initial.custom[m.custom_offset[2]*3+c]=1e-6*(1+c);
            initial.custom[(m.custom_offset[2]+1)*3+c]=1e6*(1+c);
        }
        chem_calculate_rates(&initial);
        for (size_t i=0;i<12;++i) CHECK(isfinite(initial.rates[i]) && initial.rates[i]>0);
        for (unsigned bits=0;bits<16;++bits) {
            chem_state full={0},partial={0}; unsigned char changed[4]; size_t indices[4];
            CHECK(chem_state_copy(&full,&initial)==CHEM_OK);
            CHECK(chem_state_copy(&partial,&initial)==CHEM_OK);
            for (size_t j=0;j<4;++j) {
                changed[j]=(bits>>j)&1;
                if (changed[j]) for (size_t c=0;c<3;++c) {
                    full.custom[j*3+c]*=2; partial.custom[j*3+c]*=2;
                }
            }
            size_t count=chem_rate_update_indices(&m,changed,indices);
            size_t expected=parameter_callback ? 4 : (size_t)changed[m.custom_offset[0]]+
                (size_t)(changed[m.custom_offset[2]] || changed[m.custom_offset[2]+1])+changed[m.custom_offset[3]];
            CHECK(count==expected);
            chem_calculate_rates(&full);
            for (size_t i=0;i<12;++i) CHECK(isfinite(full.rates[i]) && full.rates[i]>0);
            CHECK(chem_calculate_selected_rates(&partial,count,indices)==CHEM_OK);
            CHECK(!memcmp(full.rates,partial.rates,12*sizeof(double)));
            size_t invalid=4;
            CHECK(chem_calculate_selected_rates(&partial,1,&invalid)==CHEM_INVALID);
            CHECK(!memcmp(full.rates,partial.rates,12*sizeof(double)));
            chem_state_destroy(&full); chem_state_destroy(&partial);
        }
        chem_state_destroy(&initial); chem_model_destroy(&m);
    }
    species[2].type=CHEM_SPECIES_THIRD_BODY; species[2].parameter=NULL;
    reactions[3].rate=chem_rate_default(CHEM_LAMBDA);
    reactions[3].rate.label="fixed"; reactions[3].rate.p.lambda.function=callback;
    chem_model m={0};
    CHECK(chem_model_create(&m,&mechanism,CHEM_LU_MOZART)==CHEM_OK);
    CHECK(chem_model_set_lambda(&m,"fixed",callback,NULL)==CHEM_OK);
    unsigned char unchanged[4]={0}; size_t indices[4];
    CHECK(chem_rate_update_indices(&m,unchanged,indices)==4);
    CHECK(chem_rate_update_indices(NULL,unchanged,indices)==CHEM_NONE);
    chem_model_destroy(&m);
    puts("Rate dependencies: all custom masks, third bodies, surface slots, callbacks and invalid indices pass");
    return 0;
}
