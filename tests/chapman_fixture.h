#ifndef CHAPMAN_FIXTURE_H
#define CHAPMAN_FIXTURE_H
#include "chem_process.h"
#include "test_backend.h"
extern const chem_mechanism chapman_mechanism;
static const char *const chapman_names[]={"M","Ar","CO2","H2O","N2","O1D","O","O2","O3"};
typedef struct {
    chem_species species[9]; size_t reactants[7][3]; chem_product products[7][2]; chem_reaction reactions[7];
    chem_mechanism mechanism; chem_model model; chem_state state; integration_backend backend;
} chapman_fixture;
static inline chem_status chapman_create(chapman_fixture *f,size_t cells,size_t width,int generated,int preserve,int integration)
{
    *f=(chapman_fixture){0};
    for (size_t i=0;i<9;++i) f->species[i].name=chapman_names[i];
    const size_t reactants[7][3]={{5,4,0},{5,7,0},{6,8,0},{6,7,0},{7,0,0},{8,0,0},{8,0,0}};
    const chem_product products[7][2]={{{6,1},{4,1}},{{6,1},{7,1}},{{7,2},{0,0}},{{8,1},{0,1}},{{6,2},{0,0}},{{5,1},{7,1}},{{6,1},{7,1}}};
    memcpy(f->reactants,reactants,sizeof(reactants)); memcpy(f->products,products,sizeof(products));
    const double A[]={2.15e-11,3.3e-11,8e-12,6e-34},C[]={110,55,-2060,0};
    const char *labels[]={"jO2","jO3a","jO3b"};
    for (size_t i=0;i<7;++i) {
        chem_rate rate=chem_rate_default(i<4 ? CHEM_ARRHENIUS : CHEM_USER);
        if (i<4) { rate.p.arrhenius.A=A[i]; rate.p.arrhenius.C=C[i]; if (i==3) { rate.p.arrhenius.B=integration ? 0 : 2.4; rate.p.arrhenius.C=integration ? 2.4 : 0; } }
        else rate.label=labels[i-4];
        f->reactions[i]=(chem_reaction){i<3 ? 2 : i==3 ? 3 : 1,i==2 || i==4 ? 1 : 2,f->reactants[i],f->products[i],rate};
    }
    f->mechanism=(chem_mechanism){.n_species=9,.n_reactions=7,.species=f->species,.reactions=f->reactions,.preserve_order=preserve};
    if (generated) f->mechanism=chapman_mechanism;
    chem_status status=chem_model_create(&f->model,&f->mechanism,CHEM_LU_DOOLITTLE);
    if (status==CHEM_OK) status=chem_state_create(&f->state,&f->model,cells);
    if (status==CHEM_OK) status=integration_backend_create(&f->backend,&f->state,width,CHEM_CSR,0);
    return status;
}
static inline void chapman_destroy(chapman_fixture *f)
{ integration_backend_destroy(&f->backend); chem_state_destroy(&f->state); chem_model_destroy(&f->model); }
static inline void chapman_conditions(chem_state *s)
{
    const double photo[3][3]={{1e-4,1e-5,1e-6},{3.2e-4,7.3e-5,3.2e-6},{5.2e-4,8.2e-5,4.6e-6}};
    const double T[]={284.19,215.02,299.31},P[]={101245,100789.2,101398};
    for (size_t c=0;c<3;++c) {
        s->conditions[c]=(chem_conditions){T[c],P[c],0};
        for (size_t j=0;j<3;++j) s->custom[j*3+c]=photo[c][j];
    }
    chem_calculate_rates(s);
}
#endif
