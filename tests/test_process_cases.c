#include "state.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int update_rates(size_t width)
{
    const char *names[]={"foo","bar"}, *params[]={"bar_user","radius","count"};
    ChemTerm reactants[]={{0,1},{-2,1},{-2,1},{0,1}};
    ChemParameterization parameterization[]={ {0,0,0,.82} };
    ChemReaction reactions[]={
        {.type=CHEM_ARRHENIUS,.reactant_begin=0,.reactant_count=2,.parameter=-1,.p={12.2,0,300,300,0}},
        {.type=CHEM_USER,.reactant_begin=2,.reactant_count=1,.parameter=0,.p={1}},
        {.type=CHEM_SURFACE,.reactant_begin=3,.reactant_count=1,.parameter=1,.p={230,.025,1}}
    };
    ChemMechanism m={.species_count=2,.reaction_count=3,.parameter_count=3,.species_names=names,.parameter_names=params,
        .reactions=reactions,.reactants=reactants,.parameterized_species_count=1,.parameterized_species=parameterization};
    ChemStateParameters p={3,2,3,names,params,1e-6,NULL}; ChemState state={0};
    if(chem_state_init(&state,&p,5,width)) return 0;
    for(size_t i=0;i<5;++i) {
        state.conditions[i]=(ChemConditions){285+i,101100+10*i,10+.1*i,1};
        *chem_matrix_at(&state.custom_parameters,i,0)=.01+.001*i;
        *chem_matrix_at(&state.custom_parameters,i,1)=1e-8*(1+.1*i);
        *chem_matrix_at(&state.custom_parameters,i,2)=1e5*(1+.05*i);
    }
    int ok=!chem_state_update_rate_constants(&state,&m,NULL);
    for(size_t i=0;i<5;++i) {
        ChemConditions c=state.conditions[i]; double user=*chem_matrix_at(&state.custom_parameters,i,0);
        double radius=*chem_matrix_at(&state.custom_parameters,i,1),count=*chem_matrix_at(&state.custom_parameters,i,2);
        double gas=1.380649e-23*6.02214076e23,pi=3.14159265358979323846;
        double expected[]={12.2*exp(300/c.temperature)*c.air_density*.82, user*c.air_density*.82,
            4*count*pi*radius*radius/(radius/230+4/sqrt(8*gas*c.temperature/(pi*.025)))};
        for(size_t j=0;j<3;++j) if(fabs(*chem_matrix_at(&state.rate_constants,i,j)-expected[j])>(j==2?1e-10:1e-12)*expected[j]) ok=0;
    }
    chem_state_free(&state); return ok;
}
static int reaction_construction(int copy)
{
    ChemSpecies species[4]; const char *names[]={"O3","NO","NO2","O2"};
    ChemPhaseSpecies phase_species[4];
    for(int i=0;i<4;++i) { chem_species_init(&species[i],names[i]); phase_species[i]=(ChemPhaseSpecies){&species[i],0,0}; }
    ChemTerm reactants[]={{0,1},{1,1}},products[]={{2,1},{3,1}};
    ChemProcess reaction={.reaction={.type=CHEM_ARRHENIUS,.reactant_begin=0,.reactant_count=2,.product_begin=0,.product_count=2,.p={12.2,0,300,300}},.phase={"gas",phase_species,4}};
    ChemProcess assigned=reaction;
    int ok=reaction.reaction.reactant_count==2 && !strcmp(reaction.phase.species[products[0].species].species->name,"NO2") && !strcmp(reaction.phase.name,"gas");
    if(copy) ok &= assigned.reaction.reactant_count==reaction.reaction.reactant_count && assigned.reaction.product_count==reaction.reaction.product_count &&
        !strcmp(assigned.phase.name,reaction.phase.name) && assigned.reaction.type==CHEM_ARRHENIUS && assigned.reaction.p[0]==12.2 &&
        !strcmp(assigned.phase.species[reactants[0].species].species->name,reaction.phase.species[reactants[0].species].species->name);
    for(int i=0;i<4;++i) chem_species_free(&species[i]);
    return ok;
}
int main(void)
{
    int vector=1; for(size_t width=1;width<=4;++width) vector &= update_rates(width);
    printf("process:Process.Matrix %d\n",update_rates(1));
    printf("process:Process.VectorMatrix %d\n",vector);
    printf("process:Process.BuildsChemicalReaction %d\n",reaction_construction(0));
    printf("process:Process.ChemicalReactionCopyAssignmentSucceeds %d\n",reaction_construction(1));
    return 0;
}
