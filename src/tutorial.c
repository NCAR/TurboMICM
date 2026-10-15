/* Tutorial mechanism assembled through the public C API, independently of the
   configuration reader. Numerical inputs follow the reference tutorial. */
#include "chem.h"
#include <stdio.h>
#include <string.h>

static void print_state(ChemState *state, double time)
{
    printf("%6g,",time);
    const double *y=chem_concentrations(state);
    for (size_t i=0;i<7;++i) printf("%s%11.2e",i?",":"",y[i]);
    putchar('\n');
}
int main(int argc,char **argv)
{
    int user=argc==2 && !strcmp(argv[1],"--user");
    if(argc!=1 && !user) return 2;
    const char *names[]={"A","B","C","D","E","F","G"};
    const char *params[]={"C.effective radius [m]","C.particle number concentration [# m-3]",
        "my photolysis rate","my emission rate","my loss rate"};
    const int reactants[]={0,1,1,2,3,4,4,5,2,1};
    const ChemProduct products[]={{1,1},{2,1},{3,1},{4,1},{5,2},{6,1},{6,1},{6,1},{0,1}};
    ChemReaction reactions[10]={0};
    const size_t offsets[]={0,1,2,3,4,5,7,8,0,9};
    for(size_t i=0;i<10;++i) {
        reactions[i].reactant_count=i==5?2:i==8?0:1;
        reactions[i].reactants=reactants+offsets[i];
        reactions[i].product_count=i==9?0:1;
        reactions[i].products=i<9?products+i:NULL;
    }
    /* The final loss has no product; its pointer is unused. */
    reactions[9].products=NULL;
    reactions[0].rate=chem_rate_default(CHEM_ARRHENIUS);
    reactions[0].rate.p[0]=2.15e-4;reactions[0].rate.p[2]=110;
    reactions[1].rate=(ChemRate){.type=CHEM_BRANCHED_ALKOXY,.p={1.2,204.3,1e-3,2}};
    reactions[2].rate=reactions[1].rate;reactions[2].rate.type=CHEM_BRANCHED_NITRATE;
    reactions[3].rate=(ChemRate){.type=CHEM_SURFACE,.p={2.3e2,.025,.9},.parameter=0};
    const double conversion=1e-6*6.02214076e23;
    reactions[4].rate=(ChemRate){.type=CHEM_TERNARY,.p={1.2,2.3,302.3,2.6/conversion,-3.1,402.1,.9,1.2}};
    reactions[5].rate=(ChemRate){.type=CHEM_TROE,.p={1.8e-30*conversion*conversion,-3,0,2.8e-11*conversion,0,0,.6,1}};
    reactions[6].rate=(ChemRate){.type=CHEM_TUNNELING,.p={1.2,2.3,302.3}};
    for(size_t i=7;i<10;++i) {
        reactions[i].rate=chem_rate_default(CHEM_USER);
        reactions[i].rate.parameter=(int)i-5;
    }
    double initial[]={1,0,0,0,0,0,0};
    ChemMechanismData data={.species_count=7,.reaction_count=user?10:7,.parameter_count=user?5:2,
        .species_names=names,.parameter_names=params,.initial_concentrations=initial,
        .temperature=287.45,.pressure=101319.9,.reactions=reactions};
    ChemMechanism *m=NULL;ChemState *state=NULL;
    if(chem_mechanism_create(&data,&m)!=CHEM_OK || chem_state_create(m,1,&state)!=CHEM_OK) return 1;
    chem_parameters(state)[0]=1e-7;chem_parameters(state)[1]=2.5e6;
    if(user) {chem_parameters(state)[3]=1e-20;chem_parameters(state)[4]=1e-20*1e-3;}
    if(chem_calculate_air_density(state)!=CHEM_OK) return 1;
    printf("%6s","time");
    for(size_t i=0;i<7;++i)printf(",%11s",names[i]);
    putchar('\n');print_state(state,0);
    ChemRosenbrock options=chem_rosenbrock_parameters(3);
    double photo=1e-10;
    for(size_t step=0;step<10;++step) {
        if(user) chem_parameters(state)[2]=photo;
        if(chem_update_rates(state)!=CHEM_OK) return 1;
        double elapsed=0;
        while(elapsed<500) {
            ChemStats stats;
            ChemStatus status=chem_solve(state,&options,500-elapsed,&stats);
            if(stats.final_time<=0 || (status!=CHEM_OK && status!=CHEM_MAX_STEPS)) {
                fprintf(stderr,"%s\n",chem_status_string(status));return 1;
            }
            elapsed+=stats.final_time;
        }
        print_state(state,500*(double)(step+1));photo*=1.5;
    }
    chem_state_free(state);chem_mechanism_free(m);
    return 0;
}
