#include "chem.h"

#include <stdio.h>
#include <string.h>

extern const ChemMechanism tutorial_no_user;
extern const ChemMechanism tutorial_user;

static int set(ChemState *s, const ChemMechanism *m, const char *name, double value)
{
    int i=chem_find_parameter(m,name);
    if (i<0) return 0;
    chem_parameter_values(s,(size_t)i)[0]=value;
    return 1;
}

static void print_state(ChemState *s, const ChemMechanism *m, int time)
{
    printf("%6d",time);
    for (size_t i=0;i<m->species_count;++i) printf(",%11.2e",chem_species_values(s,i)[0]);
    putchar('\n');
}

int main(int argc, char **argv)
{
    if (argc!=2 || (strcmp(argv[1],"user") && strcmp(argv[1],"no_user"))) return 2;
    int user=!strcmp(argv[1],"user");
    const ChemMechanism *m=user ? &tutorial_user : &tutorial_no_user;
    ChemState *s=chem_state_create(m,1); if (!s) return 1;
    chem_temperatures(s)[0]=287.45; chem_pressures(s)[0]=101319.9;
    chem_species_values(s,0)[0]=1;
    if (!set(s,m,"C surface.effective radius [m]",1e-7) ||
        !set(s,m,"C surface.particle number concentration [# m-3]",2.5e6)) return 1;
    if (user && (!set(s,m,"my emission rate",1e-20) || !set(s,m,"my loss rate",1e-23))) return 1;
    printf("%6s","time");
    for (size_t i=0;i<m->species_count;++i) printf(",%11s",m->species[i].name);
    putchar('\n');
    print_state(s,m,0);
    ChemRosenbrock p=chem_rosenbrock_three_stage();
    double photo=1e-10;
    for (int step=1;step<=10;++step) {
        ChemStats stats;
        if (user && !set(s,m,"my photolysis rate",photo)) return 1;
        if (chem_update_rates(s)!=CHEM_OK || chem_solve(s,&p,500,&stats)!=CHEM_OK) return 1;
        print_state(s,m,step*500);
        photo*=1.5;
    }
    chem_state_free(s);
    return 0;
}
