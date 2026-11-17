#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int flow_tube(int backward)
{
    const char *names[]={"SOA1","SOA2","a-pinene","O3"},*parameters[]={"r1","r2"};
    const ChemTerm reactants[]={{0,1},{1,1},{2,1},{3,1}},products[]={{0,.18},{1,.09}};
    const ChemReaction reactions[]={
        {CHEM_USER,0,1,0,0,0,{1}},{CHEM_USER,1,1,0,0,1,{1}},
        {CHEM_ARRHENIUS,2,2,0,2,-1,{8.8e-17*1e-6*6.02214076e23,0,0,300}}
    };
    const ChemMechanism m={4,3,2,names,parameters,NULL,NULL,reactions,reactants,products,0,NULL};
    ChemConditions c={298.15,101325,0,0};double custom[]={.01,.05},rates[3],y[]={0,0,8e-8,.00002};
    const char *path=backward?"tests/regression_expected/backward_euler_wall_loss.csv":"tests/regression_expected/rosenbrock_wall_loss.csv";
    FILE *f=fopen(path,"r");if(!f)return 0;char line[256];int ok=fgets(line,sizeof line,f)!=NULL && !strcmp(line,"time,SOA1,SOA2,a-pinene,O3\n");
    ChemWorkspace *w=chem_workspace_create(&m);if(!w){fclose(f);return 0;}
    ChemRosenbrockParameters rp=chem_rosenbrock_two_stage();ChemBackwardEulerParameters bp=chem_backward_euler_defaults();
    if(chem_rate_constants(&m,c,custom,rates))ok=0;
    for(int step=0;step<=3600 && ok;++step){double t,expected[4];
        if(!fgets(line,sizeof line,f) || sscanf(line,"%lf,%lf,%lf,%lf,%lf",&t,&expected[0],&expected[1],&expected[2],&expected[3])!=5){ok=0;break;}
        if(step<10)for(int j=0;j<4;++j)if(fabs(y[j]-expected[j])>1e-10)ok=0;
        if(step<3600){ChemStatistics stats;int status=backward?chem_backward_euler_solve(&m,w,&bp,c,rates,1,y,&stats):chem_rosenbrock_solve(&m,w,&rp,c,rates,1,y,&stats);if(status)ok=0;}
    }
    if(fgets(line,sizeof line,f))ok=0;
    fclose(f);chem_workspace_free(w);return ok;
}
int main(void)
{
    printf("regression:RosenbrockRegression.FlowTube %d\n",flow_tube(0));
    printf("regression:BackwardEulerRegression.FlowTube %d\n",flow_tube(1));
    return 0;
}
