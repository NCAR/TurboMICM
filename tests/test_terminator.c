#include "chem.h"
#include <math.h>
#include <stdio.h>

static const char *names[]={"Cl2","Cl"},*parameter_names[]={"toy_k1"};
static const ChemTerm reactants[]={{0,1},{1,1},{1,1}},products[]={{1,2},{0,1}};
static const ChemReaction reactions[]={
    {CHEM_USER,0,1,0,1,0,{1}},{CHEM_ARRHENIUS,1,2,1,1,-1,{1,0,0,300}}
};
static const ChemMechanism mechanism={2,2,1,names,parameter_names,NULL,NULL,reactions,reactants,products,0,NULL};

static int terminator(size_t cells)
{
    const double pi=acos(-1), initial[8]={.2e-6,.3e-10,1e-6,1e-10,4e-6,3e-10,2e-6,.07e-10};
    ChemWorkspace *work[4]={0};ChemConditions conditions[4];double rates[8],y[8];
    ChemRosenbrockParameters p=chem_rosenbrock_three_stage();p.relative_tolerance=1e-8;p.absolute_tolerance=1e-20;p.max_steps=100000;
    int ok=1;for(size_t i=0;i<cells;++i){work[i]=chem_workspace_create(&mechanism);conditions[i]=(ChemConditions){298,101300,42,1};if(!work[i])ok=0;}
    for(int step=0;step<(int)floor(2*pi/.3) && ok;++step){double lon=step*.3;
        for(size_t cell=0;cell<cells;++cell){double lat=pi/180*(cell*(90.0/cells));double custom=fmax(0,sin(lat)*sin(pi*20/180)+cos(lat)*cos(pi*20/180)*cos(lon-pi*300/180));
            y[2*cell]=initial[2*cell];y[2*cell+1]=initial[2*cell+1];if(chem_rate_constants(&mechanism,conditions[cell],&custom,rates+2*cell))ok=0;}
        ChemStatistics stats;if(ok && chem_rosenbrock_solve_group(&mechanism,work,cells,&p,conditions,rates,30,y,&stats))ok=0;
        for(size_t cell=0;cell<cells;++cell){double r=rates[2*cell]/4,cl_i=initial[2*cell+1],cl2_i=initial[2*cell],cly=cl_i+2*cl2_i;
            double det=sqrt(r*r+2*r*cly),e=exp(-4*det*30),l=det*30>1e-16?(1-e)/det/30:4;
            double cl_f=-l*(cl_i-det+r)*(cl_i+det+r)/(1+e+30*l*(cl_i+r));
            double cl_ref=cl_i+30*cl_f,cl2_ref=cl2_i-15*cl_f;
            if(fabs(y[2*cell+1]-cl_ref)>fabs(cl_ref)*1e-6+1e-14 || fabs(y[2*cell]-cl2_ref)>fabs(cl2_ref)*1e-6+1e-14)ok=0;
        }
    }
    for(size_t i=0;i<cells;++i)chem_workspace_free(work[i]);
    return ok;
}
int main(void)
{
    int ok=1;for(size_t cells=1;cells<=4;++cells)ok &= terminator(cells);
    printf("terminator:RosenbrockSolver.Terminator %d\n",ok);
    printf("terminator:RosenbrockSolver.VectorTerminator %d\n",ok);
    return 0;
}
