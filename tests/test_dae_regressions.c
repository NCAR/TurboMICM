#include "chem.h"
#include <math.h>
#include <stdio.h>

typedef struct {int status; double y[4]; ChemStatistics stats;} Result;
static Result run(int system,int method,double balance_atol,double k,double k1,double k2)
{
    const char *simple_names[]={"A","B","C"},*pool_names[]={"gas","aq","P"},*cascade_names[]={"gas","aq","bq","P"};
    const char **names=system==0?simple_names:system==1?pool_names:cascade_names;
    size_t n=system==2?4:3; double gas=1.380649e-23*6.02214076e23;
    ChemTerm rs[]={{system==0?0:system==1?1:2,1}},ps[]={{system==0?1:system==1?2:3,1}};
    ChemReaction reaction={.type=CHEM_ARRHENIUS,.reactant_count=1,.product_count=1,.parameter=-1,.p={k,0,0,300}};
    ChemMechanism m={.species_count=n,.reaction_count=1,.species_names=names,.reactions=&reaction,.reactants=rs,.products=ps};
    ChemTerm mass[]={{0,1},{1,1},{2,1},{3,1}}, r1[]={{0,1}},p1[]={{1,1}},r2[]={{1,1}},p2[]={{2,1}};
    ChemConstraint c[3]={0};size_t count;
    if(system==0) {c[0]=(ChemConstraint){"mass",CHEM_LINEAR_CONSTRAINT,2,mass,NULL,3,0,1e-6,0,0,0};count=1;}
    else {
        c[0]=(ChemConstraint){"eq1",CHEM_EQUILIBRIUM_CONSTRAINT,1,r1,p1,1,1,k1,0,gas,298.15};
        if(system==2) {c[1]=(ChemConstraint){"eq2",CHEM_EQUILIBRIUM_CONSTRAINT,2,r2,p2,1,1,k2,0,gas,298.15};count=3;}
        else count=2;
        c[count-1]=(ChemConstraint){"mass",CHEM_LINEAR_CONSTRAINT,0,mass,NULL,n,0,1e-6,0,0,0};
    }
    ChemWorkspace *w=chem_workspace_create_constrained(&m,c,count);Result result={0};
    if(!w) {result.status=1;return result;}
    double atols[]={1e-12,1e-12,1e-12,1e-12};
    if(system>0) {atols[0]=balance_atol;atols[1]=1e-8;if(system==2) atols[2]=1e-8;}
    ChemRosenbrockParameters p=method==0?chem_rosenbrock_three_stage():method==1?chem_rosenbrock_four_stage():method==2?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
    p.absolute_tolerances=atols;p.relative_tolerance=1e-6;
    if(system==0) {result.y[0]=.9e-6;result.y[2]=.1e-6;}
    else if(system==1) {result.y[0]=1e-6/(1+k1);result.y[1]=k1*result.y[0];}
    else {result.y[0]=1e-6/(1+k1+k1*k2);result.y[1]=k1*result.y[0];result.y[2]=k2*result.y[1];}
    double rates[]={k};ChemConditions conditions={298,101325,0,0};
    result.status=chem_rosenbrock_solve(&m,w,&p,conditions,rates,30,result.y,&result.stats);
    chem_workspace_free(w);return result;
}
static int check_simple(Result r)
{
    double *y=r.y;
    return !r.status && fabs(y[0]+y[1]+y[2]-1e-6)<=1e-12 && y[2]>=-1e-18 &&
           y[0]<1e-12 && fabs(y[1]-.9e-6)<=1e-10 && fabs(y[2]-.1e-6)<=1e-10;
}
static int check_pool(Result r)
{
    double *y=r.y;
    return !r.status && fabs(y[0]+y[1]+y[2]-1e-6)<=1e-12 && y[0]>=-1e-18 && y[1]>=-1e-18 && fabs(y[2]-1e-6)<=1e-8;
}
static int terminator_robertson(void)
{
    const char *names[]={"Cl2","Cl","A","B","C"},*parameters[]={"terminator_k1","robertson_r1","robertson_r2","robertson_r3"};
    const ChemTerm r[]={{0,1},{1,1},{1,1},{2,1},{3,1},{3,1},{3,1},{4,1}};
    const ChemTerm q[]={{1,2},{0,1},{3,1},{3,1},{4,1},{2,1},{4,1}};
    const ChemReaction rx[]={
        {CHEM_USER,0,1,0,1,0,{1}},{CHEM_ARRHENIUS,1,2,1,1,-1,{1,0,0,300}},
        {CHEM_USER,3,1,2,1,1,{1}},{CHEM_USER,4,2,3,2,2,{1}},{CHEM_USER,6,2,5,2,3,{1}}
    };
    const ChemMechanism m={5,5,4,names,parameters,NULL,NULL,rx,r,q,0,NULL};
    const ChemTerm mass[]={{2,1},{3,1},{4,1}};
    const ChemConstraint constraint={"mass",CHEM_LINEAR_CONSTRAINT,4,mass,NULL,3,0,1,0,0,0};
    ChemWorkspace *w=chem_workspace_create_constrained(&m,&constraint,1);if(!w)return 0;
    ChemRosenbrockParameters p=chem_rosenbrock_four_stage_dae();double atol[]={1e-20,1e-20,1e-3,1e-3,1e-3};
    p.absolute_tolerances=atol;p.relative_tolerance=1e-8;
    ChemConditions conditions={298,101300,42,1};double custom[]={0,.04,3e7,1e4},rates[5],y[]={1.8e-10,1.2e-6,1,0,0};
    double total=2*y[0]+y[1],dt=1;int ok=!chem_rate_constants(&m,conditions,custom,rates);
    for(int i=0;i<12 && ok;++i){double advanced=0;while(advanced<dt){ChemStatistics stats;int status=chem_rosenbrock_solve(&m,w,&p,conditions,rates,dt-advanced,y,&stats);
            if((status && status!=2) || !(stats.final_time>0)){ok=0;break;}advanced+=stats.final_time;}
        if(fabs(y[2]+y[3]+y[4]-1)>1e-10 || fabs(2*y[0]+y[1]-total)>1e-6*total)ok=0;
        for(int j=0;j<5;++j)if(!isfinite(y[j]) || y[j]<0)ok=0;
        if(i==0 && !(y[2]<.99 && y[3]>0 && y[4]>1e-3))ok=0;
        dt*=10;
    }
    if(!(y[4]>.95 && y[2]<.05 && y[0]>.4*total && y[1]<.1*total))ok=0;
    chem_workspace_free(w);return ok;
}
int main(void)
{
    int simple=check_simple(run(0,2,0,1e4,0,0));
    int pool=check_pool(run(1,2,1e-8,1e3,10,0));
    int all=1;for(int method=0;method<4;++method) all &= check_simple(run(0,method,0,1e4,0,0));
    Result loose=run(2,2,1e-3,100,2,2),tight=run(2,2,1e-8,100,2,2),ultra=run(2,2,1e-8,1e4,100,50);
    int sensitive=!loose.status && !tight.status && tight.stats.accepted>loose.stats.accepted;
    int overshoot=!ultra.status && ultra.y[0]>=-1e-8 && fabs(ultra.y[3]-1e-6)<=1e-8 && fabs(ultra.y[0])<=1e-10;
    printf("dae_constraint_overshoot:DAEConstraintOvershoot.AlgebraicVariableStaysNonNegative %d\n",simple);
    printf("dae_constraint_overshoot:DAEConstraintOvershoot.EquilibriumPlusConservation %d\n",pool);
    printf("dae_constraint_overshoot:DAEConstraintOvershoot.AllRosenbrockOrdersConstrained %d\n",all);
    printf("dae_algebraic_error_insensitivity:DAEAlgebraicError.ErrorSensitiveToBalanceAtol %d\n",sensitive);
    printf("dae_algebraic_error_insensitivity:DAEAlgebraicError.AlgebraicVariableDoesNotOvershootDeeply %d\n",overshoot);
    printf("linear_constraint:DAESolveWithConstraint.TerminatorAndRobertson %d\n",terminator_robertson());
    return 0;
}
