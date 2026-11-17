#include "chem.h"
#include "sparse.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static const char *species[]={"M","Ar","CO2","H2O","N2","O1D","O","O2","O3"};
static const char *reaction_names[]={"O2_1","O3_1","O3_2","N2_O1D_1","O1D_O2_1","O_O3_1","M_O_O2_1"};
static const char *photolysis[]={"O2_1","O3_1","O3_2"};
static const ChemTerm reactants[]={{4,1},{5,1},{5,1},{7,1},{6,1},{8,1},{0,1},{6,1},{7,1},{7,1},{8,1},{8,1}};
static const ChemTerm products[]={{6,1},{4,1},{6,1},{7,1},{7,2},{8,1},{0,1},{6,2},{5,1},{7,1},{6,1},{7,1}};
static const ChemReaction reactions[]={
    {CHEM_ARRHENIUS,0,2,0,2,-1,{2.15e-11,0,110,300,0}},
    {CHEM_ARRHENIUS,2,2,2,2,-1,{3.3e-11,0,55,300,0}},
    {CHEM_ARRHENIUS,4,2,4,1,-1,{8e-12,0,-2060,300,0}},
    {CHEM_ARRHENIUS,6,3,5,2,-1,{6e-34,2.4,0,300,0}},
    {CHEM_USER,9,1,7,1,0,{1}},{CHEM_USER,10,1,8,2,1,{1}},{CHEM_USER,11,1,10,2,2,{1}}
};
static const ChemMechanism mechanism={9,7,3,species,photolysis,NULL,NULL,reactions,reactants,products,0,NULL};
/* The reference sparse column ordering, including the diagonal and fill. */
static const size_t row[]={0,6,7,8,1,2,3,4,5,6,5,6,6,7,8,5,6,7,8,5,6,7,8};
static const size_t col[]={0,0,0,0,1,2,3,4,4,4,5,5,6,6,6,7,7,7,7,8,8,8,8};
static int matrix_solve(int factor_case)
{
    double a[81]={0},l[81]={0},u[81]={0};unsigned char pattern[81]={0};int ok=1;
    for(size_t i=0;i<23;++i)pattern[row[i]*9+col[i]]=1;
    if(factor_case) {
        for(size_t i=0;i<23;++i)a[row[i]*9+col[i]]=-1;
        for(size_t i=0;i<9;++i)a[i*9+i]+=2;
    } else {
        for(size_t i=0;i<9;++i)l[i*9+i]=1;
        for(size_t i=0;i<23;++i) {
            size_t r=row[i],c=col[i];if(r>c)l[r*9+c]=1;else u[r*9+c]=1;
        }
        for(size_t i=0;i<9;++i)for(size_t j=0;j<9;++j)for(size_t k=0;k<9;++k)a[i*9+j]+=l[i*9+k]*u[k*9+j];
        memset(pattern,1,sizeof pattern);
    }
    ChemWorkspace *w=chem_linear_workspace_create(9,pattern);if(!w)return 0;
    if(chem_linear_factor(w,a,CHEM_LU_MOZART))ok=0;
    if(factor_case) {
        double expected[]={1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-2,-1,3,-1,-2,-3,.125};
        if(chem_linear_factors(w,l,u))ok=0;
        for(size_t i=0;i<23;++i) {
            size_t r=row[i],c=col[i];double value=r==c?1/u[r*9+c]:r>c?l[r*9+c]:u[r*9+c];
            if(fabs(value-expected[i])>.01)ok=0;
        }
    } else {
        double rhs[]={.5,.5,.5,.5,.5,.5,.5,.5,.5},expected[]={.5,.5,.5,.5,.5,-.5,-1,.5,0};
        if(chem_linear_solve(w,rhs))ok=0;
        for(size_t i=0;i<9;++i)if(rhs[i]!=expected[i])ok=0;
    }
    chem_workspace_free(w);return ok;
}
static int force(int small)
{
    double rates[7],y[9],f[9],air=small?6e-14:0;for(size_t i=0;i<7;++i)rates[i]=small?3e-8:1;
    for(size_t i=0;i<9;++i)y[i]=small?5e-6:1;
    chem_forcing(&mechanism,(ChemConditions){273.15,1e5,air,1},rates,y,f);
    double expected[]={0,0,0,0,0,small?1.49e-13:-1,small?4.55e-13:3,small?1.5e-13:2,small?-3e-13:-2};
    for(size_t i=0;i<9;++i)if(small?fabs(f[i]-expected[i])>.01:f[i]!=expected[i])return 0;
    return 1;
}
static int matrix_vector(void)
{
    ChemSparseBuilder builder={0};ChemSparseMatrix matrix={0};if(chem_sparse_builder_init(&builder,9))return 0;
    for(size_t i=0;i<23;++i)if(chem_sparse_builder_add(&builder,row[i],col[i]))return 0;
    if(chem_sparse_init(&matrix,&builder,1,1,1,0))return 0;
    chem_sparse_builder_free(&builder);
    for(size_t i=0;i<23;++i)if(i!=0 && i!=4 && i!=5 && i!=6 && i!=7)*chem_sparse_at(&matrix,0,row[i],col[i])=1;
    double x[23],result[23]={0};for(size_t i=0;i<23;++i)x[i]=.5;
    int ok=!chem_sparse_matvec(&matrix,0,x,23,result,23);
    double expected[]={0,0,0,0,0,2,3,2,2};
    for(size_t i=0;i<23;++i)if(fabs(result[i]-(i<9?expected[i]:0))>.01)ok=0;
    chem_sparse_free(&matrix);return ok;
}
static int solve(double scale,int random)
{
    ChemWorkspace *w=chem_workspace_create(&mechanism);if(!w)return 0;
    double y[]={1,.392,.0169,0,32.9,0,0,8.84,0},params[]={1e-4,1e-5,1e-6},rates[7];
    for(size_t i=0;i<9;++i)y[i]=random?((i+1)*.091):y[i]*scale;
    ChemConditions cond={273.15,1e5,2.7e19,1};int ok=!chem_rate_constants(&mechanism,cond,params,rates);
    ChemRosenbrockParameters p=chem_rosenbrock_three_stage();p.relative_tolerance=1e-4;p.h_start=.005;p.h_max=.5;p.max_steps=100;
    ChemStatistics stats;if(chem_rosenbrock_solve(&mechanism,w,&p,cond,rates,1,y,&stats))ok=0;
    if(!random) {
        double expected[]={1,.392,.0169,0,32.9,1.8039e-41,.00176789,8.83912,4.5031e-36};
        for(size_t i=0;i<9;++i)if(fabs(y[i]-expected[i]*scale)>1e-4)ok=0;
    }
    chem_workspace_free(w);return ok;
}
int main(void)
{
    ChemRosenbrockParameters p=chem_rosenbrock_three_stage();
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.DefaultConstructor %d\n",p.stages==3);
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.lin_solve %d\n",matrix_solve(0));
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.ReactionNames %d\n",sizeof reaction_names/sizeof *reaction_names==7 && !strcmp(reaction_names[0],"O2_1"));
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.PhotolysisNames %d\n",sizeof photolysis/sizeof *photolysis==3);
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.SpeciesNames %d\n",mechanism.species_count==9);
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.simple_force %d\n",force(0));
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.smaller_force %d\n",force(1));
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.factored_alpha_minus_jac %d\n",matrix_solve(1));
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.dforce_dy_times_vector %d\n",matrix_vector());
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.Solve %d\n",solve(1,0));
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.solve_10_times_larger %d\n",solve(10,0));
    printf("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.solve_with_random_number_densities %d\n",solve(1,1));
    printf("chapman_ode_solver:RegressionChapmanODESolver.solve_10_times_smaller %d\n",solve(.1,0));
    return 0;
}
