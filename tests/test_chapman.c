#include "chem.h"
#include "sparse.h"
#include "state.h"
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
/* Independent net stoichiometry used as the reference for the reaction-set
 * forcing and derivative checks. Rows for the five carrier species are zero. */
static const double net[4][7]={{-1,-1,0,0,0,1,0},{1,1,-1,-1,2,0,1},{0,0,2,-1,-1,1,1},{0,0,-1,1,0,-1,-1}};
static void reference_force_jacobian(const double *k,const double *y,double *f,double *j)
{
    double flux[]={k[0]*y[4]*y[5],k[1]*y[5]*y[7],k[2]*y[6]*y[8],k[3]*y[0]*y[6]*y[7],k[4]*y[7],k[5]*y[8],k[6]*y[8]};
    double d[7][9]={{0}};
    d[0][4]=k[0]*y[5];d[0][5]=k[0]*y[4];d[1][5]=k[1]*y[7];d[1][7]=k[1]*y[5];
    d[2][6]=k[2]*y[8];d[2][8]=k[2]*y[6];d[3][0]=k[3]*y[6]*y[7];d[3][6]=k[3]*y[0]*y[7];d[3][7]=k[3]*y[0]*y[6];
    d[4][7]=k[4];d[5][8]=k[5];d[6][8]=k[6];
    memset(f,0,9*sizeof(double));memset(j,0,81*sizeof(double));
    for(size_t i=0;i<4;++i) for(size_t r=0;r<7;++r) {
        f[5+i]+=net[i][r]*flux[r];
        for(size_t c=0;c<9;++c)j[(5+i)*9+c]+=net[i][r]*d[r][c];
    }
}
static int regression_rates(size_t width)
{
    const double temperature[]={284.19,215.02,299.31},pressure[]={101245,100789.2,101398};
    const double photo[3][3]={{1e-4,1e-5,1e-6},{3.2e-4,7.3e-5,3.2e-6},{5.2e-4,8.2e-5,4.6e-6}};
    ChemState state={0};ChemStateParameters p={7,9,3,species,photolysis,1e-6,NULL};
    if(chem_state_init(&state,&p,3,width))return 0;
    for(size_t i=0;i<3;++i) {state.conditions[i]=(ChemConditions){temperature[i],pressure[i],0,0};chem_matrix_set_row(&state.custom_parameters,i,photo[i],3);}
    int ok=!chem_state_update_rate_constants(&state,&mechanism,NULL);
    for(size_t i=0;i<3;++i) {
        double t=temperature[i],expected[]={2.15e-11*exp(110/t),3.3e-11*exp(55/t),8e-12*exp(-2060/t),6e-34*pow(t/300,2.4),photo[i][0],photo[i][1],photo[i][2]};
        for(size_t j=0;j<7;++j)if(*chem_matrix_at(&state.rate_constants,i,j)!=expected[j])ok=0;
    }
    chem_state_free(&state);return ok;
}
static int regression_derivatives(size_t width,int derivative)
{
    ChemMatrix variables={0},rate_matrix={0};if(chem_matrix_init(&variables,3,9,width,width!=1,0)||chem_matrix_init(&rate_matrix,3,7,width,width!=1,0))return 0;
    unsigned seed=739;int ok=1;
    for(size_t cell=0;cell<3;++cell) {
        double y[9],rates[7],force[9],expected[9],jac[81],expected_jac[81];
        for(size_t i=0;i<9;++i){seed=1664525u*seed+1013904223u;y[i]=.01+(seed>>8)/16777216.0*3;}
        for(size_t i=0;i<7;++i){seed=1664525u*seed+1013904223u;rates[i]=.01+(seed>>8)/16777216.0*3;}
        chem_matrix_set_row(&variables,cell,y,9);chem_matrix_set_row(&rate_matrix,cell,rates,7);
        chem_forcing(&mechanism,(ChemConditions){0},rates,y,force);chem_jacobian(&mechanism,(ChemConditions){0},rates,y,jac);
        reference_force_jacobian(rates,y,expected,expected_jac);
        if(derivative) for(size_t i=0;i<23;++i) {
            double a=jac[row[i]*9+col[i]],b=expected_jac[row[i]*9+col[i]];
            if(fabs(a-b)>(fabs(a)+fabs(b))*1e-8+1e-12)ok=0;
        } else for(size_t i=0;i<9;++i)if(fabs(force[i]-expected[i])>(fabs(force[i])+fabs(expected[i]))*1e-8+1e-12)ok=0;
    }
    chem_matrix_free(&variables);chem_matrix_free(&rate_matrix);return ok;
}
static void reference_rk4(const double *rates,double *state)
{
    double a[9],b[9],c[9],d[9],trial[9],scratch[81];
    for(int step=0;step<1000;++step) {
        reference_force_jacobian(rates,state,a,scratch);
        for(size_t i=0;i<9;++i)trial[i]=state[i]+.25*a[i];
        reference_force_jacobian(rates,trial,b,scratch);
        for(size_t i=0;i<9;++i)trial[i]=state[i]+.25*b[i];
        reference_force_jacobian(rates,trial,c,scratch);
        for(size_t i=0;i<9;++i)trial[i]=state[i]+.5*c[i];
        reference_force_jacobian(rates,trial,d,scratch);
        for(size_t i=0;i<9;++i)state[i]+=(a[i]+2*b[i]+2*c[i]+d[i])/12;
    }
}
static int regression_solve(size_t width,int method)
{
    const double temperature[]={284.19,215.02,299.31},pressure[]={101245,100789.2,101398};
    const double photo[3][3]={{1e-4,1e-5,1e-6},{3.2e-4,7.3e-5,3.2e-6},{5.2e-4,8.2e-5,4.6e-6}};
    ChemState state={0};ChemStateParameters sp={7,9,3,species,photolysis,1e-6,NULL};if(chem_state_init(&state,&sp,3,width))return 0;
    ChemWorkspace *work[3]={0};double y[27],reference[27],independent[27],rates[21];unsigned seed=2221;int ok=1;
    for(size_t cell=0;cell<3;++cell) {
        state.conditions[cell]=(ChemConditions){temperature[cell],pressure[cell],0,0};
        chem_matrix_set_row(&state.custom_parameters,cell,photo[cell],3);
        for(size_t i=0;i<9;++i){seed=1664525u*seed+1013904223u;y[cell*9+i]=.01+(seed>>8)/16777216.0;}
        chem_matrix_set_row(&state.variables,cell,y+cell*9,9);work[cell]=chem_workspace_create(&mechanism);if(!work[cell])ok=0;
    }
    if(chem_state_update_rate_constants(&state,&mechanism,NULL))ok=0;
    for(size_t i=0;i<3;++i)chem_matrix_get_row(&state.rate_constants,i,rates+7*i,7);
    memcpy(reference,y,sizeof y);memcpy(independent,y,sizeof y);
    for(size_t cell=0;cell<3;++cell)reference_rk4(rates+7*cell,independent+9*cell);
    ChemRosenbrockParameters legacy=chem_rosenbrock_three_stage();legacy.relative_tolerance=1e-4;legacy.h_start=1e-6;legacy.max_steps=100;
    for(size_t cell=0;cell<3 && ok;++cell) {ChemStatistics stats;if(chem_rosenbrock_solve(&mechanism,work[cell],&legacy,state.conditions[cell],rates+7*cell,500,reference+9*cell,&stats))ok=0;}
    ChemRosenbrockParameters p=method==0?chem_rosenbrock_two_stage():method==1?chem_rosenbrock_three_stage():method==2?chem_rosenbrock_four_stage():method==3?chem_rosenbrock_four_stage_dae():chem_rosenbrock_six_stage_dae();
    ChemStatistics stats;if(ok && chem_rosenbrock_solve_group(&mechanism,work,3,&p,state.conditions,rates,500,y,&stats))ok=0;
    double max_initial=0;for(size_t i=0;i<27;++i)if(*chem_matrix_at(&state.variables,i/9,i%9)>max_initial)max_initial=*chem_matrix_at(&state.variables,i/9,i%9);
    double rel=method==0?1e-2:1e-4;
    for(size_t i=0;i<27;++i) {
        if(fabs(y[i]-reference[i])>(fabs(y[i])+fabs(reference[i]))*rel+max_initial*1e-12)ok=0;
        if(fabs(y[i]-independent[i])>(fabs(y[i])+fabs(independent[i]))*rel+max_initial*1e-12)ok=0;
    }
    for(size_t i=0;i<3;++i)chem_workspace_free(work[i]);
    chem_state_free(&state);return ok;
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
    int vr=1,vf=1,vj=1,vs=1;
    for(size_t width=1;width<=4;++width){vr &= regression_rates(width);vf &= regression_derivatives(width,0);vj &= regression_derivatives(width,1);vs &= regression_solve(width,1);}
    printf("regression_test_rosenbrock_p_force:RegressionRosenbrock.RateConstants %d\n",regression_rates(1));
    printf("regression_test_rosenbrock_p_force:RegressionRosenbrock.VectorRateConstants %d\n",vr);
    printf("regression_test_rosenbrock_p_force:RegressionRosenbrock.Forcing %d\n",regression_derivatives(1,0));
    printf("regression_test_rosenbrock_p_force:RegressionRosenbrock.VectorForcing %d\n",vf);
    printf("regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.Jacobian %d\n",regression_derivatives(1,1));
    printf("regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.VectorJacobian %d\n",vj);
    const char *methods[]={"TwoStageSolve","ThreeStageSolve","FourStageSolve","FourStageDASolve","SixStageDASolve"};
    for(int method=0;method<5;++method)printf("regression_test_solve:RegressionRosenbrock.%s %d\n",methods[method],regression_solve(1,method));
    printf("regression_test_solve:RegressionRosenbrock.VectorSolve %d\n",vs);
    return 0;
}
