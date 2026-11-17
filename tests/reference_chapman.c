#include "chem.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static const ChemSpecies species[]={
    {"M",0,1e-3,0,0},{"Ar",0,1e-3,0,0},{"CO2",0,1e-3,0,0},{"H2O",0,1e-3,0,0},
    {"N2",0,1e-3,0,0},{"O1D",0,1e-3,0,0},{"O",0,1e-3,0,0},{"O2",0,1e-3,0,0},{"O3",0,1e-3,0,0}};
static const int reactants[]={4,5,5,7,6,8,0,6,7,7,8,8};
static const ChemProduct products[]={{6,1},{4,1},{6,1},{7,1},{7,2},{8,1},{0,1},{6,2},{5,1},{7,1},{6,1},{7,1}};
static const char *const names[]={"jO2","jO3a","jO3b","k3","k4","k5","k6"};
static const char *const reaction_names[]={"O2_1","O3_1","O3_2","N2_O1D_1","O1D_O2_1","O_O3_1","M_O_O2_1"};
static const char *const photo_names[]={"O2_1","O3_1","O3_2"};
static const double defaults[]={0,0,0,0,0,0,0};
static const int starts[]={0,2,4,6,9,10,11}, counts[]={2,2,2,3,1,1,1};
static const int pstarts[]={0,2,4,5,7,8,10}, pcounts[]={2,2,1,2,1,2,2};
static unsigned long long random_state=1;
static double uniform_random(void)
{
    random_state=(random_state*16807u)%2147483647u;
    return (double)random_state/2147483647.0;
}
static double lognormal_random(double sigma)
{
    double a=uniform_random(),b=uniform_random();
    return exp(-2+sigma*sqrt(-2*log(a))*cos(6.2831853071795864769*b));
}
typedef struct { ChemReaction reactions[7]; ChemMechanism mechanism; ChemState *state; } Chapman;
static int init(Chapman *c,size_t cells,int arbitrary_rates)
{
    memset(c,0,sizeof(*c));
    for (int i=0;i<7;++i) {
        ChemRate rate={0};
        if (arbitrary_rates) { rate.kind=CHEM_USER; rate.parameter=i; rate.a=1; }
        else if (i<4) {
            const double a[]={2.15e-11,3.3e-11,8e-12,6e-34}, b[]={0,0,0,2.4}, q[]={110,55,-2060,0};
            rate.kind=CHEM_ARRHENIUS; rate.a=a[i]; rate.b=b[i]; rate.c=q[i]; rate.d=300;
        } else { rate.kind=CHEM_USER; rate.parameter=i-4; rate.a=1; }
        c->reactions[i]=(ChemReaction){(size_t)starts[i],(size_t)counts[i],(size_t)pstarts[i],(size_t)pcounts[i],0,rate};
    }
    c->mechanism=(ChemMechanism){9,7,arbitrary_rates ? 7 : 3,species,c->reactions,reactants,products,names,defaults,298,101325};
    c->state=chem_state_create(&c->mechanism,cells); return c->state!=NULL;
}
static void reference_rates(double temperature,const double photo[3],double k[7])
{
    k[0]=2.15e-11*exp(110/temperature); k[1]=3.3e-11*exp(55/temperature);
    k[2]=8e-12*exp(-2060/temperature); k[3]=6e-34*pow(temperature/300,2.4);
    for (int i=0;i<3;++i) k[4+i]=photo[i];
}
static void reference_force(const double k[7],const double y[9],double f[9])
{
    memset(f,0,9*sizeof(*f));
    f[5]=k[5]*y[8]-k[0]*y[4]*y[5]-k[1]*y[5]*y[7];
    f[6]=2*k[4]*y[7]+k[6]*y[8]+k[0]*y[4]*y[5]+k[1]*y[5]*y[7]-k[2]*y[6]*y[8]-k[3]*y[0]*y[6]*y[7];
    f[7]=-k[4]*y[7]+k[5]*y[8]+k[6]*y[8]+2*k[2]*y[6]*y[8]-k[3]*y[0]*y[6]*y[7];
    f[8]=-k[5]*y[8]-k[6]*y[8]-k[2]*y[6]*y[8]+k[3]*y[0]*y[6]*y[7];
}
static void reference_jacobian(const double k[7],const double y[9],double j[81])
{
    memset(j,0,81*sizeof(*j));
    j[5*9+4]=-k[0]*y[5]; j[5*9+5]=-k[0]*y[4]-k[1]*y[7]; j[5*9+7]=-k[1]*y[5]; j[5*9+8]=k[5];
    j[6*9+0]=-k[3]*y[6]*y[7]; j[6*9+4]=k[0]*y[5]; j[6*9+5]=k[0]*y[4]+k[1]*y[7];
    j[6*9+6]=-k[2]*y[8]-k[3]*y[0]*y[7]; j[6*9+7]=2*k[4]+k[1]*y[5]-k[3]*y[0]*y[6]; j[6*9+8]=k[6]-k[2]*y[6];
    j[7*9+0]=-k[3]*y[6]*y[7]; j[7*9+6]=2*k[2]*y[8]-k[3]*y[0]*y[7];
    j[7*9+7]=-k[4]-k[3]*y[0]*y[6]; j[7*9+8]=k[5]+k[6]+2*k[2]*y[6];
    j[8*9+0]=k[3]*y[6]*y[7]; j[8*9+6]=-k[2]*y[8]+k[3]*y[0]*y[7];
    j[8*9+7]=k[3]*y[0]*y[6]; j[8*9+8]=-k[5]-k[6]-k[2]*y[6];
}
static int rate_check(void)
{
    Chapman c; if (!init(&c,3,0)) return 0;
    const double t[]={284.19,215.02,299.31},p[]={101245,100789.2,101398};
    const double photo[][3]={{1e-4,1e-5,1e-6},{3.2e-4,7.3e-5,3.2e-6},{5.2e-4,8.2e-5,4.6e-6}};
    for (size_t cell=0;cell<3;++cell) {
        chem_temperatures(c.state)[cell]=t[cell]; chem_pressures(c.state)[cell]=p[cell];
        for (size_t j=0;j<3;++j) chem_parameter_values(c.state,j)[cell]=photo[cell][j];
    }
    int ok=chem_update_rates(c.state)==CHEM_OK;
    for (size_t cell=0;cell<3;++cell) {
        double expected[7]; reference_rates(t[cell],photo[cell],expected);
        for (size_t r=0;r<7;++r) {
            double actual=chem_rate_constants(c.state,r)[cell];
            ok &= actual==expected[r];
        }
    }
    chem_state_free(c.state); return ok;
}
static int forcing_or_jacobian(int jacobian)
{
    Chapman c; if (!init(&c,3,1)) return 0;
    random_state=739;
    for (size_t cell=0;cell<3;++cell) {
        for (size_t i=0;i<9;++i) chem_species_values(c.state,i)[cell]=lognormal_random(jacobian ? 4 : 2);
        for (size_t r=0;r<7;++r) chem_parameter_values(c.state,r)[cell]=lognormal_random(jacobian ? 4 : 2);
    }
    int ok=chem_update_rates(c.state)==CHEM_OK;
    for (size_t cell=0;cell<3;++cell) {
        double y[9],k[7],f[9],j[81],expected[81];
        for (size_t i=0;i<9;++i) y[i]=chem_species_values(c.state,i)[cell];
        for (size_t r=0;r<7;++r) k[r]=chem_rate_constants(c.state,r)[cell];
        ok &= chem_evaluate(c.state,cell,f,j)==CHEM_OK;
        if (jacobian) {
            reference_jacobian(k,y,expected);
            for (size_t i=0;i<81;++i) ok &= fabs(j[i]+expected[i])<=(fabs(j[i])+fabs(expected[i]))*1e-8+1e-12;
        } else {
            reference_force(k,y,expected);
            for (size_t i=0;i<9;++i) ok &= fabs(f[i]-expected[i])<=(fabs(f[i])+fabs(expected[i]))*1e-8+1e-12;
        }
    }
    chem_state_free(c.state); return ok;
}
static int simple_force(double rate,double concentration,int small)
{
    Chapman c; if (!init(&c,1,1)) return 0;
    for (size_t i=0;i<9;++i) chem_species_values(c.state,i)[0]=concentration;
    for (size_t i=0;i<7;++i) chem_parameter_values(c.state,i)[0]=rate;
    chem_set_air_density(c.state,0,small ? 6e-14 : 0);
    double f[9],j[81]; int ok=chem_update_rates(c.state)==CHEM_OK && chem_evaluate(c.state,0,f,j)==CHEM_OK;
    const double expected1[]={0,0,0,0,0,-1,3,2,-2},expected2[]={0,0,0,0,0,1.49e-13,4.55e-13,1.5e-13,-3e-13};
    const double *expected=small ? expected2 : expected1;
    for (size_t i=0;i<9;++i) ok &= small ? fabs(f[i]-expected[i])<=.01 : f[i]==expected[i];
    chem_state_free(c.state); return ok;
}
static const int positions[23][2]={{0,0},{6,0},{7,0},{8,0},{1,1},{2,2},{3,3},{4,4},{5,4},{6,4},
    {5,5},{6,5},{6,6},{7,6},{8,6},{5,7},{6,7},{7,7},{8,7},{5,8},{6,8},{7,8},{8,8}};
typedef struct { size_t row_start[10]; int columns[23],lookup[81]; unsigned char algebraic[9]; double values[23]; ChemLinearSystem view; } Sparse;
static void sparse_init(Sparse *s,double value)
{
    memset(s,0,sizeof(*s)); for (size_t i=0;i<81;++i) s->lookup[i]=-1;
    size_t n=0;
    for (size_t row=0;row<9;++row) {
        s->row_start[row]=n;
        for (size_t col=0;col<9;++col) for (size_t k=0;k<23;++k)
            if (positions[k][0]==(int)row && positions[k][1]==(int)col) {
                s->columns[n]=(int)col; s->lookup[row*9+col]=(int)n; s->values[n++]=value;
            }
    }
    s->row_start[9]=n; s->view=(ChemLinearSystem){9,1,n,s->row_start,s->columns,s->lookup,s->algebraic,s->values};
}
static int lin_solve(void)
{
    Sparse s; sparse_init(&s,1); double rhs[9]; for (size_t i=0;i<9;++i) rhs[i]=.5;
    ChemLinearSolver lu=chem_mozart_linear_solver(); int ok=lu.solve(&s.view,rhs,NULL)==CHEM_OK;
    const double expected[]={.5,.5,.5,.5,.5,-.5,-1,.5,0};
    for (size_t i=0;i<9;++i) ok &= rhs[i]==expected[i];
    return ok;
}
static int factor(void)
{
    Sparse s; sparse_init(&s,-1); ChemLinearSolver lu=chem_mozart_linear_solver();
    int ok=lu.factor(&s.view,2,NULL)==CHEM_OK;
    const double expected[]={1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-2,-1,3,-1,-2,-3,.125};
    for (size_t i=0;i<23;++i) {
        int row=positions[i][0],col=positions[i][1],at=s.lookup[row*9+col]; double v=s.values[at];
        if (row==col) v=1/v; /* reference stores reciprocals of U's diagonal */
        ok &= fabs(v-expected[i])<=.01;
    }
    return ok;
}
static int multiply(void)
{
    Sparse s; sparse_init(&s,1); double input[9],output[23]={0};
    for (size_t i=0;i<9;++i) input[i]=.5;
    const unsigned char mask[]={0,0,0,0,0,1,1,1,1};
    int ok=chem_linear_multiply(&s.view,input,mask,output)==CHEM_OK;
    const double expected[]={0,0,0,0,0,2,3,2,2};
    for (size_t i=0;i<23;++i) ok &= fabs(output[i]-(i<9 ? expected[i] : 0))<=.01;
    return ok;
}
static int old_solve(double scale,int random)
{
    Chapman c; if (!init(&c,1,0)) return 0;
    const double y[]={1,.392,.0169,0,32.9,0,0,8.84,0};
    random_state=1;
    for (size_t i=0;i<9;++i) chem_species_values(c.state,i)[0]=random ? uniform_random() : y[i]*scale;
    chem_temperatures(c.state)[0]=273.15; chem_pressures(c.state)[0]=100000; chem_set_air_density(c.state,0,2.7e19);
    chem_parameter_values(c.state,0)[0]=1e-4; chem_parameter_values(c.state,1)[0]=1e-5; chem_parameter_values(c.state,2)[0]=1e-6;
    ChemRosenbrock p=chem_rosenbrock_three_stage(); ChemStats stats;
    int ok=chem_update_rates(c.state)==CHEM_OK && chem_solve(c.state,&p,1,&stats)==CHEM_OK;
    if (!random) {
        const double expected[]={1,.392,.0169,0,32.9,1.8039e-41,.00176789,8.83912,4.5031e-36};
        for (size_t i=0;i<9;++i) ok &= fabs(chem_species_values(c.state,i)[0]-expected[i]*scale)<=1e-4;
    }
    chem_state_free(c.state); return ok;
}
static void rk4(const double k[7],double y[9],double dt,int steps)
{
    double h=dt/steps;
    for (int step=0;step<steps;++step) {
        double a[9],b[9],c[9],d[9],tmp[9]; reference_force(k,y,a);
        for (size_t i=0;i<9;++i) tmp[i]=y[i]+h*a[i]/2;
        reference_force(k,tmp,b); for (size_t i=0;i<9;++i) tmp[i]=y[i]+h*b[i]/2;
        reference_force(k,tmp,c); for (size_t i=0;i<9;++i) tmp[i]=y[i]+h*c[i];
        reference_force(k,tmp,d); for (size_t i=0;i<9;++i) y[i]+=h*(a[i]+2*b[i]+2*c[i]+d[i])/6;
    }
}
static int regression_solve(int method,size_t cells)
{
    Chapman c; if (!init(&c,cells,0)) return 0;
    const double t[]={284.19,215.02,299.31},p[]={101245,100789.2,101398};
    const double photo[][3]={{1e-4,1e-5,1e-6},{3.2e-4,7.3e-5,3.2e-6},{5.2e-4,8.2e-5,4.6e-6}};
    double expected[3][9],maximum=0;
    random_state=1;
    for (size_t cell=0;cell<cells;++cell) {
        for (size_t i=0;i<9;++i) {
            double value=lognormal_random(2);
            chem_species_values(c.state,i)[cell]=expected[cell][i]=value; if (value>maximum) maximum=value;
        }
        chem_temperatures(c.state)[cell]=t[cell]; chem_pressures(c.state)[cell]=p[cell];
        for (size_t j=0;j<3;++j) chem_parameter_values(c.state,j)[cell]=photo[cell][j];
        double k[7]; reference_rates(t[cell],photo[cell],k); rk4(k,expected[cell],500,500);
    }
    ChemRosenbrock pset=method==0 ? chem_rosenbrock_two_stage() : method==1 ? chem_rosenbrock_three_stage() :
                        method==2 ? chem_rosenbrock_four_stage() : method==3 ? chem_rosenbrock_four_stage_dae() :
                        chem_rosenbrock_six_stage_dae();
    ChemStats stats; int ok=chem_update_rates(c.state)==CHEM_OK && chem_solve(c.state,&pset,500,&stats)==CHEM_OK;
    double tolerance=method==0 ? 1e-2 : 1e-4;
    for (size_t cell=0;cell<cells;++cell) for (size_t i=0;i<9;++i) {
        double actual=chem_species_values(c.state,i)[cell], ref=expected[cell][i];
        if (fabs(actual-ref)>(fabs(actual)+fabs(ref))*tolerance+maximum*1e-12) ok=0;
    }
    if (!ok) fprintf(stderr,"Chapman method=%d cells=%zu final_time=%g steps=%lu\n",method,cells,stats.final_time,stats.internal_steps);
    chem_state_free(c.state); return ok;
}
static int integration(void)
{
    Chapman c; if (!init(&c,1,0)) return 0;
    /* The integration fixture uses C=2.4 for the fourth Arrhenius rate. */
    c.reactions[3].rate.b=0; c.reactions[3].rate.c=2.4;
    const double initial[]={.2,.2,.3,.3,.3,.1,.1,.1,.2};
    for (size_t i=0;i<9;++i) chem_species_values(c.state,i)[0]=initial[i];
    const double photo[]={.1,.2,.3};
    for (size_t i=0;i<3;++i) chem_parameter_values(c.state,i)[0]=photo[i];
    chem_temperatures(c.state)[0]=2; chem_pressures(c.state)[0]=3;
    const double oxygen=initial[5]+initial[6]+2*initial[7]+3*initial[8];
    ChemRosenbrock p=chem_rosenbrock_three_stage(); int ok=1;
    for (size_t step=0;step<100 && ok;++step) {
        ChemStats stats; ok=chem_update_rates(c.state)==CHEM_OK && chem_solve(c.state,&p,30,&stats)==CHEM_OK;
        double total=chem_species_values(c.state,5)[0]+chem_species_values(c.state,6)[0]+
                     2*chem_species_values(c.state,7)[0]+3*chem_species_values(c.state,8)[0];
        ok &= fabs(total-oxygen)<=1e-8*oxygen;
        for (size_t i=0;i<9;++i) ok &= isfinite(chem_species_values(c.state,i)[0]) && chem_species_values(c.state,i)[0]>=0;
        for (size_t i=0;i<5;++i) ok &= fabs(chem_species_values(c.state,i)[0]-initial[i])<=1e-8*initial[i];
    }
    ok &= chem_species_values(c.state,7)[0]<.5*initial[7] && chem_species_values(c.state,6)[0]>initial[6] &&
          chem_species_values(c.state,6)[0]>.5*oxygen;
    chem_state_free(c.state); return ok;
}
static int failures;
static void report(const char *id,int ok)
{ printf("%s %s %s\n",id,ok ? "PASS" : "FAIL",ok ? "reference Chapman quantities matched" : "Chapman quantity differed"); failures+=!ok; }
int main(void)
{
    ChemRosenbrock defaults3=chem_rosenbrock_three_stage();
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.DefaultConstructor",defaults3.stages==3);
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.lin_solve",lin_solve());
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.ReactionNames",sizeof(reaction_names)/sizeof(*reaction_names)==7 && !strcmp(reaction_names[6],"M_O_O2_1"));
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.PhotolysisNames",sizeof(photo_names)/sizeof(*photo_names)==3 && !strcmp(photo_names[2],"O3_2"));
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.SpeciesNames",sizeof(species)/sizeof(*species)==9 && !strcmp(species[8].name,"O3"));
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.simple_force",simple_force(1,1,0));
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.smaller_force",simple_force(3e-8,5e-6,1));
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.factored_alpha_minus_jac",factor());
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.dforce_dy_times_vector",multiply());
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.Solve",old_solve(1,0));
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.solve_10_times_larger",old_solve(10,0));
    report("chapman_ode_solver:RegressionChapmanODESolver.solve_10_times_smaller",old_solve(.1,0));
    report("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral.solve_with_random_number_densities",old_solve(1,1));
    report("regression_test_rosenbrock_p_force:RegressionRosenbrock.RateConstants",rate_check());
    report("regression_test_rosenbrock_p_force:RegressionRosenbrock.VectorRateConstants",rate_check());
    report("regression_test_rosenbrock_p_force:RegressionRosenbrock.Forcing",forcing_or_jacobian(0));
    report("regression_test_rosenbrock_p_force:RegressionRosenbrock.VectorForcing",forcing_or_jacobian(0));
    report("regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.Jacobian",forcing_or_jacobian(1));
    report("regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.VectorJacobian",forcing_or_jacobian(1));
    report("regression_test_solve:RegressionRosenbrock.TwoStageSolve",regression_solve(0,3));
    report("regression_test_solve:RegressionRosenbrock.ThreeStageSolve",regression_solve(1,3));
    report("regression_test_solve:RegressionRosenbrock.FourStageSolve",regression_solve(2,3));
    report("regression_test_solve:RegressionRosenbrock.FourStageDASolve",regression_solve(3,3));
    report("regression_test_solve:RegressionRosenbrock.SixStageDASolve",regression_solve(4,3));
    report("regression_test_solve:RegressionRosenbrock.VectorSolve",regression_solve(1,1) && regression_solve(1,2) && regression_solve(1,3));
    report("chapman_integration:ChapmanIntegration.CanBuildChapmanSystem",integration());
    return failures ? 1 : 0;
}
