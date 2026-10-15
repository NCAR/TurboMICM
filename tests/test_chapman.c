/* SPDX-License-Identifier: Apache-2.0
 * Chapman unit and integration fixtures, evaluated by the generic C library. */
#include "test.h"
#include "chem_solver.h"
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <limits.h>
#include "chapman_oracle.h"

static const char *const names[]={"M","Ar","CO2","H2O","N2","O1D","O","O2","O3"};
static const char *const labels[]={"N2_O1D_1","O1D_O2_1","O_O3_1","M_O_O2_1","O2_1","O3_1","O3_2"};
static ChemMechanism *mechanism(int integration)
{
    ChemSpecies species[9]={0};for (size_t i=0;i<9;++i) species[i].name=names[i];
    size_t reactants[7][3]={{4,5},{5,7},{6,8},{0,6,7},{7},{8},{8}};
    ChemProduct products[7][2]={{{6,1},{4,1}},{{6,1},{7,1}},{{7,2}},{{8,1},{0,1}},{{6,2}},{{5,1},{7,1}},{{6,1},{7,1}}};
    size_t nr[]={2,2,2,3,1,1,1},np[]={2,2,1,2,1,2,2};
    double a[]={2.15e-11,3.3e-11,8e-12,6e-34},c[]={110,55,-2060,0};
    ChemReaction reactions[7];
    for (size_t i=0;i<7;++i) {
        ChemRate rate=chem_rate_default(i<4?CHEM_ARRHENIUS:CHEM_USER);rate.label=labels[i];
        if (i<4) { rate.data.arrhenius.A=a[i];rate.data.arrhenius.C=c[i]; }
        if (i==3) { rate.data.arrhenius.B=integration?0:2.4;rate.data.arrhenius.C=integration?2.4:0; }
        reactions[i]=(ChemReaction){nr[i],np[i],reactants[i],products[i],rate};
    }
    if (integration) {
        const size_t order[]={6,5,7,8,0,1,4,3,2};size_t inverse[9];
        for (size_t i=0;i<9;++i) { species[i].name=names[order[i]];inverse[order[i]]=i; }
        for (size_t i=0;i<7;++i) {
            for (size_t j=0;j<nr[i];++j) reactants[i][j]=inverse[reactants[i][j]];
            for (size_t j=0;j<np[i];++j) products[i][j].species=inverse[products[i][j].species];
        }
    }
    ChemMechanism *m=NULL;CHECK(chem_mechanism_create(9,species,7,reactions,&m)==CHEM_OK);return m;
}
static void metadata(int which)
{
    ChemMechanism *m=mechanism(0);
    if (which==0) { ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS3);CHECK(p.stages==3); }
    if (which==1) CHECK(chem_reaction_count(m)==7);
    if (which==2) CHECK(chem_custom_count(m)==3);
    if (which==3) CHECK(chem_species_count(m)==9);
    chem_mechanism_free(m);
}
static void force(int small)
{
    ChemMechanism *m=mechanism(0);ChemState *s=NULL;CHECK(chem_state_create(m,1,&s)==CHEM_OK);
    for (size_t i=0;i<9;++i) s->values[i]=small?5e-6:1;
    for (size_t i=0;i<7;++i) s->rates[i]=small?3e-8:1;
    s->conditions[0].air_density=small?6e-14:0;
    double f[9];chem_forcing(s,s->values,f);
    for (size_t i=0;i<5;++i) CHECK(f[i]==0);
    double expected[]={-1,3,2,-2},tiny[]={1.49e-13,4.55e-13,1.5e-13,-3e-13};
    for (size_t i=0;i<4;++i) { if (small) NEAR(f[i+5],tiny[i],.01);else CHECK(f[i+5]==expected[i]); }
    chem_state_free(s);chem_mechanism_free(m);
}
/* The original 23-entry CSC fixture includes five inert diagonals. Their
 * placeholders belong to alpha*I and LU, not to the forcing derivative. */
static const size_t row[]={0,6,7,8,1,2,3,4,5,6,5,6,6,7,8,5,6,7,8,5,6,7,8};
static const size_t col[]={0,0,0,0,1,2,3,4,4,4,5,5,6,6,6,7,7,7,7,8,8,8,8};
static ChemSparse *pattern(void)
{
    unsigned char p[81]={0};for (size_t i=0;i<23;++i) p[row[i]*9+col[i]]=1;
    ChemSparse *s=NULL;CHECK(chem_sparse_create(9,p,&s)==CHEM_OK);CHECK(chem_sparse_elements(s)==23);return s;
}
static void linear_algebra(int which)
{
    ChemSparse *s=pattern();double a[23],x[23],b[23]={0};
    for (size_t i=0;i<23;++i) { a[i]=1;x[i]=.5; }
    if (which==0) {
        chem_lu_solve(s,1,a,x);double expected[]={.5,.5,.5,.5,.5,-.5,-1,.5,0};
        for (size_t i=0;i<9;++i) CHECK(x[i]==expected[i]);
    } else if (which==1) {
        for (size_t i=0;i<23;++i) a[chem_sparse_index(s,row[i],col[i])]=row[i]==col[i]?1:-1;
        CHECK(chem_lu_factor(s,1,a,CHEM_LU_MOZART)==CHEM_OK);
        const double expected[]={1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-2,-1,3,-1,-2,-3,.125};
        for (size_t i=0;i<23;++i) {
            double value=a[chem_sparse_index(s,row[i],col[i])];
            if (row[i]==col[i]) value=1/value;
            NEAR(value,expected[i],.01);
        }
    } else {
        for (size_t i=0;i<5;++i) a[chem_sparse_index(s,i,i)]=0;
        chem_sparse_multiply(s,1,a,x,b);
        for (size_t i=0;i<23;++i) NEAR(b[i],i==6?3:i>=5 && i<=8?2:0,.01);
    }
    chem_sparse_free(s);
}
static void solve_unit(int which)
{
    ChemMechanism *m=mechanism(0);ChemState *s=NULL;ChemSolver *solver=NULL;
    CHECK(chem_state_create(m,1,&s)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    double initial[]={1,.392,.0169,0,32.9,0,0,8.84,0};
    double scale=which==1?10:which==2?.1:1;
    uint64_t seed=1;
    for (size_t i=0;i<9;++i) {
        seed=seed*48271%2147483647;
        s->values[i]=which==3?(double)(float)((double)seed/2147483647):initial[i]*scale;
    }
    s->conditions[0]=(ChemConditions){273.15,100000,2.7e19};s->relative_tolerance=1e-4;
    s->custom[0]=1e-4;s->custom[1]=1e-5;s->custom[2]=1e-6;chem_update_rates(s);
    ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS3);p.max_steps=100;p.h_start=1e-6;
    ChemSolverStats stats;CHECK(chem_solve(solver,s,1,&p,&stats)==CHEM_OK);
    if (which!=3) {
        const double expected[3][9]={{1,.392,.0169,0,32.9,1.8039e-41,.00176789,8.83912,4.5031e-36},
            {10,3.92,.169,0,329,1.8039e-38,.0176789,88.3912,4.5031e-33},
            {.1,.0392,.00169,0,3.29,1.8039e-44,.000176789,.883912,4.5031e-39}};
        for (size_t i=0;i<9;++i) NEAR(s->values[i],expected[which][i],1e-4);
    }
    chem_solver_free(solver);chem_state_free(s);chem_mechanism_free(m);
}
static double oxygen(const ChemState *s,const size_t *indices)
{ return s->values[indices[6]]+s->values[indices[5]]+2*s->values[indices[7]]+3*s->values[indices[8]]; }
static void integration(void)
{
    ChemMechanism *m=mechanism(1);ChemState *s=NULL;ChemSolver *solver=NULL;
    CHECK(chem_mechanism_reorder(m)==CHEM_OK);
    CHECK(chem_state_create(m,1,&s)==CHEM_OK);CHECK(chem_solver_create(m,1,&solver)==CHEM_OK);
    size_t indices[9];for (size_t i=0;i<9;++i) CHECK(chem_species_index(m,names[i],indices+i)==CHEM_OK);
    const double initial[]={.1,.1,.1,.2,.2,.2,.3,.3,.3};memcpy(s->values,initial,sizeof(initial));
    s->custom[0]=.1;s->custom[1]=.2;s->custom[2]=.3;s->conditions[0]=(ChemConditions){2,3,0};
    double atoms=oxygen(s,indices),o_initial=s->values[indices[6]],o2_initial=s->values[indices[7]];
    ChemRosenbrockParameters p=chem_rosenbrock_parameters(CHEM_ROS3);ChemSolverStats stats;
    for (size_t step=0;step<100;++step) {
        chem_update_rates(s);ChemStatus status=chem_solve(solver,s,30,&p,&stats);CHECK(status==CHEM_OK);
        if (status) break;
        for (size_t i=0;i<9;++i) { CHECK(isfinite(s->values[i]));CHECK(s->values[i]>=0); }
        REL(oxygen(s,indices),atoms,1e-8);
        for (size_t i=0;i<5;++i) REL(s->values[indices[i]],initial[indices[i]],1e-8);
    }
    CHECK(s->values[indices[7]]<.5*o2_initial);CHECK(s->values[indices[6]]>o_initial);
    CHECK(s->values[indices[6]]>.5*atoms);
    chem_solver_free(solver);chem_state_free(s);chem_mechanism_free(m);
}
static double lognormal_draw(uint64_t *seed,double sigma)
{
    *seed=*seed*48271%2147483647;double u=(double)*seed/2147483647;
    *seed=*seed*48271%2147483647;double v=(double)*seed/2147483647;
    return exp(-2+sigma*sqrt(-2*log(u))*cos(2*CHEM_PI*v));
}
static void independent_rates(double temperature,const double *photo,double *k)
{
    k[0]=2.15e-11*exp(110/temperature)*pow(temperature/300,0);
    k[1]=3.3e-11*exp(55/temperature)*pow(temperature/300,0);
    k[2]=8e-12*exp(-2060/temperature)*pow(temperature/300,0);
    k[3]=6e-34*exp(0/temperature)*pow(temperature/300,2.4);
    for (size_t i=0;i<3;++i) k[4+i]=photo[i];
}
static void regression(int which,size_t group,ChemRosenbrockMethod method)
{
    ChemMechanism *m=mechanism(0);CHECK(chem_mechanism_reorder(m)==CHEM_OK);
    size_t indices[9];for (size_t i=0;i<9;++i) CHECK(chem_species_index(m,names[i],indices+i)==CHEM_OK);
    double initial[27],rates[21],abs_tol=0;uint64_t seed=which==3?1:739;
    if (which==1 || which==2) {
        const char *env=getenv("MICM_TEST_SEED");
        if (env) {
            char *end;errno=0;unsigned long value=strtoul(env,&end,10);
            int valid=end!=env && !*end && *env!='-' && errno!=ERANGE && value<=UINT32_MAX;
            CHECK(valid);if (!valid) { chem_mechanism_free(m);return; }
            seed=value%2147483647;if (!seed) seed=1;
        }
    }
    for (size_t cell=0;cell<3;++cell) for (size_t v=0;v<9;++v) {
        double draw=lognormal_draw(&seed,which==2?4:2);initial[v*3+cell]=draw;abs_tol=fmax(abs_tol,draw);
    }
    if (which==2) for (size_t cell=0;cell<3;++cell) for (size_t v=0;v<9;++v)
        initial[v*3+cell]=lognormal_draw(&seed,4); /* the source replaces this array twice */
    for (size_t i=0;i<21;++i) rates[i]=lognormal_draw(&seed,which==2?4:2);
    abs_tol*=1e-12;
    const double temperature[]={284.19,215.02,299.31},pressure[]={101245,100789.2,101398};
    const double photo[3][3]={{1e-4,1e-5,1e-6},{3.2e-4,7.3e-5,3.2e-6},{5.2e-4,8.2e-5,4.6e-6}};
    size_t width=group?group:3;
    for (size_t offset=0;offset<3;offset+=width) {
        size_t cells=3-offset<width?3-offset:width;ChemState *s=NULL;
        CHECK(chem_state_create(m,cells,&s)==CHEM_OK);
        for (size_t cell=0;cell<cells;++cell) {
            s->conditions[cell]=(ChemConditions){temperature[offset+cell],pressure[offset+cell],0};
            for (size_t v=0;v<9;++v) s->values[v*cells+cell]=initial[v*3+offset+cell];
            for (size_t r=0;r<3;++r) s->custom[r*cells+cell]=photo[offset+cell][r];
            for (size_t r=0;r<7;++r) s->rates[r*cells+cell]=rates[r*3+offset+cell];
        }
        if (which==0 || which==3) chem_update_rates(s);
        double computed[243]={0};
        if (which==1) chem_add_kinetic_forcing(s,s->values,computed);
        if (which==2) chem_subtract_kinetic_jacobian(s,s->values,computed);
        if (which==3) {
            ChemSolver *solver=NULL;CHECK(chem_solver_create(m,cells,&solver)==CHEM_OK);
            ChemRosenbrockParameters p=chem_rosenbrock_parameters(method);ChemSolverStats stats;
            CHECK(chem_solve(solver,s,500,&p,&stats)==CHEM_OK);chem_solver_free(solver);
        }
        for (size_t cell=0;cell<cells;++cell) {
            double k[7],y[9],expected[23];
            if (which==0 || which==3) independent_rates(temperature[offset+cell],photo[offset+cell],k);
            else for (size_t r=0;r<7;++r) k[r]=s->rates[r*cells+cell];
            for (size_t v=0;v<9;++v) y[v]=initial[indices[v]*3+offset+cell];
            if (which==0) {
                CHECK(chem_reaction_count(m)==7);
                for (size_t r=0;r<7;++r) CHECK(s->rates[r*cells+cell]==k[r]);
            } else if (which==1) {
                chapman_force(k,y,expected);CHECK(chem_species_count(m)==9);
                for (size_t v=0;v<9;++v) {
                    double a=computed[indices[v]*cells+cell],b=expected[v];NEAR(a,b,(fabs(a)+fabs(b))*1e-8+1e-12);
                }
            } else if (which==2) {
                chapman_jacobian(k,y,expected);
                for (size_t v=0;v<23;++v) {
                    double a=computed[(indices[row[v]]*9+indices[col[v]])*cells+cell],b=expected[v];
                    NEAR(a,-b,(fabs(a)+fabs(b))*1e-8+1e-12);
                }
            } else {
                CHECK(chapman_solve(k,y,500));double relative=method==CHEM_ROS2?1e-2:1e-4;
                for (size_t v=0;v<9;++v) {
                    double a=s->values[indices[v]*cells+cell],b=y[v];NEAR(a,b,(fabs(a)+fabs(b))*relative+abs_tol);
                }
            }
        }
        chem_state_free(s);
    }
    chem_mechanism_free(m);
}
static void regression_vector(int which)
{ for (size_t group=1;group<=4;++group) regression(which,group,CHEM_ROS3); }
void test_chapman(void)
{
    RUN("regression_test_rosenbrock_p_force:RegressionRosenbrock.RateConstants",regression(0,0,CHEM_ROS3));
    RUN("regression_test_rosenbrock_p_force:RegressionRosenbrock.VectorRateConstants",regression_vector(0));
    RUN("regression_test_rosenbrock_p_force:RegressionRosenbrock.Forcing",regression(1,0,CHEM_ROS3));
    RUN("regression_test_rosenbrock_p_force:RegressionRosenbrock.VectorForcing",regression_vector(1));
    RUN("regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.Jacobian",regression(2,0,CHEM_ROS3));
    RUN("regression_test_rosenbrock_dforce_dy:RegressionRosenbrock.VectorJacobian",regression_vector(2));
    const char *methods[]={"TwoStageSolve","ThreeStageSolve","FourStageSolve","FourStageDASolve","SixStageDASolve"};
    for (int i=0;i<5;++i) {
        char id[120];snprintf(id,sizeof(id),"regression_test_solve:RegressionRosenbrock.%s",methods[i]);
        RUN(id,regression(3,0,(ChemRosenbrockMethod)i));
    }
    RUN("regression_test_solve:RegressionRosenbrock.VectorSolve",regression_vector(3));
#define C(name,call) RUN("chapman_ode_solver:ChapmanMechanismHardCodedAndGeneral." name,call)
    C("DefaultConstructor",metadata(0));C("ReactionNames",metadata(1));C("PhotolysisNames",metadata(2));C("SpeciesNames",metadata(3));
    C("lin_solve",linear_algebra(0));C("simple_force",force(0));C("smaller_force",force(1));
    C("factored_alpha_minus_jac",linear_algebra(1));C("dforce_dy_times_vector",linear_algebra(2));
    C("Solve",solve_unit(0));C("solve_10_times_larger",solve_unit(1));C("solve_with_random_number_densities",solve_unit(3));
#undef C
    RUN("chapman_ode_solver:RegressionChapmanODESolver.solve_10_times_smaller",solve_unit(2));
    RUN("chapman_integration:ChapmanIntegration.CanBuildChapmanSystem",integration());
}
