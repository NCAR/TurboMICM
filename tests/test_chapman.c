#include "test.h"
#include "chapman_fixture.h"
#include "chapman_reference.h"
#include "test_random.h"
#include <errno.h>

TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,DefaultConstructor)
{ CHECK(chem_ros_parameters_default(CHEM_ROS3).stages==3); }
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,ReactionNames)
{ chapman_fixture f; REQUIRE(chapman_create(&f,1,0,0,1,0)==CHEM_OK); CHECK(f.model.n_reactions==7); chapman_destroy(&f); }
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,PhotolysisNames)
{ chapman_fixture f; REQUIRE(chapman_create(&f,1,0,0,1,0)==CHEM_OK); CHECK(f.model.n_custom==3); chapman_destroy(&f); }
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,SpeciesNames)
{ chapman_fixture f; REQUIRE(chapman_create(&f,1,0,0,1,0)==CHEM_OK); CHECK(f.model.n_species==9); chapman_destroy(&f); }
static void chapman_force_simple(int small)
{
    chapman_fixture f; REQUIRE(chapman_create(&f,1,0,0,1,0)==CHEM_OK);
    for (size_t j=0;j<9;++j) f.state.y[j]=small ? 5e-6 : 1;
    for (size_t j=0;j<7;++j) f.state.rates[j]=small ? 3e-8 : 1;
    double out[9]={0},oracle[9]; f.model.forcing(&f.model,1,f.state.rates,f.state.y,out);
    chapman_reference_force(f.state.rates,f.state.y,oracle);
    const double expected[]={0,0,0,0,0,-1,3,2,-2},little[]={0,0,0,0,0,1.49e-13,4.55e-13,1.5e-13,-3e-13};
    for (size_t j=0;j<9;++j) {
        NEAR(out[j],small ? little[j] : expected[j],small && j>=5 ? .01 : 0);
        NEAR(oracle[j],small ? little[j] : expected[j],small && j>=5 ? .01 : 0);
    }
    chapman_destroy(&f);
}
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,simple_force) { chapman_force_simple(0); }
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,smaller_force) { chapman_force_simple(1); }
static chem_status chapman_plan(chem_lu_plan *p)
{
    unsigned char pattern[81]={0}; for (size_t i=0;i<23;++i) pattern[chapman_positions[i][0]*9+chapman_positions[i][1]]=1;
    return chem_lu_create(p,9,pattern,CHEM_LU_MOZART);
}
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,lin_solve)
{
    chem_lu_plan plan={0}; REQUIRE(chapman_plan(&plan)==CHEM_OK); CHECK(plan.nnz==23);
    double jac[23],rhs[9],oracle[9]; for (size_t i=0;i<23;++i) jac[i]=1;
    for (size_t i=0;i<9;++i) rhs[i]=.5;
    chapman_reference_linear(jac,rhs,oracle); chem_lu_solve(&plan,1,jac,rhs);
    const double expected[]={.5,.5,.5,.5,.5,-.5,-1,.5,0};
    for (size_t i=0;i<9;++i) { NEAR(rhs[i],expected[i],0); NEAR(oracle[i],expected[i],0); }
    chem_lu_destroy(&plan);
}
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,factored_alpha_minus_jac)
{
    chem_lu_plan plan={0}; REQUIRE(chapman_plan(&plan)==CHEM_OK); double jac[23],oracle[23],scratch;
    for (size_t i=0;i<23;++i) {
        size_t row=chapman_positions[i][0],col=chapman_positions[i][1];
        oracle[i]=row==col ? 1 : -1; jac[plan.map[row*9+col]]=oracle[i];
    }
    chem_lu_factor(&plan,1,jac,&scratch); chapman_reference_factor(oracle);
    const double expected[]={1,-1,-1,-1,1,1,1,1,-1,-1,1,-1,1,-1,-1,-1,-2,-1,3,-1,-2,-3,.125};
    for (size_t i=0;i<23;++i) {
        size_t row=chapman_positions[i][0],col=chapman_positions[i][1]; double actual=jac[plan.map[row*9+col]];
        if (row==col) actual=1/actual; /* Legacy storage keeps reciprocal U diagonals. */
        NEAR(actual,expected[i],.01); NEAR(oracle[i],expected[i],.01);
    }
    chem_lu_destroy(&plan);
}
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,dforce_dy_times_vector)
{
    for (size_t width=0;width<=4;++width) for (int order=CHEM_CSR;order<=CHEM_CSC;++order) {
        chem_sparse_builder builder={0}; REQUIRE(chem_sparse_builder_create(&builder,9)==CHEM_OK); builder.blocks=3;
        for (size_t j=0;j<23;++j) if (chapman_positions[j][0]>=5)
            REQUIRE(chem_sparse_builder_add(&builder,chapman_positions[j][0],chapman_positions[j][1])==CHEM_OK);
        builder.initial=1; chem_sparse matrix={0}; chem_array x={0},y={0};
        REQUIRE(chem_sparse_create(&matrix,&builder,(chem_sparse_order)order,width,CHEM_DOUBLE)==CHEM_OK);
        REQUIRE(chem_array_create(&x,3,9,width,CHEM_DOUBLE,.5)==CHEM_OK);
        REQUIRE(chem_array_create(&y,3,9,width ? 0 : 4,CHEM_DOUBLE,0)==CHEM_OK);
        REQUIRE(chem_sparse_matvec(&matrix,&x,&y)==CHEM_OK);
        for (size_t cell=0;cell<3;++cell) {
            double result[23]={0},expected[]={0,0,0,0,0,2,3,2,2};
            REQUIRE(chem_array_export_row(&y,cell,result,9)==CHEM_OK);
            for (size_t i=0;i<23;++i) NEAR(result[i],i<9 ? expected[i] : 0,.01);
        }
        CHECK(chem_sparse_matvec(&matrix,&x,&x)==CHEM_INVALID);
        chem_sparse_builder_destroy(&builder); chem_sparse_destroy(&matrix); chem_array_destroy(&x); chem_array_destroy(&y);
    }
}
static void chapman_unit_solve(double scale,int random)
{
    chapman_fixture f; REQUIRE(chapman_create(&f,1,0,0,1,0)==CHEM_OK);
    double initial[]={1,.392,.0169,0,32.9,0,0,8.84,0}; random_state rng; random_init(&rng,1,0);
    for (size_t i=0;i<9;++i) f.state.y[i]=random ? (double)((float)(draw(&rng)-1)/2147483646.0f) : initial[i]*scale;
    f.state.conditions[0]=(chem_conditions){273.15,100000,2.7e19}; memcpy(f.state.custom,(double[]){1e-4,1e-5,1e-6},3*sizeof(double));
    chem_calculate_rates(&f.state); f.state.rtol=1e-4;
    double oracle[9]; memcpy(oracle,f.state.y,sizeof(oracle)); double advanced=chapman_reference_solve(oracle,f.state.rates,1);
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3); p.h_start=1e-6; p.max_steps=100;
    chem_result r=chem_ros_solve(&f.state,1,&p,NULL,NULL); REQUIRE(r.status==CHEM_OK); NEAR(advanced,1,0);
    if (!random) {
        const double expected[]={1,.392,.0169,0,32.9,1.8039e-41,.00176789,8.83912,4.5031e-36};
        for (size_t i=0;i<9;++i) {
            double ref=expected[i]*(i==5 || i==8 ? scale*scale*scale : scale);
            NEAR(f.state.y[i],ref,1e-4); NEAR(oracle[i],ref,1e-4);
        }
    }
    chapman_destroy(&f);
}
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,Solve) { chapman_unit_solve(1,0); }
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,solve_10_times_larger) { chapman_unit_solve(10,0); }
TEST(chapman_ode_solver,RegressionChapmanODESolver,solve_10_times_smaller) { chapman_unit_solve(.1,0); }
TEST(chapman_ode_solver,ChapmanMechanismHardCodedAndGeneral,solve_with_random_number_densities) { chapman_unit_solve(1,1); }

static void chapman_rates_case(size_t width)
{
    for (int generated=0;generated<2;++generated) {
        chapman_fixture f; REQUIRE(chapman_create(&f,3,width,generated,0,0)==CHEM_OK); chapman_conditions(&f.state);
        CHECK(f.model.n_reactions==7);
        for (size_t c=0;c<3;++c) {
            double photo[3],rates[7]; for (size_t j=0;j<3;++j) photo[j]=f.state.custom[j*3+c];
            chapman_reference_rates(f.state.conditions[c].temperature,photo,rates);
            for (size_t j=0;j<7;++j) NEAR(f.state.rates[j*3+c],rates[j],0);
        }
        chapman_destroy(&f);
    }
}
TEST(regression_test_rosenbrock_p_force,RegressionRosenbrock,RateConstants) { chapman_rates_case(0); }
TEST(regression_test_rosenbrock_p_force,RegressionRosenbrock,VectorRateConstants) { for (size_t width=1;width<=4;++width) chapman_rates_case(width); }
static int chapman_seed(uint32_t *seed)
{
    *seed=739; const char *env=getenv("MICM_TEST_SEED"); if (!env) return 1;
    char *end; errno=0; unsigned long value=strtoul(env,&end,10);
    if (end==env || *end || *env=='-' || errno==ERANGE || value>UINT32_MAX) return 0;
    *seed=(uint32_t)value; return 1;
}
static void random_physical(random_state *rng,size_t columns,size_t width,double sigma,int passes,double *soa)
{
    chem_array array={0}; REQUIRE(chem_array_create(&array,3,columns,width,CHEM_DOUBLE,0)==CHEM_OK);
    for (int p=0;p<passes;++p) for (size_t j=0;j<array.size;++j) ((double *)array.data)[j]=lognormal(rng,sigma);
    for (size_t c=0;c<3;++c) for (size_t j=0;j<columns;++j) soa[j*3+c]=chem_array_get(&array,c,j);
    chem_array_destroy(&array);
}
static void chapman_derivative_case(size_t width,int jacobian)
{
    uint32_t seed; REQUIRE(chapman_seed(&seed));
    for (int generated=0;generated<2;++generated) {
        chapman_fixture f; REQUIRE(chapman_create(&f,3,width,generated,0,0)==CHEM_OK);
        random_state rng; random_init(&rng,seed,1); double sigma=jacobian ? 4 : 2;
        /* The source fills the state twice in its Jacobian case, intentionally
         * preserving both draws and all padded lanes in this reproduction. */
        random_physical(&rng,9,width,sigma,jacobian ? 2 : 1,f.state.y);
        random_physical(&rng,7,width,sigma,1,f.state.rates);
        double out[27]={0};
        if (jacobian) f.model.jacobian(&f.model,3,f.state.rates,f.state.y,f.state.jacobian);
        else f.model.forcing(&f.model,3,f.state.rates,f.state.y,out);
        for (size_t cell=0;cell<3;++cell) {
            double rates[7],y[9],ref[23];
            for (size_t i=0;i<7;++i) rates[i]=f.state.rates[i*3+cell];
            for (size_t i=0;i<9;++i) y[i]=f.state.y[f.model.spec_to_var[i]*3+cell];
            if (jacobian) chapman_reference_jacobian(rates,y,ref); else chapman_reference_force(rates,y,ref);
            for (size_t i=0;i<(jacobian ? 23u : 9u);++i) {
                double a,b=ref[i];
                if (jacobian) {
                    size_t row=f.model.spec_to_var[chapman_positions[i][0]],col=f.model.spec_to_var[chapman_positions[i][1]];
                    a=f.state.jacobian[f.model.lu.map[row*9+col]*3+cell]; b=-b;
                } else a=out[f.model.spec_to_var[i]*3+cell];
                NEAR(a,b,(fabs(a)+fabs(b))*1e-8+1e-12);
            }
        }
        chapman_destroy(&f);
    }
}
TEST(regression_test_rosenbrock_p_force,RegressionRosenbrock,Forcing) { chapman_derivative_case(0,0); }
TEST(regression_test_rosenbrock_p_force,RegressionRosenbrock,VectorForcing) { for (size_t width=1;width<=4;++width) chapman_derivative_case(width,0); }
TEST(regression_test_rosenbrock_dforce_dy,RegressionRosenbrock,Jacobian) { chapman_derivative_case(0,1); }
TEST(regression_test_rosenbrock_dforce_dy,RegressionRosenbrock,VectorJacobian) { for (size_t width=1;width<=4;++width) chapman_derivative_case(width,1); }
static void chapman_regression_solve(chem_ros_method method,size_t width,double tolerance)
{
    for (int generated=0;generated<2;++generated) {
        chapman_fixture f; REQUIRE(chapman_create(&f,3,width,generated,0,0)==CHEM_OK); chapman_conditions(&f.state);
        random_state rng; random_init(&rng,1,0); double initial[3][9],atol=0;
        for (size_t c=0;c<3;++c) for (size_t i=0;i<9;++i) {
            initial[c][i]=lognormal(&rng,2); f.state.y[i*3+c]=initial[c][i]; atol=fmax(atol,initial[c][i]);
        }
        atol*=1e-12; chem_ros_parameters p=chem_ros_parameters_default(method);
        chem_result result=chem_ros_solve(&f.state,500,&p,NULL,NULL); chem_state_clamp(&f.state); CHECK(result.status==CHEM_OK);
        for (size_t c=0;c<3;++c) {
            double ref[9],rates[7];
            for (size_t j=0;j<9;++j) ref[j]=initial[c][f.model.spec_to_var[j]];
            for (size_t j=0;j<7;++j) rates[j]=f.state.rates[j*3+c];
            NEAR(chapman_reference_solve(ref,rates,500),500,0);
            for (size_t j=0;j<9;++j) {
                double a=f.state.y[f.model.spec_to_var[j]*3+c],b=ref[j];
                NEAR(a,b,(fabs(a)+fabs(b))*tolerance+atol);
            }
        }
        chapman_destroy(&f);
    }
}
TEST(regression_test_solve,RegressionRosenbrock,TwoStageSolve) { chapman_regression_solve(CHEM_ROS2,0,1e-2); }
TEST(regression_test_solve,RegressionRosenbrock,ThreeStageSolve) { chapman_regression_solve(CHEM_ROS3,0,1e-4); }
TEST(regression_test_solve,RegressionRosenbrock,FourStageSolve) { chapman_regression_solve(CHEM_ROS4,0,1e-4); }
TEST(regression_test_solve,RegressionRosenbrock,FourStageDASolve) { chapman_regression_solve(CHEM_RODA3,0,1e-4); }
TEST(regression_test_solve,RegressionRosenbrock,SixStageDASolve) { chapman_regression_solve(CHEM_RODA4,0,1e-4); }
TEST(regression_test_solve,RegressionRosenbrock,VectorSolve) { for (size_t width=1;width<=4;++width) chapman_regression_solve(CHEM_ROS3,width,1e-4); }
TEST(chapman_integration,ChapmanIntegration,CanBuildChapmanSystem)
{
    /* This integration fixture declares a different species order and gives
     * reaction four C=2.4, whereas the regression fixture uses B=2.4. */
    chapman_fixture f; REQUIRE(chapman_create(&f,1,0,0,0,1)==CHEM_OK); chapman_destroy(&f);
    const size_t old_to_new[]={4,5,8,7,6,1,0,2,3};
    chem_species reordered[9]; for (size_t i=0;i<9;++i) reordered[old_to_new[i]]=f.species[i];
    memcpy(f.species,reordered,sizeof(reordered));
    for (size_t r=0;r<7;++r) {
        for (size_t j=0;j<f.reactions[r].n_reactants;++j) f.reactants[r][j]=old_to_new[f.reactants[r][j]];
        for (size_t j=0;j<f.reactions[r].n_products;++j) f.products[r][j].species=old_to_new[f.products[r][j].species];
    }
    REQUIRE(chem_model_create(&f.model,&f.mechanism,CHEM_LU_DOOLITTLE)==CHEM_OK);
    REQUIRE(chem_state_create(&f.state,&f.model,1)==CHEM_OK); REQUIRE(integration_backend_create(&f.backend,&f.state,0,CHEM_CSR,0)==CHEM_OK);
    memcpy(f.state.y,(double[]){.1,.1,.1,.2,.2,.2,.3,.3,.3},9*sizeof(double));
    memcpy(f.state.custom,(double[]){.1,.2,.3},3*sizeof(double)); f.state.conditions[0]=(chem_conditions){2,3,0};
    size_t ids[9]; for (size_t i=0;i<9;++i) ids[i]=chem_species_index(&f.model,chapman_names[i]);
    double initial[9]; for (size_t i=0;i<9;++i) initial[i]=f.state.y[ids[i]];
    double oxygen=initial[6]+initial[5]+2*initial[7]+3*initial[8]; chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);
    for (size_t t=0;t<100;++t) {
        chem_calculate_rates(&f.state); chem_result r=chem_ros_solve(&f.state,30,&p,NULL,NULL); chem_state_clamp(&f.state); REQUIRE(r.status==CHEM_OK);
        for (size_t i=0;i<9;++i) { CHECK(isfinite(f.state.y[i])); CHECK(f.state.y[i]>=0); }
        REL(f.state.y[ids[6]]+f.state.y[ids[5]]+2*f.state.y[ids[7]]+3*f.state.y[ids[8]],oxygen,1e-8);
        for (size_t i=0;i<5;++i) REL(f.state.y[ids[i]],initial[i],1e-8);
    }
    CHECK(f.state.y[ids[7]]<.5*initial[7]); CHECK(f.state.y[ids[6]]>initial[6]); CHECK(f.state.y[ids[6]]>.5*oxygen);
    chapman_destroy(&f);
}
