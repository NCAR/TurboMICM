#include "test.h"
#include "chem_factor.h"
typedef struct {
    chem_species species[5]; size_t reactants[3][2];chem_product products[2][2];chem_reaction reactions[3];
    chem_mechanism mechanism;chem_model model;chem_state state;
} core_fixture;
static chem_status core_create(core_fixture *f,size_t cells)
{
    memset(f,0,sizeof(*f));const char *names[]={"foo","bar","baz","quz","quuz"};
    for(size_t i=0;i<5;++i)f->species[i]=(chem_species){.name=names[i]};
    f->reactants[0][0]=0;f->reactants[0][1]=2;f->reactants[1][0]=1;f->reactants[2][0]=3;
    f->products[0][0]=(chem_product){1,1};f->products[0][1]=(chem_product){4,2.4};
    f->products[1][0]=(chem_product){0,1};f->products[1][1]=(chem_product){3,1.4};
    for(size_t i=0;i<3;++i){
        f->reactions[i]=(chem_reaction){i ? 1 : 2,i<2 ? 2 : 0,f->reactants[i],i<2 ? f->products[i] : NULL,chem_rate_default(CHEM_ARRHENIUS)};
        f->reactions[i].rate.p.arrhenius.A=(double[]){2e-11,1e-6,3.5e-6}[i];f->reactions[i].rate.p.arrhenius.C=i ? 0 : 110;
    }
    f->mechanism=(chem_mechanism){.n_species=5,.species=f->species,.n_reactions=3,.reactions=f->reactions,.preserve_order=1};
    chem_status status=chem_model_create(&f->model,&f->mechanism,CHEM_LU_MOZART);if(status)return status;
    return chem_state_create(&f->state,&f->model,cells);
}
static void core_destroy(core_fixture *f){chem_state_destroy(&f->state);chem_model_destroy(&f->model);}
static void alpha_case(size_t cells)
{
    core_fixture f;REQUIRE(core_create(&f,cells)==CHEM_OK);CHECK(f.model.n_species==5);CHECK(f.state.cells==cells);CHECK(f.model.lu.nnz*cells>=13*cells);
    const size_t rows[]={0,0,0,1,1,1,2,2,3,3,4,4,4},cols[]={0,1,2,0,1,2,0,2,1,3,0,2,4};
    const double x[]={-12.2,-24.3,-42.3,-.43,-23.4,-83.4,-4.74,-6.91,-59.1,-83.4,-78.5,-53.6,-1};
    for(size_t c=0;c<cells;++c)for(size_t j=0;j<13;++j){size_t k=f.model.lu.map[rows[j]*5+cols[j]];REQUIRE(k!=CHEM_NONE);f.state.jacobian[k*cells+c]=j==1 ? x[j]*(c+2) : j==5 ? x[j]/(c+3) : x[j];}
    chem_alpha_minus_jacobian(&f.state,42.042);
    for(size_t c=0;c<cells;++c)for(size_t j=0;j<13;++j){double expected=j==1 ? x[j]*(c+2) : j==5 ? x[j]/(c+3) : x[j];if(rows[j]==cols[j])expected+=42.042;NEAR(f.state.jacobian[f.model.lu.map[rows[j]*5+cols[j]]*cells+c],expected,0);}
    core_destroy(&f);
}
TEST(rosenbrock,RosenbrockSolver,StandardAlphaMinusJacobian){for(size_t c=1;c<=4;++c)alpha_case(c);}
TEST(rosenbrock,RosenbrockSolver,VectorAlphaMinusJacobian){alpha_case(1);alpha_case(4);alpha_case(3);alpha_case(2);}
static void normalized(size_t cells,int constrained,chem_ros_method method)
{
    core_fixture f;REQUIRE(core_create(&f,cells)==CHEM_OK);
    chem_species species[]={{.name="A"},{.name="B"},{.name="C"}};size_t reactant=0;chem_product product={1,1},b={1,1},c={2,1};
    chem_reaction reaction={1,1,&reactant,&product,chem_rate_default(CHEM_ARRHENIUS)};reaction.rate.p.arrhenius.A=.5;
    chem_constraint constraint={.name="B_C_eq",.type=CHEM_EQUILIBRIUM_CONSTRAINT,.algebraic_species=2,.n_reactants=1,.n_products=1,.reactants=&b,.products=&c,.equilibrium={.K=10,.delta_H=-2400}};
    chem_mechanism mechanism={.n_species=3,.species=species,.n_reactions=1,.reactions=&reaction,.n_constraints=1,.constraints=&constraint,.preserve_order=1};
    if(constrained){core_destroy(&f);REQUIRE(chem_model_create(&f.model,&mechanism,CHEM_LU_MOZART)==CHEM_OK);REQUIRE(chem_state_create(&f.state,&f.model,cells)==CHEM_OK);CHECK(f.model.mechanism->n_constraints==1);}
    chem_ros_parameters p=chem_ros_parameters_default(method);CHECK(p.stages>=2 && p.stages<=6);
    size_t n=f.model.n_species;double *y=f.state.y,*yn=f.state.work,*err=yn+n*cells,expected=0;
    for(size_t i=0;i<cells;++i)for(size_t j=0;j<n;++j){
        size_t k=j*cells+i;y[k]=constrained ? 1+i+.1*j : 7.7*i+j;yn[k]=constrained ? .8+.5*i+.2*j : -13.9/(j+1)-i;err[k]=constrained ? .01*(1+i+j) : 81.57/(i+7)/(j+3);
        double scale=f.state.atol[j]+f.state.rtol*fmax(fabs(y[k]),fabs(yn[k]));expected+=err[k]*err[k]/(scale*scale);
    }
    expected=fmax(sqrt(expected/(cells*n)),1e-10);double actual=chem_normalized_error(&f.state,y,yn,err);
    if(constrained)NEAR(actual,expected,1e-12);else CHECK(fabs(actual-expected)/fmax(fabs(actual),fabs(expected))<=1e-11);
    core_destroy(&f);
}
TEST(rosenbrock,RosenbrockSolver,StandardNormalizedError){for(size_t c=1;c<=4;++c)normalized(c,0,CHEM_ROS3);}
TEST(rosenbrock,RosenbrockSolver,VectorNormalizedError){size_t cells[]={1,2,3,4,1,2,3,5,3};for(size_t i=0;i<9;++i)normalized(cells[i],0,CHEM_ROS3);}
TEST(rosenbrock,RosenbrockSolver,StandardNormalizedErrorWithConstraints){for(int m=0;m<5;++m)normalized(2,1,(chem_ros_method)m);normalized(5,1,CHEM_ROS3);}
TEST(rosenbrock,RosenbrockSolver,VectorNormalizedErrorWithConstraints){for(int m=0;m<5;++m)normalized(3,1,(chem_ros_method)m);normalized(5,1,CHEM_ROS3);}
TEST(rosenbrock,RosenbrockSolver,CanSetTolerances)
{
    chem_species species[]={{.name="foo",.absolute_tolerance=1e-7},{.name="bar",.absolute_tolerance=1e-8}};
    size_t reactant=0;chem_product product={1,1};chem_reaction reaction={1,1,&reactant,&product,chem_rate_default(CHEM_ARRHENIUS)};reaction.rate.p.arrhenius.A=2e-11;reaction.rate.p.arrhenius.C=110;
    chem_mechanism mechanism={.n_species=2,.species=species,.n_reactions=1,.reactions=&reaction};chem_model model;REQUIRE(chem_model_create(&model,&mechanism,CHEM_LU_MOZART)==CHEM_OK);
    for(size_t cells=1;cells<=10;++cells){chem_state s;REQUIRE(chem_state_create(&s,&model,cells)==CHEM_OK);CHECK(model.n_species==2);NEAR(s.atol[0],1e-7,0);NEAR(s.atol[1],1e-8,0);chem_state_destroy(&s);}
    chem_model_destroy(&model);
}
/* Alternate backend owns separate sparse L and U, demonstrating replacement
 * without changing Rosenbrock, backward Euler, or the model's symbolic plan. */
typedef struct {chem_sparse a,l,u;chem_array rhs;size_t factors,solves;} split_backend;
static chem_status split_create(split_backend *b,const chem_state *s)
{
    memset(b,0,sizeof(*b));chem_sparse_builder builder;chem_status status=chem_sparse_builder_create(&builder,s->model->n_species);if(status)return status;
    builder.blocks=s->cells;
    for(size_t i=0;i<builder.n;++i)for(size_t j=0;j<builder.n;++j)if(s->model->lu.map[i*builder.n+j]!=CHEM_NONE)chem_sparse_builder_add(&builder,i,j);
    status=chem_sparse_create(&b->a,&builder,CHEM_CSR,3,CHEM_DOUBLE);chem_sparse_builder_destroy(&builder);if(status)return status;
    status=chem_sparse_lu_matrices(&s->model->lu,&b->a,&b->l,&b->u,0);if(status)return status;
    return chem_array_create(&b->rhs,s->cells,s->model->n_species,3,CHEM_DOUBLE,0);
}
static void split_destroy(split_backend *b){chem_sparse_destroy(&b->a);chem_sparse_destroy(&b->l);chem_sparse_destroy(&b->u);chem_array_destroy(&b->rhs);}
static chem_status split_factor(chem_state *s,void *context)
{
    split_backend *b=context;++b->factors;size_t n=s->model->n_species;
    for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j){size_t k=s->model->lu.map[i*n+j];if(k!=CHEM_NONE)for(size_t c=0;c<s->cells;++c)chem_sparse_set(&b->a,c,i,j,s->jacobian[k*s->cells+c]);}
    return chem_sparse_lu_factor(&s->model->lu,&b->a,&b->l,&b->u);
}
static chem_status split_solve(const chem_state *s,double *rhs,void *context)
{
    split_backend *b=context;++b->solves;size_t n=s->model->n_species;
    for(size_t i=0;i<n;++i)for(size_t c=0;c<s->cells;++c)chem_array_set(&b->rhs,c,i,rhs[i*s->cells+c]);
    chem_status status=chem_sparse_lu_solve(&s->model->lu,&b->l,&b->u,&b->rhs);if(status)return status;
    for(size_t i=0;i<n;++i)for(size_t c=0;c<s->cells;++c)rhs[i*s->cells+c]=chem_array_get(&b->rhs,c,i);
    return CHEM_OK;
}
TEST(rosenbrock,RosenbrockSolver,RejectedStepAlphaMatchesInPlaceSolver)
{
    core_fixture a,b;REQUIRE(core_create(&a,1)==CHEM_OK);REQUIRE(core_create(&b,1)==CHEM_OK);
    split_backend split;REQUIRE(split_create(&split,&b.state)==CHEM_OK);CHECK(chem_state_set_linear_backend(&b.state,split_factor,split_solve,&split)==CHEM_OK);
    for(size_t i=0;i<5;++i){a.state.y[i]=1e12;b.state.y[i]=1e12;}a.state.conditions[0].temperature=298.15;b.state.conditions[0].temperature=298.15;chem_calculate_rates(&a.state);chem_calculate_rates(&b.state);
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);p.h_start=1e3;
    chem_result ra=chem_ros_solve(&a.state,1e3,&p,NULL,NULL),rb=chem_ros_solve(&b.state,1e3,&p,NULL,NULL);
    CHECK(ra.status==CHEM_OK);CHECK(rb.status==CHEM_OK);CHECK(ra.stats.internal_steps-ra.stats.accepted>=2);CHECK(rb.stats.internal_steps-rb.stats.accepted>=2);
    CHECK(ra.stats.internal_steps==rb.stats.internal_steps);CHECK(ra.stats.accepted==rb.stats.accepted);CHECK(ra.stats.rejected==rb.stats.rejected);
    CHECK(split.factors==rb.stats.decompositions);CHECK(split.solves==rb.stats.solves);for(size_t i=0;i<5;++i)REL(a.state.y[i],b.state.y[i],1e-14);
    CHECK(chem_state_set_linear_backend(&b.state,NULL,split_solve,&split)==CHEM_INVALID);CHECK(chem_state_set_linear_backend(&b.state,NULL,NULL,NULL)==CHEM_OK);
    split_destroy(&split);core_destroy(&a);core_destroy(&b);
}
TEST(backward_euler,BackwardEuler,CanCallSolve)
{
    chem_species species[]={{.name="A",.absolute_tolerance=1e-6},{.name="B",.absolute_tolerance=1e-6},{.name="C",.absolute_tolerance=1e-6}};size_t reactants[]={0,1};chem_product products[]={{1,1},{2,1}};chem_reaction reactions[2];
    for(size_t i=0;i<2;++i){reactions[i]=(chem_reaction){1,1,&reactants[i],&products[i],chem_rate_default(CHEM_ARRHENIUS)};reactions[i].rate.p.arrhenius.A=i ? 3.3e-11 : 2.15e-11;reactions[i].rate.p.arrhenius.C=i ? 55 : 110;}
    chem_mechanism mechanism={.n_species=3,.species=species,.n_reactions=2,.reactions=reactions};chem_model m;chem_state s;
    REQUIRE(chem_model_create(&m,&mechanism,CHEM_LU_MOZART)==CHEM_OK);REQUIRE(chem_state_create(&s,&m,1)==CHEM_OK);
    s.y[0]=1;s.conditions[0]=(chem_conditions){272.5,101253.3,1e6};chem_calculate_rates(&s);chem_be_parameters p=chem_be_parameters_default();chem_result result=chem_be_solve(&s,1,&p);CHECK(result.status==CHEM_OK);
    chem_state_destroy(&s);chem_model_destroy(&m);
}
TEST(backward_euler,BackwardEuler,IsConverged)
{
    chem_species species[]={{.name="A",.absolute_tolerance=1e-6},{.name="B",.absolute_tolerance=1e-6},{.name="C",.absolute_tolerance=1e-6}};
    chem_mechanism mechanism={.n_species=3,.species=species,.preserve_order=1};chem_model m;chem_state s;
    REQUIRE(chem_model_create(&m,&mechanism,CHEM_LU_MOZART)==CHEM_OK);REQUIRE(chem_state_create(&s,&m,4)==CHEM_OK);s.rtol=1e-3;
    for(size_t v=0;v<=5;++v){ /* C SoA state replaces each original storage policy. */
        double residual[12]={0},y[12]={0};s.atol[2]=1e-6;chem_be_parameters p=chem_be_parameters_default();p.small=1e-6;
        CHECK(chem_be_converged(&s,&p,residual,y));residual[1*4]=1e-5;CHECK(!chem_be_converged(&s,&p,residual,y));p.small=1e-4;CHECK(chem_be_converged(&s,&p,residual,y));
        residual[2*4+3]=1e-3;CHECK(!chem_be_converged(&s,&p,residual,y));y[2*4+3]=10;CHECK(chem_be_converged(&s,&p,residual,y));residual[2*4+3]=1e-1;CHECK(!chem_be_converged(&s,&p,residual,y));s.atol[2]=1;CHECK(chem_be_converged(&s,&p,residual,y));
    }
    chem_state_destroy(&s);chem_model_destroy(&m);
}
