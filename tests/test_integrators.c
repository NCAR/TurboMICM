#include "test.h"
static const ck_species five_species[]={{.name="foo"},{.name="bar"},{.name="baz"},{.name="quz"},{.name="quuz"}};
static const size_t react1[]={0,2},react2[]={1},react3[]={3};
static const ck_product prod1[]={{1,1},{4,2.4}},prod2[]={{0,1},{3,1.4}};
static int model_five(ck_model **model)
{
  static ck_reaction reactions[3];
  reactions[0]=(ck_reaction){.nreactants=2,.reactants=react1,.nproducts=2,.products=prod1,.rate=ck_rate_defaults(CK_ARRHENIUS)};
  reactions[0].rate.A=2e-11;reactions[0].rate.C=110;
  reactions[1]=(ck_reaction){.nreactants=1,.reactants=react2,.nproducts=2,.products=prod2,.rate=ck_rate_defaults(CK_ARRHENIUS)};reactions[1].rate.A=1e-6;
  reactions[2]=(ck_reaction){.nreactants=1,.reactants=react3,.rate=ck_rate_defaults(CK_ARRHENIUS)};reactions[2].rate.A=3.5e-6;
  static ck_mechanism mechanism;
  mechanism=(ck_mechanism){.name="five",.nspecies=5,.species=five_species,.nreactions=3,.reactions=reactions,.temperature=298.15};
  return ck_model_create_with_order(&mechanism,0,model);
}
static int normalized(size_t cells,int constraints)
{
  ck_model *model;ck_state *s;
  if(!constraints)CHECK(model_five(&model)==CK_OK);
  else {
    static const ck_species species[]={{.name="A"},{.name="B"},{.name="C"}};
    static const ck_product reactants[]={{1,1}},products[]={{2,1}};
    static const ck_constraint constraint={.type=CK_EQUILIBRIUM_CONSTRAINT,.name="B_C_eq",.algebraic_species=2,
      .nreactants=1,.reactants=reactants,.nproducts=1,.products=products,.equilibrium_ref=10,.enthalpy=-2400};
    static const ck_mechanism mechanism={.name="normalized",.nspecies=3,.species=species,.nconstraints=1,.constraints=&constraint};
    CHECK(ck_model_create_with_order(&mechanism,0,&model)==CK_OK);CHECK(ck_model_nconstraints(model)==1);
  }
  CHECK(ck_state_create(model,cells,&s)==CK_OK);
  size_t n=ck_model_nspecies(model);double *next=calloc(n*cells,sizeof(double)),*errors=calloc(n*cells,sizeof(double));CHECK(next&&errors);
  double expected=0;
  for(size_t c=0;c<cells;++c)for(size_t j=0;j<n;++j){
    size_t i=j*cells+c;
    s->y[i]=constraints?1.0+c+0.1*j:7.7*c+j;
    next[i]=constraints?0.8+0.5*c+0.2*j:-13.9/(j+1)-c;
    errors[i]=constraints?0.01*(1+c+j):81.57/(c+7)/(j+3);
    double ymax=fmax(fabs(s->y[i]),fabs(next[i]));double scale=s->absolute_tolerance[j]+s->relative_tolerance*ymax;
    expected+=errors[i]*errors[i]/(scale*scale);
  }
  expected=fmax(sqrt(expected/(cells*n)),1e-10);
  NEAR(ck_normalized_error(s,s->y,next,errors),expected,constraints?1e-12:1e-11*expected);
  free(next);free(errors);ck_state_free(s);ck_model_free(model);return 0;
}
int test_integrator(const char *suite,const char *test)
{
  if(!strcmp(suite,"RosenbrockSolver")) {
    if(strstr(test,"NormalizedError")) {
      int constraints=strstr(test,"WithConstraints")!=NULL;
      size_t counts[]={1,2,3,4,5};
      for(size_t j=0;j<5;++j)CHECK(normalized(counts[j],constraints)==0);
      return 0;
    }
    if(IS("CanSetTolerances")) {
      const ck_species species[]={{.name="foo",.absolute_tolerance=1e-7},{.name="bar",.absolute_tolerance=1e-8}};
      const ck_mechanism mechanism={.name="tolerances",.nspecies=2,.species=species};ck_model *m;
      CHECK(ck_model_create_with_order(&mechanism,0,&m)==CK_OK);
      for(size_t c=1;c<=10;++c){ck_state *s;CHECK(ck_state_create(m,c,&s)==CK_OK);NEAR(s->absolute_tolerance[0],1e-7,0);NEAR(s->absolute_tolerance[1],1e-8,0);ck_state_free(s);}ck_model_free(m);return 0;
    }
    if(IS("StandardAlphaMinusJacobian")||IS("VectorAlphaMinusJacobian")) {
      ck_model *m;CHECK(model_five(&m)==CK_OK);
      const size_t rows[]={0,0,0,1,1,1,2,2,3,3,4,4,4},cols[]={0,1,2,0,1,2,0,2,1,3,0,2,4};
      const double values[]={-12.2,-24.3,-42.3,-0.43,-23.4,-83.4,-4.74,-6.91,-59.1,-83.4,-78.5,-53.6,-1};
      for(size_t cells=1;cells<=4;++cells){
        double *a=calloc(ck_jacobian_size(m)*cells,sizeof(double));CHECK(a);
        for(size_t c=0;c<cells;++c)for(size_t j=0;j<13;++j){size_t off=ck_jacobian_index(m,rows[j],cols[j]);CHECK(off!=SIZE_MAX);a[off*cells+c]=values[j]*(j==1?c+2:1)/(j==5?c+3:1);}
        ck_alpha_minus_jacobian(m,cells,a,42.042);
        for(size_t c=0;c<cells;++c)for(size_t j=0;j<13;++j){size_t off=ck_jacobian_index(m,rows[j],cols[j]);double expected=values[j]*(j==1?c+2:1)/(j==5?c+3:1)+(rows[j]==cols[j]?42.042:0);NEAR(a[off*cells+c],expected,0);}
        free(a);
      }ck_model_free(m);return 0;
    }
    if(IS("RejectedStepAlphaMatchesInPlaceSolver")) {
      ck_model *m;CHECK(model_five(&m)==CK_OK);ck_state *a,*b;CHECK(ck_state_create(m,1,&a)==CK_OK);
      for(size_t j=0;j<5;++j)a->y[j]=1e12;
      CHECK(ck_update_rates(a)==CK_OK);CHECK(ck_state_clone(a,&b)==CK_OK);
      ck_rosenbrock_parameters p;CHECK(ck_rosenbrock_parameters_init(3,&p)==CK_OK);p.h_start=1e3;
      ck_stats sa,sb;CHECK(ck_solve(a,1e3,&p,&sa)==CK_OK);p.lu=CK_LU_DOOLITTLE;CHECK(ck_solve(b,1e3,&p,&sb)==CK_OK);
      CHECK(sa.internal_steps-sa.accepted>=2);CHECK(sa.internal_steps==sb.internal_steps&&sa.accepted==sb.accepted&&sa.rejected==sb.rejected);
      ck_state_free(a);ck_state_free(b);ck_model_free(m);return 0;
    }
    return 2;
  }
  if(!strcmp(suite,"BackwardEuler")) {
    const ck_species species[]={{.name="A"},{.name="B"},{.name="C"}};
    const size_t r1[]={0},r2[]={1};const ck_product p1[]={{1,1}},p2[]={{2,1}};
    ck_reaction reactions[]={{.nreactants=1,.reactants=r1,.nproducts=1,.products=p1,.rate=ck_rate_defaults(CK_ARRHENIUS)},
      {.nreactants=1,.reactants=r2,.nproducts=1,.products=p2,.rate=ck_rate_defaults(CK_ARRHENIUS)}};
    reactions[0].rate.A=2.15e-11;reactions[0].rate.C=110;reactions[1].rate.A=3.3e-11;reactions[1].rate.C=55;
    ck_mechanism mechanism={.name="backward",.nspecies=3,.species=species,.nreactions=2,.reactions=reactions};
    ck_model *m;ck_state *s;CHECK(ck_model_create_with_order(&mechanism,0,&m)==CK_OK);
    CHECK(ck_state_create(m,IS("IsConverged")?4:1,&s)==CK_OK);
    double atol[]={1e-6,1e-6,1e-6};CHECK(ck_state_tolerances(s,IS("IsConverged")?1e-3:1e-6,atol,3)==CK_OK);
    if(IS("CanCallSolve")){
      s->y[0]=1;s->conditions[0]=(ck_conditions){272.5,101253.3,1e6};CHECK(ck_update_rates(s)==CK_OK);
      ck_backward_euler_parameters p=ck_backward_euler_defaults();ck_stats stats;CHECK(ck_backward_euler_solve(s,1,&p,&stats)==CK_OK);
    } else if(IS("IsConverged")){
      double residual[12]={0},small=1e-6;
      CHECK(ck_backward_euler_converged(s,residual,s->y,small));residual[4]=1e-5;CHECK(!ck_backward_euler_converged(s,residual,s->y,small));
      small=1e-4;CHECK(ck_backward_euler_converged(s,residual,s->y,small));residual[11]=1e-3;CHECK(!ck_backward_euler_converged(s,residual,s->y,small));
      s->y[11]=10;CHECK(ck_backward_euler_converged(s,residual,s->y,small));residual[11]=1e-1;CHECK(!ck_backward_euler_converged(s,residual,s->y,small));
      s->absolute_tolerance[2]=1;CHECK(ck_backward_euler_converged(s,residual,s->y,small));
    }else{ck_state_free(s);ck_model_free(m);return 2;}
    ck_state_free(s);ck_model_free(m);return 0;
  }
  return 2;
}
