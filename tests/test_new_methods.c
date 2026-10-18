#include "chem.h"
#include "chem_integrator.h"
#include <math.h>
#include <stdio.h>
static int check(int condition,const char *message)
{
  if (!condition) fprintf(stderr,"FAIL: %s\n",message);
  return condition;
}
int main(void)
{
  static const ck_species sp[]={{.name="A",.initial=1},{.name="B"}};
  static const size_t a[]={0};
  static const ck_product b[]={{1,1}};
  ck_reaction rx[2]={{.nreactants=1,.reactants=a,.nproducts=1,.products=b,.rate=ck_rate_defaults(CK_ARRHENIUS)}};
  rx[0].rate.A=3;
  ck_mechanism mech={.name="linear-transfer",.nspecies=2,.species=sp,.nreactions=1,.reactions=rx,.temperature=298};
  for (int nonlinear=0;nonlinear<2;++nonlinear)
  {
  double exact=exp(-3), mass=1;
  if (nonlinear)
  {
    static const size_t aa[]={0,0}, bb[]={1};
    static const ck_product pa[]={{0,2}};
    rx[0].nreactants=2;rx[0].reactants=aa;rx[0].rate.A=2;
    rx[1]=(ck_reaction){.nreactants=1,.reactants=bb,.nproducts=1,.products=pa,
      .rate=ck_rate_defaults(CK_ARRHENIUS)};
    rx[1].rate.A=.5;mech.nreactions=2;mass=1.2;
    double d=sqrt(.25+8*2*.5*mass), rp=(-.5+d)/8, rm=(-.5-d)/8;
    double c=(1-rp)/(1-rm)*exp(-d);
    exact=(rp-c*rm)/(1-c);
  }
  ck_model *m=NULL;
  if (ck_model_create_with_order(&mech,0,&m)) return 1;
  const int methods[]={2,3,4,43,6,12,13}, orders[]={2,3,4,3,4,2,3};
  for (size_t k=0;k<7;++k)
  {
    double error[2];
    for (size_t trial=0;trial<2;++trial)
    {
      ck_state *s=NULL;ck_integrator *v=NULL;
      if (ck_state_create(m,1,&s)) return 1;
      s->relative_tolerance=0;
      if (nonlinear) s->y[1]=.1;
      s->absolute_tolerance[0]=s->absolute_tolerance[1]=1e6;
      ck_integrator_options o=ck_integrator_defaults();
      o.method=methods[k];o.growth=1;o.start=trial ? .0125 : .025;o.positivity=0;
      if (ck_update_rates(s) || ck_integrator_create(s,&o,&v)) return 1;
      ck_stats stats;
      if (ck_integrator_advance(v,1,&stats)) return 1;
      error[trial]=fabs(s->y[0]-exact);
      if (!check(fabs(s->y[0]+(nonlinear ? 2 : 1)*s->y[1]-mass)<1e-13,"stoichiometric invariant") ||
          !check(stats.rejected==0,"fixed-step run accepted")) return 1;
      ck_integrator_free(v);ck_state_free(s);
    }
    double observed=log(error[0]/error[1])/log(2);
    printf("nonlinear=%d method=%d observed_order=%.4f expected_order=%d\n",nonlinear,methods[k],observed,orders[k]);
    if (!check(fabs(observed-orders[k])<.35,"analytic convergence order")) return 1;
  }
  ck_model_free(m);
  }
  puts("Analytic order and invariant checks PASS");
  return 0;
}
