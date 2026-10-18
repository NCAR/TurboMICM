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
  ck_reaction rx={.nreactants=1,.reactants=a,.nproducts=1,.products=b,.rate=ck_rate_defaults(CK_ARRHENIUS)};
  rx.rate.A=3;
  ck_mechanism mech={.name="linear-transfer",.nspecies=2,.species=sp,.nreactions=1,.reactions=&rx,.temperature=298};
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
      s->absolute_tolerance[0]=s->absolute_tolerance[1]=1e6;
      ck_integrator_options o=ck_integrator_defaults();
      o.method=methods[k];o.growth=1;o.start=trial ? .025 : .05;o.positivity=0;
      if (ck_update_rates(s) || ck_integrator_create(s,&o,&v)) return 1;
      ck_stats stats;
      if (ck_integrator_advance(v,1,&stats)) return 1;
      error[trial]=fabs(s->y[0]-exp(-3));
      if (!check(fabs(s->y[0]+s->y[1]-1)<1e-13,"linear invariant") ||
          !check(stats.rejected==0,"fixed-step run accepted")) return 1;
      ck_integrator_free(v);ck_state_free(s);
    }
    double observed=log(error[0]/error[1])/log(2);
    printf("method=%d observed_order=%.4f expected_order=%d\n",methods[k],observed,orders[k]);
    if (!check(fabs(observed-orders[k])<.35,"analytic convergence order")) return 1;
  }
  ck_model_free(m);
  puts("Analytic order and invariant checks PASS");
  return 0;
}
