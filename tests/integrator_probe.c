/* Small nonlinear, mechanism-generic prototype comparison and API checks. */
#include "chem.h"
#include "chem_integrator.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
static void array(const char *key,const double *a,size_t n)
{
  printf("%s",key);
  for (size_t i=0;i<n;++i) printf(" %.17g",a[i]);
  putchar('\n');
}
int main(int argc,char **argv)
{
  if (argc!=7 && argc!=8) return 2;
  ck_integrator_options o=ck_integrator_defaults();
  o.method=atoi(argv[1]); o.controller=atoi(argv[2]); o.jacobian_steps=(size_t)atoi(argv[3]);
  o.deadband=atof(argv[4]); o.norm=atoi(argv[5]); o.warm_start=atoi(argv[6]);
  if (argc==8) o.reuse_rhs=atoi(argv[7]);
  static const ck_species species[]={{.name="A"},{.name="B"},{.name="C"}};
  static const size_t aa[]={0,0}, b[]={1}, c[]={2};
  static const ck_product pb[]={{1,1}}, pa[]={{0,2}}, pc[]={{2,1}};
  ck_reaction rx[4]={
    {.nreactants=2,.reactants=aa,.nproducts=1,.products=pb},
    {.nreactants=1,.reactants=b,.nproducts=1,.products=pa},
    {.nreactants=1,.reactants=b,.nproducts=1,.products=pc},
    {.nreactants=1,.reactants=c,.nproducts=1,.products=pa}};
  const double rates[]={20,2,.4,.1};
  for (size_t r=0;r<4;++r) { rx[r].rate=ck_rate_defaults(CK_ARRHENIUS); rx[r].rate.A=rates[r]; }
  const ck_mechanism mech={.name="reversible-dimer",.nspecies=3,.species=species,
    .nreactions=4,.reactions=rx,.temperature=298,.pressure=100000};
  ck_model *m=NULL; ck_state *s=NULL; ck_integrator *v=NULL;
  if (ck_model_create_with_order(&mech,0,&m) || ck_state_create(m,2,&s)) return 1;
  const double initial[]={1,.5,.1,.2,.01,.03};
  for (size_t i=0;i<6;++i) s->y[i]=initial[i];
  for (size_t i=0;i<3;++i) s->absolute_tolerance[i]=1e-11;
  s->relative_tolerance=1e-5;
  if (ck_update_rates(s) || ck_integrator_create(s,&o,&v)) return 1;
  ck_rosenbrock_parameters p;
  if (!ck_rosenbrock_parameters_init(o.method,&p))
  {
    printf("parameters %zu %.17g %.17g\n",p.stages,p.order,p.gamma[0]);
    array("a",p.a,15); array("c",p.c,15); array("m",p.m,6); array("e",p.e,6);
    printf("new"); for (size_t i=0;i<6;++i) printf(" %d",p.new_function[i]); putchar('\n');
  }
  ck_stats sum={0};
  for (size_t call=0;call<2;++call)
  {
    ck_stats st;
    if (ck_integrator_advance(v,.2,&st)) return 1;
    ck_stats_add(&sum,&st);
    array("y",s->y,6);
  }
  printf("counts %llu %llu %llu %llu %llu %llu %llu\n",
    (unsigned long long)sum.internal_steps,(unsigned long long)sum.accepted,
    (unsigned long long)sum.rejected,(unsigned long long)sum.function_calls,
    (unsigned long long)sum.jacobian_updates,(unsigned long long)sum.decompositions,
    (unsigned long long)sum.solves);
  /* Mass conservation is independent of the prototype agreement test. */
  for (size_t cell=0;cell<2;++cell)
  {
    double mass=s->y[cell]+2*s->y[2+cell]+2*s->y[4+cell];
    double start=initial[cell]+2*initial[2+cell]+2*initial[4+cell];
    if (fabs(mass-start)>2e-12) return 1;
  }
  ck_stats zero;
  if (ck_integrator_advance(v,0,&zero) || zero.internal_steps ||
      ck_integrator_advance(v,-1,&zero)!=CK_INVALID) return 1;
  ck_integrator_free(v);
  o.start=NAN;
  if (ck_integrator_create(s,&o,&v)!=CK_INVALID || v) return 1;
  ck_state_free(s); ck_model_free(m);
  return 0;
}
