#include "test.h"
#include "tutorials.h"
#include "internal.h"
extern const ck_mechanism tutorial_no,tutorial_user,tutorial_robertson;
static int same_reaction(const ck_reaction *a,const ck_reaction *b)
{
  if(a->rate.type!=b->rate.type||a->nreactants!=b->nreactants||a->nproducts!=b->nproducts)return 0;
  for(size_t j=0;j<a->nreactants;++j)if(a->reactants[j]!=b->reactants[j])return 0;
  for(size_t j=0;j<a->nproducts;++j)if(a->products[j].species!=b->products[j].species||a->products[j].yield!=b->products[j].yield)return 0;
  return 1;
}
static int set_parameters(ck_state *s,int configuration,int robertson)
{
  for(size_t j=0;j<s->model->nc;++j){
    const char *label=ck_custom_name(s->model,j);double value=0;
    if(strstr(label,"effective radius"))value=1e-7;
    else if(strstr(label,"particle number"))value=2.5e6;
    else if(!strcmp(label,"my photolysis rate"))value=1e-10;
    else if(!strcmp(label,"my emission rate"))value=1e-20*(configuration?1e-6*CK_AVOGADRO:1);
    else if(!strcmp(label,"my loss rate"))value=1e-23;
    else if(!strcmp(label,"r1"))value=0.04;
    else if(!strcmp(label,"r2"))value=3e7/(configuration?1e-6*CK_AVOGADRO:1);
    else if(!strcmp(label,"r3"))value=1e4/(configuration?1e-6*CK_AVOGADRO:1);
    else return 1;
    for(size_t c=0;c<s->cells;++c)s->custom[j*s->cells+c]=value*(1+0.01*c);
  }
  for(size_t c=0;c<s->cells;++c){s->conditions[c]=(ck_conditions){270+15.0*c,1e5+100.0*c,0};ck_ideal_air_density(&s->conditions[c]);}
  for(size_t j=0;j<s->model->n;++j)for(size_t c=0;c<s->cells;++c)s->y[j*s->cells+c]=(robertson?1e-5:0.02)*(j+1)/(c+1);
  return ck_update_rates(s);
}
int test_configuration_equivalence(int which)
{
  const ck_mechanism *built=which==2?tutorial_hand_robertson():tutorial_hand_rates(which);
  const ck_mechanism *parsed=which==2?&tutorial_robertson:which?&tutorial_user:&tutorial_no;
  ck_model *a,*b;ck_state *sa,*sb;CHECK(ck_model_create(built,&a)==CK_OK);CHECK(ck_model_create(parsed,&b)==CK_OK);
  CHECK(a->n==b->n&&a->nr==b->nr);CHECK(ck_state_create(a,3,&sa)==CK_OK);CHECK(ck_state_create(b,3,&sb)==CK_OK);
  CHECK(set_parameters(sa,0,which==2)==0);CHECK(set_parameters(sb,1,which==2)==0);
  for(size_t r=0;r<a->nr;++r){
    const ck_reaction *rx=&built->reactions[a->reaction_order[r]];size_t match=b->nr;
    for(size_t j=0;j<b->nr;++j)if(same_reaction(rx,&parsed->reactions[b->reaction_order[j]])){CHECK(match==b->nr);match=j;}
    CHECK(match<b->nr);
    for(size_t c=0;c<3;++c)NEAR(sa->rates[r*3+c],sb->rates[match*3+c],1e-30+1e-13*fabs(sa->rates[r*3+c]));
  }
  size_t count=a->n*3;double *fa=calloc(count,sizeof(double)),*fb=calloc(count,sizeof(double));CHECK(fa&&fb);
  ck_forcing(sa,sa->y,fa);ck_forcing(sb,sb->y,fb);
  for(size_t j=0;j<a->n;++j){size_t index;CHECK(ck_species_index(b,ck_species_name(a,j),&index)==CK_OK);for(size_t c=0;c<3;++c)NEAR(fa[j*3+c],fb[index*3+c],1e-30+1e-13*fabs(fa[j*3+c]));}
  free(fa);free(fb);
  double *ja=calloc(a->nnz*3,sizeof(double)),*jb=calloc(b->nnz*3,sizeof(double));CHECK(ja&&jb);
  ck_jacobian(sa,sa->y,ja);ck_jacobian(sb,sb->y,jb);
  for(size_t i=0;i<a->n;++i)for(size_t j=0;j<a->n;++j){size_t ri,cj;CHECK(ck_species_index(b,ck_species_name(a,i),&ri)==CK_OK);CHECK(ck_species_index(b,ck_species_name(a,j),&cj)==CK_OK);
    size_t ia=ck_jacobian_index(a,i,j),ib=ck_jacobian_index(b,ri,cj);
    for(size_t c=0;c<3;++c){double va=ia==SIZE_MAX?0:ja[ia*3+c],vb=ib==SIZE_MAX?0:jb[ib*3+c];NEAR(va,vb,1e-30+1e-13*fabs(va));}}
  free(ja);free(jb);
  ck_rosenbrock_parameters p;CHECK(ck_rosenbrock_parameters_init(3,&p)==CK_OK);
  for(size_t step=0;step<10;++step){ck_stats stats;CHECK(ck_solve(sa,0.1,&p,&stats)==CK_OK);CHECK(ck_solve(sb,0.1,&p,&stats)==CK_OK);
    for(size_t j=0;j<a->n;++j){size_t index;CHECK(ck_species_index(b,ck_species_name(a,j),&index)==CK_OK);for(size_t c=0;c<3;++c)NEAR(sa->y[j*3+c],sb->y[index*3+c],1e-30+1e-12*fabs(sa->y[j*3+c]));}}
  ck_state_free(sa);ck_state_free(sb);ck_model_free(a);ck_model_free(b);return 0;
}
