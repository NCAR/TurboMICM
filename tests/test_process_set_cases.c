#include "chem.h"
#include "jacobian_verify.h"
#include "matrix.h"
#include "sparse.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int close_value(double a,double b) { return fabs(a-b)<=1e-12*fmax(1,fabs(b)); }
static int process_set(size_t width)
{
    const char *names[]={"foo","bar","baz","quz","quuz","corge"};
    ChemTerm rs[]={{0,1},{2,1},{1,1},{-2,1},{3,1},{2,1},{-2,1}};
    ChemTerm ps[]={{1,1},{4,2.4},{0,1},{3,1.4},{1,1},{3,2.5}};
    ChemReaction r[]={
        {.reactant_begin=0,.reactant_count=2,.product_begin=0,.product_count=2},
        {.reactant_begin=2,.reactant_count=2,.product_begin=2,.product_count=2},
        {.reactant_begin=4,.reactant_count=1,.product_count=0},
        {.reactant_begin=5,.reactant_count=2,.product_begin=4,.product_count=2}};
    ChemParameterization parameterized[]={{0,0,0,.72}};
    ChemMechanism m={.species_count=6,.reaction_count=4,.species_names=names,.reactions=r,.reactants=rs,.products=ps,.parameterized_species_count=1,.parameterized_species=parameterized};
    ChemMatrix vars={0}; if(chem_matrix_init(&vars,2,6,width,width!=1,0)) return 0;
    const double inputs[2][6]={{.1,.2,.3,.4,.5,0},{1.1,1.2,1.3,1.4,1.5,0}};
    const double density[]={70,80}; unsigned char pattern[36],expected_pattern[36]={0}; int ok=1;
    const int positions[][2]={{0,0},{0,1},{0,2},{1,0},{1,1},{1,2},{2,0},{2,2},{3,1},{3,2},{3,3},{4,0},{4,2}};
    for(size_t i=0;i<13;++i) expected_pattern[positions[i][0]*6+positions[i][1]]=1;
    if(chem_process_jacobian_pattern(&m,pattern) || memcmp(pattern,expected_pattern,sizeof pattern)) ok=0;
    for(size_t cell=0;cell<2;++cell) {
        for(size_t i=0;i<6;++i) *chem_matrix_at(&vars,cell,i)=inputs[cell][i];
        double y[6]; chem_matrix_get_row(&vars,cell,y,6);
        double d=density[cell]*.72, k[]={cell?110:10,(cell?120:20)*d,cell?130:30,(cell?140:40)*d};
        double a=k[0]*y[0]*y[2], b=k[1]*y[1], c=k[2]*y[3], e=k[3]*y[2];
        double expected[]={1000-a+b,1000+a-b+e,1000-a-e,1000+1.4*b-c+2.5*e,1000+2.4*a,1000};
        double f[]={1000,1000,1000,1000,1000,1000};
        if(chem_process_add_forcing(&m,(ChemConditions){0},k,y,NULL,f)) ok=0;
        for(size_t i=0;i<6;++i) if(!close_value(f[i],expected[i])) ok=0;
        double j[36]; for(size_t i=0;i<36;++i) j[i]=100;
        if(chem_process_subtract_jacobian(&m,(ChemConditions){0},k,y,NULL,j)) ok=0;
        double target[]={100+k[0]*y[2],100-k[1],100+k[0]*y[0],100-k[0]*y[2],100+k[1],100-k[0]*y[0]-k[3],100+k[0]*y[2],100+k[0]*y[0]+k[3],100-1.4*k[1],100-2.5*k[3],100+k[2],100-2.4*k[0]*y[2],100-2.4*k[0]*y[0]};
        for(size_t i=0;i<13;++i) if(!close_value(j[positions[i][0]*6+positions[i][1]],target[i])) ok=0;
    }
    chem_matrix_free(&vars); return ok;
}
static int algebraic_mask(size_t width)
{
    const char *names[]={"A","B","C","D"}; ChemTerm rs[]={{0,1},{1,1}},ps[]={{2,1},{3,1}};
    ChemReaction r={.reactant_count=2,.product_count=2}; ChemMechanism m={.species_count=4,.reaction_count=1,.species_names=names,.reactions=&r,.reactants=rs,.products=ps};
    unsigned char mask[]={0,1,0,1}; ChemMatrix state={0}; if(chem_matrix_init(&state,2,4,width,width!=1,0)) return 0;
    int ok=1;
    for(size_t cell=0;cell<2;++cell) {
        double y[]={1+4*cell,2+4*cell,3+4*cell,4+4*cell}, k[]={cell?20:10}, f[]={1000,1000,1000,1000},j[16];
        chem_matrix_set_row(&state,cell,y,4); for(size_t i=0;i<16;++i) j[i]=500;
        if(chem_process_add_forcing(&m,(ChemConditions){0},k,y,mask,f) || chem_process_subtract_jacobian(&m,(ChemConditions){0},k,y,mask,j)) ok=0;
        double v=k[0]*y[0]*y[1], dx=k[0]*y[1],dy=k[0]*y[0];
        double expected[]={1000-v,1000,1000+v,1000};
        for(size_t i=0;i<4;++i) if(f[i]!=expected[i]) ok=0;
        for(size_t i=0;i<4;++i) for(size_t col=0;col<2;++col) {
            double derivative=col?dy:dx, target=500+(i==0?derivative:i==2?-derivative:0);
            if(j[i*4+col]!=target) ok=0;
        }
    }
    chem_matrix_free(&state); return ok;
}
typedef struct {const ChemMechanism *m; double rates[2][2];} Context;
static void forcing(size_t block,const double *variables,double *out,void *ctx)
{
    Context *c=ctx; chem_process_add_forcing(c->m,(ChemConditions){0},c->rates[block],variables,NULL,out);
}
static int finite_difference(void)
{
    const char *names[]={"A","B","C"}; ChemTerm rs[]={{0,1},{1,1},{2,1}},ps[]={{2,2},{0,1}};
    ChemReaction rx[]={{.reactant_begin=0,.reactant_count=2,.product_begin=0,.product_count=1},{.reactant_begin=2,.reactant_count=1,.product_begin=1,.product_count=1}};
    ChemMechanism m={.species_count=3,.reaction_count=2,.species_names=names,.reactions=rx,.reactants=rs,.products=ps};
    Context ctx={&m,{{5,10},{8,12}}}; ChemMatrix vars={0},fd={0}; ChemSparseBuilder builder={0}; ChemSparseMatrix sparse={0};
    double inputs[][3]={{1,2,3},{.5,4,1.5}}; unsigned char pattern[9]; int ok=1;
    if(chem_matrix_init(&vars,2,3,1,0,0) || chem_matrix_init(&fd,2,9,1,0,0) || chem_process_jacobian_pattern(&m,pattern) || chem_sparse_builder_init(&builder,3)) return 0;
    for(size_t i=0;i<3;++i) for(size_t j=0;j<3;++j) if(pattern[i*3+j]) chem_sparse_builder_add(&builder,i,j);
    if(chem_sparse_init(&sparse,&builder,2,1,0,0)) ok=0;
    chem_sparse_builder_free(&builder);
    for(size_t b=0;b<2;++b) {
        double j[9]={0}; chem_matrix_set_row(&vars,b,inputs[b],3);
        if(chem_process_subtract_jacobian(&m,(ChemConditions){0},ctx.rates[b],inputs[b],NULL,j)) ok=0;
        for(size_t i=0;i<3;++i) for(size_t k=0;k<3;++k) if(pattern[i*3+k]) *chem_sparse_at(&sparse,b,i,k)=j[i*3+k];
    }
    ChemJacobianComparison comparison,complete;
    if(chem_finite_difference_jacobian(&vars,1e-8,forcing,&ctx,&fd) || chem_compare_jacobian(&sparse,&fd,1e-7,1e-7,&comparison) || chem_check_jacobian_sparsity(&sparse,&fd,1e-6,&complete) || !comparison.passed || !complete.passed) ok=0;
    chem_sparse_free(&sparse);chem_matrix_free(&vars);chem_matrix_free(&fd);return ok;
}
static uint32_t seed=1729;
static double random_value(void) {seed=1664525*seed+1013904223;return (seed>>8)/16777216.0;}
static int random_system(size_t cells,size_t reactions,size_t species)
{
    ChemReaction *r=calloc(reactions,sizeof(*r)); ChemTerm *reactants=calloc(reactions*3,sizeof(*reactants)), *products=calloc(reactions*10,sizeof(*products));
    double *y=malloc(species*sizeof(double)),*f=malloc(species*sizeof(double)),*expected=malloc(species*sizeof(double)),*rates=malloc(reactions*sizeof(double));
    const char **names=calloc(species,sizeof(*names)); int ok=r && reactants && products && y && f && expected && rates && names;
    if(ok) {
        size_t nr=0,np=0;
        for(size_t i=0;i<reactions;++i) {
            r[i].reactant_begin=nr;r[i].product_begin=np;r[i].reactant_count=(size_t)(4*random_value());r[i].product_count=(size_t)(11*random_value());
            for(size_t j=0;j<r[i].reactant_count;++j) reactants[nr++]=(ChemTerm){(int)(species*random_value()),1};
            for(size_t j=0;j<r[i].product_count;++j) products[np++]=(ChemTerm){(int)(species*random_value()),1.2};
        }
        ChemMechanism m={.species_count=species,.reaction_count=reactions,.species_names=names,.reactions=r,.reactants=reactants,.products=products};
        for(size_t cell=0;cell<cells && ok;++cell) {
            for(size_t i=0;i<species;++i) { y[i]=.1+random_value();f[i]=expected[i]=1000; }
            for(size_t i=0;i<reactions;++i) rates[i]=.1+random_value();
            for(size_t i=0;i<reactions;++i) {
                double rate=rates[i];
                for(size_t j=0;j<r[i].reactant_count;++j) rate*=y[reactants[r[i].reactant_begin+j].species];
                for(size_t j=0;j<r[i].reactant_count;++j) expected[reactants[r[i].reactant_begin+j].species]-=rate;
                for(size_t j=0;j<r[i].product_count;++j) expected[products[r[i].product_begin+j].species]+=1.2*rate;
            }
            if(chem_process_add_forcing(&m,(ChemConditions){0},rates,y,NULL,f)) ok=0;
            for(size_t i=0;i<species;++i) if(!close_value(f[i],expected[i])) ok=0;
        }
    }
    free(r);free(reactants);free(products);free(y);free(f);free(expected);free(rates);free(names);return ok;
}
int main(void)
{
    int vector=1,mask=1;
    for(size_t w=1;w<=4;++w) { vector &= process_set(w); mask &= algebraic_mask(w); }
    printf("process_set:ProcessSet.Matrix %d\n",process_set(1));
    printf("process_set:ProcessSet.VectorMatrix %d\n",vector);
    printf("process_set:RandomProcessSet.Matrix %d\n",random_system(200,50,40) && random_system(300,30,20) && random_system(400,100,80));
    printf("process_set:ProcessSetAlgebraicVariables.CudaMatrix %d\n",mask);
    printf("process_set:ProcessSetFiniteDifferenceJacobian.Matrix %d\n",finite_difference());
    return 0;
}
