#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <float.h>
#include "matrix.h"
static void sum4(size_t row,double *const *v,void *ctx)
{ (void)row; (void)ctx; *v[2]=4*(*v[0]+*v[1]+*v[2]); }
static void plus(size_t row,double *const *v,void *ctx)
{ (void)row; double factor=ctx ? *(double *)ctx : 1; *v[2]=factor*(*v[0]+*v[1]); }
static void multiply(size_t row,double *const *v,void *ctx)
{ (void)row; (void)ctx; *v[2]=*v[0]**v[1]; }
static void minus(size_t row,double *const *v,void *ctx)
{ (void)row; (void)ctx; *v[2]=*v[0]-*v[1]; }
static void weighted_sum(size_t row,double *const *v,void *ctx)
{ (void)row; double factor=*(double *)ctx; *v[2]=factor*(*v[0]+*v[1]+*v[2]); }
static void multiply_sum(size_t row,double *const *v,void *ctx)
{ (void)row; (void)ctx; *v[2]=*v[3]*(*v[0]+*v[1]+*v[2]); }
static void sumsq(double x,void *acc,void *ctx) { (void)ctx; *(double *)acc+=x*x; }
static void maximum(double x,void *acc,void *ctx) { (void)ctx; if (x>*(double *)acc) *(double *)acc=x; }
static void any_large(double x,void *acc,void *ctx) { (void)ctx; *(int *)acc|=x>10; }
static void all_finite(double x,void *acc,void *ctx) { (void)ctx; *(int *)acc&=isfinite(x); }
static void increment(double x,void *acc,void *ctx) { (void)x; (void)ctx; ++*(size_t *)acc; }
static ChemColumnView col(ChemMatrix *m,size_t i)
{ ChemColumnView c; assert(chem_matrix_column(m,i,&c)==CHEM_MATRIX_OK); return c; }
static void run(size_t width)
{
    ChemMatrix a={0},b={0};
    assert(chem_matrix_init(&a,5,3,width,-1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) for (size_t j=0;j<3;++j) *chem_matrix_at(&a,i,j)=(double)i-2+10*j;
    ChemColumnView views[]={col(&a,0),col(&a,1),col(&a,2)};
    assert(chem_columns_apply(views,3,sum4,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) {
        assert(*chem_matrix_at(&a,i,0)==(double)i-2 && *chem_matrix_at(&a,i,1)==(double)i+8);
        assert(*chem_matrix_at(&a,i,2)==4*(3*(double)i+24));
    }
    assert(chem_matrix_init(&a,3,3,width,-1)==CHEM_MATRIX_OK);
    views[0]=col(&a,0); views[1]=col(&a,1); views[2]=col(&a,2);
    assert(chem_columns_apply(views,3,sum4,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&a,i,2)==-12);
    assert(chem_matrix_init(&a,3,2,width,0)==CHEM_MATRIX_OK && chem_matrix_init(&b,3,3,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) {
        *chem_matrix_at(&a,i,0)=i; *chem_matrix_at(&a,i,1)=i+10;
        *chem_matrix_at(&b,i,0)=i*2; *chem_matrix_at(&b,i,2)=i*4;
    }
    views[0]=col(&a,0); views[1]=col(&b,2); views[2]=col(&a,1);
    assert(chem_columns_apply(views,3,plus,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&a,i,1)==i*5 && *chem_matrix_at(&b,i,0)==i*2);
    double scale=2;
    for (size_t n=2;n<=5;n+=3) {
        assert(chem_matrix_init(&a,n,2,width,0)==CHEM_MATRIX_OK && chem_matrix_init(&b,n,3,width,0)==CHEM_MATRIX_OK);
        for (size_t i=0;i<n;++i) { *chem_matrix_at(&a,i,0)=i+1; *chem_matrix_at(&b,i,2)=i*10; }
        views[0]=col(&a,0); views[1]=col(&b,2); views[2]=col(&a,1);
        assert(chem_columns_apply(views,3,plus,&scale)==CHEM_MATRIX_OK);
        for (size_t i=0;i<n;++i) assert(*chem_matrix_at(&a,i,1)==2*(i+1+i*10));
    }
    assert(chem_matrix_init(&a,2,2,width,0)==CHEM_MATRIX_OK && chem_matrix_init(&b,2,3,width,0)==CHEM_MATRIX_OK);
    *chem_matrix_at(&a,0,0)=10; *chem_matrix_at(&a,1,0)=20;
    *chem_matrix_at(&b,0,2)=5; *chem_matrix_at(&b,1,2)=15;
    views[0]=col(&a,0); views[1]=col(&b,2); views[2]=col(&a,1); scale=2;
    assert(chem_columns_apply(views,3,plus,&scale)==CHEM_MATRIX_OK);
    assert(*chem_matrix_at(&a,0,1)==30 && *chem_matrix_at(&a,1,1)==70);
    assert(chem_matrix_init(&a,5,2,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) *chem_matrix_at(&a,i,0)=i+1;
    double vector[]={0,10,20,30,40}; scale=3;
    views[0]=col(&a,0); views[1]=chem_vector_column(vector,5); views[2]=col(&a,1);
    assert(chem_columns_apply(views,3,plus,&scale)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) assert(*chem_matrix_at(&a,i,1)==3*(i+1+i*10));
    views[1]=chem_vector_column(vector,3);
    assert(chem_columns_apply(views,3,plus,&scale)==CHEM_MATRIX_SIZE);
    views[1]=col(&b,0); assert(chem_columns_apply(views,3,plus,&scale)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_init(&a,3,2,width,0)==CHEM_MATRIX_OK);
    views[0]=col(&a,0); views[1]=chem_vector_column(vector,5); views[2]=col(&a,1);
    assert(chem_columns_apply(views,3,plus,&scale)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_column(&a,4,views)==CHEM_MATRIX_INVALID);
    assert(chem_matrix_init(&a,0,3,width,1)==CHEM_MATRIX_OK);
    views[0]=col(&a,0); views[1]=col(&a,1); views[2]=col(&a,2);
    assert(chem_columns_apply(views,3,plus,NULL)==CHEM_MATRIX_OK);
    assert(chem_matrix_init(&a,3,0,width,1)==CHEM_MATRIX_OK);
    assert(chem_columns_apply(NULL,0,plus,NULL)==CHEM_MATRIX_OK);
    assert(chem_matrix_init(&a,4,5,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) { *chem_matrix_at(&a,i,0)=i+1; *chem_matrix_at(&a,i,1)=(i+1)*10; }
    double temp1[4],temp2[4];
    views[0]=col(&a,0); views[1]=col(&a,1); views[2]=chem_vector_column(temp1,4);
    assert(chem_columns_apply(views,3,multiply,NULL)==CHEM_MATRIX_OK);
    views[2]=chem_vector_column(temp2,4); assert(chem_columns_apply(views,3,plus,NULL)==CHEM_MATRIX_OK);
    views[0]=chem_vector_column(temp1,4); views[1]=chem_vector_column(temp2,4); views[2]=col(&a,2);
    assert(chem_columns_apply(views,3,plus,NULL)==CHEM_MATRIX_OK);
    views[2]=col(&a,3); assert(chem_columns_apply(views,3,minus,NULL)==CHEM_MATRIX_OK);
    views[2]=col(&a,4); assert(chem_columns_apply(views,3,multiply,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) { double x=i+1; assert(*chem_matrix_at(&a,i,2)==10*x*x+11*x);
        assert(*chem_matrix_at(&a,i,3)==10*x*x-11*x); assert(*chem_matrix_at(&a,i,4)==110*x*x*x); }
    assert(chem_matrix_init(&a,3,4,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) *chem_matrix_at(&a,i,0)=i+1;
    double two[]={2,2,2};
    views[0]=col(&a,0); views[1]=chem_vector_column(two,3); views[2]=col(&a,1);
    assert(chem_columns_apply(views,3,multiply,NULL)==CHEM_MATRIX_OK);
    views[1]=col(&a,1); views[2]=col(&a,2); assert(chem_columns_apply(views,3,plus,NULL)==CHEM_MATRIX_OK);
    views[0]=col(&a,2); views[1]=col(&a,1); views[2]=col(&a,3); assert(chem_columns_apply(views,3,multiply,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) { double x=i+1; assert(*chem_matrix_at(&a,i,1)==2*x && *chem_matrix_at(&a,i,2)==3*x && *chem_matrix_at(&a,i,3)==6*x*x); }
    assert(chem_matrix_init(&a,2,3,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<2;++i) for (size_t j=0;j<3;++j) *chem_matrix_at(&a,i,j)=i+j;
    views[0]=col(&a,0); views[1]=col(&a,1); views[2]=col(&a,2); scale=2;
    assert(chem_columns_apply(views,3,weighted_sum,&scale)==CHEM_MATRIX_OK);
    assert(*chem_matrix_at(&a,0,2)==6 && *chem_matrix_at(&a,1,2)==12);
    assert(chem_matrix_init(&a,2,3,width,5)==CHEM_MATRIX_OK);
    views[0]=col(&a,0); views[1]=col(&a,1); views[2]=col(&a,2);
    assert(chem_columns_apply(views,3,weighted_sum,&scale)==CHEM_MATRIX_OK);
    assert(*chem_matrix_at(&a,0,2)==30 && *chem_matrix_at(&a,1,2)==30);
    assert(chem_matrix_init(&a,2,3,width,0)==CHEM_MATRIX_OK);
    *chem_matrix_at(&a,1,0)=10; views[0]=col(&a,0); views[1]=col(&a,1); views[2]=col(&a,2);
    assert(chem_columns_apply(views,3,weighted_sum,&scale)==CHEM_MATRIX_OK);
    assert(*chem_matrix_at(&a,0,2)==0 && *chem_matrix_at(&a,1,2)==20);
    assert(chem_matrix_init(&a,5,3,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) for (size_t j=0;j<3;++j) *chem_matrix_at(&a,i,j)=(double)i-2+10*j;
    double weights[]={100,200,300,400,500};
    ChemColumnView weighted[]={col(&a,0),col(&a,1),col(&a,2),chem_vector_column(weights,5)};
    assert(chem_columns_apply(weighted,4,multiply_sum,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) assert(*chem_matrix_at(&a,i,2)==weights[i]*(3*(double)i+24));
    assert(chem_matrix_init(&a,3,2,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) for (size_t j=0;j<2;++j) *chem_matrix_at(&a,i,j)=2*i+j+1;
    views[0]=col(&a,0); views[1]=col(&a,1);
    double total=0; assert(chem_column_reduce(&views[0],1,sumsq,&total,NULL)==CHEM_MATRIX_OK);
    assert(chem_column_reduce(&views[1],1,sumsq,&total,NULL)==CHEM_MATRIX_OK && total==91);
    *chem_matrix_at(&a,0,1)=8; *chem_matrix_at(&a,1,0)=5; *chem_matrix_at(&a,1,1)=2;
    *chem_matrix_at(&a,2,0)=3; *chem_matrix_at(&a,2,1)=4; double max=-DBL_MAX;
    assert(chem_column_reduce(&views[0],1,maximum,&max,NULL)==CHEM_MATRIX_OK);
    assert(chem_column_reduce(&views[1],1,maximum,&max,NULL)==CHEM_MATRIX_OK && max==8);
    chem_matrix_fill(&a,1); int any=0,all=1;
    assert(chem_column_reduce(&views[0],1,any_large,&any,NULL)==CHEM_MATRIX_OK);
    assert(chem_column_reduce(&views[1],1,any_large,&any,NULL)==CHEM_MATRIX_OK && !any);
    *chem_matrix_at(&a,2,1)=42;
    assert(chem_column_reduce(&views[1],1,any_large,&any,NULL)==CHEM_MATRIX_OK && any);
    assert(chem_column_reduce(&views[0],1,all_finite,&all,NULL)==CHEM_MATRIX_OK && all);
    *chem_matrix_at(&a,1,0)=NAN;
    assert(chem_column_reduce(&views[0],1,all_finite,&all,NULL)==CHEM_MATRIX_OK && !all);
    size_t count=0; assert(chem_column_reduce(&views[0],0,increment,&count,NULL)==CHEM_MATRIX_OK && count==3);
    chem_matrix_destroy(&a); chem_matrix_destroy(&b);
}
int main(void) { for (size_t width=0;width<=4;++width) run(width); return 0; }
