#undef NDEBUG
#include <assert.h>
#include "sparse_matrix.h"
typedef struct { size_t e01,e12,e23; int mode; size_t calls; } Context;
static void add(size_t row,double *const *v,void *ctx) { (void)row; (void)ctx; *v[2]=*v[0]+*v[1]; }
static void add3(size_t row,double *const *v,void *ctx) { (void)row; (void)ctx; *v[3]=*v[0]+*v[1]+*v[2]; }
static void multiply(size_t row,double *const *v,void *ctx) { (void)row; *v[1]=*v[0]**(double *)ctx; }
static ChemColumnView col(ChemMatrix *m,size_t i)
{ ChemColumnView v; assert(chem_matrix_column(m,i,&v)==CHEM_MATRIX_OK); return v; }
static void callback(ChemMatrix *const *m,const ChemVector *vectors,size_t rows,void *ctx)
{
    Context *c=ctx; ++c->calls; ChemColumnView v[4];
    if (c->mode==0) {
        v[0]=col(m[0],c->e01); v[1]=chem_vector_column(vectors[0].data,rows); v[2]=col(m[0],c->e12);
        assert(chem_columns_apply(v,3,add,NULL)==CHEM_MATRIX_OK);
    } else if (c->mode==1) {
        v[0]=col(m[0],c->e01); v[1]=col(m[0],c->e12); v[2]=col(m[1],c->e23); v[3]=col(m[0],c->e23);
        assert(chem_columns_apply(v,4,add3,NULL)==CHEM_MATRIX_OK);
    } else if (c->mode==2) {
        v[0]=chem_vector_column(vectors[0].data,rows); v[1]=chem_vector_column(vectors[1].data,rows); v[2]=col(m[0],c->e01);
        assert(chem_columns_apply(v,3,add,NULL)==CHEM_MATRIX_OK);
    } else if (c->mode==3) {
        v[0]=col(m[0],c->e01); v[1]=chem_vector_column(vectors[0].data,rows); v[2]=col(m[1],c->e01);
        assert(chem_columns_apply(v,3,add,NULL)==CHEM_MATRIX_OK);
    } else if (c->mode==5) {
        double scale=3;
        v[0]=col(m[0],c->e01); v[1]=chem_vector_column(vectors[0].data,rows);
        assert(chem_columns_apply(v,2,multiply,&scale)==CHEM_MATRIX_OK);
    } else {
        double scale=2;
        v[0]=chem_vector_column(vectors[0].data,rows); v[1]=col(m[0],c->e01);
        assert(chem_columns_apply(v,2,multiply,&scale)==CHEM_MATRIX_OK);
    }
}
static void run(size_t width,ChemSparseOrdering order)
{
    const ChemSparseEntry entries[]={{0,1},{1,2},{2,3}};
    ChemSparseMatrix a={0},b={0},c={0}; ChemMatrixFunction function={0}; Context context={0};
    assert(chem_sparse_matrix_init(&a,4,3,width,order,entries,3,0)==CHEM_MATRIX_OK);
    assert(chem_sparse_matrix_init(&b,4,3,width,order,entries,3,0)==CHEM_MATRIX_OK);
    context.e01=chem_sparse_matrix_entry(&a,0,1); context.e12=chem_sparse_matrix_entry(&a,1,2);
    context.e23=chem_sparse_matrix_entry(&a,2,3);
    ChemMatrix *args[]={&a.values,&b.values,&c.values};
    double x[]={1,2,3,4,5},z[]={4,5,6,7,8}; ChemVector vec[]={{x,3},{z,2}};
    assert(chem_matrix_function_init(&function,args,1,1,callback,&context)==CHEM_MATRIX_OK);
    vec[0].length=2; assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_SIZE);
    vec[0].length=4; assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_SIZE);
    vec[0].length=0; assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_SIZE);
    assert(context.calls==0);
    assert(chem_sparse_matrix_init(&a,4,0,width,order,entries,3,0)==CHEM_MATRIX_OK);
    vec[0].length=3; assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_SIZE);
    vec[0]=(ChemVector){NULL,0}; assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_OK && context.calls==0);
    assert(chem_sparse_matrix_init(&a,4,5,width,order,entries,3,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) *chem_sparse_matrix_at(&a,i,0,1)=i+1;
    double vv[]={10,20,30,40,50}; vec[0]=(ChemVector){vv,5};
    assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) assert(*chem_sparse_matrix_at(&a,i,1,2)==11*(i+1));
    assert(chem_sparse_matrix_init(&a,4,3,width,order,entries,3,0)==CHEM_MATRIX_OK);
    context.mode=1; assert(chem_matrix_function_init(&function,args,2,0,callback,&context)==CHEM_MATRIX_OK);
    assert(chem_sparse_matrix_init(&b,4,4,width,order,entries,3,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,2,NULL,0)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_function_init(&function,args,2,0,callback,&context)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,2,NULL,0)==CHEM_MATRIX_SIZE);
    assert(chem_sparse_matrix_init(&a,4,4,width,order,entries,3,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) { *chem_sparse_matrix_at(&a,i,0,1)=i+1; *chem_sparse_matrix_at(&a,i,1,2)=i+10; *chem_sparse_matrix_at(&b,i,2,3)=i+100; }
    assert(chem_matrix_function_call(&function,args,2,NULL,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) assert(*chem_sparse_matrix_at(&a,i,2,3)==111+3*i);
    assert(chem_sparse_matrix_init(&b,4,3,width,order,entries,3,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,2,NULL,0)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_function_call(&function,(ChemMatrix *[]){&b.values,&a.values},2,NULL,0)==CHEM_MATRIX_SIZE);
    assert(chem_sparse_matrix_init(&c,4,3,width,order,entries,3,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_init(&function,(ChemMatrix *[]){&b.values,&b.values,&c.values},3,0,callback,&context)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,(ChemMatrix *[]){&b.values,&b.values,&c.values},3,NULL,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,(ChemMatrix *[]){&b.values,&c.values,&a.values},3,NULL,0)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_function_call(&function,(ChemMatrix *[]){&a.values,&b.values,&c.values},3,NULL,0)==CHEM_MATRIX_SIZE);
    assert(chem_sparse_matrix_init(&a,4,3,width,order,entries,3,0)==CHEM_MATRIX_OK);
    vec[0]=(ChemVector){x,3}; vec[1]=(ChemVector){z,2}; context.mode=2;
    assert(chem_matrix_function_init(&function,args,1,2,callback,&context)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,1,vec,2)==CHEM_MATRIX_SIZE);
    vec[1].length=3; assert(chem_matrix_function_call(&function,args,1,vec,2)==CHEM_MATRIX_OK);
    assert(*chem_sparse_matrix_at(&a,0,0,1)==5 && *chem_sparse_matrix_at(&a,1,0,1)==7 && *chem_sparse_matrix_at(&a,2,0,1)==9);
    context.mode=3; vec[0]=(ChemVector){vv,3};
    for (size_t i=0;i<3;++i) *chem_sparse_matrix_at(&a,i,0,1)=i+1;
    assert(chem_matrix_function_init(&function,args,2,1,callback,&context)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,2,vec,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&b,i,0,1)==11*(i+1));
    assert(chem_sparse_matrix_init(&b,4,4,width,order,entries,3,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,2,vec,1)==CHEM_MATRIX_SIZE);
    assert(chem_sparse_matrix_init(&b,4,3,width,order,entries,3,0)==CHEM_MATRIX_OK);
    vec[0].length=2; assert(chem_matrix_function_call(&function,args,2,vec,1)==CHEM_MATRIX_SIZE);
    context.mode=4; assert(chem_matrix_function_init(&function,args,1,1,callback,&context)==CHEM_MATRIX_OK);
    vec[0]=(ChemVector){x,3}; assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,0,1)==2*(i+1));
    vec[0]=(ChemVector){vv,3}; assert(chem_matrix_function_call(&function,(ChemMatrix *[]){&b.values},1,vec,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&b,i,0,1)==20*(i+1));
    double five[]={5,10,15}; vec[0]=(ChemVector){five,3};
    assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,0,1)==10*(i+1));
    for (size_t i=0;i<3;++i) *chem_sparse_matrix_at(&a,i,0,1)=i+1;
    context.mode=5; assert(chem_matrix_function_call(&function,args,1,vec,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(five[i]==3*(i+1));
    ChemColumnView missing; assert(chem_sparse_matrix_block(&a,3,3,&missing)==CHEM_MATRIX_INVALID);
    if (width>1) {
        ChemMatrix dense={0}; assert(chem_matrix_init(&dense,3,4,0,0)==CHEM_MATRIX_OK);
        assert(chem_matrix_function_init(&function,(ChemMatrix *[]){&a.values,&dense},2,0,callback,&context)==CHEM_MATRIX_SIZE);
        assert(chem_matrix_init(&dense,3,4,width==2 ? 3 : 2,0)==CHEM_MATRIX_OK);
        assert(chem_matrix_function_init(&function,(ChemMatrix *[]){&a.values,&dense},2,0,callback,&context)==CHEM_MATRIX_SIZE);
        assert(chem_sparse_matrix_init(&c,4,3,width==2 ? 3 : 2,order,entries,3,0)==CHEM_MATRIX_OK);
        assert(chem_matrix_function_init(&function,(ChemMatrix *[]){&a.values,&c.values},2,0,callback,&context)==CHEM_MATRIX_SIZE);
        chem_matrix_destroy(&dense);
    }
    chem_matrix_function_destroy(&function); chem_sparse_matrix_destroy(&a); chem_sparse_matrix_destroy(&b); chem_sparse_matrix_destroy(&c);
}
int main(void) { for (size_t width=0;width<=4;++width) for (int order=0;order<2;++order) run(width,(ChemSparseOrdering)order); return 0; }
