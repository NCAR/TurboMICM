#undef NDEBUG
#include <assert.h>
#include "matrix.h"
static void copy_scaled(size_t row,double *const *v,void *ctx)
{ (void)row; *v[1]=*v[0]**(double *)ctx; }
static void add_three(size_t row,double *const *v,void *ctx)
{ (void)row; (void)ctx; *v[2]=*v[0]+*v[1]; }
typedef struct { int mode; size_t calls; double scale; } Context;
static void invoke(ChemMatrix *const *mat,const ChemVector *vec,size_t rows,void *opaque)
{
    Context *c=opaque; ++c->calls;
    ChemColumnView v[3];
    if (c->mode==0) {
        assert(chem_matrix_column(mat[0],0,&v[0])==CHEM_MATRIX_OK);
        v[1]=chem_vector_column(vec[0].data,rows);
        assert(chem_matrix_column(mat[0],1,&v[2])==CHEM_MATRIX_OK);
        assert(chem_columns_apply(v,3,add_three,NULL)==CHEM_MATRIX_OK);
    } else if (c->mode==1) {
        assert(chem_matrix_column(mat[0],0,&v[0])==CHEM_MATRIX_OK);
        v[1]=chem_vector_column(vec[0].data,rows);
        assert(chem_columns_apply(v,2,copy_scaled,&c->scale)==CHEM_MATRIX_OK);
    } else if (c->mode==2) {
        v[0]=chem_vector_column(vec[0].data,rows); v[1]=chem_vector_column(vec[1].data,rows);
        assert(chem_matrix_column(mat[0],0,&v[2])==CHEM_MATRIX_OK);
        assert(chem_columns_apply(v,3,add_three,NULL)==CHEM_MATRIX_OK);
    } else if (c->mode==3) {
        for (size_t i=0;i<2;++i) {
            assert(chem_matrix_column(mat[i],0,&v[0])==CHEM_MATRIX_OK);
            v[1]=chem_vector_column(vec[0].data,rows);
            assert(chem_matrix_column(mat[i],1,&v[2])==CHEM_MATRIX_OK);
            assert(chem_columns_apply(v,3,add_three,NULL)==CHEM_MATRIX_OK);
        }
    } else if (c->mode==5) {
        assert(chem_matrix_column(mat[0],0,&v[0])==CHEM_MATRIX_OK);
        assert(chem_matrix_column(mat[1],0,&v[1])==CHEM_MATRIX_OK);
        assert(chem_matrix_column(mat[0],1,&v[2])==CHEM_MATRIX_OK);
        assert(chem_columns_apply(v,3,add_three,NULL)==CHEM_MATRIX_OK);
    } else {
        v[0]=chem_vector_column(vec[0].data,rows);
        assert(chem_matrix_column(mat[0],0,&v[1])==CHEM_MATRIX_OK);
        assert(chem_columns_apply(v,2,copy_scaled,&c->scale)==CHEM_MATRIX_OK);
    }
}
static void noop(ChemMatrix *const *mat,const ChemVector *vec,size_t rows,void *context)
{ (void)mat; (void)vec; (void)rows; (void)context; }
static void run(size_t width)
{
    ChemMatrix a={0},b={0}; ChemMatrixFunction function={0}; Context context={0,0,1};
    assert(chem_matrix_init(&a,5,3,width,1)==CHEM_MATRIX_OK);
    ChemMatrix *args[]={&a,&b}; double data[]={1,2,3,4,5,6,7,8,9,10};
    ChemVector vectors[]={{data,3},{data,5}};
    assert(chem_matrix_function_init(&function,args,1,1,invoke,&context)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_SIZE && context.calls==0);
    vectors[0].length=10; assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_SIZE);
    vectors[0].length=0; assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_init(&a,0,3,width,1)==CHEM_MATRIX_OK); vectors[0].length=5;
    assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_SIZE);
    vectors[0].length=0; vectors[0].data=NULL;
    assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_OK && context.calls==1);
    assert(chem_matrix_init(&a,5,4,width,1)==CHEM_MATRIX_OK);
    vectors[0]=(ChemVector){data,5}; assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_init(&a,5,2,width,1)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_init(&a,7,3,width,1)==CHEM_MATRIX_OK); vectors[0].length=7;
    assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_OK);
    assert(*chem_matrix_at(&a,0,1)==2 && *chem_matrix_at(&a,6,1)==8);
    assert(chem_matrix_init(&a,5,3,width,0)==CHEM_MATRIX_OK);
    double x[]={1,2,3,4,5},y[]={10,20,30,40,50};
    vectors[0]=(ChemVector){x,5}; vectors[1]=(ChemVector){y,3};
    context.mode=2; assert(chem_matrix_function_init(&function,args,1,2,invoke,&context)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,1,vectors,2)==CHEM_MATRIX_SIZE);
    vectors[1].length=5; assert(chem_matrix_function_call(&function,args,1,vectors,2)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) assert(*chem_matrix_at(&a,i,0)==11*(i+1));
    assert(chem_matrix_init(&a,4,2,width,0)==CHEM_MATRIX_OK && chem_matrix_init(&b,5,3,width,0)==CHEM_MATRIX_OK);
    context.mode=3; assert(chem_matrix_function_init(&function,args,2,1,invoke,&context)==CHEM_MATRIX_OK);
    vectors[0]=(ChemVector){x,4}; assert(chem_matrix_function_call(&function,args,2,vectors,1)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_init(&b,4,3,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) { *chem_matrix_at(&a,i,0)=i+1; *chem_matrix_at(&b,i,0)=i*10; x[i]=i*2; }
    assert(chem_matrix_function_call(&function,args,2,vectors,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) assert(*chem_matrix_at(&a,i,1)==3*i+1 && *chem_matrix_at(&b,i,1)==12*i);
    assert(chem_matrix_init(&a,3,2,width,0)==CHEM_MATRIX_OK);
    double v1[]={1,2,3},v2[]={5,6,7}; context.mode=4; context.scale=10;
    assert(chem_matrix_function_init(&function,args,1,1,invoke,&context)==CHEM_MATRIX_OK);
    vectors[0]=(ChemVector){v1,3}; assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&a,i,0)==10*(i+1));
    vectors[0]=(ChemVector){v2,3}; assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&a,i,0)==10*(i+5));
    vectors[0]=(ChemVector){data,5}; assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_SIZE);
    context.mode=1; context.scale=3;
    for (size_t i=0;i<3;++i) *chem_matrix_at(&a,i,0)=i+1;
    vectors[0]=(ChemVector){v2,3}; assert(chem_matrix_function_call(&function,args,1,vectors,1)==CHEM_MATRIX_OK);
    assert(v2[0]==3 && v2[1]==6 && v2[2]==9);
    assert(chem_matrix_init(&a,5,2,width,0)==CHEM_MATRIX_OK && chem_matrix_init(&b,5,3,width,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_init(&function,args,2,1,noop,NULL)==CHEM_MATRIX_OK);
    vectors[0]=(ChemVector){data,4}; assert(chem_matrix_function_call(&function,args,2,vectors,1)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_init(&a,4,4,width,1)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_init(&function,args,1,0,noop,NULL)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,1,NULL,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_column(&a,5,(ChemColumnView[1]){{0}})==CHEM_MATRIX_INVALID);
    assert(chem_matrix_init(&a,3,5,width,2)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,1,NULL,0)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_init(&a,7,4,width,1)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,1,NULL,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_init(&a,3,3,width,1)==CHEM_MATRIX_OK && chem_matrix_init(&b,4,3,width,2)==CHEM_MATRIX_OK);
    context.mode=5;
    assert(chem_matrix_function_init(&function,args,2,0,invoke,&context)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,2,NULL,0)==CHEM_MATRIX_SIZE);
    assert(chem_matrix_init(&a,5,3,width,0)==CHEM_MATRIX_OK && chem_matrix_init(&b,5,3,width,1)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,args,2,NULL,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<5;++i) assert(*chem_matrix_at(&a,i,1)==1);
    chem_matrix_function_destroy(&function); chem_matrix_destroy(&a); chem_matrix_destroy(&b);
}
int main(void) { for (size_t width=0;width<=4;++width) run(width); return 0; }
