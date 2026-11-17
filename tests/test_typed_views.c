#undef NDEBUG
#include <assert.h>
#include <stddef.h>
#include "typed_matrix.h"
#include "sparse_matrix.h"

static void sum2(size_t row,const void *const *in,void *const *out,void *context)
{ (void)row; (void)context; *(double *)out[0]=*(const double *)in[0]+*(const double *)in[1]; }
static void sum3(size_t row,const void *const *in,void *const *out,void *context)
{ (void)row; (void)context; *(double *)out[0]=*(const double *)in[0]+*(const double *)in[1]+*(const double *)in[2]; }
static void scale(size_t row,const void *const *in,void *const *out,void *context)
{ (void)row; *(double *)out[0]=*(const double *)in[0]**(const double *)context; }
static void int_scale(size_t row,const void *const *in,void *const *out,void *context)
{ (void)row; *(int *)out[0]=*(const int *)in[0]**(const int *)context; }
static void int_to_double(size_t row,const void *const *in,void *const *out,void *context)
{ (void)row; *(double *)out[0]=*(const int *)in[0]**(const double *)context; }
static void accumulate(size_t row,const void *const *in,void *const *out,void *context)
{ (void)row; (void)out; *(double *)context+=*(const double *)in[0]; }
typedef struct { const size_t *indices,*values; } Capture;
static void captured(size_t row,const void *const *in,void *const *out,void *context)
{
    (void)row; (void)in; Capture *c=context;
    for (size_t i=0;i<3;++i) *(double *)out[0]+=c->values[c->indices[i]];
    *(double *)out[1]+=c->values[c->indices[1]];
}
static void empty_callback(ChemMatrix *const *m,const ChemVector *v,size_t rows,void *context)
{ (void)m; (void)v; (void)rows; ++*(size_t *)context; }

static void dense(size_t width)
{
    ChemMatrix matrix={0}; assert(chem_matrix_init(&matrix,3,4,width,0)==CHEM_MATRIX_OK);
    for (size_t r=0;r<3;++r) for (size_t c=0;c<4;++c) *chem_matrix_at(&matrix,r,c)=10*r+c;
    const ChemMatrix *read_only=&matrix;
    ChemConstRowView inputs[3]; ChemRowView outputs[2]; double temporary[4]={0},sum=0;
    assert(chem_double_const_column(read_only,0,inputs)==CHEM_MATRIX_OK);
    assert(chem_double_const_column(read_only,1,inputs+1)==CHEM_MATRIX_OK);
    outputs[0]=chem_vector_view(temporary,3,sizeof(double));
    assert(chem_typed_rows_apply(3,inputs,2,outputs,1,sum2,NULL)==CHEM_MATRIX_OK);
    assert(temporary[0]==1 && temporary[1]==21 && temporary[2]==41);
    assert(chem_double_const_column(read_only,2,inputs)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(3,inputs,1,NULL,0,accumulate,&sum)==CHEM_MATRIX_OK && sum==36);
    assert(*chem_matrix_at(&matrix,0,0)==0 && *chem_matrix_at(&matrix,1,2)==12);

    assert(chem_matrix_init(&matrix,3,2,width,0)==CHEM_MATRIX_OK);
    const double vector[]={10,20,30},factor=2;
    inputs[0]=chem_const_vector_view(vector,3,sizeof(*vector));
    assert(chem_double_column(&matrix,0,outputs)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(3,inputs,1,outputs,1,scale,(void *)&factor)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&matrix,i,0)==20*(i+1));
    assert(chem_typed_rows_apply(2,inputs,1,outputs,1,scale,(void *)&factor)==CHEM_MATRIX_SIZE);

    assert(chem_matrix_init(&matrix,4,3,width,0)==CHEM_MATRIX_OK);
    const double v[]={100,200,300,400};
    for (size_t i=0;i<4;++i) { *chem_matrix_at(&matrix,i,0)=i+1; *chem_matrix_at(&matrix,i,1)=(i+1)*10; }
    assert(chem_double_const_column(&matrix,0,inputs)==CHEM_MATRIX_OK && chem_double_const_column(&matrix,1,inputs+1)==CHEM_MATRIX_OK);
    outputs[0]=chem_vector_view(temporary,4,sizeof(double));
    assert(chem_typed_rows_apply(4,inputs,2,outputs,1,sum2,NULL)==CHEM_MATRIX_OK);
    inputs[0]=chem_const_vector_view(temporary,4,sizeof(double)); inputs[1]=chem_const_vector_view(v,4,sizeof(double));
    assert(chem_double_column(&matrix,2,outputs)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(4,inputs,2,outputs,1,sum2,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) assert(*chem_matrix_at(&matrix,i,2)==111*(i+1));

    ChemTypedMatrix integers={0}; int zero=0;
    assert(chem_typed_matrix_init(&integers,3,2,width,sizeof(int),&zero)==CHEM_MATRIX_OK);
    const int iv[]={10,20,30},two=2;
    inputs[0]=chem_const_vector_view(iv,3,sizeof(int)); assert(chem_typed_column(&integers,0,outputs)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(3,inputs,1,outputs,1,int_scale,(void *)&two)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*(const int *)chem_typed_matrix_const_at(&integers,i,0)==20*(int)(i+1));
    chem_typed_matrix_destroy(&integers);

    assert(chem_matrix_init(&matrix,3,3,width,0)==CHEM_MATRIX_OK);
    const size_t indices[]={2,0,1},values[]={10,20,30}; Capture capture={indices,values};
    assert(chem_double_column(&matrix,0,outputs)==CHEM_MATRIX_OK && chem_double_column(&matrix,1,outputs+1)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(3,NULL,0,outputs,2,captured,&capture)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&matrix,i,0)==60 && *chem_matrix_at(&matrix,i,1)==10 && *chem_matrix_at(&matrix,i,2)==0);
    chem_matrix_destroy(&matrix);
}
static void sparse(size_t width,ChemSparseOrdering ordering)
{
    const ChemSparseEntry entries[]={{0,1},{1,2},{2,3},{3,3}};
    ChemSparseMatrix matrix={0}; assert(chem_sparse_matrix_init(&matrix,4,3,width,ordering,entries,4,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) { *chem_sparse_matrix_at(&matrix,i,0,1)=i+1; *chem_sparse_matrix_at(&matrix,i,1,2)=i+10; }
    size_t a=chem_sparse_matrix_entry(&matrix,0,1),b=chem_sparse_matrix_entry(&matrix,1,2),c=chem_sparse_matrix_entry(&matrix,2,3);
    const ChemSparseMatrix *read_only=&matrix;
    ChemConstRowView in[3]; ChemRowView out[1]; double tmp[3]={0},sum=0;
    assert(chem_double_const_column(&read_only->values,a,in)==CHEM_MATRIX_OK && chem_double_const_column(&read_only->values,b,in+1)==CHEM_MATRIX_OK);
    out[0]=chem_vector_view(tmp,3,sizeof(double));
    assert(chem_typed_rows_apply(3,in,2,out,1,sum2,NULL)==CHEM_MATRIX_OK);
    assert(tmp[0]==11 && tmp[1]==13 && tmp[2]==15);
    assert(chem_double_const_column(&read_only->values,c,in)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(3,in,1,NULL,0,accumulate,&sum)==CHEM_MATRIX_OK && sum==0);
    assert(*chem_sparse_matrix_at(&matrix,0,0,1)==1 && *chem_sparse_matrix_at(&matrix,1,1,2)==11);
    const double vector[]={5,10,15},three=3;
    in[0]=chem_const_vector_view(vector,3,sizeof(double)); assert(chem_double_column(&matrix.values,a,out)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(3,in,1,out,1,scale,(void *)&three)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&matrix,i,0,1)==15*(i+1));
    for (size_t i=0;i<3;++i) { *chem_sparse_matrix_at(&matrix,i,0,1)=i+1; *chem_sparse_matrix_at(&matrix,i,1,2)=i+2; }
    assert(chem_double_const_column(&matrix.values,a,in)==CHEM_MATRIX_OK && chem_double_const_column(&matrix.values,b,in+1)==CHEM_MATRIX_OK);
    const double more[]={100,200,300},two=2;
    in[2]=chem_const_vector_view(more,3,sizeof(double)); out[0]=chem_vector_view(tmp,3,sizeof(double));
    assert(chem_typed_rows_apply(3,in,3,out,1,sum3,NULL)==CHEM_MATRIX_OK);
    in[0]=chem_const_vector_view(tmp,3,sizeof(double)); assert(chem_double_column(&matrix.values,c,out)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(3,in,1,out,1,scale,(void *)&two)==CHEM_MATRIX_OK);
    const double expected[]={206,410,614};
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&matrix,i,2,3)==expected[i]);
    const int int_vector[]={5,10,15}; const double factor=1.5;
    in[0]=chem_const_vector_view(int_vector,3,sizeof(int)); assert(chem_double_column(&matrix.values,a,out)==CHEM_MATRIX_OK);
    assert(chem_typed_rows_apply(3,in,1,out,1,int_to_double,(void *)&factor)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&matrix,i,0,1)==7.5*(i+1));
    const double v[]={1,2,3}; in[0]=chem_const_vector_view(v,3,sizeof(double));
    assert(chem_typed_rows_apply(3,in,1,out,1,scale,(void *)&two)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&matrix,i,0,1)==2*(i+1));
    chem_sparse_matrix_destroy(&matrix);

    ChemMatrixFunction function={0}; size_t calls=0;
    assert(chem_sparse_matrix_init(&matrix,4,3,width,ordering,NULL,0,0)==CHEM_MATRIX_OK);
    ChemMatrix *arguments[]={&matrix.values};
    assert(chem_matrix_function_init(&function,arguments,1,0,empty_callback,&calls)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,arguments,1,NULL,0)==CHEM_MATRIX_OK);
    chem_matrix_function_destroy(&function); chem_sparse_matrix_destroy(&matrix);
    assert(chem_sparse_matrix_init(&matrix,4,0,width,ordering,entries,1,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_init(&function,arguments,1,0,empty_callback,&calls)==CHEM_MATRIX_OK);
    assert(chem_matrix_function_call(&function,arguments,1,NULL,0)==CHEM_MATRIX_OK);
    chem_matrix_function_destroy(&function); chem_sparse_matrix_destroy(&matrix);
}
int main(void)
{
    ChemTypedMatrix matrix={0}; int zero=0;
    assert(chem_typed_matrix_init(&matrix,3,4,0,sizeof(int),&zero)==CHEM_MATRIX_OK);
    for (size_t col=0;col<4;++col) *(int *)chem_typed_matrix_at(&matrix,1,col)=42+(int)col;
    for (size_t col=0;col<4;++col) assert(*(const int *)chem_typed_matrix_const_at(&matrix,1,col)==42+(int)col);
    assert(*(int *)chem_typed_matrix_at(&matrix,0,3)==0 && *(int *)chem_typed_matrix_at(&matrix,2,0)==0);
    chem_typed_matrix_destroy(&matrix);
    for (size_t width=0;width<=4;++width) {
        dense(width);
        for (int order=0;order<2;++order) sparse(width,(ChemSparseOrdering)order);
    }
    return 0;
}
