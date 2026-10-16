#include "test.h"
#include "chem_sparse.h"

static chem_array dense(size_t rows,size_t cols,size_t v)
{chem_array a;CHECK(chem_array_create(&a,rows,cols,v,CHEM_DOUBLE,0)==CHEM_OK);return a;}
static chem_sparse smatrix(size_t n,size_t blocks,size_t v,chem_sparse_order order,size_t count,const size_t (*positions)[2])
{
    chem_sparse_builder b;chem_sparse s;CHECK(chem_sparse_builder_create(&b,n)==CHEM_OK);b.blocks=blocks;
    for(size_t i=0;i<count;++i) CHECK(chem_sparse_builder_add(&b,positions[i][0],positions[i][1])==CHEM_OK);
    CHECK(chem_sparse_create(&s,&b,order,v,CHEM_DOUBLE)==CHEM_OK);chem_sparse_builder_destroy(&b);return s;
}
static chem_array_group group(const chem_array *a,size_t index,int writable)
{chem_array_group g;CHECK(chem_array_get_group(a,index,writable,&g)==CHEM_OK);return g;}
static size_t index0(const chem_sparse *s,size_t r,size_t c)
{size_t id;CHECK(chem_sparse_index(s,0,r,c,&id)==CHEM_OK);return id;}
static double sget(const chem_sparse *s,size_t b,size_t r,size_t c)
{double x;CHECK(chem_sparse_get(s,b,r,c,&x)==CHEM_OK);return x;}
static double times2(const double *x,size_t n,void *ctx) {(void)n;(void)ctx;return x[0]*2;}
static double plus(const double *x,size_t n,void *ctx) {(void)n;(void)ctx;return x[0]+x[1];}
static double plus5(const double *x,size_t n,void *ctx) {(void)n;(void)ctx;return x[0]+5;}
static double affine(const double *x,size_t n,void *ctx) {(void)n;(void)ctx;return x[0]*3+1;}
static double minus(const double *x,size_t n,void *ctx) {(void)n;(void)ctx;return x[0]-x[1];}
static double accum(double x,double sum,void *ctx) {(void)ctx;return sum+x;}
static double fifty_five(const double *x,size_t n,void *ctx) {(void)x;(void)n;(void)ctx;return 55.5;}
TEST(grouped_views,GroupedView,CategoryTagsAreDistinct)
{
    CHECK(CHEM_COLUMN_VIEW!=CHEM_RAW_COLUMN_VIEW);CHECK(CHEM_ELEMENT_VIEW!=CHEM_RAW_ELEMENT_VIEW);CHECK(CHEM_COLUMN_VIEW!=CHEM_ELEMENT_VIEW);
    chem_array a=dense(1,1,0);chem_array_group g=group(&a,0,1);chem_span c=chem_group_column(g,0);
    CHECK(c.kind==CHEM_COLUMN_VIEW);CHECK(c.kind!=CHEM_RAW_COLUMN_VIEW);CHECK(c.kind!=CHEM_ELEMENT_VIEW);CHECK(c.kind!=CHEM_RAW_ELEMENT_VIEW);CHECK(c.kind!=CHEM_VECTOR_VIEW);
    chem_sparse s=smatrix(1,1,0,CHEM_CSR,1,(const size_t[][2]){{0,0}});c=chem_sparse_block_span(group(&s.values,0,1),0);
    CHECK(c.kind==CHEM_ELEMENT_VIEW);CHECK(c.kind!=CHEM_RAW_ELEMENT_VIEW);CHECK(c.kind!=CHEM_COLUMN_VIEW);CHECK(c.kind!=CHEM_RAW_COLUMN_VIEW);CHECK(c.kind!=CHEM_VECTOR_VIEW);
    chem_array_destroy(&a);chem_sparse_destroy(&s);
}
TEST(grouped_views,GroupedView,MatrixColumnViewRoundTrip)
{
    chem_array a=dense(3,4,0);for(size_t r=0;r<3;++r)for(size_t c=0;c<4;++c)CHECK(chem_array_set(&a,r,c,10*r+c)==CHEM_OK);
    chem_array_group g=group(&a,1,1);chem_span x=chem_group_column(g,2),xc=chem_group_const_column(g,2),y=chem_group_const_column(g,3);
    CHECK(x.kind==CHEM_COLUMN_VIEW);CHECK(xc.kind==CHEM_COLUMN_VIEW);CHECK(y.kind==CHEM_COLUMN_VIEW);
    NEAR(*(double *)x.writable,12,0);NEAR(*(const double *)xc.data,12,0);NEAR(*(const double *)y.data,13,0);CHECK(!xc.writable);
    CHECK(chem_span_map(x,1,&y,times2,NULL)==CHEM_OK);NEAR(chem_array_get(&a,1,2),26,0);NEAR(chem_array_get(&a,0,2),2,0);NEAR(chem_array_get(&a,2,2),22,0);
    chem_array_destroy(&a);
}
TEST(grouped_views,GroupedView,VectorMatrixColumnViewContiguousBlock)
{
    chem_array a=dense(8,3,4);for(size_t r=0;r<8;++r)for(size_t c=0;c<3;++c)CHECK(chem_array_set(&a,r,c,100+10*r+c)==CHEM_OK);
    chem_array_group g=group(&a,1,1);chem_span x=chem_group_column(g,1),y=chem_group_const_column(g,2);
    CHECK(x.kind==CHEM_COLUMN_VIEW);CHECK(y.kind==CHEM_COLUMN_VIEW);
    for(size_t i=0;i<4;++i){NEAR(((double *)x.writable)[i],100+10*(4+i)+1,0);NEAR(((const double *)y.data)[i],100+10*(4+i)+2,0);}
    CHECK(chem_span_map(x,1,&y,times2,NULL)==CHEM_OK);
    for(size_t r=0;r<8;++r)NEAR(chem_array_get(&a,r,1),r<4 ? 100+10*r+1 : 2*(100+10*r+2),0);
    chem_array_destroy(&a);
}
static chem_status grouped_dense_callback(size_t n,const chem_array_group *g,void *ctx)
{
    (void)n;int plain=*(int *)ctx;chem_span in[]={chem_group_const_column(g[0],0),chem_group_const_column(g[0],1)};
    return chem_span_map(chem_group_column(g[0],plain ? 2 : 3),plain ? 2 : 1,in,plain ? minus : affine,NULL);
}
static void grouped_dense(size_t v,int plain)
{
    size_t rows=plain ? 4 : 3*v,cols=plain ? 3 : 4;chem_array a=dense(rows,cols,v),b=dense(rows,cols,v);
    for(size_t r=0;r<rows;++r)for(size_t c=0;c<cols;++c){double x=plain ? r+5*c : 10*r+c;CHECK(chem_array_set(&a,r,c,x)==CHEM_OK);CHECK(chem_array_set(&b,r,c,x)==CHEM_OK);}
    chem_array_argument arg=chem_array_mutable(&b);chem_array_function f;CHECK(chem_array_function_create(&f,1,&arg,grouped_dense_callback,&plain)==CHEM_OK);CHECK(chem_array_function_apply(&f,1,&arg)==CHEM_OK);
    for(size_t r=0;r<rows;++r)CHECK(chem_array_set(&a,r,plain ? 2 : 3,plain ? chem_array_get(&a,r,0)-chem_array_get(&a,r,1) : chem_array_get(&a,r,0)*3+1)==CHEM_OK);
    for(size_t r=0;r<rows;++r)for(size_t c=0;c<cols;++c)NEAR(chem_array_get(&b,r,c),chem_array_get(&a,r,c),0);
    /* Also bind raw views to the groups and compare their exact storage. */
    chem_raw_view raw=chem_array_column(&b,0,0);for(size_t i=0;i<chem_array_groups(&b);++i){chem_array_group g=group(&b,i,0);chem_span bound=chem_bind_view(raw,g);CHECK(bound.data==chem_group_const_column(g,0).data);}
    chem_array_function_destroy(&f);chem_array_destroy(&a);chem_array_destroy(&b);
}
TEST(grouped_views,GroupedView,DenseGroupedMatchesRawL1){grouped_dense(1,0);}
TEST(grouped_views,GroupedView,DenseGroupedMatchesRawL2){grouped_dense(2,0);}
TEST(grouped_views,GroupedView,DenseGroupedMatchesRawL4){grouped_dense(4,0);}
TEST(grouped_views,GroupedView,DenseGroupedMatchesRawL8){grouped_dense(8,0);}
TEST(grouped_views,GroupedView,PlainMatrixGroupedRoundTrip){grouped_dense(0,1);}
TEST(grouped_views,GroupedView,SparseStandardBlockViewRoundTrip)
{
    chem_sparse s=smatrix(3,2,0,CHEM_CSR,3,(const size_t[][2]){{0,1},{1,2},{2,2}});
    for(size_t b=0;b<2;++b){CHECK(chem_sparse_set(&s,b,0,1,b ? 10 : 1)==CHEM_OK);CHECK(chem_sparse_set(&s,b,1,2,b ? 20 : 2)==CHEM_OK);CHECK(chem_sparse_set(&s,b,2,2,b ? 30 : 3)==CHEM_OK);}
    chem_array_group g=group(&s.values,1,1);chem_span x=chem_sparse_block_span(g,index0(&s,0,1)),y=chem_sparse_const_block_span(g,index0(&s,1,2));
    CHECK(x.kind==CHEM_ELEMENT_VIEW);CHECK(y.kind==CHEM_ELEMENT_VIEW);NEAR(*(double *)x.writable,10,0);NEAR(*(const double *)y.data,20,0);
    CHECK(chem_span_map(x,1,&y,plus5,NULL)==CHEM_OK);NEAR(sget(&s,1,0,1),25,0);NEAR(sget(&s,0,0,1),1,0);
    chem_sparse_destroy(&s);
}
TEST(grouped_views,GroupedView,SparseVectorBlockViewContiguousBlock)
{
    chem_sparse s=smatrix(3,8,4,CHEM_CSR,2,(const size_t[][2]){{0,1},{1,2}});
    for(size_t b=0;b<8;++b){CHECK(chem_sparse_set(&s,b,0,1,100+b)==CHEM_OK);CHECK(chem_sparse_set(&s,b,1,2,200+b)==CHEM_OK);}
    chem_array_group g=group(&s.values,1,1);chem_span x=chem_sparse_block_span(g,index0(&s,0,1)),y=chem_sparse_const_block_span(g,index0(&s,1,2));
    CHECK(x.kind==CHEM_ELEMENT_VIEW);CHECK(y.kind==CHEM_ELEMENT_VIEW);
    for(size_t i=0;i<4;++i){NEAR(((double *)x.writable)[i],104+i,0);NEAR(((const double *)y.data)[i],204+i,0);}
    CHECK(chem_span_map(x,2,(chem_span[]){x,y},plus,NULL)==CHEM_OK);
    for(size_t b=0;b<8;++b)NEAR(sget(&s,b,0,1),b<4 ? 100+b : 300+2*b,0);
    chem_sparse_destroy(&s);
}
static chem_status mixed_callback(size_t n,const chem_array_group *g,void *ctx)
{(void)n;size_t k=*(size_t *)ctx;chem_span in[]={chem_group_const_column(g[1],1),chem_sparse_const_block_span(g[0],k)};return chem_span_map(chem_group_column(g[1],0),2,in,plus,NULL);}
static void mixed_grouped(size_t v)
{
    chem_sparse s=smatrix(3,2*v,v,CHEM_CSR,1,(const size_t[][2]){{0,1}});chem_array d=dense(2*v,3,v);
    for(size_t b=0;b<2*v;++b){CHECK(chem_sparse_set(&s,b,0,1,1+b)==CHEM_OK);for(size_t c=0;c<3;++c)CHECK(chem_array_set(&d,b,c,10+b+c)==CHEM_OK);}
    chem_array_argument args[]={chem_array_const(&s.values),chem_array_mutable(&d)};chem_array_function f;size_t k=index0(&s,0,1);
    CHECK(chem_array_function_create(&f,2,args,mixed_callback,&k)==CHEM_OK);CHECK(chem_array_function_apply(&f,2,args)==CHEM_OK);
    for(size_t b=0;b<2*v;++b)for(size_t c=0;c<3;++c)NEAR(chem_array_get(&d,b,c),c ? 10+b+c : 12+2*b,0);
    chem_array_function_destroy(&f);chem_sparse_destroy(&s);chem_array_destroy(&d);
}
TEST(grouped_views,GroupedView,MixedSparseAndDenseGroupedL1){mixed_grouped(1);}
TEST(grouped_views,GroupedView,MixedSparseAndDenseGroupedL2){mixed_grouped(2);}
TEST(grouped_views,GroupedView,MixedSparseAndDenseGroupedL4){mixed_grouped(4);}
TEST(grouped_views,GroupedView,MixedSparseAndDenseGroupedL8){mixed_grouped(8);}
TEST(grouped_views,GroupedView,SparseVectorRowColOverloadStillFunctional)
{
    chem_sparse s=smatrix(3,8,4,CHEM_CSR,1,(const size_t[][2]){{0,1}});for(size_t b=0;b<8;++b)CHECK(chem_sparse_set(&s,b,0,1,b+1)==CHEM_OK);
    chem_array_group g=group(&s.values,1,0);chem_span x=chem_sparse_coordinate_span(&s,g,0,1);double sum=0;
    CHECK(chem_span_reduce(x,0,&sum,accum,NULL)==CHEM_OK);NEAR(sum,26,0);
    chem_raw_view raw=chem_sparse_block(&s,index0(&s,0,1),0);CHECK(chem_bind_view(raw,g).data==x.data);chem_sparse_destroy(&s);
}
TEST(grouped_views,GroupedView,ConstGroupViewProducesGroupedViews)
{
    chem_sparse s=smatrix(3,4,4,CHEM_CSR,1,(const size_t[][2]){{0,1}});chem_array d=dense(4,2,4);
    for(size_t b=0;b<4;++b){CHECK(chem_sparse_set(&s,b,0,1,b+1)==CHEM_OK);CHECK(chem_array_set(&d,b,1,100+b)==CHEM_OK);}
    chem_span x=chem_sparse_const_block_span(group(&s.values,0,0),index0(&s,0,1)),y=chem_group_const_column(group(&d,0,0),1);
    CHECK(x.kind==CHEM_ELEMENT_VIEW);CHECK(y.kind==CHEM_COLUMN_VIEW);CHECK(!x.writable && !y.writable);
    for(size_t i=0;i<4;++i){NEAR(((const double *)x.data)[i],i+1,0);NEAR(((const double *)y.data)[i],100+i,0);}
    CHECK(chem_span_fill(x,3)==CHEM_INVALID);CHECK(chem_span_fill(y,3)==CHEM_INVALID);chem_sparse_destroy(&s);chem_array_destroy(&d);
}
static void fill_copy_group(size_t v,int order)
{
    size_t width=v ? v : 1;chem_sparse s=smatrix(3,width+1,v,order,2,(const size_t[][2]){{0,1},{2,0}});
    for(size_t b=0;b<width+1;++b){CHECK(chem_sparse_set(&s,b,0,1,b+1)==CHEM_OK);CHECK(chem_sparse_set(&s,b,2,0,100+b)==CHEM_OK);}
    size_t which=!v && order ? 1 : 0;double value=which ? 9 : 42;chem_array_group g=group(&s.values,which,1);
    chem_span x=chem_sparse_block_span(g,index0(&s,0,1)),y=chem_sparse_block_span(g,index0(&s,2,0));CHECK(chem_span_fill(x,value)==CHEM_OK);
    for(size_t b=0;b<width+1;++b){NEAR(sget(&s,b,0,1),b>=which*width && b<(which+1)*width ? value : b+1,0);NEAR(sget(&s,b,2,0),100+b,0);}
    CHECK(chem_span_copy(y,chem_span_const(x))==CHEM_OK);for(size_t b=which*width;b<(which+1)*width;++b){NEAR(sget(&s,b,0,1),value,0);NEAR(sget(&s,b,2,0),value,0);}
    if(v){g=group(&s.values,1,1);CHECK(chem_span_fill(chem_sparse_block_span(g,index0(&s,0,1)),-7)==CHEM_OK);NEAR(sget(&s,width,0,1),-7,0);}
    chem_sparse_destroy(&s);
}
TEST(grouped_views,GroupedView,FillAndCopySparseStandardCSR){fill_copy_group(0,0);}
TEST(grouped_views,GroupedView,FillAndCopySparseStandardCSC){fill_copy_group(0,1);}
TEST(grouped_views,GroupedView,FillAndCopySparseVectorCSR_L4){fill_copy_group(4,0);}
TEST(grouped_views,GroupedView,FillAndCopySparseVectorCSC_L4){fill_copy_group(4,1);}
TEST(grouped_views,GroupedView,FillAndCopySparseVectorCSR_L8){fill_copy_group(8,0);}
TEST(grouped_views,GroupedView,FillCopyCrossGroupViewSources)
{
    chem_sparse s=smatrix(2,4,4,CHEM_CSR,1,(const size_t[][2]){{0,1}}),d={0},r={0};CHECK(chem_sparse_clone(&d,&s)==CHEM_OK);CHECK(chem_sparse_clone(&r,&s)==CHEM_OK);
    for(size_t b=0;b<4;++b){CHECK(chem_sparse_set(&s,b,0,1,b+1)==CHEM_OK);CHECK(chem_sparse_set(&r,b,0,1,999)==CHEM_OK);}
    chem_span src=chem_sparse_const_block_span(group(&s.values,0,0),index0(&s,0,1)),dst=chem_sparse_block_span(group(&d.values,0,1),index0(&d,0,1));
    CHECK(chem_span_copy(dst,src)==CHEM_OK);for(size_t b=0;b<4;++b)NEAR(sget(&d,b,0,1),b+1,0);
    chem_span ref=chem_sparse_block_span(group(&r.values,0,1),index0(&r,0,1));CHECK(chem_span_map(ref,0,NULL,fifty_five,NULL)==CHEM_OK);CHECK(chem_span_fill(dst,55.5)==CHEM_OK);
    for(size_t b=0;b<4;++b)NEAR(sget(&d,b,0,1),sget(&r,b,0,1),0);
    chem_sparse_destroy(&s);chem_sparse_destroy(&d);chem_sparse_destroy(&r);
}
TEST(view_category_sfinae,ViewCategorySFINAE,VectorTypes)
{
    double x[10]={0};int y[10]={0};chem_span a=chem_span_vector(x,x,10,CHEM_DOUBLE),b=chem_span_vector(y,y,10,CHEM_INT);
    CHECK(a.kind==CHEM_VECTOR_VIEW);CHECK(b.kind==CHEM_VECTOR_VIEW);CHECK(a.type==CHEM_DOUBLE);CHECK(b.type==CHEM_INT);
}
TEST(view_category_sfinae,ViewCategorySFINAE,MatrixViews)
{
    for(size_t v=0;v<=4;v+=4){chem_array a=dense(4,2,v);chem_raw_view r=chem_array_column(&a,0,1);CHECK(r.kind==CHEM_RAW_COLUMN_VIEW);CHECK(r.kind!=CHEM_VECTOR_VIEW);chem_array_destroy(&a);}
}
TEST(view_category_sfinae,ViewCategorySFINAE,SparseMatrixViews)
{
    chem_sparse s=smatrix(2,1,0,CHEM_CSR,1,(const size_t[][2]){{0,1}});chem_raw_view a=chem_sparse_block(&s,0,1),b=chem_sparse_block(&s,0,0);
    CHECK(a.kind==CHEM_RAW_ELEMENT_VIEW);CHECK(b.kind==CHEM_RAW_ELEMENT_VIEW);CHECK(a.kind!=CHEM_VECTOR_VIEW);CHECK(b.kind!=CHEM_VECTOR_VIEW);chem_sparse_destroy(&s);
}
TEST(view_category_sfinae,ViewCategorySFINAE,MatrixTypesExcluded)
{
    chem_array a=dense(1,1,0),b=dense(4,1,4);chem_sparse s=smatrix(2,1,0,CHEM_CSR,1,(const size_t[][2]){{0,1}});
    CHECK(chem_group_vector(group(&a,0,1)).status==CHEM_INVALID);CHECK(chem_group_vector(group(&b,0,1)).status==CHEM_INVALID);CHECK(chem_group_vector(group(&s.values,0,1)).status==CHEM_INVALID);
    chem_array_destroy(&a);chem_array_destroy(&b);chem_sparse_destroy(&s);
}
TEST(view_category_sfinae,ViewCategorySFINAE,RuntimeVectorUse)
{const double x[]={1,2,3,4,5};double sum=0;CHECK(chem_span_reduce(chem_span_vector(x,NULL,5,CHEM_DOUBLE),1,&sum,accum,NULL)==CHEM_OK);NEAR(sum,15,0);}
