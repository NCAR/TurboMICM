/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "test.h"
#include "chem_sparse_matrix.h"
#include <stdlib.h>

static ChemMatrix dense(size_t rows,size_t columns,size_t group)
{ChemMatrix m={0};CHECK(chem_matrix_create(rows,columns,group,CHEM_MATRIX_REAL,0,&m)==CHEM_OK);return m;}
static ChemSparseMatrix sparse(size_t dim,size_t blocks,size_t group,ChemSparseOrdering order,size_t n,const ChemSparseElement *e)
{
    ChemSparseBuilder *b=NULL;ChemSparseMatrix m={0};CHECK(chem_sparse_builder_create(dim,&b)==CHEM_OK);
    CHECK(chem_sparse_builder_blocks(b,blocks)==CHEM_OK);for(size_t i=0;i<n;++i)CHECK(chem_sparse_builder_add(b,e[i].row,e[i].column)==CHEM_OK);
    CHECK(chem_sparse_matrix_build(b,group,order,CHEM_MATRIX_REAL,&m)==CHEM_OK);chem_sparse_builder_free(b);return m;
}
static ChemView col(ChemArray a,size_t c)
{ChemView v={0};CHECK(chem_array_column(a,c,&v)==CHEM_OK);return v;}
static ChemView group(ChemView raw,size_t g)
{ChemView v={0};CHECK(chem_view_group(raw,g,&v)==CHEM_OK);return v;}
static ChemView block(ChemSparseMatrix *m,size_t r,size_t c,int constant)
{
    size_t index=0;CHECK(chem_sparse_matrix_index(m,0,r,c,&index)==CHEM_OK);ChemView v={0};
    CHECK(chem_sparse_block_view(constant?chem_array_const_sparse_matrix(m):chem_array_sparse_matrix(m),index,&v)==CHEM_OK);return v;
}
static void set(ChemSparseMatrix *m,size_t b,size_t r,size_t c,double v)
{CHECK(chem_sparse_matrix_set(m,b,r,c,v)==CHEM_OK);}
/* Mutable destination is first; read-only sources follow. */
static void transform(double *v,size_t n,void *ctx)
{
    int mode=*(int *)ctx;
    if(mode==0)v[0]=2*v[1];
    else if(mode==1)v[0]=3*v[1]+1;
    else if(mode==2)v[0]=v[1]-v[2];
    else if(mode==3)v[0]=v[1]+5;
    else if(mode==4)v[0]+=v[1];
    else if(mode==5)v[0]=v[1]+v[2];
    else v[0]=55.5;
    (void)n;
}
static void run(size_t n,ChemView *v,int mode)
{CHECK(chem_views_for_each(n,v,transform,&mode)==CHEM_OK);}
static void categories(void)
{
    _Static_assert(CHEM_VIEW_GROUPED_DENSE!=CHEM_VIEW_DENSE,"distinct dense categories");
    _Static_assert(CHEM_VIEW_GROUPED_SPARSE!=CHEM_VIEW_SPARSE,"distinct sparse categories");
    _Static_assert(CHEM_VIEW_GROUPED_DENSE!=CHEM_VIEW_GROUPED_SPARSE,"distinct grouped categories");
    ChemMatrix m=dense(1,1,0);ChemView d=group(col(chem_array_matrix(&m),0),0);
    const ChemSparseElement e={0,0};ChemSparseMatrix s=sparse(1,1,0,CHEM_SPARSE_COLUMN,1,&e);ChemView v=group(block(&s,0,0,0),0);
    CHECK(d.kind==CHEM_VIEW_GROUPED_DENSE);CHECK(d.kind!=CHEM_VIEW_DENSE);CHECK(d.kind!=CHEM_VIEW_SPARSE);CHECK(d.kind!=CHEM_VIEW_GROUPED_SPARSE);CHECK(d.kind!=CHEM_VIEW_VECTOR);
    CHECK(v.kind==CHEM_VIEW_GROUPED_SPARSE);CHECK(v.kind!=CHEM_VIEW_SPARSE);CHECK(v.kind!=CHEM_VIEW_DENSE);CHECK(v.kind!=CHEM_VIEW_GROUPED_DENSE);CHECK(v.kind!=CHEM_VIEW_VECTOR);
    chem_matrix_free(&m);chem_sparse_matrix_free(&s);
}
static void dense_roundtrip(int vector)
{
    size_t rows=vector?8:3,columns=vector?3:4;ChemMatrix m=dense(rows,columns,vector?4:0);
    for(size_t r=0;r<rows;++r)for(size_t c=0;c<columns;++c)CHECK(chem_matrix_set(&m,r,c,(vector?100:0)+10*r+c)==CHEM_OK);
    ChemArray a=chem_array_matrix(&m);size_t target=vector?1:2,source=vector?2:3;
    ChemView dst=group(col(a,target),1),src=group(chem_view_const(col(a,source)),1),same=group(chem_view_const(col(a,target)),1);
    CHECK(dst.kind==CHEM_VIEW_GROUPED_DENSE);CHECK(src.kind==CHEM_VIEW_GROUPED_DENSE);CHECK(same.kind==CHEM_VIEW_GROUPED_DENSE);
    size_t width=vector?4:1,start=vector?4:1;
    for(size_t i=0;i<width;++i){CHECK(((double *)dst.writable)[i]==(vector?100:0)+10*(start+i)+target);
        CHECK(((const double *)src.data)[i]==(vector?100:0)+10*(start+i)+source);}
    if(!vector)CHECK(*(const double *)same.data==12);
    ChemView v[]={dst,src};run(2,v,0);
    for(size_t r=0;r<rows;++r)CHECK(chem_matrix_get(&m,r,target)==(r>=start&&r<start+width?2*((vector?100:0)+10*r+source):(vector?100:0)+10*r+target));
    chem_matrix_free(&m);
}
static ChemStatus dense_group_callback(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;int plain=*(int *)ctx;size_t width=a[0].group?a[0].group:1;
    for(size_t g=0;g<a[0].rows/width+(a[0].rows%width!=0);++g){
        ChemView v[]={group(col(a[0],plain?2:3),g),group(chem_view_const(col(a[0],0)),g),{0}};
        if(plain){v[2]=group(chem_view_const(col(a[0],1)),g);run(3,v,2);}else run(2,v,1);
    }
    return CHEM_OK;
}
static void dense_equivalence(size_t width)
{
    int plain=width==0;size_t rows=plain?4:3*width,columns=plain?3:4;ChemMatrix raw=dense(rows,columns,width),m=dense(rows,columns,width);
    for(size_t r=0;r<rows;++r)for(size_t c=0;c<columns;++c){double value=plain?r+5*c:10*r+c;CHECK(chem_matrix_set(&raw,r,c,value)==CHEM_OK);CHECK(chem_matrix_set(&m,r,c,value)==CHEM_OK);}
    ChemArray a=chem_array_matrix(&m);ChemArrayFunction *f=NULL;CHECK(chem_array_function_create(1,&a,dense_group_callback,&plain,&f)==CHEM_OK);
    CHECK(chem_array_function_call(f,1,&a)==CHEM_OK);
    for(size_t r=0;r<rows;++r)CHECK(chem_matrix_set(&raw,r,plain?2:3,plain?chem_matrix_get(&raw,r,0)-chem_matrix_get(&raw,r,1):3*chem_matrix_get(&raw,r,0)+1)==CHEM_OK);
    for(size_t r=0;r<rows;++r)for(size_t c=0;c<columns;++c)CHECK(chem_matrix_get(&m,r,c)==chem_matrix_get(&raw,r,c));
    chem_array_function_free(f);chem_matrix_free(&raw);chem_matrix_free(&m);
}
static void sparse_roundtrip(int vector)
{
    const ChemSparseElement e[]={{0,1},{1,2},{2,2}};size_t blocks=vector?8:2;
    ChemSparseMatrix m=sparse(3,blocks,vector?4:0,vector?CHEM_SPARSE_ROW:CHEM_SPARSE_COLUMN,vector?2:3,e);
    for(size_t b=0;b<blocks;++b){set(&m,b,0,1,vector?100+b:b?10:1);set(&m,b,1,2,vector?200+b:b?20:2);if(!vector)set(&m,b,2,2,b?30:3);}
    ChemView dst=group(block(&m,0,1,0),1),src=group(block(&m,1,2,1),1);CHECK(dst.kind==CHEM_VIEW_GROUPED_SPARSE);CHECK(src.kind==CHEM_VIEW_GROUPED_SPARSE);
    for(size_t i=0;i<(vector?4:1);++i){CHECK(((double *)dst.writable)[i]==(vector?104+i:10));CHECK(((const double *)src.data)[i]==(vector?204+i:20));}
    ChemView v[]={dst,src};run(2,v,vector?4:3);
    if(vector){for(size_t b=4;b<8;++b)CHECK(chem_sparse_matrix_get(&m,b,0,1)==300+2*b);for(size_t b=0;b<4;++b)CHECK(chem_sparse_matrix_get(&m,b,0,1)==100+b);}
    else {CHECK(chem_sparse_matrix_get(&m,1,0,1)==25);CHECK(chem_sparse_matrix_get(&m,0,0,1)==1);}
    chem_sparse_matrix_free(&m);
}
static ChemStatus mixed_callback(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;size_t index=*(size_t *)ctx;ChemView s={0};CHECK(chem_sparse_block_view(a[0],index,&s)==CHEM_OK);s=chem_view_const(s);
    size_t width=a[0].group?a[0].group:1;
    for(size_t g=0;g<a[0].rows/width;++g){ChemView v[]={group(col(a[1],0),g),group(chem_view_const(col(a[1],1)),g),group(s,g)};run(3,v,5);}
    return CHEM_OK;
}
static void mixed(size_t width)
{
    const ChemSparseElement e={0,1};size_t blocks=2*width;
    ChemSparseMatrix ref=sparse(3,blocks,width,CHEM_SPARSE_ROW,1,&e),m=sparse(3,blocks,width,CHEM_SPARSE_ROW,1,&e);
    ChemMatrix dr=dense(blocks,3,width),d=dense(blocks,3,width);
    for(size_t b=0;b<blocks;++b){set(&ref,b,0,1,1+b);set(&m,b,0,1,1+b);for(size_t c=0;c<3;++c){CHECK(chem_matrix_set(&dr,b,c,10+b+c)==CHEM_OK);CHECK(chem_matrix_set(&d,b,c,10+b+c)==CHEM_OK);}}
    size_t index;CHECK(chem_sparse_matrix_index(&ref,0,0,1,&index)==CHEM_OK);
    ChemArray a[]={chem_array_sparse_matrix(&m),chem_array_matrix(&d)};ChemArrayFunction *f=NULL;
    CHECK(chem_array_function_create(2,a,mixed_callback,&index,&f)==CHEM_OK);CHECK(chem_array_function_call(f,2,a)==CHEM_OK);
    for(size_t b=0;b<blocks;++b){CHECK(chem_matrix_set(&dr,b,0,chem_matrix_get(&dr,b,1)+chem_sparse_matrix_get(&ref,b,0,1))==CHEM_OK);
        for(size_t c=0;c<3;++c)CHECK(chem_matrix_get(&d,b,c)==chem_matrix_get(&dr,b,c));}
    chem_array_function_free(f);chem_sparse_matrix_free(&ref);chem_sparse_matrix_free(&m);chem_matrix_free(&dr);chem_matrix_free(&d);
}
static void add(double *v,size_t n,void *ctx){(void)n;*(double *)ctx+=v[0];}
static void row_column(void)
{
    const ChemSparseElement e={0,1};ChemSparseMatrix m=sparse(3,8,4,CHEM_SPARSE_ROW,1,&e);
    for(size_t b=0;b<8;++b)set(&m,b,0,1,b+1);
    ChemView v=group(block(&m,0,1,1),1);double sum=0;CHECK(chem_views_for_each(1,&v,add,&sum)==CHEM_OK);CHECK(sum==26);
    chem_sparse_matrix_free(&m);
}
static void const_group(void)
{
    const ChemSparseElement e={0,1};ChemSparseMatrix m=sparse(3,4,4,CHEM_SPARSE_ROW,1,&e);ChemMatrix d=dense(4,2,4);
    for(size_t b=0;b<4;++b){set(&m,b,0,1,b+1);CHECK(chem_matrix_set(&d,b,1,100+b)==CHEM_OK);}
    ChemView v=group(block(&m,0,1,1),0),dv=group(col(chem_array_const_matrix(&d),1),0);
    CHECK(v.kind==CHEM_VIEW_GROUPED_SPARSE);CHECK(dv.kind==CHEM_VIEW_GROUPED_DENSE);CHECK(!v.writable);CHECK(!dv.writable);
    for(size_t i=0;i<4;++i){CHECK(((const double *)v.data)[i]==i+1);CHECK(((const double *)dv.data)[i]==100+i);}
    chem_sparse_matrix_free(&m);chem_matrix_free(&d);
}
static void fill_copy(size_t width,ChemSparseOrdering order)
{
    const ChemSparseElement e[]={{0,1},{2,0}};size_t l=width?width:1;ChemSparseMatrix m=sparse(3,l+1,width,order,2,e);
    for(size_t b=0;b<l+1;++b){set(&m,b,0,1,b+1);set(&m,b,2,0,100+b);}
    size_t target=!width&&order==CHEM_SPARSE_COLUMN?1:0;double value=target?9:42;
    ChemView first=group(block(&m,0,1,0),target),second=group(block(&m,2,0,0),target);
    CHECK(chem_view_fill(first,value)==CHEM_OK);
    for(size_t b=target*l;b<(target+1)*l;++b){CHECK(chem_sparse_matrix_get(&m,b,0,1)==value);CHECK(chem_sparse_matrix_get(&m,b,2,0)==100+b);}
    CHECK(chem_sparse_matrix_get(&m,target?0:l,0,1)==(target?1:l+1));CHECK(chem_view_copy(second,chem_view_const(first))==CHEM_OK);
    for(size_t b=target*l;b<(target+1)*l;++b){CHECK(chem_sparse_matrix_get(&m,b,2,0)==value);CHECK(chem_sparse_matrix_get(&m,b,0,1)==value);}
    if(width){ChemView tail=group(block(&m,0,1,0),(l+1)/l);CHECK(chem_view_fill(tail,-7)==CHEM_OK);CHECK(chem_sparse_matrix_get(&m,l,0,1)==-7);
        /* Additionally verify allocated padding and strict/non-strict visits. */
        for(size_t i=0;i<tail.padded_length;++i)CHECK(((double *)tail.writable)[i]==-7);
        double total=0;ChemView read=chem_view_const(tail);CHECK(chem_views_for_each_strict(1,&read,add,&total)==CHEM_OK);CHECK(total==-7);
        total=0;CHECK(chem_views_for_each(1,&read,add,&total)==CHEM_OK);CHECK(total==-7*(double)l);
    }
    chem_sparse_matrix_free(&m);
}
static void cross_group(void)
{
    const ChemSparseElement e={0,1};ChemSparseMatrix src=sparse(2,4,4,CHEM_SPARSE_ROW,1,&e),dst=sparse(2,4,4,CHEM_SPARSE_ROW,1,&e),ref=sparse(2,4,4,CHEM_SPARSE_ROW,1,&e);
    for(size_t b=0;b<4;++b){set(&src,b,0,1,b+1);set(&ref,b,0,1,999);}
    ChemView d=group(block(&dst,0,1,0),0),s=group(block(&src,0,1,1),0),r=group(block(&ref,0,1,0),0);
    CHECK(chem_view_copy(d,s)==CHEM_OK);for(size_t b=0;b<4;++b)CHECK(chem_sparse_matrix_get(&dst,b,0,1)==b+1);
    run(1,&r,6);CHECK(chem_view_fill(d,55.5)==CHEM_OK);
    for(size_t b=0;b<4;++b)CHECK(chem_sparse_matrix_get(&dst,b,0,1)==chem_sparse_matrix_get(&ref,b,0,1));
    chem_sparse_matrix_free(&src);chem_sparse_matrix_free(&dst);chem_sparse_matrix_free(&ref);
}
static void vector_types(void)
{
    double *v=calloc(5,sizeof(double));int *iv=calloc(5,sizeof(int));double a[10]={0};CHECK(v&&iv);
    ChemArray arrays[]={chem_array_vector(v,5,CHEM_MATRIX_REAL),chem_array_vector(iv,5,CHEM_MATRIX_INT),chem_array_vector(a,10,CHEM_MATRIX_REAL)};
    for(size_t i=0;i<3;++i){CHECK(arrays[i].kind==CHEM_ARRAY_VECTOR);CHECK(col(arrays[i],0).kind==CHEM_VIEW_VECTOR);}free(v);free(iv);
}
static void matrix_types(void)
{
    for(size_t g=0;g<=4;g+=4){ChemMatrix m=dense(3,2,g);ChemArray a=chem_array_matrix(&m);ChemView v=col(a,0);
        CHECK(v.kind==CHEM_VIEW_DENSE);CHECK(v.kind!=CHEM_VIEW_VECTOR);chem_matrix_free(&m);}
}
static void sparse_types(void)
{
    const ChemSparseElement e={0,0};ChemSparseMatrix m=sparse(1,1,0,CHEM_SPARSE_COLUMN,1,&e);
    ChemView v=block(&m,0,0,0),c=block(&m,0,0,1);CHECK(v.kind==CHEM_VIEW_SPARSE);CHECK(c.kind==CHEM_VIEW_SPARSE);
    CHECK(v.kind!=CHEM_VIEW_VECTOR);CHECK(c.kind!=CHEM_VIEW_VECTOR);chem_sparse_matrix_free(&m);
}
static void excluded_types(void)
{
    ChemMatrix m=dense(1,1,0),v=dense(1,1,4);const ChemSparseElement e={0,0};ChemSparseMatrix s=sparse(1,1,0,CHEM_SPARSE_COLUMN,1,&e);
    CHECK(chem_array_matrix(&m).kind!=CHEM_ARRAY_VECTOR);CHECK(chem_array_matrix(&v).kind!=CHEM_ARRAY_VECTOR);CHECK(chem_array_sparse_matrix(&s).kind!=CHEM_ARRAY_VECTOR);
    chem_matrix_free(&m);chem_matrix_free(&v);chem_sparse_matrix_free(&s);
}
static void vector_use(void)
{
    double values[]={1,2,3,4,5},sum=0;ChemView v=col(chem_array_const_vector(values,5,CHEM_MATRIX_REAL),0);
    CHECK(chem_views_for_each(1,&v,add,&sum)==CHEM_OK);CHECK(sum==15);
}
void test_grouped_views(void)
{
    RUN("grouped_views:GroupedView.CategoryTagsAreDistinct",categories());
    RUN("grouped_views:GroupedView.MatrixColumnViewRoundTrip",dense_roundtrip(0));
    RUN("grouped_views:GroupedView.VectorMatrixColumnViewContiguousBlock",dense_roundtrip(1));
    RUN("grouped_views:GroupedView.DenseGroupedMatchesRawL1",dense_equivalence(1));
    RUN("grouped_views:GroupedView.DenseGroupedMatchesRawL2",dense_equivalence(2));
    RUN("grouped_views:GroupedView.DenseGroupedMatchesRawL4",dense_equivalence(4));
    RUN("grouped_views:GroupedView.DenseGroupedMatchesRawL8",dense_equivalence(8));
    RUN("grouped_views:GroupedView.PlainMatrixGroupedRoundTrip",dense_equivalence(0));
    RUN("grouped_views:GroupedView.SparseStandardBlockViewRoundTrip",sparse_roundtrip(0));
    RUN("grouped_views:GroupedView.SparseVectorBlockViewContiguousBlock",sparse_roundtrip(1));
    RUN("grouped_views:GroupedView.MixedSparseAndDenseGroupedL1",mixed(1));
    RUN("grouped_views:GroupedView.MixedSparseAndDenseGroupedL2",mixed(2));
    RUN("grouped_views:GroupedView.MixedSparseAndDenseGroupedL4",mixed(4));
    RUN("grouped_views:GroupedView.MixedSparseAndDenseGroupedL8",mixed(8));
    RUN("grouped_views:GroupedView.SparseVectorRowColOverloadStillFunctional",row_column());
    RUN("grouped_views:GroupedView.ConstGroupViewProducesGroupedViews",const_group());
    RUN("grouped_views:GroupedView.FillAndCopySparseStandardCSR",fill_copy(0,CHEM_SPARSE_ROW));
    RUN("grouped_views:GroupedView.FillAndCopySparseStandardCSC",fill_copy(0,CHEM_SPARSE_COLUMN));
    RUN("grouped_views:GroupedView.FillAndCopySparseVectorCSR_L4",fill_copy(4,CHEM_SPARSE_ROW));
    RUN("grouped_views:GroupedView.FillAndCopySparseVectorCSC_L4",fill_copy(4,CHEM_SPARSE_COLUMN));
    RUN("grouped_views:GroupedView.FillAndCopySparseVectorCSR_L8",fill_copy(8,CHEM_SPARSE_ROW));
    RUN("grouped_views:GroupedView.FillCopyCrossGroupViewSources",cross_group());
    RUN("view_category_sfinae:ViewCategorySFINAE.VectorTypes",vector_types());
    RUN("view_category_sfinae:ViewCategorySFINAE.MatrixViews",matrix_types());
    RUN("view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews",sparse_types());
    RUN("view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded",excluded_types());
    RUN("view_category_sfinae:ViewCategorySFINAE.RuntimeVectorUse",vector_use());
}
