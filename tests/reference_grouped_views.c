#include "chem.h"
#include "chem_matrix.h"
#include "chem_sparse.h"

#include <stdio.h>

#define D(m,r,c) (*(double *)chem_matrix_at((m),(r),(c)))
#define S(m,b,r,c) (*(double *)chem_sparse_at((m),(b),(r),(c)))
_Static_assert(CHEM_DENSE_COLUMN_VIEW != CHEM_GROUPED_DENSE_VIEW,"dense view categories differ");
_Static_assert(CHEM_SPARSE_BLOCK_VIEW != CHEM_GROUPED_SPARSE_VIEW,"sparse view categories differ");
_Static_assert(CHEM_GROUPED_DENSE_VIEW != CHEM_GROUPED_SPARSE_VIEW,"grouped view categories differ");
static int failures;
static void report(const char *name, int ok)
{
    printf("grouped_views:GroupedView.%s %s %s\n",name,ok ? "PASS" : "FAIL",ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static int dense_roundtrip(size_t width)
{
    size_t rows=width ? 8 : 3, group=1; double zero=0;
    ChemMatrix m={0};
    if (chem_matrix_init(&m,rows,4,width,sizeof zero,&zero)!=CHEM_OK) return 0;
    for (size_t r=0;r<rows;++r) for (size_t c=0;c<4;++c) D(&m,r,c)=width ? 100+10*r+c : 10*r+c;
    ChemGroupView target; ChemConstGroupView source,read_target;
    size_t dst=width ? 1 : 2, src=width ? 2 : 3;
    int ok=chem_matrix_group_column(&m,group,dst,&target)==CHEM_OK &&
           chem_matrix_const_group_column(&m,group,dst,&read_target)==CHEM_OK &&
           chem_matrix_const_group_column(&m,group,src,&source)==CHEM_OK && target.category==CHEM_GROUPED_DENSE_VIEW;
    size_t len=width ? 4 : 1, first=width ? 4 : 1;
    ok &= target.length==len && source.length==len;
    for (size_t i=0;i<len && ok;++i) {
        double *out=(double *)(void *)target.base; const double *in=(const double *)(const void *)source.base;
        const double *check=(const double *)(const void *)read_target.base;
        ok &= out[i]==(width?100+10*(first+i)+dst:10*(first+i)+dst) && check[i]==out[i];
        out[i]=2*in[i];
    }
    for (size_t r=0;r<rows && ok;++r) {
        double expected=(width?100+10*r+dst:10*r+dst);
        if (r>=first && r<first+len) expected=2*(width?100+10*r+src:10*r+src);
        ok &= D(&m,r,dst)==expected;
    }
    chem_matrix_destroy(&m); return ok;
}
static int dense_equivalence(size_t width, int plain)
{
    size_t rows=plain?4:3*width, cols=plain?3:4, group=plain?0:width; double z=0;
    ChemMatrix a={0},b={0};
    int ok=chem_matrix_init(&a,rows,cols,group,sizeof z,&z)==CHEM_OK &&
           chem_matrix_init(&b,rows,cols,group,sizeof z,&z)==CHEM_OK;
    if (!ok) goto done;
    for (size_t r=0;r<rows;++r) for (size_t c=0;c<cols;++c)
        D(&a,r,c)=D(&b,r,c)=plain ? r+5*c : r*10+c;
    size_t n_groups=plain?rows:3;
    for (size_t i=0;i<n_groups && ok;++i) {
        ChemConstGroupView input,second; ChemGroupView out;
        ok=chem_matrix_const_group_column(&b,i,0,&input)==CHEM_OK &&
           chem_matrix_group_column(&b,i,plain?2:3,&out)==CHEM_OK;
        if (plain) ok &= chem_matrix_const_group_column(&b,i,1,&second)==CHEM_OK;
        for (size_t lane=0;lane<out.length && ok;++lane) {
            const double *x=(const double *)(const void *)input.base;
            const double *y=plain?(const double *)(const void *)second.base:NULL;
            ((double *)(void *)out.base)[lane]=plain ? x[lane]-y[lane] : x[lane]*3+1;
        }
    }
    for (size_t r=0;r<rows;++r) D(&a,r,plain?2:3)=plain?D(&a,r,0)-D(&a,r,1):D(&a,r,0)*3+1;
    for (size_t r=0;r<rows;++r) for (size_t c=0;c<cols;++c) ok &= D(&a,r,c)==D(&b,r,c);
done:
    chem_matrix_destroy(&a); chem_matrix_destroy(&b); return ok;
}
static int sparse_roundtrip(size_t group)
{
    const ChemSparseCoordinate coords[]={{0,1},{1,2},{2,2}}; double z=0; size_t blocks=group?8:2;
    ChemSparseMatrix m={0};
    if (chem_sparse_init(&m,3,blocks,coords,group?2:3,CHEM_SPARSE_ROW,group,sizeof z,&z)!=CHEM_OK) return 0;
    for (size_t b=0;b<blocks;++b) {
        S(&m,b,0,1)=group?100+b:(b?10:1); S(&m,b,1,2)=group?200+b:(b?20:2);
        if (!group) S(&m,b,2,2)=b?30:3;
    }
    ChemGroupView target; ChemConstGroupView source;
    int ok=chem_sparse_group_block(&m,1,0,1,&target)==CHEM_OK &&
           chem_sparse_const_group_block(&m,1,1,2,&source)==CHEM_OK && target.category==CHEM_GROUPED_SPARSE_VIEW;
    size_t first=group?4:1, len=group?4:1;
    for (size_t lane=0;lane<len && ok;++lane) {
        double *a=(double *)(void *)target.base; const double *b=(const double *)(const void *)source.base;
        ok &= a[lane]==(group?100+first+lane:10) && b[lane]==(group?200+first+lane:20);
        a[lane]=group?a[lane]+b[lane]:b[lane]+5;
    }
    for (size_t b=0;b<blocks && ok;++b) {
        double expected=group?(b<4?100+b:300+2*b):(b?25:1);
        ok &= S(&m,b,0,1)==expected;
    }
    chem_sparse_destroy(&m); return ok;
}
static int mixed(size_t width)
{
    const ChemSparseCoordinate c[]={{0,1}}; double z=0; size_t blocks=2*width;
    ChemSparseMatrix sparse={0}; ChemMatrix dense={0};
    int ok=chem_sparse_init(&sparse,3,blocks,c,1,CHEM_SPARSE_ROW,width,sizeof z,&z)==CHEM_OK &&
           chem_matrix_init(&dense,blocks,3,width,sizeof z,&z)==CHEM_OK;
    if (!ok) goto done;
    for (size_t b=0;b<blocks;++b) {
        S(&sparse,b,0,1)=1+b;
        for (size_t j=0;j<3;++j) D(&dense,b,j)=10+b+j;
    }
    for (size_t group=0;group<2;++group) {
        ChemGroupView out; ChemConstGroupView d,s;
        ok &= chem_matrix_group_column(&dense,group,0,&out)==CHEM_OK &&
              chem_matrix_const_group_column(&dense,group,1,&d)==CHEM_OK &&
              chem_sparse_const_group_block(&sparse,group,0,1,&s)==CHEM_OK;
        for (size_t lane=0;lane<out.length && ok;++lane)
            ((double *)(void *)out.base)[lane]=((const double *)(const void *)d.base)[lane]+((const double *)(const void *)s.base)[lane];
    }
    for (size_t b=0;b<blocks;++b) for (size_t j=0;j<3;++j)
        ok &= D(&dense,b,j)==(j==0 ? 12+2*b : 10+b+j);
done:
    chem_sparse_destroy(&sparse); chem_matrix_destroy(&dense); return ok;
}
static int raw_row_col(void)
{
    const ChemSparseCoordinate coords[]={{0,1}}; double z=0; ChemSparseMatrix m={0};
    if (chem_sparse_init(&m,3,8,coords,1,CHEM_SPARSE_ROW,4,sizeof z,&z)!=CHEM_OK) return 0;
    for (size_t b=0;b<8;++b) S(&m,b,0,1)=b+1;
    ChemConstGroupView v; double sum=0;
    int ok=chem_sparse_const_group_block(&m,1,0,1,&v)==CHEM_OK;
    for (size_t i=0;i<v.length && ok;++i) sum+=((const double *)(const void *)v.base)[i];
    ok &= sum==26;
    chem_sparse_destroy(&m); return ok;
}
static int const_group(void)
{
    const ChemSparseCoordinate c[]={{0,1}}; double z=0; ChemSparseMatrix s={0}; ChemMatrix d={0};
    int ok=chem_sparse_init(&s,3,4,c,1,CHEM_SPARSE_ROW,4,sizeof z,&z)==CHEM_OK &&
           chem_matrix_init(&d,4,2,4,sizeof z,&z)==CHEM_OK;
    if (!ok) goto done;
    for (size_t i=0;i<4;++i) { S(&s,i,0,1)=i+1; D(&d,i,1)=100+i; }
    ChemConstGroupView sv,dv;
    ok=chem_sparse_const_group_block(&s,0,0,1,&sv)==CHEM_OK && chem_matrix_const_group_column(&d,0,1,&dv)==CHEM_OK &&
       sv.category==CHEM_GROUPED_SPARSE_VIEW && dv.category==CHEM_GROUPED_DENSE_VIEW;
    for (size_t i=0;i<4;++i)
        ok &= ((const double *)(const void *)sv.base)[i]==i+1 && ((const double *)(const void *)dv.base)[i]==100+i;
done:
    chem_sparse_destroy(&s); chem_matrix_destroy(&d); return ok;
}
static int fill_copy(size_t width, int col, int standard_group)
{
    const ChemSparseCoordinate coords[]={{0,1},{2,0}}; double z=0;
    size_t L=width?width:1, blocks=L+1;
    ChemSparseMatrix m={0};
    if (chem_sparse_init(&m,3,blocks,coords,2,col?CHEM_SPARSE_COLUMN:CHEM_SPARSE_ROW,width,sizeof z,&z)!=CHEM_OK) return 0;
    for (size_t b=0;b<blocks;++b) { S(&m,b,0,1)=b+1; S(&m,b,2,0)=100+b; }
    size_t group=standard_group?1:0; double value=standard_group?9:42;
    ChemGroupView target,destination; ChemConstGroupView source;
    int ok=chem_sparse_group_block(&m,group,0,1,&target)==CHEM_OK && chem_group_fill(target,&value)==CHEM_OK;
    for (size_t b=0;b<blocks;++b) {
        int inside=b>=group*L && b<group*L+target.length;
        ok &= S(&m,b,0,1)==(inside ? value : b+1) && S(&m,b,2,0)==100+b;
    }
    ok &= chem_sparse_group_block(&m,group,2,0,&destination)==CHEM_OK &&
          chem_sparse_const_group_block(&m,group,0,1,&source)==CHEM_OK && chem_group_copy(destination,source)==CHEM_OK;
    for (size_t b=group*L;b<group*L+target.length;++b) ok &= S(&m,b,2,0)==value;
    if (width) {
        double last=-7; ok &= chem_sparse_group_block(&m,1,0,1,&target)==CHEM_OK && chem_group_fill(target,&last)==CHEM_OK;
        ok &= S(&m,L,0,1)==-7;
    }
    chem_sparse_destroy(&m); return ok;
}
static int cross_copy(void)
{
    const ChemSparseCoordinate coords[]={{0,1}}; double z=0; ChemSparseMatrix src={0},dst={0};
    int ok=chem_sparse_init(&src,2,4,coords,1,CHEM_SPARSE_ROW,4,sizeof z,&z)==CHEM_OK &&
           chem_sparse_init(&dst,2,4,coords,1,CHEM_SPARSE_ROW,4,sizeof z,&z)==CHEM_OK;
    if (!ok) goto done;
    for (size_t i=0;i<4;++i) S(&src,i,0,1)=i+1;
    ChemGroupView d; ChemConstGroupView s;
    ok=chem_sparse_group_block(&dst,0,0,1,&d)==CHEM_OK && chem_sparse_const_group_block(&src,0,0,1,&s)==CHEM_OK &&
       chem_group_copy(d,s)==CHEM_OK;
    for (size_t i=0;i<4;++i) ok &= S(&dst,i,0,1)==i+1;
    double value=55.5; ok &= chem_group_fill(d,&value)==CHEM_OK;
    for (size_t i=0;i<4;++i) ok &= S(&dst,i,0,1)==55.5;
done:
    chem_sparse_destroy(&src); chem_sparse_destroy(&dst); return ok;
}
int main(void)
{
    report("CategoryTagsAreDistinct",CHEM_DENSE_COLUMN_VIEW!=CHEM_GROUPED_DENSE_VIEW &&
           CHEM_SPARSE_BLOCK_VIEW!=CHEM_GROUPED_SPARSE_VIEW && CHEM_GROUPED_DENSE_VIEW!=CHEM_GROUPED_SPARSE_VIEW &&
           CHEM_VECTOR_VIEW!=CHEM_GROUPED_DENSE_VIEW && CHEM_VECTOR_VIEW!=CHEM_GROUPED_SPARSE_VIEW);
    report("MatrixColumnViewRoundTrip",dense_roundtrip(0));
    report("VectorMatrixColumnViewContiguousBlock",dense_roundtrip(4));
    const size_t widths[]={1,2,4,8}; const char *dense_names[]={"DenseGroupedMatchesRawL1","DenseGroupedMatchesRawL2",
        "DenseGroupedMatchesRawL4","DenseGroupedMatchesRawL8"};
    const char *mixed_names[]={"MixedSparseAndDenseGroupedL1","MixedSparseAndDenseGroupedL2",
        "MixedSparseAndDenseGroupedL4","MixedSparseAndDenseGroupedL8"};
    for (size_t i=0;i<4;++i) report(dense_names[i],dense_equivalence(widths[i],0));
    report("PlainMatrixGroupedRoundTrip",dense_equivalence(1,1));
    report("SparseStandardBlockViewRoundTrip",sparse_roundtrip(0));
    report("SparseVectorBlockViewContiguousBlock",sparse_roundtrip(4));
    for (size_t i=0;i<4;++i) report(mixed_names[i],mixed(widths[i]));
    report("SparseVectorRowColOverloadStillFunctional",raw_row_col());
    report("ConstGroupViewProducesGroupedViews",const_group());
    report("FillAndCopySparseStandardCSR",fill_copy(0,0,0));
    report("FillAndCopySparseStandardCSC",fill_copy(0,1,1));
    report("FillAndCopySparseVectorCSR_L4",fill_copy(4,0,0));
    report("FillAndCopySparseVectorCSC_L4",fill_copy(4,1,0));
    report("FillAndCopySparseVectorCSR_L8",fill_copy(8,0,0));
    report("FillCopyCrossGroupViewSources",cross_copy());
    return failures ? 1 : 0;
}
