#include "matrix.h"
#include "sparse.h"
#include <stdio.h>

_Static_assert(_Generic((ChemDenseGroup){0}, ChemSparseGroup: 0, default: 1), "dense and sparse group views must be distinct");
_Static_assert(_Generic((ChemConstDenseGroup){0}, ChemConstSparseGroup: 0, default: 1), "const group views must be distinct");

static int dense_roundtrip(size_t width)
{
    size_t rows=width==1 ? 3 : 8, cols=width==1 ? 4 : 3, group=1;
    ChemMatrix m={0}; if(chem_matrix_init(&m,rows,cols,width,width!=1,0)) return 0;
    for(size_t r=0;r<rows;++r) for(size_t c=0;c<cols;++c) *chem_matrix_at(&m,r,c)=(width==1 ? 0 : 100)+10*r+c;
    size_t outcol=width==1 ? 2 : 1, incol=width==1 ? 3 : 2;
    ChemDenseGroup out; ChemConstDenseGroup in; int ok=1;
    if(chem_matrix_group_column(&m,group,outcol,&out) || chem_matrix_const_group_column(&m,group,incol,&in)) ok=0;
    else {
        for(size_t i=0;i<width;++i) {
            size_t r=group*width+i; double base=width==1?0:100;
            if(out.base[i]!=base+10*r+outcol || in.base[i]!=base+10*r+incol) ok=0;
            out.base[i]=2*in.base[i];
        }
        for(size_t r=0;r<rows;++r) {
            double expected=(width==1?0:100)+10*r+outcol;
            if(r>=group*width && r<(group+1)*width) expected=2*((width==1?0:100)+10*r+incol);
            if(*chem_matrix_at(&m,r,outcol)!=expected) ok=0;
        }
    }
    chem_matrix_free(&m); return ok;
}
static int dense_equivalence(size_t width,int plain)
{
    size_t rows=plain?4:3*width, cols=plain?3:4;
    ChemMatrix m={0}; if(chem_matrix_init(&m,rows,cols,width,!plain,0)) return 0;
    for(size_t r=0;r<rows;++r) for(size_t c=0;c<cols;++c) *chem_matrix_at(&m,r,c)=plain?r+5*c:10*r+c;
    int ok=1;
    for(size_t g=0;g<m.groups;++g) {
        ChemDenseGroup out; ChemConstDenseGroup a,b;
        if(chem_matrix_group_column(&m,g,cols-1,&out) || chem_matrix_const_group_column(&m,g,0,&a) || chem_matrix_const_group_column(&m,g,1,&b)) {ok=0;break;}
        for(size_t i=0;i<out.count;++i) out.base[i]=plain?a.base[i]-b.base[i]:3*a.base[i]+1;
    }
    for(size_t r=0;r<rows;++r) for(size_t c=0;c<cols;++c) {
        double expected=plain?r+5*c:10*r+c;
        if(c==cols-1) expected=plain ? -5.0 : (double)(30*r+1);
        if(*chem_matrix_at(&m,r,c)!=expected) ok=0;
    }
    chem_matrix_free(&m); return ok;
}
static int make_sparse(ChemSparseMatrix *m,size_t width,size_t blocks,int csc,int two)
{
    ChemSparseBuilder b={0}; if(chem_sparse_builder_init(&b,3)) return 1;
    int err=chem_sparse_builder_add(&b,0,1);
    if(two) err|=chem_sparse_builder_add(&b,2,0);
    err|=chem_sparse_init(m,&b,blocks,width,csc,0); chem_sparse_builder_free(&b); return err;
}
static int sparse_roundtrip(size_t width)
{
    size_t blocks=2*width; ChemSparseBuilder b={0}; ChemSparseMatrix m={0};
    if(chem_sparse_builder_init(&b,3)) return 0;
    chem_sparse_builder_add(&b,0,1); chem_sparse_builder_add(&b,1,2);
    if(chem_sparse_init(&m,&b,blocks,width,0,0)) return 0;
    chem_sparse_builder_free(&b);
    for(size_t i=0;i<blocks;++i) { *chem_sparse_at(&m,i,0,1)=width==1?(i?10:1):100+i; *chem_sparse_at(&m,i,1,2)=width==1?(i?20:2):200+i; }
    ChemSparseGroup a; ChemConstSparseGroup c; int ok=1;
    if(chem_sparse_group_block(&m,1,0,1,&a) || chem_sparse_const_group_block(&m,1,1,2,&c)) ok=0;
    else for(size_t i=0;i<a.count;++i) {
        double x=width==1?10:100+width+i, y=width==1?20:200+width+i;
        if(a.base[i]!=x || c.base[i]!=y) ok=0;
        a.base[i]=width==1 ? y+5 : x+y;
    }
    for(size_t i=0;i<blocks;++i) {
        double expected=width==1?(i?25:1):i<width?100+i:300+2*i;
        if(*chem_sparse_at(&m,i,0,1)!=expected) ok=0;
    }
    chem_sparse_free(&m); return ok;
}
static int mixed(size_t width)
{
    size_t blocks=2*width; ChemSparseMatrix s={0}; ChemMatrix d={0};
    if(make_sparse(&s,width,blocks,0,0) || chem_matrix_init(&d,blocks,3,width,1,0)) return 0;
    for(size_t b=0;b<blocks;++b) { *chem_sparse_at(&s,b,0,1)=1+b; for(size_t c=0;c<3;++c) *chem_matrix_at(&d,b,c)=10+b+c; }
    int ok=1;
    for(size_t g=0;g<2;++g) {
        ChemSparseGroup sview; ChemDenseGroup out; ChemConstDenseGroup in;
        if(chem_sparse_group_block(&s,g,0,1,&sview) || chem_matrix_group_column(&d,g,0,&out) || chem_matrix_const_group_column(&d,g,1,&in)) { ok=0;break; }
        for(size_t i=0;i<out.count;++i) out.base[i]=in.base[i]+sview.base[i];
    }
    for(size_t b=0;b<blocks;++b) for(size_t c=0;c<3;++c) if(*chem_matrix_at(&d,b,c)!=(c?10+b+c:12+2*b)) ok=0;
    chem_sparse_free(&s); chem_matrix_free(&d); return ok;
}
static int sparse_const(int mode)
{
    ChemSparseMatrix s={0}; ChemMatrix d={0}; if(make_sparse(&s,4,mode?4:8,0,0)) return 0;
    for(size_t b=0;b<s.blocks;++b) *chem_sparse_at(&s,b,0,1)=b+1;
    ChemConstSparseGroup v; int ok=!chem_sparse_const_group_block(&s,mode?0:1,0,1,&v);
    if(mode) {
        if(chem_matrix_init(&d,4,2,4,1,0)) ok=0;
        for(size_t b=0;b<4;++b) *chem_matrix_at(&d,b,1)=100+b;
        ChemConstDenseGroup dv; if(chem_matrix_const_group_column(&d,0,1,&dv)) ok=0;
        for(size_t i=0;i<4;++i) if(v.base[i]!=i+1 || dv.base[i]!=100+i) ok=0;
    } else { double sum=0; for(size_t i=0;i<v.count;++i) sum+=v.base[i]; if(sum!=26) ok=0; }
    chem_matrix_free(&d); chem_sparse_free(&s); return ok;
}
static int fill_copy(size_t width,int csc,int group)
{
    size_t blocks=width+1; ChemSparseMatrix m={0}; if(make_sparse(&m,width,blocks,csc,1)) return 0;
    for(size_t b=0;b<blocks;++b) { *chem_sparse_at(&m,b,0,1)=b+1; *chem_sparse_at(&m,b,2,0)=100+b; }
    ChemSparseGroup out,dest; ChemConstSparseGroup source; int ok=1; double value=group?9:42;
    if(chem_sparse_group_block(&m,group,0,1,&out) || chem_sparse_group_block(&m,group,2,0,&dest)) ok=0;
    else {
        chem_sparse_group_fill(out,value);
        for(size_t b=0;b<blocks;++b) if(*chem_sparse_at(&m,b,0,1)!=(b/width==(size_t)group?value:b+1) || *chem_sparse_at(&m,b,2,0)!=100+b) ok=0;
        if(chem_sparse_const_group_block(&m,group,0,1,&source) || chem_sparse_group_copy(dest,source)) ok=0;
        for(size_t i=0;i<dest.count;++i) if(dest.base[i]!=value) ok=0;
        if(width>1) { if(chem_sparse_group_block(&m,1,0,1,&out)) ok=0; else { chem_sparse_group_fill(out,-7); if(*chem_sparse_at(&m,width,0,1)!=-7) ok=0; } }
    }
    chem_sparse_free(&m); return ok;
}
static int cross_group(void)
{
    ChemSparseMatrix a={0},b={0}; if(make_sparse(&a,4,4,0,0)||make_sparse(&b,4,4,0,0)) return 0;
    for(size_t i=0;i<4;++i) *chem_sparse_at(&a,i,0,1)=i+1;
    ChemConstSparseGroup in; ChemSparseGroup out; int ok=!chem_sparse_const_group_block(&a,0,0,1,&in) && !chem_sparse_group_block(&b,0,0,1,&out);
    if(ok && chem_sparse_group_copy(out,in)) ok=0;
    for(size_t i=0;i<4;++i) if(*chem_sparse_at(&b,i,0,1)!=i+1) ok=0;
    chem_sparse_group_fill(out,55.5); for(size_t i=0;i<4;++i) if(*chem_sparse_at(&b,i,0,1)!=55.5) ok=0;
    chem_sparse_free(&a); chem_sparse_free(&b); return ok;
}
int main(void)
{
    const char *names[]={"CategoryTagsAreDistinct","MatrixColumnViewRoundTrip","VectorMatrixColumnViewContiguousBlock","DenseGroupedMatchesRawL1","DenseGroupedMatchesRawL2","DenseGroupedMatchesRawL4","DenseGroupedMatchesRawL8","PlainMatrixGroupedRoundTrip","SparseStandardBlockViewRoundTrip","SparseVectorBlockViewContiguousBlock","MixedSparseAndDenseGroupedL1","MixedSparseAndDenseGroupedL2","MixedSparseAndDenseGroupedL4","MixedSparseAndDenseGroupedL8","SparseVectorRowColOverloadStillFunctional","ConstGroupViewProducesGroupedViews","FillAndCopySparseStandardCSR","FillAndCopySparseStandardCSC","FillAndCopySparseVectorCSR_L4","FillAndCopySparseVectorCSC_L4","FillAndCopySparseVectorCSR_L8","FillCopyCrossGroupViewSources"};
    int results[]={1,dense_roundtrip(1),dense_roundtrip(4),dense_equivalence(1,0),dense_equivalence(2,0),dense_equivalence(4,0),dense_equivalence(8,0),dense_equivalence(1,1),sparse_roundtrip(1),sparse_roundtrip(4),mixed(1),mixed(2),mixed(4),mixed(8),sparse_const(0),sparse_const(1),fill_copy(1,0,0),fill_copy(1,1,1),fill_copy(4,0,0),fill_copy(4,1,0),fill_copy(8,0,0),cross_group()};
    for(size_t i=0;i<sizeof results/sizeof *results;++i) printf("grouped_views:GroupedView.%s %d\n",names[i],results[i]);
    return 0;
}
