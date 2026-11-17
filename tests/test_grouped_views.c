#undef NDEBUG
#include <assert.h>
#include "sparse_matrix.h"
static void run(size_t width,ChemSparseOrdering order)
{
    size_t L=width ? width : 1;
    ChemMatrix dense={0},ref={0}; ChemSparseMatrix sparse={0},other={0};
    assert(chem_matrix_init(&dense,3*L,4,width,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_init(&ref,3*L,4,width,0)==CHEM_MATRIX_OK);
    for (size_t r=0;r<3*L;++r) for (size_t c=0;c<4;++c) {
        *chem_matrix_at(&dense,r,c)=r*10+c; *chem_matrix_at(&ref,r,c)=r*10+c;
    }
    for (size_t group=0;group<3;++group) {
        ChemGroupSpan output; ChemConstGroupSpan input;
        assert(chem_matrix_group_column(&dense,group,3,&output)==CHEM_MATRIX_OK);
        assert(chem_matrix_const_group_column(&dense,group,0,&input)==CHEM_MATRIX_OK);
        assert(output.length==L && input.length==L);
        for (size_t i=0;i<L;++i) output.base[i]=input.base[i]*3+1;
    }
    for (size_t r=0;r<3*L;++r) {
        *chem_matrix_at(&ref,r,3)=3**chem_matrix_at(&ref,r,0)+1;
        for (size_t c=0;c<4;++c) assert(*chem_matrix_at(&dense,r,c)==*chem_matrix_at(&ref,r,c));
    }
    assert(chem_matrix_init(&dense,3,4,0,0)==CHEM_MATRIX_OK);
    for (size_t r=0;r<3;++r) for (size_t c=0;c<4;++c) *chem_matrix_at(&dense,r,c)=10*r+c;
    ChemGroupSpan g; ChemConstGroupSpan src;
    assert(chem_matrix_group_column(&dense,1,2,&g)==CHEM_MATRIX_OK && g.base[0]==12);
    assert(chem_matrix_const_group_column(&dense,1,3,&src)==CHEM_MATRIX_OK && src.base[0]==13);
    g.base[0]=2*src.base[0];
    assert(*chem_matrix_at(&dense,1,2)==26 && *chem_matrix_at(&dense,0,2)==2 && *chem_matrix_at(&dense,2,2)==22);
    assert(chem_matrix_init(&dense,4,3,0,0)==CHEM_MATRIX_OK);
    for (size_t r=0;r<4;++r) for (size_t c=0;c<3;++c) *chem_matrix_at(&dense,r,c)=r+5*c;
    for (size_t r=0;r<4;++r) {
        ChemConstGroupSpan x,y; ChemGroupSpan out;
        assert(chem_matrix_const_group_column(&dense,r,0,&x)==CHEM_MATRIX_OK);
        assert(chem_matrix_const_group_column(&dense,r,1,&y)==CHEM_MATRIX_OK);
        assert(chem_matrix_group_column(&dense,r,2,&out)==CHEM_MATRIX_OK);
        out.base[0]=x.base[0]-y.base[0];
    }
    for (size_t r=0;r<4;++r) assert(*chem_matrix_at(&dense,r,0)==r && *chem_matrix_at(&dense,r,1)==r+5 && *chem_matrix_at(&dense,r,2)==-5);
    assert(chem_matrix_init(&dense,8,3,4,0)==CHEM_MATRIX_OK);
    for (size_t r=0;r<8;++r) for (size_t c=0;c<3;++c) *chem_matrix_at(&dense,r,c)=100+10*r+c;
    assert(chem_matrix_group_column(&dense,1,1,&g)==CHEM_MATRIX_OK);
    assert(chem_matrix_const_group_column(&dense,1,2,&src)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) { assert(g.base[i]==100+10*(4+i)+1 && src.base[i]==100+10*(4+i)+2); g.base[i]=2*src.base[i]; }
    for (size_t r=0;r<8;++r) assert(*chem_matrix_at(&dense,r,1)==(r<4 ? 100+10*r+1 : 2*(100+10*r+2)));
    const ChemSparseEntry entries[]={{0,1},{2,0}};
    assert(chem_sparse_matrix_init(&sparse,3,L+1,width,order,entries,2,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<L+1;++i) { *chem_sparse_matrix_at(&sparse,i,0,1)=i+1; *chem_sparse_matrix_at(&sparse,i,2,0)=100+i; }
    assert(chem_sparse_matrix_group_block(&sparse,0,0,1,&g)==CHEM_MATRIX_OK);
    assert(chem_group_fill(g,42)==CHEM_MATRIX_OK);
    for (size_t i=0;i<L;++i) assert(*chem_sparse_matrix_at(&sparse,i,0,1)==42 && *chem_sparse_matrix_at(&sparse,i,2,0)==100+i);
    assert(*chem_sparse_matrix_at(&sparse,L,0,1)==L+1);
    assert(chem_sparse_matrix_const_group_block(&sparse,0,0,1,&src)==CHEM_MATRIX_OK);
    assert(chem_sparse_matrix_group_block(&sparse,0,2,0,&g)==CHEM_MATRIX_OK);
    assert(chem_group_copy(g,src)==CHEM_MATRIX_OK);
    for (size_t i=0;i<L;++i) assert(*chem_sparse_matrix_at(&sparse,i,2,0)==42);
    assert(chem_sparse_matrix_group_block(&sparse,1,0,1,&g)==CHEM_MATRIX_OK && g.length==1 && g.capacity==L);
    assert(chem_group_fill(g,-7)==CHEM_MATRIX_OK && *chem_sparse_matrix_at(&sparse,L,0,1)==-7);
    if (!width && order==CHEM_SPARSE_CSC) {
        assert(chem_sparse_matrix_init(&sparse,3,2,0,order,entries,2,0)==CHEM_MATRIX_OK);
        *chem_sparse_matrix_at(&sparse,0,0,1)=1; *chem_sparse_matrix_at(&sparse,1,0,1)=2;
        *chem_sparse_matrix_at(&sparse,1,2,0)=101;
        assert(chem_sparse_matrix_group_block(&sparse,1,0,1,&g)==CHEM_MATRIX_OK && chem_group_fill(g,9)==CHEM_MATRIX_OK);
        assert(*chem_sparse_matrix_at(&sparse,0,0,1)==1 && *chem_sparse_matrix_at(&sparse,1,2,0)==101);
        assert(chem_sparse_matrix_const_group_block(&sparse,1,0,1,&src)==CHEM_MATRIX_OK);
        assert(chem_sparse_matrix_group_block(&sparse,1,2,0,&g)==CHEM_MATRIX_OK && chem_group_copy(g,src)==CHEM_MATRIX_OK);
        assert(*chem_sparse_matrix_at(&sparse,1,2,0)==9);
    }
    const ChemSparseEntry standard_entries[]={{0,1},{1,2},{2,2}};
    assert(chem_sparse_matrix_init(&sparse,3,2,0,order,standard_entries,3,0)==CHEM_MATRIX_OK);
    *chem_sparse_matrix_at(&sparse,0,0,1)=1; *chem_sparse_matrix_at(&sparse,1,0,1)=10;
    *chem_sparse_matrix_at(&sparse,1,1,2)=20;
    assert(chem_sparse_matrix_group_block(&sparse,1,0,1,&g)==CHEM_MATRIX_OK && g.base[0]==10);
    assert(chem_sparse_matrix_const_group_block(&sparse,1,1,2,&src)==CHEM_MATRIX_OK && src.base[0]==20);
    g.base[0]=src.base[0]+5;
    assert(*chem_sparse_matrix_at(&sparse,1,0,1)==25 && *chem_sparse_matrix_at(&sparse,0,0,1)==1);
    const ChemSparseEntry vector_entries[]={{0,1},{1,2}};
    assert(chem_sparse_matrix_init(&sparse,3,8,4,order,vector_entries,2,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<8;++i) { *chem_sparse_matrix_at(&sparse,i,0,1)=100+i; *chem_sparse_matrix_at(&sparse,i,1,2)=200+i; }
    assert(chem_sparse_matrix_group_block(&sparse,1,0,1,&g)==CHEM_MATRIX_OK);
    assert(chem_sparse_matrix_const_group_block(&sparse,1,1,2,&src)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) { assert(g.base[i]==104+i && src.base[i]==204+i); g.base[i]+=src.base[i]; }
    for (size_t i=0;i<8;++i) assert(*chem_sparse_matrix_at(&sparse,i,0,1)==(i<4 ? 100+i : 300+2*i));
    assert(chem_matrix_init(&dense,4,2,4,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) *chem_matrix_at(&dense,i,1)=100+i;
    assert(chem_matrix_const_group_column(&dense,0,1,&src)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) assert(src.base[i]==100+i);
    const ChemSparseEntry e[]={{0,1}};
    assert(chem_sparse_matrix_init(&sparse,3,2*L,width,order,e,1,0)==CHEM_MATRIX_OK);
    assert(chem_matrix_init(&dense,2*L,3,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<2*L;++i) {
        *chem_sparse_matrix_at(&sparse,i,0,1)=1+i;
        for (size_t j=0;j<3;++j) *chem_matrix_at(&dense,i,j)=10+i+j;
    }
    for (size_t group=0;group<2;++group) {
        ChemConstGroupSpan d,s; ChemGroupSpan out;
        assert(chem_matrix_group_column(&dense,group,0,&out)==CHEM_MATRIX_OK);
        assert(chem_matrix_const_group_column(&dense,group,1,&d)==CHEM_MATRIX_OK);
        assert(chem_sparse_matrix_const_group_block(&sparse,group,0,1,&s)==CHEM_MATRIX_OK);
        for (size_t i=0;i<L;++i) out.base[i]=d.base[i]+s.base[i];
    }
    for (size_t i=0;i<2*L;++i) assert(*chem_matrix_at(&dense,i,0)==12+2*i && *chem_matrix_at(&dense,i,1)==11+i);
    if (width==4) {
        assert(chem_sparse_matrix_const_group_block(&sparse,1,0,1,&src)==CHEM_MATRIX_OK);
        double sum=0; for (size_t i=0;i<4;++i) sum+=src.base[i]; assert(sum==26);
        assert(chem_sparse_matrix_const_group_block(&sparse,0,0,1,&src)==CHEM_MATRIX_OK);
        for (size_t i=0;i<4;++i) assert(src.base[i]==i+1);
    }
    assert(chem_sparse_matrix_init(&other,3,2*L,width,order,e,1,0)==CHEM_MATRIX_OK);
    assert(chem_sparse_matrix_group_block(&other,0,0,1,&g)==CHEM_MATRIX_OK);
    assert(chem_sparse_matrix_const_group_block(&sparse,0,0,1,&src)==CHEM_MATRIX_OK);
    assert(chem_group_copy(g,src)==CHEM_MATRIX_OK);
    for (size_t i=0;i<L;++i) assert(*chem_sparse_matrix_at(&other,i,0,1)==i+1);
    assert(chem_group_fill(g,55.5)==CHEM_MATRIX_OK);
    for (size_t i=0;i<L;++i) assert(*chem_sparse_matrix_at(&other,i,0,1)==55.5);
    chem_matrix_destroy(&dense); chem_matrix_destroy(&ref); chem_sparse_matrix_destroy(&sparse); chem_sparse_matrix_destroy(&other);
}
int main(void)
{
    const size_t widths[]={0,1,2,4,8};
    for (size_t i=0;i<5;++i) for (int order=0;order<2;++order) run(widths[i],(ChemSparseOrdering)order);
    return 0;
}
