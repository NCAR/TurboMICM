#undef NDEBUG
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sparse_matrix.h"
static void check_print(FILE *f,const char *expected)
{
    size_t n=strlen(expected); char *buf=malloc(n+1); assert(buf); rewind(f);
    assert(fread(buf,1,n,f)==n); buf[n]=0; assert(!strcmp(buf,expected)); assert(fgetc(f)==EOF); free(buf);
}
static void run(size_t width,ChemSparseOrdering ordering)
{
    ChemSparseMatrix m={0},copy={0};
    assert(chem_sparse_matrix_init(&m,3,1,width,ordering,NULL,0,0)==CHEM_MATRIX_OK);
    assert(!m.nonzeros && m.values.storage_size==0 && chem_sparse_matrix_entry(&m,0,0)==SIZE_MAX);
    assert(!chem_sparse_matrix_at(&m,0,0,0) && !chem_sparse_matrix_const_at(&m,1,0,0));
    chem_sparse_matrix_fill(&m,2); assert(chem_sparse_matrix_clone(&copy,&m)==CHEM_MATRIX_OK);
    assert(!copy.nonzeros && !chem_sparse_matrix_const_at(&copy,0,1,1));
    const ChemSparseEntry entries[]={{0,1},{3,2},{0,1},{2,3},{1,1},{2,1}};
    assert(chem_sparse_matrix_init(&m,4,1,width,ordering,entries,6,0)==CHEM_MATRIX_OK);
    assert(m.nonzeros==5 && m.dimension==4 && m.values.rows==1);
    assert(chem_sparse_matrix_entry(&m,0,1)!=SIZE_MAX && chem_sparse_matrix_entry(&m,3,2)!=SIZE_MAX);
    assert(chem_sparse_matrix_entry(&m,2,1)!=SIZE_MAX && chem_sparse_matrix_entry(&m,2,2)==SIZE_MAX);
    assert(chem_sparse_matrix_entry(&m,0,0)==SIZE_MAX && chem_sparse_matrix_entry(&m,3,3)==SIZE_MAX);
    assert(chem_sparse_matrix_entry(&m,4,2)==SIZE_MAX && chem_sparse_matrix_entry(&m,1,5)==SIZE_MAX);
    assert(chem_sparse_matrix_index(&m,1,0,2)==SIZE_MAX && !chem_sparse_matrix_at(&m,0,2,2));
    assert(*chem_sparse_matrix_at(&m,0,2,1)==0);
    *chem_sparse_matrix_at(&m,0,2,1)=45; *chem_sparse_matrix_at(&m,0,3,2)=42; *chem_sparse_matrix_at(&m,0,2,3)=21;
    assert(chem_sparse_matrix_clone(&copy,&m)==CHEM_MATRIX_OK);
    assert(*chem_sparse_matrix_const_at(&copy,0,2,1)==45 && *chem_sparse_matrix_const_at(&copy,0,3,2)==42 &&
           *chem_sparse_matrix_const_at(&copy,0,2,3)==21 && copy.values.data!=m.values.data);
    size_t index32=chem_sparse_matrix_index(&m,0,3,2), index23=chem_sparse_matrix_index(&m,0,2,3);
    size_t scale=width ? width : 1;
    assert(index32==(ordering==CHEM_SPARSE_CSR ? 4 : 3)*scale);
    assert(index23==(ordering==CHEM_SPARSE_CSR ? 3 : 4)*scale);
    assert(m.values.data[index32]==42 && m.values.data[index23]==21);
    assert(chem_sparse_matrix_init(&m,4,53,width,ordering,entries,6,24)==CHEM_MATRIX_OK);
    assert(m.nonzeros==5 && *chem_sparse_matrix_at(&m,0,2,1)==24);
    *chem_sparse_matrix_at(&m,0,2,1)=45; *chem_sparse_matrix_at(&m,43,3,2)=29;
    *chem_sparse_matrix_at(&m,33,2,3)=79; *chem_sparse_matrix_at(&m,2,2,3)=64;
    assert(*chem_sparse_matrix_at(&m,0,2,1)==45 && *chem_sparse_matrix_at(&m,43,3,2)==29 &&
           *chem_sparse_matrix_at(&m,33,2,3)==79 && *chem_sparse_matrix_at(&m,2,2,3)==64);
    assert(chem_sparse_matrix_index(&m,53,0,1)==SIZE_MAX);
    for (size_t block=0;block<3;++block) {
        size_t diagonal[1], count=0;
        assert(chem_sparse_matrix_diagonal_indices(&m,block,diagonal,1,&count)==CHEM_MATRIX_OK && count==1);
        assert(diagonal[0]==chem_sparse_matrix_index(&m,block,1,1));
    }
    if (!width) {
        assert(chem_sparse_matrix_index(&m,0,2,3)==(ordering==CHEM_SPARSE_CSR ? 3 : 4));
        assert(chem_sparse_matrix_index(&m,2,2,1)==12);
    }
    if (width==2) {
        assert(chem_sparse_matrix_index(&m,0,2,3)==(ordering==CHEM_SPARSE_CSR ? 6 : 8));
        assert(chem_sparse_matrix_index(&m,2,2,1)==14);
    }
    assert(chem_sparse_matrix_init(&m,4,3,width,ordering,entries,6,24)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) *chem_sparse_matrix_at(&m,i,1,1)=1+2*i;
    chem_sparse_matrix_add_diagonal(&m,7);
    for (size_t i=0;i<3;++i) {
        assert(*chem_sparse_matrix_at(&m,i,1,1)==8+2*i);
        assert(*chem_sparse_matrix_at(&m,i,0,1)==24 && *chem_sparse_matrix_at(&m,i,2,1)==24 &&
               *chem_sparse_matrix_at(&m,i,2,3)==24 && *chem_sparse_matrix_at(&m,i,3,2)==24);
    }
    FILE *f=fopen("build/test_sparse_print.tmp","w+"); assert(f && chem_sparse_matrix_print(&m,f,0)==CHEM_MATRIX_OK);
    const char *expected="Block 0\n0,24,0,0\n0,8,0,0\n0,24,0,24\n0,0,24,0\n"
                         "Block 1\n0,24,0,0\n0,10,0,0\n0,24,0,24\n0,0,24,0\n"
                         "Block 2\n0,24,0,0\n0,12,0,0\n0,24,0,24\n0,0,24,0\n";
    check_print(f,expected); assert(fclose(f)==0 && remove("build/test_sparse_print.tmp")==0);
    chem_sparse_matrix_fill(&m,32); for (size_t i=0;i<3;++i) *chem_sparse_matrix_at(&m,i,1,1)=2+2*i;
    chem_sparse_matrix_add_diagonal(&m,5); f=fopen("build/test_sparse_print.tmp","w+"); assert(f && chem_sparse_matrix_print(&m,f,1)==CHEM_MATRIX_OK);
    expected="Block 0\n0, 1, 32\n1, 1, 7\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
             "Block 1\n0, 1, 32\n1, 1, 9\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
             "Block 2\n0, 1, 32\n1, 1, 11\n2, 1, 32\n2, 3, 32\n3, 2, 32\n";
    check_print(f,expected); assert(fclose(f)==0 && remove("build/test_sparse_print.tmp")==0);
    const ChemSparseEntry invalid[]={{3,0},{2,4},{6,7}};
    for (size_t i=0;i<3;++i) assert(chem_sparse_matrix_init(&m,3,1,width,ordering,invalid+i,1,0)==CHEM_MATRIX_INVALID);
    assert(m.dimension==4 && m.nonzeros==5);
    chem_sparse_matrix_destroy(&m); chem_sparse_matrix_destroy(&copy);
}
int main(void)
{
    for (size_t width=0;width<=4;++width) for (int order=0;order<2;++order) run(width,(ChemSparseOrdering)order);
    return 0;
}
