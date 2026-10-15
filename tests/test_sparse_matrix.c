/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Sparse matrix policy fixtures, with original offsets and output checks. */
#include "test.h"
#include "chem_sparse_matrix.h"
static const ChemSparseElement pattern[]={{0,1},{3,2},{2,3},{1,1},{2,1}};
static ChemSparseMatrix create(size_t order,size_t blocks,size_t group,ChemSparseOrdering ordering,ChemMatrixType type,size_t n,const ChemSparseElement *entries,double value)
{
    ChemSparseBuilder *b=NULL;ChemSparseMatrix m={0};CHECK(chem_sparse_builder_create(order,&b)==CHEM_OK);
    for(size_t i=0;i<n;++i)CHECK(chem_sparse_builder_add(b,entries[i].row,entries[i].column)==CHEM_OK);
    CHECK(chem_sparse_builder_blocks(b,blocks)==CHEM_OK);chem_sparse_builder_initial(b,value);
    CHECK(chem_sparse_builder_elements(b)==n*blocks);
    CHECK(chem_sparse_matrix_build(b,group,ordering,type,&m)==CHEM_OK);chem_sparse_builder_free(b);return m;
}
static void set(ChemSparseMatrix *m,size_t b,size_t r,size_t c,double v)
{CHECK(chem_sparse_matrix_set(m,b,r,c,v)==CHEM_OK);}
static size_t index_of(const ChemSparseMatrix *m,size_t b,size_t r,size_t c)
{size_t index=0;CHECK(chem_sparse_matrix_index(m,b,r,c,&index)==CHEM_OK);return index;}
static void zero(size_t group,ChemSparseOrdering order,int constant)
{
    ChemSparseMatrix m=create(3,1,group,order,CHEM_MATRIX_REAL,0,NULL,0);const ChemSparseMatrix *a=&m;
    CHECK(a->values.columns==0);size_t index;int is_zero;
    CHECK(chem_sparse_matrix_single_index(a,0,0,&index)==CHEM_NOT_FOUND);
    const size_t bad[][2]={{6,0},{1,3},{6,3}};
    for(size_t i=0;i<3;++i){CHECK(chem_sparse_matrix_single_index(a,bad[i][0],bad[i][1],&index)==CHEM_DIMENSION);
        CHECK(chem_sparse_matrix_is_zero(a,bad[i][0],bad[i][1],&is_zero)==CHEM_DIMENSION);}
    if(!constant){
        CHECK(chem_sparse_matrix_set(&m,0,0,4,2)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_set(&m,1,0,0,2)==CHEM_DIMENSION);
        CHECK(chem_sparse_matrix_set(&m,0,3,0,2)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_set(&m,0,1,1,2)==CHEM_NOT_FOUND);
    }
    chem_sparse_matrix_free(&m);
}
static void scalar(size_t group,ChemSparseOrdering order)
{
    ChemSparseMatrix m=create(3,1,group,order,CHEM_MATRIX_REAL,0,NULL,0);CHECK(chem_matrix_fill(&m.values,2)==CHEM_OK);
    for(size_t i=0;i<m.values.size;++i)CHECK(((double *)m.values.data)[i]==2);
    chem_sparse_matrix_free(&m);
    /* The reference empty-pattern loop is vacuous; also exercise a stored entry. */
    m=create(4,3,group,order,CHEM_MATRIX_REAL,5,pattern,0);CHECK(chem_matrix_fill(&m.values,2)==CHEM_OK);
    for(size_t i=0;i<m.values.size;++i)CHECK(((double *)m.values.data)[i]==2);
    chem_sparse_matrix_free(&m);
}
static ChemSparseMatrix diagonal_fixture(size_t group,ChemSparseOrdering order)
{
    ChemSparseBuilder *b=NULL;ChemSparseMatrix m={0};CHECK(chem_sparse_builder_create(4,&b)==CHEM_OK);
    for(size_t i=0;i<5;++i)CHECK(chem_sparse_builder_add(b,pattern[i].row,pattern[i].column)==CHEM_OK);
    CHECK(chem_sparse_builder_add(b,0,1)==CHEM_OK);CHECK(chem_sparse_builder_elements(b)==5);
    CHECK(chem_sparse_builder_blocks(b,3)==CHEM_OK);chem_sparse_builder_initial(b,24);
    CHECK(chem_sparse_matrix_build(b,group,order,CHEM_MATRIX_INT,&m)==CHEM_OK);chem_sparse_builder_free(b);
    set(&m,0,1,1,1);set(&m,1,1,1,3);set(&m,2,1,1,5);CHECK(chem_sparse_matrix_add_diagonal(&m,7)==CHEM_OK);return m;
}
static void diagonal(size_t group,ChemSparseOrdering order)
{
    ChemSparseMatrix m=diagonal_fixture(group,order);
    for(size_t b=0;b<3;++b){CHECK(chem_sparse_matrix_get(&m,b,1,1)==8+2*b);
        for(size_t i=0;i<5;++i)if(pattern[i].row!=pattern[i].column)CHECK(chem_sparse_matrix_get(&m,b,pattern[i].row,pattern[i].column)==24);}
    chem_sparse_matrix_free(&m);
}
static void single(size_t group,ChemSparseOrdering order,int constant)
{
    ChemSparseMatrix m=create(4,1,group,order,constant?CHEM_MATRIX_INT:CHEM_MATRIX_REAL,5,pattern,0),copy={0};
    if(constant){set(&m,0,2,1,45);set(&m,0,3,2,42);set(&m,0,2,3,21);CHECK(chem_sparse_matrix_copy(&m,&copy)==CHEM_OK);}
    const ChemSparseMatrix *a=constant?&copy:&m;
    size_t diag[4],count=0;CHECK(chem_sparse_matrix_diagonal_indices(a,0,diag,4,&count)==CHEM_OK);
    CHECK(count==1);CHECK(diag[0]==index_of(a,0,1,1));CHECK(a->values.columns==5);
    const size_t coords[][2]={{0,1},{3,2},{2,1},{2,2},{0,0},{3,3}};
    for(size_t i=0;i<6;++i){int is_zero=-1;CHECK(chem_sparse_matrix_is_zero(a,coords[i][0],coords[i][1],&is_zero)==CHEM_OK);CHECK(is_zero==(i>=3));}
    CHECK(chem_sparse_matrix_get(a,0,2,1)==(constant?45:0));
    if(!constant){set(&m,0,2,1,45);CHECK(chem_sparse_matrix_get(&m,0,2,1)==45);}
    size_t index;
    CHECK(chem_sparse_matrix_single_index(a,4,2,&index)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_single_index(a,1,5,&index)==CHEM_DIMENSION);
    CHECK(chem_sparse_matrix_index(a,1,0,2,&index)==CHEM_DIMENSION);
    CHECK(chem_sparse_matrix_single_index(a,2,2,&index)==CHEM_NOT_FOUND);CHECK(chem_sparse_matrix_single_index(a,0,2,&index)==CHEM_NOT_FOUND);
    if(!constant){
        CHECK(chem_sparse_matrix_set(&m,0,0,4,2)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_set(&m,1,0,0,2)==CHEM_DIMENSION);
        CHECK(chem_sparse_matrix_set(&m,0,5,0,2)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_set(&m,0,3,3,2)==CHEM_NOT_FOUND);
    }
    size_t width=group?group:1,i32=(order==CHEM_SPARSE_COLUMN?3:4)*width,i23=(order==CHEM_SPARSE_COLUMN?4:3)*width;
    CHECK(chem_sparse_matrix_single_index(a,3,2,&index)==CHEM_OK);CHECK(index==i32);
    CHECK(chem_sparse_matrix_single_index(a,2,3,&index)==CHEM_OK);CHECK(index==i23);
    if(constant){CHECK(((const int *)a->values.data)[i32]==42);CHECK(((const int *)a->values.data)[i23]==21);}
    else {((double *)m.values.data)[i32]=42;CHECK(((double *)m.values.data)[i32]==42);((double *)m.values.data)[i23]=21;CHECK(((double *)m.values.data)[i23]==21);}
    if(group){CHECK(a->values.group==group);CHECK(chem_matrix_group_size(&a->values)==5*group);CHECK(chem_matrix_group_count(&a->values)==1);}
    chem_sparse_matrix_free(&copy);chem_sparse_matrix_free(&m);
}
static void multi(size_t group,ChemSparseOrdering order)
{
    ChemSparseMatrix m=create(4,53,group,order,CHEM_MATRIX_REAL,5,pattern,24);CHECK(m.values.columns==5);
    CHECK(chem_sparse_matrix_get(&m,0,2,1)==24);set(&m,0,2,1,45);CHECK(chem_sparse_matrix_get(&m,0,2,1)==45);
    set(&m,43,3,2,29);CHECK(chem_sparse_matrix_get(&m,43,3,2)==29);set(&m,33,2,3,79);CHECK(chem_sparse_matrix_get(&m,33,2,3)==79);
    set(&m,2,2,3,64);CHECK(chem_sparse_matrix_get(&m,2,2,3)==64);
    for(size_t b=0;b<3;++b){size_t diag[4],count=0;CHECK(chem_sparse_matrix_diagonal_indices(&m,b,diag,4,&count)==CHEM_OK);CHECK(count==1);CHECK(diag[0]==index_of(&m,b,1,1));}
    size_t index;CHECK(chem_sparse_matrix_index(&m,0,4,2,&index)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_index(&m,2,1,5,&index)==CHEM_DIMENSION);
    CHECK(chem_sparse_matrix_index(&m,54,0,2,&index)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_index(&m,1,2,2,&index)==CHEM_NOT_FOUND);
    CHECK(chem_sparse_matrix_index(&m,2,0,2,&index)==CHEM_NOT_FOUND);CHECK(chem_sparse_matrix_single_index(&m,0,1,&index)==CHEM_MISSING_BLOCK);
    CHECK(chem_sparse_matrix_set(&m,0,0,4,2)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_set(&m,53,0,0,2)==CHEM_DIMENSION);
    CHECK(chem_sparse_matrix_set(&m,0,5,0,2)==CHEM_DIMENSION);CHECK(chem_sparse_matrix_set(&m,0,3,3,2)==CHEM_NOT_FOUND);
    size_t i23=(order==CHEM_SPARSE_COLUMN?4:3)*(group?group:1),i21=group?14:12;
    CHECK(index_of(&m,0,2,3)==i23);((double *)m.values.data)[i23]=21;CHECK(((double *)m.values.data)[i23]==21);
    CHECK(index_of(&m,2,2,1)==i21);((double *)m.values.data)[i21]=31;CHECK(((double *)m.values.data)[i21]==31);
    if(group){CHECK(m.values.group==2);CHECK(chem_matrix_group_size(&m.values)==10);
        ChemMatrix shape={.rows=4,.group=group};CHECK(chem_matrix_group_count(&shape)==2);}
    chem_sparse_matrix_free(&m);
}
static void print_sparse(size_t group,ChemSparseOrdering order,int nonzero)
{
    ChemSparseMatrix m;
    if(nonzero){m=create(4,3,group,order,CHEM_MATRIX_INT,5,pattern,32);set(&m,0,1,1,2);set(&m,1,1,1,4);set(&m,2,1,1,6);CHECK(chem_sparse_matrix_add_diagonal(&m,5)==CHEM_OK);}
    else m=diagonal_fixture(group,order);
    const char *expected=nonzero?
        "Block 0\n0, 1, 32\n1, 1, 7\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
        "Block 1\n0, 1, 32\n1, 1, 9\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
        "Block 2\n0, 1, 32\n1, 1, 11\n2, 1, 32\n2, 3, 32\n3, 2, 32\n":
        "Block 0\n0,24,0,0\n0,8,0,0\n0,24,0,24\n0,0,24,0\n"
        "Block 1\n0,24,0,0\n0,10,0,0\n0,24,0,24\n0,0,24,0\n"
        "Block 2\n0,24,0,0\n0,12,0,0\n0,24,0,24\n0,0,24,0\n";
    FILE *f=fopen("build/test-sparse-print.txt","w+");CHECK(f!=NULL);
    if(f){CHECK(chem_sparse_matrix_print(&m,f,nonzero)==CHEM_OK);rewind(f);char buffer[512]={0};
        CHECK(fread(buffer,1,sizeof(buffer)-1,f)==strlen(expected));CHECK(!strcmp(buffer,expected));fclose(f);remove("build/test-sparse-print.txt");}
    chem_sparse_matrix_free(&m);
}
static void bad_builder(void)
{
    ChemSparseBuilder *b=NULL;CHECK(chem_sparse_builder_create(3,&b)==CHEM_OK);
    CHECK(chem_sparse_builder_add(b,3,0)==CHEM_DIMENSION);CHECK(chem_sparse_builder_add(b,2,4)==CHEM_DIMENSION);CHECK(chem_sparse_builder_add(b,6,7)==CHEM_DIMENSION);
    chem_sparse_builder_free(b);
}
#define STANDARD(prefix,ordering) do { ChemSparseOrdering o=ordering; \
    RUN(prefix "ZeroMatrix",zero(0,o,0));RUN(prefix "ConstZeroMatrix",zero(0,o,1)); \
    RUN(prefix "SetScalar",scalar(0,o));RUN(prefix "AddToDiagonal",diagonal(0,o)); \
    RUN(prefix "SingleBlockMatrix",single(0,o,0));RUN(prefix "ConstSingleBlockMatrix",single(0,o,1)); \
    RUN(prefix "MultiBlockMatrix",multi(0,o));RUN(prefix "Print",print_sparse(0,o,0)); \
    RUN(prefix "PrintNonZero",print_sparse(0,o,1)); \
} while(0)
#define VECTOR(prefix,ordering) do { ChemSparseOrdering o=ordering; \
    RUN(prefix "ZeroMatrix",zero(2,o,0));RUN(prefix "ConstZeroMatrix",zero(3,o,1)); \
    RUN(prefix "SetScalar",scalar(3,o));RUN(prefix "AddToDiagonal",{diagonal(1,o);diagonal(3,o);diagonal(4,o);}); \
    RUN(prefix "SingleBlockMatrix",single(4,o,0));RUN(prefix "ConstSingleBlockMatrix",single(2,o,1)); \
    RUN(prefix "MultiBlockMatrix",multi(2,o));RUN(prefix "Print",{for(size_t g=1;g<=4;++g)print_sparse(g,o,0);}); \
    RUN(prefix "PrintNonZero",{for(size_t g=1;g<=4;++g)print_sparse(g,o,1);}); \
} while(0)
void test_sparse_matrix(void)
{
    STANDARD("sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.",CHEM_SPARSE_COLUMN);
    STANDARD("sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.",CHEM_SPARSE_ROW);
    VECTOR("sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.",CHEM_SPARSE_COLUMN);
    VECTOR("sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.",CHEM_SPARSE_ROW);
    RUN("sparse_matrix_standard_compressed_column_ordering:SparseMatrixBuilder.BadConfiguration",bad_builder());
    RUN("sparse_matrix_standard_compressed_row_ordering:SparseMatrixBuilder.BadConfiguration",bad_builder());
}
