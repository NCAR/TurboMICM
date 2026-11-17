#include "sparse.h"
#include <stdio.h>
#include <string.h>

static int fill_pattern(ChemSparseBuilder *builder)
{
    return chem_sparse_builder_add(builder,0,1) || chem_sparse_builder_add(builder,3,2) ||
           chem_sparse_builder_add(builder,0,1) || chem_sparse_builder_add(builder,2,3) ||
           chem_sparse_builder_add(builder,1,1) || chem_sparse_builder_add(builder,2,1);
}
static int check_print(ChemSparseMatrix *m,int nonzero)
{
    static const char *dense="Block 0\n0,24,0,0\n0,8,0,0\n0,24,0,24\n0,0,24,0\n"
        "Block 1\n0,24,0,0\n0,10,0,0\n0,24,0,24\n0,0,24,0\n"
        "Block 2\n0,24,0,0\n0,12,0,0\n0,24,0,24\n0,0,24,0\n";
    static const char *nz="Block 0\n0, 1, 32\n1, 1, 7\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
        "Block 1\n0, 1, 32\n1, 1, 9\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
        "Block 2\n0, 1, 32\n1, 1, 11\n2, 1, 32\n2, 3, 32\n3, 2, 32\n";
    FILE *f=tmpfile(); if(!f) return 0;
    int ok=!chem_sparse_print(m,f,nonzero) && !fseek(f,0,SEEK_SET);
    char text[1024]; size_t n=fread(text,1,sizeof(text)-1,f); text[n]=0;
    ok=ok && !strcmp(text,nonzero ? nz : dense); fclose(f); return ok;
}
static int sparse_case(int kind,size_t width,int column_major)
{
    ChemSparseBuilder b={0}; ChemSparseMatrix m={0}; size_t i=0; int zero=0,ok=1;
    if (chem_sparse_builder_init(&b,kind<=2 ? 3 : 4)) return 0;
    if (kind>=3 && fill_pattern(&b)) { chem_sparse_builder_free(&b); return 0; }
    size_t blocks=kind==6 ? 53 : kind==3 || kind==7 || kind==8 ? 3 : 1;
    double initial=kind==6 || kind==3 || kind==7 ? 24 : kind==8 ? 32 : 0;
    if (chem_sparse_init(&m,&b,blocks,width,column_major,initial)) { chem_sparse_builder_free(&b); return 0; }
    if (kind<=2) {
        if (b.count!=0 || m.nonzeros!=0 || !chem_sparse_index_single(&m,0,0,&i) || (chem_sparse_is_zero(&m,0,0,&zero) || !zero)) ok=0;
        if (kind==2) { chem_sparse_fill(&m,2); for(size_t k=0;k<m.storage_size;++k) if(m.data[k]!=2) ok=0; }
    } else if (kind==3 || kind==7 || kind==8) {
        for(size_t block=0;block<3;++block) *chem_sparse_at(&m,block,1,1)=kind==8 ? 2+2*block : 1+2*block;
        chem_sparse_add_diagonal(&m,kind==8 ? 5 : 7);
        for(size_t block=0;block<3;++block) {
            if (*chem_sparse_at(&m,block,1,1)!=(kind==8 ? 7+2*block : 8+2*block)) ok=0;
            if (kind==3 && (*chem_sparse_at(&m,block,0,1)!=24 || *chem_sparse_at(&m,block,2,1)!=24 || *chem_sparse_at(&m,block,2,3)!=24 || *chem_sparse_at(&m,block,3,2)!=24)) ok=0;
        }
        if ((kind==7 || kind==8) && !check_print(&m,kind==8)) ok=0;
    } else {
        if (b.count!=5 || m.nonzeros!=5) ok=0;
        size_t points[][2]={{0,1},{3,2},{2,1},{2,2},{0,0},{3,3}};
        for(size_t k=0;k<6;++k) if (chem_sparse_is_zero(&m,points[k][0],points[k][1],&zero) || zero!=(k>=3)) ok=0;
        if (!chem_sparse_at(&m,0,2,1) || *chem_sparse_at(&m,0,2,1)!=(kind==6 ? 24 : 0)) ok=0;
        *chem_sparse_at(&m,0,2,1)=45;
        if (*chem_sparse_at(&m,0,2,1)!=45) ok=0;
        size_t diag,diagonal[4],diagonal_count;
        if (chem_sparse_index(&m,0,1,1,&diag) || chem_sparse_diagonal_indices(&m,0,diagonal,4,&diagonal_count) || diagonal_count!=1 || diagonal[0]!=diag) ok=0;
        if (kind==5) {
            *chem_sparse_at(&m,0,3,2)=42; *chem_sparse_at(&m,0,2,3)=21;
            const ChemSparseMatrix *const_matrix=&m;
            if (*chem_sparse_const_at(const_matrix,0,3,2)!=42 || *chem_sparse_const_at(const_matrix,0,2,3)!=21 || *chem_sparse_const_at(const_matrix,0,2,1)!=45) ok=0;
        }
        if (!chem_sparse_index(&m,0,2,2,&i) || !chem_sparse_index(&m,0,4,2,&i)) ok=0;
        if (kind==6) {
            for (size_t block=0;block<3;++block) {
                if (chem_sparse_index(&m,block,1,1,&diag) || chem_sparse_diagonal_indices(&m,block,diagonal,4,&diagonal_count) || diagonal_count!=1 || diagonal[0]!=diag) ok=0;
            }
            *chem_sparse_at(&m,43,3,2)=29; *chem_sparse_at(&m,33,2,3)=79; *chem_sparse_at(&m,2,2,3)=64;
            if (*chem_sparse_at(&m,43,3,2)!=29 || *chem_sparse_at(&m,33,2,3)!=79 || *chem_sparse_at(&m,2,2,3)!=64 || !chem_sparse_index_single(&m,0,1,&i)) ok=0;
            if (chem_sparse_index(&m,0,2,3,&i) || i!=(column_major ? 4 : 3)*width) ok=0;
            if (chem_sparse_index(&m,2,2,1,&i) || i!=(2/width)*5*width+2*width+2%width) ok=0;
        } else {
            if (chem_sparse_index_single(&m,3,2,&i) || i!=(column_major ? 3 : 4)*width) ok=0;
            m.data[i]=42; if (m.data[i]!=42) ok=0;
            if (chem_sparse_index_single(&m,2,3,&i) || i!=(column_major ? 4 : 3)*width) ok=0;
            m.data[i]=21; if (m.data[i]!=21) ok=0;
        }
    }
    chem_sparse_free(&m); chem_sparse_builder_free(&b); return ok;
}
static int all_widths(int kind,int grouped,int csc)
{
    if (!grouped) return sparse_case(kind,1,csc);
    for(size_t width=1;width<=4;++width) if (!sparse_case(kind,width,csc)) return 0;
    return 1;
}
static int bad_builder(void)
{
    ChemSparseBuilder b={0}; if (chem_sparse_builder_init(&b,3)) return 0;
    int ok=chem_sparse_builder_add(&b,3,0) && chem_sparse_builder_add(&b,2,4) && chem_sparse_builder_add(&b,6,7);
    chem_sparse_builder_free(&b); return ok;
}
int main(void)
{
    const char *prefix[]={"sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix", "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix", "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix", "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix"};
    const char *names[]={"ZeroMatrix","ConstZeroMatrix","SetScalar","AddToDiagonal","SingleBlockMatrix","ConstSingleBlockMatrix","MultiBlockMatrix","Print","PrintNonZero"};
    for(int layout=0;layout<4;++layout) for(int kind=0;kind<9;++kind)
        printf("%s.%s %d\n",prefix[layout],names[kind],all_widths(kind,layout>=2,layout%2));
    printf("sparse_matrix_standard_compressed_row_ordering:SparseMatrixBuilder.BadConfiguration %d\n",bad_builder());
    printf("sparse_matrix_standard_compressed_column_ordering:SparseMatrixBuilder.BadConfiguration %d\n",bad_builder());
    return 0;
}
