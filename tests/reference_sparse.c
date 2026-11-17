#include "chem.h"
#include "chem_sparse.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define D(m,b,r,c) (*(double *)chem_sparse_at((m),(b),(r),(c)))
#define I(m,b,r,c) (*(int *)chem_sparse_at((m),(b),(r),(c)))
static int failures;
static void report(int grouped, int column, const char *name, int ok)
{
    printf("sparse_matrix_%s_compressed_%s_ordering:Sparse%sCompressed%sMatrix.%s %s %s\n",
           grouped ? "vector" : "standard",column ? "column" : "row",grouped ? "Vector" : "",
           column ? "Column" : "Row",name,ok ? "PASS" : "FAIL",ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static const ChemSparseCoordinate pattern[]={{0,1},{3,2},{0,1},{2,3},{1,1},{2,1}};
static int init(ChemSparseMatrix *m, size_t blocks, size_t group, int column, size_t element_size, const void *value)
{
    return chem_sparse_init(m,4,blocks,pattern,sizeof pattern/sizeof pattern[0],column ? CHEM_SPARSE_COLUMN : CHEM_SPARSE_ROW,
                            group,element_size,value)==CHEM_OK;
}
static int zero(size_t group, int column)
{
    ChemSparseMatrix m={0}; double z=0;
    int ok=chem_sparse_init(&m,3,1,NULL,0,column ? CHEM_SPARSE_COLUMN : CHEM_SPARSE_ROW,group,sizeof z,&z)==CHEM_OK;
    ok &= m.nonzeros==0 && chem_sparse_slot(&m,0,0,0)==SIZE_MAX && !chem_sparse_at(&m,0,0,0);
    chem_sparse_destroy(&m); return ok;
}
static int scalar(size_t group, int column)
{
    ChemSparseMatrix m={0}; double two=2;
    int ok=chem_sparse_init(&m,3,1,NULL,0,column ? CHEM_SPARSE_COLUMN : CHEM_SPARSE_ROW,group,sizeof two,NULL)==CHEM_OK &&
           chem_sparse_fill(&m,&two)==CHEM_OK;
    for (size_t i=0;i<m.storage.slots;++i) ok &= ((double *)(void *)m.storage.data)[i]==2;
    chem_sparse_destroy(&m); return ok;
}
static int diagonal(size_t group, int column)
{
    ChemSparseMatrix m={0}; int initial=24,inc=7;
    if (!init(&m,3,group,column,sizeof initial,&initial)) return 0;
    I(&m,0,1,1)=1; I(&m,1,1,1)=3; I(&m,2,1,1)=5;
    int ok=chem_sparse_add_diagonal(&m,&inc,CHEM_NUMBER_INT)==CHEM_OK;
    for (size_t b=0;b<3;++b) {
        ok &= I(&m,b,1,1)==(int)(8+2*b);
        ok &= I(&m,b,0,1)==24 && I(&m,b,2,1)==24 && I(&m,b,2,3)==24 && I(&m,b,3,2)==24;
    }
    chem_sparse_destroy(&m); return ok;
}
static int single(size_t group, int column, int use_const)
{
    ChemSparseMatrix m={0}; int initial=0;
    if (!init(&m,1,group,column,sizeof initial,&initial)) return 0;
    I(&m,0,2,1)=45; I(&m,0,3,2)=42; I(&m,0,2,3)=21;
    size_t scale=group ? group : 1;
    const ChemSparseMatrix *cm=&m;
    int ok=m.nonzeros==5 && chem_sparse_entry(cm,0,1)!=SIZE_MAX && chem_sparse_entry(cm,3,2)!=SIZE_MAX &&
        chem_sparse_entry(cm,2,1)!=SIZE_MAX && chem_sparse_entry(cm,2,2)==SIZE_MAX &&
        chem_sparse_entry(cm,0,0)==SIZE_MAX && chem_sparse_entry(cm,3,3)==SIZE_MAX &&
        *(const int *)chem_sparse_at_const(cm,0,2,1)==45 && chem_sparse_slot(cm,0,1,1)==scale &&
        chem_sparse_slot(cm,0,3,2)==(column ? 3 : 4)*scale && chem_sparse_slot(cm,0,2,3)==(column ? 4 : 3)*scale;
    if (use_const) ok &= *(const int *)chem_sparse_at_const(cm,0,3,2)==42 && *(const int *)chem_sparse_at_const(cm,0,2,3)==21;
    else {
        int *raw=(int *)(void *)m.storage.data;
        ok &= raw[chem_sparse_slot(cm,0,3,2)]==42 && raw[chem_sparse_slot(cm,0,2,3)]==21;
    }
    chem_sparse_destroy(&m); return ok;
}
static int multi(size_t group, int column)
{
    ChemSparseMatrix m={0}; double initial=24;
    if (!init(&m,53,group,column,sizeof initial,&initial)) return 0;
    D(&m,0,2,1)=45; D(&m,43,3,2)=29; D(&m,33,2,3)=79; D(&m,2,2,3)=64;
    int ok=m.nonzeros==5 && D(&m,0,2,1)==45 && D(&m,43,3,2)==29 && D(&m,33,2,3)==79 && D(&m,2,2,3)==64;
    size_t scale=group ? group : 1;
    ok &= chem_sparse_slot(&m,0,2,3)==(column ? 4 : 3)*scale;
    if (group==2) ok &= chem_sparse_slot(&m,2,2,1)==14;
    if (!group) ok &= chem_sparse_slot(&m,2,2,1)==12;
    for (size_t b=0;b<3;++b) ok &= chem_sparse_slot(&m,b,1,1)!=SIZE_MAX;
    ok &= chem_sparse_slot(&m,53,1,1)==SIZE_MAX && chem_sparse_slot(&m,0,2,2)==SIZE_MAX;
    chem_sparse_destroy(&m); return ok;
}
static int printing(size_t group, int column, int nonzero)
{
    ChemSparseMatrix m={0}; int initial=nonzero ? 32 : 24, inc=nonzero ? 5 : 7;
    if (!init(&m,3,group,column,sizeof initial,&initial)) return 0;
    I(&m,0,1,1)=nonzero ? 2 : 1; I(&m,1,1,1)=nonzero ? 4 : 3; I(&m,2,1,1)=nonzero ? 6 : 5;
    char buffer[1024]; int ok=chem_sparse_add_diagonal(&m,&inc,CHEM_NUMBER_INT)==CHEM_OK &&
        chem_sparse_format(&m,buffer,sizeof buffer,CHEM_NUMBER_INT,nonzero)==CHEM_OK;
    const char *all="Block 0\n0,24,0,0\n0,8,0,0\n0,24,0,24\n0,0,24,0\n"
                    "Block 1\n0,24,0,0\n0,10,0,0\n0,24,0,24\n0,0,24,0\n"
                    "Block 2\n0,24,0,0\n0,12,0,0\n0,24,0,24\n0,0,24,0\n";
    const char *nz="Block 0\n0, 1, 32\n1, 1, 7\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
                   "Block 1\n0, 1, 32\n1, 1, 9\n2, 1, 32\n2, 3, 32\n3, 2, 32\n"
                   "Block 2\n0, 1, 32\n1, 1, 11\n2, 1, 32\n2, 3, 32\n3, 2, 32\n";
    ok &= !strcmp(buffer,nonzero ? nz : all);
    chem_sparse_destroy(&m); return ok;
}
static int bad_config(int column)
{
    const ChemSparseCoordinate invalid[]={{3,0},{2,4},{6,7}}; double zero_value=0; int ok=1;
    for (size_t i=0;i<3;++i) {
        ChemSparseMatrix m={0};
        ok &= chem_sparse_init(&m,3,1,&invalid[i],1,column ? CHEM_SPARSE_COLUMN : CHEM_SPARSE_ROW,0,sizeof zero_value,&zero_value)==CHEM_INVALID;
        chem_sparse_destroy(&m);
    }
    return ok;
}
int main(void)
{
    for (int grouped=0;grouped<2;++grouped) for (int col=0;col<2;++col) {
        report(grouped,col,"ZeroMatrix",zero(grouped?2:0,col));
        report(grouped,col,"ConstZeroMatrix",zero(grouped?3:0,col));
        report(grouped,col,"SetScalar",scalar(grouped?3:0,col));
        int diag=diagonal(grouped?1:0,col); if (grouped) diag &= diagonal(3,col) && diagonal(4,col);
        report(grouped,col,"AddToDiagonal",diag);
        report(grouped,col,"SingleBlockMatrix",single(grouped?4:0,col,0));
        report(grouped,col,"ConstSingleBlockMatrix",single(grouped?2:0,col,1));
        report(grouped,col,"MultiBlockMatrix",multi(grouped?2:0,col));
        int all=1,nz=1;
        if (grouped) for (size_t g=1;g<=4;++g) { all &= printing(g,col,0); nz &= printing(g,col,1); }
        else { all=printing(0,col,0); nz=printing(0,col,1); }
        report(grouped,col,"Print",all); report(grouped,col,"PrintNonZero",nz);
        if (!grouped) {
            printf("sparse_matrix_standard_compressed_%s_ordering:SparseMatrixBuilder.BadConfiguration %s %s\n",
                   col?"column":"row",bad_config(col)?"PASS":"FAIL","reference quantities matched");
            failures+=!bad_config(col);
        }
    }
    return failures ? 1 : 0;
}
