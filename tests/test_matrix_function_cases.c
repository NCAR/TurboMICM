#include "matrix.h"
#include <stdio.h>

typedef struct { int mode; } Context;
static int row_function(size_t row, ChemMatrix *const *m, ChemVector *const *v, void *ctx)
{
    int mode=((Context *)ctx)->mode;
    double *a=chem_matrix_at(m[0],row,0), *b=chem_matrix_at(m[0],row,1);
    if (mode==0) return 0;
    if (!a || !b) return 1;
    if (mode==1) {
        double *c=chem_matrix_at(m[0],row,2); if (!c) return 1;
        *c=4*(*a+*b+*c);
    } else if (mode==2 || mode==3) {
        const double *other=chem_matrix_const_at(m[1],row,2); if (!other) return 1;
        *b=(*a+*other)*(mode==3 ? 2 : 1);
    } else if (mode==4) {
        *b=(*a+((double *)v[0]->data)[row])*3;
    } else if (mode==5) {
        *a=((double *)v[0]->data)[row]+((double *)v[1]->data)[row];
    } else if (mode==6) {
        double *other_a=chem_matrix_at(m[1],row,0), *other_b=chem_matrix_at(m[1],row,1);
        if (!other_a || !other_b) return 1;
        double value=((double *)v[0]->data)[row]; *b=*a+value; *other_b=*other_a+value;
    } else if (mode==7) {
        *a=((double *)v[0]->data)[row]*2;
    } else if (mode==8) {
        ((double *)v[0]->data)[row]=*a*3;
    } else if (mode==9) {
        *a=((double *)v[0]->data)[row]*10;
    } else if (mode==10) {
        if (!chem_matrix_at(m[0],row,5)) return 1;
    }
    return 0;
}

static int function_case(int kind,size_t width,int grouped)
{
    ChemMatrix p[2]={{0}},m[2]={{0}},m2[2]={{0}}; ChemMatrixFunction fn={0}; Context ctx={0};
    ChemMatrix *protos[]={&p[0],&p[1]},*matrices[]={&m[0],&m[1]},*matrices2[]={&m2[0],&m2[1]};
    double d1[10]={0},d2[10]={0}; ChemVector vec1={d1,0},vec2={d2,0},*vectors[]={&vec1,&vec2};
    size_t pc=1,vc=0,pr=3,pc0=3,pc1=3; int ok=1;
    if (kind==0) { pr=5; ctx.mode=1; }
    if (kind==1) { pc=2; pc0=2; ctx.mode=2; }
    if (kind==2) { pc=2; pc0=3; pc1=3; }
    if (kind==3 || kind==4) { pc0=4; ctx.mode=kind==3 ? 10 : 0; }
    if (kind==5) { pc=2; pc0=2; ctx.mode=3; }
    if (kind==6) { vc=1; ctx.mode=4; }
    if (kind==7 || kind==8) { pc0=2; vc=kind==7 ? 1 : 0; pc=kind==8 ? 2 : 1; }
    if (kind==9) { pr=4; pc0=3; }
    if (kind>=10 && kind<=15) { pr=5; vc=kind==15 ? 2 : 1; }
    if (kind==16) { pr=5; vc=2; ctx.mode=5; }
    if (kind==17) { pr=4; pc=2; pc0=2; vc=1; ctx.mode=6; }
    if (kind==18 || kind==19) { pc=2; pc0=2; vc=1; }
    if (kind>=20 && kind<=23) { pc0=2; vc=1; ctx.mode=kind==20 ? 7 : kind==21 ? 8 : kind==22 ? 9 : 0; }
    if (chem_matrix_init(&p[0],pr,pc0,width,grouped,0) || (pc==2 && chem_matrix_init(&p[1],kind==2 ? 4 : pr,pc1,width,grouped,0)) ||
        chem_matrix_function_init(&fn,protos,pc,vc,row_function,&ctx)) return 0;
    if (kind==0) {
        chem_matrix_init(&m[0],5,3,width,grouped,-1);
        for(int r=0;r<5;++r) for(int c=0;c<3;++c) *chem_matrix_at(&m[0],r,c)=r-2+10*c;
        if (chem_matrix_function_call(&fn,matrices,1,NULL,0)) ok=0;
        for(int r=0;r<5;++r) if (*chem_matrix_at(&m[0],r,0)!=r-2 || *chem_matrix_at(&m[0],r,1)!=r+8 || *chem_matrix_at(&m[0],r,2)!=4*(3*r+24)) ok=0;
        chem_matrix_init(&m2[0],3,3,width,grouped,-1);
        if (chem_matrix_function_call(&fn,matrices2,1,NULL,0)) ok=0;
        for(int r=0;r<3;++r) if (*chem_matrix_at(&m2[0],r,0)!=-1 || *chem_matrix_at(&m2[0],r,1)!=-1 || *chem_matrix_at(&m2[0],r,2)!=-12) ok=0;
    } else if (kind==1 || kind==5) {
        int rows=kind==1 ? 3 : 5;
        chem_matrix_init(&m[0],rows,2,width,grouped,0); chem_matrix_init(&m[1],rows,3,width,grouped,0);
        for(int r=0;r<rows;++r) { *chem_matrix_at(&m[0],r,0)=kind==1 ? r : r+1; *chem_matrix_at(&m[1],r,0)=r*2; *chem_matrix_at(&m[1],r,2)=kind==1 ? r*4 : r*10; }
        if (chem_matrix_function_call(&fn,matrices,2,NULL,0)) ok=0;
        for(int r=0;r<rows;++r) if (*chem_matrix_at(&m[0],r,1)!=(kind==1 ? r*5 : (r+1+r*10)*2) || *chem_matrix_at(&m[0],r,0)!=(kind==1 ? r : r+1) || *chem_matrix_at(&m[1],r,0)!=r*2) ok=0;
        if (kind==5) {
            chem_matrix_init(&m2[0],2,2,width,grouped,0); chem_matrix_init(&m2[1],2,3,width,grouped,0);
            *chem_matrix_at(&m2[0],0,0)=10; *chem_matrix_at(&m2[0],1,0)=20; *chem_matrix_at(&m2[1],0,2)=5; *chem_matrix_at(&m2[1],1,2)=15;
            if (chem_matrix_function_call(&fn,matrices2,2,NULL,0) || *chem_matrix_at(&m2[0],0,1)!=30 || *chem_matrix_at(&m2[0],1,1)!=70) ok=0;
        }
    } else if (kind==2 || kind==8 || kind==18) {
        size_t ar=kind==18 ? 4 : kind==2 ? 3 : 5, br=kind==18 ? 5 : kind==2 ? 4 : 3;
        chem_matrix_init(&m[0],ar,pc0,width,grouped,0); chem_matrix_init(&m[1],br,pc1,width,grouped,0);
        vec1.size=ar;
        if (!chem_matrix_function_call(&fn,matrices,2,vc ? vectors : NULL,vc)) ok=0;
        if (kind==2) {
            chem_matrix_init(&m2[0],5,3,width,grouped,0); chem_matrix_init(&m2[1],5,3,width,grouped,1);
            if (chem_matrix_function_call(&fn,matrices2,2,NULL,0)) ok=0;
        }
    } else if (kind==3) {
        chem_matrix_init(&m[0],3,4,width,grouped,1);
        if (!chem_matrix_function_call(&fn,matrices,1,NULL,0) || chem_matrix_at(&m[0],0,5)) ok=0;
    } else if (kind==4 || kind==9) {
        chem_matrix_init(&m[0],kind==4 ? 3 : 4,pc0,width,grouped,1); chem_matrix_init(&m2[0],kind==4 ? 3 : 4,pc0+1,width,grouped,0);
        if (chem_matrix_function_call(&fn,matrices,1,NULL,0) || !chem_matrix_function_call(&fn,matrices2,1,NULL,0)) ok=0;
        chem_matrix_free(&m2[0]); chem_matrix_init(&m2[0],7,pc0,width,grouped,0);
        if (chem_matrix_function_call(&fn,matrices2,1,NULL,0)) ok=0;
        if (kind==9) { chem_matrix_free(&m2[0]); chem_matrix_init(&m2[0],4,2,width,grouped,0); if (!chem_matrix_function_call(&fn,matrices2,1,NULL,0)) ok=0; }
    } else if (kind==6) {
        chem_matrix_init(&m[0],5,3,width,grouped,0); vec1.size=5;
        for(int r=0;r<5;++r) { *chem_matrix_at(&m[0],r,0)=r+1; d1[r]=r*10; }
        if (chem_matrix_function_call(&fn,matrices,1,vectors,1)) ok=0;
        for(int r=0;r<5;++r) if (*chem_matrix_at(&m[0],r,1)!=(r+1+r*10)*3) ok=0;
    } else if (kind==7) {
        chem_matrix_init(&m[0],5,2,width,grouped,0); vec1.size=3;
        if (!chem_matrix_function_call(&fn,matrices,1,vectors,1)) ok=0;
        chem_matrix_init(&m2[0],3,2,width,grouped,0); vec1.size=5;
        if (!chem_matrix_function_call(&fn,matrices2,1,vectors,1)) ok=0;
    } else if (kind>=10 && kind<=15) {
        size_t rows=kind==13 || kind==14 ? 0 : 5;
        chem_matrix_init(&m[0],rows,3,width,grouped,1);
        vec1.size=kind==10 ? 3 : kind==11 ? 10 : kind==12 || kind==14 ? 0 : 5;
        vec2.size=3;
        int status=chem_matrix_function_call(&fn,matrices,1,vectors,vc);
        if (kind==14 ? status!=0 : status==0) ok=0;
    } else if (kind==16) {
        chem_matrix_init(&m[0],5,3,width,grouped,0); vec1.size=vec2.size=5;
        for(int r=0;r<5;++r) { d1[r]=r+1; d2[r]=(r+1)*10; }
        if (chem_matrix_function_call(&fn,matrices,1,vectors,2)) ok=0;
        for(int r=0;r<5;++r) if (*chem_matrix_at(&m[0],r,0)!=(r+1)*11) ok=0;
    } else if (kind==17) {
        chem_matrix_init(&m[0],4,2,width,grouped,0); chem_matrix_init(&m[1],4,3,width,grouped,0); vec1.size=4;
        for(int r=0;r<4;++r) { d1[r]=r*2; *chem_matrix_at(&m[0],r,0)=r+1; *chem_matrix_at(&m[1],r,0)=r*10; }
        if (chem_matrix_function_call(&fn,matrices,2,vectors,1)) ok=0;
        for(int r=0;r<4;++r) if (*chem_matrix_at(&m[0],r,1)!=r+1+r*2 || *chem_matrix_at(&m[1],r,1)!=r*12) ok=0;
    } else if (kind==19) {
        chem_matrix_init(&m[0],5,2,width,grouped,0); chem_matrix_init(&m[1],5,3,width,grouped,0); vec1.size=4;
        if (!chem_matrix_function_call(&fn,matrices,2,vectors,1)) ok=0;
    } else if (kind>=20 && kind<=23) {
        chem_matrix_init(&m[0],3,2,width,grouped,0); vec1.size=3;
        for(int r=0;r<3;++r) { d1[r]=kind==20 ? (r+1)*10 : kind==21 ? (r+1)*5 : r+1; *chem_matrix_at(&m[0],r,0)=r+1; }
        if (chem_matrix_function_call(&fn,matrices,1,vectors,1)) ok=0;
        if (kind==20) for(int r=0;r<3;++r) if (*chem_matrix_at(&m[0],r,0)!=(r+1)*20) ok=0;
        if (kind==21) for(int r=0;r<3;++r) if (d1[r]!=(r+1)*3) ok=0;
        if (kind==22) { for(int r=0;r<3;++r) if (*chem_matrix_at(&m[0],r,0)!=(r+1)*10) ok=0;
            for(int r=0;r<3;++r) d1[r]=r+5;
            if (chem_matrix_function_call(&fn,matrices,1,vectors,1)) ok=0;
            for(int r=0;r<3;++r) if (*chem_matrix_at(&m[0],r,0)!=(r+5)*10) ok=0;
        }
    }
    chem_matrix_function_free(&fn); for(int i=0;i<2;++i) { chem_matrix_free(&p[i]); chem_matrix_free(&m[i]); chem_matrix_free(&m2[i]); }
    return ok;
}

static int all_widths(int kind,int grouped)
{
    if (!grouped) return function_case(kind,1,0);
    for(size_t width=1;width<=4;++width) if (!function_case(kind,width,1)) return 0;
    return 1;
}
int main(void)
{
    const char *names[]={"ArrayFunction","MultiMatrixArrayFunction","MismatchedRowDimensions","MismatchedColumnDimensions","WrongMatrixDimensions","MultiMatrixDifferentRowsFromCreation","MatrixVectorDifferentRowsFromCreation","MismatchedRowsAtInvocation","MultipleMatricesMismatchedRowsAtInvocation","WrongColumnCountAtInvocation","VectorTooSmall","VectorTooLarge","EmptyVectorNonEmptyMatrix","NonEmptyVectorEmptyMatrix","EmptyVectorEmptyMatrix","MultipleVectorsDifferentSizes","MultipleVectorsSameSize","MultipleMatricesOneVector","MultipleMatricesDifferentRowsVector","VectorSizeMatchesOneMatrixOnly","ConstVector","MutableVector","FunctionReusabilityWithVectors","FunctionInvocationWithWrongSizedVector"};
    for(int i=0;i<24;++i) {
        printf("matrix:Matrix.%s %d\n",names[i],all_widths(i,0));
        printf("vector_matrix:VectorMatrix.%s %d\n",names[i],all_widths(i,1));
    }
    return 0;
}
