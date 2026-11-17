#include "matrix.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

static void sum_square(double v,void *x) { *(double *)x+=v*v; }
static void maximum(double v,void *x) { if(v>*(double *)x) *(double *)x=v; }
static void any_large(double v,void *x) { *(int *)x=*(int *)x || v>10; }
static void all_finite(double v,void *x) { *(int *)x=*(int *)x && isfinite(v); }
static void count_row(double v,void *x) { (void)v; *(double *)x+=1; }

typedef struct { int mode; int *indices; int *values; } Context;
static int row_fn(size_t r,ChemMatrix *const *m,ChemVector *const *v,void *opaque)
{
    Context *ctx=opaque; ChemMatrix *a=m[0];
    if (ctx->mode==4) return 0;
    double *c0=chem_matrix_at(a,r,0), *c1=chem_matrix_at(a,r,1);
    if (!c0 || !c1) return 1;
    if (ctx->mode==0) {
        double *c2=chem_matrix_at(a,r,2),*c3=chem_matrix_at(a,r,3),*c4=chem_matrix_at(a,r,4);
        if (!c2 || !c3 || !c4) return 1;
        double t1=*c0**c1,t2=*c0+*c1; *c2=t1+t2; *c3=t1-t2; *c4=t1*t2;
    } else if(ctx->mode==1) {
        double *c2=chem_matrix_at(a,r,2),*c3=chem_matrix_at(a,r,3); if (!c2 || !c3) return 1;
        *c1=*c0*2; *c2=*c0+*c1; *c3=*c2**c1;
    } else if(ctx->mode==2 || ctx->mode==5 || ctx->mode==6) {
        double *c2=chem_matrix_at(a,r,2); if (!c2) return 1;
        double sum=*c0+*c1+*c2;
        if (ctx->mode==2) *c2=2*sum;
        else if (ctx->mode==5) *c2=((double *)v[0]->data)[r]*sum;
        else *c2=*c0+*c1+((double *)v[0]->data)[r];
    } else if(ctx->mode==3) {
        const double *c2=chem_matrix_const_at(a,r,2); if (!c2) return 1;
        double temp=*c0+*c1+*c2; (void)temp;
    } else if(ctx->mode==7) {
        *c0=((double *)v[0]->data)[r]*2;
    } else if(ctx->mode==8) {
        for(int i=0;i<3;++i) *c0+=ctx->values[ctx->indices[i]];
        *c1+=ctx->values[ctx->indices[1]];
    }
    return 0;
}

static int test_case(int kind,size_t width,int grouped)
{
    ChemMatrix m={0}, second={0}, third={0}; ChemMatrixFunction fn={0}; Context ctx={0};
    ChemMatrix *matrices[]={&m}; double data[]={100,200,300,400,500}; ChemVector vector={data,5},*vectors[]={&vector};
    size_t rows=kind==0 ? 4 : kind==1 ? 3 : kind==2 ? 2 : kind==3 ? 3 : kind==4 ? 0 : kind==5 ? 5 : kind==6 ? 4 : kind==7 ? 3 : kind==8 || kind==9 ? 3 : kind==15 ? 3 : 3;
    size_t cols=kind==0 ? 5 : kind==1 || kind==3 ? 4 : kind==7 ? 2 : kind==8 ? 2 : kind==16 ? 1 : kind>=10 && kind<=15 ? 2 : 3;
    if (chem_matrix_init(&m,rows,cols,width,grouped,0)) return 0;
    int ok=1;
    if (kind<=9) {
        ctx.mode=kind;
        if (kind==8) ctx.mode=7; /* const-signature callback reads the vector */
        if (kind==9) ctx.mode=8; /* vector capture */
        int indices[]={2,0,1},values[]={10,20,30}; ctx.indices=indices;ctx.values=values;
        size_t vc=(kind==5 || kind==6 || kind==7 || kind==8) ? 1 : 0;
        vector.size=rows;
        if (kind==0 || kind==1) for(size_t r=0;r<rows;++r) { *chem_matrix_at(&m,r,0)=r+1; if(kind==0) *chem_matrix_at(&m,r,1)=(r+1)*10; }
        if (kind==2) for(size_t r=0;r<rows;++r) for(size_t c=0;c<cols;++c) *chem_matrix_at(&m,r,c)=r+c;
        if (kind==3) for(size_t r=0;r<rows;++r) for(size_t c=0;c<cols;++c) *chem_matrix_at(&m,r,c)=r*10+c;
        if (kind==5) for(size_t r=0;r<rows;++r) for(size_t c=0;c<cols;++c) *chem_matrix_at(&m,r,c)=(double)r-2+10*c;
        if (kind==6) for(size_t r=0;r<rows;++r) { *chem_matrix_at(&m,r,0)=r+1; *chem_matrix_at(&m,r,1)=(r+1)*10; data[r]=(r+1)*100; }
        if (chem_matrix_function_init(&fn,matrices,1,vc,row_fn,&ctx) || chem_matrix_function_call(&fn,matrices,1,vc ? vectors : NULL,vc)) ok=0;
        if (kind==0) { int checks[]={0,1,3}; for(int i=0;i<3;++i) { int r=checks[i],x=r+1; if(*chem_matrix_at(&m,r,2)!=10*x*x+11*x || *chem_matrix_at(&m,r,3)!=10*x*x-11*x || *chem_matrix_at(&m,r,4)!=110*x*x*x) ok=0; } }
        if (kind==1) for(int r=0;r<3;++r) { double x=r+1; if(*chem_matrix_at(&m,r,1)!=2*x || *chem_matrix_at(&m,r,2)!=3*x || *chem_matrix_at(&m,r,3)!=6*x*x) ok=0; }
        if (kind==2) {
            if (*chem_matrix_at(&m,0,2)!=6 || *chem_matrix_at(&m,1,2)!=12) ok=0;
            chem_matrix_init(&second,2,3,width,grouped,5); matrices[0]=&second;
            if (chem_matrix_function_call(&fn,matrices,1,NULL,0) || *chem_matrix_at(&second,0,2)!=30 || *chem_matrix_at(&second,1,2)!=30) ok=0;
            chem_matrix_init(&third,2,3,width,grouped,0); *chem_matrix_at(&third,1,0)=10; matrices[0]=&third;
            if (chem_matrix_function_call(&fn,matrices,1,NULL,0) || *chem_matrix_at(&third,0,2)!=0 || *chem_matrix_at(&third,1,2)!=20) ok=0;
        }
        if (kind==3 && (*chem_matrix_at(&m,0,0)!=0 || *chem_matrix_at(&m,1,2)!=12)) ok=0;
        if (kind==4) {
            chem_matrix_function_free(&fn); chem_matrix_free(&m); chem_matrix_init(&m,3,0,width,grouped,1); matrices[0]=&m;
            if (chem_matrix_function_init(&fn,matrices,1,0,row_fn,&ctx) || chem_matrix_function_call(&fn,matrices,1,NULL,0)) ok=0;
        }
        if (kind==5) for(int r=0;r<5;++r) if (*chem_matrix_at(&m,r,2)!=(r+1)*100*(3*r+24) || *chem_matrix_at(&m,r,0)!=r-2 || *chem_matrix_at(&m,r,1)!=r+8) ok=0;
        if (kind==6) for(int r=0;r<4;++r) if (*chem_matrix_at(&m,r,2)!=(r+1)*111) ok=0;
        if (kind==7 || kind==8) for(int r=0;r<3;++r) if (*chem_matrix_at(&m,r,0)!=data[r]*2) ok=0;
        if (kind==9) for(int r=0;r<3;++r) if (*chem_matrix_at(&m,r,0)!=60 || *chem_matrix_at(&m,r,1)!=10 || *chem_matrix_at(&m,r,2)!=0) ok=0;
    } else if (kind==10) {
        if (chem_matrix_fill_column(&m,1,3.2)) ok=0;
        for(int r=0;r<3;++r) if (*chem_matrix_at(&m,r,0)!=0 || *chem_matrix_at(&m,r,1)!=3.2) ok=0;
        double vec[3]; for(int i=0;i<3;++i) vec[i]=3.2;
        if (vec[0]!=3.2 || vec[1]!=3.2 || vec[2]!=3.2 || chem_matrix_fill_column(&m,0,9.9)) ok=0;
        for(int r=0;r<3;++r) if (*chem_matrix_at(&m,r,0)!=9.9) ok=0;
    } else if (kind==11) {
        double vec[]={3.2,4.2,1.3}, out[3]={0}, temp[3]={0};
        if (chem_matrix_vector_to_column(&m,1,vec,3)) ok=0;
        for(int r=0;r<3;++r) if (*chem_matrix_at(&m,r,0)!=0 || *chem_matrix_at(&m,r,1)!=vec[r]) ok=0;
        if (chem_matrix_column_to_vector(&m,1,out,3)) ok=0;
        for(int r=0;r<3;++r) if (out[r]!=vec[r]) ok=0;
        if (chem_matrix_copy_column(&m,0,&m,1)) ok=0;
        for(int r=0;r<3;++r) if (*chem_matrix_at(&m,r,0)!=vec[r]) ok=0;
        chem_matrix_fill_column(&m,0,0); if (chem_matrix_copy_column(&m,0,&m,1)) ok=0;
        if (chem_matrix_column_to_vector(&m,1,temp,3) || chem_matrix_fill_column(&m,1,0) || chem_matrix_vector_to_column(&m,1,temp,3)) ok=0;
        for(int r=0;r<3;++r) if (*chem_matrix_at(&m,r,0)!=vec[r] || *chem_matrix_at(&m,r,1)!=vec[r]) ok=0;
    } else if (kind>=12 && kind<=16) {
        if (kind==12) { double values[]={1,2,3,4,5,6},result=0; for(int r=0;r<3;++r) for(int c=0;c<2;++c) *chem_matrix_at(&m,r,c)=values[r*2+c];
            if (chem_matrix_reduce_column(&m,0,0,sum_square,&result) || chem_matrix_reduce_column(&m,1,0,sum_square,&result) || result!=91) ok=0;
        } else if (kind==13) { double values[]={1,8,5,2,3,4},result=-DBL_MAX; for(int r=0;r<3;++r) for(int c=0;c<2;++c) *chem_matrix_at(&m,r,c)=values[r*2+c];
            if (chem_matrix_reduce_column(&m,0,0,maximum,&result) || chem_matrix_reduce_column(&m,1,0,maximum,&result) || result!=8) ok=0;
        } else if (kind==14) { int result=0; chem_matrix_fill(&m,1);
            if (chem_matrix_reduce_column(&m,0,0,any_large,&result) || chem_matrix_reduce_column(&m,1,0,any_large,&result) || result) ok=0;
            *chem_matrix_at(&m,2,1)=42; if (chem_matrix_reduce_column(&m,0,0,any_large,&result) || chem_matrix_reduce_column(&m,1,0,any_large,&result) || !result) ok=0;
        } else if (kind==15) { int result=1; chem_matrix_fill(&m,1);
            if (chem_matrix_reduce_column(&m,0,0,all_finite,&result) || chem_matrix_reduce_column(&m,1,0,all_finite,&result) || !result) ok=0;
            *chem_matrix_at(&m,1,0)=NAN; if (chem_matrix_reduce_column(&m,0,0,all_finite,&result) || chem_matrix_reduce_column(&m,1,0,all_finite,&result) || result) ok=0;
        } else { double result=0; if (chem_matrix_reduce_column(&m,0,1,count_row,&result) || result!=3) ok=0; }
    }
    chem_matrix_function_free(&fn); chem_matrix_free(&m); chem_matrix_free(&second); chem_matrix_free(&third); return ok;
}

static int all_widths(int kind,int grouped)
{
    if (!grouped) return test_case(kind,1,0);
    for(size_t width=1;width<=4;++width) if (!test_case(kind,width,1)) return 0;
    return 1;
}
static void double_int(int *out,int input,void *context) { (void)context; *out=input*2; }
static int integer_case(size_t width)
{
    ChemIntMatrix m={0}; int vector[]={10,20,30};
    if (chem_int_matrix_init(&m,3,2,width,0)) return 0;
    int ok=!chem_int_matrix_apply_vector(&m,0,vector,3,double_int,NULL) &&
        *chem_int_matrix_at(&m,0,0)==20 && *chem_int_matrix_at(&m,1,0)==40 && *chem_int_matrix_at(&m,2,0)==60;
    chem_int_matrix_free(&m); return ok;
}
int main(void)
{
    const char *names[]={"MultipleTemporaries","ColumnViewReuse","FunctionReusability","ConstMatrixFunction","EmptyMatrixFunction","VectorInMatrixFunction","MixedVectorColumnViewRowVariable",NULL,"FunctionWithConstSignature","VectorCapture","TestFill","TestCopy","ReduceSum","ReduceMax","ReduceLOr","ReduceLAnd","ReduceStrict"};
    for(int i=0;i<17;++i) if(names[i]) {
        printf("matrix:Matrix.%s %d\n",names[i],all_widths(i,0));
        printf("vector_matrix:VectorMatrix.%s %d\n",names[i],all_widths(i,1));
    }
    printf("matrix:Matrix.IntegerVector %d\n",integer_case(1));
    int integers=1; for(size_t width=1;width<=4;++width) integers=integers && integer_case(width);
    printf("vector_matrix:VectorMatrix.IntegerVector %d\n",integers);
    return 0;
}
