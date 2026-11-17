#include "chem.h"
#include "chem_matrix.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define D(m,r,c) (*(double *)chem_matrix_at((m),(r),(c)))
#define CD(m,r,c) (*(const double *)chem_matrix_at_const((m),(r),(c)))
#define I(m,r,c) (*(int *)chem_matrix_at((m),(r),(c)))
static int failures;
static void report(int grouped, const char *name, int ok)
{
    printf("%s:%s.%s %s %s\n",grouped ? "vector_matrix" : "matrix",grouped ? "VectorMatrix" : "Matrix",name,
           ok ? "PASS" : "FAIL",ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static int initd(ChemMatrix *m, size_t rows, size_t cols, size_t group, double value)
{ return chem_matrix_init(m,rows,cols,group,sizeof(value),&value)==CHEM_OK; }
static int initi(ChemMatrix *m, size_t rows, size_t cols, size_t group, int value)
{ return chem_matrix_init(m,rows,cols,group,sizeof(value),&value)==CHEM_OK; }

enum { SUM3, MULTI, SUM3_VECTOR, MATRIX_VECTOR, TEMPS, REUSE_COLS, READ_ONLY, NOOP,
       VECTOR_SUM, TWO_MATRICES_VECTOR, CONST_VECTOR, MUTABLE_VECTOR, VECTOR_TEN,
       VECTOR_COPY, MIXED, INTEGER_VECTOR, CAPTURE, BAD_COLUMN, FILL_VECTOR };
typedef struct { int mode, count; double scale, sum; const int *index, *lookup; } Context;
static int body(ChemMatrix *const *m, void *const *vectors, void *data)
{
    Context *ctx=data; ++ctx->count;
    if (ctx->mode==BAD_COLUMN) return chem_matrix_at(m[0],0,5) ? CHEM_OK : CHEM_INVALID;
    for (size_t r=0;r<m[0]->rows;++r) {
        double *v=vectors ? vectors[0] : NULL;
        switch (ctx->mode) {
        case SUM3: D(m[0],r,2)=ctx->scale*(CD(m[0],r,0)+CD(m[0],r,1)+CD(m[0],r,2)); break;
        case MULTI: D(m[0],r,1)=ctx->scale*(CD(m[0],r,0)+CD(m[1],r,2)); break;
        case SUM3_VECTOR: D(m[0],r,2)=v[r]*(CD(m[0],r,0)+CD(m[0],r,1)+CD(m[0],r,2)); break;
        case MATRIX_VECTOR: D(m[0],r,1)=3*(CD(m[0],r,0)+v[r]); break;
        case TEMPS: {
            double product=CD(m[0],r,0)*CD(m[0],r,1), sum=CD(m[0],r,0)+CD(m[0],r,1);
            D(m[0],r,2)=product+sum; D(m[0],r,3)=product-sum; D(m[0],r,4)=product*sum; break;
        }
        case REUSE_COLS:
            D(m[0],r,1)=2*CD(m[0],r,0);
            D(m[0],r,2)=CD(m[0],r,0)+CD(m[0],r,1);
            D(m[0],r,3)=CD(m[0],r,2)*CD(m[0],r,1); break;
        case READ_ONLY: ctx->sum+=CD(m[0],r,2); break;
        case NOOP: break;
        case VECTOR_SUM: D(m[0],r,0)=v[r]+((double *)vectors[1])[r]; break;
        case TWO_MATRICES_VECTOR: D(m[0],r,1)=CD(m[0],r,0)+v[r]; D(m[1],r,1)=CD(m[1],r,0)+v[r]; break;
        case CONST_VECTOR: D(m[0],r,0)=v[r]*2; break;
        case MUTABLE_VECTOR: v[r]=CD(m[0],r,0)*3; break;
        case VECTOR_TEN: D(m[0],r,0)=v[r]*10; break;
        case VECTOR_COPY: D(m[0],r,0)=v[r]; break;
        case MIXED: D(m[0],r,2)=CD(m[0],r,0)+CD(m[0],r,1)+v[r]; break;
        case INTEGER_VECTOR: I(m[0],r,0)=((int *)vectors[0])[r]*2; break;
        case CAPTURE:
            for (int i=0;i<3;++i) D(m[0],r,0)+=ctx->lookup[ctx->index[i]];
            D(m[0],r,1)+=ctx->lookup[ctx->index[1]]; break;
        case FILL_VECTOR: v[r]=3.2; break;
        default: return CHEM_INVALID;
        }
    }
    return CHEM_OK;
}
static int bind(ChemMatrixFunction *f, ChemMatrix **m, size_t count, size_t vectors, Context *ctx)
{ return chem_matrix_function_init(f,m,count,vectors,body,ctx)==CHEM_OK; }

static int array_function(size_t group)
{
    ChemMatrix m={0},other={0}; ChemMatrix *args[]={&m}; ChemMatrixFunction f={0}; Context ctx={.mode=SUM3,.scale=4};
    int ok=initd(&m,5,3,group,-1) && initd(&other,3,3,group,-1);
    if (!ok) goto done;
    for (size_t r=0;r<5;++r) for (size_t c=0;c<3;++c) D(&m,r,c)=(double)r-2+10*c;
    ok=bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    for (size_t r=0;r<5;++r) {
        double a=(double)r-2,b=(double)r+8,c=(double)r+18;
        ok &= D(&m,r,0)==a && D(&m,r,1)==b && D(&m,r,2)==4*(a+b+c);
    }
    args[0]=&other; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    for (size_t r=0;r<3;++r) ok &= D(&other,r,0)==-1 && D(&other,r,1)==-1 && D(&other,r,2)==-12;
done:
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&m); chem_matrix_destroy(&other); return ok;
}
static int multi_array(size_t group, int new_rows)
{
    ChemMatrix a={0},b={0},a2={0},b2={0}; ChemMatrix *args[]={&a,&b}; ChemMatrixFunction f={0};
    Context ctx={.mode=MULTI,.scale=new_rows?2:1}; size_t n=new_rows?5:3;
    int ok=initd(&a,3,2,group,0) && initd(&b,3,3,group,0);
    if (!ok || !bind(&f,args,2,0,&ctx)) { ok=0; goto done; }
    if (new_rows) {
        ok=initd(&a2,n,2,group,0) && initd(&b2,n,3,group,0); args[0]=&a2; args[1]=&b2;
    }
    if (!ok) goto done;
    for (size_t r=0;r<n;++r) {
        D(args[0],r,0)=new_rows ? r+1 : r;
        D(args[1],r,0)=2*r; D(args[1],r,2)=new_rows ? r*10 : r*4;
    }
    ok=chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    for (size_t r=0;r<n;++r)
        ok &= D(args[0],r,1)==ctx.scale*(D(args[0],r,0)+D(args[1],r,2)) && D(args[1],r,0)==2*r;
    if (new_rows) {
        chem_matrix_destroy(&a2); chem_matrix_destroy(&b2);
        ok &= initd(&a2,2,2,group,0) && initd(&b2,2,3,group,0);
        D(&a2,0,0)=10; D(&a2,1,0)=20; D(&b2,0,2)=5; D(&b2,1,2)=15;
        ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK && D(&a2,0,1)==30 && D(&a2,1,1)==70;
    }
done:
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&a); chem_matrix_destroy(&b);
    chem_matrix_destroy(&a2); chem_matrix_destroy(&b2); return ok;
}
static int matrix_vector(size_t group, int sum_three, int changed_rows)
{
    ChemMatrix initial={0},m={0}; ChemMatrix *args[]={&initial}; ChemMatrixFunction f={0};
    Context ctx={.mode=sum_three ? SUM3_VECTOR : MATRIX_VECTOR};
    size_t n=changed_rows ? 5 : 5, lengths[]={n}; double vec[5]; void *vectors[]={vec};
    int ok=initd(&initial,changed_rows?3:n,3,group,0) && bind(&f,args,1,1,&ctx);
    if (!ok) goto done;
    if (changed_rows) { ok=initd(&m,n,3,group,0); args[0]=&m; }
    if (!ok) goto done;
    for (size_t r=0;r<n;++r) {
        if (sum_three) {
            for (size_t c=0;c<3;++c) D(args[0],r,c)=(double)r-2+10*c;
            vec[r]=100*(r+1);
        } else { D(args[0],r,0)=r+1; vec[r]=10*r; }
    }
    ok=chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
    for (size_t r=0;r<n;++r) {
        if (sum_three) ok &= D(args[0],r,2)==vec[r]*((double)r-2+(double)r+8+(double)r+18) &&
                             D(args[0],r,0)==(double)r-2 && D(args[0],r,1)==(double)r+8;
        else ok &= D(args[0],r,1)==3*((double)r+1+vec[r]);
    }
done:
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&initial); chem_matrix_destroy(&m); return ok;
}
static int simple_operation(size_t group, int mode)
{
    size_t rows=mode==TEMPS ? 4 : 3, cols=mode==TEMPS ? 5 : 4;
    ChemMatrix m={0}; ChemMatrix *args[]={&m}; ChemMatrixFunction f={0}; Context ctx={.mode=mode};
    int ok=initd(&m,rows,cols,group,0);
    if (!ok) return 0;
    for (size_t r=0;r<rows;++r) {
        D(&m,r,0)=r+1;
        if (mode==TEMPS) D(&m,r,1)=(r+1)*10;
        if (mode==READ_ONLY) for (size_t c=0;c<cols;++c) D(&m,r,c)=r*10+c;
    }
    ok=bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    for (size_t r=0;r<rows && ok;++r) {
        double a=r+1;
        if (mode==TEMPS) ok &= D(&m,r,2)==10*a*a+11*a && D(&m,r,3)==10*a*a-11*a && D(&m,r,4)==110*a*a*a;
        if (mode==REUSE_COLS) ok &= D(&m,r,1)==2*a && D(&m,r,2)==3*a && D(&m,r,3)==6*a*a;
    }
    if (mode==READ_ONLY) ok &= D(&m,0,0)==0 && D(&m,1,2)==12 && ctx.sum==36;
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&m); return ok;
}
static int reusable(size_t group)
{
    ChemMatrix a={0},b={0},c={0}; ChemMatrix *args[]={&a}; ChemMatrixFunction f={0};
    Context ctx={.mode=SUM3,.scale=2};
    int ok=initd(&a,2,3,group,1) && initd(&b,2,3,group,5) && initd(&c,2,3,group,0);
    if (!ok) goto done;
    for (size_t r=0;r<2;++r) for (size_t col=0;col<3;++col) D(&a,r,col)=r+col;
    D(&c,1,0)=10;
    ok=bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    ok &= D(&a,0,2)==6 && D(&a,1,2)==12;
    args[0]=&b; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK && D(&b,0,2)==30 && D(&b,1,2)==30;
    args[0]=&c; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK && D(&c,0,2)==0 && D(&c,1,2)==20;
done:
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&a); chem_matrix_destroy(&b); chem_matrix_destroy(&c); return ok;
}
static int empty_function(size_t group)
{
    ChemMatrix a={0},b={0}; ChemMatrix *args[]={&a}; ChemMatrixFunction f={0}; Context ctx={.mode=NOOP};
    int ok=initd(&a,0,3,group,1) && initd(&b,3,0,group,1) && bind(&f,args,1,0,&ctx);
    if (ok) ok=chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    chem_matrix_function_destroy(&f); args[0]=&b;
    if (ok) ok=bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&a); chem_matrix_destroy(&b); return ok;
}

/* Dimension errors are validated at invocation, not at function creation. */
enum { BAD_ROW_CREATE, BAD_COLUMN_ACCESS, BAD_MATRIX_COLUMNS, BAD_VECTOR_ROWS,
       BAD_TWO_MATRIX_ROWS, BAD_COLUMNS_AT_CALL, SMALL_VECTOR, LARGE_VECTOR,
       EMPTY_VECTOR, EMPTY_MATRIX, BOTH_EMPTY, TWO_VECTOR_SIZES, TWO_MATRICES_VECTOR_ROWS,
       VECTOR_MATCHES_ONE, WRONG_VECTOR_AT_CALL };
static int dimension_case(size_t group, int test)
{
    ChemMatrix a={0},b={0},c={0},d={0}; ChemMatrix *creation[]={&a,&b}, *args[]={&a,&b};
    ChemMatrixFunction f={0}; Context ctx={.mode=test==BAD_COLUMN_ACCESS ? BAD_COLUMN : NOOP};
    double v1[10]={0},v2[10]={0}; void *vectors[]={v1,v2}; size_t lengths[]={5,5};
    size_t rows=5, columns=3, nmat=1, nvec=0;
    if (test==BAD_ROW_CREATE) { rows=3; nmat=2; }
    if (test==BAD_COLUMN_ACCESS) { rows=3; columns=4; }
    if (test==BAD_MATRIX_COLUMNS) { rows=3; columns=4; }
    if (test==BAD_VECTOR_ROWS) { rows=3; columns=2; nvec=1; lengths[0]=3; }
    if (test==BAD_TWO_MATRIX_ROWS) { rows=3; columns=2; nmat=2; }
    if (test==BAD_COLUMNS_AT_CALL) { rows=4; columns=3; }
    if (test>=SMALL_VECTOR && test<=BOTH_EMPTY) { nvec=1; if (test==EMPTY_MATRIX || test==BOTH_EMPTY) rows=0; }
    if (test==TWO_VECTOR_SIZES) { nvec=2; lengths[1]=3; }
    if (test==TWO_MATRICES_VECTOR_ROWS || test==VECTOR_MATCHES_ONE) { nmat=2; nvec=1; }
    if (test==WRONG_VECTOR_AT_CALL) { rows=3; columns=2; nvec=1; lengths[0]=3; }
    int ok=initd(&a,rows,columns,group,0);
    if (nmat==2) ok &= initd(&b,test==BAD_ROW_CREATE ? 4 : rows,test==BAD_ROW_CREATE ? 3 : 3,group,0);
    if (!ok || !bind(&f,creation,nmat,nvec,&ctx)) { ok=0; goto done; }
    if (test==BAD_ROW_CREATE) {
        ok=initd(&c,5,3,group,0) && initd(&d,5,3,group,0); args[0]=&c; args[1]=&d;
        ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
        args[0]=&a; args[1]=&b; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
    } else if (test==BAD_COLUMN_ACCESS) ok=chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
    else if (test==BAD_MATRIX_COLUMNS) {
        ok=initd(&c,3,5,group,0) && initd(&d,7,4,group,0);
        ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
        args[0]=&c; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
        args[0]=&d; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    } else if (test==BAD_VECTOR_ROWS) {
        ok=initd(&c,5,2,group,0); args[0]=&c;
        ok &= chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_INVALID;
        args[0]=&a; lengths[0]=5; ok &= chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_INVALID;
    } else if (test==BAD_TWO_MATRIX_ROWS) {
        ok=initd(&c,5,2,group,0); args[0]=&c;
        ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
    } else if (test==BAD_COLUMNS_AT_CALL) {
        ok=initd(&b,4,4,group,0) && initd(&c,4,2,group,0) && initd(&d,7,3,group,0);
        args[0]=&b; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
        args[0]=&c; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
        args[0]=&d; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    } else if (test==WRONG_VECTOR_AT_CALL) {
        ok=chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
        lengths[0]=5; ok &= chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_INVALID;
    } else {
        if (test==SMALL_VECTOR) lengths[0]=3;
        if (test==LARGE_VECTOR) lengths[0]=10;
        if (test==EMPTY_VECTOR || test==BOTH_EMPTY) { lengths[0]=0; vectors[0]=NULL; }
        if (test==TWO_MATRICES_VECTOR_ROWS) {
            chem_matrix_destroy(&b); ok &= initd(&b,4,3,group,0);
        }
        if (test==VECTOR_MATCHES_ONE) lengths[0]=4;
        int expected=test==BOTH_EMPTY ? CHEM_OK : CHEM_INVALID;
        ok &= chem_matrix_function_call(&f,args,vectors,lengths)==expected;
    }
done:
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&a); chem_matrix_destroy(&b);
    chem_matrix_destroy(&c); chem_matrix_destroy(&d); return ok;
}

static int vector_operation(size_t group, int mode)
{
    size_t n=mode==VECTOR_SUM ? 5 : (mode==TWO_MATRICES_VECTOR || mode==MIXED ? 4 : 3);
    size_t cols=mode==MIXED || mode==CAPTURE ? 3 : 2, nmat=mode==TWO_MATRICES_VECTOR ? 2 : 1;
    size_t nvec=mode==CAPTURE ? 0 : (mode==VECTOR_SUM ? 2 : 1);
    ChemMatrix a={0},b={0}; ChemMatrix *args[]={&a,&b}; ChemMatrixFunction f={0}; Context ctx={.mode=mode};
    double v1[5]={0},v2[5]={0}; int integers[]={10,20,30}, index[]={2,0,1},lookup[]={10,20,30};
    void *vectors[]={mode==INTEGER_VECTOR ? (void *)integers : (void *)v1,v2}; size_t lengths[]={n,n};
    int ok=(mode==INTEGER_VECTOR ? initi(&a,n,cols,group,0) : initd(&a,n,cols,group,0));
    if (nmat==2) ok &= initd(&b,n,3,group,0);
    if (!ok) goto done;
    ctx.index=index; ctx.lookup=lookup;
    for (size_t r=0;r<n;++r) {
        if (mode==VECTOR_SUM) { v1[r]=r+1; v2[r]=(r+1)*10; }
        else if (mode==TWO_MATRICES_VECTOR) { D(&a,r,0)=r+1; D(&b,r,0)=10*r; v1[r]=2*r; }
        else if (mode==CONST_VECTOR) v1[r]=10*(r+1);
        else if (mode==MUTABLE_VECTOR) { D(&a,r,0)=r+1; v1[r]=5*(r+1); }
        else if (mode==VECTOR_TEN || mode==VECTOR_COPY) v1[r]=r+1;
        else if (mode==MIXED) { D(&a,r,0)=r+1; D(&a,r,1)=10*(r+1); v1[r]=100*(r+1); }
    }
    ok=bind(&f,args,nmat,nvec,&ctx) && chem_matrix_function_call(&f,args,nvec?vectors:NULL,nvec?lengths:NULL)==CHEM_OK;
    for (size_t r=0;r<n && ok;++r) {
        if (mode==VECTOR_SUM) ok &= D(&a,r,0)==11*(r+1);
        if (mode==TWO_MATRICES_VECTOR) ok &= D(&a,r,1)==r+1+v1[r] && D(&b,r,1)==10*r+v1[r];
        if (mode==CONST_VECTOR) ok &= D(&a,r,0)==20*(r+1) && v1[r]==10*(r+1);
        if (mode==MUTABLE_VECTOR) ok &= v1[r]==3*(r+1) && D(&a,r,0)==r+1;
        if (mode==VECTOR_TEN) ok &= D(&a,r,0)==10*(r+1);
        if (mode==VECTOR_COPY) ok &= D(&a,r,0)==r+1;
        if (mode==MIXED) ok &= D(&a,r,2)==111*(r+1);
        if (mode==INTEGER_VECTOR) ok &= I(&a,r,0)==integers[r]*2;
        if (mode==CAPTURE) ok &= D(&a,r,0)==60 && D(&a,r,1)==10 && D(&a,r,2)==0;
    }
    if (mode==VECTOR_TEN && ok) {
        for (size_t r=0;r<n;++r) v1[r]=r+5;
        ok=chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
        for (size_t r=0;r<n;++r) ok &= D(&a,r,0)==10*(r+5);
    }
done:
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&a); chem_matrix_destroy(&b); return ok;
}
static int fill_case(size_t group)
{
    ChemMatrix m={0}; if (!initd(&m,3,2,group,0)) return 0;
    double x=3.2,vec[3]={0},temporary[3]={9.9,9.9,9.9};
    ChemMatrix *args[]={&m}; void *vectors[]={vec}; size_t lengths[]={3};
    Context ctx={.mode=FILL_VECTOR}; ChemMatrixFunction f={0};
    int ok=chem_matrix_fill_column(&m,1,&x)==CHEM_OK && bind(&f,args,1,1,&ctx) &&
           chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
    for (size_t r=0;r<3;++r) ok &= D(&m,r,0)==0 && D(&m,r,1)==3.2 && vec[r]==3.2;
    ok &= chem_matrix_set_column(&m,0,temporary,3)==CHEM_OK;
    for (size_t r=0;r<3;++r) ok &= D(&m,r,0)==9.9;
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&m); return ok;
}
static int copy_case(size_t group)
{
    ChemMatrix m={0}; if (!initd(&m,3,2,group,0)) return 0;
    double vec[]={3.2,4.2,1.3},out[3]={0},temporary[3]={0},zero=0;
    int ok=chem_matrix_set_column(&m,1,vec,3)==CHEM_OK;
    for (size_t r=0;r<3;++r) ok &= D(&m,r,0)==0 && D(&m,r,1)==vec[r];
    ok &= chem_matrix_get_column(&m,1,out,3)==CHEM_OK;
    for (size_t r=0;r<3;++r) ok &= out[r]==vec[r];
    ok &= chem_matrix_copy_column(&m,0,&m,1)==CHEM_OK;
    for (size_t r=0;r<3;++r) ok &= D(&m,r,0)==vec[r] && D(&m,r,1)==vec[r];
    ok &= chem_matrix_fill_column(&m,0,&zero)==CHEM_OK && chem_matrix_copy_column(&m,0,&m,1)==CHEM_OK;
    ok &= chem_matrix_get_column(&m,1,temporary,3)==CHEM_OK && chem_matrix_fill_column(&m,1,&zero)==CHEM_OK &&
          chem_matrix_set_column(&m,1,temporary,3)==CHEM_OK;
    for (size_t r=0;r<3;++r) ok &= D(&m,r,0)==vec[r] && D(&m,r,1)==vec[r];
    chem_matrix_destroy(&m); return ok;
}

enum { REDUCE_SUM, REDUCE_MAX, REDUCE_OR, REDUCE_AND, REDUCE_STRICT };
typedef struct { int mode; double result; } Reduction;
static void reducer(const void *value, void *accumulator, void *context)
{
    double x=*(const double *)value, *result=accumulator; int mode=*(int *)context;
    if (mode==REDUCE_SUM) *result+=x*x;
    else if (mode==REDUCE_MAX) { if (x>*result) *result=x; }
    else if (mode==REDUCE_OR) *result=(*result!=0 || x>10);
    else if (mode==REDUCE_AND) *result=(*result!=0 && isfinite(x));
    else *result+=1;
}
static int reduce_body(ChemMatrix *const *matrices, void *const *vectors, void *context)
{
    (void)vectors; Reduction *r=context;
    int status=chem_matrix_reduce_column(matrices[0],0,r->mode==REDUCE_STRICT,reducer,&r->result,&r->mode);
    if (status || r->mode==REDUCE_STRICT) return status;
    return chem_matrix_reduce_column(matrices[0],1,0,reducer,&r->result,&r->mode);
}
static int reduction(size_t group, int mode)
{
    ChemMatrix m={0}; ChemMatrix *args[]={&m}; ChemMatrixFunction f={0};
    int ok=initd(&m,3,mode==REDUCE_STRICT ? 1 : 2,group,mode==REDUCE_OR || mode==REDUCE_AND ? 1 : 0);
    if (!ok) return 0;
    if (mode==REDUCE_SUM) for (size_t r=0;r<3;++r) for (size_t c=0;c<2;++c) D(&m,r,c)=2*r+c+1;
    if (mode==REDUCE_MAX) {
        D(&m,0,0)=1; D(&m,0,1)=8; D(&m,1,0)=5; D(&m,1,1)=2; D(&m,2,0)=3; D(&m,2,1)=4;
    }
    Reduction result={mode,mode==REDUCE_MAX ? -DBL_MAX : mode==REDUCE_AND ? 1 : 0};
    ok=chem_matrix_function_init(&f,args,1,0,reduce_body,&result)==CHEM_OK &&
       chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    ok &= result.result==(mode==REDUCE_SUM ? 91 : mode==REDUCE_MAX ? 8 : mode==REDUCE_AND ? 1 : mode==REDUCE_STRICT ? 3 : 0);
    if (mode==REDUCE_OR || mode==REDUCE_AND) {
        D(&m,mode==REDUCE_OR ? 2 : 1,mode==REDUCE_OR ? 1 : 0)=mode==REDUCE_OR ? 42 : NAN;
        result.result=mode==REDUCE_AND ? 1 : 0;
        ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK && result.result==(mode==REDUCE_OR ? 1 : 0);
    }
    chem_matrix_function_destroy(&f); chem_matrix_destroy(&m); return ok;
}

static int apply_case(int grouped, int (*fn)(size_t))
{
    if (!grouped) return fn(0);
    int ok=1; for (size_t g=1;g<=4;++g) ok &= fn(g); return ok;
}
static int apply_mode(int grouped, int (*fn)(size_t,int), int mode)
{
    if (!grouped) return fn(0,mode);
    int ok=1; for (size_t g=1;g<=4;++g) ok &= fn(g,mode); return ok;
}
static int plain_multi(size_t g) { return multi_array(g,0); }
static int resized_multi(size_t g) { return multi_array(g,1); }
static int use_vector(size_t g) { return matrix_vector(g,1,0); }
static int resized_vector(size_t g) { return matrix_vector(g,0,1); }

int main(void)
{
    static const struct { const char *name; int mode; } dims[]={
        {"MismatchedRowDimensions",BAD_ROW_CREATE},{"MismatchedColumnDimensions",BAD_COLUMN_ACCESS},
        {"WrongMatrixDimensions",BAD_MATRIX_COLUMNS},{"MismatchedRowsAtInvocation",BAD_VECTOR_ROWS},
        {"MultipleMatricesMismatchedRowsAtInvocation",BAD_TWO_MATRIX_ROWS},{"WrongColumnCountAtInvocation",BAD_COLUMNS_AT_CALL},
        {"VectorTooSmall",SMALL_VECTOR},{"VectorTooLarge",LARGE_VECTOR},{"EmptyVectorNonEmptyMatrix",EMPTY_VECTOR},
        {"NonEmptyVectorEmptyMatrix",EMPTY_MATRIX},{"EmptyVectorEmptyMatrix",BOTH_EMPTY},
        {"MultipleVectorsDifferentSizes",TWO_VECTOR_SIZES},{"MultipleMatricesDifferentRowsVector",TWO_MATRICES_VECTOR_ROWS},
        {"VectorSizeMatchesOneMatrixOnly",VECTOR_MATCHES_ONE},{"FunctionInvocationWithWrongSizedVector",WRONG_VECTOR_AT_CALL}};
    static const struct { const char *name; int mode; } vectors[]={
        {"MultipleVectorsSameSize",VECTOR_SUM},{"MultipleMatricesOneVector",TWO_MATRICES_VECTOR},
        {"ConstVector",CONST_VECTOR},{"MutableVector",MUTABLE_VECTOR},
        {"FunctionReusabilityWithVectors",VECTOR_TEN},{"MixedVectorColumnViewRowVariable",MIXED},
        {"IntegerVector",INTEGER_VECTOR},{"FunctionWithConstSignature",VECTOR_COPY},{"VectorCapture",CAPTURE}};
    static const struct { const char *name; int mode; } reductions[]={
        {"ReduceSum",REDUCE_SUM},{"ReduceMax",REDUCE_MAX},{"ReduceLOr",REDUCE_OR},
        {"ReduceLAnd",REDUCE_AND},{"ReduceStrict",REDUCE_STRICT}};
    for (int g=0;g<2;++g) {
        report(g,"ArrayFunction",apply_case(g,array_function));
        report(g,"MultiMatrixArrayFunction",apply_case(g,plain_multi));
        report(g,"MultipleTemporaries",apply_mode(g,simple_operation,TEMPS));
        report(g,"ColumnViewReuse",apply_mode(g,simple_operation,REUSE_COLS));
        report(g,"FunctionReusability",apply_case(g,reusable));
        report(g,"ConstMatrixFunction",apply_mode(g,simple_operation,READ_ONLY));
        report(g,"EmptyMatrixFunction",apply_case(g,empty_function));
        report(g,"MultiMatrixDifferentRowsFromCreation",apply_case(g,resized_multi));
        report(g,"MatrixVectorDifferentRowsFromCreation",apply_case(g,resized_vector));
        report(g,"VectorInMatrixFunction",apply_case(g,use_vector));
        for (size_t i=0;i<sizeof dims/sizeof dims[0];++i) report(g,dims[i].name,apply_mode(g,dimension_case,dims[i].mode));
        for (size_t i=0;i<sizeof vectors/sizeof vectors[0];++i) report(g,vectors[i].name,apply_mode(g,vector_operation,vectors[i].mode));
        report(g,"TestFill",apply_case(g,fill_case)); report(g,"TestCopy",apply_case(g,copy_case));
        for (size_t i=0;i<sizeof reductions/sizeof reductions[0];++i)
            report(g,reductions[i].name,apply_mode(g,reduction,reductions[i].mode));
    }
    return failures ? 1 : 0;
}
