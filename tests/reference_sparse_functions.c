#include "chem.h"
#include "chem_sparse.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define D(m,b,e) (*(double *)chem_matrix_at((m),(b),(e)))
#define CD(m,b,e) (*(const double *)chem_matrix_at_const((m),(b),(e)))
#define S(m,b,r,c) (*(double *)chem_sparse_at((m),(b),(r),(c)))
static int failures;
static void report(int grouped, int column, const char *name, int ok)
{
    printf("sparse_matrix_%s_compressed_%s_ordering:Sparse%sCompressed%sMatrix.%s %s %s\n",
           grouped ? "vector" : "standard",column ? "column" : "row",grouped ? "Vector" : "",
           column ? "Column" : "Row",name,ok ? "PASS" : "FAIL",ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static int make(ChemSparseMatrix *m, size_t n, size_t blocks, const ChemSparseCoordinate *coords, size_t count,
                size_t group, int column, double initial)
{
    return chem_sparse_init(m,n,blocks,coords,count,column ? CHEM_SPARSE_COLUMN : CHEM_SPARSE_ROW,
                            group,sizeof initial,&initial)==CHEM_OK;
}
enum { ARRAY, MULTI, TEMPS, REUSE, DIAGONAL_SUM, DIFFERENT, MIX_DENSE, RESIZED_MULTI,
       RESIZED_VECTOR, READ_ONLY, NOOP, VECTOR_DOUBLE, VECTOR_SUM, TWO_MATRIX_VECTOR,
       VECTOR_TRIPLE, VECTOR_FROM_MATRIX, MIXED, INTEGER_VECTOR, GET_SLOTS, VECTOR_FILL };
typedef struct { int mode, calls; size_t e[8]; const size_t *entries; double sum; } Context;
static int body(ChemMatrix *const *m, void *const *vectors, void *context)
{
    Context *p=context; ++p->calls;
    for (size_t b=0;b<m[0]->rows;++b) {
        double *v=vectors ? vectors[0] : NULL;
        double a;
        switch (p->mode) {
        case ARRAY:
            D(m[0],b,p->e[3])=2*(CD(m[0],b,p->e[0])+CD(m[0],b,p->e[1])+CD(m[0],b,p->e[2])+CD(m[0],b,p->e[3])); break;
        case MULTI: D(m[0],b,p->e[1])=CD(m[0],b,p->e[0])+CD(m[1],b,p->e[2]); break;
        case TEMPS: {
            a=CD(m[0],b,p->e[0]); double x=CD(m[0],b,p->e[1]), product=a*x, sum=a+x;
            D(m[0],b,p->e[2])=product+sum; D(m[0],b,p->e[3])=product-sum; D(m[0],b,p->e[4])=product*sum; break;
        }
        case REUSE:
            a=CD(m[0],b,p->e[0]); D(m[0],b,p->e[1])=2*a;
            D(m[0],b,p->e[2])=a+CD(m[0],b,p->e[1]);
            D(m[0],b,p->e[3])=CD(m[0],b,p->e[2])*CD(m[0],b,p->e[1]); break;
        case DIAGONAL_SUM:
            D(m[0],b,p->e[2])=2*(CD(m[0],b,p->e[0])+CD(m[0],b,p->e[1])+CD(m[0],b,p->e[2])); break;
        case DIFFERENT:
            D(m[0],b,p->e[2])=CD(m[0],b,p->e[0])+CD(m[0],b,p->e[1])+CD(m[1],b,p->e[3]); break;
        case MIX_DENSE:
            D(m[0],b,p->e[2])=CD(m[0],b,p->e[0])+CD(m[0],b,p->e[1])+CD(m[1],b,0)+CD(m[1],b,1); break;
        case RESIZED_MULTI:
            D(m[0],b,p->e[2])=CD(m[0],b,p->e[0])+CD(m[0],b,p->e[1])+CD(m[1],b,p->e[2]); break;
        case RESIZED_VECTOR: D(m[0],b,p->e[1])=CD(m[0],b,p->e[0])+v[b]; break;
        case READ_ONLY:
            if (m[0]->columns) p->sum+=CD(m[0],b,p->e[0])+CD(m[0],b,p->e[1])+CD(m[0],b,p->e[2]);
            break;
        case NOOP: break;
        case VECTOR_DOUBLE: D(m[0],b,p->e[0])=2*v[b]; break;
        case VECTOR_SUM: D(m[0],b,p->e[0])=v[b]+((double *)vectors[1])[b]; break;
        case TWO_MATRIX_VECTOR: D(m[1],b,p->e[0])=CD(m[0],b,p->e[0])+v[b]; break;
        case VECTOR_TRIPLE: D(m[0],b,p->e[0])=3*v[b]; break;
        case VECTOR_FROM_MATRIX: v[b]=3*CD(m[0],b,p->e[0]); break;
        case MIXED: D(m[0],b,p->e[2])=2*(CD(m[0],b,p->e[0])+CD(m[0],b,p->e[1])+v[b]); break;
        case INTEGER_VECTOR: D(m[0],b,p->e[0])=1.5*((int *)vectors[0])[b]; break;
        case GET_SLOTS:
            for (size_t i=0;i<6;++i) D(m[0],b,p->entries[i])=((int *)vectors[0])[b]+10*i;
            break;
        case VECTOR_FILL: v[b]=3.2; break;
        default: return CHEM_INVALID;
        }
    }
    return CHEM_OK;
}
static int bind(ChemMatrixFunction *f, ChemMatrix **args, size_t count, size_t vectors, Context *ctx)
{ return chem_matrix_function_init(f,args,count,vectors,body,ctx)==CHEM_OK; }
static size_t entry(const ChemSparseMatrix *m, size_t r, size_t c) { return chem_sparse_entry(m,r,c); }

static int array_case(size_t g, int col)
{
    const ChemSparseCoordinate coords[]={{0,0},{1,1},{2,2},{2,3},{3,3}};
    ChemSparseMatrix a={0},b={0}; ChemMatrix *args[]={&a.storage}; ChemMatrixFunction f={0}; Context ctx={.mode=ARRAY};
    int ok=make(&a,4,3,coords,5,g,col,1) && make(&b,4,3,coords,5,g,col,-1);
    if (!ok) goto done;
    S(&a,1,1,1)=2; S(&a,2,2,2)=3; S(&a,2,2,3)=4; S(&a,2,3,3)=5;
    ctx.e[0]=entry(&a,0,0); ctx.e[1]=entry(&a,1,1); ctx.e[2]=entry(&a,2,2); ctx.e[3]=entry(&a,2,3);
    ok=bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    ok &= S(&a,0,2,3)==8 && S(&a,1,2,3)==10 && S(&a,2,2,3)==18 && S(&a,2,3,3)==5;
    for (size_t i=0;i<3;++i) ok &= S(&a,i,0,0)==1 && S(&a,i,1,1)==(i==1?2:1) && S(&a,i,2,2)==(i==2?3:1);
    args[0]=&b.storage; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    for (size_t i=0;i<3;++i) ok &= S(&b,i,2,3)==-8 && S(&b,i,0,0)==-1 && S(&b,i,1,1)==-1 && S(&b,i,2,2)==-1 && S(&b,i,3,3)==-1;
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); return ok;
}
static int multi_case(size_t g, int col)
{
    const ChemSparseCoordinate ac[]={{0,1},{1,2}},bc[]={{0,0},{1,1},{2,2}};
    ChemSparseMatrix a={0},b={0}; ChemMatrix *args[]={&a.storage,&b.storage}; ChemMatrixFunction f={0}; Context ctx={.mode=MULTI};
    int ok=make(&a,3,3,ac,2,g,col,0) && make(&b,3,3,bc,3,g,col,0);
    if (!ok) goto done;
    for (size_t i=0;i<3;++i) {
        S(&a,i,0,1)=i+10; S(&a,i,1,2)=i*2+20;
        S(&b,i,0,0)=i*4; S(&b,i,1,1)=i*3; S(&b,i,2,2)=i*5;
    }
    ctx.e[0]=entry(&a,0,1); ctx.e[1]=entry(&a,1,2); ctx.e[2]=entry(&b,2,2);
    ok=bind(&f,args,2,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    for (size_t i=0;i<3;++i) ok &= S(&a,i,1,2)==i+10+i*5 && S(&a,i,0,1)==i+10 && S(&b,i,0,0)==i*4;
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); return ok;
}
static int temp_or_reuse(size_t g, int col, int mode)
{
    const ChemSparseCoordinate temps[]={{0,1},{1,2},{2,3},{3,4},{4,4}}, reuse[]={{0,1},{1,2},{2,3},{3,3}};
    size_t n=mode==TEMPS?5:4, blocks=mode==TEMPS?4:3;
    ChemSparseMatrix a={0}; ChemMatrix *args[]={&a.storage}; ChemMatrixFunction f={0}; Context ctx={.mode=mode};
    int ok=make(&a,n,blocks,mode==TEMPS?temps:reuse,mode==TEMPS?5:4,g,col,0);
    if (!ok) goto done;
    ctx.e[0]=entry(&a,0,1); ctx.e[1]=entry(&a,1,2); ctx.e[2]=entry(&a,2,3);
    ctx.e[3]=entry(&a,mode==TEMPS?3:3,mode==TEMPS?4:3);
    if (mode==TEMPS) ctx.e[4]=entry(&a,4,4);
    for (size_t i=0;i<blocks;++i) { S(&a,i,0,1)=i+1; if (mode==TEMPS) S(&a,i,1,2)=10*(i+1); }
    ok=bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    for (size_t i=0;i<blocks;++i) {
        double x=i+1;
        if (mode==TEMPS) ok &= S(&a,i,2,3)==10*x*x+11*x && S(&a,i,3,4)==10*x*x-11*x && S(&a,i,4,4)==110*x*x*x;
        else ok &= S(&a,i,1,2)==2*x && S(&a,i,2,3)==3*x && S(&a,i,3,3)==6*x*x;
    }
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); return ok;
}
static int reusable(size_t g, int col)
{
    const ChemSparseCoordinate coords[]={{0,0},{1,1},{2,2}};
    ChemSparseMatrix a={0},b={0},c={0}; ChemMatrix *args[]={&a.storage}; ChemMatrixFunction f={0}; Context ctx={.mode=DIAGONAL_SUM};
    int ok=make(&a,3,2,coords,3,g,col,0) && make(&b,3,2,coords,3,g,col,5) && make(&c,3,2,coords,3,g,col,0);
    if (!ok) goto done;
    for (size_t i=0;i<2;++i) { S(&a,i,0,0)=i; S(&a,i,1,1)=i+1; S(&a,i,2,2)=i+2; }
    S(&c,1,0,0)=10;
    ctx.e[0]=entry(&a,0,0); ctx.e[1]=entry(&a,1,1); ctx.e[2]=entry(&a,2,2);
    ok=bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    ok &= S(&a,0,2,2)==6 && S(&a,1,2,2)==12;
    args[0]=&b.storage; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK && S(&b,0,2,2)==30 && S(&b,1,2,2)==30;
    args[0]=&c.storage; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK && S(&c,0,2,2)==0 && S(&c,1,2,2)==20;
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); chem_sparse_destroy(&c); return ok;
}
static int different(size_t g, int col)
{
    const ChemSparseCoordinate ac[]={{0,1},{1,1},{2,2}},bc[]={{0,0},{1,2},{2,3},{3,4}};
    ChemSparseMatrix a={0},b={0}; ChemMatrix *args[]={&a.storage,&b.storage}; ChemMatrixFunction f={0}; Context ctx={.mode=DIFFERENT};
    int ok=make(&a,3,4,ac,3,g,col,0) && make(&b,5,4,bc,4,g,col,0);
    if (!ok) goto done;
    for (size_t i=0;i<4;++i) {
        S(&a,i,0,1)=2*i+1; S(&a,i,1,1)=3*i+2; S(&a,i,2,2)=5*i+3;
        S(&b,i,0,0)=i+10; S(&b,i,1,2)=i+20; S(&b,i,2,3)=i+30; S(&b,i,3,4)=i+40;
    }
    ctx.e[0]=entry(&a,0,1); ctx.e[1]=entry(&a,1,1); ctx.e[2]=entry(&a,2,2); ctx.e[3]=entry(&b,0,0);
    ok=bind(&f,args,2,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    for (size_t i=0;i<4;++i) ok &= S(&a,i,2,2)==6*i+13 && S(&a,i,0,1)==2*i+1 && S(&b,i,0,0)==i+10;
    ok &= S(&b,1,1,2)==21 && S(&b,2,2,3)==32 && S(&b,3,3,4)==43;
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); return ok;
}
static int mixed_dense(size_t g, int col, size_t dense_group)
{
    const ChemSparseCoordinate coords[]={{0,1},{1,2},{2,3},{3,3}};
    ChemSparseMatrix a={0}; ChemMatrix dense={0},*args[]={&a.storage,&dense}; ChemMatrixFunction f={0}; Context ctx={.mode=MIX_DENSE};
    double zero=0;
    int ok=make(&a,4,3,coords,4,g,col,0) && chem_matrix_init(&dense,3,4,dense_group,sizeof zero,&zero)==CHEM_OK;
    if (!ok) goto done;
    for (size_t i=0;i<3;++i) {
        S(&a,i,0,1)=i+1; S(&a,i,1,2)=i+10; S(&a,i,2,3)=i+20; S(&a,i,3,3)=i+30;
        for (size_t j=0;j<4;++j) D(&dense,i,j)=100*i+10*j;
    }
    ctx.e[0]=entry(&a,0,1); ctx.e[1]=entry(&a,1,2); ctx.e[2]=entry(&a,3,3);
    ok=bind(&f,args,2,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    ok &= S(&a,0,3,3)==21 && S(&a,1,3,3)==223 && S(&a,2,3,3)==425 &&
          S(&a,0,0,1)==1 && S(&a,1,1,2)==11 && S(&a,2,2,3)==22 &&
          D(&dense,0,0)==0 && D(&dense,1,2)==120 && D(&dense,2,3)==230;
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_matrix_destroy(&dense); return ok;
}
static int read_empty(size_t g, int col, int empty)
{
    const ChemSparseCoordinate coords[]={{0,1},{1,2},{2,3},{3,3}};
    ChemSparseMatrix a={0},b={0}; ChemMatrix *args[]={&a.storage}; ChemMatrixFunction f={0}; Context ctx={.mode=READ_ONLY};
    int ok=make(&a,4,3,empty?NULL:coords,empty?0:4,g,col,0);
    if (!ok) goto done;
    if (!empty) {
        ctx.e[0]=entry(&a,0,1); ctx.e[1]=entry(&a,1,2); ctx.e[2]=entry(&a,2,3);
        for (size_t i=0;i<3;++i) { S(&a,i,0,1)=i+1; S(&a,i,1,2)=i+10; }
    }
    ok=bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    if (!empty) ok &= S(&a,0,0,1)==1 && S(&a,1,1,2)==11 && ctx.sum==39;
    if (empty) {
        chem_matrix_function_destroy(&f); ctx.mode=NOOP; ok &= make(&b,4,0,coords,1,g,col,0); args[0]=&b.storage;
        ok &= bind(&f,args,1,0,&ctx) && chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
    }
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); return ok;
}
static int resized(size_t g, int col, int use_vector)
{
    const ChemSparseCoordinate coords[]={{0,1},{1,2},{2,3}};
    ChemSparseMatrix a={0},b={0},c={0},d={0}; ChemMatrix *args[]={&a.storage,&b.storage}; ChemMatrixFunction f={0};
    Context ctx={.mode=use_vector?RESIZED_VECTOR:RESIZED_MULTI}; size_t n=use_vector?5:4;
    double vector[]={10,20,30,40,50}; void *vectors[]={vector}; size_t lengths[]={5};
    int ok=make(&a,use_vector?3:4,3,coords,use_vector?2:3,g,col,0) &&
           (use_vector || make(&b,4,3,coords,3,g,col,0));
    if (!ok) goto done;
    ctx.e[0]=entry(&a,0,1); ctx.e[1]=entry(&a,1,2); if (!use_vector) ctx.e[2]=entry(&a,2,3);
    ok=bind(&f,args,use_vector?1:2,use_vector?1:0,&ctx) && make(&c,use_vector?3:4,n,coords,use_vector?2:3,g,col,0) &&
       (use_vector || make(&d,4,n,coords,3,g,col,0));
    if (!ok) goto done;
    args[0]=&c.storage; if (!use_vector) args[1]=&d.storage;
    for (size_t i=0;i<n;++i) {
        S(&c,i,0,1)=i+1;
        if (!use_vector) { S(&c,i,1,2)=i+10; S(&d,i,2,3)=i+100; }
    }
    ok=chem_matrix_function_call(&f,args,use_vector?vectors:NULL,use_vector?lengths:NULL)==CHEM_OK;
    for (size_t i=0;i<n;++i) ok &= use_vector ? S(&c,i,1,2)==11*(i+1) : S(&c,i,2,3)==111+3*i;
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b);
    chem_sparse_destroy(&c); chem_sparse_destroy(&d); return ok;
}

enum { BLOCK_DIMENSIONS, ELEMENT_DIMENSIONS, BLOCKS_AT_CALL, THREE_BLOCKS_AT_CALL,
       VECTOR_SMALL, VECTOR_LARGE, VECTOR_EMPTY, BLOCKS_EMPTY, BOTH_EMPTY,
       TWO_VECTORS_BAD, TWO_MATRICES_BAD, VECTOR_MATCHES_ONE, VECTOR_WRONG_CALL };
static int dimensions(size_t g, int col, int mode)
{
    const ChemSparseCoordinate coords[]={{0,1},{1,2},{2,3}};
    ChemSparseMatrix a={0},b={0},c={0}; ChemMatrix *args[]={&a.storage,&b.storage,&c.storage};
    ChemMatrixFunction f={0}; Context ctx={.mode=NOOP}; double x[5]={0},y[5]={0};
    void *vectors[]={x,y}; size_t lengths[]={3,3};
    size_t blocks=mode==BLOCKS_EMPTY || mode==BOTH_EMPTY ? 0 : 3;
    size_t nmat=1,nvec=0;
    if (mode==BLOCK_DIMENSIONS || mode==BLOCKS_AT_CALL || mode==TWO_MATRICES_BAD || mode==VECTOR_MATCHES_ONE) nmat=2;
    if (mode==THREE_BLOCKS_AT_CALL) nmat=3;
    if (mode>=VECTOR_SMALL && mode<=BOTH_EMPTY) nvec=1;
    if (mode==TWO_VECTORS_BAD || mode==VECTOR_MATCHES_ONE) nvec=2;
    if (mode==TWO_MATRICES_BAD || mode==VECTOR_WRONG_CALL) nvec=1;
    int ok=make(&a,4,blocks,coords,3,g,col,0);
    if (nmat>=2) ok &= make(&b,4,mode==TWO_MATRICES_BAD ? 4 : 3,coords,3,g,col,0);
    if (nmat>=3) ok &= make(&c,4,3,coords,3,g,col,0);
    if (!ok || !bind(&f,args,nmat,nvec,&ctx)) { ok=0; goto done; }
    if (mode==ELEMENT_DIMENSIONS) ok=chem_sparse_entry(&a,3,3)==SIZE_MAX && !chem_sparse_at(&a,0,3,3);
    else if (mode==BLOCK_DIMENSIONS || mode==BLOCKS_AT_CALL) {
        args[1]=&a.storage; ok=chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
        chem_sparse_destroy(&b); ok &= make(&b,4,4,coords,3,g,col,0); args[1]=&b.storage;
        ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
        if (mode==BLOCKS_AT_CALL) { args[0]=&b.storage; args[1]=&a.storage; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID; }
    } else if (mode==THREE_BLOCKS_AT_CALL) {
        ok=chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_OK;
        chem_sparse_destroy(&c); ok &= make(&c,4,4,coords,3,g,col,0); args[2]=&c.storage;
        ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
        args[0]=&c.storage; args[2]=&a.storage; ok &= chem_matrix_function_call(&f,args,NULL,NULL)==CHEM_INVALID;
    } else if (mode==VECTOR_WRONG_CALL) {
        lengths[0]=2; ok=chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_INVALID;
        lengths[0]=3; ok &= chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
    } else {
        if (mode==VECTOR_SMALL) lengths[0]=2;
        if (mode==VECTOR_LARGE) lengths[0]=4;
        if (mode==VECTOR_EMPTY || mode==BOTH_EMPTY) { lengths[0]=0; vectors[0]=NULL; }
        if (mode==TWO_VECTORS_BAD || mode==VECTOR_MATCHES_ONE) lengths[1]=2;
        int expected=mode==BOTH_EMPTY ? CHEM_OK : CHEM_INVALID;
        ok=chem_matrix_function_call(&f,args,nvec?vectors:NULL,nvec?lengths:NULL)==expected;
    }
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); chem_sparse_destroy(&c); return ok;
}
static int incompatible(size_t sparse_group, size_t other_group, int col, int sparse_other)
{
    const ChemSparseCoordinate coords[]={{0,1},{1,2}}; double z=0;
    ChemSparseMatrix a={0},b={0}; ChemMatrix dense={0},*args[]={&a.storage,&dense}; ChemMatrixFunction f={0}; Context ctx={.mode=NOOP};
    int ok=make(&a,4,3,coords,2,sparse_group,col,0);
    if (sparse_other) { ok &= make(&b,4,3,coords,2,other_group,col,0); args[1]=&b.storage; }
    else ok &= chem_matrix_init(&dense,3,4,other_group,sizeof z,&z)==CHEM_OK;
    if (ok) ok=chem_matrix_function_init(&f,args,2,0,body,&ctx)==CHEM_INVALID;
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); chem_matrix_destroy(&dense); return ok;
}
static int vector_case(size_t g, int col, int mode)
{
    const ChemSparseCoordinate one[]={{0,1}}, three[]={{0,1},{1,2},{2,3}}, all[]={{0,0},{0,2},{1,0},{1,1},{2,1},{2,2}};
    ChemSparseMatrix a={0},b={0}; ChemMatrix *args[]={&a.storage,&b.storage}; ChemMatrixFunction f={0}; Context ctx={.mode=mode};
    double v1[]={5,10,15},v2[]={4,5,6}; int iv[]={5,10,15},block_values[]={0,1,2}; size_t lengths[]={3,3},indices[6];
    void *vectors[]={mode==INTEGER_VECTOR?(void *)iv:mode==GET_SLOTS?(void *)block_values:(void *)v1,v2};
    size_t count=mode==GET_SLOTS?6:mode==MIXED?3:1;
    const ChemSparseCoordinate *pattern=mode==GET_SLOTS?all:mode==MIXED?three:one;
    int ok=make(&a,mode==MIXED?4:3,3,pattern,count,g,col,0);
    size_t nmat=mode==TWO_MATRIX_VECTOR?2:1, nvec=mode==VECTOR_SUM?2:1;
    if (nmat==2 || mode==VECTOR_DOUBLE) ok &= make(&b,3,3,one,1,g,col,0);
    if (!ok) goto done;
    if (mode==GET_SLOTS) {
        for (size_t i=0;i<6;++i) indices[i]=entry(&a,all[i].row,all[i].column);
        ctx.entries=indices;
    } else ctx.e[0]=entry(&a,0,1);
    if (mode==VECTOR_SUM) { v1[0]=1; v1[1]=2; v1[2]=3; }
    if (mode==TWO_MATRIX_VECTOR) for (size_t i=0;i<3;++i) { S(&a,i,0,1)=i+1; v1[i]=10*(i+1); }
    if (mode==VECTOR_FROM_MATRIX) for (size_t i=0;i<3;++i) S(&a,i,0,1)=i+1;
    if (mode==MIXED) {
        ctx.e[1]=entry(&a,1,2); ctx.e[2]=entry(&a,2,3);
        for (size_t i=0;i<3;++i) { S(&a,i,0,1)=i+1; S(&a,i,1,2)=i+2; S(&a,i,2,3)=i+3; v1[i]=100*(i+1); }
    }
    ok=bind(&f,args,nmat,nvec,&ctx) && chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
    for (size_t i=0;i<3 && ok;++i) {
        if (mode==VECTOR_DOUBLE) ok &= S(&a,i,0,1)==2*v1[i];
        if (mode==VECTOR_SUM) ok &= S(&a,i,0,1)==v1[i]+v2[i];
        if (mode==TWO_MATRIX_VECTOR) ok &= S(&b,i,0,1)==11*(i+1);
        if (mode==VECTOR_TRIPLE) ok &= S(&a,i,0,1)==15*(i+1);
        if (mode==VECTOR_FROM_MATRIX) ok &= v1[i]==3*(i+1);
        if (mode==MIXED) ok &= S(&a,i,2,3)==206+204*i;
        if (mode==INTEGER_VECTOR) ok &= S(&a,i,0,1)==1.5*iv[i];
        if (mode==GET_SLOTS) for (size_t j=0;j<6;++j) {
            size_t slot=chem_sparse_slot(&a,0,all[j].row,all[j].column);
            double *by_slot=chem_sparse_block_slot(&a,i,slot);
            ok &= by_slot && *by_slot==i+10*j && S(&a,i,all[j].row,all[j].column)==i+10*j;
        }
    }
    if (mode==VECTOR_DOUBLE && ok) {
        v1[0]=1; v1[1]=2; v1[2]=3;
        ok=chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
        for (size_t i=0;i<3;++i) ok &= S(&a,i,0,1)==2*(i+1);
        args[0]=&b.storage; v1[0]=10; v1[1]=20; v1[2]=30;
        ok &= chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
        for (size_t i=0;i<3;++i) ok &= S(&b,i,0,1)==20*(i+1);
    }
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); return ok;
}
static int fill_or_copy(size_t g, int col, int copy)
{
    const ChemSparseCoordinate coords[]={{0,1},{2,0}},other[]={{1,2}};
    ChemSparseMatrix a={0},b={0}; ChemMatrix dense={0},*args[]={&a.storage}; ChemMatrixFunction f={0};
    Context ctx={.mode=VECTOR_FILL}; double z=0,vec[]={3.2,4.2,1.3},out[3]={0},temp[3]={0};
    size_t lengths[]={3}; void *vectors[]={out};
    int ok=make(&a,3,3,coords,2,g,col,0); if (!ok) goto done;
    size_t e01=entry(&a,0,1),e20=entry(&a,2,0);
    if (!copy) {
        double val=3.2; ok=chem_matrix_fill_column(&a.storage,e01,&val)==CHEM_OK && bind(&f,args,1,1,&ctx) &&
               chem_matrix_function_call(&f,args,vectors,lengths)==CHEM_OK;
        for (size_t i=0;i<3;++i) ok &= S(&a,i,0,1)==3.2 && S(&a,i,2,0)==0 && out[i]==3.2;
        for (size_t i=0;i<3;++i) temp[i]=9.9;
        ok &= chem_matrix_set_column(&a.storage,e20,temp,3)==CHEM_OK;
        for (size_t i=0;i<3;++i) ok &= S(&a,i,2,0)==9.9;
    } else {
        ok=chem_matrix_set_column(&a.storage,e01,vec,3)==CHEM_OK && chem_matrix_get_column(&a.storage,e01,out,3)==CHEM_OK &&
           chem_matrix_copy_column(&a.storage,e20,&a.storage,e01)==CHEM_OK;
        for (size_t i=0;i<3;++i) ok &= out[i]==vec[i] && S(&a,i,0,1)==vec[i] && S(&a,i,2,0)==vec[i];
        ok &= chem_matrix_get_column(&a.storage,e01,temp,3)==CHEM_OK && chem_matrix_fill_column(&a.storage,e01,&z)==CHEM_OK &&
              chem_matrix_set_column(&a.storage,e01,temp,3)==CHEM_OK;
        ok &= chem_matrix_init(&dense,3,2,g,sizeof z,&z)==CHEM_OK;
        for (size_t i=0;i<3;++i) D(&dense,i,1)=100*(i+1);
        ok &= chem_matrix_copy_column(&a.storage,e20,&dense,1)==CHEM_OK && chem_matrix_copy_column(&dense,0,&a.storage,e01)==CHEM_OK;
        ok &= make(&b,3,3,other,1,g,col,0);
        ok &= chem_matrix_copy_column(&b.storage,entry(&b,1,2),&a.storage,e01)==CHEM_OK;
        for (size_t i=0;i<3;++i) ok &= S(&a,i,2,0)==100*(i+1) && D(&dense,i,0)==vec[i] && S(&b,i,1,2)==vec[i];
    }
done:
    chem_matrix_function_destroy(&f); chem_sparse_destroy(&a); chem_sparse_destroy(&b); chem_matrix_destroy(&dense); return ok;
}

static int all2(int grouped, int col, int (*fn)(size_t,int))
{
    if (!grouped) return fn(0,col);
    int ok=1; for (size_t g=1;g<=4;++g) ok &= fn(g,col); return ok;
}
static int all3(int grouped, int col, int (*fn)(size_t,int,int), int mode)
{
    if (!grouped) return fn(0,col,mode);
    int ok=1; for (size_t g=1;g<=4;++g) ok &= fn(g,col,mode); return ok;
}

int main(void)
{
    static const struct { const char *name; int mode; } dims[]={
        {"MismatchedBlockDimensions",BLOCK_DIMENSIONS},{"MismatchedElementDimensions",ELEMENT_DIMENSIONS},
        {"MismatchedBlocksAtInvocation",BLOCKS_AT_CALL},{"MultipleSparseMatricesMismatchedBlocksAtInvocation",THREE_BLOCKS_AT_CALL},
        {"VectorTooSmall",VECTOR_SMALL},{"VectorTooLarge",VECTOR_LARGE},
        {"EmptyVectorNonEmptySparseMatrix",VECTOR_EMPTY},{"NonEmptyVectorEmptySparseMatrix",BLOCKS_EMPTY},
        {"EmptyVectorEmptySparseMatrix",BOTH_EMPTY},{"MultipleVectorsDifferentSizes",TWO_VECTORS_BAD},
        {"MultipleSparseMatricesDifferentBlocksVector",TWO_MATRICES_BAD},
        {"VectorSizeMatchesOneSparseMatrixOnly",VECTOR_MATCHES_ONE},
        {"FunctionInvocationWithWrongSizedVectorSparse",VECTOR_WRONG_CALL}};
    static const struct { const char *name; int mode; } vec[]={
        {"VectorInSparseMatrixFunction",VECTOR_DOUBLE},{"MultipleVectorsSameSize",VECTOR_SUM},
        {"MultipleSparseMatricesOneVector",TWO_MATRIX_VECTOR},{"ConstVectorSparse",VECTOR_TRIPLE},
        {"MutableVectorSparse",VECTOR_FROM_MATRIX},{"FunctionReusabilityWithVectorsSparse",VECTOR_DOUBLE},
        {"MixedVectorBlockViewBlockVariable",MIXED},{"IntegerVectorSparse",INTEGER_VECTOR},
        {"FunctionWithConstSignatureSparse",VECTOR_DOUBLE},{"GetBlockViewByVectorIndex",GET_SLOTS}};
    for (int grouped=0;grouped<2;++grouped) for (int col=0;col<2;++col) {
        report(grouped,col,"ArrayFunction",all2(grouped,col,array_case));
        report(grouped,col,"MultiMatrixArrayFunction",all2(grouped,col,multi_case));
        report(grouped,col,"MultipleTemporaries",all3(grouped,col,temp_or_reuse,TEMPS));
        report(grouped,col,"BlockViewReuse",all3(grouped,col,temp_or_reuse,REUSE));
        report(grouped,col,"FunctionReusability",all2(grouped,col,reusable));
        report(grouped,col,"TwoSparseMatricesDifferentStructure",all2(grouped,col,different));
        report(grouped,col,"SparseAndDenseMatrixFunction",mixed_dense(grouped?1:0,col,0));
        report(grouped,col,"ConstSparseMatrixFunction",all3(grouped,col,read_empty,0));
        report(grouped,col,"EmptySparseMatrixFunction",all3(grouped,col,read_empty,1));
        report(grouped,col,"MultipleSparseMatricesDifferentBlocksFromCreation",all3(grouped,col,resized,0));
        report(grouped,col,"SparseMatrixVectorDifferentBlocksFromCreation",all3(grouped,col,resized,1));
        for (size_t i=0;i<sizeof dims/sizeof dims[0];++i)
            report(grouped,col,dims[i].name,all3(grouped,col,dimensions,dims[i].mode));
        for (size_t i=0;i<sizeof vec/sizeof vec[0];++i)
            report(grouped,col,vec[i].name,all3(grouped,col,vector_case,vec[i].mode));
        report(grouped,col,"TestFill",all3(grouped,col,fill_or_copy,0));
        report(grouped,col,"TestCopy",all3(grouped,col,fill_or_copy,1));
        if (grouped) {
            int ok=1; for (size_t g=1;g<=4;++g) ok &= mixed_dense(g,col,g);
            ok &= mixed_dense(10,col,10); report(grouped,col,"SparseAndVectorMatrixFunction",ok);
            ok=incompatible(2,0,col,0) && incompatible(3,0,col,0) && incompatible(4,0,col,0);
            report(grouped,col,"IncompatibleOrdering",ok);
            ok=incompatible(2,1,col,0) && incompatible(2,3,col,0) && incompatible(4,2,col,0) && incompatible(4,10,col,0);
            report(grouped,col,"IncompatibleVectorOrdering",ok);
            ok=incompatible(1,2,col,1) && incompatible(2,3,col,1) && incompatible(3,4,col,1) && incompatible(4,10,col,1);
            report(grouped,col,"IncompatibleSparseOrdering",ok);
        }
    }
    return failures ? 1 : 0;
}
