/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * C ports of the function/view tests in test_matrix_policy.hpp. */
#include "test.h"
#include "chem_matrix.h"
#include <float.h>

static ChemMatrix mat(size_t r,size_t c,size_t g,double v)
{
    ChemMatrix m={0};CHECK(chem_matrix_create(r,c,g,CHEM_MATRIX_REAL,v,&m)==CHEM_OK);return m;
}
static void set(ChemMatrix *m,size_t r,size_t c,double v)
{ CHECK(chem_matrix_set(m,r,c,v)==CHEM_OK); }
static ChemView column(ChemArray a,size_t c)
{ ChemView v={0};CHECK(chem_array_column(a,c,&v)==CHEM_OK);return v; }
static ChemView read_column(ChemArray a,size_t c)
{ return chem_view_const(column(a,c)); }
static ChemArrayFunction *function(size_t n,const ChemArray *a,ChemArrayCallback cb,void *ctx)
{
    ChemArrayFunction *f=NULL;CHECK(chem_array_function_create(n,a,cb,ctx,&f)==CHEM_OK);return f;
}
static void call(ChemArrayFunction *f,size_t n,const ChemArray *a)
{ CHECK(chem_array_function_call(f,n,a)==CHEM_OK); }
static ChemView temporary(ChemArray a,ChemMatrix *owner)
{ *owner=mat(a.rows,1,a.group,0);return column(chem_array_matrix(owner),0); }
/* Each elementary operation writes the last argument, with read-only inputs. */
typedef struct { int operation;double scale; } Arithmetic;
static void arithmetic(double *v,size_t n,void *context)
{
    Arithmetic *p=context;double value=v[0];
    if(p->operation==0)for(size_t i=1;i+1<n;++i)value+=v[i];
    else if(p->operation==1)for(size_t i=1;i+1<n;++i)value*=v[i];
    else if(p->operation==2)for(size_t i=1;i+1<n;++i)value-=v[i];
    v[n-1]=value*p->scale;
}
static void apply(ChemView dst,size_t n,const ChemView *src,int operation,double scale)
{
    ChemView views[5];CHECK(n<=4);if(n>4)return;
    for(size_t i=0;i<n;++i) views[i]=chem_view_const(src[i]);
    views[n]=dst;
    Arithmetic a={operation,scale};CHECK(chem_views_for_each(n+1,views,arithmetic,&a)==CHEM_OK);
}
static void unary(ChemView dst,ChemView src,double scale)
{ apply(dst,1,&src,0,scale); }
static void binary(ChemView dst,ChemView a,ChemView b,int op,double scale)
{ ChemView inputs[]={a,b};apply(dst,2,inputs,op,scale); }
static ChemStatus sum_three(const ChemArray *a,size_t n,void *ctx)
{
    ChemMatrix owner={0};ChemView t=temporary(a[0],&owner);
    ChemView inputs[]={read_column(a[0],0),read_column(a[0],1),read_column(a[0],2)};
    apply(t,3,inputs,0,1);
    if(n==2)binary(column(a[0],2),t,read_column(a[1],0),1,1);
    else unary(column(a[0],2),t,*(double *)ctx);
    chem_matrix_free(&owner);return CHEM_OK;
}
static ChemStatus multi_matrix(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;ChemMatrix owner={0};ChemView t=temporary(a[0],&owner);
    binary(t,read_column(a[0],0),read_column(a[1],2),0,1);
    unary(column(a[0],1),t,*(double *)ctx);chem_matrix_free(&owner);return CHEM_OK;
}
static ChemStatus matrix_vector(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;ChemMatrix owner={0};ChemView t=temporary(a[0],&owner);
    binary(t,read_column(a[0],0),read_column(a[1],0),0,1);
    unary(column(a[0],1),t,*(double *)ctx);chem_matrix_free(&owner);return CHEM_OK;
}
static ChemStatus add_matrices(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;binary(column(a[0],1),read_column(a[0],0),read_column(a[1],0),0,1);return CHEM_OK;
}
static ChemStatus double_to_last(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;unary(column(a[0],a[0].columns-1),read_column(a[0],0),2);return CHEM_OK;
}
static ChemStatus invalid_column(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;ChemView v;return chem_array_column(a[0],5,&v);
}
static ChemStatus temporaries(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;ChemMatrix p={0},s={0};ChemView prod=temporary(a[0],&p),sum=temporary(a[0],&s);
    binary(prod,read_column(a[0],0),read_column(a[0],1),1,1);
    binary(sum,read_column(a[0],0),read_column(a[0],1),0,1);
    binary(column(a[0],2),prod,sum,0,1);binary(column(a[0],3),prod,sum,2,1);
    binary(column(a[0],4),prod,sum,1,1);
    chem_matrix_free(&p);chem_matrix_free(&s);return CHEM_OK;
}
static ChemStatus reused_columns(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;ChemView c0=read_column(a[0],0),c1=column(a[0],1),c2=column(a[0],2),c3=column(a[0],3);
    unary(c1,c0,2);binary(c2,c0,c1,0,1);binary(c3,c2,c1,1,1);return CHEM_OK;
}
static void accumulate(double *v,size_t n,void *ctx)
{ (void)n;*(double *)ctx+=v[0]; }
static ChemStatus read_matrix(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;ChemMatrix owner={0};ChemView t=temporary(a[0],&owner);
    binary(t,read_column(a[0],0),read_column(a[0],1),0,1);
    double sum=0;ChemView c=read_column(a[0],2);
    CHECK(chem_views_for_each(1,&c,accumulate,&sum)==CHEM_OK);
    chem_matrix_free(&owner);return CHEM_OK;
}
static ChemStatus empty_function(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;if(a[0].columns)CHECK(chem_view_fill(column(a[0],0),99)==CHEM_OK);return CHEM_OK;
}
static ChemStatus vector_to_matrix(const ChemArray *a,size_t n,void *ctx)
{ (void)n;unary(column(a[0],0),read_column(a[1],0),*(double *)ctx);return CHEM_OK; }
static ChemStatus matrix_to_vector(const ChemArray *a,size_t n,void *ctx)
{ (void)n;unary(column(a[1],0),read_column(a[0],0),*(double *)ctx);return CHEM_OK; }
static ChemStatus vectors_sum(const ChemArray *a,size_t n,void *ctx)
{ (void)n;(void)ctx;binary(column(a[0],0),read_column(a[1],0),read_column(a[2],0),0,1);return CHEM_OK; }
static ChemStatus matrices_vector(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;
    for(size_t i=0;i<2;++i)binary(column(a[i],1),read_column(a[i],0),read_column(a[2],0),0,1);
    return CHEM_OK;
}
static ChemStatus copy_vector_first(const ChemArray *a,size_t n,void *ctx)
{ (void)ctx;unary(column(a[0],0),read_column(a[n-1],0),1);return CHEM_OK; }
static ChemStatus mixed(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;(void)ctx;ChemMatrix owner={0};ChemView t=temporary(a[0],&owner);
    binary(t,read_column(a[0],0),read_column(a[0],1),0,1);
    binary(column(a[0],2),t,read_column(a[1],0),0,1);chem_matrix_free(&owner);return CHEM_OK;
}
static void array_test(size_t g,int vector)
{
    ChemMatrix m=mat(5,3,g,0);double v[]={100,200,300,400,500},scale=4;
    for(size_t r=0;r<5;++r)for(size_t c=0;c<3;++c)set(&m,r,c,(double)r-2+10*c);
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_vector(v,5,CHEM_MATRIX_REAL)};
    size_t n=vector?2:1;ChemArrayFunction *f=function(n,a,sum_three,&scale);call(f,n,a);
    for(size_t r=0;r<5;++r){
        CHECK(chem_matrix_get(&m,r,2)==(vector?v[r]:4)*(24+3*r));
        CHECK(chem_matrix_get(&m,r,0)==(double)r-2);CHECK(chem_matrix_get(&m,r,1)==r+8);
    }
    if(!vector){
        ChemMatrix other=mat(3,3,g,-1);a[0]=chem_array_matrix(&other);call(f,1,a);
        for(size_t r=0;r<3;++r){CHECK(chem_matrix_get(&other,r,2)==-12);CHECK(chem_matrix_get(&other,r,0)==-1);CHECK(chem_matrix_get(&other,r,1)==-1);}
        chem_matrix_free(&other);
    }
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void multi_test(size_t g)
{
    ChemMatrix m=mat(3,2,g,0),b=mat(3,3,g,0);double scale=1;
    for(size_t r=0;r<3;++r){
        for(size_t c=0;c<2;++c)set(&m,r,c,r+10*c);
        for(size_t c=0;c<3;++c)set(&b,r,c,2*r+20*c);
        set(&b,r,2,4*r);
    }
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_matrix(&b)};
    ChemArrayFunction *f=function(2,a,multi_matrix,&scale);call(f,2,a);
    for(size_t r=0;r<3;++r){CHECK(chem_matrix_get(&m,r,1)==5*r);CHECK(chem_matrix_get(&m,r,0)==r);CHECK(chem_matrix_get(&b,r,0)==2*r);}
    chem_array_function_free(f);chem_matrix_free(&m);chem_matrix_free(&b);
}
static void flexible_multi(size_t g)
{
    ChemMatrix m=mat(3,2,g,0),b=mat(3,3,g,0);double scale=2;
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_matrix(&b)};
    ChemArrayFunction *f=function(2,a,multi_matrix,&scale);
    chem_matrix_free(&m);chem_matrix_free(&b);m=mat(5,2,g,0);b=mat(5,3,g,0);
    for(size_t r=0;r<5;++r){set(&m,r,0,r+1);set(&b,r,2,10*r);}
    a[0]=chem_array_matrix(&m);a[1]=chem_array_const_matrix(&b);call(f,2,a);
    for(size_t r=0;r<5;++r)CHECK(chem_matrix_get(&m,r,1)==2+22*r);
    chem_matrix_free(&m);chem_matrix_free(&b);m=mat(2,2,g,0);b=mat(2,3,g,0);
    set(&m,0,0,10);set(&m,1,0,20);set(&b,0,2,5);set(&b,1,2,15);
    a[0]=chem_array_matrix(&m);a[1]=chem_array_const_matrix(&b);call(f,2,a);
    CHECK(chem_matrix_get(&m,0,1)==30);CHECK(chem_matrix_get(&m,1,1)==70);
    chem_array_function_free(f);chem_matrix_free(&m);chem_matrix_free(&b);
}
static void flexible_vector(size_t g)
{
    ChemMatrix m=mat(3,3,g,0);double v[5]={0},scale=3;
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(2,a,matrix_vector,&scale);
    chem_matrix_free(&m);m=mat(5,3,g,0);
    for(size_t r=0;r<5;++r){set(&m,r,0,r+1);v[r]=r*10;}
    a[0]=chem_array_matrix(&m);a[1]=chem_array_const_vector(v,5,CHEM_MATRIX_REAL);call(f,2,a);
    for(size_t r=0;r<5;++r)CHECK(chem_matrix_get(&m,r,1)==3+33*r);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void mismatch_rows(size_t g)
{
    ChemMatrix m=mat(3,3,g,1),b=mat(4,3,g,2);
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_matrix(&b)};
    ChemArrayFunction *f=function(2,a,add_matrices,NULL);
    chem_matrix_free(&m);chem_matrix_free(&b);m=mat(5,3,g,0);b=mat(5,3,g,1);
    a[0]=chem_array_matrix(&m);a[1]=chem_array_const_matrix(&b);call(f,2,a);
    chem_matrix_free(&m);chem_matrix_free(&b);m=mat(3,3,g,0);b=mat(4,3,g,1);
    a[0]=chem_array_matrix(&m);a[1]=chem_array_const_matrix(&b);
    CHECK(chem_array_function_call(f,2,a)==CHEM_DIMENSION);
    chem_array_function_free(f);chem_matrix_free(&m);chem_matrix_free(&b);
}
static void mismatch_columns(size_t g)
{
    ChemMatrix m=mat(3,4,g,1);ChemArray a=chem_array_matrix(&m);
    ChemArrayFunction *f=function(1,&a,invalid_column,NULL);
    /* Original only asserts this in debug; the C error check also runs in Release. */
    CHECK(chem_array_function_call(f,1,&a)==CHEM_DIMENSION);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void wrong_dimensions(size_t g,int invocation)
{
    size_t cols=invocation?3:4;double scale=1;
    ChemMatrix m=mat(invocation?4:3,cols,g,1);ChemArray a=chem_array_matrix(&m);
    ChemArrayFunction *f=function(1,&a,invocation?sum_three:double_to_last,&scale);
    if(!invocation){call(f,1,&a);call(f,1,&a);}
    chem_matrix_free(&m);m=mat(invocation?4:3,cols+1,g,2);a=chem_array_matrix(&m);
    CHECK(chem_array_function_call(f,1,&a)==CHEM_DIMENSION);
    if(invocation){
        chem_matrix_free(&m);m=mat(4,2,g,0);a=chem_array_matrix(&m);
        CHECK(chem_array_function_call(f,1,&a)==CHEM_DIMENSION);
    }
    chem_matrix_free(&m);m=mat(7,cols,g,1);a=chem_array_matrix(&m);call(f,1,&a);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void row_invocation(size_t g,int matrices)
{
    ChemMatrix m=mat(3,2,g,0),b=mat(3,3,g,0);double v[5]={0},scale=2;
    ChemArray a[]={chem_array_matrix(&m),matrices?chem_array_const_matrix(&b):chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};
    /* Both kernels read argument 1 and write the first matrix; validation precedes either. */
    ChemArrayFunction *f=function(2,a,vector_to_matrix,&scale);
    chem_matrix_free(&m);m=mat(5,2,g,0);a[0]=chem_array_matrix(&m);
    CHECK(chem_array_function_call(f,2,a)==CHEM_DIMENSION);
    if(!matrices){
        chem_matrix_free(&m);m=mat(3,2,g,0);a[0]=chem_array_matrix(&m);a[1]=chem_array_const_vector(v,5,CHEM_MATRIX_REAL);
        CHECK(chem_array_function_call(f,2,a)==CHEM_DIMENSION);
    }
    chem_array_function_free(f);chem_matrix_free(&m);chem_matrix_free(&b);
}
static void temporary_test(size_t g)
{
    ChemMatrix m=mat(4,5,g,0);
    for(size_t r=0;r<4;++r){set(&m,r,0,r+1);set(&m,r,1,10*(r+1));}
    ChemArray a=chem_array_matrix(&m);ChemArrayFunction *f=function(1,&a,temporaries,NULL);call(f,1,&a);
    for(size_t r=0;r<4;++r){double p=10*(r+1)*(r+1),s=11*(r+1);CHECK(chem_matrix_get(&m,r,2)==p+s);CHECK(chem_matrix_get(&m,r,3)==p-s);CHECK(chem_matrix_get(&m,r,4)==p*s);}
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void reuse_columns_test(size_t g)
{
    ChemMatrix m=mat(3,4,g,0);for(size_t r=0;r<3;++r)set(&m,r,0,r+1);
    ChemArray a=chem_array_matrix(&m);ChemArrayFunction *f=function(1,&a,reused_columns,NULL);call(f,1,&a);
    for(size_t r=0;r<3;++r){CHECK(chem_matrix_get(&m,r,1)==2*(r+1));CHECK(chem_matrix_get(&m,r,2)==3*(r+1));CHECK(chem_matrix_get(&m,r,3)==6*(r+1)*(r+1));}
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void reuse_test(size_t g)
{
    ChemMatrix m=mat(2,3,g,1);double scale=2;ChemArray a=chem_array_matrix(&m);
    ChemArrayFunction *f=function(1,&a,sum_three,&scale);
    for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)set(&m,r,c,r+c);
    call(f,1,&a);
    CHECK(chem_matrix_get(&m,0,2)==6);CHECK(chem_matrix_get(&m,1,2)==12);
    chem_matrix_free(&m);m=mat(2,3,g,5);a=chem_array_matrix(&m);call(f,1,&a);
    CHECK(chem_matrix_get(&m,0,2)==30);CHECK(chem_matrix_get(&m,1,2)==30);
    chem_matrix_free(&m);m=mat(2,3,g,0);for(size_t r=0;r<2;++r)set(&m,r,0,10*r);
    a=chem_array_matrix(&m);call(f,1,&a);CHECK(chem_matrix_get(&m,0,2)==0);CHECK(chem_matrix_get(&m,1,2)==20);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void const_test(size_t g)
{
    ChemMatrix m=mat(3,4,g,0);for(size_t r=0;r<3;++r)for(size_t c=0;c<4;++c)set(&m,r,c,r*10+c);
    ChemArray a=chem_array_const_matrix(&m);ChemArrayFunction *f=function(1,&a,read_matrix,NULL);call(f,1,&a);
    CHECK(chem_matrix_get(&m,0,0)==0);CHECK(chem_matrix_get(&m,1,2)==12);
    CHECK(chem_view_set(column(a,0),0,1)==CHEM_INVALID);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void empty_test(size_t g)
{
    ChemMatrix m=mat(0,3,g,1);ChemArray a=chem_array_matrix(&m);
    ChemArrayFunction *f=function(1,&a,empty_function,NULL);call(f,1,&a);chem_array_function_free(f);chem_matrix_free(&m);
    m=mat(3,0,g,1);a=chem_array_matrix(&m);f=function(1,&a,empty_function,NULL);call(f,1,&a);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void vector_shape_test(size_t g,int which)
{
    /* Preserve each original creation/invocation shape, including the two
     * tests whose descriptive names differ from their actual fixture. */
    size_t rows=5,vrows=which==0?3:which==1?10:which==3?5:0;
    if(which==3||which==4)rows=0;
    if(which>=5)vrows=which==5?5:4;
    ChemMatrix m=mat(rows,which>=6?2:3,g,which<5?1:0);
    ChemMatrix b=mat(which==6?5:rows,3,g,0);
    if(which==6){chem_matrix_free(&m);m=mat(4,2,g,0);}
    double vec[10]={0},vec2[3]={0},scale=1;
    ChemArray a[3]={chem_array_matrix(&m),chem_array_const_vector(vec,vrows,CHEM_MATRIX_REAL),{0}};
    size_t n=2;ChemArrayCallback cb=vector_to_matrix;
    if(which==5){a[2]=chem_array_const_vector(vec2,3,CHEM_MATRIX_REAL);n=3;cb=vectors_sum;}
    if(which>=6){a[1]=chem_array_const_matrix(&b);a[2]=chem_array_const_vector(vec,4,CHEM_MATRIX_REAL);n=3;cb=copy_vector_first;}
    ChemArrayFunction *f=function(n,a,cb,&scale);
    CHECK(chem_array_function_call(f,n,a)==(which==4?CHEM_OK:CHEM_DIMENSION));
    chem_array_function_free(f);chem_matrix_free(&m);chem_matrix_free(&b);
}
static void same_vectors(size_t g)
{
    ChemMatrix m=mat(5,3,g,0);double v[5],w[5];
    for(size_t r=0;r<5;++r){v[r]=r+1;w[r]=10*(r+1);}
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_vector(v,5,CHEM_MATRIX_REAL),chem_array_const_vector(w,5,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(3,a,vectors_sum,NULL);call(f,3,a);
    for(size_t r=0;r<5;++r)CHECK(chem_matrix_get(&m,r,0)==11*(r+1));
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void matrices_one_vector(size_t g)
{
    ChemMatrix m=mat(4,2,g,0),b=mat(4,3,g,0);double v[4];
    for(size_t r=0;r<4;++r){v[r]=r*2;set(&m,r,0,r+1);set(&b,r,0,r*10);}
    ChemArray a[]={chem_array_matrix(&m),chem_array_matrix(&b),chem_array_const_vector(v,4,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(3,a,matrices_vector,NULL);call(f,3,a);
    for(size_t r=0;r<4;++r){CHECK(chem_matrix_get(&m,r,1)==3*r+1);CHECK(chem_matrix_get(&b,r,1)==12*r);}
    chem_array_function_free(f);chem_matrix_free(&m);chem_matrix_free(&b);
}
static void const_vector(size_t g,int integer)
{
    ChemMatrix m={0};CHECK(chem_matrix_create(3,2,g,integer?CHEM_MATRIX_INT:CHEM_MATRIX_REAL,0,&m)==CHEM_OK);
    const double v[]={10,20,30};const int iv[]={10,20,30};double scale=2;
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_vector(integer?(const void *)iv:(const void *)v,3,m.type)};
    ChemArrayFunction *f=function(2,a,vector_to_matrix,&scale);call(f,2,a);
    for(size_t r=0;r<3;++r)CHECK(chem_matrix_get(&m,r,0)==20*(r+1));
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void mutable_vector(size_t g)
{
    ChemMatrix m=mat(3,2,g,0);double v[]={5,10,15},scale=3;
    for(size_t r=0;r<3;++r)set(&m,r,0,r+1);
    ChemArray a[]={chem_array_const_matrix(&m),chem_array_vector(v,3,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(2,a,matrix_to_vector,&scale);call(f,2,a);
    CHECK(v[0]==3);CHECK(v[1]==6);CHECK(v[2]==9);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void reuse_vectors(size_t g)
{
    ChemMatrix m=mat(3,2,g,0);double v[]={1,2,3},w[]={5,6,7},scale=10;
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(2,a,vector_to_matrix,&scale);call(f,2,a);
    for(size_t r=0;r<3;++r)CHECK(chem_matrix_get(&m,r,0)==10*(r+1));
    a[1]=chem_array_const_vector(w,3,CHEM_MATRIX_REAL);call(f,2,a);
    for(size_t r=0;r<3;++r)CHECK(chem_matrix_get(&m,r,0)==50+10*r);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void wrong_vector(size_t g)
{
    ChemMatrix m=mat(3,2,g,0);double v[3]={0},w[5]={0},scale=1;
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(2,a,vector_to_matrix,&scale);call(f,2,a);
    /* Original wrong-size invocation is commented out; additionally check it. */
    a[1]=chem_array_const_vector(w,5,CHEM_MATRIX_REAL);CHECK(chem_array_function_call(f,2,a)==CHEM_DIMENSION);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void mixed_test(size_t g)
{
    ChemMatrix m=mat(4,3,g,0);double v[4];
    for(size_t r=0;r<4;++r){set(&m,r,0,r+1);set(&m,r,1,10*(r+1));v[r]=100*(r+1);}
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_vector(v,4,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(2,a,mixed,NULL);call(f,2,a);
    for(size_t r=0;r<4;++r)CHECK(chem_matrix_get(&m,r,2)==111*(r+1));
    chem_array_function_free(f);chem_matrix_free(&m);
}
static void const_signature(size_t g)
{
    ChemMatrix m=mat(3,2,g,0);const double v[3]={0};double scale=2;
    ChemArray a[]={chem_array_matrix(&m),chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(2,a,vector_to_matrix,&scale);
    ChemStatus (*wrapped)(const ChemArrayFunction *,size_t,const ChemArray *)=chem_array_function_call;
    CHECK(wrapped(f,2,a)==CHEM_OK);chem_array_function_free(f);chem_matrix_free(&m);
}
static ChemStatus fill_views(const ChemArray *a,size_t n,void *ctx)
{
    int mode=*(int *)ctx;
    if(mode==0)return chem_view_fill(column(a[0],1),3.2);
    if(n==2)return chem_view_fill(column(a[1],0),3.2);
    ChemMatrix owner={0};ChemView t=temporary(a[0],&owner);
    CHECK(chem_view_fill(t,9.9)==CHEM_OK);unary(column(a[0],0),t,1);
    chem_matrix_free(&owner);return CHEM_OK;
}
static void fill_test(size_t g)
{
    ChemMatrix m=mat(3,2,g,0);double vec[3]={0};int mode=0;
    ChemArray a[]={chem_array_matrix(&m),chem_array_vector(vec,3,CHEM_MATRIX_REAL)};
    ChemArrayFunction *f=function(1,a,fill_views,&mode);call(f,1,a);
    for(size_t r=0;r<3;++r){CHECK(chem_matrix_get(&m,r,0)==0);CHECK(chem_matrix_get(&m,r,1)==3.2);}
    chem_array_function_free(f);mode=1;a[0]=chem_array_const_matrix(&m);
    f=function(2,a,fill_views,&mode);call(f,2,a);for(size_t r=0;r<3;++r)CHECK(vec[r]==3.2);
    chem_array_function_free(f);mode=2;a[0]=chem_array_matrix(&m);
    f=function(1,a,fill_views,&mode);call(f,1,a);for(size_t r=0;r<3;++r)CHECK(chem_matrix_get(&m,r,0)==9.9);
    chem_array_function_free(f);chem_matrix_free(&m);
}
static ChemStatus copy_views(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;int mode=*(int *)ctx;
    if(mode==0)return chem_view_copy(column(a[0],1),read_column(a[1],0));
    if(mode==1)return chem_view_copy(column(a[1],0),read_column(a[0],1));
    if(mode==2)return chem_view_copy(column(a[0],0),column(a[0],1));
    if(mode==3)return chem_view_copy(column(a[0],0),read_column(a[0],1));
    ChemMatrix owner={0};ChemView t=temporary(a[0],&owner);
    CHECK(chem_view_copy(t,read_column(a[0],1))==CHEM_OK);CHECK(chem_view_fill(column(a[0],1),0)==CHEM_OK);
    unary(column(a[0],1),t,1);chem_matrix_free(&owner);return CHEM_OK;
}
static void copy_test(size_t g)
{
    ChemMatrix m=mat(3,2,g,0);double vec[]={3.2,4.2,1.3},other[3]={0};
    for(int mode=0;mode<5;++mode){
        ChemArray a[]={mode==1?chem_array_const_matrix(&m):chem_array_matrix(&m),
            mode==1?chem_array_vector(other,3,CHEM_MATRIX_REAL):chem_array_const_vector(vec,3,CHEM_MATRIX_REAL)};
        if(mode==3)for(size_t r=0;r<3;++r)set(&m,r,0,0);
        size_t n=mode<2?2:1;ChemArrayFunction *f=function(n,a,copy_views,&mode);call(f,n,a);
        for(size_t r=0;r<3;++r){
            if(mode==0){CHECK(chem_matrix_get(&m,r,0)==0);CHECK(chem_matrix_get(&m,r,1)==vec[r]);}
            if(mode==1)CHECK(other[r]==vec[r]);
            if(mode==2||mode==3)CHECK(chem_matrix_get(&m,r,0)==vec[r]);
            if(mode==2||mode==4)CHECK(chem_matrix_get(&m,r,1)==vec[r]);
        }
        chem_array_function_free(f);
    }
    chem_matrix_free(&m);
}
typedef struct {int mode;double accumulator;} Reduction;
static double reduce_value(double value,double acc,void *ctx)
{
    int mode=*(int *)ctx;
    if(mode==0)return acc+value*value;
    if(mode==1)return value>acc?value:acc;
    if(mode==2)return acc!=0||value>10;
    if(mode==3)return acc!=0&&isfinite(value);
    return acc+1;
}
static ChemStatus reductions(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;Reduction *r=ctx;
    for(size_t c=0;c<a[0].columns;++c)
        CHECK(chem_view_reduce(read_column(a[0],c),r->mode==4,reduce_value,&r->mode,&r->accumulator)==CHEM_OK);
    return CHEM_OK;
}
static void reduce_test(size_t g,int mode)
{
    ChemMatrix m=mat(3,mode==4?1:2,g,mode==2||mode==3?1:0);
    if(mode==0)for(size_t r=0;r<3;++r)for(size_t c=0;c<2;++c)set(&m,r,c,r*2+c+1);
    if(mode==1){double vals[3][2]={{1,8},{5,2},{3,4}};for(size_t r=0;r<3;++r)for(size_t c=0;c<2;++c)set(&m,r,c,vals[r][c]);}
    Reduction reduction={mode,mode==1?-DBL_MAX:mode==3?1:0};
    ChemArray a=chem_array_const_matrix(&m);ChemArrayFunction *f=function(1,&a,reductions,&reduction);call(f,1,&a);
    CHECK(reduction.accumulator==(mode==0?91:mode==1?8:mode==2?0:mode==3?1:3));
    if(mode==2||mode==3){
        set(&m,mode==2?2:1,mode==2?1:0,mode==2?42:NAN);reduction.accumulator=mode==3?1:0;
        call(f,1,&a);CHECK(reduction.accumulator==(mode==2?1:0));
    }
    chem_array_function_free(f);chem_matrix_free(&m);
}
typedef struct { ChemView indices,values; } Captured;
static void add_capture(double *v,size_t n,void *ctx)
{ (void)n;v[0]+=*(double *)ctx; }
static ChemStatus capture_vector(const ChemArray *a,size_t n,void *ctx)
{
    (void)n;Captured *capture=ctx;ChemView col0=column(a[0],0),col1=column(a[0],1);
    for(size_t i=0;i<capture->indices.length;++i){
        size_t index=(size_t)chem_view_get(capture->indices,i);double value=chem_view_get(capture->values,index);
        CHECK(chem_views_for_each(1,&col0,add_capture,&value)==CHEM_OK);
    }
    double value=chem_view_get(capture->values,(size_t)chem_view_get(capture->indices,1));
    CHECK(chem_views_for_each(1,&col1,add_capture,&value)==CHEM_OK);return CHEM_OK;
}
static void capture_test(size_t g)
{
    ChemMatrix m=mat(3,3,g,0);int indices[]={2,0,1};const int values[]={10,20,30};
    Captured ctx={column(chem_array_vector(indices,3,CHEM_MATRIX_INT),0),column(chem_array_const_vector(values,3,CHEM_MATRIX_INT),0)};
    ChemArray a=chem_array_matrix(&m);ChemArrayFunction *f=function(1,&a,capture_vector,&ctx);call(f,1,&a);
    for(size_t r=0;r<3;++r){CHECK(chem_matrix_get(&m,r,0)==60);CHECK(chem_matrix_get(&m,r,1)==10);CHECK(chem_matrix_get(&m,r,2)==0);}
    chem_array_function_free(f);chem_matrix_free(&m);
}
#define BOTH(name,expr) do { \
    RUN("matrix:Matrix." name,{size_t g=0;expr;}); \
    RUN("vector_matrix:VectorMatrix." name,{for(size_t g=1;g<=4;++g){expr;}}); \
} while(0)
void test_matrix_functions(void)
{
    BOTH("ArrayFunction",array_test(g,0));
    BOTH("MultiMatrixArrayFunction",multi_test(g));
    BOTH("MismatchedRowDimensions",mismatch_rows(g));
    BOTH("MismatchedColumnDimensions",mismatch_columns(g));
    BOTH("WrongMatrixDimensions",wrong_dimensions(g,0));
    BOTH("MultipleTemporaries",temporary_test(g));
    BOTH("ColumnViewReuse",reuse_columns_test(g));
    BOTH("FunctionReusability",reuse_test(g));
    BOTH("ConstMatrixFunction",const_test(g));
    BOTH("EmptyMatrixFunction",empty_test(g));
    BOTH("MultiMatrixDifferentRowsFromCreation",flexible_multi(g));
    BOTH("MatrixVectorDifferentRowsFromCreation",flexible_vector(g));
    BOTH("MismatchedRowsAtInvocation",row_invocation(g,0));
    BOTH("MultipleMatricesMismatchedRowsAtInvocation",row_invocation(g,1));
    BOTH("WrongColumnCountAtInvocation",wrong_dimensions(g,1));
    BOTH("VectorInMatrixFunction",array_test(g,1));
    BOTH("VectorTooSmall",vector_shape_test(g,0));
    BOTH("VectorTooLarge",vector_shape_test(g,1));
    BOTH("EmptyVectorNonEmptyMatrix",vector_shape_test(g,2));
    BOTH("NonEmptyVectorEmptyMatrix",vector_shape_test(g,3));
    BOTH("EmptyVectorEmptyMatrix",vector_shape_test(g,4));
    BOTH("MultipleVectorsDifferentSizes",vector_shape_test(g,5));
    BOTH("MultipleVectorsSameSize",same_vectors(g));
    BOTH("MultipleMatricesOneVector",matrices_one_vector(g));
    BOTH("MultipleMatricesDifferentRowsVector",vector_shape_test(g,6));
    BOTH("VectorSizeMatchesOneMatrixOnly",vector_shape_test(g,7));
    BOTH("ConstVector",const_vector(g,0));
    BOTH("MutableVector",mutable_vector(g));
    BOTH("FunctionReusabilityWithVectors",reuse_vectors(g));
    BOTH("FunctionInvocationWithWrongSizedVector",wrong_vector(g));
    BOTH("MixedVectorColumnViewRowVariable",mixed_test(g));
    BOTH("IntegerVector",const_vector(g,1));
    BOTH("FunctionWithConstSignature",const_signature(g));
    BOTH("TestFill",fill_test(g));
    BOTH("TestCopy",copy_test(g));
    BOTH("ReduceSum",reduce_test(g,0));
    BOTH("ReduceMax",reduce_test(g,1));
    BOTH("ReduceLOr",reduce_test(g,2));
    BOTH("ReduceLAnd",reduce_test(g,3));
    BOTH("ReduceStrict",reduce_test(g,4));
    BOTH("VectorCapture",capture_test(g));
}
