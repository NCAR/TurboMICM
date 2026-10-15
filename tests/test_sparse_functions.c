/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Ports the function, vector and view fixtures in test_sparse_matrix_policy.hpp.
 * Small operation sequences below are C closures for the original lambdas. */
#include "test.h"
#include "chem_sparse_matrix.h"
static const ChemSparseElement chain[]={{0,1},{1,2},{2,3},{3,3}};
static const ChemSparseElement diagonal3[]={{0,0},{1,1},{2,2}};
static ChemSparseMatrix sparse(size_t dim,size_t blocks,size_t group,ChemSparseOrdering order,size_t count,const ChemSparseElement *e,double initial)
{
    ChemSparseBuilder *b=NULL;ChemSparseMatrix m={0};CHECK(chem_sparse_builder_create(dim,&b)==CHEM_OK);
    for(size_t i=0;i<count;++i)CHECK(chem_sparse_builder_add(b,e[i].row,e[i].column)==CHEM_OK);
    CHECK(chem_sparse_builder_blocks(b,blocks)==CHEM_OK);chem_sparse_builder_initial(b,initial);
    CHECK(chem_sparse_matrix_build(b,group,order,CHEM_MATRIX_REAL,&m)==CHEM_OK);chem_sparse_builder_free(b);return m;
}
static void set(ChemSparseMatrix *m,size_t b,size_t r,size_t c,double v)
{CHECK(chem_sparse_matrix_set(m,b,r,c,v)==CHEM_OK);}
static size_t idx(const ChemSparseMatrix *m,size_t r,size_t c)
{size_t i=0;CHECK(chem_sparse_matrix_index(m,0,r,c,&i)==CHEM_OK);return i;}
static ChemArray mut(ChemSparseMatrix *m){return chem_array_sparse_matrix(m);}
static ChemArray ro(const ChemSparseMatrix *m){return chem_array_const_sparse_matrix(m);}
/* Arguments 4 and 5 are block temporaries; all other references are arguments
 * to the reusable function. Sparse references carry original raw vector indices. */
typedef struct {size_t argument,index;} Ref;
#define R(a,i) ((Ref){a,i})
typedef struct {int op;size_t n;Ref dst,src[4];double scale,offset;} Operation;
typedef struct {size_t count;Operation operations[8];} Program;
enum {SUM,PRODUCT,DIFFERENCE,COPY,FILL,CONST_COPY};
static ChemView view(ChemArray a,Ref ref)
{
    ChemView v={0};CHECK((a.kind==CHEM_ARRAY_SPARSE?chem_sparse_block_view(a,ref.index,&v):chem_array_column(a,ref.index,&v))==CHEM_OK);return v;
}
static ChemView cached_view(const ChemArray *a,Ref ref,Ref *keys,ChemView *views,size_t *count)
{
    for(size_t i=0;i<*count;++i)if(keys[i].argument==ref.argument&&keys[i].index==ref.index)return views[i];
    CHECK(*count<48);keys[*count]=ref;views[*count]=view(a[ref.argument],ref);return views[(*count)++];
}
static void operation_row(double *v,size_t n,void *context)
{
    Operation *op=context;double result=v[0];
    for(size_t i=1;i+1<n;++i){if(op->op==PRODUCT)result*=v[i];else if(op->op==DIFFERENCE)result-=v[i];else result+=v[i];}
    v[n-1]=result*op->scale+op->offset;
}
static ChemStatus evaluate(const ChemArray *args,size_t count,void *context)
{
    Program *p=context;ChemArray a[6]={{0}};CHECK(count<=4);for(size_t i=0;i<count;++i)a[i]=args[i];
    ChemMatrix tmp[2]={{0}};
    for(size_t i=0;i<2;++i){CHECK(chem_matrix_create(a[0].rows,1,a[0].group,CHEM_MATRIX_REAL,0,tmp+i)==CHEM_OK);a[4+i]=chem_array_matrix(tmp+i);}
    Ref keys[48];ChemView cached[48];size_t cache_count=0;
    for(size_t i=0;i<p->count;++i){
        Operation *op=p->operations+i;ChemView dst=cached_view(a,op->dst,keys,cached,&cache_count);
        if(op->op==FILL)CHECK(chem_view_fill(dst,op->offset)==CHEM_OK);
        else if(op->op==COPY||op->op==CONST_COPY){
            ChemView source=cached_view(a,op->src[0],keys,cached,&cache_count);
            CHECK(chem_view_copy(dst,op->op==CONST_COPY?chem_view_const(source):source)==CHEM_OK);
        }
        else {ChemView v[5];for(size_t k=0;k<op->n;++k)v[k]=chem_view_const(cached_view(a,op->src[k],keys,cached,&cache_count));v[op->n]=dst;
            CHECK(chem_views_for_each(op->n+1,v,operation_row,op)==CHEM_OK);}
    }
    for(size_t i=0;i<2;++i)chem_matrix_free(tmp+i);
    return CHEM_OK;
}
static void append(Program *p,int op,Ref dst,size_t n,const Ref *src,double scale,double offset)
{
    CHECK(p->count<8);Operation *o=p->operations+p->count++;*o=(Operation){.op=op,.n=n,.dst=dst,.scale=scale,.offset=offset};
    for(size_t i=0;i<n;++i)o->src[i]=src[i];
}
static void sum(Program *p,Ref dst,size_t n,const Ref *src,double scale){append(p,SUM,dst,n,src,scale,0);}
static void copy(Program *p,Ref dst,Ref src){append(p,COPY,dst,1,&src,1,0);}
static void fill(Program *p,Ref dst,double v){append(p,FILL,dst,0,NULL,1,v);}
static ChemArrayFunction *function(size_t n,const ChemArray *a,Program *p)
{ChemArrayFunction *f=NULL;CHECK(chem_array_function_create(n,a,evaluate,p,&f)==CHEM_OK);return f;}
static void call(ChemArrayFunction *f,size_t n,const ChemArray *a){CHECK(chem_array_function_call(f,n,a)==CHEM_OK);}
static void array_function(size_t g,ChemSparseOrdering order)
{
    const ChemSparseElement e[]={{0,0},{1,1},{2,2},{2,3},{3,3}};
    ChemSparseMatrix m=sparse(4,3,g,order,5,e,1),other=sparse(4,3,g,order,5,e,-1);
    set(&m,0,0,0,1);set(&m,1,1,1,2);set(&m,2,2,2,3);set(&m,2,2,3,4);set(&m,2,3,3,5);
    Ref inputs[]={R(0,idx(&m,0,0)),R(0,idx(&m,1,1)),R(0,idx(&m,2,2)),R(0,idx(&m,2,3))};
    Program p={0};sum(&p,R(4,0),4,inputs,1);Ref t=R(4,0);sum(&p,inputs[3],1,&t,2);
    ChemArray a=mut(&m);ChemArrayFunction *f=function(1,&a,&p);call(f,1,&a);
    const double expected[]={8,10,18};
    for(size_t b=0;b<3;++b){CHECK(chem_sparse_matrix_get(&m,b,2,3)==expected[b]);CHECK(chem_sparse_matrix_get(&m,b,0,0)==1);
        CHECK(chem_sparse_matrix_get(&m,b,1,1)==(b==1?2:1));CHECK(chem_sparse_matrix_get(&m,b,2,2)==(b==2?3:1));CHECK(chem_sparse_matrix_get(&m,b,3,3)==(b==2?5:1));}
    a=mut(&other);call(f,1,&a);
    for(size_t b=0;b<3;++b)for(size_t i=0;i<5;++i)CHECK(chem_sparse_matrix_get(&other,b,e[i].row,e[i].column)==(i==3?-8:-1));
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&other);
}
static void multi_function(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(3,3,g,order,2,chain,0),b=sparse(3,3,g,order,3,diagonal3,0);
    for(size_t i=0;i<3;++i){set(&m,i,0,1,i+10);set(&m,i,1,2,2*i+20);set(&b,i,0,0,4*i);set(&b,i,1,1,3*i);set(&b,i,2,2,5*i);}
    Program p={0};Ref inputs[]={R(0,idx(&m,0,1)),R(1,idx(&b,2,2))};sum(&p,R(4,0),2,inputs,1);copy(&p,R(0,idx(&m,1,2)),R(4,0));
    ChemArray a[]={mut(&m),ro(&b)};ChemArrayFunction *f=function(2,a,&p);call(f,2,a);
    for(size_t i=0;i<3;++i){CHECK(chem_sparse_matrix_get(&m,i,1,2)==10+6*i);CHECK(chem_sparse_matrix_get(&m,i,0,1)==10+i);CHECK(chem_sparse_matrix_get(&b,i,0,0)==4*i);}
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);
}
static void mismatch_blocks(size_t g,ChemSparseOrdering order,int which)
{
    const ChemSparseElement special[]={{0,1},{1,1}};const ChemSparseElement *e=which?chain:special;size_t count=which==2?3:2,dim=which==2?4:3;
    ChemSparseMatrix m=sparse(dim,3,g,order,count,e,0),b=sparse(dim,3,g,order,count,e,0),bad=sparse(dim,4,g,order,count,e,0);
    Program p={0};Ref inputs[]={R(0,idx(&m,e[0].row,e[0].column)),R(1,idx(&m,e[1].row,e[1].column)),R(2,0)};
    if(which==0){inputs[1]=R(1,idx(&m,0,1));sum(&p,R(0,idx(&m,1,1)),2,inputs,1);}
    if(which==1){inputs[0]=R(1,idx(&m,0,1));sum(&p,R(0,idx(&m,1,2)),1,inputs,2);}
    if(which==2){inputs[2]=R(2,idx(&m,2,3));sum(&p,R(0,idx(&m,2,3)),3,inputs,1);}
    ChemArray a[]={mut(&m),ro(which==2?&b:&m),ro(&b)};size_t n=which==2?3:2;
    ChemArrayFunction *f=function(n,a,&p);call(f,n,a);
    a[n-1]=ro(&bad);CHECK(chem_array_function_call(f,n,a)==CHEM_DIMENSION);
    if(which){a[n-1]=ro(&b);a[0]=mut(&bad);CHECK(chem_array_function_call(f,n,a)==CHEM_DIMENSION);}
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);chem_sparse_matrix_free(&bad);
}
static void missing_element(size_t g,ChemSparseOrdering order)
{
    const ChemSparseElement e[]={{0,1},{1,1},{2,2}};ChemSparseMatrix m=sparse(4,3,g,order,3,e,0);size_t i;
    CHECK(chem_sparse_matrix_index(&m,0,0,1,&i)==CHEM_OK);CHECK(chem_sparse_matrix_index(&m,0,3,3,&i)==CHEM_NOT_FOUND);
    Program p={0};ChemArray a=mut(&m);ChemArrayFunction *f=function(1,&a,&p);chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void flexible_blocks(size_t g,ChemSparseOrdering order,int vector)
{
    size_t dim=vector?3:4,count=vector?2:3;
    ChemSparseMatrix m=sparse(dim,3,g,order,count,chain,0),b=sparse(dim,3,g,order,count,chain,0);
    double v[]={1,2,3},v5[]={10,20,30,40,50};
    for(size_t i=0;i<3;++i){set(&m,i,0,1,i+1);set(&m,i,1,2,i+10);if(!vector)set(&b,i,2,3,i+100);}
    Ref inputs[]={R(0,idx(&m,0,1)),vector?R(1,0):R(0,idx(&m,1,2)),R(1,vector?0:idx(&b,2,3))};
    Program p={0};sum(&p,R(4,0),vector?2:3,inputs,1);copy(&p,R(0,idx(&m,vector?1:2,vector?2:3)),R(4,0));
    ChemArray a[]={mut(&m),vector?chem_array_const_vector(v,3,CHEM_MATRIX_REAL):ro(&b)};
    ChemArrayFunction *f=function(2,a,&p);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);
    size_t blocks=vector?5:4;m=sparse(dim,blocks,g,order,count,chain,0);b=sparse(dim,blocks,g,order,count,chain,0);
    for(size_t i=0;i<blocks;++i){set(&m,i,0,1,i+1);if(!vector){set(&m,i,1,2,i+10);set(&b,i,2,3,i+100);}}
    a[0]=mut(&m);a[1]=vector?chem_array_const_vector(v5,5,CHEM_MATRIX_REAL):ro(&b);call(f,2,a);
    for(size_t i=0;i<blocks;++i)CHECK(chem_sparse_matrix_get(&m,i,vector?1:2,vector?2:3)==(vector?11*(i+1):111+3*i));
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);
}
static void temporaries(size_t g,ChemSparseOrdering order)
{
    const ChemSparseElement e[]={{0,1},{1,2},{2,3},{3,4},{4,4}};ChemSparseMatrix m=sparse(5,4,g,order,5,e,0);
    for(size_t i=0;i<4;++i){set(&m,i,0,1,i+1);set(&m,i,1,2,10*(i+1));}
    Ref input[]={R(0,idx(&m,0,1)),R(0,idx(&m,1,2))},temps[]={R(4,0),R(5,0)};Program p={0};
    append(&p,PRODUCT,temps[0],2,input,1,0);sum(&p,temps[1],2,input,1);
    sum(&p,R(0,idx(&m,2,3)),2,temps,1);append(&p,DIFFERENCE,R(0,idx(&m,3,4)),2,temps,1,0);append(&p,PRODUCT,R(0,idx(&m,4,4)),2,temps,1,0);
    ChemArray a=mut(&m);ChemArrayFunction *f=function(1,&a,&p);call(f,1,&a);
    for(size_t i=0;i<4;++i){double product=10*(i+1)*(i+1),s=11*(i+1);CHECK(chem_sparse_matrix_get(&m,i,2,3)==product+s);
        CHECK(chem_sparse_matrix_get(&m,i,3,4)==product-s);CHECK(chem_sparse_matrix_get(&m,i,4,4)==product*s);}
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void block_reuse(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(4,3,g,order,4,chain,0);for(size_t i=0;i<3;++i)set(&m,i,0,1,i+1);
    Ref v[4];for(size_t i=0;i<4;++i)v[i]=R(0,idx(&m,chain[i].row,chain[i].column));
    Program p={0};sum(&p,v[1],1,v,2);sum(&p,v[2],2,v,1);Ref in[]={v[2],v[1]};append(&p,PRODUCT,v[3],2,in,1,0);
    ChemArray a=mut(&m);ChemArrayFunction *f=function(1,&a,&p);call(f,1,&a);
    for(size_t i=0;i<3;++i){CHECK(chem_sparse_matrix_get(&m,i,1,2)==2*(i+1));CHECK(chem_sparse_matrix_get(&m,i,2,3)==3*(i+1));CHECK(chem_sparse_matrix_get(&m,i,3,3)==6*(i+1)*(i+1));}
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void reuse(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(3,2,g,order,3,diagonal3,0);Ref in[3];for(size_t i=0;i<3;++i)in[i]=R(0,idx(&m,i,i));
    Program p={0};sum(&p,R(4,0),3,in,1);Ref t=R(4,0);sum(&p,in[2],1,&t,2);ChemArray a=mut(&m);ChemArrayFunction *f=function(1,&a,&p);
    for(size_t b=0;b<2;++b)for(size_t i=0;i<3;++i)set(&m,b,i,i,b+i);
    call(f,1,&a);CHECK(chem_sparse_matrix_get(&m,0,2,2)==6);CHECK(chem_sparse_matrix_get(&m,1,2,2)==12);
    chem_sparse_matrix_free(&m);m=sparse(3,2,g,order,3,diagonal3,5);a=mut(&m);call(f,1,&a);
    CHECK(chem_sparse_matrix_get(&m,0,2,2)==30);CHECK(chem_sparse_matrix_get(&m,1,2,2)==30);
    chem_sparse_matrix_free(&m);m=sparse(3,2,g,order,3,diagonal3,0);for(size_t b=0;b<2;++b)set(&m,b,0,0,10*b);a=mut(&m);call(f,1,&a);
    CHECK(chem_sparse_matrix_get(&m,0,2,2)==0);CHECK(chem_sparse_matrix_get(&m,1,2,2)==20);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void different_structure(size_t g,ChemSparseOrdering order)
{
    const ChemSparseElement ea[]={{0,1},{1,1},{2,2}},eb[]={{0,0},{1,2},{2,3},{3,4}};
    ChemSparseMatrix m=sparse(3,4,g,order,3,ea,0),b=sparse(5,4,g,order,4,eb,0);
    for(size_t i=0;i<4;++i){set(&m,i,0,1,2*i+1);set(&m,i,1,1,3*i+2);set(&m,i,2,2,5*i+3);for(size_t k=0;k<4;++k)set(&b,i,eb[k].row,eb[k].column,i+10*(k+1));}
    Program p={0};Ref in[]={R(0,idx(&m,0,1)),R(0,idx(&m,1,1)),R(1,idx(&b,0,0))};sum(&p,R(4,0),3,in,1);copy(&p,R(0,idx(&m,2,2)),R(4,0));
    ChemArray a[]={mut(&m),ro(&b)};ChemArrayFunction *f=function(2,a,&p);call(f,2,a);
    for(size_t i=0;i<4;++i){CHECK(chem_sparse_matrix_get(&m,i,2,2)==13+6*i);CHECK(chem_sparse_matrix_get(&m,i,0,1)==1+2*i);CHECK(chem_sparse_matrix_get(&b,i,eb[i].row,eb[i].column)==10+11*i);}
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);
}
static void dense_function(size_t g,ChemSparseOrdering order,size_t dense_group)
{
    ChemSparseMatrix m=sparse(4,3,g,order,4,chain,0);ChemMatrix d={0};CHECK(chem_matrix_create(3,4,dense_group,CHEM_MATRIX_REAL,0,&d)==CHEM_OK);
    for(size_t b=0;b<3;++b){set(&m,b,0,1,b+1);set(&m,b,1,2,b+10);set(&m,b,2,3,b+20);set(&m,b,3,3,b+30);
        for(size_t c=0;c<4;++c)CHECK(chem_matrix_set(&d,b,c,100*b+10*c)==CHEM_OK);}
    Program p={0};Ref in[]={R(0,idx(&m,0,1)),R(0,idx(&m,1,2)),R(1,0),R(1,1)};sum(&p,R(4,0),4,in,1);copy(&p,R(0,idx(&m,3,3)),R(4,0));
    ChemArray a[]={mut(&m),chem_array_const_matrix(&d)};ChemArrayFunction *f=function(2,a,&p);call(f,2,a);
    for(size_t b=0;b<3;++b)CHECK(chem_sparse_matrix_get(&m,b,3,3)==21+202*b);
    CHECK(chem_sparse_matrix_get(&m,0,0,1)==1);CHECK(chem_sparse_matrix_get(&m,1,1,2)==11);CHECK(chem_sparse_matrix_get(&m,2,2,3)==22);
    CHECK(chem_matrix_get(&d,0,0)==0);CHECK(chem_matrix_get(&d,1,2)==120);CHECK(chem_matrix_get(&d,2,3)==230);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_matrix_free(&d);
}
static void incompatible(size_t g,ChemSparseOrdering order,size_t other_group,int other_sparse)
{
    ChemSparseMatrix m=sparse(4,3,g,order,2,chain,0),b={0};ChemMatrix d={0};
    if(other_sparse)b=sparse(4,3,other_group,order,2,chain,0);else CHECK(chem_matrix_create(3,4,other_group,CHEM_MATRIX_REAL,0,&d)==CHEM_OK);
    ChemArray a[]={ro(&m),other_sparse?ro(&b):chem_array_const_matrix(&d)};Program p={0};ChemArrayFunction *f=NULL;
    CHECK(chem_array_function_create(2,a,evaluate,&p,&f)==CHEM_DIMENSION);chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);chem_matrix_free(&d);
}
static void const_function(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(4,3,g,order,4,chain,0);for(size_t b=0;b<3;++b){set(&m,b,0,1,b+1);set(&m,b,1,2,b+10);}
    Ref in[]={R(0,idx(&m,0,1)),R(0,idx(&m,1,2))};Program p={0};sum(&p,R(4,0),2,in,1);
    ChemArray a=ro(&m);ChemArrayFunction *f=function(1,&a,&p);call(f,1,&a);
    ChemView v=view(a,R(0,idx(&m,2,3)));double total=0;for(size_t i=0;i<v.length;++i)total+=chem_view_get(v,i);CHECK(total==0);
    CHECK(chem_sparse_matrix_get(&m,0,0,1)==1);CHECK(chem_sparse_matrix_get(&m,1,1,2)==11);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void empty_function(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(4,3,g,order,0,NULL,0);Program p={0};ChemArray a=ro(&m);ChemArrayFunction *f=function(1,&a,&p);call(f,1,&a);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);m=sparse(4,0,g,order,1,chain,0);a=ro(&m);f=function(1,&a,&p);call(f,1,&a);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static ChemStatus should_not_execute(const ChemArray *a,size_t n,void *ctx)
{(void)a;(void)n;(void)ctx;return CHEM_INVALID;}
static void vector_shape(size_t g,ChemSparseOrdering order,int which)
{
    size_t blocks=which==3||which==4?0:3,len=which==0?2:which==1?4:which==3?3:0;
    if(which>=5)len=3;
    ChemSparseMatrix m=sparse(3,blocks,g,order,1,chain,0),b=sparse(3,which==6?4:3,g,order,1,chain,0);
    double v[]={1,2,3,4},w[]={4,5};ChemArray a[4]={mut(&m),chem_array_const_vector(v,len,CHEM_MATRIX_REAL),{0},{0}};
    Program p={0};size_t n=2;Ref src=R(1,0);
    if(which==5){a[2]=chem_array_const_vector(w,2,CHEM_MATRIX_REAL);n=3;}
    if(which==6||which==7){a[1]=ro(&b);a[2]=chem_array_const_vector(v,3,CHEM_MATRIX_REAL);src=R(2,0);n=3;}
    if(which==7){a[3]=chem_array_const_vector(w,2,CHEM_MATRIX_REAL);n=4;}
    copy(&p,R(0,0),src);ChemArrayFunction *f=NULL;
    if(which==4)CHECK(chem_array_function_create(n,a,should_not_execute,NULL,&f)==CHEM_OK);
    else f=function(n,a,&p);
    if(which==3){size_t index;CHECK(chem_sparse_matrix_index(&m,0,0,1,&index)==CHEM_DIMENSION);}
    CHECK(chem_array_function_call(f,n,a)==(which==4?CHEM_OK:CHEM_DIMENSION));
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);
}
static void vector_values(size_t g,ChemSparseOrdering order,int which)
{
    /* 0 vector-in-function; 1 const; 2 mutable; 3 integer; 4 const signature. */
    ChemSparseElement elem=which==4?(ChemSparseElement){0,0}:chain[0];
    ChemSparseMatrix m=sparse(which==0?4:which==4?2:3,3,g,order,which==0?3:1,which==0?chain:&elem,0);
    double v[]={5,10,15};int iv[]={5,10,15};if(which==4)for(size_t b=0;b<3;++b)v[b]=b+1;
    if(which==0)for(size_t b=0;b<3;++b)for(size_t i=0;i<3;++i)set(&m,b,chain[i].row,chain[i].column,b+i+1);
    if(which==2)for(size_t b=0;b<3;++b)set(&m,b,0,1,b+1);
    ChemArray a[]={which==2?ro(&m):mut(&m),which==2?chem_array_vector(v,3,CHEM_MATRIX_REAL):
        chem_array_const_vector(which==3?(const void *)iv:(const void *)v,3,which==3?CHEM_MATRIX_INT:CHEM_MATRIX_REAL)};
    Ref src=which==2?R(0,idx(&m,0,1)):R(1,0),dst=which==2?R(1,0):R(0,idx(&m,elem.row,elem.column));
    Program p={0};sum(&p,dst,1,&src,which==1||which==2?3:which==3?1.5:2);
    ChemArrayFunction *f=function(2,a,&p);
    if(which==4){ChemStatus (*wrapped)(const ChemArrayFunction *,size_t,const ChemArray *)=chem_array_function_call;CHECK(wrapped(f,2,a)==CHEM_OK);}
    else call(f,2,a);
    for(size_t b=0;b<3;++b){double expected=(which==0?10:which==1?15:which==2?3:which==3?7.5:2)*(b+1);
        CHECK((which==2?v[b]:chem_sparse_matrix_get(&m,b,elem.row,elem.column))==expected);}
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void same_vectors(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(3,3,g,order,1,chain,0);double v[]={1,2,3},w[]={4,5,6};
    ChemArray a[]={mut(&m),chem_array_const_vector(v,3,CHEM_MATRIX_REAL),chem_array_const_vector(w,3,CHEM_MATRIX_REAL)};
    Program p={0};Ref in[]={R(1,0),R(2,0)};sum(&p,R(0,idx(&m,0,1)),2,in,1);ChemArrayFunction *f=function(3,a,&p);call(f,3,a);
    for(size_t b=0;b<3;++b)CHECK(chem_sparse_matrix_get(&m,b,0,1)==5+2*b);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void matrices_vector(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(3,3,g,order,1,chain,0),b=sparse(3,3,g,order,1,chain,0);double v[]={10,20,30};
    for(size_t i=0;i<3;++i)set(&m,i,0,1,i+1);
    Program p={0};Ref in[]={R(0,idx(&m,0,1)),R(2,0)};sum(&p,R(1,idx(&b,0,1)),2,in,1);
    ChemArray a[]={ro(&m),mut(&b),chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};ChemArrayFunction *f=function(3,a,&p);call(f,3,a);
    for(size_t i=0;i<3;++i)CHECK(chem_sparse_matrix_get(&b,i,0,1)==11*(i+1));
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);
}
static void reuse_vectors(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(3,3,g,order,1,chain,0),b=sparse(3,3,g,order,1,chain,0);double v[]={1,2,3},w[]={10,20,30};
    Program p={0};Ref src=R(1,0);sum(&p,R(0,idx(&m,0,1)),1,&src,2);
    ChemArray a[]={mut(&m),chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};ChemArrayFunction *f=function(2,a,&p);call(f,2,a);
    for(size_t i=0;i<3;++i)CHECK(chem_sparse_matrix_get(&m,i,0,1)==2*(i+1));
    a[0]=mut(&b);a[1]=chem_array_const_vector(w,3,CHEM_MATRIX_REAL);call(f,2,a);
    for(size_t i=0;i<3;++i)CHECK(chem_sparse_matrix_get(&b,i,0,1)==20*(i+1));
    chem_array_function_free(f);chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&b);
}
static void wrong_vector(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(3,3,g,order,1,chain,0);double v[]={1,2,3},w[]={1,2};
    ChemArray a[]={ro(&m),chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};Program p={0};ChemArrayFunction *f=function(2,a,&p);
    a[1]=chem_array_const_vector(w,2,CHEM_MATRIX_REAL);CHECK(chem_array_function_call(f,2,a)==CHEM_DIMENSION);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void mixed_vector(size_t g,ChemSparseOrdering order)
{
    ChemSparseMatrix m=sparse(4,3,g,order,3,chain,0);double v[]={100,200,300};
    for(size_t b=0;b<3;++b)for(size_t i=0;i<3;++i)set(&m,b,chain[i].row,chain[i].column,b+i+1);
    Program p={0};Ref in[]={R(0,idx(&m,0,1)),R(0,idx(&m,1,2)),R(1,0)};sum(&p,R(4,0),3,in,1);Ref t=R(4,0);sum(&p,R(0,idx(&m,2,3)),1,&t,2);
    ChemArray a[]={mut(&m),chem_array_const_vector(v,3,CHEM_MATRIX_REAL)};ChemArrayFunction *f=function(2,a,&p);call(f,2,a);
    for(size_t b=0;b<3;++b)CHECK(chem_sparse_matrix_get(&m,b,2,3)==206+204*b);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void indices(size_t g,ChemSparseOrdering order)
{
    const ChemSparseElement e[]={{0,0},{0,2},{1,0},{1,1},{2,1},{2,2}};ChemSparseMatrix m=sparse(3,3,g,order,6,e,0);
    size_t captured[6];for(size_t i=0;i<6;++i)captured[i]=idx(&m,e[i].row,e[i].column);
    const int v[]={0,1,2};Program p={0};Ref src=R(1,0);
    for(size_t i=0;i<6;++i)append(&p,SUM,R(0,captured[i]),1,&src,1,10*i);
    ChemArray a[]={mut(&m),chem_array_const_vector(v,3,CHEM_MATRIX_INT)};ChemArrayFunction *f=function(2,a,&p);call(f,2,a);
    for(size_t b=0;b<3;++b)for(size_t i=0;i<6;++i)CHECK(chem_sparse_matrix_get(&m,b,e[i].row,e[i].column)==b+10*i);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void fills(size_t g,ChemSparseOrdering order)
{
    const ChemSparseElement e[]={{0,1},{2,0}};ChemSparseMatrix m=sparse(3,3,g,order,2,e,0);size_t first=idx(&m,0,1),second=idx(&m,2,0);
    double vec[3]={0};ChemArray a[]={mut(&m),chem_array_vector(vec,3,CHEM_MATRIX_REAL)};
    Program p={0};fill(&p,R(0,first),3.2);ChemArrayFunction *f=function(1,a,&p);call(f,1,a);
    for(size_t b=0;b<3;++b){CHECK(chem_sparse_matrix_get(&m,b,0,1)==3.2);CHECK(chem_sparse_matrix_get(&m,b,2,0)==0);}chem_array_function_free(f);
    p=(Program){0};fill(&p,R(1,0),3.2);a[0]=ro(&m);f=function(2,a,&p);call(f,2,a);
    for(size_t b=0;b<3;++b)CHECK(vec[b]==3.2);
    chem_array_function_free(f);p=(Program){0};fill(&p,R(4,0),9.9);Ref t=R(4,0);sum(&p,R(0,second),1,&t,1);
    a[0]=mut(&m);f=function(1,a,&p);call(f,1,a);for(size_t b=0;b<3;++b)CHECK(chem_sparse_matrix_get(&m,b,2,0)==9.9);
    chem_array_function_free(f);chem_sparse_matrix_free(&m);
}
static void copies(size_t g,ChemSparseOrdering order)
{
    const ChemSparseElement e[]={{0,1},{2,0}},other_e={1,2};const double original[]={3.2,4.2,1.3};
    for(int mode=0;mode<8;++mode){
        ChemSparseMatrix m=sparse(3,3,g,order,2,e,0),dst={0};for(size_t b=0;b<3;++b)set(&m,b,0,1,original[b]);
        size_t first=idx(&m,0,1),second=idx(&m,2,0);ChemMatrix dense={0};double vec[]={10,20,30};Program p={0};
        ChemArray a[3]={mut(&m),{0},{0}};size_t n=1;
        if(mode==0){for(size_t b=0;b<3;++b)set(&m,b,0,1,0);a[1]=chem_array_const_vector(vec,3,CHEM_MATRIX_REAL);n=2;copy(&p,R(0,first),R(1,0));}
        if(mode==1){for(size_t b=0;b<3;++b)vec[b]=-1;a[0]=ro(&m);a[1]=chem_array_vector(vec,3,CHEM_MATRIX_REAL);n=2;copy(&p,R(1,0),R(0,first));}
        if(mode==2||mode==3){copy(&p,R(0,second),R(0,first));if(mode==3)p.operations[0].op=CONST_COPY;}
        if(mode==4){copy(&p,R(4,0),R(0,first));fill(&p,R(0,first),0);Ref t=R(4,0);sum(&p,R(0,first),1,&t,1);}
        if(mode==5||mode==6){CHECK(chem_matrix_create(3,2,g,CHEM_MATRIX_REAL,0,&dense)==CHEM_OK);n=2;
            if(mode==5){for(size_t b=0;b<3;++b)CHECK(chem_matrix_set(&dense,b,1,(b+1)*100)==CHEM_OK);
                a[1]=chem_array_const_matrix(&dense);Ref in=R(1,1);sum(&p,R(0,second),1,&in,1);}
            else {for(size_t b=0;b<3;++b)vec[b]=0;a[0]=ro(&m);a[1]=chem_array_vector(vec,3,CHEM_MATRIX_REAL);copy(&p,R(1,0),R(0,first));}}
        if(mode==7){dst=sparse(3,3,g,order,1,&other_e,0);a[0]=mut(&dst);a[1]=ro(&m);n=2;copy(&p,R(0,idx(&dst,1,2)),R(1,first));}
        ChemArrayFunction *f=function(n,a,&p);call(f,n,a);chem_array_function_free(f);
        if(mode==6){p=(Program){0};Ref in=R(1,0);sum(&p,R(0,0),1,&in,1);a[0]=chem_array_matrix(&dense);a[1]=chem_array_const_vector(vec,3,CHEM_MATRIX_REAL);
            f=function(2,a,&p);call(f,2,a);chem_array_function_free(f);}
        for(size_t b=0;b<3;++b){
            if(mode==0)CHECK(chem_sparse_matrix_get(&m,b,0,1)==10*(b+1));
            else if(mode==1)CHECK(vec[b]==original[b]);
            else if(mode==2||mode==3)CHECK(chem_sparse_matrix_get(&m,b,2,0)==original[b]);
            else if(mode==4)CHECK(chem_sparse_matrix_get(&m,b,0,1)==original[b]);
            else if(mode==5)CHECK(chem_sparse_matrix_get(&m,b,2,0)==100*(b+1));
            else if(mode==6)CHECK(chem_matrix_get(&dense,b,0)==original[b]);
            else CHECK(chem_sparse_matrix_get(&dst,b,1,2)==original[b]);
        }
        if(mode==0)CHECK(chem_sparse_matrix_get(&m,0,2,0)==0);
        if(mode==2)CHECK(chem_sparse_matrix_get(&m,0,0,1)==3.2);
        chem_sparse_matrix_free(&m);chem_sparse_matrix_free(&dst);chem_matrix_free(&dense);
    }
}
#define CASE(prefix,name,expr) RUN(prefix name,{if(vector){for(size_t g=1;g<=4;++g){expr;}}else{size_t g=0;expr;}})
#define COMMON(prefix) do { \
    CASE(prefix,"ArrayFunction",array_function(g,o)); \
    CASE(prefix,"MultiMatrixArrayFunction",multi_function(g,o)); \
    CASE(prefix,"MismatchedBlockDimensions",mismatch_blocks(g,o,0)); \
    CASE(prefix,"MismatchedElementDimensions",missing_element(g,o)); \
    CASE(prefix,"MultipleTemporaries",temporaries(g,o)); \
    CASE(prefix,"BlockViewReuse",block_reuse(g,o)); \
    CASE(prefix,"FunctionReusability",reuse(g,o)); \
    CASE(prefix,"TwoSparseMatricesDifferentStructure",different_structure(g,o)); \
    RUN(prefix "SparseAndDenseMatrixFunction",dense_function(vector?1:0,o,0)); \
    CASE(prefix,"ConstSparseMatrixFunction",const_function(g,o)); \
    CASE(prefix,"EmptySparseMatrixFunction",empty_function(g,o)); \
    CASE(prefix,"MultipleSparseMatricesDifferentBlocksFromCreation",flexible_blocks(g,o,0)); \
    CASE(prefix,"SparseMatrixVectorDifferentBlocksFromCreation",flexible_blocks(g,o,1)); \
    CASE(prefix,"MismatchedBlocksAtInvocation",mismatch_blocks(g,o,1)); \
    CASE(prefix,"MultipleSparseMatricesMismatchedBlocksAtInvocation",mismatch_blocks(g,o,2)); \
    CASE(prefix,"VectorInSparseMatrixFunction",vector_values(g,o,0)); \
    CASE(prefix,"VectorTooSmall",vector_shape(g,o,0)); \
    CASE(prefix,"VectorTooLarge",vector_shape(g,o,1)); \
    CASE(prefix,"EmptyVectorNonEmptySparseMatrix",vector_shape(g,o,2)); \
    CASE(prefix,"NonEmptyVectorEmptySparseMatrix",vector_shape(g,o,3)); \
    CASE(prefix,"EmptyVectorEmptySparseMatrix",vector_shape(g,o,4)); \
    CASE(prefix,"MultipleVectorsDifferentSizes",vector_shape(g,o,5)); \
    CASE(prefix,"MultipleVectorsSameSize",same_vectors(g,o)); \
    CASE(prefix,"MultipleSparseMatricesOneVector",matrices_vector(g,o)); \
    CASE(prefix,"MultipleSparseMatricesDifferentBlocksVector",vector_shape(g,o,6)); \
    CASE(prefix,"VectorSizeMatchesOneSparseMatrixOnly",vector_shape(g,o,7)); \
    CASE(prefix,"ConstVectorSparse",vector_values(g,o,1)); \
    CASE(prefix,"MutableVectorSparse",vector_values(g,o,2)); \
    CASE(prefix,"FunctionReusabilityWithVectorsSparse",reuse_vectors(g,o)); \
    CASE(prefix,"FunctionInvocationWithWrongSizedVectorSparse",wrong_vector(g,o)); \
    CASE(prefix,"MixedVectorBlockViewBlockVariable",mixed_vector(g,o)); \
    CASE(prefix,"IntegerVectorSparse",vector_values(g,o,3)); \
    CASE(prefix,"FunctionWithConstSignatureSparse",vector_values(g,o,4)); \
    CASE(prefix,"GetBlockViewByVectorIndex",indices(g,o)); \
    CASE(prefix,"TestFill",fills(g,o)); \
    CASE(prefix,"TestCopy",copies(g,o)); \
}while(0)
#define VECTOR_EXTRA(prefix) do { \
    RUN(prefix "SparseAndVectorMatrixFunction",{for(size_t g=1;g<=4;++g)dense_function(g,o,g);dense_function(10,o,10);}); \
    RUN(prefix "IncompatibleOrdering",{for(size_t g=2;g<=4;++g)incompatible(g,o,0,0);}); \
    RUN(prefix "IncompatibleVectorOrdering",{incompatible(2,o,1,0);incompatible(2,o,3,0);incompatible(4,o,2,0);incompatible(4,o,10,0);}); \
    RUN(prefix "IncompatibleSparseOrdering",{incompatible(1,o,2,1);incompatible(2,o,3,1);incompatible(3,o,4,1);incompatible(4,o,10,1);}); \
}while(0)
void test_sparse_functions(void)
{
    {int vector=0;ChemSparseOrdering o=CHEM_SPARSE_COLUMN;COMMON("sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix.");}
    {int vector=0;ChemSparseOrdering o=CHEM_SPARSE_ROW;COMMON("sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix.");}
    {int vector=1;ChemSparseOrdering o=CHEM_SPARSE_COLUMN;COMMON("sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.");VECTOR_EXTRA("sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix.");}
    {int vector=1;ChemSparseOrdering o=CHEM_SPARSE_ROW;COMMON("sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.");VECTOR_EXTRA("sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix.");}
}
