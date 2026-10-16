#include "test.h"
#include "chem_sparse.h"

typedef struct { size_t r,c; } pair;
static chem_sparse sparse(size_t n,size_t blocks,size_t lanes,int order,size_t count,const pair *elements,chem_element_type type,double initial)
{
    chem_sparse_builder b; chem_sparse s={0}; CHECK(chem_sparse_builder_create(&b,n)==CHEM_OK);
    b.blocks=blocks;b.initial=initial;
    for(size_t i=0;i<count;++i) CHECK(chem_sparse_builder_add(&b,elements[i].r,elements[i].c)==CHEM_OK);
    CHECK(chem_sparse_create(&s,&b,order ? CHEM_CSC : CHEM_CSR,lanes,type)==CHEM_OK);
    chem_sparse_builder_destroy(&b);return s;
}
static const pair basic[]={{0,1},{3,2},{2,3},{1,1},{2,1}};
static const pair chain[]={{0,1},{1,2},{2,3},{3,3}};
static const pair diagonal[]={{0,0},{1,1},{2,2}};
static void put(chem_sparse *s,size_t b,size_t r,size_t c,double x)
{ CHECK(chem_sparse_set(s,b,r,c,x)==CHEM_OK); }
static double get(const chem_sparse *s,size_t b,size_t r,size_t c)
{ double x=0;CHECK(chem_sparse_get(s,b,r,c,&x)==CHEM_OK);return x; }
static size_t id(const chem_sparse *s,size_t r,size_t c)
{ size_t k=CHEM_NONE;CHECK(chem_sparse_index(s,0,r,c,&k)==CHEM_OK);return k; }
static chem_array vec(size_t n,size_t v,chem_element_type type,const double *x)
{
    chem_array a;CHECK(chem_vector_create(&a,n,v,type,0)==CHEM_OK);
    for(size_t i=0;i<n;++i) CHECK(chem_array_set(&a,i,0,x ? x[i] : 0)==CHEM_OK);
    return a;
}
static void zero_matrix(size_t v,int order,int constant,int fill)
{
    chem_sparse_builder b;REQUIRE(chem_sparse_builder_create(&b,3)==CHEM_OK);CHECK(b.nnz*b.blocks==0);
    chem_sparse s;REQUIRE(chem_sparse_create(&s,&b,order,v,CHEM_DOUBLE)==CHEM_OK);CHECK(s.nnz==0);
    chem_sparse_builder_destroy(&b);
    if(constant) { chem_sparse copy;REQUIRE(chem_sparse_clone(&copy,&s)==CHEM_OK);chem_sparse_destroy(&s);s=copy; }
    if(fill) { chem_array_fill(&s.values,2);for(size_t i=0;i<s.values.size;++i) NEAR(((double *)s.values.data)[i],2,0); }
    size_t k;int z=0;CHECK(chem_sparse_index(&s,0,0,0,&k)==CHEM_NOT_FOUND);
    CHECK(chem_sparse_index(&s,0,6,0,&k)==CHEM_SIZE);CHECK(chem_sparse_index(&s,0,1,3,&k)==CHEM_SIZE);CHECK(chem_sparse_index(&s,0,6,3,&k)==CHEM_SIZE);
    CHECK(chem_sparse_is_zero(&s,6,0,&z)==CHEM_SIZE);CHECK(chem_sparse_is_zero(&s,1,3,&z)==CHEM_SIZE);CHECK(chem_sparse_is_zero(&s,6,3,&z)==CHEM_SIZE);
    CHECK(chem_sparse_set(&s,0,0,4,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,1,0,0,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,0,3,0,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,0,1,1,2)==CHEM_NOT_FOUND);
    chem_sparse_destroy(&s);
}
static void single_matrix(size_t v,int order,int constant)
{
    chem_sparse s=sparse(4,1,v,order,5,basic,constant ? CHEM_INT : CHEM_DOUBLE,0);
    CHECK(s.nnz==5);CHECK(s.nnz*s.values.rows==5);
    int z;for(size_t i=0;i<5;++i) { CHECK(chem_sparse_is_zero(&s,basic[i].r,basic[i].c,&z)==CHEM_OK);CHECK(!z); }
    for(size_t i=0;i<4;i+=2) { CHECK(chem_sparse_is_zero(&s,i,i,&z)==CHEM_OK);CHECK(z); }
    CHECK(chem_sparse_is_zero(&s,3,3,&z)==CHEM_OK);CHECK(z);
    NEAR(get(&s,0,2,1),0,0);put(&s,0,2,1,45);NEAR(get(&s,0,2,1),45,0);
    put(&s,0,3,2,42);put(&s,0,2,3,21);
    if(constant) { chem_sparse copy;REQUIRE(chem_sparse_clone(&copy,&s)==CHEM_OK);chem_sparse_destroy(&s);s=copy; }
    size_t d[4],n,k;REQUIRE(chem_sparse_diagonal_indices(&s,0,d,4,&n)==CHEM_OK);CHECK(n==1);CHECK(d[0]==id(&s,1,1));
    size_t width=v ? v : 1;CHECK(id(&s,3,2)==(order ? 3 : 4)*width);CHECK(id(&s,2,3)==(order ? 4 : 3)*width);
    if(constant) { int *data=s.values.data;CHECK(data[id(&s,3,2)]==42);CHECK(data[id(&s,2,3)]==21); }
    else { double *data=s.values.data;data[id(&s,3,2)]=42;data[id(&s,2,3)]=21;NEAR(data[(order ? 3 : 4)*width],42,0);NEAR(data[(order ? 4 : 3)*width],21,0); }
    CHECK(chem_array_groups(&s.values)==1);CHECK(s.values.size==5*width);
    CHECK(chem_sparse_index(&s,0,4,2,&k)==CHEM_SIZE);CHECK(chem_sparse_index(&s,0,1,5,&k)==CHEM_SIZE);CHECK(chem_sparse_index(&s,1,0,2,&k)==CHEM_SIZE);
    CHECK(chem_sparse_index(&s,0,2,2,&k)==CHEM_NOT_FOUND);CHECK(chem_sparse_index(&s,0,0,2,&k)==CHEM_NOT_FOUND);
    CHECK(chem_sparse_set(&s,0,0,4,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,1,0,0,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,0,5,0,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,0,3,3,2)==CHEM_NOT_FOUND);
    REQUIRE(chem_sparse_index(&s,CHEM_NONE,3,2,&k)==CHEM_OK);CHECK(k==id(&s,3,2));chem_sparse_destroy(&s);
}
static void multi_matrix_basic(size_t v,int order)
{
    chem_sparse s=sparse(4,53,v,order,5,basic,CHEM_DOUBLE,24);CHECK(s.nnz==5);CHECK(s.nnz*s.values.rows==265);
    NEAR(get(&s,0,2,1),24,0);put(&s,0,2,1,45);NEAR(get(&s,0,2,1),45,0);
    put(&s,43,3,2,29);NEAR(get(&s,43,3,2),29,0);put(&s,33,2,3,79);NEAR(get(&s,33,2,3),79,0);put(&s,2,2,3,64);NEAR(get(&s,2,2,3),64,0);
    size_t d[4],n,k;for(size_t b=0;b<3;++b) { REQUIRE(chem_sparse_diagonal_indices(&s,b,d,4,&n)==CHEM_OK);CHECK(n==1);REQUIRE(chem_sparse_index(&s,b,1,1,&k)==CHEM_OK);CHECK(d[0]==k); }
    CHECK(id(&s,2,3)==(order ? 4 : 3)*(v ? v : 1));
    REQUIRE(chem_sparse_index(&s,2,2,1,&k)==CHEM_OK);CHECK(k==(v ? 14 : 12));((double *)s.values.data)[k]=31;NEAR(((double *)s.values.data)[v ? 14 : 12],31,0);
    ((double *)s.values.data)[id(&s,2,3)]=21;NEAR(get(&s,0,2,3),21,0);
    CHECK(chem_sparse_index(&s,0,4,2,&k)==CHEM_SIZE);CHECK(chem_sparse_index(&s,2,1,5,&k)==CHEM_SIZE);CHECK(chem_sparse_index(&s,54,0,2,&k)==CHEM_SIZE);
    CHECK(chem_sparse_index(&s,1,2,2,&k)==CHEM_NOT_FOUND);CHECK(chem_sparse_index(&s,2,0,2,&k)==CHEM_NOT_FOUND);CHECK(chem_sparse_index(&s,CHEM_NONE,0,1,&k)==CHEM_INVALID);
    CHECK(chem_sparse_set(&s,0,0,4,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,53,0,0,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,0,5,0,2)==CHEM_SIZE);CHECK(chem_sparse_set(&s,0,3,3,2)==CHEM_NOT_FOUND);
    chem_sparse_destroy(&s);
}
static void diagonal_add(size_t v,int order)
{
    const pair duplicate[]={{0,1},{3,2},{0,1},{2,3},{1,1},{2,1}};
    chem_sparse s=sparse(4,3,v,order,6,duplicate,CHEM_INT,24);CHECK(s.nnz==5);
    for(size_t b=0;b<3;++b) put(&s,b,1,1,(double)(1+2*b));
    chem_sparse_add_diagonal(&s,7);
    for(size_t b=0;b<3;++b) for(size_t j=0;j<5;++j) NEAR(get(&s,b,basic[j].r,basic[j].c),basic[j].r==1 ? 8+2*b : 24,0);
    chem_sparse_destroy(&s);
}
static void sparse_print(size_t v,int order,int nonzero)
{
    chem_sparse s=sparse(4,3,v,order,5,basic,CHEM_INT,nonzero ? 32 : 24);
    for(size_t b=0;b<3;++b) put(&s,b,1,1,(double)(nonzero ? 2+2*b : 1+2*b));
    chem_sparse_add_diagonal(&s,nonzero ? 5 : 7);
    FILE *f=tmpfile();REQUIRE(f!=NULL);CHECK(chem_sparse_print(f,&s,nonzero)==CHEM_OK);rewind(f);char output[1024]={0};fread(output,1,sizeof(output)-1,f);fclose(f);
    const char *expected=nonzero ?
        "Block 0\n0, 1, 32\n1, 1, 7\n2, 1, 32\n2, 3, 32\n3, 2, 32\nBlock 1\n0, 1, 32\n1, 1, 9\n2, 1, 32\n2, 3, 32\n3, 2, 32\nBlock 2\n0, 1, 32\n1, 1, 11\n2, 1, 32\n2, 3, 32\n3, 2, 32\n" :
        "Block 0\n0,24,0,0\n0,8,0,0\n0,24,0,24\n0,0,24,0\nBlock 1\n0,24,0,0\n0,10,0,0\n0,24,0,24\n0,0,24,0\nBlock 2\n0,24,0,0\n0,12,0,0\n0,24,0,24\n0,0,24,0\n";
    CHECK(!strcmp(output,expected));chem_sparse_destroy(&s);
}
static void bad_builder(void)
{
    chem_sparse_builder b;REQUIRE(chem_sparse_builder_create(&b,3)==CHEM_OK);
    CHECK(chem_sparse_builder_add(&b,3,0)==CHEM_SIZE);CHECK(chem_sparse_builder_add(&b,2,4)==CHEM_SIZE);CHECK(chem_sparse_builder_add(&b,6,7)==CHEM_SIZE);chem_sparse_builder_destroy(&b);
}
/* A short sequence of row lambdas, represented as C callback data. Slot numbers
 * are cached physical sparse indices, dense columns, or temporary numbers. */
typedef struct { size_t arg,slot; } ref;
typedef struct { ref out,in[4];size_t count;int op;double factor; } command;
typedef struct { command commands[8];size_t count; } block_program;
static chem_span resolve(ref r,const chem_array_group *g,double *temps)
{
    if(r.arg==CHEM_NONE) return chem_span_temporary(temps+r.slot*g[0].width,g[0].width,g[0].active);
    if(g[r.arg].array->kind==CHEM_SPARSE_ARRAY) return chem_sparse_block_span(g[r.arg],r.slot);
    if(g[r.arg].array->kind==CHEM_VECTOR_ARRAY) return chem_group_vector(g[r.arg]);
    return chem_group_column(g[r.arg],r.slot);
}
static double formula(const double *x,size_t n,void *opaque)
{
    command *c=opaque;
    if(c->op==1) return x[0]*x[1];
    if(c->op==2) return x[0]-x[1];
    double s=0;for(size_t i=0;i<n;++i) s+=x[i];return s*c->factor;
}
static chem_status execute(size_t n,const chem_array_group *g,void *opaque)
{
    (void)n;block_program *p=opaque;double *temps=calloc(2*g[0].width,sizeof(double));if(!temps) return CHEM_MEMORY;
    chem_status status=CHEM_OK;
    for(size_t i=0;i<p->count;++i) {
        command *c=&p->commands[i];chem_span out=resolve(c->out,g,temps),inputs[4];
        for(size_t j=0;j<c->count;++j) inputs[j]=resolve(c->in[j],g,temps);
        status=c->op==3 ? chem_span_fill(out,c->factor) : c->op==4 ? chem_span_copy(out,inputs[0]) : chem_span_map(out,c->count,inputs,formula,c);
        if(status) break;
    }
    free(temps);return status;
}
static ref element(const chem_sparse *s,size_t arg,size_t r,size_t c) { return (ref){arg,id(s,r,c)}; }
static ref temp(size_t slot) { return (ref){CHEM_NONE,slot}; }
static ref column(size_t arg,size_t col) { return (ref){arg,col}; }
static void cmd(block_program *p,ref out,int op,double factor,size_t n,const ref *in)
{ command *c=&p->commands[p->count++];*c=(command){.out=out,.count=n,.op=op,.factor=factor};for(size_t j=0;j<n;++j)c->in[j]=in[j]; }
#define SUM(p,out,factor,...) do {ref inputs_[]={__VA_ARGS__};cmd(p,out,0,factor,sizeof(inputs_)/sizeof(*inputs_),inputs_);}while(0)
#define COPY(p,out,in) do {ref input_=(in);cmd(p,out,4,1,1,&input_);}while(0)
static void start(chem_array_function *f,block_program *p,size_t n,chem_array_argument *args)
{CHECK(chem_array_function_create(f,n,args,execute,p)==CHEM_OK);CHECK(chem_array_function_apply(f,n,args)==CHEM_OK);}
static void sparse_array_function(size_t v,int order)
{
    const pair coords[]={{0,0},{1,1},{2,2},{2,3},{3,3}};chem_sparse a=sparse(4,3,v,order,5,coords,CHEM_DOUBLE,1);
    put(&a,0,0,0,1);put(&a,1,1,1,2);put(&a,2,2,2,3);put(&a,2,2,3,4);put(&a,2,3,3,5);
    block_program p={0};SUM(&p,temp(0),1,element(&a,0,0,0),element(&a,0,1,1),element(&a,0,2,2),element(&a,0,2,3));SUM(&p,element(&a,0,2,3),2,temp(0));
    chem_array_argument arg=chem_array_mutable(&a.values);chem_array_function f;start(&f,&p,1,&arg);
    for(size_t b=0;b<3;++b) { NEAR(get(&a,b,2,3),((double[]){8,10,18})[b],0);NEAR(get(&a,b,0,0),1,0);NEAR(get(&a,b,1,1),b==1 ? 2 : 1,0);NEAR(get(&a,b,2,2),b==2 ? 3 : 1,0);NEAR(get(&a,b,3,3),b==2 ? 5 : 1,0); }
    chem_sparse b=sparse(4,3,v,order,5,coords,CHEM_DOUBLE,-1);arg=chem_array_mutable(&b.values);CHECK(chem_array_function_apply(&f,1,&arg)==CHEM_OK);
    for(size_t k=0;k<3;++k) for(size_t j=0;j<5;++j) NEAR(get(&b,k,coords[j].r,coords[j].c),j==3 ? -8 : -1,0);
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);chem_sparse_destroy(&b);
}
static void sparse_multi_function(size_t v,int order,int different_structure)
{
    size_t blocks=different_structure ? 4 : 3;
    const pair ac[]={{0,1},{1,1},{2,2}},bc[]={{0,0},{1,2},{2,3},{3,4}};
    chem_sparse a=sparse(3,blocks,v,order,different_structure ? 3 : 2,different_structure ? ac : chain,CHEM_DOUBLE,0);
    chem_sparse b=sparse(different_structure ? 5 : 3,blocks,v,order,different_structure ? 4 : 3,different_structure ? bc : diagonal,CHEM_DOUBLE,0);
    for(size_t i=0;i<blocks;++i) {
        if(different_structure) { for(size_t j=0;j<3;++j) put(&a,i,ac[j].r,ac[j].c,i*(j==0 ? 2 : j==1 ? 3 : 5)+j+1);for(size_t j=0;j<4;++j) put(&b,i,bc[j].r,bc[j].c,i+10*(j+1)); }
        else { put(&a,i,0,1,i+10);put(&a,i,1,2,2*i+20);put(&b,i,0,0,4*i);put(&b,i,1,1,3*i);put(&b,i,2,2,5*i); }
    }
    block_program p={0};
    if(different_structure) { SUM(&p,temp(0),1,element(&a,0,0,1),element(&a,0,1,1),element(&b,1,0,0));COPY(&p,element(&a,0,2,2),temp(0)); }
    else {SUM(&p,temp(0),1,element(&a,0,0,1),element(&b,1,2,2));COPY(&p,element(&a,0,1,2),temp(0));}
    chem_array_argument args[]={chem_array_mutable(&a.values),chem_array_const(&b.values)};chem_array_function f;start(&f,&p,2,args);
    for(size_t i=0;i<blocks;++i) {
        NEAR(get(&a,i,different_structure ? 2 : 1,2),different_structure ? 13+6*i : 10+6*i,0);NEAR(get(&a,i,0,1),different_structure ? 1+2*i : 10+i,0);
        if(different_structure) for(size_t j=0;j<4;++j) NEAR(get(&b,i,bc[j].r,bc[j].c),i+10*(j+1),0);
        else NEAR(get(&b,i,0,0),4*i,0);
    }
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);chem_sparse_destroy(&b);
}
static void sparse_temporaries(size_t v,int order,int reuse_views)
{
    const pair coords[]={{0,1},{1,2},{2,3},{3,4},{4,4}};size_t blocks=reuse_views ? 3 : 4;
    chem_sparse a=sparse(reuse_views ? 4 : 5,blocks,v,order,reuse_views ? 4 : 5,reuse_views ? chain : coords,CHEM_DOUBLE,0);
    for(size_t b=0;b<blocks;++b) {put(&a,b,0,1,b+1);if(!reuse_views)put(&a,b,1,2,10*(b+1));}
    ref e[5];for(size_t j=0;j<(reuse_views ? 4u : 5u);++j) e[j]=element(&a,0,(reuse_views ? chain : coords)[j].r,(reuse_views ? chain : coords)[j].c);
    block_program p={0};
    if(reuse_views) {SUM(&p,e[1],2,e[0]);SUM(&p,e[2],1,e[0],e[1]);cmd(&p,e[3],1,1,2,(ref[]){e[2],e[1]});}
    else {cmd(&p,temp(0),1,1,2,(ref[]){e[0],e[1]});SUM(&p,temp(1),1,e[0],e[1]);SUM(&p,e[2],1,temp(0),temp(1));cmd(&p,e[3],2,1,2,(ref[]){temp(0),temp(1)});cmd(&p,e[4],1,1,2,(ref[]){temp(0),temp(1)});}
    chem_array_argument arg=chem_array_mutable(&a.values);chem_array_function f;start(&f,&p,1,&arg);
    for(size_t b=0;b<blocks;++b) {
        double x=b+1,prod=10*x*x,total=11*x;
        if(reuse_views) {NEAR(get(&a,b,1,2),2*x,0);NEAR(get(&a,b,2,3),3*x,0);NEAR(get(&a,b,3,3),6*x*x,0);}
        else {NEAR(get(&a,b,2,3),prod+total,0);NEAR(get(&a,b,3,4),prod-total,0);NEAR(get(&a,b,4,4),prod*total,0);}
    }
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);
}
static void sparse_reuse(size_t v,int order)
{
    chem_sparse a=sparse(3,2,v,order,3,diagonal,CHEM_DOUBLE,0);for(size_t b=0;b<2;++b)for(size_t j=0;j<3;++j)put(&a,b,j,j,b+j);
    block_program p={0};SUM(&p,temp(0),1,element(&a,0,0,0),element(&a,0,1,1),element(&a,0,2,2));SUM(&p,element(&a,0,2,2),2,temp(0));
    chem_array_function f;chem_array_argument arg=chem_array_mutable(&a.values);start(&f,&p,1,&arg);NEAR(get(&a,0,2,2),6,0);NEAR(get(&a,1,2,2),12,0);
    chem_sparse b=sparse(3,2,v,order,3,diagonal,CHEM_DOUBLE,5);arg=chem_array_mutable(&b.values);CHECK(chem_array_function_apply(&f,1,&arg)==CHEM_OK);NEAR(get(&b,0,2,2),30,0);NEAR(get(&b,1,2,2),30,0);
    chem_sparse c=sparse(3,2,v,order,3,diagonal,CHEM_DOUBLE,0);put(&c,1,0,0,10);arg=chem_array_mutable(&c.values);CHECK(chem_array_function_apply(&f,1,&arg)==CHEM_OK);NEAR(get(&c,0,2,2),0,0);NEAR(get(&c,1,2,2),20,0);
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);chem_sparse_destroy(&b);chem_sparse_destroy(&c);
}
static void sparse_dense(size_t v,int order,size_t dense_v)
{
    chem_sparse a=sparse(4,3,v,order,4,chain,CHEM_DOUBLE,0);chem_array b;CHECK(chem_array_create(&b,3,4,dense_v,CHEM_DOUBLE,0)==CHEM_OK);
    for(size_t i=0;i<3;++i) {for(size_t j=0;j<4;++j) {put(&a,i,chain[j].r,chain[j].c,i+(j ? 10*j : 1));CHECK(chem_array_set(&b,i,j,100*i+10*j)==CHEM_OK);}}
    block_program p={0};SUM(&p,temp(0),1,element(&a,0,0,1),element(&a,0,1,2),column(1,0),column(1,1));COPY(&p,element(&a,0,3,3),temp(0));
    chem_array_function f;chem_array_argument args[]={chem_array_mutable(&a.values),chem_array_const(&b)};start(&f,&p,2,args);
    for(size_t i=0;i<3;++i) {NEAR(get(&a,i,3,3),21+202*i,0);for(size_t j=0;j<3;++j) NEAR(get(&a,i,chain[j].r,chain[j].c),i+(j ? 10*j : 1),0);for(size_t j=0;j<4;++j) NEAR(chem_array_get(&b,i,j),100*i+10*j,0);}
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);chem_array_destroy(&b);
}
static chem_status idle(size_t n,const chem_array_group *g,void *ctx)
{(void)n;(void)g;if(ctx) ++*(size_t *)ctx;return CHEM_OK;}
static void sparse_const(size_t v,int order)
{
    chem_sparse a=sparse(4,3,v,order,4,chain,CHEM_DOUBLE,0);
    for(size_t i=0;i<3;++i) {put(&a,i,0,1,i+1);put(&a,i,1,2,i+10);}
    block_program p={0};SUM(&p,temp(0),1,element(&a,0,0,1),element(&a,0,1,2));COPY(&p,temp(1),element(&a,0,2,3));
    chem_array_argument arg=chem_array_const(&a.values);chem_array_function f;start(&f,&p,1,&arg);
    NEAR(get(&a,0,0,1),1,0);NEAR(get(&a,1,1,2),11,0);
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);
}
static void sparse_empty(size_t v,int order)
{
    for(int which=0;which<2;++which) {
        chem_sparse a=sparse(4,which ? 0 : 3,v,order,which ? 1 : 0,chain,CHEM_DOUBLE,0);
        chem_array_argument arg=chem_array_const(&a.values);chem_array_function f;size_t calls=0;
        CHECK(chem_array_function_create(&f,1,&arg,idle,&calls)==CHEM_OK);CHECK(chem_array_function_apply(&f,1,&arg)==CHEM_OK);
        if(which) CHECK(calls==0);
        chem_array_function_destroy(&f);chem_sparse_destroy(&a);
    }
}
static void sparse_incompatible(size_t v,int order,size_t other_v,int other_sparse)
{
    chem_sparse a=sparse(4,3,v,order,2,chain,CHEM_DOUBLE,0),b={0};chem_array dense={0};
    if(other_sparse) b=sparse(4,3,other_v,order,2,chain,CHEM_DOUBLE,0);
    else CHECK(chem_array_create(&dense,3,4,other_v,CHEM_DOUBLE,0)==CHEM_OK);
    chem_array_argument args[]={chem_array_const(&a.values),chem_array_const(other_sparse ? &b.values : &dense)};chem_array_function f;
    CHECK(chem_array_function_create(&f,2,args,idle,NULL)==CHEM_SIZE);
    chem_sparse_destroy(&a);chem_sparse_destroy(&b);chem_array_destroy(&dense);
}
static void sparse_dimensions(size_t v,int order,int mode)
{
    if(mode==1) {
        const pair coords[]={{0,1},{1,1},{2,2}};chem_sparse a=sparse(4,3,v,order,3,coords,CHEM_DOUBLE,0);size_t k;
        CHECK(chem_sparse_index(&a,0,3,3,&k)==CHEM_NOT_FOUND);chem_sparse_destroy(&a);return;
    }
    chem_sparse a=sparse(mode==3 ? 4 : 3,3,v,order,mode==3 ? 3 : 2,chain,CHEM_DOUBLE,0);
    chem_sparse b=sparse(mode==3 ? 4 : 3,4,v,order,mode==3 ? 3 : 2,chain,CHEM_DOUBLE,0);
    chem_array x=vec(3,v,CHEM_DOUBLE,NULL),y=vec(2,v,CHEM_DOUBLE,NULL);
    chem_array_argument args[4]={chem_array_mutable(&a.values),chem_array_const(&a.values),chem_array_const(&a.values),chem_array_const(&y)};
    size_t n=2,calls=0;chem_array_function f;
    if(mode==0) args[1]=chem_array_const(&b.values);
    if(mode==3) n=3;
    if(mode>=4 && mode<=8) {
        chem_array_destroy(&x);x=vec(mode==4 ? 2 : mode==5 ? 4 : mode==7 ? 3 : 0,v,CHEM_DOUBLE,NULL);args[1]=chem_array_const(&x);
        if(mode==7 || mode==8) {chem_sparse_destroy(&a);a=sparse(3,0,v,order,1,chain,CHEM_DOUBLE,0);}
    }
    if(mode==9) {n=3;args[1]=chem_array_const(&x);args[2]=chem_array_const(&y);}
    if(mode==10) {n=3;args[1]=chem_array_const(&b.values);args[2]=chem_array_const(&x);}
    if(mode==11) {n=4;args[2]=chem_array_const(&x);}
    if(mode==12) args[1]=chem_array_const(&x);
    CHECK(chem_array_function_create(&f,n,args,idle,&calls)==CHEM_OK);
    if(mode==0 || mode==2) {
        args[1]=chem_array_const(&a.values);CHECK(chem_array_function_apply(&f,n,args)==CHEM_OK);
        args[1]=chem_array_const(&b.values);CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE);
        if(mode==2) {args[0]=chem_array_mutable(&b.values);args[1]=chem_array_const(&a.values);CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE);}
    } else if(mode==3) {
        CHECK(chem_array_function_apply(&f,n,args)==CHEM_OK);args[2]=chem_array_const(&b.values);CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE);
        args[2]=chem_array_const(&a.values);args[0]=chem_array_mutable(&b.values);CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE);
    } else if(mode==8) {CHECK(chem_array_function_apply(&f,n,args)==CHEM_OK);CHECK(calls==0);}
    else {
        if(mode==12) args[1]=chem_array_const(&y);
        CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE);CHECK(calls==0);
        if(mode==7) {size_t k;CHECK(chem_sparse_index(&a,0,0,1,&k)==CHEM_SIZE);}
    }
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);chem_sparse_destroy(&b);chem_array_destroy(&x);chem_array_destroy(&y);
}
static void sparse_flexible(size_t v,int order,int with_vector)
{
    chem_sparse a=sparse(4,3,v,order,with_vector ? 2 : 3,chain,CHEM_DOUBLE,0),b=sparse(4,3,v,order,3,chain,CHEM_DOUBLE,0);
    chem_array x=vec(3,v,CHEM_DOUBLE,(double[]){1,2,3});block_program p={0};
    if(with_vector) {SUM(&p,temp(0),1,element(&a,0,0,1),column(1,0));COPY(&p,element(&a,0,1,2),temp(0));}
    else {SUM(&p,temp(0),1,element(&a,0,0,1),element(&a,0,1,2),element(&b,1,2,3));COPY(&p,element(&a,0,2,3),temp(0));}
    chem_array_argument args[]={chem_array_mutable(&a.values),chem_array_const(with_vector ? &x : &b.values)};chem_array_function f;
    CHECK(chem_array_function_create(&f,2,args,execute,&p)==CHEM_OK);
    chem_sparse_destroy(&a);a=sparse(4,with_vector ? 5 : 4,v,order,with_vector ? 2 : 3,chain,CHEM_DOUBLE,0);
    if(with_vector) {chem_array_destroy(&x);x=vec(5,v,CHEM_DOUBLE,(double[]){10,20,30,40,50});}
    else {chem_sparse_destroy(&b);b=sparse(4,4,v,order,3,chain,CHEM_DOUBLE,0);}
    for(size_t i=0;i<a.values.rows;++i) {put(&a,i,0,1,i+1);if(!with_vector) {put(&a,i,1,2,i+10);put(&b,i,2,3,i+100);}}
    CHECK(chem_array_function_apply(&f,2,args)==CHEM_OK);
    for(size_t i=0;i<a.values.rows;++i) NEAR(get(&a,i,with_vector ? 1 : 2,with_vector ? 2 : 3),with_vector ? 11*(i+1) : 111+3*i,0);
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);chem_sparse_destroy(&b);chem_array_destroy(&x);
}
static void sparse_vectors(size_t v,int order,int mode)
{
    /* Modes: input, two inputs, two matrices, const input, mutable output,
     * reusable input, mixed temporary, integer input, const callback signature. */
    size_t count=mode==0 || mode==6 ? 3 : 1;
    const pair first=mode==8 ? (pair){0,0} : chain[0];
    chem_sparse a=sparse(mode==8 ? 2 : 4,3,v,order,count,count==1 ? &first : chain,CHEM_DOUBLE,0);
    chem_sparse b=sparse(4,3,v,order,1,&first,CHEM_DOUBLE,0);
    chem_array x=vec(3,v,mode==7 ? CHEM_INT : CHEM_DOUBLE,(double[]){5,10,15}),y=vec(3,v,CHEM_DOUBLE,(double[]){4,5,6});
    for(size_t i=0;i<3;++i) for(size_t j=0;j<count;++j) put(&a,i,count==1 ? first.r : chain[j].r,count==1 ? first.c : chain[j].c,i+j+1);
    if(mode==1 || mode==5 || mode==8) for(size_t i=0;i<3;++i) CHECK(chem_array_set(&x,i,0,i+1)==CHEM_OK);
    if(mode==2) for(size_t i=0;i<3;++i) CHECK(chem_array_set(&x,i,0,10*(i+1))==CHEM_OK);
    if(mode==6) for(size_t i=0;i<3;++i) CHECK(chem_array_set(&x,i,0,100*(i+1))==CHEM_OK);
    block_program p={0};ref out=element(&a,0,first.r,first.c);size_t nargs=mode==1 || mode==2 ? 3 : 2;
    chem_array_argument args[]={mode==4 || mode==2 ? chem_array_const(&a.values) : chem_array_mutable(&a.values),mode==4 ? chem_array_mutable(&x) : chem_array_const(&x),chem_array_const(&y)};
    if(mode==1) {SUM(&p,out,1,column(1,0),column(2,0));}
    else if(mode==2) {args[1]=chem_array_mutable(&b.values);args[2]=chem_array_const(&x);SUM(&p,element(&b,1,0,1),1,out,column(2,0));}
    else if(mode==4) {SUM(&p,column(1,0),3,out);}
    else if(mode==6) {SUM(&p,temp(0),1,element(&a,0,0,1),element(&a,0,1,2),column(1,0));SUM(&p,element(&a,0,2,3),2,temp(0));}
    else {SUM(&p,out,mode==3 ? 3 : mode==7 ? 1.5 : 2,column(1,0));}
    chem_array_function f;start(&f,&p,nargs,args);
    for(size_t i=0;i<3;++i) {
        double expected=mode==0 ? 10*(i+1) : mode==1 ? 5+2*i : mode==2 ? 11*(i+1) : mode==3 ? 15*(i+1) : mode==4 ? 3*(i+1) : mode==5 || mode==8 ? 2*(i+1) : mode==6 ? 206+204*i : 7.5*(i+1);
        double actual=mode==4 ? chem_array_get(&x,i,0) : mode==2 ? get(&b,i,0,1) : get(&a,i,mode==6 ? 2 : first.r,mode==6 ? 3 : first.c);
        NEAR(actual,expected,0);
    }
    if(mode==5) {
        args[0]=chem_array_mutable(&b.values);for(size_t i=0;i<3;++i) CHECK(chem_array_set(&y,i,0,10*(i+1))==CHEM_OK);args[1]=chem_array_const(&y);
        CHECK(chem_array_function_apply(&f,nargs,args)==CHEM_OK);for(size_t i=0;i<3;++i) NEAR(get(&b,i,0,1),20*(i+1),0);
    }
    chem_array_function_destroy(&f);chem_sparse_destroy(&a);chem_sparse_destroy(&b);chem_array_destroy(&x);chem_array_destroy(&y);
}
static void sparse_indices(size_t v,int order)
{
    const pair coords[]={{0,0},{0,2},{1,0},{1,1},{2,1},{2,2}};
    chem_sparse a=sparse(3,3,v,order,6,coords,CHEM_DOUBLE,0);chem_array x=vec(3,v,CHEM_INDEX,(double[]){0,1,2});block_program p={0};
    /* Cached index capture, independently of CSR versus CSC ordinal order. */
    for(size_t j=0;j<6;++j) {
        cmd(&p,temp(0),3,10*j,0,NULL);
        SUM(&p,element(&a,0,coords[j].r,coords[j].c),1,column(1,0),temp(0));
        /* At most eight commands fit; apply each lambda to all groups. */
        chem_array_function f;chem_array_argument args[]={chem_array_mutable(&a.values),chem_array_const(&x)};start(&f,&p,2,args);chem_array_function_destroy(&f);p.count=0;
    }
    for(size_t b=0;b<3;++b) for(size_t j=0;j<6;++j) NEAR(get(&a,b,coords[j].r,coords[j].c),b+10*j,0);
    chem_sparse_destroy(&a);chem_array_destroy(&x);
}
static void sparse_fill(size_t v,int order)
{
    const pair coords[]={{0,1},{2,0}};chem_sparse a=sparse(3,3,v,order,2,coords,CHEM_DOUBLE,0);chem_array x=vec(3,v,CHEM_DOUBLE,NULL);
    for(int mode=0;mode<3;++mode) {
        block_program p={0};ref out=mode==0 ? element(&a,0,0,1) : mode==1 ? column(1,0) : temp(0);
        cmd(&p,out,3,mode==2 ? 9.9 : 3.2,0,NULL);if(mode==2) COPY(&p,element(&a,0,2,0),temp(0));
        chem_array_argument args[]={mode==1 ? chem_array_const(&a.values) : chem_array_mutable(&a.values),chem_array_mutable(&x)};
        chem_array_function f;start(&f,&p,mode==1 ? 2 : 1,args);
        for(size_t b=0;b<3;++b) {
            if(mode==0) {NEAR(get(&a,b,0,1),3.2,0);NEAR(get(&a,b,2,0),0,0);}
            if(mode==1) NEAR(chem_array_get(&x,b,0),3.2,0);
            if(mode==2) NEAR(get(&a,b,2,0),9.9,0);
        }
        chem_array_function_destroy(&f);
    }
    chem_sparse_destroy(&a);chem_array_destroy(&x);
}
static void sparse_copy(size_t v,int order)
{
    const pair coords[]={{0,1},{2,0}};
    for(int mode=0;mode<8;++mode) {
        chem_sparse a=sparse(3,3,v,order,2,coords,CHEM_DOUBLE,0),b=sparse(3,3,v,order,1,(pair[]){{1,2}},CHEM_DOUBLE,0);
        const double expected[]={3.2,4.2,1.3};for(size_t i=0;i<3;++i) put(&a,i,0,1,mode==0 ? 0 : expected[i]);
        chem_array x=vec(3,v,CHEM_DOUBLE,mode==0 ? (double[]){10,20,30} : NULL),dense;
        CHECK(chem_array_create(&dense,3,2,v,CHEM_DOUBLE,0)==CHEM_OK);for(size_t i=0;i<3;++i) CHECK(chem_array_set(&dense,i,1,100*(i+1))==CHEM_OK);
        block_program p={0};chem_array_argument args[]={chem_array_mutable(&a.values),chem_array_const(&x)};size_t n=1;
        if(mode==0) {n=2;COPY(&p,element(&a,0,0,1),column(1,0));}
        if(mode==1 || mode==6) {n=2;args[0]=chem_array_const(&a.values);args[1]=chem_array_mutable(&x);COPY(&p,column(1,0),element(&a,0,0,1));}
        if(mode==2 || mode==3) {COPY(&p,element(&a,0,2,0),element(&a,0,0,1));}
        if(mode==4) {COPY(&p,temp(0),element(&a,0,0,1));cmd(&p,element(&a,0,0,1),3,0,0,NULL);COPY(&p,element(&a,0,0,1),temp(0));}
        if(mode==5) {n=2;args[1]=chem_array_const(&dense);COPY(&p,element(&a,0,2,0),column(1,1));}
        if(mode==7) {n=2;args[0]=chem_array_mutable(&b.values);args[1]=chem_array_const(&a.values);COPY(&p,element(&b,0,1,2),element(&a,1,0,1));}
        chem_array_function f;start(&f,&p,n,args);chem_array_function_destroy(&f);
        if(mode==6) {
            p.count=0;args[0]=chem_array_mutable(&dense);args[1]=chem_array_const(&x);COPY(&p,column(0,0),column(1,0));start(&f,&p,2,args);chem_array_function_destroy(&f);
        }
        for(size_t i=0;i<3;++i) {
            if(mode==0) {NEAR(get(&a,i,0,1),10*(i+1),0);NEAR(get(&a,i,2,0),0,0);}
            if(mode==1) NEAR(chem_array_get(&x,i,0),expected[i],0);
            if(mode==2 || mode==3) {NEAR(get(&a,i,2,0),expected[i],0);NEAR(get(&a,i,0,1),expected[i],0);}
            if(mode==4) NEAR(get(&a,i,0,1),expected[i],0);
            if(mode==5) NEAR(get(&a,i,2,0),100*(i+1),0);
            if(mode==6) NEAR(chem_array_get(&dense,i,0),expected[i],0);
            if(mode==7) NEAR(get(&b,i,1,2),expected[i],0);
        }
        chem_sparse_destroy(&a);chem_sparse_destroy(&b);chem_array_destroy(&x);chem_array_destroy(&dense);
    }
}
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,ZeroMatrix) { zero_matrix(0,0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,ConstZeroMatrix) { zero_matrix(0,0,1,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,SetScalar) { zero_matrix(0,0,0,1); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,AddToDiagonal) { diagonal_add(0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,SingleBlockMatrix) { single_matrix(0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,ConstSingleBlockMatrix) { single_matrix(0,0,1); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultiBlockMatrix) { multi_matrix_basic(0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,Print) { sparse_print(0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,PrintNonZero) { sparse_print(0,0,1); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseMatrixBuilder,BadConfiguration) { bad_builder(); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,ArrayFunction) { sparse_array_function(0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultiMatrixArrayFunction) { sparse_multi_function(0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MismatchedBlockDimensions) { sparse_dimensions(0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MismatchedElementDimensions) { sparse_dimensions(0,0,1); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultipleTemporaries) { sparse_temporaries(0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,BlockViewReuse) { sparse_temporaries(0,0,1); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,FunctionReusability) { sparse_reuse(0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,TwoSparseMatricesDifferentStructure) { sparse_multi_function(0,0,1); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,SparseAndDenseMatrixFunction) { sparse_dense(0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,ConstSparseMatrixFunction) { sparse_const(0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,EmptySparseMatrixFunction) { sparse_empty(0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultipleSparseMatricesDifferentBlocksFromCreation) { sparse_flexible(0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,SparseMatrixVectorDifferentBlocksFromCreation) { sparse_flexible(0,0,1); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MismatchedBlocksAtInvocation) { sparse_dimensions(0,0,2); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultipleSparseMatricesMismatchedBlocksAtInvocation) { sparse_dimensions(0,0,3); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,VectorInSparseMatrixFunction) { sparse_vectors(0,0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,VectorTooSmall) { sparse_dimensions(0,0,4); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,VectorTooLarge) { sparse_dimensions(0,0,5); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,EmptyVectorNonEmptySparseMatrix) { sparse_dimensions(0,0,6); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,NonEmptyVectorEmptySparseMatrix) { sparse_dimensions(0,0,7); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,EmptyVectorEmptySparseMatrix) { sparse_dimensions(0,0,8); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultipleVectorsDifferentSizes) { sparse_dimensions(0,0,9); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultipleVectorsSameSize) { sparse_vectors(0,0,1); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultipleSparseMatricesOneVector) { sparse_vectors(0,0,2); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MultipleSparseMatricesDifferentBlocksVector) { sparse_dimensions(0,0,10); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,VectorSizeMatchesOneSparseMatrixOnly) { sparse_dimensions(0,0,11); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,ConstVectorSparse) { sparse_vectors(0,0,3); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MutableVectorSparse) { sparse_vectors(0,0,4); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,FunctionWithConstSignatureSparse) { sparse_vectors(0,0,8); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,GetBlockViewByVectorIndex) { sparse_indices(0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,FunctionReusabilityWithVectorsSparse) { sparse_vectors(0,0,5); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,FunctionInvocationWithWrongSizedVectorSparse) { sparse_dimensions(0,0,12); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,MixedVectorBlockViewBlockVariable) { sparse_vectors(0,0,6); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,IntegerVectorSparse) { sparse_vectors(0,0,7); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,TestFill) { sparse_fill(0,0); }
TEST(sparse_matrix_standard_compressed_row_ordering,SparseCompressedRowMatrix,TestCopy) { sparse_copy(0,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,ZeroMatrix) { zero_matrix(0,1,0,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,ConstZeroMatrix) { zero_matrix(0,1,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,SetScalar) { zero_matrix(0,1,0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,AddToDiagonal) { diagonal_add(0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,SingleBlockMatrix) { single_matrix(0,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,ConstSingleBlockMatrix) { single_matrix(0,1,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultiBlockMatrix) { multi_matrix_basic(0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,Print) { sparse_print(0,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,PrintNonZero) { sparse_print(0,1,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseMatrixBuilder,BadConfiguration) { bad_builder(); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,ArrayFunction) { sparse_array_function(0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultiMatrixArrayFunction) { sparse_multi_function(0,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MismatchedBlockDimensions) { sparse_dimensions(0,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MismatchedElementDimensions) { sparse_dimensions(0,1,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultipleTemporaries) { sparse_temporaries(0,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,BlockViewReuse) { sparse_temporaries(0,1,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,FunctionReusability) { sparse_reuse(0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,TwoSparseMatricesDifferentStructure) { sparse_multi_function(0,1,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,SparseAndDenseMatrixFunction) { sparse_dense(0,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,ConstSparseMatrixFunction) { sparse_const(0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,EmptySparseMatrixFunction) { sparse_empty(0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultipleSparseMatricesDifferentBlocksFromCreation) { sparse_flexible(0,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,SparseMatrixVectorDifferentBlocksFromCreation) { sparse_flexible(0,1,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MismatchedBlocksAtInvocation) { sparse_dimensions(0,1,2); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultipleSparseMatricesMismatchedBlocksAtInvocation) { sparse_dimensions(0,1,3); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,VectorInSparseMatrixFunction) { sparse_vectors(0,1,0); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,VectorTooSmall) { sparse_dimensions(0,1,4); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,VectorTooLarge) { sparse_dimensions(0,1,5); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,EmptyVectorNonEmptySparseMatrix) { sparse_dimensions(0,1,6); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,NonEmptyVectorEmptySparseMatrix) { sparse_dimensions(0,1,7); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,EmptyVectorEmptySparseMatrix) { sparse_dimensions(0,1,8); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultipleVectorsDifferentSizes) { sparse_dimensions(0,1,9); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultipleVectorsSameSize) { sparse_vectors(0,1,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultipleSparseMatricesOneVector) { sparse_vectors(0,1,2); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MultipleSparseMatricesDifferentBlocksVector) { sparse_dimensions(0,1,10); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,VectorSizeMatchesOneSparseMatrixOnly) { sparse_dimensions(0,1,11); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,ConstVectorSparse) { sparse_vectors(0,1,3); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MutableVectorSparse) { sparse_vectors(0,1,4); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,FunctionReusabilityWithVectorsSparse) { sparse_vectors(0,1,5); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,FunctionInvocationWithWrongSizedVectorSparse) { sparse_dimensions(0,1,12); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,MixedVectorBlockViewBlockVariable) { sparse_vectors(0,1,6); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,IntegerVectorSparse) { sparse_vectors(0,1,7); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,FunctionWithConstSignatureSparse) { sparse_vectors(0,1,8); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,GetBlockViewByVectorIndex) { sparse_indices(0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,TestFill) { sparse_fill(0,1); }
TEST(sparse_matrix_standard_compressed_column_ordering,SparseCompressedColumnMatrix,TestCopy) { sparse_copy(0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,ZeroMatrix) { zero_matrix(2,0,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,ConstZeroMatrix) { zero_matrix(3,0,1,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,SetScalar) { zero_matrix(3,0,0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,AddToDiagonal) { diagonal_add(1,0); diagonal_add(3,0); diagonal_add(4,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,SingleBlockMatrix) { single_matrix(4,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,ConstSingleBlockMatrix) { single_matrix(2,0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultiBlockMatrix) { multi_matrix_basic(2,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,Print) { sparse_print(1,0,0); sparse_print(2,0,0); sparse_print(3,0,0); sparse_print(4,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,PrintNonZero) { sparse_print(1,0,1); sparse_print(2,0,1); sparse_print(3,0,1); sparse_print(4,0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,ArrayFunction) { sparse_array_function(1,0); sparse_array_function(2,0); sparse_array_function(3,0); sparse_array_function(4,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultiMatrixArrayFunction) { sparse_multi_function(1,0,0); sparse_multi_function(2,0,0); sparse_multi_function(3,0,0); sparse_multi_function(4,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultipleTemporaries) { sparse_temporaries(1,0,0); sparse_temporaries(2,0,0); sparse_temporaries(3,0,0); sparse_temporaries(4,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,BlockViewReuse) { sparse_temporaries(1,0,1); sparse_temporaries(2,0,1); sparse_temporaries(3,0,1); sparse_temporaries(4,0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,FunctionReusability) { sparse_reuse(1,0); sparse_reuse(2,0); sparse_reuse(3,0); sparse_reuse(4,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,TwoSparseMatricesDifferentStructure) { sparse_multi_function(1,0,1); sparse_multi_function(2,0,1); sparse_multi_function(3,0,1); sparse_multi_function(4,0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,SparseAndDenseMatrixFunction) { sparse_dense(1,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,SparseAndVectorMatrixFunction) { sparse_dense(1,0,1); sparse_dense(2,0,2); sparse_dense(3,0,3); sparse_dense(4,0,4); sparse_dense(10,0,10); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,IncompatibleOrdering) { sparse_incompatible(2,0,0,0); sparse_incompatible(3,0,0,0); sparse_incompatible(4,0,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,IncompatibleVectorOrdering) { sparse_incompatible(2,0,1,0); sparse_incompatible(2,0,3,0); sparse_incompatible(4,0,2,0); sparse_incompatible(4,0,10,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,IncompatibleSparseOrdering) { sparse_incompatible(1,0,2,1); sparse_incompatible(2,0,3,1); sparse_incompatible(3,0,4,1); sparse_incompatible(4,0,10,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MismatchedBlockDimensions) { sparse_dimensions(1,0,0); sparse_dimensions(2,0,0); sparse_dimensions(3,0,0); sparse_dimensions(4,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MismatchedElementDimensions) { sparse_dimensions(1,0,1); sparse_dimensions(2,0,1); sparse_dimensions(3,0,1); sparse_dimensions(4,0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,ConstSparseMatrixFunction) { sparse_const(1,0); sparse_const(2,0); sparse_const(3,0); sparse_const(4,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,EmptySparseMatrixFunction) { sparse_empty(1,0); sparse_empty(2,0); sparse_empty(3,0); sparse_empty(4,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultipleSparseMatricesDifferentBlocksFromCreation) { sparse_flexible(1,0,0); sparse_flexible(2,0,0); sparse_flexible(3,0,0); sparse_flexible(4,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,SparseMatrixVectorDifferentBlocksFromCreation) { sparse_flexible(1,0,1); sparse_flexible(2,0,1); sparse_flexible(3,0,1); sparse_flexible(4,0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MismatchedBlocksAtInvocation) { sparse_dimensions(1,0,2); sparse_dimensions(2,0,2); sparse_dimensions(3,0,2); sparse_dimensions(4,0,2); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultipleSparseMatricesMismatchedBlocksAtInvocation) { sparse_dimensions(1,0,3); sparse_dimensions(2,0,3); sparse_dimensions(3,0,3); sparse_dimensions(4,0,3); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,VectorInSparseMatrixFunction) { sparse_vectors(1,0,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,VectorTooSmall) { sparse_dimensions(1,0,4); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,VectorTooLarge) { sparse_dimensions(1,0,5); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,EmptyVectorNonEmptySparseMatrix) { sparse_dimensions(1,0,6); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,NonEmptyVectorEmptySparseMatrix) { sparse_dimensions(1,0,7); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,EmptyVectorEmptySparseMatrix) { sparse_dimensions(1,0,8); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultipleVectorsDifferentSizes) { sparse_dimensions(1,0,9); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultipleVectorsSameSize) { sparse_vectors(1,0,1); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultipleSparseMatricesOneVector) { sparse_vectors(1,0,2); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MultipleSparseMatricesDifferentBlocksVector) { sparse_dimensions(1,0,10); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,VectorSizeMatchesOneSparseMatrixOnly) { sparse_dimensions(1,0,11); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,ConstVectorSparse) { sparse_vectors(1,0,3); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MutableVectorSparse) { sparse_vectors(1,0,4); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,FunctionReusabilityWithVectorsSparse) { sparse_vectors(1,0,5); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,FunctionInvocationWithWrongSizedVectorSparse) { sparse_dimensions(1,0,12); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,MixedVectorBlockViewBlockVariable) { sparse_vectors(1,0,6); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,IntegerVectorSparse) { sparse_vectors(1,0,7); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,FunctionWithConstSignatureSparse) { sparse_vectors(1,0,8); sparse_vectors(2,0,8); sparse_vectors(3,0,8); sparse_vectors(4,0,8); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,GetBlockViewByVectorIndex) { sparse_indices(1,0); sparse_indices(2,0); sparse_indices(3,0); sparse_indices(4,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,TestFill) { sparse_fill(1,0); sparse_fill(2,0); sparse_fill(3,0); sparse_fill(4,0); }
TEST(sparse_matrix_vector_compressed_row_ordering,SparseVectorCompressedRowMatrix,TestCopy) { sparse_copy(1,0); sparse_copy(2,0); sparse_copy(3,0); sparse_copy(4,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,ZeroMatrix) { zero_matrix(2,1,0,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,ConstZeroMatrix) { zero_matrix(3,1,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,SetScalar) { zero_matrix(3,1,0,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,AddToDiagonal) { diagonal_add(1,1); diagonal_add(3,1); diagonal_add(4,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,SingleBlockMatrix) { single_matrix(4,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,ConstSingleBlockMatrix) { single_matrix(2,1,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultiBlockMatrix) { multi_matrix_basic(2,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,Print) { sparse_print(1,1,0); sparse_print(2,1,0); sparse_print(3,1,0); sparse_print(4,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,PrintNonZero) { sparse_print(1,1,1); sparse_print(2,1,1); sparse_print(3,1,1); sparse_print(4,1,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,ArrayFunction) { sparse_array_function(1,1); sparse_array_function(2,1); sparse_array_function(3,1); sparse_array_function(4,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultiMatrixArrayFunction) { sparse_multi_function(1,1,0); sparse_multi_function(2,1,0); sparse_multi_function(3,1,0); sparse_multi_function(4,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultipleTemporaries) { sparse_temporaries(1,1,0); sparse_temporaries(2,1,0); sparse_temporaries(3,1,0); sparse_temporaries(4,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,BlockViewReuse) { sparse_temporaries(1,1,1); sparse_temporaries(2,1,1); sparse_temporaries(3,1,1); sparse_temporaries(4,1,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,FunctionReusability) { sparse_reuse(1,1); sparse_reuse(2,1); sparse_reuse(3,1); sparse_reuse(4,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,TwoSparseMatricesDifferentStructure) { sparse_multi_function(1,1,1); sparse_multi_function(2,1,1); sparse_multi_function(3,1,1); sparse_multi_function(4,1,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,SparseAndDenseMatrixFunction) { sparse_dense(1,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,SparseAndVectorMatrixFunction) { sparse_dense(1,1,1); sparse_dense(2,1,2); sparse_dense(3,1,3); sparse_dense(4,1,4); sparse_dense(10,1,10); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,IncompatibleOrdering) { sparse_incompatible(2,1,0,0); sparse_incompatible(3,1,0,0); sparse_incompatible(4,1,0,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,IncompatibleVectorOrdering) { sparse_incompatible(2,1,1,0); sparse_incompatible(2,1,3,0); sparse_incompatible(4,1,2,0); sparse_incompatible(4,1,10,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,IncompatibleSparseOrdering) { sparse_incompatible(1,1,2,1); sparse_incompatible(2,1,3,1); sparse_incompatible(3,1,4,1); sparse_incompatible(4,1,10,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MismatchedBlockDimensions) { sparse_dimensions(1,1,0); sparse_dimensions(2,1,0); sparse_dimensions(3,1,0); sparse_dimensions(4,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MismatchedElementDimensions) { sparse_dimensions(1,1,1); sparse_dimensions(2,1,1); sparse_dimensions(3,1,1); sparse_dimensions(4,1,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,ConstSparseMatrixFunction) { sparse_const(1,1); sparse_const(2,1); sparse_const(3,1); sparse_const(4,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,EmptySparseMatrixFunction) { sparse_empty(1,1); sparse_empty(2,1); sparse_empty(3,1); sparse_empty(4,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultipleSparseMatricesDifferentBlocksFromCreation) { sparse_flexible(1,1,0); sparse_flexible(2,1,0); sparse_flexible(3,1,0); sparse_flexible(4,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,SparseMatrixVectorDifferentBlocksFromCreation) { sparse_flexible(1,1,1); sparse_flexible(2,1,1); sparse_flexible(3,1,1); sparse_flexible(4,1,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MismatchedBlocksAtInvocation) { sparse_dimensions(1,1,2); sparse_dimensions(2,1,2); sparse_dimensions(3,1,2); sparse_dimensions(4,1,2); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultipleSparseMatricesMismatchedBlocksAtInvocation) { sparse_dimensions(1,1,3); sparse_dimensions(2,1,3); sparse_dimensions(3,1,3); sparse_dimensions(4,1,3); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,VectorInSparseMatrixFunction) { sparse_vectors(1,1,0); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,VectorTooSmall) { sparse_dimensions(1,1,4); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,VectorTooLarge) { sparse_dimensions(1,1,5); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,EmptyVectorNonEmptySparseMatrix) { sparse_dimensions(1,1,6); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,NonEmptyVectorEmptySparseMatrix) { sparse_dimensions(1,1,7); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,EmptyVectorEmptySparseMatrix) { sparse_dimensions(1,1,8); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultipleVectorsDifferentSizes) { sparse_dimensions(1,1,9); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultipleVectorsSameSize) { sparse_vectors(1,1,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultipleSparseMatricesOneVector) { sparse_vectors(1,1,2); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MultipleSparseMatricesDifferentBlocksVector) { sparse_dimensions(1,1,10); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,VectorSizeMatchesOneSparseMatrixOnly) { sparse_dimensions(1,1,11); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,ConstVectorSparse) { sparse_vectors(1,1,3); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MutableVectorSparse) { sparse_vectors(1,1,4); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,FunctionReusabilityWithVectorsSparse) { sparse_vectors(1,1,5); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,FunctionInvocationWithWrongSizedVectorSparse) { sparse_dimensions(1,1,12); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,MixedVectorBlockViewBlockVariable) { sparse_vectors(1,1,6); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,IntegerVectorSparse) { sparse_vectors(1,1,7); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,FunctionWithConstSignatureSparse) { sparse_vectors(1,1,8); sparse_vectors(2,1,8); sparse_vectors(3,1,8); sparse_vectors(4,1,8); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,GetBlockViewByVectorIndex) { sparse_indices(1,1); sparse_indices(2,1); sparse_indices(3,1); sparse_indices(4,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,TestFill) { sparse_fill(1,1); sparse_fill(2,1); sparse_fill(3,1); sparse_fill(4,1); }
TEST(sparse_matrix_vector_compressed_column_ordering,SparseVectorCompressedColumnMatrix,TestCopy) { sparse_copy(1,1); sparse_copy(2,1); sparse_copy(3,1); sparse_copy(4,1); }
