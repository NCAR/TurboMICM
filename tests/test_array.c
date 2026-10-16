#include "test.h"
#include "chem_array.h"
#include <float.h>

/* C ports of test_matrix_policy.hpp. Every VectorMatrix wrapper exercises all
 * four original lane widths unless the reference checks a particular layout. */
static chem_array array(size_t r,size_t c,size_t v,double x)
{ chem_array a; CHECK(chem_array_create(&a,r,c,v,CHEM_DOUBLE,x)==CHEM_OK); return a; }
static chem_array vector(size_t n,size_t v,const double *values)
{
    chem_array a; CHECK(chem_vector_create(&a,n,v,CHEM_DOUBLE,0)==CHEM_OK);
    if (values) for (size_t i=0;i<n;++i) CHECK(chem_array_set(&a,i,0,values[i])==CHEM_OK);
    return a;
}
static void values(const chem_array *a,const double *x)
{ for (size_t i=0;i<a->rows;++i) for (size_t j=0;j<a->columns;++j) NEAR(chem_array_get(a,i,j),x[i*a->columns+j],0); }
static void set(chem_array *a,size_t i,size_t j,double x) { CHECK(chem_array_set(a,i,j,x)==CHEM_OK); }
static void dispose(chem_array *a,chem_array *b,chem_array *c)
{ chem_array_destroy(a); if (b) chem_array_destroy(b); if (c) chem_array_destroy(c); }
static void small(size_t v,int constant)
{
    chem_array a=array(3,5,v,0), copy={0};
    set(&a,1,3,64.7); set(&a,0,0,41.2); set(&a,2,4,102.3);
    if (constant) { REQUIRE(chem_array_clone(&copy,&a)==CHEM_OK); chem_array_destroy(&a); a=copy; }
    const chem_array *p=&a;
    NEAR(chem_array_get(p,1,3),64.7,0); NEAR(chem_array_get(p,0,0),41.2,0); NEAR(chem_array_get(p,2,4),102.3,0);
    CHECK(a.size==(v ? 20 : 15));
    const double *data=a.data; NEAR(data[0],41.2,0);
    NEAR(data[v ? 18 : 14],102.3,0); NEAR(data[v ? 1+v*3 : 8],64.7,0);
    if (v) { CHECK(chem_array_groups(&a)==(v==2 ? 2u : 1u)); CHECK(a.vector_length*a.columns==v*5); }
    chem_array_destroy(&a);
}
static void initialize(size_t v)
{ chem_array a=array(2,3,v,12.4); NEAR(chem_array_get(&a,0,0),12.4,0); NEAR(chem_array_get(&a,1,0),12.4,0); NEAR(chem_array_get(&a,1,2),12.4,0); chem_array_destroy(&a); }
static void loop(size_t v,int clone)
{
    chem_array a; REQUIRE(chem_array_create(&a,3,4,v,CHEM_INT,0)==CHEM_OK);
    for (size_t i=0;i<3;++i) for (size_t j=0;j<4;++j) set(&a,i,j,(double)(100*i+j));
    chem_array copy={0}; if (clone) { REQUIRE(chem_array_clone(&copy,&a)==CHEM_OK); chem_array_destroy(&a); a=copy; }
    NEAR(chem_array_get(&a,0,0),0,0); NEAR(chem_array_get(&a,1,2),102,0); NEAR(chem_array_get(&a,2,3),203,0); NEAR(chem_array_get(&a,0,3),3,0);
    chem_array_destroy(&a);
}
static void strides(size_t v)
{
    chem_array a; REQUIRE(chem_array_create(&a,3,4,v,CHEM_INT,0)==CHEM_OK);
    CHECK(chem_array_row_stride(&a)==(v ? 1 : 4)); CHECK(chem_array_column_stride(&a)==(v ? v : 1));
    int *p=a.data;
    for (size_t i=0;i<3;++i) for (size_t j=0;j<4;++j) p[i*chem_array_row_stride(&a)+j*chem_array_column_stride(&a)]=(int)(100*i+j);
    for (size_t i=0;i<3;++i) for (size_t j=0;j<4;++j) NEAR(chem_array_get(&a,i,j),(double)(100*i+j),0);
    chem_array_destroy(&a);
}
static void to_vector(size_t v,int constant)
{
    chem_array a=array(2,3,v,0),copy={0}; const double x[]={13.2,31.2,314.2}; double y[3];
    REQUIRE(chem_array_assign_row(&a,1,x,3)==CHEM_OK);
    if (constant) { REQUIRE(chem_array_clone(&copy,&a)==CHEM_OK); chem_array_destroy(&a); a=copy; }
    REQUIRE(chem_array_export_row(&a,1,y,3)==CHEM_OK); for(size_t i=0;i<3;++i) NEAR(y[i],x[i],0);
    chem_array_destroy(&a);
}
static void from_vector(size_t v)
{
    chem_array zero={0},a={0}; REQUIRE(chem_array_import(&zero,0,NULL,NULL,v)==CHEM_OK); CHECK(zero.rows==0); chem_array_destroy(&zero);
    const double x[]={412.3,32.4,41.3,5.33,-.3,31.2}; const double *rows[]={x,x+3}; size_t lengths[]={3,3};
    REQUIRE(chem_array_import(&a,2,lengths,rows,v)==CHEM_OK); CHECK(a.rows==2); CHECK(a.columns==3); values(&a,x); chem_array_destroy(&a);
    const double a0[]={3},a1[]={4,5},a2[]={5}; const double *bad[]={a0,a1,a2}; size_t bad_lengths[]={1,2,1};
    CHECK(chem_array_import(&a,3,bad_lengths,bad,v)==CHEM_SIZE);
}
static void assign_vector(size_t v)
{
    chem_array a=array(4,3,v,0); const double x[]={12.3,15.1,24.3},big[]={14.3,52.3,65.7,16.34},little[]={13.2,52.8};
    REQUIRE(chem_array_assign_row(&a,2,x,3)==CHEM_OK);
    for(size_t j=0;j<3;++j) NEAR(chem_array_get(&a,2,j),x[j],0);
    NEAR(chem_array_get(&a,0,0),0,0); NEAR(chem_array_get(&a,3,0),0,0);
    REQUIRE(chem_array_assign_row(&a,2,big,4)==CHEM_OK);
    for(size_t j=0;j<3;++j) NEAR(chem_array_get(&a,2,j),big[j],0);
    NEAR(chem_array_get(&a,0,0),0,0); NEAR(chem_array_get(&a,3,0),0,0);
    CHECK(chem_array_assign_row(&a,2,little,2)==CHEM_SIZE); chem_array_destroy(&a);
}
static void axpy(size_t v)
{
    chem_array a=array(4,3,v,100),b=array(4,3,v,200); double expected=0,result=0;
    for(size_t i=0;i<4;++i) for(size_t j=0;j<3;++j) {
        double y=i*10.3+j*100.5,x=i*1.7+j*10.2; set(&a,i,j,y); set(&b,i,j,x); expected+=y+1.39*x;
    }
    REQUIRE(chem_array_axpy(&a,1.39,&b)==CHEM_OK);
    for(size_t i=0;i<4;++i) for(size_t j=0;j<3;++j) result+=chem_array_get(&a,i,j);
    NEAR(result,expected,1e-5); dispose(&a,&b,NULL);
}
static void scalar(size_t v,int mode)
{
    chem_array a=array(2,3,v,0);
    if (mode==0) chem_array_fill(&a,2); else if (mode==1) chem_array_max(&a,2); else chem_array_min(&a,2);
    for(size_t i=0;i<a.size;++i) NEAR(((double *)a.data)[i],mode==2 ? 0 : 2,0);
    if (mode) {
        chem_array_fill(&a,1); set(&a,1,1,3);
        if(mode==1) chem_array_max(&a,2); else chem_array_min(&a,2);
        const double mx[]={2,2,2,2,3,2},mn[]={1,1,1,1,2,1}; values(&a,mode==1 ? mx : mn);
    }
    chem_array_destroy(&a);
}
static void print_array(size_t v)
{
    chem_array a=array(2,3,v,0); set(&a,1,1,3); FILE *f=tmpfile(); REQUIRE(f!=NULL);
    REQUIRE(chem_array_print(f,&a)==CHEM_OK); rewind(f); char result[100]={0}; CHECK(fread(result,1,sizeof(result)-1,f)==12);
    CHECK(!strcmp(result,"0,0,0\n0,3,0\n")); fclose(f); chem_array_destroy(&a);
}
/* Row callbacks are the C equivalent of the reference's local lambdas. */
static double sum(const double *x,size_t n,void *ctx)
{ double s=0; for(size_t i=0;i<n;++i) s+=x[i]; return s*(ctx ? *(double *)ctx : 1); }
static double multiply(const double *x,size_t n,void *ctx)
{ (void)n;(void)ctx; return x[0]*x[1]; }
static double subtract(const double *x,size_t n,void *ctx)
{ (void)n;(void)ctx; return x[0]-x[1]; }
static double add_subtract(const double *x,size_t n,void *ctx)
{ (void)n;(void)ctx; return x[0]+x[1]-x[2]; }
static double reduce_op(double x,double acc,void *ctx)
{
    switch(*(int *)ctx) {
        case 0:return acc+x*x;
        case 1:return x>acc ? x : acc;
        case 2:return acc || x>10;
        case 3:return acc && isfinite(x);
        default:return acc+1;
    }
}
#define APPLY(destination,fn,context,...) do { chem_span sources_[]={__VA_ARGS__}; \
    chem_status s_=chem_span_map(destination,sizeof(sources_)/sizeof(*sources_),sources_,fn,context); if(s_) return s_; } while(0)
#define TRY(expr) do { chem_status s_=(expr); if(s_) return s_; } while(0)
typedef struct {
    int operation; double factor,accumulator; size_t calls;
    const size_t *captured_indices;
    const double *captured_values;
} program_context;
static chem_status program(size_t count,const chem_array_group *g,void *opaque)
{
    program_context *ctx=opaque; ++ctx->calls;
    chem_span a[5]; for(size_t j=0;j<5;++j) a[j]=chem_group_column(g[0],j);
    double *storage=calloc(2*g[0].width,sizeof(double)); if(!storage) return CHEM_MEMORY;
    chem_span t=chem_span_temporary(storage,g[0].width,g[0].active),u=chem_span_temporary(storage+g[0].width,g[0].width,g[0].active);
    /* Every temporary above is group-local, just as a native local array. */
    switch(ctx->operation) {
        case 0: /* summed columns, then scale */
            APPLY(t,sum,NULL,a[0],a[1],a[2]); APPLY(a[2],sum,&ctx->factor,t); break;
        case 1: /* matrix A column 0 plus matrix B column 2 */
            APPLY(t,sum,NULL,a[0],chem_group_const_column(g[1],2)); APPLY(a[1],sum,&ctx->factor,t); break;
        case 2: /* three-column sum multiplied by a vector */
            APPLY(t,sum,NULL,a[0],a[1],a[2]); APPLY(a[2],multiply,NULL,t,chem_group_vector(g[1])); break;
        case 3: /* first column plus vector, scaled */
            APPLY(t,sum,NULL,a[0],chem_group_vector(g[1])); APPLY(a[1],sum,&ctx->factor,t); break;
        case 4: /* multiple independent temporaries */
            APPLY(t,multiply,NULL,a[0],a[1]); APPLY(u,sum,NULL,a[0],a[1]);
            APPLY(a[2],sum,NULL,t,u); APPLY(a[3],subtract,NULL,t,u); APPLY(a[4],multiply,NULL,t,u); break;
        case 5: /* reuse saved column views */ {
            double two=2; APPLY(a[1],sum,&two,a[0]); APPLY(a[2],sum,NULL,a[0],a[1]); APPLY(a[3],multiply,NULL,a[2],a[1]); break;
        }
        case 6: /* const matrix read */
            APPLY(t,sum,NULL,a[0],a[1]);
            for(size_t i=0;i<a[2].active;++i) ctx->accumulator+=chem_span_get(a[2],i);
            break;
        case 7: /* vector scaling to column zero */
            APPLY(a[0],sum,&ctx->factor,chem_group_vector(g[1])); break;
        case 8: /* two vectors */
            APPLY(a[0],sum,NULL,chem_group_vector(g[1]),chem_group_vector(g[2])); break;
        case 9: /* two matrices share a vector */
            APPLY(a[1],sum,NULL,a[0],chem_group_vector(g[2]));
            APPLY(chem_group_column(g[1],1),sum,NULL,chem_group_const_column(g[1],0),chem_group_vector(g[2])); break;
        case 10: /* write a mutable vector using a const matrix */
            APPLY(chem_group_vector(g[1]),sum,&ctx->factor,a[0]); break;
        case 11:
            APPLY(t,sum,NULL,a[0],a[1]); APPLY(a[2],sum,NULL,t,chem_group_vector(g[1])); break;
        case 12: TRY(chem_span_fill(a[1],3.2)); break;
        case 13: TRY(chem_span_fill(chem_group_vector(g[1]),3.2)); break;
        case 14: TRY(chem_span_fill(t,9.9)); TRY(chem_span_copy(a[0],t)); break;
        case 15: TRY(chem_span_copy(a[1],chem_group_vector(g[1]))); break;
        case 16: TRY(chem_span_copy(chem_group_vector(g[1]),chem_span_const(a[1]))); break;
        case 17: TRY(chem_span_copy(a[0],a[1])); break;
        case 18: TRY(chem_span_copy(a[0],chem_span_const(a[1]))); break;
        case 19: TRY(chem_span_copy(t,chem_span_const(a[1]))); TRY(chem_span_fill(a[1],0)); TRY(chem_span_copy(a[1],t)); break;
        case 20: case 21: case 22: case 23: case 24: {
            int op=ctx->operation-20;
            for(size_t j=0;j<g[0].array->columns;++j) TRY(chem_span_reduce(chem_span_const(a[j]),op==4,&ctx->accumulator,reduce_op,&op));
            break;
        }
        case 25: { /* C callback context can retain read-only vector views. */
            const size_t *indexes=ctx->captured_indices; const double *captured=ctx->captured_values;
            for(size_t j=0;j<3;++j) for(size_t i=0;i<a[0].count;++i) TRY(chem_span_set(a[0],i,chem_span_get(a[0],i)+captured[indexes[j]]));
            for(size_t i=0;i<a[1].count;++i) TRY(chem_span_set(a[1],i,chem_span_get(a[1],i)+captured[indexes[1]]));
            break;
        }
        case 26: /* ForEach(a,b), then ForEach(a,b,c) */
            for(size_t j=0;j<g[0].array->columns;++j) {
                if(count==2) { APPLY(a[j],sum,NULL,a[j],chem_group_const_column(g[1],j)); }
                else { APPLY(a[j],add_subtract,NULL,a[j],chem_group_const_column(g[1],j),chem_group_const_column(g[2],j)); }
            }
            break;
        default: break;
    }
    free(storage); return CHEM_OK;
}
static void run(chem_array_function *f,program_context *ctx,size_t n,chem_array_argument *args)
{ CHECK(chem_array_function_create(f,n,args,program,ctx)==CHEM_OK); CHECK(chem_array_function_apply(f,n,args)==CHEM_OK); }
static void foreach(size_t v)
{
    chem_array a=array(4,3,v,0),b=array(4,3,v,0),c=array(4,3,v,0);
    for(size_t i=0;i<4;++i) for(size_t j=0;j<3;++j) { set(&a,i,j,i*10.3+j*100.5); set(&b,i,j,i*1.7+j*10.2); set(&c,i,j,i*19.5+j*32.2); }
    chem_array_argument args[]={chem_array_mutable(&a),chem_array_const(&b),chem_array_const(&c)}; chem_array_function f; program_context ctx={.operation=26};
    run(&f,&ctx,2,args);
    for(size_t i=0;i<4;++i) for(size_t j=0;j<3;++j) { NEAR(chem_array_get(&a,i,j),(i*10.3+j*100.5)+(i*1.7+j*10.2),1e-5); set(&a,i,j,i*10.3+j*100.5); }
    chem_array_function_destroy(&f); run(&f,&ctx,3,args);
    for(size_t i=0;i<4;++i) for(size_t j=0;j<3;++j) NEAR(chem_array_get(&a,i,j),(i*10.3+j*100.5)+(i*1.7+j*10.2)-(i*19.5+j*32.2),1e-5);
    chem_array_function_destroy(&f); dispose(&a,&b,&c);
}
static void array_function(size_t v,int with_vector)
{
    chem_array a=array(5,3,v,-1),b=vector(5,v,(double[]){100,200,300,400,500});
    for(size_t i=0;i<5;++i) for(size_t j=0;j<3;++j) set(&a,i,j,(double)i-2+10*j);
    program_context ctx={.operation=with_vector ? 2 : 0,.factor=4}; chem_array_function f;
    chem_array_argument args[]={chem_array_mutable(&a),chem_array_const(&b)}; run(&f,&ctx,with_vector ? 2 : 1,args);
    for(size_t i=0;i<5;++i) { NEAR(chem_array_get(&a,i,2),(24+3*i)*(with_vector ? 100*(i+1) : 4),0); NEAR(chem_array_get(&a,i,0),(double)i-2,0); NEAR(chem_array_get(&a,i,1),i+8,0); }
    if(!with_vector) {
        chem_array c=array(3,3,v,-1); args[0]=chem_array_mutable(&c); CHECK(chem_array_function_apply(&f,1,args)==CHEM_OK);
        values(&c,(double[]){-1,-1,-12,-1,-1,-12,-1,-1,-12}); chem_array_destroy(&c);
    }
    chem_array_function_destroy(&f); dispose(&a,&b,NULL);
}
static void multi_matrix(size_t v)
{
    chem_array a=array(3,2,v,1),b=array(3,3,v,2);
    for(size_t i=0;i<3;++i) { for(size_t j=0;j<2;++j) set(&a,i,j,(double)(i+10*j)); for(size_t j=0;j<3;++j) set(&b,i,j,(double)(2*i+20*j)); set(&b,i,2,(double)(4*i)); }
    program_context ctx={.operation=1,.factor=1}; chem_array_function f; chem_array_argument args[]={chem_array_mutable(&a),chem_array_const(&b)};
    run(&f,&ctx,2,args); values(&a,(double[]){0,0,1,5,2,10}); values(&b,(double[]){0,20,0,2,22,4,4,24,8});
    chem_array_function_destroy(&f); dispose(&a,&b,NULL);
}
static void different_rows(size_t v,int with_vector)
{
    chem_array a=array(3,with_vector ? 3 : 2,v,0),b=with_vector ? vector(3,v,NULL) : array(3,3,v,0);
    chem_array_argument args[]={chem_array_mutable(&a),chem_array_const(&b)}; chem_array_function f;
    program_context ctx={.operation=with_vector ? 3 : 1,.factor=with_vector ? 3 : 2};
    CHECK(chem_array_function_create(&f,2,args,program,&ctx)==CHEM_OK); dispose(&a,&b,NULL);
    a=array(5,with_vector ? 3 : 2,v,0); b=with_vector ? vector(5,v,NULL) : array(5,3,v,0);
    for(size_t i=0;i<5;++i) { set(&a,i,0,(double)(i+1)); set(&b,i,with_vector ? 0 : 2,(double)(10*i)); }
    CHECK(chem_array_function_apply(&f,2,args)==CHEM_OK);
    for(size_t i=0;i<5;++i) NEAR(chem_array_get(&a,i,1),(1+11*i)*ctx.factor,0);
    dispose(&a,&b,NULL);
    if(!with_vector) {
        a=array(2,2,v,0); b=array(2,3,v,0); set(&a,0,0,10);set(&a,1,0,20);set(&b,0,2,5);set(&b,1,2,15);
        CHECK(chem_array_function_apply(&f,2,args)==CHEM_OK); NEAR(chem_array_get(&a,0,1),30,0);NEAR(chem_array_get(&a,1,1),70,0); dispose(&a,&b,NULL);
    }
    chem_array_function_destroy(&f);
}
static void temporaries(size_t v)
{
    chem_array a=array(4,5,v,0); for(size_t i=0;i<4;++i) { set(&a,i,0,(double)(i+1)); set(&a,i,1,(double)(10*(i+1))); }
    chem_array_function f; chem_array_argument args[]={chem_array_mutable(&a)}; program_context ctx={.operation=4}; run(&f,&ctx,1,args);
    values(&a,(double[]){1,10,21,-1,110,2,20,62,18,880,3,30,123,57,2970,4,40,204,116,7040});
    chem_array_function_destroy(&f); chem_array_destroy(&a);
}
static void column_reuse(size_t v)
{
    chem_array a=array(3,4,v,0); for(size_t i=0;i<3;++i) set(&a,i,0,(double)(i+1));
    chem_array_function f; chem_array_argument args[]={chem_array_mutable(&a)}; program_context ctx={.operation=5}; run(&f,&ctx,1,args);
    values(&a,(double[]){1,2,3,6,2,4,6,24,3,6,9,54}); chem_array_function_destroy(&f); chem_array_destroy(&a);
}
static void function_reuse(size_t v)
{
    chem_array a=array(2,3,v,1); for(size_t i=0;i<2;++i) for(size_t j=0;j<3;++j) set(&a,i,j,(double)(i+j));
    chem_array_function f; chem_array_argument args[]={chem_array_mutable(&a)}; program_context ctx={.operation=0,.factor=2}; run(&f,&ctx,1,args);
    NEAR(chem_array_get(&a,0,2),6,0); NEAR(chem_array_get(&a,1,2),12,0); chem_array_destroy(&a);
    a=array(2,3,v,5); CHECK(chem_array_function_apply(&f,1,args)==CHEM_OK); NEAR(chem_array_get(&a,0,2),30,0); NEAR(chem_array_get(&a,1,2),30,0); chem_array_destroy(&a);
    a=array(2,3,v,0); set(&a,1,0,10); CHECK(chem_array_function_apply(&f,1,args)==CHEM_OK); NEAR(chem_array_get(&a,0,2),0,0); NEAR(chem_array_get(&a,1,2),20,0);
    chem_array_function_destroy(&f); chem_array_destroy(&a);
}
static void const_function(size_t v)
{
    chem_array a=array(3,4,v,0); for(size_t i=0;i<3;++i) for(size_t j=0;j<4;++j) set(&a,i,j,(double)(10*i+j));
    chem_array_function f; chem_array_argument args[]={chem_array_const(&a)}; program_context ctx={.operation=6}; run(&f,&ctx,1,args);
    NEAR(chem_array_get(&a,0,0),0,0); NEAR(chem_array_get(&a,1,2),12,0); NEAR(ctx.accumulator,36,0);
    chem_array_function_destroy(&f); chem_array_destroy(&a);
}
static chem_status no_op(size_t n,const chem_array_group *g,void *ctx)
{ (void)n;(void)g; if(ctx) ++*(size_t *)ctx; return CHEM_OK; }
static chem_status bad_column(size_t n,const chem_array_group *g,void *ctx)
{ (void)n;(void)ctx; return chem_span_fill(chem_group_column(g[0],5),1); }
static void empty_function(size_t v)
{
    chem_array a=array(0,3,v,1); chem_array_argument arg=chem_array_mutable(&a); chem_array_function f; size_t calls=0;
    CHECK(chem_array_function_create(&f,1,&arg,no_op,&calls)==CHEM_OK); CHECK(chem_array_function_apply(&f,1,&arg)==CHEM_OK); CHECK(calls==0);
    chem_array_function_destroy(&f); chem_array_destroy(&a); a=array(3,0,v,1);
    CHECK(chem_array_function_create(&f,1,&arg,no_op,&calls)==CHEM_OK); CHECK(chem_array_function_apply(&f,1,&arg)==CHEM_OK);
    chem_array_function_destroy(&f); chem_array_destroy(&a);
}
/* Validate creation and invocation independently. Different creation row counts
 * are intentional; all mismatches must be detected before callback execution. */
static void dimensions(size_t v,int which)
{
    chem_array a=array(3,3,v,1),b=array(4,3,v,2),c=vector(3,v,NULL);
    chem_array_argument args[]={chem_array_mutable(&a),chem_array_const(&b),chem_array_const(&c)};
    chem_array_function f; size_t n=2,calls=0;
    if(which==1) { /* bad column in callback */
        chem_array_destroy(&a); a=array(3,4,v,1); n=1;
        CHECK(chem_array_function_create(&f,n,args,bad_column,NULL)==CHEM_OK);
        CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE);
    } else if(which==2 || which==5) { /* wrong columns, then flexible rows */
        chem_array_destroy(&a); a=array(which==2 ? 3 : 4,which==2 ? 4 : 3,v,1); n=1;
        CHECK(chem_array_function_create(&f,n,args,no_op,&calls)==CHEM_OK);
        CHECK(chem_array_function_apply(&f,n,args)==CHEM_OK);
        chem_array_destroy(&a); a=array(4,which==2 ? 5 : 4,v,2); CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE);
        if(which==5) { chem_array_destroy(&a); a=array(4,2,v,2); CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE); }
        chem_array_destroy(&a); a=array(7,which==2 ? 4 : 3,v,1); CHECK(chem_array_function_apply(&f,n,args)==CHEM_OK);
    } else {
        if(which==3) { dispose(&a,&b,NULL); a=array(3,2,v,0); b=vector(3,v,NULL); }
        if(which==4) { dispose(&a,&b,NULL); a=array(3,2,v,0); b=array(3,3,v,0); }
        if(which>=6 && which<=10) {
            dispose(&a,&b,NULL);
            a=array(which==9 || which==10 ? 0 : 5,3,v,1);
            b=vector(which==6 ? 3 : which==7 ? 10 : which==9 ? 5 : 0,v,NULL);
        }
        if(which==11) { dispose(&a,&b,&c); a=array(5,3,v,1); b=vector(5,v,NULL); c=vector(3,v,NULL); n=3; }
        if(which==12) { dispose(&a,&b,&c); a=array(4,2,v,0); b=array(5,3,v,0); c=vector(4,v,NULL); n=3; }
        if(which==13) { dispose(&a,&b,&c); a=array(5,2,v,0); b=array(5,3,v,0); c=vector(4,v,NULL); n=3; }
        if(which==14) { dispose(&a,&b,NULL); a=array(3,2,v,0); b=vector(3,v,NULL); }
        CHECK(chem_array_function_create(&f,n,args,no_op,&calls)==CHEM_OK);
        if(which==0) { /* mismatch at creation, matching later */
            dispose(&a,&b,NULL); a=array(5,3,v,0);b=array(5,3,v,1); CHECK(chem_array_function_apply(&f,n,args)==CHEM_OK);
            dispose(&a,&b,NULL); a=array(3,3,v,0);b=array(4,3,v,1);
        }
        if(which==3 || which==4) { chem_array_destroy(&a); a=array(5,2,v,0); }
        if(which==10 || which==14) CHECK(chem_array_function_apply(&f,n,args)==CHEM_OK);
        else CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE);
        if(which==3) { dispose(&a,&b,NULL); a=array(3,2,v,0); b=vector(5,v,NULL); CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE); }
        if(which==14) { chem_array_destroy(&b); b=vector(5,v,NULL); CHECK(chem_array_function_apply(&f,n,args)==CHEM_SIZE); }
    }
    chem_array_function_destroy(&f); dispose(&a,&b,&c);
}
static void vector_operations(size_t v,int which)
{
    size_t rows=which==0 ? 5 : which==1 || which==5 ? 4 : 3;
    chem_array a=array(rows,which==0 || which==5 ? 3 : 2,v,0),b=vector(rows,v,NULL),c=vector(rows,v,NULL);
    program_context ctx={.operation=7,.factor=2};
    if(which==0) { ctx.operation=8; for(size_t i=0;i<rows;++i) { set(&b,i,0,(double)(i+1));set(&c,i,0,(double)(10*(i+1))); } }
    else if(which==1) { chem_array_destroy(&b); b=array(4,3,v,0); ctx.operation=9;
        for(size_t i=0;i<rows;++i) { set(&a,i,0,(double)(i+1));set(&b,i,0,(double)(10*i));set(&c,i,0,(double)(2*i)); } }
    else if(which==3) { ctx.operation=10;ctx.factor=3;for(size_t i=0;i<rows;++i) { set(&a,i,0,(double)(i+1));set(&b,i,0,(double)(5*(i+1))); } }
    else if(which==4) { ctx.factor=10;for(size_t i=0;i<rows;++i) set(&b,i,0,(double)(i+1)); }
    else if(which==5) { ctx.operation=11; for(size_t i=0;i<rows;++i) { set(&a,i,0,(double)(i+1));set(&a,i,1,(double)(10*(i+1)));set(&b,i,0,(double)(100*(i+1))); } }
    else if(which==6) { dispose(&a,&b,NULL); CHECK(chem_array_create(&a,3,2,v,CHEM_INT,0)==CHEM_OK); CHECK(chem_vector_create(&b,3,v,CHEM_INT,0)==CHEM_OK); for(size_t i=0;i<3;++i) set(&b,i,0,(double)(10*(i+1))); }
    else { for(size_t i=0;i<rows;++i) set(&b,i,0,(double)(10*(i+1))); }
    chem_array_argument args[]={which==3 ? chem_array_const(&a) : chem_array_mutable(&a),which==1 || which==3 ? chem_array_mutable(&b) : chem_array_const(&b),chem_array_const(&c)};
    size_t n=which<=1 ? 3 : 2; chem_array_function f; run(&f,&ctx,n,args);
    for(size_t i=0;i<rows;++i) {
        if(which==0) NEAR(chem_array_get(&a,i,0),11*(i+1),0);
        else if(which==1) { NEAR(chem_array_get(&a,i,1),1+3*i,0);NEAR(chem_array_get(&b,i,1),12*i,0); }
        else if(which==3) NEAR(chem_array_get(&b,i,0),3*(i+1),0);
        else if(which==4) NEAR(chem_array_get(&a,i,0),10*(i+1),0);
        else if(which==5) NEAR(chem_array_get(&a,i,2),111*(i+1),0);
        else NEAR(chem_array_get(&a,i,0),20*(i+1),0);
    }
    if(which==4) { chem_array_destroy(&b); b=vector(3,v,(double[]){5,6,7}); CHECK(chem_array_function_apply(&f,n,args)==CHEM_OK); for(size_t i=0;i<3;++i) NEAR(chem_array_get(&a,i,0),50+10*i,0); }
    chem_array_function_destroy(&f); dispose(&a,&b,&c);
}
static void fill_copy(size_t v,int copy)
{
    chem_array a=array(3,2,v,0),b=vector(3,v,(double[]){3.2,4.2,1.3});
    program_context ctx={0}; chem_array_function f;
    for(int i=0;i<(copy ? 5 : 3);++i) {
        ctx.operation=(copy ? 15 : 12)+i;
        size_t n=(copy ? i<2 : i==1) ? 2 : 1;
        chem_array_argument args[]={copy && i==1 ? chem_array_const(&a) : chem_array_mutable(&a),copy && i==0 ? chem_array_const(&b) : chem_array_mutable(&b)};
        if(copy && i==1) chem_array_fill(&b,0);
        if(copy && i==3) for(size_t row=0;row<3;++row) set(&a,row,0,0);
        run(&f,&ctx,n,args);
        for(size_t row=0;row<3;++row) {
            double expected=(double[]){3.2,4.2,1.3}[row];
            if(copy) {
                NEAR(chem_array_get(&a,row,1),expected,0);
                if(i==1) NEAR(chem_array_get(&b,row,0),expected,0);
                if(i==0) NEAR(chem_array_get(&a,row,0),0,0);
                if(i==2 || i==3) NEAR(chem_array_get(&a,row,0),expected,0);
            } else {
                if(i==0) { NEAR(chem_array_get(&a,row,0),0,0);NEAR(chem_array_get(&a,row,1),3.2,0); }
                if(i==1) NEAR(chem_array_get(&b,row,0),3.2,0);
                if(i==2) NEAR(chem_array_get(&a,row,0),9.9,0);
            }
        }
        chem_array_function_destroy(&f);
    }
    dispose(&a,&b,NULL);
}
static void reductions(size_t v,int op)
{
    for(int scenario=0;scenario<(op==2 || op==3 ? 2 : 1);++scenario) {
        chem_array a=array(3,op==4 ? 1 : 2,v,op==2 || op==3 ? 1 : 0);
        const double sums[]={1,2,3,4,5,6},maxs[]={1,8,5,2,3,4};
        if(op<=1) for(size_t i=0;i<3;++i) for(size_t j=0;j<2;++j) set(&a,i,j,(op ? maxs : sums)[2*i+j]);
        if(op==2 && scenario) set(&a,2,1,42);
        if(op==3 && scenario) set(&a,1,0,NAN);
        program_context ctx={.operation=20+op,.accumulator=op==1 ? -DBL_MAX : op==3 ? 1 : 0};
        chem_array_argument args[]={chem_array_const(&a)}; chem_array_function f;run(&f,&ctx,1,args);
        NEAR(ctx.accumulator,op==0 ? 91 : op==1 ? 8 : op==2 ? scenario : op==3 ? !scenario : 3,0);
        chem_array_function_destroy(&f);chem_array_destroy(&a);
    }
}
static void capture(size_t v)
{
    chem_array a=array(3,3,v,0); chem_array_argument args[]={chem_array_mutable(&a)}; chem_array_function f; program_context ctx={.operation=25};
    const size_t indexes[]={2,0,1};const double captured[]={10,20,30};
    ctx.captured_indices=indexes;ctx.captured_values=captured;
    run(&f,&ctx,1,args);values(&a,(double[]){60,10,0,60,10,0,60,10,0}); chem_array_function_destroy(&f);chem_array_destroy(&a);
}

static void iterate(void)
{
    chem_array a; REQUIRE(chem_array_create(&a,3,4,0,CHEM_INT,0)==CHEM_OK);
    int *data=a.data; for(size_t j=0;j<4;++j) data[4+j]=(int)(42+j);
    chem_array copy; REQUIRE(chem_array_clone(&copy,&a)==CHEM_OK);
    for(size_t j=0;j<4;++j) NEAR(chem_array_get(&copy,1,j),42+j,0);
    NEAR(chem_array_get(&a,0,3),0,0);NEAR(chem_array_get(&a,1,0),42,0);NEAR(chem_array_get(&a,1,1),43,0);NEAR(chem_array_get(&a,1,3),45,0);NEAR(chem_array_get(&a,2,0),0,0);
    dispose(&a,&copy,NULL);
}
TEST(matrix,Matrix,SmallMatrix) { small(0,0); }
TEST(matrix,Matrix,SmallConstMatrix) { small(0,1); }
TEST(matrix,Matrix,InitializeMatrix) { initialize(0); }
TEST(matrix,Matrix,InitializeConstMatrix) { initialize(0); }
TEST(matrix,Matrix,LoopOverMatrix) { loop(0,0); }
TEST(matrix,Matrix,LoopOverConstMatrix) { loop(0,1); }
TEST(matrix,Matrix,IterateOverMatrix) { iterate(); }
TEST(matrix,Matrix,Strides) { strides(0); }
TEST(matrix,Matrix,ConversionToVector) { to_vector(0,0); }
TEST(matrix,Matrix,ConstConversionToVector) { to_vector(0,1); }
TEST(matrix,Matrix,ConversionFromVector) { from_vector(0); }
TEST(matrix,Matrix,AssignmentFromVector) { assign_vector(0); }
TEST(matrix,Matrix,Axpy) { axpy(0); }
TEST(matrix,Matrix,ForEach) { foreach(0); }
TEST(matrix,Matrix,SetScaler) { scalar(0,0); }
TEST(matrix,Matrix,Max) { scalar(0,1); }
TEST(matrix,Matrix,Min) { scalar(0,2); }
TEST(matrix,Matrix,ArrayFunction) { array_function(0,0); }
TEST(matrix,Matrix,MultiMatrixArrayFunction) { multi_matrix(0); }
TEST(matrix,Matrix,MismatchedRowDimensions) { dimensions(0,0); }
TEST(matrix,Matrix,MismatchedColumnDimensions) { dimensions(0,1); }
TEST(matrix,Matrix,WrongMatrixDimensions) { dimensions(0,2); }
TEST(matrix,Matrix,MultipleTemporaries) { temporaries(0); }
TEST(matrix,Matrix,ColumnViewReuse) { column_reuse(0); }
TEST(matrix,Matrix,FunctionReusability) { function_reuse(0); }
TEST(matrix,Matrix,ConstMatrixFunction) { const_function(0); }
TEST(matrix,Matrix,EmptyMatrixFunction) { empty_function(0); }
TEST(matrix,Matrix,MultiMatrixDifferentRowsFromCreation) { different_rows(0,0); }
TEST(matrix,Matrix,MatrixVectorDifferentRowsFromCreation) { different_rows(0,1); }
TEST(matrix,Matrix,MismatchedRowsAtInvocation) { dimensions(0,3); }
TEST(matrix,Matrix,MultipleMatricesMismatchedRowsAtInvocation) { dimensions(0,4); }
TEST(matrix,Matrix,WrongColumnCountAtInvocation) { dimensions(0,5); }
TEST(matrix,Matrix,Print) { print_array(0); }
TEST(matrix,Matrix,VectorInMatrixFunction) { array_function(0,1); }
TEST(matrix,Matrix,VectorTooSmall) { dimensions(0,6); }
TEST(matrix,Matrix,VectorTooLarge) { dimensions(0,7); }
TEST(matrix,Matrix,EmptyVectorNonEmptyMatrix) { dimensions(0,8); }
TEST(matrix,Matrix,NonEmptyVectorEmptyMatrix) { dimensions(0,9); }
TEST(matrix,Matrix,EmptyVectorEmptyMatrix) { dimensions(0,10); }
TEST(matrix,Matrix,MultipleVectorsDifferentSizes) { dimensions(0,11); }
TEST(matrix,Matrix,MultipleVectorsSameSize) { vector_operations(0,0); }
TEST(matrix,Matrix,MultipleMatricesOneVector) { vector_operations(0,1); }
TEST(matrix,Matrix,MultipleMatricesDifferentRowsVector) { dimensions(0,12); }
TEST(matrix,Matrix,VectorSizeMatchesOneMatrixOnly) { dimensions(0,13); }
TEST(matrix,Matrix,ConstVector) { vector_operations(0,2); }
TEST(matrix,Matrix,MutableVector) { vector_operations(0,3); }
TEST(matrix,Matrix,FunctionReusabilityWithVectors) { vector_operations(0,4); }
TEST(matrix,Matrix,FunctionInvocationWithWrongSizedVector) { dimensions(0,14); }
TEST(matrix,Matrix,MixedVectorColumnViewRowVariable) { vector_operations(0,5); }
TEST(matrix,Matrix,IntegerVector) { vector_operations(0,6); }
TEST(matrix,Matrix,FunctionWithConstSignature) { vector_operations(0,7); }
TEST(matrix,Matrix,TestFill) { fill_copy(0,0); }
TEST(matrix,Matrix,TestCopy) { fill_copy(0,1); }
TEST(matrix,Matrix,ReduceSum) { reductions(0,0); }
TEST(matrix,Matrix,ReduceMax) { reductions(0,1); }
TEST(matrix,Matrix,ReduceLOr) { reductions(0,2); }
TEST(matrix,Matrix,ReduceLAnd) { reductions(0,3); }
TEST(matrix,Matrix,ReduceStrict) { reductions(0,4); }
TEST(matrix,Matrix,VectorCapture) { capture(0); }
TEST(vector_matrix,VectorMatrix,SmallVectorMatrix) { small(2,0); }
TEST(vector_matrix,VectorMatrix,SmallConstVectorMatrix) { small(4,1); }
TEST(vector_matrix,VectorMatrix,InitializeVectorMatrix) { initialize(1); }
TEST(vector_matrix,VectorMatrix,InitializeConstVectorMatrix) { initialize(2); }
TEST(vector_matrix,VectorMatrix,LoopOverVectorMatrix) { loop(2,0); }
TEST(vector_matrix,VectorMatrix,LoopOverConstVectorMatrix) { loop(1,1); }
TEST(vector_matrix,VectorMatrix,Strides) { strides(3); strides(4); }
TEST(vector_matrix,VectorMatrix,ConversionToVector) { to_vector(3,0); }
TEST(vector_matrix,VectorMatrix,ConstConversionToVector) { to_vector(1,1); }
TEST(vector_matrix,VectorMatrix,ConversionFromVector) { from_vector(2); }
TEST(vector_matrix,VectorMatrix,AssignmentFromVector) { assign_vector(2); }
TEST(vector_matrix,VectorMatrix,Axpy) { axpy(1); axpy(2); axpy(3); axpy(4); }
TEST(vector_matrix,VectorMatrix,ForEach) { foreach(1); foreach(2); foreach(3); foreach(4); }
TEST(vector_matrix,VectorMatrix,SetScaler) { scalar(1,0); scalar(2,0); scalar(3,0); scalar(4,0); }
TEST(vector_matrix,VectorMatrix,Max) { scalar(1,1); scalar(2,1); scalar(3,1); scalar(4,1); }
TEST(vector_matrix,VectorMatrix,Min) { scalar(1,2); scalar(2,2); scalar(3,2); scalar(4,2); }
TEST(vector_matrix,VectorMatrix,Print) { print_array(1); print_array(2); print_array(3); print_array(4); }
TEST(vector_matrix,VectorMatrix,ArrayFunction) { array_function(1,0); array_function(2,0); array_function(3,0); array_function(4,0); }
TEST(vector_matrix,VectorMatrix,MultiMatrixArrayFunction) { multi_matrix(1); multi_matrix(2); multi_matrix(3); multi_matrix(4); }
TEST(vector_matrix,VectorMatrix,MismatchedRowDimensions) { dimensions(1,0); dimensions(2,0); dimensions(3,0); dimensions(4,0); }
TEST(vector_matrix,VectorMatrix,MismatchedColumnDimensions) { dimensions(1,1); dimensions(2,1); dimensions(3,1); dimensions(4,1); }
TEST(vector_matrix,VectorMatrix,WrongMatrixDimensions) { dimensions(1,2); dimensions(2,2); dimensions(3,2); dimensions(4,2); }
TEST(vector_matrix,VectorMatrix,MultipleTemporaries) { temporaries(1); temporaries(2); temporaries(3); temporaries(4); }
TEST(vector_matrix,VectorMatrix,ColumnViewReuse) { column_reuse(1); column_reuse(2); column_reuse(3); column_reuse(4); }
TEST(vector_matrix,VectorMatrix,FunctionReusability) { function_reuse(1); function_reuse(2); function_reuse(3); function_reuse(4); }
TEST(vector_matrix,VectorMatrix,ConstMatrixFunction) { const_function(1); const_function(2); const_function(3); const_function(4); }
TEST(vector_matrix,VectorMatrix,EmptyMatrixFunction) { empty_function(1); empty_function(2); empty_function(3); empty_function(4); }
TEST(vector_matrix,VectorMatrix,MultiMatrixDifferentRowsFromCreation) { different_rows(1,0); different_rows(2,0); different_rows(3,0); different_rows(4,0); }
TEST(vector_matrix,VectorMatrix,MatrixVectorDifferentRowsFromCreation) { different_rows(1,1); different_rows(2,1); different_rows(3,1); different_rows(4,1); }
TEST(vector_matrix,VectorMatrix,MismatchedRowsAtInvocation) { dimensions(1,3); dimensions(2,3); dimensions(3,3); dimensions(4,3); }
TEST(vector_matrix,VectorMatrix,MultipleMatricesMismatchedRowsAtInvocation) { dimensions(1,4); dimensions(2,4); dimensions(3,4); dimensions(4,4); }
TEST(vector_matrix,VectorMatrix,WrongColumnCountAtInvocation) { dimensions(1,5); dimensions(2,5); dimensions(3,5); dimensions(4,5); }
TEST(vector_matrix,VectorMatrix,VectorInMatrixFunction) { array_function(1,1); array_function(2,1); array_function(3,1); array_function(4,1); }
TEST(vector_matrix,VectorMatrix,VectorTooSmall) { dimensions(1,6); dimensions(2,6); dimensions(3,6); dimensions(4,6); }
TEST(vector_matrix,VectorMatrix,VectorTooLarge) { dimensions(1,7); dimensions(2,7); dimensions(3,7); dimensions(4,7); }
TEST(vector_matrix,VectorMatrix,EmptyVectorNonEmptyMatrix) { dimensions(1,8); dimensions(2,8); dimensions(3,8); dimensions(4,8); }
TEST(vector_matrix,VectorMatrix,NonEmptyVectorEmptyMatrix) { dimensions(1,9); dimensions(2,9); dimensions(3,9); dimensions(4,9); }
TEST(vector_matrix,VectorMatrix,EmptyVectorEmptyMatrix) { dimensions(1,10); dimensions(2,10); dimensions(3,10); dimensions(4,10); }
TEST(vector_matrix,VectorMatrix,MultipleVectorsDifferentSizes) { dimensions(1,11); dimensions(2,11); dimensions(3,11); dimensions(4,11); }
TEST(vector_matrix,VectorMatrix,MultipleVectorsSameSize) { vector_operations(1,0); vector_operations(2,0); vector_operations(3,0); vector_operations(4,0); }
TEST(vector_matrix,VectorMatrix,MultipleMatricesOneVector) { vector_operations(1,1); vector_operations(2,1); vector_operations(3,1); vector_operations(4,1); }
TEST(vector_matrix,VectorMatrix,MultipleMatricesDifferentRowsVector) { dimensions(1,12); dimensions(2,12); dimensions(3,12); dimensions(4,12); }
TEST(vector_matrix,VectorMatrix,VectorSizeMatchesOneMatrixOnly) { dimensions(1,13); dimensions(2,13); dimensions(3,13); dimensions(4,13); }
TEST(vector_matrix,VectorMatrix,ConstVector) { vector_operations(1,2); vector_operations(2,2); vector_operations(3,2); vector_operations(4,2); }
TEST(vector_matrix,VectorMatrix,MutableVector) { vector_operations(1,3); vector_operations(2,3); vector_operations(3,3); vector_operations(4,3); }
TEST(vector_matrix,VectorMatrix,FunctionReusabilityWithVectors) { vector_operations(1,4); vector_operations(2,4); vector_operations(3,4); vector_operations(4,4); }
TEST(vector_matrix,VectorMatrix,FunctionInvocationWithWrongSizedVector) { dimensions(1,14); dimensions(2,14); dimensions(3,14); dimensions(4,14); }
TEST(vector_matrix,VectorMatrix,MixedVectorColumnViewRowVariable) { vector_operations(1,5); vector_operations(2,5); vector_operations(3,5); vector_operations(4,5); }
TEST(vector_matrix,VectorMatrix,IntegerVector) { vector_operations(1,6); vector_operations(2,6); vector_operations(3,6); vector_operations(4,6); }
TEST(vector_matrix,VectorMatrix,FunctionWithConstSignature) { vector_operations(1,7); vector_operations(2,7); vector_operations(3,7); vector_operations(4,7); }
TEST(vector_matrix,VectorMatrix,TestFill) { fill_copy(1,0); fill_copy(2,0); fill_copy(3,0); fill_copy(4,0); }
TEST(vector_matrix,VectorMatrix,TestCopy) { fill_copy(1,1); fill_copy(2,1); fill_copy(3,1); fill_copy(4,1); }
TEST(vector_matrix,VectorMatrix,ReduceSum) { reductions(1,0); reductions(2,0); reductions(3,0); reductions(4,0); }
TEST(vector_matrix,VectorMatrix,ReduceMax) { reductions(1,1); reductions(2,1); reductions(3,1); reductions(4,1); }
TEST(vector_matrix,VectorMatrix,ReduceLOr) { reductions(1,2); reductions(2,2); reductions(3,2); reductions(4,2); }
TEST(vector_matrix,VectorMatrix,ReduceLAnd) { reductions(1,3); reductions(2,3); reductions(3,3); reductions(4,3); }
TEST(vector_matrix,VectorMatrix,ReduceStrict) { reductions(1,4); reductions(2,4); reductions(3,4); reductions(4,4); }
TEST(vector_matrix,VectorMatrix,VectorCapture) { capture(1); capture(2); capture(3); capture(4); }
