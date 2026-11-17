#include "sparse.h"
#include <stdio.h>
#include <string.h>

typedef struct { int mode; size_t indices[6]; } Context;
static int block_fn(size_t b,ChemSparseMatrix *const *s,ChemMatrix *const *d,ChemVector *const *v,void *opaque)
{
    Context *ctx=opaque; int mode=ctx->mode;
    if (mode==7 || mode==8) return 0;
    if (mode==11) { ((double *)v[0]->data)[b]=3.2; return 0; }
    if (mode==10) {
        for(int i=0;i<6;++i) { double *p=chem_sparse_at_vector_index(s[0],b,ctx->indices[i]); if(!p) return 1; *p=b+10*i; }
        return 0;
    }
    if (mode==0 || mode==4) {
        double *a=chem_sparse_at(s[0],b,0,0), *x=chem_sparse_at(s[0],b,1,1), *y=chem_sparse_at(s[0],b,2,2);
        double *out=chem_sparse_at(s[0],b,mode==0 ? 2 : 2,mode==0 ? 3 : 2);
        if(!a||!x||!y||!out) return 1;
        double sum=*a+*x+*y+(mode==0 ? *out : 0); *out=2*sum; return 0;
    }
    double *a=chem_sparse_at(s[0],b,0,1); if(!a) return 1;
    if (mode==1) { double *out=chem_sparse_at(s[0],b,1,2), *x=chem_sparse_at(s[1],b,2,2); if(!out||!x) return 1; *out=*a+*x; }
    else if (mode==2) {
        double *x=chem_sparse_at(s[0],b,1,2), *sum=chem_sparse_at(s[0],b,2,3), *diff=chem_sparse_at(s[0],b,3,4), *product=chem_sparse_at(s[0],b,4,4);
        if(!x||!sum||!diff||!product) return 1;
        double t1=*a**x,t2=*a+*x; *sum=t1+t2; *diff=t1-t2; *product=t1*t2;
    } else if(mode==3) {
        double *x=chem_sparse_at(s[0],b,1,2), *y=chem_sparse_at(s[0],b,2,3), *z=chem_sparse_at(s[0],b,3,3);
        if(!x||!y||!z) return 1;
        *x=2**a; *y=*a+*x; *z=*y**x;
    } else if(mode==5) {
        double *x=chem_sparse_at(s[0],b,1,1), *out=chem_sparse_at(s[0],b,2,2), *other=chem_sparse_at(s[1],b,0,0);
        if(!x||!out||!other) return 1;
        *out=*a+*x+*other;
    } else if(mode==6) {
        double *x=chem_sparse_at(s[0],b,1,2), *out=chem_sparse_at(s[0],b,3,3);
        if(!x||!out) return 1;
        *out=*a+*x+*chem_matrix_at(d[0],b,0)+*chem_matrix_at(d[0],b,1);
    } else if(mode==9) {
        double *x=chem_sparse_at(s[0],b,1,2), *out=chem_sparse_at(s[0],b,2,3), *other=chem_sparse_at(s[1],b,2,3);
        if(!x||!out||!other) return 1;
        *out=*a+*x+*other;
    }
    return 0;
}
static int make(ChemSparseMatrix *m,size_t n,const size_t (*points)[2],size_t count,size_t blocks,size_t width,int csc,double initial)
{
    ChemSparseBuilder b={0}; if(chem_sparse_builder_init(&b,n)) return 1;
    int status=0; for(size_t i=0;i<count;++i) status|=chem_sparse_builder_add(&b,points[i][0],points[i][1]);
    status|=chem_sparse_init(m,&b,blocks,width,csc,initial); chem_sparse_builder_free(&b); return status;
}
static int test_case(int kind,size_t width,int csc)
{
    const size_t diagonal[][2]={{0,0},{1,1},{2,2},{2,3},{3,3}},
        pair_a[][2]={{0,1},{1,2}},pair_b[][2]={{0,0},{1,1},{2,2}},
        chain[][2]={{0,1},{1,2},{2,3},{3,3}},five[][2]={{0,1},{1,2},{2,3},{3,4},{4,4}},
        structure_a[][2]={{0,1},{1,1},{2,2}},structure_b[][2]={{0,0},{1,2},{2,3},{3,4}},
        block_views[][2]={{0,0},{0,2},{1,0},{1,1},{2,1},{2,2}},copy_points[][2]={{0,1},{2,0}};
    ChemSparseMatrix m[2]={{0}},other[2]={{0}}; ChemSparseMatrix *s[]={&m[0],&m[1]},*o[]={&other[0],&other[1]};
    ChemSparseFunction fn={0}; ChemMatrix dense={0},*ds[]={&dense}; Context ctx={kind,{0}}; int ok=1;
    if(kind==0) {
        make(&m[0],4,diagonal,5,3,width,csc,1); *chem_sparse_at(&m[0],1,1,1)=2; *chem_sparse_at(&m[0],2,2,2)=3; *chem_sparse_at(&m[0],2,2,3)=4; *chem_sparse_at(&m[0],2,3,3)=5;
        if(chem_sparse_function_init(&fn,s,1,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,NULL,0,NULL,0)) ok=0;
        double expected[]={8,10,18}; for(int i=0;i<3;++i) if(*chem_sparse_at(&m[0],i,2,3)!=expected[i] || *chem_sparse_at(&m[0],i,0,0)!=1 || *chem_sparse_at(&m[0],i,1,1)!=(i==1?2:1) || *chem_sparse_at(&m[0],i,2,2)!=(i==2?3:1) || *chem_sparse_at(&m[0],i,3,3)!=(i==2?5:1)) ok=0;
        make(&other[0],4,diagonal,5,3,width,csc,-1); if(chem_sparse_function_call(&fn,o,1,NULL,0,NULL,0)) ok=0;
        for(int i=0;i<3;++i) if(*chem_sparse_at(&other[0],i,2,3)!=-8 || *chem_sparse_at(&other[0],i,0,0)!=-1 || *chem_sparse_at(&other[0],i,1,1)!=-1 || *chem_sparse_at(&other[0],i,2,2)!=-1 || *chem_sparse_at(&other[0],i,3,3)!=-1) ok=0;
    } else if(kind==1) {
        make(&m[0],3,pair_a,2,3,width,csc,0); make(&m[1],3,pair_b,3,3,width,csc,0);
        for(int i=0;i<3;++i) { *chem_sparse_at(&m[0],i,0,1)=i+10; *chem_sparse_at(&m[0],i,1,2)=i*2+20; *chem_sparse_at(&m[1],i,0,0)=i*4; *chem_sparse_at(&m[1],i,1,1)=i*3; *chem_sparse_at(&m[1],i,2,2)=i*5; }
        if(chem_sparse_function_init(&fn,s,2,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,2,NULL,0,NULL,0)) ok=0;
        for(int i=0;i<3;++i) if(*chem_sparse_at(&m[0],i,1,2)!=10+i*6 || *chem_sparse_at(&m[0],i,0,1)!=10+i || *chem_sparse_at(&m[1],i,0,0)!=i*4) ok=0;
    } else if(kind==2 || kind==3) {
        make(&m[0],kind==2 ? 5 : 4,kind==2 ? five : chain,kind==2 ? 5 : 4,kind==2 ? 4 : 3,width,csc,0);
        for(size_t i=0;i<m[0].blocks;++i) { *chem_sparse_at(&m[0],i,0,1)=i+1; if(kind==2) *chem_sparse_at(&m[0],i,1,2)=(i+1)*10; }
        if(chem_sparse_function_init(&fn,s,1,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,NULL,0,NULL,0)) ok=0;
        for(size_t i=0;i<m[0].blocks;++i) { double x=i+1;
            if(kind==2 && (*chem_sparse_at(&m[0],i,2,3)!=10*x*x+11*x || *chem_sparse_at(&m[0],i,3,4)!=10*x*x-11*x || *chem_sparse_at(&m[0],i,4,4)!=110*x*x*x)) ok=0;
            if(kind==3 && (*chem_sparse_at(&m[0],i,1,2)!=2*x || *chem_sparse_at(&m[0],i,2,3)!=3*x || *chem_sparse_at(&m[0],i,3,3)!=6*x*x)) ok=0;
        }
    } else if(kind==4) {
        make(&m[0],3,pair_b,3,2,width,csc,1); make(&other[0],3,pair_b,3,2,width,csc,5); make(&other[1],3,pair_b,3,2,width,csc,0);
        for(int b=0;b<2;++b) for(int i=0;i<3;++i) *chem_sparse_at(&m[0],b,i,i)=b+i;
        if(chem_sparse_function_init(&fn,s,1,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,NULL,0,NULL,0)) ok=0;
        if(*chem_sparse_at(&m[0],0,2,2)!=6 || *chem_sparse_at(&m[0],1,2,2)!=12) ok=0;
        if(chem_sparse_function_call(&fn,o,1,NULL,0,NULL,0) || *chem_sparse_at(&other[0],0,2,2)!=30 || *chem_sparse_at(&other[0],1,2,2)!=30) ok=0;
        *chem_sparse_at(&other[1],1,0,0)=10; o[0]=&other[1];
        if(chem_sparse_function_call(&fn,o,1,NULL,0,NULL,0) || *chem_sparse_at(&other[1],0,2,2)!=0 || *chem_sparse_at(&other[1],1,2,2)!=20) ok=0;
    } else if(kind==5) {
        make(&m[0],3,structure_a,3,4,width,csc,0); make(&m[1],5,structure_b,4,4,width,csc,0);
        for(int i=0;i<4;++i) { *chem_sparse_at(&m[0],i,0,1)=i*2+1; *chem_sparse_at(&m[0],i,1,1)=i*3+2; *chem_sparse_at(&m[0],i,2,2)=i*5+3;
            *chem_sparse_at(&m[1],i,0,0)=i+10; *chem_sparse_at(&m[1],i,1,2)=i+20; *chem_sparse_at(&m[1],i,2,3)=i+30; *chem_sparse_at(&m[1],i,3,4)=i+40; }
        if(chem_sparse_function_init(&fn,s,2,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,2,NULL,0,NULL,0)) ok=0;
        for(int i=0;i<4;++i) if(*chem_sparse_at(&m[0],i,2,2)!=13+6*i || *chem_sparse_at(&m[0],i,0,1)!=1+2*i) ok=0;
        if(*chem_sparse_at(&m[1],0,0,0)!=10 || *chem_sparse_at(&m[1],1,1,2)!=21 || *chem_sparse_at(&m[1],2,2,3)!=32 || *chem_sparse_at(&m[1],3,3,4)!=43) ok=0;
    } else if(kind==6) {
        make(&m[0],4,chain,4,3,width,csc,0); chem_matrix_init(&dense,3,4,width,width!=1,0);
        for(int b=0;b<3;++b) { *chem_sparse_at(&m[0],b,0,1)=b+1; *chem_sparse_at(&m[0],b,1,2)=b+10; *chem_sparse_at(&m[0],b,2,3)=b+20; *chem_sparse_at(&m[0],b,3,3)=b+30;
            for(int c=0;c<4;++c) *chem_matrix_at(&dense,b,c)=b*100+c*10; }
        if(chem_sparse_function_init(&fn,s,1,ds,1,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,ds,1,NULL,0)) ok=0;
        for(int b=0;b<3;++b) if(*chem_sparse_at(&m[0],b,3,3)!=21+b*202) ok=0;
        if(*chem_sparse_at(&m[0],0,0,1)!=1 || *chem_sparse_at(&m[0],1,1,2)!=11 || *chem_sparse_at(&m[0],2,2,3)!=22 || *chem_matrix_at(&dense,0,0)!=0 || *chem_matrix_at(&dense,1,2)!=120 || *chem_matrix_at(&dense,2,3)!=230) ok=0;
    } else if(kind==7) {
        make(&m[0],4,chain,4,3,width,csc,0); for(int b=0;b<3;++b) { *chem_sparse_at(&m[0],b,0,1)=b+1; *chem_sparse_at(&m[0],b,1,2)=b+10; }
        if(chem_sparse_function_init(&fn,s,1,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,NULL,0,NULL,0) || *chem_sparse_at(&m[0],0,0,1)!=1 || *chem_sparse_at(&m[0],1,1,2)!=11) ok=0;
    } else if(kind==8) {
        make(&m[0],4,NULL,0,3,width,csc,0);
        if(chem_sparse_function_init(&fn,s,1,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,NULL,0,NULL,0)) ok=0;
        chem_sparse_free(&m[0]); make(&m[0],4,chain,1,0,width,csc,0);
        if(chem_sparse_function_init(&fn,s,1,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,NULL,0,NULL,0)) ok=0;
    } else if(kind==9) {
        make(&m[0],4,chain,3,3,width,csc,0); make(&m[1],4,chain,3,3,width,csc,0);
        if(chem_sparse_function_init(&fn,s,2,NULL,0,0,block_fn,&ctx)) ok=0;
        make(&other[0],4,chain,3,4,width,csc,0); make(&other[1],4,chain,3,4,width,csc,0);
        for(int i=0;i<4;++i) { *chem_sparse_at(&other[0],i,0,1)=i+1; *chem_sparse_at(&other[0],i,1,2)=i+10; *chem_sparse_at(&other[1],i,2,3)=i+100; }
        if(chem_sparse_function_call(&fn,o,2,NULL,0,NULL,0)) ok=0;
        for(int i=0;i<4;++i) if(*chem_sparse_at(&other[0],i,2,3)!=111+3*i) ok=0;
    } else if(kind==10) {
        make(&m[0],3,block_views,6,3,width,csc,0);
        for(int i=0;i<6;++i) if(chem_sparse_index(&m[0],0,block_views[i][0],block_views[i][1],&ctx.indices[i])) ok=0;
        if(chem_sparse_function_init(&fn,s,1,NULL,0,0,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,NULL,0,NULL,0)) ok=0;
        for(int b=0;b<3;++b) for(int i=0;i<6;++i) if(*chem_sparse_at(&m[0],b,block_views[i][0],block_views[i][1])!=b+10*i) ok=0;
    } else if(kind==11 || kind==12) {
        make(&m[0],3,copy_points,2,3,width,csc,0); double values[]={3.2,4.2,1.3},output[3]={0};
        if(kind==11) {
            if(chem_sparse_fill_block_view(&m[0],0,1,3.2)) ok=0;
            for(int b=0;b<3;++b) if(*chem_sparse_at(&m[0],b,0,1)!=3.2 || *chem_sparse_at(&m[0],b,2,0)!=0) ok=0;
            ChemVector vec={output,3}, *vectors[]={&vec}; ctx.mode=11;
            if(chem_sparse_function_init(&fn,s,1,NULL,0,1,block_fn,&ctx)||chem_sparse_function_call(&fn,s,1,NULL,0,vectors,1)) ok=0;
            for(int b=0;b<3;++b) if(output[b]!=3.2) ok=0;
            if(chem_sparse_fill_block_view(&m[0],2,0,9.9)) ok=0;
            for(int b=0;b<3;++b) if(*chem_sparse_at(&m[0],b,2,0)!=9.9) ok=0;
        } else {
            if(chem_sparse_vector_to_block(&m[0],0,1,values,3)||chem_sparse_block_to_vector(&m[0],0,1,output,3)||chem_sparse_copy_block_view(&m[0],2,0,&m[0],0,1)) ok=0;
            for(int b=0;b<3;++b) if(output[b]!=values[b] || *chem_sparse_at(&m[0],b,2,0)!=values[b] || *chem_sparse_at(&m[0],b,0,1)!=values[b]) ok=0;
            if(chem_sparse_fill_block_view(&m[0],0,1,0)||chem_sparse_vector_to_block(&m[0],0,1,output,3)) ok=0;
            chem_matrix_init(&dense,3,2,width,width!=1,0); for(int b=0;b<3;++b) *chem_matrix_at(&dense,b,1)=(b+1)*100;
            for(int b=0;b<3;++b) output[b]=*chem_matrix_at(&dense,b,1);
            if(chem_sparse_vector_to_block(&m[0],2,0,output,3)) ok=0;
            for(int b=0;b<3;++b) if(*chem_sparse_at(&m[0],b,2,0)!=(b+1)*100) ok=0;
            if(chem_sparse_block_to_vector(&m[0],0,1,output,3)) ok=0;
            for(int b=0;b<3;++b) *chem_matrix_at(&dense,b,0)=output[b];
            for(int b=0;b<3;++b) if(*chem_matrix_at(&dense,b,0)!=values[b]) ok=0;
            make(&m[1],3,pair_a+1,1,3,width,csc,0);
            if(chem_sparse_copy_block_view(&m[1],1,2,&m[0],0,1)) ok=0;
            for(int b=0;b<3;++b) if(*chem_sparse_at(&m[1],b,1,2)!=values[b]) ok=0;
        }
    }
    for(int i=0;i<2;++i) { chem_sparse_free(&m[i]); chem_sparse_free(&other[i]); } chem_matrix_free(&dense); return ok;
}
static int all_widths(int kind,int grouped,int csc)
{
    if(!grouped || kind==6) return test_case(kind,1,csc);
    for(size_t width=1;width<=4;++width) if(!test_case(kind,width,csc)) return 0;
    return 1;
}
static int incompatible(int mode,int csc)
{
    const size_t points[][2]={{0,1},{1,2}};
    size_t first[]={2,2,4,4}, second[]={1,3,2,10};
    if(mode==0) { first[0]=2; first[1]=3; first[2]=4; }
    if(mode==2) { first[0]=1; first[1]=2; first[2]=3; first[3]=4; second[0]=2; second[1]=3; second[2]=4; second[3]=10; }
    int count=mode==0?3:4,ok=1;
    for(int i=0;i<count;++i) {
        ChemSparseMatrix a={0},b={0},*s[]={&a,&b}; ChemMatrix d={0},*dense[]={&d}; ChemSparseFunction f={0}; Context ctx={8,{0}};
        make(&a,4,points,2,3,first[i],csc,0);
        if(mode==2) { make(&b,4,points,2,3,second[i],csc,0); if(!chem_sparse_function_init(&f,s,2,NULL,0,0,block_fn,&ctx)) ok=0; }
        else { size_t dw=mode==0?1:second[i]; chem_matrix_init(&d,3,4,dw,dw!=1,0); if(!chem_sparse_function_init(&f,s,1,dense,1,0,block_fn,&ctx)) ok=0; }
        chem_sparse_free(&a);chem_sparse_free(&b);chem_matrix_free(&d);
    }
    return ok;
}
int main(void)
{
    const char *prefix[]={"sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix", "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix", "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix", "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix"};
    const char *names[]={"ArrayFunction","MultiMatrixArrayFunction","MultipleTemporaries","BlockViewReuse","FunctionReusability","TwoSparseMatricesDifferentStructure","SparseAndDenseMatrixFunction","ConstSparseMatrixFunction","EmptySparseMatrixFunction","MultipleSparseMatricesDifferentBlocksFromCreation","GetBlockViewByVectorIndex","TestFill","TestCopy"};
    for(int layout=0;layout<4;++layout) for(int kind=0;kind<13;++kind) printf("%s.%s %d\n",prefix[layout],names[kind],all_widths(kind,layout>=2,layout%2));
    for(int layout=2;layout<4;++layout) {
        int ok=1; size_t widths[]={1,2,3,4,10}; for(int i=0;i<5;++i) ok=ok && test_case(6,widths[i],layout%2);
        printf("%s.SparseAndVectorMatrixFunction %d\n",prefix[layout],ok);
        printf("%s.IncompatibleOrdering %d\n",prefix[layout],incompatible(0,layout%2));
        printf("%s.IncompatibleVectorOrdering %d\n",prefix[layout],incompatible(1,layout%2));
        printf("%s.IncompatibleSparseOrdering %d\n",prefix[layout],incompatible(2,layout%2));
    }
    return 0;
}
