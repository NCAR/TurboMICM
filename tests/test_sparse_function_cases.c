#include "sparse.h"
#include <stdio.h>

typedef struct { int mode; } Context;
static int block_fn(size_t block,ChemSparseMatrix *const *s,ChemMatrix *const *d,ChemVector *const *v,void *opaque)
{
    (void)d;
    int mode=((Context *)opaque)->mode;
    if (!mode) return 0;
    if (mode==12) {
        double *diagonal=chem_sparse_at(s[0],block,0,0);
        if (!diagonal) return 1;
        *diagonal=2*((double *)v[0]->data)[block];
        return 0;
    }
    double *a=chem_sparse_at(s[0],block,0,1); if (!a) return 1;
    if (mode==1) *a=2*((double *)v[0]->data)[block];
    else if (mode==2) *a=((double *)v[0]->data)[block]+((double *)v[1]->data)[block];
    else if (mode==3) *chem_sparse_at(s[1],block,0,1)=*a+((double *)v[0]->data)[block];
    else if (mode==4) ((double *)v[0]->data)[block]=*a*3;
    else if (mode==5) *chem_sparse_at(s[0],block,2,3)=(*a+*chem_sparse_at(s[0],block,1,2)+((double *)v[0]->data)[block])*2;
    else if (mode==6) *a=((int *)v[0]->data)[block]*1.5;
    else if (mode==7) *chem_sparse_at(s[0],block,1,2)=*a+((double *)v[0]->data)[block];
    else if (mode==8) *chem_sparse_at(s[0],block,1,1)=*a+*chem_sparse_at(s[1],block,0,1);
    else if (mode==9) *chem_sparse_at(s[0],block,1,2)=2**chem_sparse_at(s[1],block,0,1);
    else if (mode==10) *chem_sparse_at(s[0],block,2,3)=*a+*chem_sparse_at(s[1],block,1,2)+*chem_sparse_at(s[2],block,2,3);
    else if (mode==11) *a=3*((double *)v[0]->data)[block];
    return 0;
}
static int make_matrix(ChemSparseMatrix *m,size_t blocks,size_t width,int csc)
{
    ChemSparseBuilder b={0}; if(chem_sparse_builder_init(&b,4)) return 1;
    int error=chem_sparse_builder_add(&b,0,1)||chem_sparse_builder_add(&b,1,1)||chem_sparse_builder_add(&b,1,2)||chem_sparse_builder_add(&b,2,3)||chem_sparse_builder_add(&b,0,0)||chem_sparse_init(m,&b,blocks,width,csc,0);
    chem_sparse_builder_free(&b); return error;
}
static int test_case(int kind,size_t width,int csc)
{
    ChemSparseMatrix m[3]={{0}},other[3]={{0}}; ChemSparseFunction fn={0}; Context ctx={0};
    ChemSparseMatrix *s[]={&m[0],&m[1],&m[2]},*o[]={&other[0],&other[1],&other[2]};
    double data[]={1,2,3,4,5}, data2[]={4,5,6,7,8}; int integers[]={5,10,15};
    ChemVector x={data,3},y={data2,3},*vectors[]={&x,&y}; size_t ns=1,nv=1; int ok=1;
    if (kind==0 || kind==13 || kind==14) ctx.mode=1;
    if (kind==17) ctx.mode=12;
    if (kind==6 || kind==7) { nv=2; ctx.mode=kind==7 ? 2 : 0; }
    if (kind==8 || kind==9) { ns=2; ctx.mode=kind==8 ? 3 : 0; }
    if (kind==10) { ns=2; nv=2; }
    if (kind==11) ctx.mode=11;
    if (kind==12) ctx.mode=4;
    if (kind==15) ctx.mode=5;
    if (kind==16) { ctx.mode=6; x.data=integers; }
    if (kind==18) ctx.mode=7;
    if (kind==19) { ns=2; nv=0; ctx.mode=8; }
    if (kind==20) nv=0;
    if (kind==21) { ns=2; nv=0; ctx.mode=9; }
    if (kind==22) { ns=3; nv=0; ctx.mode=10; }
    for(size_t i=0;i<ns;++i) if(make_matrix(&m[i],3,width,csc)) return 0;
    if(chem_sparse_function_init(&fn,s,ns,NULL,0,nv,block_fn,&ctx)) return 0;
    if (kind==0 || kind==11 || kind==12 || kind==15) {
        for(int i=0;i<3;++i) { data[i]=kind==0 || kind==11 || kind==12 ? (i+1)*5 : (i+1)*100;
            *chem_sparse_at(&m[0],i,0,1)=i+1; *chem_sparse_at(&m[0],i,1,2)=i+2; *chem_sparse_at(&m[0],i,2,3)=i+3; }
        if (chem_sparse_function_call(&fn,s,1,NULL,0,vectors,1)) ok=0;
        for(int i=0;i<3;++i) {
            if (kind==0 && *chem_sparse_at(&m[0],i,0,1)!=(i+1)*10) ok=0;
            if (kind==11 && *chem_sparse_at(&m[0],i,0,1)!=(i+1)*15) ok=0;
            if (kind==12 && data[i]!=(i+1)*3) ok=0;
            if (kind==15 && *chem_sparse_at(&m[0],i,2,3)!=(2*i+3+(i+1)*100)*2) ok=0;
        }
    } else if (kind>=1 && kind<=6) {
        if (kind==4 || kind==5) { chem_sparse_free(&m[0]); make_matrix(&m[0],0,width,csc); }
        x.size=kind==1 ? 2 : kind==2 ? 4 : kind==3 || kind==5 ? 0 : 3; y.size=2;
        int status=chem_sparse_function_call(&fn,s,1,NULL,0,vectors,nv);
        if (kind==5 ? status!=0 : status==0) ok=0;
        if (kind==4) { size_t index; if (!chem_sparse_index(&m[0],0,0,1,&index)) ok=0; }
    } else if (kind==7) {
        if (chem_sparse_function_call(&fn,s,1,NULL,0,vectors,2)) ok=0;
        for(int i=0;i<3;++i) if(*chem_sparse_at(&m[0],i,0,1)!=5+2*i) ok=0;
    } else if (kind==8) {
        for(int i=0;i<3;++i) { *chem_sparse_at(&m[0],i,0,1)=i+1; data[i]=(i+1)*10; }
        if (chem_sparse_function_call(&fn,s,2,NULL,0,vectors,1)) ok=0;
        for(int i=0;i<3;++i) if(*chem_sparse_at(&m[1],i,0,1)!=(i+1)*11) ok=0;
    } else if (kind==9) {
        chem_sparse_free(&m[1]); make_matrix(&m[1],4,width,csc);
        if (!chem_sparse_function_call(&fn,s,2,NULL,0,vectors,1)) ok=0;
    } else if (kind==10) {
        y.size=2; if (!chem_sparse_function_call(&fn,s,2,NULL,0,vectors,2)) ok=0;
    } else if (kind==13) {
        if (chem_sparse_function_call(&fn,s,1,NULL,0,vectors,1)) ok=0;
        for(int i=0;i<3;++i) if(*chem_sparse_at(&m[0],i,0,1)!=(i+1)*2) ok=0;
        make_matrix(&other[0],3,width,csc); for(int i=0;i<3;++i) data[i]=(i+1)*10;
        if (chem_sparse_function_call(&fn,o,1,NULL,0,vectors,1)) ok=0;
        for(int i=0;i<3;++i) if(*chem_sparse_at(&other[0],i,0,1)!=(i+1)*20) ok=0;
    } else if (kind==14) {
        x.size=2; if (!chem_sparse_function_call(&fn,s,1,NULL,0,vectors,1)) ok=0;
    } else if (kind==16) {
        if (chem_sparse_function_call(&fn,s,1,NULL,0,vectors,1)) ok=0;
        for(int i=0;i<3;++i) if(*chem_sparse_at(&m[0],i,0,1)!=(i+1)*7.5) ok=0;
    } else if (kind==17) {
        if (chem_sparse_function_call(&fn,s,1,NULL,0,vectors,1)) ok=0;
        for(int i=0;i<3;++i) if(*chem_sparse_at(&m[0],i,0,0)!=(i+1)*2) ok=0;
    } else if (kind==18) {
        make_matrix(&other[0],5,width,csc); x.size=5;
        for(int i=0;i<5;++i) { data[i]=(i+1)*10; *chem_sparse_at(&other[0],i,0,1)=i+1; }
        if (chem_sparse_function_call(&fn,o,1,NULL,0,vectors,1)) ok=0;
        for(int i=0;i<5;++i) if(*chem_sparse_at(&other[0],i,1,2)!=(i+1)*11) ok=0;
    } else if (kind==19) {
        if (chem_sparse_function_call(&fn,s,2,NULL,0,NULL,0)) ok=0;
        chem_sparse_free(&m[1]); make_matrix(&m[1],4,width,csc);
        if (!chem_sparse_function_call(&fn,s,2,NULL,0,NULL,0)) ok=0;
    } else if (kind==20) {
        size_t index;
        if (!chem_sparse_index(&m[0],0,3,3,&index)) ok=0;
    } else if (kind==21) {
        if (chem_sparse_function_call(&fn,s,2,NULL,0,NULL,0)) ok=0;
        make_matrix(&other[0],4,width,csc);
        ChemSparseMatrix *first[]={&m[0],&other[0]},*second[]={&other[0],&m[0]};
        if (!chem_sparse_function_call(&fn,first,2,NULL,0,NULL,0) || !chem_sparse_function_call(&fn,second,2,NULL,0,NULL,0)) ok=0;
    } else if (kind==22) {
        if (chem_sparse_function_call(&fn,s,3,NULL,0,NULL,0)) ok=0;
        make_matrix(&other[0],4,width,csc);
        ChemSparseMatrix *first[]={&m[0],&m[1],&other[0]},*second[]={&other[0],&m[1],&m[2]};
        if (!chem_sparse_function_call(&fn,first,3,NULL,0,NULL,0) || !chem_sparse_function_call(&fn,second,3,NULL,0,NULL,0)) ok=0;
    }
    for(int i=0;i<3;++i) { chem_sparse_free(&m[i]); chem_sparse_free(&other[i]); }
    return ok;
}
static int all_widths(int kind,int grouped,int csc)
{
    if (!grouped) return test_case(kind,1,csc);
    for(size_t width=1;width<=4;++width) if(!test_case(kind,width,csc)) return 0;
    return 1;
}
int main(void)
{
    const char *prefix[]={"sparse_matrix_standard_compressed_row_ordering:SparseCompressedRowMatrix", "sparse_matrix_standard_compressed_column_ordering:SparseCompressedColumnMatrix", "sparse_matrix_vector_compressed_row_ordering:SparseVectorCompressedRowMatrix", "sparse_matrix_vector_compressed_column_ordering:SparseVectorCompressedColumnMatrix"};
    const char *names[]={"VectorInSparseMatrixFunction","VectorTooSmall","VectorTooLarge","EmptyVectorNonEmptySparseMatrix","NonEmptyVectorEmptySparseMatrix","EmptyVectorEmptySparseMatrix","MultipleVectorsDifferentSizes","MultipleVectorsSameSize","MultipleSparseMatricesOneVector","MultipleSparseMatricesDifferentBlocksVector","VectorSizeMatchesOneSparseMatrixOnly","ConstVectorSparse","MutableVectorSparse","FunctionReusabilityWithVectorsSparse","FunctionInvocationWithWrongSizedVectorSparse","MixedVectorBlockViewBlockVariable","IntegerVectorSparse","FunctionWithConstSignatureSparse","SparseMatrixVectorDifferentBlocksFromCreation","MismatchedBlockDimensions","MismatchedElementDimensions","MismatchedBlocksAtInvocation","MultipleSparseMatricesMismatchedBlocksAtInvocation"};
    for(int layout=0;layout<4;++layout) for(int kind=0;kind<23;++kind) {
        /* These two names exist only for standard layouts. */
        if (layout>=2 && (kind==19 || kind==20)) continue;
        printf("%s.%s %d\n",prefix[layout],names[kind],all_widths(kind,layout>=2,layout%2));
    }
    return 0;
}
