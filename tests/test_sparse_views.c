#undef NDEBUG
#include <assert.h>
#include "sparse_matrix.h"
static ChemColumnView block(ChemSparseMatrix *m,size_t r,size_t c)
{ ChemColumnView v; assert(chem_sparse_matrix_block(m,r,c,&v)==CHEM_MATRIX_OK); return v; }
static ChemColumnView column(ChemMatrix *m,size_t c)
{ ChemColumnView v; assert(chem_matrix_column(m,c,&v)==CHEM_MATRIX_OK); return v; }
static void calc(size_t row,double *const *v,void *context)
{
    (void)row; int mode=*(int *)context;
    switch (mode) {
    case 0: *v[3]=*v[0]+*v[1]+*v[2]; break;
    case 1: *v[2]=*v[0]+*v[1]; break;
    case 2: *v[2]=*v[0]**v[1]; break;
    case 3: *v[2]=*v[0]-*v[1]; break;
    case 4: *v[1]=2**v[0]; break;
    case 5: *v[3]=2*(*v[0]+*v[1]+*v[2]+*v[3]); break;
    case 6: *v[4]=*v[0]+*v[1]+*v[2]+*v[3]; break;
    }
}
static void apply(ChemColumnView *views,size_t count,int mode)
{ assert(chem_columns_apply(views,count,calc,&mode)==CHEM_MATRIX_OK); }
static void run(size_t width,ChemSparseOrdering order)
{
    ChemSparseMatrix a={0},b={0}; ChemMatrix dense={0};
    const ChemSparseEntry diag[]={{0,0},{1,1},{2,2},{2,3},{3,3}};
    assert(chem_sparse_matrix_init(&a,4,3,width,order,diag,5,1)==CHEM_MATRIX_OK);
    *chem_sparse_matrix_at(&a,1,1,1)=2; *chem_sparse_matrix_at(&a,2,2,2)=3;
    *chem_sparse_matrix_at(&a,2,2,3)=4; *chem_sparse_matrix_at(&a,2,3,3)=5;
    ChemColumnView v[5]={block(&a,0,0),block(&a,1,1),block(&a,2,2),block(&a,2,3)};
    ChemColumnView by_index; size_t idx=chem_sparse_matrix_index(&a,0,2,3);
    assert(chem_sparse_matrix_block_by_index(&a,idx,&by_index)==CHEM_MATRIX_OK);
    assert(by_index.column==v[3].column);
    apply(v,4,5);
    assert(*chem_sparse_matrix_at(&a,0,2,3)==8 && *chem_sparse_matrix_at(&a,1,2,3)==10 && *chem_sparse_matrix_at(&a,2,2,3)==18);
    for (size_t i=0;i<3;++i) {
        assert(*chem_sparse_matrix_at(&a,i,0,0)==1 && *chem_sparse_matrix_at(&a,i,1,1)==(i==1 ? 2 : 1));
        assert(*chem_sparse_matrix_at(&a,i,2,2)==(i==2 ? 3 : 1) && *chem_sparse_matrix_at(&a,i,3,3)==(i==2 ? 5 : 1));
    }
    chem_sparse_matrix_fill(&a,-1); apply(v,4,5);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,2,3)==-8 && *chem_sparse_matrix_at(&a,i,0,0)==-1 &&
        *chem_sparse_matrix_at(&a,i,1,1)==-1 && *chem_sparse_matrix_at(&a,i,2,2)==-1 && *chem_sparse_matrix_at(&a,i,3,3)==-1);
    const ChemSparseEntry aa[]={{0,1},{1,2}},bb[]={{0,0},{1,1},{2,2}};
    assert(chem_sparse_matrix_init(&a,3,3,width,order,aa,2,0)==CHEM_MATRIX_OK);
    assert(chem_sparse_matrix_init(&b,3,3,width,order,bb,3,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) { *chem_sparse_matrix_at(&a,i,0,1)=i+10; *chem_sparse_matrix_at(&a,i,1,2)=i*2+20;
        *chem_sparse_matrix_at(&b,i,0,0)=i*4; *chem_sparse_matrix_at(&b,i,2,2)=i*5; }
    v[0]=block(&a,0,1); v[1]=block(&b,2,2); v[2]=block(&a,1,2); apply(v,3,1);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,1,2)==10+6*i && *chem_sparse_matrix_at(&a,i,0,1)==10+i && *chem_sparse_matrix_at(&b,i,0,0)==4*i);
    const ChemSparseEntry five[]={{0,1},{1,2},{2,3},{3,4},{4,4}};
    assert(chem_sparse_matrix_init(&a,5,4,width,order,five,5,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) { *chem_sparse_matrix_at(&a,i,0,1)=i+1; *chem_sparse_matrix_at(&a,i,1,2)=10*(i+1); }
    double t1[4],t2[4];
    v[0]=block(&a,0,1); v[1]=block(&a,1,2); v[2]=chem_vector_column(t1,4); apply(v,3,2);
    v[2]=chem_vector_column(t2,4); apply(v,3,1);
    v[0]=chem_vector_column(t1,4); v[1]=chem_vector_column(t2,4); v[2]=block(&a,2,3); apply(v,3,1);
    v[2]=block(&a,3,4); apply(v,3,3); v[2]=block(&a,4,4); apply(v,3,2);
    for (size_t i=0;i<4;++i) { double x=i+1; assert(*chem_sparse_matrix_at(&a,i,2,3)==10*x*x+11*x);
        assert(*chem_sparse_matrix_at(&a,i,3,4)==10*x*x-11*x && *chem_sparse_matrix_at(&a,i,4,4)==110*x*x*x); }
    const ChemSparseEntry four[]={{0,1},{1,2},{2,3},{3,3}};
    assert(chem_sparse_matrix_init(&a,4,3,width,order,four,4,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) *chem_sparse_matrix_at(&a,i,0,1)=i+1;
    v[0]=block(&a,0,1); v[1]=block(&a,1,2); apply(v,2,4);
    v[1]=block(&a,1,2); v[2]=block(&a,2,3); apply(v,3,1);
    v[0]=block(&a,2,3); v[1]=block(&a,1,2); v[2]=block(&a,3,3); apply(v,3,2);
    for (size_t i=0;i<3;++i) { double x=i+1; assert(*chem_sparse_matrix_at(&a,i,1,2)==2*x && *chem_sparse_matrix_at(&a,i,2,3)==3*x && *chem_sparse_matrix_at(&a,i,3,3)==6*x*x); }
    const ChemSparseEntry ds[]={{0,0},{1,1},{2,2}};
    for (size_t variant=0;variant<3;++variant) {
        assert(chem_sparse_matrix_init(&a,3,2,width,order,ds,3,variant==1 ? 5 : 0)==CHEM_MATRIX_OK);
        if (variant==0) for (size_t i=0;i<2;++i) { *chem_sparse_matrix_at(&a,i,0,0)=i; *chem_sparse_matrix_at(&a,i,1,1)=i+1; *chem_sparse_matrix_at(&a,i,2,2)=i+2; }
        if (variant==2) *chem_sparse_matrix_at(&a,1,0,0)=10;
        v[0]=block(&a,0,0); v[1]=block(&a,1,1); v[2]=block(&a,2,2);
        double temp[2]; v[3]=chem_vector_column(temp,2); apply(v,4,0);
        v[0]=v[3]; v[1]=block(&a,2,2); apply(v,2,4);
        assert(*chem_sparse_matrix_at(&a,0,2,2)==(variant==0 ? 6 : variant==1 ? 30 : 0));
        assert(*chem_sparse_matrix_at(&a,1,2,2)==(variant==0 ? 12 : variant==1 ? 30 : 20));
    }
    const ChemSparseEntry p1[]={{0,1},{1,1},{2,2}},p2[]={{0,0},{1,2},{2,3},{3,4}};
    assert(chem_sparse_matrix_init(&a,3,4,width,order,p1,3,0)==CHEM_MATRIX_OK &&
           chem_sparse_matrix_init(&b,5,4,width,order,p2,4,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) {
        *chem_sparse_matrix_at(&a,i,0,1)=2*i+1; *chem_sparse_matrix_at(&a,i,1,1)=3*i+2;
        *chem_sparse_matrix_at(&b,i,0,0)=i+10; *chem_sparse_matrix_at(&b,i,1,2)=i+20;
        *chem_sparse_matrix_at(&b,i,2,3)=i+30; *chem_sparse_matrix_at(&b,i,3,4)=i+40;
    }
    v[0]=block(&a,0,1);v[1]=block(&a,1,1);v[2]=block(&b,0,0);v[3]=block(&a,2,2); apply(v,4,0);
    for (size_t i=0;i<4;++i) {
        assert(*chem_sparse_matrix_at(&a,i,2,2)==13+6*i && *chem_sparse_matrix_at(&a,i,0,1)==2*i+1);
        assert(*chem_sparse_matrix_at(&b,i,0,0)==i+10 && *chem_sparse_matrix_at(&b,i,1,2)==i+20);
        assert(*chem_sparse_matrix_at(&b,i,2,3)==i+30 && *chem_sparse_matrix_at(&b,i,3,4)==i+40);
    }
    assert(chem_sparse_matrix_init(&a,4,3,width,order,four,4,0)==CHEM_MATRIX_OK && chem_matrix_init(&dense,3,4,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) {
        *chem_sparse_matrix_at(&a,i,0,1)=i+1; *chem_sparse_matrix_at(&a,i,1,2)=i+10;
        *chem_sparse_matrix_at(&a,i,2,3)=i+20;
        for (size_t j=0;j<4;++j) *chem_matrix_at(&dense,i,j)=i*100+j*10;
    }
    v[0]=block(&a,0,1);v[1]=block(&a,1,2);v[2]=column(&dense,0);v[3]=column(&dense,1);v[4]=block(&a,3,3); apply(v,5,6);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,3,3)==21+202*i && *chem_sparse_matrix_at(&a,i,2,3)==20+i);
    assert(*chem_matrix_at(&dense,0,0)==0 && *chem_matrix_at(&dense,1,2)==120 && *chem_matrix_at(&dense,2,3)==230);
    const ChemSparseEntry by_entry[]={{0,0},{0,2},{1,0},{1,1},{2,1},{2,2}};
    assert(chem_sparse_matrix_init(&a,3,3,width,order,by_entry,6,0)==CHEM_MATRIX_OK);
    double offsets[]={0,1,2}; ChemColumnView source=chem_vector_column(offsets,3);
    for (size_t i=0;i<6;++i) {
        size_t physical=chem_sparse_matrix_index(&a,0,by_entry[i].row,by_entry[i].column);
        assert(chem_sparse_matrix_block_by_index(&a,physical,&v[1])==CHEM_MATRIX_OK);
        v[0]=source; double increment=10*i;
        for (size_t j=0;j<3;++j) t1[j]=increment;
        v[2]=chem_vector_column(t1,3); /* output = vector + constant */
        ChemColumnView operation[]={v[0],v[2],v[1]}; apply(operation,3,1);
    }
    for (size_t i=0;i<6;++i) for (size_t j=0;j<3;++j)
        assert(*chem_sparse_matrix_at(&a,j,by_entry[i].row,by_entry[i].column)==10*i+j);
    const ChemSparseEntry copy_entries[]={{0,1},{2,0}},other_entry[]={{1,2}};
    assert(chem_sparse_matrix_init(&a,3,3,width,order,copy_entries,2,0)==CHEM_MATRIX_OK);
    ChemColumnView first=block(&a,0,1),second=block(&a,2,0);
    assert(chem_column_fill(&first,3.2)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,0,1)==3.2 && *chem_sparse_matrix_at(&a,i,2,0)==0);
    assert(chem_column_fill(&second,9.9)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,2,0)==9.9);
    double values[]={10,20,30},scratch[3]; ChemColumnView input=chem_vector_column(values,3),tmp=chem_vector_column(scratch,3);
    assert(chem_column_fill(&tmp,3.2)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(scratch[i]==3.2);
    assert(chem_column_fill(&tmp,9.9)==CHEM_MATRIX_OK && chem_column_copy(&second,&tmp)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,2,0)==9.9);
    assert(chem_column_fill(&second,0)==CHEM_MATRIX_OK);
    assert(chem_column_copy(&first,&input)==CHEM_MATRIX_OK && chem_column_copy(&tmp,&first)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(scratch[i]==values[i] && *chem_sparse_matrix_at(&a,i,2,0)==0);
    double expected[]={3.2,4.2,1.3}; input=chem_vector_column(expected,3);
    assert(chem_column_copy(&first,&input)==CHEM_MATRIX_OK && chem_column_copy(&second,&first)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,2,0)==expected[i]);
    assert(chem_column_copy(&tmp,&first)==CHEM_MATRIX_OK && chem_column_fill(&first,0)==CHEM_MATRIX_OK && chem_column_copy(&first,&tmp)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,0,1)==expected[i]);
    for (size_t i=0;i<3;++i) *chem_matrix_at(&dense,i,1)=100*(i+1);
    ChemColumnView dense_col=column(&dense,1); assert(chem_column_copy(&second,&dense_col)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&a,i,2,0)==100*(i+1));
    dense_col=column(&dense,0); assert(chem_column_copy(&tmp,&first)==CHEM_MATRIX_OK && chem_column_copy(&dense_col,&tmp)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&dense,i,0)==expected[i]);
    assert(chem_sparse_matrix_init(&b,3,3,width,order,other_entry,1,0)==CHEM_MATRIX_OK);
    ChemColumnView other=block(&b,1,2); assert(chem_column_copy(&other,&first)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_sparse_matrix_at(&b,i,1,2)==expected[i]);
    chem_sparse_matrix_destroy(&a); chem_sparse_matrix_destroy(&b); chem_matrix_destroy(&dense);
}
int main(void) { for (size_t width=0;width<=4;++width) for (int o=0;o<2;++o) run(width,(ChemSparseOrdering)o); return 0; }
