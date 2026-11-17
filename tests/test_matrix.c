#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <string.h>
#include "matrix.h"
static void add(double *y,const double *x,size_t n,void *ctx)
{ (void)ctx; assert(n==1 || n==2); *y+=x[0]; if (n==2) *y-=x[1]; }
static void run(size_t width)
{
    ChemMatrix a={0},b={0},c={0};
    assert(chem_matrix_init(&a,3,5,width,0)==CHEM_MATRIX_OK);
    *chem_matrix_at(&a,1,3)=64.7; *chem_matrix_at(&a,0,0)=41.2; *chem_matrix_at(&a,2,4)=102.3;
    assert(*chem_matrix_at(&a,1,3)==64.7 && *chem_matrix_at(&a,0,0)==41.2 && *chem_matrix_at(&a,2,4)==102.3);
    if (!width) assert(a.storage_size==15 && a.data[0]==41.2 && a.data[8]==64.7 && a.data[14]==102.3);
    if (width==2) assert(a.storage_size==20 && a.groups==2 && a.data[7]==64.7 && a.data[18]==102.3);
    if (width==4) assert(a.storage_size==20 && a.groups==1 && a.data[13]==64.7 && a.data[18]==102.3);
    assert(chem_matrix_clone(&b,&a)==CHEM_MATRIX_OK);
    assert(*chem_matrix_const_at(&b,1,3)==64.7 && *chem_matrix_const_at(&b,0,0)==41.2 &&
           *chem_matrix_const_at(&b,2,4)==102.3 && b.data!=a.data);
    assert(chem_matrix_init(&a,2,3,width,12.4)==CHEM_MATRIX_OK);
    assert(*chem_matrix_at(&a,0,0)==12.4 && *chem_matrix_at(&a,1,0)==12.4 && *chem_matrix_at(&a,1,2)==12.4);
    assert(chem_matrix_init(&a,3,4,width,0)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) for (size_t j=0;j<4;++j) *chem_matrix_at(&a,i,j)=i*100+j;
    assert(*chem_matrix_const_at(&a,0,0)==0 && *chem_matrix_const_at(&a,1,2)==102 &&
           *chem_matrix_const_at(&a,2,3)==203 && *chem_matrix_const_at(&a,0,3)==3);
    for (size_t i=0;i<3;++i) for (size_t j=0;j<4;++j) assert(a.data[chem_matrix_index(&a,i,j)]==i*100+j);
    if (width) assert(chem_matrix_index(&a,0,1)-chem_matrix_index(&a,0,0)==width);
    else assert(chem_matrix_index(&a,1,0)-chem_matrix_index(&a,0,0)==4);
    assert(chem_matrix_init(&a,2,3,width,0)==CHEM_MATRIX_OK);
    const double slice[]={13.2,31.2,314.2}; double sliced[3];
    assert(chem_matrix_set_row(&a,1,slice,3)==CHEM_MATRIX_OK);
    assert(chem_matrix_get_row(&a,1,sliced,3)==CHEM_MATRIX_OK);
    for (size_t j=0;j<3;++j) assert(sliced[j]==slice[j]);
    assert(chem_matrix_from_rows(&c,0,width,NULL,NULL)==CHEM_MATRIX_OK && c.rows==0);
    const double row0[]={412.3,32.4,41.3},row1[]={5.33,-.3,31.2};
    const double *rows[]={row0,row1}; size_t lengths[]={3,3};
    assert(chem_matrix_from_rows(&a,2,width,rows,lengths)==CHEM_MATRIX_OK);
    for (size_t i=0;i<2;++i) for (size_t j=0;j<3;++j) assert(*chem_matrix_at(&a,i,j)==rows[i][j]);
    lengths[1]=2; assert(chem_matrix_from_rows(&a,2,width,rows,lengths)==CHEM_MATRIX_SIZE);
    assert(a.rows==2 && *chem_matrix_at(&a,1,2)==31.2);
    assert(chem_matrix_init(&a,4,3,width,0)==CHEM_MATRIX_OK);
    const double one[]={12.3,15.1,24.3},big[]={14.3,52.3,65.7,16.34}; double output[4];
    assert(chem_matrix_set_row(&a,2,one,3)==CHEM_MATRIX_OK);
    assert(chem_matrix_get_row(&a,2,output,3)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(output[i]==one[i]);
    assert(chem_matrix_set_row(&a,2,big,4)==CHEM_MATRIX_OK);
    assert(chem_matrix_set_row(&a,2,one,2)==CHEM_MATRIX_SIZE);
    assert(*chem_matrix_at(&a,0,0)==0 && *chem_matrix_at(&a,2,0)==14.3 && *chem_matrix_at(&a,2,1)==52.3 &&
           *chem_matrix_at(&a,2,2)==65.7 && *chem_matrix_at(&a,3,0)==0);
    assert(chem_matrix_init(&b,4,3,width,0)==CHEM_MATRIX_OK && chem_matrix_init(&c,4,3,width,0)==CHEM_MATRIX_OK);
    double sum=0,observed=0;
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j) {
        double y=i*10.3+j*100.5,x=i*1.7+j*10.2,z=i*19.5+j*32.2;
        *chem_matrix_at(&a,i,j)=y; *chem_matrix_at(&b,i,j)=x; *chem_matrix_at(&c,i,j)=z;
        sum+=y+1.39*x;
    }
    assert(chem_matrix_axpy(&a,1.39,&b)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j) observed+=*chem_matrix_at(&a,i,j);
    assert(fabs(sum-observed)<1e-5);
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j) *chem_matrix_at(&a,i,j)=i*10.3+j*100.5;
    const ChemMatrix *inputs[]={&b,&c};
    assert(chem_matrix_apply(&a,inputs,1,add,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j)
        assert(fabs(*chem_matrix_at(&a,i,j)-((i*10.3+j*100.5)+(i*1.7+j*10.2)))<1e-5);
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j) *chem_matrix_at(&a,i,j)=i*10.3+j*100.5;
    assert(chem_matrix_apply(&a,inputs,2,add,NULL)==CHEM_MATRIX_OK);
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j)
        assert(fabs(*chem_matrix_at(&a,i,j)-((i*10.3+j*100.5)+(i*1.7+j*10.2)-(i*19.5+j*32.2)))<1e-5);
    assert(chem_matrix_init(&a,2,3,width,0)==CHEM_MATRIX_OK);
    chem_matrix_max(&a,2); for (size_t i=0;i<a.storage_size;++i) assert(a.data[i]==2);
    chem_matrix_fill(&a,0); chem_matrix_min(&a,2); for (size_t i=0;i<a.storage_size;++i) assert(a.data[i]==0);
    chem_matrix_fill(&a,2); for (size_t i=0;i<a.storage_size;++i) assert(a.data[i]==2);
    chem_matrix_fill(&a,1); *chem_matrix_at(&a,1,1)=3; chem_matrix_max(&a,2);
    for (size_t i=0;i<2;++i) for (size_t j=0;j<3;++j) assert(*chem_matrix_at(&a,i,j)==(i==1 && j==1 ? 3 : 2));
    chem_matrix_fill(&a,1); *chem_matrix_at(&a,1,1)=3; chem_matrix_min(&a,2);
    for (size_t i=0;i<2;++i) for (size_t j=0;j<3;++j) assert(*chem_matrix_at(&a,i,j)==(i==1 && j==1 ? 2 : 1));
    chem_matrix_fill(&a,0); *chem_matrix_at(&a,1,1)=3;
    FILE *file=fopen("build/test_matrix_print.tmp","w+"); assert(file && chem_matrix_print(&a,file)==CHEM_MATRIX_OK); rewind(file);
    char printed[32]={0}; assert(fread(printed,1,12,file)==12); assert(!strcmp(printed,"0,0,0\n0,3,0\n"));
    assert(fclose(file)==0 && remove("build/test_matrix_print.tmp")==0);
    assert(chem_matrix_init(&a,3,2,width,0)==CHEM_MATRIX_OK);
    double vector[]={3.2,4.2,1.3};
    assert(chem_matrix_set_column(&a,1,vector,3)==CHEM_MATRIX_OK);
    assert(chem_matrix_get_column(&a,1,output,3)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(output[i]==vector[i]);
    assert(chem_matrix_copy_column(&a,0,&a,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&a,i,0)==vector[i]);
    assert(chem_matrix_set_column(&a,1,vector,2)==CHEM_MATRIX_SIZE);
    chem_matrix_fill(&a,0);
    double fill[]={3.2,3.2,3.2},temp[]={9.9,9.9,9.9};
    assert(chem_matrix_set_column(&a,1,fill,3)==CHEM_MATRIX_OK);
    assert(chem_matrix_set_column(&a,0,temp,3)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&a,i,0)==9.9 && *chem_matrix_at(&a,i,1)==3.2);
    assert(chem_matrix_set_column(&a,1,vector,3)==CHEM_MATRIX_OK);
    assert(chem_matrix_copy_column(&a,0,&a,1)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&a,i,0)==vector[i] && *chem_matrix_at(&a,i,1)==vector[i]);
    assert(chem_matrix_get_column(&a,1,output,3)==CHEM_MATRIX_OK);
    chem_matrix_fill(&a,0);
    assert(chem_matrix_set_column(&a,1,output,3)==CHEM_MATRIX_OK);
    for (size_t i=0;i<3;++i) assert(*chem_matrix_at(&a,i,1)==vector[i]);
    chem_matrix_destroy(&a);chem_matrix_destroy(&b);chem_matrix_destroy(&c);
}
int main(void) { for (size_t width=0;width<=4;++width) run(width); return 0; }
