/* Copyright (C) 2023-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0
 * Ports test_matrix_policy.hpp, test_matrix.cpp, and test_vector_matrix.cpp. */
#include "test.h"
#include "chem_matrix.h"

static ChemMatrix matrix(size_t rows,size_t columns,size_t group,ChemMatrixType type,double value)
{
    ChemMatrix m={0};
    CHECK(chem_matrix_create(rows,columns,group,type,value,&m)==CHEM_OK);
    return m;
}
static void put(ChemMatrix *m,size_t r,size_t c,double v)
{ CHECK(chem_matrix_set(m,r,c,v)==CHEM_OK); }
static void small(size_t group,int constant)
{
    ChemMatrix m=matrix(3,5,group,CHEM_MATRIX_REAL,0),copy={0};
    put(&m,1,3,64.7);put(&m,0,0,41.2);put(&m,2,4,102.3);
    if (constant) CHECK(chem_matrix_copy(&m,&copy)==CHEM_OK);
    const ChemMatrix *a=constant?&copy:&m;
    CHECK(chem_matrix_get(a,1,3)==64.7);CHECK(chem_matrix_get(a,0,0)==41.2);
    CHECK(chem_matrix_get(a,2,4)==102.3);CHECK(a->size>=3);
    const double *data=a->data;
    CHECK(a->size==(group?20:15));CHECK(data[0]==41.2);
    if (!group) { CHECK(data[14]==102.3);CHECK(data[8]==64.7); }
    else {
        CHECK(chem_matrix_group_size(a)==group*5);
        CHECK(chem_matrix_group_count(a)==(constant?1:2));CHECK(a->group==group);
        CHECK(data[18]==102.3);CHECK(data[constant?13:7]==64.7);
    }
    chem_matrix_free(&copy);chem_matrix_free(&m);
}
static void initialize(size_t group,int constant)
{
    ChemMatrix m=matrix(2,3,group,CHEM_MATRIX_REAL,12.4);
    const ChemMatrix *a=&m;
    CHECK(chem_matrix_get(a,0,0)==12.4);CHECK(chem_matrix_get(a,1,0)==12.4);
    CHECK(chem_matrix_get(a,1,2)==12.4);
    if (!constant) CHECK(*(double *)chem_matrix_at(&m,1,2)==12.4);
    chem_matrix_free(&m);
}
static void loop(size_t group,int constant)
{
    ChemMatrix m=matrix(3,4,group,CHEM_MATRIX_INT,0),copy={0};
    for(size_t i=0;i<m.rows;++i)for(size_t j=0;j<m.columns;++j)
        *(int *)chem_matrix_at(&m,i,j)=(int)(i*100+j);
    if(constant)CHECK(chem_matrix_copy(&m,&copy)==CHEM_OK);
    const ChemMatrix *a=constant?&copy:&m;
    CHECK(*(const int *)chem_matrix_const_at(a,0,0)==0);
    CHECK(*(const int *)chem_matrix_const_at(a,1,2)==102);
    CHECK(*(const int *)chem_matrix_const_at(a,2,3)==203);
    CHECK(*(const int *)chem_matrix_const_at(a,0,3)==3);
    chem_matrix_free(&copy);chem_matrix_free(&m);
}
static void iterate(void)
{
    ChemMatrix m=matrix(3,4,0,CHEM_MATRIX_INT,0),copy={0};
    int *row=chem_matrix_at(&m,1,0);
    for(size_t j=0;j<m.columns;++j) row[j]=(int)(42+j);
    CHECK(chem_matrix_copy(&m,&copy)==CHEM_OK);
    for(size_t j=0;j<m.columns;++j) CHECK(row[j]==(int)(42+j));
    CHECK(chem_matrix_get(&m,0,3)==0);CHECK(chem_matrix_get(&m,1,0)==42);
    CHECK(chem_matrix_get(&m,1,1)==43);CHECK(chem_matrix_get(&m,1,3)==45);
    CHECK(chem_matrix_get(&m,2,0)==0);
    chem_matrix_free(&copy);chem_matrix_free(&m);
}
static void strides(size_t group)
{
    ChemMatrix m=matrix(3,4,group,CHEM_MATRIX_INT,0);
    size_t rs=chem_matrix_row_stride(&m),cs=chem_matrix_column_stride(&m);
    CHECK(rs==(group?1:4));CHECK(cs==(group?group:1));
    for(size_t i=0;i<3;++i)for(size_t j=0;j<4;++j)
        ((int *)m.data)[i*rs+j*cs]=(int)(i*100+j);
    for(size_t i=0;i<3;++i)for(size_t j=0;j<4;++j)CHECK(chem_matrix_get(&m,i,j)==i*100+j);
    chem_matrix_free(&m);
}
static void to_vector(size_t group,int constant)
{
    ChemMatrix m=matrix(2,3,group,CHEM_MATRIX_REAL,0),copy={0};
    put(&m,1,0,13.2);put(&m,1,1,31.2);put(&m,1,2,314.2);
    if(constant)CHECK(chem_matrix_copy(&m,&copy)==CHEM_OK);
    double row[3];
    CHECK(chem_matrix_copy_row(constant?&copy:&m,1,row,3)==CHEM_OK);
    CHECK(row[0]==13.2);CHECK(row[1]==31.2);CHECK(row[2]==314.2);
    chem_matrix_free(&copy);chem_matrix_free(&m);
}
static void from_vector(size_t group)
{
    ChemMatrix m={0},empty={0},bad={0};
    CHECK(chem_matrix_from_rows(0,NULL,NULL,group,CHEM_MATRIX_REAL,&empty)==CHEM_OK);
    CHECK(empty.rows==0);
    double a[]={412.3,32.4,41.3},b[]={5.33,-.3,31.2};
    const void *rows[]={a,b};size_t sizes[]={3,3};
    CHECK(chem_matrix_from_rows(2,sizes,rows,group,CHEM_MATRIX_REAL,&m)==CHEM_OK);
    CHECK(m.rows==2);CHECK(m.columns==3);
    double slice[3];
    for(size_t r=0;r<2;++r){
        CHECK(chem_matrix_copy_row(&m,r,slice,3)==CHEM_OK);
        for(size_t c=0;c<3;++c){CHECK(slice[c]==((const double *)rows[r])[c]);CHECK(chem_matrix_get(&m,r,c)==slice[c]);}
    }
    int ia[]={3},ib[]={4,5},ic[]={5};const void *invalid[]={ia,ib,ic};size_t lengths[]={1,2,1};
    CHECK(chem_matrix_from_rows(3,lengths,invalid,group,CHEM_MATRIX_INT,&bad)==CHEM_DIMENSION);
    chem_matrix_free(&empty);chem_matrix_free(&bad);chem_matrix_free(&m);
}
static void assign_vector(size_t group)
{
    ChemMatrix m=matrix(4,3,group,CHEM_MATRIX_REAL,0);
    const double a[]={12.3,15.1,24.3},b[]={14.3,52.3,65.7,16.34},c[]={13.2,52.8};
    CHECK(chem_matrix_assign_row(&m,2,a,3)==CHEM_OK);
    CHECK(chem_matrix_get(&m,0,0)==0);CHECK(chem_matrix_get(&m,3,0)==0);
    for(size_t j=0;j<3;++j)CHECK(chem_matrix_get(&m,2,j)==a[j]);
    CHECK(chem_matrix_assign_row(&m,2,b,4)==CHEM_OK);
    CHECK(chem_matrix_get(&m,0,0)==0);CHECK(chem_matrix_get(&m,3,0)==0);
    for(size_t j=0;j<3;++j)CHECK(chem_matrix_get(&m,2,j)==b[j]);
    CHECK(chem_matrix_assign_row(&m,2,c,2)==CHEM_DIMENSION);
    chem_matrix_free(&m);
}
static void axpy(size_t group)
{
    ChemMatrix y=matrix(4,3,group,CHEM_MATRIX_REAL,100),x=matrix(4,3,group,CHEM_MATRIX_REAL,200);
    double sum=0,result=0;
    for(size_t i=0;i<4;++i)for(size_t j=0;j<3;++j){
        double a=i*10.3+j*100.5,b=i*1.7+j*10.2;
        put(&y,i,j,a);put(&x,i,j,b);sum+=a+1.39*b;
    }
    CHECK(chem_matrix_axpy(&y,1.39,&x)==CHEM_OK);
    for(size_t i=0;i<4;++i)for(size_t j=0;j<3;++j)result+=chem_matrix_get(&y,i,j);
    NEAR(sum,result,1e-5);chem_matrix_free(&x);chem_matrix_free(&y);
}
static double add_sub(double a,const double *b,size_t count,void *ctx)
{ (void)ctx;return count==1?a+b[0]:a+b[0]-b[1]; }
static void foreach(size_t group)
{
    ChemMatrix a=matrix(4,3,group,CHEM_MATRIX_REAL,0),b=matrix(4,3,group,CHEM_MATRIX_REAL,0),c=matrix(4,3,group,CHEM_MATRIX_REAL,0);
    const ChemMatrix *inputs[]={&b,&c};
    for(size_t i=0;i<4;++i)for(size_t j=0;j<3;++j){put(&b,i,j,i*1.7+j*10.2);put(&c,i,j,i*19.5+j*32.2);}
    for(size_t n=1;n<=2;++n){
        for(size_t i=0;i<4;++i)for(size_t j=0;j<3;++j)put(&a,i,j,i*10.3+j*100.5);
        CHECK(chem_matrix_for_each(&a,n,inputs,add_sub,NULL)==CHEM_OK);
        for(size_t i=0;i<4;++i)for(size_t j=0;j<3;++j)
            NEAR(chem_matrix_get(&a,i,j),(i*10.3+j*100.5)+(i*1.7+j*10.2)-(n==2?(i*19.5+j*32.2):0),1e-5);
    }
    chem_matrix_free(&a);chem_matrix_free(&b);chem_matrix_free(&c);
}
static void scalar(size_t group,int mode)
{
    ChemMatrix m=matrix(2,3,group,CHEM_MATRIX_REAL,0);
    CHECK((mode==0?chem_matrix_fill(&m,2):mode==1?chem_matrix_max(&m,2):chem_matrix_min(&m,2))==CHEM_OK);
    for(size_t i=0;i<m.size;++i)CHECK(((double *)m.data)[i]==(mode==2?0:2));
    if(mode){
        CHECK(chem_matrix_fill(&m,1)==CHEM_OK);put(&m,1,1,3);
        CHECK((mode==1?chem_matrix_max(&m,2):chem_matrix_min(&m,2))==CHEM_OK);
        for(size_t r=0;r<2;++r)for(size_t c=0;c<3;++c)
            CHECK(chem_matrix_get(&m,r,c)==(r==1&&c==1?(mode==1?3:2):(mode==1?2:1)));
    }
    chem_matrix_free(&m);
}
static void print_matrix(size_t group)
{
    ChemMatrix m=matrix(2,3,group,CHEM_MATRIX_REAL,0);put(&m,1,1,3);
    FILE *f=fopen("build/test-matrix-print.txt","w+");CHECK(f!=NULL);
    if(f){CHECK(chem_matrix_print(&m,f)==CHEM_OK);rewind(f);char s[100]={0};size_t n=fread(s,1,sizeof(s)-1,f);CHECK(n==12);CHECK(!strcmp(s,"0,0,0\n0,3,0\n"));fclose(f);remove("build/test-matrix-print.txt");}
    chem_matrix_free(&m);
}
#define GROUPS(call) do { for(size_t g=1;g<=4;++g){call;} } while(0)
void test_matrix(void)
{
    RUN("matrix:Matrix.SmallMatrix",small(0,0));RUN("matrix:Matrix.SmallConstMatrix",small(0,1));
    RUN("matrix:Matrix.InitializeMatrix",initialize(0,0));RUN("matrix:Matrix.InitializeConstMatrix",initialize(0,1));
    RUN("matrix:Matrix.LoopOverMatrix",loop(0,0));RUN("matrix:Matrix.LoopOverConstMatrix",loop(0,1));
    RUN("matrix:Matrix.IterateOverMatrix",iterate());RUN("matrix:Matrix.Strides",strides(0));
    RUN("matrix:Matrix.ConversionToVector",to_vector(0,0));RUN("matrix:Matrix.ConstConversionToVector",to_vector(0,1));
    RUN("matrix:Matrix.ConversionFromVector",from_vector(0));RUN("matrix:Matrix.AssignmentFromVector",assign_vector(0));
    RUN("matrix:Matrix.Axpy",axpy(0));RUN("matrix:Matrix.ForEach",foreach(0));
    RUN("matrix:Matrix.SetScaler",scalar(0,0));RUN("matrix:Matrix.Max",scalar(0,1));RUN("matrix:Matrix.Min",scalar(0,2));
    RUN("matrix:Matrix.Print",print_matrix(0));
    RUN("vector_matrix:VectorMatrix.SmallVectorMatrix",small(2,0));RUN("vector_matrix:VectorMatrix.SmallConstVectorMatrix",small(4,1));
    RUN("vector_matrix:VectorMatrix.InitializeVectorMatrix",initialize(1,0));RUN("vector_matrix:VectorMatrix.InitializeConstVectorMatrix",initialize(2,1));
    RUN("vector_matrix:VectorMatrix.LoopOverVectorMatrix",loop(2,0));RUN("vector_matrix:VectorMatrix.LoopOverConstVectorMatrix",loop(1,1));
    RUN("vector_matrix:VectorMatrix.Strides",{strides(3);strides(4);});
    RUN("vector_matrix:VectorMatrix.ConversionToVector",to_vector(3,0));RUN("vector_matrix:VectorMatrix.ConstConversionToVector",to_vector(1,1));
    RUN("vector_matrix:VectorMatrix.ConversionFromVector",from_vector(2));RUN("vector_matrix:VectorMatrix.AssignmentFromVector",assign_vector(2));
    RUN("vector_matrix:VectorMatrix.Axpy",GROUPS(axpy(g)));RUN("vector_matrix:VectorMatrix.ForEach",GROUPS(foreach(g)));
    RUN("vector_matrix:VectorMatrix.SetScaler",GROUPS(scalar(g,0)));RUN("vector_matrix:VectorMatrix.Max",GROUPS(scalar(g,1)));
    RUN("vector_matrix:VectorMatrix.Min",GROUPS(scalar(g,2)));RUN("vector_matrix:VectorMatrix.Print",GROUPS(print_matrix(g)));
}
