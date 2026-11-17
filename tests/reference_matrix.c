#include "chem.h"
#include "chem_matrix.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define D(m,r,c) (*(double *)chem_matrix_at((m),(r),(c)))
#define CD(m,r,c) (*(const double *)chem_matrix_at_const((m),(r),(c)))
#define I(m,r,c) (*(int *)chem_matrix_at((m),(r),(c)))
static int failures;
static void report(int grouped, const char *name, int ok)
{
    printf("%s:%s.%s %s %s\n",grouped ? "vector_matrix" : "matrix",grouped ? "VectorMatrix" : "Matrix",name,
           ok ? "PASS" : "FAIL",ok ? "reference quantities matched" : "reference quantity differed");
    failures+=!ok;
}
static int initd(ChemMatrix *m, size_t r, size_t c, size_t g, double v)
{ return chem_matrix_init(m,r,c,g,sizeof(v),&v)==CHEM_OK; }
static int initi(ChemMatrix *m, size_t r, size_t c, size_t g, int v)
{ return chem_matrix_init(m,r,c,g,sizeof(v),&v)==CHEM_OK; }

static int small(size_t group)
{
    ChemMatrix m={0}; if (!initd(&m,3,5,group,0)) return 0;
    D(&m,1,3)=64.7; D(&m,0,0)=41.2; D(&m,2,4)=102.3;
    const ChemMatrix *cm=&m; const double *raw=(const double *)(const void *)m.data;
    int ok=CD(cm,1,3)==64.7 && CD(cm,0,0)==41.2 && CD(cm,2,4)==102.3;
    if (!group) ok &= m.slots==15 && raw[0]==41.2 && raw[14]==102.3 && raw[8]==64.7;
    if (group==2) ok &= m.slots==20 && raw[0]==41.2 && raw[18]==102.3 && raw[7]==64.7;
    if (group==4) ok &= m.slots==20 && raw[0]==41.2 && raw[18]==102.3 && raw[13]==64.7;
    chem_matrix_destroy(&m); return ok;
}
static int initialize(size_t group)
{
    ChemMatrix m={0}; if (!initd(&m,2,3,group,12.4)) return 0;
    const ChemMatrix *cm=&m;
    int ok=CD(cm,0,0)==12.4 && CD(cm,1,0)==12.4 && CD(cm,1,2)==12.4;
    chem_matrix_destroy(&m); return ok;
}
static int loop(size_t group)
{
    ChemMatrix m={0}; if (!initi(&m,3,4,group,0)) return 0;
    for (size_t i=0;i<3;++i) for (size_t j=0;j<4;++j) I(&m,i,j)=(int)(100*i+j);
    const ChemMatrix *cm=&m;
    int ok=*(const int *)chem_matrix_at_const(cm,0,0)==0 && *(const int *)chem_matrix_at_const(cm,1,2)==102 &&
           *(const int *)chem_matrix_at_const(cm,2,3)==203 && *(const int *)chem_matrix_at_const(cm,0,3)==3;
    chem_matrix_destroy(&m); return ok;
}
static int strides(size_t group)
{
    ChemMatrix m={0}; if (!initi(&m,3,4,group,0)) return 0;
    size_t row=chem_matrix_row_stride(&m), column=chem_matrix_column_stride(&m);
    int *data=(int *)(void *)m.data;
    for (size_t i=0;i<3;++i) for (size_t j=0;j<4;++j) data[i*row+j*column]=(int)(100*i+j);
    int ok=row==(group ? 1 : 4) && column==(group ? group : 1);
    for (size_t i=0;i<3;++i) for (size_t j=0;j<4;++j) ok &= I(&m,i,j)==(int)(100*i+j);
    chem_matrix_destroy(&m); return ok;
}
static int to_vector(size_t group)
{
    ChemMatrix m={0}; if (!initd(&m,2,3,group,0)) return 0;
    double expected[]={13.2,31.2,314.2}, actual[3];
    for (size_t i=0;i<3;++i) D(&m,1,i)=expected[i];
    const ChemMatrix *cm=&m;
    int ok=chem_matrix_get_row(cm,1,actual,3)==CHEM_OK;
    for (size_t i=0;i<3;++i) ok &= actual[i]==expected[i];
    chem_matrix_destroy(&m); return ok;
}
static int from_vector(size_t group)
{
    ChemMatrix m={0}, empty={0}, bad={0};
    const double row0[]={412.3,32.4,41.3},row1[]={5.33,-.3,31.2};
    const void *rows[]={row0,row1}; const size_t sizes[]={3,3};
    int ok=chem_matrix_from_rows(&empty,NULL,NULL,0,group,sizeof(double))==CHEM_OK && empty.rows==0;
    ok &= chem_matrix_from_rows(&m,rows,sizes,2,group,sizeof(double))==CHEM_OK;
    if (ok) for (size_t r=0;r<2;++r) for (size_t c=0;c<3;++c) ok &= CD(&m,r,c)==(r ? row1[c] : row0[c]);
    const int b0[]={3},b1[]={4,5},b2[]={5}; const void *bad_rows[]={b0,b1,b2}; const size_t bad_sizes[]={1,2,1};
    ok &= chem_matrix_from_rows(&bad,bad_rows,bad_sizes,3,group,sizeof(int))==CHEM_INVALID;
    chem_matrix_destroy(&m); chem_matrix_destroy(&empty); chem_matrix_destroy(&bad); return ok;
}
static int assignment(size_t group)
{
    ChemMatrix m={0}; if (!initd(&m,4,3,group,0)) return 0;
    const double a[]={12.3,15.1,24.3}, b[]={14.3,52.3,65.7,16.34}, shortrow[]={13.2,52.8};
    int ok=chem_matrix_set_row(&m,2,a,3)==CHEM_OK && D(&m,0,0)==0 && D(&m,3,0)==0;
    for (size_t i=0;i<3;++i) ok &= D(&m,2,i)==a[i];
    ok &= chem_matrix_set_row(&m,2,b,4)==CHEM_OK && chem_matrix_set_row(&m,2,shortrow,2)==CHEM_INVALID;
    for (size_t i=0;i<3;++i) ok &= D(&m,2,i)==b[i];
    chem_matrix_destroy(&m); return ok;
}
static int axpy(size_t group)
{
    ChemMatrix m={0}, x={0}; if (!initd(&m,4,3,group,100) || !initd(&x,4,3,group,200)) goto fail;
    double expected=0, actual=0;
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j) {
        D(&m,i,j)=i*10.3+j*100.5; D(&x,i,j)=i*1.7+j*10.2;
        expected+=D(&m,i,j)+1.39*D(&x,i,j);
    }
    int ok=chem_matrix_axpy(&m,1.39,&x)==CHEM_OK;
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j) actual+=D(&m,i,j);
    ok &= fabs(actual-expected)<=1e-5;
    chem_matrix_destroy(&m); chem_matrix_destroy(&x); return ok;
fail:
    chem_matrix_destroy(&m); chem_matrix_destroy(&x); return 0;
}
static int foreach_body(ChemMatrix *const *m, void *const *v, void *context)
{
    (void)v; int count=*(int *)context;
    for (size_t i=0;i<m[0]->rows;++i) for (size_t j=0;j<m[0]->columns;++j)
        D(m[0],i,j)+=CD(m[1],i,j)-(count==3 ? CD(m[2],i,j) : 0);
    return CHEM_OK;
}
static int foreach_case(size_t group)
{
    ChemMatrix m[3]={{0}}; ChemMatrix *args[]={&m[0],&m[1],&m[2]}; ChemMatrixFunction fn={0}; int ok=1;
    for (int i=0;i<3;++i) ok &= initd(&m[i],4,3,group,0);
    if (!ok) goto done;
    for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j) {
        D(&m[0],i,j)=i*10.3+j*100.5; D(&m[1],i,j)=i*1.7+j*10.2; D(&m[2],i,j)=i*19.5+j*32.2;
    }
    for (int count=2;count<=3 && ok;++count) {
        for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j) D(&m[0],i,j)=i*10.3+j*100.5;
        ok=chem_matrix_function_init(&fn,args,(size_t)count,0,foreach_body,&count)==CHEM_OK &&
           chem_matrix_function_call(&fn,args,NULL,NULL)==CHEM_OK;
        for (size_t i=0;i<4;++i) for (size_t j=0;j<3;++j)
            ok &= fabs(D(&m[0],i,j)-(i*10.3+j*100.5+i*1.7+j*10.2-(count==3 ? i*19.5+j*32.2 : 0)))<=1e-5;
        chem_matrix_function_destroy(&fn);
    }
done:
    chem_matrix_function_destroy(&fn); for (int i=0;i<3;++i) chem_matrix_destroy(&m[i]); return ok;
}
static int fill_or_bound(size_t group, int kind)
{
    ChemMatrix m={0}; if (!initd(&m,2,3,group,0)) return 0;
    double two=2,one=1; int ok=1;
    if (kind==0) ok=chem_matrix_fill(&m,&two)==CHEM_OK;
    else if (kind==1) ok=chem_matrix_bound(&m,2,1)==CHEM_OK;
    else ok=chem_matrix_bound(&m,2,0)==CHEM_OK;
    const double *data=(const double *)(const void *)m.data;
    for (size_t i=0;i<m.slots;++i) ok &= data[i]==(kind==2 ? 0 : 2);
    if (kind) {
        ok &= chem_matrix_fill(&m,&one)==CHEM_OK; D(&m,1,1)=3;
        ok &= chem_matrix_bound(&m,2,kind==1)==CHEM_OK;
        for (size_t i=0;i<2;++i) for (size_t j=0;j<3;++j)
            ok &= D(&m,i,j)==(kind==1 ? (i==1 && j==1 ? 3 : 2) : (i==1 && j==1 ? 2 : 1));
    }
    chem_matrix_destroy(&m); return ok;
}
static int printing(size_t group)
{
    ChemMatrix m={0}; if (!initd(&m,2,3,group,0)) return 0;
    char buffer[128]; D(&m,1,1)=3;
    int ok=chem_matrix_format(&m,buffer,sizeof buffer)==CHEM_OK && strcmp(buffer,"0,0,0\n0,3,0\n")==0;
    chem_matrix_destroy(&m); return ok;
}
static int iterate(void)
{
    ChemMatrix m={0}; if (!initi(&m,3,4,0,0)) return 0;
    int *row=chem_matrix_at(&m,1,0);
    for (size_t i=0;i<4;++i) row[i]=42+(int)i;
    int ok=I(&m,0,3)==0 && I(&m,1,0)==42 && I(&m,1,1)==43 && I(&m,1,3)==45 && I(&m,2,0)==0;
    chem_matrix_destroy(&m); return ok;
}

static int all_groups(int grouped, int (*fn)(size_t))
{
    if (!grouped) return fn(0);
    int ok=1; for (size_t g=1;g<=4;++g) ok &= fn(g); return ok;
}

int main(void)
{
    for (int g=0;g<2;++g) {
        report(g,g ? "SmallVectorMatrix" : "SmallMatrix",small(g ? 2 : 0));
        report(g,g ? "SmallConstVectorMatrix" : "SmallConstMatrix",small(g ? 4 : 0));
        report(g,g ? "InitializeVectorMatrix" : "InitializeMatrix",initialize(g ? 1 : 0));
        report(g,g ? "InitializeConstVectorMatrix" : "InitializeConstMatrix",initialize(g ? 2 : 0));
        report(g,g ? "LoopOverVectorMatrix" : "LoopOverMatrix",loop(g ? 2 : 0));
        report(g,g ? "LoopOverConstVectorMatrix" : "LoopOverConstMatrix",loop(g ? 1 : 0));
        report(g,"Strides",g ? strides(3) && strides(4) : strides(0));
        report(g,"ConversionToVector",to_vector(g ? 3 : 0));
        report(g,"ConstConversionToVector",to_vector(g ? 1 : 0));
        report(g,"ConversionFromVector",from_vector(g ? 2 : 0));
        report(g,"AssignmentFromVector",assignment(g ? 2 : 0));
        report(g,"Axpy",all_groups(g,axpy)); report(g,"ForEach",all_groups(g,foreach_case));
        for (int kind=0;kind<3;++kind) {
            int ok=1; if (g) for (size_t group=1;group<=4;++group) ok &= fill_or_bound(group,kind);
            else ok=fill_or_bound(0,kind);
            report(g,kind==0 ? "SetScaler" : kind==1 ? "Max" : "Min",ok);
        }
        report(g,"Print",all_groups(g,printing));
    }
    report(0,"IterateOverMatrix",iterate());
    return failures ? 1 : 0;
}
