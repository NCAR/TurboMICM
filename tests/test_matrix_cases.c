#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static void add(double *out, const double *a, const double *b, void *context)
{
    (void)context;
    *out += *a - (b ? *b : 0);
}

static int matrix_case(int kind, size_t width, int grouped)
{
    ChemMatrix m={0}, x={0}, z={0}; int ok=1;
    if (kind==0 || kind==1) {
        if (chem_matrix_init(&m,3,5,width,grouped,0)) return 0;
        *chem_matrix_at(&m,1,3)=64.7; *chem_matrix_at(&m,0,0)=41.2; *chem_matrix_at(&m,2,4)=102.3;
        if (*chem_matrix_const_at(&m,1,3)!=64.7 || *chem_matrix_const_at(&m,0,0)!=41.2 || *chem_matrix_const_at(&m,2,4)!=102.3) ok=0;
        size_t groups=(3+width-1)/width;
        if (m.storage_size!=groups*5*width || m.groups!=groups || m.data[0]!=41.2 ||
            m.data[(2/width)*5*width+4*width+2%width]!=102.3 || m.data[(1/width)*5*width+3*width+1%width]!=64.7) ok=0;
    } else if (kind==2 || kind==3) {
        if (chem_matrix_init(&m,2,3,width,grouped,12.4)) return 0;
        if (*chem_matrix_const_at(&m,0,0)!=12.4 || *chem_matrix_const_at(&m,1,0)!=12.4 || *chem_matrix_const_at(&m,1,2)!=12.4) ok=0;
    } else if (kind==4 || kind==5 || kind==6 || kind==7) {
        if (chem_matrix_init(&m,3,4,width,grouped,0)) return 0;
        if (kind==6) {
            for(size_t c=0;c<4;++c) *chem_matrix_at(&m,1,c)=42+c;
            if (*chem_matrix_at(&m,0,3)!=0 || *chem_matrix_at(&m,1,0)!=42 || *chem_matrix_at(&m,1,1)!=43 || *chem_matrix_at(&m,1,3)!=45 || *chem_matrix_at(&m,2,0)!=0) ok=0;
        } else {
            for(size_t r=0;r<3;++r) for(size_t c=0;c<4;++c) *chem_matrix_at(&m,r,c)=r*100+c;
            if (*chem_matrix_const_at(&m,0,0)!=0 || *chem_matrix_const_at(&m,1,2)!=102 || *chem_matrix_const_at(&m,2,3)!=203 || *chem_matrix_const_at(&m,0,3)!=3) ok=0;
            if (kind==7 && (chem_matrix_row_stride(&m)!=(grouped ? 1 : 4) || chem_matrix_column_stride(&m)!=(grouped ? width : 1))) ok=0;
        }
    } else if (kind==8 || kind==9) {
        if (chem_matrix_init(&m,2,3,width,grouped,0)) return 0;
        double row[]={13.2,31.2,314.2}, output[3];
        if (chem_matrix_set_row(&m,1,row,3) || chem_matrix_get_row(&m,1,output,3) || output[0]!=13.2 || output[1]!=31.2 || output[2]!=314.2) ok=0;
    } else if (kind==10) {
        if (chem_matrix_from_rows(&m,0,NULL,NULL,width,grouped) || m.rows!=0) return 0;
        chem_matrix_free(&m);
        double a[]={412.3,32.4,41.3}, b[]={5.33,-.3,31.2}; const double *rows[]={a,b}; size_t lengths[]={3,3};
        if (chem_matrix_from_rows(&m,2,rows,lengths,width,grouped) || m.rows!=2 || m.columns!=3) return 0;
        for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c) if (*chem_matrix_at(&m,r,c)!=rows[r][c]) ok=0;
        ChemMatrix bad={0}; size_t bad_lengths[]={1,2};
        if (!chem_matrix_from_rows(&bad,2,rows,bad_lengths,width,grouped)) ok=0;
        chem_matrix_free(&bad);
    } else if (kind==11) {
        if (chem_matrix_init(&m,4,3,width,grouped,0)) return 0;
        double other[]={12.3,15.1,24.3},big[]={14.3,52.3,65.7,16.34},small[]={13.2,52.8};
        if (chem_matrix_set_row(&m,2,other,3) || *chem_matrix_at(&m,0,0)!=0 || *chem_matrix_at(&m,2,0)!=12.3 || *chem_matrix_at(&m,2,1)!=15.1 || *chem_matrix_at(&m,2,2)!=24.3 || *chem_matrix_at(&m,3,0)!=0) ok=0;
        if (chem_matrix_set_row(&m,2,big,4) || *chem_matrix_at(&m,2,0)!=14.3 || *chem_matrix_at(&m,2,1)!=52.3 || *chem_matrix_at(&m,2,2)!=65.7 || !chem_matrix_set_row(&m,2,small,2)) ok=0;
    } else if (kind==12 || kind==13) {
        if (chem_matrix_init(&m,4,3,width,grouped,0) || chem_matrix_init(&x,4,3,width,grouped,0) || chem_matrix_init(&z,4,3,width,grouped,0)) return 0;
        double expected=0, actual=0;
        for(size_t r=0;r<4;++r) for(size_t c=0;c<3;++c) {
            double y=r*10.3+c*100.5, xx=r*1.7+c*10.2, zz=r*19.5+c*32.2;
            *chem_matrix_at(&m,r,c)=y; *chem_matrix_at(&x,r,c)=xx; *chem_matrix_at(&z,r,c)=zz; expected+=y+1.39*xx;
        }
        if (kind==12) {
            if (chem_matrix_axpy(&m,1.39,&x)) ok=0;
            for(size_t r=0;r<4;++r) for(size_t c=0;c<3;++c) actual+=*chem_matrix_at(&m,r,c);
            if (fabs(actual-expected)>1e-5) ok=0;
        } else {
            if (chem_matrix_foreach(&m,&x,NULL,add,NULL)) ok=0;
            for(size_t r=0;r<4;++r) for(size_t c=0;c<3;++c) {
                if (fabs(*chem_matrix_at(&m,r,c)-((r*10.3+c*100.5)+(r*1.7+c*10.2)))>1e-5) ok=0;
                *chem_matrix_at(&m,r,c)=r*10.3+c*100.5;
            }
            if (chem_matrix_foreach(&m,&x,&z,add,NULL)) ok=0;
            for(size_t r=0;r<4;++r) for(size_t c=0;c<3;++c)
                if (fabs(*chem_matrix_at(&m,r,c)-((r*10.3+c*100.5)+(r*1.7+c*10.2)-(r*19.5+c*32.2)))>1e-5) ok=0;
        }
    } else if (kind>=14 && kind<=16) {
        if (chem_matrix_init(&m,2,3,width,grouped,0)) return 0;
        if (kind==14) chem_matrix_fill(&m,2);
        if (kind==15) chem_matrix_max(&m,2);
        if (kind==16) chem_matrix_min(&m,2);
        for(size_t i=0;i<m.storage_size;++i) if (m.data[i]!=(kind==16 ? 0 : 2)) ok=0;
        if (kind!=14) {
            chem_matrix_fill(&m,1); *chem_matrix_at(&m,1,1)=3;
            if(kind==15) chem_matrix_max(&m,2); else chem_matrix_min(&m,2);
            for(size_t r=0;r<2;++r) for(size_t c=0;c<3;++c)
                if (*chem_matrix_at(&m,r,c)!=(kind==15 ? (r==1 && c==1 ? 3 : 2) : (r==1 && c==1 ? 2 : 1))) ok=0;
        }
    } else if (kind==17) {
        if (chem_matrix_init(&m,2,3,width,grouped,0)) return 0;
        *chem_matrix_at(&m,1,1)=3;
        FILE *f=tmpfile(); if (!f) { chem_matrix_free(&m); return 0; }
        if (chem_matrix_print(&m,f) || fseek(f,0,SEEK_SET)) ok=0;
        char buf[64]={0}; size_t read=fread(buf,1,sizeof(buf)-1,f); buf[read]=0;
        if (strcmp(buf,"0,0,0\n0,3,0\n")) ok=0;
        fclose(f);
    }
    chem_matrix_free(&m); chem_matrix_free(&x); chem_matrix_free(&z); return ok;
}

static int many(int kind, int grouped)
{
    if (!grouped) return matrix_case(kind,1,0);
    for(size_t width=1;width<=4;++width) if (!matrix_case(kind,width,1)) return 0;
    return 1;
}

int main(void)
{
    const char *standard[]={"SmallMatrix","SmallConstMatrix","InitializeMatrix","InitializeConstMatrix","LoopOverMatrix","LoopOverConstMatrix","IterateOverMatrix","Strides","ConversionToVector","ConstConversionToVector","ConversionFromVector","AssignmentFromVector","Axpy","ForEach","SetScaler","Max","Min","Print"};
    const char *grouped[]={"SmallVectorMatrix","SmallConstVectorMatrix","InitializeVectorMatrix","InitializeConstVectorMatrix","LoopOverVectorMatrix","LoopOverConstVectorMatrix",NULL,"Strides","ConversionToVector","ConstConversionToVector","ConversionFromVector","AssignmentFromVector","Axpy","ForEach","SetScaler","Max","Min","Print"};
    for(int i=0;i<18;++i) {
        printf("matrix:Matrix.%s %d\n",standard[i],many(i,0));
        if (grouped[i]) printf("vector_matrix:VectorMatrix.%s %d\n",grouped[i],many(i,1));
    }
    return 0;
}
