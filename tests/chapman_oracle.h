/* SPDX-License-Identifier: Apache-2.0
 * Test-only translation of MICM's fixed Chapman regression oracle. This is
 * intentionally independent of the mechanism, sparse LU and integrator code. */
#ifndef CHEM_CHAPMAN_ORACLE_H
#define CHEM_CHAPMAN_ORACLE_H
#include <float.h>
static void chapman_force(const double *k,const double *y,double *f)
{
    for (size_t i=0;i<9;++i) f[i]=0;
    f[5]=k[5]*y[8]-k[0]*y[4]*y[5]-k[1]*y[5]*y[7];
    f[6]=2*k[4]*y[7]+k[6]*y[8]+k[0]*y[4]*y[5]+k[1]*y[5]*y[7]-k[2]*y[6]*y[8]-k[3]*y[0]*y[6]*y[7];
    f[7]=-k[4]*y[7]+k[5]*y[8]+k[6]*y[8]+2*k[2]*y[6]*y[8]-k[3]*y[0]*y[6]*y[7];
    f[8]=-k[5]*y[8]-k[6]*y[8]-k[2]*y[6]*y[8]+k[3]*y[0]*y[6]*y[7];
}
static void chapman_jacobian(const double *k,const double *y,double *j)
{
    for (size_t i=0;i<23;++i) j[i]=0;
    j[1]=-k[3]*y[6]*y[7];j[2]=-k[3]*y[6]*y[7];j[3]=k[3]*y[6]*y[7];
    j[8]=-k[0]*y[5];j[9]=k[0]*y[5];
    j[10]=-k[0]*y[4]-k[1]*y[7];j[11]=k[0]*y[4]+k[1]*y[7];
    j[12]=-k[2]*y[8]-k[3]*y[0]*y[7];j[13]=2*k[2]*y[8]-k[3]*y[0]*y[7];
    j[14]=-k[2]*y[8]+k[3]*y[0]*y[7];j[15]=-k[1]*y[5];
    j[16]=2*k[4]+k[1]*y[5]-k[3]*y[0]*y[6];j[17]=-k[4]-k[3]*y[0]*y[6];j[18]=k[3]*y[0]*y[6];
    j[19]=k[5];j[20]=k[6]-k[2]*y[6];j[21]=k[5]+k[6]+2*k[2]*y[6];j[22]=-k[5]-k[6]-k[2]*y[6];
}
static void chapman_factor(double *j,double alpha)
{
    for (size_t i=0;i<23;++i) j[i]=-j[i];
    const size_t diagonal[]={0,4,5,6,7,10,12,17,22};
    for (size_t i=0;i<9;++i) j[diagonal[i]]+=alpha;
    j[0]=1/j[0];j[1]*=j[0];j[2]*=j[0];j[3]*=j[0];j[4]=1/j[4];j[5]=1/j[5];j[6]=1/j[6];
    j[7]=1/j[7];j[8]*=j[7];j[9]*=j[7];j[10]=1/j[10];j[11]*=j[10];
    j[16]-=j[11]*j[15];j[20]-=j[11]*j[19];j[12]=1/j[12];j[13]*=j[12];j[14]*=j[12];
    j[17]-=j[13]*j[16];j[18]-=j[14]*j[16];j[21]-=j[13]*j[20];j[22]-=j[14]*j[20];
    j[17]=1/j[17];j[18]*=j[17];j[22]-=j[18]*j[21];j[22]=1/j[22];
}
static void chapman_linear_solve(const double *j,double *x)
{
    x[5]-=j[8]*x[4];x[6]-=j[1]*x[0];x[6]-=j[9]*x[4];x[6]-=j[11]*x[5];
    x[7]-=j[2]*x[0];x[7]-=j[13]*x[6];x[8]-=j[3]*x[0];x[8]-=j[14]*x[6];x[8]-=j[18]*x[7];
    x[8]*=j[22];x[7]=(x[7]-j[21]*x[8])*j[17];
    x[6]=(x[6]-j[16]*x[7]-j[20]*x[8])*j[12];x[5]=(x[5]-j[15]*x[7]-j[19]*x[8])*j[10];
    x[4]*=j[7];x[3]*=j[6];x[2]*=j[5];x[1]*=j[4];x[0]*=j[0];
}
static int chapman_solve(const double *rates,double *y,double end)
{
    const double c[]={-0.10156171083877702091975600115545e+01,0.40759956452537699824805835358067e+01,
        0.92076794298330791242156818474003e+01};
    const double m[]={1,0.61697947043828245592553615689730e+01,-0.42772256543218573326238373806514};
    const double e[]={.5,-0.29079558716805469821718236208017e+01,0.22354069897811569627360909276199};
    const double gamma=0.43586652150845899941601945119356;
    double time=0,h=1e-6;int reject_last=0,reject_more=0;size_t steps=0;
    while (time-end+DBL_EPSILON<=0) {
        if (steps>100 || time+.1*h==time || h<=DBL_EPSILON) return 0;
        h=fmin(h,fabs(end-time));double initial[9];chapman_force(rates,y,initial);int accepted=0;
        while (!accepted) {
            if (steps>100) return 0;
            double j[23],stages[3][9],f[9],next[9],error[9]={0};
            chapman_jacobian(rates,y,j);chapman_factor(j,1/(h*gamma));
            memcpy(stages[0],initial,sizeof(initial));chapman_linear_solve(j,stages[0]);
            for (size_t v=0;v<9;++v) next[v]=y[v]+stages[0][v];
            chapman_force(rates,next,f);
            for (size_t v=0;v<9;++v) stages[1][v]=f[v]+(c[0]/h)*stages[0][v];
            chapman_linear_solve(j,stages[1]);
            for (size_t v=0;v<9;++v) stages[2][v]=f[v]+(c[1]/h)*stages[0][v]+(c[2]/h)*stages[1][v];
            chapman_linear_solve(j,stages[2]);memcpy(next,y,sizeof(next));
            for (size_t stage=0;stage<3;++stage) for (size_t v=0;v<9;++v) {
                next[v]+=m[stage]*stages[stage][v];error[v]+=e[stage]*stages[stage][v];
            }
            double sum=0;
            for (size_t v=0;v<9;++v) sum+=pow(error[v]/(1e-3+1e-4*fmax(fabs(y[v]),fabs(next[v]))),2);
            double norm=fmax(sqrt(sum/9),1e-10),hn=h*fmin(6,fmax(.2,.9/pow(norm,1.0/3)));
            ++steps;
            if (norm<1) {
                time+=h;memcpy(y,next,sizeof(next));hn=fmin(hn,end);
                if (reject_last) hn=fmin(hn,h);
                reject_last=reject_more=0;h=hn;accepted=1;
            } else {
                if (reject_more) hn=h*.1;
                reject_more=reject_last;reject_last=1;h=hn;
            }
        }
    }
    return 1;
}
#endif
