/* Test-only scalar oracle translated from MICM/UCAR ChapmanODESolver,
 * Apache-2.0. This independent legacy implementation is never linked into
 * the chemistry library or benchmark. */
#ifndef CHAPMAN_REFERENCE_H
#define CHAPMAN_REFERENCE_H
#include <float.h>

static const size_t chapman_positions[23][2]={{0,0},{6,0},{7,0},{8,0},{1,1},{2,2},{3,3},{4,4},{5,4},{6,4},{5,5},{6,5},{6,6},{7,6},{8,6},{5,7},{6,7},{7,7},{8,7},{5,8},{6,8},{7,8},{8,8}};

static inline void chapman_reference_force(const double *rate_constants,const double *number_densities,double *force)
{
    memset(force,0,9*sizeof(double));
    force[5] = force[5] + rate_constants[5] * number_densities[8];
    force[5] = force[5] - rate_constants[0] * number_densities[4] * number_densities[5];
    force[5] = force[5] - rate_constants[1] * number_densities[5] * number_densities[7];
    force[6] = force[6] + 2 * rate_constants[4] * number_densities[7];
    force[6] = force[6] + rate_constants[6] * number_densities[8];
    force[6] = force[6] + rate_constants[0] * number_densities[4] * number_densities[5];
    force[6] = force[6] + rate_constants[1] * number_densities[5] * number_densities[7];
    force[6] = force[6] - rate_constants[2] * number_densities[6] * number_densities[8];
    force[6] = force[6] - rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
    force[7] = force[7] - rate_constants[4] * number_densities[7];
    force[7] = force[7] + rate_constants[5] * number_densities[8];
    force[7] = force[7] + rate_constants[6] * number_densities[8];
    force[7] = force[7] + 2 * rate_constants[2] * number_densities[6] * number_densities[8];
    force[7] = force[7] - rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
    force[8] = force[8] - rate_constants[5] * number_densities[8];
    force[8] = force[8] - rate_constants[6] * number_densities[8];
    force[8] = force[8] - rate_constants[2] * number_densities[6] * number_densities[8];
    force[8] = force[8] + rate_constants[3] * number_densities[0] * number_densities[6] * number_densities[7];
}

static inline void chapman_reference_jacobian(const double *rate_constants,const double *number_densities,double *jacobian)
{
    memset(jacobian,0,23*sizeof(double));
    jacobian[1] = jacobian[1] - rate_constants[3] * number_densities[6] * number_densities[7];
    jacobian[2] = jacobian[2] - rate_constants[3] * number_densities[6] * number_densities[7];
    jacobian[3] = jacobian[3] + rate_constants[3] * number_densities[6] * number_densities[7];
    jacobian[8] = jacobian[8] - rate_constants[0] * number_densities[5];
    jacobian[9] = jacobian[9] + rate_constants[0] * number_densities[5];
    jacobian[10] = jacobian[10] - rate_constants[0] * number_densities[4];
    jacobian[10] = jacobian[10] - rate_constants[1] * number_densities[7];
    jacobian[11] = jacobian[11] + rate_constants[0] * number_densities[4];
    jacobian[11] = jacobian[11] + rate_constants[1] * number_densities[7];
    jacobian[12] = jacobian[12] - rate_constants[2] * number_densities[8];
    jacobian[12] = jacobian[12] - rate_constants[3] * number_densities[0] * number_densities[7];
    jacobian[13] = jacobian[13] + 2 * rate_constants[2] * number_densities[8];
    jacobian[13] = jacobian[13] - rate_constants[3] * number_densities[0] * number_densities[7];
    jacobian[14] = jacobian[14] - rate_constants[2] * number_densities[8];
    jacobian[14] = jacobian[14] + rate_constants[3] * number_densities[0] * number_densities[7];
    jacobian[15] = jacobian[15] - rate_constants[1] * number_densities[5];
    jacobian[16] = jacobian[16] + 2 * rate_constants[4];
    jacobian[16] = jacobian[16] + rate_constants[1] * number_densities[5];
    jacobian[16] = jacobian[16] - rate_constants[3] * number_densities[0] * number_densities[6];
    jacobian[17] = jacobian[17] - rate_constants[4];
    jacobian[17] = jacobian[17] - rate_constants[3] * number_densities[0] * number_densities[6];
    jacobian[18] = jacobian[18] + rate_constants[3] * number_densities[0] * number_densities[6];
    jacobian[19] = jacobian[19] + rate_constants[5];
    jacobian[20] = jacobian[20] + rate_constants[6];
    jacobian[20] = jacobian[20] - rate_constants[2] * number_densities[6];
    jacobian[21] = jacobian[21] + rate_constants[5];
    jacobian[21] = jacobian[21] + rate_constants[6];
    jacobian[21] = jacobian[21] + 2 * rate_constants[2] * number_densities[6];
    jacobian[22] = jacobian[22] - rate_constants[5];
    jacobian[22] = jacobian[22] - rate_constants[6];
    jacobian[22] = jacobian[22] - rate_constants[2] * number_densities[6];
}

static inline void chapman_reference_factor(double *jacobian)
{
    jacobian[0] = 1. / jacobian[0];
    jacobian[1] = jacobian[1] * jacobian[0];
    jacobian[2] = jacobian[2] * jacobian[0];
    jacobian[3] = jacobian[3] * jacobian[0];
    jacobian[4] = 1. / jacobian[4];
    jacobian[5] = 1. / jacobian[5];
    jacobian[6] = 1. / jacobian[6];
    jacobian[7] = 1. / jacobian[7];
    jacobian[8] = jacobian[8] * jacobian[7];
    jacobian[9] = jacobian[9] * jacobian[7];
    jacobian[10] = 1. / jacobian[10];
    jacobian[11] = jacobian[11] * jacobian[10];
    jacobian[16] = jacobian[16] - jacobian[11] * jacobian[15];
    jacobian[20] = jacobian[20] - jacobian[11] * jacobian[19];
    jacobian[12] = 1. / jacobian[12];
    jacobian[13] = jacobian[13] * jacobian[12];
    jacobian[14] = jacobian[14] * jacobian[12];
    jacobian[17] = jacobian[17] - jacobian[13] * jacobian[16];
    jacobian[18] = jacobian[18] - jacobian[14] * jacobian[16];
    jacobian[21] = jacobian[21] - jacobian[13] * jacobian[20];
    jacobian[22] = jacobian[22] - jacobian[14] * jacobian[20];
    jacobian[17] = 1. / jacobian[17];
    jacobian[18] = jacobian[18] * jacobian[17];
    jacobian[22] = jacobian[22] - jacobian[18] * jacobian[21];
    jacobian[22] = 1. / jacobian[22];
}

static inline void chapman_reference_lower(const double *jacobian,const double *b,double *y)
{
    memset(y,0,9*sizeof(double));
    y[0] = b[0];
    y[1] = b[1];
    y[2] = b[2];
    y[3] = b[3];
    y[4] = b[4];
    y[5] = b[5];
    y[5] = y[5] - jacobian[8] * y[4];
    y[6] = b[6];
    y[6] = y[6] - jacobian[1] * y[0];
    y[6] = y[6] - jacobian[9] * y[4];
    y[6] = y[6] - jacobian[11] * y[5];
    y[7] = b[7];
    y[7] = y[7] - jacobian[2] * y[0];
    y[7] = y[7] - jacobian[13] * y[6];
    y[8] = b[8];
    y[8] = y[8] - jacobian[3] * y[0];
    y[8] = y[8] - jacobian[14] * y[6];
    y[8] = y[8] - jacobian[18] * y[7];
}

static inline void chapman_reference_upper(const double *jacobian,const double *y,double *x)
{
    memset(x,0,9*sizeof(double));
    double temporary;
    temporary = y[8];
    x[8] = jacobian[22] * temporary;
    temporary = y[7];
    temporary = temporary - jacobian[21] * x[8];
    x[7] = jacobian[17] * temporary;
    temporary = y[6];
    temporary = temporary - jacobian[16] * x[7];
    temporary = temporary - jacobian[20] * x[8];
    x[6] = jacobian[12] * temporary;
    temporary = y[5];
    temporary = temporary - jacobian[15] * x[7];
    temporary = temporary - jacobian[19] * x[8];
    x[5] = jacobian[10] * temporary;
    temporary = y[4];
    x[4] = jacobian[7] * temporary;
    temporary = y[3];
    x[3] = jacobian[6] * temporary;
    temporary = y[2];
    x[2] = jacobian[5] * temporary;
    temporary = y[1];
    x[1] = jacobian[4] * temporary;
    temporary = y[0];
    x[0] = jacobian[0] * temporary;
}

static inline void chapman_reference_linear(const double *jacobian,const double *b,double *x)
{ double y[9]; chapman_reference_lower(jacobian,b,y); chapman_reference_upper(jacobian,y,x); }

static inline void chapman_reference_rates(double T,const double *photo,double *k)
{
    k[0]=2.15e-11*exp(110/T); k[1]=3.3e-11*exp(55/T); k[2]=8e-12*exp(-2060/T);
    k[3]=6e-34*pow(T/300,2.4); memcpy(k+4,photo,3*sizeof(double));
}
/* Legacy single-cell ROS3 oracle: scalar arrays, reciprocal-diagonal Mozart
 * factors, absolute initial step 1e-6 and the original 1e-4 relative tolerance. */
static inline double chapman_reference_solve(double *y,const double *rates,double dt)
{
    const double a[]={1,1,0},c[]={-0.10156171083877702091975600115545e+01,0.40759956452537699824805835358067e+01,0.92076794298330791242156818474003e+01};
    const double m[]={1,0.61697947043828245592553615689730e+01,-0.42772256543218573326238373806514};
    const double e[]={.5,-0.29079558716805469821718236208017e+01,0.22354069897811569627360909276199};
    const double gamma=0.43586652150845899941601945119356;
    double h=fmin(1e-6,dt),time=0; size_t steps=0; int reject_last=0,reject_more=0;
    while (time-dt+DBL_EPSILON<=0) {
        if (steps>100 || time+.1*h==time || h<=DBL_EPSILON) break;
        h=fmin(h,fabs(dt-time)); double initial[9]; chapman_reference_force(rates,y,initial);
        int accepted=0;
        while (!accepted && steps<=100) {
            double jac[23],k[3][9],forcing[9],ynew[9],err[9]={0};
            chapman_reference_jacobian(rates,y,jac);
            for (size_t j=0;j<23;++j) jac[j]=-jac[j]+(chapman_positions[j][0]==chapman_positions[j][1] ? 1/(h*gamma) : 0);
            chapman_reference_factor(jac); memcpy(forcing,initial,sizeof(forcing));
            for (size_t stage=0;stage<3;++stage) {
                size_t base=stage ? stage*(stage-1)/2 : 0;
                if (stage==1) {
                    memcpy(ynew,y,sizeof(ynew));
                    for (size_t j=0;j<stage;++j) for (size_t i=0;i<9;++i) ynew[i]+=a[base+j]*k[j][i];
                    chapman_reference_force(rates,ynew,forcing);
                }
                double rhs[9]; memcpy(rhs,forcing,sizeof(rhs));
                for (size_t j=0;j<stage;++j) for (size_t i=0;i<9;++i) rhs[i]+=(c[base+j]/h)*k[j][i];
                chapman_reference_linear(jac,rhs,k[stage]);
            }
            memcpy(ynew,y,sizeof(ynew));
            for (size_t j=0;j<3;++j) for (size_t i=0;i<9;++i) { ynew[i]+=m[j]*k[j][i]; err[i]+=e[j]*k[j][i]; }
            double sum=0;
            for (size_t i=0;i<9;++i) sum+=pow(err[i]/(1e-3+1e-4*fmax(fabs(y[i]),fabs(ynew[i]))),2);
            double error=fmax(sqrt(sum/9),1e-10),next=h*fmin(6,fmax(.2,.9/pow(error,1.0/3))); ++steps;
            if (error<1) {
                time+=h; memcpy(y,ynew,sizeof(ynew)); next=fmax(0,fmin(next,dt));
                if (reject_last) next=fmin(next,h);
                reject_last=reject_more=0; h=next; accepted=1;
            } else {
                if (reject_more) next=h*.1;
                reject_more=reject_last; reject_last=1; h=next;
            }
        }
    }
    return time;
}
#endif
