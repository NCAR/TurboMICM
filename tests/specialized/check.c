#include "chem_specialized.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern const chem_specialized ts1_specialized, chapman_specialized, robertson_specialized;
extern const chem_specialized rates_specialized, user_rates_specialized, edges_specialized, empty_specialized;

static int close_values(const char *what, const double *got, const double *expected, size_t n)
{
    for (size_t i=0; i<n; ++i) {
        if (!isfinite(got[i]) || !isfinite(expected[i]) || fabs(got[i]-expected[i])>1e-30+1e-10*fabs(expected[i])) {
            fprintf(stderr,"%s[%zu]: %.17g != %.17g\n",what,i,got[i],expected[i]); return 0;
        }
    }
    return 1;
}
static int check(const chem_specialized *sp, const char *name, size_t cells, int strict)
{
    chem_model m={0}; chem_state a={0},b={0}; int ok=0;
    double *x=NULL,*z=NULL;
    if (chem_model_create(&m,sp->mechanism,CHEM_LU_MOZART)!=CHEM_OK || chem_state_create(&a,&m,cells)!=CHEM_OK) goto done;
    size_t n=m.n_species*cells, nz=m.lu.nnz*cells;
    for (size_t c=0; c<cells; ++c) {
        a.conditions[c]=(chem_conditions){275+1.1*c,81000+40*c,32+0.1*c};
        for (size_t i=0; i<m.n_species; ++i) {
            double initial=a.y[i*cells+c];
            a.y[i*cells+c]=(initial ? initial : (m.n_species>20 ? 1e-12 : .01))*(1+0.01*c);
        }
        for (size_t j=0; j<m.n_custom; ++j) {
            double value=a.custom[j*cells+c];
            if (strstr(m.custom_names[j],"radius")) value=1e-7;
            else if (strstr(m.custom_names[j],"number concentration")) value=1e8;
            else value=value ? value : 1e-8;
            a.custom[j*cells+c]=value*(1+0.02*c);
        }
    }
    if (strict) { a.rtol=1e-6; for (size_t i=0; i<m.n_species; ++i) a.atol[i]=1e-15; }
    if (chem_state_copy(&b,&a)!=CHEM_OK) goto done;
    sp->rates(&a); chem_calculate_rates(&b);
    if (!close_values("rates",a.rates,b.rates,m.n_reactions*cells)) goto done;
    memset(a.work,0,n*sizeof(double)); memset(b.work,0,n*sizeof(double));
    sp->forcing(&m,cells,a.rates,a.y,a.work); chem_forcing(&m,cells,b.rates,b.y,b.work);
    if (!close_values("forcing",a.work,b.work,n)) goto done;
    sp->jacobian(&m,cells,a.rates,a.y,a.jacobian); chem_jacobian(&m,cells,b.rates,b.y,b.jacobian);
    if (!close_values("jacobian",a.jacobian,b.jacobian,nz)) goto done;
    /* A diagonally dominant matrix tests every generated fill/solve operation,
     * without relying on integration to happen to excite the sparse entries. */
    for (size_t i=0; i<m.n_species; ++i) for (size_t c=0; c<cells; ++c) {
        double total=1;
        for (size_t e=m.lu.row_start[i]; e<m.lu.row_start[i+1]; ++e) {
            double value=(double)((e*13+c*7)%19+1)*.001;
            a.jacobian[e*cells+c]=b.jacobian[e*cells+c]=value;
            total+=fabs(value);
        }
        a.jacobian[m.lu.diagonal[i]*cells+c]=b.jacobian[m.lu.diagonal[i]*cells+c]=total;
    }
    sp->factor(cells,a.jacobian,a.scratch); chem_lu_factor(&m.lu,cells,b.jacobian,b.scratch);
    if (!close_values("factor",a.jacobian,b.jacobian,nz)) goto done;
    x=malloc(n*sizeof(double)); z=malloc(n*sizeof(double)); if (!x || !z) goto done;
    for (size_t i=0; i<n; ++i) x[i]=z[i]=1+(double)(i%13)*.1;
    sp->linear(cells,a.jacobian,x); chem_lu_solve(&m.lu,cells,b.jacobian,z);
    if (!close_values("linear",x,z,n)) goto done;
    chem_ros_parameters p=chem_ros_parameters_default(CHEM_ROS3);
    for (size_t step=0; step<6; ++step) {
        /* Change runtime parameters between solves, including night and dawn. */
        double factor=step==2 ? 0 : 0.3+(double)step*.1;
        for (size_t j=0; j<m.n_custom; ++j) if (!strstr(m.custom_names[j],"."))
            for (size_t c=0; c<cells; ++c) a.custom[j*cells+c]=b.custom[j*cells+c]=1e-5*factor*(1+.01*c);
        for (size_t c=0; c<cells; ++c) {
            a.conditions[c].temperature+=.3;
            b.conditions[c]=a.conditions[c];
        }
        sp->rates(&a); chem_calculate_rates(&b);
        double dt=m.n_species>20 ? 30 : .01;
        chem_result ra=sp->solve(&a,dt), rb=chem_ros_solve(&b,dt,&p,NULL,NULL);
        if (ra.status!=CHEM_OK || rb.status!=CHEM_OK) {
            fprintf(stderr,"solve status %d/%d\n",ra.status,rb.status); goto done;
        }
        if (!close_values("solution",a.y,b.y,n)) goto done;
        if (ra.stats.internal_steps!=rb.stats.internal_steps || ra.stats.rejected!=rb.stats.rejected) {
            fprintf(stderr,"Step-control statistics differ\n"); goto done;
        }
        chem_state_clamp(&a); chem_state_clamp(&b);
    }
    ok=1;
done:
    printf("%s cells=%zu strict=%d %s\n",name,cells,strict,ok ? "PASS" : "FAIL");
    free(x); free(z); chem_state_destroy(&a); chem_state_destroy(&b); chem_model_destroy(&m); return ok;
}
int main(void)
{
    const chem_specialized *solvers[]={&ts1_specialized,&chapman_specialized,&robertson_specialized,&rates_specialized,&user_rates_specialized,&edges_specialized,&empty_specialized};
    const char *names[]={"ts1","chapman","robertson","rates","user_rates","edges","empty"};
    const size_t sizes[]={1,3,7,32,33};
    int ok=1;
    for (size_t i=0; i<sizeof(solvers)/sizeof(*solvers); ++i)
        for (size_t j=0; j<sizeof(sizes)/sizeof(*sizes); ++j)
            for (int strict=0; strict<2; ++strict) ok=check(solvers[i],names[i],sizes[j],strict) && ok;
    return ok ? 0 : 1;
}
