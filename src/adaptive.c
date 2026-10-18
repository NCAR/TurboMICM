/* Autonomous adaptive ODE methods, implemented from the formulas documented
 * in DESIGN.md. The original MICM-compatible integration path is unchanged. */
#include "chem_adaptive.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double bound(double x, double lo, double hi)
{ return fmin(hi,fmax(lo,x)); }

chem_adaptive_parameters chem_adaptive_defaults(void)
{
    return (chem_adaptive_parameters){CHEM_ADAPT_ROS,CHEM_RODA3,
        100000,1,3,0.9,5.0,0.0,0.0,0.0,0,1,0,1,8,0.03,0};
}
chem_status chem_adaptive_create(chem_adaptive_workspace *w, const chem_state *s)
{
    if (!w || !s || !s->model) return CHEM_INVALID;
    *w=(chem_adaptive_workspace){0};
    size_t n=s->model->n_species*s->cells;
    if (n>SIZE_MAX/12/sizeof(double)) return CHEM_MEMORY;
    w->model=s->model; w->cells=s->cells;
    w->jacobian=calloc(s->model->lu.nnz*s->cells,sizeof(double));
    w->vectors=calloc(12*n,sizeof(double));
    if (!w->jacobian || !w->vectors) { chem_adaptive_destroy(w); return CHEM_MEMORY; }
    return CHEM_OK;
}
void chem_adaptive_destroy(chem_adaptive_workspace *w)
{
    if (!w) return;
    free(w->jacobian); free(w->vectors); *w=(chem_adaptive_workspace){0};
}

static void add(size_t n,double a,const double *restrict x,double *restrict y)
{
    #pragma omp simd
    for (size_t i=0;i<n;++i) y[i]+=a*x[i];
}
static void forcing(chem_state *s,const double *y,double *f,chem_stats *stats)
{
    memset(f,0,s->model->n_species*s->cells*sizeof(double));
    chem_state_forcing(s,y,f); ++stats->function_calls;
}
static chem_status factor(chem_state *s,chem_adaptive_workspace *w,double alpha,chem_stats *stats)
{
    memcpy(s->jacobian,w->jacobian,s->model->lu.nnz*s->cells*sizeof(double));
    chem_alpha_minus_jacobian(s,alpha); ++stats->decompositions;
    return chem_state_factor(s);
}
static chem_status solve(chem_state *s,double *rhs,chem_stats *stats)
{ ++stats->solves; return chem_state_solve_linear(s,rhs); }

/* mode 0: all-cell RMS; 1: maximum cell RMS; 2: component maximum.
 * Projection is included in the error budget before accepting a step. */
static double error_norm(const chem_state *s,const double *yn,const double *e,int mode)
{
    double total=0,largest=0;
    for (size_t c=0;c<s->cells;++c) {
        double sum=0;
        for (size_t i=0;i<s->model->n_species;++i) {
            size_t k=i*s->cells+c;
            if (!isfinite(yn[k]) || !isfinite(e[k])) return INFINITY;
            double scale=s->atol[i]+s->rtol*fmax(fabs(s->y[k]),fabs(yn[k]));
            double z=fmax(fabs(e[k]),fmax(0,-yn[k]))/scale;
            sum+=z*z;
            if (mode==2) largest=fmax(largest,z);
        }
        total+=sum;
        if (mode==1) largest=fmax(largest,sqrt(sum/(double)s->model->n_species));
    }
    return mode ? largest : sqrt(total/(double)(s->cells*s->model->n_species));
}

static chem_status ros_step(chem_state *s,double h,const chem_ros_parameters *p,
    chem_adaptive_workspace *w,int clip,chem_stats *stats)
{
    size_t n=s->model->n_species*s->cells;
    double *yn=w->vectors,*err=yn+n,*initial=err+n,*k[6];
    for (size_t j=0;j<6;++j) k[j]=initial+(j+1)*n;
    for (size_t stage=0;stage<p->stages;++stage) {
        size_t base=stage*(stage-1)/2;
        if (!stage) memcpy(k[0],initial,n*sizeof(double));
        else if (p->new_function[stage]) {
            memcpy(yn,s->y,n*sizeof(double));
            for (size_t j=0;j<stage;++j) add(n,p->a[base+j],k[j],yn);
            if (clip) for (size_t i=0;i<n;++i) if (yn[i]<0) yn[i]=0;
            forcing(s,yn,k[stage],stats);
        }
        if (stage+1<p->stages && !p->new_function[stage+1])
            memcpy(k[stage+1],k[stage],n*sizeof(double));
        for (size_t j=0;j<stage;++j) add(n,p->c[base+j]/h,k[j],k[stage]);
        chem_status rc=solve(s,k[stage],stats); if (rc!=CHEM_OK) return rc;
    }
    memcpy(yn,s->y,n*sizeof(double)); memset(err,0,n*sizeof(double));
    for (size_t j=0;j<p->stages;++j) {
        add(n,p->m[j],k[j],yn); add(n,p->e[j],k[j],err);
    }
    return CHEM_OK;
}

/* Richardson extrapolation of linearly implicit Euler, fixed Jacobian over
 * one macro step. Sequence 1,2,3,4; polynomial extrapolation in substep size.
 * Every row starts at exactly the same initial state. */
static chem_status extrap_step(chem_state *s,double h,size_t order,
    chem_adaptive_workspace *w,chem_stats *stats,const chem_adaptive_parameters *p,size_t *used)
{
    size_t n=s->model->n_species*s->cells;
    double *yn=w->vectors,*err=yn+n,*rhs=yn+2*n,*table=yn+3*n;
    double *low=yn+8*n;
    for (size_t row=0;row<order;++row) {
        if (row) memcpy(low,table,n*sizeof(double));
        double sub=h/(double)(row+1),*v=table+row*n;
        chem_status rc=factor(s,w,1/sub,stats); if (rc!=CHEM_OK) return rc;
        memcpy(v,s->y,n*sizeof(double));
        for (size_t j=0;j<=row;++j) {
            if (!j) memcpy(rhs,w->vectors+9*n,n*sizeof(double));
            else forcing(s,v,rhs,stats);
            rc=solve(s,rhs,stats); if (rc!=CHEM_OK) return rc;
            add(n,1,rhs,v);
        }
        /* In-place Neville tableau: table[0] becomes the highest order. */
        for (size_t j=row;j>0;--j) {
            double ratio=(double)(row+1)/(double)j-1;
            for (size_t k=0;k<n;++k)
                table[(j-1)*n+k]=table[j*n+k]+(table[j*n+k]-table[(j-1)*n+k])/ratio;
        }
        if (row) {
            memcpy(yn,table,n*sizeof(double));
            for (size_t i=0;i<n;++i) err[i]=yn[i]-low[i];
            if (p->adaptive_order && error_norm(s,yn,err,p->error_norm)<=1) {
                *used=row+1; ++w->extrap_orders[row-1]; return CHEM_OK;
            }
        }
    }
    *used=order; ++w->extrap_orders[order-2];
    return CHEM_OK;
}

/* Stiffly accurate SDIRK of order two or three. Solve each stage with
 * simplified Newton and the shared diagonal factorization. The lower-order
 * embedded formula is determined by its first two quadrature moments. */
static double dirk_table(size_t order,double a[3][3],double b[3],double e[3])
{
    memset(a,0,9*sizeof(double)); memset(b,0,3*sizeof(double)); memset(e,0,3*sizeof(double));
    if (order==2) {
        double g=1-1/sqrt(2.0);
        a[0][0]=a[1][1]=g; a[1][0]=1-g;
        b[0]=1-g; b[1]=g; e[0]=-g; e[1]=g;
        return g;
    }
    /* Root in (0.43,0.44) of 6*g^3-18*g^2+9*g-1=0. */
    double g=0.43586652150845899941601945119356;
    a[0][0]=a[1][1]=a[2][2]=g;
    a[1][0]=(1-g)/2;
    a[2][0]=-1.5*g*g+4*g-.25; a[2][1]=1.5*g*g-5*g+1.25;
    for (size_t j=0;j<3;++j) b[j]=a[2][j];
    double low1=(.5-g)/((1-g)/2);
    e[0]=b[0]-(1-low1); e[1]=b[1]-low1; e[2]=b[2];
    return g;
}

static chem_status dirk_step(chem_state *s,double h,const chem_adaptive_parameters *p,
    chem_adaptive_workspace *w,chem_stats *stats)
{
    size_t n=s->model->n_species*s->cells,stages=p->extrap_order;
    double a[3][3],b[3],e[3],g=dirk_table(stages,a,b,e),alpha=1/(h*g);
    double *yn=w->vectors,*err=yn+n,*f=yn+2*n,*k=yn+3*n;
    double *base=yn+6*n,*z=yn+7*n,*delta=yn+8*n;
    for (size_t stage=0;stage<stages;++stage) {
        memcpy(base,s->y,n*sizeof(double));
        for (size_t j=0;j<stage;++j) add(n,h*a[stage][j],k+j*n,base);
        memcpy(z,base,n*sizeof(double));
        double previous=INFINITY; int converged=0;
        for (size_t it=0;it<p->newton_iterations;++it) {
            forcing(s,z,f,stats);
            for (size_t i=0;i<n;++i) delta[i]=f[i]+alpha*(base[i]-z[i]);
            chem_status rc=solve(s,delta,stats); if (rc!=CHEM_OK) return rc;
            ++w->nonlinear_iterations;
            add(n,1,delta,z);
            /* Newton corrections are measured without projection: an internal
             * negative stage is not a convergence failure in itself. */
            double norm=0;
            for (size_t c=0;c<s->cells;++c) {
                double sum=0;
                for (size_t i=0;i<s->model->n_species;++i) {
                    size_t q=i*s->cells+c;
                    double v=delta[q]/(s->atol[i]+s->rtol*fmax(fabs(s->y[q]),fabs(z[q])));
                    if (!isfinite(v)) { sum=INFINITY; break; }
                    sum+=v*v;
                }
                norm=fmax(norm,sqrt(sum/(double)s->model->n_species));
            }
            if (norm<=p->newton_tolerance) { converged=1; break; }
            if (!isfinite(norm) || (it>=2 && norm>.9*previous)) break;
            previous=norm;
        }
        if (!converged) { ++w->nonlinear_rejections; return CHEM_UNCONVERGED; }
        forcing(s,z,k+stage*n,stats);
    }
    memcpy(yn,z,n*sizeof(double)); memset(err,0,n*sizeof(double));
    for (size_t j=0;j<stages;++j) add(n,h*e[j],k+j*n,err);
    return CHEM_OK;
}

chem_result chem_adaptive_solve(chem_state *s,double dt,
    const chem_adaptive_parameters *p,chem_adaptive_workspace *w)
{
    chem_result r={CHEM_OK,{0}};
    if (!s || !p || !w || !s->model || !s->cells || !w->jacobian || !w->vectors ||
        w->model!=s->model || w->cells!=s->cells ||
        !isfinite(dt) || dt<0 || s->model->mechanism->n_constraints ||
        p->method<CHEM_ADAPT_ROS || p->method>CHEM_ADAPT_SDIRK || p->ros<CHEM_ROS2 || p->ros>CHEM_RODA4 ||
        !p->max_steps || !p->jacobian_age || p->extrap_order<2 || p->extrap_order>4 ||
        !(p->safety>0 && p->safety<1) || !(p->growth>=1 && isfinite(p->growth)) ||
        p->h_start<0 || !isfinite(p->h_start) || p->h_max<0 || !isfinite(p->h_max) ||
        !(p->reuse_window>=0 && p->reuse_window<1) || p->error_norm<0 || p->error_norm>2 ||
        p->controller<0 || p->controller>1 || !(s->rtol>0 && isfinite(s->rtol)) ||
        (p->method==CHEM_ADAPT_SDIRK && p->extrap_order>3) ||
        !p->newton_iterations || !(p->newton_tolerance>0 && p->newton_tolerance<1)) {
        r.status=CHEM_INVALID; return r;
    }
    for (size_t i=0;i<s->model->n_species;++i) if (!(s->atol[i]>0 && isfinite(s->atol[i]))) {
        r.status=CHEM_INVALID; return r;
    }
    if (!dt) return r;
    size_t n=s->model->n_species*s->cells,age=p->jacobian_age;
    chem_ros_parameters ros=chem_ros_parameters_default(p->ros);
    double order=p->method==CHEM_ADAPT_ROS ? ros.order : (double)p->extrap_order;
    double hmax=p->h_max ? fmin(dt,p->h_max) : dt;
    double h=p->warm_start && w->next_h>0 ? w->next_h : (p->h_start ? p->h_start : dt*1e-4);
    h=fmin(h,hmax);
    double t=0,last_error=1,factor_h=0;
    int initial_valid=0,rejected=0;
    while (t<dt) {
        if (r.stats.internal_steps>=p->max_steps) { r.status=CHEM_MAX_STEPS; break; }
        h=fmin(h,dt-t);
        if (!(h>0) || t+h==t) { r.status=CHEM_SMALL_STEP; break; }
        int fresh=age==0 || age>=p->jacobian_age;
        if (age>=p->jacobian_age) {
            memset(w->jacobian,0,s->model->lu.nnz*s->cells*sizeof(double));
            chem_state_jacobian(s,s->y,w->jacobian); ++r.stats.jacobian_updates;
            age=0; factor_h=0;
        } else ++w->jacobian_reuses;
        if (p->method!=CHEM_ADAPT_EXTRAP) {
            double gamma=ros.gamma[0];
            if (p->method==CHEM_ADAPT_ROS) {
                if (!initial_valid) forcing(s,s->y,w->vectors+2*n,&r.stats);
                initial_valid=1;
            } else {
                double a[3][3],b[3],e[3]; gamma=dirk_table(p->extrap_order,a,b,e);
            }
            if (h!=factor_h) {
                r.status=factor(s,w,1/(h*gamma),&r.stats); factor_h=h;
            } else ++w->factor_reuses;
            if (r.status==CHEM_OK) r.status=p->method==CHEM_ADAPT_ROS ?
                ros_step(s,h,&ros,w,p->clip_stages,&r.stats) : dirk_step(s,h,p,w,&r.stats);
        } else {
            if (!initial_valid) forcing(s,s->y,w->vectors+9*n,&r.stats);
            initial_valid=1;
            size_t used=p->extrap_order;
            r.status=extrap_step(s,h,p->extrap_order,w,&r.stats,p,&used); order=(double)used;
        }
        ++r.stats.internal_steps;
        if (r.status!=CHEM_OK && r.status!=CHEM_SINGULAR && r.status!=CHEM_UNCONVERGED &&
            r.status!=CHEM_NAN && r.status!=CHEM_INF) break;
        double norm=r.status==CHEM_OK ? error_norm(s,w->vectors,w->vectors+n,p->error_norm) : INFINITY;
        double proposal=p->safety*pow(fmax(norm,1e-12),-1/order);
        if (p->controller && !rejected && norm<=1)
            proposal=p->safety*pow(fmax(norm,1e-12),-0.7/order)*pow(last_error,0.4/order);
        if (norm<=1) {
            for (size_t i=0;i<n;++i) {
                double v=w->vectors[i];
                if (v<0) { ++w->clipped_values; w->clipped_mass-=v; v=0; }
                s->y[i]=v;
            }
            ++r.stats.accepted; ++age; t+=h; initial_valid=0;
            last_error=fmax(norm,1e-12);
            proposal=bound(proposal,0.2,rejected ? 1 : p->growth);
            if (fabs(proposal-1)<=p->reuse_window) proposal=1;
            h=fmin(hmax,h*proposal); rejected=0;
        } else {
            ++r.stats.rejected;
            if (!isfinite(norm) && r.status==CHEM_OK) ++w->nonfinite_rejections;
            else {
                for (size_t i=0;i<n;++i) if (w->vectors[i]<0) { ++w->positivity_rejections; break; }
            }
            /* Rejections from a stale Jacobian retry at the same step with a
             * fresh Jacobian before reducing the step. */
            if (!fresh) age=p->jacobian_age;
            else h*=isfinite(norm) ? bound(proposal,0.1,0.5) : 0.1;
            factor_h=0; rejected=1; r.status=CHEM_OK;
        }
    }
    w->next_h=h; r.stats.final_time=t; return r;
}
