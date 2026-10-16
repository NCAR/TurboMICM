/* Central/forward-difference checks ported from MICM/UCAR, Apache-2.0. */
#include "chem_verify.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

chem_status chem_finite_difference(size_t n, size_t cells, const double *y,
                                   chem_forcing_function function, void *context, double perturbation, double *result)
{
    if (!function || !y || !result || perturbation<=0 || !n || cells>SIZE_MAX/n/5/sizeof(double)) return CHEM_INVALID;
    size_t count=n*cells;
    double *work=calloc(count*5+1,sizeof(double));
    if (!work) return CHEM_MEMORY;
    double *plus=work,*minus=work+count,*fp=work+2*count,*fm=work+3*count,*fb=work+4*count;
    for (size_t col=0; col<n; ++col) {
        memcpy(plus,y,count*sizeof(double)); memcpy(minus,y,count*sizeof(double));
        for (size_t c=0; c<cells; ++c) {
            double x=y[col*cells+c], h=perturbation*fmax(1,fabs(x));
            plus[col*cells+c]=x+h; minus[col*cells+c]=x-h;
        }
        memset(fp,0,count*sizeof(double)); memset(fm,0,count*sizeof(double));
        function(context,plus,fp); function(context,minus,fm);
        for (size_t c=0; c<cells; ++c) {
            double x=y[col*cells+c], h=perturbation*fmax(1,fabs(x));
            int central=x-h>=0;
            if (!central) { memset(fb,0,count*sizeof(double)); function(context,y,fb); }
            for (size_t row=0; row<n; ++row) {
                result[(row*n+col)*cells+c]=central ? (fp[row*cells+c]-fm[row*cells+c])/(2*h)
                    : (fp[row*cells+c]-fb[row*cells+c])/h;
            }
        }
    }
    free(work); return CHEM_OK;
}

chem_jacobian_comparison chem_compare_jacobian(size_t n,size_t cells,const size_t *map,
    const double *a,const double *fd,double atol,double rtol)
{
    chem_jacobian_comparison result={.passed=1};
    for (size_t c=0; c<cells; ++c) for (size_t row=0; row<n; ++row) for (size_t col=0; col<n; ++col) {
        size_t flat=map[row*n+col];
        if (flat==CHEM_NONE) continue;
        double analytical=-a[flat*cells+c], numerical=fd[(row*n+col)*cells+c];
        double error=fabs(analytical-numerical), scale=fmax(fabs(analytical),fabs(numerical));
        double relative=scale>0 ? error/scale : 0;
        if (error>result.max_absolute_error) {
            result.max_absolute_error=error; result.worst_cell=c; result.worst_row=row; result.worst_column=col;
            result.worst_analytical=analytical; result.worst_finite_difference=numerical;
        }
        if (relative>result.max_relative_error) result.max_relative_error=relative;
        if (!isfinite(error) || error>atol+rtol*scale) result.passed=0;
    }
    return result;
}

chem_jacobian_comparison chem_check_sparsity(size_t n,size_t cells,const size_t *map,const double *fd,double threshold)
{
    chem_jacobian_comparison result={.passed=1};
    for (size_t c=0; c<cells; ++c) for (size_t row=0; row<n; ++row) for (size_t col=0; col<n; ++col) {
        if (map[row*n+col]!=CHEM_NONE) continue;
        double value=fabs(fd[(row*n+col)*cells+c]);
        if (!isfinite(value) || value>threshold) {
            result.passed=0;
            if (value>result.max_absolute_error) {
                result.max_absolute_error=value; result.worst_cell=c; result.worst_row=row; result.worst_column=col;
                result.worst_finite_difference=value;
            }
        }
    }
    return result;
}
