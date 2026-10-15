/* Copyright (C) 2024-2026 University Corporation for Atmospheric Research
 * SPDX-License-Identifier: Apache-2.0 */
#include "chem_jacobian.h"
#include "internal.h"
#include <math.h>
#include <string.h>
ChemStatus chem_finite_difference(size_t n,size_t cells,ChemForcingFunction forcing,void *context,
                                 const double *values,double perturbation,double *jacobian)
{
    if (!forcing || !values || !jacobian || perturbation<=0) return CHEM_INVALID;
    if (n && cells>SIZE_MAX/n/5) return CHEM_MEMORY;
    size_t count=n*cells,bytes=count*sizeof(double);
    double *work=chem_alloc(count,5*sizeof(double));
    if (!work) return CHEM_MEMORY;
    double *plus=work,*minus=plus+count,*fp=minus+count,*fm=fp+count,*base=fm+count;
    for (size_t col=0;col<n;++col) {
        memcpy(plus,values,bytes);memcpy(minus,values,bytes);
        for (size_t c=0;c<cells;++c) {
            size_t i=col*cells+c;double h=perturbation*fmax(1,fabs(values[i]));
            plus[i]+=h;minus[i]-=h;
        }
        memset(fp,0,bytes);memset(fm,0,bytes);
        forcing(plus,fp,context);forcing(minus,fm,context);
        for (size_t c=0;c<cells;++c) {
            double x=values[col*cells+c],h=perturbation*fmax(1,fabs(x));
            if (x-h<0) { memset(base,0,bytes);forcing(values,base,context); }
            for (size_t row=0;row<n;++row)
                jacobian[(row*n+col)*cells+c]=(fp[row*cells+c]-(x-h>=0?fm:base)[row*cells+c])/(x-h>=0?2*h:h);
        }
    }
    free(work);return CHEM_OK;
}
ChemJacobianComparison chem_compare_jacobian(const ChemSparse *pattern,size_t cells,const double *analytic,
                                             const double *fd,double atol,double rtol)
{
    ChemJacobianComparison result={.passed=1};size_t n=chem_sparse_order(pattern);
    for (size_t c=0;c<cells;++c) for (size_t row=0;row<n;++row) for (size_t col=0;col<n;++col) {
        size_t i=chem_sparse_index(pattern,row,col);if (i==SIZE_MAX) continue;
        double a=-analytic[i*cells+c],f=fd[(row*n+col)*cells+c];
        double error=fabs(a-f),scale=fmax(fabs(a),fabs(f)),relative=scale>0?error/scale:0;
        if (error>result.max_absolute_error) {
            result.max_absolute_error=error;result.worst_cell=c;result.worst_row=row;result.worst_column=col;
            result.worst_analytical=a;result.worst_fd=f;
        }
        if (relative>result.max_relative_error) result.max_relative_error=relative;
        if (error>atol+rtol*scale) result.passed=0;
    }
    return result;
}
ChemJacobianComparison chem_check_sparsity(const ChemSparse *pattern,size_t cells,const double *fd,double threshold)
{
    ChemJacobianComparison result={.passed=1};size_t n=chem_sparse_order(pattern);
    for (size_t c=0;c<cells;++c) for (size_t row=0;row<n;++row) for (size_t col=0;col<n;++col) {
        if (chem_sparse_index(pattern,row,col)!=SIZE_MAX) continue;
        double value=fabs(fd[(row*n+col)*cells+c]);
        if (value>threshold) {
            result.passed=0;
            if (value>result.max_absolute_error) {
                result.max_absolute_error=value;result.worst_cell=c;result.worst_row=row;result.worst_column=col;
                result.worst_fd=value;
            }
        }
    }
    return result;
}
