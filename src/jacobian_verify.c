#include "jacobian_verify.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int chem_finite_difference_jacobian(const ChemMatrix *variables,double perturbation,ChemForcingFunction fn,void *context,ChemMatrix *fd)
{
    if(!variables || !fn || !fd || !(perturbation>0)) return 1;
    size_t n=variables->columns, blocks=variables->rows;
    if((n && n>SIZE_MAX/n) || n>SIZE_MAX/(4*sizeof(double)) || fd->rows!=blocks || fd->columns!=n*n) return 1;
    if(!n) return 0;
    double *buffer=n?malloc(4*n*sizeof(double)):NULL;
    if(n && !buffer) return 1;
    double *x=buffer,*plus=buffer+n,*minus=buffer+2*n,*base=buffer+3*n;
    for(size_t b=0;b<blocks;++b) {
        for(size_t i=0;i<n;++i) x[i]=*chem_matrix_const_at(variables,b,i);
        for(size_t col=0;col<n;++col) {
            double value=x[col], h=perturbation*fmax(1,fabs(value));
            memset(plus,0,n*sizeof(double)); memset(minus,0,n*sizeof(double));
            x[col]=value+h; fn(b,x,plus,context);
            x[col]=value-h; fn(b,x,minus,context);
            x[col]=value;
            if(value-h<0) { memset(base,0,n*sizeof(double)); fn(b,x,base,context); }
            for(size_t row=0;row<n;++row)
                *chem_matrix_at(fd,b,row*n+col)=value-h>=0?(plus[row]-minus[row])/(2*h):(plus[row]-base[row])/h;
        }
    }
    free(buffer); return 0;
}
static int verify(const ChemSparseMatrix *a,const ChemMatrix *fd,double atol,double rtol,int sparsity,ChemJacobianComparison *result)
{
    if(!a || !fd || !result || atol<0 || rtol<0 || fd->rows!=a->blocks || fd->columns!=a->order*a->order) return 1;
    *result=(ChemJacobianComparison){.passed=1}; size_t n=a->order;
    for(size_t b=0;b<a->blocks;++b) for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) {
        const double *entry=chem_sparse_const_at(a,b,i,j); double numeric=*chem_matrix_const_at(fd,b,i*n+j);
        if(sparsity) {
            if(!entry && fabs(numeric)>atol) {
                result->passed=0;
                if(fabs(numeric)>result->max_abs_error) { result->max_abs_error=fabs(numeric); result->worst_block=b; result->worst_row=i; result->worst_column=j; result->worst_fd=fabs(numeric); }
            }
        } else if(entry) {
            double analytical=-*entry, err=fabs(analytical-numeric), scale=fmax(fabs(analytical),fabs(numeric));
            if(err>result->max_abs_error) { result->max_abs_error=err; result->worst_block=b; result->worst_row=i; result->worst_column=j; result->worst_analytical=analytical; result->worst_fd=numeric; }
            if(scale>0 && err/scale>result->max_rel_error) result->max_rel_error=err/scale;
            if(!isfinite(err) || err>atol+rtol*scale) result->passed=0;
        }
    }
    return 0;
}
int chem_compare_jacobian(const ChemSparseMatrix *a,const ChemMatrix *fd,double atol,double rtol,ChemJacobianComparison *r)
{ return verify(a,fd,atol,rtol,0,r); }
int chem_check_jacobian_sparsity(const ChemSparseMatrix *a,const ChemMatrix *fd,double threshold,ChemJacobianComparison *r)
{ return verify(a,fd,threshold,0,1,r); }
