#include "chem.h"
#include "chem_jacobian.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#define VALUE(m,r,c) (*(double *)chem_matrix_at((m),(r),(c)))
#define CONST_VALUE(m,r,c) (*(const double *)chem_matrix_at_const((m),(r),(c)))

int chem_finite_difference_jacobian(ChemForcingFunction forcing, void *context,
                                    const ChemMatrix *variables, double perturbation, ChemMatrix *result)
{
    if (!forcing || !variables || !variables->data || !result || result==variables || result->data ||
        variables->element_size!=sizeof(double) || !variables->columns ||
        variables->columns>SIZE_MAX/variables->columns || !isfinite(perturbation) || perturbation<=0) return CHEM_INVALID;
    size_t cells=variables->rows,n=variables->columns,g=variables->group_size;
    double zero=0; ChemMatrix plus={0},minus={0},fplus={0},fminus={0},fbase={0},output={0};
    int status=chem_matrix_init(&plus,cells,n,g,sizeof(double),&zero);
    if (!status) status=chem_matrix_init(&minus,cells,n,g,sizeof(double),&zero);
    if (!status) status=chem_matrix_init(&fplus,cells,n,g,sizeof(double),&zero);
    if (!status) status=chem_matrix_init(&fminus,cells,n,g,sizeof(double),&zero);
    if (!status) status=chem_matrix_init(&fbase,cells,n,g,sizeof(double),&zero);
    if (!status) status=chem_matrix_init(&output,cells,n*n,g,sizeof(double),&zero);
    if (status) goto done;
    status=forcing(variables,&fbase,context);
    for (size_t col=0;col<n && !status;++col) {
        status=chem_matrix_copy(&plus,variables);
        if (!status) status=chem_matrix_copy(&minus,variables);
        if (status) break;
        for (size_t b=0;b<cells;++b) {
            double x=CONST_VALUE(variables,b,col),h=perturbation*fmax(1,fabs(x));
            VALUE(&plus,b,col)=x+h; VALUE(&minus,b,col)=x-h;
        }
        status=chem_matrix_fill(&fplus,&zero);
        if (!status) status=chem_matrix_fill(&fminus,&zero);
        if (!status) status=forcing(&plus,&fplus,context);
        if (!status) status=forcing(&minus,&fminus,context);
        if (status) break;
        for (size_t b=0;b<cells;++b) {
            double x=CONST_VALUE(variables,b,col), h=perturbation*fmax(1,fabs(x));
            int central=x-h>=0;
            for (size_t row=0;row<n;++row) {
                double high=VALUE(&fplus,b,row), low=central ? VALUE(&fminus,b,row) : VALUE(&fbase,b,row);
                VALUE(&output,b,row*n+col)=(high-low)/(central ? 2*h : h);
            }
        }
    }
    if (!status) { *result=output; memset(&output,0,sizeof output); }
done:
    chem_matrix_destroy(&plus); chem_matrix_destroy(&minus); chem_matrix_destroy(&fplus);
    chem_matrix_destroy(&fminus); chem_matrix_destroy(&fbase); chem_matrix_destroy(&output);
    return status;
}
static int check_arguments(const ChemSparseMatrix *a, const ChemMatrix *fd, ChemJacobianComparison *result)
{
    if (!a || !fd || !result || !a->storage.data || !fd->data || a->storage.element_size!=sizeof(double) ||
        fd->element_size!=sizeof(double) || (a->order && a->order>SIZE_MAX/a->order) ||
        fd->columns!=a->order*a->order || fd->rows!=a->storage.rows) return CHEM_INVALID;
    memset(result,0,sizeof(*result)); result->passed=1;
    return CHEM_OK;
}
int chem_compare_jacobian(const ChemSparseMatrix *a, const ChemMatrix *fd,
                          double atol, double rtol, ChemJacobianComparison *result)
{
    int status=check_arguments(a,fd,result);
    if (status || !isfinite(atol) || !isfinite(rtol) || atol<0 || rtol<0) return CHEM_INVALID;
    for (size_t b=0;b<fd->rows;++b) for (size_t r=0;r<a->order;++r) for (size_t c=0;c<a->order;++c) {
        const double *entry=chem_sparse_at_const(a,b,r,c); if (!entry) continue;
        double actual=-*entry, expected=CONST_VALUE(fd,b,r*a->order+c), error=fabs(actual-expected);
        double scale=fmax(fabs(actual),fabs(expected)), relative=scale>0 ? error/scale : 0;
        if (error>result->max_absolute_error) {
            result->max_absolute_error=error; result->worst_block=b; result->worst_row=r; result->worst_column=c;
            result->worst_analytical=actual; result->worst_finite_difference=expected;
        }
        if (relative>result->max_relative_error) result->max_relative_error=relative;
        if (!isfinite(error) || error>atol+rtol*scale) result->passed=0;
    }
    return CHEM_OK;
}
int chem_check_jacobian_sparsity(const ChemSparseMatrix *a, const ChemMatrix *fd,
                                 double threshold, ChemJacobianComparison *result)
{
    int status=check_arguments(a,fd,result);
    if (status || !isfinite(threshold) || threshold<0) return CHEM_INVALID;
    for (size_t b=0;b<fd->rows;++b) for (size_t r=0;r<a->order;++r) for (size_t c=0;c<a->order;++c) {
        if (chem_sparse_entry(a,r,c)!=SIZE_MAX) continue;
        double value=fabs(CONST_VALUE(fd,b,r*a->order+c));
        if (!isfinite(value) || value>threshold) {
            result->passed=0;
            if (!isfinite(value) || value>result->max_absolute_error) {
                result->max_absolute_error=value; result->worst_block=b; result->worst_row=r; result->worst_column=c;
                result->worst_analytical=0; result->worst_finite_difference=value;
            }
        }
    }
    return CHEM_OK;
}
