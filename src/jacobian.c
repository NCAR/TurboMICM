/* Diagnostic finite differences and comparisons follow the reference's
   perturbations, nonnegative forward-difference fallback, and error metrics. */
#include "chem_jacobian.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

ChemStatus chem_finite_difference(size_t nc, size_t n, ChemForcingFunction forcing,
    void *context, const double *y, double perturbation, double *jacobian)
{
    if (!nc || !n || !forcing || !y || !jacobian || !(perturbation > 0)
        || n > SIZE_MAX / n || nc > SIZE_MAX / n / n / sizeof(double)
        || nc > SIZE_MAX / n / 5 / sizeof(double)) return CHEM_INVALID;
    size_t length = n * nc;
    double *work = malloc(5 * length * sizeof(*work));
    if (!work) return CHEM_MEMORY;
    double *plus = work, *minus = plus + length, *fplus = minus + length;
    double *fminus = fplus + length, *fbase = fminus + length;
    for (size_t col = 0; col < n; ++col) {
        memcpy(plus, y, length * sizeof(*y));
        memcpy(minus, y, length * sizeof(*y));
        for (size_t c = 0; c < nc; ++c) {
            size_t index = col * nc + c;
            double h = perturbation * fmax(1, fabs(y[index]));
            plus[index] += h; minus[index] -= h;
        }
        memset(fplus, 0, length * sizeof(*fplus));
        memset(fminus, 0, length * sizeof(*fminus));
        forcing(nc, n, plus, fplus, context);
        forcing(nc, n, minus, fminus, context);
        for (size_t c = 0; c < nc; ++c) {
            double x = y[col * nc + c], h = perturbation * fmax(1, fabs(x));
            if (x - h >= 0) {
                for (size_t row = 0; row < n; ++row)
                    jacobian[(row * n + col) * nc + c] = (fplus[row * nc + c] - fminus[row * nc + c]) / (2 * h);
            } else {
                memset(fbase, 0, length * sizeof(*fbase));
                forcing(nc, n, y, fbase, context);
                for (size_t row = 0; row < n; ++row)
                    jacobian[(row * n + col) * nc + c] = (fplus[row * nc + c] - fbase[row * nc + c]) / h;
            }
        }
    }
    free(work);
    return CHEM_OK;
}
static int valid_comparison(const ChemSparse *a, const double *fd, ChemJacobianComparison *out)
{
    if (!a || !fd || !out) return 0;
    ChemView values;
    ChemArray *array = chem_sparse_values((ChemSparse *)a);
    if (chem_sparse_nonzeros(a) && (chem_array_column(array, 0, 1, &values) != CHEM_OK || values.type != CHEM_REAL)) return 0;
    *out = (ChemJacobianComparison){.passed = 1};
    return 1;
}
ChemStatus chem_compare_jacobian(const ChemSparse *a, const double *fd,
    double atol, double rtol, ChemJacobianComparison *out)
{
    if (!(atol >= 0) || !(rtol >= 0) || !valid_comparison(a, fd, out)) return CHEM_INVALID;
    size_t n = chem_sparse_order(a), nc = chem_array_rows(chem_sparse_values((ChemSparse *)a));
    for (size_t c = 0; c < nc; ++c)
        for (size_t row = 0; row < n; ++row)
            for (size_t col = 0; col < n; ++col) {
                if (chem_sparse_flat_index(a, row, col) == SIZE_MAX) continue;
                double analytical;
                chem_sparse_get(a, c, row, col, &analytical);
                analytical = -analytical;
                double value = fd[(row * n + col) * nc + c];
                double error = fabs(analytical - value), scale = fmax(fabs(analytical), fabs(value));
                double relative = scale > 0 ? error / scale : 0;
                if (error > out->max_absolute_error) {
                    out->max_absolute_error = error;
                    out->worst_block = c; out->worst_row = row; out->worst_column = col;
                    out->worst_analytical = analytical; out->worst_finite_difference = value;
                }
                if (relative > out->max_relative_error) out->max_relative_error = relative;
                if (error > atol + rtol * scale) out->passed = 0;
            }
    return CHEM_OK;
}
ChemStatus chem_check_jacobian_sparsity(const ChemSparse *a, const double *fd,
    double threshold, ChemJacobianComparison *out)
{
    if (!(threshold >= 0) || !valid_comparison(a, fd, out)) return CHEM_INVALID;
    size_t n = chem_sparse_order(a), nc = chem_array_rows(chem_sparse_values((ChemSparse *)a));
    for (size_t c = 0; c < nc; ++c)
        for (size_t row = 0; row < n; ++row)
            for (size_t col = 0; col < n; ++col) {
                if (chem_sparse_flat_index(a, row, col) != SIZE_MAX) continue;
                double value = fabs(fd[(row * n + col) * nc + c]);
                if (value > threshold) {
                    out->passed = 0;
                    if (value > out->max_absolute_error) {
                        out->max_absolute_error = value;
                        out->worst_block = c; out->worst_row = row; out->worst_column = col;
                        out->worst_finite_difference = value;
                    }
                }
            }
    return CHEM_OK;
}
