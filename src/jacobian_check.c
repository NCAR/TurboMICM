#include "jacobian_check.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int chem_finite_difference_jacobian(size_t cells, size_t n, const double *states,
        ChemForcingCallback forcing, void *context, double perturbation, double *output)
{
    if (!n || !states || !forcing || !output || !(perturbation > 0.0) ||
        n > SIZE_MAX / (3 * sizeof(double)) || cells > SIZE_MAX / n ||
        cells * n > SIZE_MAX / n) return -1;
    double *work = malloc(3 * n * sizeof(*work));
    if (!work) return -1;
    double *state = work, *plus = work + n, *minus = work + 2 * n;
    for (size_t cell = 0; cell < cells; ++cell) {
        const double *base = states + cell * n;
        for (size_t col = 0; col < n; ++col) {
            double h = perturbation * fmax(1.0, fabs(base[col]));
            memcpy(state, base, n * sizeof(*state));
            state[col] = base[col] + h;
            memset(plus, 0, n * sizeof(*plus));
            forcing(context, cell, state, plus);
            int central = base[col] - h >= 0.0;
            state[col] = central ? base[col] - h : base[col];
            memset(minus, 0, n * sizeof(*minus));
            forcing(context, cell, state, minus);
            for (size_t row = 0; row < n; ++row)
                output[(cell * n + row) * n + col] = (plus[row] - minus[row]) / (central ? 2 * h : h);
        }
    }
    free(work);
    return 0;
}

ChemJacobianComparison chem_compare_jacobian(size_t cells, size_t n, const unsigned char *pattern,
        const double *analytical, const double *fd, double atol, double rtol)
{
    ChemJacobianComparison result = { .passed = 1 };
    for (size_t cell = 0; cell < cells; ++cell)
        for (size_t row = 0; row < n; ++row)
            for (size_t col = 0; col < n; ++col) {
                if (!pattern[row * n + col]) continue;
                size_t i = (cell * n + row) * n + col;
                double a = -analytical[i], b = fd[i], error = fabs(a - b);
                double scale = fmax(fabs(a), fabs(b));
                double relative = scale > 0 ? error / scale : 0;
                if (error > result.max_abs_error) {
                    result.max_abs_error = error; result.worst_cell = cell;
                    result.worst_row = row; result.worst_col = col;
                    result.worst_analytical = a; result.worst_fd = b;
                }
                if (relative > result.max_rel_error) result.max_rel_error = relative;
                if (error > atol + rtol * scale) result.passed = 0;
            }
    return result;
}

ChemJacobianComparison chem_check_jacobian_sparsity(size_t cells, size_t n, const unsigned char *pattern,
        const double *fd, double threshold)
{
    ChemJacobianComparison result = { .passed = 1 };
    for (size_t cell = 0; cell < cells; ++cell)
        for (size_t row = 0; row < n; ++row)
            for (size_t col = 0; col < n; ++col) {
                if (pattern[row * n + col]) continue;
                double value = fabs(fd[(cell * n + row) * n + col]);
                if (value > threshold) {
                    result.passed = 0;
                    if (value > result.max_abs_error) {
                        result.max_abs_error = value; result.worst_fd = value;
                        result.worst_cell = cell; result.worst_row = row; result.worst_col = col;
                    }
                }
            }
    return result;
}
