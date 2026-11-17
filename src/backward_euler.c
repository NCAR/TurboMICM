#include "backward_euler.h"
#include "linear.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

ChemBackwardEulerOptions chem_backward_euler_defaults(const double *atol)
{
    ChemBackwardEulerOptions p = {1e-40, 0, {0.5,0.5,0.5,0.5,0.1}, 11, 1e-6, atol, NULL, CHEM_LU_MOZART};
    return p;
}

int chem_backward_euler_converged(size_t cells, size_t variables, const double *residual,
    const double *state, double small, const double *atol, double rtol)
{
    for (size_t cell = 0; cell < cells; ++cell)
        for (size_t i = 0; i < variables; ++i) {
            size_t k = cell*variables+i;
            double r = fabs(residual[k]);
            if (!isfinite(r) || (r > small && r > atol[i] && r > rtol*fabs(state[k]))) return 0;
        }
    return 1;
}

ChemSolveStatus chem_backward_euler_solve(size_t n, double *y, double duration,
    ChemOdeFunction rhs, ChemOdeFunction jacobian, void *context,
    const ChemBackwardEulerOptions *o, ChemSolverStats *stats)
{
    ChemSolverStats local = {0};
    if (!stats) stats = &local;
    memset(stats, 0, sizeof(*stats));
    if (!n || !y || !rhs || !jacobian || !o || !o->absolute_tolerance ||
        !isfinite(duration) || duration <= 0 || !isfinite(o->h_start) || o->h_start < 0 ||
        !isfinite(o->small) || o->small < 0 || !isfinite(o->relative_tolerance) ||
        o->relative_tolerance < 0 || !o->max_iterations ||
        (o->sparse_plan && chem_sparse_plan_size(o->sparse_plan) != n) ||
        (o->lu_algorithm != CHEM_LU_MOZART && o->lu_algorithm != CHEM_LU_DOOLITTLE)) return CHEM_SOLVE_INVALID;
    for (size_t i = 0; i < n; ++i)
        if (!isfinite(o->absolute_tolerance[i]) || o->absolute_tolerance[i] < 0) return CHEM_SOLVE_INVALID;
    for (size_t i = 0; i < 5; ++i)
        if (!isfinite(o->time_step_reductions[i]) || o->time_step_reductions[i] <= 0 ||
            o->time_step_reductions[i] >= 1) return CHEM_SOLVE_INVALID;
    if (n > SIZE_MAX/n || n*n > SIZE_MAX/sizeof(double) ||
        n > (SIZE_MAX/sizeof(double)-n*n)/2) return CHEM_SOLVE_ALLOCATION;
    double *work = malloc((n*n+2*n)*sizeof(*work));
    if (!work) return CHEM_SOLVE_ALLOCATION;
    double *matrix = work, *yn = matrix+n*n, *residual = yn+n;
    memcpy(yn,y,n*sizeof(*y));
    double h = o->h_start == 0 ? duration : o->h_start, t = 0;
    size_t successes = 0, failures = 0;
    ChemSolveStatus status = CHEM_SOLVE_OK;
    while (t < duration) {
        int converged = 0;
        for (size_t iteration = 0; iteration < o->max_iterations; ++iteration) {
            ++stats->steps;
            rhs(context,y,residual); ++stats->function_calls;
            jacobian(context,y,matrix); ++stats->jacobian_updates;
            for (size_t i = 0; i < n*n; ++i) matrix[i] = -matrix[i];
            for (size_t i = 0; i < n; ++i) matrix[i*n+i] += 1/h;
            ++stats->decompositions;
            int factored;
            if (o->sparse_plan)
                factored = o->lu_algorithm == CHEM_LU_DOOLITTLE ? chem_sparse_factor_doolittle(o->sparse_plan,matrix) : chem_sparse_factor(o->sparse_plan,matrix);
            else
                factored = o->lu_algorithm == CHEM_LU_DOOLITTLE ? chem_dense_factor_doolittle(matrix,n) : chem_dense_factor(matrix,n);
            if (!factored) { status = CHEM_SOLVE_SINGULAR; goto done; }
            for (size_t i = 0; i < n; ++i) residual[i] -= (y[i]-yn[i])/h;
            if (o->sparse_plan) chem_sparse_solve(o->sparse_plan,matrix,residual);
            else chem_dense_solve(matrix,n,residual);
            ++stats->solves;
            for (size_t i = 0; i < n; ++i) {
                double next = y[i]+residual[i];
                y[i] = next < 0 ? 0 : next;
            }
            if (iteration > 0 && chem_backward_euler_converged(1,n,residual,y,o->small,o->absolute_tolerance,o->relative_tolerance)) {
                converged = 1; break;
            }
        }
        if (!converged) {
            ++stats->rejected; successes = 0;
            if (failures >= 5) {
                t += h; status = CHEM_SOLVE_UNCONVERGED;
                for (size_t i = 0; i < n; ++i) {
                    if (isnan(y[i])) { status = CHEM_SOLVE_NAN; break; }
                    if (isinf(y[i])) status = CHEM_SOLVE_INFINITY;
                }
                break;
            }
            memcpy(y,yn,n*sizeof(*y)); h *= o->time_step_reductions[failures++];
        } else {
            ++stats->accepted; t += h; memcpy(yn,y,n*sizeof(*yn));
            if (++successes >= 2) { successes = 0; h *= 2; }
        }
        if (t >= duration) break;
        h = fmin(h, duration-t);
        if (h <= 0 || t+h == t) { status = CHEM_SOLVE_STEP_TOO_SMALL; break; }
    }
done:
    stats->final_time = t;
    for (size_t i = 0; i < n; ++i) y[i] = y[i] > 0 ? y[i] : 0;
    free(work); return status;
}
