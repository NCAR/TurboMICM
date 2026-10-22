#include "chem_specialized.h"

void ck_run_factor_kernels(size_t cells, double *matrix, size_t count, const ck_factor_kernel *kernels)
{
  for (size_t i = 0; i < count; ++i)
    kernels[i](cells, matrix);
}
void ck_run_solve_kernels(size_t cells, const double *matrix, double *rhs,
                          size_t count, const ck_solve_kernel *kernels)
{
  for (size_t i = 0; i < count; ++i)
    kernels[i](cells, matrix, rhs);
}
