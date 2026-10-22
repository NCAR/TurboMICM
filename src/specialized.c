#include "chem_specialized.h"

/* Cells is the full row stride; lanes is the current independent-cell tile.
   Only the helper execution order changes, never the order within a cell. */
void ck_run_factor_kernels(size_t cells, double *matrix, size_t count, const ck_factor_kernel *kernels,
                           size_t tile)
{
  if (!tile || tile > cells)
    tile = cells;
  for (size_t first = 0; first < cells;)
  {
    size_t lanes = cells - first < tile ? cells - first : tile;
    for (size_t i = 0; i < count; ++i)
      kernels[i](cells, lanes, matrix + first);
    first += lanes;
  }
}
void ck_run_solve_kernels(size_t cells, const double *matrix, double *rhs, size_t count,
                          const ck_solve_kernel *kernels, size_t tile)
{
  if (!tile || tile > cells)
    tile = cells;
  for (size_t first = 0; first < cells;)
  {
    size_t lanes = cells - first < tile ? cells - first : tile;
    for (size_t i = 0; i < count; ++i)
      kernels[i](cells, lanes, matrix + first, rhs + first);
    first += lanes;
  }
}
void ck_run_chemistry_kernels(size_t cells, const double *y, const double *rates, double *out, size_t count,
                              const ck_chemistry_kernel *kernels, size_t tile)
{
  if (!tile || tile > cells)
    tile = cells;
  for (size_t first = 0; first < cells;)
  {
    size_t lanes = cells - first < tile ? cells - first : tile;
    for (size_t i = 0; i < count; ++i)
      kernels[i](cells, lanes, y + first, rates + first, out + first);
    first += lanes;
  }
}

void ck_run_joint_kernels(size_t cells, const double *y, const double *rates, double *forcing,
                          double *jacobian, size_t count, const ck_joint_kernel *kernels, size_t tile)
{
  if (!tile || tile > cells)
    tile = cells;
  for (size_t first = 0; first < cells;)
  {
    size_t lanes = cells - first < tile ? cells - first : tile;
    for (size_t i = 0; i < count; ++i)
      kernels[i](cells, lanes, y + first, rates + first, forcing + first, jacobian + first);
    first += lanes;
  }
}
