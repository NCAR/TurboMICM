#ifndef CHEM_SPECIALIZED_H
#define CHEM_SPECIALIZED_H
#include "chem.h"

/* Generated kernels use the default ordered model of this exact descriptor.
   Construct a model from mechanism and use its states with these operations. */
typedef struct
{
  const ck_mechanism *mechanism;
  ck_status (*update_rates)(ck_state *);
  void (*forcing)(const ck_state *, const double *, double *);
  void (*jacobian)(const ck_state *, const double *, double *);
  void (*factor)(size_t, double *);
  void (*linear_solve)(size_t, const double *, double *);
  ck_status (*solve)(ck_state *, double, ck_stats *);
} ck_specialized;

/* Separate translation-unit dispatch keeps C compilers from expanding all LU
   helpers into one enormous optimization unit. Arithmetic remains generated. */
typedef void (*ck_factor_kernel)(size_t, size_t, double *);
typedef void (*ck_solve_kernel)(size_t, size_t, const double *, double *);
typedef void (*ck_chemistry_kernel)(size_t, size_t, const double *, const double *, double *);
void ck_run_factor_kernels(size_t cells, double *matrix, size_t count, const ck_factor_kernel *kernels,
                           size_t tile);
void ck_run_solve_kernels(size_t cells, const double *matrix, double *rhs, size_t count,
                          const ck_solve_kernel *kernels, size_t tile);
void ck_run_chemistry_kernels(size_t cells, const double *y, const double *rates, double *out, size_t count,
                              const ck_chemistry_kernel *kernels, size_t tile);

typedef void (*ck_joint_kernel)(size_t, size_t, const double *, const double *, double *, double *);
void ck_run_joint_kernels(size_t cells, const double *y, const double *rates, double *forcing,
                          double *jacobian, size_t count, const ck_joint_kernel *kernels, size_t tile);

#endif
