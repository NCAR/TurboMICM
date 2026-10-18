#ifndef CHEM_ADAPTIVE_H
#define CHEM_ADAPTIVE_H
#include "chem.h"

/* New methods are opt-in. The inherited ChemSolver API is unchanged. */
typedef struct {
  int method;                    /* 2..6: Rosenbrock; 10: BDF */
  int controller, norm;           /* controller: 0 I, 1 PI, 2 H211b; norm: 0 RMS, 1 cell RMS, 2 max */
  int max_order;                  /* BDF maximum order, 1..5 */
  size_t jacobian_age, max_steps;
  double h_start, h_max, safety, growth, shrink, positivity;
  double reuse_ratio;            /* BDF maximum relative change of matrix shift */
  int reuse;                     /* Permit factored matrix reuse */
  double step_hold;               /* >=1: step-size deadband (1 disables it) */
  double controller_b, controller_k; /* H211b filter b and error exponent k */
} ChemAdaptive;
typedef struct {
  ChemStats base;
  uint64_t nonlinear_iterations, matrix_reuses, positivity_rejections, failed_iterations;
  uint64_t clipped_values;
  double largest_clip;
} ChemAdaptiveStats;
typedef struct ChemWorkspace ChemWorkspace;
ChemAdaptive chem_adaptive_parameters(int method);
ChemWorkspace *chem_workspace_create(const ChemState *state);
void chem_workspace_free(ChemWorkspace *workspace);
void chem_workspace_reset(ChemWorkspace *workspace);
ChemStatus chem_adaptive_solve(const ChemAdaptive *parameters, double dt,
                               ChemState *state, ChemWorkspace *workspace,
                               ChemAdaptiveStats *stats);
#endif
