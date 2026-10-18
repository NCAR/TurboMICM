#ifndef CHEM_ADAPTIVE_H
#define CHEM_ADAPTIVE_H
#include "chem.h"

/* New methods are opt-in. The inherited ChemSolver API is unchanged. */
typedef struct {
  int method;                    /* 2..6: Rosenbrock; 10: BDF */
  int controller, norm;           /* 0: I / group RMS; 1: PI / max cell RMS; 2: max */
  int max_order;                  /* BDF maximum order, 1..5 */
  size_t jacobian_age, max_steps;
  double h_start, h_max, safety, growth, shrink, positivity;
  double reuse_ratio;            /* BDF maximum relative change of matrix shift */
  int reuse;                     /* Permit factored matrix reuse */
} ChemAdaptive;
typedef struct {
  ChemStats base;
  uint64_t nonlinear_iterations, matrix_reuses, positivity_rejections, failed_iterations;
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
