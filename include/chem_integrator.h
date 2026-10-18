#ifndef CK_INTEGRATOR_H
#define CK_INTEGRATOR_H
#include "chem.h"

/* Experimental ODE integrators. The original ck_solve interface is unchanged.
   One context per state, never concurrently shared. Rates may change between
   advances: Jacobian/LU caches are invalidated at every advance boundary. */
typedef struct ck_integrator ck_integrator;
typedef struct
{
  int method;                 /* 2,3,4,43,6: Rosenbrock; 12,13: Euler extrapolation */
  int controller;             /* 0: I, 1: PI */
  int norm;                   /* 0: group RMS, 1: max of per-cell RMS, 2: max */
  size_t jacobian_steps;      /* Maximum accepted steps per Jacobian, >=1 */
  size_t max_steps;
  double start, safety, growth, shrink, deadband;
  int warm_start;             /* Carry proposed step across advance calls */
  int positivity;             /* Include negative candidate values in error norm */
} ck_integrator_options;
ck_integrator_options ck_integrator_defaults(void);
ck_status ck_integrator_create(ck_state *state, const ck_integrator_options *options,
                               ck_integrator **integrator);
void ck_integrator_free(ck_integrator *integrator);
ck_status ck_integrator_advance(ck_integrator *integrator, double dt, ck_stats *stats);
#endif
