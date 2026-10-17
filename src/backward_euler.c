/* Backward Euler integrator: placeholder until the port of MICM's
 * BackwardEuler::Solve lands. */
#include "internal.h"

tmc_result tmc_backward_euler_solve(tmc_solver *s, double time_step, tmc_state *st)
{
  (void)s;
  (void)time_step;
  (void)st;
  tmc_result r;
  memset(&r, 0, sizeof r);
  r.state = TMC_NOT_YET_CALLED;
  return r;
}
