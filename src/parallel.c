/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"
#include <omp.h>
#include <limits.h>

int ch_solve_many(ch_state *const *states, size_t count, double dt,
                  const ch_rosenbrock_parameters *parameters, int threads,
                  int *status, ch_stats *statistics) {
  if ((!states && count) || (!status && count) || (!statistics && count) || count > INT_MAX)
    return CH_INVALID;
  if (!count) return CH_OK;
  for (size_t i = 0; i < count; ++i) if (!states[i]) return CH_INVALID;
  if (threads <= 0) threads = omp_get_max_threads();
  if ((size_t)threads > count) threads = (int)count;
  int failed = 0;
  #pragma omp parallel for num_threads(threads) schedule(static) reduction(+:failed)
  for (size_t i = 0; i < count; ++i) {
    status[i] = ch_solve(states[i],dt,parameters,&statistics[i]);
    failed += status[i] != CH_OK;
  }
  return failed;
}
