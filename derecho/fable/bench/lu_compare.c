/* Compare the two factorization paths of the generated TS1 kernels bit for bit:
 * jacobian followed by lu_decompose, against jacobian_lu (the Jacobian
 * assembled inside the LU), on random states. Prints the number of LU elements
 * and pivot reciprocals that differ, and the first few.
 *
 *   lu_compare [BLOCKS]
 */
#define _POSIX_C_SOURCE 200809L
#include "ts1.h"
#include "tmc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
  const int nb = argc > 1 ? atoi(argv[1]) : 3;
  const tmc_kernels *k = ts1_kernels();
  if (!k->jacobian_lu)
  {
    printf("no jacobian_lu kernel\n");
    return 1;
  }
  const int W = TMC_WIDTH, n = TS1_NUM_VARIABLES, nr = TS1_NUM_REACTIONS, nnz = TS1_NNZ;
  double *y = aligned_alloc(64, (size_t)nb * n * W * sizeof(double));
  double *rc = aligned_alloc(64, (size_t)nb * nr * W * sizeof(double));
  double *a1 = aligned_alloc(64, (size_t)nb * nnz * W * sizeof(double));
  double *a2 = aligned_alloc(64, (size_t)nb * nnz * W * sizeof(double));
  double *r1 = aligned_alloc(64, (size_t)nb * n * W * sizeof(double));
  double *r2 = aligned_alloc(64, (size_t)nb * n * W * sizeof(double));
  srand(2026);
  for (size_t i = 0; i < (size_t)nb * n * W; ++i)
    y[i] = 1e-12 * pow(10.0, 6.0 * rand() / RAND_MAX);
  for (size_t i = 0; i < (size_t)nb * nr * W; ++i)
    rc[i] = 1e-20 * pow(10.0, 10.0 * rand() / RAND_MAX);
  const double alpha = 1.0 / (0.3 * 0.4358665215);
  memset(a1, 0x7f, (size_t)nb * nnz * W * sizeof(double)); /* NaN-like garbage: anything not written shows */
  memset(a2, 0x7f, (size_t)nb * nnz * W * sizeof(double));
  k->jacobian(y, rc, alpha, a1, nb);
  k->lu_decompose(a1, r1, nb);
  k->jacobian_lu(y, rc, alpha, a2, r2, nb);
  long diff = 0, shown = 0, nan2 = 0, zeros = 0;
  for (size_t i = 0; i < (size_t)nb * nnz * W; ++i)
  {
    if (isnan(a2[i]))
      ++nan2;
    if (a1[i] == a2[i] && memcmp(&a1[i], &a2[i], sizeof(double)) != 0)
      ++zeros; /* +0 against -0 */
    if (a1[i] != a2[i] || (isnan(a1[i]) != isnan(a2[i])))
    {
      ++diff;
      if (shown++ < 8)
        printf("element %zu lane %zu: %.17g vs %.17g (rel %.2e)\n", i / W % nnz, i % W, a1[i], a2[i], fabs(a1[i] - a2[i]) / fabs(a1[i]));
    }
  }
  long rdiff = 0;
  for (size_t i = 0; i < (size_t)nb * n * W; ++i)
    if (memcmp(&r1[i], &r2[i], sizeof(double)) != 0)
      ++rdiff;
  printf("%d blocks: %ld of %ld LU values differ (%ld NaN in jacobian_lu, %ld signed-zero differences), %ld of %ld reciprocals differ\n",
         nb, diff, (long)nb * nnz * W, nan2, zeros, rdiff, (long)nb * n * W);
  return diff || rdiff ? 1 : 0;
}
