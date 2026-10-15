#ifndef CK_TEST_RANDOM_H
#define CK_TEST_RANDOM_H
#include <stdint.h>
#include <math.h>
/* Fixed 32-bit recurrence and Box-Muller transform: C++ engine algorithms and
   distribution output are not specified consistently across standard libraries. */
static inline double test_uniform(uint32_t *state)
{
  *state = *state * UINT32_C(1664525) + UINT32_C(1013904223);
  return (*state + .5) / 4294967296.0;
}
static inline double test_lognormal(uint32_t *state, double mean, double deviation)
{
  double u = test_uniform(state), v = test_uniform(state);
  return exp(mean + deviation * sqrt(-2 * log(u)) * cos(2 * CK_PI * v));
}
/* Frozen reference-engine stream for comparisons against the legacy Chapman
   solver. This mathematical algorithm is expressed in C and is independent of
   the compiler's C++ standard library: modulus 2^31-1, multiplier 48271,
   53-bit canonical uniforms and the polar normal transform. */
typedef struct
{
  uint32_t value;
  int cached;
  double next;
} test_reference_rng;
static inline double test_reference_uniform(test_reference_rng *rng)
{
  const double range = 2147483646.0;
  rng->value = (uint32_t)((uint64_t)rng->value * 48271 % 2147483647);
  double first = (double)rng->value - 1;
  rng->value = (uint32_t)((uint64_t)rng->value * 48271 % 2147483647);
  return (first + ((double)rng->value - 1) * range) / (range * range);
}
static inline double test_reference_lognormal(test_reference_rng *rng)
{
  double normal;
  if (rng->cached)
  {
    normal = rng->next;
    rng->cached = 0;
  }
  else
  {
    double u, v, r;
    do
    {
      u = 2 * test_reference_uniform(rng) - 1;
      v = 2 * test_reference_uniform(rng) - 1;
      r = u * u + v * v;
    } while (r > 1 || r == 0);
    double scale = sqrt(-2 * log(r) / r);
    normal = u * scale;
    rng->next = v * scale;
    rng->cached = 1;
  }
  return exp(normal * 2 - 2);
}
#endif
