/* The random sequences MICM's solver tests draw, reproduced in C.
 *
 * MICM's LU, linear-solver, and Chapman tests fill their random matrices from
 * std::default_random_engine and std::mt19937 through uniform_int,
 * bernoulli, lognormal, and uniform_real distributions. Those engines are
 * fully specified by the C++ standard; the distributions are not, so these
 * follow libc++ (the library MICM's reference build uses on this machine):
 * default_random_engine is minstd_rand, uniform_int_distribution<int>(0, 1)
 * keeps the low bit of one draw, generate_canonical<double, 53> combines two
 * draws, normal_distribution is Marsaglia's polar method with the second
 * value cached, and lognormal is exp(normal). The multiply-adds are written as
 * fma() because clang contracts them when it compiles libc++'s templates; so
 * written, every compiler reproduces the libc++ values bit for bit (checked
 * against a libc++ build for the first 200 to 1000 draws of each stream). With
 * them the C tests factor and solve the same matrices MICM's tests do. (Every
 * check still compares with an independent computation, so other values would
 * test the same properties.)
 *
 * MICM binds each distribution to its own copy of a default-seeded engine
 * (std::bind copies the engine), so each stream here owns its engine. */
#ifndef TMC_TEST_LIBCXX_RANDOM_H
#define TMC_TEST_LIBCXX_RANDOM_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>

/* ---- engines -------------------------------------------------------------- */

/* std::minstd_rand: x <- 48271 x mod (2^31 - 1), default seed 1, range [1, 2^31 - 2]. */
typedef struct tmc_minstd
{
  uint32_t x;
} tmc_minstd;

static inline tmc_minstd tmc_minstd_create(void)
{
  tmc_minstd g = { 1u };
  return g;
}

static inline uint32_t tmc_minstd_next(tmc_minstd *g)
{
  g->x = (uint32_t)(((uint64_t)g->x * 48271u) % 2147483647u);
  return g->x;
}

/* std::mt19937 */
typedef struct tmc_mt19937
{
  uint32_t state[624];
  int index;
} tmc_mt19937;

static inline tmc_mt19937 tmc_mt19937_create(uint32_t seed)
{
  tmc_mt19937 g;
  g.state[0] = seed;
  for (int i = 1; i < 624; ++i)
    g.state[i] = 1812433253u * (g.state[i - 1] ^ (g.state[i - 1] >> 30)) + (uint32_t)i;
  g.index = 624;
  return g;
}

static inline uint32_t tmc_mt19937_next(tmc_mt19937 *g)
{
  if (g->index >= 624)
  {
    for (int i = 0; i < 624; ++i)
    {
      uint32_t y = (g->state[i] & 0x80000000u) | (g->state[(i + 1) % 624] & 0x7fffffffu);
      uint32_t next = g->state[(i + 397) % 624] ^ (y >> 1);
      if (y & 1u)
        next ^= 0x9908b0dfu;
      g->state[i] = next;
    }
    g->index = 0;
  }
  uint32_t y = g->state[g->index++];
  y ^= y >> 11;
  y ^= (y << 7) & 0x9d2c5680u;
  y ^= (y << 15) & 0xefc60000u;
  y ^= y >> 18;
  return y;
}

/* ---- generate_canonical<double, 53> ---------------------------------------- */

/* minstd_rand: range R = 2^31 - 2 holds 30 whole bits, so two draws. */
static inline double tmc_minstd_canonical(tmc_minstd *g)
{
  const double r = 2147483646.0;
  double sum = (double)(tmc_minstd_next(g) - 1u);
  sum = fma((double)(tmc_minstd_next(g) - 1u), r, sum);
  return sum / (r * r);
}

/* mt19937: range 2^32, two draws. */
static inline double tmc_mt19937_canonical(tmc_mt19937 *g)
{
  const double r = 4294967296.0;
  double sum = (double)tmc_mt19937_next(g);
  sum = fma((double)tmc_mt19937_next(g), r, sum);
  return sum / (r * r);
}

/* ---- distributions ------------------------------------------------------- */

/* uniform_int_distribution<int>(0, 1) on minstd_rand: the low bit of one draw. */
static inline int tmc_minstd_bool(tmc_minstd *g)
{
  return (int)((tmc_minstd_next(g) - 1u) & 1u);
}

/* bernoulli_distribution(p) on mt19937 */
static inline bool tmc_mt19937_bernoulli(tmc_mt19937 *g, double p)
{
  return tmc_mt19937_canonical(g) < p;
}

/* uniform_real_distribution<float>(0, 1) on minstd_rand: generate_canonical<float, 24>
 * takes one draw (30 bits cover 24) and forms its range in float, where
 * 2^31 - 2 rounds to 2^31. */
static inline float tmc_minstd_uniform_float(tmc_minstd *g)
{
  const float r = (float)2147483645u + 1.0f;
  float u = (float)(tmc_minstd_next(g) - 1u) / r;
  return (1.0f - 0.0f) * u + 0.0f;
}

/* lognormal_distribution(m, s): exp of a polar-method normal; one stream per engine. */
typedef struct tmc_lognormal
{
  double m, s;
  bool hot;
  double cached;
} tmc_lognormal;

static inline tmc_lognormal tmc_lognormal_create(double m, double s)
{
  tmc_lognormal d = { m, s, false, 0.0 };
  return d;
}

/* One polar-method normal from a canonical generator (u, v uniform on [-1, 1)). */
#define TMC_POLAR_NORMAL(d, canonical, g, out)                                          \
  do                                                                                    \
  {                                                                                     \
    double up_;                                                                         \
    if ((d)->hot)                                                                       \
    {                                                                                   \
      (d)->hot = false;                                                                 \
      up_ = (d)->cached;                                                                \
    }                                                                                   \
    else                                                                                \
    {                                                                                   \
      double u_, v_, s_;                                                                \
      do                                                                                \
      {                                                                                 \
        u_ = fma(1.0 - -1.0, canonical(g), -1.0);                                       \
        v_ = fma(1.0 - -1.0, canonical(g), -1.0);                                       \
        s_ = fma(u_, u_, v_ * v_);                                                      \
      } while (s_ > 1 || s_ == 0);                                                      \
      double fp_ = sqrt(-2 * log(s_) / s_);                                             \
      (d)->cached = v_ * fp_;                                                           \
      (d)->hot = true;                                                                  \
      up_ = u_ * fp_;                                                                   \
    }                                                                                   \
    (out) = fma(up_, (d)->s, (d)->m);                                                   \
  } while (0)

static inline double tmc_minstd_lognormal(tmc_lognormal *d, tmc_minstd *g)
{
  double normal;
  TMC_POLAR_NORMAL(d, tmc_minstd_canonical, g, normal);
  return exp(normal);
}

static inline double tmc_mt19937_lognormal(tmc_lognormal *d, tmc_mt19937 *g)
{
  double normal;
  TMC_POLAR_NORMAL(d, tmc_mt19937_canonical, g, normal);
  return exp(normal);
}

/* The bound generators MICM's tests build with std::bind: a distribution and
 * its own default-seeded engine. */
typedef struct tmc_bool_stream
{
  tmc_minstd engine;
} tmc_bool_stream;

typedef struct tmc_lognormal_stream
{
  tmc_lognormal dist;
  tmc_minstd engine;
} tmc_lognormal_stream;

static inline tmc_bool_stream tmc_bool_stream_create(void)
{
  tmc_bool_stream s = { tmc_minstd_create() };
  return s;
}

static inline int tmc_bool_stream_next(tmc_bool_stream *s)
{
  return tmc_minstd_bool(&s->engine);
}

static inline tmc_lognormal_stream tmc_lognormal_stream_create(double m, double sigma)
{
  tmc_lognormal_stream s = { tmc_lognormal_create(m, sigma), tmc_minstd_create() };
  return s;
}

static inline double tmc_lognormal_stream_next(tmc_lognormal_stream *s)
{
  return tmc_minstd_lognormal(&s->dist, &s->engine);
}

#endif
