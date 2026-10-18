/* Test harness: each MICM case is a C function registered under MICM's case
 * id ("<ctest name>:<Suite.Test>" for GoogleTest cases, "<ctest name>" for
 * programs). `make test-cpu` prints one line per case of the reference list:
 *   <id> PASS|FAIL|N/A <reason>
 *
 * A test returns TMC_TEST_PASS, TMC_TEST_FAIL, or TMC_TEST_NA (after
 * tmc_test_set_reason). The EXPECT macros record a failure with its location
 * and keep going; the ASSERT macros return TMC_TEST_FAIL at once. */
#ifndef TMC_TEST_HARNESS_H
#define TMC_TEST_HARNESS_H

#include "tmc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum
{
  TMC_TEST_PASS = 0,
  TMC_TEST_FAIL = 1,
  TMC_TEST_NA = 2
};

typedef struct tmc_test_case
{
  const char *id;
  int (*run)(void);
} tmc_test_case;

/* Set by tests; printed after N/A. */
void tmc_test_set_reason(const char *reason);
/* Records a failure (the case fails even if the function returns PASS). */
void tmc_test_fail(const char *file, int line, const char *format, ...);
/* Path of a reference file: MICM's source tree (TMC_TEST_MICM), the
 * experiment's reference directory (TMC_TEST_REFERENCE), or this entry's
 * tests/data directory. */
const char *tmc_test_micm_path(const char *relative);
const char *tmc_test_reference_path(const char *relative);
const char *tmc_test_data_path(const char *relative);

#define TMC_TEST(name) static int name(void)
#define TMC_CASE(id, fn) { id, fn }

#define TMC_EXPECT_TRUE(cond)                                                      \
  do                                                                               \
  {                                                                                \
    if (!(cond))                                                                   \
      tmc_test_fail(__FILE__, __LINE__, "expected true: %s", #cond);               \
  } while (0)
#define TMC_EXPECT_FALSE(cond) TMC_EXPECT_TRUE(!(cond))
#define TMC_EXPECT_EQ(a, b)                                                        \
  do                                                                               \
  {                                                                                \
    if (!((a) == (b)))                                                             \
      tmc_test_fail(__FILE__, __LINE__, "expected %s == %s", #a, #b);              \
  } while (0)
#define TMC_EXPECT_NE(a, b)                                                        \
  do                                                                               \
  {                                                                                \
    if ((a) == (b))                                                                \
      tmc_test_fail(__FILE__, __LINE__, "expected %s != %s", #a, #b);              \
  } while (0)
#define TMC_EXPECT_STREQ(a, b)                                                     \
  do                                                                               \
  {                                                                                \
    if (strcmp((a), (b)) != 0)                                                     \
      tmc_test_fail(__FILE__, __LINE__, "expected \"%s\" == \"%s\"", (a), (b));    \
  } while (0)
#define TMC_EXPECT_DOUBLE_EQ(a, b)                                                 \
  do                                                                               \
  {                                                                                \
    double tmc_a_ = (a), tmc_b_ = (b);                                             \
    if (!tmc_test_double_eq(tmc_a_, tmc_b_))                                       \
      tmc_test_fail(__FILE__, __LINE__, "expected %s (%.17g) == %s (%.17g)", #a, tmc_a_, #b, tmc_b_); \
  } while (0)
#define TMC_EXPECT_NEAR(a, b, tol)                                                 \
  do                                                                               \
  {                                                                                \
    double tmc_a_ = (a), tmc_b_ = (b), tmc_t_ = (tol);                             \
    if (!(fabs(tmc_a_ - tmc_b_) <= tmc_t_))                                        \
      tmc_test_fail(__FILE__, __LINE__, "expected |%s - %s| <= %g: %.17g vs %.17g", #a, #b, tmc_t_, tmc_a_, tmc_b_); \
  } while (0)
/* |a - b| <= tol * max(|a|, |b|) (MICM's relative comparison) */
#define TMC_EXPECT_REL(a, b, tol)                                                  \
  do                                                                               \
  {                                                                                \
    double tmc_a_ = (a), tmc_b_ = (b), tmc_t_ = (tol);                             \
    double tmc_s_ = fabs(tmc_a_) > fabs(tmc_b_) ? fabs(tmc_a_) : fabs(tmc_b_);     \
    if (!(fabs(tmc_a_ - tmc_b_) <= tmc_t_ * tmc_s_))                               \
      tmc_test_fail(__FILE__, __LINE__, "expected %s ~ %s within %g: %.17g vs %.17g", #a, #b, tmc_t_, tmc_a_, tmc_b_); \
  } while (0)
#define TMC_ASSERT_TRUE(cond)                                                      \
  do                                                                               \
  {                                                                                \
    if (!(cond))                                                                   \
    {                                                                              \
      tmc_test_fail(__FILE__, __LINE__, "required: %s", #cond);                    \
      return TMC_TEST_FAIL;                                                        \
    }                                                                              \
  } while (0)
#define TMC_NA(reason)                                                             \
  do                                                                               \
  {                                                                                \
    tmc_test_set_reason(reason);                                                   \
    return TMC_TEST_NA;                                                            \
  } while (0)

/* GoogleTest's EXPECT_DOUBLE_EQ: within 4 units in the last place. */
bool tmc_test_double_eq(double a, double b);

/* Every test file exports its cases through one of these. */
const tmc_test_case *tmc_tests_version(int *n);

#endif
