#ifndef CHEM_TEST_H
#define CHEM_TEST_H
#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int test_failures;
#define TEST(prefix,suite,name) void test_##prefix##_##suite##_##name(void)
#define CHECK(condition) do { if (!(condition)) { \
    fprintf(stderr,"%s:%d: %s\n",__FILE__,__LINE__,#condition); ++test_failures; } } while (0)
#define REQUIRE(condition) do { if (!(condition)) { CHECK(condition); return; } } while (0)
#define NEAR(actual,expected,tolerance) do { double actual_=(actual), expected_=(expected), tolerance_=(tolerance); \
    if (!(fabs(actual_-expected_)<=tolerance_)) { \
        fprintf(stderr,"%s:%d: %.17g != %.17g (tolerance %.17g)\n",__FILE__,__LINE__,actual_,expected_,tolerance_); \
        ++test_failures; } } while (0)
#define REL(actual,expected,tolerance) NEAR(actual,expected,(tolerance)*fabs((double)(expected)))
struct test_case { const char *id; void (*function)(void); const char *na_reason; };
#endif
