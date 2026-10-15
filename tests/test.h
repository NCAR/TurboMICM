#ifndef CK_TEST_H
#define CK_TEST_H
#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* Assertions remain enabled when the production library uses -DNDEBUG. */
#define CHECK(condition) do { if (!(condition)) { fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); return 1; } } while (0)
#define NEAR(actual, expected, tolerance) do { double ck_a_=(actual),ck_b_=(expected),ck_t_=(tolerance); if (!(fabs(ck_a_-ck_b_) <= ck_t_)) { fprintf(stderr,"%s:%d: %.17g != %.17g (tol %.3g)\n",__FILE__,__LINE__,ck_a_,ck_b_,ck_t_);return 1;} } while(0)
#define IS(name) (!strcmp(test, (name)))
int test_initialization(const char *suite, const char *test);
int test_integrator(const char *suite, const char *test);
int test_constraint(const char *suite, const char *test);
int test_linear(const char *suite, const char *test);
int test_rate(const char *suite, const char *test);
#endif
