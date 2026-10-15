#ifndef CHEM_TEST_H
#define CHEM_TEST_H
#include "chem.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
extern int test_failures;
extern int total_failures;
#define CHECK(x) do { if (!(x)) { ++test_failures; fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #x); } } while (0)
#define NEAR(a,b,t) do { double va_=(a), vb_=(b), vt_=(t); if (!(fabs(va_-vb_) <= vt_)) { ++test_failures; fprintf(stderr,"%s:%d: %.17g != %.17g (tolerance %.17g)\n",__FILE__,__LINE__,va_,vb_,vt_); } } while(0)
#define REL(a,b,t) do { double expected_=(b); NEAR(a,expected_,fabs(expected_)*(t)); } while(0)
#define RUN(id, call) do { test_failures=0; call; printf("%s %s %s\n",id,test_failures?"FAIL":"PASS",test_failures?"assertion failed":"original checks passed"); total_failures+=test_failures; } while(0)
void test_flowtube(void);
void test_chapman(void);
void test_aerosol(void);
void test_external(void);
void test_constraint_initialization(void);
void test_dae(void);
void test_rosenbrock(void);
void test_jacobian_verification(void);
void test_process(void);
void test_rate_store(void);
void test_state(void);
void test_matrix(void);
void test_matrix_functions(void);
void test_system(void);
void test_runtime(void);
void test_sparse_matrix(void);
void test_sparse_functions(void);
void test_grouped_views(void);
void test_constraints(void);
void test_rates(void);
void test_sparse(void);
void test_backward_euler(void);
#endif
