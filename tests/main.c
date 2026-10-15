#include "test.h"
int test_failures;
int total_failures;
int main(void)
{
    test_constraint_initialization();
    test_dae();
    test_rosenbrock();
    test_jacobian_verification();
    test_process();
    test_constraints();
    test_rates();
    test_sparse();
    test_backward_euler();
    return total_failures ? 1 : 0;
}
