#include "test.h"
int test_failures;
int total_failures;
int main(void)
{
    test_aerosol();
    test_external();
    test_constraint_initialization();
    test_dae();
    test_rosenbrock();
    test_jacobian_verification();
    test_process();
    test_rate_store();
    test_state();
    test_matrix();
    test_matrix_functions();
    test_system();
    test_runtime();
    test_sparse_matrix();
    test_sparse_functions();
    test_grouped_views();
    test_constraints();
    test_rates();
    test_sparse();
    test_backward_euler();
    return total_failures ? 1 : 0;
}
