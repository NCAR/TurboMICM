#include "test.h"
int test_failures;
int total_failures;
int main(void)
{
    test_rates();
    test_sparse();
    test_backward_euler();
    return total_failures ? 1 : 0;
}
