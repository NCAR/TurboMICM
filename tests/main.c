#include "test.h"
int test_failures;
int total_failures;
int main(void)
{
    test_rates();
    return total_failures ? 1 : 0;
}
