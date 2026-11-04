#include "test_support.h"
#include <stdio.h>

int main(void)
{
    test_state();
    test_solver();
    puts("All unit tests passed.");
    return 0;
}
