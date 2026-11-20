#include "laplace_internal.h"
#include "generated/laplace_order4_internal.h"

void laplace_translate(const double *moments,
                       double source_half_width, double target_half_width,
                       const double unit[3], double r, double *local)
{
    double kernel[LAPLACE_TERMS];
    double source_scale[LAPLACE_ORDER + 1] = {1}, target_scale[LAPLACE_ORDER + 1] = {1};
    laplace_kernel_order4(kernel, unit);
    for (int d = 1; d <= LAPLACE_ORDER; ++d) {
        source_scale[d] = source_scale[d - 1] * (-source_half_width / r);
        target_scale[d] = target_scale[d - 1] * (target_half_width / r);
    }
    laplace_contract_order4(local, moments, kernel, source_scale, target_scale, r);
}
