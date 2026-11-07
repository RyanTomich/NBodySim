#include "laplace_internal.h"
#include "laplace_order4_internal.h"

static void translation_kernel(const CartesianTerms *terms, const double unit[3], double *kernel)
{
    if (terms->order == 4) {
        laplace_kernel_order4(kernel, unit);
        return;
    }
    /* Taylor coefficients of 1/|unit+x| from the homogeneous Legendre
     * recurrence: d H_d = -(2d-1)(unit.x)H_(d-1) -(d-1)|x|²H_(d-2). */
    kernel[0] = 1;
    for (int a = 1; a < terms->terms; ++a) {
        int d = terms->degree[a];
        double linear = 0, quadratic = 0;
        for (int axis = 0; axis < 3; ++axis) {
            int b = terms->lower[a][axis];
            if (b >= 0) {
                linear += unit[axis] * kernel[b];
                int c = terms->lower[b][axis];
                if (c >= 0) {
                    quadratic += kernel[c];
                }
            }
        }
        kernel[a] = -((2 * d - 1) * linear + (d - 1) * quadratic) / d;
    }
}

void laplace_translate(const CartesianTerms *terms, const double *moments,
                       double source_half_width, double target_half_width,
                       const double unit[3], double r, double *kernel, double *local)
{
    translation_kernel(terms, unit, kernel);
    double source_scale[CARTESIAN_MAX_ORDER + 1] = {1}, target_scale[CARTESIAN_MAX_ORDER + 1] = {1};
    for (size_t d = 1; d <= terms->order; ++d) {
        source_scale[d] = source_scale[d - 1] * (-source_half_width / r);
        target_scale[d] = target_scale[d - 1] * (target_half_width / r);
    }
    if (terms->order == 4) {
        laplace_contract_order4(local, moments, kernel, source_scale, target_scale, r);
    } else {
        for (int i = 0; i < terms->splits; ++i) {
            CartesianSplit s = terms->split[i];
            local[s.b] += (moments[s.c] / r) * source_scale[terms->degree[s.c]] *
                          target_scale[terms->degree[s.b]] * s.factor * kernel[s.a];
        }
    }
}

