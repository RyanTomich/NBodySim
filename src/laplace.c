#include "laplace_internal.h"
#include "generated/laplace_order4_internal.h"

static int choose(int n, int k)
{
    int result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * (n - i + 1) / i;
    }
    return result;
}

void laplace_prepare_terms(CartesianTerms *terms)
{
    int index[CARTESIAN_MAX_ORDER + 1][CARTESIAN_MAX_ORDER + 1][CARTESIAN_MAX_ORDER + 1] = {0};
    for (int d = 0; d <= (int)terms->order; ++d) {
        for (int x = 0; x <= d; ++x) {
            for (int y = 0; y <= d - x; ++y) {
                int z = d - x - y, a = terms->terms++;
                index[x][y][z] = a;
                terms->degree[a] = d;
                terms->exponent[a][0] = x;
                terms->exponent[a][1] = y;
                terms->exponent[a][2] = z;
                terms->lower[a][0] = x ? index[x - 1][y][z] : -1;
                terms->lower[a][1] = y ? index[x][y - 1][z] : -1;
                terms->lower[a][2] = z ? index[x][y][z - 1] : -1;
                for (int i = 0; i <= x; ++i) {
                    for (int j = 0; j <= y; ++j) {
                        for (int k = 0; k <= z; ++k) {
                            terms->split[terms->splits++] =
                                (CartesianSplit){a, index[i][j][k], index[x - i][y - j][z - k],
                                                 choose(x, i) * choose(y, j) * choose(z, k)};
                        }
                    }
                }
            }
        }
    }
}

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

