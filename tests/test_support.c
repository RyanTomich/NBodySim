#include "test_support.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

Simulation *from_bodies(const ReferenceBody *bodies, size_t count)
{
    Simulation *simulation;
    assert(simulation_create(&simulation, count));
    MutableSimulationView v = simulation_mutable_view(simulation);
    for (size_t i = 0; i < count; ++i) {
        v.w[i] = bodies[i].w;
        v.x[i] = bodies[i].position.x;
        v.y[i] = bodies[i].position.y;
        v.z[i] = bodies[i].position.z;
        v.vx[i] = bodies[i].velocity.x;
        v.vy[i] = bodies[i].velocity.y;
        v.vz[i] = bodies[i].velocity.z;
        v.ax[i] = bodies[i].acceleration.x;
        v.ay[i] = bodies[i].acceleration.y;
        v.az[i] = bodies[i].acceleration.z;
    }
    return simulation;
}

ReferenceBody body_at(const Simulation *simulation, size_t i)
{
    SimulationView v = simulation_view(simulation);
    return (ReferenceBody){v.w[i], {v.x[i], v.y[i], v.z[i]},
                          {v.vx[i], v.vy[i], v.vz[i]}, {v.ax[i], v.ay[i], v.az[i]}};
}

void assert_close(double actual, double expected, double tolerance)
{
    if (!isfinite(actual) || !isfinite(expected)
        || fabs(actual - expected) > tolerance * (1 + fabs(expected))) {
        fprintf(stderr, "actual %.17g, expected %.17g, tolerance %.3g\n",
                actual, expected, tolerance);
        assert(false);
    }
}

void assert_body(ReferenceBody a, ReferenceBody b, double tolerance)
{
    assert(a.w == b.w);
    assert_close(a.position.x, b.position.x, tolerance);
    assert_close(a.position.y, b.position.y, tolerance);
    assert_close(a.position.z, b.position.z, tolerance);
    assert_close(a.velocity.x, b.velocity.x, tolerance);
    assert_close(a.velocity.y, b.velocity.y, tolerance);
    assert_close(a.velocity.z, b.velocity.z, tolerance);
    assert_close(a.acceleration.x, b.acceleration.x, tolerance);
    assert_close(a.acceleration.y, b.acceleration.y, tolerance);
    assert_close(a.acceleration.z, b.acceleration.z, tolerance);
}

/* Independent ordered-pair Newtonian reference. No production gravity helpers,
 * tree traversal, shared-pair arithmetic, or reduction-order assumptions. */
void reference_accelerations(const ReferenceBody *bodies, size_t count, Vec3 *out)
{
    for (size_t i = 0; i < count; ++i) {
        long double ax = 0, ay = 0, az = 0;
        for (size_t j = 0; j < count; ++j) {
            if (i == j) {
                continue;
            }
            long double dx = (long double)bodies[j].position.x - bodies[i].position.x;
            long double dy = (long double)bodies[j].position.y - bodies[i].position.y;
            long double dz = (long double)bodies[j].position.z - bodies[i].position.z;
            long double distance = hypotl(hypotl(dx, dy), dz);
            if (distance == 0) {
                continue;
            }
            long double scale = bodies[j].w / distance / distance;
            ax += scale * (dx / distance);
            ay += scale * (dy / distance);
            az += scale * (dz / distance);
        }
        out[i] = (Vec3){(double)ax, (double)ay, (double)az};
    }
}
