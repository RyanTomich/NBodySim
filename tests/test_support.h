#ifndef TEST_SUPPORT_H
#define TEST_SUPPORT_H
#include "nbody.h"

typedef struct {
    double w;
    Vec3 position, velocity, acceleration;
} ReferenceBody;

Simulation *from_bodies(const ReferenceBody *bodies, size_t count);
ReferenceBody body_at(const Simulation *simulation, size_t i);
void assert_close(double actual, double expected, double tolerance);
void assert_body(ReferenceBody actual, ReferenceBody expected, double tolerance);
void reference_accelerations(const ReferenceBody *bodies, size_t count, Vec3 *out);
void test_state(void);
void test_solver(void);
#endif
