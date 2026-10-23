#include "gravity.h"
#include "io.h"
#include "nbody.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void assert_vec(Vec3 actual, Vec3 expected)
{
    assert(fabs(actual.x - expected.x) < 1e-12);
    assert(fabs(actual.y - expected.y) < 1e-12);
    assert(fabs(actual.z - expected.z) < 1e-12);
}

static void assert_body_equal(Body actual, Body expected)
{
    assert(actual.w == expected.w);
    assert(actual.position.x == expected.position.x);
    assert(actual.position.y == expected.position.y);
    assert(actual.position.z == expected.position.z);
    assert(actual.velocity.x == expected.velocity.x);
    assert(actual.velocity.y == expected.velocity.y);
    assert(actual.velocity.z == expected.velocity.z);
    assert(actual.acceleration.x == expected.acceleration.x);
    assert(actual.acceleration.y == expected.acceleration.y);
    assert(actual.acceleration.z == expected.acceleration.z);
}

static void test_gravity(void)
{
    Body target = {.w = 7.0, .position = {1.0, 1.0, 1.0}};
    Body source = {.w = 27.0, .position = {2.0, 3.0, 3.0}};
    assert_vec(gravity(&target, &source), (Vec3){1.0, 2.0, 2.0});

    target.w = 99.0; /* Acceleration does not depend on target mass. */
    assert_vec(gravity(&target, &source), (Vec3){1.0, 2.0, 2.0});
    source.w = 54.0;
    assert_vec(gravity(&target, &source), (Vec3){2.0, 4.0, 4.0});

    source.w = 27.0;
    source.position = (Vec3){3.0, 5.0, 5.0}; /* Twice the distance. */
    assert_vec(gravity(&target, &source), (Vec3){0.25, 0.5, 0.5});
    source.position = (Vec3){0.0, -1.0, -1.0};
    assert_vec(gravity(&target, &source), (Vec3){-1.0, -2.0, -2.0});

    source.position = target.position;
    assert_vec(gravity(&target, &source), (Vec3){0.0, 0.0, 0.0});
}

static void test_euler(void)
{
    Body body = {
        .w = 1.0,
        .position = {1.0, 2.0, 3.0},
        .velocity = {4.0, 5.0, 6.0},
        .acceleration = {2.0, 4.0, 6.0}
    };
    assert(simulate(&body, 1, 0.25, 1, gravity));
    assert_vec(body.position, (Vec3){2.0, 3.25, 4.5});
    assert_vec(body.velocity, (Vec3){4.5, 6.0, 7.5});
    assert_vec(body.acceleration, (Vec3){0.0, 0.0, 0.0});

    assert(simulate(&body, 1, 0.25, 2, gravity));
    assert_vec(body.position, (Vec3){4.25, 6.25, 8.25});
    assert_vec(body.velocity, (Vec3){4.5, 6.0, 7.5});
}

static void test_input_acceleration(void)
{
    Body bodies[] = {
        {.w = 1.0, .position = {0.0, 0.0, 0.0}, .acceleration = {1.0, 0.0, 0.0}},
        {.w = 1.0, .position = {2.0, 0.0, 0.0}, .acceleration = {-1.0, 0.0, 0.0}}
    };
    Body together[] = {bodies[0], bodies[1]};
    assert(simulate(bodies, 2, 0.5, 1, gravity));
    assert_vec(bodies[0].position, (Vec3){0.0, 0.0, 0.0});
    assert_vec(bodies[1].position, (Vec3){2.0, 0.0, 0.0});
    assert_vec(bodies[0].velocity, (Vec3){0.5, 0.0, 0.0});
    assert_vec(bodies[1].velocity, (Vec3){-0.5, 0.0, 0.0});
    assert_vec(bodies[0].acceleration, (Vec3){0.25, 0.0, 0.0});
    assert_vec(bodies[1].acceleration, (Vec3){-0.25, 0.0, 0.0});

    assert(simulate(bodies, 2, 0.5, 1, gravity));
    assert_vec(bodies[0].position, (Vec3){0.25, 0.0, 0.0});
    assert_vec(bodies[1].position, (Vec3){1.75, 0.0, 0.0});
    assert_vec(bodies[0].velocity, (Vec3){0.625, 0.0, 0.0});
    assert_vec(bodies[1].velocity, (Vec3){-0.625, 0.0, 0.0});
    assert_vec(bodies[0].acceleration, (Vec3){4.0 / 9.0, 0.0, 0.0});
    assert_vec(bodies[1].acceleration, (Vec3){-4.0 / 9.0, 0.0, 0.0});

    assert(simulate(together, 2, 0.5, 2, gravity));
    assert_body_equal(together[0], bodies[0]);
    assert_body_equal(together[1], bodies[1]);
}

static const Body *observed_bodies;
static const Body snapshot[] = {
    {1.0, {0.5, 2.0, 3.5}, {6.0, 12.0, 18.0}, {10.0, 20.0, 30.0}},
    {2.0, {5.0, 6.5, 8.0}, {24.0, 30.0, 36.0}, {40.0, 50.0, 60.0}},
    {3.0, {9.5, 11.0, 12.5}, {42.0, 48.0, 54.0}, {70.0, 80.0, 90.0}}
};
static size_t calls;
static unsigned pair_visits[3][3];

static Vec3 custom_interaction(const Body *target, const Body *source)
{
    assert(target != source);
    ++calls;
    ++pair_visits[target - observed_bodies][source - observed_bodies];
    for (size_t i = 0; i < 3; ++i) {
        assert_body_equal(observed_bodies[i], snapshot[i]);
    }
    return (Vec3){source->w, source->velocity.y, source->acceleration.z};
}

static void test_callback_snapshot(void)
{
    Body bodies[] = {
        {1.0, {0.0, 1.0, 2.0}, {1.0, 2.0, 3.0}, {10.0, 20.0, 30.0}},
        {2.0, {3.0, 4.0, 5.0}, {4.0, 5.0, 6.0}, {40.0, 50.0, 60.0}},
        {3.0, {6.0, 7.0, 8.0}, {7.0, 8.0, 9.0}, {70.0, 80.0, 90.0}}
    };
    observed_bodies = bodies;
    assert(simulate(bodies, 3, 0.5, 1, custom_interaction));
    assert(calls == 3 * 2);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            assert(pair_visits[i][j] == (i == j ? 0u : 1u));
        }
    }
    assert_vec(bodies[0].acceleration, (Vec3){5.0, 78.0, 150.0});
    assert_vec(bodies[1].acceleration, (Vec3){4.0, 60.0, 120.0});
    assert_vec(bodies[2].acceleration, (Vec3){3.0, 42.0, 90.0});
}

static void test_zero_steps_and_arguments(void)
{
    Body body = {2.0, {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    Body original = body;
    size_t previous_calls = calls;
    assert(simulate(&body, 1, 0.1, 0, custom_interaction));
    assert_body_equal(body, original);
    assert(simulate(NULL, 0, 0.1, 10, custom_interaction));
    assert(calls == previous_calls);

    assert(!simulate(NULL, 1, 0.1, 1, gravity));
    assert(!simulate(&body, 1, 0.0, 1, gravity));
    assert(!simulate(&body, 1, -0.1, 1, gravity));
    assert(!simulate(&body, 1, NAN, 1, gravity));
    assert(!simulate(&body, 1, INFINITY, 1, gravity));
    assert(!simulate(&body, 1, 0.1, 1, NULL));
    assert(!simulate(&body, SIZE_MAX, 0.1, 1, gravity));
    assert_body_equal(body, original);
}

static FILE *text_file(const char *text)
{
    FILE *file = tmpfile();
    assert(file);
    assert(fputs(text, file) >= 0);
    rewind(file);
    return file;
}

static void test_file_roundtrip(void)
{
    Body original[] = {
        {0.1, {1e-180, -1e180, 3.141592653589793}, {0.0, -0.5, 0.75}, {2.0, 3.0, -4.0}},
        {2.5, {nextafter(1.0, 2.0), 8.0, 9.0}, {10.0, 11.0, 12.0}, {13.0, 14.0, 15.0}}
    };
    FILE *file = tmpfile();
    assert(file);
    assert(write_bodies(file, original, 2));
    rewind(file);
    Body *loaded;
    size_t count;
    assert(read_bodies(file, &loaded, &count));
    assert(count == 2);
    assert_body_equal(loaded[0], original[0]);
    assert_body_equal(loaded[1], original[1]);
    free(loaded);
    assert(fclose(file) == 0);

    file = text_file(" \n0\t\n");
    assert(read_bodies(file, &loaded, &count));
    assert(count == 0 && loaded == NULL);
    assert(fclose(file) == 0);

    file = text_file("1\n1 2 3 4 5 6 7 8 9 10"); /* No trailing newline. */
    assert(read_bodies(file, &loaded, &count));
    assert(count == 1);
    assert_body_equal(loaded[0], (Body){1.0, {2.0, 3.0, 4.0}, {5.0, 6.0, 7.0}, {8.0, 9.0, 10.0}});
    free(loaded);
    assert(fclose(file) == 0);

    file = fopen("tests/data/two_body.txt", "r");
    assert(file);
    assert(read_bodies(file, &loaded, &count));
    assert(count == 2);
    assert_vec(loaded[0].acceleration, gravity(&loaded[0], &loaded[1]));
    assert_vec(loaded[1].acceleration, gravity(&loaded[1], &loaded[0]));
    free(loaded);
    assert(fclose(file) == 0);
}

static void test_invalid_input(void)
{
    const char *invalid[] = {
        "", " \t\n", "-1", "1.5", "no", "0 extra", "0 1",
        "9999999999999999999999999999999999999999999999999999999999",
        "1\n1 0 0 0 0 0 0 0 0", /* Missing field. */
        "1\n1 0 0 0 0 0 0 0 0 0 7", /* Extra field. */
        "2\n1 0 0 0 0 0 0 0 0 0", /* Missing body. */
        "1\n0 0 0 0 0 0 0 0 0 0", /* Zero mass. */
        "1\n-1 0 0 0 0 0 0 0 0 0",
        "1\ninf 0 0 0 0 0 0 0 0 0",
        "1\n1 nan 0 0 0 0 0 0 0 0",
        "1\n1 0 0 0 1e999 0 0 0 0 0",
        "1\n1 0 0 0 0 0 0 0 0 text",
        "1\n1 0-1 0 0 0 0 0 0 0"
    };
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
        FILE *file = text_file(invalid[i]);
        Body *bodies;
        size_t count;
        assert(!read_bodies(file, &bodies, &count));
        assert(bodies == NULL && count == 0);
        assert(fclose(file) == 0);
    }

    char maximum[32];
    size_t count;
    assert(snprintf(maximum, sizeof(maximum), "%zu", SIZE_MAX) > 0);
    assert(parse_size(maximum, &count) && count == SIZE_MAX);
    FILE *file = text_file(maximum); /* Allocation size overflow. */
    Body *bodies;
    assert(!read_bodies(file, &bodies, &count));
    assert(bodies == NULL && count == 0);
    assert(fclose(file) == 0);

    file = tmpfile();
    assert(file);
    for (int i = 0; i < 200; ++i) {
        assert(fputc('0', file) != EOF);
    }
    rewind(file);
    assert(!read_bodies(file, &bodies, &count));
    assert(fclose(file) == 0);
}

int main(void)
{
    test_gravity();
    test_euler();
    test_input_acceleration();
    test_callback_snapshot();
    test_zero_steps_and_arguments();
    test_file_roundtrip();
    test_invalid_input();
    puts("All unit tests passed.");
    return 0;
}
