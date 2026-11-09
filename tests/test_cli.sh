#!/bin/sh
set -eu

expect_failure() {
    if ./build/app "$@" > build/cli-output.txt 2> build/cli-error.txt; then
        printf 'Expected failure for: %s\n' "$*" >&2
        exit 1
    fi
    test -s build/cli-error.txt
    test ! -s build/cli-output.txt
}

expect_failure
expect_failure tests/data/two_body.txt 0.1
expect_failure tests/data/two_body.txt 0.1 1 extra
expect_failure tests/data/no-such-instance.txt 0.1 1
for dt in 0 -0.1 nan inf 1e999 bad ''; do
    expect_failure tests/data/two_body.txt "$dt" 1
done
for steps in -1 1.5 nan 1x 999999999999999999999999999999999 ''; do
    expect_failure tests/data/two_body.txt 0.1 "$steps"
done
printf '1\n1 0 0\n' > build/cli-invalid.txt
expect_failure build/cli-invalid.txt 0.1 1

./build/app tests/data/two_body.txt 0.1 0 > build/cli-output.txt
cmp tests/data/two_body.txt build/cli-output.txt

printf '0\n' > build/cli-empty.txt
./build/app build/cli-empty.txt 0.1 10 > build/cli-output.txt
cmp build/cli-empty.txt build/cli-output.txt

printf '1\n1 1 2 3 4 5 6 2 4 6\n' > build/cli-one-body.txt
printf '1\n1 2 3.25 4.5 4.5 6 7.5 0 0 0\n' > build/cli-expected.txt
./build/app build/cli-one-body.txt 0.25 1 > build/cli-output.txt
cmp build/cli-expected.txt build/cli-output.txt

# Saved final states can be resumed without changing the result.
./build/app tests/data/two_body.txt 0.01 5 > build/cli-first.txt
./build/app build/cli-first.txt 0.01 5 > build/cli-resumed.txt
./build/app tests/data/two_body.txt 0.01 10 > build/cli-n2.txt
cmp build/cli-resumed.txt build/cli-n2.txt

# Each solver is selectable before or after positional arguments.
for solver in n2 barnes-hut serial-barnes-hut parallel-barnes-hut fmm; do
    for steps in 0 1 2 3; do
        ./build/app tests/data/two_body.txt 0.01 "$steps" > build/cli-default.txt
        ./build/app --solver "$solver" tests/data/two_body.txt 0.01 "$steps" > build/cli-output.txt
        cmp build/cli-default.txt build/cli-output.txt
        ./build/app tests/data/two_body.txt 0.01 "$steps" --solver "$solver" > build/cli-output.txt
        cmp build/cli-default.txt build/cli-output.txt
    done
done
expect_failure tests/data/two_body.txt 0.01 1 --solver
expect_failure tests/data/two_body.txt 0.01 1 --solver unknown
expect_failure tests/data/two_body.txt 0.01 1 --solver direct
expect_failure --solver n2 tests/data/two_body.txt 0.01 1 --solver barnes-hut

# Theta is valid only for Barnes–Hut runs, including full opening at zero.
for theta in 0 0.25 1; do
    ./build/app --theta "$theta" tests/data/two_body.txt 0.01 2 --solver barnes-hut > build/cli-output.txt
    ./build/app --solver barnes-hut --theta "$theta" tests/data/two_body.txt 0.01 2 > build/cli-default.txt
    cmp build/cli-default.txt build/cli-output.txt
done
for theta in -1 nan inf 1e999 bad ''; do
    expect_failure --solver barnes-hut --theta "$theta" tests/data/two_body.txt 0.01 1
done
expect_failure --theta 0.5 tests/data/two_body.txt 0.01 1
expect_failure --solver n2 --theta 0 tests/data/two_body.txt 0.01 1
expect_failure --solver n2 --order 4 tests/data/two_body.txt 0.01 1
expect_failure --solver barnes-hut tests/data/two_body.txt 0.01 1 --theta
expect_failure --solver barnes-hut --theta 0 --theta 1 tests/data/two_body.txt 0.01 1

expect_failure --solver approximate tests/data/two_body.txt 0.01 0
expect_failure --solver fmm --theta 0 tests/data/two_body.txt 0.01 0
expect_failure --solver barnes-hut --order 4 tests/data/two_body.txt 0.01 0
for order in 0 7 1.5; do
    expect_failure --solver fmm --order "$order" tests/data/two_body.txt 0.01 0
done

./build/app --solver barnes-hut tests/data/two_body.txt 0.01 3 > build/cli-alias.txt
./build/app --solver serial-barnes-hut tests/data/two_body.txt 0.01 3 > build/cli-canonical.txt
cmp build/cli-alias.txt build/cli-canonical.txt
expect_failure --solver parallel-barnes-hut --theta nan tests/data/two_body.txt 0.01 1
expect_failure --solver parallel-barnes-hut --order 4 tests/data/two_body.txt 0.01 1

printf 'All CLI tests passed.\n'
