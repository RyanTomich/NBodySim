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
./build/app tests/data/two_body.txt 0.01 10 > build/cli-direct.txt
cmp build/cli-resumed.txt build/cli-direct.txt

printf 'All CLI tests passed.\n'
