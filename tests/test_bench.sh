#!/bin/sh
set -eu

# Nonzero input accelerations must be restored, including after odd step counts.
printf '1\n1 1 2 3 4 5 6 2 4 6\n' > build/bench-initial.txt
printf '0\n' > build/bench-empty.txt
for input in build/bench-initial.txt build/bench-empty.txt tests/data/two_body.txt; do
    for steps in 1 2 3; do
        ./build/app "$input" 0.25 "$steps" > build/bench-expected.txt
        for warmups in 0 1 2; do
            ./build/bench "$input" 0.25 "$steps" 3 "$warmups" timing build/bench-final.txt > build/bench-samples.csv 2> build/bench-error.txt
            cmp build/bench-expected.txt build/bench-final.txt
            test "$(wc -l < build/bench-samples.csv)" -eq 4
        done
        ./build/bench "$input" 0.25 "$steps" 1 0 reference build/bench-final.txt > build/bench-samples.csv 2> build/bench-error.txt
        cmp build/bench-expected.txt build/bench-final.txt
        test ! -s build/bench-samples.csv
    done
done

# Finite input can overflow during integration. Reject both warmups and samples.
printf '1\n1 0 0 0 1e308 0 0 0 0 0\n' > build/bench-overflow.txt
for mode in timing reference; do
    for warmups in 0 1; do
        if test "$mode" = reference && test "$warmups" -eq 1; then
            continue
        fi
        rm -f build/bench-invalid-final.txt
        if ./build/bench build/bench-overflow.txt 2 1 1 "$warmups" "$mode" build/bench-invalid-final.txt > build/bench-samples.csv 2> build/bench-error.txt; then
            printf 'Expected nonfinite benchmark state to fail.\n' >&2
            exit 1
        fi
        if test "$warmups" -eq 1; then
            grep -q 'Warmup produced a nonfinite state' build/bench-error.txt
        else
            grep -q 'Sample 1 produced a nonfinite state' build/bench-error.txt
        fi
        test ! -s build/bench-samples.csv
        test ! -e build/bench-invalid-final.txt
    done
done

printf 'All benchmark tests passed.\n'
