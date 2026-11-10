#!/bin/sh
set -eu

# Reset the logical initial state before every call. Leapfrog recomputes its
# acceleration, including after either parity of acceleration-buffer swaps.
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

# Exercise saved timing/reference metadata with a short existing-driver run.
python3 - <<'PY'
import json
from pathlib import Path
import runpy
from tempfile import TemporaryDirectory
from types import SimpleNamespace

runner = runpy.run_path("bench/run.py")
binary = Path("build/bench").resolve()
instance = Path("tests/data/two_body.txt").resolve()
for solver in ("n2", "barnes-hut", "fmm"):
    args = SimpleNamespace(solver=solver, theta=None, order=None,
                           dt=0.0001, steps=2, samples=2, warmups=1)
    with TemporaryDirectory(dir="build") as temporary:
        folder = Path(temporary).resolve()
        status, summary = runner["collect_timing"](binary, instance, folder / "timing", args)
        assert status["status"] == "ok", status
        assert summary["samples"] == args.samples
        if solver == "barnes-hut":
            assert status["theta"] == 0.5
            assert set(runner["solver_settings"](binary, instance, solver, None)) == {
                "theta", "bucket_capacity", "max_depth", "node_capacity"
            }
            assert status["diagnostics"]["calls"] == (args.steps + 1) * (args.samples + args.warmups)
            assert set(status["diagnostics"]) == {
                "calls", "capacity_limited_leaves",
                "source_node_visits", "direct_evaluations", "aggregate_evaluations", "scope"
            }
        assert runner["compare_n2"](binary, instance, folder, args)
        for mode in ("timing", "reference"):
            saved = json.loads((folder / mode / "status.json").read_text())
            assert saved["status"] == "ok", saved
            assert saved["integrator"] == "leapfrog"
            assert saved["steps"] == args.steps and saved["dt"] == args.dt

try:
    runner["solver_settings"](binary, instance, "unknown", None)
except RuntimeError as error:
    assert "Unknown solver." in str(error), error
else:
    raise AssertionError("Unknown solver must fail")
PY

printf 'All benchmark tests passed.\n'
