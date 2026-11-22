#!/bin/sh
set -eu
unset NBODY_INTEGRATOR

# Reset the logical initial state before every call, including after either
# parity of acceleration-buffer swaps and for Euler's first velocity update.
printf '1\n1 1 2 3 4 5 6 2 4 6\n' > build/bench-initial.txt
printf '0\n' > build/bench-empty.txt
for integrator in leapfrog euler; do
    for input in build/bench-initial.txt build/bench-empty.txt tests/data/two_body.txt; do
        for steps in 1 2 3; do
            ./build/app --integrator "$integrator" "$input" 0.25 "$steps" > build/bench-expected.txt
            for warmups in 0 1 2; do
                NBODY_INTEGRATOR="$integrator" ./build/bench "$input" 0.25 "$steps" 3 "$warmups" timing build/bench-final.txt > build/bench-samples.csv 2> build/bench-error.txt
                cmp build/bench-expected.txt build/bench-final.txt
                test "$(wc -l < build/bench-samples.csv)" -eq 4
            done
            NBODY_INTEGRATOR="$integrator" ./build/bench "$input" 0.25 "$steps" 1 0 reference build/bench-final.txt > build/bench-samples.csv 2> build/bench-error.txt
            cmp build/bench-expected.txt build/bench-final.txt
            test ! -s build/bench-samples.csv
        done
    done
done

if NBODY_INTEGRATOR=unknown ./build/bench build/bench-initial.txt 0.25 1 1 0 timing build/bench-final.txt > build/bench-samples.csv 2> build/bench-error.txt; then
    printf 'Expected invalid benchmark integrator to fail.\n' >&2
    exit 1
fi
grep -q 'NBODY_INTEGRATOR must be leapfrog or euler' build/bench-error.txt

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
import csv
import json
import os
from pathlib import Path
import runpy
import subprocess
from tempfile import TemporaryDirectory
from types import SimpleNamespace

runner = runpy.run_path("bench/run.py")
binary = Path("build/bench").resolve()
instance = Path("tests/data/two_body.txt").resolve()
names = json.loads(subprocess.check_output([binary, "--metadata"]))
openmp = json.loads(subprocess.check_output([binary, "--build-info"]))["openmp"]
expected = ["n2", "serial-barnes-hut"]
if openmp:
    expected.append("parallel-barnes-hut")
assert names == expected + ["grouped-barnes-hut", "fmm"]
for solver, integrator in [*((name, "leapfrog") for name in names), ("grouped-barnes-hut", "euler")]:
    calls = 3 if integrator == "leapfrog" else 2
    os.environ["NBODY_INTEGRATOR"] = "euler" if integrator == "leapfrog" else "leapfrog"
    args = SimpleNamespace(solver=solver, theta=None, integrator=integrator,
                           dt=0.0001, steps=2, samples=2, warmups=1)
    with TemporaryDirectory(dir="build") as temporary:
        folder = Path(temporary).resolve()
        status, summary = runner["collect_timing"](binary, instance, folder / "timing", args)
        assert status["status"] == "ok", status
        assert summary["samples"] == args.samples
        if solver == "grouped-barnes-hut":
            assert status["theta"] == 0.5
            assert set(runner["solver_settings"](binary, instance, solver, None)) == {
                "theta", "bucket_capacity", "max_depth", "node_capacity"
            }
            assert status["diagnostics"]["calls"] == calls * (args.samples + args.warmups)
            expected = {
                "calls", "capacity_limited_leaves",
                "source_node_visits", "direct_evaluations", "aggregate_evaluations", "scope"
            }
            if "membership_seconds" in status["diagnostics"]:
                import sys
                sys.path.insert(0, "bench")
                from profile_barnes_hut import phase_samples, summarize, COUNTERS, TIMERS, WORKER_TIMERS
                expected.update(COUNTERS)
                expected.update(TIMERS)
                expected.update(WORKER_TIMERS)
                expected.update(("evaluation_seconds", "min_threads", "max_threads"))
                rows = phase_samples(folder / "timing", args.steps, args.samples, integrator=integrator)
                assert len(rows) == args.samples
                assert all(row["tree_reuses"] == calls for row in rows)
                assert all(row["tree_builds"] == 0 for row in rows)
                assert abs(sum(row["percent_of_total"] for row in summarize(rows, args.steps)) - 100) < 1e-9
            assert set(status["diagnostics"]) == expected
        if solver == "fmm":
            assert status["order"] == summary["order"] == 4
            assert "direct_pair_cutoff" not in status
        assert runner["compare_n2"](binary, instance, folder, args)
        with (folder / "timing/distribution-error.csv").open() as source:
            distances = list(csv.DictReader(source))
        actual = runner["read_state"](folder / "timing/final.txt")
        reference = runner["read_state"](folder / "reference/final.txt")
        assert [float(row["ks_distance"]) for row in distances] == [
            row["ks_distance"] for row in runner["distribution_errors"](actual, reference)]
        for mode in ("timing", "reference"):
            saved = json.loads((folder / mode / "status.json").read_text())
            assert saved["status"] == "ok", saved
            assert saved["integrator"] == integrator
            assert saved["force_evaluations_per_run"] == calls
            assert saved["steps"] == args.steps and saved["dt"] == args.dt

for warmups in (0, 1):
    args = SimpleNamespace(solver="fmm", theta=None, integrator="leapfrog",
                           dt=0.0001, steps=2, samples=2, warmups=warmups)
    with TemporaryDirectory(dir="build") as temporary:
        timing = Path(temporary).resolve() / "timing"
        overflow = Path("tests/data/aggregate_overflow.txt").resolve()
        status, summary = runner["collect_timing"](binary, overflow, timing, args)
        assert status["status"] == "failed" and summary is None
        assert json.loads((timing / "status.json").read_text()) == status
        assert (timing / "stderr.log").read_text().count("fmm: nonfinite multipole") == 1
        assert not (timing / "final.txt").exists()
        assert not (timing / "samples.csv").read_text()

try:
    runner["solver_settings"](binary, instance, "unknown", None)
except RuntimeError as error:
    assert "Unknown solver." in str(error), error
else:
    raise AssertionError("Unknown solver must fail")
PY

printf 'All benchmark tests passed.\n'
