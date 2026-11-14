#!/usr/bin/env python3
"""Serial solver comparison or --parallel scaling sweep for all three solvers.

Each point averages three 20-step trajectories by default, resetting the saved
initial state before each run. Error bars are sample SD across run averages.
"""

import argparse
import csv
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import run as runner

SOLVERS = ("n2", "barnes-hut", "fmm")
PARALLEL_SERIES = tuple((solver, threads) for solver in SOLVERS for threads in (1, 4, 8))


def configurations(count, parallel):
    series = PARALLEL_SERIES if parallel else tuple((solver, 1) for solver in SOLVERS)
    return [(solver, threads) for solver, threads in series
            if not (parallel and count > 50000 and solver == "n2" and threads != 8)]


@contextmanager
def thread_environment(threads):
    settings = {"OMP_NUM_THREADS": str(threads), "OMP_DYNAMIC": "FALSE",
                "OMP_THREAD_LIMIT": str(threads)}
    previous = {key: os.environ.get(key) for key in settings}
    os.environ.update(settings)
    try:
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def check_threads(log, threads, parallel):
    records = [json.loads(line.split(": ", 1)[1]) for line in log.splitlines()
               if line.startswith("thread_check: ")]
    if parallel:
        expected = {"openmp": True, "requested_threads": threads, "effective_threads": threads}
        if records != [expected]:
            raise ValueError("Missing or mismatched OpenMP thread check")
        return expected
    if records:
        raise ValueError("Expected a serial executable")
    return {}


def summarize(durations, steps):
    statistics = runner.timing_summary(durations, steps)
    return {"mean_ms_per_step": 1000 * statistics["mean_seconds_per_step"],
            "stddev_ms_per_step": 1000 * statistics["stddev_seconds_per_step"]}


def build(directory, parallel=False):
    target = "build/bench_steps_openmp" if parallel else "build/bench_steps"
    option = "OPENMP=" + str(int(parallel))
    if "NBODY_PROFILE_BH" in os.environ.get("CPPFLAGS", ""):
        raise ValueError("Solver sweeps require phase profiling to be disabled")
    with (directory / "build.stdout.log").open("w") as log:
        for command in (["make", "--no-print-directory", "clean"],
                        ["make", "--no-print-directory", option, "print-config"],
                        ["make", "--no-print-directory", option, target]):
            subprocess.run(command, cwd=runner.ROOT, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
    binary = directory / "bench"
    shutil.copy2(runner.ROOT / target, binary)
    return binary


def collect(binary, instance, case, solver, args, threads=1):
    case.mkdir()
    status = {"status": "failed", "solver": solver, "threads": threads,
              "steps": args.steps, "dt": args.dt, "samples": args.samples, "warmups": 0,
              "integrator": "leapfrog", "force_evaluations_per_run": args.steps + 1,
              "started": runner.now()}
    try:
        command = [binary, instance, str(args.dt), str(args.steps), str(args.samples), "0",
                   "timing", case / "final.txt", solver]
        if solver == "barnes-hut":
            command.extend(["theta", str(args.theta)])
        elif solver == "fmm":
            command.extend(["order", str(args.order)])
        with thread_environment(threads):
            result = runner.run_logged(command, case, stdout="samples.csv")
        if result:
            raise RuntimeError("Driver failed; see " + str(case / "stderr.log"))
        with (case / "samples.csv").open(newline="") as source:
            samples = list(csv.DictReader(source))
        if [int(row["sample"]) for row in samples] != list(range(1, args.samples + 1)):
            raise ValueError("Missing or duplicate whole-run samples")
        log = (case / "stderr.log").read_text()
        status.update(check_threads(log, threads, getattr(args, "parallel", False)))
        steps = [json.loads(line.split(": ", 1)[1])
                 for line in log.splitlines()
                 if line.startswith("step_timing: ")]
        if [(row["sample"], row["step"]) for row in steps] != [
                (sample, step) for sample in range(1, args.samples + 1)
                for step in range(1, args.steps + 1)]:
            raise ValueError("Missing or duplicate step timings")
        durations = [float(row["seconds"]) for row in samples]
        for sample, seconds in enumerate(durations, 1):
            values = [row["seconds"] for row in steps if row["sample"] == sample]
            if (any(not math.isfinite(x) or x < 0 for x in values)
                    or math.fsum(values) > seconds + 1e-9):
                raise ValueError("Invalid individual-step timings")
        summary = summarize(durations, args.steps)
        final = runner.read_state(case / "final.txt")
        runner.validate_bodies(final, runner.read_state(instance))
        runner.driver_metadata(case, status)
        runner.write_csv(case / "steps.csv", ["sample", "step", "seconds"], steps)
        status["status"] = "ok"
        return summary, final
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        status["error"] = str(error)
        raise
    finally:
        status["finished"] = runner.now()
        runner.write_json(case / "status.json", status)


def experiment(args):
    directory = runner.create_experiment_directory(args.results_dir)
    suite = "parallel-sweep" if args.parallel else "solver-sweep"
    max_speed = 0 if args.distribution == "clusters" else 3
    settings = {"bodies": args.bodies, "steps": args.steps, "dt": args.dt,
                "theta": args.theta, "order": args.order,
                "configurations": [{"bodies": count, "solver": solver, "threads": threads}
                                   for count in args.bodies
                                   for solver, threads in configurations(count, args.parallel)],
                "threads": [1, 4, 8] if args.parallel else 1,
                "samples": args.samples, "warmups": 0, "seed": args.seed,
                "generator": {"distribution": args.distribution,
                              "L": 10 if args.distribution == "uniform" else None,
                              "w_min": 1, "w_max": 100, "v_max": max_speed, "a_max": 0}}
    status = {"status": "running", "suite": suite, "settings": settings,
              "started": runner.now(), "error_bars": "whole_run_sample_sd"}
    runner.write_json(directory / "status.json", status)
    rows = []
    try:
        binary = build(directory, args.parallel)
        for count in args.bodies:
            folder = directory / str(count)
            folder.mkdir()
            instance = folder / "input.txt"
            command = [sys.executable, runner.ROOT / "generate_instance.py", "--n", str(count),
                       "--distribution", args.distribution, "--L", "10",
                       "--w-min", "1", "--w-max", "100", "--v-max", str(max_speed),
                       "--a-max", "0", "--seed", str(args.seed), "--output", instance]
            if runner.run_logged(command, folder, stdout="generation.stdout.log",
                                 stderr="generation.stderr.log"):
                raise RuntimeError("Input generation failed; see " + str(folder))
            reference = None
            for solver, threads in configurations(count, args.parallel):
                print("N={} {}: {} runs × {} steps, {} threads".format(
                    count, solver, args.samples, args.steps, threads), flush=True)
                case = folder / ("{}-{}threads".format(solver, threads) if args.parallel else solver)
                summary, final = collect(binary, instance, case, solver, args, threads)
                if solver == "n2" and reference is None:
                    reference = final
                    reference_threads = threads
                distribution = runner.save_distribution_errors(case, final, reference)
                errors = runner.trajectory_errors(final, reference)
                runner.write_csv(case / "error.csv", list(errors[0]), errors)
                row = {"bodies": count, "solver": solver, "threads": threads, "steps": args.steps,
                       "dt": args.dt, "samples": args.samples, "reference_threads": reference_threads, "theta": args.theta if solver == "barnes-hut" else None,
                       "order": args.order if solver == "fmm" else None, **summary, **distribution}
                if args.parallel:
                    row.update(requested_threads=threads, effective_threads=threads)
                for error in errors:
                    row[error["field"] + "_relative_l2"] = error["relative_l2_error"]
                    row[error["field"] + "_rms"] = error["rms_error"]
                rows.append(row)
                runner.write_csv(directory / (suite + ".csv"), list(rows[0]), rows)
                print("  {:.6g} ms/step; run SD {:.3g} ms".format(
                    row["mean_ms_per_step"], row["stddev_ms_per_step"]), flush=True)
        status["status"] = "ok"
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError, KeyboardInterrupt) as error:
        status.update(status="failed", error=str(error) or "Interrupted by user")
        (directory / "error.txt").write_text(status["error"] + "\n")
    status["finished"] = runner.now()
    runner.write_json(directory / "status.json", status)
    print("Results: {} ({})".format(directory, status["status"]), flush=True)
    return 0 if status["status"] == "ok" else 1


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parallel", action="store_true", help="all solvers at 1/4/8 threads; above 50k use Direct at 8 threads only")
    parser.add_argument("--bodies", nargs="+", type=int, help="body counts (parallel default extends through 100000)")
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--samples", type=int, default=3, help="full runs per configuration (at least two)")
    parser.add_argument("--dt", type=float, default=0.0001)
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--order", type=int, choices=range(1, 7), default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--distribution", choices=("uniform", "clusters"), default="uniform",
                        help="initial positions; clusters start from rest")
    parser.add_argument("--results-dir", type=Path)
    args = parser.parse_args(argv)
    if args.bodies is None:
        args.bodies = [1000, 5000, 10000, 20000, 50000] + ([100000] if args.parallel else [])
    if args.results_dir is None:
        args.results_dir = runner.ROOT / "results" / ("parallel-sweep" if args.parallel else "solver-sweep")
    if any(count <= 0 for count in args.bodies) or len(set(args.bodies)) != len(args.bodies):
        parser.error("body counts must be positive and unique")
    if args.samples < 2:
        parser.error("at least two samples are required for between-run SD")
    if args.steps < 2 or not math.isfinite(args.dt) or args.dt <= 0:
        parser.error("at least two steps and a finite positive dt are required")
    if not math.isfinite(args.theta) or args.theta < 0:
        parser.error("theta must be finite and nonnegative")
    args.bodies.sort()
    return args


def main():
    args = parse_args()
    try:
        return experiment(args)
    except (OSError, RuntimeError) as error:
        print("Cannot create experiment: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
