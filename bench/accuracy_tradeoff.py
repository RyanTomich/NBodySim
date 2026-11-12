#!/usr/bin/env python3
"""Compare serial Barnes–Hut/FMM accuracy and cost against a shared direct reference."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import run as runner


def configurations(thetas):
    return ([("barnes-hut", theta, None) for theta in thetas]
            + [("fmm", None, order) for order in range(1, 7)])


def case_name(solver, theta, order):
    if solver == "barnes-hut":
        return "barnes-hut-theta-" + str(theta)
    return "fmm-order-" + str(order) if solver == "fmm" else "n2"


def build(directory, openmp=False):
    setting = "OPENMP={}".format(int(openmp))
    prefix = "reference-" if openmp else ""
    with (directory / (prefix + "build.stdout.log")).open("w") as log:
        for command in (["make", "--no-print-directory", "clean"],
                        ["make", "--no-print-directory", setting, "print-config"],
                        ["make", "--no-print-directory", setting, "build/bench"]):
            subprocess.run(command, cwd=runner.ROOT, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
    binary = directory / (prefix + "bench")
    shutil.copy2(runner.ROOT / "build/bench", binary)
    info = subprocess.run([str(binary), "--build-info"], capture_output=True,
                          text=True, check=True)
    if json.loads(info.stdout) != {"openmp": openmp, "instrumented": False}:
        raise ValueError("Accuracy executable has unexpected build settings")
    return binary


def collect_reference(binary, instance, directory, args):
    """Produce one untimed direct trajectory, optionally using OpenMP."""
    directory.mkdir()
    status = {"status": "failed", "mode": "reference", "solver": "n2",
              "integrator": "leapfrog", "dt": args.dt, "steps": args.steps,
              "samples": 1, "warmups": 0, "threads": args.reference_threads,
              "openmp": args.reference_threads > 1, "started": runner.now()}
    try:
        command = [str(binary), str(instance), str(args.dt), str(args.steps),
                   "1", "0", "reference", str(directory / "final.txt"), "n2"]
        environment = {**os.environ, "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
                       "OMP_NUM_THREADS": str(args.reference_threads),
                       "OMP_THREAD_LIMIT": str(args.reference_threads)}
        with (directory / "stdout.log").open("w") as out, (directory / "stderr.log").open("w") as err:
            result = subprocess.run(command, cwd=runner.ROOT, env=environment,
                                    stdin=subprocess.DEVNULL, stdout=out, stderr=err)
        if result.returncode:
            raise RuntimeError("Direct reference failed; see reference/stderr.log")
        final = runner.read_state(directory / "final.txt")
        runner.validate_bodies(final, runner.read_state(instance))
        status["status"] = "ok"
        return final
    except (OSError, ValueError, RuntimeError, KeyboardInterrupt) as error:
        status["error"] = str(error) or "Interrupted"
        raise
    finally:
        status["finished"] = runner.now()
        runner.write_json(directory / "status.json", status)


def read_samples(path, samples):
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source))
    if [int(row["sample"]) for row in rows] != list(range(1, samples + 1)):
        raise ValueError("Missing or duplicate timing samples")
    durations = [float(row["seconds"]) for row in rows]
    if any(not math.isfinite(value) or value <= 0 for value in durations):
        raise ValueError("Expected finite positive timing samples")
    return durations


def experiment(args):
    directory = runner.create_experiment_directory(args.results_dir)
    settings = {key: getattr(args, key) for key in
                ("bodies", "dt", "steps", "samples", "warmups", "seed", "thetas", "reference_threads")}
    settings.update(threads=1, openmp=False, orders=list(range(1, 7)),
                    generator={"L": 10, "w_min": 1, "w_max": 100, "v_max": 3, "a_max": 0})
    status = {"status": "running", "suite": "accuracy-tradeoff", "settings": settings,
              "started": runner.now()}
    runner.write_json(directory / "status.json", status)
    rows = []
    try:
        binary = build(directory)
        reference_binary = build(directory, openmp=True) if args.reference_threads > 1 else binary
        instance = directory / "input.txt"
        command = [sys.executable, runner.ROOT / "generate_instance.py", "--n", str(args.bodies),
                   "--L", "10", "--w-min", "1", "--w-max", "100", "--v-max", "3",
                   "--a-max", "0", "--seed", str(args.seed), "--output", instance]
        if runner.run_logged(command, directory, stdout="generation.stdout.log",
                             stderr="generation.stderr.log"):
            raise RuntimeError("Input generation failed; see generation.stderr.log")
        initial = runner.read_state(instance)
        print("N={}: n2 reference (one untimed run, {} thread(s))".format(
            args.bodies, args.reference_threads), flush=True)
        trajectory_reference = collect_reference(reference_binary, instance, directory / "reference", args)
        for solver, theta, order in configurations(args.thetas):
            name = case_name(solver, theta, order)
            print("N={}: {} ({} measured runs)".format(args.bodies, name, args.samples), flush=True)
            case = directory / name
            case.mkdir()
            options = argparse.Namespace(**vars(args), solver=solver, theta=theta, order=order)
            timing_status, timing_summary = runner.collect_timing(binary, instance, case / "timing", options)
            if timing_status["status"] != "ok":
                raise RuntimeError(timing_status["error"])
            read_samples(case / "timing/samples.csv", args.samples)
            final = runner.read_state(case / "timing/final.txt")
            runner.validate_bodies(final, initial)
            position_error = runner.trajectory_errors(final, trajectory_reference)[0]
            runner.write_csv(case / "timing/error.csv", list(position_error), [position_error])
            row = {"case": name, "bodies": args.bodies, "solver": solver, "theta": theta,
                   "order": order, "threads": 1, "dt": args.dt, "steps": args.steps,
                   "samples": args.samples, "warmups": args.warmups,
                   "mean_seconds_per_step": timing_summary["mean_seconds_per_step"],
                   "stddev_seconds_per_step": timing_summary["stddev_seconds_per_step"],
                   "position_relative_l2": position_error["relative_l2_error"]}
            rows.append(row)
            runner.write_csv(directory / "accuracy-tradeoff.csv", list(row), rows)
            print("  trajectory {:.4g} ms/step; position error {}".format(
                1000 * row["mean_seconds_per_step"], row["position_relative_l2"]), flush=True)
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
    parser.add_argument("--bodies", type=int, default=1000)
    parser.add_argument("--dt", type=float, default=0.0001)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--reference-threads", type=int, default=1,
                        help="threads for the single untimed n2 reference; >1 requires an OpenMP-capable compiler")
    parser.add_argument("--thetas", nargs="+", type=float, default=[0.2, 0.3, 0.5, 0.7])
    parser.add_argument("--results-dir", type=Path, default=runner.ROOT / "results/accuracy-tradeoff")
    args = parser.parse_args(argv)
    if args.bodies <= 0 or args.steps <= 0 or not math.isfinite(args.dt) or args.dt <= 0:
        parser.error("bodies, steps, and finite dt must be positive")
    if args.samples < 2 or args.warmups < 0:
        parser.error("at least two samples and nonnegative warmups are required")
    if args.reference_threads < 1:
        parser.error("reference threads must be positive")
    if (any(not math.isfinite(theta) or theta < 0 for theta in args.thetas)
            or len(set(args.thetas)) != len(args.thetas)):
        parser.error("thetas must be unique, finite, and nonnegative")
    args.thetas.sort()
    return args


def main():
    try:
        return experiment(parse_args())
    except (OSError, RuntimeError) as error:
        print("Cannot create experiment: {}".format(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
