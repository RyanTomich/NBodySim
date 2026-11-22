#!/usr/bin/env python3
"""Compare grouped Barnes–Hut and REBOUND on one saved initial state."""

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys

from run import (ROOT, create_experiment_directory, driver_metadata, read_state,
                 run_logged, timing_summary, trajectory_errors, validate_bodies,
                 write_csv, write_json)


GROUPED = "grouped-barnes-hut"
REBOUND = "rebound"
THREADS = 4
DT, THETA = 0.0001, 0.5
GENERATOR = {"distribution": "uniform", "seed": 42, "L": 10, "w-min": 1,
             "w-max": 100, "v-max": 3, "a-max": 0}


def prepare_rebound(directory):
    submodule = "third_party/rebound"
    command = ["git", "submodule", "status", "--", submodule]
    pin = subprocess.check_output(command, cwd=ROOT, text=True)
    if pin.startswith("-"):
        if run_logged(["git", "submodule", "update", "--init", "--", submodule],
                      directory, stdout="submodule.log", stderr=None):
            raise RuntimeError("Could not initialize REBOUND. See submodule.log")
        pin = subprocess.check_output(command, cwd=ROOT, text=True)
    if not pin.startswith(" ") or not (ROOT / submodule / "src/rebound.h").is_file():
        raise RuntimeError("Expected the pinned REBOUND checkout. Check git submodule status -- " + submodule)
    return pin.split()[0]


def build(directory):
    if run_logged(["make", "-B", "OPENMP=1", "print-config", "rebound-benchmark"],
                  directory, stdout="build.stdout.log", stderr=None):
        raise RuntimeError("Build failed. See build.stdout.log")
    for name in ("bench_grouped_compare", "bench_rebound", "librebound.so"):
        shutil.copy2(ROOT / "build" / name, directory / name)


def sample(directory, name, instance, dt, steps, folder):
    folder.mkdir()
    final = folder / "final.txt"
    if name == GROUPED:
        command = [directory / "bench_grouped_compare", instance, dt, steps, 1, 0,
                   "timing", final, GROUPED, "theta", THETA]
    else:
        command = [directory / "bench_rebound", instance, dt, steps, final]
    env = {"OMP_NUM_THREADS": str(THREADS), "OMP_DYNAMIC": "FALSE",
           "OMP_THREAD_LIMIT": str(THREADS), "NBODY_INTEGRATOR": "leapfrog"}
    if run_logged(command, folder, stdout="sample.csv", env=env):
        raise RuntimeError("Driver failed. See " + str(folder / "stderr.log"))
    with (folder / "sample.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1 or rows[0]["sample"] != "1":
        raise ValueError("Expected one sample in " + str(folder))
    seconds = float(rows[0]["seconds"])
    if not math.isfinite(seconds) or seconds <= 0:
        raise ValueError("Invalid duration in " + str(folder))

    lines = (folder / "stderr.log").read_text().splitlines()
    checks = [json.loads(line.split(": ", 1)[1]) for line in lines if line.startswith("thread_check: ")]
    if len(checks) != 1 or checks[0]["requested_threads"] != THREADS or checks[0]["effective_threads"] != THREADS:
        raise ValueError("Expected four OpenMP workers in " + str(folder))
    metadata = {}
    driver_metadata(folder, metadata)
    calls = metadata["diagnostics"]["calls"]
    if calls != steps + (name == GROUPED):
        raise ValueError("Unexpected force evaluation count in " + str(folder))
    return seconds, calls, read_state(final)


def validate_adapter(directory):
    fixture = ROOT / "tests/data/two_body.txt"
    dt, steps = 0.01, 3
    _, _, actual = sample(directory, REBOUND, fixture, dt, steps, directory / "two-body")
    expected = read_state(fixture)
    left, right = expected
    for _ in range(steps):
        for body in expected:
            for axis in range(3):
                body[1 + axis] += 0.5 * dt * body[4 + axis]
        delta = [right[1 + axis] - left[1 + axis] for axis in range(3)]
        inverse_cube = math.hypot(*delta) ** -3
        for axis, displacement in enumerate(delta):
            left[7 + axis] = right[0] * displacement * inverse_cube
            right[7 + axis] = -left[0] * displacement * inverse_cube
        for body in expected:
            for axis in range(3):
                body[4 + axis] += dt * body[7 + axis]
                body[1 + axis] += 0.5 * dt * body[4 + axis]
    validate_bodies(actual, expected)
    if any(not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
           for row, ref in zip(actual, expected) for a, b in zip(row, ref)):
        raise ValueError("REBOUND failed the independent two-body drift–kick–drift check")


def compare(directory, bodies, steps, pairs):
    instance = directory / "input.txt"
    command = [sys.executable, ROOT / "generate_instance.py", "--n", bodies, "--output", instance]
    for key, value in GENERATOR.items():
        command.extend(["--" + key, value])
    if run_logged(command, directory, stdout="generation.log", stderr=None):
        raise RuntimeError("Input generation failed. See generation.log")
    initial = read_state(instance)
    if len(initial) != bodies:
        raise ValueError("Generated input has the wrong number of bodies")
    samples, differences = [], []
    for pair in range(1, pairs + 1):
        states = {}
        order = (GROUPED, REBOUND) if pair % 2 else (REBOUND, GROUPED)
        for name in order:
            folder = directory / f"pair-{pair}-{name}"
            seconds, calls, states[name] = sample(directory, name, instance, DT, steps, folder)
            validate_bodies(states[name], initial)
            samples.append({"pair": pair, "solver": name, "seconds": seconds,
                            "steps": steps, "threads": THREADS, "force_evaluations": calls})
        metrics = trajectory_errors(states[REBOUND], states[GROUPED])
        difference = {"pair": pair}
        for row in metrics:
            if row["field"] in ("position", "velocity"):
                difference[row["field"] + "_rms"] = row["rms_error"]
        differences.append(difference)

    summary = []
    for name in (GROUPED, REBOUND):
        durations = [row["seconds"] for row in samples if row["solver"] == name]
        timing = timing_summary(durations, steps)
        summary.append({"solver": name, "samples": timing["samples"],
                        "mean_ms_per_step": 1000 * timing["mean_seconds_per_step"],
                        "stddev_ms_per_step": 1000 * timing["stddev_seconds_per_step"]})
    for filename, rows in (("samples.csv", samples), ("summary.csv", summary),
                           ("trajectory-differences.csv", differences)):
        write_csv(directory / filename, list(rows[0]), rows)
    for row in summary:
        print("{solver}: {mean_ms_per_step:.4f} ± {stddev_ms_per_step:.4f} ms/step (sample SD)".format(**row))
    ratio = summary[1]["mean_ms_per_step"] / summary[0]["mean_ms_per_step"]
    print(f"REBOUND mean time / grouped Barnes–Hut mean time = {ratio:.4f}")
    for row in differences:
        print("Pair {pair} trajectory differences: position RMS {position_rms:.6g}, "
              "velocity RMS {velocity_rms:.6g}".format(**row))
    return ratio


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true",
                        help="run two pairs of two-step trajectories on 64 bodies")
    args = parser.parse_args()
    bodies, steps, pairs = (64, 2, 2) if args.smoke else (10000, 20, 3)
    directory = create_experiment_directory(ROOT / "results/rebound")
    print(directory, flush=True)
    status = {"status": "running", "bodies": bodies, "steps": steps, "dt": DT,
              "samples_per_solver": pairs, "threads": THREADS, "warmups": 0,
              "generator": GENERATOR, "theta": THETA, "opening_angle2": 0.25,
              "G": 1, "softening": 0, "rebound_root_size": 32, "collisions": False,
              "integrators": {GROUPED: "kick-drift-kick", REBOUND: "drift-kick-drift"}}
    write_json(directory / "status.json", status)
    try:
        status["local_revision"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        status["local_dirty"] = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT))
        status["rebound_revision"] = prepare_rebound(directory)
        build(directory)
        validate_adapter(directory)
        status["two_body_validation"] = "ok"
        status["rebound_over_grouped_mean_time"] = compare(directory, bodies, steps, pairs)
        status["status"] = "ok"
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as error:
        status.update(status="failed", error=str(error))
        print(error, file=sys.stderr)
    write_json(directory / "status.json", status)
    return int(status["status"] != "ok")


if __name__ == "__main__":
    sys.exit(main())
