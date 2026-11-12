#!/usr/bin/env python3
"""Profile Barnes-Hut with one or more thread counts using the simulation timing driver."""

import argparse
from contextlib import contextmanager
import csv
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess

import run as runner

PHASES = {
    "membership_seconds": "Reuse membership checks",
    "construction_seconds": "Tree construction",
    "summaries_seconds": "Mass summaries",
    "bounds_seconds": "Target leaf bounds",
    "traversal_seconds": "Traversal and opening-distance checks",
    "force_seconds": "Force distances and accumulation",
    "evaluation_other_seconds": "Scheduling, idle time, and timer overhead",
    "simulation_other_seconds": "Integration and other overhead",
}
COUNTERS = ("calls", "membership_checks", "tree_builds", "tree_reuses",
            "source_node_visits", "direct_evaluations", "aggregate_evaluations", "worker_slots")
TIMERS = tuple(key for key in PHASES if not key.endswith("other_seconds"))


WORKER_TIMERS = ("bounds_worker_seconds", "traversal_worker_seconds", "force_worker_seconds")
PHASE_TIME_BASIS = "serial_wall_and_mean_worker_elapsed"


def phase_samples(folder, steps, samples, threads=None):
    snapshots = {}
    for line in (folder / "stderr.log").read_text().splitlines():
        if line.startswith("phase_sample: "):
            value = json.loads(line.split(": ", 1)[1])
            key = (value.pop("sample"), value.pop("boundary"))
            if key in snapshots:
                raise ValueError("Duplicate phase snapshot")
            snapshots[key] = value
    with (folder / "samples.csv").open() as source:
        durations = list(csv.DictReader(source))
    if len(durations) != samples or len(snapshots) != 2 * samples:
        raise ValueError("Wrong sample/snapshot count")
    rows = []
    for index, duration in enumerate(durations, 1):
        before, after = snapshots[index, "before"], snapshots[index, "after"]
        # Old serial archives have no worker counters; interpret them as 1T.
        legacy = "worker_slots" not in after
        counters = tuple(key for key in COUNTERS if not legacy or key != "worker_slots")
        delta = {key: after[key] - before[key]
                 for key in (*counters, *TIMERS, "evaluation_seconds")}
        actual = after.get("min_threads", 1)
        if actual < 1 or actual != after.get("max_threads", 1):
            raise ValueError("Runtime changed the OpenMP team size")
        if threads is not None and actual != threads:
            raise ValueError("Requested {} threads but measured {}".format(threads, actual))
        if legacy:
            delta["worker_slots"] = delta["calls"]
            delta.update({key: delta[key.replace("_worker", "")] for key in WORKER_TIMERS})
        else:
            delta.update({key: after[key] - before[key] for key in WORKER_TIMERS})
        if delta["worker_slots"] != actual * delta["calls"]:
            raise ValueError("Unexpected team size or missing profiled force calls")
        for key in WORKER_TIMERS:
            if not math.isclose(delta[key] / actual, delta[key.replace("_worker", "")],
                                rel_tol=1e-8, abs_tol=1e-9):
                raise ValueError("Worker time was not normalized by team size")
        if any(not math.isfinite(x) or x < 0 for x in delta.values()):
            raise ValueError("Invalid phase measurement")
        if delta["calls"] != steps + 1:
            raise ValueError("Expected initial acceleration plus one force call per step")
        if delta["tree_builds"] + delta["tree_reuses"] != delta["calls"]:
            raise ValueError("Missing topology outcome")
        seconds = float(duration["seconds"])
        if int(duration["sample"]) != index or not math.isfinite(seconds) or seconds <= 0:
            raise ValueError("Invalid or unresolvable duration; increase --steps")
        delta["evaluation_other_seconds"] = delta["evaluation_seconds"] - sum(
            delta[key] for key in ("bounds_seconds", "traversal_seconds", "force_seconds"))
        delta["simulation_other_seconds"] = seconds - sum(delta[key] for key in (
            "membership_seconds", "construction_seconds", "summaries_seconds", "evaluation_seconds"))
        # Allow sub-nanosecond rounding from subtracting cumulative floating totals.
        if min(delta["evaluation_other_seconds"], delta["simulation_other_seconds"]) < -1e-9:
            raise ValueError("Phase times exceed enclosing wall time")
        rows.append({"sample": index, "threads": actual, "seconds": seconds, **delta})
    return rows


def summarize(rows, steps):
    total = statistics.mean(row["seconds"] for row in rows)
    summaries = []
    for key, label in PHASES.items():
        values = [row[key] * 1000 / steps for row in rows]
        mean, sd = statistics.mean(values), statistics.stdev(values)
        summaries.append({"phase": label, "field": key, "samples": len(rows),
                          "mean_ms_per_step": mean, "sd_ms_per_step": sd,
                          "approx_95ci_half_width_ms": 1.96 * sd / math.sqrt(len(rows)),
                          "percent_of_total": 100 * statistics.mean(row[key] for row in rows) / total})
    return sorted(summaries, key=lambda row: -row["mean_ms_per_step"])


def build(directory, profile, openmp):
    directory.mkdir()
    # All thread counts share the same two executables and compiler/LTO flags.
    # OMPFLAGS/OMPLIBS are inherited (Apple Clang requires libomp overrides).
    flags = os.environ.get("CPPFLAGS", "")
    if "NBODY_PROFILE_BH" in flags:
        raise ValueError("Let this runner manage NBODY_PROFILE_BH")
    if profile:
        flags += " -DNBODY_PROFILE_BH"
    options = ["OPENMP=" + str(int(openmp)), "CPPFLAGS=" + flags]
    with (directory / "build.stdout.log").open("w") as log:
        for command in (["make", "--no-print-directory", "print-config", *options],
                        ["make", "--no-print-directory", "-B", "all", *options]):
            subprocess.run(command, cwd=runner.ROOT, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
    binary = directory / "bench"
    shutil.copy2(runner.ROOT / "build/bench", binary)
    return binary


@contextmanager
def thread_environment(threads):
    """Apply the same fixed team settings to metadata queries and driver runs."""
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


def measure_case(folder, instance, baseline, profiled, args, threads):
    folder.mkdir()
    initial = runner.read_state(instance)
    timings = {}
    for name, binary in (("baseline-before", baseline), ("profile", profiled),
                         ("baseline-after", baseline)):
        print("{} [{}T]: {}".format(instance.parent.name, threads, name), flush=True)
        measured, summary = runner.collect_timing(binary, instance, folder / name, args)
        if measured["status"] != "ok":
            raise RuntimeError(str(measured))
        runner.validate_bodies(runner.read_state(folder / name / "final.txt"), initial)
        timings[name] = summary["mean_seconds_per_step"]
    reference = (folder / "baseline-before/final.txt").read_bytes()
    if any((folder / name / "final.txt").read_bytes() != reference
           for name in ("profile", "baseline-after")):
        raise ValueError("Instrumentation changed the saved trajectory")
    rows = phase_samples(folder / "profile", args.steps, args.samples, threads)
    runner.write_csv(folder / "phase-samples.csv", list(rows[0]), rows)
    summary = [{"instance": instance.parent.name, "bodies": len(initial),
                "threads": threads, "steps": args.steps, "dt": args.dt, "theta": args.theta,
                "phase_time_basis": PHASE_TIME_BASIS, **row}
               for row in summarize(rows, args.steps)]
    runner.write_csv(folder / "phase-summary.csv", list(summary[0]), summary)
    first, second = summary[:2]
    gaps = [(row[first["field"]] - row[second["field"]]) * 1000 / args.steps
            for row in rows]
    baseline_mean = (timings["baseline-before"] + timings["baseline-after"]) / 2
    if min(timings.values()) <= 0:
        raise ValueError("Unresolvable timing; increase --steps")
    ratio = (first["mean_ms_per_step"] / second["mean_ms_per_step"]
             if second["mean_ms_per_step"] > 0 else None)
    checks = sum(row["membership_checks"] for row in rows)
    builds = sum(row["tree_builds"] for row in rows)
    reuses = sum(row["tree_reuses"] for row in rows)
    details = {"status": "ok", "bodies": len(initial), "threads": threads,
               "phase_time_basis": PHASE_TIME_BASIS, "bitwise_equal_final_states": True,
               "timings_seconds_per_step": timings,
               "profile_vs_baseline_percent": 100 * (timings["profile"] / baseline_mean - 1),
               "baseline_after_vs_before_percent": 100 * (timings["baseline-after"] / timings["baseline-before"] - 1),
               "measured_force_calls": sum(row["calls"] for row in rows),
               "membership_checks": checks, "tree_builds": builds, "tree_reuses": reuses,
               "mean_ms_per_check": 1000 * sum(row["membership_seconds"] for row in rows) / checks if checks else None,
               "mean_ms_per_build": 1000 * sum(row["construction_seconds"] for row in rows) / builds if builds else None,
               "largest_phase": first["phase"], "second_phase": second["phase"],
               "largest_to_second_ratio": ratio,
               "paired_gap_ms_per_step": statistics.mean(gaps),
               "paired_gap_approx_95ci_half_width_ms": 1.96 * statistics.stdev(gaps) / math.sqrt(len(gaps))}
    runner.write_json(folder / "status.json", details)
    print("N={} [{}T]: {} {:.2f}%, {} {:.2f}%; ratio {}; rebuilds {}/{}".format(
        len(initial), threads, first["phase"], first["percent_of_total"], second["phase"],
        second["percent_of_total"], "{:.2f}x".format(ratio) if ratio is not None else "unresolved",
        builds, builds + reuses), flush=True)
    return summary, details


def experiment(args):
    directory = runner.create_experiment_directory(args.results_dir)
    print("Results: {}".format(directory), flush=True)
    status = {"status": "failed", "started": runner.now(), "threads": args.threads,
              "openmp": not args.serial, "platform": platform.platform(), "machine": platform.machine(),
              "dt": args.dt, "steps": args.steps, "samples": args.samples,
              "warmups": args.warmups, "theta": args.theta, "phase_time_basis": PHASE_TIME_BASIS,
              "normalization": "whole simulation / steps; includes initial force evaluation",
              "parallel_phase_attribution": "sum of worker elapsed times / actual team size; residual includes idle and synchronization"}
    try:
        baseline = build(directory / "baseline-build", False, not args.serial)
        profiled = build(directory / "profile-build", True, not args.serial)
        all_summary, speedups = [], []
        for index, path in enumerate(args.inputs, 1):
            folder = directory / ("{:04d}-".format(index) + Path(path).stem)
            folder.mkdir()
            instance = folder / "input.txt"
            shutil.copyfile(path, instance)
            initial = runner.read_state(instance)
            if not initial:
                raise ValueError("Phase profiling needs a nonempty input")
            reference = None
            cases = {}
            for threads in args.threads:
                case = folder / ("threads-" + str(threads))
                with thread_environment(threads):
                    summary, details = measure_case(case, instance, baseline, profiled, args, threads)
                final = (case / "profile/final.txt").read_bytes()
                if reference is not None and final != reference:
                    raise ValueError("Thread counts changed the saved Barnes-Hut trajectory")
                reference = final
                all_summary.extend(summary)
                cases[threads] = details
            reference_threads = min(args.threads)
            def baseline_time(details):
                times = details["timings_seconds_per_step"]
                return (times["baseline-before"] + times["baseline-after"]) / 2
            reference_time = baseline_time(cases[reference_threads])
            for threads, details in cases.items():
                speedups.append({"instance": folder.name, "bodies": len(initial), "threads": threads,
                                 "reference_threads": reference_threads,
                                 "baseline_ms_per_step": 1000 * baseline_time(details),
                                 "speedup": reference_time / baseline_time(details)})
            runner.write_json(folder / "status.json", {"status": "ok", "threads": args.threads,
                              "bitwise_equal_across_threads": True})
        runner.write_csv(directory / "phase-summary.csv", list(all_summary[0]), all_summary)
        runner.write_csv(directory / "thread-summary.csv", list(speedups[0]), speedups)
        status["status"] = "ok"
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.CalledProcessError) as error:
        status["error"] = str(error)
        (directory / "error.txt").write_text(str(error) + "\n")
    finally:
        status["finished"] = runner.now()
        runner.write_json(directory / "status.json", status)
    return 0 if status["status"] == "ok" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--samples", type=int, default=50)
    parser.add_argument("--warmups", type=int, default=3)
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4],
                        help="OpenMP team sizes to compare (default: 1 4)")
    parser.add_argument("--serial", action="store_true",
                        help="build without OpenMP; requires --threads 1")
    parser.add_argument("--results-dir", type=Path, default=runner.ROOT / "results/bh-phases")
    args = parser.parse_args()
    if (not math.isfinite(args.dt) or args.dt <= 0 or args.steps < 1 or args.samples < 2
            or args.warmups < 0 or not math.isfinite(args.theta) or args.theta < 0):
        parser.error("Need positive finite dt, positive steps, >=2 samples, >=0 warmups and finite nonnegative theta")
    if (any(value < 1 for value in args.threads) or len(set(args.threads)) != len(args.threads)
            or (args.serial and args.threads != [1])):
        parser.error("threads must be unique positive counts; --serial requires --threads 1")
    args.threads.sort()
    args.solver, args.order = "barnes-hut", None
    raise SystemExit(experiment(args))


if __name__ == "__main__":
    main()
