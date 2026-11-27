#!/usr/bin/env python3
"""Serial solver comparison or OpenMP scaling sweep.

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

SOLVERS = ("n2", "serial-barnes-hut", "grouped-barnes-hut", "fmm")
PARALLEL_SOLVERS = ("n2", *runner.BARNES_HUT_SOLVERS, "fmm")
PARALLEL_SERIES = tuple(
    (solver, threads) for solver in PARALLEL_SOLVERS
    for threads in ((1,) if solver == "serial-barnes-hut" else (1, 4, 8))
)


def configurations(count, parallel, solvers=None, teams=None, engaging=False):
    legacy_selection = parallel and solvers is None and teams is None and not engaging
    selected = solvers or (PARALLEL_SOLVERS if parallel else SOLVERS)
    teams = sorted(teams or ((1, 4, 8) if parallel else (1,)))
    skipped_direct = {50000: (1,), 100000: (1, 4)}.get(count, ()) if engaging else ()
    cases = []
    # Visit direct first so its smallest eligible team supplies the reference.
    for solver in PARALLEL_SOLVERS:
        if solver not in selected:
            continue
        for team in teams:
            if solver == "serial-barnes-hut" and team != 1:
                continue
            if solver == "n2" and team in skipped_direct:
                continue
            if legacy_selection and count > 50000 and solver == "n2" and team != 8:
                continue
            cases.append((solver, team))
    return cases


def check_placement(required):
    """Inspect Linux job-step affinity before the build or any timed run."""
    if not hasattr(os, "sched_getaffinity"):
        raise ValueError("Job-step affinity validation requires Linux")
    cpus = sorted(os.sched_getaffinity(0))
    cores = set()
    sockets = set()
    governors = set()
    for cpu in cpus:
        root = Path("/sys/devices/system/cpu/cpu{}".format(cpu))
        socket = (root / "topology/physical_package_id").read_text().strip()
        sockets.add(socket)
        cores.add((socket, (root / "topology/core_id").read_text().strip()))
        governor = root / "cpufreq/scaling_governor"
        if governor.is_file():
            governors.add(governor.read_text().strip())
    print("Step CPUs={} sockets={} physical cores={} governors={}".format(
        cpus, sorted(sockets), len(cores), sorted(governors) or "unavailable"), flush=True)
    if len(sockets) != 1 or len(cores) < required or len(cores) != len(cpus):
        raise ValueError("Unsuitable job-step affinity {} for {} physical cores on one socket".format(
            cpus, required))
    return {"cpus": cpus, "socket": next(iter(sockets)), "governors": sorted(governors)}


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
    integrator = getattr(args, "integrator", "leapfrog")
    warmups = getattr(args, "warmups", 0)
    status = {"status": "failed", "solver": solver, "threads": threads,
              "steps": args.steps, "dt": args.dt, "samples": args.samples, "warmups": warmups,
              "integrator": integrator,
              "started": runner.now()}
    try:
        command = [binary, instance, str(args.dt), str(args.steps), str(args.samples), str(warmups),
                   "timing", case / "final.txt", solver]
        if solver in runner.BARNES_HUT_SOLVERS:
            command.extend(["theta", str(args.theta)])
        with thread_environment(threads):
            result = runner.run_logged(command, case, stdout="samples.csv",
                                       env={"NBODY_INTEGRATOR": integrator})
        if result:
            raise RuntimeError("Driver failed; see " + str(case / "stderr.log"))
        with (case / "samples.csv").open(newline="") as source:
            samples = list(csv.DictReader(source))
        if [int(row["sample"]) for row in samples] != list(range(1, args.samples + 1)):
            raise ValueError("Missing or duplicate whole-run samples")
        durations = [float(row["seconds"]) for row in samples]
        if any(not math.isfinite(value) or value <= 0 for value in durations):
            raise ValueError("Expected finite positive whole-run timings; increase --steps if a duration is zero")
        log = (case / "stderr.log").read_text()
        status.update(check_threads(log, threads, getattr(args, "parallel", False)))
        steps = [json.loads(line.split(": ", 1)[1])
                 for line in log.splitlines()
                 if line.startswith("step_timing: ")]
        if [(row["sample"], row["step"]) for row in steps] != [
                (sample, step) for sample in range(1, args.samples + 1)
                for step in range(1, args.steps + 1)]:
            raise ValueError("Missing or duplicate step timings")
        for sample, seconds in enumerate(durations, 1):
            values = [row["seconds"] for row in steps if row["sample"] == sample]
            if (any(not math.isfinite(x) or x < 0 for x in values)
                    or math.fsum(values) > seconds + 1e-9):
                raise ValueError("Invalid individual-step timings")
        summary = summarize(durations, args.steps)
        final = runner.read_state(case / "final.txt")
        runner.validate_bodies(final, runner.read_state(instance))
        status["force_evaluations_per_run"] = runner.force_evaluations(integrator, args.steps, len(final))
        runner.driver_metadata(case, status)
        runner.write_csv(case / "steps.csv", ["sample", "step", "seconds"], steps)
        status.update(summary, status="ok")
        return summary, final
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        status["error"] = str(error)
        raise
    finally:
        status["finished"] = runner.now()
        runner.write_json(case / "status.json", status)


def sweep_status(args):
    suite = "parallel-sweep" if args.parallel else "solver-sweep"
    max_speed = 0 if args.distribution == "clusters" else 3
    integrator = getattr(args, "integrator", "leapfrog")
    cases = {count: configurations(count, args.parallel, args.solvers, args.threads, args.engaging)
             for count in args.bodies}
    settings = {"bodies": args.bodies, "steps": args.steps, "dt": args.dt,
                "integrator": integrator,
                "theta": args.theta, "order": runner.FMM_ORDER,
                "configurations": [{"bodies": count, "solver": solver, "threads": threads}
                                   for count in args.bodies
                                   for solver, threads in cases[count]],
                "threads": (args.threads or [1, 4, 8]) if args.parallel else 1,
                "samples": args.samples, "warmups": args.warmups, "seed": args.seed,
                "generator": {"distribution": args.distribution,
                              "L": 10 if args.distribution == "uniform" else None,
                              "w_min": 1, "w_max": 100, "v_max": max_speed, "a_max": 0}}
    return {"status": "running", "suite": suite, "settings": settings,
            "started": runner.now(), "error_bars": "whole_run_sample_sd"}


def generate_input(directory, count, args):
    folder = directory / str(count)
    folder.mkdir()
    instance = folder / "input.txt"
    command = [sys.executable, runner.ROOT / "generate_instance.py", "--n", str(count),
               "--distribution", args.distribution, "--L", "10",
               "--w-min", "1", "--w-max", "100", "--v-max", str(0 if args.distribution == "clusters" else 3),
               "--a-max", "0", "--seed", str(args.seed), "--output", instance]
    if runner.run_logged(command, folder, stdout="generation.stdout.log",
                         stderr="generation.stderr.log"):
        raise RuntimeError("Input generation failed; see " + str(folder))
    return instance


def result_row(args, count, solver, threads, case, summary, final, reference, reference_threads):
    distribution = {}
    errors = []
    if reference is not None:
        distribution = runner.save_distribution_errors(case, final, reference)
        errors = runner.trajectory_errors(final, reference)
        runner.write_csv(case / "error.csv", list(errors[0]), errors)
    row = {"bodies": count, "solver": solver, "threads": threads, "steps": args.steps,
           "dt": args.dt, "samples": args.samples, "integrator": args.integrator,
           "force_evaluations_per_run": runner.force_evaluations(args.integrator, args.steps, count),
           "reference_threads": reference_threads,
           "theta": args.theta if solver in runner.BARNES_HUT_SOLVERS else None,
           "order": runner.FMM_ORDER if solver == "fmm" else None, **summary, **distribution}
    if args.parallel:
        row.update(requested_threads=threads, effective_threads=threads)
    for error in errors:
        row[error["field"] + "_relative_l2"] = error["relative_l2_error"]
        row[error["field"] + "_rms"] = error["rms_error"]
    return row


def experiment(args):
    directory = runner.create_experiment_directory(args.results_dir)
    status = sweep_status(args)
    suite = status["suite"]
    runner.write_json(directory / "status.json", status)
    rows = []
    try:
        if args.check_placement:
            check_placement(max(case["threads"] for case in status["settings"]["configurations"]))
        binary = build(directory, args.parallel)
        for count in args.bodies:
            instance = generate_input(directory, count, args)
            folder = instance.parent
            reference = None
            reference_threads = None
            for solver, threads in configurations(count, args.parallel, args.solvers, args.threads, args.engaging):
                print("N={} {}: {} runs × {} steps, {} threads".format(
                    count, solver, args.samples, args.steps, threads), flush=True)
                case = folder / ("{}-{}threads".format(solver, threads) if args.parallel else solver)
                summary, final = collect(binary, instance, case, solver, args, threads)
                if solver == "n2" and reference is None:
                    reference = final
                    reference_threads = threads
                row = result_row(args, count, solver, threads, case, summary, final,
                                 reference, reference_threads)
                rows.append(row)
                fields = dict.fromkeys(key for result in rows for key in result)
                runner.write_csv(directory / (suite + ".csv"), list(fields), rows)
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
    parser.add_argument("--parallel", action="store_true",
                        help="compare OpenMP solvers at 1/4/8 threads and the serial baseline at 1")
    parser.add_argument("--engaging", action="store_true", help="use the Engaging scaling workload and large-body skips")
    parser.add_argument("--check-placement", action="store_true", help="validate Linux job-step CPU affinity")
    parser.add_argument("--threads", type=int, nargs="+", help="OpenMP team sizes")
    parser.add_argument("--solvers", nargs="+", choices=PARALLEL_SOLVERS)
    parser.add_argument("--warmups", type=int, help="discarded trajectories per case")
    parser.add_argument("--integrator", choices=runner.INTEGRATORS, default="leapfrog")
    parser.add_argument("--bodies", nargs="+", type=int, help="body counts (parallel default extends through 100000)")
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--samples", type=int, default=3, help="full runs per configuration (at least two)")
    parser.add_argument("--dt", type=float, default=0.0001)
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--distribution", choices=("uniform", "clusters"), default="uniform",
                        help="initial positions; clusters start from rest")
    parser.add_argument("--results-dir", type=Path)
    args = parser.parse_args(argv)
    if args.engaging:
        args.parallel = True
        args.threads = args.threads or [1, 4, 8, 16]
        args.solvers = args.solvers or ["n2", "grouped-barnes-hut", "fmm"]
        args.bodies = args.bodies or [1000, 5000, 10000, 50000, 100000]
    if args.warmups is None:
        args.warmups = int(args.engaging)
    if args.warmups < 0:
        parser.error("warmups must be nonnegative")
    if args.threads:
        if not args.parallel or min(args.threads) <= 0 or len(set(args.threads)) != len(args.threads):
            parser.error("--threads requires a parallel sweep and unique positive team sizes")
        args.threads.sort()
    if args.solvers and (len(set(args.solvers)) != len(args.solvers)
                         or (not args.parallel and "parallel-barnes-hut" in args.solvers)):
        parser.error("duplicate or incompatible solver selection")
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
    for count in args.bodies:
        if not configurations(count, args.parallel, args.solvers, args.threads, args.engaging):
            parser.error("no eligible solver/thread configurations at {} bodies".format(count))
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
