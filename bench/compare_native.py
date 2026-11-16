#!/usr/bin/env python3
"""Compare native historical Barnes-Hut or current n2/src/Rushil implementations."""

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import random
import shlex
import statistics
import subprocess
import sys
from tempfile import TemporaryDirectory

from run import read_state, write_csv, write_json

ROOT = Path(__file__).resolve().parents[1]
NAMES = ("original-src", "merged", "rushil")
CURRENT_NAMES = ("n2", "rushil", "src")


def configurations_for(args):
    names = CURRENT_NAMES if args.suite == "current" else NAMES
    return [(n, name, t) for n in args.bodies for name in names
            for t in ([1] if name == "n2" else args.threads)
            if not (getattr(args, "skip_large_n2", False) and name == "n2" and n == max(args.bodies))]


def samples_for(args, key):
    count, name, _ = key
    if (args.suite == "current" and name == "n2" and count == max(args.bodies)
            and args.large_n2_samples is not None):
        return args.large_n2_samples
    return args.batch_size


def resume_current(parser, requested):
    """Reuse saved binaries and successful points; retain interrupted attempts."""
    out = Path(requested.resume).resolve()
    try:
        status = json.loads((out / "status.json").read_text())
        if status["settings"]["suite"] != "current" or status["status"] == "ok":
            raise ValueError("Resume requires an incomplete current-suite run")
        if json.loads((out / "validation/status.json").read_text())["status"] != "ok":
            raise ValueError("Resume requires successful adapter validation")
        args = argparse.Namespace(**status["settings"])
        args.resume = str(out)
        args.large_n2_samples = (requested.large_n2_samples if requested.large_n2_samples is not None
                                 else getattr(args, "large_n2_samples", None))
        args.skip_large_n2 = requested.skip_large_n2 or getattr(args, "skip_large_n2", False)
        if args.large_n2_samples is not None and not 2 <= args.large_n2_samples <= args.samples:
            raise ValueError("--large-n2-samples must be between 2 and --samples")
        if not all((out / name).is_file() for name in CURRENT_NAMES):
            raise ValueError("Saved comparison executables are missing")
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    status.setdefault("interruptions", []).append({"finished": status.get("finished"),
                                                  "error": status.pop("error", "interrupted")})
    status.update(status="running", settings=vars(args), completed_batches=0)
    write_json(out / "status.json", status)
    try:
        execute(args, out, status)
        status["status"] = "ok"
    except Exception as error:
        status.update(status="failed", error=str(error))
        print(str(error), file=sys.stderr)
    finally:
        status["finished"] = datetime.now(timezone.utc).isoformat()
        write_json(out / "status.json", status)
    return 0 if status["status"] == "ok" else 1


def build_current(out, flags, libs):
    driver = ROOT / "bench/comparison/native_driver.c"
    source = ROOT / "src"
    sources = sorted(p for p in source.glob("*.c") if p.name != "main.c")
    for name in ("n2", "src"):
        defines = ["-DN2_BACKEND"] if name == "n2" else []
        logged([*flags, *defines, "-I" + str(source), driver, *sources,
                *libs, "-lm", "-o", out / name], out / ("build-" + name + ".log"))
    rushil = ROOT / "rushil"
    # The public state views used by I/O are identical; allocation and native
    # stepping still belong to Rushil, including its ten-array representation.
    logged([*flags, '-DRUSHIL_SOURCE="' + str(rushil / "barnes_hut_nbody.c") + '"',
            "-I" + str(source), driver, rushil / "state.c", rushil / "physics.c",
            source / "io.c", source / "parse.c", *libs, "-lm", "-o", out / "rushil"],
           out / "build-rushil.log")


def logged(command, path, env=None, expected_success=True):
    with path.open("w") as out, path.with_suffix(path.suffix + ".stderr").open("w") as err:
        process = subprocess.run([str(p) for p in command], cwd=ROOT, env=env,
                                 stdout=out, stderr=err)
    if expected_success and process.returncode:
        raise RuntimeError("Command failed; see " + str(path))
    return process.returncode


def environment(threads):
    # Prevent inherited team limits or affinity policies from silently changing
    # the requested comparison. Let the OS schedule the heterogeneous CPU.
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("OMP_", "KMP_", "GOMP_"))}
    env.update(OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE", LC_ALL="C")
    return env


def build(out, revision, flags, libs):
    driver = ROOT / "bench/comparison/native_driver.c"
    with TemporaryDirectory(prefix="nbody-original-") as directory:
        original = Path(directory)
        archive = subprocess.check_output(["git", "archive", revision, "src", "rushil"], cwd=ROOT)
        subprocess.run(["tar", "-x", "-C", str(original)], input=archive, check=True)
        for name, source in (("original-src", original / "src"), ("merged", ROOT / "src")):
            sources = sorted(p for p in source.glob("*.c") if p.name != "main.c")
            defines = ["-DORIGINAL_SRC"] if name == "original-src" else []
            logged([*flags, *defines, "-I" + str(source), driver, *sources,
                    *libs, "-lm", "-o", out / name], out / ("build-" + name + ".log"))
        source = original / "rushil"
        # Shared I/O only uses the identical public state views; Rushil retains
        # its own ten-array state allocator and all native force/integration code.
        logged([*flags, '-DRUSHIL_SOURCE="' + str(source / "barnes_hut_nbody.c") + '"',
                "-I" + str(original / "src"), driver, source / "state.c", source / "physics.c",
                original / "src/io.c", original / "src/parse.c", *libs, "-lm", "-o", out / "rushil"],
               out / "build-rushil.log")


def scalar_reference(rows, dt, steps, leapfrog):
    state = [row[:] for row in rows]

    def forces():
        for i, target in enumerate(state):
            total = [0.0] * 3
            for j, source in enumerate(state):
                if i == j:
                    continue
                d = [source[k] - target[k] for k in range(1, 4)]
                distance = math.hypot(*d)
                if distance:
                    for axis in range(3):
                        total[axis] += source[0] / distance / distance * (d[axis] / distance)
            target[7:10] = total

    if leapfrog:
        forces()
    for _ in range(steps):
        for row in state:
            for axis in range(3):
                if leapfrog:
                    row[4 + axis] += 0.5 * dt * row[7 + axis]
                row[1 + axis] += dt * row[4 + axis]
                if not leapfrog:
                    row[4 + axis] += dt * row[7 + axis]
        forces()
        if leapfrog:
            for row in state:
                for axis in range(3):
                    row[4 + axis] += 0.5 * dt * row[7 + axis]
    return state


def metadata(path):
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith("threads: "):
            result.update(json.loads(line.split(": ", 1)[1]))
        elif line.startswith("force_error: "):
            result["force_error"] = json.loads(line.split(": ", 1)[1])
        elif line.startswith("rushil_task_cutoff: "):
            result["rushil_task_cutoff"] = int(line.split(": ", 1)[1])
    return result


def validate_adapters(out, args):
    folder = out / "validation"
    folder.mkdir()
    fixture = ROOT / "tests/data/two_body.txt"
    original = read_state(fixture)
    names = CURRENT_NAMES if args.suite == "current" else NAMES
    for name in names:
        for threads in ([1] if name == "n2" else args.threads):
            label = "{}-{}".format(name, threads)
            final = folder / (label + "-final.txt")
            log = folder / (label + ".csv")
            logged([out / name, fixture, 0.01, 3, 2, 1, 0.5, final], log, environment(threads))
            actual = read_state(final)
            expected = scalar_reference(original, 0.01, 3, name != "original-src")
            if len(actual) != len(expected) or not all(
                    math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-11)
                    for row_a, row_b in zip(actual, expected) for a, b in zip(row_a, row_b)):
                raise RuntimeError("Scalar integration regression failed: " + label)
            info = metadata(log.with_suffix(".csv.stderr"))
            if info.get("effective") != (1 if name == "original-src" else threads):
                raise RuntimeError("Unexpected team size: " + label)
            if info["force_error"]["relative_l2"] > 1e-11:
                raise RuntimeError("Direct-force regression failed: " + label)
        bad = folder / (name + "-invalid.log")
        if logged([out / name, fixture, "nan", 3, 2, 1, 0.5, folder / "invalid.txt"],
                  bad, environment(1), expected_success=False) == 0:
            raise RuntimeError("Invalid timestep accepted: " + name)
    write_json(folder / "status.json", {"status": "ok", "tests": "scalar Euler/leapfrog, forces, teams, invalid dt"})


def summarize(durations, steps):
    per_step = [seconds * 1000 / steps for seconds in durations]
    mean = statistics.mean(per_step)
    sd = statistics.stdev(per_step)
    return {"samples": len(per_step), "mean_ms_per_step": mean, "stddev_ms_per_step": sd,
            "median_ms_per_step": statistics.median(per_step),
            "relative_standard_error": sd / math.sqrt(len(per_step)) / mean}


def execute(args, out, status):
    configurations = configurations_for(args)
    status["skipped_configurations"] = ([{"bodies": max(args.bodies), "implementation": "n2",
                                          "threads": 1, "reason": "explicitly skipped"}]
                                        if args.skip_large_n2 else [])
    binaries = {name: out / name for _, name, _ in configurations}
    inputs = {}
    for count in args.bodies:
        folder = out / str(count)
        folder.mkdir(exist_ok=bool(args.resume))
        inputs[count] = folder / "input.txt"
        if args.resume:
            continue
        logged([sys.executable, ROOT / "generate_instance.py", "--n", count, "--L", 1,
                "--w-min", 1 / count, "--w-max", 1 / count, "--v-max", 0, "--a-max", 0,
                "--seed", 42, "--output", inputs[count]], folder / "generate.log")
    samples = {key: [] for key in configurations}
    information = {}
    baselines = {}
    batch_counts = {key: 0 for key in configurations}
    rows = []

    def save_summary():
        rows.clear()
        for key in configurations:
            if len(samples[key]) < 2:
                continue
            count, name, threads = key
            summary = summarize(samples[key], args.steps)
            info = information[key]
            row = {"bodies": count, "implementation": name, "requested_threads": threads,
                   "effective_threads": info["effective"], "integrator": "euler" if name == "original-src" else "leapfrog",
                   "steps": args.steps, "dt": args.dt, "theta": "" if name == "n2" else args.theta,
                   "force_evaluations_per_run": args.steps + (name != "original-src"),
                   **summary, "stable": summary["relative_standard_error"] <= args.rse,
                   "force_relative_l2": info["force_error"]["relative_l2"],
                   "force_rms_absolute": info["force_error"]["rms_absolute"],
                   "force_targets": info["force_error"]["targets"]}
            rows.append(row)
        if rows:
            write_csv(out / "comparison.csv", list(rows[0]), rows)

    def run_batches(keys, rounds):
        schedule = [key for key in keys for _ in range(rounds)]
        random.Random(42 + status["completed_batches"]).shuffle(schedule)
        for key in schedule:
            count, name, threads = key
            case = out / str(count) / (name + "-" + str(threads) + "threads")
            case.mkdir(exist_ok=True)
            batch_counts[key] += 1
            batch = case / ("batch-{:02d}".format(batch_counts[key]))
            reused = False
            if args.resume and batch.exists():
                saved_status = json.loads((batch / "status.json").read_text())
                reused = saved_status.get("status") == "ok"
                if not reused:
                    batch.rename(case / ("interrupted-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")))
            target_samples = samples_for(args, key)
            if not reused:
                batch.mkdir()
                write_json(batch / "status.json", {"status": "running"})
            try:
                if not reused:
                    logged([binaries[name], inputs[count], args.dt, args.steps, target_samples,
                            args.warmups, args.theta, batch / "final.txt"], batch / "samples.csv", environment(threads))
                with (batch / "samples.csv").open() as source:
                    seconds = [float(row["seconds"]) for row in csv.DictReader(source)]
                if len(seconds) != target_samples or any(not math.isfinite(s) or s <= 0 for s in seconds):
                    raise ValueError("Invalid sample data")
                info = metadata(batch / "samples.csv.stderr")
                if info.get("requested") != threads or info.get("effective") != (1 if name == "original-src" else threads):
                    raise ValueError("Unexpected effective thread count")
                if name == "rushil" and info.get("rushil_task_cutoff") != 256:
                    raise ValueError("Rushil default task cutoff changed")
                actual = read_state(batch / "final.txt")
                initial = read_state(inputs[count])
                if len(actual) != count or any(a[0] != b[0] for a, b in zip(actual, initial)):
                    raise ValueError("Body count/masses changed")
                baseline_key = (count, name)
                if baseline_key in baselines:
                    if actual != baselines[baseline_key]:
                        raise ValueError("Thread settings or batches changed final state")
                else:
                    baselines[baseline_key] = actual
                if not all(math.isfinite(info["force_error"][k]) for k in ("relative_l2", "rms_absolute")):
                    raise ValueError("Invalid force error")
                information[key] = info
                samples[key].extend(seconds)
                write_json(batch / "status.json", {"status": "ok", **info, "samples": len(seconds)})
            except Exception as error:
                write_json(batch / "status.json", {"status": "failed", "error": str(error)})
                raise
            status["completed_batches"] += 1
            status["latest"] = {"bodies": count, "implementation": name, "threads": threads,
                                "samples": len(samples[key])}
            write_json(out / "status.json", status)
            save_summary()
            print("Batch {}: {} bodies, {}, {} threads, {} samples".format(
                status["completed_batches"], count, name, threads, len(samples[key])), flush=True)

    run_batches(configurations, args.samples // args.batch_size)
    extend = [key for key in configurations
              if summarize(samples[key], args.steps)["relative_standard_error"] > args.rse]
    status["extended_configurations"] = ([list(key) for key in extend]
                                         if args.max_samples > args.samples else [])
    if extend and args.max_samples > args.samples:
        run_batches(extend, (args.max_samples - args.samples) // args.batch_size)
    status["unstable_configurations"] = [list(key) for key in configurations
        if summarize(samples[key], args.steps)["relative_standard_error"] > args.rse]
    status["configurations"] = len(configurations)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=("historical", "current"), default="historical",
                        help="current: working-tree n2 (1 thread), Rushil and src (1/4 threads)")
    parser.add_argument("--revision", default="25f7df2")
    parser.add_argument("--large-n2-samples", type=int,
                        help="current suite: fewer samples only for n2 at the largest body count (minimum 2)")
    parser.add_argument("--skip-large-n2", action="store_true",
                        help="current suite: explicitly omit n2 at the largest body count")
    parser.add_argument("--resume", help="resume an incomplete current run using its saved settings and binaries; "
                        "only --large-n2-samples and --skip-large-n2 override its settings")
    parser.add_argument("--bodies", type=int, nargs="+", default=None)
    parser.add_argument("--threads", type=int, nargs="+", default=None)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--dt", type=float, default=0.0001)
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--samples", type=int, default=None)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--rse", type=float, default=0.02)
    parser.add_argument("--machine-description", default="")
    args = parser.parse_args()
    if args.resume:
        return resume_current(parser, args)
    current = args.suite == "current"
    if args.skip_large_n2 and not current:
        parser.error("--skip-large-n2 requires current suite")
    if args.bodies is None:
        args.bodies = [1000, 10000, 100000] if current else [10000, 20000, 50000]
    if args.threads is None:
        args.threads = [1, 4] if current else [1, 4, 8]
    if args.samples is None:
        args.samples = 10 if current else 30
    if args.max_samples is None:
        args.max_samples = args.samples if current else 60
    if args.batch_size is None:
        args.batch_size = args.samples if current else 5
    if args.large_n2_samples is not None and (not current or not 2 <= args.large_n2_samples <= args.samples):
        parser.error("--large-n2-samples requires current suite and a value between 2 and --samples")
    if current and (args.threads != [1, 4] or args.max_samples != args.samples
                    or args.batch_size != args.samples):
        parser.error("current suite requires --threads 1 4 and one fixed-size batch per configuration")
    if (min(args.bodies) < 2 or min(args.threads) < 1 or args.steps < 1 or args.samples < 2
        or args.batch_size < 2 or args.max_samples < args.samples or args.warmups < 0
        or args.samples % args.batch_size or args.max_samples % args.batch_size
        or any(not math.isfinite(v) or v <= 0 for v in (args.dt, args.theta, args.rse))
        or len(set(args.bodies)) != len(args.bodies) or len(set(args.threads)) != len(args.threads)):
        parser.error("Invalid benchmark settings")
    revision = (None if current else subprocess.check_output(
        ["git", "rev-parse", args.revision + "^{commit}"], cwd=ROOT, text=True).strip())
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    out = ROOT / "results" / (("current-comparison-" if current else "native-barnes-hut-") + stamp)
    out.mkdir(parents=True)
    compiler = shlex.split(os.environ.get("CC", "clang"))
    omp = shlex.split(os.environ.get("OMPFLAGS", "-fopenmp"))
    libs = shlex.split(os.environ.get("OMPLIBS", "-fopenmp"))
    flags = [*compiler, "-std=c17", "-Wall", "-Wextra", "-Wpedantic", "-O3", "-g",
             "-fno-omit-frame-pointer", "-fno-fast-math", "-UNDEBUG", "-flto", *omp]
    status = {"status": "running", "settings": vars(args), "original_revision": revision,
              "merged_base_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "merged_working_tree_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
              "compiler": subprocess.check_output([*compiler, "--version"], text=True),
              "compiler_flags": flags, "link_flags": libs, "platform": platform.platform(),
              "logical_cpus": os.cpu_count(), "machine_description": args.machine_description,
              "scope": "native integration and forces; fresh workspace per run; no allocation/reset/I/O",
              "completed_batches": 0}
    write_json(out / "status.json", status)
    print("Results: " + str(out), flush=True)
    try:
        if current:
            build_current(out, flags, libs)
        else:
            build(out, revision, flags, libs)
        validate_adapters(out, args)
        execute(args, out, status)
        status["status"] = "ok"
    except Exception as error:
        status.update(status="failed", error=str(error))
        print(str(error), file=sys.stderr)
    finally:
        status["finished"] = datetime.now(timezone.utc).isoformat()
        write_json(out / "status.json", status)
    return 0 if status["status"] == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
