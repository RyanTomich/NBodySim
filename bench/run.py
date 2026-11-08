#!/usr/bin/env python3
"""Build and run timing experiments, optionally compare accuracy, and save results.

The C driver times only simulate_with_solver(); Python summarizes its samples.
"""

import argparse
from contextlib import ExitStack
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_logged(command, directory, stdout="stdout.log", stderr="stderr.log"):
    """Save process output; stderr=None combines both streams in stdout's file."""
    command = [str(part) for part in command]
    environment = {**os.environ, "LC_ALL": "C"}
    with ExitStack() as files:
        out = files.enter_context((directory / stdout).open("w"))
        err = files.enter_context((directory / stderr).open("w")) if stderr else out
        try:
            process = subprocess.run(command, cwd=ROOT, env=environment,
                                     stdin=subprocess.DEVNULL, stdout=out, stderr=err)
            return process.returncode
        except OSError as error:
            err.write(str(error) + "\n")
            return 127


def timing_summary(durations, steps):
    """Summarize whole-run averages per integration step, in seconds per step."""
    if type(steps) is not int or steps <= 0:
        raise ValueError("timing requires a positive integer step count")
    if len(durations) < 2 or any(not math.isfinite(x) or x < 0 for x in durations):
        raise ValueError("timing requires at least two finite nonnegative samples")
    per_step = [duration / steps for duration in durations]
    return {"samples": len(durations), "mean_seconds_per_step": statistics.mean(per_step),
            "variance_seconds_per_step_squared": statistics.variance(per_step),
            "stddev_seconds_per_step": statistics.stdev(per_step)}


def read_state(path):
    """Read the saved instance format, rejecting malformed/nonfinite states."""
    tokens = path.read_text().split()
    if not tokens or not re.fullmatch(r"[0-9]+", tokens[0]):
        raise ValueError("Invalid body count in {}".format(path))
    count = int(tokens[0])
    if len(tokens) != 1 + 10 * count:
        raise ValueError("Incorrect state length in {}".format(path))
    if any(len(token) > 127 or not token.isascii() or "_" in token for token in tokens):
        raise ValueError("Invalid numeric token in {}".format(path))
    try:
        # strtod in the C reader also accepts hexadecimal floating-point tokens.
        values = [float.fromhex(token) if token.lstrip("+-").lower().startswith("0x")
                  else float(token) for token in tokens[1:]]
    except OverflowError as error:
        raise ValueError("Nonfinite state in {}".format(path)) from error
    if any(not math.isfinite(value) for value in values):
        raise ValueError("Nonfinite state in {}".format(path))
    rows = [values[i:i + 10] for i in range(0, len(values), 10)]
    if any(row[0] <= 0 for row in rows):
        raise ValueError("Nonpositive mass in {}".format(path))
    return rows


def validate_bodies(actual, reference):
    if len(actual) != len(reference):
        raise ValueError("Body counts differ")
    if any(a[0] != b[0] for a, b in zip(actual, reference)):
        raise ValueError("Body masses changed")


def trajectory_errors(actual, reference):
    """Errors over ordered bodies; positions, velocities and accelerations separate."""
    validate_bodies(actual, reference)
    if any(not math.isfinite(value) for state in (actual, reference)
           for row in state for value in row):
        raise ValueError("Nonfinite state")
    rows = []
    for field, start in (("position", 1), ("velocity", 4), ("acceleration", 7)):
        error_norm = reference_norm = maximum = 0.0
        for a, b in zip(actual, reference):
            distance = math.hypot(*(a[j] - b[j] for j in range(start, start + 3)))
            error_norm = math.hypot(error_norm, distance)
            reference_norm = math.hypot(reference_norm, *b[start:start + 3])
            maximum = max(maximum, distance)
        rms = error_norm / math.sqrt(len(actual)) if actual else 0.0
        relative = error_norm / reference_norm if reference_norm else (0.0 if not error_norm else None)
        metrics = [error_norm, reference_norm, rms, maximum]
        if relative is not None:
            metrics.append(relative)
        if any(not math.isfinite(value) for value in metrics):
            raise ValueError("Nonfinite error metric for " + field)
        rows.append({"field": field, "rms_error": rms, "max_error": maximum,
                     "relative_l2_error": relative,
                     "relative_status": "ok" if relative is not None else "undefined"})
    return rows


def compare_n2(binary, instance, folder, args):
    """Compare the timed final state with one unmeasured n2 reference run."""
    directory = folder / "reference"
    directory.mkdir()
    status = {"solver": "n2", "steps": args.steps, "dt": args.dt,
              "started": now(), "status": "failed"}
    reference = None
    try:
        command = [binary, instance, str(args.dt), str(args.steps), "1", "0",
                   "reference", directory / "final.txt", "n2"]
        if run_logged(command, directory):
            raise RuntimeError("n2 reference failed; see reference/stderr.log.")
        reference = read_state(directory / "final.txt")
        validate_bodies(reference, read_state(instance))
        status["status"] = "ok"
    except (OSError, ValueError, RuntimeError) as error:
        status["error"] = str(error)
    status["finished"] = now()
    write_json(directory / "status.json", status)
    comparison = {"solver": args.solver, "reference_solver": "n2",
                  "started": now(), "status": "failed"}
    timing = folder / "timing"
    try:
        if status["status"] != "ok":
            raise RuntimeError(status["error"])
        rows = trajectory_errors(read_state(timing / "final.txt"), reference)
        write_csv(timing / "error.csv", list(rows[0]), rows)
        comparison["status"] = "ok"
    except (OSError, ValueError, RuntimeError) as error:
        comparison["error"] = str(error)
    comparison["finished"] = now()
    write_json(timing / "error-status.json", comparison)
    return comparison["status"] == "ok"


def solver_settings(binary, instance, solver, theta, order=None):
    """Ask the measured executable for its effective configuration, unmeasured."""
    with instance.open() as source:
        first = next((line.split()[0] for line in source if line.split()), "")
    command = [str(binary), "--metadata", solver, first]
    if theta is not None:
        command.extend(["theta", str(theta)])
    if order is not None:
        command.extend(["order", str(order)])
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "Could not resolve solver settings")
    return json.loads(result.stdout)


def driver_metadata(directory, status):
    path = directory / "stderr.log"
    if path.exists():
        for line in path.read_text().splitlines():
            if line.startswith(("approximate_config: ", "solver_config: ")):
                status.update(json.loads(line.split(": ", 1)[1]))
            elif line.startswith(("approximate_totals: ", "solver_totals: ")):
                status["diagnostics"] = json.loads(line.split(": ", 1)[1])


def collect_timing(binary, instance, directory, args):
    """Run warmups and timing samples, then save per-step statistics and status."""
    directory.mkdir()
    status = {"mode": "timing", "solver": args.solver, "steps": args.steps, "dt": args.dt,
              "started": now(), "status": "failed"}
    summary = None
    theta = args.theta
    try:
        status.update(solver_settings(binary, instance, args.solver, theta, args.order))
        theta = status.get("theta")
        driver = [binary, instance, str(args.dt), str(args.steps), str(args.samples),
                  str(args.warmups), "timing", directory / "final.txt", args.solver]
        if theta is not None:
            driver.extend(["theta", str(theta)])
        if "order" in status:
            driver.extend(["order", str(status["order"])])

        if run_logged(driver, directory, stdout="samples.csv"):
            raise RuntimeError("Benchmark failed; see stderr.log.")
        with (directory / "samples.csv").open(newline="") as source:
            durations = [float(row["seconds"]) for row in csv.DictReader(source)]
        if len(durations) != args.samples:
            raise ValueError("Unexpected number of timing samples.")
        summary = {"solver": args.solver, "theta": theta, "order": status.get("order"),
                   "steps": args.steps, "dt": args.dt,
                   **timing_summary(durations, args.steps)}
        write_csv(directory / "summary.csv", list(summary), [summary])
        status["status"] = "ok"
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        status["error"] = str(error)
        summary = None
    try:
        driver_metadata(directory, status)
    except (OSError, ValueError) as error:
        status.update(status="failed", error="Invalid driver metadata: " + str(error))
        summary = None
    status["finished"] = now()
    write_json(directory / "status.json", status)
    return status, summary


def create_experiment_directory(results_dir):
    try:
        revision = subprocess.run(["git", "rev-parse", "--short=7", "--verify", "HEAD^{commit}"],
                                  cwd=ROOT, capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RuntimeError(f"Cannot identify Git commit: {error}") from error
    short_hash = revision.stdout.strip()
    if revision.returncode != 0 or not re.fullmatch(r"[0-9a-f]{7,64}", short_hash):
        raise RuntimeError("Cannot identify Git commit: " + (revision.stderr.strip() or "invalid Git output"))
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%fZ")
    directory = results_dir.resolve() / f"{short_hash}_{timestamp}_timing"
    directory.mkdir(parents=True, exist_ok=False)
    return directory


def build_and_archive(directory):
    if run_logged(["make", "--no-print-directory", "-B", "all"], directory,
                  "build.stdout.log", stderr=None):
        raise RuntimeError("Build failed; see build.stdout.log.")
    binary = directory / "bench"
    shutil.copy2(ROOT / "build" / "bench", binary)
    return binary


def experiment(args):
    directory = create_experiment_directory(args.results_dir)
    summaries = []
    failed = False
    try:
        binary = build_and_archive(directory)
        for index, path in enumerate(args.inputs, 1):
            path = Path(path).resolve()
            name = f"{index:04d}-{re.sub(r'[^A-Za-z0-9_.-]', '_', path.stem)[:80]}"
            folder = directory / name
            folder.mkdir()
            try:
                instance = folder / "input.txt"
                shutil.copyfile(path, instance)
                status, summary = collect_timing(binary, instance, folder / "timing", args)
                failed |= status["status"] != "ok"
                if summary:
                    summaries.append({"instance": name, **summary})
                if args.compare_n2 and status["status"] == "ok":
                    failed |= not compare_n2(binary, instance, folder, args)
            except OSError as error:
                write_json(folder / "status.json", {"status": "failed", "error": str(error)})
                failed = True
        fields = ["instance", "solver", "theta", "order", "steps", "dt", "samples",
                  "mean_seconds_per_step", "variance_seconds_per_step_squared",
                  "stddev_seconds_per_step"]
        write_csv(directory / "timing-summary.csv", fields, summaries)
    except (OSError, ValueError, RuntimeError) as error:
        (directory / "error.txt").write_text(str(error) + "\n", encoding="utf-8")
        failed = True
    except KeyboardInterrupt:
        (directory / "error.txt").write_text("Interrupted by user.\n", encoding="utf-8")
        failed = True
    for summary in summaries:
        print("Timing {instance} [{solver}{option_label}]: {mean:.6g} ± {stddev:.6g} ms/step "
              "(mean ± sample SD; {steps} steps, {samples} samples)".format(
                  **summary, option_label="".join(
                      " {}={:g}".format(key, summary[key]) for key in ("theta", "order")
                      if summary.get(key) is not None),
                  mean=summary["mean_seconds_per_step"] * 1000,
                  stddev=summary["stddev_seconds_per_step"] * 1000))
    print(f"Results: {directory} ({'failed' if failed else 'ok'})")
    return 1 if failed else 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("timing",), help="collect runtime statistics")
    parser.add_argument("inputs", nargs="+", help="instance files (shell globs supported)")
    parser.add_argument("--dt", type=float, required=True)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--samples", type=int, default=10)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--solver", default="n2")
    parser.add_argument("--order", help="FMM expansion order (1 through 6)")
    parser.add_argument("--theta", type=float,
                        help="solver opening parameter (default supplied by the solver)")
    parser.add_argument("--compare-n2", action="store_true",
                        help="compare final trajectory state with an unmeasured n2 run")
    parser.add_argument("--results-dir", type=Path, default=ROOT / "results",
                        help="parent directory for automatically named experiments (default: results/)")
    args = parser.parse_args(argv)
    if args.solver == "direct":
        parser.error("solver direct was renamed to n2")
    if args.solver == "approximate":
        parser.error("solver approximate was renamed to barnes-hut")
    if args.theta is not None and args.order is not None:
        parser.error("--theta and --order apply to different solvers")
    if args.theta is not None:
        if not math.isfinite(args.theta) or args.theta < 0:
            parser.error("--theta must be finite and nonnegative")
    if not math.isfinite(args.dt) or args.dt <= 0 or args.steps <= 0 or args.warmups < 0:
        parser.error("dt and steps must be positive, dt finite, and warmups nonnegative")
    if args.samples < 2:
        parser.error("timing requires --samples >= 2")
    return args


def main():
    args = parse_args()
    try:
        return experiment(args)
    except (OSError, RuntimeError) as error:
        print(f"Cannot create experiment: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
