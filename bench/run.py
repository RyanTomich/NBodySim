#!/usr/bin/env python3
"""Measure simulate_with_solver() and save results with the exact benchmark executable."""

import argparse
from contextlib import ExitStack
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVENTS = "cycles:u,instructions:u,branches:u,branch-misses:u,cache-references:u,cache-misses:u"
COUNTER_FIELDS = ["event", "value", "unit", "status", "event_runtime_ns", "running_percent"]


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


def compare_direct(binary, instance, folder, modes, args):
    """One ungated reference call per input, after all measurement modes finish."""
    directory = folder / "reference"
    directory.mkdir()
    status = {"solver": "direct", "steps": args.steps, "dt": args.dt,
              "started": now(), "status": "failed"}
    reference = None
    try:
        command = [binary, instance, str(args.dt), str(args.steps), "1", "0",
                   "reference", directory / "final.txt", "direct"]
        if run_logged(command, directory):
            raise RuntimeError("Direct reference failed; see reference/stderr.log.")
        reference = read_state(directory / "final.txt")
        validate_bodies(reference, read_state(instance))
        status["status"] = "ok"
    except (OSError, ValueError, RuntimeError) as error:
        status["error"] = str(error)
    status["finished"] = now()
    write_json(directory / "status.json", status)
    failed = status["status"] != "ok"
    for mode in modes:
        comparison = {"solver": args.solver, "reference_solver": "direct",
                      "started": now(), "status": "failed"}
        try:
            if status["status"] != "ok":
                raise RuntimeError(status["error"])
            rows = trajectory_errors(read_state(folder / mode / "final.txt"), reference)
            write_csv(folder / mode / "error.csv", list(rows[0]), rows)
            comparison["status"] = "ok"
        except (OSError, ValueError, RuntimeError) as error:
            comparison["error"] = str(error)
            failed = True
        comparison["finished"] = now()
        write_json(folder / mode / "error-status.json", comparison)
    return not failed


def parse_counters(text):
    """Parse perf stat -x ';' without mistaking process elapsed time for a counter."""
    rows = []
    for fields in csv.reader(text.splitlines(), delimiter=";"):
        if not fields or fields[0].lstrip().startswith("#"):
            continue
        fields = [field.strip() for field in fields]
        if len(fields) < 3 or not fields[2]:
            continue
        if any("seconds time elapsed" in field or field in ("seconds user", "seconds sys")
               for field in fields):
            continue
        raw, unit, event = fields[:3]
        try:
            value = float(raw)
            if not math.isfinite(value) or value < 0:
                raise ValueError
            status = "ok"
        except ValueError:
            value = None
            status = raw or "missing"
        rows.append({"event": event, "value": value, "unit": unit, "status": status,
                     "event_runtime_ns": fields[3] if len(fields) > 3 else "",
                     "running_percent": fields[4] if len(fields) > 4 else ""})
    return rows


def counter_ratios(rows):
    # Avoid merging different PMUs or deriving ratios from duplicate event rows.
    def value(name):
        matching = [row for row in rows if row["event"].split(":")[0] == name]
        if len(matching) == 1 and matching[0]["status"] == "ok":
            return matching[0]["value"]
        return None

    ratios = []
    for name, numerator, denominator in (
        ("ipc", "instructions", "cycles"),
        ("branch_miss_rate", "branch-misses", "branches"),
        ("cache_miss_rate", "cache-misses", "cache-references"),
    ):
        a, b = value(numerator), value(denominator)
        result = a / b if a is not None and b is not None and b > 0 else None
        ratios.append({"metric": name, "value": result,
                       "unit": "instructions/cycle" if name == "ipc" else "fraction",
                       "status": "ok" if result is not None else "unavailable"})
    return ratios


def perf_available():
    if platform.system() != "Linux":
        return "Linux perf is required for this mode; portable timing remains available."
    if not shutil.which("perf"):
        return "perf was not found. Install the Linux perf tools for this host's kernel."
    return None


def collect(mode, binary, instance, directory, args):
    directory.mkdir()
    status = {"mode": mode, "solver": args.solver, "steps": args.steps, "dt": args.dt,
              "started": now(), "status": "failed"}
    summary = None
    try:
        if mode != "timing":
            unavailable = perf_available()
            if unavailable:
                raise RuntimeError(unavailable)
        samples = args.samples if mode == "timing" else 1
        driver = [binary, instance, str(args.dt), str(args.steps), str(samples),
                  str(args.warmups), "timing" if mode == "timing" else "perf",
                  directory / "final.txt", args.solver]

        if mode == "timing":
            if run_logged(driver, directory, stdout="samples.csv"):
                raise RuntimeError("Benchmark failed; see stderr.log.")
            with (directory / "samples.csv").open(newline="") as source:
                durations = [float(row["seconds"]) for row in csv.DictReader(source)]
            if len(durations) != samples:
                raise ValueError("Unexpected number of timing samples.")
            summary = {"solver": args.solver, "steps": args.steps, "dt": args.dt,
                       **timing_summary(durations, args.steps)}
            write_csv(directory / "summary.csv", list(summary), [summary])
        elif mode == "profile":
            data = directory / "perf.data"
            command = ["perf", "record", "--delay=-1", "-e", "cpu-clock:u", "-F", "99",
                       "--call-graph", "fp", "-o", data, "--", *driver]
            if run_logged(command, directory):
                raise RuntimeError("perf record failed; see stderr.log for permissions/tool errors.")
            if run_logged(["perf", "report", "--stdio", "--no-children", "-i", data],
                          directory, "report.txt", "report.stderr.log"):
                raise RuntimeError("perf report failed; see report.stderr.log.")
            if run_logged(["perf", "script", "-G", "-i", data, "-F", "event"],
                          directory, "sample-events.txt", "script.stderr.log"):
                raise RuntimeError("Could not inspect profile samples; see script.stderr.log.")
            events = (directory / "sample-events.txt").read_text().splitlines()
            status["sample_count"] = sum(bool(line.strip()) and not line.startswith("#") for line in events)
            if not status["sample_count"]:
                raise RuntimeError("Profile contains no samples. Increase --steps and run again.")
        else:
            raw = directory / "perf-stat.txt"
            command = ["perf", "stat", "--delay=-1", "-x", ";", "--no-big-num",
                       "-e", args.events, "-o", raw, "--", *driver]
            code = run_logged(command, directory)
            rows = parse_counters(raw.read_text()) if raw.exists() else []
            write_csv(directory / "counters.csv", COUNTER_FIELDS, rows)
            write_csv(directory / "ratios.csv", ["metric", "value", "unit", "status"],
                      counter_ratios(rows))
            if code:
                raise RuntimeError("perf stat failed; see stderr.log and perf-stat.txt.")
            if not rows or any(row["status"] != "ok" for row in rows):
                raise RuntimeError("Some counters are unavailable/not counted; see counters.csv. Use --events to select supported events.")
        status["status"] = "ok"
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        status["error"] = str(error)
        summary = None
    status["finished"] = now()
    write_json(directory / "status.json", status)
    return status, summary


def create_experiment_directory(results_dir, mode):
    try:
        revision = subprocess.run(["git", "rev-parse", "--short=7", "--verify", "HEAD^{commit}"],
                                  cwd=ROOT, capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RuntimeError(f"Cannot identify Git commit: {error}") from error
    short_hash = revision.stdout.strip()
    if revision.returncode != 0 or not re.fullmatch(r"[0-9a-f]{7,64}", short_hash):
        raise RuntimeError("Cannot identify Git commit: " + (revision.stderr.strip() or "invalid Git output"))
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%fZ")
    directory = results_dir.resolve() / f"{short_hash}_{timestamp}_{mode}"
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
    directory = create_experiment_directory(args.results_dir, args.mode)
    summaries = []
    failed = False
    try:
        binary = build_and_archive(directory)
        modes = ("timing", "profile", "counters") if args.mode == "all" else (args.mode,)
        for index, path in enumerate(args.inputs, 1):
            path = Path(path).resolve()
            name = f"{index:04d}-{re.sub(r'[^A-Za-z0-9_.-]', '_', path.stem)[:80]}"
            folder = directory / name
            folder.mkdir()
            try:
                instance = folder / "input.txt"
                shutil.copyfile(path, instance)
                successful_modes = []
                for mode in modes:
                    status, summary = collect(mode, binary, instance, folder / mode, args)
                    failed |= status["status"] != "ok"
                    if status["status"] == "ok":
                        successful_modes.append(mode)
                    if summary:
                        summaries.append({"instance": name, **summary})
                if args.compare_direct:
                    failed |= not compare_direct(binary, instance, folder, successful_modes, args)
            except OSError as error:
                write_json(folder / "status.json", {"status": "failed", "error": str(error)})
                failed = True
        if args.mode in ("timing", "all"):
            fields = ["instance", "solver", "steps", "dt", "samples", "mean_seconds_per_step",
                      "variance_seconds_per_step_squared", "stddev_seconds_per_step"]
            write_csv(directory / "timing-summary.csv", fields, summaries)
    except (OSError, ValueError, RuntimeError) as error:
        (directory / "error.txt").write_text(str(error) + "\n", encoding="utf-8")
        failed = True
    except KeyboardInterrupt:
        (directory / "error.txt").write_text("Interrupted by user.\n", encoding="utf-8")
        failed = True
    for summary in summaries:
        print("Timing {instance} [{solver}]: {mean:.6g} ± {stddev:.6g} ms/step "
              "(mean ± sample SD; {steps} steps, {samples} samples)".format(
                  **summary, mean=summary["mean_seconds_per_step"] * 1000,
                  stddev=summary["stddev_seconds_per_step"] * 1000))
    print(f"Results: {directory} ({'failed' if failed else 'ok'})")
    return 1 if failed else 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    for mode in ("timing", "profile", "counters", "all"):
        sub = modes.add_parser(mode)
        sub.add_argument("inputs", nargs="+", help="instance files (shell globs supported)")
        sub.add_argument("--dt", type=float, required=True)
        sub.add_argument("--steps", type=int, required=True)
        sub.add_argument("--warmups", type=int, default=1)
        sub.add_argument("--solver", choices=("direct", "approximate"), default="direct")
        sub.add_argument("--compare-direct", action="store_true",
                         help="compare final trajectory state with an unmeasured direct run")
        sub.add_argument("--results-dir", type=Path, default=ROOT / "results",
                         help="parent directory for automatically named experiments (default: results/)")
        if mode in ("timing", "all"):
            sub.add_argument("--samples", type=int, default=10)
        else:
            sub.set_defaults(samples=1)
        if mode in ("counters", "all"):
            sub.add_argument("--events", default=DEFAULT_EVENTS, help="perf event list")
        else:
            sub.set_defaults(events=DEFAULT_EVENTS)
    args = parser.parse_args(argv)
    if not math.isfinite(args.dt) or args.dt <= 0 or args.steps <= 0 or args.warmups < 0:
        parser.error("dt and steps must be positive, dt finite, and warmups nonnegative")
    if args.mode in ("timing", "all") and args.samples < 2:
        parser.error("timing requires --samples >= 2")
    if not args.events.strip():
        parser.error("--events must not be empty")
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
