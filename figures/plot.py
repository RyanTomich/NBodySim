#!/usr/bin/env python3
"""Generate all supported paper figures from saved benchmark data.

Usage: python3 figures/plot.py DATA [DATA ...] [--baseline RUN]
DATA may be a phase-summary.csv, comparison.csv, solver-sweep.csv, parallel-sweep.csv,
accuracy-tradeoff.csv, timing-summary.csv, its run
folder, or a parent containing ordinary timing runs. Outputs are fixed-name 600-dpi PNGs
saved beside this script; each run replaces the previous figure of that type. Requires Python 3.8+ and Matplotlib; no TeX installation needed.
"""

import argparse
import csv
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import random
import re
import statistics
import sys
from tempfile import TemporaryDirectory
import textwrap

OUTPUT = Path(__file__).resolve().parent

# Shared NeurIPS paper style: 5.5-inch text width and serif typography.
# Titles and experimental captions belong in the paper.
# https://media.neurips.cc/Conferences/NeurIPS2023/Styles/neurips_2023.pdf
STYLE = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
    "font.size": 9, "axes.labelsize": 10, "legend.fontsize": 9,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": "white",
}


# Colorblind-friendly palette; no boundaries artificially widen tiny segments.
GROUPS = (
    ("Force evaluation", ("force_seconds",), "#0072B2", "white"),
    ("Traversal + opening tests", ("traversal_seconds",), "#E69F00", "black"),
    ("Tree construction", ("construction_seconds",), "#CC79A7", "black"),
    ("Mass summaries", ("summaries_seconds",), "#009E73", "black"),
    ("Reuse verification", ("membership_seconds",), "#D55E00", "black"),
    ("Bounds + other overhead", ("bounds_seconds", "evaluation_other_seconds",
                                "simulation_other_seconds"), "#B8BEC7", "black"),
)


@dataclass
class PhaseShares:
    shares: dict  # (body count, actual threads) -> phase percentages.
    total_ms: dict  # Same keys; mean profiled whole-call milliseconds / step.


def read_shares(path):
    expected = {field for _, fields, _, _ in GROUPS for field in fields}
    instances, workloads = {}, {}
    with path.open(newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            count = int(row["bodies"])
            threads = int(row.get("threads") or 1)  # Historical serial data.
            field = row["field"]
            value = float(row["mean_ms_per_step"])
            if (count <= 0 or threads <= 0 or field not in expected
                    or not math.isfinite(value) or value < 0):
                raise ValueError("Invalid body count, threads, phase, or time in {}".format(path))
            if threads > 1 and row.get("phase_time_basis") != "serial_wall_and_mean_worker_elapsed":
                raise ValueError("Parallel phase times must be normalized by actual team size")
            workload = tuple(row.get(key) for key in ("instance", "steps", "dt", "theta"))
            if count in workloads and workloads[count] != workload:
                raise ValueError("Thread profiles use different workloads for N={}".format(count))
            workloads[count] = workload
            phases = instances.setdefault((count, threads), {})
            if field in phases:
                raise ValueError("Multiple instances or duplicate phases for N={}, threads={}".format(count, threads))
            phases[field] = value
    if not instances:
        raise ValueError("No phase data in {}".format(path))
    teams = {threads for _, threads in instances}
    counts = {count for count, _ in instances}
    if set(instances) != {(count, threads) for count in counts for threads in teams}:
        raise ValueError("Missing thread configuration for one or more body counts")
    shares, totals = {}, {}
    for key, phases in sorted(instances.items()):
        if set(phases) != expected:
            raise ValueError("Missing phase measurements for N={}, threads={}".format(*key))
        total = math.fsum(phases.values())
        if total <= 0:
            raise ValueError("Nonpositive total time for N={}, threads={}".format(*key))
        totals[key] = total
        shares[key] = [100 * math.fsum(phases[field] for field in fields) / total
                       for _, fields, _, _ in GROUPS]
        if not math.isclose(math.fsum(shares[key]), 100, abs_tol=1e-8):
            raise ValueError("Phase shares do not sum to 100%")
    return PhaseShares(shares, totals)


SERIES = (("n2", 1), ("rushil", 1), ("rushil", 4), ("src", 1), ("src", 4))


def load_comparison(directory):
    status = json.loads((directory / "status.json").read_text())
    if status.get("status") != "ok" or status.get("settings", {}).get("suite") != "current":
        raise ValueError("Expected a completed, successful current-suite run")
    settings = status["settings"]
    expected = {(count, name, threads) for count in settings["bodies"] for name, threads in SERIES}
    if settings.get("skip_large_n2"):
        expected.remove((max(settings["bodies"]), "n2", 1))
    with (directory / "comparison.csv").open(newline="") as source:
        rows = list(csv.DictReader(source))
    seen = set()
    for row in rows:
        count, name, threads = int(row["bodies"]), row["implementation"], int(row["requested_threads"])
        key = (count, name, threads)
        if key not in expected or key in seen:
            raise ValueError("Unexpected or duplicate configuration: {}".format(key))
        seen.add(key)
        steps = int(row["steps"])
        if (steps <= 0 or steps != settings["steps"] or float(row["dt"]) != settings["dt"]
                or row["integrator"] != "leapfrog" or int(row["force_evaluations_per_run"]) != steps + 1
                or int(row["effective_threads"]) != threads
                or row["theta"] != ("" if name == "n2" else str(settings["theta"]))):
            raise ValueError("Inconsistent workload metadata: {}".format(key))
        case = directory / str(count) / (name + "-" + str(threads) + "threads")
        batches = sorted(case.glob("batch-*"))
        if len(batches) != 1:
            raise ValueError("Expected one batch per current-suite configuration")
        batch_status = json.loads((batches[0] / "status.json").read_text())
        if (batch_status.get("status") != "ok" or batch_status.get("requested") != threads
                or batch_status.get("effective") != threads):
            raise ValueError("Incomplete batch or incorrect thread count: {}".format(key))
        with (batches[0] / "samples.csv").open(newline="") as source:
            raw = list(csv.DictReader(source))
        seconds = [float(sample["seconds"]) for sample in raw]
        expected_samples = (settings.get("large_n2_samples")
                            if name == "n2" and count == max(settings["bodies"]) else None) or settings["samples"]
        if (len(seconds) != expected_samples or len(seconds) < 2
                or int(row["samples"]) != len(seconds) or batch_status.get("samples") != len(seconds)
                or [int(sample["sample"]) for sample in raw] != list(range(1, len(seconds) + 1))
                or any(not math.isfinite(value) or value <= 0 for value in seconds)):
            raise ValueError("Invalid raw timing samples: {}".format(key))
        values = [1000 * value / steps for value in seconds]
        for field, actual in (("mean_ms_per_step", statistics.mean(values)),
                              ("stddev_ms_per_step", statistics.stdev(values)),
                              ("median_ms_per_step", statistics.median(values))):
            saved = float(row[field])
            if not math.isfinite(saved) or not math.isclose(saved, actual, rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("Summary does not match raw samples: {} {}".format(key, field))
            row[field] = actual
        for field in ("force_relative_l2", "force_rms_absolute"):
            value = float(row[field])
            if not math.isfinite(value) or value < 0:
                raise ValueError("Invalid force error: {}".format(key))
        row.update(bodies=count, requested_threads=threads)
    if seen != expected:
        raise ValueError("Missing configurations: {}".format(sorted(expected - seen)))
    return status, rows


SWEEP_SOLVERS = (("n2", "Direct", "#E69F00"),
                 ("barnes-hut", "Barnes–Hut", "#0072B2"),
                 ("fmm", "FMM", "#009E73"))


def load_accuracy_tradeoff(directory):
    status = json.loads((directory / "status.json").read_text())
    if status.get("status") != "ok" or status.get("suite") != "accuracy-tradeoff":
        raise ValueError("Expected a completed accuracy tradeoff run")
    settings = status["settings"]
    if (settings["bodies"] <= 0 or settings["steps"] <= 0 or settings["samples"] < 2
            or settings["warmups"] < 0 or not math.isfinite(settings["dt"]) or settings["dt"] <= 0
            or settings["threads"] != 1 or settings["openmp"] is not False
            or settings["orders"] != list(range(1, 7)) or not settings["thetas"]
            or len(set(settings["thetas"])) != len(settings["thetas"])
            or any(not math.isfinite(theta) or theta < 0 for theta in settings["thetas"])):
        raise ValueError("Invalid accuracy tradeoff settings")
    # New runs save a single untimed direct reference outside the timing table.
    # Keep historical runs with a timed n2 row readable without modifying them.
    expected = {}
    if "reference_threads" in settings:
        threads = settings["reference_threads"]
        reference = json.loads((directory / "reference/status.json").read_text())
        if (type(threads) is not int or threads < 1
                or any(reference.get(key) != value for key, value in {
                    "status": "ok", "mode": "reference", "solver": "n2", "integrator": "leapfrog",
                    "dt": settings["dt"], "steps": settings["steps"], "samples": 1, "warmups": 0,
                    "threads": threads, "openmp": threads > 1}.items())
                or not (directory / "reference/final.txt").is_file()):
            raise ValueError("Invalid or incomplete direct reference")
    else:
        expected[("n2", None, None)] = "n2"
    expected.update({("barnes-hut", theta, None): "barnes-hut-theta-" + str(theta)
                     for theta in settings["thetas"]})
    expected.update({("fmm", None, order): "fmm-order-" + str(order) for order in range(1, 7)})
    with (directory / "accuracy-tradeoff.csv").open(newline="") as source:
        rows = list(csv.DictReader(source))
    seen = set()
    for row in rows:
        theta = float(row["theta"]) if row["theta"] else None
        order = int(row["order"]) if row["order"] else None
        key = row["solver"], theta, order
        if key not in expected or key in seen or row["case"] != expected[key]:
            raise ValueError("Unexpected or duplicate accuracy configuration")
        seen.add(key)
        for field in ("bodies", "steps", "samples", "warmups", "threads"):
            if int(row[field]) != settings[field]:
                raise ValueError("Inconsistent accuracy workload: " + field)
        if float(row["dt"]) != settings["dt"]:
            raise ValueError("Inconsistent accuracy timestep")
        folder = directory / row["case"] / "timing"
        saved = json.loads((folder / "status.json").read_text())
        if saved.get("status") != "ok" or saved.get("solver") != row["solver"]:
            raise ValueError("Incomplete accuracy measurement")
        if saved.get("theta") != theta or saved.get("order") != order:
            raise ValueError("Inconsistent solver parameter")
        if (saved.get("dt") != settings["dt"] or saved.get("steps") != settings["steps"]
                or saved.get("integrator") != "leapfrog"):
            raise ValueError("Inconsistent trajectory measurement")
        with (folder / "samples.csv").open(newline="") as source:
            samples = list(csv.DictReader(source))
        if [int(sample["sample"]) for sample in samples] != list(range(1, settings["samples"] + 1)):
            raise ValueError("Missing or duplicate accuracy samples")
        durations = [float(sample["seconds"]) / settings["steps"] for sample in samples]
        if any(not math.isfinite(value) or value <= 0 for value in durations):
            raise ValueError("Invalid accuracy timing sample")
        for prefix, value in (("mean", statistics.mean(durations)),
                              ("stddev", statistics.stdev(durations))):
            field = prefix + "_seconds_per_step"
            if not math.isclose(float(row[field]), value, rel_tol=1e-12, abs_tol=1e-15):
                raise ValueError("Accuracy summary does not match raw samples")
            row[field] = value
        with (folder / "error.csv").open(newline="") as source:
            # Older runs also saved velocity and acceleration errors. Only
            # position is needed; their other files and columns stay untouched.
            errors = [error for error in csv.DictReader(source) if error["field"] == "position"]
        if len(errors) != 1:
            raise ValueError("Missing or duplicate position error metric")
        error = errors[0]
        value = float(error["relative_l2_error"]) if error["relative_l2_error"] else None
        aggregate = float(row["position_relative_l2"]) if row["position_relative_l2"] else None
        if (value != aggregate or (value is not None and (not math.isfinite(value) or value < 0))
                or error["relative_status"] != ("undefined" if value is None else "ok")
                or (row["solver"] == "n2" and value != 0)):
            raise ValueError("Invalid position error metric")
        row["position_relative_l2"] = value
        row.update(theta=theta, order=order)
    if seen != set(expected):
        raise ValueError("Missing accuracy configurations")
    return status, rows


def pareto_frontier(points):
    """Indices of nondominated (error, runtime) points; exact ties both survive."""
    return [i for i, (error, runtime) in enumerate(points)
            if not any(other_error <= error and other_time <= runtime
                       and (other_error < error or other_time < runtime)
                       for other_error, other_time in points)]


def draw_accuracy_trajectory(plt, data, min_theta=None):
    status, rows = data
    metric = "position_relative_l2"
    mean, sd = "mean_seconds_per_step", "stddev_seconds_per_step"
    points = []
    for row in rows:
        if row["solver"] == "n2":
            continue
        if row["solver"] == "barnes-hut" and min_theta is not None and row["theta"] < min_theta:
            continue
        if row[metric] is None or row[metric] <= 0:
            print("Omitting {} from {}: zero or undefined error".format(row["case"], metric), file=sys.stderr)
        else:
            points.append(row)
    if not points:
        raise ValueError("No positive errors available for " + metric)
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    for solver, label, color in SWEEP_SOLVERS[1:]:
        series = sorted((row for row in points if row["solver"] == solver),
                        key=lambda row: row["theta"] if solver == "barnes-hut" else row["order"])
        if not series:
            continue
        ax.errorbar([row[metric] for row in series], [1000 * row[mean] for row in series],
                    yerr=[1000 * row[sd] for row in series], color=color, label=label,
                    marker="s" if solver == "barnes-hut" else "^", markersize=4,
                    linewidth=1, capsize=3, elinewidth=0.8)
        for row in series:
            label = r"$\theta={:g}$".format(row["theta"]) if solver == "barnes-hut" else "$p={}$".format(row["order"])
            # Separate labels where the two curves cluster at low runtimes.
            crossing = solver == "barnes-hut" and row["theta"] == 0.3
            order_four = solver == "fmm" and row["order"] == 4
            theta_five = solver == "barnes-hut" and row["theta"] == 0.5
            lower_theta = solver == "barnes-hut" and row["theta"] == 0.7
            ax.annotate(label, (row[metric], 1000 * row[mean]),
                        xytext=(0, -16) if crossing else ((0, -28) if order_four else ((0, 10) if theta_five else ((8, -8) if lower_theta else (5, 6)))),
                        ha="center" if crossing or order_four or theta_five else "left",
                        va="top" if crossing or order_four or lower_theta else "bottom",
                        arrowprops={"arrowstyle": "-", "color": "0.5", "linewidth": 0.5} if crossing or order_four else None,
                        textcoords="offset points", fontsize=7)
    frontier = [points[i] for i in pareto_frontier([(row[metric], row[mean]) for row in points])]
    frontier.sort(key=lambda row: row[metric])
    ax.plot([row[metric] for row in frontier], [1000 * row[mean] for row in frontier],
            color="0.25", linestyle="--", linewidth=0.8, marker="o", markersize=8,
            markerfacecolor="none", markeredgewidth=0.8, label="Pareto frontier", zorder=4)
    ax.set_xscale("log")
    ax.set_yscale("log", nonpositive="clip")
    ax.set_xlabel("Relative L2 final-position error")
    ax.set_ylabel("Simulation runtime (ms / step)")
    ax.margins(x=0.18, y=0.2)
    ax.grid(which="major", color="0.9", linewidth=0.4)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8)
    settings = status["settings"]
    caption = "N={:,}; serial; {} runs; bars: between-run SD".format(settings["bodies"], settings["samples"])
    fig.text(0.55, 0.025, caption, ha="center", fontsize=8)
    fig.subplots_adjust(left=0.15, right=0.97, top=0.97, bottom=0.22)
    if any(row[sd] >= row[mean] for row in points):
        print("Note: nonpositive lower SD endpoints extend below the log axis.", file=sys.stderr)
    return fig


def load_parallel_sweep(directory):
    return load_solver_sweep(directory, parallel=True)


def load_solver_sweep(directory, parallel=False):
    suite = "parallel-sweep" if parallel else "solver-sweep"
    solvers = SWEEP_SOLVERS[:2] if parallel else SWEEP_SOLVERS
    teams = (1, 4, 8) if parallel else (1,)
    status = json.loads((directory / "status.json").read_text())
    if status.get("status") != "ok" or status.get("suite") != suite:
        raise ValueError("Expected a successful solver sweep")
    settings = status["settings"]
    counts, steps = settings["bodies"], settings["steps"]
    samples_count = settings["samples"]
    legacy = status.get("error_bars") == "individual_step_sample_sd"
    if (steps < 2 or not counts or len(set(counts)) != len(counts)
            or any(count <= 0 for count in counts)
            or settings["threads"] != ([1, 4, 8] if parallel else 1)
            or type(samples_count) is not int
            or (samples_count != 1 if legacy else samples_count < 2)
            or settings["warmups"] != 0
            or status.get("error_bars") not in ("individual_step_sample_sd", "whole_run_sample_sd")):
        raise ValueError("Invalid sweep sample count or error-bar definition")
    expected = {(count, solver, threads) for count in counts
                for solver, _, _ in solvers for threads in teams}
    if "configurations" in settings:
        configured = [(point["bodies"], point["solver"], point["threads"])
                      for point in settings["configurations"]]
        allowed = {(count, solver, threads) for count in counts
                   for solver, _, _ in SWEEP_SOLVERS
                   for threads in teams}
        expected = set(configured)
        if (not expected or len(expected) != len(configured) or not expected <= allowed
                or {point[0] for point in expected} != set(counts)):
            raise ValueError("Invalid solver-sweep configurations")
    with (directory / (suite + ".csv")).open(newline="") as source:
        rows = list(csv.DictReader(source))
    seen = set()
    for row in rows:
        count, solver = int(row["bodies"]), row["solver"]
        threads = int(row["threads"])
        key = count, solver, threads
        if key not in expected or key in seen:
            raise ValueError("Unexpected or duplicate solver-sweep point: {}".format(key))
        seen.add(key)
        if (int(row["samples"]) != samples_count
                or int(row["steps"]) != steps or float(row["dt"]) != settings["dt"]
                or row["theta"] != (str(settings["theta"]) if solver == "barnes-hut" else "")
                or row["order"] != (str(settings["order"]) if solver == "fmm" else "")):
            raise ValueError("Inconsistent solver-sweep workload")
        case = directory / str(count) / ("{}-{}threads".format(solver, threads) if parallel else solver)
        case_status = json.loads((case / "status.json").read_text())
        if (case_status.get("status") != "ok" or case_status.get("threads") != threads
                or case_status.get("solver") != solver or case_status.get("steps") != steps
                or case_status.get("dt") != settings["dt"] or case_status.get("samples") != samples_count
                or case_status.get("warmups") != 0 or case_status.get("integrator") != "leapfrog"
                or case_status.get("force_evaluations_per_run") != steps + 1):
            raise ValueError("Incomplete or inconsistent solver-sweep case")
        if parallel and (case_status.get("openmp") is not True
                         or any(int(row[field]) != threads or case_status.get(field) != threads
                                for field in ("requested_threads", "effective_threads"))):
            raise ValueError("Inconsistent parallel thread metadata")
        with (case / "samples.csv").open(newline="") as source:
            samples = list(csv.DictReader(source))
        with (case / "steps.csv").open(newline="") as source:
            raw_steps = list(csv.DictReader(source))
        if [int(x["sample"]) for x in samples] != list(range(1, samples_count + 1)):
            raise ValueError("Missing or duplicate whole-run samples")
        if [(int(x.get("sample", 1)) if legacy else int(x["sample"]), int(x["step"]))
                for x in raw_steps] != [(sample, step) for sample in range(1, samples_count + 1)
                                      for step in range(1, steps + 1)]:
            raise ValueError("Expected every step timing for every run")
        durations = [float(x["seconds"]) for x in samples]
        if any(not math.isfinite(x) or x <= 0 for x in durations):
            raise ValueError("Invalid solver-sweep durations")
        for sample, seconds in enumerate(durations, 1):
            values = [float(x["seconds"]) for x in raw_steps if int(x.get("sample", 1)) == sample]
            if (any(not math.isfinite(x) or x < 0 for x in values)
                    or math.fsum(values) > seconds + 1e-9):
                raise ValueError("Invalid solver-sweep step durations")
        if legacy:
            calculated = {"seconds": durations[0], "mean_ms_per_step": 1000 * durations[0] / steps,
                          "step_mean_ms": 1000 * statistics.mean(values),
                          "step_sd_ms": 1000 * statistics.stdev(values)}
        else:
            averages = [1000 * seconds / steps for seconds in durations]
            calculated = {"mean_ms_per_step": statistics.mean(averages),
                          "stddev_ms_per_step": statistics.stdev(averages)}
        for field, value in calculated.items():
            if not math.isclose(float(row[field]), value, rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("Solver-sweep summary does not match raw timings")
            row[field] = value
        for field in ("position", "velocity", "acceleration"):
            for suffix in ("_relative_l2", "_rms"):
                text = row[field + suffix]
                if text == "" and suffix == "_relative_l2":
                    continue  # Zero reference norm can make relative error undefined.
                value = float(text)
                if not math.isfinite(value) or value < 0:
                    raise ValueError("Invalid trajectory error")
        row.update(bodies=count, threads=threads)
    if seen != expected:
        raise ValueError("Missing solver-sweep configurations")
    return status, rows


def draw_parallel_sweep(plt, data):
    return draw_solver_sweep(plt, data, parallel=True)


def draw_parallel_sweep_normalized(plt, data):
    """Divide runtime and SD by each solver's asymptotic growth factor.

    Separate linear y axes expose departures from constant normalized cost;
    these factors are not operation counts or predictions of absolute runtime.
    """
    status, rows = data
    legacy = status.get("error_bars") == "individual_step_sample_sd"
    error_field = "step_sd_ms" if legacy else "stddev_ms_per_step"
    factors = {"n2": lambda n: n ** 2,
               "barnes-hut": lambda n: n * math.log2(n),
               "fmm": lambda n: n}
    formulas = {"n2": r"T/n^2", "barnes-hut": r"T/(n\log_2 n)", "fmm": r"T/n"}
    solvers = [(solver, label, color) for solver, label, color in SWEEP_SOLVERS
               if any(row["solver"] == solver for row in rows)]
    fig, axes = plt.subplots(len(solvers), 1, sharex=True, squeeze=False,
                             figsize=(5.5, 2 * len(solvers) + 0.5))
    fig.suptitle("Runtime normalized by expected scaling", fontsize=11, y=0.985)
    fig.text(0.55, 0.95, "Flat curves indicate runtime growing at the expected rate.",
             ha="center", fontsize=8)
    for ax, (solver, label, color) in zip(axes[:, 0], solvers):
        for threads, line, marker in ((1, "-", "o"), (4, "--", "s"), (8, ":", "^")):
            points = sorted((row for row in rows if row["solver"] == solver and row["threads"] == threads),
                            key=lambda row: row["bodies"])
            # n log(n) is zero at n=1, so that point has no normalized value.
            points = [row for row in points if factors[solver](row["bodies"]) > 0]
            if not points:
                continue
            scales = [1e6 / factors[solver](row["bodies"]) for row in points]  # ms -> ns
            ax.errorbar([row["bodies"] for row in points],
                        [row["mean_ms_per_step"] * scale for row, scale in zip(points, scales)],
                        yerr=[row[error_field] * scale for row, scale in zip(points, scales)],
                        label="{}T".format(threads), color=color, linestyle=line, marker=marker,
                        markersize=3.5, linewidth=1, capsize=3, elinewidth=0.8)
        ax.set_title(label, loc="left", fontsize=10)
        ax.set_ylabel("Normalized runtime\n${}$ (ns / step)".format(formulas[solver]), fontsize=9)
        ax.set_xscale("log")
        ax.set_ylim(bottom=0)
        ax.grid(which="major", color="0.9", linewidth=0.4)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, ncol=3, loc="lower right", bbox_to_anchor=(1, 1),
                  borderaxespad=0, columnspacing=1.2, handletextpad=0.5, fontsize=8)
    counts = sorted(status["settings"]["bodies"])
    axes[-1, 0].set_xticks(counts)
    axes[-1, 0].set_xticklabels(["{:g}k".format(count / 1000) if count >= 1000 else str(count)
                                for count in counts])
    axes[-1, 0].set_xlabel("Number of bodies, $n$")
    caption = ("One {}-step run; bars: normalized step-to-step SD".format(status["settings"]["steps"])
               if legacy else "{} runs × {} steps; bars: normalized between-run SD".format(
                   status["settings"]["samples"], status["settings"]["steps"]))
    fig.text(0.55, 0.025, caption, ha="center", fontsize=8)
    fig.subplots_adjust(left=0.18, right=0.97, top=0.89, bottom=0.13, hspace=0.45)
    return fig


def draw_solver_sweep(plt, data, parallel=False):
    status, rows = data
    legacy = status.get("error_bars") == "individual_step_sample_sd"
    error_field = "step_sd_ms" if legacy else "stddev_ms_per_step"
    fig, ax = plt.subplots(figsize=(5.5, 4 if parallel else 3.5))
    styles = ((1, "-"), (4, "--"), (8, ":")) if parallel else ((1, "-"),)
    markers = {"n2": "o", "barnes-hut": "s", "fmm": "^"}
    for solver, label, color in SWEEP_SOLVERS:
        for threads, line in styles:
            points = sorted((row for row in rows if row["solver"] == solver and row["threads"] == threads),
                            key=lambda row: row["bodies"])
            if not points:
                continue
            ax.errorbar([row["bodies"] for row in points],
                        [row["mean_ms_per_step"] for row in points],
                        yerr=[row[error_field] for row in points],
                        label="{} ({}T)".format(label, threads) if parallel else label,
                        color=color, marker=markers[solver], linestyle=line, markersize=3.5,
                        linewidth=1, capsize=3, elinewidth=0.8)
    counts = sorted(status["settings"]["bodies"])
    ax.set_xscale("log")
    ax.set_yscale("log", nonpositive="clip")
    ax.set_xticks(counts)
    ax.set_xticklabels(["{:g}k".format(count / 1000) if count >= 1000 else str(count) for count in counts])
    ax.set_xlabel("Number of bodies" if parallel else "Number of bodies (single thread)")
    ax.set_ylabel("Runtime (ms / step)")
    ax.grid(which="major", color="0.9", linewidth=0.4)
    ax.set_axisbelow(True)
    if parallel:
        ax.legend(frameon=False, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.02),
                  borderaxespad=0, columnspacing=1.2, handletextpad=0.5, fontsize=8)
    else:
        ax.legend(frameon=False)
    caption = ("One {}-step run per point; bars: step-to-step SD".format(status["settings"]["steps"])
               if legacy else "{} runs × {} steps; bars: between-run SD".format(
                   status["settings"]["samples"], status["settings"]["steps"]))
    fig.text(0.55, 0.025, caption, ha="center", fontsize=8)
    fig.subplots_adjust(left=0.14, right=0.98, top=0.8 if parallel else 0.98,
                        bottom=0.2 if parallel else 0.23)
    if any(row[error_field] >= row["mean_ms_per_step"] for row in rows):
        print("Note: nonpositive lower SD endpoints extend below the log axis.", file=sys.stderr)
    return fig


RUN_NAME = re.compile(
    r"^[^_]+_(\d{8}-\d{6}\.\d+Z)_(?:timing|all)(?:_(.+))?$"
)


@dataclass
class Timing:
    run: Path
    label: str
    instance: str
    samples: list  # Whole-run averages in seconds per step.
    input_path: Path


def run_details(run):
    match = RUN_NAME.fullmatch(run.name)
    if match:
        return match.group(1), match.group(2) or run.name
    return run.name, run.name.rsplit("_", 1)[-1]


def load_results(results_dir, selected_runs=None):
    """Normalize each raw duration by its run's step count; keep runs separate."""
    groups = {}
    if selected_runs is None:
        directories = sorted(results_dir.glob("*/*/timing"))
    else:
        available = {run.name: run for run in results_dir.iterdir() if run.is_dir()}
        missing = set(selected_runs) - available.keys()
        if missing:
            raise ValueError("Unknown run(s) in {}: {}".format(results_dir, ", ".join(sorted(missing))))
        directories = sorted(directory for name in set(selected_runs)
                             for directory in available[name].glob("*/timing"))
    for directory in directories:
        if not (directory / "status.json").exists():
            print("Skipping {}: no timing status yet".format(directory), file=sys.stderr)
            continue
        try:
            status = json.loads((directory / "status.json").read_text())
            if status.get("status") != "ok":
                reason = status.get("error", "timing did not finish successfully")
                print("Skipping {}: {}".format(directory, reason), file=sys.stderr)
                continue
            if "steps" not in status:
                print("Skipping {}: no recorded step count; cannot calculate time per step".format(
                    directory), file=sys.stderr)
                continue
            steps = status["steps"]
            if type(steps) is not int or steps <= 0:
                raise ValueError("expected a positive integer step count")
            with (directory / "samples.csv").open(newline="") as source:
                samples = [float(row["seconds"]) for row in csv.DictReader(source)]
            if len(samples) < 2 or any(not math.isfinite(x) or x < 0 for x in samples):
                raise ValueError("expected at least two finite, nonnegative timings")
            samples = [sample / steps for sample in samples]
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
            # Do not silently discard observations from a successful but corrupt run.
            raise ValueError("{}: {}".format(directory, error)) from error

        run = directory.parent.parent
        _, label = run_details(run)
        instance = re.sub(r"^\d+-", "", directory.parent.name)
        if status.get("solver"):
            settings = status["solver"]
            if status.get("order") is not None:
                settings += " order={}".format(status["order"])
            if status.get("theta") is not None:
                settings += " theta={:g}".format(status["theta"])
            label += " [" + settings + "]"
        timing = Timing(run, label, instance, samples, directory.parent / "input.txt")
        groups.setdefault(instance, []).append(timing)

    if not groups:
        raise ValueError("No successful timing results with recorded step counts found in {}".format(results_dir))
    for runs in groups.values():
        runs.sort(key=lambda timing: (run_details(timing.run)[0], timing.run.name))
    return dict(sorted(groups.items()))


def choose_baseline(runs, requested=None):
    if requested:
        matches = [run for run in runs if requested in (run.label, run.run.name, run_details(run.run)[1])]
        if len(matches) != 1:
            raise ValueError(
                "Input {}: baseline {!r} matches {} runs; use a unique label or "
                "the full result directory name".format(runs[0].instance, requested, len(matches))
            )
        return matches[0]
    return next((run for run in runs if run_details(run.run)[1].lstrip("+").casefold() == "baseline"), runs[0])


def speedup_samples(runs, baseline):
    """Compare each time per step against the fixed baseline mean time per step."""
    reference = statistics.mean(baseline.samples)
    if any(sample <= 0 for run in runs for sample in run.samples):
        raise ValueError("Input {}: speedup requires positive timings".format(baseline.instance))
    baseline_input = baseline.input_path.read_bytes()
    for run in runs:
        if run.input_path.read_bytes() != baseline_input:
            raise ValueError("Input {} differs between {} and {}; cannot compare speedup".format(
                run.instance, baseline.run.name, run.run.name))
    return [[reference / sample for sample in run.samples] for run in runs]



def draw_phase_shares(plt, profile):
    keys = list(profile.shares)
    teams = sorted({threads for _, threads in keys})
    paired = len(teams) > 1
    positions = [i + 0.45 * (i // len(teams)) for i in range(len(keys))]
    height = max(3.1, 1.7 + 0.34 * len(keys)) if paired else max(3.1, 1.4 + 0.425 * len(keys))
    fig, ax = plt.subplots(figsize=(5.5, height))
    fig.subplots_adjust(left=0.23 if paired else 0.16, right=0.84 if paired else 0.98,
                        bottom=1.25 / height if paired else 0.37, top=0.94 if paired else 0.98)
    left = [0.0] * len(keys)
    for group, (label, _, color, text_color) in enumerate(GROUPS):
        values = [profile.shares[key][group] for key in keys]
        ax.barh(positions, values, left=left, height=0.72 if paired else 0.62,
                color=color, edgecolor="none", label=label, zorder=2)
        for position, value, start in zip(positions, values, left):
            if value >= (15 if paired else 10):
                ax.text(start + value / 2, position, "{:.1f}%".format(value),
                        ha="center", va="center", color=text_color, fontsize=9)
        left = [start + value for start, value in zip(left, values)]
    ax.set_yticks(positions)
    ax.set_yticklabels(["{:,} / {}T".format(*key) if paired else "{:,}".format(key[0]) for key in keys])
    ax.set_ylim(positions[-1] + 0.55, -0.55)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_xlabel("Attributed share of runtime (%)" if max(teams) > 1 else "Share of simulation runtime (%)", labelpad=5)
    ax.set_ylabel("Bodies / threads" if paired else "Number of bodies", labelpad=6)
    if paired:
        ax.text(1.04, 1.03, "ms/step", transform=ax.transAxes, ha="left", va="bottom", fontsize=8)
        for position, key in zip(positions, keys):
            ax.text(104, position, "{:.2f}".format(profile.total_ms[key]),
                    ha="left", va="center", fontsize=8, clip_on=False)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="0.9", linewidth=0.4)
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center",
               bbox_to_anchor=(0.5, 0.015), ncol=2, frameon=False,
               handlelength=1.3, handleheight=0.8, columnspacing=1.6,
               labelspacing=0.55)
    return fig


def draw_comparison(plt, data):
    status, rows = data
    colors = {"n2": "#E69F00", "rushil": "#0072B2", "src": "#009E73"}
    names = {"n2": "Direct", "rushil": "Rushil Barnes-Hut", "src": "src Barnes-Hut"}
    markers = {"n2": "o", "rushil": "s", "src": "^"}
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    for name, threads in SERIES:
        points = sorted((r for r in rows if (r["implementation"], r["requested_threads"]) == (name, threads)),
                        key=lambda row: row["bodies"])
        if not points:
            continue
        ax.errorbar([r["bodies"] for r in points], [r["mean_ms_per_step"] for r in points],
                    yerr=[r["stddev_ms_per_step"] for r in points], color=colors[name],
                    marker=markers[name], markersize=3.5, linewidth=1,
                    linestyle="-" if threads == 1 else "--", capsize=2, elinewidth=0.7,
                    label="{} ({}T)".format(names[name], threads))
    ax.set_xscale("log")
    ax.set_yscale("log", nonpositive="clip")
    counts = sorted(status["settings"]["bodies"])
    ax.set_xticks(counts)
    ax.set_xticklabels([format(count, ",") for count in counts])
    ax.set_xlabel("Number of bodies")
    ax.set_ylabel("Runtime (ms / step)")
    ax.grid(which="major", color="0.9", linewidth=0.4)
    ax.set_axisbelow(True)
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center",
               bbox_to_anchor=(0.5, 0.01), ncol=2, frameon=False,
               columnspacing=1.5, handlelength=2)
    fig.subplots_adjust(left=0.14, right=0.98, top=0.98, bottom=0.36)
    if any(row["mean_ms_per_step"] <= row["stddev_ms_per_step"] for row in rows):
        print("Note: nonpositive lower SD endpoints extend below the log axis.", file=sys.stderr)
    if status["settings"].get("skip_large_n2"):
        print("Note: largest direct-solver point was explicitly skipped.", file=sys.stderr)
    if len({int(row["samples"]) for row in rows}) > 1:
        print("Note: sample counts differ; retain this detail in the paper caption.", file=sys.stderr)
    return fig


def draw_timing(plt, data, speedup=False):
    runs, observations, baseline = data
    observations = observations if speedup else [[1000 * x for x in samples] for samples in observations]
    # Short series labels belong inside the figure; full source paths stay in
    # the data. Duplicate solver labels receive local indices, printed below.
    labels = [run.label.rsplit(" [", 1)[-1].rstrip("]") if " [" in run.label
              else run_details(run.run)[1] for run in runs]
    labels = [label if labels.count(label) == 1 else "{} ({})".format(label, i + 1)
              for i, label in enumerate(labels)]
    labels = [textwrap.fill(label, 19, break_long_words=False) for label in labels]
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    randomizer = random.Random(0)
    for position, samples in enumerate(observations):
        offsets = [0.36 * (i / (len(samples) - 1) - 0.5) for i in range(len(samples))]
        randomizer.shuffle(offsets)
        ax.scatter([position + x for x in offsets], samples, s=8, color="#0072B2",
                   alpha=0.25, edgecolors="none", zorder=2,
                   label="Whole-run averages" if position == 0 else None)
    ax.errorbar(range(len(runs)), [statistics.mean(x) for x in observations],
                yerr=[statistics.stdev(x) for x in observations], fmt="o", linestyle="none",
                color="#0072B2", markersize=3.5, elinewidth=0.8, capsize=2,
                zorder=3, label="Mean ± sample SD")
    if speedup:
        ax.axhline(1, color="0.6", linewidth=0.6, linestyle="--", zorder=1)
    ax.set_xticks(range(len(runs)))
    ax.set_xticklabels(labels, rotation=30 if len(runs) > 3 else 0,
                       ha="right" if len(runs) > 3 else "center")
    ax.set_ylabel("Speedup relative to baseline" if speedup else "Runtime (ms / step)")
    ax.set_xlim(-0.5, len(runs) - 0.5)
    ax.margins(y=0.15)
    ax.set_ylim(bottom=0)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="0.9", linewidth=0.4)
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, 0.01), frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 1), pad=0.7)
    print("Input {}: baseline {}; series: {}".format(
        runs[0].instance, baseline.run, ", ".join(run.label for run in runs)))
    return fig


def prepare_figures(paths, baseline=None):
    """Validate all data and output names before writing any figures.

    Each supported dataset has one dispatch branch here. New plot types add
    a loader and a draw function in this file and reuse STYLE/save_figure.
    """
    jobs, seen_data, seen_names = [], set(), set()
    for path in paths:
        path = path.expanduser().resolve()
        if not path.exists():
            raise ValueError("Data path does not exist: {}".format(path))
        directory = path.parent if path.is_file() else path
        if directory in seen_data:
            raise ValueError("Dataset supplied more than once: {}".format(directory))
        seen_data.add(directory)
        if path.is_file():
            kind = path.name
        else:
            candidates = [name for name in ("phase-summary.csv", "comparison.csv", "solver-sweep.csv", "parallel-sweep.csv", "accuracy-tradeoff.csv", "timing-summary.csv")
                          if (directory / name).is_file()]
            if len(candidates) > 1:
                raise ValueError("Ambiguous dataset {}; pass a specific summary CSV".format(directory))
            kind = candidates[0] if candidates else "timing-parent"
        if kind == "phase-summary.csv":
            jobs.append(("phase-shares", draw_phase_shares,
                         read_shares(directory / kind)))
        elif kind == "accuracy-tradeoff.csv":
            data = load_accuracy_tradeoff(directory)
            jobs.append(("accuracy-trajectory", draw_accuracy_trajectory, data))
        elif kind == "solver-sweep.csv":
            jobs.append(("solver-sweep", draw_solver_sweep, load_solver_sweep(directory)))
        elif kind == "parallel-sweep.csv":
            data = load_parallel_sweep(directory)
            jobs.append(("parallel-sweep", draw_parallel_sweep, data))
            jobs.append(("parallel-sweep-normalized", draw_parallel_sweep_normalized, data))
        elif kind == "comparison.csv":
            jobs.append(("performance", draw_comparison, load_comparison(directory)))
        elif kind in ("timing-summary.csv", "timing-parent"):
            groups = (load_results(directory.parent, [directory.name]) if kind == "timing-summary.csv"
                      else load_results(directory))
            for instance, runs in groups.items():
                reference = choose_baseline(runs, baseline)
                ratios = speedup_samples(runs, reference)
                jobs.append(("runtime", draw_timing, (runs, [run.samples for run in runs], reference)))
                jobs.append(("speedup", draw_speedup, (runs, ratios, reference)))
        else:
            raise ValueError("Unsupported data file: {}".format(path))
    for name, _, _ in jobs:
        if name in seen_names:
            raise ValueError("Output name collision: {}; supply only one dataset or input group per figure type".format(name))
        seen_names.add(name)
    return jobs


def draw_speedup(plt, data):
    return draw_timing(plt, data, speedup=True)


def save_figure(fig, name):
    # Preserve the chosen physical size; do not silently resize with a tight crop.
    path = OUTPUT / (name + ".png")
    fig.savefig(path, dpi=600)
    print("Saved " + str(path))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", nargs="+", type=Path, help="summary CSVs, run folders, or timing-run parent folders")
    parser.add_argument("--baseline", help="baseline run name/label for ordinary timing and speedup figures")
    parser.add_argument("--accuracy-min-theta", type=float,
                        help="omit smaller Barnes–Hut theta values from accuracy figures without changing saved data")
    args = parser.parse_args(argv)
    if args.accuracy_min_theta is not None and (not math.isfinite(args.accuracy_min_theta) or args.accuracy_min_theta < 0):
        parser.error("accuracy minimum theta must be finite and nonnegative")
    try:
        jobs = prepare_figures(args.data, args.baseline)
        # Keep imports lazy: data validation and repository tests need only the
        # standard library; drawing alone requires Matplotlib.
        with TemporaryDirectory(prefix="nbody-figures-") as cache:
            os.environ.setdefault("MPLCONFIGDIR", cache)
            os.environ.setdefault("XDG_CACHE_HOME", cache)
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            with plt.rc_context(STYLE):
                for name, draw, data in jobs:
                    fig = (draw(plt, data, min_theta=args.accuracy_min_theta)
                           if name == "accuracy-trajectory" else draw(plt, data))
                    try:
                        save_figure(fig, name)
                    finally:
                        plt.close(fig)
    except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
        parser.exit(1, "error: {}\n".format(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
