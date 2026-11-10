#!/usr/bin/env python3
"""Plot a completed --suite current benchmark, checking summaries against raw data."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import statistics
import sys
from tempfile import TemporaryDirectory
import textwrap

from run import write_csv

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


def draw(plt, status, rows, output):
    colors = {"n2": "#e68613", "rushil": "#2677b8", "src": "#269447"}
    names = {"n2": "n²", "rushil": "Rushil Barnes–Hut", "src": "src Barnes–Hut"}
    fig, ax = plt.subplots(figsize=(12, 7.5))
    clipped = False
    for name, threads in SERIES:
        points = sorted((r for r in rows if (r["implementation"], r["requested_threads"]) == (name, threads)),
                        key=lambda row: row["bodies"])
        if not points:
            continue
        means = [r["mean_ms_per_step"] for r in points]
        deviations = [r["stddev_ms_per_step"] for r in points]
        clipped |= any(mean <= sd for mean, sd in zip(means, deviations))
        ax.errorbar([r["bodies"] for r in points], means, yerr=deviations,
                    color=colors[name], marker="o", markersize=7, linewidth=2,
                    linestyle="-" if threads == 1 else "--", capsize=4,
                    label="{} ({} {})".format(names[name], threads, "thread" if threads == 1 else "threads"))
    ax.set_xscale("log")
    ax.set_yscale("log", nonpositive="clip")
    counts = sorted(status["settings"]["bodies"])
    ax.set_xticks(counts)
    ax.set_xticklabels([format(count, ",") for count in counts])
    ax.set_xlabel("Number of bodies")
    ax.set_ylabel("Mean runtime (ms / leapfrog step)")
    ax.set_title("Solver performance; error bars = sample standard deviation", pad=14)
    ax.grid(True, which="major", color="#dddddd", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=10)
    settings = status["settings"]
    machine = status.get("machine_description") or "{}; {} logical CPUs".format(
        status["platform"], status["logical_cpus"])
    sampling = "{} runs/point".format(settings["samples"])
    if settings.get("skip_large_n2"):
        sampling += "; n² at {:,} skipped".format(max(settings["bodies"]))
    elif settings.get("large_n2_samples") is not None:
        sampling += "; n² at {:,}: {} runs".format(max(settings["bodies"]), settings["large_n2_samples"])
    caption = ("{} · {} warmups per point · {} steps/run · dt={:g} · θ={:g} · OpenMP, -O3, LTO\n"
               "Mean ± 1 sample SD of whole-run averages; initial force evaluation included. n² force kernel is serial.\n{}"
               ).format(sampling, settings["warmups"], settings["steps"], settings["dt"],
                        settings["theta"], textwrap.fill(machine, width=125))
    if clipped:
        caption += "\nNonpositive lower error-bar endpoints extend below the logarithmic axis."
    fig.text(0.5, 0.02, caption, ha="center", va="bottom", fontsize=9, color="#444444")
    fig.tight_layout(rect=(0, 0.16 if clipped else 0.14, 1, 1))
    for extension in ("png", "svg"):
        fig.savefig(output / ("performance." + extension), dpi=220, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="completed current-comparison result directory")
    parser.add_argument("--output-dir", type=Path, help="default: RUN/plots")
    args = parser.parse_args()
    try:
        status, rows = load_comparison(args.run)
        output = args.output_dir or args.run / "plots"
        with TemporaryDirectory(prefix="nbody-plot-") as cache:
            os.environ.setdefault("MPLCONFIGDIR", cache)
            os.environ.setdefault("XDG_CACHE_HOME", cache)
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            output.mkdir(parents=True, exist_ok=True)
            with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 12}):
                draw(plt, status, rows, output)
        write_csv(output / "comparison.csv", list(rows[0]), rows)
        print("Saved PNG, SVG, and verified summary CSV to " + str(output))
    except (OSError, ValueError, KeyError, TypeError, ImportError) as error:
        parser.exit(1, "error: {}\n".format(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
