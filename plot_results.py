#!/usr/bin/env python3
"""Plot time per step and speedup, with all run averages, means, and sample SDs."""

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
import tempfile
import textwrap


ROOT = Path(__file__).resolve().parent
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


def draw_figure(plt, runs, observations, baseline, output, speedup=False):
    # Keep calculations in seconds per step; display runtime in milliseconds.
    if not speedup:
        observations = [[sample * 1000 for sample in samples] for samples in observations]
    # Keep each result about an inch apart; long labels should not widen every gap.
    width = max(5.2, 1.2 + len(runs))
    labels = ["\n".join(textwrap.fill(part, width=14, break_on_hyphens=False)
                        for part in run.label.split("_")) for run in runs]
    fig, ax = plt.subplots(figsize=(width, 4.6))
    color = "#3978a8"
    randomizer = random.Random(0)

    means = [statistics.mean(samples) for samples in observations]
    deviations = [statistics.stdev(samples) for samples in observations]
    for position, samples in enumerate(observations):
        # Spread every sample horizontally so even identical values remain visible.
        offsets = [0.36 * (i / (len(samples) - 1) - 0.5) for i in range(len(samples))]
        randomizer.shuffle(offsets)
        ax.scatter([position + offset for offset in offsets], samples,
                   s=24, color=color, alpha=0.22, edgecolors="none", zorder=2)
    ax.errorbar(range(len(runs)), means, yerr=deviations, fmt="o", linestyle="none",
                color=color, markersize=4.5,
                elinewidth=1.3, capsize=4, capthick=1.3, zorder=3)
    if speedup:
        ax.axhline(1, color="#999999", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_xticks(range(len(runs)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Speedup (×)" if speedup else "Time per step (ms)")
    title = "Speedup vs {}".format(baseline.label) if speedup else "Time per simulation step"
    title = textwrap.fill(title.replace("_", " "), width=int(width * 7))
    if speedup:
        title += "\nBased on time per simulation step"
    ax.set_title("{}\nInput {}".format(title, runs[0].instance), loc="left", pad=12)
    ax.set_xlim(-0.5, len(runs) - 0.5)
    ax.margins(y=0.12)
    ax.set_ylim(bottom=0)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e9e9e9", linewidth=0.7)
    ax.tick_params(axis="both", length=0, pad=8)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#cccccc")

    note = "Faint dots: whole-run averages per step\nMarkers + error bars: mean ± 1 sample SD across runs"
    if speedup:
        note += "\nSpeedup = baseline mean time per step / measured time per step"
    fig.text(0.5, 0.025, note, ha="center", va="bottom", fontsize=9, color="#666666")
    fig.tight_layout(rect=(0, 0.14 if speedup else 0.11, 1, 1), pad=1.2)
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=ROOT / "results",
                        help="saved benchmark runs (default: results/ beside this script)")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "plots",
                        help="image directory (default: plots/ beside this script)")
    parser.add_argument("--runs", nargs="+", metavar="RUN",
                        help="plot only these exact run directory names under --results-dir; "
                        "default: all runs")
    parser.add_argument("--baseline", help="baseline label or full run directory name; "
                        "defaults to +Baseline, then the earliest run for each input")
    args = parser.parse_args(argv)
    try:
        groups = load_results(args.results_dir, args.runs)
        baselines = [choose_baseline(runs, args.baseline) for runs in groups.values()]
        speedups = [speedup_samples(runs, baseline)
                    for runs, baseline in zip(groups.values(), baselines)]
        figures = []
        for (instance, runs), baseline, ratios in zip(groups.items(), baselines, speedups):
            for kind, values, speedup in (
                    ("raw_time", [run.samples for run in runs], False),
                    ("speedup", ratios, True)):
                output = args.output_dir / "{}_{}.png".format(kind, instance)
                figures.append((runs, values, baseline, output, speedup))
        args.output_dir.mkdir(parents=True, exist_ok=True)
        # Avoid requiring a writable home directory or a graphical desktop.
        with tempfile.TemporaryDirectory(prefix="nbody-matplotlib-") as config_dir:
            os.environ.setdefault("MPLCONFIGDIR", config_dir)
            os.environ.setdefault("XDG_CACHE_HOME", config_dir)
            try:
                import matplotlib
            except ImportError as error:
                raise ValueError("Matplotlib is required. Install it with: "
                                 "python3 -m pip install matplotlib") from error
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 11,
                                 "axes.titlesize": 13}):
                for runs, values, baseline, output, speedup in figures:
                    draw_figure(plt, runs, values, baseline, output,
                                speedup=speedup)
        for (instance, runs), baseline in zip(groups.items(), baselines):
            print("Input {}: {} runs, {} samples; baseline {}".format(
                instance, len(runs), sum(len(run.samples) for run in runs), baseline.label))
        for _, _, _, output, _ in figures:
            print("Saved {}".format(output))
    except (OSError, ValueError) as error:
        parser.exit(1, "error: {}\n".format(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
