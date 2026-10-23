#!/usr/bin/env python3
"""Plot every successful timing sample, with means and sample standard deviations."""

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


ROOT = Path(__file__).resolve().parent
RUN_NAME = re.compile(
    r"^[^_]+_(\d{8}-\d{6}\.\d+Z)_(?:timing|all)(?:_(.+))?$"
)


@dataclass
class Timing:
    run: Path
    label: str
    instance: str
    samples: list
    input_path: Path


def run_details(run):
    match = RUN_NAME.fullmatch(run.name)
    if match:
        return match.group(1), match.group(2) or run.name
    return run.name, run.name.rsplit("_", 1)[-1]


def load_results(results_dir):
    """Keep runs separate and read raw samples, never averages of summaries."""
    groups = {}
    for directory in sorted(results_dir.glob("*/*/timing")):
        try:
            status = json.loads((directory / "status.json").read_text())
            if status.get("status") != "ok":
                reason = status.get("error", "timing did not finish successfully")
                print("Skipping {}: {}".format(directory, reason), file=sys.stderr)
                continue
            with (directory / "samples.csv").open(newline="") as source:
                samples = [float(row["seconds"]) for row in csv.DictReader(source)]
            if len(samples) < 2 or any(not math.isfinite(x) or x < 0 for x in samples):
                raise ValueError("expected at least two finite, nonnegative timings")
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
            # Do not silently discard observations from a successful but corrupt run.
            raise ValueError("{}: {}".format(directory, error)) from error

        run = directory.parent.parent
        _, label = run_details(run)
        instance = re.sub(r"^\d+-", "", directory.parent.name)
        timing = Timing(run, label, instance, samples, directory.parent / "input.txt")
        groups.setdefault(instance, []).append(timing)

    if not groups:
        raise ValueError("No successful timing results found in {}".format(results_dir))
    for runs in groups.values():
        runs.sort(key=lambda timing: (run_details(timing.run)[0], timing.run.name))
    return dict(sorted(groups.items()))


def choose_baseline(runs, requested=None):
    if requested:
        matches = [run for run in runs if requested in (run.label, run.run.name)]
        if len(matches) != 1:
            raise ValueError(
                "Input {}: baseline {!r} matches {} runs; use a unique label or "
                "the full result directory name".format(runs[0].instance, requested, len(matches))
            )
        return matches[0]
    return next((run for run in runs if run.label.lstrip("+").casefold() == "baseline"), runs[0])


def speedup_samples(runs, baseline):
    """Normalize each observation by the fixed baseline arithmetic mean."""
    reference = statistics.mean(baseline.samples)
    if any(sample <= 0 for run in runs for sample in run.samples):
        raise ValueError("Input {}: speedup requires positive timings".format(baseline.instance))
    baseline_input = baseline.input_path.read_bytes()
    for run in runs:
        if run.input_path.read_bytes() != baseline_input:
            raise ValueError("Input {} differs between {} and {}; cannot compare speedup".format(
                run.instance, baseline.run.name, run.run.name))
    return [[reference / sample for sample in run.samples] for run in runs]


def draw_figure(plt, groups, values, baselines, output, speedup=False):
    columns = min(2, len(groups))
    rows = math.ceil(len(groups) / columns)
    longest = max(len(run.label) for runs in groups.values() for run in runs)
    widest = max(len(runs) for runs in groups.values())
    panel_width = max(7.2, widest * max(1.25, longest * 0.082))
    fig, axes = plt.subplots(rows, columns, squeeze=False,
                             figsize=(panel_width * columns, 4.6 * rows))
    color = "#3978a8"
    randomizer = random.Random(0)

    for ax, (instance, runs), observations, baseline in zip(
            axes.flat, groups.items(), values, baselines):
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
        ax.set_xticklabels([run.label for run in runs], fontsize=10)
        ax.set_ylabel("Speedup (×)" if speedup else "Time (s)")
        title = "Speedup vs {}".format(baseline.label) if speedup else "Runtime"
        ax.set_title("{} · Input {}".format(title, instance), loc="left", pad=16)
        ax.set_xlim(-0.5, len(runs) - 0.5)
        ax.margins(y=0.12)
        ax.set_ylim(bottom=0)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#e9e9e9", linewidth=0.7)
        ax.tick_params(axis="both", length=0, pad=8)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#cccccc")

    for ax in list(axes.flat)[len(groups):]:
        ax.set_visible(False)
    note = "Faint dots: all samples   •   Mean markers + error bars: mean ± 1 sample SD"
    if speedup:
        note += "\nEach sample = baseline mean time / measured time"
    fig.text(0.5, 0.025, note, ha="center", va="bottom", fontsize=9, color="#666666")
    fig.tight_layout(rect=(0, 0.11 if speedup else 0.08, 1, 1), pad=1.5)
    fig.savefig(output, dpi=220, facecolor="white")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=ROOT / "results",
                        help="saved benchmark runs (default: results/ beside this script)")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "plots",
                        help="image directory (default: plots/ beside this script)")
    parser.add_argument("--baseline", help="baseline label or full run directory name; "
                        "defaults to +Baseline, then the earliest run for each input")
    args = parser.parse_args(argv)
    try:
        groups = load_results(args.results_dir)
        baselines = [choose_baseline(runs, args.baseline) for runs in groups.values()]
        speedups = [speedup_samples(runs, baseline)
                    for runs, baseline in zip(groups.values(), baselines)]
        raw_times = [[run.samples for run in runs] for runs in groups.values()]
        figures = (
            ("raw_time.png", raw_times, False),
            ("speedup.png", speedups, True),
        )
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
                for filename, values, speedup in figures:
                    draw_figure(plt, groups, values, baselines, args.output_dir / filename,
                                speedup=speedup)
        for (instance, runs), baseline in zip(groups.items(), baselines):
            print("Input {}: {} runs, {} samples; baseline {}".format(
                instance, len(runs), sum(len(run.samples) for run in runs), baseline.label))
        for filename, _, _ in figures:
            print("Saved {}".format(args.output_dir / filename))
    except (OSError, ValueError) as error:
        parser.exit(1, "error: {}\n".format(error))
    return 0


if __name__ == "__main__":
    sys.exit(main())
