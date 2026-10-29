import contextlib
import importlib.util
import io
import json
from pathlib import Path
import statistics
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("plot_results", ROOT / "plot_results.py")
plotter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plotter)


class PlotResultsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / "build")
        self.addCleanup(self.temp.cleanup)
        self.results = Path(self.temp.name)

    def save_run(self, name, samples, instance="0001-1000", status="ok", steps=1):
        timing = self.results / name / instance / "timing"
        timing.mkdir(parents=True)
        (timing.parent / "input.txt").write_text("1\n1 0 0 0 0 0 0 0 0 0\n")
        (timing / "status.json").write_text(json.dumps({"status": status, "steps": steps}))
        (timing / "samples.csv").write_text("sample,seconds\n" + "".join(
            "{},{}\n".format(index, sample) for index, sample in enumerate(samples, 1)))
        return timing

    def test_reads_all_samples_and_orders_by_timestamp_not_commit_hash(self):
        self.save_run("aaaaaaa_20261023-160502.425884Z_timing_+Linking", [2, 3, 4])
        self.save_run("fffffff_20261023-154013.688489Z_timing_+Baseline", [4, 6])
        self.save_run("bbbbbbb_20261023-170502.425884Z_all_+Local_Sum", [1, 2, 3, 4],
                      instance="0002-1000")
        runs = plotter.load_results(self.results)["1000"]
        self.assertEqual([run.label for run in runs], ["+Baseline", "+Linking", "+Local_Sum"])
        self.assertEqual([run.samples for run in runs], [[4, 6], [2, 3, 4], [1, 2, 3, 4]])

    def test_solver_labels_preserve_baseline_selection_and_optional_solver_metadata(self):
        self.save_run("aaaaaaa_20261023-150000.000000Z_timing_+Earlier", [1, 2])
        second = self.save_run("bbbbbbb_20261023-160000.000000Z_timing_+Baseline", [4, 8])
        third = self.save_run("ccccccc_20261023-170000.000000Z_timing_+Changed", [2, 3])
        (second / "status.json").write_text(json.dumps({"status": "ok", "steps": 1, "solver": "direct"}))
        (third / "status.json").write_text(json.dumps({"status": "ok", "steps": 1, "solver": "approximate"}))
        runs = plotter.load_results(self.results)["1000"]
        self.assertEqual([run.label for run in runs],
                         ["+Earlier", "+Baseline [direct]", "+Changed [approximate]"])
        self.assertIs(plotter.choose_baseline(runs), runs[1])
        self.assertIs(plotter.choose_baseline(runs, "+Baseline"), runs[1])
        self.assertIs(plotter.choose_baseline(runs, "+Changed [approximate]"), runs[2])

    def test_baseline_selection_and_per_sample_speedup(self):
        self.save_run("aaaaaaa_20261023-150000.000000Z_timing_+Earlier", [1, 2])
        self.save_run("bbbbbbb_20261023-160000.000000Z_timing_+Baseline", [4, 8])
        runs = plotter.load_results(self.results)["1000"]
        baseline = plotter.choose_baseline(runs)
        self.assertEqual(baseline.label, "+Baseline")
        ratios = plotter.speedup_samples(runs, baseline)
        self.assertEqual(ratios, [[6, 3], [1.5, 0.75]])
        self.assertEqual(statistics.mean(ratios[0]), 4.5)
        self.assertAlmostEqual(statistics.stdev(ratios[0]), 4.5 ** 0.5)
        self.assertIs(plotter.choose_baseline(runs, "+Earlier"), runs[0])
        self.assertIs(plotter.choose_baseline(runs, runs[0].run.name), runs[0])
        self.assertIs(plotter.choose_baseline(runs[:1]), runs[0])
        with self.assertRaisesRegex(ValueError, "matches 0 runs"):
            plotter.choose_baseline(runs, "+Missing")

    def test_step_counts_normalize_samples_and_speedups(self):
        self.save_run("aaaaaaa_20261023-150000.000000Z_timing_+Baseline", [2, 2], steps=10)
        self.save_run("bbbbbbb_20261023-160000.000000Z_timing_+Longer", [20, 20], steps=100)
        self.save_run("ccccccc_20261023-170000.000000Z_timing_+Faster", [10, 10], steps=100)
        runs = plotter.load_results(self.results)["1000"]
        self.assertEqual([run.samples for run in runs], [[0.2, 0.2], [0.2, 0.2], [0.1, 0.1]])
        self.assertEqual(plotter.speedup_samples(runs, runs[0]), [[1, 1], [1, 1], [2, 2]])

    def test_skips_legacy_results_without_changing_files(self):
        self.save_run("aaaaaaa_20261023-160000.000000Z_timing", [1, 2], steps=10)
        legacy = self.save_run("bbbbbbb_20261023-170000.000000Z_timing", [1, 2])
        (legacy / "status.json").write_text(json.dumps({"status": "ok"}))
        before = {path: path.read_bytes() for path in legacy.iterdir()}
        errors = io.StringIO()
        with contextlib.redirect_stderr(errors):
            runs = plotter.load_results(self.results)["1000"]
        self.assertEqual(len(runs), 1)
        self.assertIn("no recorded step count; cannot calculate time per step", errors.getvalue())
        self.assertEqual(before, {path: path.read_bytes() for path in legacy.iterdir()})

    def test_rejects_present_but_invalid_step_counts(self):
        timing = self.save_run("aaaaaaa_20261023-160000.000000Z_timing", [1, 2])
        for steps in (0, -1, 1.5, 2.0, "2", None, True):
            (timing / "status.json").write_text(json.dumps({"status": "ok", "steps": steps}))
            with self.subTest(steps=steps), self.assertRaisesRegex(ValueError, "positive integer step count"):
                plotter.load_results(self.results)

    def test_plot_displays_milliseconds_per_step_and_run_variation(self):
        self.save_run("aaaaaaa_20261023-160000.000000Z_timing", [0.002, 0.004, 0.006], steps=2)
        runs = plotter.load_results(self.results)["1000"]
        plt = mock.MagicMock()
        fig, ax = mock.MagicMock(), mock.MagicMock()
        plt.subplots.return_value = fig, ax
        plotter.draw_figure(plt, runs, [runs[0].samples], runs[0], self.results / "runtime.png")
        ax.set_ylabel.assert_called_once_with("Time per step (ms)")
        self.assertEqual(ax.scatter.call_args.args[1], [1, 2, 3])
        self.assertEqual(ax.errorbar.call_args.args[1], [2])
        self.assertEqual(ax.errorbar.call_args.kwargs["yerr"], [1])
        self.assertIn("SD across runs", fig.text.call_args.args[2])
        self.assertEqual(runs[0].samples, [0.001, 0.002, 0.003])  # Conversion is display-only.

        ax.reset_mock()
        ratios = plotter.speedup_samples(runs, runs[0])
        plotter.draw_figure(plt, runs, ratios, runs[0], self.results / "speedup.png", speedup=True)
        ax.set_ylabel.assert_called_once_with("Speedup (×)")
        self.assertEqual(ax.scatter.call_args.args[1], ratios[0])

    def test_failed_runs_are_reported_and_inputs_are_kept_separate(self):
        name = "aaaaaaa_20261023-160000.000000Z_timing_+Baseline"
        self.save_run(name, [1, 2])
        self.save_run(name, [3, 4], instance="0002-two_body")
        self.save_run("bbbbbbb_20261023-170000.000000Z_timing_+Failed", [100], status="failed")
        errors = io.StringIO()
        with contextlib.redirect_stderr(errors):
            groups = plotter.load_results(self.results)
        self.assertEqual(list(groups), ["1000", "two_body"])
        self.assertEqual([len(runs) for runs in groups.values()], [1, 1])
        self.assertIn("+Failed", errors.getvalue())

    def test_creates_two_figures_per_instance_with_its_own_runs_and_baseline(self):
        baseline_name = "aaaaaaa_20261023-160000.000000Z_timing_+Baseline"
        changed_name = "bbbbbbb_20261023-170000.000000Z_timing_+Changed"
        self.save_run(baseline_name, [4, 6])
        self.save_run(changed_name, [2, 3], instance="0002-1000")
        self.save_run(baseline_name, [40, 60], instance="0002-100000")
        self.save_run(changed_name, [20, 30], instance="0001-100000")
        output = self.results / "plots"
        matplotlib = mock.MagicMock()
        messages = io.StringIO()
        # Exercise CLI grouping without making Matplotlib a test-suite dependency.
        with mock.patch.dict("sys.modules", {
                "matplotlib": matplotlib, "matplotlib.pyplot": matplotlib.pyplot}), \
                mock.patch.object(plotter, "draw_figure") as draw, \
                contextlib.redirect_stdout(messages):
            self.assertEqual(plotter.main([
                "--results-dir", str(self.results), "--output-dir", str(output)]), 0)

        self.assertEqual(draw.call_count, 4)
        calls = {call.args[4].name: call for call in draw.call_args_list}
        self.assertEqual(set(calls), {
            "raw_time_1000.png", "speedup_1000.png",
            "raw_time_100000.png", "speedup_100000.png"})
        for instance, raw, reference in (
                ("1000", [[4, 6], [2, 3]], 5),
                ("100000", [[40, 60], [20, 30]], 50)):
            for kind in ("raw_time", "speedup"):
                call = calls["{}_{}.png".format(kind, instance)]
                _, runs, observations, baseline, path = call.args
                self.assertEqual([run.instance for run in runs], [instance, instance])
                self.assertEqual([run.label for run in runs], ["+Baseline", "+Changed"])
                self.assertIs(baseline, runs[0])
                self.assertEqual(statistics.mean(baseline.samples), reference)
                expected = raw if kind == "raw_time" else [
                    [reference / sample for sample in samples] for samples in raw]
                self.assertEqual(observations, expected)
                self.assertEqual(call.kwargs["speedup"], kind == "speedup")
                self.assertEqual(path.parent, output)
                self.assertIn("Saved {}".format(path), messages.getvalue())

    def test_unfinished_runs_are_skipped_but_successful_runs_require_samples(self):
        self.save_run("aaaaaaa_20261023-160000.000000Z_timing_+Baseline", [2, 4])
        unfinished = self.save_run("bbbbbbb_20261023-170000.000000Z_timing", [1, 2])
        (unfinished / "status.json").unlink()
        errors = io.StringIO()
        with contextlib.redirect_stderr(errors):
            runs = plotter.load_results(self.results)["1000"]
        self.assertEqual(len(runs), 1)
        self.assertIn("no timing status yet", errors.getvalue())
        (unfinished / "status.json").write_text(json.dumps({"status": "ok", "steps": 1}))
        (unfinished / "samples.csv").unlink()
        with self.assertRaisesRegex(ValueError, "samples.csv"):
            plotter.load_results(self.results)

    def test_rejects_invalid_samples_and_incomparable_inputs(self):
        first = self.save_run("aaaaaaa_20261023-160000.000000Z_timing_+Baseline", [2, 4])
        second = self.save_run("bbbbbbb_20261023-170000.000000Z_timing_+Changed", [1, 2])
        runs = plotter.load_results(self.results)["1000"]
        (second.parent / "input.txt").write_text("0\n")
        with self.assertRaisesRegex(ValueError, "differs"):
            plotter.speedup_samples(runs, runs[0])
        runs[0].samples = [0, 1]
        with self.assertRaisesRegex(ValueError, "positive timings"):
            plotter.speedup_samples(runs, runs[0])
        (first / "samples.csv").write_text("sample,seconds\n1,nan\n2,1\n")
        with self.assertRaisesRegex(ValueError, "finite"):
            plotter.load_results(self.results)


if __name__ == "__main__":
    unittest.main()
