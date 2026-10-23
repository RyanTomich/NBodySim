import contextlib
import importlib.util
import io
import json
from pathlib import Path
import statistics
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("plot_results", ROOT / "plot_results.py")
plotter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plotter)


class PlotResultsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / "build")
        self.addCleanup(self.temp.cleanup)
        self.results = Path(self.temp.name)

    def save_run(self, name, samples, instance="0001-1000", status="ok"):
        timing = self.results / name / instance / "timing"
        timing.mkdir(parents=True)
        (timing.parent / "input.txt").write_text("1\n1 0 0 0 0 0 0 0 0 0\n")
        (timing / "status.json").write_text(json.dumps({"status": status}))
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
