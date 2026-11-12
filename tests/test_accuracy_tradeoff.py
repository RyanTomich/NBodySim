"""Final-position accuracy-tradeoff data and plotting validation."""

from contextlib import redirect_stdout, redirect_stderr
import csv
import io
import json
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bench"))
import accuracy_tradeoff as tradeoff
import run as runner
sys.path.insert(0, str(ROOT / "figures"))
import plot


class TradeoffTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name)
        args = tradeoff.parse_args(["--bodies", "64", "--steps", "2", "--thetas", "0.1", "0.7",
                                    "--results-dir", str(cls.root / "source")])
        with patch.object(tradeoff, "build", return_value=ROOT / "build/bench"), redirect_stdout(io.StringIO()):
            if tradeoff.experiment(args) != 0:
                raise RuntimeError("Small accuracy experiment failed")
        cls.source = next((cls.root / "source").iterdir())

    def setUp(self):
        self.case_tmp = TemporaryDirectory()
        self.addCleanup(self.case_tmp.cleanup)
        self.directory = Path(self.case_tmp.name) / "run"
        shutil.copytree(self.source, self.directory)

    def rows(self, path):
        with path.open(newline="") as source:
            return list(csv.DictReader(source))

    def test_defaults_and_validation(self):
        args = tradeoff.parse_args([])
        self.assertEqual((args.bodies, args.steps, args.samples, args.warmups), (1000, 500, 3, 1))
        self.assertEqual(args.dt, 0.0001)
        self.assertEqual(args.reference_threads, 1)
        self.assertEqual(args.thetas, [0.2, 0.3, 0.5, 0.7])
        self.assertEqual(len(tradeoff.configurations(args.thetas)), 10)
        for argv in (["--samples", "1"], ["--dt", "nan"], ["--bodies", "0"],
                     ["--thetas", "0.1", "0.1"], ["--thetas", "inf"], ["--warmups", "-1"], ["--reference-threads", "0"]):
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                tradeoff.parse_args(argv)

    def test_failed_experiment_retains_status(self):
        args = tradeoff.parse_args(["--results-dir", str(self.directory / "failed")])
        with patch.object(tradeoff, "build", side_effect=RuntimeError("test build failure")), redirect_stdout(io.StringIO()):
            self.assertEqual(tradeoff.experiment(args), 1)
        failed = next((self.directory / "failed").iterdir())
        self.assertEqual(json.loads((failed / "status.json").read_text())["status"], "failed")
        self.assertIn("test build failure", (failed / "error.txt").read_text())

    def test_complete_position_only_data_and_one_figure(self):
        status, rows = plot.load_accuracy_tradeoff(self.directory)
        self.assertEqual(len(rows), 8)
        self.assertFalse((self.directory / "n2").exists())
        reference = json.loads((self.directory / "reference/status.json").read_text())
        self.assertEqual((reference["mode"], reference["samples"], reference["warmups"]),
                         ("reference", 1, 0))
        self.assertFalse((self.directory / "reference/samples.csv").exists())
        for path in (self.directory, self.directory / "accuracy-tradeoff.csv"):
            self.assertEqual([job[0] for job in plot.prepare_figures([path])],
                             ["accuracy-trajectory"])
        for row in rows:
            case = self.directory / row["case"]
            self.assertFalse((case / "force").exists())
            self.assertEqual([error["field"] for error in self.rows(case / "timing/error.csv")], ["position"])
            self.assertNotIn("force_relative_l2", row)
            self.assertNotIn("acceleration_relative_l2", row)
            self.assertNotIn("velocity_relative_l2", row)
            final = runner.read_state(case / "timing/final.txt")
            direct_final = runner.read_state(self.directory / "reference/final.txt")
            self.assertEqual(row["position_relative_l2"], runner.trajectory_errors(final, direct_final)[0]["relative_l2_error"])
            self.assertEqual(len(self.rows(case / "timing/samples.csv")), status["settings"]["samples"])

    def test_parallel_reference_is_one_untimed_call_with_its_own_environment(self):
        args = tradeoff.parse_args(["--steps", "2", "--reference-threads", "4"])
        folder = self.directory / "parallel-reference"
        def produce_reference(*unused, **kwargs):
            shutil.copy2(self.directory / "reference/final.txt", folder / "final.txt")
            return Mock(returncode=0)
        with patch.object(tradeoff.subprocess, "run", side_effect=produce_reference) as call:
            tradeoff.collect_reference(ROOT / "build/bench", self.directory / "input.txt", folder, args)
        call.assert_called_once()
        self.assertEqual(call.call_args.args[0][4:7], ["1", "0", "reference"])
        environment = call.call_args.kwargs["env"]
        self.assertEqual((environment["OMP_NUM_THREADS"], environment["OMP_DYNAMIC"]), ("4", "FALSE"))
        status = json.loads((folder / "status.json").read_text())
        self.assertEqual((status["status"], status["threads"], status["openmp"]), ("ok", 4, True))

    def test_reference_failure_stops_before_timing(self):
        args = tradeoff.parse_args(["--bodies", "8", "--results-dir", str(self.directory / "failed")])
        with patch.object(tradeoff, "build", return_value=ROOT / "build/bench"), \
                patch.object(tradeoff, "collect_reference", side_effect=RuntimeError("reference failed")), \
                patch.object(runner, "collect_timing") as timing, redirect_stdout(io.StringIO()):
            self.assertEqual(tradeoff.experiment(args), 1)
        timing.assert_not_called()

    def test_reject_failed_or_mismatched_reference(self):
        path = self.directory / "reference/status.json"
        original = json.loads(path.read_text())
        for changes in ({"status": "failed"}, {"steps": 999}, {"threads": 4}, {"samples": 3}):
            runner.write_json(path, {**original, **changes})
            with self.assertRaisesRegex(ValueError, "direct reference"):
                plot.load_accuracy_tradeoff(self.directory)

    def test_pareto_ties_and_dominance(self):
        self.assertEqual(plot.pareto_frontier([(1, 4), (2, 2), (4, 1), (2, 3), (2, 2), (1, 5)]), [0, 1, 2, 4])
        self.assertEqual(plot.pareto_frontier([]), [])

    def test_log_figure_uses_correct_units_and_omits_zero(self):
        status, rows = plot.load_accuracy_tradeoff(self.directory)
        metric = "position_relative_l2"
        selected = [dict(row) for row in rows if row["solver"] != "n2"]
        for i, row in enumerate(selected):
            row[metric] = 0.01 * (i + 1)
        selected[0][metric], selected[1][metric] = 0.0, None
        plt, fig, ax = Mock(), Mock(), Mock()
        plt.subplots.return_value = fig, ax
        with redirect_stderr(io.StringIO()) as warnings:
            self.assertIs(plot.draw_accuracy_trajectory(plt, (status, selected)), fig)
        self.assertEqual(warnings.getvalue().count("Omitting"), 2)
        ax.set_xscale.assert_called_once_with("log")
        self.assertEqual(ax.set_yscale.call_args.args, ("log",))
        self.assertEqual(ax.annotate.call_count, len(selected) - 2)
        for call in ax.errorbar.call_args_list:
            for error, milliseconds in zip(*call.args):
                row = next(row for row in selected if row[metric] == error)
                self.assertEqual(milliseconds, 1000 * row["mean_seconds_per_step"])

    def test_minimum_theta_filters_points_and_frontier(self):
        status, rows = plot.load_accuracy_tradeoff(self.directory)
        for index, row in enumerate(rows):
            row["position_relative_l2"] = 0.01 * (index + 1)
        expected = [row for row in rows if row["solver"] != "barnes-hut" or row["theta"] >= 0.2]
        plt, fig, ax = Mock(), Mock(), Mock()
        plt.subplots.return_value = fig, ax
        with patch.object(plot, "pareto_frontier", wraps=plot.pareto_frontier) as frontier:
            plot.draw_accuracy_trajectory(plt, (status, rows), min_theta=0.2)
        frontier.assert_called_once_with([(row["position_relative_l2"], row["mean_seconds_per_step"])
                                          for row in expected])
        self.assertEqual(ax.annotate.call_count, len(expected))
        self.assertEqual(sum(len(call.args[0]) for call in ax.errorbar.call_args_list), len(expected))

    def test_legacy_extra_errors_do_not_require_force_data(self):
        path = self.directory / "accuracy-tradeoff.csv"
        rows = self.rows(path)
        for row in rows:
            row.update(force_relative_l2="unused", velocity_relative_l2="unused", acceleration_relative_l2="unused")
            errors_path = self.directory / row["case"] / "timing/error.csv"
            errors = self.rows(errors_path)
            errors.extend([{**errors[0], "field": field} for field in ("velocity", "acceleration")])
            runner.write_csv(errors_path, list(errors[0]), errors)
        runner.write_csv(path, list(rows[0]), rows)
        self.assertEqual([job[0] for job in plot.prepare_figures([self.directory])], ["accuracy-trajectory"])

    def test_reject_missing_duplicate_and_bad_summary(self):
        path = self.directory / "accuracy-tradeoff.csv"
        rows = self.rows(path)
        for corrupt in (rows[:-1], rows + [rows[-1]],
                        [{**rows[0], "mean_seconds_per_step": "100"}] + rows[1:],
                        [{**rows[0], "position_relative_l2": "nan"}] + rows[1:]):
            runner.write_csv(path, list(rows[0]), corrupt)
            with self.assertRaises(ValueError):
                plot.load_accuracy_tradeoff(self.directory)

    def test_reject_missing_and_nonfinite_samples(self):
        path = self.directory / "barnes-hut-theta-0.1/timing/samples.csv"
        rows = self.rows(path)
        for corrupt in (rows[:-1], [{**rows[0], "seconds": "nan"}] + rows[1:]):
            runner.write_csv(path, list(rows[0]), corrupt)
            with self.assertRaises(ValueError):
                plot.load_accuracy_tradeoff(self.directory)

    def test_reject_failed_or_wrong_timing_policy(self):
        path = self.directory / "barnes-hut-theta-0.1/timing/status.json"
        original = json.loads(path.read_text())
        for changes in ({"status": "failed"}, {"integrator": "other"}, {"steps": 999}):
            runner.write_json(path, {**original, **changes})
            with self.assertRaises(ValueError):
                plot.load_accuracy_tradeoff(self.directory)

    def test_zero_and_undefined_errors_are_retained(self):
        # A single body at the origin has zero position norm; the existing metric convention
        # returns zero for zero/zero, undefined for a nonzero error/zero norm.
        reference = [[1.0] + [0.0] * 9]
        actual = [[1.0, 1.0] + [0.0] * 8]
        self.assertEqual(runner.trajectory_errors(reference, reference)[0]["relative_l2_error"], 0)
        self.assertIsNone(runner.trajectory_errors(actual, reference)[0]["relative_l2_error"])


if __name__ == "__main__":
    unittest.main()
