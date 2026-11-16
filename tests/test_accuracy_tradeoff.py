"""Accuracy-tradeoff benchmark data validation."""

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
        self.assertEqual(args.distribution, "uniform")
        self.assertEqual(args.thetas, [0.2, 0.3, 0.5, 0.7])
        self.assertEqual(len(tradeoff.configurations(args.thetas)), 10)
        for argv in (["--samples", "1"], ["--dt", "nan"], ["--bodies", "0"],
                     ["--thetas", "0.1", "0.1"], ["--thetas", "inf"], ["--warmups", "-1"], ["--reference-threads", "0"]):
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                tradeoff.parse_args(argv)

    def test_cluster_experiment_saves_distribution_and_rest_state(self):
        args = tradeoff.parse_args(["--distribution", "clusters", "--bodies", "64",
                                    "--steps", "2", "--samples", "2", "--dt", "0.000001",
                                    "--results-dir", str(self.directory / "clusters")])
        with patch.object(tradeoff, "build", return_value=ROOT / "build/bench"), redirect_stdout(io.StringIO()):
            self.assertEqual(tradeoff.experiment(args), 0)
        directory = next((self.directory / "clusters").iterdir())
        status = json.loads((directory / "status.json").read_text())
        self.assertEqual(status["status"], "ok")
        self.assertEqual(status["settings"]["generator"],
                         {"distribution": "clusters", "L": None, "w_min": 1,
                          "w_max": 100, "v_max": 0, "a_max": 0})
        initial = runner.read_state(directory / "input.txt")
        self.assertTrue(all(row[4:] == [0] * 6 for row in initial))
        self.assertTrue(all(abs(abs(value) - 1) < 2 for row in initial for value in row[1:4]))
        self.assertEqual(len(self.rows(directory / "accuracy-tradeoff.csv")), 10)

    def test_failed_experiment_retains_status(self):
        args = tradeoff.parse_args(["--results-dir", str(self.directory / "failed")])
        with patch.object(tradeoff, "build", side_effect=RuntimeError("test build failure")), redirect_stdout(io.StringIO()):
            self.assertEqual(tradeoff.experiment(args), 1)
        failed = next((self.directory / "failed").iterdir())
        self.assertEqual(json.loads((failed / "status.json").read_text())["status"], "failed")
        self.assertIn("test build failure", (failed / "error.txt").read_text())

    def test_complete_distribution_data(self):
        status = json.loads((self.directory / "status.json").read_text())
        self.assertEqual(status["status"], "ok")
        rows = self.rows(self.directory / "accuracy-tradeoff.csv")
        self.assertEqual(len(rows), 8)
        self.assertFalse((self.directory / "n2").exists())
        reference = json.loads((self.directory / "reference/status.json").read_text())
        self.assertEqual((reference["mode"], reference["samples"], reference["warmups"]),
                         ("reference", 1, 0))
        self.assertFalse((self.directory / "reference/samples.csv").exists())
        for path in (self.directory, self.directory / "accuracy-tradeoff.csv"):
            with self.assertRaisesRegex(ValueError, "use --compare-states"):
                plot.prepare_figures([path])
        for row in rows:
            case = self.directory / row["case"]
            self.assertFalse((case / "force").exists())
            self.assertEqual([error["field"] for error in self.rows(case / "timing/error.csv")], ["position"])
            self.assertNotIn("force_relative_l2", row)
            self.assertNotIn("acceleration_relative_l2", row)
            self.assertNotIn("velocity_relative_l2", row)
            final = runner.read_state(case / "timing/final.txt")
            direct_final = runner.read_state(self.directory / "reference/final.txt")
            self.assertEqual(float(row["position_relative_l2"]), runner.trajectory_errors(final, direct_final)[0]["relative_l2_error"])
            distances = runner.distribution_errors(final, direct_final)
            saved = self.rows(case / "timing/distribution-error.csv")
            self.assertEqual([(item["field"], item["component"], float(item["ks_distance"])) for item in saved],
                             [(item["field"], item["component"], item["ks_distance"]) for item in distances])
            for key, value in runner.distribution_summary(distances).items():
                self.assertEqual(float(row[key]), value)
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

    def test_zero_and_undefined_errors_are_retained(self):
        # A single body at the origin has zero position norm; the existing metric convention
        # returns zero for zero/zero, undefined for a nonzero error/zero norm.
        reference = [[1.0] + [0.0] * 9]
        actual = [[1.0, 1.0] + [0.0] * 8]
        self.assertEqual(runner.trajectory_errors(reference, reference)[0]["relative_l2_error"], 0)
        self.assertIsNone(runner.trajectory_errors(actual, reference)[0]["relative_l2_error"])


if __name__ == "__main__":
    unittest.main()
