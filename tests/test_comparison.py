"""Regression checks for the current benchmark matrix and saved-data validation."""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))
from compare_native import configurations_for, samples_for, summarize
from plot_comparison import load_comparison, SERIES
from run import write_csv, write_json


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        write_json(self.root / "status.json", {
            "status": "ok", "settings": {"suite": "current", "bodies": [1000],
            "steps": 10, "dt": 0.0001, "theta": 0.5, "samples": 2}})
        self.rows = []
        for name, threads in SERIES:
            batch = self.root / "1000" / (name + "-" + str(threads) + "threads") / "batch-01"
            batch.mkdir(parents=True)
            write_json(batch / "status.json", {"status": "ok", "requested": threads,
                                               "effective": threads, "samples": 2})
            write_csv(batch / "samples.csv", ["sample", "seconds"],
                      [{"sample": 1, "seconds": 0.01}, {"sample": 2, "seconds": 0.03}])
            self.rows.append({"bodies": 1000, "implementation": name, "requested_threads": threads,
                              "effective_threads": threads, "integrator": "leapfrog", "steps": 10,
                              "dt": 0.0001, "theta": "" if name == "n2" else 0.5,
                              "force_evaluations_per_run": 11, **summarize([0.01, 0.03], 10),
                              "force_relative_l2": 0.001, "force_rms_absolute": 0.0001})
        self.save_rows()

    def save_rows(self):
        write_csv(self.root / "comparison.csv", list(self.rows[0]), self.rows)

    def test_current_matrix_has_only_five_series(self):
        args = SimpleNamespace(suite="current", bodies=[1000, 10000, 100000], threads=[1, 4])
        cases = configurations_for(args)
        self.assertEqual(len(cases), 15)
        self.assertEqual({(name, threads) for _, name, threads in cases}, set(SERIES))

    def test_historical_matrix_unchanged(self):
        args = SimpleNamespace(suite="historical", bodies=[10000], threads=[1, 4, 8])
        self.assertEqual(len(configurations_for(args)), 9)

    def test_explicit_skip_omits_only_largest_n2(self):
        args = SimpleNamespace(suite="current", bodies=[1000, 10000, 100000],
                               threads=[1, 4], skip_large_n2=True)
        cases = configurations_for(args)
        self.assertEqual(len(cases), 14)
        self.assertNotIn((100000, "n2", 1), cases)
        self.assertIn((10000, "n2", 1), cases)

    def test_plot_accepts_only_explicitly_skipped_point(self):
        self.rows.pop(0)
        self.save_rows()
        path = self.root / "status.json"
        status = json.loads(path.read_text())
        status["settings"]["skip_large_n2"] = True
        write_json(path, status)
        self.assertEqual(len(load_comparison(self.root)[1]), 4)

    def test_reduced_samples_apply_only_to_largest_n2(self):
        args = SimpleNamespace(suite="current", bodies=[1000, 10000, 100000],
                               threads=[1, 4], batch_size=10, large_n2_samples=2)
        for key in configurations_for(args):
            self.assertEqual(samples_for(args, key), 2 if key == (100000, "n2", 1) else 10)

    def test_plot_accepts_declared_mixed_sample_counts(self):
        path = self.root / "status.json"
        status = json.loads(path.read_text())
        status["settings"].update(samples=3, large_n2_samples=2)
        write_json(path, status)
        for row in self.rows:
            if row["implementation"] == "n2":
                continue
            case = self.root / "1000" / (row["implementation"] + "-" + str(row["requested_threads"]) + "threads") / "batch-01"
            with (case / "samples.csv").open("a") as output:
                output.write("3,0.02\n")
            info = json.loads((case / "status.json").read_text())
            info["samples"] = 3
            write_json(case / "status.json", info)
            row.update(summarize([0.01, 0.03, 0.02], 10))
        self.save_rows()
        _, rows = load_comparison(self.root)
        self.assertEqual([int(row["samples"]) for row in rows], [2, 3, 3, 3, 3])

    def test_mean_and_sample_standard_deviation(self):
        _, rows = load_comparison(self.root)
        self.assertAlmostEqual(rows[0]["mean_ms_per_step"], 2)
        self.assertAlmostEqual(rows[0]["stddev_ms_per_step"], 2 ** 0.5)

    def test_missing_point_rejected(self):
        self.rows.pop()
        self.save_rows()
        with self.assertRaisesRegex(ValueError, "Missing configurations"):
            load_comparison(self.root)

    def test_duplicate_point_rejected(self):
        self.rows.append(self.rows[0])
        self.save_rows()
        with self.assertRaisesRegex(ValueError, "duplicate"):
            load_comparison(self.root)

    def test_wrong_normalization_rejected(self):
        self.rows[0]["mean_ms_per_step"] *= 10
        self.save_rows()
        with self.assertRaisesRegex(ValueError, "Summary does not match"):
            load_comparison(self.root)

    def test_incomplete_run_rejected(self):
        path = self.root / "status.json"
        status = json.loads(path.read_text())
        status["status"] = "running"
        write_json(path, status)
        with self.assertRaisesRegex(ValueError, "completed"):
            load_comparison(self.root)

    def test_invalid_or_missing_samples_rejected(self):
        path = self.root / "1000/n2-1threads/batch-01/samples.csv"
        for contents in ("sample,seconds\n1,nan\n2,0.03\n", "sample,seconds\n1,0.01\n"):
            with self.subTest(contents=contents):
                path.write_text(contents)
                with self.assertRaisesRegex(ValueError, "Invalid raw"):
                    load_comparison(self.root)


if __name__ == "__main__":
    unittest.main()
