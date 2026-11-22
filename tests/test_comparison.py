"""Regression checks for REBOUND checkout and saved comparison data."""

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))
import compare_rebound
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
from plot import load_comparison, SERIES
from run import write_csv, write_json


class ReboundCheckoutTests(unittest.TestCase):
    def test_initializes_missing_checkout_and_verifies_pin(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            header = root / "third_party/rebound/src/rebound.h"
            header.parent.mkdir(parents=True)
            header.touch()
            with patch.object(compare_rebound, "ROOT", root), \
                 patch.object(compare_rebound.subprocess, "check_output",
                              side_effect=["-abc third_party/rebound", " abc third_party/rebound"]), \
                 patch.object(compare_rebound, "run_logged", return_value=0) as run:
                self.assertEqual(compare_rebound.prepare_rebound(root), "abc")
                run.assert_called_once_with(
                    ["git", "submodule", "update", "--init", "--", "third_party/rebound"],
                    root, stdout="submodule.log", stderr=None)

    def test_leaves_different_checkout_in_place(self):
        with patch.object(compare_rebound.subprocess, "check_output",
                          return_value="+def third_party/rebound"), \
             patch.object(compare_rebound, "run_logged") as run:
            with self.assertRaisesRegex(RuntimeError, "pinned REBOUND"):
                compare_rebound.prepare_rebound(Path("."))
            run.assert_not_called()


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
                              "force_evaluations_per_run": 11, "samples": 2,
                              "mean_ms_per_step": 2, "median_ms_per_step": 2,
                              "stddev_ms_per_step": 2 ** 0.5,
                              "force_relative_l2": 0.001, "force_rms_absolute": 0.0001})
        self.save_rows()

    def save_rows(self):
        write_csv(self.root / "comparison.csv", list(self.rows[0]), self.rows)

    def test_plot_accepts_only_explicitly_skipped_point(self):
        self.rows.pop(0)
        self.save_rows()
        path = self.root / "status.json"
        status = json.loads(path.read_text())
        status["settings"]["skip_large_n2"] = True
        write_json(path, status)
        self.assertEqual(len(load_comparison(self.root)[1]), 4)

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
            row.update(samples=3, stddev_ms_per_step=1)
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
