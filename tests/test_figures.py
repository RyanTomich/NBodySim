"""Saved-data routing and phase accounting for the shared figure generator."""

import csv
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
from plot import GROUPS, prepare_figures, read_shares


class FigureTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.run = self.root / "phase-experiment"
        self.run.mkdir()
        self.path = self.run / "phase-summary.csv"
        self.rows = [{"bodies": 100, "field": field, "mean_ms_per_step": i + 1}
                     for i, field in enumerate(field for _, fields, _, _ in GROUPS for field in fields)]
        self.save()

    def save(self):
        with self.path.open("w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(self.rows[0]))
            writer.writeheader()
            writer.writerows(self.rows)

    def test_phase_shares_preserve_total_and_group_overhead(self):
        shares = read_shares(self.path).shares[(100, 1)]
        self.assertAlmostEqual(sum(shares), 100)
        self.assertAlmostEqual(shares[0], 100 / 36)
        self.assertAlmostEqual(shares[-1], 100 * (6 + 7 + 8) / 36)

    def test_reject_incomplete_duplicate_and_nonfinite_phases(self):
        original = list(self.rows)
        for rows in (original[:-1], original + original[:1],
                     [{**original[0], "mean_ms_per_step": "nan"}] + original[1:]):
            with self.subTest(rows=rows):
                self.rows = rows
                self.save()
                with self.assertRaises(ValueError):
                    read_shares(self.path)

    def test_one_and_four_threads_share_one_figure(self):
        original = list(self.rows)
        self.rows = [{**row, "threads": threads, "phase_time_basis": "serial_wall_and_mean_worker_elapsed",
                      "mean_ms_per_step": row["mean_ms_per_step"] / threads}
                     for threads in (1, 4) for row in original]
        self.save()
        profile = read_shares(self.path)
        self.assertEqual(set(profile.shares), {(100, 1), (100, 4)})
        self.assertEqual(profile.total_ms[(100, 1)] / profile.total_ms[(100, 4)], 4)
        self.assertEqual(len(prepare_figures([self.path])), 1)
        self.rows[0]["phase_time_basis"] = "wrong"
        # Only parallel data requires explicit worker normalization.
        self.rows[-1]["phase_time_basis"] = "wrong"
        self.save()
        with self.assertRaisesRegex(ValueError, "normalized"):
            read_shares(self.path)

    def test_reject_mismatched_thread_workloads(self):
        original = list(self.rows)
        self.rows = [{**row, "threads": threads, "phase_time_basis": "serial_wall_and_mean_worker_elapsed",
                      "dt": 0.001 * threads}
                     for threads in (1, 4) for row in original]
        self.save()
        with self.assertRaisesRegex(ValueError, "different workloads"):
            read_shares(self.path)

    def test_file_and_directory_route_to_same_figure(self):
        for path in (self.path, self.run):
            jobs = prepare_figures([path])
            self.assertEqual(jobs[0][0], "phase-shares")
            self.assertEqual(len(jobs), 1)

    def test_output_names_are_independent_of_dataset_folder(self):
        other = self.root / "another-experiment"
        other.mkdir()
        (other / self.path.name).write_bytes(self.path.read_bytes())
        self.assertEqual(prepare_figures([self.run])[0][0], prepare_figures([other])[0][0])
        with self.assertRaisesRegex(ValueError, "one dataset"):
            prepare_figures([self.run, other])

    def test_reject_duplicate_data_and_output_name_collisions(self):
        with self.assertRaisesRegex(ValueError, "more than once"):
            prepare_figures([self.run, self.path])
        other = self.root / "nested" / self.run.name
        other.mkdir(parents=True)
        (other / self.path.name).write_bytes(self.path.read_bytes())
        with self.assertRaisesRegex(ValueError, "collision"):
            prepare_figures([self.run, other])

    def test_reject_unsupported_and_missing_data(self):
        unsupported = self.run / "other.csv"
        unsupported.write_text("unrelated data")
        for path in (unsupported, self.root / "missing"):
            with self.subTest(path=path), self.assertRaises(ValueError):
                prepare_figures([path])


if __name__ == "__main__":
    unittest.main()
