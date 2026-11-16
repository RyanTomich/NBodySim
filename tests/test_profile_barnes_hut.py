"""Check parallel profiling accounting independently of the measured machine."""

import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))
from profile_barnes_hut import COUNTERS, TIMERS, WORKER_TIMERS, phase_samples, thread_environment


class PhaseAccountingTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.delta = {key: 0 for key in (*COUNTERS, *TIMERS, *WORKER_TIMERS, "evaluation_seconds")}
        self.delta.update(calls=3, tree_builds=1, tree_reuses=2, membership_checks=2, worker_slots=12,
                          membership_seconds=0.05, construction_seconds=0.1, summaries_seconds=0.1,
                          bounds_seconds=0.05, traversal_seconds=0.15, force_seconds=0.3,
                          bounds_worker_seconds=0.2, traversal_worker_seconds=0.6, force_worker_seconds=1.2,
                          evaluation_seconds=0.6)

    def save(self, legacy=False):
        snapshots = []
        delta = {key: value for key, value in self.delta.items()
                 if not legacy or key not in (*WORKER_TIMERS, "worker_slots")}
        for sample in (1, 2):
            for boundary, factor in (("before", sample - 1), ("after", sample)):
                value = {key: factor * amount for key, amount in delta.items()}
                if not legacy:
                    value.update(min_threads=4 if factor else 0, max_threads=4 if factor else 0)
                value.update(sample=sample, boundary=boundary)
                snapshots.append("phase_sample: " + json.dumps(value))
        (self.folder / "stderr.log").write_text("\n".join(snapshots) + "\n")
        (self.folder / "samples.csv").write_text("sample,seconds\n1,1\n2,1\n")

    def test_parallel_phases_count_overlap_once(self):
        self.save()
        rows = phase_samples(self.folder, 2, 2, threads=4)
        for row in rows:
            self.assertEqual(row["threads"], 4)
            self.assertAlmostEqual(row["evaluation_other_seconds"], 0.1)
            self.assertAlmostEqual(row["simulation_other_seconds"], 0.15)
            self.assertAlmostEqual(sum(row[key] for key in TIMERS)
                                   + row["evaluation_other_seconds"] + row["simulation_other_seconds"], 1)

    def test_reject_requested_team_mismatch(self):
        self.save()
        with self.assertRaisesRegex(ValueError, "Requested 1 threads but measured 4"):
            phase_samples(self.folder, 2, 2, threads=1)

    def test_reject_unscaled_worker_time(self):
        self.delta["force_seconds"] = self.delta["force_worker_seconds"]
        self.save()
        with self.assertRaisesRegex(ValueError, "normalized"):
            phase_samples(self.folder, 2, 2, threads=4)

    def test_reject_missing_team_slots(self):
        self.delta["worker_slots"] = 9
        self.save()
        with self.assertRaisesRegex(ValueError, "team size"):
            phase_samples(self.folder, 2, 2, threads=4)

    def test_historical_serial_samples_remain_readable(self):
        self.save(legacy=True)
        rows = phase_samples(self.folder, 2, 2, threads=1)
        self.assertEqual(rows[0]["worker_slots"], 3)
        self.assertEqual(rows[0]["force_worker_seconds"], rows[0]["force_seconds"])

    def test_thread_environment_restored_after_failure(self):
        before = dict(os.environ)
        with self.assertRaises(RuntimeError):
            with thread_environment(4):
                self.assertEqual(os.environ["OMP_NUM_THREADS"], "4")
                self.assertEqual(os.environ["OMP_DYNAMIC"], "FALSE")
                raise RuntimeError("deliberate")
        self.assertEqual(dict(os.environ), before)


if __name__ == "__main__":
    unittest.main()
