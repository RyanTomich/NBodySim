"""Public CLI trajectory regressions; optionally compare a saved baseline binary."""

import math
import os
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[2]


class TrajectoryTests(unittest.TestCase):
    def read_state(self, text):
        lines = text.splitlines()
        rows = [tuple(map(float, line.split())) for line in lines[1:]]
        self.assertEqual(len(rows), int(lines[0]))
        for row in rows:
            self.assertEqual(len(row), 10)
            self.assertTrue(all(math.isfinite(value) for value in row))
        return rows

    def run_solver(self, binary, options):
        result = subprocess.run(
            [str(binary), *options, str(ROOT / "tests/data/1000.txt"), "0.0001", "10"],
            check=True, capture_output=True, text=True)
        return self.read_state(result.stdout)

    def errors(self, actual, reference):
        self.assertEqual([row[0] for row in actual], [row[0] for row in reference])
        for start in (1, 4, 7):
            error = math.sqrt(math.fsum(
                (a[k] - b[k]) ** 2 for a, b in zip(actual, reference)
                for k in range(start, start + 3)))
            norm = math.sqrt(math.fsum(
                b[k] ** 2 for b in reference for k in range(start, start + 3)))
            yield error, norm

    def test_backend_trajectories(self):
        binary = ROOT / "build/app"
        original = self.read_state((ROOT / "tests/data/1000.txt").read_text())
        direct = self.run_solver(binary, ["--solver", "direct"])
        # These are regression budgets for this fixture and short trajectory.
        # Order four has different tree settings: do not assume monotonic error.
        variants = [
            ("direct", ["--solver", "direct"], (0, 0, 0)),
            ("barnes-hut exact", ["--solver", "barnes-hut", "--theta", "0"],
             (1e-11, 1e-11, 1e-11)),
            ("barnes-hut default", ["--solver", "barnes-hut"],
             (5e-8, 7e-4, 3.5e-3)),
            ("fmm default", ["--solver", "fmm"], (2e-9, 2e-5, 1e-4)),
        ]
        acceleration_budgets = (0.04, 0.005, 0.001, 1e-4, 1e-4, 5e-5)
        for order, budget in enumerate(acceleration_budgets, 1):
            limits = (2e-9, 2e-5, 1e-4) if order == 4 else (1e-6, 0.01, budget)
            variants.append((f"fmm order {order}", ["--solver", "fmm", "--order", str(order)], limits))

        baseline = os.environ.get("NBODY_BASELINE_APP")
        for name, options, limits in variants:
            with self.subTest(solver=name):
                actual = direct if name == "direct" else self.run_solver(binary, options)
                self.assertEqual([row[0] for row in actual], [row[0] for row in original])
                for (error, norm), limit in zip(self.errors(actual, direct), limits):
                    self.assertLessEqual(error, limit * norm)
                # A supplied pre-refactor binary validates each approximation
                # against its own previous result, not just against direct.
                if baseline:
                    previous = self.run_solver(Path(baseline).resolve(), options)
                    for error, norm in self.errors(actual, previous):
                        self.assertLessEqual(error, 1e-12 * (1 + norm))


if __name__ == "__main__":
    unittest.main()
