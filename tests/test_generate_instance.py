"""Seeded uniform and clustered initial conditions."""

from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
import io
import math
from pathlib import Path
import statistics
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import generate_instance as generator


class GeneratorTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.output = self.root / "input.txt"

    def run_generator(self, *arguments):
        argv = ["generate_instance.py", "--n", "1003", "--w-min", "1", "--w-max", "100",
                "--v-max", "0", "--a-max", "0", "--seed", "42", *arguments]
        with patch.object(sys, "argv", argv), \
                patch.object(generator, "__file__", str(self.root / "generate_instance.py")), \
                redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            generator.main()

    def generate(self, *arguments):
        self.run_generator("--output", str(self.output), *arguments)
        lines = self.output.read_text().splitlines()
        self.assertEqual(int(lines[0]), len(lines) - 1)
        return [list(map(float, line.split())) for line in lines[1:]]

    def test_clusters_are_balanced_gaussians_at_cube_corners(self):
        rows = self.generate("--distribution", "clusters")
        centers = [(-1, -1, -1), (1, -1, -1), (-1, 1, -1), (1, 1, -1),
                   (-1, -1, 1), (1, -1, 1), (-1, 1, 1), (1, 1, 1)]
        assignments = [tuple(1 if coordinate > 0 else -1 for coordinate in row[1:4])
                       for row in rows]
        self.assertEqual(assignments, [centers[i % 8] for i in range(1003)])
        self.assertEqual([Counter(assignments)[center] for center in centers],
                         [126, 126, 126, 125, 125, 125, 125, 125])
        self.assertEqual(len({tuple(row[1:4]) for row in rows}), len(rows))
        for row in rows:
            self.assertEqual(len(row), 10)
            self.assertTrue(all(math.isfinite(value) for value in row))
            self.assertTrue(1 <= row[0] <= 100)
            self.assertEqual(row[4:], [0] * 6)
        for axis in range(3):
            offsets = [row[axis + 1] - center[axis] for row, center in zip(rows, assignments)]
            self.assertAlmostEqual(statistics.mean(offsets), 0, delta=0.03)
            self.assertAlmostEqual(statistics.stdev(offsets), 0.25, delta=0.03)
        self.assertTrue(any(abs(value) > 1 for row in rows for value in row[1:4]))

    def test_seed_reproducibility_and_cluster_bounds(self):
        first = self.generate("--distribution", "clusters")
        self.assertEqual(first, self.generate("--distribution", "clusters", "--L", "0"))
        self.assertNotEqual(first, self.generate("--distribution", "clusters", "--seed", "43"))

    def test_uniform_default_and_vector_options(self):
        options = ("--L", "10", "--v-max", "3", "--a-max", "2")
        rows = self.generate(*options)
        self.assertEqual(rows, self.generate(*options, "--distribution", "uniform"))
        self.assertEqual(rows[0][:4], [64.30325304733049, -9.49978489554666,
                                      -4.499413632617614, -5.5357852370235445])
        for row in rows:
            self.assertTrue(all(-10 <= value <= 10 for value in row[1:4]))
            self.assertLessEqual(math.hypot(*row[4:7]), 3)
            self.assertLessEqual(math.hypot(*row[7:10]), 2)

    def test_default_filenames_keep_distributions_separate(self):
        self.run_generator("--L", "10")
        uniform = self.root / "bench/data/instance-n1003-seed42.txt"
        before = uniform.read_bytes()
        self.run_generator("--distribution", "clusters")
        self.assertEqual(uniform.read_bytes(), before)
        self.assertTrue((self.root / "bench/data/instance-clusters-n1003-seed42.txt").is_file())

    def test_invalid_arguments_and_empty_cluster_input(self):
        for arguments in ((), ("--L", "0"), ("--L", "nan"),
                          ("--distribution", "unknown"), ("--distribution", "clusters", "--n", "-1")):
            with self.subTest(arguments=arguments), self.assertRaises(SystemExit):
                self.generate(*arguments)
        self.assertEqual(self.generate("--distribution", "clusters", "--n", "0"), [])

    def test_clusters_still_reject_nonfinite_and_duplicate_positions(self):
        args = SimpleNamespace(output=self.output, n=2, distribution="clusters",
                               w_min=1, w_max=100, v_max=0, a_max=0)
        for body in ("1 nan 1 1 0 0 0 0 0 0\n", "1 1 1 1 0 0 0 0 0 0\n"):
            self.output.write_text("2\n" + body * 2)
            with self.assertRaises(ValueError):
                generator.verify_instance(args)


if __name__ == "__main__":
    unittest.main()
