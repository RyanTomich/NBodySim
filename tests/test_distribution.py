"""Final-state marginal CDF comparisons, independent of particle identity."""

import math
from pathlib import Path
import random
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import MagicMock, Mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bench"))
import run as runner
sys.path.insert(0, str(ROOT / "figures"))
import plot


def state(values):
    return [[float(i + 1), value, value, value, 0, 0, 0, value, value, value]
            for i, value in enumerate(values)]


class DistributionTests(unittest.TestCase):
    def test_known_distances_and_ties(self):
        for a, b, expected in (([1, 1, 2], [1, 1, 2], 0),
                               ([0, 0], [1, 1], 1),
                               ([0, 0, 1, 2], [0, 1, 1, 3], 0.25),
                               ([0], [0], 0), ([0], [1], 1)):
            with self.subTest(a=a, b=b):
                rows = runner.distribution_errors(state(a), state(b))
                self.assertEqual(len(rows), 6)
                self.assertEqual([row["ks_distance"] for row in rows], [expected] * 6)

    def test_order_and_mass_do_not_weight_particle_distributions(self):
        reference = state([0, 1, 2])
        actual = list(reversed(reference))
        actual = [[row[0] * 100] + row[1:] for row in actual]
        self.assertEqual(runner.distribution_summary(runner.distribution_errors(actual, reference)),
                         {"position_ks_max": 0, "acceleration_ks_max": 0})

    def test_components_are_independent_and_velocity_is_unused(self):
        reference = state([0, 0])
        for field, component, column in runner.DISTRIBUTION_COMPONENTS:
            actual = [row[:] for row in reference]
            actual[0][column] = 1
            actual[0][4:7] = [100, 200, 300]
            rows = runner.distribution_errors(actual, reference)
            self.assertEqual([(row["component"], row["ks_distance"]) for row in rows],
                             [(name, 0.5 if name == component else 0)
                              for _, name, _ in runner.DISTRIBUTION_COMPONENTS])
            summary = runner.distribution_summary(rows)
            self.assertEqual(summary[field + "_ks_max"], 0.5)
            self.assertEqual(sum(summary.values()), 0.5)

    def test_matches_independent_cdf_counts(self):
        rng = random.Random(42)
        for n in (1, 2, 7, 32):
            for _ in range(20):
                a, b = [[rng.randrange(-4, 5) for _ in range(n)] for _ in range(2)]
                expected = max(abs(sum(x <= v for x in a) - sum(x <= v for x in b))
                               for v in set(a + b)) / n
                self.assertEqual(runner.distribution_errors(state(a), state(b))[0]["ks_distance"], expected)
                self.assertEqual(runner.distribution_errors(state(b), state(a))[0]["ks_distance"], expected)

    def test_rejects_empty_mismatched_malformed_and_nonfinite_states(self):
        for a, b in (([], []), (state([1]), state([1, 2])),
                     ([[1, 2]], state([1])), (state([math.nan]), state([1])),
                     (state([math.inf]), state([1])), ([[0] + [0] * 9], state([1]))):
            with self.subTest(a=a, b=b), self.assertRaises(ValueError):
                runner.distribution_errors(a, b)

    def test_plot_loads_only_two_state_files_and_draws_six_overlays(self):
        actual, reference = state([1, 2, 3]), state([0, 1, 2])
        with TemporaryDirectory() as tmp:
            paths = [Path(tmp) / name for name in ("actual.txt", "reference.txt")]
            for path, values in zip(paths, (actual, reference)):
                path.write_text(str(len(values)) + "\n" + "\n".join(
                    " ".join(map(str, row)) for row in values) + "\n")
            data = plot.load_distribution_comparison(*paths, ("Candidate", "Reference"))
        # Drawing still works after the files have been removed.
        axes = [Mock() for _ in range(6)]
        axes[0].get_legend_handles_labels.return_value = ([], [])
        grid = MagicMock()
        grid.flat = axes
        grid.__getitem__.side_effect = lambda key: axes[0] if key == (0, 0) else [axes[0], axes[3]]
        plt, fig = Mock(), Mock()
        plt.subplots.return_value = fig, grid
        self.assertIs(plot.draw_distribution_comparison(plt, data), fig)
        for ax in axes:
            self.assertEqual(ax.step.call_count, 2)
            for call in ax.step.call_args_list:
                x, y = call.args
                self.assertEqual(x, sorted(x))
                self.assertEqual((y[0], y[-1]), (0, 1))
                self.assertEqual(call.kwargs["where"], "post")
            self.assertIn("D = 0.3333", ax.text.call_args.args)

    def test_zoom_changes_limits_without_filtering_cdf_samples(self):
        actual, reference = state(range(100)), state(range(10, 110))
        data = actual, reference, runner.distribution_errors(actual, reference), ("A", "B")
        axes = [Mock() for _ in range(6)]
        axes[0].get_legend_handles_labels.return_value = ([], [])
        grid = MagicMock()
        grid.flat = axes
        grid.__getitem__.side_effect = lambda key: axes[0] if key == (0, 0) else [axes[0], axes[3]]
        plt, fig = Mock(), Mock()
        plt.subplots.return_value = fig, grid
        plot.draw_distribution_comparison(plt, data, central_fraction=0.5)
        for ax in axes:
            ax.set_xlim.assert_called_once_with(25, 84)
            for call in ax.step.call_args_list:
                self.assertEqual(len(call.args[0]), 102)
                self.assertEqual(call.args[1][-1], 1)
            self.assertIn("D = 0.1", ax.text.call_args.args)
        for fraction in (0, -1, 2, math.nan):
            with self.assertRaises(ValueError):
                plot.draw_distribution_comparison(plt, data, central_fraction=fraction)


if __name__ == "__main__":
    unittest.main()
