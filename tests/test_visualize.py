"""Standalone visualization checks; run with the optional NumPy environment."""

from contextlib import redirect_stdout, redirect_stderr
import io
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import visualize as visual


class VisualizationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.input = self.directory / "initial.txt"
        self.input.write_text((ROOT / "tests/data/two_body.txt").read_text())
        self.app = self.directory / "app"
        self.app.write_text("frozen executable")
        self.app.chmod(0o755)

    def options(self, *extra):
        return visual.parser().parse_args([
            "capture", str(self.input), "--app", str(self.app), "--dt", "0.01",
            "--steps", "5", "--every", "2", "--output-dir", str(self.directory / "output"),
            *extra])

    def quiet_capture(self, args):
        with redirect_stdout(io.StringIO()):
            return visual.capture(args)

    def test_frame_schedule_includes_initial_and_exact_endpoint(self):
        for steps, every, expected in ((5, 2, [0, 2, 4, 5]), (4, 2, [0, 2, 4]),
                                       (0, 2, [0]), (2, 20, [0, 2])):
            np.testing.assert_array_equal(visual.frame_steps(steps, every), expected)
        for steps, every in ((-1, 2), (2, 0)):
            with self.assertRaises(ValueError):
                visual.frame_steps(steps, every)

    def test_parser_checks_complete_finite_positive_mass_state(self):
        valid = "1\n1 2 3 4 5 6 7 8 9 10\n"
        self.assertEqual(visual.parse_state(valid).shape, (1, 10))
        for text in ("", "0\n", "1.0\n", valid + "0", valid.replace("10", "nan"),
                     valid.replace("10", "inf"), "1\n-1 2 3 4 5 6 7 8 9 10", "1\n1 2"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                visual.parse_state(text)

    def test_option_validation(self):
        for extra in (("--dt", "0"), ("--theta", "0.5"), ("--order", "4"),
                      ("--solver", "barnes-hut", "--theta", "-1")):
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                visual.validate_options(self.options(*extra))
        for extra in (("--dt", "nan"), ("--steps", "-1"), ("--every", "0"),
                      ("--solver", "fmm", "--order", "7")):
            with self.subTest(extra=extra), redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                self.options(*extra)
        args = visual.parser().parse_args(["render", "saved.npz", "--width", "101"])
        with self.assertRaises(ValueError):
            visual.validate_options(args)

    def test_frozen_replay_and_archive_round_trip(self):
        original = self.input.read_text()
        calls = []
        parent_environment = os.environ.copy()

        def run(command, **kwargs):
            self.assertEqual(Path(command[0]).read_text(), "frozen executable")
            self.assertEqual(Path(command[-3]).read_text(), original)
            self.assertEqual(kwargs["env"]["OMP_NUM_THREADS"], "1")
            self.assertEqual(kwargs["env"]["OMP_DYNAMIC"], "FALSE")
            self.assertEqual(command[1:3], ["--solver", "n2"])
            calls.append(int(command[-1]))
            # Changes to external files must not affect subsequent frames.
            self.input.write_text("changed")
            self.app.write_text("changed")
            state = visual.parse_state(original)
            state[:, 1] += int(command[-1])
            text = "2\n" + "\n".join(" ".join(map(str, row)) for row in state) + "\n"
            return subprocess.CompletedProcess(command, 0, text, "diagnostic\n")

        with patch.object(visual.subprocess, "run", side_effect=run):
            path = self.quiet_capture(self.options())
        self.assertEqual(calls, [2, 4, 5])
        self.assertEqual(os.environ, parent_environment)
        data = visual.load_trajectory(path)
        np.testing.assert_array_equal(data["steps"], [0, 2, 4, 5])
        np.testing.assert_array_equal(data["times"], [0, .02, .04, .05])
        np.testing.assert_array_equal(data["positions"][:, 0, 0], [-1, 1, 3, 4])
        self.assertEqual(data["positions"].dtype, np.float64)
        self.assertEqual(data["metadata"]["capture"], "faithful-replay")
        self.assertEqual((path.parent / "input.txt").read_text(), original)
        np.testing.assert_array_equal(visual.parse_state((path.parent / "final.txt").read_text())[:, 1:4],
                                      data["positions"][-1])
        with self.assertRaises(FileExistsError):
            visual.save_trajectory(path, data)

    def test_capture_failures_do_not_publish_trajectory(self):
        original = self.input.read_text()
        cases = [(1, "", "solver failure"), (0, original.replace("0.25", "nan"), "nonfinite"),
                 (0, original.replace("1 -1", "2 -1"), "mass changed"),
                 (0, "1\n1 0 0 0 0 0 0 0 0 0\n", "count changed")]
        for index, (code, stdout, stderr) in enumerate(cases):
            with self.subTest(index=index):
                args = self.options("--output-dir", str(self.directory / str(index)))
                result = subprocess.CompletedProcess([], code, stdout, stderr)
                with patch.object(visual.subprocess, "run", return_value=result):
                    with self.assertRaisesRegex(ValueError, "frame 1, step 2"):
                        self.quiet_capture(args)
                self.assertFalse((args.output_dir / "trajectory.npz").exists())
                self.assertIn(stderr, (args.output_dir / "capture.stderr.log").read_text())

    def test_zero_steps_and_output_collision(self):
        args = self.options("--steps", "0")
        with patch.object(visual.subprocess, "run", side_effect=AssertionError("unexpected execution")):
            path = self.quiet_capture(args)
        self.assertEqual(visual.load_trajectory(path)["positions"].shape, (1, 2, 3))
        with self.assertRaises(FileExistsError):
            self.quiet_capture(args)

    def test_launch_error_identifies_frame(self):
        with patch.object(visual.subprocess, "run", side_effect=OSError("launch failed")):
            with self.assertRaisesRegex(ValueError, "frame 1, step 2: launch failed"):
                self.quiet_capture(self.options())

    def test_solver_parameters_are_forwarded(self):
        for solver, flag, value in (("barnes-hut", "--theta", "0.3"), ("fmm", "--order", "6")):
            args = self.options("--solver", solver, flag, value, "--steps", "1",
                                "--output-dir", str(self.directory / solver))
            result = subprocess.CompletedProcess([], 0, self.input.read_text(), "")
            with patch.object(visual.subprocess, "run", return_value=result) as run:
                path = self.quiet_capture(args)
            self.assertEqual(run.call_args.args[0][1:5], ["--solver", solver, flag, value])
            self.assertEqual(str(visual.load_trajectory(path)["metadata"][flag[2:]]), value)

    def test_archive_validation_and_no_pickle(self):
        path = self.quiet_capture(self.options("--steps", "0"))
        data = visual.load_trajectory(path)
        data["positions"][0, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            visual.save_trajectory(self.directory / "invalid.npz", data)
        bad = self.directory / "object.npz"
        np.savez(bad, positions=np.array([{}], dtype=object))
        with self.assertRaises(ValueError):
            visual.load_trajectory(bad)

    def test_bounds_for_planar_and_stationary_states(self):
        for positions in (np.zeros((3, 1, 3)), np.array([[[-1., 0, 0], [1., 0, 0]]])):
            center, half = visual.scene_bounds(positions)
            self.assertTrue(np.isfinite(center).all())
            self.assertGreater(half, 0)
            self.assertTrue((np.abs(positions - center) <= half).all())

    def test_render_command_does_not_capture_or_require_app(self):
        with patch.object(visual, "capture", side_effect=AssertionError("unexpected capture")):
            with patch.object(visual, "render") as render:
                self.assertEqual(visual.main(["render", str(self.directory / "saved.npz")]), 0)
                self.assertFalse(hasattr(render.call_args.args[2], "app"))


@unittest.skipUnless((ROOT / "build/app").is_file(), "run make for real solver checks")
class SolverCaptureTests(unittest.TestCase):
    def test_replay_matches_each_solver_on_nontrivial_input(self):
        with TemporaryDirectory() as temporary:
            for solver in ("n2", "barnes-hut", "fmm"):
                with self.subTest(solver=solver):
                    input_path = ROOT / "bench/data/1000.txt"
                    args = visual.parser().parse_args([
                        "capture", str(input_path), "--solver", solver, "--dt", "0.0001",
                        "--steps", "3", "--every", "2", "--output-dir", str(Path(temporary) / solver)])
                    with redirect_stdout(io.StringIO()):
                        trajectory = visual.load_trajectory(visual.capture(args))
                    for frame, step in enumerate(trajectory["steps"]):
                        completed = subprocess.run(
                            [str(ROOT / "build/app"), "--solver", solver, str(input_path),
                             "0.0001", str(step)], text=True, capture_output=True, check=True,
                            env=dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE"))
                        reference = visual.parse_state(completed.stdout)
                        np.testing.assert_allclose(trajectory["positions"][frame], reference[:, 1:4],
                                                   rtol=1e-13, atol=1e-13)

    def test_two_body_orbit_stays_circular_and_planar(self):
        with TemporaryDirectory() as temporary:
            args = visual.parser().parse_args([
                "capture", str(ROOT / "tests/data/two_body.txt"), "--dt", "0.01",
                "--steps", "1260", "--every", "126", "--output-dir", str(Path(temporary) / "orbit")])
            with redirect_stdout(io.StringIO()):
                data = visual.load_trajectory(visual.capture(args))
            np.testing.assert_allclose(np.linalg.norm(data["positions"], axis=2), 1, atol=2e-5, rtol=0)
            np.testing.assert_array_equal(data["positions"][:, :, 2], 0)
            np.testing.assert_allclose(data["positions"].sum(axis=1), 0, atol=1e-13)


@unittest.skipUnless(os.environ.get("NBODY_RENDER_TESTS") == "1", "opt-in graphics/encoding checks")
class RenderTests(unittest.TestCase):
    def test_render_without_executable_and_failed_encoding_publication(self):
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            trajectory = directory / "saved.npz"
            visual.save_trajectory(trajectory, {
                "positions": np.array([[[-1., 0, 0], [1., 0, 0]],
                                       [[-.9, -.1, 0], [.9, .1, 0]],
                                       [[-.8, -.2, 0], [.8, .2, 0]]]),
                "masses": np.ones(2), "steps": np.arange(3, dtype=np.int64),
                "times": np.arange(3, dtype=np.float64) * .01,
                "metadata": {"version": 1, "solver": "n2", "dt": .01,
                             "app": str(directory / "nonexistent-app")},
            })
            args = visual.parser().parse_args([
                "render", str(trajectory), "--width", "640", "--height", "360", "--camera", "orbit"])
            output = directory / "video.mp4"
            with redirect_stdout(io.StringIO()):
                visual.render(trajectory, output, args)
            visual.verify_video(output, 3, 640, 360)
            with self.assertRaisesRegex(ValueError, "dimensions"):
                visual.verify_video(output, 3, 640, 362)
            with self.assertRaisesRegex(ValueError, "expected 4"):
                visual.verify_video(output, 4, 640, 360)
            with self.assertRaises(FileExistsError):
                visual.render(trajectory, output, args)
            failed_output = directory / "failed.mp4"
            with patch.object(visual, "verify_video", side_effect=ValueError("encoding failed")):
                with redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, "encoding failed"):
                    visual.render(trajectory, failed_output, args)
            self.assertFalse(failed_output.exists())
            self.assertTrue(trajectory.exists())


if __name__ == "__main__":
    unittest.main()
