from contextlib import redirect_stdout, redirect_stderr
import csv
from datetime import datetime, timezone
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("benchmark", ROOT / "bench" / "run.py")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


class BenchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / "build")
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.instance = self.directory / "one.txt"
        self.instance.write_text("1\n1 0 0 0 0 0 0 0 0 0\n")

    def driver(self, *, fake=False, mode="timing", samples=3, warmups=2, steps=3,
               instance=None, env=None, solver="direct"):
        return subprocess.run(
            [str(ROOT / "build" / ("bench-test" if fake else "bench")),
             str(instance or self.instance), "0.01", str(steps), str(samples), str(warmups), mode,
             str(self.directory / "final.txt"), solver],
            capture_output=True, text=True, env=env)

    def args(self, mode="timing"):
        return SimpleNamespace(mode=mode, dt=0.01, steps=3, samples=3, warmups=1,
                               events=benchmark.DEFAULT_EVENTS, inputs=[str(self.instance)],
                               results_dir=self.directory / "results", solver="direct", compare_direct=False)

    def test_solver_selection_and_reference_boundaries(self):
        for solver in ("direct", "approximate"):
            for mode in ("timing", "perf"):
                result = self.driver(fake=True, mode=mode, samples=1, solver=solver,
                                     env={**os.environ, "BENCH_TEST_SOLVER": solver})
                self.assertEqual(result.returncode, 0, result.stderr)
        result = self.driver(fake=True, mode="reference", samples=1, warmups=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr.splitlines(), ["simulate:ungated"])
        self.assertEqual(result.stdout, "")
        for options in ({"solver": "unknown"}, {"mode": "reference"},
                        {"mode": "reference", "samples": 1, "warmups": 0, "solver": "approximate"}):
            self.assertNotEqual(self.driver(**options).returncode, 0)

    def test_trajectory_error_metrics(self):
        reference = [[1, 3, 4, 0, 0, 0, 0, 0, 0, 0],
                     [2, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        actual = [[1, 6, 8, 0, 1, 0, 0, 0, 0, 0], reference[1]]
        rows = benchmark.trajectory_errors(actual, reference)
        self.assertAlmostEqual(rows[0]["rms_error"], 5 / 2 ** 0.5)
        self.assertEqual(rows[0]["max_error"], 5)
        self.assertEqual(rows[0]["relative_l2_error"], 1)
        self.assertIsNone(rows[1]["relative_l2_error"])
        self.assertEqual(rows[1]["relative_status"], "undefined")
        self.assertEqual(rows[2]["relative_l2_error"], 0)
        for rows in (benchmark.trajectory_errors([], []),
                     benchmark.trajectory_errors(reference, reference)):
            self.assertTrue(all(row["rms_error"] == row["max_error"] == row["relative_l2_error"] == 0
                                for row in rows))
        for bad in ([], [[2] + reference[0][1:], reference[1]],
                    [[1, float("nan")] + [0] * 8, reference[1]]):
            with self.assertRaises(ValueError):
                benchmark.trajectory_errors(bad, reference)
        with self.assertRaisesRegex(ValueError, "Nonfinite error metric"):
            benchmark.trajectory_errors([[1, 1e308] + [0] * 8], [[1, -1e308] + [0] * 8])
        # Tiny and large representable norms must not underflow/overflow by squaring.
        for value in (1e-200, 1e200):
            rows = benchmark.trajectory_errors([[1, value] + [0] * 8], [[1] + [0] * 9])
            self.assertEqual(rows[0]["rms_error"], value)

    def test_state_reader_rejects_corrupt_outputs(self):
        for content in ("", "-1", "1\n1 0", "0 extra", "1\n1 nan 0 0 0 0 0 0 0 0",
                        "1\n0 0 0 0 0 0 0 0 0 0", "1\n1 1_0 0 0 0 0 0 0 0 0",
                        "1\n1 0x1p9999 0 0 0 0 0 0 0 0"):
            self.instance.write_text(content)
            with self.assertRaises(ValueError):
                benchmark.read_state(self.instance)

    def test_state_reader_matches_c_hexadecimal_input(self):
        self.instance.write_text("1\n0x1p0 0x1p2 0 0 0 0 0 0 0 0\n")
        self.assertEqual(benchmark.read_state(self.instance), [[1, 4] + [0] * 8])
        result = self.driver(mode="reference", samples=1, warmups=0)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(benchmark.read_state(self.directory / "final.txt"),
                         benchmark.read_state(self.instance))

    def test_comparison_failure_does_not_prevent_other_modes(self):
        for mode in ("timing", "profile"):
            (self.directory / mode).mkdir()
        (self.directory / "timing/final.txt").write_text("0\n")
        (self.directory / "profile/final.txt").write_text(self.instance.read_text())
        self.assertFalse(benchmark.compare_direct(ROOT / "build" / "bench", self.instance,
                                                  self.directory, ("timing", "profile"), self.args()))
        self.assertEqual(json.loads((self.directory / "reference/status.json").read_text())["status"], "ok")
        self.assertEqual(json.loads((self.directory / "timing/error-status.json").read_text())["status"], "failed")
        self.assertEqual(json.loads((self.directory / "profile/error-status.json").read_text())["status"], "ok")

    def test_reference_once_and_comparison_artifacts(self):
        args = self.args()
        args.solver = "approximate"
        modes = ("timing", "profile", "counters")
        for mode in modes:
            folder = self.directory / mode
            folder.mkdir()
            (folder / "final.txt").write_text(self.instance.read_text())
        with patch.object(benchmark, "run_logged", wraps=benchmark.run_logged) as run:
            self.assertTrue(benchmark.compare_direct(ROOT / "build" / "bench", self.instance,
                                                     self.directory, modes, args))
            self.assertEqual(run.call_count, 1)
            self.assertEqual(run.call_args.args[0][4:7], ["1", "0", "reference"])
        for mode in modes:
            folder = self.directory / mode
            self.assertEqual(json.loads((folder / "error-status.json").read_text())["status"], "ok")
            with (folder / "error.csv").open() as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 3)
            self.assertTrue(all(float(row["rms_error"]) == 0 for row in rows))

    def test_evaluation_failure_preserves_measurements(self):
        args = self.args()
        args.compare_direct = True
        real_run = benchmark.run_logged

        def fail_reference(command, directory, *a, **kw):
            if "reference" in command:
                return 1
            return real_run(command, directory, *a, **kw)

        with patch.object(benchmark, "build_and_archive", return_value=ROOT / "build" / "bench"), \
             patch.object(benchmark, "run_logged", side_effect=fail_reference):
            self.assertEqual(benchmark.experiment(args), 1)
        folder = next(args.results_dir.iterdir()) / "0001-one"
        self.assertEqual(json.loads((folder / "timing/status.json").read_text())["status"], "ok")
        self.assertEqual(json.loads((folder / "timing/error-status.json").read_text())["status"], "failed")
        self.assertEqual(json.loads((folder / "reference/status.json").read_text())["status"], "failed")
        self.assertTrue((folder / "timing/samples.csv").exists())
        self.assertFalse((folder / "timing/error.csv").exists())

    def test_approximate_experiment_with_reference(self):
        args = self.args()
        args.solver = "approximate"
        args.compare_direct = True
        args.inputs.append(str(ROOT / "tests/data/two_body.txt"))
        with patch.object(benchmark, "build_and_archive", return_value=ROOT / "build" / "bench"):
            self.assertEqual(benchmark.experiment(args), 0)
        experiment = next(args.results_dir.iterdir())
        for folder in sorted(experiment.glob("000*")):
            self.assertEqual((folder / "timing/final.txt").read_bytes(),
                             (folder / "reference/final.txt").read_bytes())
            status = json.loads((folder / "timing/status.json").read_text())
            self.assertEqual(status["solver"], "approximate")
            self.assertEqual(status["theta"], 0.5)
        with (experiment / "timing-summary.csv").open() as source:
            self.assertTrue(all(row["solver"] == "approximate" and float(row["theta"]) == 0.5
                                for row in csv.DictReader(source)))

    def test_theta_driver_and_recorded_settings(self):
        args = self.args()
        args.solver = "approximate"
        args.theta = 0.0
        instance = ROOT / "tests/data/1000.txt"
        status, summary = benchmark.collect("timing", ROOT / "build/bench", instance,
                                            self.directory / "full", args)
        self.assertEqual(status["status"], "ok", status)
        for key, value in benchmark.approximate_settings(instance, 0.0).items():
            self.assertEqual(status[key], value)
        totals = status["diagnostics"]
        self.assertEqual(totals["scope"], "warmups_and_samples")
        self.assertEqual(totals["calls"], (args.samples + args.warmups) * args.steps)
        self.assertEqual(totals["accepted_nodes"], 0)
        self.assertEqual(summary["theta"], 0)
        # Actual forwarding: a looser theta accepts aggregates for this input.
        args.theta = 0.8
        status, _ = benchmark.collect("timing", ROOT / "build/bench", instance,
                                      self.directory / "loose", args)
        self.assertEqual(status["status"], "ok", status)
        self.assertEqual(status["theta"], 0.8)
        self.assertGreater(status["diagnostics"]["accepted_nodes"], 0)
        with patch.object(benchmark, "perf_available", return_value="unavailable"):
            status, _ = benchmark.collect("profile", ROOT / "build/bench", instance,
                                          self.directory / "unavailable", args)
        self.assertEqual(status["status"], "failed")
        self.assertEqual(status["node_capacity"], 4000)
        self.assertEqual(status["theta"], 0.8)
        command = [str(ROOT / "build/bench"), str(self.instance), "0.01", "1", "2", "0",
                   "timing", str(self.directory / "final.txt")]
        # Historical internal-driver syntax continues to work.
        self.assertEqual(subprocess.run(command, capture_output=True).returncode, 0)
        for solver, theta in (("direct", "0"), ("approximate", "nan"),
                              ("approximate", "-1"), ("approximate", "inf")):
            self.assertNotEqual(subprocess.run(command + [solver, theta],
                                               capture_output=True).returncode, 0)

    def test_theta_cli_validation(self):
        command = ["timing", str(self.instance), "--dt", "0.01", "--steps", "2"]
        self.assertIsNone(benchmark.parse_args(command).theta)
        self.assertEqual(benchmark.parse_args(command + ["--solver", "approximate"]).theta, 0.5)
        self.assertEqual(benchmark.parse_args(command + ["--solver", "approximate", "--theta", "0"]).theta, 0)
        for options in (["--theta", "0"], ["--solver", "direct", "--theta", "0.5"],
                        *[["--solver", "approximate", "--theta", value]
                          for value in ("-1", "nan", "inf", "bad")]):
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                benchmark.parse_args(command + options)

    def test_statistics(self):
        summary = benchmark.timing_summary([2.0, 4.0, 6.0], 2)
        self.assertEqual(summary["samples"], 3)
        self.assertEqual(summary["mean_seconds_per_step"], 2.0)
        self.assertEqual(summary["variance_seconds_per_step_squared"], 1.0)
        self.assertEqual(summary["stddev_seconds_per_step"], 1.0)
        for bad in ([1.0], [float("nan"), 1.0], [-1.0, 2.0]):
            with self.assertRaises(ValueError):
                benchmark.timing_summary(bad, 2)
        for steps in (0, -1, 1.5, 2.0, "2", True, None):
            with self.subTest(steps=steps), self.assertRaisesRegex(ValueError, "positive integer"):
                benchmark.timing_summary([2.0, 4.0, 6.0], steps)

    def test_per_step_summary_and_display_preserve_raw_durations(self):
        args = self.args()
        args.steps = 2
        args.dt = 0.125  # Normalization uses the number of steps, not dt.
        raw = "sample,seconds\n1,0.002\n2,0.004\n3,0.006\n"

        def fake_run(command, folder, stdout="stdout.log", stderr="stderr.log"):
            self.assertEqual(command[2:7], ["0.125", "2", "3", "1", "timing"])
            self.assertEqual(stdout, "samples.csv")
            (folder / stdout).write_text(raw)
            return 0

        messages = io.StringIO()
        with patch.object(benchmark, "build_and_archive", return_value=Path("driver")), \
             patch.object(benchmark, "run_logged", side_effect=fake_run), redirect_stdout(messages):
            self.assertEqual(benchmark.experiment(args), 0)
        experiment = next(args.results_dir.iterdir())
        timing = experiment / "0001-one" / "timing"
        self.assertEqual((timing / "samples.csv").read_text(), raw)
        status = json.loads((timing / "status.json").read_text())
        self.assertEqual((status["steps"], status["dt"]), (2, 0.125))
        for path in (timing / "summary.csv", experiment / "timing-summary.csv"):
            with path.open() as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual((row["solver"], row["steps"], row["dt"], row["samples"]),
                             ("direct", "2", "0.125", "3"))
            self.assertEqual(float(row["mean_seconds_per_step"]), 0.002)
            self.assertAlmostEqual(float(row["variance_seconds_per_step_squared"]), 1e-6)
            self.assertEqual(float(row["stddev_seconds_per_step"]), 0.001)
            self.assertNotIn("mean_seconds", row)
        self.assertIn("Timing 0001-one [direct]: 2 ± 1 ms/step "
                      "(mean ± sample SD; 2 steps, 3 samples)", messages.getvalue())

    def test_warmup_reset_and_timing_sample_count(self):
        for steps in (1, 2, 3, 4):
            for warmups in (0, 1, 2):
                with self.subTest(steps=steps, warmups=warmups):
                    result = self.driver(fake=True, steps=steps, warmups=warmups)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(result.stderr.splitlines(), ["simulate:ungated"] * (warmups + 3))
                    rows = list(csv.DictReader(result.stdout.splitlines()))
                    self.assertEqual(len(rows), 3)  # Warmups are excluded.
                    self.assertTrue(all(float(row["seconds"]) >= 0 for row in rows))
                    self.assertEqual((self.directory / "final.txt").read_text(),
                                     "1\n1 2 3 4 5 6 7 8 9 10\n")

    def test_perf_gate_surrounds_only_measured_call(self):
        result = self.driver(fake=True, mode="perf", samples=1)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertEqual(result.stderr.splitlines(),
                         ["simulate:ungated", "simulate:ungated", "enable", "simulate:measured", "disable"])

    def test_nonfinite_simulation_still_disables_perf(self):
        for steps in (1, 2):
            for axis in ("1", "2", "3"):
                with self.subTest(steps=steps, axis=axis):
                    result = self.driver(fake=True, mode="perf", samples=1, warmups=0, steps=steps,
                                         env={**os.environ, "BENCH_TEST_FAIL": axis})
                    self.assertNotEqual(result.returncode, 0)
                    self.assertEqual(result.stderr.splitlines()[:3],
                                     ["enable", "simulate:measured", "disable"])
                    self.assertIn("Sample 1 produced a nonfinite state.", result.stderr)
                    self.assertFalse((self.directory / "final.txt").exists())

    def test_nonfinite_warmup_is_rejected(self):
        for steps in (1, 2):
            for axis in ("1", "2", "3"):
                with self.subTest(steps=steps, axis=axis):
                    result = self.driver(fake=True, steps=steps,
                                         env={**os.environ, "BENCH_TEST_FAIL": axis})
                    self.assertNotEqual(result.returncode, 0)
                    self.assertEqual(result.stdout, "")
                    self.assertIn("Warmup produced a nonfinite state.", result.stderr)
                    self.assertFalse((self.directory / "final.txt").exists())

    def test_matches_application(self):
        empty = self.directory / "empty.txt"
        empty.write_text("0\n")
        for instance in (empty, self.instance, ROOT / "tests" / "data" / "two_body.txt"):
            for steps in (1, 2, 3, 4):
                with self.subTest(instance=instance.name, steps=steps):
                    normal = subprocess.run(
                        [str(ROOT / "build" / "app"), str(instance), "0.01", str(steps)],
                        capture_output=True, text=True, check=True)
                    measured = self.driver(instance=instance, steps=steps)
                    self.assertEqual(measured.returncode, 0, measured.stderr)
                    self.assertEqual((self.directory / "final.txt").read_text(), normal.stdout)

    def test_nonfinite_state_fails_without_summary(self):
        self.instance.write_text("1\n1 1e308 0 0 1e308 0 0 0 0 0\n")
        args = self.args()
        args.dt = 10.0
        args.warmups = 0
        folder = self.directory / "timing"
        status, summary = benchmark.collect("timing", ROOT / "build" / "bench", self.instance, folder, args)
        self.assertEqual(status["status"], "failed")
        self.assertIsNone(summary)
        self.assertIn("nonfinite", (folder / "stderr.log").read_text())
        self.assertFalse((folder / "summary.csv").exists())

    def test_generator_compatible_and_reproducible(self):
        command = [sys.executable, str(ROOT / "generate_instance.py"), "--n", "12", "--L", "10",
                   "--w-min", "1", "--w-max", "2", "--v-max", "0", "--a-max", "0", "--seed", "42"]
        outputs = [self.directory / "generated1.txt", self.directory / "generated2.txt"]
        for output in outputs:
            subprocess.run([*command, "--output", str(output)], check=True, capture_output=True)
        self.assertEqual(outputs[0].read_bytes(), outputs[1].read_bytes())
        self.assertEqual(outputs[0].read_text().splitlines()[0], "12")
        subprocess.run([str(ROOT / "build" / "app"), str(outputs[0]), "0.01", "0"],
                       check=True, capture_output=True)

    def test_counter_parsing_and_unavailable_values(self):
        rows = benchmark.parse_counters(
            "# started\n100;;cycles:u;10000;50.00;;\n250;;instructions:u;10000;50.00;;\n"
            "20;;branches:u;10000;100.00;;\n2;;branch-misses:u;10000;100.00;;\n"
            "<not supported>;;cache-references:u;0;100.00;;\n"
            "<not counted>;;cache-misses:u;0;0.00;;\n1.0;seconds time elapsed\n")
        self.assertEqual(len(rows), 6)
        self.assertEqual(rows[0]["running_percent"], "50.00")
        self.assertIsNone(rows[4]["value"])
        self.assertIsNone(rows[5]["value"])
        ratios = benchmark.counter_ratios(rows)
        self.assertEqual(ratios[0]["value"], 2.5)
        self.assertEqual(ratios[1]["value"], 0.1)
        self.assertIsNone(ratios[2]["value"])
        self.assertEqual(benchmark.counter_ratios(rows + [rows[0]])[0]["status"], "unavailable")

    def test_missing_tools_and_platforms_saved(self):
        for platform_name in ("Darwin", "Linux"):
            with patch.object(benchmark.platform, "system", return_value=platform_name), \
                 patch.object(benchmark.shutil, "which", return_value=None):
                folder = self.directory / platform_name
                status, summary = benchmark.collect("counters", Path("unused"), self.instance, folder, self.args())
                self.assertEqual(status["status"], "failed")
                self.assertIsNone(summary)
                self.assertEqual(json.loads((folder / "status.json").read_text()), status)

    def test_perf_command_boundaries_and_zero_sample_profile(self):
        commands = []

        def fake_run(command, folder, stdout="stdout.log", stderr="stderr.log"):
            commands.append([str(item) for item in command])
            (folder / stdout).write_text("")
            return 0

        with patch.object(benchmark, "perf_available", return_value=None), \
             patch.object(benchmark, "run_logged", side_effect=fake_run):
            status, _ = benchmark.collect("profile", Path("driver"), self.instance,
                                          self.directory / "profile", self.args())
        self.assertEqual(status["status"], "failed")
        self.assertIn("no samples", status["error"])
        self.assertIn("--delay=-1", commands[0])
        self.assertIn("cpu-clock:u", commands[0])
        self.assertIn("perf", commands[0][commands[0].index("--") + 1:])

    def test_partial_counters_and_permission_failure_saved(self):
        def fake_run(command, folder, stdout="stdout.log", stderr="stderr.log"):
            self.assertIn("--delay=-1", command)
            (folder / "perf-stat.txt").write_text(
                "10;;cycles:u;50;100.00;;\n<not supported>;;cache-misses:u;0;0.00;;\n")
            (folder / stderr).write_text("permission denied\n")
            return 1

        with patch.object(benchmark, "perf_available", return_value=None), \
             patch.object(benchmark, "run_logged", side_effect=fake_run):
            folder = self.directory / "counters"
            status, _ = benchmark.collect("counters", Path("driver"), self.instance, folder, self.args())
        self.assertEqual(status["status"], "failed")
        with (folder / "counters.csv").open() as source:
            rows = list(csv.DictReader(source))
        self.assertEqual(rows[0]["value"], "10.0")
        self.assertEqual(rows[1]["value"], "")
        self.assertIn("permission denied", (folder / "stderr.log").read_text())

    def test_multi_instance_experiment_and_overwrite_protection(self):
        args = self.args("all")
        args.inputs.append(str(ROOT / "tests" / "data" / "two_body.txt"))
        with patch.object(benchmark, "datetime") as clock:
            clock.now.return_value = datetime(2026, 10, 22, 5, 34, 19, 808566, tzinfo=timezone.utc)
            with patch.object(benchmark, "perf_available", return_value="Linux perf unavailable in test"):
                self.assertEqual(benchmark.experiment(args), 1)
            experiment = next(args.results_dir.iterdir())
            self.assertRegex(experiment.name, r"^[0-9a-f]{7,64}_20261022-053419\.808566Z_all$")
            self.assertEqual({path.name for path in experiment.iterdir()},
                             {"bench", "build.stdout.log", "timing-summary.csv", "0001-one", "0002-two_body"})
            for name in ("0001-one", "0002-two_body"):
                instance = experiment / name
                for mode, expected in (("timing", "ok"), ("profile", "failed"), ("counters", "failed")):
                    status = json.loads((instance / mode / "status.json").read_text())
                    self.assertEqual(status["status"], expected)
                    self.assertEqual((status["steps"], status["dt"]), (args.steps, args.dt))
                self.assertTrue((instance / "input.txt").is_file())
                self.assertTrue((instance / "timing" / "final.txt").is_file())
            self.assertEqual((experiment / "bench").read_bytes(), (ROOT / "build" / "bench").read_bytes())
            self.assertTrue(os.access(experiment / "bench", os.X_OK))
            self.assertIn("-O3", (experiment / "build.stdout.log").read_text())
            with (experiment / "timing-summary.csv").open() as source:
                self.assertEqual(len(list(csv.DictReader(source))), 2)
            forbidden = {"source", "source-hashes.json", "manifest.json", "commands.jsonl", "build.stderr.log"}
            self.assertFalse(any(path.name in forbidden for path in experiment.rglob("*")))

            before = {path.relative_to(experiment): path.read_bytes()
                      for path in experiment.rglob("*") if path.is_file()}
            with patch.object(benchmark, "build_and_archive") as build:
                with self.assertRaises(FileExistsError):
                    benchmark.experiment(args)
                build.assert_not_called()
            self.assertEqual(before, {path.relative_to(experiment): path.read_bytes()
                                      for path in experiment.rglob("*") if path.is_file()})

    def test_directory_naming_for_each_mode(self):
        revision = subprocess.CompletedProcess([], 0, stdout="abc1234\n", stderr="")
        with patch.object(benchmark.subprocess, "run", return_value=revision) as git, \
             patch.object(benchmark, "datetime") as clock:
            clock.now.return_value = datetime(2026, 10, 22, 5, 34, 19, 808566, tzinfo=timezone.utc)
            for mode in ("timing", "profile", "counters", "all"):
                directory = benchmark.create_experiment_directory(self.directory / "custom", mode)
                self.assertEqual(directory.name, f"abc1234_20261022-053419.808566Z_{mode}")
                self.assertEqual(directory.parent, self.directory / "custom")
                self.assertEqual(list(directory.iterdir()), [])
            clock.now.assert_called_with(timezone.utc)
            self.assertEqual(git.call_args.args[0],
                             ["git", "rev-parse", "--short=7", "--verify", "HEAD^{commit}"])

    def test_missing_git_commit_does_not_create_misleading_directory(self):
        revision = subprocess.CompletedProcess([], 128, stdout="", stderr="no commit\n")
        parent = self.directory / "results"
        with patch.object(benchmark.subprocess, "run", return_value=revision):
            with self.assertRaisesRegex(RuntimeError, "no commit"):
                benchmark.create_experiment_directory(parent, "timing")
        self.assertFalse(parent.exists())

    def test_combined_build_log_preserves_errors(self):
        command = [sys.executable, "-c",
                   "import sys; print('build command', flush=True); print('compiler error', file=sys.stderr); sys.exit(2)"]
        code = benchmark.run_logged(command, self.directory, "build.stdout.log", stderr=None)
        self.assertEqual(code, 2)
        log = (self.directory / "build.stdout.log").read_text()
        self.assertIn("build command", log)
        self.assertIn("compiler error", log)
        self.assertFalse((self.directory / "build.stderr.log").exists())
        self.assertFalse((self.directory / "commands.jsonl").exists())

    def test_build_failure_saved_in_error_file(self):
        args = self.args()

        def failed_build(command, folder, stdout="stdout.log", stderr="stderr.log"):
            self.assertEqual(stdout, "build.stdout.log")
            self.assertIsNone(stderr)
            (folder / stdout).write_text("compiler command\ncompiler error\n")
            return 2

        with patch.object(benchmark, "run_logged", side_effect=failed_build):
            self.assertEqual(benchmark.experiment(args), 1)
        experiment = next(args.results_dir.iterdir())
        self.assertEqual({path.name for path in experiment.iterdir()}, {"error.txt", "build.stdout.log"})
        self.assertIn("Build failed; see build.stdout.log", (experiment / "error.txt").read_text())
        self.assertIn("compiler error", (experiment / "build.stdout.log").read_text())

    def test_setup_failure_saved_in_error_file(self):
        args = self.args()
        with patch.object(benchmark, "build_and_archive", side_effect=OSError("cannot save binary")):
            self.assertEqual(benchmark.experiment(args), 1)
        experiment = next(args.results_dir.iterdir())
        self.assertEqual((experiment / "error.txt").read_text(), "cannot save binary\n")

    def test_results_directory_cli(self):
        required = ["timing", str(self.instance), "--dt", "0.01", "--steps", "10"]
        self.assertEqual(benchmark.parse_args(required).results_dir, ROOT / "results")
        args = benchmark.parse_args([*required, "--results-dir", str(self.directory)])
        self.assertEqual(args.results_dir, self.directory)
        self.assertEqual(args.solver, "direct")
        self.assertFalse(args.compare_direct)
        args = benchmark.parse_args([*required, "--solver", "approximate", "--compare-direct"])
        self.assertEqual(args.solver, "approximate")
        self.assertTrue(args.compare_direct)


@unittest.skipUnless(os.environ.get("NBODY_TEST_PERF") == "1", "set NBODY_TEST_PERF=1 on a perf-enabled Linux host")
class LinuxPerfTests(unittest.TestCase):
    def test_real_profile_and_counters(self):
        self.assertIsNone(benchmark.perf_available())
        with tempfile.TemporaryDirectory(dir=ROOT / "build") as name:
            directory = Path(name)
            instance = directory / "instance.txt"
            subprocess.run([sys.executable, str(ROOT / "generate_instance.py"), "--n", "1000",
                            "--L", "100", "--w-min", "1", "--w-max", "2", "--v-max", "0",
                            "--a-max", "0", "--seed", "42", "--output", str(instance)],
                           check=True, capture_output=True)
            args = SimpleNamespace(dt=1e-5, steps=100, samples=2, warmups=1,
                                   events="cycles:u,instructions:u", solver="direct")
            for mode in ("profile", "counters"):
                status, _ = benchmark.collect(mode, ROOT / "build" / "bench", instance,
                                              directory / mode, args)
                self.assertEqual(status["status"], "ok", status)


if __name__ == "__main__":
    unittest.main()
