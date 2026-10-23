import csv
from datetime import datetime, timezone
import importlib.util
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

    def driver(self, *, fake=False, mode="timing", samples=3, warmups=2, instance=None, env=None):
        return subprocess.run(
            [str(ROOT / "build" / ("bench-test" if fake else "bench")),
             str(instance or self.instance), "0.01", "3", str(samples), str(warmups), mode,
             str(self.directory / "final.txt")],
            capture_output=True, text=True, env=env)

    def args(self, mode="timing"):
        return SimpleNamespace(mode=mode, dt=0.01, steps=3, samples=3, warmups=1,
                               events=benchmark.DEFAULT_EVENTS, inputs=[str(self.instance)],
                               results_dir=self.directory / "results")

    def test_statistics(self):
        summary = benchmark.timing_summary([1.0, 2.0, 3.0])
        self.assertEqual(summary["mean_seconds"], 2.0)
        self.assertEqual(summary["variance_seconds_squared"], 1.0)
        self.assertEqual(summary["stddev_seconds"], 1.0)
        for bad in ([1.0], [float("nan"), 1.0], [-1.0, 2.0]):
            with self.assertRaises(ValueError):
                benchmark.timing_summary(bad)

    def test_warmup_reset_and_timing_sample_count(self):
        result = self.driver(fake=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr.splitlines(), ["simulate:ungated"] * 5)
        rows = list(csv.DictReader(result.stdout.splitlines()))
        self.assertEqual(len(rows), 3)  # Two warmups are excluded.
        self.assertTrue(all(float(row["seconds"]) >= 0 for row in rows))
        self.assertEqual((self.directory / "final.txt").read_text(), "1\n1 1 0 0 2 0 0 3 0 0\n")

    def test_perf_gate_surrounds_only_measured_call(self):
        result = self.driver(fake=True, mode="perf", samples=1)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertEqual(result.stderr.splitlines(),
                         ["simulate:ungated", "simulate:ungated", "enable", "simulate:measured", "disable"])

    def test_failed_simulation_still_disables_perf(self):
        result = self.driver(fake=True, mode="perf", samples=1, warmups=0,
                             env={**os.environ, "BENCH_TEST_FAIL": "1"})
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stderr.splitlines()[:3], ["enable", "simulate:measured", "disable"])
        self.assertFalse((self.directory / "final.txt").exists())

    def test_matches_application(self):
        instance = ROOT / "tests" / "data" / "two_body.txt"
        normal = subprocess.run([str(ROOT / "build" / "app"), str(instance), "0.01", "3"],
                                capture_output=True, text=True, check=True)
        measured = self.driver(instance=instance)
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
                                   events="cycles:u,instructions:u")
            for mode in ("profile", "counters"):
                status, _ = benchmark.collect(mode, ROOT / "build" / "bench", instance,
                                              directory / mode, args)
                self.assertEqual(status["status"], "ok", status)


if __name__ == "__main__":
    unittest.main()
