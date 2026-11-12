"""Repeated-run averages, reset trajectories, and saved-data validation."""

import csv
import json
import statistics
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bench"))
sys.path.insert(0, str(ROOT / "figures"))
import compare_solvers as sweep
import run as runner
from plot import draw_parallel_sweep, load_solver_sweep, load_parallel_sweep, prepare_figures


class SolverSweepTests(unittest.TestCase):
    def test_default_workload(self):
        args = sweep.parse_args([])
        self.assertEqual(args.bodies, [1000, 5000, 10000, 20000, 50000])
        self.assertEqual(args.steps, 20)
        self.assertEqual(args.samples, 3)
        self.assertEqual((args.theta, args.order), (0.5, 4))

    def test_parallel_defaults_and_overrides(self):
        args = sweep.parse_args(["--parallel"])
        self.assertEqual(args.bodies, [1000, 5000, 10000, 20000, 50000, 100000])
        self.assertEqual(args.steps, 20)
        self.assertEqual(args.results_dir, ROOT / "results/parallel-sweep")
        self.assertEqual(sweep.PARALLEL_SERIES,
                         (("n2", 1), ("n2", 4), ("n2", 8),
                          ("barnes-hut", 1), ("barnes-hut", 4), ("barnes-hut", 8),
                          ("fmm", 1), ("fmm", 4), ("fmm", 8)))
        self.assertEqual(sweep.parse_args(["--parallel", "--steps", "5"]).steps, 5)

    def test_large_parallel_configurations(self):
        self.assertEqual(sweep.configurations(50000, True), list(sweep.PARALLEL_SERIES))
        for count in (50001, 100000, 200000):
            self.assertEqual(sweep.configurations(count, True),
                             [("n2", 8), ("barnes-hut", 1), ("barnes-hut", 4),
                              ("barnes-hut", 8), ("fmm", 1), ("fmm", 4), ("fmm", 8)])
        self.assertEqual(sweep.configurations(100000, False),
                         [("n2", 1), ("barnes-hut", 1), ("fmm", 1)])

    def test_large_parallel_run_uses_eight_thread_reference(self):
        with TemporaryDirectory() as temporary:
            directory = Path(temporary) / "run"
            directory.mkdir()
            args = sweep.parse_args(["--parallel", "--bodies", "100000"])
            final = runner.read_state(ROOT / "tests/data/two_body.txt")
            summary = {"mean_ms_per_step": 1.0, "stddev_ms_per_step": 0.1}
            with patch.object(runner, "create_experiment_directory", return_value=directory), \
                    patch.object(sweep, "build", return_value=Path("unused")), \
                    patch.object(runner, "run_logged", return_value=0), \
                    patch.object(sweep, "collect") as collect, patch("builtins.print"):
                def collect_case(binary, instance, case, solver, args, threads):
                    case.mkdir()
                    return summary, final
                collect.side_effect = collect_case
                self.assertEqual(sweep.experiment(args), 0)
            self.assertEqual([(call.args[3], call.args[5]) for call in collect.call_args_list],
                             sweep.configurations(100000, True))
            with (directory / "parallel-sweep.csv").open() as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 7)
            self.assertTrue(all(row["reference_threads"] == "8" for row in rows))
            self.assertEqual(json.loads((directory / "status.json").read_text())["status"], "ok")

    def test_sample_count_override_and_validation(self):
        self.assertEqual(sweep.parse_args(["--samples", "5"]).samples, 5)
        with patch("sys.stderr"), self.assertRaises(SystemExit):
            sweep.parse_args(["--samples", "1"])

    def test_thread_check_and_environment(self):
        record = {"openmp": True, "requested_threads": 4, "effective_threads": 4}
        self.assertEqual(sweep.check_threads("thread_check: " + json.dumps(record), 4, True), record)
        for log in ("", "thread_check: " + json.dumps({**record, "effective_threads": 1})):
            with self.assertRaisesRegex(ValueError, "thread check"):
                sweep.check_threads(log, 4, True)
        with self.assertRaisesRegex(ValueError, "serial"):
            sweep.check_threads("thread_check: " + json.dumps(record), 1, False)
        before = dict(os.environ)
        with sweep.thread_environment(8):
            self.assertEqual(os.environ["OMP_NUM_THREADS"], "8")
            self.assertEqual(os.environ["OMP_THREAD_LIMIT"], "8")
            self.assertEqual(os.environ["OMP_DYNAMIC"], "FALSE")
        self.assertEqual(dict(os.environ), before)

    def test_mean_and_sd_use_whole_run_averages(self):
        result = sweep.summarize([1, 2, 3], 20)
        self.assertEqual(result["mean_ms_per_step"], 100)
        self.assertAlmostEqual(result["stddev_ms_per_step"], 50)
        for durations, steps in (([1], 20), ([1, float("nan")], 20), ([-1, 2], 20), ([1, 2], 0)):
            with self.subTest(durations=durations, steps=steps), self.assertRaises(ValueError):
                sweep.summarize(durations, steps)

    def test_parallel_figure_has_nine_lines_with_between_run_error_bars(self):
        plt, fig, ax = Mock(), Mock(), Mock()
        plt.subplots.return_value = fig, ax
        rows = [{"solver": solver, "threads": threads, "bodies": 1000,
                 "mean_ms_per_step": 10 / threads, "stddev_ms_per_step": 0.1 / threads}
                for solver, threads in sweep.PARALLEL_SERIES]
        draw_parallel_sweep(plt, ({"settings": {"bodies": [1000], "steps": 20, "samples": 3},
                                   "error_bars": "whole_run_sample_sd"}, rows))
        self.assertEqual(ax.errorbar.call_count, 9)
        for call, row in zip(ax.errorbar.call_args_list, rows):
            self.assertEqual(call.args, ([1000], [row["mean_ms_per_step"]]))
            self.assertEqual(call.kwargs["yerr"], [row["stddev_ms_per_step"]])
        self.assertEqual({call.kwargs["linestyle"] for call in ax.errorbar.call_args_list},
                         {"-", "--", ":"})
        calls = ax.errorbar.call_args_list
        for group in (calls[:3], calls[3:6], calls[6:9]):
            self.assertEqual(len({call.kwargs["color"] for call in group}), 1)
            self.assertEqual(len({call.kwargs["marker"] for call in group}), 1)
            self.assertEqual([call.kwargs["linestyle"] for call in group], ["-", "--", ":"])
        self.assertNotEqual(calls[0].kwargs["marker"], calls[3].kwargs["marker"])
        self.assertNotEqual(calls[0].kwargs["color"], calls[3].kwargs["color"])

    def test_parallel_figure_preserves_series_endpoints(self):
        plt, fig, ax = Mock(), Mock(), Mock()
        plt.subplots.return_value = fig, ax
        rows = [{"solver": solver, "threads": threads, "bodies": count,
                 "mean_ms_per_step": 1, "stddev_ms_per_step": 0.1}
                for count in (50000, 100000)
                for solver, threads in sweep.configurations(count, True)]
        draw_parallel_sweep(plt, ({"settings": {"bodies": [50000, 100000], "steps": 20,
                                               "samples": 3}}, rows))
        self.assertEqual(ax.errorbar.call_count, 9)
        for index, call in enumerate(ax.errorbar.call_args_list):
            self.assertEqual(call.args[0], [50000] if index < 2 else [50000, 100000])

    def test_instrumented_trajectory_and_plot_validation(self):
        self.check_trajectory_and_plot(False)

    @unittest.skipUnless((ROOT / "build/bench_steps_openmp").exists(), "requires OpenMP build")
    def test_parallel_trajectory_and_plot_validation(self):
        self.check_trajectory_and_plot(True)

    @unittest.skipUnless((ROOT / "build/bench_steps_openmp").exists(), "requires OpenMP build")
    def test_driver_rejects_reduced_team(self):
        with TemporaryDirectory() as temporary:
            result = subprocess.run([str(ROOT / "build/bench_steps_openmp"),
                                     str(ROOT / "tests/data/two_body.txt"), "0.001", "2", "1", "0",
                                     "timing", str(Path(temporary) / "final.txt"), "n2"],
                                    capture_output=True, text=True,
                                    env={**os.environ, "OMP_NUM_THREADS": "4", "OMP_THREAD_LIMIT": "1",
                                         "OMP_DYNAMIC": "FALSE"})
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Expected a fixed OpenMP team", result.stderr)

    def test_legacy_sweep_keeps_within_run_error_bars(self):
        self.check_trajectory_and_plot(False, legacy=True)

    def check_trajectory_and_plot(self, parallel, legacy=False):
        suite = "parallel-sweep" if parallel else "solver-sweep"
        binary = ROOT / ("build/bench_steps_openmp" if parallel else "build/bench_steps")
        series = sweep.PARALLEL_SERIES if parallel else tuple((solver, 1) for solver in sweep.SOLVERS)
        load = load_parallel_sweep if parallel else load_solver_sweep
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            folder = directory / "129"
            folder.mkdir()
            instance = folder / "input.txt"
            subprocess.run([sys.executable, str(ROOT / "generate_instance.py"), "--n", "129",
                            "--L", "10", "--w-min", "1", "--w-max", "100", "--v-max", "3",
                            "--a-max", "0", "--seed", "42", "--output", str(instance)],
                           capture_output=True, check=True)
            args = SimpleNamespace(dt=0.0001, steps=3, theta=0.5, order=4, parallel=parallel, samples=3)
            rows = []
            reference = None
            for solver, threads in series:
                case = folder / ("{}-{}threads".format(solver, threads) if parallel else solver)
                summary, final = sweep.collect(binary, instance, case, solver, args, threads)
                command = [str(ROOT / "build/app"), "--solver", solver, str(instance), "0.0001", "3"]
                ordinary = subprocess.run(command, capture_output=True, check=True,
                                          env={**os.environ, "OMP_NUM_THREADS": str(threads),
                                               "OMP_THREAD_LIMIT": str(threads), "OMP_DYNAMIC": "FALSE"})
                self.assertEqual((case / "final.txt").read_bytes(), ordinary.stdout)
                if solver == "n2" and threads == 1:
                    reference = final
                if solver == "n2":
                    for metric in runner.trajectory_errors(final, reference):
                        self.assertLess(metric["rms_error"], 1e-10)
                row = {"bodies": 129, "solver": solver, "threads": threads, "steps": 3,
                       "dt": args.dt, "samples": 1 if legacy else 3, "theta": 0.5 if solver == "barnes-hut" else None,
                       "order": 4 if solver == "fmm" else None, **summary}
                if legacy:
                    with (case / "samples.csv").open() as source:
                        raw = list(csv.DictReader(source))[:1]
                    with (case / "steps.csv").open() as source:
                        raw_steps = [r for r in csv.DictReader(source) if r["sample"] == "1"]
                    runner.write_csv(case / "samples.csv", ["sample", "seconds"], raw)
                    runner.write_csv(case / "steps.csv", ["step", "seconds"],
                                     [{"step": r["step"], "seconds": r["seconds"]} for r in raw_steps])
                    values = [float(r["seconds"]) for r in raw_steps]
                    row.pop("stddev_ms_per_step")
                    row.update(seconds=float(raw[0]["seconds"]),
                               mean_ms_per_step=1000 * float(raw[0]["seconds"]) / args.steps,
                               step_mean_ms=1000 * statistics.mean(values),
                               step_sd_ms=1000 * statistics.stdev(values))
                    case_status = json.loads((case / "status.json").read_text())
                    runner.write_json(case / "status.json", {**case_status, "samples": 1})
                for error in runner.trajectory_errors(final, reference):
                    row[error["field"] + "_relative_l2"] = error["relative_l2_error"]
                    row[error["field"] + "_rms"] = error["rms_error"]
                if parallel:
                    row.update(requested_threads=threads, effective_threads=threads)
                rows.append(row)
            status = {"status": "ok", "suite": suite, "error_bars": "individual_step_sample_sd" if legacy else "whole_run_sample_sd",
                      "settings": {"bodies": [129], "threads": [1, 4, 8] if parallel else 1, "steps": 3, "dt": args.dt,
                                   "samples": 1 if legacy else 3, "warmups": 0, "theta": 0.5, "order": 4}}
            if parallel:
                status["settings"]["configurations"] = [
                    {"bodies": 129, "solver": solver, "threads": threads}
                    for solver, threads in series]
            runner.write_json(directory / "status.json", status)
            path = directory / (suite + ".csv")
            runner.write_csv(path, list(rows[0]), rows)
            self.assertEqual(len(load(directory)[1]), len(series))
            self.assertEqual(prepare_figures([directory])[0][0], prepare_figures([path])[0][0])
            if parallel:
                # Historical sweeps may contain serial FMM only, or no FMM.
                previous = [row for row in rows if row["solver"] != "fmm" or row["threads"] == 1]
                settings = {**status["settings"], "configurations": [
                    point for point in status["settings"]["configurations"]
                    if point["solver"] != "fmm" or point["threads"] == 1]}
                runner.write_json(directory / "status.json", {**status, "settings": settings})
                runner.write_csv(path, list(rows[0]), previous)
                self.assertEqual(len(load(directory)[1]), 7)
                settings.pop("configurations")
                runner.write_json(directory / "status.json", {**status, "settings": settings})
                runner.write_csv(path, list(rows[0]), [row for row in previous if row["solver"] != "fmm"])
                self.assertEqual(len(load(directory)[1]), 6)
                runner.write_json(directory / "status.json", status)
            for bad in (rows[:-1], rows + rows[:1],
                        [{**rows[0], "mean_ms_per_step": 1e9}] + rows[1:],
                        [{**rows[0], ("step_sd_ms" if legacy else "stddev_ms_per_step"): 1e9}] + rows[1:],
                        [{**rows[0], "threads": 2}] + rows[1:]):
                runner.write_csv(path, list(rows[0]), bad)
                with self.assertRaises(ValueError):
                    load(directory)
            runner.write_csv(path, list(rows[0]), rows)
            if parallel:
                original = json.loads((case / "status.json").read_text())
                for field, bad_value in (("effective_threads", 4), ("status", "failed")):
                    runner.write_json(case / "status.json", {**original, field: bad_value})
                    with self.assertRaises(ValueError):
                        load(directory)
                runner.write_json(case / "status.json", original)
                bad = [{**rows[0], "effective_threads": 8}] + rows[1:]
                runner.write_csv(path, list(rows[0]), bad)
                with self.assertRaises(ValueError):
                    load(directory)
                runner.write_csv(path, list(rows[0]), rows)
            samples_path = case / "samples.csv"
            saved_samples = samples_path.read_text()
            samples_path.write_text("sample,seconds\n")
            with self.assertRaisesRegex(ValueError, "whole-run"):
                load(directory)
            samples_path.write_text(saved_samples)
            steps_path = case / "steps.csv"
            steps_path.write_text("sample,step,seconds\n1,1,0.001\n")
            with self.assertRaisesRegex(ValueError, "every step"):
                load(directory)
            status["status"] = "failed"
            runner.write_json(directory / "status.json", status)
            with self.assertRaisesRegex(ValueError, "successful"):
                load(directory)

    def test_failed_driver_records_failure(self):
        with TemporaryDirectory() as temporary, patch.object(runner, "run_logged", return_value=1):
            case = Path(temporary) / "failed"
            with self.assertRaisesRegex(RuntimeError, "Driver failed"):
                sweep.collect(Path("unused"), Path("unused"), case, "n2",
                              SimpleNamespace(dt=0.001, steps=3, samples=3))
            self.assertEqual(json.loads((case / "status.json").read_text())["status"], "failed")


if __name__ == "__main__":
    unittest.main()
