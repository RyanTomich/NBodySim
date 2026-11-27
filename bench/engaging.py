"""Prepare, time, and aggregate one Engaging solver sweep on a shared filesystem."""

import argparse
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import tempfile

import compare_solvers as sweep
import run as runner


SMALL_SWEEP_LIMIT = 10000


def assign_tasks(points):
    """Keep small cases together and each larger solver's scaling curve intact."""
    groups = {}
    for point in points:
        key = (point["bodies"], point["solver"]) if point["bodies"] > SMALL_SWEEP_LIMIT else "small"
        point["array_task"] = groups.setdefault(key, len(groups))
    return len(groups)


def worker_nodes(partition):
    # ORCD lists the common 9474F nodes on mit_normal as dual-socket, 48 cores/socket.
    inventory = subprocess.check_output(
        ["sinfo", "-N", "-h", "-p", partition, "-o", "%N %X %Y"], text=True)
    excluded = []
    eligible = False
    for line in inventory.splitlines():
        nodes, sockets, cores = line.split()
        if (sockets, cores) == ("2", "48"):
            eligible = True
        else:
            excluded.append(nodes)
    if not eligible:
        raise ValueError("No dual-socket 48-core/socket nodes listed on " + partition)
    return ["--exclude=" + ",".join(excluded)] if excluded else []


def submit_job(options):
    response = subprocess.check_output(["sbatch", "--parsable", *options], text=True).strip()
    job_id = response.split(";")[0]
    if not job_id.isdigit():
        raise RuntimeError("Unexpected sbatch response: " + response)
    return job_id


def submit(argv):
    args = sweep.parse_args(["--engaging", *argv])
    directory = runner.create_experiment_directory(args.results_dir)
    print("Results:", directory, flush=True)
    status = sweep.sweep_status(args)
    points = status["settings"]["configurations"]
    tasks = assign_tasks(points)
    runner.write_json(directory / "status.json", status)
    try:
        partition = os.environ["SLURM_JOB_PARTITION"]
        selection = worker_nodes(partition)
        shared = os.environ.get("NBODY_SHARED_NODES", "0")
        if shared not in ("0", "1"):
            raise ValueError("NBODY_SHARED_NODES must be 0 or 1")
        status["exclusive_nodes"] = shared == "0"
        if shared == "0":
            selection.append("--exclusive")

        if "NBODY_PROFILE_BH" in os.environ.get("CPPFLAGS", ""):
            raise ValueError("Solver sweeps require phase profiling to be disabled")
        with (directory / "build.stdout.log").open("w") as log:
            subprocess.run(["make", "--no-print-directory", "OPENMP=1", "LTO=1",
                            "print-config", "sweep-binary"],
                           env={**os.environ, "SWEEP_BINARY": str(directory / "bench")},
                           cwd=runner.ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
        for count in args.bodies:
            sweep.generate_input(directory, count, args)

        common = ["--partition=" + partition, "--chdir=" + str(runner.ROOT)]
        workers = submit_job([*common, *selection, "--array=0-{}%4".format(tasks - 1),
                              "--cpus-per-task=" + str(max(p["threads"] for p in points)),
                              "--dependency=afterok:" + os.environ["SLURM_JOB_ID"],
                              "--kill-on-invalid-dep=yes", "--output=" + str(directory / "worker-%A_%a.out"),
                              str(runner.ROOT / "bench/engaging_worker.sbatch"), str(directory)])
        print("Worker array:", workers, flush=True)
        # Run even after a worker times out, so incomplete results get a failed status.
        command = shlex.join([sys.executable, str(runner.ROOT / "bench/engaging.py"), "merge", str(directory)])
        merge = submit_job([*common, "--dependency=afterany:" + workers, "--ntasks=1",
                            "--cpus-per-task=1", "--mem=4G", "--time=01:00:00", "--job-name=nbody-merge",
                            "--output=" + str(directory / "merge-%j.out"), "--wrap=" + command])
        status.update(array_job_id=workers, merge_job_id=merge)
        print("Merge job:", merge, flush=True)
    except (ValueError, OSError, KeyError, RuntimeError, subprocess.CalledProcessError) as error:
        status.update(status="failed", error=str(error), finished=runner.now())
        raise
    finally:
        runner.write_json(directory / "status.json", status)


def case_path(directory, point):
    return directory / str(point["bodies"]) / "{}-{}threads".format(point["solver"], point["threads"])


def worker_placement():
    required = int(os.environ["SLURM_CPUS_PER_TASK"])
    topology = sweep.check_placement(required)
    if len(topology["cpus"]) != required:
        raise ValueError("Expected exactly {} CPUs in the timed step".format(required))
    cpuinfo = Path("/proc/cpuinfo").read_text()
    model = next((line.split(":", 1)[1].strip() for line in cpuinfo.splitlines()
                  if line.startswith("model name")), "unknown")
    if "AMD EPYC 9474F " not in model:
        raise ValueError("Expected an AMD EPYC 9474F worker, found " + model)
    return {**topology, "cpu_model": model}


def worker(directory, task):
    settings = json.loads((directory / "status.json").read_text())["settings"]
    points = [p for p in settings["configurations"] if p["array_task"] == task]
    if not points:
        raise ValueError("Unknown array task {}".format(task))
    args = argparse.Namespace(**settings, parallel=True)
    marker = directory / "task-{}.json".format(task)
    status = {"status": "running", "node": socket.gethostname(), "started": runner.now()}
    runner.write_json(marker, status)
    try:
        status.update(worker_placement())
        previous = None
        for point in points:
            case = case_path(directory, point)
            if case.exists():
                # A requeued element repeats its whole scaling curve on the new node.
                if previous is None:
                    previous = Path(tempfile.mkdtemp(prefix="previous-task-{}-".format(task), dir=directory))
                case.rename(previous / (case.parent.name + "-" + case.name))
            print("N={bodies} {solver}, {threads} threads".format(**point), flush=True)
            summary, _ = sweep.collect(directory / "bench", case.parent / "input.txt",
                                       case, point["solver"], args, point["threads"])
            print("  {mean_ms_per_step:.6g} ms/step".format(**summary), flush=True)
        status["status"] = "ok"
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        status.update(status="failed", error=str(error))
        raise
    finally:
        runner.write_json(marker, {**status, "finished": runner.now()})


def merge(directory):
    status = json.loads((directory / "status.json").read_text())
    settings = status["settings"]
    args = argparse.Namespace(**settings, parallel=True)
    rows = []
    try:
        points = settings["configurations"]
        for task in sorted({p["array_task"] for p in points}):
            marker = directory / "task-{}.json".format(task)
            if not marker.exists() or json.loads(marker.read_text())["status"] != "ok":
                raise ValueError("Incomplete task {}. See its worker log in {}".format(task, directory))
        for count in args.bodies:
            cases = [p for p in points if p["bodies"] == count]
            direct = [p for p in cases if p["solver"] == "n2"]
            reference = reference_threads = None
            if direct:
                baseline = min(direct, key=lambda p: p["threads"])
                reference_threads = baseline["threads"]
                reference = runner.read_state(case_path(directory, baseline) / "final.txt")
            for point in cases:
                case = case_path(directory, point)
                saved = json.loads((case / "status.json").read_text())
                if saved["status"] != "ok":
                    raise ValueError("Incomplete case " + str(case))
                summary = {key: saved[key] for key in ("mean_ms_per_step", "stddev_ms_per_step")}
                final = runner.read_state(case / "final.txt")
                rows.append(sweep.result_row(args, count, point["solver"], point["threads"],
                                             case, summary, final, reference, reference_threads))
        output = directory / "parallel-sweep.csv"
        staging = output.with_suffix(".csv.tmp")
        runner.write_csv(staging, list(dict.fromkeys(key for row in rows for key in row)), rows)
        staging.replace(output)
        status["status"] = "ok"
        status.pop("error", None)
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        status.update(status="failed", error=str(error))
        raise
    finally:
        runner.write_json(directory / "status.json", {**status, "finished": runner.now()})


if __name__ == "__main__":
    action = sys.argv[1]
    if action == "submit":
        submit(sys.argv[2:])
    elif action == "worker":
        worker(Path(sys.argv[2]), int(sys.argv[3]))
    elif action == "merge":
        merge(Path(sys.argv[2]))
    else:
        sys.exit("Expected submit, worker, or merge")
