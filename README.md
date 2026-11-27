# 3D N-body simulator

Simulate gravity in three dimensions and compare five solvers on the same initial conditions.

## Run a simulation

You need Python 3.8+, Make, and a C17 toolchain with OpenMP support.
Run this from the root of a Git checkout to build the simulator and time a two-body orbit.

```sh
python3 bench/run.py timing tests/data/two_body.txt \
  --solver n2 --dt 0.01 --steps 100 --samples 3
```

- **Initial conditions.** Use `tests/data/two_body.txt` for the orbit, `bench/data/1000.txt` for a 1,000-body system, or your own file in the final-state format below.
- **Duration.** The example advances 100 steps of 0.01 time units each. Each of the three timing samples starts from the same input.
- **Integrator.** Leapfrog is the default. Use `--integrator euler` for forward Euler, which uses the input accelerations for its first velocity update.

| `--solver` | How it computes gravity |
| --- | --- |
| `n2` (default) | Evaluates every pair directly. Use it as the numerical reference. |
| `serial-barnes-hut` | Approximates distant groups of bodies using a serial tree walk. |
| `parallel-barnes-hut` | Runs independent tree walks for individual bodies in parallel. |
| `grouped-barnes-hut` | Shares one tree walk among nearby target bodies and processes groups in parallel. |
| `fmm` | Uses fourth-order fast multipole expansions for distant interactions. |

For Barnes–Hut, `--theta` defaults to `0.5`. Smaller values generally improve accuracy by evaluating more interactions in detail.

## Read the results

Open the directory printed after `Results:` to find the files from the two-body run.

| Result | Where to find it |
| --- | --- |
| Timing | The terminal reports the mean ± sample standard deviation in milliseconds per step. `timing-summary.csv` records both in seconds per step. |
| Final state | `0001-two_body/timing/final.txt` contains the state after step 100. The first line is the body count. Each remaining row lists mass, position, velocity, and acceleration in the order `w x y z vx vy vz ax ay az`. |

## Run the Engaging scaling benchmark

From the repository root on the shared filesystem, submit the preparation job:

```sh
sbatch bench/engaging.sbatch
```

Preparation builds one GCC/OpenMP binary in the result directory, generates one
input per size, and submits the worker array and merge job. The default array has
seven elements. Element 0 runs all three solvers at 1k, 5k and 10k bodies.
Elements 1–3 run `n2`, grouped Barnes–Hut and FMM respectively at 50k bodies.
Elements 4–6 repeat that split at 100k. Each scaling curve runs on one node, in
ascending thread count. The teams are 1, 4, 8 and 16, with the 1-thread direct
case skipped at 50k and the 1- and 4-thread direct cases skipped at 100k.
Each point uses one warmup and three measured 20-step leapfrog trajectories.

Workers request exclusive EPYC 9474F nodes on `mit_normal`, with up to four
workers running concurrently as resources permit. Each timed step requests a
compact allocation on physical cores, local memory and the performance CPU
governor. Workers verify the CPU model, single-socket affinity and effective
OpenMP team sizes. The governor is a request to the site. Check the observed
governors in `task-<index>.json`. Exclusive jobs can consume the full 96-core
base allowance on one node. For more concurrency under that allowance,
`NBODY_SHARED_NODES=1 sbatch bench/engaging.sbatch` requests only the cores
needed per worker. Jobs on shared nodes can contend for cache and memory bandwidth.
[Engaging node inventory](https://orcd-docs.mit.edu/running-jobs/available-resources/) ·
[Partition limits](https://orcd-docs.mit.edu/running-jobs/overview/)

Pass workload options after the script name, for example
`sbatch bench/engaging.sbatch --samples 10 --steps 40`. The preparation log
`slurm-nbody-prepare-<jobid>.out` prints the result path and both downstream job
IDs. The merge runs after all workers finish and checks every expected task. On
success, `parallel-sweep.csv` contains all timings and comparisons to the
smallest eligible `n2` team for each size. For timing-only solver selections,
the reference fields are empty. Check the run's `status.json`, the per-task
status files, and `worker-<arrayid>_<index>.out` for failures.

## Contribute to a solver

Edit the gravity algorithm in [src/](src/). To add a solver, declare its descriptor in [solver.h](src/solver.h) and register it in [solver.c](src/solver.c).

Add a behavioral regression case in [test_solver.c](tests/test_solver.c) and run `make test`. For an approximate solver, add `--compare-n2` to the timing command to measure final-state differences. To compare runtime, time `n2` separately with the same input, integrator, timestep, and step count. See the [repository guidelines](AGENTS.md) for contribution conventions.
