# 3D N-body simulator

A C17 simulator for comparing direct (`n2`), Barnes–Hut, and fast multipole
(FMM) gravity solvers. All backends use the same kick-drift-kick leapfrog loop
and input format. Compare runtime and numerical differences using the same
saved input, timestep, and step count; approximate solvers are not necessarily
faster on small inputs.

Requires Make and a C compiler with the standard math library. Benchmarking and
instance generation use Python 3.8+ with no extra packages. Run the commands
below from the repository root.

## Build and run

```sh
make

# Direct all-pairs reference (also the default backend).
./build/app --solver n2 tests/data/two_body.txt 0.01 100 > build/final-n2.txt

# Barnes–Hut: smaller theta generally improves accuracy at greater cost.
./build/app --solver barnes-hut --theta 0.5 tests/data/two_body.txt 0.01 100 > build/final-bh.txt

# FMM: higher expansion orders retain more terms.
./build/app --solver fmm --order 4 tests/data/two_body.txt 0.01 100 > build/final-fmm.txt
```

The positional arguments are `INPUT DT STEPS`. The timestep must be positive and
finite; steps must be a nonnegative integer. Barnes–Hut defaults to theta `0.5`
and accepts finite values of zero or greater. FMM defaults to order `4` and
accepts integers from `1` through `6`. Use `--theta` only with Barnes–Hut and
`--order` only with FMM.

The final state is written to stdout in the input format, so it can be used to
resume a simulation. Direct-force restarts reproduce uninterrupted runs; a newly
prepared approximate tree can differ from a reused tree.
**Always use a different output filename from the input.**

Builds enable link-time optimization by default. If your compiler does not
support it, run `make clean && make LTO=0`. Clean before changing compiler
settings; `make print-config` shows the active configuration.

Builds are serial by default (`OPENMP=0`). To parallelize leapfrog body updates
for every solver and Barnes–Hut target-leaf force evaluation with an
OpenMP-capable compiler and runtime:

```sh
make clean && make OPENMP=1
OMP_NUM_THREADS=4 ./build/app --solver barnes-hut tests/data/two_body.txt 0.01 100 > build/final-bh.txt

# Apple Clang with Homebrew libomp; adjust paths to your installation.
make clean
make OPENMP=1 \
  OMPFLAGS='-Xclang -fopenmp -isystem /opt/homebrew/opt/libomp/include' \
  OMPLIBS='-L/opt/homebrew/opt/libomp/lib -lomp'
```

`OMPFLAGS` and `OMPLIBS` override the compiler and linker OpenMP flags. An
explicit OpenMP build fails if support is unavailable. Clean before switching
between serial and OpenMP builds, including tests. Set `OMP_NUM_THREADS` before
starting the program (or before solver creation when using the API). Barnes–Hut
reserves one count-element source buffer per worker and fixes the maximum team
size at creation; smaller runtime teams are safe. A prepared solver must not be
called concurrently. Tree construction, mass aggregation, and the `n2` and FMM
force evaluations remain serial. Leapfrog uses two static parallel loops per
timestep: half-kick and drift, then the final half-kick. Every position update
finishes before the force solver is called once outside those parallel regions;
the final half-kick starts after force evaluation completes. Parallel scheduling
preserves each body's arithmetic and each target's force accumulation order.
The two integration parallel-region launches can add overhead for small inputs;
measure complete simulation time to assess benefit for your workload and thread
count.

## Benchmarking

The runner rebuilds the executables, resets the input before every run, and times
only the simulation call. Parsing, setup, resets, and output are excluded.
Initial acceleration is recomputed inside each measured call: a nonempty run
uses `steps + 1` force evaluations. Whole-call time is still divided by `steps`.

```sh
# Direct reference timing.
python3 bench/run.py timing bench/data/1000.txt \
  --dt 0.0001 --steps 10 --samples 30

# Approximate backends, with final-state comparisons against n2.
python3 bench/run.py timing bench/data/1000.txt \
  --dt 0.0001 --steps 10 --samples 30 \
  --solver barnes-hut --theta 0.5 --compare-n2

python3 bench/run.py timing bench/data/1000.txt \
  --dt 0.0001 --steps 10 --samples 30 \
  --solver fmm --order 4 --compare-n2
```

Pass multiple input files or `bench/data/*.txt` to benchmark several instances.
Both `--dt` and `--steps` are required and positive. Defaults are 10 samples
(minimum 2) and one warmup; set `--warmups` to change the latter. For a quick
workflow check, use `--samples 3`; use more samples and enough steps for
performance comparisons. Run benchmarks sequentially because they share the
build directory. Prefix a command with `LTO=0` for a non-LTO benchmark.

The runner requires a checkout with a Git commit and prints its result directory
under `results/<short-hash>_<UTC-timestamp>_timing/`. Use `--results-dir PATH` to
choose another location. Start with these saved files:

| File | Contents |
| --- | --- |
| `timing-summary.csv` | Per-step timing summaries for all inputs |
| `<input>/timing/samples.csv` | Raw whole-run durations in seconds |
| `<input>/timing/status.json` | Success/failure and run settings |
| `<input>/timing/final.txt` | Final state from the timed run |
| `<input>/timing/error.csv` | Position, velocity, and acceleration differences with `--compare-n2` |

Here `<input>` is a numbered folder such as `0001-1000`. Each run also saves the
input files, measured executable, build log, and per-input logs. Terminal timings
are in milliseconds per step; summaries divide whole-run durations by the step
count. Check status files and logs before using results. With `--compare-n2`, also
check `reference/status.json` and `timing/error-status.json` inside each input
folder. Comparisons measure final trajectory differences, including accumulated
drift; success means metrics were computed, not that an accuracy threshold passed.
New timing and reference status files identify `"integrator": "leapfrog"`;
historical Euler results are preserved unchanged.

### Compare the original and merged implementations

`bench/compare_native.py` compares the current Barnes–Hut implementation with
original `src` and Rushil sources from commit `25f7df2`, using each one's native
integrator. It builds isolated executables without changing the working tree.
Original `src` remains serial and uses Euler; merged and Rushil use leapfrog.

```sh
# Supply OMPFLAGS/OMPLIBS as above on Apple Clang.
python3 bench/compare_native.py
```

Defaults are 10,000/20,000/50,000 bodies, requested thread counts 1/4/8,
theta 0.5, dt 0.0001, and ten steps. The same seed-42 saved inputs contain
uniform positions in [-1,1], total mass one, and zero velocities/accelerations.
Euler uses the supplied acceleration and makes ten force calls; leapfrog
reinitializes acceleration and makes eleven. Equal theta does not mean equal
accuracy between the different trees.

Each configuration gets 30 samples in shuffled five-sample batches, with two
warmups per batch. Configurations with relative standard error above 2% extend
to 60 samples. Timing covers native integration and forces, including tree work;
fresh workspace preparation, state resets, validation, and I/O stay outside it.
Every run validates finite state, masses, and repeatability. Final states must
also match across thread settings within each implementation. Independent direct
forces at 128 evenly spaced targets provide an untimed accuracy check.

Results go into a new `results/native-barnes-hut-<timestamp>/` directory with
executables, build logs, saved inputs, per-batch samples/final states/status,
`comparison.csv`, and root `status.json`. Samples are whole-run seconds;
summaries are mean, standard deviation, and median milliseconds per step.
Use `--machine-description` to attach CPU/core/memory details. A short harness
check uses `--bodies 129 --steps 3 --samples 2 --max-samples 2 --batch-size 2`.

### Plot current n², Rushil, and src performance

```sh
# Supply the same OMPFLAGS/OMPLIBS overrides as above on Apple Clang.
python3 bench/compare_native.py --suite current \
  --machine-description 'CPU model, core count, memory'
python3 bench/plot_comparison.py results/current-comparison-<timestamp>

# Quick adapter/workflow check before the full experiment:
python3 bench/compare_native.py --suite current --bodies 129 --steps 3 --samples 2
```

The current suite builds the actual working-tree sources, including uncommitted
changes. It measures n² with one thread and both Rushil and src Barnes–Hut with
one and four threads, at 1,000, 10,000, and 100,000 bodies: 15 configurations.
Defaults are ten measured runs and two warmups per configuration, ten leapfrog
steps, dt 0.0001, and theta 0.5. All builds use the same compiler, OpenMP, `-O3`,
LTO, assertions, and disabled fast-math. The n² force kernel stays serial.
Rushil retains its native loop and default task cutoff of 256.

Each body count uses one saved seed-42 input, with uniform positions in [-1,1],
equal masses totaling one, and zero velocities/accelerations. Fresh workspace
preparation and state reset happen outside timing; native integration, initial
acceleration, and all tree work happen inside it (eleven force calls per run).
Configurations run sequentially in a seed-42 shuffled order. Each has one batch;
there is no automatic sample extension or outlier removal. `--samples`, `--bodies`,
`--steps`, `--dt`, `--theta`, and `--warmups` can customize the workload; current
mode requires `--threads 1 4`, with batch size and maximum samples equal to samples.
The historical suite remains the default when `--suite` is omitted.

For an expensive largest n² point, add `--large-n2-samples 2` to reduce only that
point to two measured runs, retaining two warmups and ten samples elsewhere.
An interrupted current-suite run can be continued with
`python3 bench/compare_native.py --resume results/current-comparison-<timestamp> --large-n2-samples 2`.
Resume uses saved settings and binaries, revalidates completed points, and keeps
interrupted attempt logs. Only `--large-n2-samples` and `--skip-large-n2` override
saved settings; they should target an unfinished point. Use `--skip-large-n2` to omit that
point entirely; this also works with `--resume`. The plot caption identifies
mixed sample counts or the explicitly skipped point. Missing points without an
explicit skip remain errors. Two samples are the minimum for a sample standard deviation.

Results are saved in a new `results/current-comparison-<timestamp>/` directory,
with the same binaries, logs, inputs, raw seconds, final states, and status files
as the native comparison. `comparison.csv` also includes untimed direct-force
error at up to 128 targets. Equal theta does not guarantee equal accuracy.

The plotter requires Matplotlib and rejects incomplete runs, missing points,
invalid samples, and summaries inconsistent with raw data. It writes
`plots/performance.png`, `plots/performance.svg`, and a verified `plots/comparison.csv`
inside the run directory. Both axes are logarithmic; markers show mean
milliseconds per leapfrog step and error bars show one sample standard deviation
across whole-run averages. Colors identify implementations; solid/dashed lines
identify one/four threads. Nonpositive lower error-bar endpoints, if any, extend
below the logarithmic axis and are noted in the caption.

## Run tests

```sh
make test

# Rebuild everything from scratch and run all tests.
make clean && make && make test

# Test four OpenMP threads, then one (supply platform overrides if needed).
make clean && OMP_NUM_THREADS=4 make OPENMP=1 test
OMP_NUM_THREADS=1 ./build/tests
```

This runs the C state and solver tests, command-line tests, and benchmark-driver
tests. `make clean` removes `build/` and preserves saved benchmark results.
On Apple Clang, pass the same `OMPFLAGS` and `OMPLIBS` overrides from the build
example to `make OPENMP=1 test`.

## Make test instances

```sh
python3 generate_instance.py --n 1000 --L 10 \
  --w-min 1 --w-max 100 --v-max 3 --a-max 2 \
  --seed 42

./build/app --solver fmm bench/data/instance-n1000-seed42.txt \
  0.0001 10 > build/generated-final.txt
```

This creates `bench/data/instance-n1000-seed42.txt` with 1,000 bodies, position
coordinates in `[-10, 10]`, masses in `[1, 100]`, velocity magnitudes up to `3`,
and acceleration magnitudes up to `2`. Reusing the arguments and seed reproduces
the input. Use `--output PATH` to choose a destination; existing output files are
overwritten. Generated files work directly with both the app and benchmark runner.

The instance format is a body count followed by one row per body:

```text
w x y z vx vy vz ax ay az
```

Mass `w` must be positive, and all values must be finite. The remaining values
are position, velocity, and acceleration components. Saved acceleration fields
round-trip through I/O, but **every positive-step run recomputes initial
accelerations from current masses and positions**. Generated random acceleration
values therefore do not affect integration. Zero-step runs preserve all fields.

## Architecture and project map

The command-line app loads a state, prepares the selected gravity backend, and
runs a shared leapfrog loop. After initializing acceleration, each step applies a
half-kick to velocity, drifts positions, evaluates new acceleration for the whole
system, and completes the half-kick. Returned velocities are full-step values and
accelerations correspond to final positions. State ownership, integration, and
gravity calculations are separate, so changing backends preserves the simulation
workflow. Barnes–Hut and FMM share an octree for grouping nearby bodies.

Barnes–Hut reuses tree topology while every body remains inside its assigned
leaf, checking every member of each bucket. Uncertain boundary checks trigger a
rebuild. Mass summaries, target bounds, and source lists are refreshed on every
force call, including calls with changed masses or a different state allocation.
Topology reuse is automatic; rebuilding occurs only when the topology is
unavailable or membership validation fails. At nonzero theta, reused and rebuilt
trees may produce different approximations. FMM rebuilds on every force call.
Barnes–Hut's cumulative `capacity_limited_leaves` diagnostic counts the current
topology's limited leaves on each evaluation, including reuse. Barnes–Hut uses
direct particle calculations at leaves, but never switches to the all-pairs
solver. If tree bounds or mass summaries cannot be represented, it fills all
acceleration outputs with `NaN`. A later call with valid data can recover.

Gravity uses normalized units (`G = 1`). Leapfrog is second-order for smooth
forces; it does not conserve energy exactly, and approximate tree forces do not
guarantee symplectic behavior. There is no adaptive timestep, softening, or
collision handling. Solver callbacks must compute acceleration from masses and
positions only: velocities seen after a drift are half-step values, and previous
accelerations must not affect the force. Empty states and zero-step calls perform
no force evaluations. Every later positive-step invocation reinitializes forces
to account for changes through mutable state views.

| Location | Purpose |
| --- | --- |
| `src/main.c`, `src/io.*` | Command-line entry point and instance I/O |
| `src/nbody.h`, `src/state.c`, `src/leapfrog.c` | Public simulation interface, state ownership, and time stepping |
| `src/solver.*` | Backend selection and lifecycle |
| `src/n2.c`, `src/barnes_hut.c`, `src/fmm.c` | Gravity backends |
| `src/gravity.c`, `src/octree.c`, `src/laplace.*` | Shared interactions, spatial tree, and FMM expansion calculations |
| `bench/run.py`, `bench/driver.c` | Benchmark orchestration and measurement |
| `tests/`, `tests/data/` | Regression tests and small fixtures |
| `generate_instance.py`, `bench/data/` | Input generation and benchmark instances |
| `build/`, `results/` | Generated executables and persistent benchmark output |
| `plot_results.py` | Plot saved timings (requires Matplotlib) |

See the [architecture diagram](docs/architecture.mmd) and
[infrastructure diagram](docs/infrastructure.mmd) for module relationships.
Public headers document interfaces; `*_internal.h` files contain implementation
details. [AGENTS.md](AGENTS.md) contains contributor workflow and design constraints.
