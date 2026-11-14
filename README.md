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
for every solver, direct `n2` pair evaluation, Barnes–Hut target-leaf evaluation,
and FMM root-target interactions with an OpenMP-capable compiler and runtime:

```sh
make clean && make OPENMP=1
OMP_NUM_THREADS=4 ./build/app --solver barnes-hut tests/data/two_body.txt 0.01 100 > build/final-bh.txt
OMP_NUM_THREADS=4 ./build/app --solver n2 tests/data/two_body.txt 0.01 100 > build/final-n2.txt

# Apple Clang with Homebrew libomp; adjust paths to your installation.
make clean
make OPENMP=1 \
  OMPFLAGS='-Xclang -fopenmp -isystem /opt/homebrew/opt/libomp/include' \
  OMPLIBS='-L/opt/homebrew/opt/libomp/lib -lomp'
```

`OMPFLAGS` and `OMPLIBS` override the compiler and linker OpenMP flags. An
explicit OpenMP build fails if support is unavailable. Clean before switching
between serial and OpenMP builds, including tests. Set `OMP_NUM_THREADS` before
starting the program (or before solver creation when using the API). Both `n2`
and Barnes–Hut fix the maximum team size at creation; smaller runtime teams are
safe. Barnes–Hut reserves one count-element source buffer per worker. For at
least two bodies and multiple workers, `n2` reserves three count-element double
arrays per worker (approximately `24 × N × T` bytes for N bodies and T workers).
Each unordered pair is evaluated once, updating both bodies in private arrays;
the arrays are combined in worker-index order after all pairs finish. Buffers
are reused without allocating during force evaluation. Serial builds and
one-thread execution use the original serial kernel. Parallel `n2` changes
addition order, so results can differ slightly across thread counts through
floating-point rounding; bitwise equality is not guaranteed.

FMM reads the current OpenMP thread setting on each force evaluation. With at
least 1,024 bodies and two occupied root children, it distributes those target
subtrees over at most eight workers. Each job uses small stack-local scratch;
source traversal and accumulation order remain unchanged. Smaller inputs,
single-worker runs, and roots with fewer than two children use the serial
traversal. Uneven subtree populations can limit speedup.

A prepared solver must not be called concurrently. Tree construction, mass
aggregation, FMM upward/downward passes, and FMM's direct fallback remain serial.
Leapfrog uses two static parallel loops per timestep: half-kick and drift,
then the final half-kick. Every position update
finishes before the force solver is called once outside those parallel regions;
the final half-kick starts after force evaluation completes. Parallel scheduling
preserves each body's integration arithmetic and Barnes–Hut/FMM per-target force
accumulation order. Parallel-region launches and `n2` buffer clearing and
combining can add overhead for small inputs; measure complete simulation time
to assess benefit for your workload and thread count.

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
| `<input>/timing/distribution-error.csv` | Six final-state position/acceleration CDF distances with `--compare-n2` |
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

### Barnes–Hut phase experiment: one and four threads

```sh
# Supply OMPFLAGS/OMPLIBS as in the OpenMP build instructions on Apple Clang.
python3 bench/profile_barnes_hut.py bench/data/1000.txt \
  --dt 0.0001 --steps 20 --samples 50 --warmups 3 --theta 0.5 --threads 1 4

# Pure serial profiling without an OpenMP dependency remains available.
python3 bench/profile_barnes_hut.py bench/data/1000.txt \
  --dt 0.0001 --threads 1 --serial

# One aggregate file produces a combined thread-comparison figure.
python3 figures/plot.py results/bh-phases/<run>/phase-summary.csv
```

Pass additional saved inputs to compare body counts/distributions. The default
thread counts are `1 4`; `--threads 4` measures only four workers. Ordinary and
instrumented executables are built once with identical OpenMP/compiler/LTO
settings and reused for all requested thread counts. `OMP_NUM_THREADS` and
`OMP_THREAD_LIMIT` select each team, with `OMP_DYNAMIC=FALSE`. The instrumented
backend records actual minimum/maximum team sizes and total worker slots; the
runner rejects any mismatch. Four threads means four OpenMP workers, not a
promise of physical-core affinity. All runs are sequential.

For every input and thread count, the runner measures baseline, profiled, and
baseline again. State resets, driver behavior, and workloads are identical.
Final states must agree exactly both with and without instrumentation and across
thread counts. Normal builds have no profiling timers. Instrumentation uses
`CPPFLAGS=-DNBODY_PROFILE_BH`; rebuild normally with `make clean && make` afterwards.

Each input contains `threads-1/`, `threads-4/`, etc., with `phase-samples.csv`,
`phase-summary.csv`, and `status.json`. Raw whole-call seconds, raw worker-seconds,
actual thread counts, residuals, baseline drift, and profiling overhead estimates
are retained. The root `phase-summary.csv` includes every input/thread combination;
`thread-summary.csv` reports speedup using the mean of the two uninstrumented
baseline measurements. The figure pairs the thread-count bars for each body
count and annotates total *profiled* ms/step. Old summaries without a `threads`
column remain readable as one-thread profiles.

Membership checks, tree construction, and mass summaries are measured as serial
wall time. Inside parallel force evaluation, each worker privately times target
bounds, traversal/opening-distance checks, and force distances/accumulation.
For each force call, the backend sums each phase's worker elapsed times and
divides by the **actual team size**. These attributed elapsed times measure
average worker occupancy, not separate sequential phase latencies. Adding raw
worker-seconds directly would double-count overlap. The residual of the enclosing
force-evaluation wall time includes scheduling, synchronization, idle/imbalanced
work, startup, and timer overhead. The outer residual includes leapfrog body
updates and other overhead. These mutually exclusive contributions sum to total
runtime; raw worker-seconds are also saved for independent analysis.

Timers surround whole target-leaf operations, never individual interactions.
Warmups are excluded by differencing snapshots outside timing. Per-step values
include initial force evaluation: every sample makes `steps + 1` force calls.
Reuse remains automatic, including across resets; per-build cost and actual
rebuild/reuse counts are saved. Direct and aggregate interaction counts are
retained without separate timers. Confidence intervals describe repeated
whole-run means of the fixed trajectory; they exclude thermal drift, differences
between inputs, and instrumentation bias.

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
python3 figures/plot.py results/current-comparison-<timestamp>

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
point entirely; this also works with `--resume`. The figure generator reports
mixed sample counts or the explicitly skipped point for inclusion in the paper caption. Missing points without an
explicit skip remain errors. Two samples are the minimum for a sample standard deviation.

Results are saved in a new `results/current-comparison-<timestamp>/` directory,
with the same binaries, logs, inputs, raw seconds, final states, and status files
as the native comparison. `comparison.csv` also includes untimed direct-force
error at up to 128 targets. Equal theta does not guarantee equal accuracy.

The plotter requires Matplotlib and rejects incomplete runs, missing points,
invalid samples, and summaries inconsistent with raw data. It writes
`figures/<run-name>-performance.png` using the shared
paper style. Both axes are logarithmic; markers show mean
milliseconds per leapfrog step and error bars show one sample standard deviation
across whole-run averages. Colors identify implementations; solid/dashed lines
identify one/four threads. Nonpositive lower error-bar endpoints, if any, extend
below the logarithmic axis and are reported in the terminal for the paper caption.

## Single-threaded three-solver sweep

```sh
# Three 20-step runs per solver and body count; each resets the saved initial state.
python3 bench/compare_solvers.py
# Use the exact result directory printed above.
python3 figures/plot.py results/solver-sweep/<run>/solver-sweep.csv
```

The default sweep compares **Direct (`n2`), Barnes–Hut (theta 0.5), and FMM
(order 4)** at **1,000, 5,000, 10,000, 20,000, and 50,000 bodies**, with `dt=0.0001`.
Defaults are `--samples 3 --steps 20`, with no warmups. Override `--bodies`,
`--samples` (at least two), `--steps` (at least two), `--dt`, `--theta`, `--order`,
`--seed`, `--distribution`, or `--results-dir` as needed. For a quick workflow check use
`--bodies 16 32 --steps 3`. The runner cleans and builds a dedicated serial
`build/bench_steps` executable; it inherits the compiler and LTO settings.
Run it sequentially with other benchmarks because they share `build/`.

For each body count, the existing generator creates one saved uniform-cube
input (`L=10`, masses 1–100, maximum speed 3, initial acceleration zero,
seed 42 by default). All three solvers use that same file, timestep, and step
count. The timed direct trajectory is also the accuracy reference, so no extra
direct run is needed. Errors are final-trajectory position, velocity, and
acceleration differences, not force errors at identical positions.

Use `--distribution clusters` for eight Gaussian clusters centered at the cube
corners `(±1, ±1, ±1)`, with standard deviation `0.25` per coordinate and zero
initial velocities. Masses and the other experiment settings stay the same.
The saved generator settings record the distribution. This option also works
with `--parallel`.

The plot contains three lines. Each full-run duration is divided by its step
count and converted to milliseconds: `ms_per_step = seconds * 1000 / steps`.
The point is the **mean of the three run averages**, and the bars are **± one
sample standard deviation of those three averages** (denominator `samples - 1`).
The bars describe between-run variation, not variation between the 60 individual
steps, and are not confidence intervals. Every run includes its initial force
evaluation (`steps + 1` force calls). State resets before each measured call;
the prepared solver and its workspace are reused between runs.

Per-step timings are retained as diagnostics, with `sample,step,seconds` columns,
but do not determine the plotted error bars. Step observation adds a clock
reading at each integration boundary in the dedicated binary; its small overhead
is included in total time. Ordinary application and benchmark builds have no
step observer. Allocation, state resets, input generation, validation, and output
are untimed. Saved final states and accuracy metrics describe the last run.

Each run saves its executable, build log, settings and completion status,
`solver-sweep.csv`, and a shared `<N>/input.txt`. Each `<N>/<solver>/` contains
`samples.csv` (three raw whole-run durations), `steps.csv`, `final.txt`, `error.csv`,
`status.json`, and `stderr.log` with solver settings and diagnostics. The plotter
rejects failed/incomplete sweeps and verifies summaries against raw durations.
It writes `figures/solver-sweep.png` at 600 dpi. Historical single-run
sweeps remain readable and retain their explicitly labeled step-to-step SD bars.

## Accuracy–runtime tradeoff

```sh
python3 bench/accuracy_tradeoff.py
# Use four threads for the untimed direct reference (requires OpenMP support):
python3 bench/accuracy_tradeoff.py --reference-threads 4
# Small end-to-end check of the same parameter sweep:
python3 bench/accuracy_tradeoff.py --bodies 128 --steps 2
```

The default experiment uses one saved **1,000-body** uniform-cube instance,
seed 42, `L=10`, masses 1–100, maximum speed 3, and initial acceleration zero.
It compares Barnes–Hut **theta = 0.2, 0.3, 0.5, 0.7** and FMM **orders 1–6**
against one shared, untimed direct (`n2`) reference. Every approximate configuration gets **three measured
runs and one warmup**, with **dt=0.0001 and 500 steps** for trajectories
(total simulated duration **0.05**).
Override `--bodies`, `--dt`, `--steps`, `--samples`, `--warmups`, `--seed`,
`--thetas`, `--distribution`, or `--results-dir` as needed. Use
`--distribution clusters` for the same eight clusters starting from rest described
above. The measured executable is built with
`OPENMP=0` and checked to be serial and uninstrumented. The direct reference runs
once, with no warmup or timing samples. It defaults to one thread; use
`--reference-threads N` to build a separate OpenMP reference executable and use
N threads, while all measured approximate runs remain serial. Both builds inherit
the compiler/LTO settings and the reference build honors `OMPFLAGS`/`OMPLIBS`.
Parallel direct summation changes rounding order, so record the reference thread
count when comparing results. Run benchmarks sequentially because they share `build/`.

The benchmark times the entire `simulate_with_solver` call, including its initial
force evaluation (`steps + 1` calls). It resets the saved state before every run
and retains ordinary prepared-solver reuse. Preparation, allocation, state resets,
validation, and output remain untimed.

Accuracy is measured from the two saved final states using the maximum absolute
empirical CDF difference for each of `x`, `y`, `z`, `ax`, `ay`, and `az`.
Each particle has equal weight, and particle ordering is ignored. Distances range
from zero (identical empirical distributions) to one. `position_ks_max` and
`acceleration_ks_max` are the maximum distances across their respective three
components. The comparison takes O(N log N) time outside the measured simulation.
It requires equal, nonzero body counts and finite states. No initial state,
intermediate trajectory, or force recomputation is used by this calculation.

These are marginal distribution distances, not a statistical equivalence test.
Matching marginals does not establish agreement of joint distributions or
correlations. Acceleration differences include both different final positions
and different force approximations. No p-values or acceptance thresholds are
assigned. The runner uses the same timestep for the direct reference; to study
integration accuracy, separately produce a finer-step reference at the same
final time. State files contain no time metadata, so the caller must ensure
that the two states represent the intended comparison time.

Results live under `results/accuracy-tradeoff/<run>/`: shared `input.txt`, binary,
build/generation logs, overall status, and `accuracy-tradeoff.csv` (approximate
solvers only). The `reference/` folder contains the single direct final state,
status including its thread count, and logs; parallel runs also retain
`reference-bench` and `reference-build.stdout.log`. Each approximate solver
configuration has a `timing/` folder containing raw `samples.csv`, `summary.csv`,
`final.txt`, six-row `distribution-error.csv` (`field,component,ks_distance`),
legacy position-only `error.csv`, `status.json`, and logs. The aggregate CSV
includes both maximum CDF distances. Raw times are seconds per complete measured
call; summaries are seconds/step. Failed experiments retain their error details.
Ordinary `--compare-n2` runs and solver sweeps also save `distribution-error.csv`;
the maxima appear in comparison status and sweep summaries, respectively.

Accuracy sweeps retain CSV summaries for analysis. Use the six-panel final-state
comparison below to visualize distributions; no accuracy–runtime figures are generated.

To overlay all six empirical distributions from any two saved final states,
without rerunning a simulation:

```sh
python3 figures/plot.py --compare-states \
  results/accuracy-tradeoff/<run>/barnes-hut-theta-0.5/timing/final.txt \
  results/accuracy-tradeoff/<run>/reference/final.txt \
  --state-labels "Barnes–Hut" "Direct reference"
```

This writes **`figures/distribution-comparison.png`** at 600 dpi. The top row
compares position components and the bottom row compares acceleration components.
Each panel overlays both ECDFs and reports its distance D; the six distances and
two maxima are also printed. This command needs only the two state files and
works with historical results. Existing result files are not modified.
Add `--cdf-central-fraction 0.98` to write `distribution-comparison-zoom.png`,
zooming each x axis to the combined central 98% ranges of the two samples.
Only the visible range changes; ECDF normalization and distances still use every
particle, and the figure labels the zoom. The default shows the full range.

## Parallel scaling sweep

```sh
# Three 20-step runs per solver, thread count, and body count.
# Supply OMPFLAGS/OMPLIBS as in the OpenMP build instructions on Apple Clang.
python3 bench/compare_solvers.py --parallel
python3 figures/plot.py results/parallel-sweep/<run>/parallel-sweep.csv
# Short workflow check:
python3 bench/compare_solvers.py --parallel --bodies 16 32 --steps 3
```

This mode compares **Direct, Barnes–Hut, and FMM (order 4) at 1, 4, and 8
threads** using one OpenMP executable.
Default body counts are **1k, 5k, 10k, 20k, 50k, and 100k**. Above 50k, Direct
runs only at 8 threads; Barnes–Hut and FMM keep all three thread counts across
every body count. Generator settings, timestep, and default theta match the serial sweep.
Each body count has one saved input shared by its configurations. The Direct
one-thread trajectory supplies the accuracy reference through 50k; above 50k,
the 8-thread Direct trajectory supplies it without an extra run. The summary
records `reference_threads`, and run settings record the expected configurations.
Each configuration defaults to three 20-step runs, resetting state before each,
with no trajectory warmups. Thread counts are fixed, dynamic teams are disabled, and an untimed OpenMP team-formation check
must match the requested count before measurement. This check initializes the
thread runtime but does not execute a simulation warmup.

The runner builds `build/bench_steps_openmp` once, inherits compiler, LTO, and
OpenMP flags, and runs all cases sequentially. Ordinary binaries and the serial
sweep remain unchanged. The nine-line log–log figure plots whole simulation
milliseconds per step averaged over three runs with **± between-run sample SD**,
using the same timing definition as the serial sweep. Error bars are not confidence
intervals; three runs do not guarantee a precise variance estimate. Solver colors
match the serial figure; markers identify the solver across thread counts. Solid,
dashed, and dotted lines identify one, four, and eight threads respectively.
All three FMM lines span every body count; Direct 1T and 4T end at 50k.

Results are saved under `results/parallel-sweep/<run>/`, with
`parallel-sweep.csv`, build log, executable, status, and shared `<N>/input.txt`.
Cases use `<N>/<solver>-<threads>threads/` and retain raw whole-run and step
timings, final states, trajectory errors, logs, and requested/effective thread
metadata. The plotter accepts the CSV or run directory and validates every case
and its raw durations before writing `figures/parallel-sweep.png` at 600 dpi.
The same command also writes `figures/parallel-sweep.csv`, with one row per
plotted point: body count, solver, threads, mean runtime and sample SD in ms/step,
plus steps, timestep, sample count, theta, and FMM order. Historical single-run
data exports its plotted step-to-step SD as `step_sd_ms`; newer runs use
`stddev_ms_per_step` for between-run SD.

## Paper figures

All figure creation lives in **`figures/plot.py`**. Supply one or more saved data
paths; all figures are written directly to `figures/`, regardless of the current
working directory. Benchmark data stays in `results/`. The generator never runs
benchmarks or changes the source data. It replaces the former `plot_results.py`,
`bench/plot_comparison.py`, and the experiment-local phase plotting script.
Figure filenames contain only their type, without commit hashes or timestamps.
Regenerating a type overwrites its previous PNG; raw experiment archives remain
unchanged in `results/`.

```sh
# Recreate the Barnes-Hut phase-share figure.
python3 figures/plot.py \
  results/bh-phases/2b79744_20261111-005459.630660Z_timing/phase-summary.csv

# Generate multiple figure families in one invocation.
python3 figures/plot.py \
  results/bh-phases/2b79744_20261111-005459.630660Z_timing \
  results/current-comparison-<timestamp>

# Compare ordinary timing runs in one parent folder, choosing a baseline.
python3 figures/plot.py results/my-experiment --baseline <baseline-run-name>
```

Install Matplotlib in the Python environment used to run the script. Loaders and
validation tests use only the standard library; drawing uses a headless backend
and a temporary font cache. No TeX installation is required.

| Input | Figures generated |
| --- | --- |
| Aggregate `phase-summary.csv` or its run folder | `phase-shares`: stacked attributed runtime shares, paired by thread count when available |
| Current-suite `comparison.csv` or its run folder | `performance`: solver runtime with sample SD |
| `solver-sweep.csv` or its run folder | `solver-sweep`: three serial solver lines with between-run sample SD |
| `parallel-sweep.csv` or its run folder | `parallel-sweep`: Direct/Barnes–Hut/FMM at 1/4/8 threads; Direct only 8T above 50k |
| `--compare-states ACTUAL REFERENCE` | `distribution-comparison`: six final-state ECDF overlays |
| `timing-summary.csv` or its run folder | `runtime` and `speedup` |
| Parent folder containing ordinary timing runs | `runtime` and `speedup`, comparing runs |

Every figure uses the same NeurIPS-style 5.5-inch width, Times-style serif font
(with portable serif fallbacks), compact typography, and colorblind-friendly
palette. Each figure is saved as one 600-dpi PNG. Figure titles, numbers, workload details, and
captions belong in the manuscript. Include PNGs at their natural text width,
for example `\includegraphics[width=\linewidth]{figures/phase-shares.png}`.
The source experiment settings remain in its `status.json` and saved inputs.

Phase shares divide each mean attributed phase time by the mean total time.
Parallel phases use worker elapsed time divided by team size; their meaning and
residuals are described in the profiling section above. Bounds and remaining
overhead are grouped as in the original chart. Paired profiles must use matching
workloads and include all thread counts for every body count. Ordinary runtime
figures show every whole-run average, the mean, and one sample SD; speedup divides
the fixed baseline mean time per step by each measured time per step. Identical
saved inputs are required for speedup. Legacy runs without recorded step counts
are skipped with an explanation. Comparison figures retain validation against
raw timing samples, including explicitly skipped configurations.

The script validates every supplied dataset before drawing. Duplicate datasets
and conflicting output names are rejected. Supply one dataset per figure type
and one input group for ordinary runtime/speedup plots per invocation. Each
invocation replaces the corresponding fixed-name figures and parallel-sweep CSV.
Generated PNG files and the exported CSV are ignored by Git; the script is tracked.

To add a figure, add its data loader and drawing function to `figures/plot.py`,
then register its input format in `prepare_figures`. Drawing functions return a
Matplotlib figure and use the shared `STYLE`; `save_figure` writes PNGs and
`save_parallel_sweep_table` writes the parallel-sweep CSV.
Add data-validation coverage in `tests/test_figures.py`, update the input table
above, and inspect the generated PNG at paper size. Do not create separate
plotting scripts in benchmark or result directories.

## 3D demonstration videos

`visualize.py` is an optional presentation tool, separate from benchmarking and
paper figures. It uses the existing application without modifying the C solver
or rebuilding automatically. Create a separate environment (Python 3.12 is a
tested choice) and install the optional rendering packages:

```sh
make
python3.12 -m venv results/.visualization-venv
source results/.visualization-venv/bin/activate
python -m pip install -r requirements-visualization.txt

# Approximately one orbit: 421 frames, about 14 seconds at 30 fps.
python visualize.py movie tests/data/two_body.txt \
    --solver n2 --dt 0.01 --steps 1260 --every 3 \
    --output-dir results/two-body-demo
```

The resulting `simulation.mp4` shows planar two-body motion from an oblique 3D
view, with fixed framing, colored particle markers, fading trails, and simulation
time. Marker sizes are illustrative, not physical radii. Defaults are 1920×1080
and 30 fps; use `--width`, `--height` (both even), `--fps`, `--trail-frames`
(`0` disables trails), and `--camera fixed|orbit` to adjust presentation.

Capture and rendering can also be run separately:

```sh
python visualize.py capture tests/data/two_body.txt \
    --solver n2 --dt 0.01 --steps 1260 --every 3 \
    --output-dir results/orbit-capture

python visualize.py render results/orbit-capture/trajectory.npz \
    --camera orbit --output results/orbit-capture/orbit-camera.mp4
```

Each frame is a **faithful replay** from the original input to its cumulative
simulation step. This repeats computation, but preserves uninterrupted solver
behavior, including Barnes–Hut topology reuse within each replay. The script
freezes copies of the executable and input, runs subprocesses sequentially with
one OpenMP thread, and includes the initial and exact final states. It accepts
`--solver n2|barnes-hut|fmm`, `--theta` for Barnes–Hut, `--order` for FMM, and
`--app PATH` for another prebuilt executable. `--dt` and `--steps` are required;
`--every` defaults to 10. Empty inputs are rejected; zero steps captures the
initial state only.

The output directory contains `trajectory.npz`, the original `input.txt`, the
complete `final.txt`, and `capture.stderr.log`. The archive stores float64
positions and masses, int64 step indices, simulation times, and JSON capture
settings; it loads without pickle. Omitting `--output-dir` creates a timestamped
directory under `results/`. Existing directories and output files are never
overwritten. Capture errors do not publish a trajectory; rendering errors leave
completed capture data available for another attempt.

Rendering a saved trajectory requires no simulator executable. It writes one
frame per saved state, without interpolating motion. Simulation time advances by
`dt × every` between ordinary frames, independently of playback fps; the final
interval can be shorter when `steps` is not divisible by `every`. Replay and
render durations are **not solver benchmark measurements**. Use the benchmark
runner for performance studies.

The renderer uses PyVista/VTK and the FFmpeg executable supplied by
`imageio-ffmpeg`; a separate system FFmpeg install is not required. A working
graphics backend is needed even for off-screen rendering. Capture needs only
NumPy beyond the Python standard library. The optional dependencies do not affect
normal builds or benchmark requirements. Run the additional tests explicitly:

```sh
python tests/test_visualize.py
# Also exercise real off-screen rendering, decoding, and publication failures:
NBODY_RENDER_TESTS=1 python tests/test_visualize.py
```

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

To generate eight clusters starting from rest:

```sh
python3 generate_instance.py --distribution clusters --n 1000 \
  --w-min 1 --w-max 100 --v-max 0 --a-max 0 --seed 42
```

This writes `bench/data/instance-clusters-n1000-seed42.txt`. Body `i` is assigned
to corner `i % 8` of `(±1, ±1, ±1)`, with independent Gaussian noise of standard
deviation `0.25` added to each coordinate. For 1,000 bodies, each cluster has
125 members. Cluster positions are not clipped; `--L` is ignored for this mode
and remains required for the default `--distribution uniform`. Mass and vector
options work with either distribution. Reusing the same seed and arguments
reproduces the instance.

To view the full clustered setup and a zoom of one cluster:

```sh
python3 figures/plot.py --cluster-state bench/data/instance-clusters-n1000-seed42.txt
```

This writes `figures/cluster-setup.png`. The zoom shows x–y offsets from the
first cluster's center. Use an initial clustered input, whose rows cycle through
the eight clusters in generator order.

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
| `figures/plot.py`, `figures/` | Central paper-figure generator and its PNG outputs (requires Matplotlib) |

See the [architecture diagram](docs/architecture.mmd) and
[infrastructure diagram](docs/infrastructure.mmd) for module relationships.
Public headers document interfaces; `*_internal.h` files contain implementation
details. [AGENTS.md](AGENTS.md) contains contributor workflow and design constraints.
