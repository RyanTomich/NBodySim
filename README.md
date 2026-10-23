# 3D N-body baseline

A minimal serial C17 simulator. Each step evaluates every ordered pair of
distinct bodies: O(n²) work per step and O(n) memory. Requires Make and a C
compiler on a Unix-like system; no external libraries or OpenMP are needed.
Benchmark orchestration, generation, and the full test suite require Python 3.8+
(standard library only).

```sh
make
./build/app tests/data/two_body.txt 0.01 100 > build/final.txt
make test
```

The command is `./build/app INPUT DT STEPS`. `DT` must be positive and finite;
`STEPS` must be a nonnegative integer. The final state goes to stdout. Errors go
to stderr and return a nonzero exit status. Zero steps returns the input state
unchanged, although numeric formatting may change.

## Instance format

The first value is the number of bodies. Each following row contains ten values:

```text
w x y z vx vy vz ax ay az
```

`w` is a positive scalar mass. Position, velocity, and acceleration each have
three components. All values must be finite. Values are separated by whitespace;
line breaks are for readability. Comments are not supported, numeric tokens are
limited to 127 characters, and missing or extra values are rejected. An empty
instance is simply `0`.

For example, `tests/data/two_body.txt` contains:

```text
2
1 -1 0 0 0 -0.5 0 0.25 0 0
1 1 0 0 0 0.5 0 -0.25 0 0
```

Output uses the same format and `%.17g` precision, so a saved final state can be
used as the next input. Use a different output path from the input path when
redirecting stdout.

## Generate random instances

`generate_instance.py` uses only Python 3's standard library. Run from the
repository root:

```sh
mkdir -p build
python3 generate_instance.py --n 1000 --L 10 \
  --w-min 1 --w-max 100 --v-max 3 --a-max 2 \
  --seed 42 --output build/instance.txt
```

This generates 1,000 bodies with properties in `[1, 100]`, position coordinates
in `[-10, 10]`, velocity magnitudes up to 3, and acceleration magnitudes up to 2.
The same arguments and seed reproduce the same file. The output file is
overwritten, then reopened to verify record counts, finite values, bounds, and
unique positions before reporting success.

The generator writes only `n` on the first line, followed by the ten values per
body shown above, using 17 significant digits. Its output works directly with
the simulator and benchmark tools; the seed is a generation option, not part of
the file format. Use positive property bounds for masses; signed bounds are
supported by the generator for charges. Random accelerations represent
arbitrary input states; physically consistent initial accelerations must be
computed from the interaction law.

## Benchmarking and profiling

All targets use `-O3 -g -fno-omit-frame-pointer`, with fast-math disabled,
assertions enabled, and link-time optimization (`-flto`) enabled by default.
`make` builds both the ordinary
application and a separate benchmark driver linked against the same simulator
sources. The runner force-rebuilds these targets before each experiment and
saves the compiler commands and flags in `build.stdout.log`.

For a comparison without LTO, set `LTO=0` in the runner's environment:

```sh
LTO=0 python3 bench/run.py timing tests/data/1000.txt --dt 0.0001 --steps 10 --samples 30
```

Run these commands from the repository root. The runner builds the executable
automatically; invoke it directly instead of wrapping `build/app` in `perf`,
which would also measure file loading and output.

| Mode | What it measures | Where it runs |
| --- | --- | --- |
| `timing` | Repeated runtime samples, mean, variance, and standard deviation | macOS or Linux |
| `profile` | Sampled CPU call stacks and a function-level report | Linux with usable `perf` |
| `counters` | Hardware events, including cache references and misses | Linux with access to the requested counters |
| `all` | Each of the above in separate executions | Linux for complete results; macOS saves timing plus failure statuses for the other modes |

Run one mode at a time, or use `all`:

```sh
# Runtime statistics across all instance files.
python3 bench/run.py timing tests/data/*.txt --dt 0.0001 --steps 10 --samples 30

# CPU profiling of the 1,000-body instance (Linux).
python3 bench/run.py profile tests/data/1000.txt --dt 0.0001 --steps 100

# General hardware counters, including cache misses (Linux).
python3 bench/run.py counters tests/data/1000.txt --dt 0.0001 --steps 100

# Timing, profiling, and counters across multiple instances (Linux).
python3 bench/run.py all tests/data/*.txt --dt 0.0001 --steps 100 --samples 30
```

`--dt` and `--steps` are required and positive. Timing defaults to 10 samples and
requires at least two. `--samples` applies only to `timing` and the timing portion
of `all`; profiling and counters each measure one simulation call. Every mode
defaults to one unmeasured warmup; change this with `--warmups` (including zero).
Pass one file, several filenames, or a shell glob. Inputs run sequentially.
Each invocation creates a directory under `results/` named
`<short-git-hash>_<UTC-timestamp>_<type>`, for example:

```text
results/d0b5cc2_20261022-053419.808566Z_timing/
```

The timestamp includes year, month, day, time, and microseconds. Types are
`timing`, `profile`, `counters`, and `all`; an `all` run keeps its individual
measurements in separate mode subdirectories. Existing directories are never
overwritten. Use `--results-dir PATH` to select a different parent directory;
the child name is always generated automatically. This replaces the previous
exact-directory `--output` option.

### Measurement boundary

The measured region is the entire `simulate(...)` call, including its workspace
allocation and deallocation. Instance parsing, copying/resetting state, result
validation, and file writing are outside the region. Every warmup and sample
starts from the same initial bodies in a warm process. A nonfinite final state
fails the run. Timing uses `CLOCK_MONOTONIC` wall-clock seconds, retains all raw
samples, and reports mean, sample variance (`n - 1` denominator, seconds²), and
standard deviation. Profiling and counters never contribute runtime samples.

Profiler/counter runs each measure one call, with the same reset and warmup
policy. Linux events start disabled (`--delay=-1`) and the driver uses
`PR_TASK_PERF_EVENTS_ENABLE`/`DISABLE` immediately around the call. The boundary
has small clock/control overhead; very short workloads can be dominated by it.
Increase `--steps` for meaningful timing or profiles. No automatic batching or
cache flushing is performed. Event gating currently targets the serial calling
task; revisit it when adding worker threads or OpenMP.

### Linux setup and interpretation

Timing works on macOS and Linux. Profiling and counters require Linux `perf`
with support for `--delay=-1`. Install the perf tools appropriate to the target
kernel using that system's package manager, and check `perf --version` and
`perf list`. The tools need permission to collect per-task events under the
host's `perf_event_paranoid` policy. The runner never uses sudo, installs
packages, or changes system settings; tool/permission failures are saved.
Run the same Python commands from a checkout inside the Linux environment.
A VM or container must expose the requested performance events; installing
`perf` alone does not guarantee hardware-counter access.

Profiling samples user-space CPU time at 99 Hz and records frame-pointer call
stacks. A profile with no samples is marked failed with instructions to increase
`--steps`. System libraries may lack frame pointers or debug symbols.

Counters default to user-space cycles, instructions, branches, branch misses,
cache references, and cache misses. Override the perf event list when needed:

```sh
python3 bench/run.py counters tests/data/1000.txt --dt 0.0001 --steps 100 \
  --events cycles:u,instructions:u,cache-references:u,cache-misses:u

# Cache-only collection still uses the "counters" mode and directory suffix.
python3 bench/run.py counters tests/data/1000.txt --dt 0.0001 --steps 100 \
  --events cache-references:u,cache-misses:u
```

With cache-only collection, IPC and branch-miss ratios are marked unavailable
because their input events were not requested; the cache-miss ratio remains
available when both cache events are supported and references are nonzero.

Available events and cache meanings depend on the CPU; generic cache events are
not a portable per-cache-level model. Inspect `perf list` for L1/LLC-specific
events. Preserve counter running percentages when comparing multiplexed/scaled
counts. IPC and miss rates are derived only from unambiguous available pairs;
rates are fractions. Missing/unsupported counters have blank values and explicit
statuses, never fabricated zeros. Whole-process elapsed time in raw perf output
is not reported as simulation runtime.

### Saved artifacts

The runner prints `Results: <directory> (ok|failed)`. Open that directory to find
the artifacts. For example, a successful `all` run on `tests/data/1000.txt`
creates this layout:

```text
results/d0b5cc2_20261022-053419.808566Z_all/
  bench
  build.stdout.log
  timing-summary.csv
  0001-1000/
    input.txt
    timing/
    profile/
    counters/
```

Each additional input gets its own numbered directory in command-line order.
A single-mode run creates only that mode's subdirectory. The main files are:

| Location relative to the experiment directory | Expected contents |
| --- | --- |
| `timing-summary.csv` | Runtime summaries across instances (`timing` and `all`) |
| `0001-1000/timing/samples.csv` | Individual durations in seconds |
| `0001-1000/timing/summary.csv` | Mean, sample variance in seconds², and standard deviation |
| `0001-1000/profile/report.txt` | Readable function/call-graph profile |
| `0001-1000/profile/perf.data` | Raw Linux perf recording |
| `0001-1000/profile/sample-events.txt` | Events used to check that the recording contains samples |
| `0001-1000/counters/perf-stat.txt` | Raw perf counter output |
| `0001-1000/counters/counters.csv` | Counts, availability, and counter running percentages |
| `0001-1000/counters/ratios.csv` | IPC, branch-miss rate, and cache-miss rate where available |

Each experiment also contains:

- `bench`: the exact executable used for measurement, saved directly in the
  experiment directory.
- `build.stdout.log`: compiler commands and combined build stdout/stderr.
- One numbered directory per input, containing its archived `input.txt` and
  mode directories. Each mode saves `status.json`, output/error logs, and the last
  successful `final.txt` when available.

Check each mode's `status.json` before interpreting measurements. On failure,
some measurement files may be absent; the status gives the reason and points
to logs such as `stderr.log`, `report.stderr.log`, or `script.stderr.log`.
An unsupported platform or missing `perf` produces a status file without
measurement data. A profile with zero samples or unavailable counters is a
failed measurement, not a successful benchmark with zero cost.

Setup/build failures or interruptions add a short `error.txt` at the experiment
root. Successful runs have no extra experiment metadata. Source snapshots,
manifests, hashes, command journals, and separate build stderr logs are not saved.

The hash in the directory name identifies the current Git commit. Uncommitted
changes are allowed; the saved binary includes the actual code measured, but
reconstructing that source state requires committing those changes. A valid Git
commit is required to name a run. Previous result directories remain unchanged.

The terminal prints the result directory and overall status; detailed evidence
stays in files. A failed mode or input produces a nonzero exit status, while
`all` preserves the results of successful modes. On macOS, `all` saves timing
results and explicit unsupported-platform statuses for Linux-only modes.
`make clean` preserves `results/`. Avoid editing sources during an experiment.

### Plot saved timings

Install Matplotlib in a virtual environment, then run the root-level plotter:

```sh
python3 -m venv build/plot-venv
build/plot-venv/bin/python -m pip install matplotlib
build/plot-venv/bin/python plot_results.py
```

This reads every successful `timing/samples.csv` in `results/` and saves
`plots/raw_time.png` and `plots/speedup.png`. Labels use the suffix after the run
type: `c7a9d7e_20261023-160502.425884Z_timing_+Linking` becomes `+Linking`.
Runs appear in chronological order. Each input gets a separate panel in the
same two figures; failed timing runs are skipped with a message.

Every sample appears as a faint dot with horizontal jitter. Separate markers show
the arithmetic means, with error bars of ± one sample standard deviation and no
lines connecting means. Speedup samples are the baseline's mean time divided by
each measured time; the markers and error bars summarize those ratios.
Baseline uncertainty is not propagated,
and the baseline's own mean ratio can be slightly above 1 because of the
reciprocal transformation. The default reference is the earliest run labelled
`+Baseline` (or `Baseline`), otherwise the earliest run for that input. Override
it with `--baseline LABEL` or a full result directory name. Use comparable `DT`
and `STEPS` across experiments; these parameters are not saved in the timing
CSVs. The script checks that each input matches the baseline's archived input.

Use `--results-dir PATH` and `--output-dir PATH` to select other directories.
Re-running replaces the two images; saved benchmark results are unchanged.

## Simulation and interaction interface

`src/nbody.h` defines `Vec3`, `Body`, and the callback interface:

```c
typedef Vec3 (*Interaction)(const Body *target, const Body *source);

bool simulate(Body *bodies, size_t count, double dt, size_t steps, Interaction f);
```

`Body` contains `w`, `position`, `velocity`, and `acceleration`. The simulator
updates the array in place. `f(target, source)` reads two bodies and returns the
source's contribution to the target's acceleration. It must not mutate state.
When `f == gravity`, the simulator uses `gravity_accelerations` to evaluate each
unordered pair once, sharing its displacement and distance between both bodies.
Pass a different callback to use the generic ordered-pair loop. Gravity
calculations remain in the gravity module. The file reader enforces positive
weights because this executable uses masses.

Each step performs these operations in order:

1. Update all positions using `position += velocity * dt`.
2. Update all velocities using `velocity += acceleration * dt`.
3. Accumulate accelerations in separate workspace using the gravity kernel or
   the supplied callback.
4. Assign the summed accelerations after all pairs have been evaluated.

Steps 1 and 2 are implemented together per body, before any interaction calls.
This is forward Euler: position uses the old velocity. **Input accelerations
drive the first step.** If the file supplies zero accelerations, gravity first
changes velocity on the second step. The example supplies gravity-consistent
initial accelerations. After each step, stored accelerations have been evaluated
at the updated positions and velocities. During that evaluation, callbacks see
the previous acceleration values consistently for every body.

Contributions are combined by vector addition in increasing source-index order.
Mathematically this reduction is associative and commutative; floating-point
roundoff can make results depend on order.
The gravity kernel weights the shared vector by each source mass separately.
Its arithmetic can round differently from two calls to `gravity`; bitwise
equivalence is not guaranteed.

`src/gravity.c` implements Newtonian gravity with normalized units, `G = 1`:

```text
displacement = source.position - target.position
acceleration = source.w * displacement / |displacement|³
```

The target mass cancels when converting force to acceleration. Self-interactions
are skipped by the simulator. Distinct bodies at exactly the same position
contribute zero by convention to avoid division by zero; this is not collision
physics. There is no softening or collision handling.

Forward Euler is first-order and does not conserve orbital energy. This is a
baseline for interface and algorithm work, not a high-accuracy orbital solver.
Large time steps, close encounters, or extreme values can cause unstable motion
or floating-point overflow; the simulator does not adapt the time step or
detect numerical overflow during evolution.

## Project layout and tests

- `src/nbody.h` and `src/nbody.c`: body types, callback, and serial simulation.
- `src/gravity.h` and `src/gravity.c`: gravitational interaction only.
- `src/io.h` and `src/io.c`: text input/output and numeric parsing.
- `src/main.c`: command-line arguments and selection of the gravity callback.
- `bench/`: separate measurement driver, Linux event gate, and Python runner.
- `tests/`: assertion-based numerical/interface tests, CLI tests, and instances.
- `build/`: ignored executables and generated test output.
- `results/`: ignored, persistent experiment artifacts.

`make test` checks gravity, Euler update order, supplied initial accelerations,
callback snapshots and pair counts, input validation, round-tripping, and
resuming saved states. Python tests also cover benchmark boundaries, warmup/state
resets, statistics, output organization, and tool failures. Actual Linux perf
integration is opt-in: `NBODY_TEST_PERF=1 make test` on a host with perf access.
An explicit skip in ordinary tests does not validate Linux event collection.
That integration test checks profiling and cycles/instructions using temporary
files under `build/`, which it removes afterward. Use the `bench/run.py` commands
above when results need to persist under `results/`. Do not describe mocked or
skipped integration tests as successful hardware measurements.

Run commands from the repository root. `make clean` removes build artifacts.
The Makefile discovers application C sources directly under `src/`, unit tests
directly under `tests/`, and benchmark sources directly under `bench/`.
Override `CC`, `CFLAGS`, `CPPFLAGS`, `LDFLAGS`, or `LDLIBS` as needed; mandatory
optimization/debug/frame-pointer flags are appended to `CFLAGS` for every target.
LTO defaults to `LTO=1` and adds `-flto` to the combined compile-and-link command
for every target. The compiler and linker must support LTO; a failure is reported
without silently disabling it. To build with an unsupported toolchain or compare
against the non-LTO configuration, use:

```sh
make clean
make LTO=0
```

Run `make clean` before manual builds with different compiler settings, including
switching `LTO` between `0` and `1`; the benchmark runner always forces a fresh
build. Use `make print-config` (or `make LTO=0 print-config`) to inspect the
selected configuration. The simulation and callback interfaces are unchanged.
