# 3D N-body baseline

A minimal serial C17 simulator. The direct gravity solver evaluates each
unordered pair of distinct bodies once: O(n²) work per step and O(n) memory. Requires Make and a C
compiler on a Unix-like system; no external libraries or OpenMP are needed.
Benchmark orchestration, generation, and the full test suite require Python 3.8+
(standard library only).

```sh
make
./build/app tests/data/two_body.txt 0.01 100 > build/final.txt
make test
```

The command is `./build/app [--solver direct|barnes-hut|fmm] [--theta VALUE | --order P] INPUT DT STEPS`.
`DT` must be positive and finite;
`STEPS` must be a nonnegative integer. The final state goes to stdout. Errors go
to stderr and return a nonzero exit status. Zero steps returns the input state
unchanged, although numeric formatting may change.

The default solver is `direct`. `barnes-hut` uses a serial Barnes–Hut octree:
rebuild the tree, aggregate masses and centers of mass, then discover shared
source lists for target leaves and evaluate each particle. All three solvers share
the Euler loop and instance format. The default
opening parameter is theta = 0.5. Set `--theta` to any finite nonnegative value
with `--solver barnes-hut`; theta zero fully opens the tree. Other Barnes–Hut tree settings can be passed as named options to `solver_create`. Typical
work is approximately O(n log n), but unfavorable distributions or full opening
can approach O(n²).

```sh
./build/app --solver barnes-hut --theta 0.5 tests/data/two_body.txt 0.01 100 > build/final-barnes-hut.txt
```

The serial fast multipole solver (`--solver fmm`) accepts `--order P` from 1
through 6 and defaults to 4. It builds Cartesian multipole expansions on the
shared octree, translates distant cells into local potential expansions, and
evaluates nearby particles directly. Increasing the order retains more expansion
terms. Separation uses a cube enclosing each cell's occupied positions about
its expansion center. At order 4, leaves hold up to 32 particles when subdivision
is possible, and cell pairs with at most 192 source–target particle pairs are
evaluated directly. Larger, well-separated pairs use a specialized 35-coefficient
translation. Other orders use a leaf capacity of 16 and the generic translation.
Use `barnes-hut` for commands that previously selected `approximate`.

```sh
./build/app --solver fmm --order 4 tests/data/1000.txt 0.0001 1 > build/final-fmm.txt
```

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
| `timing` | Whole-run samples; mean, variance, and standard deviation of time per step | macOS or Linux |
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

The measured region is the entire `simulate_with_solver(...)` call, which performs no
allocation or data conversion. State and acceleration workspace are allocated
before measurement. Instance parsing, copying/resetting state, result validation,
deallocation, and file writing are outside the region. Every warmup and sample
starts from the same initial bodies in a warm process. A nonfinite final state
fails the run. Timing uses `CLOCK_MONOTONIC` wall-clock seconds, retains all raw
whole-run samples in `samples.csv` (`sample,seconds`). Reporting divides each
sample by `--steps`, then computes the mean, sample variance (`n - 1` denominator),
and standard deviation. The normalization is by the number of integration steps,
not by `dt` or elapsed simulated time. There are no timers around individual steps.
Profiling and counters never contribute runtime samples.

Both `timing/summary.csv` and the aggregate `timing-summary.csv` record `solver`,
`steps`, `dt`, `samples`, `mean_seconds_per_step`,
`variance_seconds_per_step_squared`, and `stddev_seconds_per_step`; the aggregate
also includes `instance`. Variance has units of (seconds/step)² and scales by
`1 / steps²`; mean and standard deviation scale by `1 / steps`. These columns
replace the previous total-time summary columns. Each mode's `status.json` also
records `steps` and `dt`.

The terminal reports one line per timed input, for example:
`Timing 0001-1000 [approximate]: 0.689 ± 0.158 ms/step (mean ± sample SD; 10 steps, 30 samples)`.
The standard deviation describes variation between whole-run averages per step,
not variation between individual steps within a run.
Older revisions allocated workspace inside `simulate`; exclude that allocation
from the baseline as well when comparing layout changes alone.

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
| `timing-summary.csv` | Time-per-step summaries across instances (`timing` and `all`) |
| `0001-1000/timing/samples.csv` | Raw whole-run durations in seconds |
| `0001-1000/timing/summary.csv` | Mean and standard deviation in seconds/step, variance in (seconds/step)², and run parameters |
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

### Select a solver and evaluate error

Every runner mode accepts `--solver direct|barnes-hut|fmm` (default `direct`),
`--theta` for Barnes–Hut (default 0.5), `--order` for FMM (default 4), and optional `--compare-direct`:

```sh
python3 bench/run.py timing tests/data/two_body.txt --dt 0.0001 --steps 100 \
  --samples 30 --solver barnes-hut --theta 0.5 --compare-direct
```

Solver selection applies to warmups and measured calls. Each mode's `status.json`,
`timing/summary.csv`, and the aggregate `timing-summary.csv` identify the selected
solver. Barnes–Hut status records effective `theta`, `bucket_capacity`,
`max_depth`, and `node_capacity`, including modes unavailable on the host.
FMM status records the expansion order and effective tree settings. Timing summaries
record theta and order, leaving inapplicable fields blank. Plot and terminal labels
show the applicable setting. Older results without theta remain readable unchanged.
The generator, input files, and raw timing sample format stay the same.

The driver writes `solver_config` and `solver_totals` JSON lines to
`stderr.log` outside measured intervals. The runner copies totals into the mode
status's `diagnostics` field. Empty metadata produces no log record. The reader
also accepts the older `approximate_config` and `approximate_totals` prefixes. Their scope is explicitly `warmups_and_samples`:
`calls`, capacity-limited leaf counts, bounds/aggregate fallback counts, and
accepted-node counts accumulate across every completed solver call. `source_node_visits`
counts source nodes visited during shared-list discovery, including leaves and
the root. `direct_evaluations` counts directed particle contributions (excluding
self), and `aggregate_evaluations` counts aggregate contributions evaluated for
individual targets. `accepted_nodes` retains the same meaning as
`aggregate_evaluations`, rather than counting list entries. Direct fallback reports
`count * (count - 1)` directed contributions, even though the direct kernel shares
distance calculations between both directions. Failed modes
retain available diagnostics; unsupported modes have no measured diagnostics.

FMM records the effective `direct_pair_cutoff` in its settings (192 at order 4,
zero at other orders). It reports cumulative `translations` (directed cell-to-cell interactions),
`direct_evaluations`, `capacity_limited_leaves`, and `direct_fallbacks`.

With `--compare-direct`, the runner executes one additional direct trajectory per
input after all measurement modes finish, starting from the archived input with
identical `dt` and `steps`. The archived benchmark executable runs this reference
without timers, warmups, or perf gates. Reference work and Python error evaluation
are outside every measured interval. The internal driver CLI is
`INPUT DT STEPS SAMPLES WARMUPS timing|perf|reference FINAL_FILE [SOLVER [OPTION VALUE]]`;
For example, append `fmm order 4` or `barnes-hut theta 0.5`. `reference` requires one sample, zero
warmups, and the direct solver.

Each input gains `reference/final.txt`, `reference/status.json`, and process logs.
Each successful measurement mode gains `error.csv` and `error-status.json`.
The CSV has a row for position, velocity, and acceleration, with these columns:

- `rms_error`: square root of the mean squared Euclidean error per body.
- `max_error`: maximum Euclidean error over bodies.
- `relative_l2_error`: total difference L2 norm divided by the reference L2 norm.
- `relative_status`: `ok`, or `undefined` when the reference norm is zero and
  the error is nonzero; the relative value is then empty. Both norms zero gives zero.

Empty states have zero error. Bodies are matched by input order; count or mass
changes, malformed/nonfinite states, and nonfinite metrics fail evaluation.
The direct reference is also checked against the original input's counts and masses.
An evaluation failure returns a nonzero overall status but preserves measurement
statuses, samples, and summaries. Check `reference/status.json` and each
`error-status.json` as well as the measurement status. No accuracy threshold is
imposed. Without `--compare-direct`, there is no reference run or evaluation output.

These metrics measure **final trajectory differences**, including accumulated
drift. Final accelerations may be evaluated at different positions once the
approximate solver changes; they do not isolate force error on identical positions.
Barnes–Hut errors depend on the opening parameter and particle distribution;
zero trajectory error is not expected in general.

### Plot saved timings

Install Matplotlib in a virtual environment, then run the root-level plotter:

```sh
python3 -m venv build/plot-venv
build/plot-venv/bin/python -m pip install matplotlib
build/plot-venv/bin/python plot_results.py
```

This reads successful `timing/samples.csv` files with recorded step counts in
`results/`, divides each raw duration by its run's `steps`, and saves
two images per test instance: `plots/raw_time_<instance>.png` and
`plots/speedup_<instance>.png`. For example, inputs named `1000` and `100000`
produce four images, each accumulating all successful runs for that input.
Labels use the suffix after the run
type: `c7a9d7e_20261023-160502.425884Z_timing_+Linking` becomes `+Linking`.
New runs append the solver to the label, such as `+Linking [barnes-hut theta=0.5]`; older
results without solver metadata retain their labels. Runs appear in chronological
order. Each input gets its own baseline and axes;
failed or unfinished timing runs are skipped with a message. Historical runs
without `steps` in their timing status are also skipped with an explanation;
step counts cannot be inferred from saved durations. Present but invalid step
counts are errors. No historical results are modified. Compact figures
wrap long labels to avoid excessive spacing between results.

To compare only specific runs, pass their exact directory names with `--runs`.
Selection happens before loading timings, so unrelated runs are not read. Duplicate
names are ignored; unknown names are errors. The usual chronological order and
input-matching checks still apply. The baseline must be among the selected runs.
Use a separate output directory to retain other plots:

```sh
build/plot-venv/bin/python plot_results.py \
  --runs cf3be1a_20261029-233126.502386Z_timing_+Direct-n12272-steps10 cf3be1a_20261029-233212.934159Z_timing_+BarnesHut-PerTarget-n12272-theta0.5-steps10 \
  --baseline cf3be1a_20261029-233126.502386Z_timing_+Direct-n12272-steps10 \
  --output-dir plots/exact-vs-approximate
```

Time plots use milliseconds per simulation step (`Time per step (ms)`), while
retaining the `raw_time_<instance>.png` filenames. Every run's average per step
appears as a faint dot with horizontal jitter. Separate markers show arithmetic
means, with error bars of ± one sample standard deviation across runs and no
lines connecting means. These bars do not measure variation between individual
steps. Speedup samples are the baseline's mean time per step divided by each
run's time per step; the markers and error bars summarize those ratios.
Baseline uncertainty is not propagated,
and the baseline's own mean ratio can be slightly above 1 because of the
reciprocal transformation. The default reference is the earliest run labelled
`+Baseline` (or `Baseline`), otherwise the earliest run for that input. Override
it with `--baseline LABEL` or a full result directory name. Runs with different
step counts can be compared using normalized times; `steps` and `dt` are saved
in the summaries and mode statuses. Different `dt` values or trajectory lengths
can still change the workload, so normalization does not guarantee identical
conditions. The script checks that each input matches the baseline's archived input.

Use `--results-dir PATH` and `--output-dir PATH` to select other directories.
Re-running updates each instance's two images with all available results; saved
benchmark results are unchanged.

## State, integration, and solver interfaces

Public headers describe contracts; `*_internal.h` headers describe private
representations for implementations. Tests use public headers only.

- `nbody.h`: opaque state ownership, borrowed array views, state copying and
  validation, and the whole-system acceleration callback used by integration.
- `solver.h`: solver registration queries, preparation/destruction, named options,
  and typed settings/diagnostic visitors.
- `io.h`: instance input/output. `parse.h` supplies shared numeric parsing.

The integration layer depends only on state and the acceleration contract. Solver
implementations depend on array views, never on state ownership or scratch-buffer
representation. The executables depend on public interfaces, never on tree nodes
or solver workspace layouts. Benchmark serialization is separate from solver
metadata production.

```c
Simulation *simulation;
if (!simulation_create(&simulation, count)) { return EXIT_FAILURE; }
MutableSimulationView state = simulation_mutable_view(simulation);
/* Initialize all ten arrays through state. */
Solver solver;
const char *error;
if (!solver_create(&solver, "barnes-hut", count, NULL, 0, &error)) {
    simulation_destroy(simulation);
    fprintf(stderr, "%s\n", error);
    return EXIT_FAILURE;
}
simulate_with_solver(simulation, dt, steps, solver.acceleration, solver.context);
SimulationView final = simulation_view(simulation);
/* Consume final state before destroying it. */
solver_destroy(&solver);
simulation_destroy(simulation);
```

`simulate_with_solver` is the only production integration entry point. Once per
step it calls the prepared solver with read-only state and three disjoint output
arrays. The solver fills every output without modifying input or allocating.
Preparation validates configuration and allocates reusable workspace for a fixed
body count; no state pointers are retained. Direct accepts no options.
Barnes–Hut accepts `theta`, `bucket_capacity`, `max_depth`, and `node_capacity`
as named strings; duplicate/unknown options and invalid values fail preparation.
The command line retains its existing `--theta` spelling. All defaults live in
the backend. Zero bodies or zero steps makes no solver call.

The benchmark driver supports `--metadata` to list registered solver names and
`--metadata SOLVER COUNT [OPTION VALUE]` to report effective settings as JSON. This
unmeasured query prepares and releases a solver without integrating; it also works
on hosts without perf. Python reads these settings rather than duplicating tree
configuration rules. Settings are recorded even when profiling is unavailable.
Diagnostics accumulate across warmups and samples; logging and state resets stay
outside measurement. Historical approximate log and result fields remain readable.

### Adding a CPU solver

Implement a `SolverDescriptor` using the private descriptor contract: prepare,
destroy, acceleration, settings, and diagnostics. Only the name and acceleration
callback are required. A stateless solver omits lifecycle callbacks and accepts
no options; absent metadata callbacks produce no fields. Preparation must clean
up on failure. Export the descriptor and
add its address to the static registry in `solver.c`. Add backend-specific
numerical tests. The shared solver tests iterate the registry automatically.
No changes to integration, state ownership, or benchmark measurement loops are
needed. New CLI option spellings require explicit CLI support; the C preparation
interface already accepts named backend options.

### Numerical implementations

Gravity calculations and spatial traversal have separate internal boundaries.
`gravity.c` owns particle interactions, indexed source-range accumulation,
positive-mass monopole construction/evaluation, and the optimized whole-system
direct kernel. Both tree solvers use those source-range and fallback operations.
The direct kernel retains its symmetric unordered-pair optimization; tree source
ranges retain their directed evaluation order.

The octree stores geometry and particle ranges. Barnes–Hut owns a separate array
of `GravityMonopole` summaries, cleared and rebuilt on every evaluation. FMM owns
its moment and local-expansion arrays. All storage is allocated during solver
preparation and reused; no particle-state pointers are retained between calls.

FMM keeps Cartesian term preparation, polynomial shifts, and traversal in
`fmm.c`. The private `laplace.c` module owns the positive `1/r` Taylor recurrence,
distance scaling, and multipole-to-local contraction, including the generated
order-four specialization. Its inputs are term tables, moments, cell scales,
normalized displacement and distance, and scratch/output arrays; it does not
depend on particle state or the octree. Regenerate the specialization with
`python3 generate_fmm_order4.py > src/laplace_order4_internal.h`.

Gravity supplies one source component, mass (`w`). Its FMM field readout takes
the positive gradient of the positive `1/r` expansion; this field is already
acceleration, so target response requires no extra buffer or pass. Source
weighting and field readout are the extension points for future interactions.
There is currently no selectable physics model, charge storage, or magnetic
interaction; the public API and instance format remain gravity-specific.

Each Barnes–Hut call rebuilds a bounded octree without allocation, accumulates
mass and center of mass from leaves upward, and discovers one shared source
interaction list for each target leaf. Each particle in that leaf evaluates the
list independently. After construction, the partition scratch buffer holds the
reusable list of at most `count` source node indices; each entry represents a
disjoint, nonempty source range. Nodes cache their leaf status so evaluation does
not scan child slots for every target.
Default tree settings are theta 0.5, bucket capacity 8, maximum depth 32, and
a node pool of `max(1, 4 * count)`. Depth, capacity, and floating-point subdivision
limits retain larger leaf buckets. Unrepresentable root bounds or aggregate
mass/centers select direct fallback. Workspace diagnostics report capacity-limited
leaves, bounds/aggregate fallback, source-node visits, direct evaluations, and
aggregate evaluations for the latest call.
The public diagnostic visitor reports cumulative totals for the solver lifetime.
Zero bodies or zero steps makes no solver call.

A source node whose particle range overlaps the target leaf is always opened.
Other internal nodes are accepted only when full cube width divided by the
minimum distance from the source center of mass to the target group's actual
particle bounding box is less than theta. This test is conservative for every
target in the group geometrically. The threshold uses ordinary floating-point
arithmetic, so rounding can affect decisions very close to theta. Zero-distance or
unrepresentable distance bounds cause opening. Leaves evaluate particles directly,
skipping the target by index. Theta zero opens every internal node and agrees
with direct summation within rounding error. Grouped lists do not exactly preserve
pairwise force symmetry. Smaller theta generally improves accuracy at higher
cost; errors need not decrease monotonically for each particle.

One allocation holds the opaque `Simulation` header and thirteen contiguous arrays of `count`
doubles, initially ordered as `w`, `x`, `y`, `z`, `vx`, `vy`, `vz`, `ax`, `ay`,
`az`, followed by `next_ax`, `next_ay`, `next_az` workspace. The header stores the
count and active array pointers; the storage follows it. Each step swaps the acceleration
and workspace pointers after evaluating all interactions, eliminating the final
three array copies. Always access state through its named pointers: the first
ten slices of `storage` no longer necessarily contain the active state.
The allocation never moves. Gravity reads mass and position
directly from their four arrays. There is no packing step or production array of
body structs; the instance file format is unchanged.

`simulation_create` checks allocation-size overflow and allocation failure once.
The destination must not already own storage. Initialize all ten state arrays
before simulation; the workspace is overwritten on every step. The file reader
allocates this storage and parses directly into the arrays. Free the allocation
with `simulation_destroy(simulation)`; copying a handle does not duplicate ownership.
An empty simulation has zero count and null pointers.

`simulate_with_solver` updates arrays in place and trusts its caller: storage must come from
`simulation_create`, `dt` must be positive and finite, and the callback must be
nonnull. It performs no argument or bounds checks and cannot fail allocation.
Zero steps or zero bodies is a no-op. The CLI validates arguments before calling
it. `SimulationView` exposes the count and ten state arrays through pointers to
`const double`, excluding workspace and ownership. A view captures pointers, not
array contents. Reacquire views and cached acceleration pointers after calling
`simulate_with_solver`, because their previous arrays may now be workspace. The simulator
refreshes its callback view each step. `simulation_copy` copies only the ten
logical state arrays through their active pointers, including after odd/even
buffer swaps. `simulation_is_finite` likewise ignores inactive workspace.

The direct solver uses `gravity_accelerations` to process each unordered pair
once using scalar loops over `i` and `j > i`. The first
target's source loop initializes the output arrays with each body's first
contribution, including zero for coincident pairs. Later targets accumulate
into those outputs in increasing index order. No separate clearing pass is
needed, and prior workspace values are never read before initialization.

Each pair shares its displacement and inverse distance cubed between both
bodies. The current target's coordinates, mass,
and running acceleration sums stay in local variables throughout the source
loop. These sums start with contributions already received from lower-index
bodies. Each pair immediately updates the source's workspace acceleration;
the target's sums are stored after its source loop. The pair calculations use
scalar arithmetic, with no explicit SIMD.
Input state stays read-only throughout the force calculation. The kernel uses
C `restrict` array pointers to express that its three output arrays do not
overlap one another or any input array, as guaranteed by `simulation_create`.
Callers supplying separate output arrays must honor the same contract.

There are `n(n-1)/2` pair calculations and O(n²) work, with constant auxiliary
storage beyond the existing state and output arrays. Each pair updates both
bodies, so parallelizing this traversal would require coordinating shared
output writes.

Custom whole-system callbacks can use the same integration loop. Tests use an
independent scalar Newtonian reference, accessed through the public solver API.
The file reader enforces positive weights because this executable uses masses.

Each step performs these operations in order:

1. Update all positions using `position += velocity * dt`.
2. Update all velocities using `velocity += acceleration * dt`.
3. Accumulate accelerations in separate workspace using the gravity kernel or
   the supplied callback.
4. Swap the current and next acceleration pointers after all pairs have been
   evaluated, making the summed accelerations active without copying arrays.

Steps 1 and 2 are implemented together per body, before any interaction calls.
Their private update helper uses `restrict` array parameters to express the
storage's non-overlap guarantee. Each step passes the current acceleration
pointers, so the guarantee also holds after buffer swaps.
This is forward Euler: position uses the old velocity. **Input accelerations
drive the first step.** If the file supplies zero accelerations, gravity first
changes velocity on the second step. The example supplies gravity-consistent
initial accelerations. After each step, stored accelerations have been evaluated
at the updated positions and velocities. During that evaluation, callbacks see
the previous acceleration values consistently for every body.

Contributions are combined by vector addition in increasing source-index order.
Mathematically this reduction is associative and commutative; floating-point
roundoff can make results depend on order.
The unordered-pair kernel computes
`inv_r3 = 1 / (r2 * sqrt(r2))`, multiplies it by each mass separately, and updates
both bodies. Each body receives contributions from lower indices first, then
higher indices, preserving increasing source-index order. Each contribution
is added directly to its running sum without intermediate subtotals.
Pair-sharing arithmetic can round differently from two calls to `gravity`;
bitwise equivalence across implementations or compiler configurations is not
guaranteed.

`src/gravity.c` implements Newtonian gravity with normalized units, `G = 1`:

```text
displacement = source.position - target.position
acceleration = source.w * displacement / |displacement|³
```

The target mass cancels when converting force to acceleration. Self-interactions
are skipped by both interaction paths. Distinct bodies at exactly the same position
contribute zero by convention to avoid division by zero; this is not collision
physics. There is no softening or collision handling.

Forward Euler is first-order and does not conserve orbital energy. This is a
baseline for interface and algorithm work, not a high-accuracy orbital solver.
Large time steps, close encounters, or extreme values can cause unstable motion
or floating-point overflow; the simulator does not adapt the time step or
detect numerical overflow during evolution.

## Project layout and tests

- `src/nbody.h` and `src/nbody.c`: array types, callback, and serial simulation.
- `src/state.c`: one-time allocation and read-only array views.
- `src/gravity.c`: gravity kernels, mass monopoles, and the direct solver descriptor.
- `src/barnes_hut.c`: Barnes–Hut solver and private workspace lifecycle.
- `src/fmm.c` implements the fast multipole solver. Both tree solvers use `src/octree.c`.
- `src/laplace.c`: generic and generated order-four `1/r` translations, using private Cartesian term tables.
- `src/solver.h`, `src/solver.c`: shared lifecycle, metadata, and static registry.
- `src/io.h` and `src/io.c`: text input/output and numeric parsing.
- `src/main.c`: command-line arguments and solver preparation.
- `bench/`: separate measurement driver, Linux event gate, and Python runner.
- `tests/test_state.c`, `tests/test_solver.c`: public state, integration, I/O, and solver contract tests, called by one runner.
- `tests/test_support.*`: public-view fixtures and an independent scalar gravity reference.
- `build/`: ignored executables and generated test output.
- `results/`: ignored, persistent experiment artifacts.

The C tests include only public headers. A small table of representative inputs
covers empty/singleton states, unequal masses, coincident particles, and a 3D
cloud. It checks input immutability, complete bounded output writes, solver reuse,
Euler update order and snapshot consistency, state copying after odd/even runs,
configuration validation, fallback/recovery, and file round-tripping and rejection.
Direct and full-opening results are compared to an independent Newtonian reference
with rounding tolerance. A fixed cloud also checks approximate accuracy with a
scaled L2 budget `||error|| <= 0.02 * (1 + ||reference||)` at theta 0.8; this is
not a universal theta error guarantee. Public diagnostics verify that aggregation and fallback were exercised.
All six FMM orders are checked for small-state contracts, approximation reuse,
and fallback/recovery. A 1,000-body, ten-step CLI regression compares all backends
against direct, including separate position, velocity, and acceleration budgets
for default Barnes–Hut and FMM. It does not require monotonic accuracy across FMM
orders, whose tree settings differ. To additionally compare every backend against
a saved executable built before a refactor, run
`NBODY_BASELINE_APP=/absolute/path/to/saved-app make test`; the comparison uses
`||difference|| <= 1e-12 * (1 + ||baseline||)` for each state vector group. Keep
the baseline binary outside `build/` if cleaning between compiler settings.

Tests do not prescribe tree shape, node counts, traversal order, pointer identity,
scratch layouts, or bitwise agreement between numerical implementations. Such
choices can change without changing the public contract. Benchmark tests retain
measurement boundaries, resets, statistics, historical artifact compatibility,
plot units, and failure reporting. The measurement test double uses public state
views only; actual solver selection and results are also checked with real binaries.
Actual Linux perf integration is opt-in:
`NBODY_TEST_PERF=1 make test` on a host with perf access.
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
selected configuration. Production integration uses `simulate_with_solver`;
the test reference evaluates Newtonian gravity independently.
