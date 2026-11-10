# Repository Guidelines

## Project Context

The project compares exact and approximate gravity solvers through one shared
forward Euler simulation loop. `n2` is the direct numerical reference;
Barnes–Hut trades accuracy for speed through theta, and FMM through expansion
order. Solver speed depends on body count and spatial distribution: approximate
solvers need not beat `n2` on small inputs. Evaluate numerical differences and
runtime together using the same saved initial state, timestep, and step count.

The implementation is C17 with OpenMP for parallel Barnes-Hut and the math
library. Building requires Make and a C compiler; generation and benchmarking
use Python 3.8+ and its standard library. See `README.md` for usage and
`docs/architecture.mmd` and `docs/infrastructure.mmd` for architecture diagrams.

## Project Structure & Module Organization

This is a C17/OpenMP 3D N-body simulator with an `n2` O(n²) gravity solver and two approximate
solvers, Barnes–Hut and FMM.

- `src/nbody.h` defines public state views and the simulation interface;
  `src/euler.c` implements forward Euler time stepping.
- `src/state.c` owns private storage and constructs borrowed array views.
- `src/solver.*` defines the shared preparation, destruction, and metadata interface.
  Keep representation details in `*_internal.h`; executables use public headers.
- `src/gravity.c` and `src/gravity_internal.h` provide shared direct gravity interactions and source-range accumulation.
  `src/serial_barnes_hut.c` owns its private mass monopoles and their construction/evaluation helpers.
- `src/n2.c` implements the all-pairs solver; `src/n2_internal.h` exposes its kernel for tree-solver fallbacks.
- `src/io.*` handles instance files.
- `src/octree.c` builds the bounded octree used by `src/serial_barnes_hut.c` and `src/fmm.c`.
- `src/laplace.*` owns expansion table preparation and translations;
  `src/generated/` contains generated order-four arithmetic.
- `src/main.c` parses CLI arguments and prepares the selected solver.
- `tests/test_main.c` runs the public state and solver contract tests in
  `tests/test_state.c` and `tests/test_solver.c`; `tests/test_cli.sh` checks CLI
  behavior, and `tests/test_bench.sh` checks benchmark resets and failure handling.
- `bench/driver.c` measures simulation calls and reports backend settings;
  `bench/run.py` orchestrates timing, optional `--compare-n2` comparisons, and saved results.
- `tests/data/` holds regression fixtures; `bench/data/` holds benchmark inputs.
  Generated `bench/data/instance-*.txt` files are ignored. Keep persistent benchmark
  artifacts in ignored `results/`.

The Makefile automatically discovers C sources directly under `src/` and `tests/`.

## Editing Workflow

Inspect `git status --short` and the relevant source before editing. The working
tree may contain ongoing changes; preserve unrelated edits and existing results.
Keep changes in the owning module listed above, retain public/private header
boundaries, and follow the architecture constraints below. Update README examples
when command-line behavior changes.

For order-four arithmetic changes, edit `generate_fmm_order4.py` and regenerate
the checked-in header instead of editing generated expressions by hand:

```sh
python3 generate_fmm_order4.py > src/generated/laplace_order4_internal.h
```

Add behavioral regression coverage for numerical or interface changes, run
`make test`, and inspect the diff for unintended changes. For benchmark runner
or driver changes, also run the end-to-end smoke benchmark below and inspect its
saved artifacts. Run benchmark invocations sequentially: each rebuilds the shared
`build/` executables. Avoid source edits while a benchmark is running.

## Build, Test, and Development Commands

Run commands from the repository root:

```sh
make                 # Build build/app and build/bench
make test            # Build and run unit, CLI, and benchmark-driver tests
make clean           # Remove build/; preserve results/
./build/app tests/data/two_body.txt 0.01 100 > build/final.txt
```

The application accepts `[--solver n2|barnes-hut|fmm] [--theta VALUE | --order P] INPUT DT STEPS`
and writes the final state to stdout. Theta defaults to 0.5 for Barnes–Hut runs;
explicit theta requires Barnes–Hut and must be finite and nonnegative.
FMM accepts integer orders 1 through 6 and defaults to 4. The benchmark runner accepts the same options.
Use a different output filename from the input. Override compiler settings with
commands such as `make CC=clang`; clean before changing compiler flags.
All targets use `-O3 -g -fno-omit-frame-pointer`, with fast-math disabled and
assertions enabled. Link-time optimization (`-flto`) defaults to `LTO=1` for
every target and requires compiler/linker support; failures do not silently
disable it. Use `make LTO=0` for a non-LTO build, cleaning before switching
settings. `make print-config` reports the selected LTO setting and flags.
Benchmark runs force a rebuild and inherit `LTO` from the environment; prefix
the runner command with `LTO=0` for a non-LTO comparison.

For a clean verification of the current working tree, run
`make clean && make && make test`. This checks the actual uncommitted code as well
as committed code. Require a clean compilation and all three test groups to pass;
a previous successful run does not verify later edits.

## Running Benchmarks and Finding Results

Run from the repository root using Python 3.8+ and a checkout with a valid Git
commit. Use the timing runner so only `simulate_with_solver(...)` is measured.
State/workspace allocation, parsing, reset, validation, and output stay outside
the measured call; `simulate_with_solver` performs no allocation or conversion.

```sh
# Quick end-to-end check of every solver, including saved reference comparisons.
for solver in n2 barnes-hut fmm; do
    python3 bench/run.py timing bench/data/1000.txt \
        --dt 0.0001 --steps 10 --samples 3 --warmups 1 \
        --solver "$solver" --compare-n2 || break
done

# Timing: macOS or Linux; mean and sample variance of time per step over 30 runs.
python3 bench/run.py timing bench/data/*.txt --dt 0.0001 --steps 10 --samples 30
```

The three-sample run verifies the workflow, not a reliable performance ranking.
Use more samples and enough steps for performance studies. Select Barnes–Hut
settings with `--solver barnes-hut --theta 0.5`, or FMM settings with
`--solver fmm --order 4`.

`--dt` and `--steps` are required and positive; `dt` must be finite.
Timing defaults to 10 samples (minimum 2);
the runner defaults to `--warmups 1`; state resets before each call.
Only `timing` is supported. Use `--results-dir PATH` to change the parent directory.

Measure the entire simulation call. Keep `samples.csv` as raw whole-run seconds;
normalize each sample by `--steps` only for summaries and plotting (never by
`dt`). Summary CSVs record solver, theta and order (blank when inapplicable), steps, dt, samples, `mean_seconds_per_step`,
`variance_seconds_per_step_squared`, and `stddev_seconds_per_step`. `timing/status.json` records steps and dt. Approximate runs also record effective tree
settings and available diagnostics, whose totals span warmups and samples. Terminal output and runtime plots use ms/step;
standard deviations describe variation across whole-run averages, not individual
steps. The plotter skips older runs without recorded step counts with an
explanation; preserve those files unchanged.

The printed result path is `results/<short-hash>_<UTC-timestamp>_timing/`.
Each run retains `bench` and combined `build.stdout.log`. Each input gets a
numbered directory, such as `0001-1000/`, containing `input.txt` and a timing folder:

| Mode folder | Main output files |
| --- | --- |
| `timing/` | `samples.csv`, `summary.csv`; aggregate `timing-summary.csv` at the run root |

Read `timing/status.json` and error logs before reporting success. Setup or
build failures produce root `error.txt`. With `--compare-n2`, also check
`reference/status.json`, `timing/error-status.json`, and `timing/error.csv`.
Comparison success means the metrics were computed successfully; no accuracy
threshold is imposed. Metrics describe final trajectory differences, including
accumulated drift, rather than force error at identical positions.

Check that raw sample counts match `--samples`, summaries normalize by steps,
and saved final states are finite. `stderr.log` may contain normal solver
configuration and diagnostic records; inspect its contents rather than treating
every nonempty log as an error. Preserve prior runs; do not recreate
source snapshots, manifests, hash files, or command journals. The directory hash
identifies HEAD; the saved binary also includes any uncommitted code changes.

See [README.md](README.md#benchmarking) for the complete layout.

## Coding Style & Naming Conventions

Use four-space indentation, function opening braces on a separate line, and
control-flow opening braces on the same line. Use `snake_case` for functions and
variables, `PascalCase` for typedefs, and uppercase header guards. Keep private
helpers `static` and read-only callback parameters `const`.

No formatter or linter is configured. Compile cleanly with the existing
`-std=c17 -Wall -Wextra -Wpedantic` flags. Prefer straightforward loops and small
modules over additional abstractions.

## Architecture Constraints

Keep gravity separate from the simulation loop. `AccelerationSolver` fills three
disjoint acceleration arrays without modifying input or allocating. Prepare the
solver once through `solver_create`, reuse it for the same body count, then call
`solver_destroy`. All configuration defaults and validation belong to the backend.
Preserve forward Euler ordering: input accelerations drive the first step, and
new accelerations are published only after the whole solver call. Preserve the
instance format: body count followed by `w x y z vx vy vz ax ay az` rows.
`Simulation` is opaque and owns thirteen contiguous count-element arrays. Create
it with `simulation_create`, initialize the ten state arrays through its mutable
view, and release it with `simulation_destroy`. Reacquire borrowed views after
integration; acceleration buffers swap. The benchmark driver resets logical state through public views.
Validate input and allocation sizes at boundaries; `simulate_with_solver` trusts
storage, positive finite timestep, and callback contracts. Tests use an independent scalar Newtonian reference. Preserve optimized kernel arithmetic and
iteration order when changing interfaces.

## Testing Guidelines

Use standard C `assert` and name test functions `test_<behavior>`. Register new
tests in the appropriate `test_state()` or `test_solver()` suite; register new
suites in `tests/test_main.c`. Application sources except `src/main.c` are
linked automatically. Tests and fixtures must use public headers and views.
Prefer table-driven behavioral cases and an independent numerical reference; do
not assert private layouts, tree structure, pointer swaps, or traversal counts. Add regression tests for numerical or interface changes
and run `make test`. No coverage percentage is configured. Check update order,
callback consistency, malformed input, and saved-state round-tripping when relevant.

## Generating Additional Test Cases

Use `generate_instance.py` when more random inputs are needed; it requires only
Python 3's standard library. Run from the repository root:

```sh
python3 generate_instance.py --n 1000 --L 10 \
  --w-min 1 --w-max 100 --v-max 3 --a-max 2 \
  --seed 42

# Use the generated file directly with the application or benchmark runner.
make
./build/app bench/data/instance-n1000-seed42.txt 0.0001 10 > build/generated-final.txt
python3 bench/run.py timing bench/data/instance-n1000-seed42.txt \
  --dt 0.0001 --steps 10 --samples 3
```

Adjust the body count, bounds, and seed for the case being investigated. Keep
generated benchmark files in `bench/data/`; place only intentional regression
fixtures in `tests/data/`. The default output is
`bench/data/instance-n<N>-seed<SEED>.txt`, relative to the repository. Override
it with `--output PATH`; missing parent directories are created. Record the arguments and seed to reproduce failures, and reuse
the same saved files across implementations for timing comparisons. The output
path is overwritten, so choose distinct filenames for cases you want to retain.

The generator verifies record counts, finite values, bounds, vector magnitudes,
and unique positions before reporting success. Its header contains only `n`, so
the file works directly with the simulator and benchmark runner. `--seed`
controls generation without adding seed metadata to the file.

Use `--w-min > 0` for the gravity simulator; signed property bounds are supported
by the generator but rejected by the current C reader. Random accelerations test
arbitrary input states. For physically consistent cases, compute initial
accelerations from the interaction law. Generator verification establishes input
validity; check simulation correctness on small instances against a trusted
direct O(n²) implementation before using larger instances for timing.

## Commit & Pull Request Guidelines

Use concise, imperative subjects describing one focused change. PRs should
explain the resulting behavior, relevant numerical or interface implications,
and commands used for validation. Link related issues when applicable and update
README examples when usage changes.

## Parallel backend

`src/parallel_barnes_hut.c` owns the ported Rushil solver's private tree/workspace.
Registered names are `n2`, `serial-barnes-hut`, `parallel-barnes-hut`, `fmm`;
`barnes-hut` is a compatibility alias for the serial backend. OpenMP is required.
Keep all backends on the shared Euler integrator. Use `bench/compare_solvers.py`
for a four-way accuracy/performance and thread-scaling study. Tests should run
with `OMP_NUM_THREADS=4` to avoid excessive tiny-case thread overhead.

`src/leapfrog.c` provides optional kick-drift-kick integration. Euler's existing
API and input-acceleration semantics remain unchanged. Benchmark sweeps select
leapfrog with `--integrator leapfrog`; it computes its own starting acceleration,
uses steps+1 solver evaluations, and shares the same solver callback contract.
