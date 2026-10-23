# Repository Guidelines

## Project Structure & Module Organization

This is a serial C17 3D N-body simulator with O(n²) work per step.

- `src/nbody.h` and `src/nbody.c` define array state types and simulation logic;
  `src/state.c` allocates their storage and constructs read-only views.
- `src/gravity.*` implements gravity; `src/io.*` handles instance files.
- `src/main.c` parses CLI arguments and selects the interaction callback.
- `tests/test_main.c` contains unit tests; `tests/test_cli.sh` checks CLI behavior.
- `bench/` contains the separate measurement driver and Python runner;
  `tests/bench/` tests these integrations.
- `tests/data/` holds instances. Keep generated inputs in ignored `build/` and
  persistent benchmark artifacts in ignored `results/`.

The Makefile automatically discovers C sources directly under `src/` and `tests/`.

## Build, Test, and Development Commands

Run commands from the repository root:

```sh
make                 # Build build/app and build/bench
make test            # Run C, shell, and Python tests
make clean           # Remove build/; preserve results/
./build/app tests/data/two_body.txt 0.01 100 > build/final.txt
```

The application accepts `INPUT DT STEPS` and writes the final state to stdout.
Use a different output filename from the input. Override compiler settings with
commands such as `make CC=clang`; clean before changing compiler flags.
All targets use `-O3 -g -fno-omit-frame-pointer`, with fast-math disabled and
assertions enabled. Link-time optimization (`-flto`) defaults to `LTO=1` for
every target and requires compiler/linker support; failures do not silently
disable it. Use `make LTO=0` for a non-LTO build, cleaning before switching
settings. `make print-config` reports the selected LTO setting and flags.
Benchmark runs force a rebuild and inherit `LTO` from the environment; prefix
the runner command with `LTO=0` for a non-LTO comparison.

## Running Benchmarks and Finding Results

Run from the repository root using Python 3.8+ and a checkout with a valid Git
commit. Use the runner so only `simulate(...)` is measured; do not profile the
ordinary application's loading and output by wrapping `build/app` in perf.
State/workspace allocation, parsing, reset, validation, and output stay outside
the measured call; `simulate` performs no allocation or conversion.

```sh
# Timing: macOS or Linux; mean and sample variance over 30 runs.
python3 bench/run.py timing tests/data/*.txt --dt 0.0001 --steps 10 --samples 30

# Profiling: Linux perf CPU call stacks and function report.
python3 bench/run.py profile tests/data/1000.txt --dt 0.0001 --steps 100

# General hardware counters, including cache references/misses: Linux.
python3 bench/run.py counters tests/data/1000.txt --dt 0.0001 --steps 100

# Cache-only collection: Linux with these events available.
python3 bench/run.py counters tests/data/1000.txt --dt 0.0001 --steps 100 \
  --events cache-references:u,cache-misses:u

# All three modes, independently, across multiple inputs: Linux.
python3 bench/run.py all tests/data/*.txt --dt 0.0001 --steps 100 --samples 30
```

`--dt` and `--steps` are required. Timing defaults to 10 samples (minimum 2);
`--samples` is available only for `timing` and `all`. Profiling/counters each
measure one call. Every mode defaults to `--warmups 1`; state resets before each
call. Use `--results-dir PATH` to change the parent directory.

The printed result path is `results/<short-hash>_<UTC-timestamp>_<type>/`.
Each run retains `bench` and combined `build.stdout.log`. Each input gets a
numbered directory, such as `0001-1000/`, containing `input.txt` and mode folders:

| Mode folder | Main output files |
| --- | --- |
| `timing/` | `samples.csv`, `summary.csv`; aggregate `timing-summary.csv` at the run root |
| `profile/` | `report.txt`, `perf.data`, `sample-events.txt` |
| `counters/` | `perf-stat.txt`, `counters.csv`, `ratios.csv` |

Read each mode's `status.json` and error logs before reporting success. Setup or
build failures produce root `error.txt`. Preserve prior runs; do not recreate
source snapshots, manifests, hash files, or command journals. The directory hash
identifies HEAD; the saved binary also includes any uncommitted code changes.

Profiling/counters require Linux perf and permission to access the selected
events; a VM/container may not expose hardware counters. On macOS, `all` saves
timing and failure statuses for the other modes. Do not bypass platform checks
or claim mocked/skipped tests are real measurements. Increase `--steps` if a
profile has no samples; use supported `--events` when counters are unavailable.
See [README.md](README.md#benchmarking-and-profiling) for the complete layout.

## Coding Style & Naming Conventions

Use four-space indentation, function opening braces on a separate line, and
control-flow opening braces on the same line. Use `snake_case` for functions and
variables, `PascalCase` for typedefs, and uppercase header guards. Keep private
helpers `static` and read-only callback parameters `const`.

No formatter or linter is configured. Compile cleanly with the existing
`-std=c17 -Wall -Wextra -Wpedantic` flags. Prefer straightforward loops and small
modules over additional abstractions.

## Architecture Constraints

Keep gravity separate from the simulation loop. `Interaction` returns one
source body's acceleration contribution without mutating state. Preserve forward
Euler ordering: input accelerations drive the first step, and newly computed
accelerations are assigned only after all pairs are evaluated. Preserve the
instance format: body count followed by `w x y z vx vy vz ax ay az` rows.
`Simulation` owns thirteen contiguous arrays, each of length `count`, in one
allocation. Construct it with `simulation_create` and release `storage` with
`free`. Initialize the ten state arrays before calling `simulate`; its three
workspace arrays are overwritten. Validate input and allocation sizes at the
boundary; `simulate` trusts the storage, timestep, and callback contracts.

## Testing Guidelines

Use standard C `assert` and name test functions `test_<behavior>`. Register new
tests in the existing test `main`; application sources except `src/main.c` are
linked automatically. Add regression tests for numerical or interface changes
and run `make test`. No coverage percentage is configured. Check update order,
callback consistency, malformed input, and saved-state round-tripping when relevant.
For actual Linux profiler/counter validation, run `NBODY_TEST_PERF=1 make test`
on a host with perf access. This opt-in test uses temporary files under `build/`
and removes them afterward; use the benchmark runner for persistent results.
Ordinary `make test` skips this real-perf test unless explicitly enabled.

## Generating Additional Test Cases

Use `generate_instance.py` when more random inputs are needed; it requires only
Python 3's standard library. Run from the repository root:

```sh
mkdir -p build
python3 generate_instance.py --n 1000 --L 10 \
  --w-min 1 --w-max 100 --v-max 3 --a-max 2 \
  --seed 42 --output build/instance-seed42.txt
```

Adjust the body count, bounds, and seed for the case being investigated. Keep
generated files in ignored `build/`; place only intentional regression fixtures
in `tests/data/`. Record the arguments and seed to reproduce failures, and reuse
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

Git history contains only `init`, so no established message convention exists.
Use concise, imperative subjects describing one focused change. PRs should
explain the resulting behavior, relevant numerical or interface implications,
and commands used for validation. Link related issues when applicable and update
README examples when usage changes.
