# 3D N-body simulator

A serial C17 simulator for comparing direct (`n2`), Barnes–Hut, and fast multipole
(FMM) gravity solvers. All backends use the same forward Euler simulation loop
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
resume a simulation. **Always use a different output filename from the input.**

Builds enable link-time optimization by default. If your compiler does not
support it, run `make clean && make LTO=0`. Clean before changing compiler
settings; `make print-config` shows the active configuration.

## Benchmarking

The runner rebuilds the executables, resets the input before every run, and times
only the simulation call. Parsing, setup, resets, and output are excluded.

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

## Run tests

```sh
make test

# Rebuild everything from scratch and run all tests.
make clean && make && make test
```

This runs the C state and solver tests, command-line tests, and benchmark-driver
tests. `make clean` removes `build/` and preserves saved benchmark results.

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
are position, velocity, and acceleration components. **Input accelerations drive
the first step.** Generated random accelerations represent arbitrary states;
physically consistent initial accelerations must be calculated from gravity.

## Architecture and project map

The command-line app loads a state, prepares the selected gravity backend, and
runs a shared forward Euler loop. Each step advances positions and velocities,
then asks the backend for new accelerations. State ownership, integration, and
gravity calculations are separate, so changing backends preserves the simulation
workflow. Barnes–Hut and FMM share an octree for grouping nearby bodies.

Gravity uses normalized units (`G = 1`). Forward Euler is a simple first-order
integrator; there is no adaptive timestep, softening, or collision handling.

| Location | Purpose |
| --- | --- |
| `src/main.c`, `src/io.*` | Command-line entry point and instance I/O |
| `src/nbody.h`, `src/state.c`, `src/euler.c` | Public simulation interface, state ownership, and time stepping |
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
