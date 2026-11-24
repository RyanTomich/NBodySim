# Initial distributions and Bash runner

This folder collects **actual copies of all 88 saved initial states** from our
previous experiments: 14 families, including both galaxy initializations.
`run_all.sh` runs the three production implementations on identical inputs.
The C executable `build/bench` contains all three and selects them with its
`n2`, `barnes-hut`, and `fmm` arguments; separate solver executables are not needed.
The script runs that binary directly, without rebuilding between cases.

## Inputs

| `inputs/` subfolder | Body counts | Seeds / notes |
|---|---|---|
| `pbbs_cube`, `pbbs_shell`, `pbbs_plummer`, `pbbs_ball`, `pbbs_kuzmin` | 256, 1024, 4096, 10000 | Native PBBS distributions; equal masses, total mass 1 |
| `galaxy_merger` | 256, 1024, 4096, 10000 | Original GalaxyMerger-derived, pre-relaxed two-disk setup |
| `spiral_merger` | 4096 | Revised live bulge/disk/halo model with seeded spiral arms |
| `synthetic_uniform`, `synthetic_clusters`, `synthetic_core_halo`, `synthetic_disk`, `synthetic_filament`, `synthetic_close_pairs`, `synthetic_unequal_masses` | 256, 2048, 10000 | Seeds 42, 43, 44 for each size; normalized cold starts |

Files are named `n10000.txt` or `n10000_seed42.txt`. Their format is the repository's
body count followed by `mass x y z vx vy vz ax ay az`. Each has a sibling
`.source.json` describing its origin. They are physical copies, so the runner
does not depend on the old results folders once preparation is complete.
Generated `inputs/` files are ignored by Git; copy the whole folder when moving
the inputs to another machine.

To recollect them from the saved experiments, from NBodySim:

```bash
python3 experiments/distribution_suite/prepare.py
```

Preparation validates counts, finite state values, positive masses and unique
positions. It preserves bytes and refuses to replace a locally changed input.
It reads the two completed synthetic studies, the repository-distribution
study, and the revised Barnes–Hut movie study; paths can be overridden with
`prepare.py --help`. The earlier experiment generators remain available if the
source result folders need to be regenerated.

## Run

From NBodySim, use an existing `build/bench`, or build once:

```bash
make clean
make

# Inspect the 88 inputs without running any binaries.
./experiments/distribution_suite/run_all.sh --list

# All inputs, all three solvers, one thread: 30 samples, 3 warmups, 10 steps.
./experiments/distribution_suite/run_all.sh

# Quick check of one input.
SAMPLES=3 WARMUPS=1 STEPS=2 ./experiments/distribution_suite/run_all.sh \
  experiments/distribution_suite/inputs/pbbs_cube/n256.txt

# Only the native PBBS inputs.
./experiments/distribution_suite/run_all.sh \
  experiments/distribution_suite/inputs/pbbs_*/*.txt
```

For a thread sweep, build OpenMP first (clean before changing build flags):

```bash
make clean
# GCC/LLVM on Linux:
make OPENMP=1 CPPFLAGS=-DNBODY_BENCH_OPENMP

# Apple Clang with this machine's Homebrew libomp, instead:
make OPENMP=1 CPPFLAGS=-DNBODY_BENCH_OPENMP \
  OMPFLAGS='-Xclang -fopenmp -isystem /opt/homebrew/opt/libomp/include' \
  OMPLIBS='-L/opt/homebrew/opt/libomp/lib -lomp'

THREADS="1 2 4 8" SAMPLES=15 ./experiments/distribution_suite/run_all.sh
```

`NBODY_BENCH_OPENMP` adds the untimed team-formation check, without per-step
instrumentation. The runner rejects instrumented binaries and rejects a serial
binary when more than one thread is requested. Thread budgets apply equally to
each solver. FMM's force team can be smaller, since it is limited by occupied
root children; construction and aggregation also have serial phases.

Configuration is through environment variables:

| Variable | Default |
|---|---|
| `BENCH_BIN` | `NBodySim/build/bench`; can point to another compatible benchmark executable |
| `THREADS` | `1`; whitespace-separated budgets, e.g. `"1 2 4 8"` |
| `DT`, `STEPS` | `1e-5`, `10` |
| `SAMPLES`, `WARMUPS` | `30`, `3` |
| `THETA`, `ORDER` | `0.3`, `4` |
| `OUTPUT_DIR` | A new timestamped directory under `results/distribution-suite/` |

## Outputs and interpretation

Every run archives the executable and build-info record, copies each input once,
and writes `summary.tsv`. Each input/solver/thread case retains `samples.csv`,
`final.txt`, `stderr.log`, and `status.txt`. The root status is `ok` only after
every case succeeds. Existing output directories are refused; failed runs retain
their logs and stop immediately.

The summary records mean and sample SD in **ms/step**, speedup over n2 at the
same thread budget, and final position/velocity RMS differences from that n2
run. Timings measure only the production leapfrog simulation, including its
initial force evaluation; parsing, allocation, reset and output are excluded.
Differences are unweighted RMS values over particles at matching dt and duration.
They compare against same-timestep n2, not an exactly converged physical solution.
For identical-position force errors and timestep-refined references, use the
full runners in [initial_setups](../initial_setups/README.md) and
[repository_distributions](../repository_distributions/README.md).

**This three-solver sweep uses unsoftened G=1 throughout.** The spiral-merger
input preserves the new movie's masses, positions and velocities, but its use
here adapts it to the common unsoftened law. The softened Barnes–Hut movie is a
separate experiment; FMM does not implement that softened law. This Bash runner
saves final states, rather than animation snapshots.
