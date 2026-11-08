# Rushil solvers

Build and run from the repository root on this Mac (Homebrew libomp):

```sh
make -C rushil
./rushil/naive 40 20000
OMP_NUM_THREADS=4 ./rushil/barnes_hut 40 20000 0.6
```

Arguments are body count and step count; Barnes-Hut also accepts theta.
Both drivers use seed 42, `init_bodies`, timestep 0.002, and kick-drift-kick
leapfrog. Both use unsoftened Newtonian gravity with G = 1 and skip coincident
positions. Close encounters may require a smaller timestep. Smaller theta
improves tree accuracy; 0.000001 effectively disables grouping on ordinary inputs.

`naive.c` calls `n2_accelerations` directly and performs leapfrog updates.
The unused Euler integrator and solver registry/descriptor headers have been
removed from this folder. `physics.c` and `physics.h` now contain only the shared
body type, initialization, and CM diagnostics. The original `src/` is unchanged.

## Files and memory layout

- `n2.c`, `n2_internal.h`: direct acceleration kernel and declaration.
- `nbody.h`: state views and allocation API.
- `state.c`, `state_internal.h`: one allocation containing ten contiguous arrays
  (mass, three positions, three velocities, three accelerations).
- `naive.c`, `barnes_hut_nbody.c`: runnable leapfrog drivers.
- `physics.c`, `physics.h`: shared initialization and diagnostics.
- `test_state.c`: shared initialization, CM diagnostics, and direct-force checks.

The direct solver uses a structure of arrays: x values are adjacent in memory,
then y, z, etc. Its force loop streams positions and masses without loading
interleaved velocity/other body fields. `restrict` pointers describe disjoint
input/output storage to the compiler. Target coordinates and accumulated target
accelerations stay in local scalar variables. The first target initializes the
outputs, avoiding a separate clearing pass. Each pair's distance is computed
once and used for both bodies. These optimizations are present in the code;
the isolated runtime benefit of the layout has not been benchmarked.
There is no explicit SIMD implementation or special over-alignment.

Both drivers now use Simulation arrays throughout. Initialization and CM
measurements operate directly on views, with no Body struct or conversion.
Barnes-Hut borrows coordinate pointer tables from those same arrays; nodes
store body indices. The private Simulation layout stays in state_internal.h.

Run local regression tests with `make -C rushil test`.
On Linux, override the OpenMP options:
`make -C rushil CC=gcc OMPFLAGS=-fopenmp OMPLIBS=-fopenmp`.

## Compare the two Barnes-Hut implementations

```sh
python3 rushil/compare_performance.py
# Customize the workload:
python3 rushil/compare_performance.py --bodies 1000 4000 --steps 50 \
  --dt 0.0001 --theta 0.6 --samples 5 --warmups 1 --threads 1 4
```

The script builds separate benchmark adapters for `src/barnes_hut.c` and
`rushil/barnes_hut_nbody.c`, using the same compiler flags (no LTO/fast-math),
saved input, and leapfrog integration loop. This deliberately replaces the
original src Euler loop for the comparison; production sources are unchanged.
The src force backend is serial and is run with one thread. The rushil backend
is measured with each requested thread count. Integration scheduling is the same
in both adapters, using each run's configured thread count.

Timing includes initial acceleration evaluation, integration, and all tree
construction/reuse, summarization, and force work. It excludes build, input,
allocation, state reset, validation, and output. Each run is reset to the same
saved input; rushil topology is invalidated before each run. Prepared workspace
is reused across samples. `rushil/bench/compare_driver.c` includes the standalone
rushil source with its main renamed to exercise its actual static force routines.

Runs save binaries, build logs, input, raw whole-run samples, final states,
`comparison.csv`, and `status.json` under a new `results/barnes-hut-src-vs-rushil-*`
directory. Check status before using results. Summaries show mean/stddev ms per
step, src-relative speedup, and final-state errors against serial direct n2 using
the same leapfrog loop. Same theta does not imply equal accuracy between trees.
Inputs use equal total mass 1, uniform positions in [-1,1], zero velocity, and
Python RNG seed 42. Default results describe only this distribution and short
simulation duration; they are not a universal performance ranking.

### Sweep the summarize task cutoff

```sh
python3 rushil/compare_performance.py --theta 0.6 --threads 4 \
  --bodies 1000 4000 16000 --steps 20 --samples 5 --warmups 1 \
  --base-cases 0 64 256 512 2000 8000 32768 2147483647
```

Both implementations receive the same theta. Each cutoff is compiled into a
separate binary using `-DSUMMARIZE_BASE_CASE=VALUE`; the production default
remains 2000. At/below the cutoff, summarize still recurses but creates no new
tasks. Zero tasks every internal node; 2147483647 effectively disables task
creation. Cutoff/thread configurations must produce exactly identical final
states or the benchmark fails. The cutoff appears in comparison.csv, binary
names, and raw sample filenames. Variant order is shuffled reproducibly.

The November 8 sweep and ten-sample confirmation found no universal optimum:
64/256/2000 were close at 4000 bodies, while 256 led the confirmation at 16000.
Use 256 as a candidate for larger four-thread runs, then measure your workload.
Cutoffs above body count all use serial summarization, so timing differences
between those configurations reflect noise/code-generation effects, not more
parallel work. This sweep measures complete timestep performance, not isolated
summarization time.

To try a cutoff in the standalone executable (force rebuild when changing flags):

```sh
make -B -C rushil CFLAGS='-O3 -std=c17 -Wall -Wextra -Wpedantic -DSUMMARIZE_BASE_CASE=256'
```

The comparison's direct reference now uses the local `rushil/naive.c`
acceleration wrapper and `rushil/n2.c`, rather than `src/n2.c`.
The common benchmark loop supplies matched integration without CM printing;
`naive.c`'s normal main remains the entry point for `./rushil/naive`.
`comparison.csv` includes `n2_speedup` relative to the serial direct reference.
