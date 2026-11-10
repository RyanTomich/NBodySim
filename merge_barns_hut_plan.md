# Barnes–Hut merge plan

Use `src/` as the implementation backbone. Adapt Rushil's parallel force
evaluation, conditional tree rebuilding, and fast kick-drift-kick integration
to the existing public solver/state interfaces. This document is a plan; it does
not implement the merge or claim a speedup.

Validation is limited to correctness checks and short runnable examples. No
benchmark campaign, profiling, parameter sweep, or extensive new test framework
is required. All callers, including the benchmark driver, will use leapfrog through
the shared simulation entry point. Limit benchmark-related edits to compatibility
and identifying the new integrator; preserve existing results unchanged.

## 1. Keep and replace decisions

References below name existing files and exact symbols. New symbols are marked
as proposed; function names are the primary references as line numbers change.

| Component | Retain from `src/` | Adapt from `rushil/`; replacement |
| --- | --- | --- |
| Solver lifecycle | `solver_create`, `solver_destroy`, `AccelerationSolver`, and the private `BarnesHutWorkspace` | Keep all ownership in the workspace; do not import Rushil's file-scope simulation/tree globals. |
| Spatial tree | `Octree`, `OctreeNode`, bounded node pool, bucket leaves, and index partitioning | Replace unconditional rebuilding with validation followed by rebuilding only when needed. Do not import `Node`, `insert`, or the single-body leaf representation. |
| Source discovery | `compute_target_bounds`, `accept_source`, and `collect_sources` | Keep the group acceptance rule and child traversal order. Replace shared source-list storage with one list per worker. Do not replace this with Rushil's per-body `walk`. |
| Force evaluation | `evaluate_group`, `gravity_accumulate_sources`, and `gravity_monopole_acceleration` | Replace the serial target-leaf loop with OpenMP dynamic scheduling over target leaves. Each leaf remains one indivisible unit of work. |
| Mass summaries | `GravityMonopole`, `gravity_monopole_add`, and `aggregate_tree` | Keep the serial bottom-up pass. Do not port task-based `summarize` or the unused center-of-mass velocity `vcom`. |
| Integration | The `simulate_with_solver` signature, opaque `Simulation`, thirteen-array storage, and separate acceleration output buffers | Replace the Euler implementation with fast kick-drift-kick leapfrog based on Rushil's `step`. Recompute initial acceleration; remove Euler rather than keeping an integrator selector. |
| Reliability | Finite-value checks, bounded subdivision, direct fallback, and diagnostics | Keep these safeguards instead of Rushil's exit-on-capacity-exhaustion behavior. |

Source map:

- [src/barnes_hut.c](src/barnes_hut.c): workspace, aggregation, source discovery,
  evaluation, backend options, and diagnostics.
- [src/octree_internal.h](src/octree_internal.h) and
  [src/octree.c](src/octree.c): tree representation, `octree_create`,
  `octree_build`, `compute_root_bounds`, and `subdivide`.
- [src/solver.h](src/solver.h), [src/nbody.h](src/nbody.h), and
  [src/state_internal.h](src/state_internal.h): public contracts and private state.
- [rushil/barnes_hut_nbody.c](rushil/barnes_hut_nbody.c): `accelerations`
  contains OpenMP dynamic scheduling; `tree_is_valid` implements reuse checks;
  `summarize` illustrates the aggregation dependency; `step` implements leapfrog.

## 2. Phase A: isolate worker storage, then parallelize target leaves

First refactor `src/barnes_hut.c` without changing numerical behavior.

1. Add a proposed private `BarnesHutWorker` containing a `size_t *sources`
   buffer and the three counters `source_node_visits`, `direct_evaluations`, and
   `aggregate_evaluations`. Each source list needs capacity `count`, matching
   the existing disjoint-source-range bound in `collect_sources`.
2. Allocate workers and their source buffers in `barnes_hut_create`; free them
   in `barnes_hut_destroy`, including cleanup after partial allocation failure.
   Check size multiplications before allocation. This uses O(worker_count × n)
   additional scratch memory. Keep `tree.scratch` for serial tree partitioning.
3. Change `collect_sources` to read a const workspace and append to an explicit
   worker source list. Change `evaluate_group` to accept that list and the
   worker's counters. Neither function may update shared workspace counters.
   Preserve source ordering and each target's accumulation arithmetic.
4. Replace `evaluate`'s serial node loop with an OpenMP work-sharing loop using
   `schedule(dynamic, 1)`. Iterations for internal nodes do nothing; each leaf
   iteration computes current target bounds, collects sources, then evaluates
   every body in that leaf. No additional leaf-index array is necessary.
5. Each thread selects its own worker using `omp_get_thread_num()`. The tree,
   monopoles, and input arrays are read-only during this phase. Leaves cover
   disjoint target ranges, so acceleration outputs need neither locks nor
   floating-point reductions.
6. Zero worker counters before the loop, then combine their integer counts after
   its completion barrier. Leave cumulative totals updates in
   `gravity_barnes_hut` after `evaluate` returns. Keep direct fallback outside
   the parallel region.

For the first implementation, capture `omp_get_max_threads()` at solver creation
and reserve that many worker slots. Use an explicit `num_threads(worker_count)`
on the region so later runtime settings cannot exceed allocated capacity.
Smaller actual teams are safe. In a serial build, reserve one worker and use the
same evaluation helpers. Document that `OMP_NUM_THREADS` must be set before
solver creation and that a prepared solver is not concurrently callable.

Update the root [Makefile](Makefile) with proposed `OPENMP=0|1`, defaulting to 0,
and overridable `OMPFLAGS`/`OMPLIBS`. With `OPENMP=1`, apply flags to every target
that links the solver; guard OpenMP headers and calls with `_OPENMP`. Model the
platform overrides on [rushil/Makefile](rushil/Makefile), while retaining the
root compiler, LTO, assertion, and floating-point settings. An explicitly
requested OpenMP build should fail clearly if support is unavailable.

## 3. Phase B: reuse topology, refresh all numerical summaries

Adapt `tree_is_valid` to bucket leaves; checking a single representative body
would be incorrect.

1. Add a proposed private `octree_is_valid(state, tree)` declaration and
   implementation in `src/octree_internal.h` and `src/octree.c`. For a nonempty
   input, require an existing successfully built tree with matching count.
   For each leaf, iterate every `indices[p]` in `[begin,end)` and check all three
   finite coordinates against that leaf's cube. Reject uncertain or
   unrepresentable containment calculations rather than allowing NaN comparisons
   to pass. Preserve the builder's containment invariants; do not add a tolerance
   that accepts bodies outside a source cell. Handle an empty input explicitly.
2. In `BarnesHutWorkspace`, track whether the previous topology is usable. Add
   a backend `reuse_tree` option accepting `0` or `1`, defaulting to `1`, through
   `barnes_hut_default_config` and `barnes_hut_create`. A value of `0` retains
   rebuild-every-call behavior for focused correctness comparisons. This option
   initially needs only the existing `SolverOption` API, not a new CLI flag.
3. Replace `evaluate`'s unconditional `octree_build` with: validate when reuse
   is enabled and topology is usable; otherwise rebuild. Invalidate the cached
   topology after a failed build. Preserve the existing `n2_accelerations`
   fallback for failed bounds or aggregation; never use stale monopoles.
4. Run `aggregate_tree` on every successful evaluation, even when topology is
   reused. Recompute `TargetBounds` for every target leaf as well. Positions and
   masses can change between calls, and callers may supply a different state
   allocation of the same count. Cache no borrowed state pointers or source lists.
5. Keep reuse policy in Barnes–Hut. FMM continues calling `octree_build` as
   before; do not assume its `occupied_half_width` or expansion data can be reused.
   Keep current diagnostic meanings explicit: `capacity_limited_leaves` may
   describe the retained topology per evaluation, while fallback flags are
   reset on each call. Report `reuse_tree` through `barnes_hut_settings`.

Only topology and membership are reused. A valid older tree can yield a different
approximation from a newly rebuilt tree at nonzero theta; exact equality between
those policies is not a requirement. Self-exclusion still uses disjoint index
ranges, and leaves still evaluate their sources directly.

## 4. Phase C: replace Euler with fast leapfrog

Adapt the arithmetic in Rushil's `step` while preserving `src`'s ownership and
separate input/output acceleration buffers. Leapfrog becomes the only production
integrator for `n2`, Barnes–Hut, and FMM. This deliberately changes numerical
results and the initial-acceleration contract.

### Entry point and initialization

- Rename [src/euler.c](src/euler.c) to proposed `src/leapfrog.c` and replace its
  implementation of `simulate_with_solver`. Keep that public name and signature;
  do not add a second entry point or `--integrator` option. The Makefile discovers
  the renamed source automatically; ensure only one definition remains.
- Update the contract in `src/nbody.h`: positive-step calls recompute acceleration
  from the current masses and positions, regardless of the input acceleration
  values. Compute into `next_ax/next_ay/next_az` and publish only after the entire
  callback returns. Return immediately for zero steps or zero bodies, leaving all
  state untouched and invoking no solver.
- Initialize once per nonempty, positive-step invocation, then reuse each step's
  new acceleration for the next step. A call with `steps` steps therefore makes
  exactly `steps + 1` force evaluations. Reinitializing on a later invocation is
  intentional: public mutable views allow callers to change positions and masses
  between calls, so cached accelerations cannot be assumed valid.

### Fast timestep implementation

Precompute `half_dt = 0.5 * dt`. For each step:

1. In one body loop, fuse the first half-kick and drift: for each coordinate,
   compute `v_half = v + half_dt * a`, store `v_half`, then advance
   `x += dt * v_half`.
2. After every position has been updated, acquire a fresh `SimulationView` and
   invoke the selected solver once into the alternate acceleration buffers.
   Publish those buffers only after the callback has completed.
3. In one body loop, finish the kick with `v += half_dt * a_new`. Returned
   velocities are full-step velocities and returned accelerations correspond
   to the final positions; do not expose a staggered-velocity state to callers.

Use small private helpers with disjoint `restrict` array parameters, following
`advance_state` in the current Euler source. Keep scalar temporaries in the fused
loop, reuse the thirteen existing arrays, and swap acceleration pointers rather
than copying arrays. No per-step allocation, conversion, or extra force call is
needed. Keep integration loops serial initially to avoid extra OpenMP team
startup and synchronization; the force phase already distributes target leaves.
Retain the existing compiler optimization settings and disabled fast-math.
These choices avoid unnecessary work but do not constitute a measured speedup.

Document that solver acceleration must depend on current masses and positions,
not velocities or previous accelerations. During the force evaluation after a
drift, the visible velocities are half-step values. All current gravity backends
satisfy this restriction; replace the old acceleration-dependent test callback
instead of claiming that arbitrary callbacks retain their previous semantics.
Reacquire borrowed views after buffer publication and keep any swap helper private.

### Callers, saved state, and documentation

- [src/main.c](src/main.c) continues calling `simulate_with_solver` with the same
  CLI. Its behavior now uses leapfrog automatically for every backend.
- [bench/driver.c](bench/driver.c)'s warmup, sample, and direct-reference calls
  automatically use the same leapfrog implementation. Keep initialization inside
  the measured simulation call, so its one extra force evaluation per invocation
  is included. Keep resets outside that call; do not move initial force work into
  an untimed preparation phase. Normalize whole-call time by `steps` as before.
- Add an `integrator: "leapfrog"` field to newly saved timing/reference status
  records in [bench/run.py](bench/run.py), so new outputs can be distinguished
  from historical Euler runs. Do not rewrite old artifacts or add timing studies.
- Preserve the ten-field instance format, including saved acceleration fields.
  Those fields round-trip through I/O but are recomputed when integration starts.
  Restart correctness should use deterministic direct forces; a reused
  approximate tree and a newly prepared tree can produce different approximations.
- Update [README.md](README.md), [AGENTS.md](AGENTS.md), and the architecture
  documentation in [docs/architecture.mmd](docs/architecture.mmd) and its SVG to
  describe leapfrog, the renamed source, and the new initialization contract.
  Remove stale instructions to preserve Euler ordering and add OpenMP examples.
  Explain second-order integration for smooth forces without promising exact
  energy conservation or symplectic behavior for the approximate tree solver.

## 5. Minimal correctness and run verification

Reuse the existing public-interface helpers in
[tests/test_solver.c](tests/test_solver.c), especially `fill_cloud`, `check_force`,
`test_gravity_accuracy`, and `test_registered_solver_contract`. Extend a few
behavioral cases rather than creating a comprehensive test matrix.

1. **Forces and parallel execution:** use the existing small deterministic cloud
   to exercise multiple leaves. At theta zero, compare with the independent
   scalar reference using the existing `1e-11` normalized tolerance. At nonzero
   theta, retain the existing fixture-specific error budget. Run the same checks
   with one and four OpenMP threads. For identical topology and inputs, thread
   scheduling must preserve each target's result; compare the outputs directly.
2. **Reuse and invalidation:** evaluate a prepared solver repeatedly with small
   in-cell movements and changed masses, then with a large position change that
   requires rebuilding. Compare theta-zero forces with the scalar reference;
   include a nonzero-theta force check so stale aggregates cannot hide behind
   direct evaluation. Use public views and options, not private tree assertions.
   Existing tests already cover coincident positions and capacity/depth limits.
3. **Integration:** update `check_step`, `test_integration_contract`, and the
   callback used by `test_io_contract` in [tests/test_state.c](tests/test_state.c)
   for position-dependent leapfrog forces. Check one and a few steps against an
   independent scalar kick-drift-kick calculation, including whole-system output
   publication, zero-step/empty-state behavior, full-step final velocities, and
   the initial evaluation plus one force evaluation per step. Split invocations
   each have their own initial evaluation; account for that in callback counts.
   Verify that changing only supplied initial accelerations does not change a
   positive-step result, and retain a direct-force saved-state restart check.
4. **CLI and driver compatibility:** update the one-body expected result in
   [tests/test_cli.sh](tests/test_cli.sh). With no other masses, input acceleration
   is replaced by zero: at `dt=0.25`, the existing fixture must return
   `1 2 3.25 4.5 4 5 6 0 0 0` as its body row. Retain zero-step round-tripping
   and direct-force restart checks. Adjust the reset explanation in
   [tests/test_bench.sh](tests/test_bench.sh) to the recomputation contract;
   preserve its short app-versus-driver final-state checks. If status metadata
   is changed, check that new status objects identify leapfrog without launching
   a benchmark campaign.

Run the existing `make test` target once in a clean serial build and once in a
clean OpenMP build; this includes the existing benchmark-driver regression
checks, not a timing study. Run the OpenMP unit executable once more with the
other thread count. Do not add benchmarks, sanitizer matrices, long trajectories,
or exhaustive parameter sweeps.

Example commands after implementation on this Mac, using Rushil's existing
Homebrew OpenMP paths:

```sh
make clean && make && make test
make clean
OMP_NUM_THREADS=4 make OPENMP=1 \
    OMPFLAGS='-Xclang -fopenmp -isystem /opt/homebrew/opt/libomp/include' \
    OMPLIBS='-L/opt/homebrew/opt/libomp/lib -lomp' test
OMP_NUM_THREADS=1 ./build/tests
OMP_NUM_THREADS=4 ./build/app --solver barnes-hut --theta 0 \
    tests/data/two_body.txt 0.0001 3 > build/merged-leapfrog.txt
./build/app --solver n2 tests/data/two_body.txt 0.0001 3 \
    > build/reference-leapfrog.txt
```

Inspect both saved outputs for the expected body count, ten finite values per
body, and successful exit. The cloud checks above exercise parallel tree work;
the two-body runs establish that the merged solver and direct reference both run
with leapfrog. Compare their final states within floating-point tolerance, not
against historical Euler trajectories.

The merge is complete when both builds compile cleanly, these correctness checks
pass, and both short runs produce valid states. Parallel aggregation, parallel
tree construction, task-cutoff tuning, and performance ranking are outside this
plan. Leave `rushil/` intact as reference code.
