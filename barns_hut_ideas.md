This document proposes performance experiments for the serial Barnes–Hut baseline
described below. First, investigate target traversal order and node metadata.
Then use profiles to choose between arithmetic and tree construction experiments.
These proposals are hypotheses based on the source. This document contains no
measurements of additional speedups.

The proposals concern the earlier per-particle traversal baseline of the solver
in [src/approximate.c](src/approximate.c). That baseline rebuilds an octree each
step and aggregates mass and center of mass. It then recursively traverses the
tree for each target. Defaults are theta 0.5, bucket capacity 8, maximum depth 32,
and a pool of `max(1, 4 * count)` nodes. The baseline already uses reusable
workspace, separate state arrays, `-O3`, and link-time optimization (LTO).

Preserve the serial C17 interface and Euler update order. Preserve body/output
order. Keep the simulation call free of allocations and conversions.

First, measure time in `root_bounds`, `subdivide`, `aggregate`, and `traverse`,
including leaf gravity evaluations. The baseline diagnostics count accepted
aggregates and fallbacks, but do not identify the dominant cost. An optional
diagnostic build could count node visits, leaf interactions, maximum depth, and
leaf occupancy. It could also time the major phases. Exclude instrumentation
from the final timing comparison. Diagnostic totals include warmups and samples.

The experiment order reflects implementation cost and likely usefulness.
Use profiles to decide which experiments to pursue.

1. **Visit targets in the tree's existing spatial order.**

   In the baseline, `evaluate` visits targets in original input order.
   `workspace->indices` already groups bodies by octant. Try visiting targets as
   `target = indices[p]`. Write accelerations to `ax[target]`, `ay[target]`, and
   `az[target]`. Nearby targets may reuse recently accessed nodes and source data.
   This requires no additional sort or state conversion.

   Pass the target position and rank `p` into the traversal once to avoid
   repeated lookups. Keep the rank array initially because tree validation also
   uses it. The tree, opening decisions, and summation order within each target
   can remain unchanged.

   Compare against the baseline at the same theta, including shuffled input
   orders. Scattered output writes may offset some of the cache benefit.
   Measure the whole simulation.

2. **Make leaf detection and child iteration cheaper.**

   `is_leaf` scans up to eight child slots. Traversal scans them again before
   descending. Store a child count or occupancy mask to reduce leaf detection
   to one check. The builder already reserves occupied siblings consecutively.
   Test a representation with `first_child` and `child_count`.

   This representation could also shrink `OctreeNode`. The baseline layout
   stores eight `size_t` child indices plus geometry, ranges, and aggregates.
   Measure actual `sizeof(OctreeNode)` and cache behavior before further layout
   changes. Another experiment could separate cube centers, which construction
   alone uses, from traversal fields. This separation adds array accesses.

   Preserve octant order. Retain an occupancy mask if construction needs the
   original octant mapping. Preserve the parent-before-child index invariant
   that reverse-order aggregation uses. Update structural tests.

3. **Reduce opening-test arithmetic with a guarded fast path.**

   In the baseline, every eligible internal-node visit evaluates two nested
   `hypot` calls before testing `s / d < theta`. For ordinary coordinate
   magnitudes, test whether a squared comparison can reject a node before a
   square root. For positive theta and distance, the algebraically equivalent
   condition is
   `s*s < theta*theta*d2`, where `d2 = dx*dx + dy*dy + dz*dz`.

   Treat this expression as a starting point. Its intermediate products can
   overflow or underflow even when the baseline comparison is meaningful.
   Use proven range guards or scaled arithmetic. Retain the baseline `hypot`
   path for exceptional magnitudes and ambiguous threshold cases.

   Preserve these numerical conditions:

   - Keep the strict inequality.
   - Exclude the target.
   - Preserve the handling of zero and nonfinite distances.
   - When theta is zero, open the node immediately without the acceptance calculation.

   Separately investigate reciprocal distance reuse in `point_mass` for accepted
   nodes. Multiplication order changes rounding and representable range.
   Check these numerical effects separately. Benchmark whether guards cost less
   than the arithmetic they save. Include extreme coordinates, tiny and large
   theta, and values just on either side of the opening boundary.

4. **Reduce repeated work during tree construction.**

   `subdivide` classifies each body once to count octants and again to scatter
   indices. It then copies the partition back. Cache an octant byte during the
   counting pass in workspace allocated before simulation. Reuse that byte
   during scattering. The extra memory traffic may outweigh three coordinate
   comparisons. Test this change independently.

   If copying dominates, investigate alternating the existing index and scratch
   buffers across levels. Define which buffer owns each completed range.
   Ensure that leaves and rank construction consistently use the final
   permutation. Mixed-depth leaves make this change more involved than removing
   one `memcpy`.

   A smaller change could remove the initial `rank[i] = i` write. Successful
   construction overwrites all ranks before traversal. A bounds failure selects
   the direct solver immediately. Preserve diagnostics for empty calls.

5. **Tune leaf size and the small-system crossover.**

   Compare bucket capacities 4, 8, 16, and 32. Larger leaves reduce tree work but
   increase direct interactions. The best balance depends on body distribution
   and theta. These settings require `ApproximateConfig` changes in C.
   The command-line interface exposes only theta.

   Record effective settings for every experiment. Keep runner metadata
   consistent with the executable.

   For small systems, test dispatching to `gravity_direct` before building the
   tree. Its symmetric pair loop evaluates each unordered pair once, whereas
   leaf traversal evaluates both target directions. Start with systems that
   already fit in one leaf. Use measurements to choose any broader crossover
   threshold. Add explicit dispatch diagnostics. Update interface tests.

   Bucket size, depth, and pool capacity can change which groups are available
   for approximation. Compare their accuracy as well as runtime. Increase pool
   capacity only when diagnostics show capacity-limited leaves. Larger unused
   pools alone do not reduce traversal work.

6. **Test iterative traversal if recursion or branches dominate.**

   Replace recursive calls with an explicit depth-first stack. Push children in
   reverse order to preserve visitation order. An eight-way tree can require
   up to `1 + 7 * max_depth` pending node entries. One entry per level is
   insufficient. Allocate scratch before simulation or use a bounded local
   array that supports the maximum depth of 64.

   Keep target coordinates, rank, and acceleration sums local throughout the
   walk. Compare instructions, branch misses, and total runtime. A manual stack
   also adds loads and stores. Do not assume that removing recursion improves
   performance.

7. **Investigate leaf arithmetic and SIMD after checking compiler output.**

   Leaves call `gravity` for each source, with indirect source loads and three
   running sums. Inspect whether LTO already inlines the call and hoists target
   data. If necessary, factor a shared internal kernel that makes those values
   explicit. Preserve the baseline behavior for coincident positions.

   Clang's `-Rpass=loop-vectorize`, `-Rpass-missed=loop-vectorize`, and
   `-Rpass-analysis=loop-vectorize` explain vectorization decisions. Floating-point
   reductions and gathers can limit profitable vectorization. Support varies
   by target architecture. Keep fast-math disabled. Inspect the generated code
   before proposing intrinsics. See the primary [LLVM vectorization documentation](https://llvm.org/docs/Vectorizers.html).

   A later experiment could batch independent targets to share source loads,
   but each target still needs its own opening decision and self exclusion.
   Retain scalar handling for divergent traversal and small leaves. Any manual
   reassociation of sums requires accuracy checks.

8. **Consider tree reuse with containment checks only if rebuilding dominates.**

   For slowly moving bodies, try retaining the previous tree while every body
   remains within its assigned leaf cube. At each step:

   - Check containment.
   - Rebuild immediately when containment fails.
   - Recompute aggregates from current positions.

   Preserve the baseline fallback behavior for bounds and aggregates.

   Reusing topology can produce different grouping from rebuilding the root
   cube every step, even when containment remains valid. Measure accuracy and
   traversal degradation as well as saved build time. After each state reset
   between benchmark samples, invalidate cached topology or independently check
   all memberships before workspace reuse. Skipping force
   evaluations or reusing stale centers of mass would change the simulation.

Theta tuning is a separate accuracy/performance experiment. Sweep values such
as 0.25, 0.5, 0.75, and 1.0 only against an agreed error budget. Larger theta
usually accepts more aggregates, but provides no percentage-error guarantee.
Report a time-versus-error comparison separately from implementation speedups
at fixed settings. Even at fixed theta, changed tree topology can change error.
Multithreading and GPU execution would require revisiting the project's serial
scope and should follow the serial experiments above.

Use the existing runner for reproducible comparisons. The following commands
describe proposed experiments. The author did not run them for this document:

```sh
python3 bench/run.py timing tests/data/1000.txt \
  --dt 0.0001 --steps 100 --samples 30 --warmups 1 \
  --solver approximate --theta 0.5 --compare-direct

# Linux with perf access. Increase steps if the profile has too few samples.
python3 bench/run.py profile tests/data/1000.txt \
  --dt 0.0001 --steps 100 --solver approximate --theta 0.5

python3 bench/run.py counters tests/data/1000.txt \
  --dt 0.0001 --steps 100 --solver approximate --theta 0.5 \
  --events cycles:u,instructions:u,branches:u,branch-misses:u,cache-references:u,cache-misses:u
```

Use supported events on the measurement host. Portable timing works on macOS.
The runner's Linux profiling and counter modes do not work on macOS.
Check each mode's `status.json` and error logs before interpreting results. Also check
`reference/status.json` and `error-status.json` when comparing against direct.

Reuse identical saved inputs across variants. Include small cases, larger clouds
generated with `generate_instance.py`, several seeds, clustered distributions,
and nearly coincident bodies. Keep generated inputs in ignored `build/` and
persistent benchmark results in `results/`. Record generator arguments and
seeds. Random initial accelerations are valid arbitrary states. Use accelerations
computed from direct gravity for physically consistent cases.

Keep compiler, LTO setting, theta, timestep, steps, samples, and warmups fixed
within each comparison. The runner rebuilds and retains its executable.
Preserve earlier result directories. Uncommitted variants can share the same
HEAD hash.

Measure the complete `simulate_with_solver` call, including tree construction.
Keep raw samples as whole-run seconds. Summarize them in ms/step, with variation
across whole-run averages. Alternate baseline and candidate runs when machine
load or thermal drift could dominate.

For each implemented change, run `make test`. Add relevant numerical or
interface regression tests. Preserve the following behavior:

- Theta zero agrees with direct gravity.
- The solver excludes the target from its own force calculation.
- Subdivision remains bounded.
- The solver retains its handling of degenerate inputs.
- The solver retains its fallback behavior.

Check accelerations on identical snapshots. Also check final trajectories.
`--compare-direct` currently measures trajectory differences, not isolated force
error on identical positions. Retain an optimization only when the complete
simulation improves reproducibly at the required accuracy. Report memory use
and important regressions alongside runtime.
