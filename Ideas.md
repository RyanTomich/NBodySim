## 4. Native array layout (implemented)

[`Simulation`](src/nbody.h) now stores mass, position, velocity, acceleration,
and acceleration workspace in separate contiguous arrays. Gravity reads `w`,
`x`, `y`, and `z` directly, without stepping over velocity and acceleration
fields. Input is parsed into this layout; allocation and reset happen outside
the measured `simulate()` call, and there is no packing step.

Layout alone leaves the number of square roots, divisions, and dependent
accumulations unchanged. A small timing improvement does not establish which
of those costs dominates. The next opportunity is to use the contiguous arrays
to reuse data across interactions and calculate several interactions together.
The following items are proposed experiments, not implemented optimizations.

The current [`gravity_accelerations()`](src/gravity.c) evaluates only `i < j`
and updates both bodies. Preserve that advantage when testing blocking; compare
it explicitly against any SIMD design that computes both ordered interactions.

**## 5. Block the interaction loop for cache reuse

Blocking, or tiling, processes a small group of bodies against another group
before moving on. A source block can remain in cache while several targets use
it, reducing repeated transfers from larger caches or memory. This is the
general purpose of [loop blocking](https://www.intel.com/content/www/us/en/developer/articles/technical/loop-optimizations-where-blocks-are-required.html).

For this kernel, divide the body indices into blocks of `B` and visit the upper
triangle of the block grid:

```text
zero next_ax, next_ay, next_az
for target block I in increasing order:
    for source block J from I through the last block:
        for i in I in increasing order:
            load i's accumulated acceleration into local sums
            for j in J in increasing order, with j > i:
                compute displacement and inv_r3; skip coincident positions
                add (mass[j] * inv_r3) * displacement to i's local sums
                subtract (mass[i] * inv_r3) * displacement from j's workspace
            store i's local sums
```

An off-diagonal tile evaluates all `B * B` pairs between its blocks. A diagonal
tile evaluates only `i < j`, avoiding self-interactions and duplicate pairs.
Use smaller ranges at the end when `count` is not divisible by `B`.

The arrays make both blocks compact: each body needs four input doubles and
three workspace doubles during this kernel. Two full blocks therefore occupy
about `112 * B` bytes of relevant data, before cache-line overhead and other
working data. Try `B = 32, 64, 128, 256` as measurement candidates, not fixed
defaults. Too-small blocks add loop overhead; too-large blocks lose cache reuse.
The workspace for `j` still receives an update for each pair, but those updates
now revisit a small region.

The traversal above can preserve each body's increasing source-index order:
earlier blocks supply its lower-index contributions first, then its own block,
then later blocks. Keep the current pair arithmetic and add each contribution
directly to the running sum. Computing a separate subtotal per tile and adding
that subtotal later would change floating-point grouping.
**
Start with scalar blocking alone to measure its effect. Then use a much smaller
group of targets inside a tile, retaining several sets of coordinates and
accumulators in BBregisters. This creates independent calculations for SIMD and
can reuse a loaded source across those targets. Cache block size and the number
of targets kept in registers are separate tuning choices.

## 6. Use SIMD to compute several interactions at once

SIMD instructions apply one operation to several numbers, called lanes. With
separate arrays, adjacent `x` values can be loaded together, as can `y`, `z`,
and mass; the kernel need not extract fields from interleaved body structs.
The choice of what each lane represents determines the accumulation problem.

**Keep pair sharing and vectorize across sources.** For one target `i`, lanes
represent consecutive sources `j`. Broadcast `i`'s coordinates and mass, load
vectors of source coordinates and masses, and calculate several displacements,
squared distances, square roots, and scales together. Updates to the distinct
`next_a[j]` elements can be stored as vectors. Contributions to `i`, however,
must be combined into one acceleration. A horizontal reduction combines lanes
and typically changes addition order. An ordered reduction can preserve order
but may be slower; support depends on target and compiler.
[LLVM documents these reduction constraints](https://llvm.org/docs/Vectorizers.html#reductions).

**Alternatively, vectorize across independent targets.** Give each lane its own
target and three acceleration sums. Load target coordinates once and broadcast
one source at a time:

```text
lanes = consecutive targets i, i+1, ...
ax, ay, az = zero vectors
for source j = 0 through count-1:
    broadcast w[j], x[j], y[j], z[j]
    calculate one contribution per target lane
    add each contribution to that lane's own sums
store the completed target accelerations
```

Here the lanes never need to be reduced together, so each target can retain its
original source order. This is an ordered-pair gravity kernel: it computes both
`(i, j)` and `(j, i)`, roughly doubling the distance calculations relative to the
current kernel. It could still win through better SIMD throughput and less
workspace traffic, but that must be measured. Keeping pair sharing in this
version introduces a reduction for the shared source's acceleration instead;
the two-sided update does not become independent just by changing lane direction.

Handle self-pairs and coincident positions with masks. For inactive lanes, use
a safe denominator before division and select a zero contribution afterward;
computing infinity and multiplying by zero can produce NaN. Handle the final
partial vector with a scalar tail or safe masked loads, never an out-of-range
read. Keep the existing unsoftened gravity law.

Keep fast-math disabled initially. Matching source order alone does not ensure
bitwise equality: fused multiply-add, reciprocal approximations, and changed
pair arithmetic can also alter rounding. Keep contraction settings consistent
between variants; assess approximate reciprocal square roots as a separate
numerical tradeoff. [Clang's floating-point options](https://clang.llvm.org/docs/UsersManual.html#controlling-floating-point-behavior)
describe these controls.

First try compiler vectorization on a small, directly expressed gravity loop;
use architecture-specific intrinsics only if diagnostics and measurements
justify them. Inspect the force loop itself, not just initialization loops.
Use `-Rpass=loop-vectorize`, `-Rpass-missed=loop-vectorize`, and
`-Rpass-analysis=loop-vectorize`, then inspect generated instructions and time
the result. A successful vectorization remark does not prove a speedup.
[LLVM vectorization diagnostics](https://llvm.org/docs/Vectorizers.html#diagnostics).

## Other opportunities enabled by the layout

- **Expose non-overlap to the compiler.** Construction already guarantees that
  state and workspace slices do not overlap. Express the relevant contracts
  with C `restrict` on the actual array pointers used by the gravity kernel,
  especially its writable outputs. This may eliminate alias checks and allow
  better load scheduling. `const` alone does not promise non-overlap, and
  restricting a pointer to the container does not establish independence of
  every array it contains. These must be true contracts, not new runtime checks.
  See [restricted pointers](https://gcc.gnu.org/onlinedocs/gcc/Restricted-Pointers.html).
- **Vectorize the Euler update.** Contiguous position, velocity, and acceleration
  arrays also allow several bodies to be advanced together. Preserve old velocity
  for the position update and old acceleration for the velocity update. This
  removes no O(n²) work, so expect it to matter mainly for smaller instances.
- **Swap acceleration buffers instead of copying.** Swap the three current and
  next pointers after all interactions, saving the final O(n) copy. This needs
  changes to the current fixed-offset assumptions: refresh the callback view
  after swaps and make benchmark reset/validation follow the active state arrays
  instead of assuming state is always the first ten slices. Free the original
  `storage` allocation as before. Treat this as lower priority than the pair loop.
- **Parallelize independent target blocks later.** An ordered-pair kernel can
  assign disjoint output ranges to workers while sharing read-only source arrays.
  The current pair-sharing kernel writes both ends of each pair, so simply
  parallelizing its outer loop would race. It needs conflict-free tile scheduling
  or per-worker acceleration buffers and a later reduction, with added memory,
  synchronization, and possible rounding changes. Complete all interactions
  before publishing accelerations or starting the next Euler step.

## How to compare the experiments

Measure scalar blocking first, then each SIMD strategy, using the same saved
inputs and compiler/LTO settings. Try sizes both within and beyond cache capacity,
including counts that are not multiples of a tile or vector width. Keep setup
and reset outside timing; include all per-step work inside `simulate()`. Use
`bench/run.py` with repeated samples, and inspect each run's status and error
logs. Linux cache counters can help test the cache-reuse hypothesis where available.

Check results against the direct reference and run `make test` after code changes.
Exercise unequal masses, coincident positions, extreme separations, zero and one
body, multiple steps, and saved-state resumption. Preserve custom callback
behavior and Euler ordering. Compare runtime for the same simulation work;
report separately whether the kernel evaluates ordered or unordered pairs.
