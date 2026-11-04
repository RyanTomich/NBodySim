The serial C17 Barnes–Hut solver resides in [src/approximate.c](src/approximate.c).
The first version (`6cdfbe1`) traversed the tree recursively for each particle.
The optimized version (`0448712`) shares traversal across each target leaf.

The Euler loop calls `gravity_approximate` through a callback to compute accelerations
for all particles. Each call rebuilds and aggregates the tree from current positions.
The solver writes separate acceleration arrays, then the integrator swaps the buffers.
Structure-of-arrays storage preserves input/output body order. The caller allocates
the workspace before simulation.

The octree stores nodes contiguously and identifies children by integer indices.
Each node owns a half-open range in a particle-index permutation. The
builder counts particles by octant and scatters their indices through scratch
storage. It creates only occupied children, after their parents.
A reverse linear pass aggregates mass and center of mass. Incremental
weighted means avoid mass–position products.

Defaults are theta 0.5, bucket capacity 8, maximum depth 32, and `max(1, 4N)`
node slots. Depth and capacity limits prevent further subdivision, as do rounding
limits. The solver retains larger leaves in these cases. It uses direct summation
when root bounds or aggregates are unrepresentable. The first version already
used these storage and construction methods, reusable allocations, `-O3`, and
link-time optimization (LTO).

The optimized version applies four changes:

1. **Share source discovery across a target leaf.** `collect_sources` traverses
   the tree once per leaf and records source leaves or accepted internal nodes.
   `evaluate_group` reuses this list for every target in the leaf, reducing
   repeated recursive calls, node loads, and opening decisions. Targets evaluate
   particles in source leaves directly and exclude themselves by index. Each
   target also evaluates monopole contributions from accepted internal nodes.
   Neighboring targets can reuse source data.

2. **Use one conservative opening decision for the group.** The traversal tests
   only internal sources whose particle ranges are disjoint from the target leaf.
   It accepts a source when `source_width / d_min < theta`. Here, `d_min` is the
   distance from the source’s center of mass to the targets’ actual bounding box.
   This criterion is geometrically conservative for every target. Discovery uses
   squared coordinate gaps and one `sqrt` instead of the old per-target nested
   `hypot` calls.

   The traversal opens nodes when squared distances are zero or nonfinite.
   It also opens nodes when the farthest L1-distance bound overflows. Theta zero
   skips the acceptance calculation. Aggregate forces still use nested `hypot`,
   and threshold comparisons retain ordinary floating-point rounding.
   Conservative grouping can open more nodes and change results at the same theta.

3. **Cache leaf status.** An `is_leaf` flag replaces repeated scans of up to
   eight child slots during aggregation, discovery, and evaluation.

4. **Reuse scratch and remove inverse ranks.** After partitioning, the existing
   `N`-entry scratch array holds one group’s source list. Its entries represent
   disjoint, nonempty particle ranges, which limit the list to `N` entries.
   Comparisons between group and source ranges replace per-particle rank lookups.
   This removes the inverse-rank allocation (`N * sizeof(size_t)` bytes) and the
   writes that initialize and construct it. The solver allocates no lists during
   simulation.

Saved comparisons measure the complete simulation call, including tree construction.
Both versions use identical inputs, theta 0.5, `dt = 0.0001`, one warmup, and LTO.
All four timing statuses report success with no tree fallbacks. Times are means
in ms/step:

| Bodies | Steps × samples | First Barnes–Hut | Optimized | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 1,000 | 100 × 30 | [1.501][first-1000] | [0.906][optimized-1000] | 1.66× |
| 20,000 | 10 × 10 | [111.348][first-20000] | [75.546][optimized-20000] | 1.47× |

These measurements establish the combined implementation gain. They do not
isolate each optimization’s contribution or establish equal numerical error.

[first-1000]: results/0448712_20261030-012721.705567Z_timing_+FourWay-FirstBarnesHut-n1000-steps100/timing-summary.csv
[optimized-1000]: results/0448712_20261030-012726.631023Z_timing_+FourWay-OptimizedBarnesHut-n1000-steps100/timing-summary.csv
[first-20000]: results/0448712_20261030-013830.865026Z_timing_+FourWay-FirstBarnesHut-n20000-steps10/timing-summary.csv
[optimized-20000]: results/0448712_20261030-013843.387511Z_timing_+FourWay-OptimizedBarnesHut-n20000-steps10/timing-summary.csv
