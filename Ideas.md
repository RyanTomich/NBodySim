# $N^2$ Optimizations

These are performance experiments for the direct all-pairs implementation, not
implemented changes or measured speedups. Start with the cost of each interaction
before restructuring the whole simulator.

## Current implementation

Each step in [`simulate()`](src/nbody.c) updates all positions and velocities,
evaluates all `n * (n - 1)` ordered interactions, and then assigns the new
accelerations. Work is O(n²) per step and memory is O(n). The acceleration
workspace is allocated once per `simulate()` call, not once per step.

Inspection with Apple Clang 17 on ARM64 using the project's normal compiler
flags showed:

- An indirect callback call for every interaction.
- Two actual `hypot` library calls per [`gravity()`](src/gravity.c) evaluation.
- No loop vectorization of the inner pair loop; compiler diagnostics cite the
  callback and control flow.
- Loads and stores of the acceleration accumulator after every interaction.

These observations identify opportunities; profiling is still needed to establish
their relative runtime cost.

## 1. Make gravity visible to the optimizer

The callback currently prevents gravity from being inlined into the pair loop.
Inlining could also let the compiler retain the target's coordinates in registers
and optimize several interactions together.

First experiment with link-time optimization (`-flto`), which enables optimization
across source files. Inspect the generated code to verify whether it actually
removes the indirect call. If insufficient, consider a gravity-specific kernel
alongside the generic callback implementation. Keep gravity separate from time
integration and preserve the generic `Interaction` behavior.

Reference: [LLVM link-time optimization](https://llvm.org/docs/LinkTimeOptimization.html).

## 2. Simplify distance arithmetic

For displacement `(dx, dy, dz)`, consider computing:

```text
r2 = dx*dx + dy*dy + dz*dz
scale = source_mass / (r2 * sqrt(r2))
```

This offers one square root and straightforward arithmetic instead of two nested
`hypot` calls. It is a high-priority experiment, but not a numerically equivalent
replacement for every input.

Directly squaring components can overflow or underflow, and results will round
differently. Establish acceptable coordinate ranges and error tolerances, retain
the zero contribution for coincident positions, and test extreme separations.
Adding gravitational softening would change the physical model and should be
considered separately.

## 3. Accumulate locally and store once per target

Use local `ax`, `ay`, and `az` sums throughout the inner loop, then write the
completed vector to the acceleration workspace. The current assembly accesses
accumulator memory for every pair; local sums give the compiler a clearer
opportunity to retain them in registers. Verify the generated code and runtime
rather than assuming the compiler will eliminate all spills.

This experiment can preserve source visitation and addition order. Keep the
separate workspace: the [callback tests](tests/test_main.c) require every callback
to observe a consistent body state, including the previous accelerations, until
all pairs have been evaluated.

## 4. Improve data layout, then add blocking and SIMD

Each [`Body`](src/nbody.h) contains ten doubles, while gravity reads only mass and
position: four doubles. Separate arrays for `x`, `y`, `z`, and mass would provide
compact, contiguous inputs for a gravity-specific kernel. The public instance
format can remain unchanged. Include any packing or conversion cost in timing.

Next, process small blocks of targets against blocks of sources, reusing source
data across several targets. This can improve cache reuse and expose independent
calculations for SIMD. Tune block sizes through measurement.

Vectorizing across sources often changes floating-point reduction order.
Vectorizing across independent targets can retain each target's source order.
Keep fast-math disabled for the baseline comparison and explicitly evaluate any
proposal to relax floating-point behavior.

Use compiler vectorization reports to verify what happens instead of assuming
`-O3` vectorizes the loop. Useful Clang diagnostics are
`-Rpass=loop-vectorize`, `-Rpass-missed=loop-vectorize`, and
`-Rpass-analysis=loop-vectorize`.

Reference: [LLVM vectorization documentation](https://llvm.org/docs/Vectorizers.html).

## 5. Evaluate each unordered pair once

For gravity, the displacement and distance for `(i, j)` can serve both bodies.
For each pair with `i < j`, the mathematical updates are:

```text
q = (position[j] - position[i]) / distance(i, j)^3
acceleration[i] += mass[j] * q
acceleration[j] -= mass[i] * q
```

This halves the number of distance calculations, but does not guarantee a 2x
speedup. Updating two accumulators complicates memory access, vectorization, and
threading; compare it against an optimized ordered-pair kernel.

Accelerations are not equal and opposite when masses differ. This optimization
applies to gravity, not arbitrary `Interaction` callbacks. Check rounding changes
and accumulate into separate workspace before assigning body accelerations.

## Measurement and correctness

- Measure one change at a time with `bench/run.py`, using several body counts and
  the same saved inputs. Use repeated samples and inspect each mode's
  `status.json` and error logs.
- Track runtime and runtime divided by `n * (n - 1) * steps` for `n > 1`. For the
  symmetric variant, retain that denominator for comparison and separately
  distinguish its actual unordered-pair count.
- Include all work inside `simulate()`, including allocation, packing, and
  conversion. The existing benchmark excludes input loading, state reset, and
  final output from its measured interval.
- Timing works on macOS; the repository's profiling and hardware-counter modes
  require Linux with usable perf access. Use profiles and counters to determine
  whether math, instruction overhead, or memory traffic dominates.
- Compare numerical results against the baseline on small cases and run
  `make test` after implementation changes. Exercise unequal masses, coincident
  positions, extreme separations, custom callbacks, and multiple steps.
- Preserve forward Euler ordering: input accelerations drive the first step,
  positions use the previous velocities, and new accelerations are assigned only
  after all pair evaluations. Compute final accelerations because they are part
  of the saved state.

Allocation reuse, the self-interaction branch, and the final O(n) acceleration
copy are lower priorities for large n. Begin with compiler visibility, local
accumulation, and distance arithmetic, then measure whether layout changes,
blocking, pair symmetry, or parallel execution provide further gains.
