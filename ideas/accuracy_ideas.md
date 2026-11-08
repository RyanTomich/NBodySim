# Accuracy ideas

## Leapfrog integration

Add leapfrog (kick–drift–kick) integration to reduce the long-term numerical drift of forward Euler: advance velocity by half a timestep, advance position by a full timestep, recompute acceleration, then finish the velocity update. It provides second-order accuracy with one new force evaluation per step after initialization. Fixed timesteps and conservative forces generally give bounded, oscillating energy error, although approximate solver forces can weaken this benefit. Validate by halving the timestep over the same physical horizon and checking energy and trajectory error.

## Close-particle interactions

Handle close encounters separately, since leapfrog can still miss rapidly changing forces when the timestep is too large. For smooth collective gravity, consider softening to limit forces at small separations, applying the same physical model consistently across n2, Barnes–Hut, and FMM. When encounters or binaries matter, preserve their physics with adaptive timesteps or specialized encounter integration; individual or block timesteps can reduce the cost of resolving only the difficult regions, at the expense of implementation complexity. Choose the treatment based on whether small-scale encounters are part of the intended workload.
