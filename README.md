# 3D N-body simulator

Simulate gravity in three dimensions and compare five solvers on the same initial conditions.

## Run a simulation

You need Python 3.8+, Make, and a C17 toolchain with OpenMP support.
Run this from the root of a Git checkout to build the simulator and time a two-body orbit.

```sh
python3 bench/run.py timing tests/data/two_body.txt \
  --solver n2 --dt 0.01 --steps 100 --samples 3
```

- **Initial conditions.** Use `tests/data/two_body.txt` for the orbit, `bench/data/1000.txt` for a 1,000-body system, or your own file in the final-state format below.
- **Duration.** The example advances 100 steps of 0.01 time units each. Each of the three timing samples starts from the same input.
- **Integrator.** Leapfrog is the default. Use `--integrator euler` for forward Euler, which uses the input accelerations for its first velocity update.

| `--solver` | How it computes gravity |
| --- | --- |
| `n2` (default) | Evaluates every pair directly. Use it as the numerical reference. |
| `serial-barnes-hut` | Approximates distant groups of bodies using a serial tree walk. |
| `parallel-barnes-hut` | Runs independent tree walks for individual bodies in parallel. |
| `grouped-barnes-hut` | Shares one tree walk among nearby target bodies and processes groups in parallel. |
| `fmm` | Uses fourth-order fast multipole expansions for distant interactions. |

For Barnes–Hut, `--theta` defaults to `0.5`. Smaller values generally improve accuracy by evaluating more interactions in detail.

## Read the results

Open the directory printed after `Results:` to find the files from the two-body run.

| Result | Where to find it |
| --- | --- |
| Timing | The terminal reports the mean ± sample standard deviation in milliseconds per step. `timing-summary.csv` records both in seconds per step. |
| Final state | `0001-two_body/timing/final.txt` contains the state after step 100. The first line is the body count. Each remaining row lists mass, position, velocity, and acceleration in the order `w x y z vx vy vz ax ay az`. |

## Contribute to a solver

Edit the gravity algorithm in [src/](src/). To add a solver, declare its descriptor in [solver.h](src/solver.h) and register it in [solver.c](src/solver.c).

Add a behavioral regression case in [test_solver.c](tests/test_solver.c) and run `make test`. For an approximate solver, add `--compare-n2` to the timing command to measure final-state differences. To compare runtime, time `n2` separately with the same input, integrator, timestep, and step count. See the [repository guidelines](AGENTS.md) for contribution conventions.
