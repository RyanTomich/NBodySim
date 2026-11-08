#!/usr/bin/env python3
"""Generate and verify random body states using only the standard library.

Use positive property bounds for masses; signed bounds can represent charges.
Random accelerations describe arbitrary states, not necessarily physical ones.
The header contains only the body count; --seed controls reproducible generation.
"""

import argparse
import math
from pathlib import Path
import random


def uniform(rng, lower, upper):
    # A weighted sum avoids overflowing upper - lower for large finite bounds.
    fraction = rng.random()
    value = (1.0 - fraction) * lower + fraction * upper
    return min(upper, max(lower, value))


def random_vector(rng, maximum):
    magnitude = rng.uniform(0.0, maximum)
    while True:
        direction = [rng.gauss(0.0, 1.0) for _ in range(3)]
        length = math.hypot(*direction)
        if length != 0.0:
            return [magnitude * (component / length) for component in direction]


def verify_instance(args):
    positions = set()
    count = 0
    with open(args.output, encoding="utf-8") as source:
        if source.readline().split() != [str(args.n)]:
            raise ValueError("incorrect body-count header")
        for count, line in enumerate(source, start=1):
            fields = line.split()
            if count > args.n or len(fields) != 10:
                raise ValueError(f"invalid body record {count}")
            values = [float(field) for field in fields]
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"body {count} contains a nonfinite number")
            if not args.w_min <= values[0] <= args.w_max:
                raise ValueError(f"body {count} has an out-of-bounds property")
            position = tuple(values[1:4])
            if not all(-args.L <= coordinate <= args.L for coordinate in position):
                raise ValueError(f"body {count} has an out-of-bounds position")
            if position in positions:
                raise ValueError(f"body {count} has a duplicate position")
            positions.add(position)
            for name, vector, maximum in (
                ("velocity", values[4:7], args.v_max),
                ("acceleration", values[7:10], args.a_max),
            ):
                magnitude = math.hypot(*vector)
                if magnitude > maximum and not math.isclose(
                    magnitude, maximum, rel_tol=1e-12, abs_tol=0.0
                ):
                    raise ValueError(f"body {count} has an out-of-bounds {name}")
    if count != args.n:
        raise ValueError(f"expected {args.n} bodies, found {count}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, required=True, help="number of bodies")
    parser.add_argument("--L", type=float, required=True, help="position bound")
    parser.add_argument("--w-min", type=float, required=True)
    parser.add_argument("--w-max", type=float, required=True)
    parser.add_argument("--v-max", type=float, required=True)
    parser.add_argument("--a-max", type=float, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path,
                        help="output file (overwritten; default: bench/data/instance-nN-seedS.txt)")
    args = parser.parse_args()

    bounds = (args.L, args.w_min, args.w_max, args.v_max, args.a_max)
    if not all(math.isfinite(bound) for bound in bounds):
        parser.error("all bounds must be finite")
    if args.n < 0 or args.L < 0 or args.v_max < 0 or args.a_max < 0:
        parser.error("n, L, v-max, and a-max must be nonnegative")
    if args.w_min > args.w_max:
        parser.error("w-min must not exceed w-max")
    if args.n > 1 and args.L == 0:
        parser.error("L must be positive for multiple distinct positions")

    if args.output is None:
        args.output = (Path(__file__).resolve().parent / "bench" / "data"
                       / f"instance-n{args.n}-seed{args.seed}.txt")

    rng = random.Random(args.seed)
    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as destination:
            destination.write(f"{args.n}\n")
            for _ in range(args.n):
                values = [uniform(rng, args.w_min, args.w_max)]
                values.extend(uniform(rng, -args.L, args.L) for _ in range(3))
                values.extend(random_vector(rng, args.v_max))
                values.extend(random_vector(rng, args.a_max))
                destination.write(" ".join(format(value, ".17g") for value in values) + "\n")
        verify_instance(args)
    except (OSError, ValueError) as error:
        parser.exit(1, f"error: {error}\n")
    print(f"Verified {args.n} bodies in {args.output} (seed {args.seed}).")


if __name__ == "__main__":
    main()
