# C project

A minimal C17 project. Requires Make and a C compiler on a Unix-like system.

```sh
make          # Build build/app
./build/app   # Run the application
make test     # Build and run tests
make clean    # Remove build artifacts
```

## Layout

- `src/`: application source and headers; `main.c` is the entry point.
- `tests/`: tests using standard C `assert`; `test_main.c` runs the tests.
- `tests/data/`: committed test inputs, fixtures, and instances.
- `build/`: generated executables and test output, ignored by Git.

Add `.c` and `.h` files directly to `src/` or `tests/`. The Makefile picks them
up automatically. Tests link against all application sources except `src/main.c`.
Call additional test functions from `tests/test_main.c`; keep only one test
`main` function. The initial test is only a smoke check.

Run commands from the repository root. Tests can read fixtures using paths such
as `tests/data/example.txt` and write generated output under `build/`.

Override compiler settings as needed, for example `make CC=clang` or
`make CFLAGS='-std=c17 -Wall -Wextra -Wpedantic -g'`. Standard `CPPFLAGS`,
`LDFLAGS`, and `LDLIBS` overrides are also supported. Run `make clean` before
changing compiler settings. Test builds keep assertions enabled.
