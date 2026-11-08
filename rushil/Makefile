CC = clang
CFLAGS = -O3 -std=c17 -Wall -Wextra -Wpedantic
# macOS Homebrew libomp defaults. Override OMPFLAGS/OMPLIBS on other systems.
OMPFLAGS = -Xclang -fopenmp -isystem /opt/homebrew/opt/libomp/include
OMPLIBS = -L/opt/homebrew/opt/libomp/lib -lomp

.PHONY: all test
all: naive barnes_hut

naive: naive.c n2.c state.c physics.c nbody.h n2_internal.h state_internal.h physics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) naive.c n2.c state.c physics.c $(OMPLIBS) -lm -o $@

barnes_hut: barnes_hut_nbody.c physics.c state.c physics.h nbody.h state_internal.h
	$(CC) $(CFLAGS) $(OMPFLAGS) barnes_hut_nbody.c physics.c state.c $(OMPLIBS) -lm -o $@

test_state: test_state.c physics.c n2.c state.c physics.h n2_internal.h nbody.h state_internal.h
	$(CC) $(CFLAGS) test_state.c physics.c n2.c state.c -lm -o $@

test: all test_state
	./test_state
