CC ?= cc
LTO ?= 1
ifeq ($(LTO),1)
LTO_FLAGS := -flto
else ifeq ($(LTO),0)
LTO_FLAGS :=
else
$(error LTO must be 0 or 1)
endif

OPENMP ?= 1
ifeq ($(shell uname -s),Darwin)
ifeq ($(shell uname -m),arm64)
BREW_PREFIX := /opt/homebrew
else
BREW_PREFIX := /usr/local
endif
ifneq ($(wildcard $(BREW_PREFIX)/opt/libomp/include/omp.h),)
OMPFLAGS ?= -Xclang -fopenmp -isystem $(BREW_PREFIX)/opt/libomp/include
OMPLIBS ?= -L$(BREW_PREFIX)/opt/libomp/lib -lomp
else ifneq ($(wildcard $(BREW_PREFIX)/opt/llvm/bin/clang),)
OMPFLAGS ?= -Xclang -fopenmp -idirafter $(shell $(BREW_PREFIX)/opt/llvm/bin/clang -print-resource-dir)/include
OMPLIBS ?= -L$(BREW_PREFIX)/opt/llvm/lib -lomp
else
OMPFLAGS ?= -Xclang -fopenmp -isystem $(BREW_PREFIX)/opt/libomp/include
OMPLIBS ?= -L$(BREW_PREFIX)/opt/libomp/lib -lomp
endif
else
OMPFLAGS ?= -fopenmp
OMPLIBS ?= -fopenmp
endif
ifeq ($(OPENMP),1)
OPENMP_FLAGS := $(OMPFLAGS) -DNBODY_REQUIRE_OPENMP
OPENMP_LIBS := $(OMPLIBS)
else ifeq ($(OPENMP),0)
OPENMP_FLAGS :=
OPENMP_LIBS :=
else
$(error OPENMP must be 0 or 1)
endif

CFLAGS ?= -std=c17 -Wall -Wextra -Wpedantic
override CFLAGS += -O3 -g -fno-omit-frame-pointer -fno-fast-math -UNDEBUG $(LTO_FLAGS)
LDLIBS += -lm

SOURCES := $(wildcard src/*.c)
LIB_SOURCES := $(filter-out src/main.c,$(SOURCES))
TEST_SOURCES := $(wildcard tests/*.c) $(LIB_SOURCES)
BENCH_SOURCES := $(wildcard bench/*.c) $(LIB_SOURCES)
HEADERS := $(wildcard src/*.h src/generated/*.h tests/*.h bench/*.h)

.PHONY: all test clean print-config

all: build/app build/bench

build/app: $(SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) $(OPENMP_FLAGS) -Isrc $(SOURCES) $(LDLIBS) $(OPENMP_LIBS) -o $@

# The REBOUND comparison uses the same timer plus the existing untimed team check.
build/bench_grouped_compare: override CPPFLAGS += -DNBODY_BENCH_OPENMP
build/bench build/bench_grouped_compare: $(BENCH_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) $(OPENMP_FLAGS) -Isrc $(BENCH_SOURCES) $(LDLIBS) $(OPENMP_LIBS) -o $@

REBOUND_SRC := third_party/rebound/src
ifeq ($(shell uname -s),Darwin)
REBOUND_LIB_ID := -Wl,-install_name,@rpath/librebound.so
REBOUND_RPATH := @loader_path
else
REBOUND_RPATH := $$ORIGIN
endif

.PHONY: rebound-benchmark
rebound-benchmark: build/bench_grouped_compare build/bench_rebound

build/librebound.so: $(REBOUND_SRC)/Makefile $(REBOUND_SRC)/Makefile.defs $(wildcard $(REBOUND_SRC)/*.[chs]) Makefile | build
	$(MAKE) -B -C $(REBOUND_SRC) OPENMP=1 QUADRUPOLE=0 SERVER=0 MPI=0 OPENGL=0 FFTW=0 \
	    CC="$(CC)" OPT="$(CPPFLAGS) $(CFLAGS) -Wno-strict-prototypes -D_GNU_SOURCE -fPIC $(OMPFLAGS)" \
	    LIB="$(LDFLAGS) $(LDLIBS) $(OMPLIBS) $(REBOUND_LIB_ID)"
	cp $(REBOUND_SRC)/librebound.so $@

build/bench_rebound: bench/comparison/rebound_driver.c src/io.c src/parse.c src/state.c $(HEADERS) build/librebound.so
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) $(OPENMP_FLAGS) -DOPENMP -Isrc -isystem $(REBOUND_SRC) \
	    $(filter %.c,$^) -Lbuild -lrebound -Wl,-rpath,'$(REBOUND_RPATH)' $(LDLIBS) $(OPENMP_LIBS) -o $@

# Serial, opt-in per-step observation; ordinary benchmark binaries are unchanged.
build/bench_steps: $(BENCH_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) -DNBODY_BENCH_STEPS -Ibench -Isrc $(BENCH_SOURCES) $(LDLIBS) -o $@

.PHONY: sweep-binary
sweep-binary: output = "$${SWEEP_BINARY:?Set SWEEP_BINARY to the output path}"
build/bench_steps_openmp: output = $@
build/bench_steps_openmp: | build
sweep-binary build/bench_steps_openmp: $(BENCH_SOURCES) $(HEADERS) Makefile
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) $(OPENMP_FLAGS) -DNBODY_BENCH_STEPS -DNBODY_BENCH_OPENMP -Ibench -Isrc $(BENCH_SOURCES) $(LDLIBS) $(OPENMP_LIBS) -o $(output)

build/tests: $(TEST_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) $(OPENMP_FLAGS) -Isrc $(TEST_SOURCES) $(LDLIBS) $(OPENMP_LIBS) -o $@

build:
	mkdir -p $@

test: all build/tests build/bench_steps
	./build/tests
	sh tests/test_cli.sh
	sh tests/test_bench.sh
	python3 tests/test_generate_instance.py
	python3 tests/test_comparison.py
	python3 tests/test_figures.py
	python3 tests/test_distribution.py
	python3 tests/test_profile_barnes_hut.py
	python3 tests/test_solver_sweep.py

ifeq ($(OPENMP),1)
test: build/bench_steps_openmp
endif

print-config:
	@printf '%s\n' 'CC=$(CC)' 'LTO=$(LTO)' 'LTO_FLAGS=$(LTO_FLAGS)' \
	    'OPENMP=$(OPENMP)' 'OPENMP_FLAGS=$(OPENMP_FLAGS)' 'OPENMP_LIBS=$(OPENMP_LIBS)' \
	    'CPPFLAGS=$(CPPFLAGS)' 'CFLAGS=$(CFLAGS)' 'LDFLAGS=$(LDFLAGS)' 'LDLIBS=$(LDLIBS)'

clean:
	rm -rf build
