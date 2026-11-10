CC ?= cc
LTO ?= 1
ifeq ($(LTO),1)
LTO_FLAGS := -flto
else ifeq ($(LTO),0)
LTO_FLAGS :=
else
$(error LTO must be 0 or 1)
endif

OPENMP ?= 0
# GCC/LLVM defaults; Apple Clang needs explicit libomp paths (see README).
OMPFLAGS ?= -fopenmp
OMPLIBS ?= -fopenmp
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
HEADERS := $(wildcard src/*.h src/generated/*.h tests/*.h)

.PHONY: all test clean print-config

all: build/app build/bench

build/app: $(SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) $(OPENMP_FLAGS) -Isrc $(SOURCES) $(LDLIBS) $(OPENMP_LIBS) -o $@

build/bench: $(BENCH_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) $(OPENMP_FLAGS) -Isrc $(BENCH_SOURCES) $(LDLIBS) $(OPENMP_LIBS) -o $@

build/tests: $(TEST_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) $(OPENMP_FLAGS) -Isrc $(TEST_SOURCES) $(LDLIBS) $(OPENMP_LIBS) -o $@

build:
	mkdir -p $@

test: all build/tests
	./build/tests
	sh tests/test_cli.sh
	sh tests/test_bench.sh
	python3 tests/test_comparison.py

print-config:
	@printf '%s\n' 'CC=$(CC)' 'LTO=$(LTO)' 'LTO_FLAGS=$(LTO_FLAGS)' \
	    'OPENMP=$(OPENMP)' 'OPENMP_FLAGS=$(OPENMP_FLAGS)' 'OPENMP_LIBS=$(OPENMP_LIBS)' \
	    'CPPFLAGS=$(CPPFLAGS)' 'CFLAGS=$(CFLAGS)' 'LDFLAGS=$(LDFLAGS)' 'LDLIBS=$(LDLIBS)'

clean:
	rm -rf build
