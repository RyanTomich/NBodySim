CC ?= cc
LTO ?= 1
ifeq ($(LTO),1)
LTO_FLAGS := -flto
else ifeq ($(LTO),0)
LTO_FLAGS :=
else
$(error LTO must be 0 or 1)
endif

CFLAGS ?= -std=c17 -Wall -Wextra -Wpedantic
override CFLAGS += -O3 -g -fno-omit-frame-pointer -fno-fast-math -UNDEBUG $(LTO_FLAGS)
# OpenMP is required by parallel-barnes-hut. Override for other installations.
ifeq ($(shell uname -s),Darwin)
OMPFLAGS ?= -Xclang -fopenmp -isystem /opt/homebrew/opt/libomp/include
OMPLIBS ?= -L/opt/homebrew/opt/libomp/lib -lomp
else
OMPFLAGS ?= -fopenmp
OMPLIBS ?= -fopenmp
endif
override CFLAGS += $(OMPFLAGS)
LDLIBS += $(OMPLIBS) -lm

SOURCES := $(wildcard src/*.c)
LIB_SOURCES := $(filter-out src/main.c,$(SOURCES))
TEST_SOURCES := $(wildcard tests/*.c) $(LIB_SOURCES)
BENCH_SOURCES := $(wildcard bench/*.c) $(LIB_SOURCES)
HEADERS := $(wildcard src/*.h src/generated/*.h tests/*.h)

.PHONY: all test clean print-config

all: build/app build/bench

build/app: $(SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) -Isrc $(SOURCES) $(LDLIBS) -o $@

build/bench: $(BENCH_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) -Isrc $(BENCH_SOURCES) $(LDLIBS) -o $@

build/tests: $(TEST_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) -Isrc $(TEST_SOURCES) $(LDLIBS) -o $@

build:
	mkdir -p $@

test: all build/tests
	./build/tests
	sh tests/test_cli.sh
	sh tests/test_bench.sh

print-config:
	@printf '%s\n' 'CC=$(CC)' 'LTO=$(LTO)' 'LTO_FLAGS=$(LTO_FLAGS)' \
	    'CPPFLAGS=$(CPPFLAGS)' 'CFLAGS=$(CFLAGS)' 'LDFLAGS=$(LDFLAGS)' 'LDLIBS=$(LDLIBS)'

clean:
	rm -rf build
