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
LDLIBS += -lm
PYTHON ?= python3

SOURCES := $(wildcard src/*.c)
LIB_SOURCES := $(filter-out src/main.c,$(SOURCES))
TEST_SOURCES := $(wildcard tests/*.c) $(LIB_SOURCES)
BENCH_SOURCES := $(wildcard bench/*.c) $(LIB_SOURCES)
HEADERS := $(wildcard src/*.h tests/*.h bench/*.h)

.PHONY: all test clean print-config

all: build/app build/bench

build/app: $(SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) -Isrc $(SOURCES) $(LDLIBS) -o $@

build/bench: $(BENCH_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) -Isrc $(BENCH_SOURCES) $(LDLIBS) -o $@

build/tests: $(TEST_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) -Isrc $(TEST_SOURCES) $(LDLIBS) -o $@

build/bench-test: bench/driver.c bench/metadata.c tests/bench/fake_simulation.c $(filter-out src/nbody.c,$(LIB_SOURCES)) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(CFLAGS) -Isrc -Ibench bench/driver.c bench/metadata.c tests/bench/fake_simulation.c $(filter-out src/nbody.c,$(LIB_SOURCES)) $(LDLIBS) -o $@

build:
	mkdir -p $@

test: all build/tests build/bench-test
	./build/tests
	sh tests/test_cli.sh
	$(PYTHON) -m unittest discover -s tests/bench -v

print-config:
	@printf '%s\n' 'CC=$(CC)' 'LTO=$(LTO)' 'LTO_FLAGS=$(LTO_FLAGS)' \
	    'CPPFLAGS=$(CPPFLAGS)' 'CFLAGS=$(CFLAGS)' 'LDFLAGS=$(LDFLAGS)' 'LDLIBS=$(LDLIBS)'

clean:
	rm -rf build
