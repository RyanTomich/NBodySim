CC ?= cc
CFLAGS ?= -std=c17 -Wall -Wextra -Wpedantic

SOURCES := $(wildcard src/*.c)
TEST_SOURCES := $(wildcard tests/*.c) $(filter-out src/main.c,$(SOURCES))
HEADERS := $(wildcard src/*.h tests/*.h)

.PHONY: all test clean

all: build/app

build/app: $(SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Isrc $(LDFLAGS) $(SOURCES) $(LDLIBS) -o $@

build/tests: $(TEST_SOURCES) $(HEADERS) Makefile | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -Isrc -UNDEBUG $(LDFLAGS) $(TEST_SOURCES) $(LDLIBS) -o $@

build:
	mkdir -p $@

test: build/tests
	./build/tests

clean:
	rm -rf build
