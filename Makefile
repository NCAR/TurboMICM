# tmc: a C11 port of MICM's CPU solver.
#
#   make env        create the conda environment in .conda (Python for build and test only)
#   make cpu        build the library and bin/ts1_bench_cpu
#   make bench      cpu, plus bin/ts1_scenario (Tier 3)
#   make test-cpu   run the test suite: one line per MICM case, "<id> PASS|FAIL|N/A <reason>"
#   make micm-ref   build MICM's reference TS1 driver with clang++ and OPTFLAGS
#   make proto      the library as a shared object for the Python prototypes (prototypes/)
#   make micm-scenario  MICM's scenario driver, bin/micm_scenario, for small references
#
# The compiler comes from CC and the optimization flags from OPTFLAGS; the
# Makefile adds only language, warning, include, OpenMP, and link flags.
.DEFAULT_GOAL := cpu

ifeq ($(origin CC),default)
CC := clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
ROOT ?= ..
MICM ?= $(ROOT)/micm
REFERENCE ?= $(ROOT)/reference
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON := $(CURDIR)/.conda/bin/python

# Tuning parameters (see DESIGN.md). Every value gives correct results.
WIDTH ?= 2
GROUP_SIZE ?= 32
KERNEL_CHUNK ?= 16
KERNEL_MODE ?= fast
# the TS1 kernels' species ordering: minimum fill, or MICM's Markowitz in exact mode
ifeq ($(KERNEL_MODE),exact)
TS1_ORDERING ?= markowitz
else
TS1_ORDERING ?= minimum-fill
endif
GEN_FLAGS ?=

UNAME_S := $(shell uname -s)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -n 1)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
endif
ifneq (,$(findstring clang,$(CC_VERSION)))
ifeq ($(UNAME_S),Darwin)
OMP_CFLAGS := -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include
OMP_LIBS := -L/opt/homebrew/opt/libomp/lib -lomp
else
OMP_CFLAGS := -fopenmp
OMP_LIBS := -fopenmp
endif
else
OMP_CFLAGS := -fopenmp
OMP_LIBS := -fopenmp
endif

CFLAGS := -std=c11 -Wall -Wextra -Wpedantic $(OPTFLAGS) $(OMP_CFLAGS) -Iinclude -DTMC_WIDTH=$(WIDTH) -DTMC_DEFAULT_GROUP_SIZE=$(GROUP_SIZE)
LDLIBS := $(OMP_LIBS) -lm

LIB_SOURCES := $(wildcard src/*.c)
LIB_OBJECTS := $(patsubst src/%.c,build/%.o,$(LIB_SOURCES))

.PHONY: env cpu bench test-cpu micm-ref clean

# The environment is created once, when its Python is missing; after editing
# environment.yml, `make env-recreate` builds it again.
env: | $(PYTHON)
$(PYTHON):
	CONDA_REGISTER_ENVS=false CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes
env-recreate:
	rm -rf .conda
	$(MAKE) env
.PHONY: env-recreate

cpu: bin/ts1_bench_cpu
bench: cpu bin/ts1_scenario

# Rebuild everything when the compiler or flags change: the objects depend on
# build/flags.txt, and when it changes every object is removed as well, so that
# an object compiled in the same second as the change cannot pass as current.
build/flags.txt: FORCE
	@mkdir -p build
	@printf '%s\n' "$(CC)" "$(CFLAGS)" "$(LDLIBS)" > build/flags.tmp
	@if cmp -s build/flags.tmp $@; then rm build/flags.tmp; else mv build/flags.tmp $@; rm -f build/*.o build/proto/*.o build/*.a build/*.dylib build/*.so generated/*.o generated/tests/*.o; fi
.PHONY: FORCE
FORCE:

build/%.o: src/%.c $(wildcard include/*.h) src/internal.h build/flags.txt
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

build/libtmc.a: $(LIB_OBJECTS)
	$(AR) rcs $@ $^

# ---- generated mechanisms ----------------------------------------------
configs/ts1.json: tools/translate_ts1.py $(MICM)/test/benchmark/ts1_mechanism.hpp | $(PYTHON)
	$(PYTHON) tools/translate_ts1.py $(MICM)/test/benchmark/ts1_mechanism.hpp $@

generated/ts1.c generated/ts1.h: configs/ts1.json tools/generate.py tools/mechanism.py | $(PYTHON)
	@mkdir -p generated
	$(PYTHON) tools/generate.py configs/ts1.json generated/ts1 --name ts1 --width $(WIDTH) --chunk $(KERNEL_CHUNK) --mode $(KERNEL_MODE) --ordering $(TS1_ORDERING) $(GEN_FLAGS)

build/ts1.o: generated/ts1.c generated/ts1.h include/tmc.h build/flags.txt
	$(CC) $(CFLAGS) -Igenerated -c $< -o $@

bin/ts1_bench_cpu: bench/ts1_bench.c bench/csv_table.h build/ts1.o build/libtmc.a generated/ts1.h $(wildcard include/*.h)
	@mkdir -p bin
	$(CC) $(CFLAGS) -Igenerated bench/ts1_bench.c build/ts1.o build/libtmc.a $(LDLIBS) -o $@

bin/ts1_scenario: bench/ts1_scenario.c bench/csv_table.h build/ts1.o build/libtmc.a generated/ts1.h $(wildcard include/*.h)
	@mkdir -p bin
	$(CC) $(CFLAGS) -Igenerated bench/ts1_scenario.c build/ts1.o build/libtmc.a $(LDLIBS) -o $@

# ---- tests --------------------------------------------------------------
TEST_SOURCES := $(wildcard tests/*.c)
EXAMPLES := robertson rate_constants_no_user_defined rate_constants_user_defined
TEST_MECHANISMS := $(patsubst %,generated/%.o,$(EXAMPLES))

test-cpu: cpu bin/test_cpu
	TMC_TEST_MICM=$(MICM) TMC_TEST_REFERENCE=$(REFERENCE) bin/test_cpu $(TEST_FILTER)
	@TMC_TEST_MICM=$(MICM) TMC_TEST_REFERENCE=$(REFERENCE) bin/test_cpu --extra 1>&2
	@$(MAKE) --no-print-directory test-methods 1>&2

# Tier 3: each C method against its Python prototype (states to roundoff, equal
# operation counts) and the transformed tableaus against their standard form,
# on the mini scenario; printed to stderr by test-cpu, so the 680 lines stay the list.
test-methods: bench proto prototypes/mini/reference_snapshots.csv
	$(PYTHON) prototypes/tableau_check.py
	$(PYTHON) prototypes/agreement.py
prototypes/mini/reference_snapshots.csv: prototypes/make_mini.py bin/micm_scenario | $(PYTHON)
	$(PYTHON) prototypes/make_mini.py
.PHONY: test-methods

generated/tests/registry.h generated/tests/case_list.inc: $(TEST_SOURCES) $(REFERENCE)/micm_cases.txt tools/test_registry.py | $(PYTHON)
	@mkdir -p generated/tests
	$(PYTHON) tools/test_registry.py $(REFERENCE)/micm_cases.txt tests generated/tests

# the example configurations; their tests compare them with the same mechanisms built in code
.SECONDEXPANSION:
.SECONDARY: # keep the generated example sources between builds
generated/%.c generated/%.h: $$(wildcard configs/examples/%/*) tools/generate.py tools/mechanism.py | $(PYTHON)
	@mkdir -p generated
	$(PYTHON) tools/generate.py $(firstword $(wildcard configs/examples/$*/config.yaml configs/examples/$*/config.json)) generated/$* --name $* --width $(WIDTH) --chunk $(KERNEL_CHUNK) --mode $(KERNEL_MODE)

generated/%.o: generated/%.c generated/%.h include/tmc.h build/flags.txt
	$(CC) $(CFLAGS) -Igenerated -c $< -o $@

bin/test_cpu: $(TEST_SOURCES) $(wildcard tests/*.h include/*.h) generated/tests/registry.h generated/tests/case_list.inc build/ts1.o $(TEST_MECHANISMS) build/libtmc.a
	@mkdir -p bin
	$(CC) $(CFLAGS) -Itests -Igenerated -Igenerated/tests -Isrc $(TEST_SOURCES) build/ts1.o $(TEST_MECHANISMS) build/libtmc.a $(LDLIBS) -o $@

# ---- prototypes: the library as a shared object for ctypes (development tool) ----
# Built from its own objects with -fPIC, so the benchmark's objects are untouched.
ifeq ($(UNAME_S),Darwin)
SHARED_FLAGS := -dynamiclib
PROTO_LIB := build/libtmc_proto.dylib
else
SHARED_FLAGS := -shared
PROTO_LIB := build/libtmc_proto.so
endif
PROTO_OBJECTS := $(patsubst src/%.c,build/proto/%.o,$(LIB_SOURCES)) build/proto/ts1.o build/proto/shim.o
proto: $(PROTO_LIB)
build/proto/%.o: src/%.c $(wildcard include/*.h) src/internal.h build/flags.txt
	@mkdir -p build/proto
	$(CC) $(CFLAGS) -fPIC -c $< -o $@
build/proto/ts1.o: generated/ts1.c generated/ts1.h include/tmc.h build/flags.txt
	@mkdir -p build/proto
	$(CC) $(CFLAGS) -fPIC -Igenerated -c $< -o $@
build/proto/shim.o: prototypes/shim.c generated/ts1.h $(wildcard include/*.h) src/internal.h build/flags.txt
	@mkdir -p build/proto
	$(CC) $(CFLAGS) -fPIC -Igenerated -Isrc -c $< -o $@
$(PROTO_LIB): $(PROTO_OBJECTS)
	$(CC) $(SHARED_FLAGS) $(PROTO_OBJECTS) $(LDLIBS) -o $@
.PHONY: proto

# ---- MICM reference -----------------------------------------------------
micm-ref: bin/micm_ts1_ref
micm-scenario: bin/micm_scenario
bin/micm_scenario: $(ROOT)/scenarios/micm_scenario.cpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(MICM)/include -I$(MICM)/test/benchmark $< -o $@
.PHONY: micm-scenario
bin/micm_ts1_ref: $(REFERENCE)/micm_ts1_ref.cpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(MICM)/include -I$(MICM)/test/benchmark $< -o $@

clean:
	rm -rf build bin generated

# ---- development tools (not part of the entry's targets) -----------------
bin/kernel_bench: bench/kernel_bench.c build/ts1.o build/libtmc.a generated/ts1.h $(wildcard include/*.h)
	@mkdir -p bin
	$(CC) $(CFLAGS) -Igenerated -Isrc bench/kernel_bench.c build/ts1.o build/libtmc.a $(LDLIBS) -o $@
