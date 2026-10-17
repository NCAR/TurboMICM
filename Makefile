# tmc: a C11 port of MICM's CPU solver.
#
#   make env        create the conda environment in .conda (Python for build and test only)
#   make cpu        build the library and bin/ts1_bench_cpu
#   make bench      alias of cpu
#   make test-cpu   run the test suite: one line per MICM case, "<id> PASS|FAIL|N/A <reason>"
#   make micm-ref   build MICM's reference TS1 driver with clang++ and OPTFLAGS
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
GROUP_SIZE ?= 8
KERNEL_CHUNK ?= 16
KERNEL_MODE ?= fast

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

env: .conda/conda-meta/history
.conda/conda-meta/history: environment.yml
	CONDA_REGISTER_ENVS=false CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes
	@touch $@

cpu: bin/ts1_bench_cpu
bench: cpu

# Rebuild everything when the compiler or flags change.
build/flags.txt: FORCE
	@mkdir -p build
	@printf '%s\n' "$(CC)" "$(CFLAGS)" "$(LDLIBS)" > build/flags.tmp
	@cmp -s build/flags.tmp $@ && rm build/flags.tmp || mv build/flags.tmp $@
.PHONY: FORCE
FORCE:

build/%.o: src/%.c include/tmc.h src/internal.h build/flags.txt
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

build/libtmc.a: $(LIB_OBJECTS)
	$(AR) rcs $@ $^

# ---- generated mechanisms ----------------------------------------------
configs/ts1.json: tools/translate_ts1.py $(MICM)/test/benchmark/ts1_mechanism.hpp .conda/conda-meta/history
	$(PYTHON) tools/translate_ts1.py $(MICM)/test/benchmark/ts1_mechanism.hpp $@

generated/ts1.c generated/ts1.h: configs/ts1.json tools/generate.py tools/mechanism.py .conda/conda-meta/history
	@mkdir -p generated
	$(PYTHON) tools/generate.py configs/ts1.json generated/ts1 --name ts1 --width $(WIDTH) --chunk $(KERNEL_CHUNK) --mode $(KERNEL_MODE)

build/ts1.o: generated/ts1.c generated/ts1.h include/tmc.h build/flags.txt
	$(CC) $(CFLAGS) -Igenerated -c $< -o $@

bin/ts1_bench_cpu: bench/ts1_bench.c build/ts1.o build/libtmc.a generated/ts1.h
	@mkdir -p bin
	$(CC) $(CFLAGS) -Igenerated bench/ts1_bench.c build/ts1.o build/libtmc.a $(LDLIBS) -o $@

# ---- tests --------------------------------------------------------------
test-cpu: cpu bin/test_cpu
	bin/test_cpu $(TEST_FILTER)

bin/test_cpu: $(wildcard tests/*.c) $(wildcard tests/*.h) build/libtmc.a
	@mkdir -p bin
	$(CC) $(CFLAGS) -Itests $(wildcard tests/*.c) build/libtmc.a $(LDLIBS) -o $@

# ---- MICM reference -----------------------------------------------------
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: $(REFERENCE)/micm_ts1_ref.cpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(MICM)/include -I$(MICM)/test/benchmark $< -o $@

clean:
	rm -rf build bin generated

# ---- development tools (not part of the entry's targets) -----------------
bin/kernel_bench: bench/kernel_bench.c build/ts1.o build/libtmc.a generated/ts1.h
	@mkdir -p bin
	$(CC) $(CFLAGS) -Igenerated -Isrc bench/kernel_bench.c build/ts1.o build/libtmc.a $(LDLIBS) -o $@
