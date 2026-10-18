# tmc: a C11 port of MICM's CPU solver, as tuned for NCAR's Derecho (2x AMD
# EPYC 7763 Milan per node; gcc 14.3.0 from ncarenv/25.10, see ../env.sh).
# The defaults below (compiler, OPTFLAGS, WIDTH, paths) are Derecho's; every
# other target and rule is the entry's, see ../../fable/Makefile.
#
#   make env        create the conda environment in .conda (Python for build and test only)
#   make cpu        build the library and bin/ts1_bench_cpu
#   make bench      alias of cpu
#   make test-cpu   run the test suite: one line per MICM case, "<id> PASS|FAIL|N/A <reason>"
#   make micm-ref   build MICM's reference TS1 driver with $(CXX) and OPTFLAGS
#
# The compiler comes from CC and the optimization flags from OPTFLAGS; the
# Makefile adds only language, warning, include, OpenMP, and link flags.
.DEFAULT_GOAL := cpu

ifeq ($(origin CC),default)
CC := gcc
endif
ifeq ($(origin CXX),default)
CXX := g++
endif
# Zen 3 (AVX2, FMA): GCC contracts a*b+c into fused multiply-adds only when told
# to in -std=c11 mode, and the kernels were written for contraction.
OPTFLAGS ?= -O3 -DNDEBUG -march=znver3 -ffp-contract=fast
ROOT ?= ..
MICM ?= $(ROOT)/micm
REFERENCE ?= $(ROOT)/../reference
# conda from PATH, else Derecho's system miniforge (module conda, not loaded: see ../env.sh)
CONDA ?= $(shell command -v conda 2>/dev/null || echo /glade/u/apps/opt/miniforge/default/condabin/conda)
PYTHON := $(CURDIR)/.conda/bin/python

# Tuning parameters (see DESIGN.md). Every value gives correct results.
WIDTH ?= 4
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
bench: cpu

# Rebuild everything when the compiler or flags change.
build/flags.txt: FORCE
	@mkdir -p build
	@printf '%s\n' "$(CC)" "$(CFLAGS)" "$(LDLIBS)" > build/flags.tmp
	@cmp -s build/flags.tmp $@ && rm build/flags.tmp || mv build/flags.tmp $@
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

bin/ts1_bench_cpu: bench/ts1_bench.c build/ts1.o build/libtmc.a generated/ts1.h $(wildcard include/*.h)
	@mkdir -p bin
	$(CC) $(CFLAGS) -Igenerated bench/ts1_bench.c build/ts1.o build/libtmc.a $(LDLIBS) -o $@

# ---- tests --------------------------------------------------------------
TEST_SOURCES := $(wildcard tests/*.c)
EXAMPLES := robertson rate_constants_no_user_defined rate_constants_user_defined
TEST_MECHANISMS := $(patsubst %,generated/%.o,$(EXAMPLES))

test-cpu: cpu bin/test_cpu
	TMC_TEST_MICM=$(MICM) TMC_TEST_REFERENCE=$(REFERENCE) bin/test_cpu $(TEST_FILTER)
	@TMC_TEST_MICM=$(MICM) TMC_TEST_REFERENCE=$(REFERENCE) bin/test_cpu --extra 1>&2

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

# ---- MICM reference -----------------------------------------------------
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: $(REFERENCE)/micm_ts1_ref.cpp
	@mkdir -p bin
	$(CXX) -std=c++20 $(OPTFLAGS) -I$(MICM)/include -I$(MICM)/test/benchmark $< -o $@

clean:
	rm -rf build bin generated

# ---- development tools (not part of the entry's targets) -----------------
bin/kernel_bench: bench/kernel_bench.c build/ts1.o build/libtmc.a generated/ts1.h $(wildcard include/*.h)
	@mkdir -p bin
	$(CC) $(CFLAGS) -Igenerated -Isrc bench/kernel_bench.c build/ts1.o build/libtmc.a $(LDLIBS) -o $@
