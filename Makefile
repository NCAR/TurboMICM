.DEFAULT_GOAL := cpu
ifeq ($(origin CC),default)
CC := clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON := $(CURDIR)/.conda/bin/python
MICM ?= ../micm
REFERENCE ?= ../reference
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
CFLAGS := $(OPTFLAGS) -std=c11 -Wall -Wextra -Wpedantic $(OMP_CFLAGS) -Iinclude
LDLIBS := $(OMP_LIBS) -lm

.PHONY: env cpu test-cpu bench micm-ref clean

env: .conda/conda-meta/history
.conda/conda-meta/history: environment.yml
	CONDA_REGISTER_ENVS=false CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes

SOURCES := src/rates.c src/linear.c src/model.c src/state.c src/parameters.c src/rosenbrock.c src/backward_euler.c src/constraints.c src/verify.c src/array.c src/sparse.c src/factor.c src/version.c src/system.c src/process.c src/solver.c src/external.c src/adaptive.c
OBJECTS := $(patsubst src/%.c,build/%.o,$(SOURCES))

cpu: bin/ts1_bench_cpu bin/tutorial build/libchem.a
bench: cpu bin/ts1_scenario

bin/ts1_scenario: build/scenario.o build/ts1.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@

test-cpu: cpu bin/test_cpu
	bin/test_cpu $(TEST_FILTER)

generated/test_registry.inc: $(wildcard tests/test_*.c) tests/cases.txt tools/test_registry.py .conda/conda-meta/history
	$(PYTHON) tools/test_registry.py $@

CONFIG_OBJECTS := build/config_robertson.o build/config_rates.o build/config_user_rates.o

bin/test_cpu: $(wildcard tests/*.c) $(wildcard tests/*.h) generated/test_registry.inc build/chapman.o $(CONFIG_OBJECTS) build/tutorials.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CFLAGS) -Itests -Igenerated -Iexamples $(wildcard tests/*.c) build/tutorials.o build/chapman.o $(CONFIG_OBJECTS) build/libchem.a $(LDLIBS) -o $@

build/%.o: src/%.c $(wildcard include/*.h) $(wildcard src/*.h) build/compiler.conf
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

build/libchem.a: $(OBJECTS)
	$(AR) rcs $@ $^

generated/ts1.c: configs/ts1.json tools/configure.py tools/generate.py .conda/conda-meta/history
	$(PYTHON) tools/generate.py $< $@ --symbol ts1_mechanism

build/ts1.o: generated/ts1.c $(wildcard include/*.h) build/compiler.conf
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

generated/chapman.c: configs/chapman.json tools/configure.py tools/generate.py .conda/conda-meta/history
	$(PYTHON) tools/generate.py $< $@ --symbol chapman_mechanism

build/chapman.o: generated/chapman.c $(wildcard include/*.h) build/compiler.conf
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

bin/ts1_bench_cpu: build/bench.o build/ts1.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@

micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: $(REFERENCE)/micm_ts1_ref.cpp build/reference-compiler.conf
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(MICM)/include -I$(MICM)/test/benchmark $< -o $@

clean:
	rm -rf build bin generated

.PHONY: FORCE
FORCE:

build/compiler.conf: FORCE
	@mkdir -p build
	@printf '%s\n' "$(CC)" "$(CFLAGS)" "$(LDLIBS)" > build/compiler.conf.tmp
	@cmp -s build/compiler.conf.tmp $@ && rm build/compiler.conf.tmp || mv build/compiler.conf.tmp $@

build/tutorials.o: examples/tutorials.c examples/tutorials.h $(wildcard include/*.h) build/compiler.conf
	@mkdir -p build
	$(CC) $(CFLAGS) -Iexamples -c $< -o $@

bin/tutorial: examples/main.c examples/tutorials.h build/tutorials.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CFLAGS) -Iexamples examples/main.c build/tutorials.o build/libchem.a $(LDLIBS) -o $@

generated/robertson.c: $(wildcard configs/robertson/*) tools/configure.py tools/generate.py .conda/conda-meta/history
	$(PYTHON) tools/generate.py configs/robertson/config.json $@ --symbol robertson_mechanism

build/config_robertson.o: generated/robertson.c $(wildcard include/*.h) build/compiler.conf
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

generated/rates.c: $(wildcard configs/rate_constants_no_user_defined/*) tools/configure.py tools/generate.py .conda/conda-meta/history
	$(PYTHON) tools/generate.py configs/rate_constants_no_user_defined/config.yaml $@ --symbol rates_mechanism

build/config_rates.o: generated/rates.c $(wildcard include/*.h) build/compiler.conf
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

generated/user_rates.c: $(wildcard configs/rate_constants_user_defined/*) tools/configure.py tools/generate.py .conda/conda-meta/history
	$(PYTHON) tools/generate.py configs/rate_constants_user_defined/config.json $@ --symbol user_rates_mechanism

build/config_user_rates.o: generated/user_rates.c $(wildcard include/*.h) build/compiler.conf
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

build/reference-compiler.conf: FORCE
	@mkdir -p build
	@printf '%s\n' "clang++" "$(OPTFLAGS)" > build/reference-compiler.conf.tmp
	@cmp -s build/reference-compiler.conf.tmp $@ && rm build/reference-compiler.conf.tmp || mv build/reference-compiler.conf.tmp $@
