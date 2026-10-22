ifeq ($(origin CC),default)
CC = clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
GROUP_SIZE ?= 32
KERNEL_BLOCK ?= 32
KERNEL_FIXED ?= 1
KERNEL_LU ?= mozart
KERNEL_WIDTH ?= $(GROUP_SIZE)
KERNEL_FUSION ?= 1
GENERATOR_OPTIONS = --block $(KERNEL_BLOCK) --fixed $(KERNEL_FIXED) --lu $(KERNEL_LU) --fusion $(KERNEL_FUSION)
LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON = $(CURDIR)/.conda/bin/python
UNAME_S := $(shell uname -s)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -1)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
SHARED = -dynamiclib
else
SHARED = -shared
endif
ifneq (,$(findstring clang,$(CC_VERSION)))
ifeq ($(UNAME_S),Darwin)
OMP_CFLAGS = -Xpreprocessor -fopenmp -isystem $(LIBOMP_PREFIX)/include
OMP_LIBS = -L$(LIBOMP_PREFIX)/lib -lomp
else
OMP_CFLAGS = -fopenmp
OMP_LIBS = -fopenmp
endif
else
OMP_CFLAGS = -fopenmp
OMP_LIBS = -fopenmp
endif
CPPFLAGS += -Iinclude -Isrc -Itools -Ibuild/generated -DCH_DEFAULT_GROUP_SIZE=$(GROUP_SIZE) -DCH_KERNEL_WIDTH=$(KERNEL_WIDTH)
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic $(OPTFLAGS) $(OMP_CFLAGS) -fPIC
LDLIBS = $(OMP_LIBS) -lm
SOURCES = $(wildcard src/*.c)
OBJECTS = $(patsubst src/%.c,build/%.o,$(SOURCES))
EXAMPLE_OBJECTS = build/example_no_user.o build/example_user.o build/example_robertson.o
SPECIALIZED_EXAMPLES = no_user user robertson all_rates empty
SPECIALIZED_TESTS = bin/test_specialized_ts1 $(addprefix bin/test_specialized_,$(SPECIALIZED_EXAMPLES))
.SECONDARY:

.PHONY: all cpu env bench test-cpu test-specialized micm-ref clean FORCE
all: cpu
env:
	@if test -x .conda/bin/python; then \
	  CONDA_PKGS_DIRS='$(CURDIR)/.conda-pkgs' $(CONDA) env update --prefix '$(CURDIR)/.conda' --file environment.yml; \
	else \
	  CONDA_PKGS_DIRS='$(CURDIR)/.conda-pkgs' $(CONDA) env create --yes --prefix '$(CURDIR)/.conda' --file environment.yml; \
	fi
.conda/bin/python:
	CONDA_PKGS_DIRS='$(CURDIR)/.conda-pkgs' $(CONDA) env create --yes --prefix '$(CURDIR)/.conda' --file environment.yml

build/flags: FORCE
	@mkdir -p build
	@echo '$(CC) $(CPPFLAGS) $(CFLAGS) $(LDLIBS)' > build/flags.tmp
	@cmp -s build/flags.tmp $@ || cp build/flags.tmp $@

build/%.o: src/%.c include/chem.h include/chem_array.h src/internal.h build/flags
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

lib/libchem.so: $(OBJECTS)
	@mkdir -p lib
	$(CC) $(OPTFLAGS) $(SHARED) $^ $(LDLIBS) -o $@

cpu: lib/libchem.so bin/ts1_bench_cpu
bench: cpu bin/ts1_bench_specialized bin/ts1_scenario_specialized

build/generator-flags: FORCE
	@mkdir -p build
	@echo '$(GENERATOR_OPTIONS)' > $@.tmp
	@cmp -s $@.tmp $@ || cp $@.tmp $@

build/generated/ts1_specialized.c: configs/ts1.json tools/configure.py tools/specialize.py build/generator-flags .conda/bin/python
	$(PYTHON) tools/specialize.py $< $@ --name ts1_specialized $(GENERATOR_OPTIONS)

build/ts1_specialized.o: build/generated/ts1_specialized.c tools/specialized_solver.h src/internal.h include/chem.h build/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

bin/ts1_bench_specialized: tools/ts1_bench.c tools/ts1_io.h build/ts1_specialized.o build/ts1_specialized_kernels.o $(OBJECTS) include/chem.h build/flags
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCH_SPECIALIZED tools/ts1_bench.c build/ts1_specialized.o build/ts1_specialized_kernels.o $(OBJECTS) $(LDLIBS) -o $@

bin/ts1_scenario_specialized: tools/ts1_scenario.c tools/ts1_io.h build/ts1_specialized.o build/ts1_specialized_kernels.o $(OBJECTS) include/chem.h build/flags
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCH_SPECIALIZED tools/ts1_scenario.c build/ts1_specialized.o build/ts1_specialized_kernels.o $(OBJECTS) $(LDLIBS) -o $@

build/generated/no_user_specialized.c: $(wildcard configs/examples/rate_constants_no_user_defined/*) tools/configure.py tools/specialize.py build/generator-flags .conda/bin/python
	$(PYTHON) tools/specialize.py configs/examples/rate_constants_no_user_defined/config.yaml $@ --name no_user_specialized $(GENERATOR_OPTIONS)
build/generated/user_specialized.c: $(wildcard configs/examples/rate_constants_user_defined/*) tools/configure.py tools/specialize.py build/generator-flags .conda/bin/python
	$(PYTHON) tools/specialize.py configs/examples/rate_constants_user_defined/config.json $@ --name user_specialized $(GENERATOR_OPTIONS)
build/generated/robertson_specialized.c: $(wildcard configs/examples/robertson/*) tools/configure.py tools/specialize.py build/generator-flags .conda/bin/python
	$(PYTHON) tools/specialize.py configs/examples/robertson/config.json $@ --name robertson_specialized $(GENERATOR_OPTIONS)
build/generated/test-configs.stamp: tests/generator_cases.py .conda/bin/python
	$(PYTHON) tests/generator_cases.py build/generated
	@touch $@
build/generated/all_rates_specialized.c: build/generated/test-configs.stamp tools/configure.py tools/specialize.py build/generator-flags .conda/bin/python
	$(PYTHON) tools/specialize.py build/generated/all_rates.json $@ --name all_rates_specialized $(GENERATOR_OPTIONS)
build/generated/empty_specialized.c: build/generated/test-configs.stamp tools/configure.py tools/specialize.py build/generator-flags .conda/bin/python
	$(PYTHON) tools/specialize.py build/generated/empty.json $@ --name empty_specialized $(GENERATOR_OPTIONS)
build/%_specialized.o: build/generated/%_specialized.c tools/specialized_solver.h src/internal.h include/chem.h build/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
bin/test_specialized_%: tests/specialized.c build/%_specialized.o build/%_specialized_kernels.o $(OBJECTS) include/chem.h src/internal.h build/flags
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -DSPEC_HEADER='"$*_specialized.h"' -DSPEC_NAME=$*_specialized tests/specialized.c build/$*_specialized.o build/$*_specialized_kernels.o $(OBJECTS) $(LDLIBS) -o $@

test-specialized: $(SPECIALIZED_TESTS)
	@set -e; for test in $(SPECIALIZED_TESTS); do $$test; done

build/generated/ts1.c: configs/ts1.json tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py $< $@ --name ts1

build/ts1.o: build/generated/ts1.c include/chem.h build/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

bin/ts1_bench_cpu: tools/ts1_bench.c tools/ts1_io.h build/ts1.o $(OBJECTS) include/chem.h build/flags
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) tools/ts1_bench.c build/ts1.o $(OBJECTS) $(LDLIBS) -o $@

bin/test_native_cpu: tests/native_cases.c $(OBJECTS) include/chem.h build/flags
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) tests/native_cases.c $(OBJECTS) $(LDLIBS) -o $@

build/generated/example_no_user.c: $(wildcard configs/examples/rate_constants_no_user_defined/*) tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py configs/examples/rate_constants_no_user_defined/config.yaml $@ --name example_no_user
build/generated/example_user.c: $(wildcard configs/examples/rate_constants_user_defined/*) tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py configs/examples/rate_constants_user_defined/config.json $@ --name example_user
build/generated/example_robertson.c: $(wildcard configs/examples/robertson/*) tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py configs/examples/robertson/config.json $@ --name example_robertson
build/example_%.o: build/generated/example_%.c include/chem.h build/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
bin/tutorial_cpu: tools/tutorials.c $(EXAMPLE_OBJECTS) $(OBJECTS) include/chem.h include/chem_array.h build/flags
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) tools/tutorials.c $(EXAMPLE_OBJECTS) $(OBJECTS) $(LDLIBS) -o $@

test-cpu: cpu bin/test_native_cpu bin/tutorial_cpu
	@bin/test_native_cpu extensions
	$(PYTHON) tests/run_cases.py

micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp build/flags
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@

clean:
	rm -rf build bin lib

build/%_specialized_kernels.o: build/generated/%_specialized.c src/internal.h include/chem.h build/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -c build/generated/$*_specialized_kernels.c -o $@
