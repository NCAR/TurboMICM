ifeq ($(origin CC),default)
CC = clang
endif
.DEFAULT_GOAL := cpu
.SUFFIXES:
ifeq ($(origin CXX),default)
CXX = clang++
endif
OPTFLAGS ?= -O3 -DNDEBUG
PYTHON = $(CURDIR)/.conda/bin/python
CONDA = $(HOME)/miniforge3/bin/conda
LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
GROUP_SIZE ?= 32
UNAME_S := $(shell uname -s)
COMPILER_VERSION := $(shell $(CC) --version 2>/dev/null)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
ifneq (,$(findstring clang,$(COMPILER_VERSION)))
OMPFLAGS = -Xpreprocessor -fopenmp -isystem $(LIBOMP_PREFIX)/include
OMPLIBS = -L$(LIBOMP_PREFIX)/lib -lomp
else
OMPFLAGS = -fopenmp
OMPLIBS = -fopenmp
endif
else
OMPFLAGS = -fopenmp
OMPLIBS = -fopenmp
endif
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Iinclude -Isrc -DCHEM_GROUP_SIZE=$(GROUP_SIZE) $(OMPFLAGS) $(OPTFLAGS)
LDLIBS = $(OMPLIBS) -lm
.PHONY: env cpu bench test-cpu micm-ref clean
env: .conda/.ready
.conda/.ready: environment.yml
	@if [ -d .conda/conda-meta ]; then \
	  CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env update --prefix $(CURDIR)/.conda --file environment.yml; \
	else \
	  CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes; \
	fi
	@touch $@
cpu: bin/ts1_bench_cpu build/libchem.so
bench: cpu
micm-ref: bin/micm_ts1_ref bin/micm_scenario
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp
	@mkdir -p bin
	$(CXX) -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@
clean:
	@for d in build bin; do if [ -d "$$d" ]; then rm -r "$$d"; fi; done

LIBSRC = src/mechanism.c src/state.c src/rates.c src/process.c src/linear.c src/parameters.c src/rosenbrock.c src/backward_euler.c src/constraints.c src/array.c src/sparse.c src/system.c src/jacobian.c src/builder.c
LIBHEADERS = src/rosenbrock_impl.h include/chem.h include/chem_array.h include/chem_sparse.h include/chem_system.h include/chem_jacobian.h include/chem_builder.h src/internal.h
build/ts1.c: configs/ts1.json tools/configure.py .conda/.ready
	@mkdir -p build
	$(PYTHON) tools/configure.py $< $@ --symbol ts1_data
bin/ts1_bench_cpu: build/compiler-flags $(LIBSRC) src/bench.c src/bench_io.h $(LIBHEADERS) build/ts1.c
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LIBSRC) src/bench.c build/ts1.c -o $@ $(LDLIBS)

.PHONY: FORCE
FORCE:
build/compiler-flags: FORCE
	@mkdir -p build
	@printf '%s\n' '$(CC)' '$(CFLAGS)' '$(LDLIBS)' > build/compiler-flags.new
	@cmp -s build/compiler-flags.new $@ || cp build/compiler-flags.new $@
ifeq ($(UNAME_S),Darwin)
SHAREDFLAGS = -dynamiclib
else
SHAREDFLAGS = -shared -fPIC
endif
build/libchem.so: build/compiler-flags $(LIBSRC) $(LIBHEADERS)
	$(CC) $(CFLAGS) $(SHAREDFLAGS) $(LIBSRC) -o $@ $(LDLIBS)
test-cpu: cpu
	$(PYTHON) tests/run_cases.py
bin/test_rates: build/compiler-flags $(LIBSRC) tests/rate_cases.c $(LIBHEADERS)
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LIBSRC) tests/rate_cases.c -o $@ $(LDLIBS)

test-cpu: bin/test_rates

build/no_user.c: configs/rate_constants_no_user_defined/config.yaml configs/rate_constants_no_user_defined/species.yaml configs/rate_constants_no_user_defined/reactions.yaml tools/configure.py .conda/.ready
	$(PYTHON) tools/configure.py $< $@ --symbol no_user_data
build/user.c: configs/rate_constants_user_defined/config.json configs/rate_constants_user_defined/species.json configs/rate_constants_user_defined/reactions.json tools/configure.py .conda/.ready
	$(PYTHON) tools/configure.py $< $@ --symbol user_data
build/robertson.c: configs/robertson/config.json configs/robertson/species.json configs/robertson/reactions.json tools/configure.py .conda/.ready
	$(PYTHON) tools/configure.py $< $@ --symbol robertson_data
build/parameterized.c: configs/parameterized/config.json tools/configure.py .conda/.ready
	$(PYTHON) tools/configure.py $< $@ --symbol parameterized_data
build/configs.so: build/compiler-flags build/no_user.c build/user.c build/robertson.c build/parameterized.c include/chem.h
	$(CC) $(CFLAGS) $(SHAREDFLAGS) build/no_user.c build/user.c build/robertson.c build/parameterized.c -o $@ $(LDLIBS)
test-cpu: build/configs.so
bin/rate_tutorial: build/compiler-flags $(LIBSRC) src/tutorial.c $(LIBHEADERS)
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LIBSRC) src/tutorial.c -o $@ $(LDLIBS)
test-cpu: bin/rate_tutorial

bin/cpu_tutorials: build/compiler-flags $(LIBSRC) src/programs.c $(LIBHEADERS)
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LIBSRC) src/programs.c -o $@ $(LDLIBS)
test-cpu: bin/cpu_tutorials

SPECIALIZE_BATCH ?= 1
SPECIALIZE_WIDTH ?= $(GROUP_SIZE)
SPECIALIZE_FLAGS ?=
.PHONY: test-specialized
bench: bin/ts1_bench_specialized bin/ts1_scenario_specialized
build/generator-flags: FORCE
	@mkdir -p build
	@printf '%s\n' '$(SPECIALIZE_BATCH)' '$(SPECIALIZE_WIDTH)' '$(SPECIALIZE_FLAGS)' > $@.new
	@cmp -s $@.new $@ || cp $@.new $@
build/ts1_specialized.c: configs/ts1.json tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix ts1 --batch $(SPECIALIZE_BATCH) --fixed-width $(SPECIALIZE_WIDTH) $(SPECIALIZE_FLAGS)
build/ts1_specialized.o: build/ts1_specialized.c build/compiler-flags $(LIBHEADERS)
	$(CC) $(CFLAGS) -fPIC -c $< -o $@
bin/ts1_bench_specialized: build/compiler-flags $(LIBSRC) src/bench.c src/bench_io.h include/chem_specialized.h $(LIBHEADERS) build/ts1_specialized.o
	@mkdir -p bin
	$(CC) $(CFLAGS) -DCHEM_SPECIALIZED $(LIBSRC) src/bench.c build/ts1_specialized.o -o $@ $(LDLIBS)
bin/ts1_scenario_specialized: build/compiler-flags $(LIBSRC) src/scenario.c src/bench_io.h include/chem_specialized.h $(LIBHEADERS) build/ts1_specialized.o
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LIBSRC) src/scenario.c build/ts1_specialized.o -o $@ $(LDLIBS)

SPECIALIZED_CONFIGS = no_user user robertson parameterized checks inert odd shapes
build/no_user_specialized.c: configs/rate_constants_no_user_defined/config.yaml $(wildcard configs/rate_constants_no_user_defined/*) tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix no_user --batch $(SPECIALIZE_BATCH) --fixed-width $(SPECIALIZE_WIDTH) $(SPECIALIZE_FLAGS)
build/user_specialized.c: configs/rate_constants_user_defined/config.json $(wildcard configs/rate_constants_user_defined/*) tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix user --batch $(SPECIALIZE_BATCH) --fixed-width $(SPECIALIZE_WIDTH) $(SPECIALIZE_FLAGS)
build/robertson_specialized.c: configs/robertson/config.json $(wildcard configs/robertson/*) tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix robertson --batch $(SPECIALIZE_BATCH) --fixed-width $(SPECIALIZE_WIDTH) $(SPECIALIZE_FLAGS)
build/parameterized_specialized.c: configs/parameterized/config.json tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix parameterized --batch $(SPECIALIZE_BATCH) --fixed-width $(SPECIALIZE_WIDTH) $(SPECIALIZE_FLAGS)
build/%_specialized.o: build/%_specialized.c build/compiler-flags $(LIBHEADERS)
	$(CC) $(CFLAGS) -fPIC -c $< -o $@
build/specialized_tests.so: build/compiler-flags $(LIBSRC) $(LIBHEADERS) build/ts1_specialized.o $(addprefix build/,$(addsuffix _specialized.o,$(SPECIALIZED_CONFIGS)))
	$(CC) $(CFLAGS) $(SHAREDFLAGS) $(LIBSRC) build/ts1_specialized.o $(addprefix build/,$(addsuffix _specialized.o,$(SPECIALIZED_CONFIGS))) -o $@ $(LDLIBS)
test-specialized: build/specialized_tests.so build/libchem.so bin/ts1_bench_specialized
	$(PYTHON) tests/test_specialized.py

build/checks_specialized.c: configs/specialization_checks.json tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix checks --batch $(SPECIALIZE_BATCH) --fixed-width $(SPECIALIZE_WIDTH) $(SPECIALIZE_FLAGS)

# Developer timing variants use the same compiler and flags as release builds.
build/tuning/bench%: build/tuning/batch%_specialized.o $(LIBSRC) $(LIBHEADERS) src/bench.c src/bench_io.h include/chem_specialized.h
	$(CC) $(CFLAGS) -DCHEM_SPECIALIZED $(LIBSRC) src/bench.c $< -o $@ $(LDLIBS)

build/inert_specialized.c: configs/inert.json tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix inert --batch $(SPECIALIZE_BATCH) --fixed-width $(SPECIALIZE_WIDTH) $(SPECIALIZE_FLAGS)

.SECONDARY: $(addprefix build/,$(addsuffix _specialized.o,$(SPECIALIZED_CONFIGS)))

build/odd_specialized.c: configs/specialization_checks.json tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix odd --batch 3 --fixed-width 7 --tile 3 $(SPECIALIZE_FLAGS) --unrolled-linear

.PHONY: test-scenarios
test-scenarios: bin/ts1_scenario_specialized
	$(PYTHON) tools/evaluate.py --scenarios --directory build/scenario-gates --label '$(CC)'

build/shapes_specialized.c: configs/specialization_checks.json tools/specialize.py tools/compact_linear.py tools/compact_chemistry.py tools/configure.py build/generator-flags .conda/.ready
	$(PYTHON) tools/specialize.py $< $@ --prefix shapes --batch 3 --fixed-width 7 --tile 3 $(SPECIALIZE_FLAGS) --compact-chemistry --narrow-indices

bin/micm_scenario: ../scenarios/micm_scenario.cpp ../micm/test/benchmark/ts1_mechanism.hpp
	@mkdir -p bin
	$(CXX) -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@
