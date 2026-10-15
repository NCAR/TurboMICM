ifeq ($(origin CC),default)
CC = clang
endif
.DEFAULT_GOAL := cpu
OPTFLAGS ?= -O3 -DNDEBUG
PYTHON = $(CURDIR)/.conda/bin/python
CONDA = $(HOME)/miniforge3/bin/conda
UNAME_S := $(shell uname -s)
COMPILER_VERSION := $(shell $(CC) --version 2>/dev/null)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
ifneq (,$(findstring clang,$(COMPILER_VERSION)))
OMPFLAGS = -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include
OMPLIBS = -L/opt/homebrew/opt/libomp/lib -lomp
else
OMPFLAGS = -fopenmp
OMPLIBS = -fopenmp
endif
else
OMPFLAGS = -fopenmp
OMPLIBS = -fopenmp
endif
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Iinclude $(OMPFLAGS) $(OPTFLAGS)
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
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@
clean:
	@for d in build bin; do if [ -d "$$d" ]; then rm -r "$$d"; fi; done

LIBSRC = src/mechanism.c src/state.c src/rates.c src/process.c src/linear.c src/parameters.c src/rosenbrock.c src/backward_euler.c src/constraints.c src/array.c src/sparse.c src/system.c src/jacobian.c
LIBHEADERS = include/chem.h include/chem_array.h include/chem_sparse.h include/chem_system.h include/chem_jacobian.h src/internal.h
build/ts1.c: configs/ts1.json tools/configure.py .conda/.ready
	@mkdir -p build
	$(PYTHON) tools/configure.py $< $@ --symbol ts1_data
bin/ts1_bench_cpu: build/compiler-flags $(LIBSRC) src/bench.c $(LIBHEADERS) build/ts1.c
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
