CC = clang
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
env: .conda/bin/python
.conda/bin/python: environment.yml
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes
cpu: bin/ts1_bench_cpu
bench: cpu
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@
clean:
	rm -rf build bin

LIBSRC = src/mechanism.c src/state.c src/rates.c src/process.c src/linear.c src/parameters.c src/rosenbrock.c
build/ts1.c: configs/ts1.json tools/configure.py .conda/bin/python
	@mkdir -p build
	$(PYTHON) tools/configure.py $< $@ --symbol ts1_data
bin/ts1_bench_cpu: $(LIBSRC) src/bench.c include/chem.h src/internal.h build/ts1.c
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LIBSRC) src/bench.c build/ts1.c -o $@ $(LDLIBS)
