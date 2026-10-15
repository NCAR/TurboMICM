ifeq ($(origin CC),default)
CC = clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON := $(CURDIR)/.conda/bin/python
UNAME_S := $(shell uname -s)
COMPILER_VERSION := $(shell $(CC) --version 2>/dev/null)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
endif
ifneq (,$(findstring clang,$(COMPILER_VERSION)))
ifeq ($(UNAME_S),Darwin)
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
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic $(OPTFLAGS) $(OMPFLAGS) -Iinclude
LDLIBS = $(OMPLIBS) -lm
.PHONY: env cpu test-cpu bench micm-ref clean
env: .conda/bin/python
.conda/bin/python: environment.yml
	CONDA_PKGS_DIRS=$(CURDIR)/.conda_pkgs $(CONDA) env create --yes --prefix $(CURDIR)/.conda --file environment.yml

micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I ../micm/include -I ../micm/test/benchmark $< -o $@

clean:
	rm -rf build bin

SOURCES = src/rates.c src/mechanism.c src/state.c
OBJECTS = $(SOURCES:src/%.c=build/%.o)
cpu: build/libchem.a
build/%.o: src/%.c include/chem.h src/internal.h
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@
build/libchem.a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

bin/test_cpu: tests/main.c tests/test_rates.c tests/test.h build/libchem.a
	@mkdir -p bin
	$(CC) $(CFLAGS) tests/main.c tests/test_rates.c build/libchem.a $(LDLIBS) -o $@
test-cpu: env bin/test_cpu
	$(PYTHON) tools/run_tests.py
