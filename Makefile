ifeq ($(origin CC),default)
CC = clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON = $(CURDIR)/.conda/bin/python
UNAME := $(shell uname -s)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -1)
ifeq ($(UNAME),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
ifneq (,$(findstring clang,$(CC_VERSION)))
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
CPPFLAGS += -Iinclude -Igenerated
CFLAGS = $(OPTFLAGS) -std=c11 -Wall -Wextra -Wpedantic $(OMPFLAGS)
LDLIBS = $(OMPLIBS) -lm
.PHONY: env cpu bench test-cpu micm-ref clean
env: .conda/bin/python
.conda/bin/python: environment.yml
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml
cpu: bin/ts1_bench_cpu
bench: cpu
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@
clean:
	rm -rf build bin generated

LIB_SOURCES = src/rates.c src/model.c src/state.c src/process.c src/parameters.c src/solver.c
LIB_OBJECTS = $(LIB_SOURCES:src/%.c=build/%.o)
build/%.o: src/%.c include/chem.h src/internal.h
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
build/libchem.a: $(LIB_OBJECTS)
	ar rcs $@ $^

generated/ts1.c: mechanisms/ts1.json tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py $< $@ --symbol ts1
build/ts1.o: generated/ts1.c include/chem.h
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

GROUP_SIZE ?= 128
bin/ts1_bench_cpu: src/bench.c build/ts1.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCK_DEFAULT_GROUP_SIZE=$(GROUP_SIZE) $< build/ts1.o build/libchem.a $(LDLIBS) -o $@
