ifeq ($(origin CC),default)
CC = clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
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
OMP_CFLAGS = -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include
OMP_LIBS = -L/opt/homebrew/opt/libomp/lib -lomp
else
OMP_CFLAGS = -fopenmp
OMP_LIBS = -fopenmp
endif
else
OMP_CFLAGS = -fopenmp
OMP_LIBS = -fopenmp
endif
CPPFLAGS += -Iinclude -Ibuild/generated
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic $(OPTFLAGS) $(OMP_CFLAGS) -fPIC
LDLIBS = $(OMP_LIBS) -lm
SOURCES = $(wildcard src/*.c)
OBJECTS = $(patsubst src/%.c,build/%.o,$(SOURCES))

.PHONY: all cpu env bench test-cpu micm-ref clean FORCE
all: cpu
env: .conda/bin/python
.conda/bin/python: environment.yml
	CONDA_PKGS_DIRS='$(CURDIR)/.conda-pkgs' $(CONDA) env create --yes --prefix '$(CURDIR)/.conda' --file environment.yml

build/flags: FORCE
	@mkdir -p build
	@echo '$(CC) $(CPPFLAGS) $(CFLAGS) $(LDLIBS)' > build/flags.tmp
	@cmp -s build/flags.tmp $@ || cp build/flags.tmp $@

build/%.o: src/%.c include/chem.h src/internal.h build/flags
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

lib/libchem.so: $(OBJECTS)
	@mkdir -p lib
	$(CC) $(OPTFLAGS) $(SHARED) $^ $(LDLIBS) -o $@

cpu: lib/libchem.so bin/ts1_bench_cpu
bench: cpu

build/generated/ts1.c: configs/ts1.json tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py $< $@ --name ts1

build/ts1.o: build/generated/ts1.c include/chem.h build/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

bin/ts1_bench_cpu: tools/ts1_bench.c build/ts1.o $(OBJECTS) include/chem.h build/flags
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) tools/ts1_bench.c build/ts1.o $(OBJECTS) $(LDLIBS) -o $@

test-cpu: cpu
	$(PYTHON) tests/run_cases.py

micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp build/flags
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@

clean:
	rm -rf build bin lib
