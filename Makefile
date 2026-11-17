ifeq ($(origin CC),default)
CC = clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
PYTHON = .conda/bin/python
CONDA ?= $(HOME)/miniforge3/bin/conda
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Iinclude $(OPTFLAGS)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
endif
ifeq ($(UNAME_S),Darwin)
ifneq (,$(findstring clang,$(shell $(CC) --version | head -1)))
OMP_CFLAGS = -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include
OMP_LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp
else
OMP_CFLAGS = -fopenmp
OMP_LDFLAGS = -fopenmp
endif
else
OMP_CFLAGS = -fopenmp
OMP_LDFLAGS = -fopenmp
endif

.PHONY: env cpu test-cpu test-reader bench micm-ref clean

env:
	CONDA_PKGS_DIRS="$(CURDIR)/.conda-pkgs" $(CONDA) env create --prefix "$(CURDIR)/.conda" --file environment.yml

build bin:
	mkdir -p $@

build/chem.o: src/chem.c include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@

build/ts1.c: configs/ts1/config.json configs/ts1/species.json configs/ts1/reactions.json tools/generate.py tools/mechanism.py | build
	$(PYTHON) tools/generate.py $< $@

build/ts1.o: build/ts1.c include/chem.h
	$(CC) $(CFLAGS) -c $< -o $@

build/rosenbrock.o: src/rosenbrock.c include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@

build/libchem.a: build/chem.o build/rosenbrock.o
	ar rcs $@ $^

cpu: bin/ts1_bench_cpu

micm-ref: | bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark ../reference/micm_ts1_ref.cpp -o bin/micm_ts1_ref

bench: bin/ts1_bench_cpu

bin/ts1_bench_cpu: src/ts1_bench.c build/ts1.o build/libchem.a | bin
	$(CC) $(CFLAGS) $(OMP_CFLAGS) $^ $(OMP_LDFLAGS) -lm -o $@

build/test_chem: tests/test_chem.c build/ts1.o build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@

build/test_solver: tests/test_solver.c build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@

build/test_cases: tests/test_cases.c build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@

test-cpu:
	@$(MAKE) --no-print-directory -s build/test_cases
	@$(PYTHON) tools/test_cpu.py

test-reader: build/test_chem build/test_solver
	$(PYTHON) -m unittest discover -s tests -v
	./build/test_chem
	./build/test_solver

clean:
	rm -rf build bin
