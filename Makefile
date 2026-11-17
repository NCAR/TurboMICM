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

build/chem.o: src/chem.c src/internal.h include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@

build/ts1.c: configs/ts1/config.json configs/ts1/species.json configs/ts1/reactions.json tools/generate.py tools/mechanism.py | build
	$(PYTHON) tools/generate.py $< $@

build/ts1.o: build/ts1.c include/chem.h
	$(CC) $(CFLAGS) -c $< -o $@

build/rosenbrock.o: src/rosenbrock.c src/internal.h include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@

build/constraint.o: src/constraint.c include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@

build/libchem.a: build/chem.o build/rosenbrock.o build/constraint.o
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
	@$(MAKE) --no-print-directory -s build/test_cases build/tutorial_no_user build/tutorial_user
	@$(PYTHON) tools/test_cpu.py

build/test_dae: tests/test_dae.c build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@

test-reader: build/test_chem build/test_solver build/test_robertson_config build/test_dae build/test_parameterized
	$(PYTHON) -m unittest discover -s tests -v
	./build/test_chem
	./build/test_solver
	./build/test_robertson_config
	./build/test_dae
	./build/test_parameterized

clean:
	rm -rf build bin

build/tutorial_no_user.c: ../micm/test/tutorial/configs/rate_constants_no_user_defined/config.yaml ../micm/test/tutorial/configs/rate_constants_no_user_defined/species.yaml ../micm/test/tutorial/configs/rate_constants_no_user_defined/reactions.yaml tools/generate.py tools/mechanism.py | build
	$(PYTHON) tools/generate.py $< $@

build/tutorial_user.c: ../micm/test/tutorial/configs/rate_constants_user_defined/config.json ../micm/test/tutorial/configs/rate_constants_user_defined/species.json ../micm/test/tutorial/configs/rate_constants_user_defined/reactions.json tools/generate.py tools/mechanism.py | build
	$(PYTHON) tools/generate.py $< $@

build/tutorial_no_user: tests/test_tutorial.c build/tutorial_no_user.c build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@

build/tutorial_user: tests/test_tutorial.c build/tutorial_user.c build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@

build/robertson.c: ../micm/test/tutorial/configs/robertson/config.json ../micm/test/tutorial/configs/robertson/species.json ../micm/test/tutorial/configs/robertson/reactions.json tools/generate.py tools/mechanism.py | build
	$(PYTHON) tools/generate.py $< $@

build/test_robertson_config: tests/test_robertson_config.c build/robertson.c build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@

build/parameterized.c: configs/parameterized/config.json configs/parameterized/species.json configs/parameterized/reactions.json tools/generate.py tools/mechanism.py | build
	$(PYTHON) tools/generate.py $< $@

build/test_parameterized: tests/test_parameterized.c build/parameterized.c build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@
