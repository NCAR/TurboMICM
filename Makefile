CC = clang
OPTFLAGS ?= -O3 -DNDEBUG
CFLAGS = $(OPTFLAGS) -std=c11 -Wall -Wextra -Wpedantic -Iinclude
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
endif

.PHONY: env cpu check clean
env:
	CONDA_PKGS_DIRS="$(CURDIR)/.conda-pkgs" $(HOME)/miniforge3/bin/conda env create --prefix "$(CURDIR)/.conda" --file environment.yml
cpu: build/libchem.a
build:
	mkdir -p $@
build/chem.o: src/chem.c include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/rosenbrock.o: src/rosenbrock.c include/rosenbrock.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/rates.o: src/rates.c include/rates.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/libchem.a: build/chem.o build/rosenbrock.o build/rates.o
	ar rcs $@ $^
build/test_chem: tests/test_chem.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
clean:
	rm -rf build bin
build/robertson.c: tools/generate.py tools/config.py | build
	.conda/bin/python tools/generate.py ../micm/test/tutorial/configs/robertson/config.json $@
build/test_generated: tests/test_generated.c build/robertson.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_generated.c build/robertson.c build/libchem.a -lm -o $@
build/test_rosenbrock: tests/test_rosenbrock.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
build/test_chemical_solver: tests/test_chemical_solver.c build/robertson.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_chemical_solver.c build/robertson.c build/libchem.a -lm -o $@
check: build/test_chem build/test_generated build/test_rosenbrock build/test_chemical_solver
	./build/test_chem
	./build/test_rosenbrock
	./build/test_chemical_solver
	./build/test_tutorial
	./build/test_tutorial_custom
	./build/test_ts1_rates
	./build/test_generated
	.conda/bin/python tests/test_config.py
	.conda/bin/python tests/test_ts1_config.py
build/tutorial.c: tools/generate.py tools/config.py | build
	.conda/bin/python tools/generate.py ../micm/test/tutorial/configs/rate_constants_no_user_defined/config.yaml $@
build/tutorial_custom.c: tools/generate.py tools/config.py | build
	.conda/bin/python tools/generate.py ../micm/test/tutorial/configs/rate_constants_user_defined/config.json $@
build/test_tutorial: tests/test_tutorial.c build/tutorial.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_tutorial.c build/tutorial.c build/libchem.a -lm -o $@
build/test_tutorial_custom: tests/test_tutorial.c build/tutorial_custom.c build/libchem.a
	$(CC) $(CFLAGS) -DWITH_CUSTOM tests/test_tutorial.c build/tutorial_custom.c build/libchem.a -lm -o $@
check: build/test_tutorial build/test_tutorial_custom
build/ts1.c: tools/generate.py tools/config.py configs/ts1/config.json configs/ts1/species.json configs/ts1/reactions.json | build
	.conda/bin/python tools/generate.py configs/ts1/config.json $@
build/test_ts1_rates: tests/test_ts1_rates.c build/ts1.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_ts1_rates.c build/ts1.c build/libchem.a -lm -o $@
check: build/test_ts1_rates

# Apple clang needs Homebrew libomp; GCC and non-Apple compilers use -fopenmp.
ifeq ($(UNAME_S),Darwin)
ifneq (,$(findstring clang,$(shell $(CC) --version | head -1)))
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
bin:
	mkdir -p $@
bin/ts1_bench_cpu: bench/ts1_bench.c build/ts1.c build/libchem.a | bin
	$(CC) $(CFLAGS) $(OMP_CFLAGS) bench/ts1_bench.c build/ts1.c build/libchem.a $(OMP_LIBS) -lm -o $@
.PHONY: bench micm-ref check-ts1
bench: bin/ts1_bench_cpu
check-ts1: bin/ts1_bench_cpu
	./bin/ts1_bench_cpu --threads 10 100 100 30 ../reference/ts1_cells_100_seed1.csv
	.conda/bin/python ../reference/compare_final.py ts1_final.csv ../reference/micm_ts1_final_100_seed1.csv
cpu: bin/ts1_bench_cpu
micm-ref: | bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark ../reference/micm_ts1_ref.cpp -o bin/micm_ts1_ref
