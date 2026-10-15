.DEFAULT_GOAL := cpu
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

SOURCES = src/version.c src/system.c src/matrix.c src/rates.c src/mechanism.c src/state.c src/sparse.c src/rosenbrock.c src/rosenbrock_parameters.c src/backward_euler.c src/constraints.c src/constraint_initialization.c src/jacobian_verification.c
OBJECTS = $(SOURCES:src/%.c=build/%.o)
cpu: build/libchem.a
build/%.o: src/%.c include/chem.h include/chem_solver.h include/chem_constraints.h include/chem_jacobian.h src/internal.h include/chem_matrix.h include/chem_system.h include/chem_version.h
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@
build/libchem.a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

bin/test_cpu: tests/main.c tests/test_rates.c tests/test_sparse.c tests/test_backward_euler.c tests/test_constraints.c tests/test_constraint_initialization.c tests/test_dae.c tests/test_rosenbrock.c tests/test_jacobian_verification.c tests/test_process.c tests/test_rate_store.c tests/test_state.c tests/test_matrix.c tests/test_matrix_functions.c tests/test_system.c tests/test_runtime.c tests/test.h build/libchem.a
	@mkdir -p bin
	$(CC) $(CFLAGS) tests/main.c tests/test_rates.c tests/test_sparse.c tests/test_backward_euler.c tests/test_constraints.c tests/test_constraint_initialization.c tests/test_dae.c tests/test_rosenbrock.c tests/test_jacobian_verification.c tests/test_process.c tests/test_rate_store.c tests/test_state.c tests/test_matrix.c tests/test_matrix_functions.c tests/test_system.c tests/test_runtime.c build/libchem.a $(LDLIBS) -o $@
test-cpu: env bin/test_cpu bin/test_config bin/ts1_bench_cpu
	$(PYTHON) tools/run_tests.py

build/ts1.c: mechanisms/ts1.json tools/configure.py | env
	$(PYTHON) tools/configure.py $< $@ --symbol ts1
build/ts1.o: build/ts1.c include/chem.h
	$(CC) $(CFLAGS) -c $< -o $@
bin/ts1_bench_cpu: src/ts1_bench.c build/ts1.o build/libchem.a include/chem_solver.h
	@mkdir -p bin
	$(CC) $(CFLAGS) src/ts1_bench.c build/ts1.o build/libchem.a $(LDLIBS) -o $@
cpu: bin/ts1_bench_cpu
bench: bin/ts1_bench_cpu

.PHONY: FORCE
FORCE:
build/compiler-settings: FORCE
	@mkdir -p build
	@printf '%s\n' '$(CC)' '$(CFLAGS)' '$(LDLIBS)' > build/compiler-settings.tmp
	@cmp -s build/compiler-settings.tmp $@ || cp build/compiler-settings.tmp $@
	@rm build/compiler-settings.tmp
$(OBJECTS) build/ts1.o bin/test_cpu bin/ts1_bench_cpu: build/compiler-settings

build/no_user.c: mechanisms/tutorial/rate_constants_no_user_defined/config.yaml mechanisms/tutorial/rate_constants_no_user_defined/species.yaml mechanisms/tutorial/rate_constants_no_user_defined/reactions.yaml tools/configure.py | env
	$(PYTHON) tools/configure.py $< $@ --symbol no_user
build/with_user.c: mechanisms/tutorial/rate_constants_user_defined/config.json mechanisms/tutorial/rate_constants_user_defined/species.json mechanisms/tutorial/rate_constants_user_defined/reactions.json tools/configure.py | env
	$(PYTHON) tools/configure.py $< $@ --symbol with_user
build/robertson.c: mechanisms/tutorial/robertson/config.json mechanisms/tutorial/robertson/species.json mechanisms/tutorial/robertson/reactions.json tools/configure.py | env
	$(PYTHON) tools/configure.py $< $@ --symbol robertson
build/no_user.o build/with_user.o build/robertson.o: build/%.o: build/%.c include/chem.h build/compiler-settings
	$(CC) $(CFLAGS) -c $< -o $@
bin/test_config: tests/config_equivalence.c tests/test.h build/no_user.o build/with_user.o build/robertson.o build/ts1.o build/libchem.a build/compiler-settings
	$(CC) $(CFLAGS) tests/config_equivalence.c build/no_user.o build/with_user.o build/robertson.o build/ts1.o build/libchem.a $(LDLIBS) -o $@
.PHONY: test-config
test-config: bin/test_config
	bin/test_config

$(OBJECTS): src/solver_internal.h
