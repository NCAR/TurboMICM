ifeq ($(origin CC),default)
CC := clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
CONDA ?= conda
ENV := $(CURDIR)/.conda
PYTHON := $(ENV)/bin/python
ROOT := $(abspath ..)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
ifneq ($(findstring clang,$(shell $(CC) --version 2>/dev/null | head -1)),)
LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null || printf '%s' /opt/homebrew/opt/libomp)
OMP_CFLAGS := -Xpreprocessor -fopenmp -isystem $(LIBOMP_PREFIX)/include
OMP_LIBS := -L$(LIBOMP_PREFIX)/lib -lomp
else
OMP_CFLAGS := -fopenmp
OMP_LIBS := -fopenmp
endif
else
OMP_CFLAGS := -fopenmp
OMP_LIBS := -fopenmp
endif
CFLAGS := -std=c11 -Wall -Wextra -Wpedantic $(OPTFLAGS) $(OMP_CFLAGS) -Iinclude
LDLIBS := $(OMP_LIBS) -lm

.PHONY: env clean micm-ref cpu bench test-core test-cpu FORCE
env:
	@if [ ! -x "$(PYTHON)" ]; then CONDA_PKGS_DIRS="$(CURDIR)/.conda-pkgs" $(CONDA) env create --prefix "$(ENV)" --file environment.yml; fi

bin:
	@mkdir -p $@

build:
	@mkdir -p $@

build/.compiler: FORCE | build
	@printf '%s\n' '$(CC)|$(OPTFLAGS)|$(OMP_CFLAGS)' | cmp -s - $@ || printf '%s\n' '$(CC)|$(OPTFLAGS)|$(OMP_CFLAGS)' > $@

SOURCES := src/state.c src/rates.c src/linear.c src/rosenbrock.c src/backward_euler.c src/constraints.c src/projection.c src/matrix.c src/sparse.c
OBJECTS := $(SOURCES:src/%.c=build/%.o)

build/%.o: src/%.c src/internal.h include/chem.h include/chem_matrix.h include/chem_sparse.h build/.compiler | build
	@$(CC) $(CFLAGS) -c $< -o $@

build/libchem.a: $(OBJECTS)
	@$(AR) rcs $@ $(OBJECTS)

build/ts1.c: mechanisms/ts1.json tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --reorder

build/ts1.o: build/ts1.c include/chem.h build/.compiler
	@$(CC) $(CFLAGS) -c $< -o $@

bin/ts1_bench_cpu: src/ts1_bench.c build/ts1.o build/libchem.a | bin
	@$(CC) $(CFLAGS) $< build/ts1.o build/libchem.a $(LDLIBS) -o $@

bin/test_core: tests/test_core.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_rates: tests/reference_rates.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_rate_store: tests/reference_rate_store.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_solver: tests/reference_solver.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_constraints: tests/reference_constraints.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_constraint_set: tests/reference_constraint_set.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_equilibrium_integration: tests/reference_equilibrium_integration.c build/constrained_original.o build/constrained_ordered.o build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/constrained_original.o build/constrained_ordered.o build/libchem.a $(LDLIBS) -o $@

bin/reference_stiff_constraints: tests/reference_stiff_constraints.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_analytical: tests/reference_analytical.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_nonlinear_integration: tests/reference_nonlinear_integration.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_state: tests/reference_state.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_matrix: tests/reference_matrix.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_matrix_functions: tests/reference_matrix_functions.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_sparse: tests/reference_sparse.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_sparse_functions: tests/reference_sparse_functions.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_grouped_views: tests/reference_grouped_views.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_initialization: tests/reference_initialization.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/reference_lu: tests/reference_lu.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/test_constraints: tests/test_constraints.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

bin/test_linear: tests/test_linear.c build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/libchem.a $(LDLIBS) -o $@

build/tutorial_no_user.c: tests/fixtures/rate_constants_no_user_defined/config.yaml tests/fixtures/rate_constants_no_user_defined/species.yaml tests/fixtures/rate_constants_no_user_defined/reactions.yaml tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --symbol tutorial_no_user

build/tutorial_user.c: tests/fixtures/rate_constants_user_defined/config.json tests/fixtures/rate_constants_user_defined/species.json tests/fixtures/rate_constants_user_defined/reactions.json tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --symbol tutorial_user

build/tutorial_ordered.c: tests/fixtures/reorder.json tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --symbol tutorial_ordered --reorder

build/tutorial_unordered.c: tests/fixtures/reorder.json tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --symbol tutorial_unordered

build/tutorial_robertson.c: tests/fixtures/robertson/config.json tests/fixtures/robertson/species.json tests/fixtures/robertson/reactions.json tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --symbol tutorial_robertson

build/tutorial_parameterized.c: tests/fixtures/parameterized.json tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --symbol tutorial_parameterized

build/constrained_original.c: tests/fixtures/constrained_reorder.json tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --symbol constrained_original

build/constrained_ordered.c: tests/fixtures/constrained_reorder.json tools/config.py tools/generate.py
	@$(PYTHON) tools/generate.py $< $@ --symbol constrained_ordered --reorder

build/tutorial_%.o: build/tutorial_%.c include/chem.h build/.compiler
	@$(CC) $(CFLAGS) -c $< -o $@

build/constrained_%.o: build/constrained_%.c include/chem.h build/.compiler
	@$(CC) $(CFLAGS) -c $< -o $@

bin/test_configurations: tests/test_configurations.c build/tutorial_no_user.o build/tutorial_user.o build/tutorial_ordered.o build/tutorial_unordered.o build/tutorial_robertson.o build/tutorial_parameterized.o build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/tutorial_no_user.o build/tutorial_user.o build/tutorial_ordered.o build/tutorial_unordered.o build/tutorial_robertson.o build/tutorial_parameterized.o build/libchem.a $(LDLIBS) -o $@

bin/tutorial_program: tests/tutorial_program.c build/tutorial_no_user.o build/tutorial_user.o build/libchem.a build/.compiler | bin
	@$(CC) $(CFLAGS) $< build/tutorial_no_user.o build/tutorial_user.o build/libchem.a $(LDLIBS) -o $@

test-core: bin/test_core bin/test_configurations bin/test_constraints bin/test_linear
	./bin/test_core
	./bin/test_configurations
	./bin/test_constraints
	./bin/test_linear

test-cpu: bin/reference_rates bin/reference_rate_store bin/reference_solver bin/reference_constraints bin/reference_constraint_set bin/reference_equilibrium_integration bin/reference_stiff_constraints bin/reference_analytical bin/reference_nonlinear_integration bin/reference_state bin/reference_matrix bin/reference_matrix_functions bin/reference_sparse bin/reference_sparse_functions bin/reference_grouped_views bin/reference_initialization bin/reference_lu bin/tutorial_program bin/test_core bin/test_configurations bin/test_constraints bin/test_linear
	@$(PYTHON) tests/run_cases.py

cpu bench: bin/ts1_bench_cpu

micm-ref: | bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(ROOT)/micm/include -I$(ROOT)/micm/test/benchmark $(ROOT)/reference/micm_ts1_ref.cpp -o bin/micm_ts1_ref

clean:
	rm -rf build bin
