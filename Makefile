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
build/rosenbrock.o: src/rosenbrock.c include/rosenbrock.h include/sparse.h include/constraint.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/rosenbrock_group.o: src/rosenbrock_group.c include/rosenbrock.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/rates.o: src/rates.c include/rates.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/rate_store.o: src/rate_store.c include/rate_store.h include/rates.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/sparse.o: src/sparse.c include/sparse.h include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/linear.o: src/linear.c include/linear.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/backward_euler.o: src/backward_euler.c include/backward_euler.h include/rosenbrock.h include/linear.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/backward_euler_group.o: src/backward_euler_group.c include/backward_euler.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/constraint.o: src/constraint.c include/constraint.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/constraint_init.o: src/constraint_init.c include/rosenbrock.h include/constraint.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/jacobian_check.o: src/jacobian_check.c include/jacobian_check.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/state.o: src/state.c include/state.h include/rates.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/species.o: src/species.c include/species.h include/rates.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/process.o: src/process.c include/process.h include/species.h include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/matrix.o: src/matrix.c include/matrix.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/sparse_matrix.o: src/sparse_matrix.c include/sparse_matrix.h include/matrix.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/composite.o: src/composite.c include/composite.h include/rates.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/libchem.a: build/backward_euler_group.o build/process.o build/species.o build/rosenbrock_group.o build/rate_store.o build/composite.o build/sparse_matrix.o build/matrix.o build/state.o build/jacobian_check.o build/chem.o build/rosenbrock.o build/rates.o build/sparse.o build/linear.o build/backward_euler.o build/constraint.o build/constraint_init.o
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
	./build/test_matrix
	./build/test_matrix_views
	./build/test_grouped_views
	./build/test_matrix_function
	./build/test_jacobian_check
	./build/test_process_set
	./build/test_process_random
	./build/test_state
	./build/test_species
	./build/test_process_builder
	./build/test_flow_tube
	./build/test_openmp
	./build/test_rosenbrock
	./build/test_normalized_error
	./build/test_analytical_methods
	./build/test_chain_rates
	./build/test_nonlinear
	./build/test_backward_euler
	./build/test_constraints
	./build/test_external_constraints
	./build/test_external_composition
	./build/test_aerosol
	./build/test_aerosol_state
	./build/test_grouped_constraints
	./build/test_chapman_regression
	./build/test_constraint_set
	./build/test_constraint_init
	./build/test_dae
	./build/test_dae_overshoot
	./build/test_rates
	./build/test_rate_store
	./build/test_parameterized
	./build/test_callback
	./build/test_rate_extensions
	./build/test_generated_rate_extensions
	./build/test_falloff
	./build/test_sparse
	./build/test_sparse_matrix
	./build/test_sparse_views
	./build/test_sparse_functions
	./build/test_doolittle
	./build/test_linear_solve
	./build/test_reordered
	./build/test_chemical_solver
	./build/test_tutorial
	./build/test_tutorial_custom
	.conda/bin/python tools/test_tutorial_programs.py
	./build/test_ts1_rates
	./build/test_generated
	.conda/bin/python tests/test_config.py
	.conda/bin/python tests/test_ts1_config.py
	.conda/bin/python tests/test_markowitz.py
build/tutorial.c: tools/generate.py tools/config.py | build
	.conda/bin/python tools/generate.py ../micm/test/tutorial/configs/rate_constants_no_user_defined/config.yaml $@
build/tutorial_custom.c: tools/generate.py tools/config.py | build
	.conda/bin/python tools/generate.py ../micm/test/tutorial/configs/rate_constants_user_defined/config.json $@
build/test_tutorial: tests/test_tutorial.c build/tutorial.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_tutorial.c build/tutorial.c build/libchem.a -lm -o $@
build/test_tutorial_custom: tests/test_tutorial.c build/tutorial_custom.c build/libchem.a
	$(CC) $(CFLAGS) -DWITH_CUSTOM tests/test_tutorial.c build/tutorial_custom.c build/libchem.a -lm -o $@
check: build/test_tutorial build/test_tutorial_custom
build/rate_constants_no_user_defined_example_by_hand: examples/rate_constants.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
build/rate_constants_user_defined_example_by_hand: examples/rate_constants.c build/libchem.a
	$(CC) $(CFLAGS) -DWITH_CUSTOM $< build/libchem.a -lm -o $@
check: build/rate_constants_no_user_defined_example_by_hand build/rate_constants_user_defined_example_by_hand
build/README_example: examples/solver_examples.c build/libchem.a
	$(CC) $(CFLAGS) -DEXAMPLE=0 $< build/libchem.a -lm -o $@
build/multiple_grid_cells: examples/solver_examples.c build/libchem.a
	$(CC) $(CFLAGS) -DEXAMPLE=1 $< build/libchem.a -lm -o $@
build/solve_results: examples/solver_examples.c build/libchem.a
	$(CC) $(CFLAGS) -DEXAMPLE=2 $< build/libchem.a -lm -o $@
check: build/README_example build/multiple_grid_cells build/solve_results
build/solver_configuration: examples/solver_configuration.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/solver_configuration
build/ts1.c: tools/generate.py tools/config.py configs/ts1/config.json configs/ts1/species.json configs/ts1/reactions.json | build
	.conda/bin/python tools/generate.py --reorder configs/ts1/config.json $@
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
build/reordered.c: tools/generate.py tools/config.py tests/reorder/config.json tests/reorder/species.json tests/reorder/reactions.json | build
	.conda/bin/python tools/generate.py --reorder tests/reorder/config.json $@
build/test_reordered: tests/test_reordered.c build/reordered.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_reordered.c build/reordered.c build/libchem.a -lm -o $@
build/test_sparse: tests/test_sparse.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_sparse.c build/libchem.a -lm -o $@
check: build/test_sparse build/test_reordered
build/test_backward_euler: tests/test_backward_euler.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_backward_euler.c build/libchem.a -lm -o $@
check: build/test_backward_euler
.PHONY: test-cpu
test-cpu:
	@$(MAKE) -s build/test_flow_tube build/test_openmp build/openmp_tutorial build/test_process_builder build/test_species build/solver_configuration build/README_example build/multiple_grid_cells build/solve_results build/rate_constants_no_user_defined_example_by_hand build/rate_constants_user_defined_example_by_hand build/test_aerosol build/test_aerosol_state build/test_chapman_regression build/test_rate_store build/test_dae_overshoot build/test_grouped_views build/test_sparse_functions build/test_sparse_views build/test_sparse_matrix build/test_matrix_function build/test_matrix_views build/test_matrix build/test_external_composition build/test_external_constraints build/test_nonlinear build/test_process_random build/test_chain_rates build/test_state build/test_process_set build/test_jacobian_check build/test_backward_euler build/test_rates build/test_rate_extensions build/test_falloff build/test_normalized_error build/test_doolittle build/test_analytical_methods build/test_constraints build/test_constraint_set build/test_constraint_init build/test_dae build/test_linear_solve build/test_parameterized build/test_callback CC="$(CC)" OPTFLAGS="$(OPTFLAGS)"
	@.conda/bin/python tools/test_cpu.py
build/test_rates: tests/test_rates.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_rates.c build/libchem.a -lm -o $@
check: build/test_rates
build/test_rate_extensions: tests/test_rate_extensions.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_rate_extensions.c build/libchem.a -lm -o $@
check: build/test_rate_extensions
build/test_falloff: tests/test_falloff.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_falloff.c build/libchem.a -lm -o $@
check: build/test_falloff
build/extended_rates.c: tools/generate.py tools/config.py tests/rate_extensions/config.json tests/rate_extensions/species.json tests/rate_extensions/reactions.json | build
	.conda/bin/python tools/generate.py tests/rate_extensions/config.json $@
build/test_generated_rate_extensions: tests/test_generated_rate_extensions.c build/extended_rates.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_generated_rate_extensions.c build/extended_rates.c build/libchem.a -lm -o $@
check: build/test_generated_rate_extensions
build/test_analytical_methods: tests/test_analytical_methods.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_analytical_methods.c build/libchem.a -lm -o $@
check: build/test_analytical_methods
build/test_normalized_error: tests/test_normalized_error.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_normalized_error.c build/libchem.a -lm -o $@
check: build/test_normalized_error
build/test_doolittle: tests/test_doolittle.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_doolittle.c build/libchem.a -lm -o $@
check: build/test_doolittle
build/test_constraints: tests/test_constraints.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_constraints.c build/libchem.a -lm -o $@
check: build/test_constraints
build/test_constraint_set: tests/test_constraint_set.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_constraint_set.c build/libchem.a -lm -o $@
check: build/test_constraint_set
build/test_constraint_init: tests/test_constraint_init.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_constraint_init.c build/libchem.a -lm -o $@
check: build/test_constraint_init
build/test_dae: tests/test_dae.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_dae.c build/libchem.a -lm -o $@
check: build/test_dae
build/test_linear_solve: tests/test_linear_solve.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_linear_solve.c build/libchem.a -lm -o $@
check: build/test_linear_solve
build/parameterized.c: tools/generate.py tools/config.py tests/parameterized/config.json tests/parameterized/species.json tests/parameterized/reactions.json | build
	.conda/bin/python tools/generate.py tests/parameterized/config.json $@
build/test_parameterized: tests/test_parameterized.c build/parameterized.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_parameterized.c build/parameterized.c build/libchem.a -lm -o $@
check: build/test_parameterized
build/callback.c: tools/generate.py tools/config.py tests/callback/config.json tests/callback/species.json tests/callback/reactions.json | build
	.conda/bin/python tools/generate.py tests/callback/config.json $@
build/test_callback: tests/test_callback.c build/callback.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_callback.c build/callback.c build/libchem.a -lm -o $@
check: build/test_callback

build/test_jacobian_check: tests/test_jacobian_check.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_jacobian_check

build/process_fixture.c: tools/generate.py tools/config.py tests/process_fixture/config.json tests/process_fixture/species.json tests/process_fixture/reactions.json | build
	.conda/bin/python tools/generate.py tests/process_fixture/config.json $@
build/test_process_set: tests/test_process_set.c build/process_fixture.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_process_set.c build/process_fixture.c build/libchem.a -lm -o $@
check: build/test_process_set

build/test_state: tests/test_state.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_state

build/test_chain_rates: tests/test_chain_rates.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_chain_rates

build/test_process_random: tests/test_process_random.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_process_random

build/test_nonlinear: tests/test_nonlinear.c build/robertson.c build/libchem.a
	$(CC) $(CFLAGS) tests/test_nonlinear.c build/robertson.c build/libchem.a -lm -o $@
check: build/test_nonlinear

build/test_external_constraints: tests/test_external_constraints.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_external_constraints

build/test_external_composition: tests/test_external_composition.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_external_composition

build/test_matrix: tests/test_matrix.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_matrix

build/test_matrix_views: tests/test_matrix_views.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_matrix_views

build/test_matrix_function: tests/test_matrix_function.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_matrix_function

build/test_sparse_matrix: tests/test_sparse_matrix.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_sparse_matrix

build/test_sparse_views: tests/test_sparse_views.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_sparse_views

build/test_sparse_functions: tests/test_sparse_functions.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_sparse_functions

build/test_grouped_views: tests/test_grouped_views.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_grouped_views

build/test_dae_overshoot: tests/test_dae_overshoot.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_dae_overshoot

build/test_aerosol: tests/test_aerosol.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_aerosol
build/test_rate_store: tests/test_rate_store.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_rate_store
build/test_chapman_regression: tests/test_chapman_regression.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_chapman_regression
build/test_grouped_constraints: tests/test_grouped_constraints.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_grouped_constraints
build/test_aerosol_state: tests/test_aerosol_state.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_aerosol_state
build/test_species: tests/test_species.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_species
build/test_flow_tube: tests/test_flow_tube.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_flow_tube
build/test_process_builder: tests/test_process_builder.c build/libchem.a
	$(CC) $(CFLAGS) $< build/libchem.a -lm -o $@
check: build/test_process_builder
build/test_openmp: examples/openmp_tutorial.c build/libchem.a
	$(CC) $(CFLAGS) $(OMP_CFLAGS) -DTEST_MODE $< build/libchem.a $(OMP_LIBS) -lm -o $@
build/openmp_tutorial: examples/openmp_tutorial.c build/libchem.a
	$(CC) $(CFLAGS) $(OMP_CFLAGS) $< build/libchem.a $(OMP_LIBS) -lm -o $@
check: build/test_openmp build/openmp_tutorial
