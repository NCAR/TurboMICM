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
OMP_PREFIX ?= $(shell brew --prefix libomp)
OMPFLAGS = -Xpreprocessor -fopenmp -isystem $(OMP_PREFIX)/include
OMPLIBS = -L$(OMP_PREFIX)/lib -lomp
else
OMPFLAGS = -fopenmp
OMPLIBS = -fopenmp
endif
else
OMPFLAGS = -fopenmp
OMPLIBS = -fopenmp
endif
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic $(OPTFLAGS) $(OMPFLAGS) -Iinclude -Isrc
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

SOURCES = src/builder.c src/external.c src/sparse_matrix.c src/version.c src/system.c src/matrix.c src/rates.c src/mechanism.c src/state.c src/sparse.c src/rosenbrock.c src/rosenbrock_parameters.c src/backward_euler.c src/constraints.c src/constraint_initialization.c src/jacobian_verification.c
OBJECTS = $(SOURCES:src/%.c=build/%.o)
cpu: build/libchem.a
build/%.o: src/%.c include/chem.h include/chem_solver.h include/chem_constraints.h include/chem_jacobian.h src/internal.h include/chem_matrix.h include/chem_system.h include/chem_version.h include/chem_sparse_matrix.h include/chem_external.h include/chem_builder.h
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@
build/libchem.a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

bin/test_cpu: tests/main.c tests/test_builder.c tests/test_analytical.c tests/test_flowtube.c tests/test_chapman.c tests/test_aerosol.c tests/test_external.c tests/test_rates.c tests/test_sparse.c tests/test_backward_euler.c tests/test_constraints.c tests/test_constraint_initialization.c tests/test_dae.c tests/test_rosenbrock.c tests/test_jacobian_verification.c tests/test_process.c tests/test_rate_store.c tests/test_state.c tests/test_matrix.c tests/test_matrix_functions.c tests/test_system.c tests/test_runtime.c tests/test_sparse_matrix.c tests/test_sparse_functions.c tests/test_grouped_views.c tests/test.h tests/chapman_oracle.h tests/analytical_reference.h build/libchem.a
	@mkdir -p bin
	$(CC) $(CFLAGS) tests/main.c tests/test_builder.c tests/test_analytical.c tests/test_flowtube.c tests/test_chapman.c tests/test_aerosol.c tests/test_external.c tests/test_rates.c tests/test_sparse.c tests/test_backward_euler.c tests/test_constraints.c tests/test_constraint_initialization.c tests/test_dae.c tests/test_rosenbrock.c tests/test_jacobian_verification.c tests/test_process.c tests/test_rate_store.c tests/test_state.c tests/test_matrix.c tests/test_matrix_functions.c tests/test_system.c tests/test_runtime.c tests/test_sparse_matrix.c tests/test_sparse_functions.c tests/test_grouped_views.c build/libchem.a $(LDLIBS) -o $@
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
	bin/test_config_extra
	$(PYTHON) tools/test_configuration.py

$(OBJECTS): src/solver_internal.h

bin/tutorial_cpu: tests/tutorials.c tests/tutorial_mechanisms.h build/libchem.a build/compiler-settings
	$(CC) $(CFLAGS) tests/tutorials.c build/libchem.a $(LDLIBS) -o $@
bin/test_config: tests/tutorial_mechanisms.h
test-cpu: bin/tutorial_cpu

build/config_extra.c: mechanisms/config_extra.json tools/configure.py | env
	$(PYTHON) tools/configure.py $< $@ --symbol extra
build/config_extra.o: build/config_extra.c include/chem.h build/compiler-settings
	$(CC) $(CFLAGS) -c $< -o $@
bin/test_config_extra: tests/config_extra.c tests/test.h build/config_extra.o build/libchem.a build/compiler-settings
	$(CC) $(CFLAGS) tests/config_extra.c build/config_extra.o build/libchem.a $(LDLIBS) -o $@
test-cpu test-config: bin/test_config_extra

# Mechanism-specialized build products; never checked in or hand edited.
GROUP_SIZE ?= 32
GEN_CHUNK ?= 32
GEN_OPTIONS ?=
SPECIALIZED_NAMES = ts1 no_user with_user robertson extra
SPECIALIZED_OBJECTS = $(SPECIALIZED_NAMES:%=build/specialized/%.o)
build/specialization-settings: FORCE
	@mkdir -p build
	@printf '%s\n' '$(GROUP_SIZE)' '$(GEN_CHUNK)' '$(GEN_OPTIONS)' > $@.tmp
	@cmp -s $@.tmp $@ || cp $@.tmp $@
	@rm $@.tmp
build/specialized/ts1.c: mechanisms/ts1.json tools/configure.py tools/specialize.py build/specialization-settings | env
	$(PYTHON) tools/specialize.py $< $@ --symbol ts1 --chunk $(GEN_CHUNK) $(GEN_OPTIONS)
build/specialized/no_user.c: mechanisms/tutorial/rate_constants_no_user_defined/config.yaml mechanisms/tutorial/rate_constants_no_user_defined/species.yaml mechanisms/tutorial/rate_constants_no_user_defined/reactions.yaml tools/configure.py tools/specialize.py build/specialization-settings | env
	$(PYTHON) tools/specialize.py $< $@ --symbol no_user --chunk $(GEN_CHUNK) $(GEN_OPTIONS)
build/specialized/with_user.c: mechanisms/tutorial/rate_constants_user_defined/config.json mechanisms/tutorial/rate_constants_user_defined/species.json mechanisms/tutorial/rate_constants_user_defined/reactions.json tools/configure.py tools/specialize.py build/specialization-settings | env
	$(PYTHON) tools/specialize.py $< $@ --symbol with_user --chunk $(GEN_CHUNK) $(GEN_OPTIONS)
build/specialized/robertson.c: mechanisms/tutorial/robertson/config.json mechanisms/tutorial/robertson/species.json mechanisms/tutorial/robertson/reactions.json tools/configure.py tools/specialize.py build/specialization-settings | env
	$(PYTHON) tools/specialize.py $< $@ --symbol robertson --chunk $(GEN_CHUNK) $(GEN_OPTIONS)
build/specialized/extra.c: mechanisms/config_extra.json tools/configure.py tools/specialize.py build/specialization-settings | env
	$(PYTHON) tools/specialize.py $< $@ --symbol extra --chunk $(GEN_CHUNK) $(GEN_OPTIONS)
build/specialized/%.o: build/specialized/%.c src/rosenbrock_step.inc src/solver_internal.h src/internal.h include/chem.h include/chem_solver.h build/compiler-settings
	$(CC) $(CFLAGS) -c $< -o $@
$(SPECIALIZED_OBJECTS): build/compiler-settings
build/rosenbrock.o: src/rosenbrock_step.inc
bin/ts1_bench_cpu: src/ts1_driver.h
bin/ts1_bench_specialized: src/ts1_bench.c src/ts1_driver.h build/specialized/ts1.o build/libchem.a build/compiler-settings build/specialization-settings
	@mkdir -p bin
	$(CC) $(CFLAGS) -DCHEM_SPECIALIZED -DCHEM_GROUP_SIZE=$(GROUP_SIZE) src/ts1_bench.c build/specialized/ts1.o build/libchem.a $(LDLIBS) -o $@
bin/ts1_scenario_specialized: src/ts1_scenario.c src/ts1_driver.h build/specialized/ts1.o build/libchem.a build/compiler-settings build/specialization-settings
	@mkdir -p bin
	$(CC) $(CFLAGS) -DCHEM_SPECIALIZED -DCHEM_GROUP_SIZE=$(GROUP_SIZE) src/ts1_scenario.c build/specialized/ts1.o build/libchem.a $(LDLIBS) -o $@
bench: bin/ts1_bench_specialized bin/ts1_scenario_specialized
bin/test_specialized: tests/test_specialized.c $(SPECIALIZED_OBJECTS) build/libchem.a build/compiler-settings
	@mkdir -p bin
	$(CC) $(CFLAGS) tests/test_specialized.c $(SPECIALIZED_OBJECTS) build/libchem.a $(LDLIBS) -o $@
.PHONY: test-specialized
test-specialized: bin/test_specialized
	bin/test_specialized
