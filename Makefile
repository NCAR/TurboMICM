ifeq ($(origin CC),default)
CC := clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
PYTHON ?= ./.conda/bin/python
ROOT := $(abspath ..)
BUILD := build
BIN := bin
GEN := $(BUILD)/generated
LIBSRC := src/chem.c src/state.c src/objects.c src/constraints.c src/linalg.c src/views.c src/matrix.c src/sparse_matrix.c src/backward_euler.c
LIBOBJ := $(patsubst src/%.c,$(BUILD)/%.o,$(LIBSRC))
LIBHDR := include/chem.h include/chem_constraints.h include/chem_views.h include/chem_matrix.h include/chem_sparse_matrix.h include/jacobian_verification.h
LIBHDR += include/chem_specialized.h
CONFIG := config/ts1.json
GENSRC := $(GEN)/ts1.c
NO_USER_CONFIG := config/rate_constants_no_user_defined/config.yaml
USER_CONFIG := config/rate_constants_user_defined/config.json
NO_USER_GEN := $(GEN)/rate_constants_no_user_defined.c
USER_GEN := $(GEN)/rate_constants_user_defined.c
LIBOMP_PREFIX ?= /Users/dfillmor/homebrew/opt/libomp

UNAME_S := $(shell uname -s)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -1)
ifeq ($(findstring clang,$(CC_VERSION)),clang)
  ifeq ($(UNAME_S),Darwin)
    OPENMP_CFLAGS := -Xpreprocessor -fopenmp -isystem $(LIBOMP_PREFIX)/include
    OPENMP_LDFLAGS := -L$(LIBOMP_PREFIX)/lib -lomp
  else
    OPENMP_CFLAGS := -fopenmp
    OPENMP_LDFLAGS := -fopenmp
  endif
else
  OPENMP_CFLAGS := -fopenmp
  OPENMP_LDFLAGS := -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
  SDKROOT ?= $(shell xcrun --sdk macosx --show-sdk-path)
  export SDKROOT
endif
CFLAGS := -std=c11 $(OPTFLAGS) -Wall -Wextra -Wpedantic $(OPENMP_CFLAGS) -Iinclude
LDFLAGS := $(OPENMP_LDFLAGS) -lm

.PHONY: all env cpu test-cpu test-specialized bench micm-ref clean
all: cpu

env:
	CONDA_PKGS_DIRS=$(CURDIR)/.conda_pkgs conda env create --prefix $(CURDIR)/.conda -f environment.yml

$(GENSRC): $(CONFIG) tools/generate.py include/chem_specialized.h
	@mkdir -p $(GEN)
	$(PYTHON) tools/generate.py $(CONFIG) $@

cpu: $(BIN)/libchem.a $(BIN)/ts1_bench_cpu

$(BUILD)/%.o: src/%.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/ts1.o: src/ts1_bench.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c src/ts1_bench.c -o $@

$(BUILD)/ts1_mechanism.o: $(GENSRC) $(LIBHDR)
	$(CC) $(CFLAGS) -c $(GENSRC) -o $@

$(BIN)/libchem.a: $(LIBOBJ)
	@mkdir -p $(BIN)
	ar rcs $@ $^

$(BIN)/ts1_bench_cpu: $(BUILD)/ts1.o $(BUILD)/ts1_mechanism.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_rates.o: tests/test_rates.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_rates: $(BUILD)/test_rates.o $(BUILD)/ts1_mechanism.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(NO_USER_GEN): $(NO_USER_CONFIG) config/rate_constants_no_user_defined/species.yaml config/rate_constants_no_user_defined/reactions.yaml tools/generate.py
	@mkdir -p $(GEN)
	$(PYTHON) tools/generate.py $(NO_USER_CONFIG) $@

$(USER_GEN): $(USER_CONFIG) config/rate_constants_user_defined/species.json config/rate_constants_user_defined/reactions.json tools/generate.py
	@mkdir -p $(GEN)
	$(PYTHON) tools/generate.py $(USER_CONFIG) $@

$(BUILD)/test_example_no_user.o: tests/test_example_descriptor.c $(LIBHDR)
	@mkdir -p $(BUILD)
	$(CC) $(CFLAGS) -DEXPECT_USER_DEFINED=0 -c $< -o $@

$(BUILD)/test_example_user.o: tests/test_example_descriptor.c $(LIBHDR)
	@mkdir -p $(BUILD)
	$(CC) $(CFLAGS) -DEXPECT_USER_DEFINED=1 -c $< -o $@

$(BUILD)/no_user_mechanism.o: $(NO_USER_GEN) $(LIBHDR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/user_mechanism.o: $(USER_GEN) $(LIBHDR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_example_no_user: $(BUILD)/test_example_no_user.o $(BUILD)/no_user_mechanism.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BIN)/test_example_user: $(BUILD)/test_example_user.o $(BUILD)/user_mechanism.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_tutorial_programs.o: tests/test_tutorial_programs.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_tutorial_programs: $(BUILD)/test_tutorial_programs.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_jacobian_verification.o: tests/test_jacobian_verification.c include/jacobian_verification.h
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_jacobian_verification: $(BUILD)/test_jacobian_verification.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(OPENMP_LDFLAGS) -lm -o $@

$(BUILD)/test_solver_configuration.o: tests/test_solver_configuration.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_solver_configuration: $(BUILD)/test_solver_configuration.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_integration_equivalents.o: tests/test_integration_equivalents.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_integration_equivalents: $(BUILD)/test_integration_equivalents.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_state.o: tests/test_state.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_state: $(BUILD)/test_state.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_objects.o: tests/test_objects.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_objects: $(BUILD)/test_objects.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_processes.o: tests/test_processes.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_processes: $(BUILD)/test_processes.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_constraints.o: tests/test_constraints.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_constraints: $(BUILD)/test_constraints.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_lu.o: tests/test_lu.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_lu: $(BUILD)/test_lu.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_chapman.o: tests/test_chapman.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_chapman: $(BUILD)/test_chapman.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_grouped_views.o: tests/test_grouped_views.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_grouped_views: $(BUILD)/test_grouped_views.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_matrix.o: tests/test_matrix.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_matrix: $(BUILD)/test_matrix.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_sparse_matrix.o: tests/test_sparse_matrix.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_sparse_matrix: $(BUILD)/test_sparse_matrix.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_backward_euler.o: tests/test_backward_euler.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_backward_euler: $(BUILD)/test_backward_euler.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_analytical_simple.o: tests/test_analytical_simple.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_analytical_simple: $(BUILD)/test_analytical_simple.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_rosenbrock_utilities.o: tests/test_rosenbrock_utilities.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_rosenbrock_utilities: $(BUILD)/test_rosenbrock_utilities.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_analytical_additional.o: tests/test_analytical_additional.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_analytical_additional: $(BUILD)/test_analytical_additional.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_regression.o: tests/test_regression.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_regression: $(BUILD)/test_regression.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_aerosol_equivalents.o: tests/test_aerosol_equivalents.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_aerosol_equivalents: $(BUILD)/test_aerosol_equivalents.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/test_external_constraints.o: tests/test_external_constraints.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN)/test_external_constraints: $(BUILD)/test_external_constraints.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

test-cpu: cpu $(BIN)/test_rates $(BIN)/test_example_no_user $(BIN)/test_example_user $(BIN)/test_tutorial_programs $(BIN)/test_jacobian_verification $(BIN)/test_solver_configuration $(BIN)/test_integration_equivalents $(BIN)/test_state $(BIN)/test_objects $(BIN)/test_processes $(BIN)/test_constraints $(BIN)/test_lu $(BIN)/test_chapman $(BIN)/test_grouped_views $(BIN)/test_matrix $(BIN)/test_sparse_matrix $(BIN)/test_backward_euler $(BIN)/test_analytical_simple $(BIN)/test_rosenbrock_utilities $(BIN)/test_analytical_additional $(BIN)/test_regression $(BIN)/test_aerosol_equivalents $(BIN)/test_external_constraints
	$(PYTHON) tools/test_cases.py

$(BUILD)/ts1_specialized.o: src/ts1_bench.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -DCHEM_USE_SPECIALIZED -c src/ts1_bench.c -o $@

$(BIN)/ts1_bench_specialized: $(BUILD)/ts1_specialized.o $(BUILD)/ts1_mechanism.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BUILD)/ts1_scenario.o: src/ts1_bench.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -DCHEM_USE_SPECIALIZED -DCHEM_SCENARIO -c src/ts1_bench.c -o $@

$(BIN)/ts1_scenario_specialized: $(BUILD)/ts1_scenario.o $(BUILD)/ts1_mechanism.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

$(BIN)/test_specialized_ts1: tests/test_specialized.c $(BUILD)/ts1_mechanism.o $(BIN)/libchem.a $(LIBHDR)
	$(CC) $(CFLAGS) tests/test_specialized.c $(BUILD)/ts1_mechanism.o $(BIN)/libchem.a $(LDFLAGS) -o $@

$(BUILD)/specialized_no_user.o: $(NO_USER_GEN) $(LIBHDR)
	$(CC) $(CFLAGS) -c $(NO_USER_GEN) -o $@

$(BUILD)/specialized_user.o: $(USER_GEN) $(LIBHDR)
	$(CC) $(CFLAGS) -c $(USER_GEN) -o $@

$(BIN)/test_specialized_no_user: tests/test_specialized.c $(BUILD)/specialized_no_user.o $(BIN)/libchem.a $(LIBHDR)
	$(CC) $(CFLAGS) tests/test_specialized.c $(BUILD)/specialized_no_user.o $(BIN)/libchem.a $(LDFLAGS) -o $@

$(BIN)/test_specialized_user: tests/test_specialized.c $(BUILD)/specialized_user.o $(BIN)/libchem.a $(LIBHDR)
	$(CC) $(CFLAGS) tests/test_specialized.c $(BUILD)/specialized_user.o $(BIN)/libchem.a $(LDFLAGS) -o $@

test-specialized: cpu $(BIN)/test_specialized_ts1 $(BIN)/test_specialized_no_user $(BIN)/test_specialized_user
	$(BIN)/test_specialized_ts1
	$(BIN)/test_specialized_no_user
	$(BIN)/test_specialized_user

bench: cpu $(BIN)/ts1_bench_specialized $(BIN)/ts1_scenario_specialized

micm-ref:
	@mkdir -p $(BIN)
	SDKROOT=$$(xcrun --sdk macosx --show-sdk-path) clang++ -std=c++20 $(OPTFLAGS) \
	  -I $(ROOT)/micm/include -I $(ROOT)/micm/test/benchmark \
	  $(ROOT)/reference/micm_ts1_ref.cpp -o $(BIN)/micm_ts1_ref

clean:
	rm -rf $(BUILD) $(BIN)
