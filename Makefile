ifeq ($(origin CC),default)
CC := clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
PYTHON ?= ./.conda/bin/python
ROOT := $(abspath ..)
BUILD := build
BIN := bin
GEN := $(BUILD)/generated
LIBSRC := src/chem.c src/state.c
LIBOBJ := $(patsubst src/%.c,$(BUILD)/%.o,$(LIBSRC))
LIBHDR := include/chem.h include/jacobian_verification.h
CONFIG := config/ts1.json
GENSRC := $(GEN)/ts1.c
NO_USER_CONFIG := config/rate_constants_no_user_defined/config.yaml
USER_CONFIG := config/rate_constants_user_defined/config.json
NO_USER_GEN := $(GEN)/rate_constants_no_user_defined.c
USER_GEN := $(GEN)/rate_constants_user_defined.c

UNAME_S := $(shell uname -s)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -1)
ifeq ($(findstring clang,$(CC_VERSION)),clang)
  ifeq ($(UNAME_S),Darwin)
    OPENMP_CFLAGS := -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include
    OPENMP_LDFLAGS := -L/opt/homebrew/opt/libomp/lib -lomp
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

.PHONY: all env cpu test-cpu bench micm-ref clean
all: cpu

env:
	CONDA_PKGS_DIRS=$(CURDIR)/.conda_pkgs conda env create --prefix $(CURDIR)/.conda -f environment.yml

$(GENSRC): $(CONFIG) tools/generate.py
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

test-cpu: cpu $(BIN)/test_rates $(BIN)/test_example_no_user $(BIN)/test_example_user $(BIN)/test_tutorial_programs $(BIN)/test_jacobian_verification $(BIN)/test_solver_configuration $(BIN)/test_integration_equivalents $(BIN)/test_state
	$(PYTHON) tools/test_cases.py

bench: cpu

micm-ref:
	@mkdir -p $(BIN)
	SDKROOT=$$(xcrun --sdk macosx --show-sdk-path) clang++ -std=c++20 $(OPTFLAGS) \
	  -I $(ROOT)/micm/include -I $(ROOT)/micm/test/benchmark \
	  $(ROOT)/reference/micm_ts1_ref.cpp -o $(BIN)/micm_ts1_ref

clean:
	rm -rf $(BUILD) $(BIN)
