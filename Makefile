CC = clang
OPTFLAGS ?= -O3 -DNDEBUG
AR ?= ar
PYTHON = $(if $(wildcard $(CURDIR)/.conda/bin/python),$(CURDIR)/.conda/bin/python,python3)

UNAME_S := $(shell uname -s)
CC_BASENAME := $(notdir $(firstword $(CC)))
OPENMP_CFLAGS :=
OPENMP_LDFLAGS :=
ifeq ($(UNAME_S),Darwin)
  SDKROOT ?= $(shell xcrun --sdk macosx --show-sdk-path)
  ifneq (,$(findstring clang,$(CC_BASENAME)))
    OMP_PREFIX ?= /opt/homebrew/opt/libomp
    OPENMP_CFLAGS := -Xpreprocessor -fopenmp -isystem $(OMP_PREFIX)/include
    OPENMP_LDFLAGS := -L$(OMP_PREFIX)/lib -lomp
  else ifneq (,$(findstring gcc,$(CC_BASENAME)))
    OPENMP_CFLAGS := -fopenmp
    OPENMP_LDFLAGS := -fopenmp
  endif
else
  OPENMP_CFLAGS := -fopenmp
  OPENMP_LDFLAGS := -fopenmp
endif

CPPFLAGS := -Iinclude -Ibuild
CFLAGS := -std=c11 -Wall -Wextra -Wpedantic $(OPENMP_CFLAGS)
LDLIBS := -lm $(OPENMP_LDFLAGS)
BUILD := build
BIN := bin
LIB := lib
TS1_JSON := mechanisms/ts1.json
TS1_MODEL := $(BUILD)/ts1_model.c
TS1_OBJECTS := $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/solver.o $(BUILD)/ts1_model.o $(BUILD)/benchmark.o

.PHONY: all cpu bench test-cpu env clean micm-ref mechanism
all: cpu bench

cpu: $(LIB)/libchem.a $(BIN)/ts1_bench_cpu

$(BUILD) $(BIN) $(LIB):
	mkdir -p $@

$(BUILD)/ts1_model.c: $(TS1_JSON) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@

CONFIG ?= $(TS1_JSON)
MECHANISM_C := $(BUILD)/mechanism_model.c
MECHANISM_O := $(BUILD)/mechanism_model.o

mechanism: $(LIB)/libmechanism.a

$(MECHANISM_C): $(CONFIG) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@

$(MECHANISM_O): $(MECHANISM_C) include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(LIB)/libmechanism.a: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/state.o $(BUILD)/solver.o $(MECHANISM_O) | $(LIB)
	$(AR) rcs $@ $^

$(BUILD)/kinetics.o: src/kinetics.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/matrix.o: src/matrix.c include/matrix.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/state.o: src/state.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/solver.o: src/solver.c include/kinetics.h src/solver_internal.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/ts1_model.o: $(TS1_MODEL) include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/benchmark.o: src/benchmark.c include/kinetics.h src/solver_internal.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(LIB)/libchem.a: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/state.o $(BUILD)/solver.o | $(LIB)
	$(AR) rcs $@ $^

bench: $(BIN)/ts1_bench_cpu

$(BIN)/ts1_bench_cpu: $(TS1_OBJECTS) | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

SDKROOT_FLAG := $(if $(filter Darwin,$(UNAME_S)),-isysroot $(SDKROOT),)

env:
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs conda env create --prefix $(CURDIR)/.conda -f environment.yml

test-cpu: $(BIN)/test_cpu $(BIN)/ts1_bench_cpu $(BIN)/test_example_no_user $(BIN)/test_example_user $(BUILD)/tutorial_no_user_model.o $(BUILD)/tutorial_user_model.o $(BUILD)/robertson_model.o
	$(PYTHON) tests/run_cases.py

$(BIN)/test_cpu: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/state.o $(BUILD)/solver.o $(BUILD)/ts1_model.o $(BUILD)/test_cpu.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BUILD)/example_no_user.o: tests/test_example_runtime.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -DEXAMPLE_USER_CONFIG=0 -c $< -o $@

$(BUILD)/example_user.o: tests/test_example_runtime.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -DEXAMPLE_USER_CONFIG=1 -c $< -o $@

$(BIN)/test_example_no_user: $(BUILD)/kinetics.o $(BUILD)/solver.o $(BUILD)/tutorial_no_user_model.o $(BUILD)/example_no_user.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BIN)/test_example_user: $(BUILD)/kinetics.o $(BUILD)/solver.o $(BUILD)/tutorial_user_model.o $(BUILD)/example_user.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BUILD)/test_cpu.o: tests/test_cpu.c include/kinetics.h include/matrix.h include/version.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/tutorial_no_user_model.c: examples/rate_constants_no_user_defined/config.yaml tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@

$(BUILD)/tutorial_user_model.c: examples/rate_constants_user_defined/config.json tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@

$(BUILD)/robertson_model.c: examples/robertson/config.json tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@

$(BUILD)/tutorial_no_user_model.o: $(BUILD)/tutorial_no_user_model.c include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/tutorial_user_model.o: $(BUILD)/tutorial_user_model.c include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/robertson_model.o: $(BUILD)/robertson_model.c include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

micm-ref: | $(BIN)
	SDKROOT=$$(xcrun --sdk macosx --show-sdk-path) clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark ../reference/micm_ts1_ref.cpp -o $(BIN)/micm_ts1_ref

clean:
	rm -rf $(BUILD) $(BIN) $(LIB)
