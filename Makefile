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
    OMP_PREFIX ?= /Users/dfillmor/homebrew/opt/libomp
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
TS1_SPECIALIZED_MODEL := $(BUILD)/ts1_model_specialized.c
CORE_OBJECTS := $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/constraint.o $(BUILD)/solver.o
TS1_OBJECTS := $(CORE_OBJECTS) $(BUILD)/ts1_model.o $(BUILD)/benchmark.o
TS1_SPECIALIZED_OBJECTS := $(CORE_OBJECTS) $(BUILD)/ts1_model_specialized.o $(BUILD)/benchmark.o

.PHONY: all cpu bench bench-ablation test-cpu test-specialized env clean micm-ref mechanism
all: cpu bench

cpu: $(LIB)/libchem.a $(BIN)/ts1_bench_cpu

$(BUILD) $(BIN) $(LIB):
	mkdir -p $@

$(BUILD)/ts1_model.c: $(TS1_JSON) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@

$(TS1_SPECIALIZED_MODEL): $(TS1_JSON) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized

CONFIG ?= $(TS1_JSON)
MECHANISM_C := $(BUILD)/mechanism_model.c
MECHANISM_O := $(BUILD)/mechanism_model.o

mechanism: $(LIB)/libmechanism.a

$(MECHANISM_C): $(CONFIG) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@

$(MECHANISM_O): $(MECHANISM_C) include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(LIB)/libmechanism.a: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/state.o $(BUILD)/species.o $(BUILD)/system.o $(BUILD)/constraint.o $(BUILD)/backward_euler.o $(BUILD)/solver.o $(MECHANISM_O) | $(LIB)
	$(AR) rcs $@ $^

$(BUILD)/kinetics.o: src/kinetics.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/matrix.o: src/matrix.c include/matrix.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/state.o: src/state.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/species.o: src/species.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/system.o: src/system.c include/system.h include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/constraint.o: src/constraint.c include/constraint.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/backward_euler.o: src/backward_euler.c include/kinetics.h include/matrix.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/solver.o: src/solver.c include/kinetics.h src/solver_internal.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/ts1_model.o: $(TS1_MODEL) include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/benchmark.o: src/benchmark.c include/kinetics.h src/solver_internal.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(LIB)/libchem.a: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/state.o $(BUILD)/species.o $(BUILD)/system.o $(BUILD)/constraint.o $(BUILD)/backward_euler.o $(BUILD)/solver.o | $(LIB)
	$(AR) rcs $@ $^

bench: $(BIN)/ts1_bench_cpu $(BIN)/ts1_bench_specialized $(BIN)/ts1_scenario_specialized

$(BIN)/ts1_bench_cpu: $(TS1_OBJECTS) | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BUILD)/ts1_model_specialized.o: $(TS1_SPECIALIZED_MODEL) include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BIN)/ts1_bench_specialized: $(TS1_SPECIALIZED_OBJECTS) | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BUILD)/scenario.o: src/scenario.c include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BIN)/ts1_scenario_specialized: $(CORE_OBJECTS) $(BUILD)/ts1_model_specialized.o $(BUILD)/scenario.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BUILD)/ts1_ablate_rates.c: $(TS1_JSON) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized --without-rate-kernel

$(BUILD)/ts1_ablate_chemistry.c: $(TS1_JSON) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized --without-chemistry-kernel

$(BUILD)/ts1_ablate_factor.c: $(TS1_JSON) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized --without-factor-kernel

$(BUILD)/ts1_ablate_solve.c: $(TS1_JSON) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized --without-solve-kernel

$(BUILD)/ts1_ablate_%.o: $(BUILD)/ts1_ablate_%.c include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BIN)/ts1_bench_ablate_rates: $(CORE_OBJECTS) $(BUILD)/benchmark.o $(BUILD)/ts1_ablate_rates.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BIN)/ts1_bench_ablate_chemistry: $(CORE_OBJECTS) $(BUILD)/benchmark.o $(BUILD)/ts1_ablate_chemistry.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BIN)/ts1_bench_ablate_factor: $(CORE_OBJECTS) $(BUILD)/benchmark.o $(BUILD)/ts1_ablate_factor.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BIN)/ts1_bench_ablate_solve: $(CORE_OBJECTS) $(BUILD)/benchmark.o $(BUILD)/ts1_ablate_solve.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

bench-ablation: $(BIN)/ts1_bench_ablate_rates $(BIN)/ts1_bench_ablate_chemistry $(BIN)/ts1_bench_ablate_factor $(BIN)/ts1_bench_ablate_solve

SDKROOT_FLAG := $(if $(filter Darwin,$(UNAME_S)),-isysroot $(SDKROOT),)

env:
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs conda env create --prefix $(CURDIR)/.conda -f environment.yml

test-cpu: $(BIN)/test_cpu $(BIN)/test_programs $(BIN)/ts1_bench_cpu $(BIN)/test_example_no_user $(BIN)/test_example_user $(BUILD)/tutorial_no_user_model.o $(BUILD)/tutorial_user_model.o $(BUILD)/robertson_model.o
	$(PYTHON) tests/run_cases.py

$(BUILD)/test_specialized.o: tests/test_specialized.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BUILD)/specialized_ts1.c: $(TS1_JSON) tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized

$(BUILD)/specialized_robertson.c: examples/robertson/config.json tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized

$(BUILD)/specialized_tutorial_user.c: examples/rate_constants_user_defined/config.json tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized

$(BUILD)/specialized_tutorial_no_user.c: examples/rate_constants_no_user_defined/config.yaml tools/generate_c.py | $(BUILD)
	$(PYTHON) tools/generate_c.py $< $@ --specialized

$(BUILD)/specialized_%.o: $(BUILD)/specialized_%.c include/kinetics.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BIN)/test_specialized_ts1: $(CORE_OBJECTS) $(BUILD)/test_specialized.o $(BUILD)/specialized_ts1.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BIN)/test_specialized_robertson: $(CORE_OBJECTS) $(BUILD)/test_specialized.o $(BUILD)/specialized_robertson.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BIN)/test_specialized_tutorial_user: $(CORE_OBJECTS) $(BUILD)/test_specialized.o $(BUILD)/specialized_tutorial_user.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BIN)/test_specialized_tutorial_no_user: $(CORE_OBJECTS) $(BUILD)/test_specialized.o $(BUILD)/specialized_tutorial_no_user.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

test-specialized: $(BIN)/test_specialized_ts1 $(BIN)/test_specialized_robertson $(BIN)/test_specialized_tutorial_user $(BIN)/test_specialized_tutorial_no_user
	$(BIN)/test_specialized_ts1
	$(BIN)/test_specialized_robertson
	$(BIN)/test_specialized_tutorial_user
	$(BIN)/test_specialized_tutorial_no_user

$(BIN)/test_cpu: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/state.o $(BUILD)/species.o $(BUILD)/system.o $(BUILD)/constraint.o $(BUILD)/backward_euler.o $(BUILD)/solver.o $(BUILD)/ts1_model.o $(BUILD)/test_cpu.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BUILD)/test_programs.o: tests/test_programs.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

$(BIN)/test_programs: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/constraint.o $(BUILD)/solver.o $(BUILD)/test_programs.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BUILD)/example_no_user.o: tests/test_example_runtime.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -DEXAMPLE_USER_CONFIG=0 -c $< -o $@

$(BUILD)/example_user.o: tests/test_example_runtime.c include/kinetics.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) -DEXAMPLE_USER_CONFIG=1 -c $< -o $@

$(BIN)/test_example_no_user: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/constraint.o $(BUILD)/solver.o $(BUILD)/tutorial_no_user_model.o $(BUILD)/example_no_user.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BIN)/test_example_user: $(BUILD)/kinetics.o $(BUILD)/matrix.o $(BUILD)/constraint.o $(BUILD)/solver.o $(BUILD)/tutorial_user_model.o $(BUILD)/example_user.o | $(BIN)
	$(CC) $(OPTFLAGS) $(SDKROOT_FLAG) $^ $(LDLIBS) -o $@

$(BUILD)/test_cpu.o: tests/test_cpu.c include/kinetics.h include/matrix.h include/system.h include/constraint.h include/version.h | $(BUILD)
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
