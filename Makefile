ifeq ($(origin CC),default)
CC = clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON = $(CURDIR)/.conda/bin/python
UNAME := $(shell uname -s)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -1)
ifeq ($(UNAME),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
ifneq (,$(findstring clang,$(CC_VERSION)))
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
CPPFLAGS += -Iinclude -Igenerated
CFLAGS = $(OPTFLAGS) -std=c11 -Wall -Wextra -Wpedantic $(OMPFLAGS)
LDLIBS = $(OMPLIBS) -lm
.PHONY: env cpu bench test-cpu micm-ref clean FORCE
env: .conda/bin/python
.conda/bin/python: environment.yml
	@if test -x $@; then CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env update --prefix $(CURDIR)/.conda --file environment.yml; else CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml; fi
	@touch $@
cpu: bin/ts1_bench_cpu
bench: cpu
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@
clean:
	rm -rf build bin generated

LIB_SOURCES = src/rates.c src/model.c src/state.c src/process.c src/parameters.c src/solver.c src/constraints.c src/backward_euler.c src/matrix.c src/system.c src/external.c src/verify.c
LIB_OBJECTS = $(LIB_SOURCES:src/%.c=build/%.o)
build/%.o: src/%.c $(wildcard include/*.h) src/internal.h build/configuration.txt
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
build/libchem.a: $(LIB_OBJECTS)
	ar rcs $@ $^

generated/ts1.c: mechanisms/ts1.json tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py $< $@ --symbol ts1
build/ts1.o: generated/ts1.c include/chem.h build/configuration.txt
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

GROUP_SIZE ?= 128
bin/ts1_bench_cpu: src/bench.c build/ts1.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCK_DEFAULT_GROUP_SIZE=$(GROUP_SIZE) $< build/ts1.o build/libchem.a $(LDLIBS) -o $@

TUTORIAL_OBJECTS = build/tutorial_no.o build/tutorial_user.o build/tutorial_robertson.o
TEST_SOURCES = tests/runner.c tests/test_rates.c tests/test_linear.c tests/test_constraints.c tests/test_integrators.c tests/test_initialization.c tests/tutorials.c tests/test_configuration.c tests/test_matrices.c tests/test_sparse.c tests/test_state.c tests/test_system.c tests/test_constraint_set.c tests/test_builder.c tests/test_verification.c
bin/test_cpu: $(TEST_SOURCES) tests/test.h $(TUTORIAL_OBJECTS) build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -Itests -Isrc $(TEST_SOURCES) $(TUTORIAL_OBJECTS) build/libchem.a $(LDLIBS) -o $@

test-cpu: bin/test_cpu
	@bin/test_cpu tests/cases.txt

generated/tutorial_no.c: mechanisms/tutorial/rate_constants_no_user_defined/config.yaml mechanisms/tutorial/rate_constants_no_user_defined/species.yaml mechanisms/tutorial/rate_constants_no_user_defined/reactions.yaml tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py $< $@ --symbol tutorial_no
generated/tutorial_user.c: mechanisms/tutorial/rate_constants_user_defined/config.json mechanisms/tutorial/rate_constants_user_defined/species.json mechanisms/tutorial/rate_constants_user_defined/reactions.json tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py $< $@ --symbol tutorial_user
generated/tutorial_robertson.c: mechanisms/tutorial/robertson/config.json mechanisms/tutorial/robertson/species.json mechanisms/tutorial/robertson/reactions.json tools/configure.py .conda/bin/python
	$(PYTHON) tools/configure.py $< $@ --symbol tutorial_robertson
build/tutorial_%.o: generated/tutorial_%.c include/chem.h build/configuration.txt
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
bin/tutorial_cpu: tests/tutorial_main.c tests/tutorials.c tests/test_configuration.c $(TUTORIAL_OBJECTS) build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -Itests -Isrc tests/tutorial_main.c tests/tutorials.c tests/test_configuration.c $(TUTORIAL_OBJECTS) build/libchem.a $(LDLIBS) -o $@

# Switching compilers, optimization flags, or tuning parameters rebuilds all C.
FORCE:
build/configuration.txt: FORCE
	@mkdir -p build
	@printf '%s\n' '$(CC)' '$(OPTFLAGS)' '$(CPPFLAGS)' '$(OMPFLAGS)' '$(OMPLIBS)' '$(GROUP_SIZE)' > build/configuration.new
	@if cmp -s build/configuration.new $@; then rm build/configuration.new; else mv build/configuration.new $@; fi
