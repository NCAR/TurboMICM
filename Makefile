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
OMP_PREFIX ?= /Users/dfillmor/homebrew/opt/libomp
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
CPPFLAGS += -Iinclude -Igenerated -Isrc
CFLAGS = $(OPTFLAGS) -std=c11 -Wall -Wextra -Wpedantic $(OMPFLAGS)
LDLIBS = $(OMPLIBS) -lm
.PHONY: env cpu bench test-cpu micm-ref clean FORCE
env: .conda/bin/python
.conda/bin/python: environment.yml
	@if test -x $@; then CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env update --prefix $(CURDIR)/.conda --file environment.yml; else CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml; fi
	@touch $@
cpu: bin/ts1_bench_cpu
bench: cpu bin/ts1_bench_specialized bin/ts1_scenario_specialized
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp build/configuration.txt
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@
clean:
	rm -rf build bin generated

LIB_SOURCES = src/specialized.c src/rates.c src/model.c src/state.c src/process.c src/parameters.c src/solver.c src/constraints.c src/backward_euler.c src/matrix.c src/system.c src/external.c src/verify.c
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

GROUP_SIZE ?= 32
bin/ts1_bench_cpu: src/bench.c src/bench_common.h build/ts1.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCK_DEFAULT_GROUP_SIZE=$(GROUP_SIZE) $< build/ts1.o build/libchem.a $(LDLIBS) -o $@

TUTORIAL_OBJECTS = build/tutorial_no.o build/tutorial_user.o build/tutorial_robertson.o
TEST_SOURCES = tests/runner.c tests/test_rates.c tests/test_linear.c tests/test_constraints.c tests/test_integrators.c tests/test_initialization.c tests/tutorials.c tests/test_configuration.c tests/test_matrices.c tests/test_sparse.c tests/test_state.c tests/test_system.c tests/test_constraint_set.c tests/test_builder.c tests/test_verification.c tests/test_analytical.c tests/test_process.c tests/test_process_set.c tests/test_grouped.c tests/test_aerosol.c tests/test_external.c tests/test_dae.c tests/test_chemistry.c tests/test_chapman.c tests/test_stiff.c
bin/test_cpu: $(TEST_SOURCES) $(wildcard tests/*.h) $(TUTORIAL_OBJECTS) build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -Itests -Isrc $(TEST_SOURCES) $(TUTORIAL_OBJECTS) build/libchem.a $(LDLIBS) -o $@

test-cpu: bin/test_cpu
	@$(PYTHON) tools/test_configuration.py
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

# Mechanism specialization: generated C is always a disposable build product.
SPECIALIZE_BATCH ?= 8
SPECIALIZE_SPLIT ?= 32
SPECIALIZE_ALGEBRA_BATCH ?= 8
SPECIALIZE_WIDTH ?= $(GROUP_SIZE)
SPECIALIZE_TILE ?= 0
SPECIALIZE_FUSE ?= 0
SPECIALIZE_JOINT_BATCH ?= 1
SPECIALIZE_LU_ORDER ?= column
SPECIALIZED_OBJECTS = build/sp_ts1.o build/sp_tutorial_no.o build/sp_tutorial_user.o build/sp_tutorial_robertson.o
SPECIALIZE_DEPS = tools/specialize.py tools/configure.py .conda/bin/python build/specialization.txt
build/specialization.txt: FORCE
	@mkdir -p build
	@printf '%s\n' '$(SPECIALIZE_BATCH)' '$(SPECIALIZE_SPLIT)' '$(SPECIALIZE_ALGEBRA_BATCH)' '$(SPECIALIZE_WIDTH)' '$(SPECIALIZE_TILE)' '$(SPECIALIZE_FUSE)' '$(SPECIALIZE_JOINT_BATCH)' '$(SPECIALIZE_LU_ORDER)' > build/specialization.new
	@if cmp -s build/specialization.new $@; then rm build/specialization.new; else mv build/specialization.new $@; fi

generated/sp_ts1.c: mechanisms/ts1.json $(SPECIALIZE_DEPS)
	$(PYTHON) tools/specialize.py $< $@ --symbol ts1_specialized --batch $(SPECIALIZE_BATCH) --split $(SPECIALIZE_SPLIT) --algebra-batch $(SPECIALIZE_ALGEBRA_BATCH) --width $(SPECIALIZE_WIDTH) --tile $(SPECIALIZE_TILE) --fuse $(SPECIALIZE_FUSE) --joint-batch $(SPECIALIZE_JOINT_BATCH) --lu-order $(SPECIALIZE_LU_ORDER)
generated/sp_tutorial_no.c: $(wildcard mechanisms/tutorial/rate_constants_no_user_defined/*) $(SPECIALIZE_DEPS)
	$(PYTHON) tools/specialize.py mechanisms/tutorial/rate_constants_no_user_defined/config.yaml $@ --symbol sp_tutorial_no --batch $(SPECIALIZE_BATCH) --split $(SPECIALIZE_SPLIT) --algebra-batch $(SPECIALIZE_ALGEBRA_BATCH) --width $(SPECIALIZE_WIDTH) --tile $(SPECIALIZE_TILE) --fuse $(SPECIALIZE_FUSE) --joint-batch $(SPECIALIZE_JOINT_BATCH) --lu-order $(SPECIALIZE_LU_ORDER)
generated/sp_tutorial_user.c: $(wildcard mechanisms/tutorial/rate_constants_user_defined/*) $(SPECIALIZE_DEPS)
	$(PYTHON) tools/specialize.py mechanisms/tutorial/rate_constants_user_defined/config.json $@ --symbol sp_tutorial_user --batch $(SPECIALIZE_BATCH) --split $(SPECIALIZE_SPLIT) --algebra-batch $(SPECIALIZE_ALGEBRA_BATCH) --width $(SPECIALIZE_WIDTH) --tile $(SPECIALIZE_TILE) --fuse $(SPECIALIZE_FUSE) --joint-batch $(SPECIALIZE_JOINT_BATCH) --lu-order $(SPECIALIZE_LU_ORDER)
generated/sp_tutorial_robertson.c: $(wildcard mechanisms/tutorial/robertson/*) $(SPECIALIZE_DEPS)
	$(PYTHON) tools/specialize.py mechanisms/tutorial/robertson/config.json $@ --symbol sp_tutorial_robertson --batch $(SPECIALIZE_BATCH) --split $(SPECIALIZE_SPLIT) --algebra-batch $(SPECIALIZE_ALGEBRA_BATCH) --width $(SPECIALIZE_WIDTH) --tile $(SPECIALIZE_TILE) --fuse $(SPECIALIZE_FUSE) --joint-batch $(SPECIALIZE_JOINT_BATCH) --lu-order $(SPECIALIZE_LU_ORDER)
build/sp_%.o: generated/sp_%.c $(wildcard include/*.h) src/internal.h src/specialized_solver.h build/configuration.txt
	@mkdir -p build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
bin/ts1_bench_specialized: src/bench.c src/bench_common.h build/sp_ts1.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCK_DEFAULT_GROUP_SIZE=$(GROUP_SIZE) -DCK_SPECIALIZED $< build/sp_ts1.o build/libchem.a $(LDLIBS) -o $@

bin/ts1_scenario_specialized: src/scenario.c src/bench_common.h build/sp_ts1.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCK_DEFAULT_GROUP_SIZE=$(GROUP_SIZE) -DCK_SPECIALIZED $< build/sp_ts1.o build/libchem.a $(LDLIBS) -o $@

generated/extended.json: tools/make_specialized_fixture.py .conda/bin/python
	$(PYTHON) $<
generated/sp_extended.c: generated/extended.json $(SPECIALIZE_DEPS)
	$(PYTHON) tools/specialize.py $< $@ --symbol sp_extended --batch $(SPECIALIZE_BATCH) --split $(SPECIALIZE_SPLIT) --algebra-batch $(SPECIALIZE_ALGEBRA_BATCH) --width $(SPECIALIZE_WIDTH) --tile $(SPECIALIZE_TILE) --fuse $(SPECIALIZE_FUSE) --joint-batch $(SPECIALIZE_JOINT_BATCH) --lu-order $(SPECIALIZE_LU_ORDER)
bin/test_specialized: tests/test_specialized.c $(SPECIALIZED_OBJECTS) build/sp_extended.o build/sp_emission.o build/sp_alternate.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(SPECIALIZED_OBJECTS) build/sp_extended.o build/sp_emission.o build/sp_alternate.o build/libchem.a $(LDLIBS) -o $@
.PHONY: test-specialized
test-specialized: bin/test_specialized
	$(PYTHON) tools/test_specialize.py
	bin/test_specialized

bin/bench_kernels: tools/bench_kernels.c build/sp_ts1.o build/libchem.a
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) $< build/sp_ts1.o build/libchem.a $(LDLIBS) -o $@

generated/emission.json: generated/extended.json
	@test -f $@ || $(PYTHON) tools/make_specialized_fixture.py
generated/sp_emission.c: generated/emission.json $(SPECIALIZE_DEPS)
	$(PYTHON) tools/specialize.py $< $@ --symbol sp_emission --batch $(SPECIALIZE_BATCH) --split $(SPECIALIZE_SPLIT) --algebra-batch $(SPECIALIZE_ALGEBRA_BATCH) --width $(SPECIALIZE_WIDTH) --tile $(SPECIALIZE_TILE) --fuse $(SPECIALIZE_FUSE) --joint-batch $(SPECIALIZE_JOINT_BATCH) --lu-order $(SPECIALIZE_LU_ORDER)

# Exercise optional transformations together, including odd fixed/tile widths.
generated/sp_alternate.c: generated/extended.json $(SPECIALIZE_DEPS)
	$(PYTHON) tools/specialize.py $< $@ --symbol sp_alternate --batch 64 --split 7 --algebra-batch 64 --width 3 --tile 2 --fuse 1 --joint-batch 4 --lu-order row
