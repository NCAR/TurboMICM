ROOT := $(abspath ..)
.SILENT:
.SUFFIXES:
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON := $(CURDIR)/.conda/bin/python
export TMPDIR := $(CURDIR)/build/tmp
OPTFLAGS ?= -O3 -DNDEBUG
GROUP_SIZE ?= 32
KERNEL_CHUNK ?= 32
LINEAR_CHUNK ?= 8
FUSE_KERNELS ?= 0
FUSE_STAGES ?= 1
SPECIALIZE_RATES ?= 1
LINEAR_MODE ?= unrolled
ifeq ($(origin CC),default)
CC := clang
endif
UNAME_S := $(shell uname -s)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -1)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
SHARED_EXT := dylib
SHARED_FLAG := -dynamiclib
else
SHARED_EXT := so
SHARED_FLAG := -shared
endif
ifneq (,$(findstring clang,$(CC_VERSION)))
ifeq ($(UNAME_S),Darwin)
LIBOMP_PREFIX ?= $(shell brew --prefix libomp)
OMPFLAGS := -Xpreprocessor -fopenmp -isystem $(LIBOMP_PREFIX)/include
OMPLIBS := -L$(LIBOMP_PREFIX)/lib -lomp
else
OMPFLAGS := -fopenmp
OMPLIBS := -fopenmp
endif
else
OMPFLAGS := -fopenmp
OMPLIBS := -fopenmp
endif
CFLAGS := -std=c11 -Wall -Wextra -Wpedantic $(OPTFLAGS) $(OMPFLAGS)
CPPFLAGS := -Iinclude -Igenerated -DCHEM_GROUP_SIZE=$(GROUP_SIZE)
GENFLAGS := --kernel-chunk $(KERNEL_CHUNK)
SPECFLAGS := $(GENFLAGS) --linear-chunk $(LINEAR_CHUNK) --fuse-kernels $(FUSE_KERNELS) --fuse-stages $(FUSE_STAGES) --specialize-rates $(SPECIALIZE_RATES) --linear $(LINEAR_MODE)
LDLIBS := $(OMPLIBS) -lm
BUILD := build/$(notdir $(CC))
SOURCES := $(filter-out src/bench.c src/scenario.c,$(wildcard src/*.c))
OBJECTS := $(patsubst src/%.c,$(BUILD)/%.o,$(SOURCES))
GENERATED_OBJECTS := $(BUILD)/ts1.o $(BUILD)/example_no_user.o $(BUILD)/example_user.o $(BUILD)/example_robertson.o $(BUILD)/validation.o
SPEC_NAMES := ts1 example_no_user example_user example_robertson validation edges ts1_variant
SPECIALIZED_OBJECTS := $(addprefix $(BUILD)/,$(addsuffix _specialized.o,$(SPEC_NAMES)))

.PHONY: all env cpu test-cpu test-specialized bench micm-ref clean
.SECONDARY: $(addprefix generated/,$(addsuffix _specialized.c,$(SPEC_NAMES)))
all: cpu
env: .conda/conda-meta/history
	mkdir -p $(TMPDIR)
.conda/conda-meta/history: environment.yml
	mkdir -p $(TMPDIR)
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes >&2
cpu: bin/libchem.a bin/libchem.$(SHARED_EXT) bin/ts1_bench_cpu bin/tutorials
bench: cpu bin/ts1_bench_specialized bin/ts1_scenario_specialized
kernel-bench: bin/kernel_bench
bin/kernel_bench: tests/kernel_bench.c $(BUILD)/ts1_specialized.o $(BUILD)/libchem.a $(BUILD)/flags
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(BUILD)/ts1_specialized.o $(BUILD)/libchem.a $(LDLIBS) -o $@
test-cpu: cpu
	$(PYTHON) tests/run.py
test-specialized: bench bin/libspecialized.$(SHARED_EXT)
	$(PYTHON) tests/specialized.py
$(SPECIALIZED_OBJECTS): $(BUILD)/%_specialized.o: generated/%_specialized.c include/chem.h $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
$(addprefix generated/,$(addsuffix _specialized.c,$(SPEC_NAMES))): generated/%_specialized.c: generated/%.c tools/specialize.py tools/configure.py tools/rosenbrock3.c.in generated/specialized-flags | env
	$(PYTHON) tools/specialize.py $(CONFIG_$*) $@ --name $*_specialized $(SPECFLAGS)
CONFIG_ts1 := configs/ts1.json
CONFIG_example_no_user := configs/examples/rate_constants_no_user_defined/config.yaml
CONFIG_example_user := configs/examples/rate_constants_user_defined/config.json
CONFIG_example_robertson := configs/examples/robertson/config.json
CONFIG_validation := configs/validation.yaml
CONFIG_edges := configs/specialization_edges.yaml
CONFIG_ts1_variant := generated/ts1_variant.json
generated/edges.c: configs/specialization_edges.yaml tools/configure.py generated/flags | env
	$(PYTHON) tools/configure.py $< $@ --name edges $(GENFLAGS)
generated/ts1_variant.json: configs/ts1.json tools/specialization_fixture.py | env
	$(PYTHON) tools/specialization_fixture.py $< $@
generated/ts1_variant.c: generated/ts1_variant.json tools/configure.py generated/flags | env
	$(PYTHON) tools/configure.py $< $@ --name ts1_variant $(GENFLAGS)
$(BUILD)/edges.o: generated/edges.c include/chem.h $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
$(BUILD)/ts1_variant.o: generated/ts1_variant.c include/chem.h $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
generated/specialized-flags: force
	@mkdir -p generated
	@echo '$(SPECFLAGS)' > generated/specialized-flags.new
	@cmp -s generated/specialized-flags.new $@ || cp generated/specialized-flags.new $@
$(BUILD)/libspecialized.$(SHARED_EXT): $(SPECIALIZED_OBJECTS) $(BUILD)/edges.o $(BUILD)/ts1_variant.o $(BUILD)/libchem.a
	$(CC) $(OPTFLAGS) $(SHARED_FLAG) $^ $(LDLIBS) -o $@
bin/libspecialized.$(SHARED_EXT): $(BUILD)/libspecialized.$(SHARED_EXT) force
	@mkdir -p bin
	cp $< $@
$(BUILD)/bench_specialized.o: src/bench.c src/bench_io.h $(wildcard include/*.h) $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -DCHEM_SPECIALIZED -c $< -o $@
$(BUILD)/ts1_bench_specialized: $(BUILD)/bench_specialized.o $(BUILD)/ts1_specialized.o $(BUILD)/libchem.a
	$(CC) $(OPTFLAGS) $^ $(LDLIBS) -o $@
$(BUILD)/ts1_scenario_specialized: $(BUILD)/scenario.o $(BUILD)/ts1_specialized.o $(BUILD)/libchem.a
	$(CC) $(OPTFLAGS) $^ $(LDLIBS) -o $@
bin/ts1_bench_specialized: $(BUILD)/ts1_bench_specialized force
	@mkdir -p bin
	cp $< $@
bin/ts1_scenario_specialized: $(BUILD)/ts1_scenario_specialized force
	@mkdir -p bin
	cp $< $@
$(BUILD)/%.o: src/%.c $(wildcard include/*.h) src/internal.h src/bench_io.h $(BUILD)/flags
	@mkdir -p $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
$(BUILD)/ts1.o: generated/ts1.c include/chem.h $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
$(BUILD)/example_%.o: generated/example_%.c include/chem.h $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
$(BUILD)/validation.o: generated/validation.c include/chem.h $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
generated/validation.c: configs/validation.yaml tools/configure.py generated/flags | env
	$(PYTHON) tools/configure.py $< $@ --name validation $(GENFLAGS)
generated/example_no_user.c: $(wildcard configs/examples/rate_constants_no_user_defined/*) tools/configure.py generated/flags | env
	$(PYTHON) tools/configure.py configs/examples/rate_constants_no_user_defined/config.yaml $@ --name example_no_user $(GENFLAGS)
generated/example_user.c: $(wildcard configs/examples/rate_constants_user_defined/*) tools/configure.py generated/flags | env
	$(PYTHON) tools/configure.py configs/examples/rate_constants_user_defined/config.json $@ --name example_user $(GENFLAGS)
generated/example_robertson.c: $(wildcard configs/examples/robertson/*) tools/configure.py generated/flags | env
	$(PYTHON) tools/configure.py configs/examples/robertson/config.json $@ --name example_robertson $(GENFLAGS)
generated/ts1.c: configs/ts1.json tools/configure.py generated/flags | env
	$(PYTHON) tools/configure.py $< $@ --name ts1 $(GENFLAGS)
.PHONY: force
force:
generated/flags: force
	@mkdir -p generated
	@echo '$(GENFLAGS)' > generated/flags.new
	@cmp -s generated/flags.new $@ || cp generated/flags.new $@
build/reference-flags: force
	@mkdir -p build $(TMPDIR)
	@echo 'clang++ $(OPTFLAGS)' > build/reference-flags.new
	@cmp -s build/reference-flags.new $@ || cp build/reference-flags.new $@
$(BUILD)/flags: force
	@mkdir -p $(BUILD) $(TMPDIR)
	@echo '$(CC) $(CFLAGS) $(CPPFLAGS) $(LDLIBS)' > $(BUILD)/flags.new
	@cmp -s $(BUILD)/flags.new $@ || cp $(BUILD)/flags.new $@
$(BUILD)/libchem.a: $(OBJECTS) $(GENERATED_OBJECTS)
	ar rcs $@ $(OBJECTS) $(GENERATED_OBJECTS)
$(BUILD)/libchem.$(SHARED_EXT): $(OBJECTS) $(GENERATED_OBJECTS)
	$(CC) $(OPTFLAGS) $(SHARED_FLAG) $(OBJECTS) $(GENERATED_OBJECTS) $(LDLIBS) -o $@
$(BUILD)/ts1_bench_cpu: $(BUILD)/bench.o $(BUILD)/libchem.a
	$(CC) $(OPTFLAGS) $(BUILD)/bench.o $(BUILD)/libchem.a $(LDLIBS) -o $@
$(BUILD)/tutorials: tests/tutorials.c $(BUILD)/libchem.a include/chem.h include/chem_array.h $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(BUILD)/libchem.a $(LDLIBS) -o $@
bin/tutorials: $(BUILD)/tutorials force
	@mkdir -p bin
	cp $< $@
bin/libchem.a: $(BUILD)/libchem.a force
	@mkdir -p bin
	cp $< $@
bin/libchem.$(SHARED_EXT): $(BUILD)/libchem.$(SHARED_EXT) force
	@mkdir -p bin
	cp $< $@
bin/ts1_bench_cpu: $(BUILD)/ts1_bench_cpu force
	@mkdir -p bin
	cp $< $@
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: $(ROOT)/reference/micm_ts1_ref.cpp build/reference-flags
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(ROOT)/micm/include -I$(ROOT)/micm/test/benchmark $< -o $@
clean:
	rm -rf build bin generated
