ROOT := $(abspath ..)
.SILENT:
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON := $(CURDIR)/.conda/bin/python
OPTFLAGS ?= -O3 -DNDEBUG
GROUP_SIZE ?= 32
KERNEL_CHUNK ?= 32
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
OMPFLAGS := -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include
OMPLIBS := -L/opt/homebrew/opt/libomp/lib -lomp
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
LDLIBS := $(OMPLIBS) -lm
BUILD := build/$(notdir $(CC))
SOURCES := $(filter-out src/bench.c,$(wildcard src/*.c))
OBJECTS := $(patsubst src/%.c,$(BUILD)/%.o,$(SOURCES))
GENERATED_OBJECTS := $(BUILD)/ts1.o $(BUILD)/example_no_user.o $(BUILD)/example_user.o $(BUILD)/example_robertson.o $(BUILD)/validation.o

.PHONY: all env cpu test-cpu bench micm-ref clean
all: cpu
env: .conda/conda-meta/history
.conda/conda-meta/history: environment.yml
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes >&2
cpu: bin/libchem.a bin/libchem.$(SHARED_EXT) bin/ts1_bench_cpu bin/tutorials
bench: cpu
test-cpu: cpu
	$(PYTHON) tests/run.py
$(BUILD)/%.o: src/%.c $(wildcard include/*.h) src/internal.h $(BUILD)/flags
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
	@mkdir -p build
	@echo 'clang++ $(OPTFLAGS)' > build/reference-flags.new
	@cmp -s build/reference-flags.new $@ || cp build/reference-flags.new $@
$(BUILD)/flags: force
	@mkdir -p $(BUILD)
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
