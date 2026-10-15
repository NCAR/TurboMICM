ROOT := $(abspath ..)
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON := $(CURDIR)/.conda/bin/python
OPTFLAGS ?= -O3 -DNDEBUG
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
CPPFLAGS := -Iinclude -Igenerated
LDLIBS := $(OMPLIBS) -lm
BUILD := build/$(notdir $(CC))
SOURCES := $(filter-out src/bench.c,$(wildcard src/*.c))
OBJECTS := $(patsubst src/%.c,$(BUILD)/%.o,$(SOURCES))

.PHONY: all env cpu test-cpu bench micm-ref clean
all: cpu
env: .conda/conda-meta/history
.conda/conda-meta/history: environment.yml
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes
cpu: bin/libchem.a bin/libchem.$(SHARED_EXT) bin/ts1_bench_cpu
bench: cpu
test-cpu: cpu
	$(PYTHON) tests/run.py
$(BUILD)/%.o: src/%.c $(wildcard include/*.h) src/internal.h $(BUILD)/flags
	@mkdir -p $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
$(BUILD)/ts1.o: generated/ts1.c include/chem.h $(BUILD)/flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c $< -o $@
generated/ts1.c: configs/ts1.json tools/configure.py | env
	$(PYTHON) tools/configure.py $< $@ --name ts1
.PHONY: force
force:
$(BUILD)/flags: force
	@mkdir -p $(BUILD)
	@echo '$(CC) $(CFLAGS) $(CPPFLAGS) $(LDLIBS)' > $(BUILD)/flags.new
	@cmp -s $(BUILD)/flags.new $@ || cp $(BUILD)/flags.new $@
bin/libchem.a: $(OBJECTS)
	@mkdir -p bin
	ar rcs $@ $(OBJECTS)
bin/libchem.$(SHARED_EXT): $(OBJECTS)
	@mkdir -p bin
	$(CC) $(OPTFLAGS) $(SHARED_FLAG) $(OBJECTS) $(LDLIBS) -o $@
bin/ts1_bench_cpu: $(BUILD)/bench.o $(BUILD)/ts1.o bin/libchem.a
	$(CC) $(OPTFLAGS) $(BUILD)/bench.o $(BUILD)/ts1.o bin/libchem.a $(LDLIBS) -o $@
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: $(ROOT)/reference/micm_ts1_ref.cpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(ROOT)/micm/include -I$(ROOT)/micm/test/benchmark $< -o $@
clean:
	rm -rf build bin generated
