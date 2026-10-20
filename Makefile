ifeq ($(origin CC),default)
CC := clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
PYTHON ?= ./.conda/bin/python
ROOT := $(abspath ..)
BUILD := build
BIN := bin
GEN := $(BUILD)/generated
LIBSRC := src/chem.c
LIBHDR := include/chem.h
CONFIG := config/ts1.json
GENSRC := $(GEN)/ts1.c

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

$(BUILD)/chem.o: $(LIBSRC) $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c $(LIBSRC) -o $@

$(BUILD)/ts1.o: src/ts1_bench.c $(LIBHDR)
	@mkdir -p $(BUILD) $(BIN)
	$(CC) $(CFLAGS) -c src/ts1_bench.c -o $@

$(BUILD)/ts1_mechanism.o: $(GENSRC) $(LIBHDR)
	$(CC) $(CFLAGS) -c $(GENSRC) -o $@

$(BIN)/libchem.a: $(BUILD)/chem.o
	@mkdir -p $(BIN)
	ar rcs $@ $^

$(BIN)/ts1_bench_cpu: $(BUILD)/ts1.o $(BUILD)/ts1_mechanism.o $(BIN)/libchem.a
	$(CC) $(OPTFLAGS) $(OPENMP_CFLAGS) $^ $(LDFLAGS) -o $@

test-cpu:
	$(PYTHON) tools/test_cases.py

bench: cpu

micm-ref:
	@mkdir -p $(BIN)
	SDKROOT=$$(xcrun --sdk macosx --show-sdk-path) clang++ -std=c++20 $(OPTFLAGS) \
	  -I $(ROOT)/micm/include -I $(ROOT)/micm/test/benchmark \
	  $(ROOT)/reference/micm_ts1_ref.cpp -o $(BIN)/micm_ts1_ref

clean:
	rm -rf $(BUILD) $(BIN)
