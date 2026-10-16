.DEFAULT_GOAL := cpu
ifeq ($(origin CC),default)
CC := clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
CONDA ?= $(HOME)/miniforge3/bin/conda
PYTHON := $(CURDIR)/.conda/bin/python
MICM ?= ../micm
REFERENCE ?= ../reference
UNAME_S := $(shell uname -s)
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -n 1)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
endif
ifneq (,$(findstring clang,$(CC_VERSION)))
ifeq ($(UNAME_S),Darwin)
OMP_CFLAGS := -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include
OMP_LIBS := -L/opt/homebrew/opt/libomp/lib -lomp
else
OMP_CFLAGS := -fopenmp
OMP_LIBS := -fopenmp
endif
else
OMP_CFLAGS := -fopenmp
OMP_LIBS := -fopenmp
endif
CFLAGS := $(OPTFLAGS) -std=c11 -Wall -Wextra -Wpedantic $(OMP_CFLAGS) -Iinclude
LDLIBS := $(OMP_LIBS) -lm

.PHONY: env cpu test-cpu bench micm-ref clean

env: .conda/conda-meta/history
.conda/conda-meta/history: environment.yml
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml --yes

cpu:
	@echo 'Implementation in progress'; exit 1

test-cpu: cpu
bench: cpu

micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: $(REFERENCE)/micm_ts1_ref.cpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(MICM)/include -I$(MICM)/test/benchmark $< -o $@

clean:
	rm -rf build bin generated
