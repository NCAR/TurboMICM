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
.PHONY: env cpu bench test-cpu micm-ref clean
env: .conda/bin/python
.conda/bin/python: environment.yml
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(CURDIR)/.conda --file environment.yml
cpu: bin/ts1_bench_cpu
bench: cpu
micm-ref: bin/micm_ts1_ref
bin/micm_ts1_ref: ../reference/micm_ts1_ref.cpp ../micm/test/benchmark/ts1_mechanism.hpp
	@mkdir -p bin
	clang++ -std=c++20 $(OPTFLAGS) -I../micm/include -I../micm/test/benchmark $< -o $@
clean:
	rm -rf build bin generated
