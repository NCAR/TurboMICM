ifeq ($(origin CC),default)
CC := clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
CONDA ?= $(HOME)/miniforge3/bin/conda
ENV := $(CURDIR)/.conda
PYTHON := $(ENV)/bin/python
ROOT := $(abspath ..)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
ifneq ($(findstring clang,$(shell $(CC) --version 2>/dev/null | head -1)),)
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
CFLAGS := -std=c11 -Wall -Wextra -Wpedantic $(OPTFLAGS) $(OMP_CFLAGS) -Iinclude
LDLIBS := $(OMP_LIBS) -lm

.PHONY: env clean micm-ref
env:
	CONDA_PKGS_DIRS=$(CURDIR)/.conda-pkgs $(CONDA) env create --prefix $(ENV) --file environment.yml

bin:
	mkdir -p $@

micm-ref: | bin
	clang++ -std=c++20 $(OPTFLAGS) -I$(ROOT)/micm/include -I$(ROOT)/micm/test/benchmark $(ROOT)/reference/micm_ts1_ref.cpp -o bin/micm_ts1_ref

clean:
	rm -rf build bin
