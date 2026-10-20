CC ?= cc
OPTFLAGS ?= -O3 -DNDEBUG
PYTHON ?= python3
UNAME_S := $(shell uname -s)

CPPFLAGS += -Iinclude
CFLAGS += -std=c11 -Wall -Wextra -Wpedantic
LDFLAGS += -lm

ifeq ($(UNAME_S),Darwin)
  SDKROOT ?= $(shell xcrun --sdk macosx --show-sdk-path)
  export SDKROOT
  ifneq (,$(findstring clang,$(shell $(CC) --version 2>/dev/null | head -1)))
    OMPFLAGS = -Xpreprocessor -fopenmp -isystem /opt/homebrew/opt/libomp/include
    OMPLIBS = -L/opt/homebrew/opt/libomp/lib -lomp
  else ifneq (,$(findstring gcc,$(shell $(CC) --version 2>/dev/null | head -1)))
    OMPFLAGS = -fopenmp
    OMPLIBS = -fopenmp
  endif
else
  OMPFLAGS = -fopenmp
  OMPLIBS = -fopenmp
endif

.PHONY: all cpu test-cpu bench clean env
all: cpu
cpu: bin/ts1_bench_cpu
bench: bin/ts1_bench_cpu

generated/ts1.c: mechanisms/ts1.json tools/generate_c.py
	$(PYTHON) tools/generate_c.py $< $@

bin/ts1_bench_cpu: src/main.c src/mechanism.c generated/ts1.c include/chemistry.h
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) $(OMPFLAGS) src/main.c src/mechanism.c generated/ts1.c -o $@ $(LDFLAGS) $(OMPLIBS)

bin/test_core: tests/test_core.c src/mechanism.c include/chemistry.h
	@mkdir -p bin
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OPTFLAGS) $(OMPFLAGS) tests/test_core.c src/mechanism.c -o $@ $(LDFLAGS) $(OMPLIBS)

test-cpu: cpu bin/test_core
	@bin/test_core
	$(PYTHON) tools/test_cases.py

env:
	CONDA_PKGS_DIRS="$(CURDIR)/.conda-pkgs" conda env create --prefix "$(CURDIR)/.conda" --file environment.yml

clean:
	rm -rf bin build generated ts1_final.csv
