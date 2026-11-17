ifeq ($(origin CC),default)
CC = clang
endif
OPTFLAGS ?= -O3 -DNDEBUG
PYTHON = .conda/bin/python
CONDA ?= $(HOME)/miniforge3/bin/conda
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Iinclude $(OPTFLAGS)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
endif

.PHONY: env cpu test-reader clean

env:
	CONDA_PKGS_DIRS="$(CURDIR)/.conda-pkgs" $(CONDA) env create --prefix "$(CURDIR)/.conda" --file environment.yml

build bin:
	mkdir -p $@

build/chem.o: src/chem.c include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@

build/ts1.c: configs/ts1/config.json configs/ts1/species.json configs/ts1/reactions.json tools/generate.py tools/mechanism.py | build
	$(PYTHON) tools/generate.py $< $@

build/ts1.o: build/ts1.c include/chem.h
	$(CC) $(CFLAGS) -c $< -o $@

build/libchem.a: build/chem.o
	ar rcs $@ $^

cpu: build/libchem.a build/ts1.o

build/test_chem: tests/test_chem.c build/ts1.o build/libchem.a
	$(CC) $(CFLAGS) $^ -lm -o $@

test-reader: build/test_chem
	$(PYTHON) -m unittest discover -s tests -v
	./build/test_chem

clean:
	rm -rf build bin
