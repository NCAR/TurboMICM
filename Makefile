CC = clang
OPTFLAGS ?= -O3 -DNDEBUG
CFLAGS = $(OPTFLAGS) -std=c11 -Wall -Wextra -Wpedantic -Iinclude
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
export SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path)
endif

.PHONY: env cpu check clean
env:
	CONDA_PKGS_DIRS="$(CURDIR)/.conda-pkgs" $(HOME)/miniforge3/bin/conda env create --prefix "$(CURDIR)/.conda" --file environment.yml
cpu: build/libchem.a
build:
	mkdir -p $@
build/chem.o: src/chem.c include/chem.h | build
	$(CC) $(CFLAGS) -c $< -o $@
build/libchem.a: build/chem.o
	ar rcs $@ $^
build/test_chem: tests/test_chem.c build/libchem.a
	$(CC) $(filter-out -DNDEBUG,$(CFLAGS)) $< build/libchem.a -lm -o $@
clean:
	rm -rf build
build/robertson.c: tools/generate.py tools/config.py | build
	.conda/bin/python tools/generate.py ../micm/test/tutorial/configs/robertson/config.json $@
build/test_generated: tests/test_generated.c build/robertson.c build/libchem.a
	$(CC) $(filter-out -DNDEBUG,$(CFLAGS)) tests/test_generated.c build/robertson.c build/libchem.a -lm -o $@
check: build/test_chem build/test_generated
	./build/test_chem
	./build/test_generated
	.conda/bin/python tests/test_config.py
