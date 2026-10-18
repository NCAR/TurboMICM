#!/bin/bash
# Build a variant of the Derecho copy of the Fable entry in its own directory.
#
#   derecho/tools/build_variant.sh NAME [MAKE_VARIABLE=value ...]
#
# Copies derecho/fable (sources, tools, tests; not .conda, build, bin, or
# generated) to derecho/work/variants/NAME, links the conda environment, and
# runs `make cpu` there with the make variables (WIDTH, OPTFLAGS, GROUP_SIZE,
# KERNEL_CHUNK, KERNEL_MODE, GEN_FLAGS, CC, ...). The benchmark is then
# derecho/work/variants/NAME/bin/ts1_bench_cpu, and NAME/variant.txt records the
# variables, the compiler, the flags, the build time, and the kernel code size.
# An existing variant is rebuilt only where make finds work to do.
set -eu
name=$1; shift
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
dir=$here/work/variants/$name
mkdir -p "$dir"
rsync -a --exclude .conda --exclude .conda-pkgs --exclude build --exclude bin --exclude generated --exclude results \
    --exclude '*.o' --exclude '*.a' --exclude ts1_final.csv --exclude __pycache__ "$here/fable/" "$dir/"
[[ -e $dir/.conda ]] || ln -s "$here/fable/.conda" "$dir/.conda"
cd "$dir"
start=$(date +%s)
make cpu MICM="$here/micm" REFERENCE="$repo/reference" "$@" 2>&1 | grep -vE '^(ar |.*/gcc .*-c src/)' || true
end=$(date +%s)
[[ -x bin/ts1_bench_cpu ]] || { echo "build_variant: $name failed"; exit 1; }
{
  echo "variant: $name"
  echo "variables: $*"
  echo "compiler: $(sed -n 1p build/flags.txt) ($($(sed -n 1p build/flags.txt) --version | head -1))"
  echo "cflags: $(sed -n 2p build/flags.txt)"
  echo "build seconds: $((end - start))"
  echo "generated/ts1.c lines: $(wc -l < generated/ts1.c)"
  echo "ts1.o text bytes: $(size build/ts1.o | awk 'NR==2 {print $1}')"
  echo "built: $(date -Is) on $(hostname)"
} | tee variant.txt
