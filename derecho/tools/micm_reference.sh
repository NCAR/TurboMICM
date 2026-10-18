#!/bin/bash
# MICM's final state for a TS1 input on Derecho: the gate's reference.
#
#   [OPTFLAGS="..."] [THREADS=1] derecho/tools/micm_reference.sh OUT CELLS_CSV [CELLS] [STEPS] [DT]
#
# Builds reference/micm_ts1_ref.cpp (MICM's default configuration: standard
# ordering, in-place Mozart LU) with `g++ -std=c++20 $OPTFLAGS` into
# OUT/micm_ts1_ref, runs it with THREADS threads (1; MICM's final state does
# not depend on the thread count or grouping, as every run here confirmed) on
# CELLS (10000) STEPS (100) DT (30) and CELLS_CSV, and leaves
# OUT/micm_final.csv, OUT/micm_run.txt (its output), and OUT/flags.txt. Run
# inside a PBS job (derecho/tools/submit.sh). An existing OUT/micm_ts1_ref
# built with the same flags is reused.
#
# No -fopenmp: the reference threads with std::thread, and a linked libgomp
# pins the main thread to one core under OMP_PROC_BIND (env.sh), which every
# std::thread inherits: 128 threads on one core.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1")
cells_csv=$(realpath "$2")
cells=${3:-10000}
steps=${4:-100}
dt=${5:-30}
threads=${THREADS:-1}
optflags=${OPTFLAGS:--O3 -DNDEBUG -march=znver3 -ffp-contract=fast}
mkdir -p "$out"
echo "compiler: $(g++ --version | head -1); flags: $optflags"
if [[ ! -x $out/micm_ts1_ref || "$(cat "$out/flags.txt" 2>/dev/null)" != "$optflags" ]]; then
  printf '%s\n' "$optflags" > "$out/flags.txt"
  g++ -std=c++20 $optflags -I "$here/micm/include" -I "$here/micm/test/benchmark" \
      "$repo/reference/micm_ts1_ref.cpp" -o "$out/micm_ts1_ref"
fi
cd "$out"
./micm_ts1_ref --threads "$threads" "$cells" "$steps" "$dt" "$cells_csv" | tee micm_run.txt
mv ts1_final.csv micm_final.csv
echo "wrote $out/micm_final.csv"
