#!/bin/bash
# Finer cell groups at a whole node: 8, 12, and 16 cells at 128 threads, and
# 8 against 16 at one thread, with the final build (derecho/work/variants/final),
# MICM's fastest at 128 threads (vector8, one group per thread) in the same
# rounds. Part A: 5 rounds at 128 threads; part B: 3 rounds at one thread.
#   derecho/tools/submit.sh grouping-fine 128 01:30:00 pbs/grouping_fine.sh OUT
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1")
micm_final=$here/work/micm-reference/micm_final.csv
cells=$here/work/cells_2026.csv
bench=$here/work/variants/final/bin/ts1_bench_cpu
micm=$here/work/micm-configs/bin/micm-vector8-in-place-mozart
cd "$here"
[[ -x $bench && -x $micm ]] || { echo "missing $bench or $micm"; exit 1; }
mkdir -p "$out"
tools/environment.sh "$out"
cp work/variants/final/variant.txt "$out/variant-final.txt"
echo "== part A: 128 threads ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out/threads-128" --threads 128 --rounds 5 \
    "micm-v8=$micm" \
    "g8=$bench --group-size 8" \
    "g12=$bench --group-size 12" \
    "g16=$bench --group-size 16" \
    "g8-static=$bench --group-size 8 --schedule static"
echo "== part B: 1 thread ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out/threads-1" --threads 1 --rounds 3 \
    "g8=$bench --group-size 8" \
    "g16=$bench --group-size 16"
echo "== done ($(date -Is))"
