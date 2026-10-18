#!/bin/bash
# Session 4 of the Fable entry on Derecho: the one-thread knobs.
#
#   derecho/pbs/port_session_4.sh OUT MICM_DEFAULT=COMMAND MICM_FASTEST=COMMAND
#
# Variants (built already by derecho/tools/build_variant.sh, or built here):
#   w4-prenorm   the defaults before the fast error norm (the control)
#   w4n          the defaults: fast norm (per-lane partial sums), omp simd lanes
#   split8       solve rows of 8 or more terms split into partial sums (12)
#   split6       ... of 6 or more
#   inter8       8 accumulations side by side in the gathers (4)
#   inter2       2 accumulations side by side
#   pairs        the solves take two blocks per loop iteration, interleaved
# Part A times every variant with both MICM programs at 128 threads; part B
# times them with MICM's fastest at one thread. 3 rotating rounds each, the
# 10,000-cell seed-2026 input, every run gated against MICM's final state.
#   derecho/tools/submit.sh port-session-4 128 02:00:00 pbs/port_session_4.sh OUT micm-default=... micm-v8-g32=...
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1"); shift
micm_final=$here/work/micm-reference/micm_final.csv
cells=$here/work/cells_2026.csv
cd "$here"
mkdir -p "$out" work/variants
tools/environment.sh "$out"
micm=()
for item in "$@"; do
  name=${item%%=*}; command=${item#*=}
  [[ $command == /* ]] || command="$here/$command"
  micm+=("$name=$command")
done
variants=(w4-prenorm w4n split8 split6 inter8 inter2 pairs)

echo "== building variants ($(date -Is))"
[[ -x work/variants/w4n/bin/ts1_bench_cpu ]] || tools/build_variant.sh w4n > work/variants/w4n.build.log 2>&1 &
[[ -x work/variants/split8/bin/ts1_bench_cpu ]] || tools/build_variant.sh split8 GEN_FLAGS="--split-min 8" > work/variants/split8.build.log 2>&1 &
[[ -x work/variants/split6/bin/ts1_bench_cpu ]] || tools/build_variant.sh split6 GEN_FLAGS="--split-min 6" > work/variants/split6.build.log 2>&1 &
[[ -x work/variants/inter8/bin/ts1_bench_cpu ]] || tools/build_variant.sh inter8 GEN_FLAGS="--interleave 8" > work/variants/inter8.build.log 2>&1 &
[[ -x work/variants/inter2/bin/ts1_bench_cpu ]] || tools/build_variant.sh inter2 GEN_FLAGS="--interleave 2" > work/variants/inter2.build.log 2>&1 &
[[ -x work/variants/pairs/bin/ts1_bench_cpu ]] || tools/build_variant.sh pairs GEN_FLAGS="--solve-pairs" > work/variants/pairs.build.log 2>&1 &
wait
programs=()
for v in "${variants[@]}"; do
  [[ -x work/variants/$v/bin/ts1_bench_cpu ]] || { echo "variant $v missing"; exit 1; }
  cp work/variants/$v/variant.txt "$out/variant-$v.txt"
  programs+=("$v=$here/work/variants/$v/bin/ts1_bench_cpu")
done

echo "== part A: 128 threads ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out/threads-128" --threads 128 --rounds 3 \
    "${micm[@]}" "${programs[@]}"
echo "== part B: 1 thread ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out/threads-1" --threads 1 --rounds 3 \
    "${micm[1]}" "${programs[@]}"
echo "== done ($(date -Is))"
