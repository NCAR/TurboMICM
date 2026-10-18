#!/bin/bash
# Session 1 of the Fable entry on Derecho: which width and flags, against MICM.
#
#   derecho/pbs/port_session_1.sh OUT NAME=COMMAND...
#
# Builds five variants of the port in parallel (derecho/tools/build_variant.sh):
# the Derecho defaults (width 4, -O3 -march=znver3 -ffp-contract=fast, OpenMP
# simd lane loops, aligned arrays), widths 2 and 8, the defaults without FMA
# contraction, and plain -O3 for the generic x86-64; then times them with the
# given MICM programs (NAME=COMMAND, as admin/time_entries.py takes them) at 128
# and 1 threads for 3 rotating rounds on the 10,000-cell seed-2026 input, every
# run gated against MICM's final state from derecho/tools/micm_reference.sh.
# It waits for that reference and for the MICM programs if the baseline job is
# still producing them. Run with a whole node:
#   derecho/tools/submit.sh port-session-1 128 02:30:00 pbs/port_session_1.sh OUT micm-default=... micm-fastest=...
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1"); shift
micm_final=$here/work/micm-reference/micm_final.csv
cells=$here/work/cells_2026.csv
cd "$here"
mkdir -p "$out" work/variants
tools/environment.sh "$out"
variants=(w4 w2 w8 w4-nofma w4-generic-x86)

echo "== building variants ($(date -Is))"
tools/build_variant.sh w4 > work/variants/w4.build.log 2>&1 &
tools/build_variant.sh w2 WIDTH=2 > work/variants/w2.build.log 2>&1 &
tools/build_variant.sh w8 WIDTH=8 > work/variants/w8.build.log 2>&1 &
tools/build_variant.sh w4-nofma OPTFLAGS="-O3 -DNDEBUG -march=znver3" > work/variants/w4-nofma.build.log 2>&1 &
tools/build_variant.sh w4-generic-x86 OPTFLAGS="-O3 -DNDEBUG" > work/variants/w4-generic-x86.build.log 2>&1 &
wait
programs=()
for v in "${variants[@]}"; do
  [[ -x work/variants/$v/bin/ts1_bench_cpu ]] || { echo "variant $v failed:"; tail -20 work/variants/$v.build.log; exit 1; }
  cp work/variants/$v/variant.txt "$out/variant-$v.txt"
  programs+=("$v=$here/work/variants/$v/bin/ts1_bench_cpu")
done
grep -H 'build seconds' work/variants/{w4,w2,w8,w4-nofma,w4-generic-x86}/variant.txt

# the gate's reference and the MICM programs come from the baseline job
# admin/time_entries.py runs every program from a fresh working directory, so
# the programs' paths must be absolute: relative ones are taken from derecho/
micm=()
for item in "$@"; do
  name=${item%%=*}; command=${item#*=}
  [[ $command == /* ]] || command="$here/$command"
  micm+=("$name=$command")
done
for i in $(seq 1 120); do
  ready=1
  [[ -f $micm_final ]] || ready=0
  for item in "${micm[@]}"; do
    program=${item#*=}; program=${program%% *}
    [[ -x $program ]] || ready=0
  done
  [[ $ready == 1 ]] && break
  sleep 30
done
[[ -f $micm_final ]] || { echo "no $micm_final"; exit 1; }

echo "== timing ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out" --threads 128,1 --rounds 3 \
    "${micm[@]}" "${programs[@]}"
