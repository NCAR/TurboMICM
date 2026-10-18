#!/bin/bash
# Session 2 of the Fable entry on Derecho: code-generation and tuning knobs.
#
#   derecho/pbs/port_session_2.sh OUT NAME=COMMAND...
#
# Builds variants of the Derecho defaults (width 4, -O3 -march=znver3
# -ffp-contract=fast, OpenMP simd lane loops, 64-byte aligned arrays) in parallel:
#   w4           the defaults again (the control)
#   nosimd       plain lane loops (--no-lane-pragma): GCC's SLP vectorizer
#   fmachains    GCC's Zen tuning keeps a*b+c chains as multiply and add;
#                -mtune-ctrl=^avoid_fma_chains,^avoid_fma256_chains fuses them
#   panel16      LU panels of 16 pivots (16 vector registers) instead of 64
#   panel32      LU panels of 32 pivots
#   noalign      malloc instead of aligned_alloc (TMC_ALIGN=0)
# then times them with the given MICM programs at 128 and 1 threads for 3
# rotating rounds on the 10,000-cell seed-2026 input, gated against MICM. It
# waits for the gate's reference and the MICM programs if the baseline job is
# still producing them.
#   derecho/tools/submit.sh port-session-2 128 02:30:00 pbs/port_session_2.sh OUT micm-default=... micm-fastest=...
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1"); shift
micm_final=$here/work/micm-reference/micm_final.csv
cells=$here/work/cells_2026.csv
cd "$here"
mkdir -p "$out" work/variants
tools/environment.sh "$out"
base="-O3 -DNDEBUG -march=znver3 -ffp-contract=fast"
fma='-mtune-ctrl=^avoid_fma_chains,^avoid_fma256_chains'
variants=(w4 nosimd fmachains panel16 panel32 noalign)

echo "== building variants ($(date -Is))"
tools/build_variant.sh w4 > work/variants/w4.build.log 2>&1 &
tools/build_variant.sh nosimd GEN_FLAGS=--no-lane-pragma > work/variants/nosimd.build.log 2>&1 &
tools/build_variant.sh fmachains OPTFLAGS="$base $fma" > work/variants/fmachains.build.log 2>&1 &
tools/build_variant.sh panel16 GEN_FLAGS="--panel 16" > work/variants/panel16.build.log 2>&1 &
tools/build_variant.sh panel32 GEN_FLAGS="--panel 32" > work/variants/panel32.build.log 2>&1 &
tools/build_variant.sh noalign OPTFLAGS="$base -DTMC_ALIGN=0" > work/variants/noalign.build.log 2>&1 &
wait
programs=()
for v in "${variants[@]}"; do
  [[ -x work/variants/$v/bin/ts1_bench_cpu ]] || { echo "variant $v failed:"; tail -20 work/variants/$v.build.log; exit 1; }
  cp work/variants/$v/variant.txt "$out/variant-$v.txt"
  programs+=("$v=$here/work/variants/$v/bin/ts1_bench_cpu")
done
grep -H 'build seconds' work/variants/{w4,nosimd,fmachains,panel16,panel32,noalign}/variant.txt

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
