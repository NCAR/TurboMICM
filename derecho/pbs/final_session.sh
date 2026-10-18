#!/bin/bash
# The final timing session of the Fable entry on Derecho, as the experiment's
# procedure prescribes: 5 rotating rounds at a whole node and at one thread,
# every run gated against MICM's final state.
#
#   derecho/pbs/final_session.sh OUT VARIANT THREADS NAME=COMMAND...
#
# VARIANT is a built variant (derecho/work/variants/VARIANT) run in its default
# configuration, and THREADS the comma-separated thread counts (e.g. 128,1).
# The NAME=COMMAND programs are MICM's default and fastest configurations.
#   derecho/tools/submit.sh final 128 03:00:00 pbs/final_session.sh OUT final 128,1 micm-default=... micm-fastest=...
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1"); variant=$2; threads=$3; shift 3
micm_final=$here/work/micm-reference/micm_final.csv
cells=$here/work/cells_2026.csv
bench=$here/work/variants/$variant/bin/ts1_bench_cpu
cd "$here"
# admin/time_entries.py runs every program from a fresh working directory, so
# the programs' paths must be absolute: relative ones are taken from derecho/
micm=()
for item in "$@"; do
  name=${item%%=*}; command=${item#*=}
  [[ $command == /* ]] || command="$here/$command"
  micm+=("$name=$command")
done
[[ -x $bench ]] || { echo "no $bench"; exit 1; }
mkdir -p "$out"
tools/environment.sh "$out"
cp "work/variants/$variant/variant.txt" "$out/variant-$variant.txt"
echo "== timing ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out" --threads "$threads" --rounds 5 \
    "${micm[@]}" "fable-derecho=$bench"
echo "== done ($(date -Is))"
