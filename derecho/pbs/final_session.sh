#!/bin/bash
# The final timing session of the Fable entry on Derecho, as the experiment's
# procedure prescribes: rotating rounds at a whole node and at one thread,
# every run gated against MICM's final state.
#
#   derecho/pbs/final_session.sh OUT THREADS ROUNDS NAME=COMMAND...
#
# THREADS is the comma-separated thread counts (e.g. 128,256,1) and ROUNDS the
# number of rounds (5 for a ranking). The NAME=COMMAND programs are MICM's
# default and fastest configurations and the port's builds
# (derecho/work/variants/NAME/bin/ts1_bench_cpu, with run-time options);
# relative paths are taken from derecho/. Each port variant's variant.txt is
# copied to OUT.
#   derecho/tools/submit.sh final 128 04:00:00 pbs/final_session.sh OUT 128,1 5 micm-default=... micm-fastest=... fable-derecho=work/variants/final/bin/ts1_bench_cpu
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1"); threads=$2; rounds=$3; shift 3
micm_final=$here/work/micm-reference/micm_final.csv
cells=$here/work/cells_2026.csv
cd "$here"
mkdir -p "$out"
tools/environment.sh "$out"
programs=()
for item in "$@"; do
  name=${item%%=*}; command=${item#*=}
  [[ $command == /* ]] || command="$here/$command"
  programs+=("$name=$command")
  program=${command%% *}
  if [[ $program == */work/variants/*/bin/ts1_bench_cpu ]]; then
    variant=${program%/bin/ts1_bench_cpu}
    cp "$variant/variant.txt" "$out/variant-$(basename "$variant").txt"
  fi
done
echo "== timing ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out" --threads "$threads" --rounds "$rounds" "${programs[@]}"
echo "== done ($(date -Is))"
