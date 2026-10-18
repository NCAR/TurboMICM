#!/bin/bash
# Session 3 of the Fable entry on Derecho: grouping, scheduling, and threads.
#
#   derecho/pbs/port_session_3.sh OUT VARIANT NAME=COMMAND...
#
# Run-time choices of one built variant (derecho/work/variants/VARIANT, from a
# previous session): cell groups of 16, 32, 64, 79 (one group per thread at 128
# threads, as MICM's default grouping), and 128; the static schedule, which
# keeps each group on the thread that first touched it; and 256 threads, two
# per core. Part A times every choice with the given MICM programs at 128 and
# 256 threads; part B times groups of 32, 64, and 128 with the MICM programs
# at one thread. 3 rotating rounds each, the 10,000-cell seed-2026 input, every
# run gated against MICM's final state.
#   derecho/tools/submit.sh port-session-3 128 02:30:00 pbs/port_session_3.sh OUT w4 micm-default=... micm-fastest=...
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1"); variant=$2; shift 2
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

echo "== part A: 128 and 256 threads ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out/threads-128-256" --threads 128,256 --rounds 3 \
    "${micm[@]}" \
    "g16=$bench --group-size 16" \
    "g32=$bench --group-size 32" \
    "g64=$bench --group-size 64" \
    "g79=$bench --group-size 79" \
    "g128=$bench --group-size 128" \
    "g16-static=$bench --group-size 16 --schedule static" \
    "g32-static=$bench --group-size 32 --schedule static" \
    "g79-static=$bench --group-size 79 --schedule static"

echo "== part B: 1 thread ($(date -Is))"
python3 "$repo/admin/time_entries.py" --input "$cells" --reference "$micm_final" --out "$out/threads-1" --threads 1 --rounds 3 \
    "${micm[@]}" \
    "g32=$bench --group-size 32" \
    "g64=$bench --group-size 64" \
    "g128=$bench --group-size 128"
echo "== done ($(date -Is))"
