#!/bin/bash
# Profile a built variant with perf: where one thread's time goes, and counters.
#
#   derecho/pbs/profile.sh OUT VARIANT [GROUP_SIZE]
#
# On the 1,000-cell seed-2026 input, 100 steps of 30 s, one thread, groups of
# GROUP_SIZE (32): `perf stat` for cycles, instructions, and cache misses, then
# `perf record` at 999 Hz with the time by function (OUT/perf-report.txt).
#   derecho/tools/submit.sh profile 8 00:30:00 pbs/profile.sh OUT w4
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
out=$(realpath -m "$1"); variant=$2; group=${3:-32}
bench=$here/work/variants/$variant/bin/ts1_bench_cpu
cells=$here/work/cells_1000_2026.csv
[[ -x $bench ]] || { echo "no $bench"; exit 1; }
mkdir -p "$out"
cd "$out"
"$here/tools/environment.sh" "$out" > /dev/null
echo "== perf stat, 1 thread, groups of $group"
perf stat -e task-clock,cycles,instructions,branches,branch-misses,L1-dcache-loads,L1-dcache-load-misses,L1-icache-load-misses,cache-misses \
    -o perf-stat.txt -- "$bench" --threads 1 --group-size "$group" 1000 100 30 "$cells" | tee run-stat.txt
cat perf-stat.txt
echo "== perf record"
perf record -F 999 -g -o perf.data -- "$bench" --threads 1 --group-size "$group" 1000 100 30 "$cells" | tee run-record.txt
perf report -i perf.data --stdio --no-children --percent-limit 0.5 2>/dev/null | grep -vE '^$|^#' | grep -E '^\s+[0-9.]+%' | head -40 | tee perf-report.txt
echo "== by kernel kind"
perf report -i perf.data --stdio --no-children --percent-limit 0 2>/dev/null | grep -E '^\s+[0-9.]+%' \
  | awk '{pct=$1; sub(/%/,"",pct); name=$NF; kind=name; if (name ~ /^ts1_/) {sub(/^ts1_/,"",kind); sub(/_[0-9]+.*$/,"",kind); sub(/\..*$/,"",kind)} else kind="library:" name; s[kind]+=pct} END {for (k in s) printf "%6.1f%%  %s\n", s[k], k}' \
  | sort -rn | head -25 | tee perf-by-kind.txt
rm -f perf.data
