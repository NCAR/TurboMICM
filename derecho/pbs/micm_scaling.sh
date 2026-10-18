#!/bin/bash
# How MICM's reference scales with threads on a Derecho node, and why.
#
#   derecho/pbs/micm_scaling.sh OUT PROGRAM [STEPS]
#
# Runs PROGRAM (a MICM TS1 reference binary, e.g. work/micm-configs/bin/micm-standard-in-place-mozart)
# on the 10,000-cell seed-2026 input for STEPS (10) steps at 1, 2, 4, 8, 16,
# 32, 64, 128, and 256 threads, each under /usr/bin/time -v (user and system
# time, page faults, context switches) and perf stat (page faults, migrations,
# cycles, instructions). Then the same at 128 threads with the variables that
# change glibc's allocator: MALLOC_ARENA_MAX=1 and 128, MALLOC_TOP_PAD_ and
# MALLOC_MMAP_THRESHOLD_ large (no trimming, no mmap), and with transparent
# huge pages disabled for the process (prctl via `setarch -R`? no: via
# GLIBC_TUNABLES=glibc.malloc.hugetlb=0 and THP madvise unavailable per
# process, so recorded only). OUT/scaling.txt has the summary lines.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
out=$(realpath -m "$1"); program=$(realpath "$2"); steps=${3:-10}
cells=$here/work/cells_2026.csv
mkdir -p "$out"
cd "$out"
"$here/tools/environment.sh" "$out" > /dev/null
: > scaling.txt
run() {
  local label=$1; shift
  local threads=$1; shift
  /usr/bin/time -v -o time-$label.txt env "$@" "$program" --threads "$threads" 10000 "$steps" 30 "$cells" > run-$label.txt 2>&1 || true
  local elapsed user sys faults ctx
  elapsed=$(sed -n 's/^elapsed_ms=\([0-9.]*\).*/\1/p' run-$label.txt)
  user=$(sed -n 's/.*User time (seconds): //p' time-$label.txt)
  sys=$(sed -n 's/.*System time (seconds): //p' time-$label.txt)
  faults=$(sed -n 's/.*Minor (reclaiming a frame) page faults: //p' time-$label.txt)
  ctx=$(sed -n 's/.*Involuntary context switches: //p' time-$label.txt)
  printf '%-28s threads %3d  elapsed_ms %10s  user_s %8s  sys_s %8s  minor_faults %10s  invol_ctx %8s\n' "$label" "$threads" "$elapsed" "$user" "$sys" "$faults" "$ctx" | tee -a scaling.txt
}
echo "== $program, $steps steps" | tee -a scaling.txt
for t in 1 2 4 8 16 32 64 128 256; do
  run "t$t" "$t"
done
echo "== allocator variables at 128 threads" | tee -a scaling.txt
run "arena1" 128 MALLOC_ARENA_MAX=1
run "arena128" 128 MALLOC_ARENA_MAX=128
run "notrim-nommap" 128 MALLOC_TOP_PAD_=1073741824 MALLOC_TRIM_THRESHOLD_=4294967295 MALLOC_MMAP_THRESHOLD_=4294967295
run "tunables-nommap" 128 GLIBC_TUNABLES=glibc.malloc.mmap_threshold=4294967295:glibc.malloc.trim_threshold=4294967295:glibc.malloc.top_pad=1073741824
echo "== perf stat at 128 threads" | tee -a scaling.txt
perf stat -e task-clock,page-faults,context-switches,cpu-migrations,cycles,instructions -o perf-t128.txt -- "$program" --threads 128 10000 "$steps" 30 "$cells" > run-perf-t128.txt 2>&1 || true
cat perf-t128.txt | tee -a scaling.txt
echo "== perf record at 128 threads (top symbols)" | tee -a scaling.txt
perf record -F 499 -g -o perf.data -- "$program" --threads 128 10000 "$steps" 30 "$cells" > run-perf-record.txt 2>&1 || true
perf report -i perf.data --stdio --no-children --percent-limit 1 2>/dev/null | grep -E '^\s+[0-9.]+%' | head -30 | tee -a scaling.txt
rm -f perf.data
