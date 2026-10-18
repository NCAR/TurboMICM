#!/bin/bash
# Build the Derecho copy of the Fable entry, run its ported MICM test suite,
# and apply the TS1 gate to several run-time settings.
#
#   derecho/tools/verify_port.sh OUT MICM_FINAL CELLS_CSV [MAKE_VARIABLE=value ...]
#
# OUT receives build.txt (compiler, flags, build time, code size), tests.txt
# (the 680 `<id> PASS|FAIL|N/A` lines), tests-extra.txt (the entry's two extra
# checks), tests-summary.txt (counts and the id-set check against
# reference/micm_cases.txt), and gate.csv: for each setting of bin/ts1_bench_cpu
# on CELLS_CSV (10,000 cells, 100 steps of 30 s), elapsed_ms, the internal
# steps, and the gate against MICM_FINAL (reference/compare_final.py). The make
# variables (CC, OPTFLAGS, WIDTH, GROUP_SIZE, KERNEL_MODE, ...) go to every make.
# SETTINGS (environment) overrides the list of settings, one per line
# ("name|arguments"). Run inside a PBS job (derecho/tools/submit.sh) with a
# whole node, so the default setting's thread count is the node's.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1")
micm_final=$(realpath "$2")
cells_csv=$(realpath "$3")
shift 3
mkdir -p "$out"
cd "$here/fable"

echo "== build" | tee "$out/build.txt"
start=$(date +%s)
make cpu "$@" 2>&1 | tail -3
end=$(date +%s)
{
  echo "compiler: $(make -s -f - <<'EOF' "$@"
include Makefile
print:
	@echo $(CC) $$($(CC) --version | head -1)
EOF
)"
  echo "variables: $*"
  echo "cflags: $(sed -n 2p build/flags.txt)"
  echo "build seconds: $((end - start)) (library, generated kernels, and benchmark; nothing rebuilt that was up to date)"
  echo "generated/ts1.c lines: $(wc -l < generated/ts1.c)"
  echo "ts1.o text bytes: $(size build/ts1.o | awk 'NR==2 {print $1}')"
} | tee -a "$out/build.txt"

echo "== tests"
make test-cpu "$@" > "$out/tests.txt" 2> "$out/tests-extra.txt" || true
{
  echo "PASS $(grep -c ' PASS' "$out/tests.txt")  FAIL $(grep -c ' FAIL' "$out/tests.txt")  N/A $(grep -c ' N/A' "$out/tests.txt")  case lines $(grep -cE '^[^ ]+ (PASS|FAIL|N/A)' "$out/tests.txt")"
  if diff <(awk '{print $1}' "$repo/reference/micm_cases.txt" | sort) <(grep -E '^[^ ]+ (PASS|FAIL|N/A)' "$out/tests.txt" | awk '{print $1}' | sort) > /dev/null; then
    echo "ids: every id of reference/micm_cases.txt exactly once"
  else
    echo "ids: MISMATCH with reference/micm_cases.txt"
  fi
  grep ' FAIL' "$out/tests.txt" || true
  echo "extra checks:"; cat "$out/tests-extra.txt"
} | tee "$out/tests-summary.txt"

echo "== gate on $(basename "$cells_csv") against $(basename "$micm_final")"
settings=${SETTINGS:-$'default|\n1 thread|--threads 1\ngroups of 16|--group-size 16\ngroups of 64|--group-size 64\ngroups of 128|--group-size 128\ngeneric kernels|--generic\nstatic schedule|--schedule static'}
echo "setting,arguments,elapsed_ms,internal_steps,gate,max_relative_difference" > "$out/gate.csv"
while IFS='|' read -r name arguments; do
  [[ -z $name ]] && continue
  work=$(mktemp -d)
  if (cd "$work" && "$here/fable/bin/ts1_bench_cpu" $arguments 10000 100 30 "$cells_csv" > run.txt 2>&1); then
    elapsed=$(sed -n 's/^elapsed_ms=\([0-9.]*\).*/\1/p' "$work/run.txt")
    steps=$(sed -n 's/^internal_steps=\([0-9]*\).*/\1/p' "$work/run.txt")
    if compare=$(python3 "$repo/reference/compare_final.py" "$work/ts1_final.csv" "$micm_final" 2>&1); then gate=PASS; else gate=FAIL; fi
    largest=$(echo "$compare" | sed -n 's/.*largest relative difference \([^ ]*\).*/\1/p' | head -1)
  else
    elapsed= steps= gate="ERROR" largest="$(tail -1 "$work/run.txt")"
  fi
  echo "$name,\"$arguments\",$elapsed,$steps,$gate,$largest" | tee -a "$out/gate.csv"
  rm -rf "$work"
done <<< "$settings"
