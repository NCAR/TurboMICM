#!/bin/bash
# Exercise C addendum for luna-high: 64 distinct cells. Same build as the main
# runs (converted mechanism swapped in for mechanisms/ts1.json, `make bench`
# with clang -O3 -DNDEBUG); both programs on cells64.csv (warm-up + 60 x 60 s,
# 1 thread), compared with each other and with MICM's 64-cell final state.
# Also records which kernels the binaries' mechanism points to. Restores TS1.
WORK=$HOME/NSF-NCAR/tier2-work/exercises/c/luna-high
E=$WORK/luna-high
CANON=$HOME/NSF-NCAR/tier2-work/exercises/canonical
COMPARE=$HOME/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py
CONFIG=$E/mechanisms/ts1.json

restore() {
  cp "$WORK/backup/ts1.json" "$CONFIG"
  if git -C "$E" diff --quiet -- mechanisms/ts1.json; then
    echo "restored mechanisms/ts1.json (matches HEAD)"
  else
    echo "WARNING: mechanisms/ts1.json differs from HEAD after restore"
  fi
}
trap restore EXIT

for mech in ${@:-chapman carbon_bond_5 carbon_bond_5_no_emissions}; do
  OUT=$WORK/runs64/$mech
  rm -rf "$OUT"; mkdir -p "$OUT/generic" "$OUT/specialized"
  echo "===== $mech"
  cp "$WORK/converted/$mech.json" "$CONFIG"
  ( cd "$E" && /usr/bin/time -p make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" ) > "$OUT/build.log" 2>&1
  rc=$?
  echo "build_exit=$rc build_real_s=$(awk '/^real/{print $2}' "$OUT/build.log")"
  if [ $rc -ne 0 ]; then tail -30 "$OUT/build.log"; continue; fi
  cmp -s "$E/generated/ts1.c" "$WORK/gen/$mech.c" && echo "generated/ts1.c identical to gen/$mech.c"
  echo "-- kernels, specialized binary:"; python3 "$WORK/check_ops.py" "$E/bin/ts1_bench_specialized" | tee "$OUT/ops_specialized.txt"
  echo "-- kernels, generic binary:"; python3 "$WORK/check_ops.py" "$E/bin/ts1_bench_cpu" | tee "$OUT/ops_generic.txt"
  CELLS=$CANON/$mech/cells64.csv
  ( cd "$OUT/generic" && "$E/bin/ts1_bench_cpu" --threads 1 64 60 60 "$CELLS" ) > "$OUT/generic.log" 2>&1
  echo "generic_exit=$? $(tr '\n' ' ' < "$OUT/generic.log")"
  ( cd "$OUT/specialized" && "$E/bin/ts1_bench_specialized" --threads 1 64 60 60 "$CELLS" ) > "$OUT/specialized.log" 2>&1
  echo "specialized_exit=$? $(tr '\n' ' ' < "$OUT/specialized.log")"
  MICM=$CANON/$mech/micm-run64/ts1_final.csv
  echo -n "specialized_vs_generic: "; python3 "$COMPARE" "$OUT/specialized/ts1_final.csv" "$OUT/generic/ts1_final.csv"
  echo -n "generic_vs_micm:        "; python3 "$COMPARE" "$OUT/generic/ts1_final.csv" "$MICM"
  echo -n "specialized_vs_micm:    "; python3 "$COMPARE" "$OUT/specialized/ts1_final.csv" "$MICM"
  echo -n "specialized_vs_micm_g32: "; python3 "$COMPARE" "$OUT/specialized/ts1_final.csv" "$CANON/$mech/micm-run64-g32/ts1_final.csv"
  cmp -s "$OUT/generic/ts1_final.csv" "$OUT/specialized/ts1_final.csv" && echo "generic and specialized files byte-identical"
done
