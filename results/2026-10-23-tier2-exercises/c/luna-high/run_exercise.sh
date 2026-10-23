#!/bin/bash
# Exercise C for luna-high: swap each converted mechanism in for
# mechanisms/ts1.json, time `make bench`, run both benchmark programs on the
# mechanism's 4 cells (warm-up + 60 x 60 s), compare. Restores TS1 on exit.
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

for m in "${@:-chapman carbon_bond_5 carbon_bond_5_no_emissions}"; do
  for mech in $m; do
    OUT=$WORK/runs/$mech
    rm -rf "$OUT"; mkdir -p "$OUT/generic" "$OUT/specialized"
    echo "===== $mech"
    cp "$WORK/converted/$mech.json" "$CONFIG"
    ( cd "$E" && /usr/bin/time -p make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" ) > "$OUT/build.log" 2>&1
    rc=$?
    echo "build_exit=$rc build_real_s=$(awk '/^real/{print $2}' "$OUT/build.log")"
    if [ $rc -ne 0 ]; then tail -30 "$OUT/build.log"; continue; fi
    cmp -s "$E/generated/ts1.c" "$WORK/gen/$mech.c" && echo "generated/ts1.c identical to gen/$mech.c" \
      || echo "NOTE: generated/ts1.c differs from gen/$mech.c"
    CELLS=$CANON/$mech/cells.csv
    ( cd "$OUT/generic" && "$E/bin/ts1_bench_cpu" --threads 1 4 60 60 "$CELLS" ) > "$OUT/generic.log" 2>&1
    echo "generic_exit=$? $(tr '\n' ' ' < "$OUT/generic.log")"
    ( cd "$OUT/specialized" && "$E/bin/ts1_bench_specialized" --threads 1 4 60 60 "$CELLS" ) > "$OUT/specialized.log" 2>&1
    echo "specialized_exit=$? $(tr '\n' ' ' < "$OUT/specialized.log")"
    MICM=$CANON/$mech/micm-run/ts1_final.csv
    echo -n "specialized_vs_generic: "; python3 "$COMPARE" "$OUT/specialized/ts1_final.csv" "$OUT/generic/ts1_final.csv"
    echo -n "generic_vs_micm:        "; python3 "$COMPARE" "$OUT/generic/ts1_final.csv" "$MICM"
    echo -n "specialized_vs_micm:    "; python3 "$COMPARE" "$OUT/specialized/ts1_final.csv" "$MICM"
  done
done
