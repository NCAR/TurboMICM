#!/bin/bash
# Exercise C for astra-xhigh: swap each converted mechanism in for configs/ts1.json,
# make bench, run both benchmark programs, compare; restore TS1 on exit.
# Usage: run_exercise.sh MECH...
set -u
W=$HOME/NSF-NCAR/tier2-work/exercises/c/astra-xhigh
T=$W/astra-xhigh
CAN=$HOME/NSF-NCAR/tier2-work/exercises/canonical
CMP=$HOME/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py
PY=$T/.conda/bin/python
SPECIAL_FLAGS="--chemistry-batch 8 --algebra-batch 8 --function-lines 256 --tile-size 4"
mkdir -p $W/backup $W/logs $W/gen $W/runs

# Back up the TS1 configuration once, checked against the committed blob.
if [ ! -f $W/backup/ts1.json ]; then cp -p $T/configs/ts1.json $W/backup/ts1.json; fi
[ "$(git -C $T hash-object $W/backup/ts1.json)" = "$(git -C $T rev-parse HEAD:configs/ts1.json)" ] \
  || { echo "backup does not match HEAD's configs/ts1.json"; exit 1; }
restore() { cp $W/backup/ts1.json $T/configs/ts1.json; echo "restored configs/ts1.json"; }
trap restore EXIT

seconds() { awk '/^real/ {print $2}' "$1"; }

for m in "$@"; do
  echo "=== $m"
  G=$W/gen/$m; mkdir -p $G
  if [ "$m" = ts1 ]; then config=$W/backup/ts1.json; else config=$W/converted/$m/ts1.json; fi
  # Generators, standalone, with the Makefile's arguments (output outside the entry).
  (cd $T && /usr/bin/time -p $PY tools/specialize.py $config $G/ts1_specialized.c --symbol ts1_specialized $SPECIAL_FLAGS) \
    > $W/logs/$m-specialize.log 2>&1
  echo "specialize.py exit $? real $(seconds $W/logs/$m-specialize.log) s"
  (cd $T && /usr/bin/time -p $PY tools/generate.py $config $G/ts1.c --symbol ts1_mechanism) \
    > $W/logs/$m-generate.log 2>&1
  echo "generate.py exit $? real $(seconds $W/logs/$m-generate.log) s"
  [ "$m" = ts1 ] && continue
  # Swap, build, run.
  cp $config $T/configs/ts1.json
  (cd $T && /usr/bin/time -p make bench CC=clang OPTFLAGS="-O3 -DNDEBUG") > $W/logs/$m-build.log 2>&1
  status=$?
  echo "make bench exit $status real $(seconds $W/logs/$m-build.log) s"
  if [ $status -ne 0 ]; then tail -20 $W/logs/$m-build.log; continue; fi
  for f in ts1_specialized.c ts1.c; do
    cmp -s $T/generated/$f $G/$f && same=identical || same=DIFFERENT
    echo "generated/$f: $(wc -l < $T/generated/$f) lines, $(wc -c < $T/generated/$f) bytes; standalone copy $same"
  done
  for p in ts1_bench_cpu ts1_bench_specialized; do
    R=$W/runs/$m/$p; rm -rf $R; mkdir -p $R
    (cd $R && $T/bin/$p --threads 1 4 60 60 $CAN/$m/cells.csv) > $W/logs/$m-$p.log 2>&1
    echo "$p exit $?: $(tr '\n' ' ' < $W/logs/$m-$p.log)"
  done
  S=$W/runs/$m/ts1_bench_specialized/ts1_final.csv; C=$W/runs/$m/ts1_bench_cpu/ts1_final.csv; M=$CAN/$m/micm-run/ts1_final.csv
  echo "specialized vs generic: $(python3 $CMP $S $C 2>&1)"
  echo "generic vs MICM:        $(python3 $CMP $C $M 2>&1)"
  echo "specialized vs MICM:    $(python3 $CMP $S $M 2>&1)"
done
