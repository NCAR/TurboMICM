#!/bin/bash
# Exercise C for astra-max: swap each converted mechanism into configs/ts1.json,
# run make bench, run both benchmark programs on the mechanism's 4 cells, compare.
# The original TS1 configuration is restored on any exit.
set -u
WORK=~/NSF-NCAR/tier2-work/exercises/c/astra-max
TREE=$WORK/astra-max
CAN=~/NSF-NCAR/tier2-work/exercises/canonical
CMP=~/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py
PY=$TREE/.conda/bin/python
CONFIG=$TREE/configs/ts1.json
BACKUP=$WORK/backup/ts1.json
MECHS=${MECHS:-"chapman carbon_bond_5 carbon_bond_5_no_emissions"}

mkdir -p $WORK/backup $WORK/logs
if [ ! -f $BACKUP ]; then cp -p $CONFIG $BACKUP; fi
cmp -s $CONFIG $BACKUP || { echo "configs/ts1.json differs from the backup; stopping"; exit 1; }
shasum -a 256 $BACKUP > $WORK/backup/ts1.json.sha256

restore() {
  cp $BACKUP $CONFIG
  if cmp -s $CONFIG $BACKUP; then echo "restored configs/ts1.json"; else echo "RESTORE FAILED"; fi
}
trap restore EXIT

TIMEFORMAT='%3R'
for M in $MECHS; do
  echo "=== $M"
  sleep 2  # make 3.81 compares mtimes to the second
  cp $WORK/converted/$M.json $CONFIG
  cd $TREE
  { time make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" > $WORK/logs/make-bench-$M.log 2>&1 ; } 2> $WORK/logs/make-bench-$M.time
  status=$?
  echo "make bench exit=$status time=$(cat $WORK/logs/make-bench-$M.time) s"
  echo "$status" > $WORK/logs/make-bench-$M.exit
  mkdir -p $WORK/gen/$M/from-make $WORK/runs/$M/bin
  cp build/generated/ts1.c build/generated/ts1_specialized.c build/generated/ts1_specialized_kernels.c \
     build/generated/ts1_specialized.h build/generated/ts1_specialized.json $WORK/gen/$M/from-make/ 2>/dev/null
  if [ $status -ne 0 ]; then continue; fi
  cp bin/ts1_bench_cpu bin/ts1_bench_specialized $WORK/runs/$M/bin/
  for prog in generic specialized; do
    exe=$TREE/bin/ts1_bench_cpu; [ $prog = specialized ] && exe=$TREE/bin/ts1_bench_specialized
    rm -rf $WORK/runs/$M/$prog && mkdir -p $WORK/runs/$M/$prog && cd $WORK/runs/$M/$prog
    $exe --threads 1 4 60 60 $CAN/$M/cells.csv > stdout.txt 2> stderr.txt
    echo "$?" > exit.txt
    echo "$prog exit=$(cat exit.txt): $(tr '\n' ' ' < stdout.txt) $(cat stderr.txt)"
  done
  cd $WORK/runs/$M
  for pair in "specialized generic" "generic micm" "specialized micm"; do
    set -- $pair
    a=$WORK/runs/$M/$1/ts1_final.csv
    b=$WORK/runs/$M/$2/ts1_final.csv; [ $2 = micm ] && b=$CAN/$M/micm-run/ts1_final.csv
    if [ -f $a ] && [ -f $b ]; then
      out=$($PY $CMP $a $b); echo "$1 vs $2: $out"; echo "$out" > $WORK/runs/$M/compare-$1-vs-$2.txt
    else
      echo "$1 vs $2: missing final state"
    fi
  done
done
