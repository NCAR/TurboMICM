#!/bin/bash
# Addendum: 64 distinct cells (two full 32-cell groups) per mechanism, the entry
# built as in the main runs (configs/ts1.json swapped, make bench with the
# Makefile's defaults). Both programs, each in its own empty directory, compared
# with each other and with MICM's 64-cell final state (micm-run64).
# The original TS1 configuration is restored on any exit.
set -u
WORK=~/NSF-NCAR/tier2-work/exercises/c/astra-max
TREE=$WORK/astra-max
CAN=~/NSF-NCAR/tier2-work/exercises/canonical
CMP=~/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py
PY=$TREE/.conda/bin/python
CONFIG=$TREE/configs/ts1.json
BACKUP=$WORK/backup/ts1.json
cmp -s $CONFIG $BACKUP || { echo "configs/ts1.json differs from the backup; stopping"; exit 1; }
restore() {
  cp $BACKUP $CONFIG
  if cmp -s $CONFIG $BACKUP; then echo "restored configs/ts1.json"; else echo "RESTORE FAILED"; fi
}
trap restore EXIT
TARGETS=('configure.py configs/ts1.json build/generated/ts1.c' '-o build/ts1.o' '-o bin/ts1_bench_cpu'
  'specialize.py configs/ts1.json build/generated/ts1_specialized.c' '-o build/ts1_specialized.o'
  '-o build/ts1_specialized_kernels.o' '-o bin/ts1_bench_specialized' '-o bin/ts1_scenario_specialized')
for M in chapman carbon_bond_5 carbon_bond_5_no_emissions; do
  D=$WORK/runs/$M-cells64
  mkdir -p $D/bin
  sleep 2  # make 3.81 compares mtimes to the second
  cp $WORK/converted/$M.json $CONFIG
  cd $TREE
  make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" > $WORK/logs/make-bench-$M-cells64.log 2>&1
  status=$?
  n=0; for t in "${TARGETS[@]}"; do command grep -qF -- "$t" $WORK/logs/make-bench-$M-cells64.log && n=$((n+1)); done
  echo "=== $M: make bench exit=$status, $n/8 targets rebuilt"
  [ $status -ne 0 ] && continue
  cp bin/ts1_bench_cpu bin/ts1_bench_specialized $D/bin/
  cp build/generated/ts1_specialized.c build/generated/ts1_specialized_kernels.c $D/bin/
  for prog in generic specialized; do
    exe=$TREE/bin/ts1_bench_cpu; [ $prog = specialized ] && exe=$TREE/bin/ts1_bench_specialized
    rm -rf $D/$prog && mkdir -p $D/$prog && cd $D/$prog
    $exe --threads 1 64 60 60 $CAN/$M/cells64.csv > stdout.txt 2> stderr.txt
    echo "$?" > exit.txt
    echo "$prog exit=$(cat exit.txt): $(tr '\n' ' ' < stdout.txt) $(cat stderr.txt)"
  done
  for pair in "specialized generic" "generic micm" "specialized micm"; do
    set -- $pair
    a=$D/$1/ts1_final.csv
    b=$D/$2/ts1_final.csv; [ $2 = micm ] && b=$CAN/$M/micm-run64/ts1_final.csv
    out=$($PY $CMP $a $b); echo "$1 vs $2: $out"; echo "$out" > $D/compare-$1-vs-$2.txt
  done
done
