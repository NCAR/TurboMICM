#!/bin/zsh
# Addendum to exercise C for astra-high: 64 distinct cells per mechanism, with
# the entry built as in the main runs (swap mechanisms/ts1.json, default
# `make bench`, GROUP_SIZE 32). Also a same-binary timing check: the driver's
# CHEM_GROUP_SIZE=16 makes every group differ from the generated stride (32),
# which selects the generic kernels. Restores and rebuilds TS1 at the end.
set -u
W=~/NSF-NCAR/tier2-work/exercises/c/astra-high
E=$W/astra-high
C=~/NSF-NCAR/tier2-work/exercises/canonical
CMP=~/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py
pause() { python3 -c 'import time; time.sleep(1.1)'; }   # make 3.81: 1-s mtimes
cmp -s $W/ts1.json.orig $E/mechanisms/ts1.json || { echo "ts1.json not at its original; stop"; exit 1; }

run_bench() {  # $1 program, $2 directory, $3 cells.csv, [$4 CHEM_GROUP_SIZE]
  rm -rf $2; mkdir -p $2
  if [[ -n ${4:-} ]]; then
    (cd $2 && CHEM_GROUP_SIZE=$4 $E/bin/$1 --threads 1 64 60 60 $3 > stdout.txt 2> stderr.txt; echo $? > exit.txt)
  else
    (cd $2 && $E/bin/$1 --threads 1 64 60 60 $3 > stdout.txt 2> stderr.txt; echo $? > exit.txt)
  fi
}
compare() { python3 $CMP $1 $2 > $3 2>&1; echo "exit $?" >> $3; }

mechs=("$@"); (( $#mechs )) || mechs=(chapman carbon_bond_5 carbon_bond_5_no_emissions)
for m in $mechs; do
  R=$W/runs64/$m; rm -rf $R; mkdir -p $R
  cells=$C/$m/cells64.csv; micm=$C/$m/micm-run64/ts1_final.csv
  pause; cp $W/configs/$m/ts1.json $E/mechanisms/ts1.json   # the temporary swap
  (cd $E && make bench CC=clang 'OPTFLAGS=-O3 -DNDEBUG' > $R/build.log 2>&1); echo "exit $?" > $R/build.status
  cmp -s $E/build/specialized/ts1.c $W/runs/$m/specialized.c && echo "generated C identical to main run" >> $R/build.status
  run_bench ts1_bench_cpu $R/generic $cells
  run_bench ts1_bench_specialized $R/specialized $cells
  compare $R/specialized/ts1_final.csv $R/generic/ts1_final.csv $R/specialized_vs_generic.txt
  compare $R/generic/ts1_final.csv $micm $R/generic_vs_micm.txt
  compare $R/specialized/ts1_final.csv $micm $R/specialized_vs_micm.txt
  compare $R/specialized/ts1_final.csv $C/$m/micm-run64-g32/ts1_final.csv $R/specialized_vs_micm_g32.txt
  # Timing evidence, 5 repetitions each (groups of 32 = stride, and of 16).
  for rep in 1 2 3 4 5; do
    run_bench ts1_bench_cpu $R/t/generic-g32-$rep $cells
    run_bench ts1_bench_specialized $R/t/specialized-g32-$rep $cells
    run_bench ts1_bench_cpu $R/t/generic-g16-$rep $cells 16
    run_bench ts1_bench_specialized $R/t/specialized-g16-$rep $cells 16
  done
done

pause; cp $W/ts1.json.orig $E/mechanisms/ts1.json            # restore and rebuild TS1
(cd $E && make bench CC=clang 'OPTFLAGS=-O3 -DNDEBUG' > $W/runs64/restore-build.log 2>&1); echo "restore build exit $?"
git -C $E diff --quiet -- mechanisms/ts1.json && echo "ts1.json restored"
