#!/bin/zsh
# Exercise C for astra-high: swap each converted mechanism in for
# mechanisms/ts1.json, run `make bench`, run both benchmark programs on the
# mechanism's 4 cells (1 warm-up step, then 60 steps of 60 s), compare, and
# restore the TS1 file. Writes only under $W.
set -u
W=~/NSF-NCAR/tier2-work/exercises/c/astra-high
E=$W/astra-high
C=~/NSF-NCAR/tier2-work/exercises/canonical
CMP=~/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py
MAKEFLAGS_EX=(CC=clang 'OPTFLAGS=-O3 -DNDEBUG')
now() { python3 -c 'import time; print(f"{time.time():.3f}")'; }
elapsed() { python3 -c "print(f'{$2-$1:.2f}')"; }

# Back up the TS1 configuration once; it must match the committed file.
if [[ ! -f $W/ts1.json.orig ]]; then
  cp $E/mechanisms/ts1.json $W/ts1.json.orig
fi
git -C $E diff --quiet -- mechanisms/ts1.json || { echo "ts1.json already modified; stop"; exit 1; }
cmp -s $W/ts1.json.orig $E/mechanisms/ts1.json || { echo "backup differs; stop"; exit 1; }

run_bench() {  # $1 program, $2 directory, $3 cells.csv
  rm -rf $2; mkdir -p $2
  (cd $2 && $E/bin/$1 --threads 1 4 60 60 $3 > stdout.txt 2> stderr.txt; echo $? > exit.txt)
}
compare() {  # $1 yours, $2 other, $3 log
  if [[ -f $1 && -f $2 ]]; then python3 $CMP $1 $2 > $3 2>&1; echo "exit $?" >> $3; else echo "missing input" > $3; fi
}

mechs=("$@"); (( $#mechs )) || mechs=(chapman carbon_bond_5 carbon_bond_5_no_emissions)
for m in $mechs; do
  R=$W/runs/$m; rm -rf $R; mkdir -p $R
  cells=$C/$m/cells.csv; micm=$C/$m/micm-run/ts1_final.csv
  cp $W/configs/$m/ts1.json $E/mechanisms/ts1.json        # the temporary swap
  # 1. Default build (GROUP_SIZE 32: generated stride 32).
  t0=$(now); (cd $E && make bench $MAKEFLAGS_EX > $R/build.log 2>&1); ec=$?; t1=$(now)
  echo "exit $ec; $(elapsed $t0 $t1) s" > $R/build.status
  cp $E/build/specialized/ts1.c $R/specialized.c 2>/dev/null
  cp $E/build/ts1.c $R/builder.c 2>/dev/null
  if [[ $ec = 0 ]]; then
    run_bench ts1_bench_cpu $R/generic $cells
    run_bench ts1_bench_specialized $R/specialized $cells
    compare $R/specialized/ts1_final.csv $R/generic/ts1_final.csv $R/specialized_vs_generic.txt
    compare $R/generic/ts1_final.csv $micm $R/generic_vs_micm.txt
    compare $R/specialized/ts1_final.csv $micm $R/specialized_vs_micm.txt
  fi
  # 2. Same problem with the generated stride equal to the 4-cell group
  #    (the Makefile's documented GROUP_SIZE variable), so the generated
  #    kernels run instead of the generic fallback for groups != stride.
  t0=$(now); (cd $E && make bench $MAKEFLAGS_EX GROUP_SIZE=4 > $R/build-g4.log 2>&1); ec=$?; t1=$(now)
  echo "exit $ec; $(elapsed $t0 $t1) s" > $R/build-g4.status
  if [[ $ec = 0 ]]; then
    run_bench ts1_bench_specialized $R/specialized-g4 $cells
    compare $R/specialized-g4/ts1_final.csv $R/generic/ts1_final.csv $R/specialized_g4_vs_generic.txt
    compare $R/specialized-g4/ts1_final.csv $micm $R/specialized_g4_vs_micm.txt
  fi
done

cp $W/ts1.json.orig $E/mechanisms/ts1.json                 # restore
cmp -s $W/ts1.json.orig $E/mechanisms/ts1.json && git -C $E diff --quiet -- mechanisms/ts1.json && echo "ts1.json restored"
