#!/bin/bash
# Exercise C for sol-high: swap configs/ts1.json for a converted mechanism,
# rebuild with the entry's Makefile, run both benchmark programs on the
# mechanism's 4 cells, compare. Usage: run_exercise.sh PHASE MECH...
#   protocol  : make bench (defaults), run generic + specialized, compare
#   width4    : make bench SPECIALIZE_WIDTH=4 (generated kernels at 4 cells), run specialized, compare
#   clean     : make clean, then time make bench (defaults)
#   cells64   : make bench (defaults), run both programs on the 64 distinct cells of
#               cells64.csv (two full 32-cell groups), compare, and run the fallback probe
#   restore   : put the original TS1 configuration back, rebuild, TS1 sanity check
set -u
W=$HOME/NSF-NCAR/tier2-work/exercises/c/sol-high
E=$W/sol-high
C=$HOME/NSF-NCAR/tier2-work/exercises/canonical
CMP=$HOME/NSF-NCAR/ODE/benchmarks/tier2/reference/compare_final.py
mkdir -p "$W/backup" "$W/logs" "$W/runs"
if [ ! -f "$W/backup/ts1.json" ]; then
  cp -p "$E/configs/ts1.json" "$W/backup/ts1.json"
fi

swap() {  # put the converted mechanism in place of TS1's configuration
  cp "$W/converted/$1/ts1.json" "$E/configs/ts1.json"
}

build() {  # build LOG [extra make variables]
  local log=$1; shift
  (cd "$E" && /usr/bin/time -p make bench CC=clang OPTFLAGS="-O3 -DNDEBUG" "$@") > "$log" 2>&1
  local status=$?
  echo "build exit $status real $(awk '/^real/{print $2}' "$log")  ($log)"
  return $status
}

run() {  # run MECH LABEL PROGRAM
  local d=$W/runs/$1/$2
  rm -rf "$d" && mkdir -p "$d"
  (cd "$d" && "$E/bin/$3" --threads 1 4 60 60 "$C/$1/cells.csv") > "$d/run.log" 2>&1
  echo "$1 $2 ($3) exit $?: $(tr '\n' ' ' < "$d/run.log")"
}

compare() {  # compare MECH A B (labels under runs/MECH, or micm, micm64, micm64g32)
  local a=$W/runs/$1/$2/ts1_final.csv b=$W/runs/$1/$3/ts1_final.csv
  [ "$3" = micm ] && b=$C/$1/micm-run/ts1_final.csv
  [ "$3" = micm64 ] && b=$C/$1/micm-run64/ts1_final.csv
  [ "$3" = micm64g32 ] && b=$C/$1/micm-run64-g32/ts1_final.csv
  echo "$1 $2 vs $3: $(python3 "$CMP" "$a" "$b" 2>&1)"
}

run64() {  # run64 MECH LABEL PROGRAM: the 64 distinct cells of cells64.csv
  local d=$W/runs/$1/$2
  rm -rf "$d" && mkdir -p "$d"
  (cd "$d" && "$E/bin/$3" --threads 1 64 60 60 "$C/$1/cells64.csv") > "$d/run.log" 2>&1
  echo "$1 $2 ($3, 64 cells) exit $?: $(tr '\n' ' ' < "$d/run.log")"
}

# Evidence probe (not the measured programs): the current build/ts1_specialized.c,
# compiled with probe/probe.h forced in so that its generic width fallbacks are
# counted, linked like bin/ts1_bench_specialized, run at 4 and at 64 cells.
LIBSRC="src/mechanism.c src/state.c src/rates.c src/process.c src/linear.c src/parameters.c src/rosenbrock.c src/backward_euler.c src/constraints.c src/array.c src/sparse.c src/system.c src/jacobian.c src/builder.c"
LIBOMP=$(brew --prefix libomp 2>/dev/null)
probe() {  # probe MECH
  local p=$W/probe/$1
  rm -rf "$p" && mkdir -p "$p/run4" "$p/run64"
  local cf="-std=c11 -Wall -Wextra -Wpedantic -Iinclude -Isrc -DCHEM_GROUP_SIZE=32 -Xpreprocessor -fopenmp -isystem $LIBOMP/include -O3 -DNDEBUG"
  cp "$E/build/ts1_specialized.c" "$p/ts1_specialized.c"
  (cd "$E" && export SDKROOT=$(xcrun --sdk macosx --show-sdk-path) && \
    clang $cf -fPIC -include "$W/probe/probe.h" -c build/ts1_specialized.c -o "$p/ts1_specialized_probe.o" && \
    clang $cf -DCHEM_SPECIALIZED $LIBSRC src/bench.c "$p/ts1_specialized_probe.o" "$W/probe/probe.c" \
      -o "$p/bench_specialized_probe" -L"$LIBOMP/lib" -lomp -lm) > "$p/build.log" 2>&1
  echo "$1 probe build exit $? (warnings: $(command grep -c 'warning:' "$p/build.log"))"
  (cd "$p/run4" && "$p/bench_specialized_probe" --threads 1 4 60 60 "$C/$1/cells.csv") > "$p/run4/run.log" 2>&1
  echo "$1 probe 4 cells exit $?: $(tr '\n' ' ' < "$p/run4/run.log")"
  cmp -s "$p/run4/ts1_final.csv" "$W/runs/$1/specialized/ts1_final.csv" \
    && echo "$1 probe 4-cell output byte-identical to bin/ts1_bench_specialized's" || echo "$1 probe 4-cell output DIFFERS"
  (cd "$p/run64" && "$p/bench_specialized_probe" --threads 1 64 60 60 "$C/$1/cells64.csv") > "$p/run64/run.log" 2>&1
  echo "$1 probe 64 cells exit $?: $(tr '\n' ' ' < "$p/run64/run.log")"
  cmp -s "$p/run64/ts1_final.csv" "$W/runs/$1/specialized-64/ts1_final.csv" \
    && echo "$1 probe 64-cell output byte-identical to bin/ts1_bench_specialized's" || echo "$1 probe 64-cell output DIFFERS"
}

phase=$1; shift
for m in "$@"; do
  case $phase in
    protocol)
      swap "$m"
      build "$W/logs/$m-make.log"
      cmp "$E/build/ts1_specialized.c" "$W/gen/$m/ts1_specialized.c" && echo "$m build/ts1_specialized.c identical to gen/$m"
      echo "$m build/ts1_specialized.c: $(wc -l < "$E/build/ts1_specialized.c") lines, $(wc -c < "$E/build/ts1_specialized.c") bytes"
      run "$m" generic ts1_bench_cpu
      run "$m" specialized ts1_bench_specialized
      compare "$m" specialized generic
      compare "$m" generic micm
      compare "$m" specialized micm
      ;;
    width4)
      swap "$m"
      build "$W/logs/$m-make-width4.log" SPECIALIZE_WIDTH=4
      run "$m" specialized-width4 ts1_bench_specialized
      compare "$m" specialized-width4 generic
      compare "$m" specialized-width4 micm
      compare "$m" specialized-width4 specialized
      ;;
    clean)
      swap "$m"
      (cd "$E" && make clean) > /dev/null 2>&1
      build "$W/logs/$m-make-clean.log"
      ;;
    cells64)  # the Makefile's defaults, as in the protocol phase
      swap "$m"
      build "$W/logs/$m-make-cells64.log"
      cmp "$E/build/ts1_specialized.c" "$W/gen/$m/ts1_specialized.c" && echo "$m build/ts1_specialized.c identical to gen/$m (default width 32)"
      echo "$m width guards in build/ts1_specialized.c: $(command grep -c 'cells != 32)' "$E/build/ts1_specialized.c")"
      run64 "$m" generic-64 ts1_bench_cpu
      run64 "$m" specialized-64 ts1_bench_specialized
      compare "$m" specialized-64 generic-64
      compare "$m" generic-64 micm64
      compare "$m" specialized-64 micm64
      compare "$m" specialized-64 micm64g32
      probe "$m"
      ;;
  esac
done
if [ "$phase" = restore ]; then
  cp -p "$W/backup/ts1.json" "$E/configs/ts1.json"
  touch "$E/configs/ts1.json"   # newer than the generated files, so make regenerates TS1
  build "$W/logs/ts1-make-restore.log"
  mkdir -p "$W/runs/ts1" && rm -rf "$W/runs/ts1/generic" "$W/runs/ts1/specialized"
  for p in generic:ts1_bench_cpu specialized:ts1_bench_specialized; do
    d=$W/runs/ts1/${p%%:*}; mkdir -p "$d"
    (cd "$d" && "$E/bin/${p#*:}" --threads 1 100 100 30 "$W/reference/ts1_cells_100_seed1.csv") > "$d/run.log" 2>&1
    echo "ts1 ${p%%:*} exit $?: $(tr '\n' ' ' < "$d/run.log")"
    echo "ts1 ${p%%:*} vs MICM: $(python3 "$CMP" "$d/ts1_final.csv" "$W/reference/micm_ts1_final_100_seed1.csv")"
  done
  git -C "$E" status --porcelain --untracked-files=no
  echo "tracked changes listed above (none expected)"
fi
