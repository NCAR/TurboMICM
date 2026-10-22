#!/bin/zsh
# Make MICM's results at the gate's settings, and a reference solution, for each
# Tier 2 scenario.
#
#   admin/make_scenario_data.sh MICM_SCENARIO SCENARIO_DIR...
#
# MICM_SCENARIO is scenarios/micm_scenario.cpp built with the Tier 1 flags:
#   clang++ -std=c++20 -O3 -DNDEBUG -I micm/include -I micm/test/benchmark \
#       scenarios/micm_scenario.cpp -o micm_scenario
#
# For each scenario directory, this runs four MICM solutions:
#   gate     MICM's three-stage Rosenbrock at relative 1e-6 (its default) and
#            absolute 1e-15 mol m-3, with 32-cell groups: the gate's settings,
#            kept as micm_final.csv and micm_snapshots.csv
#   default  the same with MICM's default absolute tolerance, 1e-3 mol m-3
#   tight    the six-stage Rosenbrock at relative 1e-9, absolute 1e-22 mol m-3
#   tighter  the six-stage Rosenbrock at relative 1e-10, absolute 1e-23 mol m-3,
#            kept as reference_final.csv and reference_snapshots.csv
# The references run with 10 threads and MICM's default grouping, one group of
# 100 cells per thread, as Tier 3's did; the gate run uses all cores, which with
# fixed 32-cell groups does not change its result. It scores tight against
# tighter (the reference's convergence), and gate and default against tighter
# (MICM's accuracy), and writes the runs, the scores, and the SHA-256 of every
# file it keeps to micm.log in the scenario directory.
set -u
here=${0:A:h}
micm=${1:A}
shift
score=$here/../scenarios/score.py
for dir in "$@"; do
  dir=${dir:A}
  scenario=$dir/scenario.txt
  work=$(mktemp -d)
  for run in gate default tight tighter; do
    mkdir -p $work/$run
    case $run in
      gate) flags=(--group-size 32 --atol 1e-15) ;;
      default) flags=(--group-size 32) ;;
      tight) flags=(--threads 10 --method rosenbrock6dae --rtol 1e-9 --atol 1e-22 --max-steps 10000000) ;;
      tighter) flags=(--threads 10 --method rosenbrock6dae --rtol 1e-10 --atol 1e-23 --max-steps 10000000) ;;
    esac
    print -r -- "\$ micm_scenario $flags scenario.txt" > $work/$run/run.log
    if ! (cd $work/$run && $micm $flags $scenario >> run.log 2>&1); then
      if [ $run = default ]; then
        echo "FAILED" >> $work/$run/run.log  # MICM's default failing is a finding; the rest still gets made
      else
        echo "${dir:t} $run failed:"; cat $work/$run/run.log; continue 2
      fi
    fi
  done
  result=ts1_final.csv
  [ -f $work/tighter/ts1_snapshots.csv ] && result=ts1_snapshots.csv
  cp $work/gate/ts1_final.csv $dir/micm_final.csv
  cp $work/tighter/ts1_final.csv $dir/reference_final.csv
  if [ -f $work/tighter/ts1_snapshots.csv ]; then
    cp $work/gate/ts1_snapshots.csv $dir/micm_snapshots.csv
    cp $work/tighter/ts1_snapshots.csv $dir/reference_snapshots.csv
  fi
  {
    echo "# ${dir:t}: $(head -1 $scenario | sed 's/^# //')"
    echo "# made $(date '+%Y-%m-%d %H:%M') on $(sysctl -n machdep.cpu.brand_string), with ${micm:t}"
    for run in gate default tight tighter; do echo "$run: $(tr '\n' ' ' < $work/$run/run.log)"; done
    echo "convergence (tight against tighter): $(python3 $score $work/tight/$result $work/tighter/$result | head -1)"
    echo "gate run against tighter: $(python3 $score $work/gate/$result $work/tighter/$result | head -1)"
    if [ -f $work/default/$result ]; then
      echo "MICM default against tighter: $(python3 $score $work/default/$result $work/tighter/$result | head -1)"
    else
      echo "MICM default against tighter: no result (the default run failed)"
    fi
    echo "files:"
    (cd $dir && shasum -a 256 cells.csv photolysis.csv scenario.txt micm_*.csv reference_*.csv | sed 's/^/  /')
  } > $dir/micm.log
  cat $dir/micm.log
  rm -rf $work
done
