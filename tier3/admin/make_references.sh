#!/bin/zsh
# Make the high-accuracy reference solution for each Tier 3 scenario.
#
#   admin/make_references.sh MICM_SCENARIO THREADS SCENARIO_DIR...
#
# MICM_SCENARIO is scenarios/micm_scenario.cpp built with the Tier 1 flags:
#   clang++ -std=c++20 -O3 -DNDEBUG -I micm/include -I micm/test/benchmark \
#       scenarios/micm_scenario.cpp -o micm_scenario
#
# For each scenario directory, this runs three MICM solutions:
#   default  MICM's three-stage Rosenbrock with its default tolerances (the baseline)
#   tight    the six-stage Rosenbrock at relative 1e-9, absolute 1e-22 mol m-3
#   tighter  the six-stage Rosenbrock at relative 1e-10, absolute 1e-23 mol m-3
# It scores tight against tighter (convergence) and default against tighter
# (MICM's own accuracy), keeps tighter as reference_final.csv and, for scenarios
# with snapshots, reference_snapshots.csv, and writes the scores to
# references.log in the scenario directory.
set -u
here=${0:A:h}
micm=${1:A}
threads=$2
shift 2
for dir in "$@"; do
  dir=${dir:A}
  scenario=$dir/scenario.txt
  work=$(mktemp -d)
  for run in default tight tighter; do
    mkdir -p $work/$run
    case $run in
      default) flags=() ;;
      tight) flags=(--method rosenbrock6dae --rtol 1e-9 --atol 1e-22 --max-steps 10000000) ;;
      tighter) flags=(--method rosenbrock6dae --rtol 1e-10 --atol 1e-23 --max-steps 10000000) ;;
    esac
    if ! (cd $work/$run && $micm --threads $threads $flags $scenario > run.log 2>&1); then
      if [ $run = default ]; then
        echo "FAILED" >> $work/$run/run.log  # MICM's default failing is a finding; the reference still gets made
      else
        echo "${dir:t} $run failed:"; cat $work/$run/run.log; continue 2
      fi
    fi
  done
  result=ts1_final.csv
  [ -f $work/tighter/ts1_snapshots.csv ] && result=ts1_snapshots.csv
  {
    echo "# ${dir:t}: $(head -1 $scenario | sed 's/^# //')"
    echo "# made with ${micm:t} --threads $threads"
    for run in default tight tighter; do echo "$run: $(tr '\n' ' ' < $work/$run/run.log)"; done
    echo "convergence (tight against tighter): $(python3 $here/../scenarios/score.py $work/tight/$result $work/tighter/$result | head -1)"
    if [ -f $work/default/$result ]; then
      echo "MICM default against tighter: $(python3 $here/../scenarios/score.py $work/default/$result $work/tighter/$result | head -1)"
    else
      echo "MICM default against tighter: no result (the default run failed)"
    fi
  } > $dir/references.log
  cp $work/tighter/ts1_final.csv $dir/reference_final.csv
  [ -f $work/tighter/ts1_snapshots.csv ] && cp $work/tighter/ts1_snapshots.csv $dir/reference_snapshots.csv
  cat $dir/references.log
  rm -rf $work
done
