#!/bin/sh
# Reproduce the disclosed day-night gates without writing into reference inputs.
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
label=${1:-local}
threads=${2:-18}
for name in diurnal-30min diurnal-10min diurnal-5min polluted-10min clean-10min; do
  out="$repo/results/$label/$name"
  mkdir -p "$out"
  (
    cd "$out"
    "$repo/bin/ts1_scenario_specialized" --threads "$threads" --group-size 32 --rtol 1e-6 --atol 1e-15 "$repo/../scenarios/$name/scenario.txt" > run.log
    "$repo/.conda/bin/python" "$repo/../scenarios/score.py" ts1_snapshots.csv "$repo/../scenarios/$name/micm_snapshots.csv" --level 10 --json > micm-score.json
    "$repo/.conda/bin/python" "$repo/../scenarios/score.py" ts1_snapshots.csv "$repo/../scenarios/$name/reference_snapshots.csv" --json > reference-score.json
    "$repo/.conda/bin/python" "$repo/../scenarios/score.py" ts1_final.csv "$repo/../scenarios/$name/micm_final.csv" --level 10 --json > final-score.json
    cat run.log
    cat micm-score.json
  )
done
