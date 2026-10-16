#!/bin/zsh
# Write the disclosed Tier 3 scenarios: cells.csv, photolysis.csv, scenario.txt.
#
#   admin/make_scenarios.sh DIR [CELLS]
#
# With DIR = scenarios and the default 1,000 cells, this remakes the disclosed set
# exactly; admin/make_references.sh then adds the reference solutions. With
# CELLS = 10000 it makes the 10,000-cell versions used for timing, whose first
# 1,000 cells are the disclosed ones.
set -eu
make=${0:A:h}/../scenarios/make_scenario.py
dir=$1
cells=${2:-1000}
mkdir -p $dir
cd $dir
python3 $make tier1-30s --cells $cells --seed 2026 --dt 30 --steps 101 --constant
for spec in 5:300:288:72 10:600:144:36 30:1800:48:12; do
  IFS=: read minutes dt steps every <<< $spec
  python3 $make diurnal-${minutes}min --cells $cells --seed 101 --dt $dt --steps $steps \
      --diurnal --lat 40 --doy 80 --start 6 --snapshot-every $every
done
python3 $make polluted-10min --cells $cells --seed 102 --dt 600 --steps 144 \
    --diurnal --lat 40 --doy 80 --start 6 --nox-scale 10 --snapshot-every 36
python3 $make clean-10min --cells $cells --seed 103 --dt 600 --steps 144 \
    --diurnal --lat 40 --doy 80 --start 6 --nox-scale 0.1 --snapshot-every 36
