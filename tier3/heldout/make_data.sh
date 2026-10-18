#!/bin/zsh
# Make Tier 3's held-out scenario files and references, and the
# 10,000-cell timing versions of the disclosed and held-out scenarios, outside the
# run root, in ~/NSF-NCAR/.tier3-heldout/:
#   bin/micm_scenario   MICM's scenario driver (tier3/scenarios/micm_scenario.cpp)
#   heldout/            the five held-out scenarios, 1,000 cells, with references
#   heldout-10000/      the same at 10,000 cells, for timing
#   disclosed-10000/    the six disclosed scenarios at 10,000 cells, for timing
# Each part is skipped if it exists. On 2026-10-18 it took 27 minutes, 26 of them
# for the references; compare the result with SHA256SUMS in this folder.
set -eu
here=${0:A:h}
T3=$HOME/NSF-NCAR/ODE/benchmarks/tier3
H=$HOME/NSF-NCAR/.tier3-heldout
export SDKROOT=$(xcrun --sdk macosx --show-sdk-path)
mkdir -p $H/bin
if [[ ! -x $H/bin/micm_scenario ]]; then
  echo "build: $(date '+%H:%M:%S')"
  clang++ -std=c++20 -O3 -DNDEBUG -I $T3/micm/include -I $T3/micm/test/benchmark \
      $T3/scenarios/micm_scenario.cpp -o $H/bin/micm_scenario
fi
echo "files: $(date '+%H:%M:%S')"
[[ -d $H/heldout ]] || $here/make_heldout.sh $T3/scenarios/make_scenario.py $H/heldout 1000
[[ -d $H/heldout-10000 ]] || $here/make_heldout.sh $T3/scenarios/make_scenario.py $H/heldout-10000 10000
[[ -d $H/disclosed-10000 ]] || $T3/admin/make_scenarios.sh $H/disclosed-10000 10000
for d in $H/heldout/*(/); do
  [[ -f $d/reference_final.csv ]] && continue
  echo "reference ${d:t}: $(date '+%H:%M:%S')"
  $T3/admin/make_references.sh $H/bin/micm_scenario 10 $d
done
echo "done: $(date '+%H:%M:%S')"
