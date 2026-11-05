#!/bin/zsh
# Time the paper's retaining solver, the library's restarting loop, and MICM's
# baselines on the 10,000-cell versions of the eleven scenarios, in rotating
# round-robin order, scoring the first 64 cells against the archived references.
#   time_paper.sh OUT [ROUNDS] [THREADS]
set -eu
TM=$HOME/NSF-NCAR/TurboMICM
WORK=$HOME/NSF-NCAR/tier3-paper-work
ARCHIVE=$HOME/DavidFillmore/Fast-ODE-Atmospheric-Chem/data/scenarios
out=$1 rounds=${2:-3} threads=${3:-10}
programs=()
while IFS= read -r line; do
  [[ -z $line || $line == \#* ]] && continue
  programs+=("${line//\$HOME/$HOME}")
done < $WORK/programs.txt
for pair in tier1-30s:disclosed diurnal-30min:disclosed diurnal-10min:disclosed diurnal-5min:disclosed polluted-10min:disclosed clean-10min:disclosed summer-15min:heldout tropics-20min:heldout winter60-10min:heldout polluted-48h-30min:heldout clean-48h-5min:heldout; do
  name=${pair%%:*} set=${pair##*:}
  if [[ $set == disclosed ]]; then d=$HOME/NSF-NCAR/tier2-work/check/disclosed-10000/$name; else d=$WORK/heldout-10000/$name; fi
  ref=$ARCHIVE/$name/reference_snapshots.csv
  [[ -f $ref ]] || ref=$ARCHIVE/$name/reference_final.csv
  echo "== $set/$name: $(date '+%H:%M:%S')"
  python3 $TM/admin/time_entries.py --scenario $d/scenario.txt --reference $ref \
      --out $out/$set/$name --threads $threads --rounds $rounds --pause 5 "${programs[@]}"
done
echo "== done: $(date '+%H:%M:%S')"
