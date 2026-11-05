#!/bin/zsh
# Time MICM's backward Euler in the library: at its default tolerances on all
# eleven 10,000-cell scenarios, and with the equal sub-steps that first meet
# criterion 2 on the three scenarios of the paper's comparison, beside the
# retaining solver's criterion-2 point as the in-session control.
#   time_be.sh OUT [ROUNDS] [THREADS]
set -eu
TM=$HOME/NSF-NCAR/TurboMICM
WORK=$HOME/NSF-NCAR/tier3-paper-work
ARCHIVE=$HOME/DavidFillmore/Fast-ODE-Atmospheric-Chem/data/scenarios
BIN=$HOME/NSF-NCAR/TurboMICM-paper/bin/ts1_scenario
out=$1 rounds=${2:-3} threads=${3:-10}
be_default="be/default=$BIN --method be --option rtol=1e-6 --option atol=1e-3"
retaining="retaining/C2=$BIN --method ros3 --option rtol=0.01 --option atol=1e-16 --option quantize=1 --option controller=h211b --option safety=0.9 --option exact_end=1 --option reject_elementary=1 --option strict_min=1"
typeset -A subs
subs=(tier1-30s 4 diurnal-10min 64 diurnal-30min 256)
for pair in tier1-30s:disclosed diurnal-30min:disclosed diurnal-10min:disclosed diurnal-5min:disclosed polluted-10min:disclosed clean-10min:disclosed summer-15min:heldout tropics-20min:heldout winter60-10min:heldout polluted-48h-30min:heldout clean-48h-5min:heldout; do
  name=${pair%%:*} set=${pair##*:}
  if [[ $set == disclosed ]]; then d=$HOME/NSF-NCAR/tier2-work/check/disclosed-10000/$name; else d=$WORK/heldout-10000/$name; fi
  ref=$ARCHIVE/$name/reference_snapshots.csv
  [[ -f $ref ]] || ref=$ARCHIVE/$name/reference_final.csv
  programs=("$be_default")
  if (( ${+subs[$name]} )); then
    n=$subs[$name]
    programs+=("be/x$n=$BIN --method be --option rtol=1e-6 --option atol=1e-3 --option substeps=$n" "$retaining")
  fi
  echo "== $set/$name: $(date '+%H:%M:%S')"
  python3 $TM/admin/time_entries.py --scenario $d/scenario.txt --reference $ref \
      --out $out/$set/$name --threads $threads --rounds $rounds --pause 5 "${programs[@]}"
done
echo "== done: $(date '+%H:%M:%S')"
