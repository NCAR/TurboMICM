#!/bin/zsh
# Benchmark MICM's 24 CPU solver configurations on TS1 against its default.
#
#   [THREADS=10] [ROUNDS=3] [ONLY=pattern] [LABEL=name] \
#       admin/micm_configs.sh ROOT OUT CELLS_CSV MICM_FINAL [NAME=COMMAND ...]
#
# micm_bench offers 24 CPU configurations: the matrix layout (standard ordering,
# or vector ordering with 1, 2, 4, 8, or 128 cells interleaved), the LU storage
# (in-place, or separate L and U factors), and the LU algorithm (Mozart or
# Doolittle). The default is standard ordering with in-place Mozart LU. This
# script:
# 1. builds ROOT/reference/micm_ts1_ref.cpp in each configuration into OUT/bin,
#    from a copy whose solver types are chosen by the TS1_VECTOR,
#    TS1_LU_SEPARATE, and TS1_LU_DOOLITTLE macros (admin/micm_configs_source.py,
#    which leaves the reference itself alone);
# 2. runs each on the 100-cell example in ROOT/reference, on one thread, against
#    MICM's result there (OUT/check.csv). The default must reproduce it exactly;
# 3. times each, with one cell group per thread and with 32-cell groups, using
#    admin/time_entries.py: 48 programs on CELLS_CSV (10,000 cells, 100 steps of
#    30 s), at THREADS threads (default 10), for ROUNDS rounds (default 3), each
#    run checked against MICM_FINAL (OUT/timing). Any NAME=COMMAND programs, such
#    as the entries' bin/ts1_bench_cpu, are timed in the same rounds;
#    ONLY, a zsh pattern such as '*-in-place-mozart', times only the matching
#    configurations, plus the default as the baseline; each session's timings go
#    to OUT/timing-LABEL (LABEL defaults to all). Binaries already built from an
#    unchanged copy are reused;
# 4. writes OUT/configs-LABEL.md: each configuration's median time, relative to the
#    default with the same grouping, and its largest difference from MICM_FINAL,
#    then the extra programs, relative to both defaults.
# Confirm the fastest configurations with admin/time_entries.py, at 10 and 1
# threads, alongside the entries.
set -eu
here=${0:A:h}
root=${1:A}
out=${2:A}
cells_csv=${3:A}
micm_final=${4:A}
shift 4
threads=${THREADS:-10}
rounds=${ROUNDS:-3}
only=${ONLY:-*}
label=${LABEL:-all}
mkdir -p $out/bin $out/logs
export SDKROOT=$(xcrun --sdk macosx --show-sdk-path)

configs=()
for matrix in standard vector1 vector2 vector4 vector8 vector128; do
  for storage in in-place separate; do
    for lu in mozart doolittle; do
      configs+=($matrix-$storage-$lu)
    done
  done
done
default=standard-in-place-mozart
source=$out/micm_ts1_ref_configs.cpp
python3 $here/micm_configs_source.py $root/reference/micm_ts1_ref.cpp $source.new
if cmp -s $source.new $source; then rm $source.new; else mv $source.new $source; fi

build() {
  local config=$1
  local -a flags
  flags=()
  [[ $config == vector* ]] && flags+=(-DTS1_VECTOR=${${config%%-*}#vector})
  [[ $config == *-separate-* ]] && flags+=(-DTS1_LU_SEPARATE)
  [[ $config == *-doolittle ]] && flags+=(-DTS1_LU_DOOLITTLE)
  /usr/bin/clang++ -std=c++20 -O3 -DNDEBUG $flags -I $root/micm/include -I $root/micm/test/benchmark \
      $source -o $out/bin/micm-$config
}

echo "building ${#configs} configurations"
i=0
for config in $configs; do
  [[ -x $out/bin/micm-$config && $out/bin/micm-$config -nt $source ]] && continue
  build $config > $out/logs/build-$config.log 2>&1 &
  (( ++i % 6 == 0 )) && wait
done
wait
for config in $configs; do
  [[ -x $out/bin/micm-$config ]] || { echo "build failed: $config"; cat $out/logs/build-$config.log; exit 1; }
done

echo "checking every configuration on the 100-cell example"
example=$root/reference/ts1_cells_100_seed1.csv
expected=$root/reference/micm_ts1_final_100_seed1.csv
echo "config,internal_steps,identical_to_default,gate" > $out/check.csv
for config in $configs; do
  work=$(mktemp -d)
  if (cd $work && $out/bin/micm-$config --threads 1 100 100 30 $example > run.log 2>&1); then
    steps_taken=$(sed -n 's/^internal_steps=\([0-9]*\).*/\1/p' $work/run.log)
    identical=no
    cmp -s $work/ts1_final.csv $expected && identical=yes
    gate=$(python3 $root/reference/compare_final.py $work/ts1_final.csv $expected || true)
  else
    steps_taken= identical=no gate="run failed: $(tail -1 $work/run.log)"
  fi
  echo "$config,$steps_taken,$identical,\"$gate\"" | tee -a $out/check.csv
  rm -rf $work
done
[[ $(command grep -c "^$default,[0-9]*,yes," $out/check.csv) == 1 ]] || { echo "the default configuration does not reproduce MICM's result"; exit 1; }

programs=()
for config in $configs; do
  [[ $config == $default || $config == ${~only} ]] || continue
  programs+=("micm-$config=$out/bin/micm-$config" "micm-$config-g32=$out/bin/micm-$config --group-size 32")
done
echo "timing $(( ${#programs} / 2 )) configurations x 2 groupings, $threads threads, $rounds rounds -> $out/timing-$label"
python3 $here/time_entries.py --input $cells_csv --reference $micm_final --out $out/timing-$label \
    --threads $threads --rounds $rounds $programs "$@"

python3 $here/micm_configs_table.py $out $default $label
