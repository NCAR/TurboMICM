#!/bin/bash
# Build and time MICM's 24 CPU solver configurations on Derecho, with GCC.
#
#   [THREADS=128] [ROUNDS=3] [ONLY=regex] [LABEL=name] [OPTFLAGS="..."] \
#       derecho/tools/micm_configs.sh OUT CELLS_CSV MICM_FINAL [NAME=COMMAND ...]
#
# A bash and g++ port of admin/micm_configs.sh (zsh, Apple clang), run inside a
# PBS job (derecho/tools/submit.sh). It
# 1. writes the configurable copy of reference/micm_ts1_ref.cpp with
#    admin/micm_configs_source.py and builds the 24 configurations (matrix
#    layout: standard, vector 1, 2, 4, 8, 128; LU storage: in-place, separate;
#    LU algorithm: Mozart, Doolittle) into OUT/bin with
#    `g++ -std=c++20 $OPTFLAGS -fopenmp`; OPTFLAGS defaults to the entry's
#    Derecho flags, -O3 -DNDEBUG -march=znver3 -ffp-contract=fast;
# 2. runs each on the 100-cell example on one thread against MICM's result
#    from the experiment's machine (OUT/check.csv): the gate, and whether the
#    file is identical;
# 3. times the configurations matching ONLY (an extended regex; default all),
#    plus the default, each with one cell group per thread and with 32-cell
#    groups, with admin/time_entries.py at THREADS threads for ROUNDS rounds on
#    CELLS_CSV (10,000 cells, 100 steps of 30 s), every run checked against
#    MICM_FINAL (OUT/timing-LABEL); any NAME=COMMAND programs run in the same
#    rounds;
# 4. writes OUT/configs-LABEL.md with admin/micm_configs_table.py.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1")
cells_csv=$(realpath "$2")
micm_final=$(realpath "$3")
shift 3
threads=${THREADS:-128}
rounds=${ROUNDS:-3}
only=${ONLY:-.*}
label=${LABEL:-all}
optflags=${OPTFLAGS:--O3 -DNDEBUG -march=znver3 -ffp-contract=fast}
micm=$here/micm
reference=$repo/reference
admin=$repo/admin
mkdir -p "$out/bin" "$out/logs"

configs=()
for matrix in standard vector1 vector2 vector4 vector8 vector128; do
  for storage in in-place separate; do
    for lu in mozart doolittle; do
      configs+=("$matrix-$storage-$lu")
    done
  done
done
default=standard-in-place-mozart
source=$out/micm_ts1_ref_configs.cpp
python3 "$admin/micm_configs_source.py" "$reference/micm_ts1_ref.cpp" "$source.new"
if cmp -s "$source.new" "$source"; then rm "$source.new"; else mv "$source.new" "$source"; fi
if [[ ! -f $out/flags.txt || "$(cat "$out/flags.txt")" != "$optflags" ]]; then
  rm -f "$out"/bin/micm-*
  printf '%s\n' "$optflags" > "$out/flags.txt"
fi

build() {
  local config=$1
  local flags=()
  if [[ $config == vector* ]]; then
    local v=${config%%-*}
    flags+=("-DTS1_VECTOR=${v#vector}")
  fi
  [[ $config == *-separate-* ]] && flags+=(-DTS1_LU_SEPARATE)
  [[ $config == *-doolittle ]] && flags+=(-DTS1_LU_DOOLITTLE)
  g++ -std=c++20 $optflags -fopenmp "${flags[@]}" -I "$micm/include" -I "$micm/test/benchmark" \
      "$source" -o "$out/bin/micm-$config"
}

echo "compiler: $(g++ --version | head -1); flags: $optflags"
echo "building ${#configs[@]} configurations"
for config in "${configs[@]}"; do
  [[ -x $out/bin/micm-$config && $out/bin/micm-$config -nt $source ]] && continue
  build "$config" > "$out/logs/build-$config.log" 2>&1 &
done
wait
for config in "${configs[@]}"; do
  [[ -x $out/bin/micm-$config ]] || { echo "build failed: $config"; cat "$out/logs/build-$config.log"; exit 1; }
done

echo "checking every configuration on the 100-cell example"
example=$reference/ts1_cells_100_seed1.csv
expected=$reference/micm_ts1_final_100_seed1.csv
echo "config,internal_steps,identical_to_default,gate" > "$out/check.csv"
for config in "${configs[@]}"; do
  work=$(mktemp -d)
  if (cd "$work" && "$out/bin/micm-$config" --threads 1 100 100 30 "$example" > run.log 2>&1); then
    steps_taken=$(sed -n 's/^internal_steps=\([0-9]*\).*/\1/p' "$work/run.log")
    identical=no
    cmp -s "$work/ts1_final.csv" "$expected" && identical=yes
    gate=$(python3 "$reference/compare_final.py" "$work/ts1_final.csv" "$expected" || true)
  else
    steps_taken= identical=no gate="run failed: $(tail -1 "$work/run.log")"
  fi
  echo "$config,$steps_taken,$identical,\"$gate\"" | tee -a "$out/check.csv"
  rm -rf "$work"
done

programs=()
for config in "${configs[@]}"; do
  [[ $config == "$default" || $config =~ $only ]] || continue
  programs+=("micm-$config=$out/bin/micm-$config" "micm-$config-g32=$out/bin/micm-$config --group-size 32")
done
echo "timing $(( ${#programs[@]} / 2 )) configurations x 2 groupings, $threads threads, $rounds rounds -> $out/timing-$label"
python3 "$admin/time_entries.py" --input "$cells_csv" --reference "$micm_final" --out "$out/timing-$label" \
    --threads "$threads" --rounds "$rounds" "${programs[@]}" "$@"

python3 "$admin/micm_configs_table.py" "$out" "$default" "$label"
