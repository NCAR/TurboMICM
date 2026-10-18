#!/bin/bash
# Build some of MICM's TS1 reference configurations with g++.
#
#   [OPTFLAGS="..."] [JOBS=4] derecho/tools/micm_build.sh OUT CONFIG...
#
# CONFIG is <matrix>-<storage>-<lu> as in tools/micm_configs.sh (standard,
# vector1, vector2, vector4, vector8, vector128; in-place, separate; mozart,
# doolittle). The configurable source OUT/micm_ts1_ref_configs.cpp is written
# if missing (admin/micm_configs_source.py); binaries go to OUT/bin/micm-CONFIG,
# JOBS at a time (each compile needs about 1 GB). Existing binaries newer than
# the source are kept.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
repo=$(cd "$here/.." && pwd)
out=$(realpath -m "$1"); shift
jobs=${JOBS:-4}
optflags=${OPTFLAGS:--O3 -DNDEBUG -march=znver3 -ffp-contract=fast}
mkdir -p "$out/bin" "$out/logs"
source=$out/micm_ts1_ref_configs.cpp
[[ -f $source ]] || python3 "$repo/admin/micm_configs_source.py" "$repo/reference/micm_ts1_ref.cpp" "$source"
[[ -f $out/flags.txt ]] || printf '%s\n' "$optflags" > "$out/flags.txt"
[[ "$(cat "$out/flags.txt")" == "$optflags" ]] || { echo "OUT was built with other flags: $(cat "$out/flags.txt")"; exit 1; }

# No -fopenmp: the reference threads with std::thread, and a linked libgomp
# pins the main thread to one core under OMP_PROC_BIND (env.sh), which every
# std::thread inherits: 128 threads on one core.
build() {
  local config=$1
  local flags=()
  if [[ $config == vector* ]]; then
    local v=${config%%-*}
    flags+=("-DTS1_VECTOR=${v#vector}")
  fi
  [[ $config == *-separate-* ]] && flags+=(-DTS1_LU_SEPARATE)
  [[ $config == *-doolittle ]] && flags+=(-DTS1_LU_DOOLITTLE)
  g++ -std=c++20 $optflags "${flags[@]}" -I "$here/micm/include" -I "$here/micm/test/benchmark" \
      "$source" -o "$out/bin/micm-$config" > "$out/logs/build-$config.log" 2>&1 \
    && echo "built micm-$config" || { echo "FAILED micm-$config"; tail -3 "$out/logs/build-$config.log"; }
}

running=0
for config in "$@"; do
  if [[ -x $out/bin/micm-$config && $out/bin/micm-$config -nt $source ]]; then
    echo "kept micm-$config"
    continue
  fi
  build "$config" &
  running=$((running + 1))
  if (( running >= jobs )); then
    wait -n
    running=$((running - 1))
  fi
done
wait
