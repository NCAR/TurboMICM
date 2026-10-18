# Derecho environment for the TurboMICM Fable entry: source this file.
#   source derecho/env.sh
# The module stack MUSICA builds with on Derecho (MUSICA/machines/derecho/build.sh):
# ncarenv/25.10 and gcc/14.3.0. Lmod's functions are not set -e safe, so load
# before any `set -e`.
if ! command -v module >/dev/null 2>&1; then
  source /etc/profile.d/z00_modules.sh 2>/dev/null || source /usr/share/lmod/lmod/init/bash
fi
module --force purge >/dev/null 2>&1
module load ncarenv/25.10 >/dev/null 2>&1
module load gcc/14.3.0 >/dev/null 2>&1
# ncarenv/25.10 leaves LIBRARY_PATH pointing at gcc 12.5.0's lib64, so g++ 14
# would link MICM's C++ reference against the older libstdc++ and fail with
# undefined _M_replace_cold and __cxa_call_terminate. Without it, g++ 14 links
# and runs with its own libstdc++.
unset LIBRARY_PATH
# the conda module is not needed: the fable Makefile finds Derecho's miniforge
export TURBOMICM_PBS_ACCOUNT=${TURBOMICM_PBS_ACCOUNT:-P19010000}
export TURBOMICM_DERECHO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# one OpenMP thread per core, threads kept where they start
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_PLACES=${OMP_PLACES:-cores}
