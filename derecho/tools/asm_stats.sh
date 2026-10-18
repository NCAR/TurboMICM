#!/bin/bash
# Instruction classes of generated-kernel object files (x86-64), one line each.
#
#   derecho/tools/asm_stats.sh OBJECT...
#
# Counts, from objdump: all instructions; 256-bit (ymm) and 128-bit (xmm)
# vector FMAs and multiply/add/subtract; scalar double arithmetic; references
# to the stack (spills and stack temporaries); shuffles and inserts; and the
# object's text size.
set -eu
printf '%-34s %8s %8s %8s %8s %8s %7s %7s %8s %8s\n' object instrs ymm-fma xmm-fma ymm-mas xmm-mas scalar stack shuffle text
for o in "$@"; do
  asm=$(mktemp)
  objdump -d --no-show-raw-insn "$o" > "$asm"
  printf '%-34s %8d %8d %8d %8d %8d %7d %7d %8d %8d\n' "$(basename "$(dirname "$(dirname "$o")")")/$(basename "$o")" \
    "$(grep -cE '^\s+[0-9a-f]+:' "$asm")" \
    "$(grep -cE 'vfn?m(add|sub)[0-9]*pd.*ymm' "$asm")" \
    "$(grep -cE 'vfn?m(add|sub)[0-9]*pd.*xmm' "$asm")" \
    "$(grep -cE 'v(mul|add|sub)pd.*ymm' "$asm")" \
    "$(grep -cE '(v|\s)(mul|add|sub)pd.*xmm' "$asm")" \
    "$(grep -cE '(mul|add|sub|fn?madd|fn?msub)[0-9]*sd ' "$asm")" \
    "$(grep -cE '\(%r[sb]p' "$asm")" \
    "$(grep -cE 'v(insert|extract|perm|unpck|shuf)' "$asm")" \
    "$(size "$o" | awk 'NR==2 {print $1}')"
  rm -f "$asm"
done
