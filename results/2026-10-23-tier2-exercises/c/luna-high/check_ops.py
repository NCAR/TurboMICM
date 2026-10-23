#!/usr/bin/env python3
"""Report what chem_generated_mechanism.specialized_ops points to in a binary.

    check_ops.py BINARY

Reads symbol addresses with `nm` and the pointer fixups with `dyld_info
-fixups` (no execution). specialized_ops is the 10th 8-byte field of
struct chem_mechanism (include/chemistry.h), offset 0x48. A NULL pointer has
no rebase entry; otherwise the entry's target is named, and for generated_ops
its four function pointers are resolved too.
"""
import re
import subprocess
import sys

binary = sys.argv[1]
symbols = {}
for line in subprocess.run(["nm", binary], capture_output=True, text=True, check=True).stdout.splitlines():
    parts = line.split()
    if len(parts) == 3:
        symbols[parts[2]] = int(parts[0], 16)
by_address = {address: name for name, address in symbols.items()}
rebase = {}
for line in subprocess.run(["dyld_info", "-fixups", binary], capture_output=True, text=True,
                           check=True).stdout.splitlines():
    match = re.search(r"0x([0-9A-Fa-f]+)\s+rebase\s+0x([0-9A-Fa-f]+)", line)
    if match:
        rebase[int(match.group(1), 16)] = int(match.group(2), 16)

mechanism = symbols["_chem_generated_mechanism"]
field = mechanism + 0x48
target = rebase.get(field)
generated = sorted(name for name in symbols if name.startswith("_generated_"))
print(f"generated symbols: {' '.join(generated) if generated else 'none'}")
if target is None:
    print("specialized_ops = NULL (no fixup at _chem_generated_mechanism+0x48): generic kernels only")
    sys.exit(0)
print(f"specialized_ops -> {by_address.get(target, hex(target))}")
for slot, label in enumerate(("forcing_with_constants", "negative_jacobian_with_pattern",
                              "factor", "solve_factored")):
    pointer = rebase.get(target + 8 * slot)
    print(f"  {label} -> {by_address.get(pointer, 'NULL' if pointer is None else hex(pointer))}")
