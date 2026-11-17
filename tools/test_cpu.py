#!/usr/bin/env python3
"""Report each reference case independently; unported cases are failures."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
run = subprocess.run([str(root / 'build/test_cases')], text=True, capture_output=True, check=True)
results = {}
for line in run.stdout.splitlines():
    case, passed = line.split()
    if case in results:
        raise ValueError(f'duplicate case result: {case}')
    results[case] = passed == '1'
failed = False
for line in (root / 'tests/micm_cases.txt').read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    case, kind = line.split()
    if case in results:
        ok = results.pop(case)
        status, reason = ('PASS', 'ported assertions') if ok else ('FAIL', 'ported assertion failed')
        failed |= not ok
    else:
        status, reason = 'FAIL', 'case not yet ported'
        failed = True
    print(f'{case} {status} {reason}')
if results:
    raise ValueError(f'unknown case results: {list(results)}')
sys.exit(int(failed))
