#!/usr/bin/env python3
"""Check the Tier 1 case-ID contract against the read-only manifest."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
expected=[line.split() for line in (ROOT.parent/'reference/micm_cases.txt').read_text().splitlines() if line.strip() and not line.startswith('#')]
actual=[line.split() for line in (ROOT/'tests/cases.txt').read_text().splitlines() if line.strip() and not line.startswith('#')]
assert actual==expected, 'Tier 1 case manifest changed'
print(f'Exact reference case manifest PASS: {len(actual)} IDs')
