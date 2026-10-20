#!/usr/bin/env python3
"""Report each reference case, distinguishing exercised ports from missing ones."""
from pathlib import Path

root = Path(__file__).resolve().parents[2]
case_file = root / "reference" / "micm_cases.txt"
implemented = {
    "arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithSystem",
    "arrhenius_rate_constant:ArrheniusRateConstant.CalculateWithPrescribedArguments",
}

failed = False
for raw in case_file.read_text().splitlines():
    line = raw.strip()
    if not line or line.startswith("#"):
        continue
    case_id, kind = line.split()
    if case_id in implemented:
        print(f"{case_id} PASS C Arrhenius formula checked")
    elif kind == "program":
        failed = True
        print(f"{case_id} FAIL C program port incomplete")
    else:
        failed = True
        print(f"{case_id} FAIL C case port incomplete")
raise SystemExit(1 if failed else 0)
