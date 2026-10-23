# TS1 exercises: sol-high (tier2/sol/high)

TS1 configuration: `configs/ts1.json`

Baseline: make env exit 0; make cpu bench exit 0 (19 s); specialized against generic: {'pass': True, 'largest': '3.9e-14'}

| Exercise | Edit | make bench | Rebuilt | Specialized vs generic (1e-10) | Generic changed | Specialized changed | Verdict |
|---|---|---|---|---|---|---|---|
| add | appended a copy of reaction 342 (NO + O3 -> NO2 + O2), a second channel | exit 0 (16.1 s) | yes | PASS 4.39e-14 | yes, 2.83e+11 | yes, 2.83e+11 | PASS |
| change | multiplied A of reaction 342 (NO + O3 -> NO2 + O2) by 1.5 | exit 0 (15.8 s) | yes | PASS 4.71e-14 | yes, 1.62e+07 | yes, 1.62e+07 | PASS |
