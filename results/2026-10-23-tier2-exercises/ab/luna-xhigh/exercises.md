# TS1 exercises: luna-xhigh (tier2/luna/xhigh)

TS1 configuration: `config/ts1.json`

Baseline: make env exit 0; make cpu bench exit 0 (71 s); specialized against generic: {'pass': True, 'largest': '3.46e-14'}

| Exercise | Edit | make bench | Rebuilt | Specialized vs generic (1e-10) | Generic changed | Specialized changed | Verdict |
|---|---|---|---|---|---|---|---|
| add | appended a copy of reaction 342 (NO + O3 -> NO2 + O2), a second channel | exit 0 (67.2 s) | yes | PASS 3.07e-14 | yes, 2.83e+11 | yes, 2.83e+11 | PASS |
| change | multiplied A of reaction 342 (NO + O3 -> NO2 + O2) by 1.5 | exit 0 (62.3 s) | yes | PASS 3.38e-14 | yes, 1.62e+07 | yes, 1.62e+07 | PASS |
