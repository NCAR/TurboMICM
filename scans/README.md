# Scans

One CSV per scan, written by `tools/scan.py` (C runs: time, counts, and
accuracy on the full scenarios) or by `prototypes/scenario.py --csv-row`
(Python runs: counts and accuracy on small cases). SCANS.md names each scan,
what it varied, and what it decided. Rows are never deleted; a rerun appends
with a later `started` time.
