"""Generate 10,000-cell disclosed extensions; retain the supplied 1,000-cell reference.

score.py checks the gate on every output value, but reference-relative errors
only on the first 1,000 cells. This scope is explicit in the results report.
"""
import csv
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'scenarios'


def main():
    for config in json.loads((ROOT / 'scans/large-scenarios.json').read_text()):
        name = config['name']; out = ROOT / 'results/large' / name
        command = [sys.executable, str(SOURCE / 'make_scenario.py'), str(out)]
        for key, value in config.items():
            if key == 'name':
                continue
            command += ['--'+key.replace('_', '-')]
            if not isinstance(value, bool):
                command += [str(value)]
        subprocess.run(command, capture_output=True, text=True, check=True)
        with (out / 'cells.csv').open() as a, (SOURCE / name / 'cells.csv').open() as b:
            actual, reference = csv.reader(a), csv.reader(b)
            assert next(actual) == next(reference)
            for row in reference:
                assert [float(v) for v in next(actual)] == [float(v) for v in row], name
        for filename in ['reference_final.csv', 'reference_snapshots.csv']:
            path = SOURCE / name / filename
            if path.exists():
                (out / filename).write_bytes(path.read_bytes())
        print(name, '10000 cells; first 1000 match supplied input exactly', flush=True)


if __name__ == '__main__':
    main()
