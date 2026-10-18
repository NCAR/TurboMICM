"""Generate extra 48-hour inputs locally and independently converged MICM references.

Run after `make micm-scenario`; all output remains inside this repository. These
inputs are additional validation, not undisclosed experimenter scenarios.
"""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'scenarios'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    configurations = json.loads((ROOT / 'scans/stress-scenarios.json').read_text())
    destination = ROOT / 'results/stress'
    scratch = ROOT / 'results/stress-reference-run'; scratch.mkdir(parents=True, exist_ok=True)
    log = ROOT / 'scans/stress-reference-generation.csv'
    fields = ['scenario', 'rtol', 'atol', 'status', 'elapsed_ms', 'overall_digits', 'worst_species', 'worst_digits',
              'gate', 'internal_steps', 'accepted', 'rejected', 'function_calls', 'jacobian_updates', 'decompositions', 'solves',
              'input_sha256', 'output_sha256', 'command']
    completed = set()
    if log.exists():
        completed = {(r['scenario'], r['rtol']) for r in csv.DictReader(log.open()) if r['status'] == 'OK'}
    with log.open('a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator='\n')
        if f.tell() == 0:
            writer.writeheader()
        for config in configurations:
            name = config['name']; path = destination / name
            command = [sys.executable, str(SOURCE / 'make_scenario.py'), str(path), '--diurnal']
            for key, value in config.items():
                if key != 'name':
                    command += ['--'+key.replace('_', '-'), str(value)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            for rtol, atol in [(1e-10, 1e-23), (1e-9, 1e-22)]:
                if (name, str(rtol)) in completed and (path / 'reference_snapshots.csv').exists():
                    continue
                command = [str(ROOT / 'bin/micm_scenario'), '--threads', '8', '--group-size', '8', '--method', 'rosenbrock6dae',
                           '--rtol', str(rtol), '--atol', str(atol), '--max-steps', '100000', str(path / 'scenario.txt')]
                result = subprocess.run(command, cwd=scratch, capture_output=True, text=True, timeout=1800)
                row = dict(scenario=name, rtol=rtol, atol=atol, status='OK' if result.returncode == 0 else 'FAIL',
                           input_sha256=sha(path / 'cells.csv'), command=json.dumps(command))
                for token in result.stdout.split():
                    if '=' in token:
                        key, value = token.split('=', 1)
                        if key in fields:
                            row[key] = value
                if result.returncode:
                    writer.writerow(row); f.flush(); raise RuntimeError(result.stderr)
                if rtol == 1e-10:
                    for filename in ['final', 'snapshots']:
                        (path / f'reference_{filename}.csv').write_bytes((scratch / f'ts1_{filename}.csv').read_bytes())
                score = subprocess.run([sys.executable, str(SOURCE / 'score.py'), str(scratch / 'ts1_snapshots.csv'),
                                        str(path / 'reference_snapshots.csv'), '--json'], capture_output=True, text=True, check=True)
                info = json.loads(score.stdout)
                row.update(overall_digits=info['overall_digits'], worst_species=info['worst_species'][0][0],
                           worst_digits=info['worst_species'][0][1], gate=info['gate'], output_sha256=sha(scratch / 'ts1_snapshots.csv'))
                writer.writerow(row); f.flush(); print(name, rtol, info['overall_digits'], info['worst_species'][0], flush=True)
                if rtol == 1e-9:
                    assert info['overall_digits'] >= 8 and info['worst_species'][0][1] >= 7, info


if __name__ == '__main__':
    main()
