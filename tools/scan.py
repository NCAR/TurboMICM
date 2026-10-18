"""Reproducible native structure scans; every attempted run gets one CSV row.

Only elapsed_ms printed by C is a timing measurement. The supplied score.py is
invoked unchanged. Prefix screening uses matching prefix reference rows; final
operating-point validation uses all cells and all snapshots.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ROOT.parent / 'scenarios'
NAMES = ['tier1-30s', 'diurnal-30min', 'diurnal-10min', 'diurnal-5min',
         'polluted-10min', 'clean-10min']
COUNTERS = ('internal_steps accepted rejected function_calls jacobian_updates '
            'decompositions solves nonlinear_iterations matrix_reuses '
            'positivity_rejections failed_iterations').split()
FIELDS = ['study', 'configuration', 'method', 'options', 'scenario', 'cells',
          'repeat', 'compiler', 'threads', 'status', 'elapsed_ms',
          'overall_digits', 'worst_digits', 'worst_species', 'gate', 'values_counted'] + COUNTERS + [
          'input_sha256', 'reference_sha256', 'output_sha256', 'diagnostic']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(name, count):
    source = SCENARIOS / name
    if not count:
        return source
    out = ROOT / 'results' / f'prefix-{count}' / name
    out.mkdir(parents=True, exist_ok=True)
    for filename in ['cells.csv', 'reference_final.csv', 'reference_snapshots.csv']:
        path = source / filename
        if not path.exists():
            continue
        with path.open(newline='') as f, (out / filename).open('w', newline='') as g:
            reader = csv.reader(f); writer = csv.writer(g)
            header = next(reader); writer.writerow(header)
            for i, row in enumerate(reader):
                if (int(row[1]) if filename == 'reference_snapshots.csv' else i) < count:
                    writer.writerow(row)
    for filename in ['scenario.txt', 'photolysis.csv']:
        (out / filename).write_bytes((source / filename).read_bytes())
    return out


def run(args):
    configs = json.loads(args.configs.read_text())
    if args.select:
        configs = [c for c in configs if c['name'] in args.select.split(',')]
    if not configs:
        raise ValueError('no configurations selected')
    out = ROOT / 'scans' / f'{args.study}.csv'
    out.parent.mkdir(exist_ok=True)
    working = ROOT / 'results' / 'scan-run'; working.mkdir(parents=True, exist_ok=True)
    done = set()
    if out.exists():
        with out.open() as f:
            done = {(r['configuration'], r['scenario'], int(r['repeat'])) for r in csv.DictReader(f)}
    with out.open('a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if f.tell() == 0:
            writer.writeheader()
        for config in configs:
            options = config.get('options', {})
            for name in args.scenarios.split(','):
                scenario = prepare(name, args.cells)
                reference = scenario / ('reference_snapshots.csv' if (scenario / 'reference_snapshots.csv').exists() else 'reference_final.csv')
                for repeat in range(args.repeats):
                    if (config['name'], name, repeat) in done:
                        continue
                    for filename in ['ts1_final.csv', 'ts1_snapshots.csv']:
                        (working / filename).unlink(missing_ok=True)
                    command = [str(ROOT / 'bin/ts1_scenario'), '--threads', str(args.threads), '--method', config['method']]
                    for key, value in options.items():
                        command += ['--option', f'{key}={value}']
                    command.append(str(scenario / 'scenario.txt'))
                    row = dict(study=args.study, configuration=config['name'], method=config['method'],
                               options=json.dumps(options, sort_keys=True, separators=(',', ':')),
                               scenario=name, cells=args.cells or 1000, repeat=repeat,
                               compiler=args.compiler, threads=args.threads,
                               input_sha256=digest(scenario / 'cells.csv'), reference_sha256=digest(reference))
                    try:
                        result = subprocess.run(command, cwd=working, capture_output=True, text=True,
                                                timeout=args.timeout, env=dict(os.environ, OMP_SCHEDULE='static'))
                        row['status'] = 'OK' if result.returncode == 0 else 'FAIL'
                        row['diagnostic'] = result.stderr.strip()[:1500]
                        for token in result.stdout.split():
                            if '=' in token:
                                key, value = token.split('=', 1)
                                if key in FIELDS:
                                    row[key] = value
                        if result.returncode == 0:
                            output = working / ('ts1_snapshots.csv' if reference.name == 'reference_snapshots.csv' else 'ts1_final.csv')
                            score = subprocess.run([sys.executable, str(SCENARIOS / 'score.py'), str(output), str(reference), '--json'],
                                                   capture_output=True, text=True, check=False)
                            info = json.loads(score.stdout)
                            row.update({key: info[key] for key in ['overall_digits', 'gate', 'values_counted']})
                            row['worst_species'], row['worst_digits'] = info['worst_species'][0]
                            row['output_sha256'] = digest(output)
                            if score.returncode:
                                row['status'] = 'GATE_FAIL'
                    except subprocess.TimeoutExpired:
                        row.update(status='TIMEOUT', diagnostic=f'native process exceeded {args.timeout}s')
                    writer.writerow(row); f.flush()
                    print(config['name'], name, row['status'], row.get('elapsed_ms', ''),
                          row.get('overall_digits', ''), row.get('worst_species', ''), row.get('worst_digits', ''), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('configs', type=Path)
    p.add_argument('--study', required=True)
    p.add_argument('--select', default='')
    p.add_argument('--cells', type=int, default=0)
    p.add_argument('--threads', type=int, default=10)
    p.add_argument('--repeats', type=int, default=1)
    p.add_argument('--compiler', default='clang')
    p.add_argument('--timeout', type=float, default=300)
    p.add_argument('--scenarios', default=','.join(NAMES))
    run(p.parse_args())
