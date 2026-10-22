"""Run built benchmarks and scenario gates, retaining commands and JSON evidence.

Run after `make bench micm-ref`. This tool never builds or runs anything outside
this repository; reference scripts and input files are read-only.
"""
import argparse
import json
from pathlib import Path
import re
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
SCENARIOS = ['diurnal-30min', 'diurnal-10min', 'diurnal-5min', 'polluted-10min', 'clean-10min']


def run(args, directory):
    result = subprocess.run([str(x) for x in args], cwd=directory, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (directory / 'last-command.log').write_text(' '.join(str(x) for x in args) + '\n' + result.stdout)
    if result.returncode:
        raise RuntimeError(result.stdout)
    return result.stdout


def metrics(output):
    found = re.search(r'elapsed_ms=([0-9.]+)', output)
    assert found, output
    return {'elapsed_ms': float(found[1]),
            'statistics': dict((k, int(v)) for k, v in re.findall(r'(\w+)=(\d+)(?=\s|$)', output))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, default=ROOT / 'build/evaluation')
    p.add_argument('--scenarios', action='store_true')
    p.add_argument('--benchmarks', action='store_true')
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--label', default='local')
    a = p.parse_args()
    a.directory = a.directory.resolve()
    a.directory.mkdir(parents=True, exist_ok=True)
    evidence = {'label': a.label, 'scenarios': {}, 'benchmarks': {}}
    evidence_file = a.directory / 'results.json'

    def save():
        evidence_file.write_text(json.dumps(evidence, indent=2) + '\n')

    if a.scenarios:
        for name in SCENARIOS:
            directory = a.directory / name
            directory.mkdir(exist_ok=True)
            source = ROOT.parent / 'scenarios' / name
            output = run([ROOT / 'bin/ts1_scenario_specialized', '--threads', 18,
                          '--group-size', 32, '--rtol', '1e-6', '--atol', '1e-15', source / 'scenario.txt'], directory)
            (directory / 'run.log').write_text(output)
            record = metrics(output)
            for kind in ('micm', 'reference'):
                args = [sys.executable, ROOT.parent / 'scenarios/score.py',
                        directory / 'ts1_snapshots.csv', source / (kind + '_snapshots.csv'), '--json']
                if kind == 'micm':
                    args += ['--level', 10]
                record[kind] = json.loads(run(args, directory))
            evidence['scenarios'][name] = record
            print(name, record, flush=True)
            save()

    if a.benchmarks:
        assert a.repeats > 0
        cells = a.directory / 'cells_10000_seed2026.csv'
        run([sys.executable, ROOT.parent / 'reference/make_cells_csv.py', 10000, 2026, '--output', cells], a.directory)
        for threads in (18, 1):
            for program in ('ts1_bench_specialized', 'ts1_bench_cpu', 'micm_ts1_ref'):
                directory = a.directory / f'{program}-{threads}'
                directory.mkdir(exist_ok=True)
                records = []
                for trial in range(a.repeats):
                    output = run([ROOT / 'bin' / program, '--threads', threads, 10000, 100, 30, cells], directory)
                    (directory / f'run-{trial}.log').write_text(output)
                    records.append(metrics(output))
                    print(program, threads, trial, records[-1], flush=True)
                evidence['benchmarks'][f'{program}-{threads}'] = {
                    'runs': records, 'median_ms': statistics.median(r['elapsed_ms'] for r in records)}
                save()
            reference = a.directory / f'micm_ts1_ref-{threads}/ts1_final.csv'
            for program in ('ts1_bench_specialized', 'ts1_bench_cpu'):
                directory = a.directory / f'{program}-{threads}'
                output = run([sys.executable, ROOT.parent / 'reference/compare_final.py',
                              directory / 'ts1_final.csv', reference], directory)
                evidence['benchmarks'][f'{program}-{threads}']['accuracy'] = output.strip()
                save()
    save()


if __name__ == '__main__':
    main()
