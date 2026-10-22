"""Run Tier 2 numerical gates and diagnostic timings using only local inputs.

Build the requested compiler with make bench test-cpu test-specialized micm-ref
first. This script does not change compiler flags. Full CSV evidence remains in
validation/; a compact summary is written there for inclusion in RESULTS.md.
Python orchestrates tests only; it never executes inside either timed C driver.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run(command, directory, log):
    directory.mkdir(parents=True, exist_ok=True)
    result = subprocess.run([str(x) for x in command], cwd=directory, capture_output=True, text=True)
    (directory / log).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f'{log} failed: {result.stdout}\n{result.stderr}')
    return result.stdout


def stats(stdout):
    return {key: float(value) if key == 'elapsed_ms' else int(value)
            for key, value in re.findall(r'(\w+)=([0-9.]+)', stdout)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', required=True)
    parser.add_argument('--skip-one-thread', action='store_true')
    args = parser.parse_args()
    base = ROOT / 'validation' / args.compiler
    base.mkdir(parents=True, exist_ok=True)
    summary = {'compiler': args.compiler, 'benchmarks': [], 'scenarios': []}
    source = base / 'cells-10000-2026.csv'
    run([sys.executable, ROOT.parent / 'reference/make_cells_csv.py', '10000', '2026', '--output', source], base, 'cells.log')
    reference = None
    for threads in ([18] if args.skip_one_thread else [18, 1]):
        for executable in ('micm_ts1_ref', 'ts1_bench_cpu', 'ts1_bench_specialized'):
            label = f'{executable}-{threads}'
            directory = base / label
            print(label, flush=True)
            output = run([ROOT / 'bin' / executable, '--threads', threads, 10000, 100, 30, source], directory, 'run.log')
            print(output, flush=True)
            if reference is None:
                reference = directory / 'ts1_final.csv'
            comparison = run([sys.executable, ROOT.parent / 'reference/compare_final.py',
                              directory / 'ts1_final.csv', reference], directory, 'comparison.log')
            summary['benchmarks'].append(dict(executable=executable, threads=threads, **stats(output), comparison=comparison.strip()))
            (base / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    for scenario in sorted((ROOT.parent / 'scenarios').glob('*/scenario.txt')):
        name = scenario.parent.name
        directory = base / name
        print(name, flush=True)
        output = run([ROOT / 'bin/ts1_scenario_specialized', '--threads', 18,
                      '--group-size', 32, '--rtol', '1e-6', '--atol', '1e-15', scenario], directory, 'run.log')
        record = dict(name=name, **stats(output))
        for kind in ('micm', 'reference'):
            for data in ('snapshots', 'final'):
                command = [sys.executable, ROOT.parent / 'scenarios/score.py', directory / f'ts1_{data}.csv',
                           scenario.parent / f'{kind}_{data}.csv', '--json']
                if kind == 'micm':
                    command += ['--level', 10]
                record[f'{kind}_{data}'] = json.loads(run(command, directory, f'{kind}_{data}.log'))
        summary['scenarios'].append(record)
        print(json.dumps(record), flush=True)
        (base / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    main()
