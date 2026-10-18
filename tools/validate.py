"""Run both original gates and Tier 3 checks, retaining portable transcripts."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', choices=['clang', 'gcc-15'], required=True)
    args = parser.parse_args()
    label = 'gcc' if args.compiler == 'gcc-15' else 'clang'
    dest = ROOT / 'results' / ('validation-' + label)
    dest.mkdir(parents=True, exist_ok=True)
    assert (ROOT / 'tests/cases.txt').read_bytes() == (ROOT.parent / 'reference/micm_cases.txt').read_bytes()
    jobs = [
        ('checks.txt', ['make', 'cpu', 'bench', 'test-cpu', 'CC=' + args.compiler, '-j4']),
        ('tier3.txt', ['make', 'test-tier3', 'CC=' + args.compiler, '-j4']),
        ('ts1_gate.txt', [sys.executable, 'tools/check_ts1.py', '--output', 'build/final-gate-' + label]),
    ]
    for filename, command in jobs:
        print('Running', ' '.join(command).replace(str(ROOT), '.'), flush=True)
        run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
        transcript = (run.stdout + run.stderr).replace(str(ROOT), '.')
        (dest / filename).write_text(transcript)
        if run.returncode:
            print(transcript[-8000:])
            raise SystemExit(run.returncode)
        print(filename, 'PASS', flush=True)
    for filename in ['prototype_agreement.csv', 'convergence.csv', 'scenario_io.csv']:
        shutil.copyfile(ROOT / 'results' / filename, dest / filename)
    hashes = {name: hashlib.sha256((ROOT / 'bin' / name).read_bytes()).hexdigest()
              for name in ['ts1_scenario', 'ts1_bench_cpu', 'adaptive_probe', 'test_cpu', 'test_rate_updates']}
    (dest / 'build.json').write_text(json.dumps({
        'compiler': (ROOT / 'build/compiler.conf').read_text(),
        'compiler_version': subprocess.check_output([args.compiler, '--version'], text=True).splitlines()[0],
        'binaries': hashes, 'original_case_count': 680, 'original_case_list_identical': True,
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
