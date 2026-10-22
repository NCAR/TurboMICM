"""Run the disclosed day-night gates and retain both MICM and reference scores."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, default=Path('bin/ts1_scenario_specialized'))
p.add_argument('--scenarios', type=Path, default=Path('../scenarios'))
p.add_argument('--output', type=Path, required=True)
p.add_argument('--threads', type=int, default=18)
p.add_argument('--names', nargs='*')
a = p.parse_args()
root = Path.cwd().resolve()
output = a.output.resolve()
if not output.is_relative_to(root):
    p.error('Output must remain in the checkout')
source = a.scenarios.resolve()
names = a.names or sorted(path.parent.name for path in source.glob('*/scenario.txt'))
summary = {}
for name in names:
    work = output/name
    work.mkdir(parents=True, exist_ok=True)
    command = [str(a.binary.resolve()), '--threads', str(a.threads), '--group-size', '32',
               '--rtol', '1e-6', '--atol', '1e-15', str(source/name/'scenario.txt')]
    run = subprocess.run(command, cwd=work, capture_output=True, text=True)
    (work/'run.log').write_text(run.stdout+run.stderr)
    run.check_returncode()
    result = {'command': command, 'metrics': {k: float(v) if k == 'elapsed_ms' else int(v)
              for k, v in (word.split('=') for word in run.stdout.split() if '=' in word)}}
    for kind in ['snapshots', 'final']:
        for reference in ['micm', 'reference']:
            command = [sys.executable, str(source/'score.py'), str(work/f'ts1_{kind}.csv'),
                       str(source/name/f'{reference}_{kind}.csv'), '--json']
            if reference == 'micm':
                command += ['--level', '10']
            score = subprocess.run(command, capture_output=True, text=True)
            (work/f'{reference}_{kind}.json').write_text(score.stdout)
            score.check_returncode()
            result[f'{reference}_{kind}'] = json.loads(score.stdout)
    summary[name] = result
    print(name, result['metrics'], result['micm_snapshots'], flush=True)
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
