"""Run all disclosed scenarios and score snapshots and finals with the supplied tool.

This is test orchestration only; the native executable handles every timed step.
Outputs and machine-readable evidence stay under results/ in this repository.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAMES = ('diurnal-30min', 'diurnal-10min', 'diurnal-5min', 'polluted-10min', 'clean-10min')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', required=True)
    parser.add_argument('--threads', type=int, default=18)
    parser.add_argument('--program',type=Path,default=ROOT/'bin/ts1_scenario_specialized')
    args = parser.parse_args()
    destination = (ROOT/'results'/args.label).resolve()
    if not destination.is_relative_to(ROOT/'results'):
        parser.error('label must resolve inside results/')
    destination.mkdir(parents=True,exist_ok=True)
    summary = {}
    for name in NAMES:
        source = ROOT.parent/'scenarios'/name
        output = destination/name
        output.mkdir(exist_ok=True)
        command = [str(args.program.resolve()),'--threads',str(args.threads),
                   '--group-size','32','--rtol','1e-6','--atol','1e-15',str(source/'scenario.txt')]
        run = subprocess.run(command,cwd=output,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        (output/'run.log').write_text(run.stdout)
        run.check_returncode()
        scores = {}
        for artifact in ('snapshots','final'):
            for reference in ('micm','reference'):
                command = [sys.executable,str(ROOT.parent/'scenarios/score.py'),str(output/f'ts1_{artifact}.csv'),
                           str(source/f'{reference}_{artifact}.csv'),'--json']
                if reference == 'micm':
                    command += ['--level','10']
                score = subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
                (output/f'{reference}_{artifact}_score.json').write_text(score.stdout)
                score.check_returncode()
                scores[f'{reference}_{artifact}'] = json.loads(score.stdout)
        summary[name] = dict(command=run.args,stdout=run.stdout,scores=scores)
        print(name,run.stdout.splitlines()[0],
              'MICM digits',scores['micm_snapshots']['overall_digits'],
              'worst',scores['micm_snapshots']['worst_species'][0],flush=True)
        (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__ == '__main__':
    main()
