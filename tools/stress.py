"""Prepare three deterministic 48-hour stress cases and converged MICM references.

Inputs/references stay in build/stress. The experiment recipe, hashes, C timing,
operation counts, and supplied-scorer convergence results are retained in git.
"""
import csv
import argparse
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from scan import FIELDS, ROOT, SCENARIOS

CASES = [
    dict(name='summer48', seed=777, lat=60, doy=172, start=0, dt=1800, steps=96, nox=5, snapshot=24),
    dict(name='winter48', seed=2027, lat=55, doy=355, start=12, dt=600, steps=288, nox=.03, snapshot=72),
    dict(name='tropical48', seed=858, lat=-15, doy=270, start=18, dt=900, steps=192, nox=20, snapshot=48),
]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scan',default='stress_references')
    args=parser.parse_args()
    destination = ROOT/'build/stress'
    destination.mkdir(parents=True, exist_ok=True)
    output = ROOT/'scans'/f'{args.scan}.csv'
    if output.exists():
        raise SystemExit(f'Refusing to overwrite {output}')
    manifest=[]
    with output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS); writer.writeheader(); f.flush()
        for case in CASES:
            directory=destination/case['name']; directory.mkdir(exist_ok=True)
            command=[sys.executable,str(SCENARIOS/'make_scenario.py'),str(directory),
                     '--cells','128','--seed',str(case['seed']),'--diurnal']
            for key in ['lat','doy','start','dt','steps']:
                command += ['--'+key,str(case[key])]
            command += ['--nox-scale',str(case['nox']),'--snapshot-every',str(case['snapshot'])]
            subprocess.run(command,check=True)
            rows=[]
            for label,rtol,atol in [('tight',1e-9,1e-22),('tighter',1e-10,1e-23)]:
                work=directory/label; work.mkdir(exist_ok=True)
                binary=ROOT/'bin/micm_scenario'
                options={'rtol':rtol,'atol':atol,'max-steps':100000,'group-size':32}
                command=[str(binary),'--threads','4','--method','rosenbrock6dae']
                for key,value in options.items(): command += ['--'+key,str(value)]
                command += [str(directory/'scenario.txt')]
                print(case['name'],label,flush=True)
                run=subprocess.run(command,cwd=work,text=True,capture_output=True)
                log=run.stdout+run.stderr
                (work/'run.txt').write_text(' '.join(command)+'\n'+log)
                values=dict(re.findall(r'(\w+)=([^\s]+)',log))
                row={k:values[k] for k in FIELDS if k in values}
                row.update(scan=args.scan,config=label,method='rosenbrock6dae',
                           options=json.dumps(options,sort_keys=True),scenario=case['name'],cells=128,
                           repeat=0,threads=4,status='ok' if run.returncode==0 else f'exit-{run.returncode}',
                           utc=datetime.now(timezone.utc).isoformat(),
                           binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                           compiler='clang++ -std=c++20 -O3 -DNDEBUG',
                           run_log=str((work/'run.txt').relative_to(ROOT)))
                if run.returncode:
                    writer.writerow(row); f.flush(); raise SystemExit(log)
                rows.append(row)
            for suffix in ['final.csv','snapshots.csv']:
                (directory/('reference_'+suffix)).write_bytes((directory/'tighter'/('ts1_'+suffix)).read_bytes())
            for label,row in zip(['tight','tighter'],rows):
                score=subprocess.run([sys.executable,str(SCENARIOS/'score.py'),
                    str(directory/label/'ts1_snapshots.csv'),str(directory/'reference_snapshots.csv'),'--json'],
                    capture_output=True,text=True,check=True)
                report=json.loads(score.stdout); worst=report['worst_species'][0]
                row.update(overall_digits=report['overall_digits'],worst_species=worst[0],worst_digits=worst[1],
                           gate=report['gate'],score_json=score.stdout.strip())
                writer.writerow(row); f.flush()
                if label=='tight':
                    assert report['overall_digits']>7 and worst[1]>6, report
                    print('reference convergence',report,flush=True)
            manifest.append(case | {'cells':128,'sha256':{name:hashlib.sha256((directory/name).read_bytes()).hexdigest()
                for name in ['scenario.txt','cells.csv','photolysis.csv','reference_snapshots.csv','reference_final.csv']}})
            (ROOT/'results/stress_inputs.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    main()
