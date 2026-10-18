#!/usr/bin/env python3
"""Generate and test independent 48-hour conditions; all artifacts stay local.

Run after selecting operating points. Reference C++ runner must be built from
../scenarios/micm_scenario.cpp first. No reference source is copied or changed.
"""
import argparse
import csv
import json
import hashlib
from pathlib import Path
import re
import shutil
import subprocess
import sys
from scan import ROOT, SCENARIOS, FIELDS

CASES={
 'stress-summer':['--seed','951','--dt','1800','--steps','96','--lat','60','--doy','172','--nox-scale','0.02','--snapshot-every','24'],
 'stress-winter':['--seed','952','--dt','900','--steps','192','--lat','60','--doy','355','--nox-scale','50','--snapshot-every','48'],
 'stress-equator':['--seed','953','--dt','1200','--steps','144','--lat','0','--doy','260','--nox-scale','5','--snapshot-every','36']}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--references',action='store_true')
    p.add_argument('--method',default='rodas3')
    p.add_argument('--option',action='append')
    p.add_argument('--point',help='Use a configuration and target from operating_points.json')
    p.add_argument('--level',type=float,default=3)
    p.add_argument('--repeat',type=int,default=0)
    a=p.parse_args()
    if a.point:
        point=next(p for p in json.loads((ROOT/'operating_points.json').read_text()) if p['id']==a.point)
        a.method=point['method'];a.level=point['level']
        a.option=[f'{key}={value}' for key,value in point['options'].items()]
    if not a.option: a.option=['rtol=0.0003','atol=1e-18']
    stamp=(ROOT/'build/configuration.txt').read_text().splitlines()
    with (ROOT/'bin/ts1_scenario').open('rb') as binary: digest=hashlib.file_digest(binary,'sha256').hexdigest()
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    failed=False
    work=ROOT/'results/stress-work';work.mkdir(parents=True,exist_ok=True)
    output=ROOT/'scans/stress.csv';exists=output.exists()
    with output.open('a',newline='') as file:
        writer=csv.DictWriter(file,FIELDS,extrasaction='ignore')
        if not exists: writer.writeheader()
        for name,args in CASES.items():
            case=ROOT/'results'/name;case.mkdir(parents=True,exist_ok=True)
            subprocess.run([sys.executable,str(SCENARIOS/'make_scenario.py'),str(case),
                '--cells','128','--diurnal',*args],check=True,capture_output=True)
            if a.references:
                # Two independently solved tolerances establish reference convergence.
                for rtol,atol,label in [('1e-9','1e-22','check'),('1e-10','1e-23','reference')]:
                    cmd=[str(ROOT/'bin/micm_scenario'),'--threads','10','--group-size','16','--method',
                         'rosenbrock6dae','--rtol',rtol,'--atol',atol,'--max-steps','100000',str(case/'scenario.txt')]
                    result=subprocess.run(cmd,cwd=work,capture_output=True,text=True,check=True)
                    (case/f'{label}.log').write_text(result.stdout+result.stderr)
                    for kind in ['final','snapshots']: shutil.copyfile(work/f'ts1_{kind}.csv',case/f'{label}_{kind}.csv')
                scored=subprocess.check_output([sys.executable,str(SCENARIOS/'score.py'),
                    str(case/'check_snapshots.csv'),str(case/'reference_snapshots.csv'),'--json'],text=True)
                print(name,'reference convergence',scored.strip(),flush=True)
                (case/'convergence.json').write_text(scored)
            if not (case/'reference_snapshots.csv').exists():
                raise RuntimeError('Generate references first with --references')
            cmd=[str(ROOT/'bin/ts1_scenario'),'--threads','10','--method',a.method]
            config={'method':a.method}
            for option in a.option:
                key,value=option.split('=',1);config[key]=float(value);cmd+=['--option',option]
            result=subprocess.run(cmd+[str(case/'scenario.txt')],cwd=work,capture_output=True,text=True)
            row=dict(re.findall(r'(\w+)=([^\s]+)',result.stdout))
            row.update(suite='stress',configuration=json.dumps(config,sort_keys=True,separators=(',',':')),
                scenario=name,repeat=a.repeat,compiler=stamp[0],
                threads=row.get('threads',10),returncode=result.returncode,message=result.stderr.strip().replace('\n',' | '),
                binary_sha256=digest,optflags=stamp[1],revision=revision)
            if result.returncode==0:
                score=json.loads(subprocess.run([sys.executable,str(SCENARIOS/'score.py'),
                    str(work/'ts1_snapshots.csv'),str(case/'reference_snapshots.csv'),'--json','--level',str(a.level)],capture_output=True,text=True).stdout)
                row.update(overall_digits=score['overall_digits'],worst_species=score['worst_species'][0][0],
                    worst_digits=score['worst_species'][0][1],gate=score['gate'],
                    level=score['level'],level_result=score['level_result'])
            else: row['gate']='FAIL'
            writer.writerow(row);file.flush()
            failed |= row.get('gate')!='PASS' or row.get('level_result')!='PASS'
            print(name,row.get('overall_digits','FAIL'),row.get('worst_digits',''),flush=True)
    if failed: raise SystemExit('Stress case failed its declared gate or accuracy target')
if __name__=='__main__': main()
