#!/usr/bin/env python3
"""Large-cell grouping/scheduling study, scoring the disclosed first 1,000 cells.

All 10,000 cells are checked for the finite/nonnegative gate by score.py.
Only the disclosed first 1,000 cells have high-accuracy references; do not
misrepresent their score as a 10,000-cell accuracy measurement.
"""
import csv
import hashlib
import itertools
import json
import subprocess
import sys
from scan import ROOT,SCENARIOS,FIELDS,run

def main():
    case=ROOT/'results/scale-diurnal-10min';case.mkdir(parents=True,exist_ok=True)
    subprocess.run([sys.executable,str(SCENARIOS/'make_scenario.py'),str(case),
        '--cells','10000','--seed','101','--dt','600','--steps','144','--diurnal',
        '--snapshot-every','36'],check=True,capture_output=True)
    with (case/'cells.csv').open() as generated, (SCENARIOS/'diurnal-10min/cells.csv').open() as supplied:
        a,b=csv.reader(generated),csv.reader(supplied)
        for row in b: assert row==next(a),'Generated prefix differs from disclosed cells'
    work=ROOT/'results/scale-work';work.mkdir(parents=True,exist_ok=True)
    path=ROOT/'scans/scale.csv';exists=path.exists()
    stamp=(ROOT/'build/configuration.txt').read_text().splitlines()
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    with (ROOT/'bin/ts1_scenario').open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
    with path.open('a',newline='') as f:
        writer=csv.DictWriter(f,FIELDS,extrasaction='ignore')
        if not exists: writer.writeheader()
        for group,schedule in itertools.product([16,32,64],[0,1,2]):
            options=dict(rtol=3e-4,atol=1e-18,group_size=group,schedule=schedule,schedule_chunk=1 if schedule else 0)
            row=run('rodas31',options,'scale-diurnal-10min',10,work,3,case,SCENARIOS/'diurnal-10min')
            row.update(suite='scale',configuration=json.dumps({'method':'rodas31',**options},sort_keys=True,separators=(',',':')),
                scenario='scale-diurnal-10min',repeat=0,compiler=stamp[0],threads=10,
                binary_sha256=digest,optflags=stamp[1],revision=revision,
                message='Accuracy: disclosed first 1000 cells; finite/nonnegative gate: all 10000 cells')
            writer.writerow(row);f.flush()
            print(options,row.get('overall_digits','FAIL'),row.get('worst_digits',''),row.get('elapsed_ms',''),flush=True)
if __name__=='__main__': main()
