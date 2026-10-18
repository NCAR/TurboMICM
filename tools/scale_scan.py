#!/usr/bin/env python3
"""Large-cell grouping/scheduling study, scoring the disclosed first 1,000 cells.

All 10,000 cells are checked for the finite/nonnegative gate by score.py.
Only the disclosed first 1,000 cells have high-accuracy references; do not
misrepresent their score as a 10,000-cell accuracy measurement.
"""
import csv
import argparse
import hashlib
import itertools
import json
import subprocess
import sys
from scan import ROOT,SCENARIOS,FIELDS,run

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points',action='store_true',help='Validate all declared points at 10000 cells')
    args=parser.parse_args()
    case=ROOT/'results/scale-diurnal-10min';case.mkdir(parents=True,exist_ok=True)
    subprocess.run([sys.executable,str(SCENARIOS/'make_scenario.py'),str(case),
        '--cells','10000','--seed','101','--dt','600','--steps','144','--diurnal',
        '--snapshot-every','36'],check=True,capture_output=True)
    with (case/'cells.csv').open() as generated, (SCENARIOS/'diurnal-10min/cells.csv').open() as supplied:
        a,b=csv.reader(generated),csv.reader(supplied)
        for row in b: assert row==next(a),'Generated prefix differs from disclosed cells'
    work=ROOT/'results/scale-work';work.mkdir(parents=True,exist_ok=True)
    path=ROOT/'scans'/('scale-points.csv' if args.points else 'scale.csv');exists=path.exists()
    stamp=(ROOT/'build/configuration.txt').read_text().splitlines()
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    with (ROOT/'bin/ts1_scenario').open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
    with path.open('a',newline='') as f:
        writer=csv.DictWriter(f,FIELDS+(['point'] if args.points else []),extrasaction='ignore')
        if not exists: writer.writeheader()
        jobs=json.loads((ROOT/'operating_points.json').read_text()) if args.points else [
            dict(id='',method='rodas31',level=3,options=dict(rtol=3e-4,atol=1e-18,group_size=group,schedule=schedule,schedule_chunk=1 if schedule else 0))
            for group,schedule in itertools.product([16,32,64],[0,1,2])]
        failed=False
        for point in jobs:
            options=point['options']
            row=run(point['method'],options,'scale-diurnal-10min',10,work,point['level'],case,SCENARIOS/'diurnal-10min')
            row.update(suite='scale',point=point['id'],configuration=json.dumps({'method':point['method'],**options},sort_keys=True,separators=(',',':')),
                scenario='scale-diurnal-10min',repeat=0,compiler=stamp[0],threads=10,
                binary_sha256=digest,optflags=stamp[1],revision=revision,
                message='Accuracy: disclosed first 1000 cells; finite/nonnegative gate: all 10000 cells')
            writer.writerow(row);f.flush()
            failed |= row.get('gate')!='PASS' or row.get('level_result')!='PASS'
            print(options,row.get('overall_digits','FAIL'),row.get('worst_digits',''),row.get('elapsed_ms',''),flush=True)
    if failed: raise SystemExit('Large-cell gate or accuracy check failed')
if __name__=='__main__': main()
