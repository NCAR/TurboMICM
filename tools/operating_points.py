#!/usr/bin/env python3
"""Run all declared operating points, checking their declared scorer levels.

The manifest is the single source for the commands in OPERATING_POINTS.md.
Timings are always from C. Runs are sequential to avoid CPU contention.
"""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
from scan import ROOT, FIELDS, NAMES, run

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repeat',type=int,default=5)
    p.add_argument('--output',type=Path,default=ROOT/'scans/operating-points.csv')
    p.add_argument('--point',action='append')
    p.add_argument('--scenario',action='append',choices=NAMES)
    p.add_argument('--threads',type=int,default=10)
    a=p.parse_args()
    points=json.loads((ROOT/'operating_points.json').read_text())
    if a.point: points=[point for point in points if point['id'] in a.point]
    assert points and a.repeat>0 and a.threads>0
    stamp=(ROOT/'build/configuration.txt').read_text().splitlines()
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    with (ROOT/'bin/ts1_scenario').open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
    work=ROOT/'results/operating-points';work.mkdir(parents=True,exist_ok=True)
    exists=a.output.exists();failed=False
    with a.output.open('a',newline='') as f:
        writer=csv.DictWriter(f,FIELDS+['point'],extrasaction='ignore')
        if not exists: writer.writeheader()
        for point in points:
            for name,repeat in itertools.product(a.scenario or NAMES,range(a.repeat)):
                row=run(point['method'],point['options'],name,a.threads,work,point['level'])
                row.update(suite='operating-points',point=point['id'],
                    configuration=json.dumps({'method':point['method'],**point['options']},sort_keys=True,separators=(',',':')),
                    scenario=name,repeat=repeat,compiler=stamp[0],threads=a.threads,
                    binary_sha256=digest,optflags=stamp[1],revision=revision)
                writer.writerow(row);f.flush()
                good=row.get('level_result')=='PASS' and row.get('gate')=='PASS'
                failed |= not good
                print(point['id'],name,repeat,row.get('overall_digits','FAIL'),row.get('worst_digits',''),row.get('elapsed_ms',''),flush=True)
    if failed: raise SystemExit('At least one declared point missed its accuracy/gate target')
if __name__=='__main__': main()
