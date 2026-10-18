#!/usr/bin/env python3
"""Measure original C++ MICM scenario baseline with the same OPTFLAGS.

Group size 32 matches the C scan, and 100 matches the supplied reference runner's
1,000-cell ten-thread setup. This tool runs only the supplied reference program.
"""
import argparse
import csv
import json
import hashlib
import re
import subprocess
import sys
from scan import ROOT,SCENARIOS,NAMES,FIELDS

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',default='scans/micm.csv')
    p.add_argument('--group-size',type=int,default=32)
    p.add_argument('--repeat',type=int,default=1)
    p.add_argument('--rtol')
    p.add_argument('--atol')
    a=p.parse_args()
    work=ROOT/'results/micm-scan';work.mkdir(parents=True,exist_ok=True)
    output=ROOT/a.output;exists=output.exists()
    config={'method':'MICM-three-stage','group_size':a.group_size}
    with (ROOT/'bin/micm_scenario').open('rb') as binary: digest=hashlib.file_digest(binary,'sha256').hexdigest()
    if a.rtol: config['rtol']=float(a.rtol)
    if a.atol: config['atol']=float(a.atol)
    with output.open('a',newline='') as f:
        writer=csv.DictWriter(f,FIELDS,extrasaction='ignore')
        if not exists: writer.writeheader()
        for name in NAMES:
            for repeat in range(a.repeat):
                cmd=[str(ROOT/'bin/micm_scenario'),'--threads','10','--group-size',str(a.group_size)]
                if a.rtol: cmd+=['--rtol',a.rtol]
                if a.atol: cmd+=['--atol',a.atol]
                result=subprocess.run(cmd+[str(SCENARIOS/name/'scenario.txt')],cwd=work,capture_output=True,text=True)
                row=dict(re.findall(r'(\w+)=([^\s]+)',result.stdout))
                row.update(suite='micm-baseline',configuration=json.dumps(config,sort_keys=True,separators=(',',':')),
                    scenario=name,repeat=repeat,compiler='clang++',threads=10,returncode=result.returncode,
                    nonlinear_iterations=0,message=result.stderr.strip().replace('\n',' | '),
                    binary_sha256=digest,optflags=(ROOT/'build/configuration.txt').read_text().splitlines()[1],revision='MICM-e15898ce')
                if result.returncode==0:
                    kind='snapshots' if (SCENARIOS/name/'reference_snapshots.csv').exists() else 'final'
                    score=json.loads(subprocess.run([sys.executable,str(SCENARIOS/'score.py'),
                        str(work/f'ts1_{kind}.csv'),str(SCENARIOS/name/f'reference_{kind}.csv'),'--json','--level','3'],capture_output=True,text=True).stdout)
                    row.update(overall_digits=score['overall_digits'],worst_species=score['worst_species'][0][0],
                        worst_digits=score['worst_species'][0][1],gate=score['gate'],
                        level=score['level'],level_result=score['level_result'])
                else: row['gate']='FAIL'
                writer.writerow(row);f.flush()
                print(name,row.get('overall_digits','FAIL'),row.get('elapsed_ms',''),flush=True)
if __name__=='__main__': main()
