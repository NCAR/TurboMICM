#!/usr/bin/env python3
"""Reproducible C-only timing / accuracy scan. One durable row per run.

Examples: .conda/bin/python tools/scan.py --output scans/family.csv --suite family
All timing comes from bin/ts1_scenario; accuracy comes from the supplied score.py.
"""
import argparse
import csv
import itertools
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ROOT.parent / 'scenarios'
NAMES = ['tier1-30s', 'diurnal-30min', 'diurnal-10min', 'diurnal-5min', 'polluted-10min', 'clean-10min']
FIELDS = ['suite', 'configuration', 'scenario', 'repeat', 'compiler', 'threads', 'returncode',
          'overall_digits', 'worst_species', 'worst_digits', 'gate', 'elapsed_ms', 'internal_steps',
          'accepted', 'rejected', 'function_calls', 'jacobian_updates', 'decompositions', 'solves',
          'nonlinear_iterations', 'cells', 'groups', 'message', 'level', 'level_result',
          'binary_sha256', 'optflags', 'revision']

def configurations(suite):
    if suite == 'final-tune':
        return [('rodas31', {'atol':1e-18,**o}) for o in [
            {'rtol':.002}, {'rtol':.0015}, {'rtol':.0007},
            {'rtol':.0003,'atol':1e-17}, {'rtol':.0003,'atol':1e-16},
            {'rtol':.0015,'atol':1e-16},
            {'rtol':.0005,'controller':1}, {'rtol':.0005,'controller':2},
            {'rtol':.001,'controller':1}, {'rtol':.001,'controller':2}]]
    if suite == 'baseline':
        return [('micm3', {}), ('micm3', {'rtol': 1e-4, 'atol': 1e-18})]
    if suite == 'family':
        return [(m, {'rtol': t, 'atol': 1e-18, 'reuse_rhs': 0}) for m, t in itertools.product(
            ['ros2', 'ros3', 'ros4', 'rodas3', 'rodas4', 'extrap2', 'extrap3'], [1e-3, 1e-4, 1e-5])]
    if suite == 'loose':
        return [(m, {'rtol': t, 'atol': 1e-18}) for m,t in itertools.product(
            ['ros3','ros4','rodas3','rodas4'], [0.003,0.01,0.03])]
    if suite == 'rodas31':
        return [('rodas31', {'rtol': t, 'atol': 1e-18}) for t in [1e-2,3e-3,1e-3,5e-4,3e-4,1e-4]]
    if suite == 'h211b':
        return [('rodas3', dict(rtol=3e-4,atol=1e-18,**o)) for o in [
            {}, {'controller':1}, {'controller':2},
            {'controller':2,'controller_b':2}, {'controller':2,'controller_b':1},
            {'controller':2,'controller_b':1,'controller_k':2},
            {'controller':2,'controller_b':1,'controller_k':1.7},
            {'safety':1.1}, {'safety':1.3}, {'safety':1.5}, {'growth':12}, {'growth':30}]]
    if suite == 'mid':
        return [(m, {'rtol': t, 'atol': 1e-18}) for m,t in itertools.product(
            ['rodas3','rodas4'], [7e-4,5e-4,3e-4])]
    if suite == 'absolute':
        return [('rodas3', {'rtol': 3e-4, 'atol': t}) for t in [1e-14,1e-16,1e-18,1e-20]]
    if suite == 'group':
        return [('rodas3', {'rtol': 3e-4, 'atol': 1e-18, 'group_size': w}) for w in [4,8,16,32,64,128]]
    if suite == 'controller':
        return [('rodas3', dict(rtol=3e-4, atol=1e-18, **o)) for o in [
            {}, {'controller': 1}, {'norm': 0}, {'norm': 2}, {'warm_start': 0},
            {'start': 1e-4}, {'start': 1}, {'positivity': 0},
            {'safety': .7}, {'safety': .95}, {'growth': 2}, {'growth': 3}, {'reuse_rhs': 0}]]
    if suite == 'reuse':
        return [('rodas3', dict(rtol=3e-4, atol=1e-18, jacobian_steps=j, deadband=d))
                for j,d in [(1,0),(2,0),(4,0),(8,0),(2,0.2),(4,0.2),(8,0.2),(4,0.5)]]
    raise ValueError(suite)

def run(method, options, name, threads, work, level=3, scenario_dir=None, reference_dir=None):
    scenario_dir=Path(scenario_dir) if scenario_dir else SCENARIOS/name
    reference_dir=Path(reference_dir) if reference_dir else scenario_dir
    command = [str(ROOT/'bin/ts1_scenario'), '--threads', str(threads), '--method', method]
    for key,value in options.items(): command += ['--option', f'{key}={value}']
    command += [str(scenario_dir/'scenario.txt')]
    result = subprocess.run(command, cwd=work, capture_output=True, text=True)
    data = dict(re.findall(r'(\w+)=([^\s]+)', result.stdout))
    data.update(returncode=result.returncode, message=result.stderr.strip().replace('\n',' | '))
    if result.returncode == 0:
        filename = 'snapshots' if (reference_dir/'reference_snapshots.csv').exists() else 'final'
        scored = subprocess.run([sys.executable, str(SCENARIOS/'score.py'),
            str(work/f'ts1_{filename}.csv'), str(reference_dir/f'reference_{filename}.csv'),
            '--json', '--level', str(level)],
            capture_output=True, text=True)
        try:
            score = json.loads(scored.stdout)
            data.update(overall_digits=score['overall_digits'], worst_species=score['worst_species'][0][0],
                worst_digits=score['worst_species'][0][1], gate=score['gate'],
                level=score['level'], level_result=score['level_result'])
        except (ValueError, KeyError, IndexError):
            data.update(gate='FAIL', message='Scoring failed: '+scored.stdout.strip()+' '+scored.stderr.strip())
    else: data['gate']='FAIL'
    return data

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--suite',required=True)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--scenario',action='append',choices=NAMES)
    p.add_argument('--threads',type=int,default=10)
    p.add_argument('--repeat',type=int,default=1)
    p.add_argument('--level',type=float,default=3)
    p.add_argument('--method')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--option',action='append',default=[])
    a=p.parse_args()
    configs = [(a.method,dict(v.split('=',1) for v in a.option))] if a.method else configurations(a.suite)
    work=ROOT/'results'/('scan-'+a.output.stem)
    work.mkdir(parents=True,exist_ok=True)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    stamp=(ROOT/'build/configuration.txt').read_text().splitlines()
    compiler=stamp[0]
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    with (ROOT/'bin/ts1_scenario').open('rb') as binary:
        digest=hashlib.file_digest(binary,'sha256').hexdigest()
    exists=a.output.exists()
    done=set()
    if a.resume and exists:
        with a.output.open() as previous:
            done={(r['configuration'],r['scenario'],int(r['repeat'])) for r in csv.DictReader(previous)}
    with a.output.open('a',newline='') as file:
        writer=csv.DictWriter(file,FIELDS,extrasaction='ignore')
        if not exists: writer.writeheader()
        for method,options in configs:
            config=json.dumps({'method':method,**options},sort_keys=True,separators=(',',':'))
            for name,repeat in itertools.product(a.scenario or NAMES,range(a.repeat)):
                if (config,name,repeat) in done: continue
                row=run(method,options,name,a.threads,work,a.level)
                row.update(suite=a.suite,configuration=config,scenario=name,repeat=repeat,compiler=compiler,
                    threads=a.threads,binary_sha256=digest,optflags=stamp[1],revision=revision)
                writer.writerow(row); file.flush()
                print(config,name,row.get('overall_digits','FAIL'),row.get('worst_digits',''),row.get('elapsed_ms',''),flush=True)
if __name__=='__main__': main()
