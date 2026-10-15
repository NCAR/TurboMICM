"""Run the native benchmark and the supplied all-cell/all-species accuracy gate.

Run `make cpu CC=...` first. Outputs and run metadata are kept in the selected
directory, so measurements from different compilers and settings cannot collide.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('cells',type=int)
parser.add_argument('steps',type=int)
parser.add_argument('dt',type=float)
parser.add_argument('input',type=Path)
parser.add_argument('reference',type=Path)
parser.add_argument('output',type=Path)
parser.add_argument('--threads',type=int)
parser.add_argument('--group',type=int,default=32)
parser.add_argument('--schedule',default='static')
parser.add_argument('--compiler',required=True)
args = parser.parse_args()
out = args.output.resolve(); out.mkdir(parents=True,exist_ok=True)
command = [str(ROOT/'bin/ts1_bench_cpu'),'--group',str(args.group)]
if args.threads is not None: command += ['--threads',str(args.threads)]
command += [str(args.cells),str(args.steps),str(args.dt),str(args.input.resolve())]
environment = dict(os.environ,OMP_SCHEDULE=args.schedule)
result = subprocess.run(command,cwd=out,env=environment,capture_output=True,text=True)
(out/'run.log').write_text(result.stdout+result.stderr)
print(result.stdout,end='',flush=True)
if result.returncode: raise SystemExit(result.stderr)
gate = subprocess.run([str(ROOT/'.conda/bin/python'),str(ROOT.parent/'reference/compare_final.py'),str(out/'ts1_final.csv'),str(args.reference.resolve())],capture_output=True,text=True)
(out/'gate.log').write_text(gate.stdout+gate.stderr)
print(gate.stdout,end='',flush=True)
metadata = {'compiler':args.compiler,'command':command,'schedule':args.schedule,'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
            'reference_sha256':hashlib.sha256(args.reference.read_bytes()).hexdigest(),'gate_exit_code':gate.returncode,'stdout':result.stdout,'gate':gate.stdout}
(out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
if gate.returncode: raise SystemExit(gate.stderr or gate.stdout)
