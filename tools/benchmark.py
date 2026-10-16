"""Record sequential benchmark runs and optionally apply the supplied TS1 gate."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,default=Path('bin/ts1_bench_cpu'))
p.add_argument('--input',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--expected',type=Path)
p.add_argument('--reference',type=Path,default=Path('../reference'))
p.add_argument('--cells',type=int,default=10000)
p.add_argument('--steps',type=int,default=100)
p.add_argument('--dt',type=float,default=30)
p.add_argument('--threads',type=int)
p.add_argument('--group-size',type=int)
p.add_argument('--schedule',choices=['static','dynamic','guided'])
p.add_argument('--repeats',type=int,default=5)
a=p.parse_args()
root=Path.cwd().resolve(); output=a.output.resolve()
if not output.is_relative_to(root): p.error('Output must be inside the checkout')
if a.repeats<1: p.error('Repeats must be positive')
command=[str(a.binary.resolve())]
for key in ['threads','group_size','schedule']:
    value=getattr(a,key)
    if value is not None: command += ['--'+key.replace('_','-'),str(value)]
command += [str(a.cells),str(a.steps),str(a.dt),str(a.input.resolve())]
runs=[]
for index in range(a.repeats):
    work=output/str(index+1);work.mkdir(parents=True,exist_ok=True)
    result=subprocess.run(command,cwd=work,capture_output=True,text=True,check=True)
    (work/'run.txt').write_text(result.stdout+result.stderr)
    metrics={k:float(v) if k=='elapsed_ms' else int(v)
             for k,v in (word.split('=',1) for word in result.stdout.split() if '=' in word)}
    if 'elapsed_ms' not in metrics: raise ValueError('Benchmark did not report elapsed_ms')
    if a.expected:
        gate=subprocess.run([sys.executable,str(a.reference.resolve()/'compare_final.py'),str(work/'ts1_final.csv'),str(a.expected.resolve())],
                            capture_output=True,text=True,check=True)
        (work/'gate.txt').write_text(gate.stdout+gate.stderr)
        metrics['gate']=gate.stdout.strip()
    runs.append(metrics);print(index+1,metrics,flush=True)
summary={'command':command,'runs':runs,'median_ms':statistics.median(r['elapsed_ms'] for r in runs)}
(output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('median_ms='+str(summary['median_ms']),flush=True)
