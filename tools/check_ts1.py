"""Run the C TS1 executable and the supplied cell-by-cell gate in isolated directories.

Example: .conda/bin/python tools/check_ts1.py --binary bin/ts1_bench_cpu
All subprocess working directories remain inside this checkout.
"""
import argparse
from pathlib import Path
import subprocess
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,default=Path('bin/ts1_bench_cpu'))
p.add_argument('--reference',type=Path,default=Path('../reference'))
p.add_argument('--input',type=Path)
p.add_argument('--expected',type=Path)
p.add_argument('--cells',type=int,default=100)
p.add_argument('--steps',type=int,default=100)
p.add_argument('--dt',type=float,default=30)
p.add_argument('--output',type=Path,default=Path('build/ts1-checks'))
p.add_argument('--quick',action='store_true',help='Only the default all-core run')
a=p.parse_args()
root=Path.cwd().resolve(); output=a.output.resolve()
if not output.is_relative_to(root): p.error('Output must be inside the checkout')
binary=a.binary.resolve(); reference=a.reference.resolve()
source=(a.input or reference/'ts1_cells_100_seed1.csv').resolve()
expected=(a.expected or reference/'micm_ts1_final_100_seed1.csv').resolve()
variants=[('default',[])]
if not a.quick:
    variants += [(f't{threads}-g{group}',['--threads',str(threads),'--group-size',str(group)])
                 for threads,group in [(1,1),(4,7),(1,32),(3,64),(4,a.cells+1)]]
    variants += [('dynamic',['--threads','4','--group-size','7','--schedule','dynamic','--chunk','2']),
                 ('guided',['--threads','3','--group-size','31','--schedule','guided','--chunk','1'])]
    variants += [('generic',['--threads','2','--group-size','13','--generic'])]
for name,options in variants:
    work=output/name; work.mkdir(parents=True,exist_ok=True)
    command=[str(binary),*options,str(a.cells),str(a.steps),str(a.dt),str(source)]
    result=subprocess.run(command,cwd=work,capture_output=True,text=True,check=True)
    (work/'run.txt').write_text(result.stdout+result.stderr)
    gate=subprocess.run([sys.executable,str(reference/'compare_final.py'),str(work/'ts1_final.csv'),str(expected)],
                        capture_output=True,text=True,check=True)
    (work/'gate.txt').write_text(gate.stdout+gate.stderr)
    print(name,result.stdout.strip(),gate.stdout.strip(),sep='\n',flush=True)
