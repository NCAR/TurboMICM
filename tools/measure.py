"""Record native benchmark runs and optional full-state MICM comparison gates."""
import argparse
import hashlib
import json
import re
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('program', type=Path)
    parser.add_argument('cells_csv', type=Path)
    parser.add_argument('--label',required=True)
    parser.add_argument('--threads',type=int,default=18)
    parser.add_argument('--cells',type=int,default=10000)
    parser.add_argument('--steps',type=int,default=100)
    parser.add_argument('--dt',type=float,default=30)
    parser.add_argument('--repeats',type=int,default=1)
    parser.add_argument('--reference',type=Path)
    parser.add_argument('--micm',action='store_true')
    args = parser.parse_args()
    output = (ROOT/'results'/args.label).resolve()
    if not output.is_relative_to(ROOT/'results'):
        parser.error('label must resolve inside results/')
    output.mkdir(parents=True,exist_ok=True)
    program, cells_csv = args.program.resolve(), args.cells_csv.resolve()
    reference = args.reference.resolve() if args.reference else None
    command = [str(program),'--threads',str(args.threads),
               '--group-size' if args.micm else '--group','32',
               str(args.cells),str(args.steps),str(args.dt),str(cells_csv)]
    records = []
    for repeat in range(args.repeats):
        run_dir = output/f'run-{repeat+1}'
        run_dir.mkdir(exist_ok=True)
        run = subprocess.run(command,cwd=run_dir,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        (run_dir/'run.log').write_text(run.stdout)
        run.check_returncode()
        elapsed = float(re.search(r'elapsed_ms=([0-9.]+)',run.stdout).group(1))
        record = dict(elapsed_ms=elapsed,stdout=run.stdout,final_sha256=digest(run_dir/'ts1_final.csv'))
        if reference:
            gate = subprocess.run([sys.executable,str(ROOT.parent/'reference/compare_final.py'),
                                   str(run_dir/'ts1_final.csv'),str(reference)],
                                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            (run_dir/'gate.log').write_text(gate.stdout)
            gate.check_returncode()
            record['gate'] = gate.stdout.strip()
        records.append(record)
        summary = dict(command=command,input_sha256=digest(cells_csv),runs=records,
                       median_ms=statistics.median(r['elapsed_ms'] for r in records))
        if reference:
            summary['reference_sha256'] = digest(reference)
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(args.label,repeat+1,elapsed,record.get('gate',''),flush=True)


if __name__ == '__main__':
    main()
