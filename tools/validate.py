"""Rebuild and run the full Tier 1 plus new-method gates for one compiler."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run(command,directory=ROOT):
    result = subprocess.run(command,cwd=directory,capture_output=True,text=True)
    if result.returncode:
        raise RuntimeError(f'{command[0]} failed ({result.returncode}):\n{result.stdout}\n{result.stderr}')
    return result.stdout+result.stderr


def save(path,text):
    path.write_text(text.replace(str(ROOT),'<repository>'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',required=True,choices=['clang','gcc-15'])
    args = parser.parse_args()
    label = 'clang' if args.compiler == 'clang' else 'gcc'
    output = run(['make','cpu','bench','test-cpu',f'CC={args.compiler}'])
    save(ROOT/'validation'/f'tier3-{label}-cases.txt',output)
    expected = [line.split()[0] for line in (ROOT.parent/'reference/micm_cases.txt').read_text().splitlines()
                if line and not line.startswith('#')]
    found = [line.split()[0] for line in output.splitlines()
             if len(line.split()) > 1 and line.split()[0] in expected and line.split()[1] == 'PASS']
    assert sorted(found) == sorted(expected) and len(found) == 680
    method_output = run([sys.executable,'tests/test_integrator.py',f'validation/prototype-{label}.csv'])
    save(ROOT/'validation'/f'methods-{label}.txt',method_output)
    binary = ROOT/'bin/ts1_bench_cpu'
    result = dict(compiler=run([args.compiler,'--version']).splitlines()[0],cases_passed=len(found),
                  binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),gates={})
    large_reference = ROOT/'runs/gates/micm-large/ts1_final.csv'
    if not large_reference.exists():
        run(['make','micm-ref'])
        directory = large_reference.parent; directory.mkdir(parents=True,exist_ok=True)
        text = run([str(ROOT/'bin/micm_ts1_ref'),'--threads','10','10000','100','30',
                    str(ROOT/'runs/large/tier1-30s/cells.csv')],directory)
        save(directory/'stdout.txt',text)
    for name,n,threads,cells,reference in [
        ('small',100,1,ROOT.parent/'reference/ts1_cells_100_seed1.csv',
         ROOT.parent/'reference/micm_ts1_final_100_seed1.csv'),
        ('large',10000,10,ROOT/'runs/large/tier1-30s/cells.csv',large_reference)]:
        directory = ROOT/'runs/gates'/label/name; directory.mkdir(parents=True,exist_ok=True)
        text = run([str(binary),'--threads',str(threads),str(n),'100','30',str(cells)],directory)
        save(directory/'stdout.txt',text)
        comparison = run([sys.executable,str(ROOT.parent/'reference/compare_final.py'),
                          str(directory/'ts1_final.csv'),str(reference),'--rtol','1e-10','--atol','1e-30'])
        result['gates'][name] = dict(cells=n,threads=threads,stdout=text.strip(),comparison=comparison.strip(),
            input_sha256=hashlib.sha256(cells.read_bytes()).hexdigest(),
            reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
            result_sha256=hashlib.sha256((directory/'ts1_final.csv').read_bytes()).hexdigest())
        print(label,name,comparison.strip(),flush=True)
    path = ROOT/'validation/tier3-gates.json'
    report = json.loads(path.read_text()) if path.exists() else {}
    report[label] = result
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS: {label}: 680 named cases, new methods, small and 10,000-cell Tier 1 gates',flush=True)


if __name__ == '__main__':
    main()
