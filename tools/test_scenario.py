"""Check forcing interpolation, snapshot indexing, CSV paths, and generic kernels.

The supplied MICM C++ runner is an independent oracle for the compatibility
path. All generated inputs and outputs are confined to build/scenario-tests.
"""
import csv
import math
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
work = root / 'build' / 'scenario-tests'
work.mkdir(parents=True, exist_ok=True)
source = root.parent / 'scenarios' / 'tier1-30s' / 'cells.csv'
with source.open() as f, (work/'cells with spaces.csv').open('w', newline='') as out:
    reader, writer = csv.reader(f), csv.writer(out)
    for _ in range(4):
        writer.writerow(next(reader))
(work/'photo.csv').write_text('time_s,factor\n-10,0.1\n5,0.3\n13,0.9\n')
(work/'scenario.txt').write_text('''# inline comments and paths containing spaces
cells = cells with spaces.csv
dt = 7 # seconds
steps = 5
photolysis = photo.csv
snapshot_every = 2
''')


def run(label, binary, options, scenario_name='scenario.txt'):
    dest = work/label; dest.mkdir(exist_ok=True)
    result = subprocess.run([str(root/'bin'/binary), *options, str(work/scenario_name)],
                            cwd=dest, text=True, capture_output=True, check=True)
    (dest/'run.txt').write_text(result.stdout+result.stderr)
    return dest


def compare(a, b, name):
    with (a/name).open() as f:
        rows_a = list(csv.DictReader(f))
    with (b/name).open() as f:
        rows_b = list(csv.DictReader(f))
    assert len(rows_a) == len(rows_b)
    if name == 'ts1_snapshots.csv':
        assert [(int(r['step']), int(r['cell'])) for r in rows_a] == [(s,c) for s in [2,4,5] for c in range(3)]
    largest = 0.
    for row_a, row_b in zip(rows_a, rows_b):
        assert row_a.keys() == row_b.keys()
        for key in row_a:
            x,y = float(row_a[key]),float(row_b[key])
            assert math.isfinite(x) and math.isfinite(y)
            assert abs(x-y) <= 1e-30+1e-10*abs(y), (name,key,x,y)
            if y:
                largest=max(largest,abs(x-y)/abs(y))
    return largest


c = run('c-reference', 'ts1_scenario', ['--threads','2','--method','micm4dae','--option','group-size=3'])
ref = run('micm-reference', 'micm_scenario', ['--threads','2','--method','rosenbrock4dae','--group-size','3'])
checks=[]
for filename in ['ts1_final.csv','ts1_snapshots.csv']:
    checks.append(('MICM '+filename, compare(c,ref,filename)))
for method in ['ros2','ros3','ros4','rodas3','rodas4','extrap','sdirk']:
    options=['--threads','2','--method',method,'--option','group-size=2']
    native=run(method+'-generated','ts1_scenario',options)
    generic=run(method+'-generic','ts1_scenario',options+['--option','generic=1'])
    uncached=run(method+'-uncached','ts1_scenario',options+['--option','rate-cache=0'])
    blocked=run(method+'-blocked','ts1_scenario',options+['--option','step-block=99','--option','schedule=dynamic'])
    for filename in ['ts1_final.csv','ts1_snapshots.csv']:
        checks.append((method+' '+filename,compare(native,generic,filename)))
        assert (native/filename).read_bytes()==(uncached/filename).read_bytes(), (method,filename,'rate cache mismatch')
        checks.append((method+' cached/uncached '+filename,0.))
        assert (native/filename).read_bytes()==(blocked/filename).read_bytes(), (method,filename,'blocking mismatch')
        checks.append((method+' blocked/dynamic '+filename,0.))
    (work/'scenario-final.txt').write_text((work/'scenario.txt').read_text().replace('snapshot_every = 2\n',''))
    final_only=run(method+'-final-only','ts1_scenario',options+['--option','step-block=3','--option','schedule=guided'], 'scenario-final.txt')
    assert not (final_only/'ts1_snapshots.csv').exists()
    assert (native/'ts1_final.csv').read_bytes()==(final_only/'ts1_final.csv').read_bytes()
with (root/'results/scenario_io.csv').open('w',newline='') as f:
    writer=csv.writer(f); writer.writerow(['check','max_relative_difference']); writer.writerows(checks)
print('Scenario forcing/snapshots agree with MICM; seven methods agree with generic kernels and bit-identical cached/uncached rates')
