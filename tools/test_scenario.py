#!/usr/bin/env python3
"""Driver contract checks on small copies of disclosed data, inside workspace."""
import csv
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
SC=ROOT.parent/'scenarios'

def main():
    with tempfile.TemporaryDirectory(dir=ROOT/'results',prefix='driver-test-') as tmp:
        work=Path(tmp)
        with (SC/'diurnal-10min/cells.csv').open() as f:
            rows=list(csv.reader(f))[:4]
        with (work/'cells.csv').open('w') as f: csv.writer(f).writerows(rows)
        (work/'photo.csv').write_text('time_s,factor\n0,0\n10,0.2\n20,0.1\n')
        (work/'scenario.txt').write_text('cells = cells.csv\ndt = 10\nsteps = 5\nphotolysis = photo.csv\nsnapshot_every = 2\n')
        cmd=[str(ROOT/'bin/ts1_scenario'),'--threads','2','--method','ros4','--option','group_size=2']
        output=subprocess.check_output(cmd+[str(work/'scenario.txt')],cwd=work,text=True)
        assert 'cells=3 groups=2 group_size=2 threads=2' in output
        with (work/'ts1_snapshots.csv').open() as f: snapshots=list(csv.DictReader(f))
        with (work/'ts1_final.csv').open() as f: final=list(csv.DictReader(f))
        assert [(int(r['step']),int(r['cell'])) for r in snapshots]==[(s,c) for s in [2,4,5] for c in range(3)]
        assert len(final)==3
        for snap,last in zip(snapshots[-3:],final):
            assert all(snap[k]==v for k,v in last.items())
        reference=ROOT/'bin/micm_scenario'
        if reference.exists():
            # Independent C++ driver checks factor interpolation, constant
            # extension, rate updates, initial defaults and snapshot semantics.
            ccmd=[str(ROOT/'bin/ts1_scenario'),'--threads','2','--method','micm3',
                  '--option','group_size=2','--option','rtol=1e-4','--option','atol=1e-18']
            subprocess.run(ccmd+[str(work/'scenario.txt')],cwd=work,check=True,capture_output=True)
            with (work/'ts1_snapshots.csv').open() as f: actual=list(csv.DictReader(f))
            subprocess.run([str(reference),'--threads','2','--group-size','2',
                '--rtol','1e-4','--atol','1e-18',str(work/'scenario.txt')],
                cwd=work,check=True,capture_output=True)
            with (work/'ts1_snapshots.csv').open() as f: expected=list(csv.DictReader(f))
            assert len(actual)==len(expected)
            for x,z in zip(actual,expected):
                for key,value in z.items():
                    a,b=float(x[key]),float(value)
                    assert abs(a-b)<=1e-30+1e-10*abs(b),(key,a,b)
            print('Independent MICM scenario driver parity PASS (all cells/species/snapshots)')
        # Input failures are explicit, not a silent return to default options.
        for arg in ['group_size=0','jacobian_steps=0','rtol=nan','start=-1','unknown=3','norm=3']:
            assert subprocess.run(cmd+['--option',arg,str(work/'scenario.txt')],cwd=work,capture_output=True).returncode==2,arg
        bad=work/'bad.txt'
        bad.write_text((work/'scenario.txt').read_text()+'dt = 4\n')
        assert subprocess.run(cmd+[str(bad)],cwd=work,capture_output=True).returncode==2
        print('Scenario parser, photolysis input, final partial group, snapshot remainder and invalid options PASS')
if __name__=='__main__': main()
