"""Configuration and benchmark regressions beyond the supplied case list."""
import csv
from pathlib import Path
import subprocess
from test_configs import verify_reader, robertson_equivalence, parameterization_equivalence

ROOT=Path(__file__).resolve().parent.parent


def run():
    verify_reader()
    robertson_equivalence()
    parameterization_equivalence()
    directory=ROOT/'build/ts1-regression';directory.mkdir(exist_ok=True)
    for group in [1]:
        subprocess.run([ROOT/'bin/ts1_bench_cpu','--threads','2','--group-size',str(group),'1','100','30',
            ROOT/'tests/data/ts1_stiff_cell.csv'],cwd=directory,check=True,capture_output=True,text=True)
        with (directory/'ts1_final.csv').open() as f,(ROOT/'tests/data/ts1_stiff_final.csv').open() as ref:
            actual,expected=csv.DictReader(f),csv.DictReader(ref)
            assert set(actual.fieldnames)==set(expected.fieldnames)
            rows=list(actual);refs=list(expected);assert len(rows)==len(refs)==1
            for name in rows[0]:
                got,want=float(rows[0][name]),float(refs[0][name])
                assert abs(got-want)<=1e-30+1e-10*abs(want),(group,name,got,want)

    for group in [1,7,32,128]:
        subprocess.run([ROOT/'bin/ts1_bench_cpu','--threads','2','--group-size',str(group),'100','100','30',
            ROOT.parent/'reference/ts1_cells_100_seed1.csv'],cwd=directory,check=True,capture_output=True,text=True)
        with (directory/'ts1_final.csv').open() as f,(ROOT.parent/'reference/micm_ts1_final_100_seed1.csv').open() as ref:
            actual,expected=csv.DictReader(f),csv.DictReader(ref)
            rows=list(actual);refs=list(expected);assert len(rows)==len(refs)==100
            assert set(actual.fieldnames)==set(expected.fieldnames)
            for cell,(row,wanted) in enumerate(zip(rows,refs)):
                for name in row:
                    got,want=float(row[name]),float(wanted[name])
                    assert abs(got-want)<=1e-30+1e-10*abs(want),(group,cell,name,got,want)
