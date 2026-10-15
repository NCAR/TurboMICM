"""Plain-program stdout comparisons against the supplied double build."""
from functools import partial
from pathlib import Path
import re
import json
import subprocess
ROOT=Path(__file__).resolve().parent.parent


def normalize(text,name):
    if name=='solver_configuration':
        # Wall-clock durations are observations, not reproducible answers.
        text=re.sub(r'Total solve time: [0-9.e+]+ nanoseconds','Total solve time: <elapsed> nanoseconds',text)
    if name=='openmp_tutorial':
        # Unsynchronized reference thread log messages can interleave characters.
        text=text[text.index('Thread 1\n'):]
    return text


def program(name):
    result=subprocess.run([ROOT/'bin/cpu_tutorials',name],text=True,capture_output=True,timeout=45)
    assert result.returncode==0,(result.returncode,result.stderr)
    outputs=json.loads((ROOT/'tests/data/program_outputs.json').read_text())
    expected=''.join(outputs[name])
    actual=normalize(result.stdout,name);expected=normalize(expected,name)
    if actual!=expected:
        import difflib
        raise AssertionError('\n'.join(difflib.unified_diff(expected.splitlines(),actual.splitlines(),n=2)))


def openmp():
    result=subprocess.run([ROOT/'bin/cpu_tutorials','openmp_check'],capture_output=True,text=True,timeout=45)
    assert result.returncode==0,(result.returncode,result.stderr)


CASES={name:partial(program,name) for name in ['README_example','solve_results','multiple_grid_cells','solver_configuration','vectorized_matrix_solver','openmp_tutorial']}
CASES['openmp:OpenMP.OneSolverManyStates']=openmp
