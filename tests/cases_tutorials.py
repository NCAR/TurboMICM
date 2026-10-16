"""Native tutorial output is checked against the copied reference executables."""
import re
import subprocess
from pathlib import Path
from functools import lru_cache
from harness import case

ROOT = Path(__file__).resolve().parents[1]
IDS = '''README_example solve_results multiple_grid_cells rate_constants_no_user_defined_example_by_hand
rate_constants_user_defined_example_by_hand solver_configuration vectorized_matrix_solver openmp_tutorial'''


@lru_cache(maxsize=1)
def configuration_parity():
    result = subprocess.run([str(ROOT/'bin/tutorial_cpu'),'config_parity'],capture_output=True,text=True)
    assert result.returncode == 0, result.stdout+result.stderr


@case(*IDS.split())
def tutorial(id_):
    result = subprocess.run([str(ROOT/'bin/tutorial_cpu'),id_],capture_output=True,text=True)
    assert result.returncode == 0, result.stdout+result.stderr
    expected = (ROOT/'tests/reference'/f'{id_}.txt').read_text()
    def normalize(text):
        lines = text.splitlines()
        # Wall-clock measurements and thread start ordering are nondeterministic.
        timing = [line for line in lines if line.startswith('Total solve time:')]
        if id_ == 'solver_configuration':
            assert len(timing) == 5
            for line in timing:
                assert re.fullmatch(r'Total solve time: [0-9.e+]+ nanoseconds',line)
        return '\n'.join(line.rstrip() for line in lines if not line.startswith(('Total solve time:','Running solver on thread')))
    actual,reference = normalize(result.stdout),normalize(expected)
    assert actual == reference, '\n'.join(f'line {i}: {a!r} != {b!r}' for i,(a,b) in enumerate(zip(actual.splitlines(),reference.splitlines())) if a != b)[:1000]
    if id_ == 'openmp_tutorial':
        assert sorted(re.findall(r'Running solver on thread (\d)',result.stdout)) == ['0','1','2']
    if id_ in ('multiple_grid_cells','rate_constants_no_user_defined_example_by_hand','rate_constants_user_defined_example_by_hand'):
        configuration_parity()
