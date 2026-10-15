"""Run the C tutorial ports and compare their actual printed reference results."""
import difflib
import re
import subprocess
from api import ROOT
from harness import case

PROGRAMS = ('README_example','solve_results','multiple_grid_cells','rate_constants_no_user_defined_example_by_hand',
            'rate_constants_user_defined_example_by_hand','solver_configuration','vectorized_matrix_solver','openmp_tutorial')


def normalize(text,name):
    # Elapsed time is measured by both programs, and is necessarily variable.
    text = re.sub(r'Total solve time: [0-9.e+]+ nanoseconds','Total solve time: <elapsed> nanoseconds',text)
    if name == 'openmp_tutorial':
        # Concurrent startup diagnostics in the original interleave individual
        # characters; the three result tables begin at the following marker.
        assert text.count('Running solver on thread') == 3
        text = text[text.index('Thread 1\n'):]
    return text


@case(*PROGRAMS)
def program(identifier):
    process = subprocess.run([str(ROOT/'bin/tutorials'),identifier],capture_output=True,text=True)
    assert process.returncode == 0,process.stderr
    expected = (ROOT/'tests/reference-output'/f'{identifier}.txt').read_text()
    actual,expected = normalize(process.stdout,identifier),normalize(expected,identifier)
    assert actual == expected,''.join(difflib.unified_diff(expected.splitlines(True),actual.splitlines(True)))
    if identifier == 'README_example':
        process = subprocess.run([str(ROOT/'bin/tutorials'),'configuration'],capture_output=True,text=True)
        assert process.returncode == 0,process.stderr
    return 'C program output matches reference; generated configuration checks passed'
