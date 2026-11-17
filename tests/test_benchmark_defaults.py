"""Check omitted condition columns against an independently generated reference."""
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(dir=root / 'build') as directory:
    work = Path(directory)
    inputs = work / 'cells.csv'
    inputs.write_text('OH,jo2_a,jh2o_b,jh2o_c,jch4_a,jch4_b,jhf\n1e-10,0,0,0,0,0,0\n')
    subprocess.run([root / 'bin/ts1_bench_cpu', '--threads', '1', '1', '2', '30', inputs],
                   cwd=work, check=True, capture_output=True)
    subprocess.run([root / '.conda/bin/python', root / '../reference/compare_final.py',
                    work / 'ts1_final.csv', root / 'tests/reference/ts1_missing_conditions.csv'],
                   check=True, capture_output=True)
    # With all photolysis defaults retained, this is an extremely sensitive
    # unphysical input. The relaxed bound accommodates compiler roundoff.
    inputs.write_text('OH\n1e-10\n')
    subprocess.run([root / 'bin/ts1_bench_cpu', '--threads', '1', '1', '0', '30', inputs],
                   cwd=work, check=True, capture_output=True)
    subprocess.run([root / '.conda/bin/python', root / '../reference/compare_final.py',
                    work / 'ts1_final.csv', root / 'tests/reference/ts1_stiff_defaults.csv',
                    '--rtol', '1e-3'], check=True, capture_output=True)
print('Benchmark default conditions match the reference')
