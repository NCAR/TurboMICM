"""Independent cell solutions must not depend on thread scheduling."""
import os
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(dir=root / 'build') as directory:
    work = Path(directory)
    outputs = []
    for threads, schedule in [('1', 'static'), ('3', 'static,2'),
                              ('3', 'dynamic,1'), ('3', 'guided,2')]:
        environment = os.environ.copy()
        environment['OMP_SCHEDULE'] = schedule
        subprocess.run([root / 'bin/ts1_bench_cpu', '--threads', threads, '5', '2', '30',
                        root / '../reference/ts1_cells_100_seed1.csv'], cwd=work,
                       env=environment, check=True, capture_output=True)
        outputs.append((work / 'ts1_final.csv').read_bytes())
    assert all(output == outputs[0] for output in outputs[1:])
print('Benchmark schedules produce identical cell results')
