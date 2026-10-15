"""Run explicit C ports, preserving every reference case ID and failure status."""
from pathlib import Path
import argparse
import json
import subprocess
from test_constraints import CASES as CONSTRAINT_CASES
from test_solvers import CASES as SOLVER_CASES
from test_configs import CASES as CONFIG_CASES
from test_arrays import CASES as ARRAY_CASES
from test_sparse_arrays import CASES as SPARSE_CASES
from test_grouped import CASES as GROUPED_CASES
from test_linear import CASES as LINEAR_CASES
from test_state import CASES as STATE_CASES
from test_system import CASES as SYSTEM_CASES
from test_rate_store import CASES as STORE_CASES
from test_jacobian import CASES as JACOBIAN_CASES
from test_process_set import CASES as PROCESS_CASES
from test_analytical import CASES as ANALYTICAL_CASES
from test_terminator import CASES as TERMINATOR_CASES
from test_dae_integration import CASES as DAE_CASES
import supplemental
PYTHON_CASES = {**CONSTRAINT_CASES, **SOLVER_CASES, **CONFIG_CASES, **ARRAY_CASES, **SPARSE_CASES, **GROUPED_CASES, **LINEAR_CASES, **STATE_CASES, **SYSTEM_CASES, **STORE_CASES, **JACOBIAN_CASES, **PROCESS_CASES, **ANALYTICAL_CASES, **TERMINATOR_CASES, **DAE_CASES}

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--filter', default='')
    args = parser.parse_args()
    binary = ROOT / 'bin/test_rates'
    registered = set(subprocess.check_output([binary, '--list'], text=True).splitlines())
    counts = {'PASS': 0, 'FAIL': 0, 'N/A': 0}
    not_applicable = json.loads((ROOT / 'tests/not_applicable.json').read_text())
    for line in (ROOT / 'tests/cases.txt').read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        case, kind = line.split()
        if args.filter not in case:
            continue
        name = case.split(':', 1)[-1]
        key = case if case in PYTHON_CASES else name
        if case in not_applicable:
            status, reason = 'N/A', not_applicable[case]
        elif key in PYTHON_CASES:
            try:
                PYTHON_CASES[key]()
                status, reason = 'PASS', 'reference assertions'
            except Exception as error:
                status, reason = 'FAIL', f'{type(error).__name__}: {error}'
        elif name in registered:
            result = subprocess.run([binary, name], text=True, capture_output=True)
            status = 'PASS' if result.returncode == 0 else 'FAIL'
            reason = 'reference assertions' if status == 'PASS' else ' '.join(result.stderr.split())
        else:
            status, reason = 'FAIL', 'case port pending'
        counts[status] += 1
        print(f'{case} {status} {reason}', flush=True)
    try:
        supplemental.run()
        print('supplemental PASS configuration equivalence and TS1 clamp regression', flush=True)
    except Exception as error:
        counts['FAIL'] += 1
        print(f'supplemental FAIL {type(error).__name__}: {error}', flush=True)
    print(' '.join(f'{key}={value}' for key, value in counts.items()), flush=True)
    return int(counts['FAIL'] > 0)


if __name__ == '__main__':
    raise SystemExit(main())
