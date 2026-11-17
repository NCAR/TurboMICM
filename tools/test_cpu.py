#!/usr/bin/env python3
"""Report each reference case independently; unported cases are failures."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
results = {}
for executable in ('test_cases', 'test_matrix_cases', 'test_matrix_function_cases', 'test_matrix_remaining', 'test_sparse_cases', 'test_sparse_function_cases', 'test_sparse_remaining'):
    run = subprocess.run([str(root / 'build' / executable)], text=True, capture_output=True, check=True)
    for line in run.stdout.splitlines():
        case, passed = line.split()
        if case in results:
            raise ValueError(f'duplicate case result: {case}')
        results[case] = passed == '1'
expected = (root / 'tests/tutorial_expected.txt').read_text()
for case, executable, args in [
    ('rate_constants_no_user_defined_example_by_hand', 'tutorial_no_user', []),
    ('rate_constants_user_defined_example_by_hand', 'tutorial_user', ['user']),
]:
    result = subprocess.run([str(root / 'build' / executable), *args], text=True, capture_output=True)
    results[case] = result.returncode == 0 and result.stdout == expected
not_applicable = {
    'view_category_sfinae:ViewCategorySFINAE.VectorTypes': 'C has no C++ concepts or SFINAE',
    'view_category_sfinae:ViewCategorySFINAE.MatrixViews': 'C has no template view-category traits',
    'view_category_sfinae:ViewCategorySFINAE.SparseMatrixViews': 'C has no template view-category traits',
    'view_category_sfinae:ViewCategorySFINAE.MatrixTypesExcluded': 'C has no C++ concepts or SFINAE',
}
failed = False
for line in (root / 'tests/micm_cases.txt').read_text().splitlines():
    if not line or line.startswith('#'):
        continue
    case, kind = line.split()
    if case in not_applicable:
        status, reason = 'N/A', not_applicable.pop(case)
    elif case in results:
        ok = results.pop(case)
        status, reason = ('PASS', 'ported assertions') if ok else ('FAIL', 'ported assertion failed')
        failed |= not ok
    else:
        status, reason = 'FAIL', 'case not yet ported'
        failed = True
    print(f'{case} {status} {reason}')
if results or not_applicable:
    raise ValueError(f'unknown case results: {list(results)}, {list(not_applicable)}')
sys.exit(int(failed))
